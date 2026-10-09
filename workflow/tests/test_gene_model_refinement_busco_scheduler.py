"""Controlled scheduler/cache fixtures, not a scientific BUSCO run."""
import hashlib
import json
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from importlib import import_module
from pathlib import Path
from unittest import mock

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
MODULE = import_module("gene_model_refinement_busco")

class Batch:
    def read(self, paths):
        return {str(path): hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in paths}

    def check(self):
        pass


class RecordingPool(ThreadPoolExecutor):
    instances = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.submitted = []
        self.drain_entered = threading.Event()
        self.guard = threading.Lock()
        self.inflight = 0
        self.maximum_inflight = 0
        self.instances.append(self)

    def submit(self, fn, pair):
        with self.guard:
            self.submitted.append(pair["species"])
            self.inflight += 1
            self.maximum_inflight = max(self.maximum_inflight, self.inflight)
        future = super().submit(fn, pair)

        def complete(_):
            with self.guard:
                self.inflight -= 1
        future.add_done_callback(complete)
        return future

    def shutdown(self, *args, **kwargs):
        self.drain_entered.set()
        return super().shutdown(*args, **kwargs)


class ConcurrencyFixtures(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.lineage = self.root / "lineage"
        self.lineage.mkdir()
        (self.lineage / "dataset.cfg").write_text("toy\n")
        self.pairs = []
        for index in range(5):
            path = self.root / f"species{index}.fa"
            path.write_text(">toy\nATGTAA\n")
            self.pairs.append({"species": f"species{index}", "before": str(path), "after": str(path)})
        self.stack = []
        overrides = {
            "FreshDigestBatch": Batch,
            "digest": lambda path: hashlib.sha256(str(Path(path).name).encode()).hexdigest(),
            "staged_result": lambda pair, before, after, pre: {
                "species": pair["species"], "before_result": {"complete": 100},
                "after_result": {"complete": 100}, "delta_complete": 0, "delta_complete_pp": 0,
                "delta_duplicated": 0, "refinement_status": "analysed", "reason": ""},
            "plot_comparison": lambda *args: None,
            "ThreadPoolExecutor": RecordingPool,
        }
        for name, value in overrides.items():
            patch = mock.patch.object(MODULE, name, value)
            self.stack.append(patch)
            patch.start()
        self.stack.extend([mock.patch.object(MODULE.subprocess, "check_output", return_value="BUSCO fixture"),
                           mock.patch.object(MODULE.shutil, "which", side_effect=lambda name: "/" + name)])
        for patch in self.stack[len(overrides):]:
            patch.start()
        RecordingPool.instances.clear()

    def tearDown(self):
        for patch in reversed(self.stack):
            patch.stop()
        self.tmp.cleanup()

    def launch(self, run_one, pairs, jobs, name="report"):
        outcome = {}
        patch = mock.patch.object(MODULE, "run_one", run_one)
        patch.start()

        def work():
            try:
                outcome["rows"] = MODULE.evaluate(pairs, self.root / name, self.lineage, self.root, cpus=4, jobs=jobs)
            except BaseException as error:
                outcome["error"] = error
        thread = threading.Thread(target=work)
        thread.start()
        return thread, outcome, patch

    def test_first_observed_failure_stops_unsubmitted_and_naturally_drains(self):
        for fail_index in (0, 1):
            with self.subTest(failure_input_index=fail_index):
                entered = [threading.Event(), threading.Event()]
                release_failure = threading.Event()
                release_survivor = threading.Event()
                calls = []
                lock = threading.Lock()

                def run_one(pair, phase, report, contract, cpus, *, lock=lock, calls=calls,
                            entered=entered, fail_index=fail_index,
                            release_failure=release_failure, release_survivor=release_survivor):
                    index = int(pair["species"][-1])
                    with lock:
                        calls.append((index, phase))
                    if phase == "before":
                        entered[index].set()
                        if index == fail_index:
                            if not release_failure.wait(5):
                                raise AssertionError("Failure release timed out")
                            raise RuntimeError("deliberate fixture failure")
                        if not release_survivor.wait(5):
                            raise AssertionError("Survivor release timed out")
                    return {}
                thread, outcome, patch = self.launch(run_one, self.pairs[:3], 2, f"failure{fail_index}")
                try:
                    self.assertTrue(entered[0].wait(5))
                    self.assertTrue(entered[1].wait(5))
                    pool = RecordingPool.instances[-1]
                    release_failure.set()
                    self.assertTrue(pool.drain_entered.wait(5))
                    self.assertTrue(thread.is_alive(), "Must naturally wait for admitted survivor")
                    self.assertEqual(pool.submitted, ["species0", "species1"])
                    self.assertLessEqual(pool.maximum_inflight, 2)
                    release_survivor.set()
                    thread.join(5)
                    self.assertFalse(thread.is_alive())
                    self.assertIsInstance(outcome.get("error"), RuntimeError)
                    self.assertEqual(pool.submitted, ["species0", "species1"])
                    self.assertIn((1 - fail_index, "after"), calls, "Running whole-species job drains without abort")
                    self.assertNotIn((fail_index, "after"), calls)
                    self.assertFalse((self.root / f"failure{fail_index}/busco_comparison.json").exists())
                finally:
                    release_failure.set()
                    release_survivor.set()
                    thread.join(5)
                    patch.stop()

    def test_dynamic_refill_preserves_species_order_and_jobs_bound(self):
        entered = [threading.Event() for _ in self.pairs]
        releases = [threading.Event() for _ in self.pairs]

        def run_one(pair, phase, report, contract, cpus):
            index = int(pair["species"][-1])
            if phase == "before":
                entered[index].set()
                if not releases[index].wait(5):
                    raise AssertionError("Ordered fixture release timed out")
            return {}
        thread, outcome, patch = self.launch(run_one, self.pairs, 3)
        try:
            for event in entered[:3]:
                self.assertTrue(event.wait(5))
            self.assertFalse(entered[3].is_set())
            releases[2].set()
            self.assertTrue(entered[3].wait(5), "Freed slot must refill despite held first species")
            releases[3].set()
            self.assertTrue(entered[4].wait(5))
            for release in releases:
                release.set()
            thread.join(5)
            self.assertFalse(thread.is_alive())
            self.assertNotIn("error", outcome)
            self.assertEqual([row["species"] for row in outcome["rows"]], [p["species"] for p in self.pairs])
            pool = RecordingPool.instances[-1]
            self.assertLessEqual(pool.maximum_inflight, 3)
            self.assertEqual(pool.submitted, [p["species"] for p in self.pairs])
        finally:
            for release in releases:
                release.set()
            thread.join(5)
            patch.stop()

    def test_jobs_one_and_three_have_same_complete_contract_and_order(self):
        for jobs in (1, 3):
            with mock.patch.object(MODULE, "run_one", return_value={}):
                rows = MODULE.evaluate(self.pairs, self.root / f"jobs{jobs}", self.lineage, self.root, cpus=4, jobs=jobs)
                self.assertEqual([row["species"] for row in rows], [p["species"] for p in self.pairs])
        one = json.loads((self.root / "jobs1/contract.json").read_text())
        three = json.loads((self.root / "jobs3/contract.json").read_text())
        self.assertEqual(one, three)
        self.assertEqual(one["contract"]["cpus_per_job"], 4)
        self.assertNotIn("jobs", one["contract"])

    def test_zero_jobs_still_rejects_without_normal_comparison(self):
        with mock.patch.object(MODULE, "run_one", return_value={}):
            with self.assertRaises(ValueError):
                MODULE.evaluate(self.pairs, self.root / "jobs0", self.lineage, self.root, jobs=0)
        self.assertFalse((self.root / "jobs0/busco_comparison.json").exists())



if __name__ == "__main__":
    unittest.main()
