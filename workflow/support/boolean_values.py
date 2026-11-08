"""Strict CLI boolean conversion and nullable annotation-table conversion."""

import pandas

try:
    from distutils.util import strtobool as strtobool
except ModuleNotFoundError:
    def strtobool(val):
        val = str(val).strip().lower()
        if val in ("y", "yes", "t", "true", "on", "1"):
            return 1
        if val in ("n", "no", "f", "false", "off", "0"):
            return 0
        raise ValueError("invalid truth value {!r}".format(val))


def parse_boolish(value):
    if pandas.isna(value):
        return pandas.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return pandas.NA
