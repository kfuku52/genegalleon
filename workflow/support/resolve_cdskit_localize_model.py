#!/usr/bin/env python3
"""Select the newest published CDSKIT localization checkpoint once per cache."""

import argparse

from cdskit.latest_model import DOWNLOAD_PREFIX as DOWNLOAD_PREFIX
from cdskit.latest_model import RELEASES_URL as RELEASES_URL
from cdskit.latest_model import latest_model as _latest_model
from cdskit.latest_model import open_url as open_url
from cdskit.latest_model import resolve as _resolve
from cdskit.latest_model import verify as verify


def latest_model():
    return _latest_model(opener=open_url)


def resolve(cache_dir, allow_download=True):
    return _resolve(cache_dir, allow_download, namespace="genegalleon-latest", selector=latest_model, opener=open_url)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--model-download", choices=("yes", "no"), default="yes")
    args = parser.parse_args()
    print(resolve(args.cache_dir, args.model_download == "yes"))


if __name__ == "__main__":
    main()
