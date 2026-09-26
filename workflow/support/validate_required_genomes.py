#!/usr/bin/env python3
"""Compatibility entrypoint for the formatted-genome requirement."""

from validate_required_species_outputs import main

if __name__ == "__main__":
    main(default_required=("genome",))
