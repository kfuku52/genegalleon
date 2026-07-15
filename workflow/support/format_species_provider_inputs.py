from pathlib import Path

from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS, PROVIDERS


def manifest_declared_providers(rows, provider_filter="all"):
    normalized_filter = str(provider_filter or "all").strip().lower()
    if normalized_filter != "all":
        return [normalized_filter]

    providers = []
    seen = set()
    for row in rows:
        provider = str((row or {}).get("provider", "") or "").strip().lower()
        if provider == "" or provider not in PROVIDERS or provider in seen:
            continue
        seen.add(provider)
        providers.append(provider)
    return providers


def manifest_declared_species_keys(rows, provider):
    normalized_provider = str(provider or "").strip().lower()
    return {
        str((row or {}).get("species_key", "") or "").strip()
        for row in rows
        if str((row or {}).get("provider", "") or "").strip().lower() == normalized_provider
        and str((row or {}).get("species_key", "") or "").strip() != ""
    }


def resolve_provider_inputs(args, manifest_rows=None):
    if args.provider == "all":
        if args.input_dir == "":
            raise ValueError("--input-dir is required when --provider all is used.")
        input_root = Path(args.input_dir).expanduser().resolve()
        provider_names = list(PROVIDERS)
        if manifest_rows is not None:
            provider_names = manifest_declared_providers(manifest_rows, args.provider)
        return [
            (provider, (input_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]).resolve())
            for provider in provider_names
        ]

    provider = args.provider
    if args.input_dir != "":
        return [(provider, Path(args.input_dir).expanduser().resolve())]
    raise ValueError("Specify --input-dir.")
