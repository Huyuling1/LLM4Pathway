"""Frozen currency-molecule configuration, matched by reference identifier."""

from pathlib import Path

from .common import digest, read_json

DEFAULT_CONFIG = Path(__file__).parent / "config" / "currency.v1.json"


def load_currency(path=None):
    config = read_json(path or DEFAULT_CONFIG)
    if not config.get("version") or not isinstance(config.get("molecules"), list):
        raise ValueError("Currency configuration requires version and molecules")
    identifiers = {row["chebi_id"] for row in config["molecules"]}
    if len(identifiers) != len(config["molecules"]):
        raise ValueError("Repeated currency ChEBI identifier")
    return config, digest(config), identifiers, set(config.get("reactome_st_ids", []))
