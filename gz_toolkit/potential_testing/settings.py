"""Project-wide settings: ``<root>/settings.json`` (or ``settings.yaml``).

Example::

    {
      "potinfo_savepath": true
    }

``potinfo_savepath`` controls where tested-potential info is saved once a
potential's tests finish (``summarize-pot``):

* ``false``          -> do not save.
* ``true``           -> save into ``<root>/pot_infobank/``.
* ``"<path>"``       -> save into ``<path>/pot_infobank/`` (created if missing).
  Relative paths are taken relative to ``<root>``; ``~`` is expanded.

A missing settings file means ``potinfo_savepath = false`` (older projects).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SETTINGS_FILENAMES = ("settings.json", "settings.yaml", "settings.yml")
INFOBANK_DIRNAME = "pot_infobank"
DEFAULT_SETTINGS: dict[str, Any] = {"potinfo_savepath": True}


def settings_path(root: str | Path = ".") -> Path | None:
    for name in SETTINGS_FILENAMES:
        p = Path(root) / name
        if p.is_file():
            return p
    return None


def load_project_settings(root: str | Path = ".") -> dict[str, Any]:
    """Read the root settings file; ``{}`` if there is none."""
    p = settings_path(root)
    if p is None:
        return {}
    text = p.read_text(encoding="utf-8")
    if p.suffix == ".json":
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in '{p}': {exc.msg} at line {exc.lineno}") from exc
    else:
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover - optional dep
            raise RuntimeError(f"'{p}' requires PyYAML; install with `pip install pyyaml`.") from exc
        data = yaml.safe_load(text) or {}
    if not isinstance(data, dict):
        raise ValueError(f"'{p}' must contain a mapping (e.g. {{\"potinfo_savepath\": true}})")
    return data


def write_default_settings(root: str | Path = ".") -> Path | None:
    """Create ``<root>/settings.json`` with defaults unless a settings file exists."""
    if settings_path(root) is not None:
        return None
    p = Path(root) / "settings.json"
    p.write_text(json.dumps(DEFAULT_SETTINGS, indent=2) + "\n", encoding="utf-8")
    return p


def resolve_infobank_dir(root: str | Path = ".", settings: dict[str, Any] | None = None) -> Path | None:
    """Directory to save potential info into, or ``None`` when saving is off."""
    root = Path(root).resolve()
    if settings is None:
        settings = load_project_settings(root)
    value = settings.get("potinfo_savepath", False)
    if value is False or value is None:
        return None
    if value is True:
        return root / INFOBANK_DIRNAME
    if isinstance(value, str):
        if value.strip().lower() in ("false", "no", "none", ""):
            return None
        if value.strip().lower() in ("true", "yes"):
            return root / INFOBANK_DIRNAME
        base = Path(value).expanduser()
        if not base.is_absolute():
            base = root / base
        return base / INFOBANK_DIRNAME
    raise ValueError(f"potinfo_savepath must be true, false or a path string, got {value!r}")
