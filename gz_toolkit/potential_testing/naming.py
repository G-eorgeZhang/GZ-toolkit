"""Canonical potential naming: ``<name>_<system>_<year>`` (e.g. ``gao_FeHe_2011``).

The name part may be anything (``tao_2.31``); the year is optional and, when
missing, the name ends with a bare underscore (``tao_2.31_FeHeCr_``).

Users may start with a short name such as ``gao2011`` or ``gao_2011``. Once the
``pot_inputs/<pot_name>.json`` exists (so the element system is known),
``normalize_pot_names`` renames, consistently:

* ``pot_inputs/<old>.json``  -> ``pot_inputs/<new>.json`` (``potential.pot_name`` updated)
* ``potentials/<old>/``      -> ``potentials/<new>/``
* ``<root>/<old>/`` (run dir, if already built) -> ``<root>/<new>/``

Names already in the canonical form are left alone.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

_YEAR = r"(?:19|20)\d{2}"
# Already canonical: ends in _<element symbols>_ with an optional year, e.g. _FeHe_2011 / _FeHeCr_
_CANONICAL = re.compile(rf"_(?:[A-Z][a-z]?)+_(?:{_YEAR})?$")
# Short name with a trailing year: gao2011, gao_2011, tao_2.31_2026
_WITH_YEAR = re.compile(rf"^(.*?[A-Za-z0-9.])[_-]?({_YEAR})$")


def system_string(potential: dict) -> str:
    """Element system like 'FeCrHe', in type order (type_map), else metals then gases."""
    type_map = potential.get("type_map") or {}
    if type_map:
        elements = [el for el, _ in sorted(type_map.items(), key=lambda kv: kv[1])]
    else:
        elements = list(potential.get("single_elements", [])) + list(potential.get("gases", []))
    return "".join(elements)


def canonical_pot_name(pot_name: str, system: str) -> str:
    """Return ``<name>_<system>_<year>``; ``<name>_<system>_`` when no year was given.

    Names already in that form are returned unchanged.
    """
    if not system or _CANONICAL.search(pot_name):
        return pot_name
    m = _WITH_YEAR.match(pot_name)
    if m:
        return f"{m.group(1)}_{system}_{m.group(2)}"
    return f"{pot_name}_{system}_"


def _rename_dir(src: Path, dst: Path, log: list[str]) -> None:
    if not src.is_dir():
        return
    if dst.exists():
        log.append(f"  skipped {src} -> {dst}: target already exists, merge by hand")
        return
    src.rename(dst)


def normalize_pot_names(pot_inputs_dir: str | Path, root: str | Path | None = None,
                        print_fn=print) -> list[tuple[str, str]]:
    """Rename every non-canonical potential under ``pot_inputs_dir``. Returns [(old, new)]."""
    pot_inputs = Path(pot_inputs_dir)
    root = Path(root) if root is not None else pot_inputs.resolve().parent
    renamed: list[tuple[str, str]] = []
    warnings: list[str] = []
    if not pot_inputs.is_dir():
        return renamed

    for f in sorted(pot_inputs.glob("*.json")):
        text = f.read_text(encoding="utf-8")
        try:
            potential = json.loads(text).get("potential", {})
        except (json.JSONDecodeError, AttributeError):
            continue  # loader reports bad JSON with a better message
        old = potential.get("pot_name", f.stem)
        new = canonical_pot_name(old, system_string(potential))
        if new == old:
            continue
        new_json = pot_inputs / f"{new}.json"
        if new_json.exists():
            warnings.append(f"  skipped {old}: {new_json.name} already exists")
            continue

        new_text = re.sub(r'("pot_name"\s*:\s*)"' + re.escape(old) + '"',
                          lambda m: f'{m.group(1)}"{new}"', text, count=1)
        new_json.write_text(new_text, encoding="utf-8")
        if f != new_json:
            f.unlink()
        _rename_dir(root / "potentials" / old, root / "potentials" / new, warnings)
        _rename_dir(root / old, root / new, warnings)
        renamed.append((old, new))

    for old, new in renamed:
        print_fn(f"[potential_testing] renamed potential '{old}' -> '{new}' "
                 f"(json, potentials/ dir and run dir)")
    for w in warnings:
        print_fn(f"[potential_testing] {w.strip()}")
    return renamed
