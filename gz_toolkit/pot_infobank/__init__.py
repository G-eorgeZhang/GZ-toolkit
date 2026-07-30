"""Cross-run store of per-potential structure-building facts (lc, type_map, masses).

Once a potential's tests finish, ``gz_toolkit.potential_testing.summary.
build_pot_tag_json`` writes a minimal record to
``<run_dir>/<pot_name>/<pot_name>.json`` (e.g. ``gao2011/gao2011.json``) —
just ``pot_name``, ``type_map``, ``masses``, and ``lc`` (lattice constant per
element/alloy composition), the four things you need to go from
``pot="gao2011"`` to an actual cell, without re-parsing LAMMPS logs.
Everything else (Ecoh, elastic constants, defect formation energies, ...) is
deliberately left out of this record — see ``summary.csv`` for that.
``promote_tag_json`` optionally copies that same file (same filename) into a
shared ``pot_infobank/`` directory (keyed by
``PotentialMetadata.promote_to_infobank``), so later work — ``buildmtx``,
``defect``, and the ``analyze`` module — can pull a potential's lattice
constant by tag alone:

    pot = gz_toolkit.pot_infobank.load_tag("gao2011", infobank_dir="/path/to/pot_infobank")
    lc_fe = get_lc(pot, "Fe")
    lc_fe9cr = get_lc(pot, "Fe", solute="Cr", at_pct=9)  # nearest tested composition

Overriding a looked-up value (e.g. in ``buildmtx``) needs no mechanism here
— callers should just accept an explicit ``lc=`` argument that short-circuits
the ``get_lc``/pot_infobank lookup entirely when given, rather than reaching
into or mutating this module's state.
"""

from __future__ import annotations

import json
import re
import shutil
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gz_toolkit.potential_testing.config import PotentialConfig


def promote_tag_json(cfg: "PotentialConfig", run_dir: str | Path = ".") -> Path | None:
    """Copy ``<run_dir>/<pot_name>/<pot_name>.json`` into the shared pot_infobank, if requested.

    Controlled by ``cfg.potential.promote_to_infobank``:
      * ``False``/``None`` (default) -> no-op, returns ``None``.
      * ``True``  -> destination is resolved from the pot's cluster
        (``cfg.hpc.cluster``) via its ``path2gz_toolkit`` in
        ``gz_toolkit/jobs/cluster_info/<CLUSTER>/meta.json``, validated to
        actually point at a ``gz_toolkit`` install before appending
        ``pot_infobank/``.
      * a path string -> used verbatim as the destination directory. No
        validation — the caller is trusted.

    Raises if ``<pot_name>.json`` hasn't been built yet (run
    ``summary.build_pot_tag_json`` first) or if ``True`` is requested without
    a resolvable ``path2gz_toolkit``.
    """
    run_dir = Path(run_dir).resolve()
    tag_json_path = run_dir / cfg.potential.pot_name / f"{cfg.potential.pot_name}.json"
    promote = cfg.potential.promote_to_infobank
    if not promote:
        return None
    if not tag_json_path.is_file():
        raise FileNotFoundError(
            f"No {cfg.potential.pot_name}.json to promote at {tag_json_path} — "
            "run build_pot_tag_json() for this pot first."
        )

    if isinstance(promote, str):
        dest_dir = Path(promote)
    else:
        cluster = cfg.hpc.cluster
        if not cluster:
            raise ValueError(
                "promote_to_infobank=True requires hpc.cluster to be set "
                "(path2gz_toolkit is looked up per-cluster)."
            )
        from gz_toolkit.jobs.template_engine import load_cluster_meta
        path2gz_toolkit = load_cluster_meta(cluster).get("path2gz_toolkit")
        if not path2gz_toolkit:
            raise ValueError(
                f"Cluster '{cluster}' has no path2gz_toolkit set in "
                f"gz_toolkit/jobs/cluster_info/{cluster.upper()}/meta.json"
            )
        if Path(path2gz_toolkit).name != "gz_toolkit":
            raise ValueError(
                f"Invalid path2gz_toolkit for cluster '{cluster}': {path2gz_toolkit!r} "
                "does not point at a gz_toolkit install (expected the path to end in "
                "'gz_toolkit')."
            )
        dest_dir = Path(path2gz_toolkit) / "pot_infobank"

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{cfg.potential.pot_name}.json"
    shutil.copy2(tag_json_path, dest)
    return dest


@lru_cache(maxsize=None)
def _load_cached(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_tag(tag: str, infobank_dir: str | Path) -> dict:
    """Load one potential's standardized record by tag from ``infobank_dir``."""
    path = Path(infobank_dir) / f"{tag}.json"
    if not path.is_file():
        raise FileNotFoundError(f"No pot_infobank record for '{tag}' at {path}")
    return _load_cached(str(path.resolve()))


def _solute_fraction(comp_tag: str, solute: str) -> float | None:
    """Extract the solute's at% from a composition tag like ``"Fe95_Cr5"``.

    Returns ``None`` if ``solute`` doesn't appear in ``comp_tag`` (e.g. it's
    the majority-element part, or a different pair entirely).
    """
    for part in comp_tag.split("_"):
        m = re.match(rf"^{re.escape(solute)}(\d+(?:\.\d+)?)$", part)
        if m:
            return float(m.group(1))
    return None


def get_lc(
    pot: dict,
    element: str,
    solute: str | None = None,
    at_pct: float | None = None,
) -> float:
    """Look up a lattice constant from a ``load_tag``/``load_potential`` record.

        get_lc(pot, "Fe")                          -> pure-element lc
        get_lc(pot, "Fe", solute="Cr", at_pct=9)    -> alloy lc, nearest tested
                                                        Cr at% to 9 (exact match
                                                        if it was tested)

    Alloy pairs are stored under ``pot["lc"]["alloys"]`` keyed by the
    concatenated element symbols in whichever order the potential's
    ``two_element_suites`` used (e.g. ``"FeCr"``) — both orderings of
    ``element``/``solute`` are tried, since the caller may not know which.

    Raises ``KeyError`` if ``element`` (pure case) or the alloy pair (alloy
    case) isn't in the record at all — that means it was never tested, not
    that the composition just wasn't exact.
    """
    if solute is None:
        return pot["lc"][element]

    alloys = pot.get("lc", {}).get("alloys", {})
    pair = alloys.get(element + solute)
    if pair is None:
        pair = alloys.get(solute + element)
    if pair is None:
        raise KeyError(
            f"No alloy lc for {element}-{solute} in this pot record "
            f"(available alloy pairs: {sorted(alloys)})"
        )

    if at_pct is None:
        if len(pair) == 1:
            return next(iter(pair.values()))
        raise ValueError(
            f"Multiple {element}-{solute} compositions available "
            f"({sorted(pair)}) — pass at_pct to pick one."
        )

    best_tag, best_diff = None, None
    for comp_tag in pair:
        frac = _solute_fraction(comp_tag, solute)
        if frac is None:
            continue
        diff = abs(frac - at_pct)
        if best_diff is None or diff < best_diff:
            best_tag, best_diff = comp_tag, diff
    if best_tag is None:
        raise KeyError(
            f"Could not parse solute fraction from any composition tag for "
            f"{element}-{solute}: {sorted(pair)}"
        )
    return pair[best_tag]


__all__ = ["promote_tag_json", "load_tag", "get_lc"]
