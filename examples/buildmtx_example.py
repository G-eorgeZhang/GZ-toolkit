from gz_toolkit.buildmtx.buildstr import Gen_crystal

g = Gen_crystal(atom_style="atomic")

# 1) BCC Fe: seed → replicate to .data
g.seed_crystal("A2", ["Fe", "Fe"], "bcc_seed.data")
g.replicate("bcc_seed.data", 50, 50, 50, 2.86, "bcc_50x50x50.data")

# 2) FCC Cu: seed → replicate to POSCAR
g.seed_crystal("A1", ["Cu", "Cu", "Cu", "Cu"], "fcc_seed.data")
g.replicate("fcc_seed.data", 20, 20, 20, 3.615, "POSCAR_Cu")

# 3) B2 NiAl: seed → replicate to .vasp
g.seed_crystal("B2", ["Ni", "Al"], "NiAl_seed.data")
g.replicate("NiAl_seed.data", 30, 30, 30, 2.887, "NiAl_30x30x30.vasp")

# =====================================================
# trans_unit  —  rotate seed to a new orientation
# =====================================================

# 4) BCC Fe: reorient the cubic seed so that
#      x → [1,-1,0],  y → [1,1,1],  z → [1,1,-2]
#    This is the standard orientation for edge-dislocation studies.
g.seed_crystal("A2", ["Fe", "Fe"], "bcc_seed.data")
transformed = g.trans_unit(
    seed_file="bcc_seed.data",
    directions=[[1, -1, 0],
                [1,  1, 1],
                [1,  1, -2]],
    # filename is optional; auto-generates "bcc_110_111_112.data"
)
print(f"Transformed unit cell written to: {transformed}")

# 5) trans_unit + replicate  —  build a large oriented supercell
#    First create the rotated unit cell, then tile it.
g.replicate(transformed, 30, 40, 20, 2.87, "bcc_110_111_112_30x40x20.data")

# 6) FCC Cu: full pipeline  seed → trans_unit → replicate
g.seed_crystal("A1", ["Cu", "Cu", "Cu", "Cu"], "fcc_seed.data")
fcc_unit = g.trans_unit(
    seed_file="fcc_seed.data",
    directions=[[1, 1, 0],
                [0, 0, 1],
                [1, -1, 0]],
    filename="fcc_110_001_110.data",
)
g.replicate(fcc_unit, 25, 25, 25, 3.615, "fcc_110_001_110_25x25x25.data")

# =====================================================
# pot_infobank  —  resolve lc from a tested potential instead of hardcoding it
# =====================================================
#
# Once `gz-toolkit-potential-testing summarize-pot` has run for a potential
# (automatically, at the end of its job chain — see the potential-testing
# manual), its lc lives in <run_dir>/<pot_name>/<pot_name>.json, and gets
# copied into gz_toolkit/pot_infobank/ if that potential's JSON sets
# "promote_to_infobank". Here we fake up one such record (a tempdir standing
# in for a real pot_infobank/) to show the lookup pattern end-to-end.

import json
import tempfile
from pathlib import Path
from gz_toolkit.pot_infobank import load_tag, get_lc

infobank_dir = Path(tempfile.mkdtemp())
(infobank_dir / "bonny.json").write_text(json.dumps({
    "pot_name": "bonny",
    "type_map": {"Fe": 1, "Cr": 2},
    "masses": {"1": 55.845, "2": 51.996},
    "lc": {
        "Fe": 2.8553,
        "Cr": 2.8841,
        "alloys": {"FeCr": {"Fe95_Cr5": 2.859, "Fe90_Cr10": 2.861}},
    },
}))

pot = load_tag("bonny", infobank_dir=infobank_dir)

# 7) Pure Fe — look up lc by tag instead of hardcoding a number.
lc_fe = get_lc(pot, "Fe")
g.seed_crystal("A2", ["Fe", "Fe"], "bcc_seed.data")
g.replicate("bcc_seed.data", 50, 50, 50, lc_fe, "bcc_bonny_50x50x50.data")

# 8) Fe-9Cr — nearest tested composition to 9 at% Cr (here: Fe90_Cr10, 10%).
lc_fe9cr = get_lc(pot, "Fe", solute="Cr", at_pct=9)
g.replicate("bcc_seed.data", 50, 50, 50, lc_fe9cr, "bcc_bonny_Fe9Cr_50x50x50.data")

# 9) Overriding the looked-up value — there's no special override mechanism
#    in pot_infobank; an explicit lc always wins simply because you control
#    what gets passed to replicate() instead of calling get_lc() at all.
lc_override = 2.87  # e.g. from your own DFT reference instead of pot_infobank
g.replicate("bcc_seed.data", 50, 50, 50, lc_override, "bcc_bonny_override_50x50x50.data")
