"""Interactive project setup — plain-language questions, no JSON editing.

Designed for a materials scientist with no computing background: every
question shows its default in [brackets] (press Enter to accept), invalid
answers re-prompt with an explanation, and the wizard ends by validating
the config and printing exactly what to do next.

Fully testable: pass ``input_fn`` / ``print_fn`` to drive it from a script.
"""

from __future__ import annotations

import re
from pathlib import Path

from gz_toolkit.potential_testing.config import (
    HPCOptions,
    PotentialConfig,
    PotentialMetadata,
    WorkflowOptions,
    save_potential_config,
)


def suggest_mass(symbol: str) -> float | None:
    """Atomic mass (amu) from the periodic table; None for unknown symbols."""
    try:
        from pymatgen.core.periodic_table import Element
        return float(Element(symbol).atomic_mass)
    except Exception:  # noqa: BLE001 - unknown symbol / pymatgen quirk
        return None


# Rough lattice-constant starting guesses (Å) for common bcc/fcc metals —
# just a seed for the box relax, accuracy doesn't matter.
_LC_GUESS = {
    "Fe": 2.86, "Cr": 2.88, "W": 3.17, "Mo": 3.15, "V": 3.03, "Nb": 3.30,
    "Ta": 3.31, "Al": 4.05, "Cu": 3.61, "Ni": 3.52, "Au": 4.08, "Ag": 4.09,
    "Pt": 3.92, "Pb": 4.95, "Ti": 2.95, "Zr": 3.23, "Mg": 3.21,
}
_STRUCT_GUESS = {
    "Fe": "bcc", "Cr": "bcc", "W": "bcc", "Mo": "bcc", "V": "bcc",
    "Nb": "bcc", "Ta": "bcc",
    "Al": "fcc", "Cu": "fcc", "Ni": "fcc", "Au": "fcc", "Ag": "fcc",
    "Pt": "fcc", "Pb": "fcc",
    "Ti": "hcp", "Zr": "hcp", "Mg": "hcp",
}


def _ask(prompt, default=None, parse=str, validate=None, choices=None,
         input_fn=input, print_fn=print):
    """Ask one question; re-prompt until the answer parses and validates."""
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input_fn(f"{prompt}{suffix}: ").strip()
        if not raw:
            if default is not None:
                return default
            print_fn("  Please type an answer (there is no default for this one).")
            continue
        try:
            value = parse(raw)
        except (TypeError, ValueError):
            print_fn(f"  '{raw}' is not a valid value here — please try again.")
            continue
        if choices is not None and value not in choices:
            print_fn(f"  Please answer one of: {', '.join(str(c) for c in choices)}")
            continue
        if validate is not None:
            problem = validate(value)
            if problem:
                print_fn(f"  {problem}")
                continue
        return value


def _ask_yesno(prompt, default: bool, input_fn=input, print_fn=print) -> bool:
    d = "Y/n" if default else "y/N"
    val = _ask(prompt, default=d, parse=lambda s: s.lower(),
               input_fn=input_fn, print_fn=print_fn)
    if val in ("Y/n", "y/N"):
        return default
    if val in ("y", "yes"):
        return True
    if val in ("n", "no"):
        return False
    return default


def run_init_wizard(root: str | Path = ".", *, input_fn=input, print_fn=print) -> Path:
    """Interactive questionnaire; scaffolds the project and writes
    pot_inputs/<pot_name>.json. Returns the config file path."""
    print_fn("")
    print_fn("=== GZ-toolkit potential-testing setup ===")
    print_fn("Press Enter to accept the value shown in [brackets].")
    print_fn("")

    # ---- potential identity -------------------------------------------------
    pot_name = _ask("Short name for this potential (used for folder names)",
                    default="my_potential",
                    validate=lambda s: None if all(c.isalnum() or c in "._-" for c in s)
                    else "Use only letters, numbers, '.', '_' or '-'.",
                    input_fn=input_fn, print_fn=print_fn)

    print_fn("")
    print_fn("Paste the LAMMPS pair lines for this potential, e.g.")
    print_fn("  pair_style eam/fs")
    print_fn("  pair_coeff * * Fe.eam.fs Fe")
    print_fn("Type each line, then an empty line to finish.")
    pot_lines_list: list[str] = []
    while True:
        line = input_fn("pair line (empty = done): ").strip()
        if not line:
            if pot_lines_list:
                break
            print_fn("  At least one pair_style line is needed.")
            continue
        pot_lines_list.append(line)
    pot_lines = "\n".join(pot_lines_list)

    print_fn("")
    print_fn("Potential file(s) to copy into every run dir, e.g. Fe.eam.fs")
    print_fn("(some potentials need more than one file, e.g. a .meam + a library.meam).")
    print_fn("Type each filename, then an empty line to finish (leave blank right away "
              "if everything is in the pair lines).")
    potential_files: list[str] = []
    while True:
        fname = input_fn("potential file (empty = done): ").strip()
        if not fname:
            break
        potential_files.append(fname)

    # ---- elements ------------------------------------------------------------
    print_fn("")
    elements = _ask("Element symbols this potential covers, comma-separated "
                    "(metals first, then gases — e.g. Fe,Cr,He)",
                    parse=lambda s: [e.strip() for e in s.replace(";", ",").split(",") if e.strip()],
                    validate=lambda lst: None if lst else "Give at least one element.",
                    input_fn=input_fn, print_fn=print_fn)

    type_map: dict[str, int] = {}
    masses: dict[int, float] = {}
    crystal_structures: dict[str, dict] = {}
    gases: list[str] = []
    metals: list[str] = []

    for i, el in enumerate(elements, start=1):
        type_map[el] = i
        m = suggest_mass(el)
        mass = _ask(f"Atomic mass of {el} (amu)",
                    default=round(m, 4) if m is not None else None,
                    parse=float,
                    validate=lambda v: None if v > 0 else "Mass must be positive.",
                    input_fn=input_fn, print_fn=print_fn)
        masses[i] = float(mass)
        is_gas = _ask_yesno(f"Is {el} a gas species (He, H, ...)?",
                            default=el in ("He", "H", "Ne", "Ar", "Kr", "Xe", "D", "T"),
                            input_fn=input_fn, print_fn=print_fn)
        if is_gas:
            gases.append(el)
        else:
            metals.append(el)
            struct = _ask(f"Crystal structure of pure {el}",
                          default=_STRUCT_GUESS.get(el, "bcc"),
                          parse=lambda s: s.lower(), choices=("bcc", "fcc", "hcp"),
                          input_fn=input_fn, print_fn=print_fn)
            crystal_structures[el] = {"structure": struct}

    lc_default = _LC_GUESS.get(metals[0], 2.85) if metals else 2.85
    lc_initial = _ask("Starting lattice-constant guess in Angstrom "
                      "(just a seed; the code relaxes it)",
                      default=lc_default, parse=float,
                      validate=lambda v: None if v > 0 else "Must be positive.",
                      input_fn=input_fn, print_fn=print_fn)

    # ---- which tests -----------------------------------------------------------
    print_fn("")
    print_fn("--- Which tests do you want? ---")
    wf = WorkflowOptions()
    # point defects always run (lc + Ecoh + formation energies are the core).
    print_fn("Lattice constant, cohesive energy and point-defect formation")
    print_fn("energies always run — they are the core of the test.")

    wf.include_elastic = _ask_yesno(
        "Elastic constants (full 6x6 Cij) for each pure element?",
        default=True, input_fn=input_fn, print_fn=print_fn)

    two_element_suites: list[dict] = []
    if len(metals) >= 2:
        if _ask_yesno(f"Alloy lattice-constant scan for {metals[0]}-{metals[1]}?",
                      default=False, input_fn=input_fn, print_fn=print_fn):
            wf.include_alloy_suite = True
            fracs = _ask(f"Atomic percent of {metals[1]} to test, comma-separated",
                         default="3,5,8,10,50",
                         parse=lambda s: [float(x) for x in s.split(",") if x.strip()],
                         validate=lambda lst: None if all(0 < f < 100 for f in lst)
                         else "Each value must be between 0 and 100.",
                         input_fn=input_fn, print_fn=print_fn)
            two_element_suites.append({
                "A": metals[0], "B": metals[1],
                "fractions_atpct": fracs, "ordering": "random",
            })
            if wf.include_elastic:
                wf.elastic_for_alloys = _ask_yesno(
                    "Also compute elastic constants for each alloy composition?\n"
                    "  (one random arrangement on a small cell — treat as an estimate)",
                    default=False, input_fn=input_fn, print_fn=print_fn)

    multi_element_suites: list[dict] = []
    if len(metals) >= 3:
        if _ask_yesno(f"Also test a fixed {'/'.join(metals)} multi-element "
                      f"composition (3+ elements)?",
                      default=False, input_fn=input_fn, print_fn=print_fn):
            wf.include_alloy_suite = True
            print_fn("Give the at.% for each element; leave blank to split the "
                      "remainder evenly across every blank entry.")
            comp: dict[str, float | None] = {}
            for el in metals:
                while True:
                    raw = input_fn(f"  {el} at.% (blank = rest): ").strip()
                    if not raw:
                        comp[el] = None
                        break
                    try:
                        comp[el] = float(raw)
                        break
                    except ValueError:
                        print_fn(f"  '{raw}' is not a number — leave blank for "
                                  f"'rest' or type a number.")
            multi_element_suites.append({"composition": comp, "ordering": "random"})
            if wf.include_elastic:
                wf.elastic_for_alloys = wf.elastic_for_alloys or _ask_yesno(
                    "Also compute elastic constants for this composition?\n"
                    "  (one random arrangement on a small cell — treat as an estimate)",
                    default=False, input_fn=input_fn, print_fn=print_fn)

    if gases:
        wf.include_gas_complexes = _ask_yesno(
            f"Gas-defect complexes (e.g. {gases[0]}-vacancy clusters) with binding energies?",
            default=False, input_fn=input_fn, print_fn=print_fn)

    wf.include_dislocation_loops = _ask_yesno(
        "Dislocation loops (interstitial loops, bcc metals only)?",
        default=False, input_fn=input_fn, print_fn=print_fn)

    wf.seakmc_enabled = _ask_yesno(
        "Defect migration barriers via SEAKMC (needs seakmc_p on the cluster)?",
        default=False, input_fn=input_fn, print_fn=print_fn)

    # ---- cluster ------------------------------------------------------------------
    print_fn("")
    print_fn("--- Cluster settings ---")
    print_fn("Defaults below match the ISAAC/mtp_env setup — press Enter to accept "
              "each one, or type a new value.")
    hpc = HPCOptions()

    hpc.mode = _ask("Submission mode (scatter/batch)", default="scatter",
                    choices=("scatter", "batch"), input_fn=input_fn, print_fn=print_fn)
    cluster = _ask("Cluster name ('none' for a generic SLURM script)",
                   default="ISAAC", input_fn=input_fn, print_fn=print_fn)
    hpc.cluster = None if cluster.lower() == "none" else cluster
    partition_key = _ask("Partition key ('none' for no partition flag)",
                         default="s", input_fn=input_fn, print_fn=print_fn)
    hpc.partition_key = None if partition_key.lower() == "none" else partition_key
    hpc.walltime = _ask("Time limit per job (HH:MM:SS)", default="3:00:00",
                        validate=lambda s: None if re.match(r"^\d+:\d{2}:\d{2}$", s)
                        else "Format must be like 3:00:00.",
                        input_fn=input_fn, print_fn=print_fn)
    hpc.nodes = _ask("Nodes per job", default=1, parse=int,
                     validate=lambda v: None if v >= 1 else "At least 1.",
                     input_fn=input_fn, print_fn=print_fn)
    hpc.ntasks = _ask("CPU cores per job", default=48, parse=int,
                      validate=lambda v: None if v >= 1 else "At least 1.",
                      input_fn=input_fn, print_fn=print_fn)
    hpc.output_name = _ask("Job script filename", default="run.job",
                           input_fn=input_fn, print_fn=print_fn)
    hpc.batch_max_per_file = _ask("Max sims batched per .job file", default=50, parse=int,
                                  validate=lambda v: None if v >= 1 else "At least 1.",
                                  input_fn=input_fn, print_fn=print_fn)
    hpc.max_jobs_per_batch = _ask("Hard cap on sims per batch job (site-imposed)",
                                  default=200, parse=int,
                                  validate=lambda v: None if v >= 1 else "At least 1.",
                                  input_fn=input_fn, print_fn=print_fn)
    hpc.run_command = _ask("LAMMPS run command "
                           "({{NTASKS}}, {{INPUT}}, {{LAMMPS_BIN}} are placeholders)",
                           default="srun -n {{NTASKS}} {{LAMMPS_BIN}} -in {{INPUT}}",
                           input_fn=input_fn, print_fn=print_fn)
    hpc.device = _ask("Device", default="cpu", choices=("cpu", "gpu"),
                      input_fn=input_fn, print_fn=print_fn)
    hpc.gpus_per_node = _ask("GPUs per node", default=0, parse=int,
                             validate=lambda v: None if v >= 0 else "Cannot be negative.",
                             input_fn=input_fn, print_fn=print_fn)
    modules_raw = _ask("Modules to load, comma-separated (blank for none)",
                       default="", input_fn=input_fn, print_fn=print_fn)
    hpc.modules = [m.strip() for m in modules_raw.split(",") if m.strip()]
    conda_env = _ask("Conda env to activate ('none' to skip)",
                     default="/lustre/isaac24/scratch/qzhang55/myVenvs/mtp_env",
                     input_fn=input_fn, print_fn=print_fn)
    hpc.conda_env = None if conda_env.lower() == "none" else conda_env
    lammps_binary_path = _ask("Path to lmp_mpi binary ('none' to use PATH)",
                              default="~/bins/lmp_mpi",
                              input_fn=input_fn, print_fn=print_fn)
    hpc.lammps_binary_path = None if lammps_binary_path.lower() == "none" else lammps_binary_path
    lammps_src_path = _ask("LAMMPS src/ dir for LD_LIBRARY_PATH ('none' to skip)",
                           default="/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src",
                           input_fn=input_fn, print_fn=print_fn)
    hpc.lammps_src_path = None if lammps_src_path.lower() == "none" else lammps_src_path

    print_fn("")
    print_fn("How parallel should the runs be?")
    print_fn("  1 = everything in ONE job (slowest, simplest)")
    print_fn("  2 = one job per potential")
    print_fn("  3 = one job per test group (defects, elastic, alloys, ...)")
    print_fn("  4 = one job per simulation (fastest; recommended)")
    wf.paral_degree = _ask("Parallelization level", default=4, parse=int,
                           choices=(1, 2, 3, 4), input_fn=input_fn, print_fn=print_fn)

    # ---- assemble + validate + save -------------------------------------------------
    cfg = PotentialConfig(
        potential=PotentialMetadata(
            pot_name=pot_name,
            pot_lines=pot_lines,
            type_map=type_map,
            masses=masses,
            single_elements=metals,
            gases=gases,
            two_element_suites=two_element_suites,
            multi_element_suites=multi_element_suites,
            crystal_structures=crystal_structures,
            lc_initial=float(lc_initial),
            potential_files=potential_files,
        ),
        workflow=wf,
        hpc=hpc,
    )

    from gz_toolkit.potential_testing.validate import validate_config
    problems = validate_config(cfg)

    # scaffold the project, then overwrite the starter JSON with our answers
    from gz_toolkit.potential_testing.project import init_potential_project
    root_path = init_potential_project(root_dir=root, pot_name=pot_name)
    cfg_path = Path(root_path) / "pot_inputs" / f"{pot_name}.json"
    save_potential_config(cfg, cfg_path)

    print_fn("")
    print_fn(f"Config written to: {cfg_path}")
    if problems:
        print_fn("Heads-up — the validator flagged these:")
        for p in problems:
            print_fn(f"  - {p}")
    print_fn("")
    print_fn("Next steps:")
    if potential_files:
        print_fn(f"  1. Copy {', '.join(potential_files)} into potentials/{pot_name}/")
        print_fn("  2. Sanity-check them:    python check.py")
        print_fn("  3. On the cluster, run:  python main.py")
        print_fn("  4. Then submit:          bash submit_all.sh")
        print_fn("  5. Check progress:       python -m gz_toolkit.potential_testing.cli status")
        print_fn("  6. When done, summarize: python summarize.py")
    else:
        print_fn("  1. On the cluster, run:  python main.py")
        print_fn("  2. Then submit:          bash submit_all.sh")
        print_fn("  3. Check progress:       python -m gz_toolkit.potential_testing.cli status")
        print_fn("  4. When done, summarize: python summarize.py")
    return cfg_path
