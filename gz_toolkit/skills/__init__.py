"""Install the bundled materials-simulation skill for a local agent."""

from __future__ import annotations

import argparse
from importlib.resources import files
import os
from pathlib import Path
import tempfile

SKILL_NAME = "materials-simulation-skill"


def skill_path(destination: str | Path | None = None) -> Path:
    """Return the optional installation path without creating any files."""
    if destination is None:
        codex_home = Path(os.environ.get("CODEX_HOME") or Path.home() / ".codex")
        destination = codex_home / "skills"
    return Path(destination).expanduser().resolve() / SKILL_NAME


def _copy_resources(source, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.name == "__pycache__" or item.name.endswith((".pyc", ".pyo")):
            continue
        target = destination / item.name
        if item.is_dir():
            _copy_resources(item, target)
        else:
            target.write_bytes(item.read_bytes())


def install_skill(destination: str | Path | None = None, *, force: bool = False) -> Path:
    """Copy the complete skill to an agent's skills root without running its scripts.

    The default is $CODEX_HOME/skills, or ~/.codex/skills. Pass another skills
    root explicitly for a different agent. Existing installs require force=True.
    """
    target = skill_path(destination)
    root = target.parent
    if target.is_symlink():
        raise ValueError(f"Refusing to replace a linked skill directory: {target}")
    if target.exists() and (not force or not target.is_dir()):
        raise FileExistsError(f"Skill already exists: {target}. Use --force to replace it.")

    source = files(__package__).joinpath(SKILL_NAME)
    if not source.joinpath("SKILL.md").is_file():
        raise FileNotFoundError("Bundled SKILL.md is missing; reinstall GZ-toolkit.")
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gz-toolkit-skill-", dir=root) as staging:
        staged = Path(staging) / SKILL_NAME
        _copy_resources(source, staged)
        backup = Path(staging) / "previous"
        if target.exists():
            target.rename(backup)
        try:
            staged.rename(target)
        except OSError:
            if backup.exists():
                backup.rename(target)
            raise
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Install GZ-toolkit's materials-simulation skill.")
    parser.add_argument("--dest", type=Path, help="Agent skills root (default: $CODEX_HOME/skills or ~/.codex/skills)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--force", action="store_true", help="Replace an existing installation of this skill")
    mode.add_argument("--check", action="store_true", help="Report whether the skill is installed; do not create or change files")
    args = parser.parse_args(argv)
    if args.check:
        target = skill_path(args.dest)
        installed = (target / "SKILL.md").is_file()
        print(f"{'Installed' if installed else 'Not installed'}: {target}")
        return 0 if installed else 1
    try:
        target = install_skill(args.dest, force=args.force)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Could not install skill: {exc}\n")
    print(f"Installed {SKILL_NAME}: {target}")
    print("The skill will be available on your next Codex turn.")
    return 0
