"""
tree.py - Generate a readable directory tree for agent workflows.

The default base folder is the current working directory.
"""

import os


def _safe_listdir(path, include_hidden=False):
    try:
        names = sorted(os.listdir(path))
    except OSError:
        return []

    if include_hidden:
        return names
    return [n for n in names if not n.startswith(".")]


def _build_tree_lines(base_path, max_depth=4, include_files=True, include_hidden=False):
    base_path = os.path.abspath(base_path)
    lines = [base_path]

    def walk(path, prefix, depth):
        if depth >= max_depth:
            return

        entries = _safe_listdir(path, include_hidden=include_hidden)
        kept = []
        for name in entries:
            full = os.path.join(path, name)
            if os.path.isdir(full) or include_files:
                kept.append((name, full, os.path.isdir(full)))

        for i, (name, full, is_dir) in enumerate(kept):
            last = i == len(kept) - 1
            branch = "└── " if last else "├── "
            suffix = "/" if is_dir else ""
            lines.append(f"{prefix}{branch}{name}{suffix}")
            if is_dir:
                next_prefix = prefix + ("    " if last else "│   ")
                walk(full, next_prefix, depth + 1)

    walk(base_path, "", 0)
    return lines


def generate_folder_tree(
    base_folder=None,
    path=None,
    output_file="folder_tree.txt",
    max_depth=4,
    include_files=True,
    include_hidden=False,
):
    """
    Generate and write a folder tree text file.

    Parameters
    ----------
    base_folder : str or None
        Base folder to scan. If None, uses current working directory.
    path : str or None
        Alias of ``base_folder``. If both are given, ``path`` is used.
    output_file : str
        Output text file path. Relative paths resolve under base folder.
    max_depth : int
        Maximum directory depth to traverse.
    include_files : bool
        If True, include files in the tree; otherwise directories only.
    include_hidden : bool
        If True, include hidden entries (leading '.').

    Returns
    -------
    str
        Absolute path to the written tree file.
    """
    if path is not None:
        base_folder = path
    if base_folder is None:
        base_folder = os.getcwd()
    base_folder = os.path.abspath(base_folder)

    if not os.path.isdir(base_folder):
        raise FileNotFoundError(f"Base folder '{base_folder}' not found.")
    if max_depth < 1:
        raise ValueError("max_depth must be >= 1.")

    if os.path.isabs(output_file):
        out_path = output_file
    else:
        out_path = os.path.join(base_folder, output_file)
    out_path = os.path.abspath(out_path)

    lines = _build_tree_lines(
        base_path=base_folder,
        max_depth=int(max_depth),
        include_files=bool(include_files),
        include_hidden=bool(include_hidden),
    )

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Wrote folder tree: {out_path}")
    return out_path
