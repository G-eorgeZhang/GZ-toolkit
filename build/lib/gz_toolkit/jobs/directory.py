"""
directory.py — Navigate and filter simulation run directories.

Supports flexible directory naming conventions so users with different
folder structures can all use the same traversal logic.
"""

import os
import re
import glob as _glob


class DirectoryManager:
    """
    Traverse and filter subdirectories under a base folder.

    Supports multiple naming conventions:
      - "prefix-suffix"  e.g. "0-Fe", "12-He"   (default)
      - "flat"           e.g. any folder name, sorted alphabetically
      - custom regex     e.g. r"temp_(\\d+)" for "temp_300", "temp_600"

    Parameters
    ----------
    base_folder : str
        Path to the parent directory that contains the run subdirectories.
        Can be relative (resolved against ``base_dir``) or absolute.
    base_dir : str or None
        Root directory to resolve ``base_folder`` against.
        If None, uses the current working directory.
    naming : str or compiled regex
        How to interpret subdirectory names for sorting and filtering.
        Built-in options:

        - ``"prefix-suffix"`` — expects ``<int>-<anything>``, sorts/filters
          by the integer prefix.  This is the default.
        - ``"flat"`` — treat every subdirectory as valid, sort alphabetically.
        - A compiled ``re.Pattern`` — must have a named group ``(?P<key>...)``
          that captures the sortable part.  Example:
          ``re.compile(r"temp_(?P<key>\\d+)")`` matches ``temp_300``.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, base_folder, base_dir=None, naming="prefix-suffix"):
        if base_dir is None:
            base_dir = os.getcwd()

        self.base_dir = os.path.abspath(base_dir)
        self.folder_path = (
            base_folder
            if os.path.isabs(base_folder)
            else os.path.join(self.base_dir, base_folder)
        )

        if not os.path.isdir(self.folder_path):
            raise FileNotFoundError(
                f"Base folder '{self.folder_path}' not found."
            )

        # Store naming mode
        if isinstance(naming, re.Pattern):
            if "key" not in naming.groupindex:
                raise ValueError(
                    "Custom regex must contain a named group (?P<key>...)."
                )
            self._naming = "regex"
            self._regex = naming
        elif naming in ("prefix-suffix", "flat"):
            self._naming = naming
            self._regex = None
        else:
            raise ValueError(
                f"Unknown naming mode '{naming}'. "
                "Use 'prefix-suffix', 'flat', or a compiled re.Pattern."
            )

    # ------------------------------------------------------------------
    # Key extraction (used for sorting and range filtering)
    # ------------------------------------------------------------------

    def _extract_key(self, folder_name):
        """
        Return the sortable key for a subdirectory name.

        Returns None if the name does not match the expected pattern
        (the folder will be skipped).
        """
        if self._naming == "prefix-suffix":
            if "-" not in folder_name:
                return None
            prefix_str = folder_name.split("-", 1)[0]
            try:
                return int(prefix_str)
            except ValueError:
                return None

        elif self._naming == "flat":
            return folder_name  # every directory is valid

        else:  # regex
            m = self._regex.match(folder_name)
            if m is None:
                return None
            raw = m.group("key")
            try:
                return int(raw)
            except ValueError:
                return raw

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def get_subdirs(self, selectors=None):
        """
        Return a sorted list of subdirectory *names* that match ``selectors``.

        Parameters
        ----------
        selectors : None, str, list[str], or tuple[int, int]
            - ``None``            — include every valid subdirectory.
            - ``str``             — include only that exact name.
            - ``list[str]``       — include those exact names.
            - ``(lo, hi)``        — include subdirs whose numeric key is
              in ``[lo, hi]`` (only works with numeric keys).
            - ``"glob:<pattern>"``— include subdirs matching a glob pattern
              relative to the base folder (e.g. ``"glob:temp_*"``).

        Returns
        -------
        list[str]
            Sorted subdirectory names.
        """
        all_items = os.listdir(self.folder_path)

        # Determine mode
        if selectors is None:
            mode = "all"
        elif isinstance(selectors, str):
            if selectors.startswith("glob:"):
                mode = "glob"
                pattern = selectors[5:]
            else:
                mode = "list"
                sel_set = {selectors}
        elif isinstance(selectors, (list, set)):
            mode = "list"
            sel_set = set(selectors)
        elif (
            isinstance(selectors, tuple)
            and len(selectors) == 2
            and all(isinstance(x, (int, float)) for x in selectors)
        ):
            mode = "range"
            lo, hi = int(selectors[0]), int(selectors[1])
        else:
            raise ValueError(
                "`selectors` must be None, a str, a list of str, "
                "a tuple of two ints, or a 'glob:<pattern>' string."
            )

        if mode == "glob":
            # Use glob to match directory names
            full_pattern = os.path.join(self.folder_path, pattern)
            matched = [
                os.path.basename(p)
                for p in _glob.glob(full_pattern)
                if os.path.isdir(p)
            ]
            matched.sort(key=lambda n: (self._extract_key(n) is None, self._extract_key(n) or n))
            return matched

        filtered = []
        for item in all_items:
            item_path = os.path.join(self.folder_path, item)
            if not os.path.isdir(item_path):
                continue

            key = self._extract_key(item)
            if key is None:
                continue

            if mode == "all":
                filtered.append((key, item))
            elif mode == "list":
                if item in sel_set:
                    filtered.append((key, item))
            elif mode == "range":
                if isinstance(key, (int, float)) and lo <= key <= hi:
                    filtered.append((key, item))

        # Sort by key
        filtered.sort(key=lambda pair: pair[0])
        return [name for _, name in filtered]

    def iter_subdirs(self, selectors=None):
        """
        Yield ``(subdir_name, subdir_full_path)`` for each matched subdirectory.
        """
        for name in self.get_subdirs(selectors):
            yield name, os.path.join(self.folder_path, name)

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def subdir_path(self, name):
        """Return the full path for a single subdirectory name."""
        return os.path.join(self.folder_path, name)

    def __repr__(self):
        return (
            f"DirectoryManager(folder_path='{self.folder_path}', "
            f"naming='{self._naming}')"
        )
