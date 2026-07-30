"""Thin shim: the real CLI lives in cli.py.

Kept so the console-script entry point
(``gz-toolkit-potential-testing = gz_toolkit.potential_testing.main:main``)
and any ``python -m gz_toolkit.potential_testing.main`` invocations keep
working after the CLI consolidation.
"""

from __future__ import annotations

from gz_toolkit.potential_testing.cli import build_parser, main

__all__ = ["build_parser", "main"]


if __name__ == "__main__":
    raise SystemExit(main())
