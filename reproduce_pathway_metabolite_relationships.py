#!/usr/bin/env python3
from __future__ import annotations

import runpy
from pathlib import Path


SCRIPT = Path(__file__).with_name("make_pathway_relationship_figures.py")


if __name__ == "__main__":
    runpy.run_path(str(SCRIPT), run_name="__main__")
