"""Copies the BRENDA release that geckopy downloads into GECKO's databases folder.

Usage: python scripts/ci/stage_brenda.py repos/GECKO

ec_model_full_yeastgem reads databases/brenda-2026.1-1 of the GECKO checkout on both
sides. The folder is not part of the repository, and the Python side runs before MATLAB
could fetch it, so geckopy downloads the files (checked against the toolbox-assets
manifest) and they are copied to where GECKO's fetchBrendaAsset looks first.
"""

import shutil
import sys
from pathlib import Path

from geckopy.databases.brenda_asset import fetch_brenda_asset

FILES = ("kcat.tsv", "sa.tsv", "mw.tsv", "measurements.tsv")

source = fetch_brenda_asset(FILES)
target = Path(sys.argv[1]) / "databases" / "brenda-2026.1-1"
target.mkdir(parents=True, exist_ok=True)
for name in FILES:
    shutil.copyfile(source / name, target / name)
print(f"BRENDA files in {target}")
