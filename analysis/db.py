"""Shared SQLite connection helper for the analysis scripts."""

import sqlite3
from pathlib import Path

DB_PATH = Path(__file__).parent.parent / "cell_counts.db"
OUTPUTS_DIR = Path(__file__).parent.parent / "outputs"


def get_connection():
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"{DB_PATH.name} not found. Run `python load_data.py` first."
        )
    return sqlite3.connect(DB_PATH)


def ensure_outputs_dir():
    OUTPUTS_DIR.mkdir(exist_ok=True)
    return OUTPUTS_DIR
