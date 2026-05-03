"""
col_normalize.py
----------------
Shared column-normalization helper used by EDA, modeling, and report scripts.
Standardizes DataFrame column names to UPPER_SNAKE_CASE so downstream code
can use consistent names regardless of the raw source format.
"""

import re
import pandas as pd


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to UPPER_SNAKE_CASE.

    Rules applied in order:
    1. Strip leading/trailing whitespace.
    2. Replace spaces, hyphens, dots, and consecutive underscores with a single underscore.
    3. Convert to uppercase.
    4. Strip leading/trailing underscores.

    Returns a *new* DataFrame with renamed columns (does not mutate in place).
    """
    def _clean(name):
        if not isinstance(name, str):
            name = str(name)
        name = name.strip()
        name = re.sub(r"[\s\-\.]+", "_", name)
        name = re.sub(r"_+", "_", name)
        name = name.upper()
        name = name.strip("_")
        return name

    new_cols = [_clean(c) for c in df.columns]
    # Handle duplicate column names after normalization
    seen = {}
    deduped = []
    for c in new_cols:
        if c in seen:
            seen[c] += 1
            deduped.append(f"{c}_{seen[c]}")
        else:
            seen[c] = 0
            deduped.append(c)
    df = df.copy()
    df.columns = deduped
    return df
