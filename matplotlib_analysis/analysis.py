from pathlib import Path
import re
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from io import StringIO

import matplotlib.pyplot as plt

bp_to_inspect = 5
hexamer_sequence = "ACGTAC"

root_dir = Path(__file__).resolve().parents[1]
prev_dir = root_dir / "previous_iteration"
hexamer_csv = (
    root_dir
    / "Offset_energy_calculator"
    / "hexamers_csv"
    / "DNA"
    / "coords_grooves_DNA_hexamers_table.csv"
)

param_keys = ("rise", "twist", "shift", "slide")


def _find_column(columns, candidates):
    for col in columns:
        low = col.lower()
        if any(c in low for c in candidates):
            return col
    return None


def _read_avg_params(path, bp_index):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Unexpected format in {path}")
    header = lines[1].lstrip("#").strip()
    data = "\n".join([header] + lines[2:])
    df = pd.read_csv(
        StringIO(data),
        delim_whitespace=True,
        engine="python",
    )
    bp_col = _find_column(df.columns, ("set", "bp", "base", "step", "index"))
    if bp_col is None:
        bp_col = df.columns[0]
    row = df[df[bp_col] == bp_index]
    if row.empty:
        raise ValueError(f"bp {bp_index} not found in {path}")
    row = row.iloc[0]
    out = {}
    for key in param_keys:
        col = _find_column(df.columns, (key,))
        if col is None:
            raise ValueError(f"Missing column for {key} in {path}")
        out[key] = float(row[col])
    return out


def _read_equilibrium_params(path, sequence):
    df = pd.read_csv(path)
    seq_col = _find_column(df.columns, ("seq", "hexamer"))
    if seq_col is None:
        raise ValueError("Hexamer sequence column not found")
    row = df[df[seq_col].str.upper() == sequence.upper()]
    if row.empty:
        raise ValueError(f"Hexamer {sequence} not found in {path}")
    row = row.iloc[0]
    out = {}
    for key in param_keys:
        col = _find_column(df.columns, (key,))
        if col is None:
            raise ValueError(f"Missing column for {key} in {path}")
        out[key] = float(row[col])
    return out


def _find_md_file(iter_dir):
    candidates = [
        iter_dir / "MD_Results" / "MD_averaged_paramers_of_helix_0.dat",
        iter_dir / "MD_Results" / "MD_averaged_parameters_of_helix_0.dat",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"No averaged parameter file in {iter_dir}")


iter_dirs = sorted(
    prev_dir.glob("iteration_*"),
    key=lambda p: int(re.findall(r"\d+", p.name)[0]),
)

iterations = [int(re.findall(r"\d+", p.name)[0]) for p in iter_dirs]
avg_by_param = {k: [] for k in param_keys}

for d in iter_dirs:
    md_file = _find_md_file(d)
    params = _read_avg_params(md_file, bp_to_inspect)
    for k in param_keys:
        avg_by_param[k].append(params[k])

eq_params = _read_equilibrium_params(hexamer_csv, hexamer_sequence)
avg_helical_coords = {
    k: list(zip(avg_by_param[k], [eq_params[k]] * len(iterations)))
    for k in param_keys
}

_original_bar = Axes.bar


def _bar_with_equilibrium(self, x, height, *args, **kwargs):
    try:
        first = height[0]
    except Exception:
        return _original_bar(self, x, height, *args, **kwargs)
    if isinstance(first, (tuple, list)) and len(first) == 2:
        avg = np.array([h[0] for h in height], dtype=float)
        eq = np.array([h[1] for h in height], dtype=float)
        width = kwargs.pop("width", 0.35)
        x = np.array(x, dtype=float)
        _original_bar(self, x - width / 2, avg, width=width, *args, **kwargs)
        _original_bar(self, x + width / 2, eq, width=width, color="orange")
        return None
    return _original_bar(self, x, height, *args, **kwargs)


Axes.bar = _bar_with_equilibrium

# # Example averaged helical coordinates for iterations 0-3
# # Replace with real values as needed
# avg_helical_coords = {
#     "rise": [3.2, 3.3, 3.1, 3.4],
#     "twist": [34.5, 34.7, 34.6, 34.8],
#     "shift": [0.1, 0.2, 0.15, 0.18],
#     "slide": [0.05, 0.04, 0.06, 0.05],
# }

# iterations = [0, 1, 2, 3]

print(avg_helical_coords)
print(eq_params)
fig, axes = plt.subplots(len(avg_helical_coords), 1, figsize=(6, 8), sharex=True)

for ax, (param, values) in zip(axes, avg_helical_coords.items()):
    ax.bar(iterations, values, color="steelblue")
    ax.set_ylabel(param)
    ax.set_xticks(iterations)

axes[-1].set_xlabel("Iteration")
fig.suptitle("Averaged Helical Coordinates (Iterations 0-3)")
plt.tight_layout()
plt.show()
