from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import re
import matplotlib
from matplotlib.patches import ConnectionPatch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from io import StringIO
import matplotlib.pyplot as plt
import sys
import os, os.path
# Adds the parent directory (project/) to the path
sys.path.append(str(Path(__file__).resolve().parents[1]))
import Offset_energy_calculator.calculate as calc



# Base-pair index to read from each MD averaged parameters file.
iterations = [0,1,2,3]


# Project root (two levels up from this file).
root_dir = Path(__file__).resolve().parents[1]
# Folder that contains the iteration_* directories.
prev_dir = root_dir / "result_data"

plt.figure(figsize=(10, 6))
sequence_dirs = sorted(
    prev_dir.glob("sequence_*"),
    key=lambda p: int(re.findall(r"\d+", p.name)[0]),
)

energy_by_sequence = []
for i, sequence_dir in enumerate(sequence_dirs):

    # Collect iteration directories and sort by their numeric suffix.
    iter_dirs = sorted(
        sequence_dir.glob("iteration_*"),
        key=lambda p: int(re.findall(r"\d+", p.name)[0]),
    )

    energy_by_iteration = []
    with open ( sequence_dir / "energy_log.txt" ) as f:
        for line in f:
            energy = float(line.split()[2])
            energy_by_iteration.append(energy)
    energy_by_sequence.append(energy_by_iteration)

avg_energy_by_sequence = [0, 0, 0, 0]
for strt_seq in energy_by_sequence:
    for i in range(4):
        avg_energy_by_sequence[i] += strt_seq[i]

avg_energy_by_sequence = [x / len(energy_by_sequence) for x in avg_energy_by_sequence]

plt.plot(iterations, avg_energy_by_sequence, label="Average Energy")

# for i in range(len(energy_by_sequence)):
#     plt.plot(iterations, energy_by_sequence[i], label=f"Sequence {i}")

plt.legend()
# plt.ylim(10, 50)
plt.xlabel("Sequence Index")
plt.ylabel("Energy")
plt.title("Energy Results")
plt.savefig(root_dir / "matplotlib_analysis" / "energy_results.pdf", bbox_inches="tight")

#pdf.savefig(fig, bbox_inches="tight")
