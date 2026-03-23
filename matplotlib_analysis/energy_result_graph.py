import shutil

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

from uncertainties import ufloat



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




# Function to count entries in a file (assuming entries are lines)
def count_entries(file_path):
    with open(file_path, 'r') as f:
        return sum(1 for line in f)

# Filter out directories that contain files with fewer than n entries
filtered_sequence_dirs = []
for dir_path in sequence_dirs:
    if os.path.exists(dir_path):
        files = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
        remove_dir = False
        for file in files:
            file_path = os.path.join(dir_path, file)
            if count_entries(file_path) < 4:
                remove_dir = True
                break
        if not remove_dir:
            filtered_sequence_dirs.append(dir_path)
        else:
            # Remove the directory
            shutil.rmtree(dir_path)

sequence_dirs = filtered_sequence_dirs

# ...existing code...


energy_by_sequence = []
sequences = []

for i, sequence_dir in enumerate(sequence_dirs):

    # Collect iteration directories and sort by their numeric suffix.
    iter_dirs = sorted(
        sequence_dir.glob("iteration_*"),
        key=lambda p: int(re.findall(r"\d+", p.name)[0]),
    )
    with open ( iter_dirs[0] / "cpptraj_base_pairing.txt" ) as f:
        for line in f:
            sequences.append(line.split()[2])
            break
    energy_by_iteration = []
    intermed_energy = None
    with open ( sequence_dir / "energy_log.txt" ) as f:
        for line_idx, line in enumerate(f):
            nom_val = float(line.split()[2].split("+/-")[0])
            std_val = float(line.split()[2].split("+/-")[1])
            energy = ufloat(nom_val, std_val)
            if line_idx == 0:
                nom_val = float(line.split()[3].split("+/-")[0])
                std_val = float(line.split()[3].split("+/-")[1])
                intermed_energy = ufloat(nom_val, std_val)
            energy_by_iteration.append(energy)
    energy_by_iteration.insert(1, intermed_energy)
    energy_by_sequence.append(energy_by_iteration)


# avg_energy_by_sequence = [0, 0, 0, 0]
# for strt_seq in energy_by_sequence:
#     for i in range(4):
#         avg_energy_by_sequence[i] += strt_seq[i]

# avg_energy_by_sequence = [x / len(energy_by_sequence) for x in avg_energy_by_sequence]

# # plt.plot(iterations, avg_energy_by_sequence, label="Average Energy")


iterations.insert(1, 0.5)
for i in range(len(energy_by_sequence)):
    plt.plot(iterations[:3], [e.nominal_value for e in energy_by_sequence[i][:3]], label=f"{sequences[i]}", c = plt.cm.tab10(i))
    plt.errorbar(iterations[:3], [e.nominal_value for e in energy_by_sequence[i][:3]], yerr=[e.std_dev for e in energy_by_sequence[i][:3]], c = plt.cm.tab10(i))
    plt.plot(iterations[3:], [e.nominal_value for e in energy_by_sequence[i][3:]], c = plt.cm.tab10(i))
    plt.errorbar(iterations[3:], [e.nominal_value for e in energy_by_sequence[i][3:]], yerr=[e.std_dev for e in energy_by_sequence[i][3:]], c = plt.cm.tab10(i))

plt.legend(title="Starting Sequence:", prop={'family': 'monospace'})
plt.ylim(4, 45)
plt.xlabel("Iteration")
plt.ylabel("$E_{offset} (kcal/mol)$")
plt.savefig(root_dir / "matplotlib_analysis" / "energy_results.pdf", bbox_inches="tight")

