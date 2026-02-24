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




n_params = 3 # Select top n most impactful parameters
range_to_inspect = range(5, 12)  # Range of base pair indices to inspect





# Base-pair index to read from each MD averaged parameters file.
iterations = [0,1,2,3]


# Project root (two levels up from this file).
root_dir = Path(__file__).resolve().parents[1]
# Folder that contains the iteration_* directories.
prev_dir = root_dir / "previous_iteration"


# Collect iteration directories and sort by their numeric suffix.
iter_dirs = sorted(
    prev_dir.glob("iteration_*"),
    key=lambda p: int(re.findall(r"\d+", p.name)[0]),
)

bp_keys = ['shear', 'stretch', 'stagger', 'buckle', 'prop', 'open']
step_keys = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
heli_keys = ['xdisp', 'ydisp', 'hrise', 'incl', 'tip', 'htwist']

# Read helical parameters and initialize helices
helix_by_iteration = [] # this is for now only one helix per iteration which needs to be extended later
for d in iter_dirs:
    md_file = file = d / "MD_Results" / "MD_averaged_parameters_of_helix_0.dat"

    helix = calc.extract_data(md_file)
    helix['energys'], helix['stiffs'], helix['eq_params'], helix['differences'] = calc.calculate_displacement_energy(helix)
    helix_by_iteration.append(helix)

# Calculate total energy contribution per parameter type across all helices/iterations
energy_by_param = {k: 0.0 for k in bp_keys + step_keys + heli_keys}
for helix in helix_by_iteration:
    for i, k in enumerate(bp_keys):
        energy_by_param[k] += sum(row[i] for row in helix['energys']['bp'][range_to_inspect.start:range_to_inspect.stop])
    for i, k in enumerate(step_keys):
        energy_by_param[k] += sum(row[i] for row in helix['energys']['step'][range_to_inspect.start:range_to_inspect.stop])
    for i, k in enumerate(heli_keys):
        energy_by_param[k] += sum(row[i] for row in helix['energys']['heli'][range_to_inspect.start:range_to_inspect.stop])
# Sort parameters by total energy (descending) and select top n
sorted_params = sorted(energy_by_param.items(), key=lambda x: x[1], reverse=True)


pdf_path = "matplotlib_analysis/Equalibrium_Offset_Vectors_graph.pdf"

# for param_type, energys in helix_by_iteration[2]['energys'].items():
#     print(param_type)
#     for bp_idx, energy_values in enumerate(energys):
#         print(f"  Base pair {bp_idx}: {energy_values}")
# print()

with PdfPages(pdf_path) as pdf:

    for param_keys in [tuple(k for k, _ in sorted_params[i:i+n_params]) for i in range(0, len(sorted_params), n_params)]:
        # Stores averages per parameter across iterations.
        # Structure: {param: [iteration_data, ...]} where iteration_data is [(param, eql_param) for each bp in range]
        avg_by_param = {k: [] for k in param_keys}



        for k in param_keys:
            for helix in helix_by_iteration:
                seq = helix['strand_sequences'][0]
                iteration_data = []
                for place_to_inspect in range_to_inspect:
                    if k in bp_keys:
                        #heptamer_sequence = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(place_to_inspect-3,place_to_inspect+4)]) 
                        param = helix['bp_params'][place_to_inspect][bp_keys.index(k)]
                        eql_param = helix['eq_params']['bp'][place_to_inspect][bp_keys.index(k)]
                        stiffs = helix['stiffs']['bp'][place_to_inspect][bp_keys.index(k)]
                        iteration_data.append( (param, eql_param, stiffs) )
                    elif k in step_keys:
                        #hexamer_sequence = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(place_to_inspect-2,place_to_inspect+4)]) 
                        param = helix['step_params'][place_to_inspect][step_keys.index(k)]
                        eql_param = helix['eq_params']['step'][place_to_inspect][step_keys.index(k)]
                        stiffs = helix['stiffs']['step'][place_to_inspect][step_keys.index(k)]
                        iteration_data.append( (param, eql_param, stiffs) )
                    elif k in heli_keys:
                        #hexamer_sequence = ''.join([ seq[i] if (i in range(len(seq))) else '-' for i in range(place_to_inspect-2,place_to_inspect+4)]) 
                        param = helix['heli_params'][place_to_inspect][heli_keys.index(k)]
                        eql_param = helix['eq_params']['heli'][place_to_inspect][heli_keys.index(k)]
                        stiffs = helix['stiffs']['heli'][place_to_inspect][heli_keys.index(k)]
                        iteration_data.append( (param, eql_param, stiffs) )
                avg_by_param[k].append(iteration_data)


        # Create one subplot per parameter, stacked vertically and sharing x-axis.
        fig, axes = plt.subplots(len(avg_by_param), 1, figsize=(6, 8), sharex=True)

        num_bp_in_range = len(range_to_inspect)
        colors = plt.cm.tab10(np.linspace(0, 1, num_bp_in_range))

        for ax, (param, values) in zip(axes, avg_by_param.items()):
            # values is a list of iterations, each containing [(param, eql_param, stiffs) for each bp in range]
            
            # Create x-positions: for each iteration, place base pairs side by side
            spacing = 0.8 / num_bp_in_range  # Space allocated for all bp within one iteration
            
            arrow_width_factor = 0.1/max([x[2] for iteration_data in values for x in iteration_data])  # Adjust this factor to control arrow width
            for iter_idx, iteration_data in enumerate(values):
                base_x = iterations[iter_idx]
                for bp_offset, (param_val, eql_val, stiff_val) in enumerate(iteration_data):
                    # Position base pairs next to each other within the iteration
                    x_pos = base_x - 0.4 + spacing * bp_offset + spacing / 2
                    dy = param_val - eql_val  # vertical displacement
                    ax.arrow(x_pos, eql_val, 0, dy, head_width=arrow_width_factor*stiff_val, head_length=abs(dy)*0.25 if dy != 0 else 0.1,
                            fc=colors[bp_offset], ec=colors[bp_offset], length_includes_head=True, width=arrow_width_factor*stiff_val)
            
            # Label y-axis with the parameter name.
            ax.set_ylabel(str(param) + " " + str(energy_by_param[param]/sum(energy_by_param.values())*100)[:4] + "%")
            # Show every iteration index as a tick.
            ax.set_xticks(iterations)
            # Add horizontal grid
            ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.5)

        # Label the shared x-axis on the bottom subplot.
        axes[-1].set_xlabel("Iteration")

        # Display the sequence with the inspected range highlighted
        for i, helix in enumerate(helix_by_iteration):
            fig.text(0.19*i+0.15, 0.9, helix['strand_sequences'][0] , fontsize=6,fontname='monospace',)
            for k, j in enumerate(range_to_inspect):
                range_to_inspect_seq = ''.join([ c if l == j else " " for l,c in enumerate(helix['strand_sequences'][0])])
                fig.text(0.19*i+0.15, 0.9, range_to_inspect_seq , fontsize=6, fontname='monospace', color=colors[k])


        # Add vertical dividers between iterations
        for i in [0.5+j for j in range(3)]:
            divider_display = axes[0].transData.transform((i, 0))
            divider_fig_x = fig.transFigure.inverted().transform(divider_display)[0]
            fig.add_artist(Line2D([divider_fig_x, divider_fig_x], [0.1, 0.93], transform=fig.transFigure, color="black", linewidth=0.5))


        fig.suptitle("Restrained                     non Restrained")
        pdf.savefig(fig, bbox_inches="tight")
        #plt.close(fig)
