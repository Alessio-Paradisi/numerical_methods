#This code is the visualizer for the particle on a circle.
#
#It has the plateau plots and  the winding number evolution.
#
#
#======================================
#		Libraries
#======================================


import math
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


#======================================
#		Plotting plateaus mosaic
#======================================


#N_t_list = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]		#number of sites in the euclidean_time_lattice 
N_t_list = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325]

total_plots_count = len(N_t_list)
columns_number = math.ceil(math.sqrt(total_plots_count))
rows_number = math.ceil(total_plots_count / columns_number)


# Creating the figure and the grid of subplots
figure_mosaic, axes_matrix = plt.subplots(rows_number, columns_number,  figsize=(4 * columns_number, 3.5 * rows_number),squeeze=False)
flattened_axes_list = axes_matrix.flatten()

for plot_index, current_N_t in enumerate(N_t_list):
  current_axis = flattened_axes_list[plot_index]
  
  #loading data from file
  plateau_data = np.loadtxt(f"data/plateau_{current_N_t}.txt")
  block_size_values = plateau_data[:, 0]
  error_values = plateau_data[:, 1]

  #plotting
  current_axis.plot(block_size_values, error_values, '-o', markersize = 2, linewidth = 0.5, color = 'black')
  current_axis.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
  current_axis.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

  current_axis.set_title(f"$N_t = {current_N_t}$")
  current_axis.grid(True, linestyle='--', color = 'black',  alpha=0.3)

#two global axis for every plot
figure_mosaic.supxlabel("Block Size")
figure_mosaic.supylabel(r"Jackknife Error $\sigma_{\langle Q^2 \rangle}$")

# layout adjustment: rect says talk only to the subplots: left, bottom, right, top
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])


plt.savefig("data/plateau.pdf") 
plt.close()


#======================================
#		Plotting <Q^2> and sigma_<Q^2>
#======================================


results = np.loadtxt("data/results.txt")
N_t_values = results[:, 0]
Q2_means = results[:, 1]
Q2_errors = results[:, 2]

eta2_values = (2 / N_t_values) * (2 / N_t_values) 

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot <Q^2> vs eta^2
ax1.errorbar(eta2_values, Q2_means, yerr=Q2_errors, fmt='o', markersize=2, color='black', ecolor='red', elinewidth=0.5, capsize=3, capthick=0.5)
ax1.set_ylim(1.25, 2.25)
ax1.set_xlabel(r"$\eta^2$", fontsize=14)
ax1.set_ylabel(r"$\langle Q^2 \rangle$", fontsize=14)
#ax1.grid(False, linestyle='--', color = 'black', alpha=0.3)

#Subplot zoom on the critical zone
#Define position: [x, y, width, height] relative to the parent axes (0 to 1)
axins = ax1.inset_axes([0.55, 0.55, 0.4, 0.4])

# Data plotted in the inset
axins.errorbar(eta2_values, Q2_means, yerr=Q2_errors, fmt='o', markersize=2, color='black', ecolor='red', elinewidth=0.5, capsize=3, capthick=0.5)

#Define the ZOOM zone (Limit the view to where data is amassed)
#change these values to match the specific X and Y range of cluster
x1, x2 = 0.0, 0.0005  # <--- Change these X limits to zoom in on your cluster
y1, y2 = 1.97, 2.03  # <--- Change these Y limits to zoom in on your cluster

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Force scientific notation (like 10^-4)
axins.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Rotate labels 45 degrees and make font smaller
axins.tick_params(axis='x', labelrotation=45)

#Error vs N_t
ax2.plot(N_t_values, Q2_errors, 'o', markersize=3, color='black')
ax2.set_xlabel(r"$N_t$", fontsize=14)
ax2.set_ylabel(r"$\sigma_{\langle Q^2 \rangle}$", fontsize=14)
ax2.set_yscale('log')
#ax2.grid(False, linestyle='--', color='black', alpha=0.3)

plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])

plt.savefig("data/results_analysis.pdf")
plt.close()



















