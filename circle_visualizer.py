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

# Loop through the remaining (empty) axes and hide them
for i in range(len(N_t_list), len(flattened_axes_list)):
    flattened_axes_list[i].axis('off')

#two global axis for every plot
figure_mosaic.supxlabel("Block Size")
figure_mosaic.supylabel(r"Jackknife Error $\sigma_{\langle Q^2 \rangle}$")

# layout adjustment: rect says talk only to the subplots: left, bottom, right, top
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])


plt.savefig("results/plateau.pdf") 
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
ax1.errorbar(eta2_values, Q2_means, yerr=Q2_errors, fmt='o', markersize=4, mfc='none', mec='red', ecolor='red', elinewidth=0.5, capsize=3, capthick=0.5)
ax1.hlines(y=2, xmin=0, xmax=0.001, color='black', linestyle='--', linewidth=0.5)
ax1.set_ylim(1.25, 2.25)
ax1.set_xlabel(r"$\eta^2$", fontsize=14)
ax1.set_ylabel(r"$\langle Q^2 \rangle$", fontsize=14)
#ax1.grid(False, linestyle='--', color = 'black', alpha=0.3)

#Subplot zoom on the critical zone
#Define position: [x, y, width, height] relative to the parent axes (0 to 1)
axins = ax1.inset_axes([0.55, 0.55, 0.4, 0.4])

# Data plotted in the inset
axins.errorbar(eta2_values, Q2_means, yerr=Q2_errors, fmt='o', markersize=4, mfc='none', mec='red', ecolor='red', elinewidth=0.5, capsize=3, capthick=0.5)

#Define the ZOOM zone (Limit the view to where data is amassed)
#change these values to match the specific X and Y range of cluster
x1, x2 = 0.0, 0.0005  # <--- Change these X limits to zoom in on your cluster
y1, y2 = 1.97, 2.03  # <--- Change these Y limits to zoom in on your cluster

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins.hlines(y=2, xmin=0, xmax=0.001, color='black', linestyle='--', linewidth=0.5)

# Force scientific notation (like 10^-4)
axins.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Rotate labels 45 degrees and make font smaller
axins.tick_params(axis='x', labelrotation=45)

#Error vs N_t
ax2.plot(N_t_values, Q2_errors, 'o', markersize=3, color='red')
ax2.set_xlabel(r"$N_t$", fontsize=14)
ax2.set_ylabel(r"$\sigma_{\langle Q^2 \rangle}$", fontsize=14)
ax2.set_yscale('log')
#ax2.grid(False, linestyle='--', color='black', alpha=0.3)

plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])

plt.savefig("results/results_analysis.pdf")
plt.close()



#======================================
#   Plotting <Q^2> vs beta' 
#   high temperature scan
#======================================



#extracting data from results file: temperature, <Q^2>, error
results_high_temperature = np.loadtxt("data1/results.txt")
beta_prime_values        = results_high_temperature[:, 1]
Q2_mean_values           = results_high_temperature[:, 2]
Q2_error_values          = results_high_temperature[:, 3]


#plots
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

#temperature span for plot
beta_prime_theory = np.linspace(0.15, 0.5, 1000)

#function to evaluate the high temperature approximation of Q^2
Q2_high_T_approximation = 2.0 * np.exp(-1.0 / (2.0 * beta_prime_theory))

#plotting points
ax.errorbar(
    beta_prime_values, Q2_mean_values, yerr=Q2_error_values,
    fmt='o', markersize=4, mfc='none', mec='red',
    ecolor='red', elinewidth=0.5, capsize=3, capthick=0.5,
    label=r"$\langle Q^2 \rangle$  ($N_t = 20$)"
)
ax.plot(beta_prime_theory, Q2_high_T_approximation, 'k--', linewidth=1.0,
        label=r"$2\,e^{-1/(2\beta')}$ (high-T approx.)")
ax.set_xlabel(r"$\beta'$", fontsize=14)
ax.set_ylabel(r"$\langle Q^2 \rangle$", fontsize=14)
ax.legend(fontsize=9, loc='upper left')


plt.title(r"High-T regime: $\langle Q^2 \rangle$ vs $\beta'$", fontsize=13)


#inset zoom
# [x, y, width, height] in axes coordinates
axins = ax.inset_axes([0.55, 0.1, 0.4, 0.5])

# Plot same data in inset
axins.errorbar(
    beta_prime_values, Q2_mean_values, yerr=Q2_error_values,
    fmt='o', markersize=3, mfc='none', mec='red',
    ecolor='red', elinewidth=0.5, capsize=2, capthick=0.5
)

#plot theory function for high T and make the axis numbers small
axins.plot(beta_prime_theory, Q2_high_T_approximation, 'k--', linewidth=1.0)
axins.tick_params(labelsize=8)

#set zoom zone
axins.set_xlim(0.17, 0.30)  
axins.set_ylim(0.15, 0.25) 


plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])

plt.savefig("results/high_temp_analysis.pdf")
plt.close()

