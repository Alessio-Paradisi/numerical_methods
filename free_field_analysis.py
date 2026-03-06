# This code reads the output of free_field.c and performs the thermodynamic analysis:
#
# data_O{1,2,3}.txt contain rows:  mean, error, N_T, N_S
# rows with N_T = N_S are zero temperature references (set by free_field.py)
# rows with N_T < N_S are finite temperature measurements
#
# Steps:
#  1. Read data and separate finite temperature from zero temperature rows
#  2. Subtract zero temperature reference to remove UV divergences
#  3. Compute renormalized energy density epsilon/T^D (eq. 15.1.11):
#       epsilon/T^D = (N_T^D / 2) * (<O1+O2-O3>_Nt - <O1+O2-O3>_Nt_bar)
#  4. Compute trace anomaly (epsilon - P)/T^D (eq. 15.2.19):
#       (epsilon - P)/T^2 = N_T^2 * (<O1>_Nt - <O1>_Nt_bar)    (D=2)
#  5. Plot epsilon/T^2, (epsilon-P)/T^2



#======================================
#		Libraries
#======================================


import numpy as np
import matplotlib.pyplot as plt



#======================================
#		Parameters
#======================================


D = 2

N_T_values = [4, 5, 6, 7, 8, 10]

ratios = [4, 6, 8]



#======================================
#		Reading and separating data
#         with functions
#======================================


def read_pairs(filename):
  rows = []
  with open(filename, "r") as f:
    for line in f:
      parts = line.strip().split(",")
      rows.append([float(parts[0]), float(parts[1]), int(parts[2]), int(parts[3])])
  return [(rows[i], rows[i+1]) for i in range(0, len(rows), 2)]


pairs_O1 = read_pairs("data/data_O1.txt")
pairs_O2 = read_pairs("data/data_O2.txt")
pairs_O3 = read_pairs("data/data_O3.txt")



#======================================
#	Computing renormalized observables
#======================================



results = []

for (r1f, r1z), (r2f, r2z), (r3f, r3z) in zip(pairs_O1, pairs_O2, pairs_O3):

  m1f, e1f, N_T, N_S = r1f
  m1z, e1z = r1z[0], r1z[1]
  m2f, e2f = r2f[0], r2f[1]
  m2z, e2z = r2z[0], r2z[1]
  m3f, e3f = r3f[0], r3f[1]
  m3z, e3z = r3z[0], r3z[1]
  
  # renormalized O1+O2-O3 combination
  combination_renorm = (m1f + m2f - m3f) - (m1z + m2z - m3z)
  combination_renorm_error = np.sqrt(e1f**2 + e2f**2 + e3f**2 + e1z**2 + e2z**2 + e3z**2)

  # renormalized O1 alone for trace anomaly
  O1_renorm = m1f - m1z
  O1_renorm_error = np.sqrt(e1f**2 + e1z**2)


  energy_density = (N_T**D / 2.0) * combination_renorm
  energy_density_error = (N_T**D / 2.0) * combination_renorm_error

  trace_anomaly = (N_T**D) * O1_renorm
  trace_anomaly_error = (N_T**D) * O1_renorm_error  

  results.append((energy_density, energy_density_error, trace_anomaly, trace_anomaly_error, N_T, N_S))

results.sort(key=lambda row: (row[4], row[0]))




#======================================
#       Plotting the results
#======================================


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Group results by ratio
by_ratio = {4: [], 6: [], 8: []}

for r in results:
	ratio = r[5] // r[4]
	by_ratio[ratio].append(r)

# Markers for each ratio
markers = {4: 'o', 6: 's', 8: 'D'}
colors = {4: 'blue', 6: 'red', 8:'green'}

# --- Left Plot: Energy Density ---
for ratio, data in by_ratio.items():
	x = [r[4]**(-2) for r in data]
	y = [r[0] for r in data]
	err = [r[1] for r in data]
	ax1.errorbar(x, y, yerr=err, fmt=markers[ratio], mfc='none', mec=colors[ratio], capsize=3, markersize=6, label=f'$N_S/N_T={ratio}$')
	coeffs = np.polyfit(x, y, 1)
	x_fit = np.linspace(0, max(x), 100)
	ax1.plot(x_fit, np.polyval(coeffs, x_fit), linestyle='dotted', color=colors[ratio], alpha=0.6)

ax1.set_xlabel('$1/N_T^2$', fontsize=12)
ax1.set_ylabel('$\\epsilon/T^2$', fontsize=12)
ax1.set_title('Energy Density vs $1/N_T^2$', fontsize=13)
ax1.legend(loc='best')
ax1.axhline(y=0.406123498885, color='black', linestyle='dashed')


# --- Right Plot: Trace Anomaly ---
for ratio, data in by_ratio.items():
	x = [r[4]**(-2) for r in data]
	y = [r[2] for r in data]
	err = [r[3] for r in data]
	ax2.errorbar(x, y, yerr=err, fmt=markers[ratio], mfc='none', mec=colors[ratio], capsize=3, markersize=6, label=f'$N_S/N_T={ratio}$')
	coeffs = np.polyfit(x, y, 1)
	x_fit = np.linspace(0, max(x), 100)
	ax2.plot(x_fit, np.polyval(coeffs, x_fit), linestyle='dotted', color=colors[ratio], alpha=0.6)

ax2.set_xlabel('$1/N_T^2$', fontsize=12)
ax2.set_ylabel('$(\\epsilon-P)/T^2$', fontsize=12)
ax2.set_title('Trace Anomaly vs $1/N_T^2$', fontsize=13)
ax2.legend(loc='best')
ax2.set_ylim(0.18, 0.21)
ax2.axhline(y=0.1866576057, color='black', linestyle='dashed')

plt.tight_layout()
plt.show()



