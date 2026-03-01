import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

beta = 0.125
gamma = 1.75
beta_c = 0.4407

mag_files= [('magbeta/magbeta60.txt', 'L=60', 'blue'),('magbeta/magbeta80.txt', 'L=80', 'red'),('magbeta/magbeta100.txt', 'L=100', 'green'),('magbeta/magbeta120.txt', 'L=120', 'magenta'),('magbeta/magbeta140.txt', 'L=140', 'orange'),('magbeta/magbeta160.txt', 'L=160', 'cyan'),('magbeta/magbeta180.txt', 'L=180', 'gray')]

chi_files = [('chibeta/chibeta60.txt', 'L=60', 'blue'),('chibeta/chibeta80.txt', 'L=80', 'red'),('chibeta/chibeta100.txt', 'L=100', 'green'),('chibeta/chibeta120.txt', 'L=120', 'magenta'),('chibeta/chibeta140.txt', 'L=140', 'orange'),('chibeta/chibeta160.txt', 'L=160', 'cyan'),('chibeta/chibeta180.txt', 'L=180', 'gray')]

chip_files = [('chipbeta/chipbeta60.txt', 'L=60', 'blue'),('chipbeta/chipbeta80.txt', 'L=80', 'red'),('chipbeta/chipbeta100.txt', 'L=100', 'green'),('chipbeta/chipbeta120.txt', 'L=120', 'magenta'),('chipbeta/chipbeta140.txt', 'L=140', 'orange'),('chipbeta/chipbeta160.txt', 'L=160', 'cyan'),('chipbeta/chipbeta180.txt', 'L=180', 'gray')]

u_files = [('ubeta/ubeta60.txt', 'L=60', 'blue'),('ubeta/ubeta80.txt', 'L=80', 'red'),('ubeta/ubeta100.txt', 'L=100', 'green'),('ubeta/ubeta120.txt', 'L=120', 'magenta'),('ubeta/ubeta140.txt', 'L=140', 'orange'),('ubeta/ubeta160.txt', 'L=160', 'cyan'),('ubeta/ubeta180.txt', 'L=180', 'gray')]


#-------FIGURE 1: ENSAMBLE AVARAGE |M|---------
plt.figure(1, figsize=(12,6))

for filename, label, color in mag_files:
	L = int(label.split('=')[1])
	df = pd.read_csv(filename, sep =',', header = None, names = ['beta', 'm', 'error'] )
	df = df.sort_values('beta')
	df['m'] = df['m'] * ( L**beta)
	df['beta'] = (df['beta'] - beta_c) * L
	plt.plot(df['beta'], df['m'],'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

# styling
plt.title('Ising Model: FSS $\\langle |m| \\rangle$')
plt.xlabel(' ($(\\beta - \\beta_c)L^{1/\\nu}$)')
plt.ylabel('$\\langle |m| \\rangle L^{\\gamma/\\nu }$')
#plt.grid(True, linestyle='--', color='black', linewidth= 0.35, alpha=0.5)

#legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.subplots_adjust(right=0.8)



# ------- FIGURE 2: SUSCEPTIBILITY (CHI) --------
plt.figure(2, figsize=(12, 6))

for filename, label, color in chi_files:
	L = int(label.split('=')[1])
	df_chi = pd.read_csv(filename, sep=',', header=None, names=['beta', 'chi'])
	df_chi = df_chi.sort_values('beta')
	df_chi['chi'] = df_chi['chi'] / (L**gamma)
	df_chi['beta'] = (df_chi['beta'] - beta_c) * L  
	plt.plot(df_chi['beta'], df_chi['chi'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

#styling
plt.title('Ising Model: FSS $\\chi$')
plt.xlabel(' ($(\\beta - \\beta_c)L^{1/\\nu}$)')
plt.ylabel('$\\chi / L^{\\gamma/\\nu}$')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)
	


# --- FIGURE 3: SUBTRACTED SUSCEPTIBILITY (CHIP) ---
plt.figure(3, figsize=(12, 6))

for filename, label, color in chip_files:
	L = int(label.split('=')[1])
	df_chip = pd.read_csv(filename, sep=',', header=None, names=['beta', 'chip', 'error'])
	df_chip = df_chip.sort_values('beta')
	df_chip['chip'] = df_chip['chip'] / (L**gamma)
	df_chip['beta'] = (df_chip['beta'] - beta_c) * L 
	plt.plot(df_chip['beta'], df_chip['chip'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

#style
plt.title('Ising Model: FSS Subtracted Susceptibility $\\chi\'$')
plt.xlabel(' ($(\\beta - \\beta_c)L^{1/\\nu}$)')
plt.ylabel(' $\\chi\' /  L^{\\gamma/\\nu} $')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)



# ------- FIGURE 4: BINDER CUMULANT (U) --------
plt.figure(4, figsize=(12, 6))

for filename, label, color in u_files:
	L = int(label.split('=')[1])
	df_bin = pd.read_csv(filename, sep=',', header=None, names=['beta', 'U4'])
	df_bin = df_bin.sort_values('beta')
	df_bin['beta'] = (df_bin['beta'] - beta_c) * L 
	plt.plot(df_bin['beta'], df_bin['U4'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

# Theoretical intersection 
plt.axhline(y=1.1679227, color='black', linestyle='--', linewidth='0.5', label='$U_4^* = 1.1679$', alpha=0.8)

#style
plt.title('Ising Model: Binder Cumulant $U$')
plt.xlabel('($(\\beta - \\beta_c)L^{1/\\nu}$)')
plt.ylabel('$U$')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legendj
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)

plt.show()
