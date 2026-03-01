import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
from scipy.optimize import curve_fit


#Define the theoretical function
def onsager_m(beta):
    beta_c = 0.44068679
    if beta < beta_c:
        return 0.0
    else:
        # m = (1 - [sinh(2*beta)]^-4)^(1/8)
        term = np.sinh(2 * beta)**(-4)
        return (1.0 - term)**(0.125)
        
        
mag_files= [('magbeta/magbeta60.txt', 'L=60', 'blue'),('magbeta/magbeta80.txt', 'L=80', 'red'),('magbeta/magbeta100.txt', 'L=100', 'green'),('magbeta/magbeta120.txt', 'L=120', 'magenta'),('magbeta/magbeta140.txt', 'L=140', 'orange'),('magbeta/magbeta160.txt', 'L=160', 'cyan'),('magbeta/magbeta180.txt', 'L=180', 'gray')]

chi_files = [('chibeta/chibeta60.txt', 'L=60', 'blue'),('chibeta/chibeta80.txt', 'L=80', 'red'),('chibeta/chibeta100.txt', 'L=100', 'green'),('chibeta/chibeta120.txt', 'L=120', 'magenta'),('chibeta/chibeta140.txt', 'L=140', 'orange'),('chibeta/chibeta160.txt', 'L=160', 'cyan'),('chibeta/chibeta180.txt', 'L=180', 'gray')]

chip_files = [('chipbeta/chipbeta60.txt', 'L=60', 'blue'),('chipbeta/chipbeta80.txt', 'L=80', 'red'),('chipbeta/chipbeta100.txt', 'L=100', 'green'),('chipbeta/chipbeta120.txt', 'L=120', 'magenta'),('chipbeta/chipbeta140.txt', 'L=140', 'orange'),('chipbeta/chipbeta160.txt', 'L=160', 'cyan'),('chipbeta/chipbeta180.txt', 'L=180', 'gray')]

u_files = [('ubeta/ubeta60.txt', 'L=60', 'blue'),('ubeta/ubeta80.txt', 'L=80', 'red'),('ubeta/ubeta100.txt', 'L=100', 'green'),('ubeta/ubeta120.txt', 'L=120', 'magenta'),('ubeta/ubeta140.txt', 'L=140', 'orange'),('ubeta/ubeta160.txt', 'L=160', 'cyan'),('ubeta/ubeta180.txt', 'L=180', 'gray')]





#----------- FIGURE 1: ANSAMBLE AVARAGE OF MAGNETIZATION (<|m|>) --------------
plt.figure(1, figsize=(12,6))

for filename, label, color in mag_files:
	df = pd.read_csv(filename, sep =',', header = None, names = ['beta', 'm', 'error'] )
	df = df.sort_values('beta')
	plt.plot(df['beta'], df['m'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)


#Plot the Theoretical Function (L -> infinity)
beta_range = np.linspace(0.42, 0.47, 5000)
theory_m = [onsager_m(b) for b in beta_range]
plt.plot(beta_range, theory_m, color='black', linestyle='--', linewidth=0.5, label='Theory ($L \\to \\infty$)', zorder=10)

# styling
plt.title('Ising Model: Simulation vs. Onsager Exact Solution')
plt.xlabel('Inverse Temperature ($\\beta$)')
plt.ylabel('Average Absolute Magnetization $\\langle |m| \\rangle$')
#plt.grid(True, linestyle='--', color='black', linewidth= 0.35, alpha=0.5)

# Place the legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Manually adjust the right margin (0.8 means the plot takes 80% of the width)
plt.subplots_adjust(right=0.8)



# ------- FIGURE 2: SUSCEPTIBILITY (CHI) --------
plt.figure(2, figsize=(12, 6))

for filename, label, color in chi_files:
    df_chi = pd.read_csv(filename, sep=',', header=None, names=['beta', 'chi'])
    df_chi = df_chi.sort_values('beta')
    df_chi['chi'] = df_chi['chi'] / 1000
    
    plt.plot(df_chi['beta'], df_chi['chi'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

#styling
plt.title('Ising Model: Magnetic Susceptibility $\\chi$')
plt.xlabel('Inverse Temperature ($\\beta$)')
plt.ylabel('Susceptibility $\\chi / \\ 10^3$')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)



# --- FIGURE 3: SUBTRACTED SUSCEPTIBILITY (CHIP) ---
plt.figure(3, figsize=(12, 6))

for filename, label, color in chip_files:
    df_chip = pd.read_csv(filename, sep=',', header=None, names=['beta', 'chip', 'error'])
    df_chip = df_chip.sort_values('beta')
    df_chip['chip'] = df_chip['chip'] / 100
    plt.plot(df_chip['beta'], df_chip['chip'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

#style
plt.title('Ising Model: Subtracted Susceptibility $\\chi\'$')
plt.xlabel('Inverse Temperature ($\\beta$)')
plt.ylabel('Subtracted Susceptibility $\\chi\' /  \\ 10^2 $')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)



# ------- FIGURE 4: BINDER CUMULANT (U) --------
plt.figure(4, figsize=(12, 6))

for filename, label, color in u_files:
    df_bin = pd.read_csv(filename, sep=',', header=None, names=['beta', 'U4'])
    df_bin = df_bin.sort_values('beta')
    plt.plot(df_bin['beta'], df_bin['U4'], 'o-', markersize=3, linewidth=0.6, label=label, color=color, alpha=0.5)

# Theoretical intersection 
plt.axhline(y=1.1679227, color='black', linestyle='--', linewidth='0.5', label='$U_4^* = 1.1679$', alpha=0.8)

#style
plt.title('Ising Model: Binder Cumulant $U$')
plt.xlabel('Inverse Temperature ($\\beta$)')
plt.ylabel('Binder Cumulant $U$')
#plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

#legendj
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(right=0.80)

plt.show()




#--------------------------DATA ANALYSIS FOR CRITICAL EXPONENTS-------------------------------
peak_eights = []
peak_eights_error = []
beta_pc = []
beta_pc_error = []
L_values = []

starting_index = 3

plt.figure(5, figsize=(12, 6))

def fit_function(x, a, beta_pc, chi_prime_max):
	return a * ( (x - beta_pc)**2 ) + chi_prime_max

def reduced_Chi_squared(y_data, y_model, errors, number_parameters):
	residuals = (y_data - y_model) / errors
	chi_sq = np.sum(residuals**2)
	dof = len(y_data) - number_parameters
	return chi_sq / dof

for counter in range(starting_index, len(chip_files)):
	filename, label, color = chip_files[counter]
	current_L = int(label.split('=')[1])
	L_values.append(current_L)
	df_chip = pd.read_csv(filename, sep=',', header=None, names=['beta', 'chip', 'error'])
	df_chip = df_chip.sort_values('beta')
	
	
	chi_prime = df_chip['chip'].to_numpy()
	beta = df_chip['beta'].to_numpy()
	chi_prime_error = df_chip['error'].to_numpy()
	
	
	#mask condition, takes only array values that satisfy mask or the beta t.c. mask is satisfied for the chi_prime corrispondent
	mask = chi_prime > (0.8 * np.max(chi_prime))
	beta_fit = beta[mask]
	chi_prime_fit = chi_prime[mask]
	chi_prime_error_fit = chi_prime_error[mask]
	
	#dynamic guessing of p0
	beta_pc_tentative = beta_fit[np.argmax(chi_prime_fit)]
	chi_prime_max_tentative = np.max(chi_prime_fit)
	a_tentative = - chi_prime_max_tentative / 0.001
	
	popt, pcov = curve_fit(fit_function, beta_fit, chi_prime_fit, p0 = [a_tentative, beta_pc_tentative, chi_prime_max_tentative], sigma = chi_prime_error_fit, absolute_sigma = True)
	
	beta_pc.append(popt[1])
	beta_pc_error.append(np.sqrt(pcov[1,1]))
	
	peak_eights.append(popt[2])
	peak_eights_error.append(np.sqrt(pcov[2,2]))
	
	#Chi squared reduced
	chi_prime_model = fit_function(beta_fit, *popt)
	Chi_squared_value = reduced_Chi_squared(chi_prime_fit, chi_prime_model, chi_prime_error_fit, 3)
	
	
	#plotting fits
	plt.figure(5)
	plt.errorbar(beta_fit, chi_prime_fit, yerr = chi_prime_error_fit, fmt = 'o', color=color, markersize=4, label=label, alpha=0.6)
	
	beta_model_line = np.linspace(beta_fit.min(), beta_fit.max(), 1000)
	chi_model_line = fit_function(beta_model_line, *popt)

	plt.plot(beta_model_line, chi_model_line, '-', color='black', linewidth=0.5, label=f"{label}: $\\chi^2_{{red}} = {Chi_squared_value:.4f}$", alpha = 0.6)
	
# style
plt.figure(5)
plt.title('Ising Model: Parabolic Fits of $\\chi\'$ Peaks')
plt.xlabel('Inverse Temperature ($\\beta$)')
plt.ylabel('Subtracted Susceptibility $\\chi\'$')
plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

# Legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.subplots_adjust(right=0.80)
plt.show()




#----------------- FITTING FOR CRITICAL EXPONENT: NU ---------------------

def beta_pc_fit_function(x, beta_c, b, nu):
	return beta_c + b * (x**(-1/nu))

beta_pc = np.array(beta_pc)
beta_pc_error = np.array(beta_pc_error)
L_values = np.array(L_values)

popt, pcov = curve_fit(beta_pc_fit_function, L_values, beta_pc, p0 = [0.4406, -0.1, 1.0], sigma = beta_pc_error, absolute_sigma= True)

beta_c_fit = popt[0]
nu_fit = popt[2]
beta_c_error = np.sqrt(pcov[0,0])
nu_error = np.sqrt(pcov[2,2])

print(f"--- FSS Results: beta_pc vs L fit for nu and beta_c---")
print(f"Beta_c : {beta_c_fit:.6f} +/- {beta_c_error:.6f}")
print(f"Critical exponent nu: {nu_fit:.4f} +/- {nu_error:.4f}")

plt.figure(6, figsize=(12,6))

plt.errorbar(L_values, beta_pc, yerr= beta_pc_error, fmt='o', markersize=4, alpha=0.6)
L_plot = np.linspace(L_values.min()-10, L_values.max()+10, 1000)
plt.plot(L_plot, beta_pc_fit_function(L_plot, *popt), '-', label='fit : $\\beta_c + b (x^{ -1/\\nu })$', color ='black', linewidth=0.5)

plt.xlabel('Lattice Size $L$')
plt.ylabel('$\\beta_{pc}(L)$')

plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

# Legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.subplots_adjust(right=0.80)

plt.show()






#--------------------------DATA ANALYSIS FOR CRITICAL EXPONENTS: GAMMA/NU-------------------------------

def chi_prime_max_fit_function(x, c_0, c_1, gamma_over_nu):
	return c_0 + c_1 * x**(gamma_over_nu)

#making lists into array(L values are already array from previous fit)
peak_eights = np.array(peak_eights)
peak_eights_error = np.array(peak_eights_error)

popt, pcov = curve_fit(chi_prime_max_fit_function, L_values, peak_eights, p0 = [0.0, 1.0, 1.75], sigma = peak_eights_error, absolute_sigma= True)

gamma_over_nu = popt[2]
gamma_over_nu_error = np.sqrt(pcov[2,2])

print(f"--- FSS Results: chi_prime_max vs L fit for nu and gamma---")
print(f"gamma / nu : {gamma_over_nu:.6f} +/- {gamma_over_nu_error:.6f}")


plt.figure(7, figsize=(12,6))

plt.errorbar(L_values, peak_eights, yerr= peak_eights_error, fmt='o', markersize=4, alpha=0.6)
L_plot = np.linspace(L_values.min()-10, L_values.max()+10, 1000)
plt.plot(L_plot, chi_prime_max_fit_function(L_plot, *popt), '-', label='fit : $c_0 + c_1 (x^{ \\gamma/\\nu })$', color ='black', linewidth=0.5)

plt.xlabel('Lattice Size $L$')
plt.ylabel('$\\chi\'_{max}(L)$')

plt.grid(True, linestyle='--', color='black', linewidth=0.35, alpha=0.5)

# Legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.subplots_adjust(right=0.80)

plt.show()





#--------------------------<|M|> ERROR VS LATTICE SIZE L-------------------------------

plt.figure(8, figsize=(12,6))

for filename, label, color in mag_files:
	df = pd.read_csv(filename, sep =',', header = None, names = ['beta', 'm', 'error'] )
	df = df.sort_values('beta')
	
	current_L = int(label.split('=')[1])
	avg_error = df['error'].mean() * 10**3
	
	
	plt.plot(current_L, avg_error, 'o', markersize=5, color='black', alpha=1)



# styling
plt.title('Ising Model: mean error scaling with lattice size')
plt.xlabel('Lattice size ($L$)')
plt.ylabel('mean error rescaled $\\sigma_{<|m|>} \\cdot 10^3 $')
plt.ylim(0,1.5)
plt.grid(True, linestyle='--', color='black', linewidth= 0.35, alpha=0.5)

plt.show()


















