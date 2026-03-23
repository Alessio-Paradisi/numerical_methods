


#=============================
#      libraries
#=============================


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit



#=============================
#     extract the data
#        from file 
#=============================


#loading the results from the file in different arrays
data = np.loadtxt("results/results2.txt", delimiter=",")

lattice_size                          = data[:, 0].astype(int)
beta                                  = data[:, 1]
energy                                = data[:, 2]
energy_error                          = data[:, 3]
energy_variance                       = data[:, 4]
energy_variance_error                 = data[:, 5]
absolute_magnetization                = data[:, 6]
absolute_magnetization_error          = data[:, 7]
absolute_magnetization_variance       = data[:, 8]
absolute_magnetization_variance_error = data[:, 9]
magnetization2                        = data[:, 10]
magnetization2_error                  = data[:, 11]
binder                                = data[:, 12]
binder_error                          = data[:, 13]


lattice_sizes = sorted(set(lattice_size))


#markers list to use in plots
markers    = ['o', 'o', 'o', 's', '^', 'v', 'D', 'o']
markermap = {l: m for l, m in zip(lattice_sizes, markers)}

#colors list to use in plots
colors = [
    '#000000',   # L=20 black
    '#555555',   # L=40  dark gray
    '#5C00A3',   # L=60  purple
    '#0050A0',   # L=80  blue
    '#FF1493',   # L=100  hot pink
    '#00B300',   # L=120 pure green
    '#FFD700',   # L=140 vivid yellow
    '#FF0000',   # L=160 pure red
]
colormap = {l: c for l, c in zip(lattice_sizes, colors)}


#list of objects with associated values
observables = [
    ("energy density",           energy,                          energy_error,                        r"$\langle e \rangle$", False),
    ("energy variance",          energy_variance,                 energy_variance_error,               r"$C = L^2 (\langle e^2 \rangle - \langle e \rangle^2)$", True),
    ("absolute magnetization",   absolute_magnetization,          absolute_magnetization_error,        r"$\langle |m| \rangle$", False),
    ("magnetization variance",   absolute_magnetization_variance, absolute_magnetization_variance_error, r"$\chi' = L^2 (\langle m^2 \rangle - \langle |m| \rangle^2)$", True),
    ("magnetization squared",    magnetization2,                  magnetization2_error,                r"$\chi = L^2 (\langle m^2 \rangle)$", True),
    ("binder cumulant",          binder,                          binder_error,                        r"$U = \langle m^4 \rangle / \langle m^2 \rangle^2$", False),
]


#making the grid of the plots mosaic
n_obs  = len(observables)
n_cols = 2
n_rows = (n_obs + 1) // n_cols



#=============================
#   Plotting the mosaic
#=============================


fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, n_rows * 5))
axes = axes.flatten()


for index, (title, observable, error, label, rescale) in enumerate(observables):
    ax = axes[index]

    for l in lattice_sizes:
        mask = lattice_size == l
        scale = l**2 if rescale else 1.0
        ax.errorbar(
            beta[mask], observable[mask]*scale, yerr = error[mask]*scale,
            fmt=markermap[l],
            markerfacecolor = 'none',
            markeredgecolor = colormap[l],
            ecolor=colormap[l],
            capsize = 2,
            markersize = 5,
            zorder = l,
            label = f"L = {l}"
        )
    
    #binder cumulant special plots like the line and the insert zoomed
    if title == "binder cumulant":
        ax.axhline(y=1.1679227, color='black', linestyle='dotted', linewidth=0.8, label=r"$U^* = 1.1679$", alpha=0.7, zorder=0)
        axins = ax.inset_axes([0.67, 0.24, 0.3, 0.3])

        for l in lattice_sizes:
            mask = lattice_size == l
            axins.errorbar(
                beta[mask], binder[mask], yerr=binder_error[mask],
                fmt='-' + markermap[l],
                markerfacecolor='none',
                markeredgecolor=colormap[l],
                ecolor=colormap[l],
                color=colormap[l],
                capsize=2,
                markersize=4,
                linewidth=0.8,
                zorder=l
            )
        axins.axhline(y=1.1679227, color='black', linestyle='dotted', linewidth=0.8, alpha=0.7, zorder=0)
        axins.set_xlim(0.438, 0.444)
        axins.set_ylim(1.12, 1.23)
        axins.tick_params(labelsize=6)
    
    
    if title == "absolute magnetization":
        beta_critical    = 0.4406868
        beta_theory      = np.linspace(beta_critical, max(beta) + 0.01, 5000)
        m_theory         = (1 - np.sinh(2 * beta_theory)**(-4))**(1/8)
        ax.plot(beta_theory, m_theory, color='black', linestyle='dotted', label=r"$\langle |m| \rangle_\infty$ (Onsager)", alpha=0.7,  zorder=0)    

    
    ax.set_xlabel(r"$\beta$", fontsize=12)
    ax.set_ylabel(label, fontsize=12)
    ax.legend(fontsize=9, loc='best')


for idx in range(n_obs, len(axes)):
    axes[idx].set_visible(False)


#plotting and saving the result
plt.suptitle("2D Ising model observables", fontsize=15, y=1.01)
plt.tight_layout()


with PdfPages("results/observables_mosaic.pdf") as pdf:
    pdf.savefig(fig, bbox_inches="tight")

plt.close(fig)



#============================================================================================



#=============================
#     2D-ising critical
#        exponents  
#=============================


beta_c =  0.4406868
nu = 1.0
gamma = 7/4
b = 1/8



#=============================
#     values rescaled
#=============================


collapse_observables = [
    ("magnetization collapse",        absolute_magnetization,          absolute_magnetization_error,          r"$\langle |m| \rangle \cdot L^{\beta/\nu}$",  lambda obs, err, l: (obs * l**(b/nu),      err * l**(b/nu))),
    ("susceptibility collapse",       absolute_magnetization_variance, absolute_magnetization_variance_error, r"$\chi' \cdot L^{-\gamma/\nu}$",                    lambda obs, err, l: (obs * l**2 * l**(-gamma/nu), err * l**2 * l**(-gamma/nu))),
    ("chi collapse",                  magnetization2,                  magnetization2_error,                  r"$\chi \cdot L^{-\gamma/\nu}$",                     lambda obs, err, l: (obs * l**2 * l**(-gamma/nu), err * l**2 * l**(-gamma/nu))),
    ("binder collapse",               binder,                          binder_error,                          r"$U$",                                              lambda obs, err, l: (obs,                         err)),
]



n_obs  = len(collapse_observables)
n_cols = 2
n_rows = (n_obs + 1) // n_cols


#dictionary for the insets
# "Title": ([left, bottom, width, height], (x_min, x_max), (y_min, y_max))
inset_config = {
    "magnetization collapse":  ([0.65, 0.15, 0.3, 0.3], (-0.45, -0.14), (0.72, 0.84)), # Left plot
    "susceptibility collapse": ([0.65, 0.65, 0.3, 0.3], (-0.63, -0.18), (0.0975, 0.1150)),  # Right plot - zoom on peak
    "chi collapse":            ([0.65, 0.15, 0.3, 0.3], (-0.3, -0.12), (0.76, 0.94)),  # Left plot
    "binder collapse":         ([0.65, 0.65, 0.3, 0.3], (-0.35, -0.15), (1.275, 1.450)) # Right plot - central crossing
}



#=============================
#   plotting the FSS results
#=============================


fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, n_rows * 5))
axes = axes.flatten()


for index, (title, observable, error, label, rescale) in enumerate(collapse_observables):
    ax = axes[index]
    position, x_lim, y_lim = inset_config[title]
    axins = ax.inset_axes(position)

    for l in lattice_sizes:
        mask = lattice_size == l
        x = (beta[mask] - beta_c) * l**(1.0/nu)
        y, yerr = rescale(observable[mask], error[mask], l)

        ax.errorbar(
            x, y, yerr=yerr,
            fmt=markermap[l],
            markerfacecolor='none',
            markeredgecolor=colormap[l],
            ecolor=colormap[l],
            capsize=2,
            markersize=5,
            zorder=l,
            label=f"L = {l}"
        )
        
        axins.errorbar(
            x, y, yerr=yerr,
            fmt=markermap[l],         
            markerfacecolor='none',
            markeredgecolor=colormap[l],
            ecolor=colormap[l],
            capsize=2,
            markersize=5,
            zorder=l
        )

    axins.set_xlim(*x_lim)
    axins.set_ylim(*y_lim)
    axins.tick_params(labelsize=6)
    
    ax.set_xlabel(r"$(\beta - \beta_c) \cdot L^{1/\nu}$", fontsize=12)
    ax.set_ylabel(label, fontsize=12)
    if index % 2 == 1:
        ax.legend(fontsize=9, loc='lower left')
    else:
        ax.legend(fontsize=9, loc='upper left')

for idx in range(n_obs, len(axes)):
    axes[idx].set_visible(False)


plt.suptitle(
    rf"Finite size scaling collapse  ($\beta_c = {beta_c}$, $\nu = {nu}$, $\gamma = {gamma}$, $\beta = {b}$)",
    fontsize=13, y=1.01
)
plt.tight_layout()


 
with PdfPages("results/collapse.pdf") as pdf:
    pdf.savefig(fig, bbox_inches="tight")
 
plt.close(fig)



#============================================================================================



def chi_squared_red(y_data, y_fit, y_err, n_params):
    return np.sum(((y_data - y_fit) / y_err)**2) / (len(y_data) - n_params)



#=============================
#  parabolic fit for the
#  peaks of chi prime 
# obtaining:
# beta_pc, chi_prime_max
#=============================


#function to fit the peaks of chi' with a parabola X'(beta) = a(beta - beta_pc)^2 + X'_max
def parabola(beta, a, beta_pc, chi_prime_max):
    return a*(beta- beta_pc)**2 + chi_prime_max 


#using a dictionary to have n_points adjusted values for every single L
n_points_dict = {
    80:  15,
    100: 10,
    120: 11,
    140: 13,
    160: 13,
}

#arrays to save the results of the fits: X^2, beta_pc +- err, chi'max +- err
chi_prime_max_list = []
chi_prime_max_error_list = []
beta_pc_list = []
a_list = []
beta_pc_error_list = []
X2_list = []


#figure of the parabolic fits 
fig, ax = plt.subplots(figsize=(10, 6))


#taking only the lattice sizes with the more prominent peaks
L_fit = [l for l in lattice_sizes if l >= 80]


#extracting only the wanted values from the data bytearray
for l in L_fit:
    mask = lattice_size == l
    x = beta[mask]
    y = absolute_magnetization_variance[mask] * l**2
    yerr = absolute_magnetization_variance_error[mask] * l**2

    n_points = n_points_dict[l]
    peak_index = np.argmax(y)
    starting_index = max(0, peak_index - n_points//2)
    ending_index = min(len(x), starting_index + n_points)

    x_cut = x[starting_index:ending_index]
    y_cut = y[starting_index:ending_index]
    yerr_cut = yerr[starting_index:ending_index]

    popt,pcov = curve_fit(parabola, x_cut, y_cut, sigma= yerr_cut, absolute_sigma=True, p0=[-abs(max(y_cut)), x_cut[np.argmax(y_cut)], max(y_cut)])
    X2_red = chi_squared_red(y_cut, parabola(x_cut, *popt), yerr_cut, 3)

    #appending the fits results, used later for analysis
    beta_pc_list.append(popt[1])
    beta_pc_error_list.append(np.sqrt(pcov[1,1]))
    chi_prime_max_list.append(popt[2])
    chi_prime_max_error_list.append(np.sqrt(pcov[2,2]))
    a_list.append(popt[0])
    X2_list.append(X2_red)

    
    #plotting the fits results
    ax.errorbar(
        x_cut, y_cut, yerr=yerr_cut,
        fmt=markermap[l],
        markerfacecolor='none',
        markeredgecolor=colormap[l],
        ecolor=colormap[l],
        capsize=2,
        markersize=5,
        zorder=l,
        label=f"L = {l}"
    )
    
    x_fit = np.linspace(min(x_cut), max(x_cut), 500)
    ax.plot(x_fit, parabola(x_fit, *popt), color='black', linewidth=1, zorder=0)
    

    #printing in terminal    
    print(f"\n{'='*50}")
    print(f"  L = {l}")
    print(f"  n_points  = {len(x_cut)}")
    print(f"  beta_pc   = {popt[1]:.6f} +/- {np.sqrt(pcov[1,1]):.6f}")
    print(f"  chi'_max  = {popt[2]:.4f} +/- {np.sqrt(pcov[2,2]):.4f}")
    print(f"  X2_red    = {X2_red:.4f}")
    print(f"{'='*50}")


ax.set_xlabel(r"$\beta$", fontsize=12)
ax.set_ylabel(r"$\chi' = L^2(\langle m^2 \rangle - \langle |m| \rangle^2)$", fontsize=12)
ax.legend(fontsize=9, loc='best')
plt.tight_layout()


with PdfPages("results/parabolic_fits.pdf") as pdf:
    pdf.savefig(fig, bbox_inches="tight")

plt.close(fig)



#=============================
# Fitting the values obtained
#  with the scaling law 
#  for beta_pc and 
#  for chi_prime_max
#=============================


def scaling_law_chi_prime_max(L, c_0, c_1 , g_over_nu):
    return c_0 + c_1 * (L**(g_over_nu))

def scaling_law_beta_pc(L, beta_c, b, nu):
    return beta_c + b * (L**(-1/nu))



#from lists to arrays
beta_pc_fit             = np.array(beta_pc_list)
beta_pc_error_fit       = np.array(beta_pc_error_list)
chi_prime_max_fit       = np.array(chi_prime_max_list)
chi_prime_max_error_fit = np.array(chi_prime_max_error_list)
L_fit                   = np.array(L_fit)


#fits
popt_beta, pcov_beta = curve_fit(scaling_law_beta_pc, L_fit, beta_pc_fit, sigma=beta_pc_error_fit, absolute_sigma=True, p0=[0.4407, -0.5, 1.0])

beta_c_inf = popt_beta[0]
b_fit      = popt_beta[1]
nu_fit     = popt_beta[2]

beta_c_inf_err = np.sqrt(pcov_beta[0,0])
b_fit_err      = np.sqrt(pcov_beta[1,1])
nu_fit_err     = np.sqrt(pcov_beta[2,2])


popt_chi, pcov_chi = curve_fit(scaling_law_chi_prime_max, L_fit, chi_prime_max_fit ,sigma=chi_prime_max_error_fit, absolute_sigma=True, p0=[0.0, 1.0, 1.75])

c_0      = popt_chi[0]
c_1      = popt_chi[1]
g_nu     = popt_chi[2]

c_0_err  = np.sqrt(pcov_chi[0,0])
c_1_err  = np.sqrt(pcov_chi[1,1])
g_nu_err = np.sqrt(pcov_chi[2,2]) 


# --- χ²_red for both fits ---
y_fit_beta = scaling_law_beta_pc(L_fit, *popt_beta)
chi2_beta  = chi_squared_red(beta_pc_fit, y_fit_beta, beta_pc_error_fit, 3)

y_fit_chi  = scaling_law_chi_prime_max(L_fit, *popt_chi)
chi2_chi   = chi_squared_red(chi_prime_max_fit, y_fit_chi, chi_prime_max_error_fit, 3)


# --- print ALL parameters ---
print(f"\n{'='*75}")
print("          EXTRAPOLATION RESULTS — ALL PARAMETERS")
print(f"  β_pc fit:")
print(f"     β_c     = {beta_c_inf:.8f} ± {beta_c_inf_err:.8f}")
print(f"     b       = {b_fit:.6f} ± {b_fit_err:.6f}")
print(f"     ν       = {nu_fit:.4f} ± {nu_fit_err:.4f}    (χ²_red = {chi2_beta:.3f})")
print(f"  χ'_max fit:")
print(f"     c_0     = {c_0:.6f} ± {c_0_err:.6f}")
print(f"     c_1     = {c_1:.6f} ± {c_1_err:.6f}")
print(f"     γ/ν     = {g_nu:.4f} ± {g_nu_err:.4f}    (χ²_red = {chi2_chi:.3f})")
print(f"{'='*75}")



#=============================
# plot of the scaling law fits
#=============================


fig, axes = plt.subplots(1, 2, figsize=(14, 6))

L_smooth = np.linspace(min(L_fit), max(L_fit), 600)

# beta_pc vs L
axes[0].errorbar(
    L_fit, beta_pc_fit, yerr=beta_pc_error_fit,
    fmt='D',
    markerfacecolor = 'none',
    markeredgecolor = 'red',
    ecolor='red',
    capsize=2,
    markersize=5,
    zorder=3,
    label=r"$\beta_{pc}(L)$"
)
axes[0].plot(L_smooth, scaling_law_beta_pc(L_smooth, *popt_beta),
             color='black', linestyle='-', linewidth=1, zorder=1,
             label=r"$\beta_c + b \cdot L^{-1/\nu}$")
axes[0].set_xlabel(r"$L$", fontsize=12)
axes[0].set_ylabel(r"$\beta_{pc}(L)$", fontsize=12)
axes[0].legend(fontsize=9, loc='best')


# chi_prime_max vs L
axes[1].errorbar(
    L_fit, chi_prime_max_fit, yerr=chi_prime_max_error_fit,
    fmt='D',
    markerfacecolor='none',
    markeredgecolor = 'red',
    ecolor='red',
    capsize=2,
    markersize=5,
    zorder=3,
    label=r"$\chi'_{max}(L)$"
)
axes[1].plot(L_smooth, scaling_law_chi_prime_max(L_smooth, *popt_chi),
             color='black', linestyle='-', linewidth=1, zorder=1,
             label=r"$c_0 + c_1 \cdot L^{\gamma/\nu}$")
axes[1].set_xlabel(r"$L$", fontsize=12)
axes[1].set_ylabel(r"$\chi'_{max}(L)$", fontsize=12)
axes[1].legend(fontsize=9, loc='best')


plt.suptitle("Scaling law fits", fontsize=13, y=1.01)
plt.tight_layout()

with PdfPages("results/scaling_fits.pdf") as pdf:
    pdf.savefig(fig, bbox_inches="tight")

plt.close(fig)
