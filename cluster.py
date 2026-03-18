

#=============================
#      libraries
#=============================


import numpy as np
import os
import time



#=============================
#   Simulation parameters
#          setup      
#=============================


L = [20, 40, 60, 80, 100, 120, 140, 160]
beta_points = 40
beta_critical = 0.4406867935

#values that enter the function to estimate the pseudo-critical beta, center of the temperature scan (ideally the peak)
nu = 1.0     
a = -0.41     

#half of the window of scan for the minimum and maximum value of L
window_min_L = 0.011   
window_max_L = 0.004



#=============================
# functions to generate beta
#=============================


def beta_pseudocritical(l):
    return beta_critical + a * l**(-1.0/nu)


def window_halfwidth(l): 
    t = (l - L[0])/(L[-1] - L[0])
    return window_min_L + t * (window_max_L - window_min_L)


def beta_values(l):
    beta_pc = beta_pseudocritical(l)
    if L.index(l) < 3:
        return np.linspace(0.42, 0.46, beta_points)
    else:
        half_width = window_halfwidth(l)
        return np.linspace(beta_pc - half_width, beta_pc + half_width, beta_points)


#passes to the c code N_CLUSTER= number of cluster updates x sweep, N_SWEEPS=total number of cycles done, N_SKIP=every n_skip we sample the base obs
def get_simulation_params(l):
    params = {
        20:  (17, 58937, 3),   # 58937*17 = 1,001,929  skip=3
        40:  (17, 58937, 3),   # 1,001,929             skip=3
        60:  (17, 64000, 3),   # 1,088,000             skip=3
        80:  (17, 70000, 3),   # 1,190,000             skip=3
        100: (17, 76000, 1),   # 1,292,000             skip=2
        120: (17, 82000, 1),   # 1,394,000             skip=2
        140: (17, 88000, 1),   # 1,496,000             skip=2
        160: (17, 94000, 1),   # 1,598,000             skip=2
    }
    return params[l]



#=============================
# function for the simulation
#=============================


def compile_and_run(l, beta, index):
    n_cluster, n_sweeps, n_skip = get_simulation_params(l)
    
    with open("parameters.h", "w") as f:
        f.write(f"#define L {l}\n")
        f.write(f"#define BETA {beta:.10f}\n")
        f.write(f"#define N_CLUSTER {n_cluster}\n")
        f.write(f"#define N_SWEEPS {n_sweeps}\n")
        f.write(f"#define N_SKIP {n_skip}\n")
    
    
    start_time = time.time()
    
    os.system(f"gcc -O3 cluster.c /usr/local/lib/pcg32.c -o cluster.exe -lm")
    os.system(f"./cluster.exe")
    
    end_time = time.time()
    seconds = end_time - start_time
    minutes = int(seconds / 60)
    seconds = seconds % 60
    
    print(f"{index})  beta = {beta:.6f} time: {minutes:.0f} min  {seconds:.0f} sec\n")



#=============================
# running the simulations
#     with a cycle
#=============================


for l in L:
    betas = beta_values(l)
    beta_pc = beta_pseudocritical(l)

    print(f"\n{'='*50}")
    print(f"  SIMULATION:  L = {l}")
    print(f"{'='*50}")
    print(f"  beta_pc  = {beta_pc:.6f}")
    print(f"  beta_min = {betas[0]:.6f}")
    print(f"  beta_max = {betas[-1]:.6f}")
    print(f"{'='*50}")


    scan_start = time.time()
    
    for index, beta in enumerate(betas, start=1):
        compile_and_run(l,beta, index)

    scan_end     = time.time()
    scan_seconds = scan_end - scan_start
    scan_minutes = int(scan_seconds / 60)
    scan_seconds = scan_seconds % 60

    print(f"\n  total scan time for L={l}: {scan_minutes:.0f} min  {scan_seconds:.0f} sec")
    print(f"\n{'='*50}")

print(f"\n{'='*50}")
print("SIMULATIONS ENDED")
print(f"{'='*50}\n")

