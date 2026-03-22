#
#
# This code is born with the purpose of controlling "metro_circle.c"
# to study the High-Temperature regime of the quantum particle on a circle.
# More specifically here we modify the "parameters.h" file in the header.
#
#
# eta = (beta h^2) / m(2 pi R)^2  x (1/N_t)
#
#
#The idea is to keep fixed the value of N_t and scan in the temperature (so in the value of eta), evaluate the winding number and see if goes
#as expected with a fit. 
#


# b' = b h^2 / m(2 pi R)^2
#Codntion for low temperature is b h^2/ (m R^2)<<1 so   4 pi^2 b' << 1  
#



#======================================
#	    	 Libraries
#======================================


import math
import os
import time
import numpy as np



#======================================
#     Writing into parameters.h
#======================================


#max size of the random jump in the circle
DELTA = 0.5
#number of sweeps we do before data acquisition
N_THERMAL = 50000
#fix here the value of N_t to high(continuous limit) 
N_t = 20
#number of single updates in a sweep
N_SAMPLES = N_t
#number of times we sweep trough the lattice: from 5x10^8 to 10^8 since at high T freezing is milder
N_SWEEPS = 300000000
#minimal block size aceptable
MINIMAL_BLOCK_SIZE = 50


#defining a list of rescaled inverse temperatures:
# b' = b h^2 / m(2 pi R)^2
#BETA_PRIME_list = np.linspace(2.0, 0.35, 22)
BETA_PRIME_list = [0.30, 0.27, 0.25, 0.23, 0.22]

for beta_prime in BETA_PRIME_list:
    ETA = beta_prime/N_t
    PREFACTOR = 1/(2*ETA)
    
    with open("parameters.h", "w") as f:
        f.write(f"#define MINIMAL_BLOCK_SIZE {MINIMAL_BLOCK_SIZE}\n")
        f.write(f"#define N_SAMPLES {N_SAMPLES}\n")
        f.write(f"#define N_THERMAL {N_THERMAL}\n")
        f.write(f"#define N_t {N_t}\n")
        f.write(f"#define N_SWEEPS {N_SWEEPS}\n")
        f.write(f"#define DELTA {DELTA:.4f}\n")
        f.write(f"#define PREFACTOR {PREFACTOR:.15f}\n")
        f.write(f"#define BETA {beta_prime:.5f}\n")
    
    
    
    
    #======================================
    #       running simulations
    #======================================
    
    print("="*50)
    print(f" BETA' = {beta_prime:.6f} (ETA = {ETA:.6e}, PREFACTOR = {PREFACTOR:.1f})\n")
    
    start_time = time.time()
    
    os.system("gcc -O3 metro_circle.c /usr/local/lib/pcg32.c -o circle_high_temp.exe -lm ")
    os.system("./circle_high_temp.exe")
    
    end_time = time.time()
    total = end_time - start_time
    minutes = int(total // 60)
    seconds = total % 60
    print(f"time: {minutes:.0f} min  {seconds:.0f} sec")
    print("="*50)

print("All simulations finished\n")

