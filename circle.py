# This code is born with the purpose of controlling "metro_circle.c".
# More specifically here we we modify the "parameters.h" file in the header.
#
#
#

#======================================
#		Libraries
#======================================


import math
import os
import time
import numpy as np


#======================================
#		Writing into parameters.h
#======================================

open("data/results.txt", "w").close() #clears the file results if i want to run sims again

DELTA = 0.5

N_TERMAL = 50000   #number of sweeps  to do before data acquisition

#N_t_list = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]		#number of sites in the euclidean_time_lattice 
N_t_list = [325]

MINIMAL_BLOCK_SIZE = 50

#number of times we call the sweep function for a single lattice spacing
N_SWEEPS = 500000000  

for N_t in N_t_list:

  N_SAMPLES = N_t		#number of updates in a sweep    

  ETA = 2 / N_t
  PREFACTOR = 1 / (2 * ETA)
  
  with open("parameters.h", "w") as f:
	  f.write(f"#define MINIMAL_BLOCK_SIZE {MINIMAL_BLOCK_SIZE}\n")
	  f.write(f"#define N_SAMPLES {N_SAMPLES}\n")
	  f.write(f"#define N_TERMAL {N_TERMAL}\n")
	  f.write(f"#define N_t {N_t}\n")
	  f.write(f"#define N_SWEEPS {N_SWEEPS}\n")
	  f.write(f"#define DELTA {DELTA:.4f}\n")
	  f.write(f"#define PREFACTOR {PREFACTOR:.15f}\n")
	  f.write(f"#define ETA {ETA:.15f}\n")
  
  
  #======================================
  #       Running simulations		
  #======================================


  print(f" ---- Processing N_t = {N_t} ----\n")

  start_time = time.time()
  
  os.system("gcc -O3 metro_circle.c /usr/local/lib/pcg32.c -o circle.exe -lm ")
  os.system("./circle.exe")

  end_time = time.time()
  seconds = end_time - start_time
  seconds = (end_time - start_time) %  60
  minutes = ((end_time - start_time) - seconds) / 60
  print(f"time: {minutes:.0f} min  {seconds:.0f} sec \n")
    
  print("Cool down...\n")
  time.sleep(10)      #seconds of pause to cool down CPU
  
print("All simulations finished\n")



