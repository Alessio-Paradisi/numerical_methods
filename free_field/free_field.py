# This code controls "free_field.c" giving it in input the following values: 
# N_T, N_S, M_HAT, N_T_BAR
#
# Then runs sims for:
# N_T values: [4, 5, 6, 7, 8, 10]
#
# For each N_T, N_S takes values:
#  - 4 * N_T
#  - 6 * N_T  
#  - 8 * N_T
#
# For each (N_T, N_S) pair, two simulations are run:
#  1. Finite temperature:    N_T, N_S, M_HAT
#  2. Zero temperature ref:  N_T_bar = N_S, N_S, same M_HAT
#
#======================================
#		Libraries
#======================================


import math
import os
import time
import numpy as np


#======================================
#		Setting up the parameters
#======================================


N_T_values = [4, 5, 6, 7, 8, 10]

ratios = [4, 6, 8]

simulation_queue = []

for N_T in N_T_values:

  for ratio in ratios:
  
    N_S = N_T * ratio
    simulation_queue.append((N_T, N_S))

#======================================
#	Function to run the given simulation
#======================================


def run_simulation(N_T, N_S, M_HAT):
 
  with open("parameters.h", "w") as f:
    f.write(f"#define N_T {N_T}\n")
    f.write(f"#define N_S {N_S}\n")
    f.write(f"#define M_HAT {M_HAT:.10f}\n")

  start_time = time.time()
  
  os.system("gcc -O3 free_field.c /usr/local/lib/pcg32.c -o free_field.exe -lm")
  os.system("./free_field.exe")
  
  end_time = time.time()
  seconds = end_time - start_time
  seconds_disp = seconds % 60
  minutes_disp = (seconds - seconds_disp) / 60
      
  print(f"time: {minutes_disp:.0f} min  {seconds_disp:.0f} sec \n")
  
  print("Cool down...\n")
  time.sleep(6)


#======================================
#       Running simulations		
#======================================


print("Starting simulations...\n")

counter = 0

for N_T, N_S in simulation_queue:

  M_HAT = 1.0 / N_T
  counter = counter + 1

  print(f"\n \n \n----------- START SIMULATION NUMBER: {counter} ----------- \n" )
  
  # finite temperature simulation
  print(f" ---- Finite temperature: N_T = {N_T}, N_S = {N_S}, M_HAT = {M_HAT:.4f} ----\n")
  run_simulation(N_T, N_S, M_HAT)

  # zero temperature reference simulation with same N_S and M_HAT, N_T_bar = N_S
  N_T_bar = N_S
  print(f" ---- Zero temperature ref: N_T = {N_T_bar}, N_S = {N_S}, M_HAT = {M_HAT:.4f} ----\n")
  run_simulation(N_T_bar, N_S, M_HAT)

  print(f"----------- END SIMULATION  -----------\n \n \n" )

print("All simulations finished\n")

  
