# This script automates the execution of the C simulation 'free_field.c' for various
# lattice temporal sizes (Nt = 4,6,8,10), spatial sizes (Ns = Nt * 5), and lattice masses (mhat).
# It performs a linear scan in mhat from 2/Ns to 1 with 20 points per Nt, 
# running simulations for both finite temperature (Nt, Ns) and zero temperature (Ns, Ns)
# cases, collects the averages and errors of observables O1, O2, O3 from the output files,
# and saves all raw data to 'data1/data.txt' for further analysis or plotting.



#======================================
#		Libraries
#======================================


import os
import time
import math
import numpy as np



#======================================
#		Setting up the parameters
#======================================


#N_T_values = [4, 6, 8, 10]
N_T_values = [4, 6, 8]

ratio = 5

#number of points in the plot for every N_T value
M_HAT_points = 15

simulation_queue = []

for N_T in N_T_values:
  
  N_S = N_T * ratio
  simulation_queue.append((N_T, N_S))


#======================================
#  Function writing parameters.h 
#       and run sims
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
  time.sleep(2)



#======================================
#       Running simulations		
#======================================


counter = 0
print("Starting simulations...\n")

for N_T, N_S in simulation_queue:

  # generating points that will be linear in T/m
  T_over_m_min = 0.2
  T_over_m_max = 2.5
  T_over_m_list = np.linspace(T_over_m_min, T_over_m_max, M_HAT_points)

  # Calculate the lattice mass needed for each point
  # Formula: M_HAT = 1 / (N_T * (T/m))
  M_HAT_list = 1.0 / (N_T * T_over_m_list)

  counter = counter + 1
  print(f"\n\n\n ------------ SCAN NUMBER: {counter} ------------ \n")
  
  for M_HAT in M_HAT_list:
    #finite temperature simulation with scan in M_HAT
    print(f" ---- Finite temperature: N_T = {N_T}, N_S = {N_S}, M_HAT = {M_HAT:.4f} ----\n ")
    run_simulation(N_T, N_S, M_HAT)

    #zero temperature simulation with scan in M_HAT
    print(f" ---- Zero temperature: N_T = {N_S}, N_S = {N_S}, M_HAT = {M_HAT:.4f} ----\n")
    run_simulation(N_S, N_S, M_HAT)

  print("\n\n ------------ SCAN FINISHED ------------ ")


"""
print(f"\n\n ------------ FINITE SIZE EFFECT SIMULATION  ------------\n\n")
N_S_check = 100
N_T_check = 10

for M_HAT in M_HAT_list:
  
  print(f" ---- finite temperature check: N_T = {N_S_check}, N_S = {N_S_check}, M_HAT = {M_HAT:.4f} ----\n\n")
  run_simulation(N_T, N_S, M_HAT)
  
  print(f" ---- zero temperature check: N_T = {N_S_check}, N_S = {N_S_check}, M_HAT = {M_HAT:.4f} ----\n\n")
  run_simulation(N_S, N_S, M_HAT)
  
print("Simulations finished...\n")

"""
