================================================================================
2D Ising Model Cluster Simulation - README
================================================================================

This repository contains a complete suite of programs designed to simulate the 
2D Ising model using a Monte Carlo cluster update algorithm. The suite is 
designed to manage multiple simulation runs across different lattice sizes and 
temperatures, and subsequently analyze the resulting data to extract critical 
exponents via finite-size scaling.

The project is divided into three main files, each handling a distinct part of 
the workflow: the simulation manager, the computational engine, and the data 
analyzer.


--------------------------------------------------------------------------------
System Architecture: How the Codes Interact
--------------------------------------------------------------------------------

The three scripts form an automated pipeline that takes the user from parameter 
generation all the way to finalized scientific plots. 

1. The Orchestration:
   The workflow begins with `cluster.py`. This script is the manager. It decides 
   which lattice sizes to simulate and calculates the specific temperatures 
   (inverse temperature, beta) to test. Instead of doing the heavy math itself, 
   it acts as a compiler and execution manager. It dynamically writes the chosen 
   parameters into a C header file (`parameters.h`), compiles the C program, and 
   runs the resulting executable.

2. The Computation:
   Once called by the Python manager, the compiled `cluster.c` engine takes over. 
   It reads the parameters passed to it, initializes a random lattice, and runs 
   the Monte Carlo simulation. After completing its sweeps and calculating the 
   statistical errors, it appends a single line of summary data (averages and 
   errors for that specific lattice size and temperature) into a text file 
   (`results/results.txt`). Control is then handed back to the Python manager 
   for the next run.

3. The Analysis:
   After `cluster.py` has finished running all parameter combinations and the 
   data file is fully populated, the final script, `analysis_cluster.py`, is 
   used. This script acts completely independently of the simulation phase. It 
   reads the aggregated data from the text file, performs complex statistical 
   fitting, and generates visual representations of the physics, outputting 
   everything as PDF files.


--------------------------------------------------------------------------------
Theory: The Simulation Engine (cluster.c)
--------------------------------------------------------------------------------

The core of this project is simulating the 2D Ising model, a fundamental 
mathematical model of ferromagnetism in statistical mechanics. The system 
undergoes a continuous phase transition at a critical temperature.

The Cluster Update Algorithm:
Standard Monte Carlo methods (like the Metropolis algorithm) flip single spins 
one at a time. Near a phase transition, these methods suffer from "critical 
slowing down," where the system takes an impractically long time to evolve 
because large regions of spins become highly correlated. 

To overcome this, `cluster.c` implements a cluster update algorithm. 
- The algorithm picks a random seed spin and flips it.
- It then looks at neighboring spins that have the same original alignment.
- It adds these neighbors to a stack with a specific probability defined by 
  P_add = 1 - exp(-2 * beta).
- The algorithm proceeds to check the neighbors of the newly added spins, adding 
  them with the same probability. 
- This continues iteratively until no more spins are added. As a result, entire 
  domains of correlated spins are flipped simultaneously, vastly improving the 
  efficiency of the simulation near the critical point.

Observables:
During the simulation, the program continuously measures several physical 
properties of the lattice:
- Energy Density: The sum of nearest-neighbor interactions.
- Magnetization: The average spin direction across the entire lattice.
- Susceptibility: Measured in two ways, using the standard variance of the 
  magnetization and the subtracted absolute variance.
- Binder Cumulant: A specialized statistical measure involving the fourth moment 
  of the magnetization. The Binder cumulant is particularly useful because its 
  value at the critical point is theoretically independent of the lattice size, 
  making it an excellent tool for locating the critical temperature.

Error Analysis:
Monte Carlo data is inherently correlated (the state of the lattice at step N 
depends heavily on step N-1). To provide accurate error bars, the code uses 
advanced statistical techniques:
- Blocking Method: Used for direct observables to group data and find the true 
  standard deviation.
- Jackknife Method: Used for composite, non-linear observables (like variance 
  and the Binder cumulant) to remove bias and estimate correct statistical errors.


--------------------------------------------------------------------------------
File Breakdown: Detailed Explanations
--------------------------------------------------------------------------------

1. cluster.c (The Simulation Core)
----------------------------------
This C program is built for performance and handles the intense mathematical 
operations. 

- Initialization: It sets up a random grid of spins using the highly efficient 
  `pcg32` random number generator. It first runs a "thermalization" phase 
  (e.g., 10,000 sweeps) where no data is collected, ensuring the system reaches 
  thermal equilibrium before measurements begin.
- Sweeping and Acquisition: The code executes a set number of cluster updates. 
  Crucially, it skips a specific number of sweeps between actual measurements to 
  further guarantee that the recorded data points are statistically independent.
- Output: Once all sweeps are finished, it computes the final averages and the 
  Jackknife/blocking errors for every observable, formatting them into a neat, 
  comma-separated string appended to the results file.


2. cluster.py (The Simulation Manager)
--------------------------------------
This Python script prevents the user from having to manually compile and run 
the C code hundreds of times.

- Targeted Parameter Generation: It loops through a predefined list of lattice 
  sizes (e.g., L=20 up to L=160). For each size, it calculates an optimal window 
  of beta values. Instead of scanning the whole temperature range uniformly, it 
  uses theoretical estimations to zoom in precisely on the critical region, 
  which narrows as the lattice gets larger.
- Automation: For every temperature point, it creates a custom `parameters.h` 
  file. It then invokes the system shell to compile `cluster.c` with the `-O3` 
  optimization flag, ensuring maximum execution speed.
- Benchmarking: To help the user monitor progress, the script tracks the exact 
  time taken for each individual run and prints a formatted log to the terminal.


3. analysis_cluster.py (Data Analysis & Visualization)
------------------------------------------------------
This Python script handles all post-processing, transforming raw numbers into 
meaningful physics.

- Data Extraction: It reads the comma-separated results file into organized 
  NumPy arrays, separating the data by lattice size.
- Raw Observables Plotting: It creates a comprehensive mosaic plot showing how 
  energy, magnetization, specific heat, susceptibility, and the Binder cumulant 
  behave as a function of temperature. It even overlays theoretical exact 
  solutions (like the Onsager solution for infinite lattices) for comparison.
- Finite Size Scaling (FSS): The script applies theoretical critical exponents 
  to the axes of the plots. If the physics is correct, data from completely 
  different lattice sizes will magically collapse onto a single universal curve.
- Peak Fitting & Extrapolation: It automatically locates the peak of the 
  susceptibility curve for each lattice size by fitting a parabola. It then 
  plots these pseudo-critical temperatures against the lattice size, applying a 
  scaling law fit to extrapolate the true critical temperature for an infinite 
  lattice, printing the final extracted critical exponents.

