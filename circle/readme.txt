# Quantum Particle on a Circle — Metropolis Simulation

Path-integral Monte Carlo simulation of a quantum particle constrained to move
on a circle, using the Metropolis algorithm on a euclidean time lattice.
The observable of interest is the winding number squared <Q^2>, which is
measured as a function of the lattice spacing eta = 2/N_t and extrapolated to
the continuum limit eta -> 0.

---

## Metropolis Update — How sites are updated

The euclidean time lattice has N_t sites {x_0, ..., x_{N_t-1}}, each holding
the position of the particle (a number in [0, 1)) at that imaginary time slice.
Periodic boundary conditions are imposed: site N_t wraps back to site 0.

Each call to site_update() does the following:

1. SITE SELECTION: a random site index i is chosen uniformly in [0, N_t).

2. PROPOSAL: a new position is proposed by adding a uniform random displacement
   in [-DELTA, +DELTA] to the current value. If the result falls outside [0, 1)
   it is wrapped back with modular arithmetic, preserving the circular topology.

3. ACTION DIFFERENCE: the change in the discretised euclidean action is computed
   from only the two bonds connecting site i to its neighbours (i-1) and (i+1):

       delta_S = PREFACTOR * [ d(x_new, x_{i+1})^2 + d(x_{i-1}, x_new)^2
                              - d(x_old, x_{i+1})^2 - d(x_{i-1}, x_old)^2 ]

   where PREFACTOR = 1 / (2 * eta) and d(x, y) is the shortest-arc distance on
   the circle, wrapping differences that exceed 0.5 back into [-0.5, 0.5].

4. ACCEPT / REJECT: the move is accepted if delta_S <= 0 (lower action), or
   otherwise with probability exp(-delta_S), following the standard Metropolis
   criterion. If rejected, the lattice is unchanged.

A "sweep" (lattice_sweep) calls site_update() N_SAMPLES = N_t times, giving
on average one proposed update per site per sweep.

---

## Winding Number

After each sweep the winding number Q is computed by summing the signed
arc-distances along the full time loop:

    Q = sum_{i=0}^{N_t - 1}  d(x_i, x_{i+1 mod N_t})

Q counts how many times the particle winds around the circle over the full
euclidean time period. The quantity of physical interest is <Q^2>, the ensemble
average of the squared winding number.

---

## Files

### `metro_circle.c` — C simulation

Runs the Metropolis algorithm on an N_t-site euclidean time lattice. All
parameters are read from `parameters.h`, which is written by `circle.py` before
each compilation. Key parameters:

| Parameter        | Description                                      |
|------------------|--------------------------------------------------|
| N_t              | Number of euclidean time sites                   |
| ETA              | Lattice spacing  (= 2 / N_t)                     |
| PREFACTOR        | Action prefactor (= 1 / (2 * eta))               |
| DELTA            | Half-width of Metropolis proposal window (= 0.5) |
| N_SWEEPS         | Total number of sweeps (measurements)            |
| N_TERMAL         | Thermalisation sweeps before data acquisition    |
| N_SAMPLES        | Site updates per sweep (= N_t)                   |
| MINIMAL_BLOCK_SIZE | Minimum block size for plateau search          |

After thermalisation the code runs N_SWEEPS sweeps, writing one Q^2 value per
sweep to a binary file (data/measurements.bin).

Error estimation uses a jackknife-over-blocks method. The optimal block size is
found automatically by scanning increasing block sizes and detecting a plateau
in the jackknife error using a sliding window of size WINDOW_SIZE=8. The
coefficient of variation (std/mean) of the window must fall below a threshold
(starting at 0.001, relaxed up to 0.02 if no plateau is found). The plateau
curve is saved to data/plateau_{N_t}.txt. Final results are appended to
data/results.txt.

**Compile and run (handled automatically by circle.py):**
    gcc -O3 metro_circle.c /usr/local/lib/pcg32.c -o circle.exe -lm
    ./circle.exe

---

### `circle.py` — Simulation controller

Loops over a list of N_t values, writing parameters.h for each one, then
compiling and running metro_circle.c. A 10-second cooldown pause is inserted
between runs.

Key settings to edit at the top of the script:

| Variable         | Default     | Description                              |
|------------------|-------------|------------------------------------------|
| N_t_list         | [325]       | List of euclidean time lattice sizes     |
| DELTA            | 0.5         | Metropolis proposal half-width           |
| N_TERMAL         | 50000       | Thermalisation sweeps                    |
| N_SWEEPS         | 500000000   | Measurement sweeps per N_t              |
| MINIMAL_BLOCK_SIZE | 50        | Minimum block size for plateau detection |

ETA and PREFACTOR are computed automatically from N_t. Results are accumulated
in data/results.txt (cleared at the start of each full run).

**Run:**
    python circle.py

---

### `circle_visualizer.py` — Plots

Reads the output files and produces two figures saved as PDFs:

data/plateau.pdf
    Mosaic of jackknife error vs block size for each N_t in N_t_list. Each
    subplot shows the plateau used to select the optimal block size.

data/results_analysis.pdf
    Left panel: <Q^2> vs eta^2 with error bars. A dashed line marks the exact
    continuum value <Q^2> = 2. An inset zooms into the small-eta region.
    Right panel: statistical error sigma_{<Q^2>} vs N_t on a log scale.

**Run:**
    python circle_visualizer.py

---

## Expected directory structure

    .
    ├── metro_circle.c
    ├── circle.py
    ├── circle_visualizer.py
    ├── parameters.h          # auto-generated by circle.py, do not edit manually
    └── data/
        ├── results.txt       # N_t, <Q^2>, sigma per line
        ├── plateau_{N_t}.txt # block size vs jackknife error for each N_t
        └── measurements.bin  # temporary binary file, deleted after each run

---

## Quick workflow

    # 1. Edit N_t_list and run parameters in circle.py, then:
    python circle.py

    # 2. Once all N_t runs are complete, visualise:
    python circle_visualizer.py
