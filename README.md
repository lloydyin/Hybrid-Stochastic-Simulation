# Hybrid Stochastic Simulation
This repository contains the simulation code and results for a hybrid stochastic algorithm developed for FRP and DT polymerzation. All components are modular and can be used independently.

## Repository Contents
- `Julia_Notebook.ipynb/html/pdf`
A Jupyter notebook containing the complete simulation procedure, including algorithm implementation, simulation steps, and result generation.

- `Simulation_Code.jl`
The core **Julia** script that implements the hybrid stochastic simulation algorithm.

- `mydata.jld2`
Simulation output data saved in Julia's JLD2 format. This file includes the generated data from simulation runs.

- `Codes_for_Visualization.jl`
Contains Julia code snippets for visualizing the simulation results in `mydata.jld2`. You can use these scripts to reproduce figures from the simulation outputs.

- `Plots/`
This folder contains the generated plots and figures exported from the simulations.

- `Benchmark.txt`
This folder includes results from repeated simulation runs (10 times) under default parameter settings, useful for performance comparison and consistency checking.

- `Execution_time.txt`
Contains timing data from different simulation methods and settings.

- `Execution_time.R`
**R** script used to generate performance comparison plots from Execution_time.txt.

## How to Use
- To view or run the simulation code, use either:
  - `Julia_Notebook.ipynb` (interactive, Jupyter);
  - `Simulation_Code.jl` (script-based, Julia REPL).

- To visualize simulation results, run:
  - `Codes_for_Visualization.jl` using **Julia**. It loads `mydata.jld2` and produces plots similar to those in the Plots/ folder.

- To analyze performance, you can:
  - Read `Benchmark.txt` for default-parameter benchmarks;
  - Run `Execution_time.R` in **R** to generate execution time comparison plots from `Execution_time.txt`.

## Data Availability
This code and all simulation results are made publicly available in accordance with **MDPI**'s open data policy. See the Data Availability Statement for more information.
