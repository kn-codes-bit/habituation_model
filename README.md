# Mathematical Modeling for a Primitive Form of Habituation in an Amoeba

[![Paper](https://img.shields.io/badge/Paper-DOI-blue)](https://doi.org/10.1007/s00285-026-02438-7)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281/zenodo.20438610-blue)](https://doi.org/10.5281/zenodo.20438610)

This repository contains the official C++ implementation of the numerical simulations and the Python scripts for parameter estimation in our manuscript:  
> **"Mathematical Modeling for a Primitive Form of Habituation in an Amoeba"**
> Kota Nishi, Atsushi Tero, Yukinori Nishigami, Toshiyuki Nakagaki
> *Journal of Mathematical Biology*, 2026. DOI: [10.1007/s00285-026-02438-7](https://doi.org/10.1007/s00285-026-02438-7)

---

## Requirements
### For the C++ Simulator
- **C++ Compiler**: `g++` / `clang++` (Apple Clang version 17.0 or higher)
- **Target Architecture**: Apple Silicon (`arm64-apple-darwin` verified)
- **C++ Standard**: C++11 or higher (depending on your codebase)

### For Python Scripts (Parameter Estimation & Visualization)
- **Python**: version 3.10.12 or higher
- **Core Standard Libraries**: `multiprocessing`, `functools`, `datetime`, `time`, `os`, `ast`
- **External Libraries**: Install the necessary packages via pip:
  ```bash
  pip install numpy scipy pandas matplotlib deap numba

---

## Usage & Figure Reproduction Steps

Our computational framework is divided into two core phases: **Parameter Estimation (Python)** and **Numerical Simulations (C++)**. To guarantee the rigorous reproducibility of all results, follow the execution paths below.

---

### Phase 1: Parameter Estimation & Sensitivity Analysis (Python)

Before running the simulations, the 7 free parameters are estimated and analyzed using the automated Python framework. 

```bash
# 1. Run the Hybrid Parameter Estimation (Genetic Algorithm via DEAP + Local Search)
# Output: Results are saved in 'param_estimation_result/'
python param_search/param_search.py

# 2. Plot the loss function trajectory (Modified RMSE transition) during optimization
# Output: Results are saved in 'param_estimation_result/'
python param_search/plot_GA.py

# 3. Perform Local Sensitivity Analysis to evaluate parameter robustness
# Output: Results are saved in 'sensitivity_result/'
python param_search/local_sensitivity.py
```

### Phase 2: Core Numerical Simulations (C++)

Navigate to the respective directories, compile the C++ source code using standard g++ commands, and execute the binaries to reproduce the manuscript figures.

```bash
# 1. Main Model Simulation (Bridge-Crossing & Variable Dynamics)
# Output: Results are saved in 'result/fig/'
cd simulation/model_simulation
g++ -o simulate simulate.cpp
./simulate

# 2. Potentiation Simulation (Hallmark of habituation dynamics)
g++ -o potentiation potentiation.cpp
./potentiation

# 3. Frequency sensitivity Simulation (Hallmark of habituation dynamics)
g++ -o frequency_sensitivity frequency_sensitivity.cpp
./frequency_sensitivity

# 4. Subliminal accumulation Simulation (Hallmark of habituation dynamics)
g++ -o subliminal_accumulation subliminal_accumulation.cpp
./subliminal_accumulation

# 5. Phase Plane Analysis (Model behaviors on the phase plane)
# Output: Results are saved in 'result/fig/'
cd ../phase_diagram
g++ -o simulate simulate.cpp
./simulate

# 6. Contour Lines Plot (Analytical iso-velocity & iso-acceleration lines)
# Output: Results are saved in 'result/fig/'
g++ -o contour contour.cpp
./contour

# 7. Perturbation injection simulation (Bridge-Crossing)
# Output: Results are saved in 'result/fig/'
cd ../sensitivity
g++ -o simulate simulate.cpp
./simulate
```
