# GRAPE\_n: Generalized GRAPE Algorithm for Quantum State Control

## Overview

This repository implements a generalized version of the **GRAPE (Gradient Ascent Pulse Engineering)** algorithm to minimize infidelity in quantum state evolution. The code supports **multi-qubit systems**, allowing them to be steered towards an **n-dimensional target manifold** with high fidelity.

## Features

- **Generalized multi-qubit GRAPE**: Supports quantum systems with `n` qubits evolving under a control field.
- **Fidelity Computation**: Evaluates the accuracy of quantum state evolution.
- **Gradient-Based Optimization**: Uses `BFGS()` to optimize the control parameters.
- **Custom Evolution Functions**: Integrates quantum evolution dynamics via `evolution.jl`.

## Installation

Ensure you have Julia installed. You also need the following Julia packages:

```julia
using Pkg
Pkg.add(["LinearAlgebra", "Optimization", "OptimizationOptimJL", "Plots", "DiffEqDiffTools", "FiniteDifferences"])
```

## Usage

### 1. Load Required Modules

Include the necessary modules before running any function:

```julia
include("evolution.jl")
```

### 2. Define System Parameters

Set the number of control parameters (`M`), evolution time (`T`), and initial/target states.

```julia
initial_states = [...]  # Define initial quantum states
target_states = [...]   # Define target states
T = 8.0                 # Total evolution time
M = 500                  # Number of control parameters
initial_phi = rand(M)   # Random initial control parameters
alpha = rand(length(initial_states))  # Random phase adjustments
```

### 3. Optimize Quantum Evolution

Call the **GRAPE\_n** function to find the optimal control parameters:

```julia
minimum_infidelity, optimized_phi, optimized_alpha = GRAPE_n(
    initial_states, target_states, T, M, initial_phi, alpha
)
```

### 4. Output Results

The function prints the minimum infidelity achieved and returns the optimized control parameters.

```julia
println("Optimized Control Parameters:", optimized_phi)
println("Optimized Alpha Parameters:", optimized_alpha)
```

## Functions Explained

### `Fidelity_n(vars, T, initial_states, target_states)`

Computes the **fidelity** between the final evolved states and the target states.

### `Infidelity_n(vars, T, initial_states, target_states)`

Computes **infidelity** as `1 - Fidelity_n(...)`.

### `GRAPE_n(initial_states, target_states, T, M, initial_phi, alpha)`

Optimizes the control parameters using gradient-based optimization (`BFGS()`) to minimize the infidelity function.

## Dependencies

- **LinearAlgebra**: For matrix operations.
- **Optimization, OptimizationOptimJL**: For numerical optimization.
- **Plots**: For visualization.
- **DiffEqDiffTools, FiniteDifferences**: For computing gradients.

## Contributions

Feel free to fork, improve, and contribute to this project. Open an issue for any bugs or suggestions!
