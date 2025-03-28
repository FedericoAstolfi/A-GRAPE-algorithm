using LinearAlgebra
using Optimization, OptimizationOptimJL
using Plots
using DiffEqDiffTools
using FiniteDifferences
include("evolution.jl")

# ------------------------ Target Manifold Functions ------------------------ #
"""
Computes the target states with global phase factors.

# Arguments
- `alpha::Vector`: Phase parameters.
- `target_states::Vector`: Target states to be steered to.

# Returns
- A vector of target states with applied phases.
"""
function target_n(alpha::Vector, target_states::Vector)
    return [target_states[i] * exp(im * alpha[i]) for i in eachindex(target_states)]
end

# ------------------------ Fidelity Computation ------------------------ #
"""
Computes the Fidelity between final and initial states for n qubits.

# Arguments
- `vars`: Control parameters (phase phi and global phases alpha).
- `T`: Evolution time.
- `initial_states`: Initial states of the qubits.
- `target_states`: Target states.

# Returns
- Fidelity value (a measure of state overlap).
"""
function Fidelity_n(vars, T, initial_states, target_states)
    M = length(vars) - length(initial_states)  # Extract M from vars size
    phi = vars[1:M]                            # Control parameters
    alpha = vars[M+1:end]                      # Phase parameters
    
    # Evolve initial states with phi
    final_states = evolution(initial_states, phi, T, "psi")
    
    # Apply global phases to target states
    target_states = target_n(alpha, target_states)
    
    # Compute inner products
    N = length(initial_states)
    products = Complex{Float64}[final_states[k]' * target_states[k] for k in 1:N]
    
    # Compute fidelity
    fidelity = (1 / (N^2 + N)) * (abs(sum(products))^2 + sum(abs.(products).^2))
    
    return fidelity
end

"""
Computes the infidelity as 1 - Fidelity.
"""
function Infidelity_n(vars, T, initial_states, target_states)
    return 1 - Fidelity_n(vars, T, initial_states, target_states)
end

# ------------------------ GRAPE Optimization ------------------------ #
"""
Generalized GRAPE algorithm for minimizing the Infidelity function.
Supports n qubits steered to an n-dimensional target manifold.

# Arguments
- `initial_states`: Initial quantum states.
- `target_states`: Desired target states.
- `T`: Total evolution time.
- `M`: Number of pieces in which phi is discretized.
- `initial_phi`: Initial guess for control parameters.
- `alpha`: Initial guess for phase parameters.

# Returns
- `minimum_value`: Minimum infidelity achieved.
- `optimized_phi`: Optimized laser phase.
- `optimized_alpha`: Optimized phase parameters.
"""
function GRAPE_n(initial_states, target_states, T, M, initial_phi, alpha)
    # Ensure the input dimensions match
    if size(target_states) != size(initial_states)
        throw(ErrorException("Target and initial states have different dimensions!"))
    end
    
    N = length(initial_states)
    
    # Initialize states
    psi = initial_states
    chi = target_n(alpha, target_states)  # Target function with phases
    psi_targ = deepcopy(chi)
    
    # Placeholder for gradient
    gradI = Vector{Float64}(undef, M + N)
    
    # Normalization factor for gradient computation
    signed_normalization_factor = -2 / (N^2 + N)
    
    # Compute gradient over control parameters
    for m in 1:M
        chi = evolution(chi, initial_phi, T, "chi")
        der = [chi[i]' * (dU_dphi(i, T/M, initial_phi[m]) * psi[i]) for i in 1:N]
        psi = evolution(psi, initial_phi, T, "psi")
        gradI[m] = signed_normalization_factor * real(sum(chi[i]' * psi[i] * der[i] for i in 1:N))
    end
    
    # Compute gradient over phase parameters
    dalpha = [-sin(alpha[i]) + cos(alpha[i]) * im for i in 1:N]
    for i in 1:N
        gradI[M+i] = signed_normalization_factor * real(sum(chi[j]' * psi[j] * dalpha[j] for j in 1:N))
    end
    
    # Optimize the infidelity function
    result = optimize(vars -> Infidelity_n(vars, T, initial_states, psi_targ),
                      vcat(initial_phi, alpha), BFGS())
    
    # Extract results
    minimum_value = Optim.minimum(result)
    println("\nMinimum Infidelity:", minimum_value)
    optimized_vars = Optim.minimizer(result)
    optimized_phi = optimized_vars[1:M]
    optimized_alpha = optimized_vars[M+1:end]
    
    return minimum_value, optimized_phi, optimized_alpha
end