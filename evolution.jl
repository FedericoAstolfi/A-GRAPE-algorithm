using LinearAlgebra


function U(k, dt, phi_m)
    """
    Temporal evolution operator for a single qubit.
    
    Arguments:
    - k: Index of the qubit being evolved.
    - dt: Discretized time step.
    - phi_m: Discretized phase value at time m*dt.
    """
    cos_term = cos(sqrt(k) * dt / 2)
    sin_term = sin(sqrt(k) * dt / 2) * exp(im * phi_m)
    
    return [cos_term  sin_term;
            -sin_term' cos_term]
end

function dU_dphi(k, dt, phi_m)
    """
    Derivative of the evolution operator U with respect to phi.
    
    Arguments:
    - k: Index of the qubit being evolved.
    - dt: Discretized time step.
    - phi_m: Discretized phase value at time m*dt.
    """
    sin_term = sin(sqrt(k) * dt / 2)
    exp_term = exp(im * phi_m)
    
    return [0      -im * sin_term * exp_term;
            im * sin_term * exp_term' 0]
end

function evolution(initial_states, phi, T, key)
    """
    Evolves a system of qubits forward or backward in time.
    
    Arguments:
    - initial_states: Vector of N qubit states (each with 2 complex components).
    - phi: Discretized laser pulse sequence.
    - T: Total evolution time.
    - key: Evolution direction ('psi' for forward, 'chi' for backward).
    """
    if !(key in ["psi", "chi"])
        throw(ErrorException("Invalid key: must be either 'psi' or 'chi'"))
    end
    
    N = length(initial_states)  # Number of qubits
    M = length(phi)             # Number of time steps
    dt = T / M                 # Discretized time step
    evolved_states = deepcopy(initial_states)
    
    for k in 1:N
        # Evolve each qubit independently
        if key == "psi"
            for m in 1:M
                evolved_states[k] = U(k, dt, phi[m]) * evolved_states[k]
            end
        elseif key == "chi"
            for m in 1:M
                evolved_states[k] = U(k, dt, phi[M + 1 - m]) * evolved_states[k]
            end
            reverse!(evolved_states[k])
        end
    end
    
    return evolved_states
end
