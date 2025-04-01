function direct_deg(n0, k, Δt, T)
    """
    Simulates the process:
    
    1) Generate r ~ Uniform(0, 1).
    2) If r < A(t)*k*Δt, then A(t+Δt) = A(t) - 1, else A(t+Δt) = A(t).
    
    for t from 0 to T, with step size Δt.
    
    Returns two arrays: 
    - t_vec: the array of time points
    - A_vec: the values of A at each time point
    """
    
    # Number of steps
    nsteps = Int(floor(T / Δt))
    
    # Creates vectors for time and A(t)
    t_vec = [i*Δt for i in 0:nsteps]
    A_vec = zeros(length(t_vec))
    
    # Initial condition
    A_vec[1] = n0
    
    # Main loop
    for i in 1:nsteps
        # Current A
        A_current = A_vec[i]
        
        # Draw a random number in (0,1)
        r = rand()
        
        # Decide if A decreases by 1
        if r < A_current * k * Δt
            A_vec[i+1] = A_current - 1
        else
            A_vec[i+1] = A_current
        end
    end
    
    return t_vec, A_vec
end
