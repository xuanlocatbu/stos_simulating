using Random

function ssa_degen(n, k, T)
    """
    Perform the improved SSA for the reaction A -> ∅ with rate k*A(t).

    Arguments:
      - n: Initial number of molecules of A (A(0))
      - k:  Rate constant
      - T:  Final time to stop simulation

    Returns:
      - t_array: Vector of times at which events (or final check) occurred
      - A_array: Vector of A(t) values corresponding to t_array
    """
    
    # Initialize time and molecule count
    t = 0.0
    A = n
    
    # Store the time and A(t) in arrays for plotting or analysis
    t_vec = Float64[t]
    A_vec = [A]
    
    # Main loop: continue until no molecules left or time exceeds T
    while A > 0 && t < T
        # 1) Generate random number r from the exp(1) dis
        r = randexp()
        
        # 2) Compute tau = 1/(A*k)*ln(1/r)
        #    Note: A*k must be > 0; if A=0, we won't enter this loop
        tau = (1 / (A_vec[end]*k)) * r
        
        # 3) Advance time
        t += tau
        
        # Check if the new time is still within T
        if t <= T
            # 4) Decrement the molecule count by 1
            A -= 1
            
            # Save the new state. Each time we update t anf A, we append them to t_vec and A_vec:
            push!(t_vec, t)
            push!(A_vec, A)
        else
            push!(t_vec,T)
            push!(A_vec,A)
            # If we've passed T, we stop; 
            break
        end
    end
    
    return t_vec, A_vec
end
