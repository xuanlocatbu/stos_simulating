function ssa_second_dimerization(n, k1, k2, T)
    """
    Perform the improved SSA for the reaction A + A -> ∅ with rate k1v-1*A(t).
    Perform the improved SSA for the reaction ∅ ->A with rate k2v*A(t).

    Arguments:
      - n: Initial number of molecules of A (A(0))
      - k1v-1, k2v:  Rate constant
      - T:  Final time to stop simulation

    Returns:
      - t_vec: Vector of times at which events (or final check) occurred
      - A_vec: Vector of A(t) values corresponding to t_array
    """
    
    # Initialize time and molecule count
    t = 0
    A = n
    
    # Store the time and A(t) in arrays for plotting or analysis
    t_vec = Float64[t]
    A_vec = [A]
    
    # Main loop: continue until no molecules left or time exceeds T
    while A >= 0 && t < T
        # 1) Generate random number r
        r1 = rand()
        r2 = rand()
        
        a1 = A_vec[end]*(A_vec[end]-1)*k1 
        a2 = k2
        r = a1 + a2
        # 2) Compute tau = 1/(A*k)*ln(1/r)
        tau = (1 / r) * log(1 / r1)
        
        # 3) Advance time
        t += tau
        
        # Check if the new time is still within T
        if t <= T
            if r2 < a1/r
                A -= 2
            else
                A += 1
            end
            # Save the new state. Each time we update t anf A, we append them to t_vec and A_vec:
            push!(t_vec, t)
            push!(A_vec, A)
        else
            push!(t_vec, T) #adding the last time
            push!(A_vec, A)# If we've passed T, we stop; 
            break
        end
    end
    return t_vec, A_vec
end
