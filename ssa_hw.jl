function ssa_hw(na, nb, nc, k1, k2, k3, T)
    """
    Perform the improved SSA for the reaction ∅ -> A with rate k1.
    Perform the improved SSA for the reaction A + B -> C with rate k2.
    Perform the improved SSA for the reaction B -> ∅ with rate k3.

    Arguments:
      - na, nb, nc: Initial number of molecules of A, B, C (A(0), B(0), C(0))
      - k1, k2, k3:  Rate constant
      - T:  Final time to stop simulation

    Returns:
      - t_vec: Vector of times at which events (or final check) occurred
      - A_vec: Vector of A(t) values corresponding to t_array
      - B_vec: Vector of B(t) values corresponding to t_array
      - C_vec: Vector of C(t) values corresponding to t_array
    """
    
    # Initialize time and molecule count
    t = 0
    A = na
    B = nb
    C = nc
    
    # Store the time and A(t) in arrays for plotting or analysis
    t_vec = Float64[t]
    A_vec = [A]
    B_vec = [B]
    C_vec = [C]
    
    # Main loop: continue until no molecules left or time exceeds T
    while A >= 0 && B >= 0 && C >= 0 && t < T
        # 1) Generate random number r
        r1 = rand()
        r2 = rand()
        
        
        p = [k1]
        push!(p, p[1] + k2*A*B)
        push!(p, p[2]+ k3*B)
        
        # 2) Compute tau = 1/(A*k)*ln(1/r)
        tau = (1 / p[end]) * log(1 / r1)
        
        # 3) Advance time
        t += tau
        
        reaction = 0
        for i in 1:3
            if r2 <= p[i]/p[end]
                reaction = i
                break
            end
        end

        # Check if the new time is still within T
        if t <= T
            if reaction == 1
                A += 1
            elseif reaction == 2
                A -= 1
                B -= 1
                C += 1
            else
                B-= 1
            end
            # Save the new state. Each time we update t anf A, we append them to t_vec and A_vec:
            push!(t_vec, t)
            push!(A_vec, A)
            push!(B_vec, B)
            push!(C_vec, C)
        else
            # If we've passed T, we stop; 
            break
        end
    end
    return t_vec, A_vec, B_vec, C_vec
end
