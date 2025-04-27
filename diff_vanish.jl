function diff_vanish(x0, D, Δt, L, T)
    # Solve X(t+Δt) = X(t) + sqrt(2DΔt)·ε with
    #  • terminate at 0,
    #  • terminate if it hits 0 along the way,
    #  • reflection at L.
    X    = [x0]
    Xnew = x0
    t    = [0.0]
    tnew = 0.0
    con  = sqrt(2 * D * Δt)

    while tnew < T
        ε    = randn()
        Xnew = X[end] + con * ε
        tnew += Δt
        if tnew <= T
            push!(t, tnew)
        else
            push!(t, T)
        end
        # (c13) Absorb at 0
        if Xnew < 0
            # record absorption at 0
            break
        end

        # Interior, possible terminate
        if 0 <= Xnew <= L
            terminate_prob = exp(- X[end] * Xnew / (D * Δt))
            if rand() < terminate_prob
                # record absorption in interior
                break
            end

        else
            # Reflect at L
            Xnew = 2L - X[end] - con * ε
        end

        # record the successful step
        push!(X, Xnew)
        
    end
    return t, X
end
