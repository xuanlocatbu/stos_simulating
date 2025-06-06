function simple_sde(x0, Δt, T)
    # Solve X(t+Δt) = X(t) + \sqrt(Δt)ϵ
    X = [x0]
    tnew = 0.0
    t = [0.0]
    Xnew = x0
    sqdt = sqrt(Δt)
    while tnew <T
        ϵ = randn()
        Xnew = X[end] + sqdt*ϵ
        push!(X, Xnew)
        tnew += Δt
        if tnew  <= T
            push!(t, tnew)
        else
            push!(t, T)
        end
    end
    return t, X
end