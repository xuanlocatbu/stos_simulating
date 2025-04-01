function simple_sde(x0, Δt, T)
    # Solve X(t+Δt) = X(t) + \sqrt(Δt)ϵ
    X = [x0]
    tnew = 0.0
    t = [0.0]
    Xnew = x0
    while tnew <T
        ϵ = rand(Normal(0, 1))
        Xnew = X[end] + sqrt(Δt)*ϵ
        push!(X, Xnew)
        tnew += Δt
        push!(t, tnew)
    end
    return t, X
end