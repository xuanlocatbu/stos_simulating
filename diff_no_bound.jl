function diff_no_bound(x0, D, Δt, T)
    # Solve X(t+Δt) = X(t) + \sqrt(Δt)ϵ
    X = [x0]
    tnew = 0.0
    t = [0.0]
    Xnew = x0
    while tnew <T
        ϵ = rand(Normal(0, 1))
        Xnew = X[end] + sqrt(2*D*Δt)*ϵ
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