function diff_0L(x0, D, Δt, L, T)
    # Solve X(t+Δt) = X(t) + \sqrt(2DΔt)ϵ with reflective boundary on [0,L]
    X = [x0]
    tnew = 0.0
    t = [0.0]
    Xnew = x0
    con = sqrt(2*D*Δt)
    while tnew <T
        ϵ = randn()
        Xnew = X[end] + con*ϵ
        if Xnew <= 0
            Xnew = - X[end] - con*ϵ
        elseif Xnew >= L 
            Xnew = 2*L - X[end] - con*ϵ
        end
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