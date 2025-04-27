"""
    simulate_diffusion_compartments(x0, D, L, tmax)

Simulate pure diffusion on a 1D chain of M compartments via Gillespie SSA.

# Arguments
- `x0::Vector{Int}` : initial molecule counts in each of M compartments
- `D::Float64`      : diffusion coefficient
- `L::Float64`      : total length of the domain
- `tmax::Float64`   : final time to simulate up to

# Returns
- `ts::Vector{Float64}`       : event times, starting at t=0
- `traj::Vector{Vector{Int}}`: list of state‐vectors at each recorded time
"""
function comapartment_0L(x0, D, L, K, sim, tmax)
    # compartment width
    h = L / K
    # diffusion jump rate per particle
    d = D / h^2

    A = zeros(Int, K)
    i_exact = x0/h
    # state
    if isinteger(i_exact) 
        i = Int(i_exact)
        A[i] = div(sim,2)
        A[i + 1] = sim - A[i]
    else
        i = floor(Int,i_exact)+1
        A[i] = sim 
    
    end

    # recorders
    t = 0.0

    while t < tmax
        # 1) Generate random number r
        r1 = rand()
            #Initilize
        a = zeros(2*(K-1))
        for j in 1:K-1
            a[j] = d * A[j]      # j → j+1
            a[K-1 + j] = d * A[j+1]    # j+1 → j
        end

        a0 = sum(a)
        # draw next reaction time
        τ  = -log(r1)/a0
        t += τ

        if t <= tmax
            r2 = rand()*a0
            cumsum = 0.0
            index = findfirst(j -> (cumsum += a[j]) ≥ r2, eachindex(a))
    # eachindex(a) generates an array of indices of the array a.
    """
    findfirst(f, collection)

This scans through collection in order, applying the predicate function f to each element.

It returns the first element of collection for which f(elem) is true.

Here, collection is the list of indices 1,2,…,length(a), and f(j) is our anonymous cumul-check.
    """

            # update state
            if index ≤ K-1
                # forward jump index → μ+1
                A[index]   -= 1
                A[index+1] += 1
            else
                # backward jump: index k = index-(K-1), so (k+1)→k
                k = index - (K-1)
                A[k+1] -= 1
                A[k]   += 1
            end
        end
    end
    return A
end
