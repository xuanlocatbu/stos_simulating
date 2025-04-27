include("diff_vanish.jl")

# Initial data
x0 = 0.4
Δt = 10^-1
T = 15.0
L= 1.0
D = 10^-3
num_sim = 1000
K = 40 #number of bins/ to create compartment

using Plots
using Distributions # Must have this to sample from a normal distribution
simulations = [diff_vanish(x0, D, Δt, L, T) for i in 1:num_sim]

# Create an empty plot to build on
p = plot(
    xlabel="t",
    ylabel="x",
    title="Sample Paths of SDE",
    legend = false) #not to inlcude legend 

# Plot each sample path

for i in 1:num_sim
    t, x = simulations[i]
    if last(t) == T
        plot!(p, t, x)
    end    
end

#Test with Fokker Planck to plot the distribution at t = 1
location = [last(simulations[1][2])]
for i in 2:num_sim
    if last(simulations[i][1]) == T
        push!(location,last(simulations[i][2]))
    end    
end
using LaTeXStrings
hist = histogram(location,
    bins = K,
    normalize = :pdf,
    xlabel = L"$x$",
    ylabel = L"$p(x,1)$",
    label = "Brownian Motions",
    title = L"The probability distribution $p(x,1)$",
    legend = :outerright
)

display(hist)
display(p)