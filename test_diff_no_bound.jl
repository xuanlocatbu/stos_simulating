include("simple_sde.jl")

# Initial data
x0 = 0.0
Δt = 10^-2
T = 1.0
D = 10^-4
num_sim = 20

using Plots
using Distributions # Must have this to sample from a normal distribution
simulations = [simple_sde(x0,Δt,T) for i in 1:num_sim]

# Create an empty plot to build on
p = plot(
    xlabel="t",
    ylabel="X",
    title="Sample Paths of SDE",
    legend = false) #not to inlcude legend 

# Plot each sample path

for i in 1:num_sim
    t, x = simulations[i]
    plot!(p, t, x)
end

display(p)