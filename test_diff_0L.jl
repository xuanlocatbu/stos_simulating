include("diff_0L.jl")

# Initial data
x0 = 0.4
Δt = 10^-1
T = 240.0
L= 1.0
D = 10^-4
num_sim = 1000
K = 40 #number of bins/ to create compartment

using Plots
using Distributions # Must have this to sample from a normal distribution
simulations = [diff_0L(x0,D, Δt, L, T) for i in 1:num_sim]

# Create an empty plot to build on
p = plot(
    xlabel="x",
    ylabel="t",
    title="Sample Paths of SDE",
    legend = false) #not to inlcude legend 

# Plot each sample path

for i in 1:num_sim
    t, x = simulations[i]
    plot!(p, x, t)
end

#Test with Fokker Planck to plot the distribution at t = 1
location = [last(simulations[1][2])]
for i in 2:num_sim
    push!(location,last(simulations[i][2]))
end
using LaTeXStrings
hist = histogram(location,
    bins = K,
    normalize = :pdf,
    xlabel = L"$x$",
    ylabel = L"$p(x,1)$",
    label = "SSA",
    title = L"The probability distribution $p(x,1)$",
    legend = :outerright
)
# p(x,t) = 1/L + (2/L) * ∑_{n=1}^N cos(nπ x0/L) cos(nπ x/L) exp(−D (nπ/L)^2 t)
function dist(x, t, D, L, x0, N)
    # start with the n=0 (constant) mode
    p = 1/L
    # accumulate n=1..N
    for n in 1:N
        α = n*pi/L
        p += (2/L)*cos(α*x0)*cos(α*x)*exp(-D*α^2*t)
    end
    return p
end

x_values = range(xmin, xmax, length=300)
dist_values = [dist(x, T, D, L, x0, 1000) for x in x_values]

# Overlay the analytic PDF on the same plot
plot!(
    x_values,
    dist_values;
    linewidth = 2,
    color = :red,
    label = L"True distribution $p(x,t)$"
)

display(p)
display(hist)

## Now test the compartment model associated to the diffusion eqn 
include("compartment_0L.jl")
h = L / K
# compute the midpoint x‐coordinate of each compartment
xs = [h*(i-0.5) for i in 1:K]    
compartment = comapartment_0L(x0,D,L,K,num_sim,T)
comp_hist = bar(
  xs, compartment;
  xlabel="x (from 0 to 1)",
  ylabel="Count in compartment",
  xlim=(0,1),
  legend=false,
  bar_width = h*0.9        # a little narrower than h so bars don’t touch
)

display(comp_hist)
