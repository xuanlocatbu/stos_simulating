include("simple_sde.jl")

# Initial data
x0 = 0.0
Δt = 10^-2
T = 1.0
num_sim = 10000

# X(t+Δt) = X(t) + sqrt(2Δt)ϵ
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

for i in 1:20
    t, x = simulations[i]
    plot!(p, t, x)
end

# Include the dotted mean line
t = simulations[1][1]
y = zeros(length(t))
plot!(p, t, y, label="The mean M(t)", line=:dash, color=:black)

#The plot will not be shown without this code
display(p)

#Test with Fokker Planck to plot the distribution at t = 1
location = [last(simulations[1][2])]
for i in 2:num_sim
    push!(location,last(simulations[i][2]))
end
using LaTeXStrings
Δx = 0.1
xmin = minimum(location)
xmax = maximum(location)
bins = xmin:Δx:xmax
hist = histogram(location,
    bins = bins,
    normalize = :pdf,
    xlabel = L"$x$",
    ylabel = L"$p(x,1)$",
    label = "SSA",
    title = L"The probability distribution $p(x,1)$",
    legend = :outerright
)
# p(x,t) = (1 / sqrt(2πt)) * exp(-x^2 / (2t))
function dist(x, t)
    return 1 / sqrt(2 * pi * t) * exp(-x^2 / (2 * t))
end
x_values = range(xmin, xmax, length=300)
dist_values = [dist(x, 1) for x in x_values]

# Overlay the analytic PDF on the same plot
plot!(
    x_values,
    dist_values;
    linewidth = 2,
    color = :red,
    label = L"True distribution $p(x,t)$"
)

display(hist)
