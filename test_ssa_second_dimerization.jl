include("ssa_second_dimerization.jl") #call out the location of the file

#ssa_second_dimerization
n = 0
k1 = 0.005
k2 = 1.0
T = 100.0

t_vals_1, A_vals_1 = ssa_second_dimerization(n, k1, k2, T)
t_vals_2, A_vals_2 = ssa_second_dimerization(n, k1, k2, T)
t_vals_3, A_vals_3 = ssa_second_dimerization(n, k1, k2, T)

#ODE problem
# Let A(t) be the number at t
# Package to use (if you did not install it) 
# using Pkg # use Pkg to add, update, remove, and manage external libraries or packages.
# Pkg.add("DifferentialEquations")
using DifferentialEquations

# Define the problem dA(t)/dt = -k1A(t) +k2v
function f!(a, p, t) #Always remember to include !
    k1, k2 = p       # Unpack parameters from the tuple
    return -2*k1*a^2 + k2
end
t_ode = (0.0,T)
prob = ODEProblem(f!, n, t_ode, (k1, k2)) #input n is the initial condition

# Solve the ODE 
sol = solve(prob) 

# Plot simulation
using Plots
using LaTeXStrings
pic1 = plot(t_vals_1, A_vals_1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Second Dimerization using SSA simulation",
     label="Simulation 1", seriestype = :steppost, legend = :outerright)
plot!(t_vals_2, A_vals_2, label="Simulation 2", seriestype = :steppost)
plot!(t_vals_3, A_vals_3, label="Simulation 3", seriestype = :steppost)
plot!(sol, label = L"Solution of $\dfrac{dA(t)}{dt} = -2k1v^{-1}A^2(t) +k2$", color = :black, linewidth = 2)

# Histogram
# Simulation
num_sim = 10000
results = [last(ssa_second_dimerization(n, k1, k2, T)[2]) for i in 1:num_sim]
min_n = minimum(results) #minimum value from all results
max_n = maximum(results) #maximum value from all results
bins = (min_n - 1):1:(max_n + 1) 

# Plot the histogram of simulation outcomes
pic2 = histogram(results,
    bins = bins,
    normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of A",
    ylabel = "Distribution",
    label = "Gillespie SSA",
    title = "Histogram of 10000 Simulation Outcomes as t close to 100"
)

# Compare the true distribution function 
#Computes the normalizing constant C₁ for the distribution.
using SpecialFunctions
#Computes the stationary probability for having n molecules.
function phi(n, k1, k2)
    # Compute the piece (k2*ν^2 / k1)^n / n!
    coeff = (k2 / k1)^n / (factorial(big(n)) *sqrt(2) * besseli(1, 2 * sqrt(2 * k2 / k1)) )
    # Compute the modified Bessel function I_{n-1}
    bessel_part = besseli(n-1, 2 * sqrt(k2/k1))
    return coeff * bessel_part
end

n_values = collect(min_n:max_n)
dist_func = [phi(n, k1, k2) for n in n_values]

# Overlay the theoretical distribution using a scatter plot

plot!(n_values, dist_func,
    label = "Theoretical Distribution",
    lw = 2,           # line width
    linecolor = :red,  # line color
    legend = :outerright
)

display(pic1)
display(pic2)