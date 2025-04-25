include("ssa_prod_deg.jl") #call out the location of the file

#ssa_degen_and_prod
n = 0
k1 = 0.1
k2v = 1.0
T = 100.0

t_vals_1, A_vals_1 = ssa_prod_deg(n, k1, k2v, T)
t_vals_2, A_vals_2 = ssa_prod_deg(n, k1, k2v, T)
t_vals_3, A_vals_3 = ssa_prod_deg(n, k1, k2v, T)

#ODE problem
# Let A(t) be the number at t
# Package to use (if you did not install it) 
# using Pkg # use Pkg to add, update, remove, and manage external libraries or packages.
# Pkg.add("DifferentialEquations")
using DifferentialEquations

# Define the problem dA(t)/dt = -k1A(t) +k2v
function f(a, p, t) #Remember to include ! if return a vector
    k1, k2v = p       # Unpack parameters from the tuple
    return -k1*a + k2v
end
t_ode = (0.0,T)
prob = ODEProblem(f, n, t_ode, (k1, k2v)) #input n is the initial condition

# Solve the ODE 
sol = solve(prob) 

# Plot simulation
using Plots
using LaTeXStrings
pic1 = plot(t_vals_1, A_vals_1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Production and degration using SSA simulation",
     label="Simulation 1", seriestype = :steppost, legend = :outerright)
plot!(t_vals_2, A_vals_2, label="Simulation 2", seriestype = :steppost)
plot!(t_vals_3, A_vals_3, label="Simulation 3", seriestype = :steppost)
plot!(sol, label = L"Solution of $\dfrac{dA(t)}{dt} = -k1A(t) +k2v$", color = :black, linewidth = 2)

# Histogram
# Simulation
num_sim = 10000
results = [last(ssa_prod_deg(n, k1, k2v, T)[2]) for i in 1:num_sim]
min_n = minimum(results) #minimum value from all results
max_n = maximum(results) #maximum value from all results
bins = (min_n - 0.5):1:(max_n + 0.5) 

# Plot the histogram of simulation outcomes
pic2 = histogram(results,
    bins = bins,
    normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of A",
    ylabel = "Distribution",
    label = "Gillespie SSA",
    title = "Histogram of 1000 Simulation Outcomes as t close to 100"
)

# Compare the true distribution function f(n) = C/n!(k2v/k1)^n

rate = k2v/k1
dist_function(n, k1, k2v) = exp(-k2v/k1) * (k2v/k1)^n / factorial(big(n)) #Define a quick the function
n_values = collect(min_n:max_n)
dist_func = [dist_function(n, k1, k2v) for n in n_values]

# Overlay the theoretical distribution using a scatter plot
plot!(n_values, dist_func,
    label = "Theoretical Distribution",
    lw = 2,           # line width
    linecolor = :red,  # line color
    legend = :outerright
)

display(pic1)
display(pic2)