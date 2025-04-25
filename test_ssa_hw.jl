include("ssa_hw.jl") #call out the location of the file

#ssa_homework
na = 0
nb = 100
nc = 0
k1 = 1.0
k2 = 0.2
k3 = 0.005
T = 100.0

t_1, A_1, B_1, C_1  = ssa_hw(na, nb, nc, k1, k2, k3, T)
t_2, A_2, B_2, C_2 = ssa_hw(na, nb, nc, k1, k2, k3, T)
t_3, A_3, B_3, C_3 = ssa_hw(na, nb, nc, k1, k2, k3, T)

#ODE problem
# Let A(t), B(t), C(t) be the number of A, B, C at t
# Package to use (if you did not install it) 
# using Pkg # use Pkg to add, update, remove, and manage external libraries or packages.
# Pkg.add("DifferentialEquations")

using DifferentialEquations

# Define the problem dA(t)/dt = -k1A(t) +k2v
function f!(du, u, p, t) #Remember to include ! if return a vector
    k1, k2, k3 = p       # Unpack parameters from the tuple
    du[1] = k1 - k2*u[1]*u[2]
    du[2] = -k2*u[1]*u[2] - k3*u[2]
    du[3] = k2*u[1]*u[2]
end
t_ode = (0.0,T)
n = [na, nb, nc]
prob = ODEProblem(f!, n, t_ode, (k1, k2, k3)) #input n is the initial condition

# Solve the ODE 
sol = solve(prob) 
t_sol = sol.t
A_sol = [u[1] for u in sol.u]  # extract A values
B_sol = [u[2] for u in sol.u]  # extract B values
C_sol = [u[3] for u in sol.u]  # extract C values

# Plot simulation
using Plots
using LaTeXStrings
figA = plot(t_1, A_1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Homework using SSA simulation",
     label="Simulation 1", seriestype = :steppost, legend = :outerright)
plot!(t_2, A_2, label="Simulation 2", seriestype = :steppost)
plot!(t_3, A_3, label="Simulation 3", seriestype = :steppost)
plot!(t_sol, A_sol, label = "Solution of ODE", color = :black, linewidth = 2)

figB = plot(t_1, B_1, 
     xlabel="Time t", ylabel="B(t)", 
     title="Homework using SSA simulation",
     label="Simulation 1", seriestype = :steppost, legend = :outerright)
plot!(t_2, B_2, label="Simulation 2", seriestype = :steppost)
plot!(t_3, B_3, label="Simulation 3", seriestype = :steppost)
plot!(t_sol, B_sol, label = "Solution of ODE", color = :black, linewidth = 2)

figC = plot(t_1, C_1, 
     xlabel="Time t", ylabel="C(t)", 
     title="Homework using SSA simulation",
     label="Simulation 1", seriestype = :steppost, legend = :outerright)
plot!(t_2, C_2, label="Simulation 2", seriestype = :steppost)
plot!(t_3, C_3, label="Simulation 3", seriestype = :steppost)
plot!(t_sol, C_sol, label = "Solution of ODE", color = :black, linewidth = 2)

figABC = plot(t_1, A_1, 
     xlabel="Time t", ylabel="N(t)", 
     title="Homework using SSA simulation",
     label="Number of A", seriestype = :steppost, legend = :outerright)
plot!(t_1, B_1, label="Number of B", seriestype = :steppost)
plot!(t_1, C_1, label="Number of C", seriestype = :steppost)


# Histogram
# Simulation


num_sim = 10000
simulations = [ssa_hw(na, nb, nc, k1, k2, k3, T) for i in 1:num_sim]

# Plot the histogram of number of A
results_A = [last(simulations[1][2])]
for i in 2:num_sim
    push!(results_A, last(simulations[i][2]))
end
min_nA = minimum(results_A) #minimum value from all results
max_nA = maximum(results_A) #maximum value from all results
binsA = (min_nA - 0.5):1:(max_nA + 0.5) 

hisA = histogram(results_A,
    bins = binsA,
    normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of A",
    ylabel = "Distribution",
    label = "Gillespie SSA",
    title = "Histogram of number of A"
)

results_B = [last(simulations[1][3])]
for i in 2:num_sim
    push!(results_B, last(simulations[i][3]))
end
min_nB = minimum(results_B) #minimum value from all results
max_nB = maximum(results_B) #maximum value from all results
binsB = (min_nB - 0.5):1:(max_nB + 0.5) 

hisB = histogram(results_B,
    bins = binsB,
    normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of B",
    ylabel = "Distribution",
    label = "Gillespie SSA",
    title = "Histogram of number of B"
)

results_C = [last(simulations[1][4])]
for i in 2:num_sim
    push!(results_C, last(simulations[i][4]))
end
min_nC = minimum(results_C) #minimum value from all results
max_nC = maximum(results_C) #maximum value from all results
binsC = (min_nC - 0.5):1:(max_nC + 0.5) 

hisC = histogram(results_C,
    bins = 30,
    #bins = binsC,
    #normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of C",
    ylabel = "Distribution",
    label = "Gillespie SSA",
    title = "Histogram of number of C"
)


display(figA)
display(figB)
display(figC)
display(figABC)
display(hisA)
display(hisB)
display(hisC)


