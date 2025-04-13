include("ssa_hw.jl")

#Define ODE problem
na = 0
nb = 100
nc = 0
k1 = 1.0
k2 = 0.2
k3 = 0.005
T = 150.0

#ODE problem
# Let A(t), B(t), C(t) be the number of A, B, C at t
# Package to use (if you did not install it) 
# using Pkg # use Pkg to add, update, remove, and manage external libraries or packages.
# Pkg.add("DifferentialEquations")

using DifferentialEquations

# Define the problem dA(t)/dt = -k1A(t) +k2v
function f!(du, u, p, t) #Always remember to include !
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


# for small number of simulations, calculate the mean
small_num_sim = 10
small_simulations = [ssa_hw(na, nb, nc, k1, k2, k3, T) for i in 1:small_num_sim]

# Plot the mean of A
small_t_A = [small_simulations[1][1]]
small_A = [small_simulations[1][2]]
for i in 2:small_num_sim
    push!(small_t_A, small_simulations[i][1])
    push!(small_A, small_simulations[i][2])
end

using Plots
using LaTeXStrings
figA_sim = plot(xlabel="Time t", ylabel="A(t)", 
    title="Comparison small_A and ODE using SSA simulation", 
    seriestype = :steppost, legend = false)
plot!(t_sol, A_sol, color = :black, label = "ODE for A", linewidth = 2)

for i in 1:small_num_sim
    plot!(figA_sim, small_t_A[i], small_A[i], seriestype = :steppost)
end

using Statistics

mean_A = mean(small_A)

figA_mean = plot(xlabel="Time t", ylabel="A(t)", 
    title="Comparison mean small_A and ODE using SSA simulation", 
    seriestype = :steppost, legend = :outerright)
plot!(t_sol, A_sol, color = :black, linewidth = 2)
plot!(t_sol, small_A, color = :red,   
    linestyle = :dash, linewidth = 2,
    label = "Average of small_A using SSA")


display(figA_sim)
display(figA_mean)