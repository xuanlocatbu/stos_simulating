include("ssa_hw.jl")
include("get_value_at.jl")

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
small_num_sim = 5
small_simulations = [ssa_hw(na, nb, nc, k1, k2, k3, T) for i in 1:small_num_sim]

large_num_sim = 100
large_simulations = [ssa_hw(na, nb, nc, k1, k2, k3, T) for i in 1:large_num_sim]

# Plot the mean of A
small_t_A = [small_simulations[1][1]]
small_A = [small_simulations[1][2]]
for i in 2:small_num_sim
    push!(small_t_A,small_simulations[i][1])
    push!(small_A, small_simulations[i][2])
end

large_t_A = [large_simulations[1][1]]
large_A = [large_simulations[1][2]]
for i in 2:large_num_sim
    push!(large_t_A, large_simulations[i][1])
    push!(large_A, large_simulations[i][2])
end

using Plots
using LaTeXStrings
figA_sim = plot(xlabel="Time t", ylabel="A(t)", 
    title="Comparison between ODE and mean of SSA simulation", 
    seriestype = :steppost, legend = false)
plot!(t_sol, A_sol, color = :black, label = "ODE for A", linewidth = 2)

for i in 1:small_num_sim
    plot!(figA_sim, small_t_A[i], small_A[i], seriestype = :steppost)
end

for i in 1:large_num_sim
    plot!(figA_sim, large_t_A[i], large_A[i], seriestype = :steppost)
end

#   add to Amean the values of trajects[i] evaluated at the 
#   times in t by using the evaltraject_attimes function:
Ameans = zeros(length(t_sol))
for i in 1:small_num_sim
    traj = get_value_at(small_t_A[i],small_A[i],t_sol)
    for i in 1:length(Ameans)
        Ameans[i] += traj[i]
    end
end

Ameanl = zeros(length(t_sol))
for i in 1:large_num_sim
    traj = get_value_at(large_t_A[i],large_A[i],t_sol)
    for i in 1:length(Ameanl)
        Ameanl[i] += traj[i]
    end
end

# divide by number of sims to get the mean
# The .= takes that result and writes it back into each element of Amean, 
# rather than rebinding the variable to a fresh array.
Ameans .= Ameans ./ small_num_sim
Ameanl .= Ameanl ./ large_num_sim

figA_mean = plot(xlabel="Time t", ylabel="A(t)", 
    title="Comparison between ODE and mean of SSA simulation", 
    seriestype = :steppost, legend = :outerright)
plot!(t_sol, A_sol, color = :black, linewidth = 2)
plot!(t_sol, Ameans, color = :red,   
    linestyle = :dash, linewidth = 2,
    label = "Average of small simulation using SSA")
plot!(t_sol, Ameanl, color = :blue,   
    linestyle = :dash, linewidth = 2,
    label = "Average of large simulation using SSA")


display(figA_sim)
display(figA_mean)