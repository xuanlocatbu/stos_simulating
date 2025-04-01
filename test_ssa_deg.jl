include("ssa_deg.jl") #call out the location of the file

#direct_degen
n = 20
k = 0.1
T = 30

t_vals_1, A_vals_1 = ssa_degen(n, k, T)
t_vals_2, A_vals_2 = ssa_degen(n, k, T)
t_vals_3, A_vals_3 = ssa_degen(n, k, T)
using Plots
plot(t_vals_1, A_vals_1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Simulation of A(t) using SSA simulation",
     label="Simulation 1", seriestype = :steppost)
plot!(t_vals_2, A_vals_2, label="Simulation 2", seriestype = :steppost)
plot!(t_vals_3, A_vals_3, label="Simulation 3", seriestype = :steppost)
