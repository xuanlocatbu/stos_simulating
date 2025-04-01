include("direct_deg.jl") #call out the location of the file

#direct_degen
n = 20
k = 0.1
Δt = 0.005
T = 30

t_vals_1, A_vals_1 = direct_deg(n, k, Δt, T)
t_vals_2, A_vals_2 = direct_deg(n, k, Δt, T)
t_vals_3, A_vals_3 = direct_deg(n, k, Δt, T)
using Plots
plot(t_vals_1, A_vals_1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Simulation of A(t) using direct simulation",
     label="Simulation 1")
plot!(t_vals_2, A_vals_2, label="Simulation 2")
plot!(t_vals_3, A_vals_3, label="Simulation 3")
