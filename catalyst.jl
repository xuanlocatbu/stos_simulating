using Catalyst, OrdinaryDiffEqDefault, Plots

# Deterministic ODE simulation
# Create model.
model = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

# Create an ODE that can be simulated.
u0 = [:S => 50.0, :E => 10.0, :SE => 0.0, :P => 0.0]
tspan = (0., 200.)
ps = [:kB => 0.01, :kD => 0.1, :kP => 0.1]
ode = ODEProblem(model, u0, tspan, ps)

# Simulate ODE and plot results.
sol = solve(ode)
ode = plot(sol; lw = 5)
plot!(legend = :outerright)

display(ode)

# Stochastic jump simulations
# Create and simulate a jump process (here using Gillespie's direct algorithm).
# The initial conditions are now integers as we track exact populations for each species.
using JumpProcesses # Have to import JumpProcesses, import Pkg; Pkg.add("JumpProcesses")
u0_integers = [:S => 50, :E => 10, :SE => 0, :P => 0]
jinput = JumpInputs(model, u0_integers, tspan, ps)
jprob = JumpProblem(jinput)
jump_sol = solve(jprob)
stos = plot(jump_sol; lw = 2)
plot!(legend = :outerright)
