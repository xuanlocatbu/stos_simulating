#import Pkg; Pkg.add("Plots") #add Plots package to Julia
using Plots
A₀ = 30 #Set the initial value to 30
u₀ = zeros(A₀ + 1) #Create a vector of length 31, indicating a probability vector at t=0
u₀[end] = 1.0        # set p(30,0) = 1.0
bar(0:A₀, u₀; label = "p(a,0)", xlabel = "a")
