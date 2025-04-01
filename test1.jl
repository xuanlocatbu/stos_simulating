#import Pkg; Pkg.add("Plots") #add Plots package to Julia
using Plots
A0 = 30 #Set the initial value to 30
u0 = zeros(A0 + 1) #Create a vector of length 31, indicating a probability vector at t=0
u0[end] = 1.0        # set p(30,0) = 1.0
bar(0:A0, u0; label = "p(a,0)", xlabel = "a")
