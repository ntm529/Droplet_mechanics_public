include("myLoopMatrix.jl")
include("E:/Research Scripts and Functions/Julia Scripts/mydropletvectorj!.jl")

using .MyDropletVectorJ!
using LinearAlgebra
using DifferentialEquations
using Plots
using DataFrames  

radii = [3, 3]
osmolarities = [0.1, 1]
D = 1 
tspan = (0.0, 100.0);
m = 1 
k = 100
L = 0.8
γ = 200
locations = [0 0 0; 6 0 0]
velocities = [0 0 0; 0 0 0]
timeend = 200.0




ICmatrix = [radii osmolarities locations velocities]; # This is another u0 for testing


# First unravel ICmatrix into a usable ICvector
ICvector = reshape(ICmatrix, :, 1);
display(ICvector)
# Define some of the model parameters 
dimensions = (size(ICmatrix,2)-2)/2;
dimensions = Int(dimensions);
dropletcount = size(ICmatrix,1);

# display(radii)
dropvecj!(m, k, L , γ, D, timeend, ICmatrix)
print("done")


# Now use this as an input for the next function 



# str = ["Carlos" "Sada" ; "Ella" "Olsen"]
# permutedims(str)
# outfile = "random_script.txt"
# f = open(outfile, "w")

# for i in eachindex(str)   
#     println(f, str[i])
# end

# close(f)













# newArray = Array(Float64, size(sol))

# Message #general2
# newarray = convert(Array{Float64,2}, sol)
# a = Iterators.Stateful("sol");
# newArray = Array{Float64, size(sol)}

# for x in size(sol)
#     newArray[x] = sol[x][1]
# end 



# [t+u for (u,t) in tuples(sol)]

# myOvitoPrint!(sol, dropletcount, dimensions)


# plot(sol, vars = 11)
# plot!(sol, vars = 12)
# plot!(sol, vars = 3)
# plot!(sol, vars = 4)