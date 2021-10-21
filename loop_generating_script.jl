module XX

include("myLoopMatrix.jl")
include("E:/Research Scripts and Functions/Julia Scripts/mydropletvectorj!.jl")

using .MyDropletVectorJ!
using LinearAlgebra
using DifferentialEquations
using Plots
using DataFrames 


#Characters I might use for some initialized matrixes
# ℵ ℶ ℷ ℸ 



radii = 3
dropletcount =  38
osmolaritybase = 0.1
osmolaritytop = 1 

m = 1 
k = 10000
L = 0.8
γ = 200
timeend = 500.0
D = 0.01
boolean_osm = [0.0]



ICmatrix = myLoopMatrix(radii ,dropletcount, osmolaritybase, osmolaritytop); # This is another u0 for testing
display(ICmatrix)

@time dropvecj!(m, k, L , γ, D, timeend, ICmatrix, boolean_osm)
print("done")


end 
