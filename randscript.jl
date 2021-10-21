include("myDropvecDistances.jl")
include("myDropvecRadius_Osmolarity.jl")
include("mydropfun!.jl")
include("myOvitoPrint!.jl")
include("myLoopMatrix.jl")
include("myDropletVectorJ!.jl")

using LinearAlgebra
using DifferentialEquations
using Plots
using DataFrames 







# Note: Setting velocities [1 1 1 ; 1 1 1] causes an error
radii = 3 
dropletcount = 38
osmolaritybase = 0.1
osmolaritytop = 1 

m = 1 
k = 100
L = 0.8
γ = 200
timeend = 100.0
D = 0.01 




ICmatrix = myLoopMatrix(radii ,dropletcount, osmolaritybase, osmolaritytop); # This is another u0 for testing
display(ICmatrix)



myDropletVectorJ!(m, k, L , γ, D, timeend, ICmatrix)


