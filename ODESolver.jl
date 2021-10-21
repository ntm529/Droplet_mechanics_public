#First let's define the function we have, and test it with two dropletcount: 
using LinearAlgebra
using DifferentialEquations
using Plots
include("myDropvecDistances.jl")
include("myDropvecRadius_Osmolarity.jl")


#STEP 1: DEFINE THE PROBLEM
function mydropfun!(du,u, p, t)
    ϕ, ζ , m, k, L, γ, dropletcount, dimensions, D, N = p;
    
    
    # Extract the values for radii, osmolarities, positions and velocities
    radii = u[1:dropletcount]
    # display(radii)
    osmolarities = u[Int(dropletcount+1):Int(2*dropletcount)]
    # display(osmolarities)
    positions = u[Int(2*dropletcount+1) : Int(dimensions * dropletcount + 2 * dropletcount)];
    # display(positions)
    velocities = u[(2*dropletcount + dropletcount * dimensions + 1) : end]
    # display(velocities)

     # ================================================

    params =  myDropvecDistances(positions, dropletcount, ϕ , ζ, radii); # Returns, in this order: distance_combinations, norm_combinations, direction, radii_combinations , logic_spring
    distance_combinations = params[1];
    norm_combinations = params[2];
    direction = params[3] ;
    radii_combinations = params[4];
    logic_spring = params[5];
    bilayer_spring_distance = params[6]; 

    dvelocities = sum( (logic_spring .* direction .* -k .* (abs.(bilayer_spring_distance.*direction) .- abs.(distance_combinations)) ), dims = 2) - γ .*  velocities

    #This will become important when making osmolarity activate at the right time 
    #The gradient should be tested as overall increasing or decreasing somehow, ask prof. later if there's a good function to use? 
    #This gradient being negative is one half of the puzzle
    # In addition we will calculate the velocity_norm and have velocity_norm -------> 0 
    ∇velocity = sum(dvelocities);
    normvelocity = norm(velocities)

    # ================================================
    # Perform the osmolarity component
    params2 = myDropvecRadius_Osmolarity(radii, osmolarities, D, dropletcount, norm_combinations, ∇velocity, normvelocity, N)
    
    dradii = params2[1]
    dosmolarities = params2[2]    
    
    
    
    
    dposition = velocities
    # println("Gradient is", gradient) 


    du[1:dropletcount] = dradii
    du[Int(dropletcount+1):Int(2*dropletcount)] = dosmolarities
    du[Int(2*dropletcount + dropletcount * dimensions + 1) : end] = dvelocities
    du[Int(2*dropletcount+1) : Int(dimensions * dropletcount + 2 * dropletcount) ] = dposition
    end
# ===========================================================================================================================================================================

# sum((logic_spring .* direction .* -k.*(abs(bilayer_spring_distance.* direction) - abs(iota))),2) - gamma .* s(2*dropletcount + dropletcount*dimensions + 1:b)]

#  ICmatrix = [1 2 3 0 0 0 ; 4 5 6 0 0 0 ; 7 8 9 0 0 0]; #This is our u0
radii = [3; 3]
osmolarities = [0.5, 1]
display([radii osmolarities])



locations = [0 0 0 ; 6 0 0]
# Note: Settomg velocities [1 1 1 ; 1 1 1] causes an error
velocities = [0.1 0 0 ; 0 0 0]

display([radii osmolarities locations velocities])


ICmatrix = [radii osmolarities locations velocities]; # This is another u0 for testing


# First unravel ICmatrix into a usable ICvector
ICvector = reshape(ICmatrix, :, 1);
display(ICvector)
# Define some of the model parameters 
dimensions = (size(ICmatrix,2)-2)/2;
dimensions = Int(dimensions);
dropletcount = size(ICmatrix,1);
radii = [3, 3]
D = 1 
tspan = (0.0, 100.0);


# initialize the two logical matrixes
ϕ =  .~I(dropletcount);
ϕ = repeat(ϕ, 1, dimensions);
# display(ϕ);

ζ = I(dropletcount-1);
ζ = repeat(ζ,dimensions,1);
ζ = reshape(ζ,(dropletcount-1),:) ;
ζ = ζ';
ζ = repeat(ζ,dropletcount,1);
ζ = reshape(ζ,dimensions*(dropletcount-1),:);

m = 1 
k = 20
L = 0.8
γ = 200
N = Int64[0]

# display(radii)
p = [ϕ, ζ , m, k, L, γ, dropletcount, dimensions, D, N] 
prob = ODEProblem(mydropfun!, ICvector, tspan, p)
sol = solve(prob, reltol=1e-6, saveat = 0.1)

using Plots
plot(sol, vars = 11)
plot!(sol, vars = 12)

plot!(sol, vars = 3)

plot!(sol, vars = 4)

# ##
# w = ICvector[1:Int(dimensions*dropletcount)]
# display(w)
# W = mydropf(ICvector,ϕ,ζ, m, k, L, γ, dropletcount, dimensions, radii)

# display(W)


# #testme(u,phi)


 








































