module MyDropletVectorJ!
include("mydropletvectorjcomponents.jl")
using .MyDropletVectorJComponents
using LinearAlgebra
using DifferentialEquations
export dropvecj!

function dropvecj!(m, k, L , γ, D, timeend, ICmatrix, boolean_osm)

ICvector = reshape(ICmatrix, :, 1);
# Define some of the model parameters 
dimensions = (size(ICmatrix,2)-2)/2;
dimensions = Int(dimensions);
dropletcount = Int(size(ICmatrix,1));
tspan = (0.0, timeend);
# initialize the two logical matrixes
ϕ =  .~I(dropletcount);
ϕ = repeat(ϕ, 1, dimensions);

ζ = I(dropletcount-1);
ζ = repeat(ζ,dimensions,1);
ζ = reshape(ζ,(dropletcount-1),:) ;
ζ = ζ';
ζ = repeat(ζ,dropletcount,1);
ζ = reshape(ζ,dimensions*(dropletcount-1),:);



# initialize certain empty matrixes
# The 3 value is the number of dimensions, x, y, z 

distance_combinations = zeros(dropletcount*3, dropletcount-1)

direction = zeros(size(distance_combinations))



sphere_intersection = zeros(dropletcount, dropletcount-1)
sphere_dist_component_1 = zeros(dropletcount, dropletcount-1)
sphere_dist_component_2 = zeros(dropletcount, dropletcount-1)



circular_area = zeros(dropletcount, dropletcount-1)
comosm = zeros(dropletcount, dropletcount-1)
J = zeros(dropletcount, dropletcount-1) 
volume = zeros(dropletcount,1)
moles = zeros(dropletcount,1)
dradii = zeros(dropletcount,1)
dosmolarities = zeros(dropletcount,1)




norm_combinations = zeros(dropletcount*3,dropletcount-1)
comradii = zeros(dropletcount, dropletcount-1); 
radii_combinations = zeros(dropletcount, dropletcount-1)
bilayer_spring_distance = zeros(dropletcount*3, dropletcount-1); 
logic_spring = zeros(dropletcount*3, dropletcount-1)





# Create the current residual that makes the osmolarity activate after the droplets settle


p = [ϕ, ζ, sphere_intersection , distance_combinations, norm_combinations, direction,   m, k, L, γ, dropletcount, dimensions, D, boolean_osm, comradii, bilayer_spring_distance, circular_area, radii_combinations, comosm, J, logic_spring, volume, moles, dradii, dosmolarities, sphere_dist_component_1, sphere_dist_component_2] 
prob = ODEProblem(mydropfun!, ICvector, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, saveat = 0.1)

W = myOvitoPrint!(sol, dropletcount, dimensions)
outfile = "E:/Research Scripts and Functions/OVITO_Files/visual_model.ovito"
f = open(outfile, "w")
for i in eachindex(W)   
    println(f, W[i])
end
close(f)

end


# END module
end
