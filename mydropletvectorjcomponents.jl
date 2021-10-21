module MyDropletVectorJComponents
include("myodecomponents.jl")
using .myODEcomponents
using LinearAlgebra
export mydropfun!, myOvitoPrint!


function mydropfun!(du,u, p, t)
    ϕ, ζ , sphere_intersection, distance_combinations, norm_combinations, direction, m, k, L, γ, dropletcount, dimensions, D, boolean_osm, comradii, bilayer_spring_distance, circular_area, radii_combinations, comosm, J, logic_spring, volume, moles, dradii, dosmolarities, sphere_dist_component_1, sphere_dist_component_2 = p;
    
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

    params =  myDropvecDistances(positions, dropletcount, ϕ , ζ, radii, L, dimensions, comradii, bilayer_spring_distance, norm_combinations, distance_combinations, direction, radii_combinations, logic_spring); # Returns, in this order: distance_combinations, norm_combinations, direction, radii_combinations , logic_spring
    comradii = params[1];
    distance_combinations = params[2];
    norm_combinations = params[3];
    direction = params[4] ;
    logic_spring = params[5];
    bilayer_spring_distance = params[6]; 

    dvelocities = sum( (logic_spring .* direction .* -k .* (abs.(bilayer_spring_distance.*direction) .- abs.(distance_combinations)) ), dims = 2) - γ .*  velocities

    # Use the norm velocity to find later where the boolean goes
    normvelocity = norm(velocities)
    
    # ================================================
    # Perform the osmolarity component
    params2 = myDropvecRadius_Osmolarity(comradii, radii, osmolarities, D, dropletcount, norm_combinations, normvelocity, boolean_osm, t, sphere_intersection, circular_area, comosm, J, volume, moles, dradii, dosmolarities, sphere_dist_component_1, sphere_dist_component_2)

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





function myOvitoPrint!(sol, dropletcount, dimensions)

    # This function will recast the solution of our ordinary differential equation into a printed file that can be visualized in ovito.
    
    V = convert(Array,sol) # Convert the ODE solution into an array that can be worked on 
    V = dropdims(V, dims = 2) # Dimension 2 is a singleton and should be eliminated
    V = V'
    tsteps = size(V,1)
    omicron = .~I(dimensions);
    omicron = repeat(omicron, dropletcount, 1)
    omicron = reshape(omicron, dimensions, :)
    omicron = omicron'
    omicron = repeat(omicron, tsteps, 1)
    omicron = reshape(omicron, dimensions*dropletcount, :)
    
    # Shave off the velocities for now 
    W = V[:,2*dropletcount+1:2*dropletcount + dimensions*dropletcount]; ##ADJUST##
    # Perform a sequence of transformations to re-arrange W to what we want
    Wtrue = [W[:,1:2]; W[:,3:4] ; W[:,5:6]]
    Wtransptrue = Wtrue'
    W = W'
    W = repeat(W, 1, dimensions)
    W = W[.~omicron]
    W = reshape(W, dropletcount, :)
    W = W'
    W = reshape(W, tsteps, :)
    W = W'
    W = reshape(W, dimensions, :)
    W = W'
    
    radvector = reshape(V[:,1:dropletcount]',1,:)';
    osmvector = reshape(V[:,dropletcount+1:2*dropletcount]',1,:)'
    
    W = [W radvector osmvector]
    # Now create a unified string vector 
    α = []
    ω = size(W,1)
    
    W = string.(W)
    for i = 1:ω
        push!(α,join(W[i,:]," "))
    end
    timetitle = ["ITEM: TIMESTEP"];
    timetitle = repeat(timetitle, Int(size(W,1)/dropletcount), 1);
    
    timenumber = convert(Array, sol.t)
    timenumber = string.(timenumber)
    
    dropnumtitle = ["ITEM: NUMBER OF ATOMS"]
    dropnumtitle = repeat(dropnumtitle,Int(size(W,1)/dropletcount),1 )
    
    dropnumnumber = [string.(dropletcount)]
    dropnumnumber = repeat(dropnumnumber, Int(size(W,1)/dropletcount), 1)
    
    
    typetitle = ["ITEM: ATOMS x y z radius osmolarity"]
    typetitle = repeat(typetitle, Int(size(W,1)/dropletcount), 1)
    
    assembly = permutedims([timetitle timenumber dropnumtitle dropnumnumber typetitle])
    assembly = reshape(assembly, :, 1);
    # timenumber dropnumtitle dropnumnumber typetitle]
    
    boolean_me = repeat([["boolean_me"]; repeat(["boolean_me2"],Int(dropletcount-1) ,1)],tsteps,5);
    W = α
    W = [W boolean_me]
    W = permutedims(W) # Seems good to this point
    W = reverse(W, dims = 1)
    W = permutedims(reshape(W,1,:))
    jj = W .== "boolean_me2"
    W= W[.~jj]
    boolean_me = [repeat(["boolean_me2"],4,1); ["boolean_me"]];
    assembly = [assembly repeat(boolean_me,tsteps,dropletcount)]
    assembly = permutedims(assembly)
    assembly = reshape(assembly,1,:)
    assembly = permutedims(assembly)
    
    kk = assembly .== "boolean_me2"
    assembly = assembly[.~kk]
    W = [W assembly]
    W = permutedims(W)
    W = reshape(W,1,:)
    ii = W .== "boolean_me"
    W = W[.~ii]
    
    return W
    
    end
        

# Module end here:
end






























