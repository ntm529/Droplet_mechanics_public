module myODEcomponents
export myDropvecDistances, myDropvecRadius_Osmolarity
using LinearAlgebra

boolean_osm = false

function myDropvecDistances(droplet_locations, dropletcount, boolean_distances, boolean_norms, radii, L, dimensions, comradii, bilayer_spring_distance, norm_combinations, distance_combinations, direction, radii_combinations, logic_spring)

   

    for l in 0:dropletcount:2*dropletcount # The '2' value is the total number of dimensions (3) minus one. Needed for iteration to work
        for i in 1: dropletcount
            k = 1 
            for j in 1:(i-1)
                distance_combinations[i+ l,k] = (droplet_locations[j+l])         

            k = k +  1
        end
            for j in i+1:dropletcount
                distance_combinations[i+ l,k] = (droplet_locations[j+l])        
                k = k +  1
            end
        end  
    end
    @.distance_combinations = distance_combinations .- droplet_locations;

    # Get the norm distance between all droplets

    # display(norm_combinations)
    for i in 1:dropletcount-1
        for j in 1:dropletcount
            # zed = view(distance_combinations,j:dropletcount:size(distance_combinations,1),i)
            # norm_combinations[j,i] = norm(zed, 2)
            norm_combinations[j,i] = norm(view(distance_combinations,j:dropletcount:size(distance_combinations,1),i), 2)
            norm_combinations[j+dropletcount,i] = norm_combinations[j,i]     
            norm_combinations[j+2*dropletcount,i] = norm_combinations[j,i]   
            
        end
    end
    
  
    # This line of code might be necessary in the future
    # norm_combinations[isnan.(norm_combinations)] .= 0

    # ℶ = view(norm_combinations, 1:dropletcount, :)
    # display(ℶ)
    #Get the normalized directional vector (ranges between 0 and 1) for all droplets

    @.direction = distance_combinations./norm_combinations;

   
    


    for i in 1: length(radii)
         k = 1 
         for j in 1:(i-1)
            comradii[i,k] = (radii[j])         
            
                      
        
            bilayer_spring_distance[i,k] = L * (radii[j] + radii[i])          
            bilayer_spring_distance[i+dropletcount,k] = bilayer_spring_distance[i,k]      
            bilayer_spring_distance[i+2*dropletcount,k] = bilayer_spring_distance[i,k]  
                

            

            k = k +  1
        end
         for j in i+1:length(radii)
            comradii[i,k] = (radii[j])        

                     
            bilayer_spring_distance[i,k] = L * (radii[j] + radii[i])          
            bilayer_spring_distance[i+dropletcount,k] = bilayer_spring_distance[i,k]      
            bilayer_spring_distance[i+2*dropletcount,k] = bilayer_spring_distance[i,k]   

    
       
            k = k +  1
        end
    end  
    
    @.radii_combinations = comradii.+radii

    ℸ = view(norm_combinations, 1:dropletcount, :)

    @.logic_spring[1:dropletcount, :] = (radii_combinations .-  broadcast(abs, ℸ))  .>= 0 
    logic_spring[dropletcount+1: 2* dropletcount, :] = logic_spring[1:dropletcount, :]
    logic_spring[2*dropletcount+1: 3* dropletcount, :] = logic_spring[1:dropletcount, :]
  
    return comradii, distance_combinations, norm_combinations, direction, logic_spring , bilayer_spring_distance

end # <<<<< End of function
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function myDropvecRadius_Osmolarity(comradii, radii, osmolarities, D, dropletcount, norm_combinations, normvelocity, boolean_osm, t, sphere_intersection, circular_area, comosm, J, volume, moles, dradii, dosmolarities, sphere_dist_component_1, sphere_dist_component_2)
    if (t > 0.1 && normvelocity/dropletcount <= 1e-4)  
        boolean_osm[1] = 1.0
    end

    if norm(boolean_osm) < 1.0 
        fill!(dradii,0.0)
        fill!(dosmolarities,0.0)

    else

    normdists = view(norm_combinations, 1:dropletcount, :)
    # display(ℶ)
    # display(radii)
    # display(comradii)
    # The sphere intersection evaluation should first be turned into a logical matrix, and then be operated on to get numerical results for circular intersections
    


     #CASE 1: Simplest case where near circle is outside of the radius of,...
    # far circle
    @.sphere_intersection = (radii .+ comradii.-      normdists) .>= 0 ;

    #CASE 2: near circle is fully inside radius of the far circle
    @.sphere_intersection = sphere_intersection .* ((normdists .+ radii .- comradii) .>= 0 );

    #CASE 3: far circle is completely inside the near circle
    @.sphere_intersection = sphere_intersection .* ((normdists .+ comradii .- radii) .>= 0 );

  
   

    # Despite this, there are some times very small values below 

    # pre-allocate this matrix
    # sphere_intersection = zeros(n)
    #@. will make julia try to devectorize line
    # Initialize sphere_intersection -------> use the devectorizing @. to get it to work.

    #Now with the logical values set up, show the distance at which the sphere intersection is occuring


    @.sphere_dist_component_1 = ((normdists.^2-comradii.^2 .+ radii.^2).^2)
  

    @.sphere_dist_component_2 = ( 4*normdists.^2 .*radii.^2 )
    @.sphere_intersection = sphere_intersection .* ( sphere_dist_component_2 - sphere_dist_component_1)




    # bilayer_spring_distance = (comradii.+radii);

    @.sphere_intersection = 1 ./(2*normdists) .* sqrt.(sphere_intersection); #radius of intersection
    
    @.circular_area = pi * sphere_intersection.^2; #circular area of intersection


    # @.circular_area[isnan.(circular_area)] .= 0; #dealing with case where circular area is NaN <---- Might be useful, but I think the code works fine...
    # ... Without this line of code, IE, use it if an issue crops up but the other lines before this probably eliminate the case that could lead to NaN in code


    for i in 1: length(osmolarities)
        k = 1 
        for j in 1:(i-1)
           comosm[i,k] = (osmolarities[j])         
           
           k = k +  1
       end
        for j in i+1:length(osmolarities)
           comosm[i,k] = (osmolarities[j])       

           k = k +  1
       end
   end  
    @.comosm = -1 *(comosm .- osmolarities)

    # display(circular_area)
    # display(comosm)

    @.J = D * circular_area .* comosm;


    # display(volume)
    @.volume = 4/3 * pi * radii.^3;
    # display(volume)
    volumeflux = sum(J, dims= 2)


    # display(moles)
    @.moles = osmolarities .* volume;
    
    @.volume = volume .+  volumeflux #sum(J, dims= 2); # Julia won't let me put the @. macro here when I just use the function sum(J, dims= 2)
    # display(volume)
    
    @.dradii =  (3/(4*pi) * volume).^(1/3) - radii
    @.dosmolarities = moles./volume - osmolarities

    #add dmass here and make it one the exports 





end

    return  dradii, dosmolarities



end # <<<<< End of function


end # <<<<< End of module 













