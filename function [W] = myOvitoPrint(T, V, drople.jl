function myOvitoPrint!(T, sol, dropletcount, dimensions)

# This function will recast the solution of our ordinary differential equation into a printed file that can be visualized in ovito.

display(sol)

end




# function [W] = myOvitoPrint(T, V, dropletcount, dimensions)
#     %now print out a list that can be put into particle visualization software
    
    
    
#     tsteps = size(V,1);
#     %create the boolean matrix that will be used for this transformation
    
#     omicron = ~eye(dimensions);
#     omicron = repmat(omicron,dropletcount,1);
#     omicron = reshape(omicron, dimensions, []);
#     omicron = omicron';
#     omicron = repmat(omicron,tsteps,1);
#     omicron = reshape(omicron,dimensions*dropletcount,[]);
    
#     %shave off the velocities for now
#     W = V(:,2*dropletcount+1:2*dropletcount + dimensions*dropletcount); %##ADJUST##
#     %perform a sequence of transformations to re-arrange W to what we want
#     Wtrue = [W(:,1:2); W(:,3:4) ; W(:,5:6)];
#     Wtransptrue = Wtrue';
#     W = W';
#     W = repmat(W,1,dimensions);
#     W = W(~omicron);
#     W = reshape(W,dropletcount,[]);
#     W = W';
#     W = reshape(W,tsteps,[]);
    
#     %check this part later with pen and paper
#     W = W';
    
#     W = reshape(W, dimensions, []);
    
#     W = W';
    
#     radvector = reshape(V(:,1:dropletcount)',1,[])';
#     osmvector = reshape(V(:,dropletcount+1:2*dropletcount)',1,[])';
    
#     W = [W, radvector,osmvector];
#     %now convert W to a string vector that can be printed in a text file
#     W = string(W);
#     W = join(W);
    
    
#     %assemble the string vector that provides information for each time step
#     timetitle = "ITEM: TIMESTEP";
#     timetitle = repmat(timetitle, length(W)/dropletcount, 1);
    
    
#     timenumber = (T);
#     timenumber = string(timenumber);
    
#     dropnumtitle = "ITEM: NUMBER OF ATOMS" ;
#     dropnumtitle = repmat(dropnumtitle,length(W)/dropletcount, 1);
    
#     dropnumnumber = string(dropletcount);
#     dropnumnumber = repmat(dropnumnumber,length(W)/dropletcount, 1);
    
#     typetitle = "ITEM: ATOMS x y z radius osmolarity";
#     typetitle = repmat(typetitle, length(W)/dropletcount, 1);
    
#     assembly = [timetitle, timenumber, dropnumtitle, dropnumnumber, typetitle]';
#     assembly = reshape(assembly,[],1);
    
    
#     %create a boolean matrix that will be used to concatenate the...,
#     %time step information vector to the actual droplet data
    
#     boolean_me = repmat(["boolean_me"; repmat("boolean_me2",dropletcount-1,1)],tsteps,5);
#     W = [W, boolean_me];
#     W = W';
#     W = flip(W);
#     W = reshape(W,1,[])' ;
#     jj = W == "boolean_me2";
#     W = W(~jj);
    
#     boolean_me = [repmat("boolean_me2",4,1); "boolean_me"];
#     assembly = [assembly, repmat(boolean_me,tsteps,dropletcount)];
#     assembly = assembly';
#     assembly = reshape(assembly,1,[]);
#     assembly = assembly';
#     kk = assembly == "boolean_me2";
#     assembly = assembly(~kk);
    
#     W = [W, assembly];
#     W = W' ;
#     W = reshape(W,1,[])';
#     ii = W == "boolean_me";
#     W = W(~ii);
    
    
    
#     end