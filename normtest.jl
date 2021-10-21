using LinearAlgebra
droplets = 3 

# ζ = zeros(9,2)
# display(ζ)

# A = [ ones(9,1) 2 * ones(9,1) ]
# display(A)

# zed = view(A,1:droplets:size(A,1),1)
# C =norm(zed,2)
# display(C)
# B = sqrt(1+1+1)
# display(B)



using LinearAlgebra
droplets = 3 

ζ = zeros(droplets,droplets-1)
display(ζ)

A = [ ones(9,1) 2 * ones(9,1) ]
display(A)

# i = 1
#     for j in 1:droplets
#         zed = view(A,j:droplets:size(A,1),i)
#         C =norm(zed,2)
#         ζ[j,i] = C
#     end

# i = 2
#     for j in 1:droplets
#         zed = view(A,j:droplets:size(A,1),i)
#         C =norm(zed,2)
#         ζ[j,i] = C
#     end

# display(ζ)




for i in 1:droplets-1
    for j in 1:droplets
        zed = view(A,j:droplets:size(A,1),i)
        C =norm(zed,2)
        ζ[j,i] = C
    end

end

display(ζ)