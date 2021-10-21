# A = [1 2 3 4 ; 11 10 12 18 ; 21 22 23 24]
# string.(A)

# B = [1 2]
# C = [3 4]
# string.(B)
# join(B)
# D = [B; C]
# E = join.(D[1,:])
# E'

# str = ["Carlos" "Sada" ; "Ella" "Olsen"]
# permutedims(str)
# outfile = "random_script.txt"
# f = open(outfile, "w")

# for i in eachindex(str)   
#     println(f, str[i])
# end

# close(f)

# booleantest = false
x = 0 
boolean_test = false
function myfn(x, boolean_test)
    for i = 1:20
        x += 1 
        if x>10
            boolean_test = true
        end

        if boolean_test == true
             println(x, boolean_test)
        end
        
    end
end

myfn(x, boolean_test)