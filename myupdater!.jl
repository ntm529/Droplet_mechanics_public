using LinearAlgebra
using DifferentialEquations
using Plots

# f(u,p,t) = 1.01*u
# u0 = 1/2



# tspan = (0.0,10.0)
# prob = ODEProblem(f,u0,tspan)
# sol = solve(prob)
# display(sol)
# plot(sol, vars = 1)


##

a = Int64[]

for i = 1:7
    
    push!(a,i)
    b = a[end] + 1 
end

display(a)
display(b)

# push!(a,1)
# push!(a,2)
# display(a)