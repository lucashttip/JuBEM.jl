using Revise
using JuBEM
using Plots
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

inp_file = "meshes/dynamic/bars/bar_2_3.msh"

solve(inp_file)

mesh,material,problem,solver_var = readvars_out("output")

node = findfirst(x->x==10,mesh.nodes[:,2])

u, t, freqs = getnoderes_out("output",node)
# u,t = getfreqres_out("output",0)

writevtk(mesh,u,t,"test")

plot(freqs, abs.(real.(u[:,1])))