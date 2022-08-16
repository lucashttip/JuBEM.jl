using Revise
using JuBEM
using Plots
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

inp_file = "meshes/static/soils/soilEE.msh"

mesh, material, problem, solver_var = read_msh(inp_file)
derive_data!(material, problem, solver_var)
generate_mesh!(mesh)
calc_GH!(mesh, material, solver_var,-1.0)


solve(inp_file)

mesh,material,problem,solver_var = readvars_out("output")

node = findfirst(x->x==10,mesh.nodes[:,2])

u, t, freqs = getnoderes_out("output",node)
# u,t = getfreqres_out("output",0)

writevtk(mesh,u,t,"test")

plot(freqs, abs.(real.(u[:,1])))