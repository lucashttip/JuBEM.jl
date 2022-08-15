using Revise
using JuBEM
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

inp_file = "meshes/dynamic/bars/bar_2_3.msh"

mesh, material, problem, solver_var = read_msh(inp_file)
derive_data!(material, problem, solver_var)
generate_mesh!(mesh)
calc_GH!(mesh, material, solver_var,-1.0)
output_vars_h5("output", mesh, problem, solver_var, material)

