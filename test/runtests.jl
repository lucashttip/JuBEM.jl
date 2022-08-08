using Revise
using JuBEM
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

# inp_file = "meshes/vigas/viga_4_15.msh"
inp_file = "meshes/vigas/viga_2_5.msh"

mesh, material, problem, solver_var = read_msh(inp_file)

derive_data!(material, problem, solver_var)

generate_mesh!(mesh)

calc_GH!(mesh, material, problem, solver_var)

applyBC_nonrb3!(mesh, solver_var)

x = solver_var.ma \ mesh.zbcvalue

u,t = returnut3(mesh,x)
# maximum(u[:,1])
# minimum(u[:,2])

function solve_lin_viga(A,mesh)
    x = A\mesh.zbcvalue
    u,t = returnut3(mesh,x)

    ut = u[mesh.IEN[:,mesh.bc[:,2].==2][:],:]

    neu = sum(mesh.bc[:,2].==1)
    nnel = (mesh.eltype+1)^2
    ud = sum(ut[1:nnel*neu,2])./(nnel*neu)

    return u, t, ud
end