using Revise
using JuBEM
using Plots
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

inp_file = "meshes/static/bars/bar_2_3.msh"
# inp_file = "meshes/static/bars/bar_1_1.msh"

solve(inp_file)
mesh,material,problem,solver_var = readvars_out("output")
u,t = getfreqres_out("output",0)

solver_var.H

node = findfirst(x->x==10,mesh.nodes[:,2])
node = 1

u, t, freqs = getnoderes_out("output",node)
u,t = getfreqres_out("output",freqs[3])
u,t = getfreqres_out("output",0)


writevtk(mesh,u,t,"test")

xmin = minimum(abs.(mesh.nodes[:,2]))
nodes = mesh.nodes[findall(x->x≈xmin,mesh.nodes[:,2]),:]
idx = Int.(nodes[sortperm(nodes[:,3]),1])
plot(mesh.nodes[idx,3],u[idx,3])
plot!(mesh.nodes[idx,3],real.(u_din[idx,3]))


plot(freqs, abs.(real.(u[:,2])))
freq = freqs[3]
animate_res_freq(mesh,u,freq;frac = 2.0, filename = "anim.mp4",res = (1920, 1080))

##

mesh, material, problem, solver_var = read_msh(inp_file)
derive_data!(material, problem, solver_var)
generate_mesh!(mesh)
calc_GH!(mesh, material, solver_var,-1.0)

frequency = problem.frequencies[1]; calc_GH!(mesh, material, solver_var, frequency)
