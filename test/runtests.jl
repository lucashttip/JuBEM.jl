using Revise
using JuBEM
using Plots
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

inp_file = "meshes/dynamic/soils/soilEE_216_rb.msh"
# inp_file = "meshes/static/bars/bar_2_3.msh"

file_out = "output_rb_lin"
# JuBEM.solve_rb(inp_file;file_out=file_out)

solve(inp_file;file_out=file_out)
mesh,material,problem,solver_var = readvars_out(file_out)

using Statistics, LinearAlgebra
a = mean(diag(solver_var.H))

u,t = getfreqres_out(file_out,0)

points_int = [9 0.5 0.5
    5 0.5 0.5]

ui, si = JuBEM.calc_interior_static(mesh,material,u,t,points_int)

solver_var.H

node = findfirst(x->x==10,mesh.nodes[:,2])
node = 1

u, t, freqs = getnoderes_out(file_out,node)
u,t = getfreqres_out(file_out,freqs[3])
u,t = getfreqres_out(file_out,0)

phiyFx_rb,freqs = geturb_out(file_out,5)

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
mesh.eltype=1
derive_data!(material, problem, solver_var)
generate_mesh!(mesh)
calc_GH!(mesh, material, solver_var,-1.0)

JuBEM.remove_EE!(mesh, solver_var)
# frequency = problem.frequencies[1]
frequency = 0.0001
calc_GH!(mesh, material, solver_var, frequency)

mesh, solver_var, C = JuBEM.applyBC_rb(mesh, solver_var,solver_var.H,solver_var.G)
solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
u,t,urb = JuBEM.returnut_rb(mesh,solver_var.zvetsol, C)


mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.zH,solver_var.zG)
zma = solver_var.ma
zbc = mesh.zbcvalue
solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
zu,zt,zurb = JuBEM.returnut(mesh,solver_var.zvetsol, C)

mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.H,solver_var.G)
ma = solver_var.ma
bc = mesh.zbcvalue
solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
u,t,urb = JuBEM.returnut(mesh,solver_var.zvetsol, C)