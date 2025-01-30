using Revise
using JuBEM

# mesh_file = "./input/meshes/homogbar_ne=1x10_l=10x1x1_eo=1.msh"
# mesh_file = "./input/meshes/multibar_ne=1x(2+2)_l=(5+5)x1x1_eo=1.msh"
mesh_file = "./input/meshes/SoilLayers_a=1.0_fl=2.0_fs=8.0_h=4,nf=4_eo=1.msh"

# mesh_file = "./input/meshes/soilEE_109.msh"

problem_file="./input/problems/soilrb_layer_EE_static.prob"
# problem_file = "./input/problems/multi_bar.prob"
# problem_file = "./input/problems/bar_static.prob"
# problem_file = "./input/problems/soilrb_EE_static.prob"
# problem_file = "./input/problems/soil_EE_static.prob"


output_file = "test_static.h5"
##

# t1 = time()

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
derive_data!(mesh,materials,problem)

# plot_meshtags(mesh)

assembly = JuBEM.statics_assembly(mesh,materials,problem)
# output_vars(output_file, mesh, problem, materials,assembly)

JuBEM.remove_EE!(mesh,assembly,problem)

LHS, RHS = applyBC_simple(mesh,problem,assembly.H,assembly.G)
LHS, RHS = applyBC_rb(mesh,problem,assembly.H,assembly.G)
LHS,RHS = JuBEM.applyBC_multi_simple(mesh,problem,assembly.H,assembly.G)
LHS,RHS = JuBEM.applyBC_rb_multi(mesh,problem,assembly.H,assembly.G)


x = LHS\RHS

u2,t2 = JuBEM.returnut_simple(x,mesh,problem)
u,t,urb = JuBEM.returnut_rb(x,mesh,problem)
u,t,urb = JuBEM.returnut_rb_multi(x,mesh,problem)
u,t = JuBEM.returnut_multi_simple(x,mesh,problem)
# t2 = time()

sol = Solution(u,t)
sol = Solution(u,t,urb,0.0,0.0)


plot_disp(mesh,sol,3)

# JuBEM.output_solution(output_file,sol)

# JuBEM.output_time(output_file,t2-t1,"totaltime")

## Pos-Processing

e = findall(mesh.tag.==3)
idx = vec(mesh.IEN[:,e])
ul = u[idx,:]

plot_disp(mesh,sol,3)

writevtk(mesh,sol,"teste")