using Revise
using JuBEM

# mesh_file = "./input/meshes/homogbar_ne=1x10_l=10x1x1_eo=1.msh"
mesh_file = "./input/meshes/soilEE_109.msh"
# problem_file = "./input/problems/bar_static.prob"
problem_file = "./input/problems/soilrb_EE_static.prob"
# problem_file = "./input/problems/soil_EE_static.prob"


output_file = "test_static.h5"
##

# t1 = time()

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
derive_data!(materials,problem)

# output_vars(output_file, mesh, problem, materials)

assembly = JuBEM.statics_assembly(mesh,materials)

JuBEM.remove_EE!(mesh,assembly,problem)

LHS, RHS = JuBEM.applyBC_simple(mesh,problem,assembly.H,assembly.G)
LHS, RHS = JuBEM.applyBC_rb(mesh,problem,assembly.H,assembly.G)


x = LHS\RHS

u,t = JuBEM.returnut_simple(x,mesh,problem)
u,t,urb = JuBEM.returnut_rb(x,mesh,problem)

# t2 = time()

# sol = Solution(u,t,Float64[],t2-t1,0.0)
sol = Solution(u,t,urb,0.0,0.0)

# JuBEM.output_solution(output_file,sol)

# JuBEM.output_time(output_file,t2-t1,"totaltime")

## Pos-Processing

e = findall(mesh.tag.==3)
idx = vec(mesh.IEN[:,e])
ul = u[idx,:]

plot_disp(mesh,sol,3)

writevtk(mesh,sol,"teste")