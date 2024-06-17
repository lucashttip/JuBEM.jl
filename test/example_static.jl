using JuBEM

mesh_file = "./input/meshes/bar_1_10.msh"
problem_file = "./input/problems/bar.prob"
output_file = "test_out.h5"
##

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
assembly = derive_data!(materials,problem)

JuBEM.statics_assembly_tests!(mesh,materials,assembly)

LHS, RHS = JuBEM.applyBC_simple(mesh::Mesh,problem::Problem,assembly::Assembly,assembly.H,assembly.G)

x = LHS\RHS

u,t = JuBEM.returnut_simple(x,mesh,problem)

e = findall(mesh.tag.==3)
idx = vec(mesh.IEN[:,e])
ul = u[idx,:]

# solution = solve(Assembly)

# plot_disps(mesh,solution)