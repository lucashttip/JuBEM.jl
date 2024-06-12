using JuBEM

mesh_file = "./input/meshes/bar_1_10.msh"
problem_file = "./input/problems/barproblem.prob"
output_file = "test_out.h5"
##

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
assembly = derive_data!(mesh,problem,materials)

# Assembly = statics_assembly(mesh,problem,materials,assembly)

# solution = solve(Assembly)

# plot_disps(mesh,solution)