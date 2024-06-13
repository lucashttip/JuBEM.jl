using JuBEM

mesh_file = "./input/meshes/cube.msh"
problem_file = "./input/problems/barproblem.prob"
output_file = "test_out.h5"
##

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
assembly = derive_data!(materials,problem)

JuBEM.statics_assembly_tests!(mesh,materials,assembly)

# solution = solve(Assembly)

# plot_disps(mesh,solution)