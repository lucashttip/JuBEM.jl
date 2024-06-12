using JuBEM

mesh_file = "../input/meshes/bar_1_10.msh"
problem_file = "../input/problems/barproblem.prob"
output_file = "test_out.h5"

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

function solve(mesh,problem,materials)

    output_vars_h5(file_out, mesh, problem, material)

    static_ass = statics_assembly(mesh,problem,materials)

    remove_ee!(mesh,problem,static_ass)

    for freq in problem.frequencies

        assembly = dynamic_assembly(mesh,problem,materials,static_ass,freq)
        solution = solve(assembly)


        output_freq_h5(file_out,freq,solution)

    end

end

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
derive_data!(mesh,problem,materials)

solve(mesh,problem,materials)