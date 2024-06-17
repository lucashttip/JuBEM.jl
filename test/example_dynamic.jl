using JuBEM

mesh_file = "./input/meshes/bar_1_10.msh"
problem_file = "./input/problems/bar_dynamic.prob"
output_file = "test_out.h5"

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
assembly = derive_data!(materials, problem)

output_vars_h5(output_file, mesh, problem, assembly, materials)

static_ass = statics_assembly(mesh,problem,materials)

# remove_ee!(mesh,problem,static_ass)

for freq in problem.frequencies

    assembly = dynamic_assembly(mesh,problem,materials,static_ass,freq)
    solution = solve(assembly)


    output_freq_h5(file_out,freq,solution)

end

