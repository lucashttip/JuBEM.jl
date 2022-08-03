function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, problem::problem_type, solver_var::solver_var_type)

    if !isempty(problem.nFr)
        calc_GH_dynamic_non_const!(mesh, material, solver_var, problem.frequency)
    else
        calc_GH_static_non_const!(mesh, material, solver_var)
    end

    return solver_var
end

function solvestatic(inp_file)
    mesh, material, problem, solver_var = read_msh(inp_file)

    derive_data!(material, problem, solver_var)

    generate_mesh!(mesh)

    calc_GH!(mesh, material, problem, solver_var)

    applyBC_nonrb3!(mesh, solver_var)

    x = solver_var.ma \ mesh.zbcvalue

    u,t = returnut3(mesh,x)

    return mesh, material, problem, solver_var, u, t
end