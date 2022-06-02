function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, problem::problem_type, solver_var::solver_var_type)

    # calc_GH_dynamic_non_const!(mesh, material, solver_var, problem.frequency)

    calc_GH_static_non_const!(mesh, material, solver_var)

    return solver_var
end


