function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, problem::problem_type, solver_var::solver_var_type)

    if !isempty(problem.nFr)
        calc_GH_dynamic_non_const!(mesh, material, solver_var, problem.frequency)
    else
        calc_GH_static_non_const!(mesh, material, solver_var)
    end

    return solver_var
end


function apply_bc!(mesh,solver_var)
    
end