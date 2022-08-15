function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type,frequency=-1.0)

    if frequency >= 0.0
        if mesh.eltype >0
            calc_GH_dynamic_non_const!(mesh, material, solver_var, frequency)
        else
            calc_GH_dynamic_const!(mesh, material, solver_var, frequency)
        end
    else
        if mesh.eltype >0
            calc_GH_static_non_const!(mesh, material, solver_var)
        else
            calc_GH_static_const!(mesh, material, solver_var)
        end
    end

    return solver_var
end

function solvestatic(inp_file)
    mesh, material, problem, solver_var = read_msh(inp_file)

    derive_data!(material, problem, solver_var)

    generate_mesh!(mesh)

    calc_GH!(mesh, material, solver_var,-1.0)

    applyBC_nonrb3!(mesh, solver_var)

    x = solver_var.ma \ mesh.zbcvalue

    u,t = returnut3(mesh,x)

    return mesh, material, problem, solver_var, u, t
end

function solvedynamic(inp_file)
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)

    calc_GH!(mesh, material, solver_var,-1.0)
    output_vars_h5("output", mesh, problem, solver_var, material)


end