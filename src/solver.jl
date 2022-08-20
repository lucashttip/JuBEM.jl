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

function solvestatic(inp_file;file_out="output")
    mesh, material, problem, solver_var = read_msh(inp_file)

    derive_data!(material, problem, solver_var)

    generate_mesh!(mesh)

    solvestatic(mesh, material, problem, solver_var;file_out=file_out)
end

function solvestatic(mesh, material, problem, solver_var;file_out = "output")

    calc_GH!(mesh, material, solver_var,-1.0)

    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end

    applyBC_nonrb3!(mesh, solver_var, solver_var.H, solver_var.G)

    solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue

    u,t = returnut3(mesh,solver_var.zvetsol)

    output_vars_h5(file_out, mesh, problem, solver_var, material)
    output_freq_h5(file_out,u,t,0)

end

function solvedynamic(inp_file;file_out="output")

    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)
    solvedynamic(mesh, material, problem, solver_var;file_out=file_out)

end

function solvedynamic(mesh, material, problem, solver_var;file_out="output")

    calc_GH!(mesh, material, solver_var,-1.0)
    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end
    output_vars_h5(file_out, mesh, problem, solver_var, material)

    for frequency in problem.frequencies
    # frequency = problem.frequencies[1]
        println("Rodando para frequencia: ", frequency)
        calc_GH!(mesh, material, solver_var, frequency)
        applyBC_nonrb3!(mesh, solver_var, solver_var.zH, solver_var.zG)
        solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
        zu,zt = returnut3(mesh,solver_var.zvetsol)

        output_freq_h5(file_out,zu,zt,frequency)
    end

end

function solve(inp_file;file_out="output")
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)

    if isempty(problem.nFr)
        solvestatic(mesh, material, problem, solver_var;file_out=file_out)
    else
        solvedynamic(mesh, material, problem, solver_var;file_out=file_out)
    end
end