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

function solvestatic(mesh, material, problem, solver_var;file_out = "output", savemat = false)

    calc_GH!(mesh, material, solver_var,-1.0)

    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end

    mesh, solver_var, C = applyBC(mesh, solver_var,solver_var.H,solver_var.G)
    solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
    u,t,urb = returnut(mesh,solver_var.zvetsol, C)


    # applyBC_nonrb3!(mesh, solver_var, solver_var.H, solver_var.G)
    # solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
    # u,t = returnut3(mesh,solver_var.zvetsol)

    if savemat
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end
    output_freq_h5(file_out,0,u,t, urb)

end

function solvedynamic(inp_file;file_out="output")

    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)
    solvedynamic(mesh, material, problem, solver_var;file_out=file_out)

end

function solvedynamic(mesh, material, problem, solver_var;file_out="output", savemat = false)

    calc_GH!(mesh, material, solver_var,-1.0)
    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end
    if savemat
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end

    for frequency in problem.frequencies
    # frequency = problem.frequencies[1]
        println("Rodando para frequencia: ", frequency)
        calc_GH!(mesh, material, solver_var, frequency)
        mesh, solver_var, C = applyBC(mesh, solver_var,solver_var.zH,solver_var.zG)
        solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
        zu,zt,zurb = returnut(mesh,solver_var.zvetsol, C)

        # applyBC_nonrb3!(mesh, solver_var, solver_var.zH, solver_var.zG)
        # solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
        # zu,zt = returnut3(mesh,solver_var.zvetsol)

        output_freq_h5(file_out,frequency,zu,zt,zurb)
    end

end

function solve(inp_file;file_out="output", savemat = false)
    t1 = time()
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)
    if savemat == false
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end

    if isempty(problem.nFr)
        solvestatic(mesh, material, problem, solver_var;file_out=file_out, savemat = savemat)
    else
        solvedynamic(mesh, material, problem, solver_var;file_out=file_out, savemat = savemat)
    end

    t2 = time()
    output_time(file_out,t2-t1,"totaltime")
end

# function solve_rb(inp_file;file_out="output"; savemat = savemat)
#     mesh, material, problem, solver_var = read_msh(inp_file)
#     derive_data!(material, problem, solver_var)
#     generate_mesh!(mesh)
#     calc_GH!(mesh, material, solver_var,-1.0)
#     if 0 in mesh.bc
#         remove_EE!(mesh, solver_var)
#     end
#     output_vars_h5(file_out, mesh, problem, solver_var, material)

#     for frequency in problem.frequencies
#     # frequency = problem.frequencies[1]
#         println("Rodando para frequencia: ", frequency)
#         calc_GH!(mesh, material, solver_var, frequency)
#         mesh, solver_var,C = applyBC_rb(mesh, solver_var, solver_var.zH, solver_var.zG)
#         solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
#         output_frb_h5(file_out,solver_var.zvetsol[1:6],frequency)
#     end
# end