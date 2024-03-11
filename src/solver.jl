
# Interface to choose GH calculation type
function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type,frequency=-1.0)

    if frequency >= 0.0
        calc_GH_dynamic!(mesh, material, solver_var, frequency)
    else
        calc_GH_static!(mesh, material, solver_var)
    end

    return solver_var
end

function solvestatic(mesh, material, problem, solver_var;file_out = "output", savemat = false)

    calc_GH!(mesh, material, solver_var,-1.0)

    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end

    mesh, solver_var, C = applyBC(mesh, solver_var,solver_var.H,solver_var.G)
    solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
    u,t,urb = returnut(mesh,solver_var.zvetsol, C)

    if savemat
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end
    output_freq_h5(file_out,0,u,t, urb)

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
        println("Rodando para frequencia: ", frequency)
        calc_GH!(mesh, material, solver_var, frequency)
        mesh, solver_var, C = applyBC(mesh, solver_var,solver_var.zH,solver_var.zG)
        solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
        zu,zt,zurb = returnut(mesh,solver_var.zvetsol, C)

        output_freq_h5(file_out,frequency,zu,zt,zurb)
    end

end

# Main function
function solve(inp_file;file_out="output", savemat = false)
    
    t1 = time()

    # Read mesh and generate 
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)

    if savemat == false
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end

    # Choose solver
    if isempty(problem.nFr)
        solvestatic(mesh, material, problem, solver_var;file_out=file_out, savemat = savemat)
    else
        solvedynamic(mesh, material, problem, solver_var;file_out=file_out, savemat = savemat)
    end

    t2 = time()
    output_time(file_out,t2-t1,"totaltime")
end

function solve_flex_dyn(inp_file;file_out = "output", savemat = false)

    t1 = time()

    forces = Float64[
        0 1 0 0 0 0
        0 0 1 0 0 0
        1 0 0 0 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        0 0 0 1 0 0
    ]

    # Read mesh and generate 
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)

    if savemat == false
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end
    
    calc_GH!(mesh, material, solver_var,-1.0)
    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end
    if savemat
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end

    for frequency in problem.frequencies
        println("Rodando para frequencia: ", frequency)
        calc_GH!(mesh, material, solver_var, frequency)

        N = zeros(ComplexF64,6,6)
        @showprogress 1 "Calculating flexibilities..." for i in axes(forces,2)
            mesh.forces = zeros(6,1)
            mesh.forces[:,1] = forces[:,i]

            mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.zH,solver_var.zG)
            solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
            zu,zt,urb = JuBEM.returnut(mesh,solver_var.zvetsol, C)
        
            N[:,i] = urb[[3,1,2,6,4,5]]

        end

        output_freqflex_h5(file_out, frequency, N)
    end
    

    t2 = time()
    output_time(file_out,t2-t1,"totaltime")

end

function solve_flex_dyn2(inp_file;file_out = "output", savemat = false)

    t1 = time()

    forces = Float64[
        1 0
        0 0
        1 0
        0 1
        0 0
        0 1
    ]

    # Read mesh and generate 
    mesh, material, problem, solver_var = read_msh(inp_file)
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)

    if savemat == false
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end
    
    calc_GH!(mesh, material, solver_var,-1.0)
    if 0 in mesh.bc
        remove_EE!(mesh, solver_var)
    end
    if savemat
        output_vars_h5(file_out, mesh, problem, solver_var, material)
    end

    for frequency in problem.frequencies
        println("Rodando para frequencia: ", frequency)
        calc_GH!(mesh, material, solver_var, frequency)

        N = zeros(ComplexF64,6,6)
        @showprogress 1 "Calculating flexibilities..." for i in 1:2
            mesh.forces = zeros(6,1)
            mesh.forces[:,1] = forces[:,i]

            mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.zH,solver_var.zG)
            solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
            zu,zt,urb = JuBEM.returnut(mesh,solver_var.zvetsol, C)
        
            if i == 1
                N[1,1] = urb[3]
                N[2,2] = urb[1]
                N[3,3] = urb[1]
                N[6,2] = urb[5]
                N[5,3] = -urb[5]
            else
                N[4,4] = urb[6]
                N[5,5] = urb[4]
                N[6,6] = urb[4]
                N[3,5] = urb[2]
                N[2,6] = -urb[2]
            end

        end

        output_freqflex_h5(file_out, frequency, N)
    end
    
    t2 = time()
    output_time(file_out,t2-t1,"totaltime")

end