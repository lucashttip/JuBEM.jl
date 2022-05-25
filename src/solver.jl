function calc_GH!(mesh::mesh_type, material::Vector{material_table_type}, problem::problem_type, solver_var::solver_var_type)

    calc_GH_dynamic!(mesh, material, solver_var, problem.frequency)

    return solver_var
end

function calc_GH_dynamic!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type, frequency::Float64)

    nelem = mesh.nelem
    m = 1
    calc_static_consts!()
    calc_complex_constants!()
    delta = I(3)
    zHselem = complex(delta)./2

    max_GL = maximum(mesh.ID)
    solver_var.H = complex(zeros(max_GL,max_GL))
    solver_var.G = complex(zeros(max_GL,max_GL))


    # Field loop:
    for fe in 1:nelem
    
        N = calc_N_matriz()
        dN = calc_dN_matriz()

        normal = calc_normal()

        # source loop:
        for se in 1:nelem
       
            if fe != se
                HELEM, GELEM = calc_nonsing()
            else
                HELEM, GELEM = calc_sing()
            end

            solver_var.H[mesh.LM[:,se], mesh.LM[:,fe]] = HELEM
            solver_var.G[mesh.LM[:,se], mesh.LM[:,fe]] = GELEM

        end
    
    end

end
