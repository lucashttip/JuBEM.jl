function calc_integration_rules!(rules::Rules,mesh::Mesh)
    
    # Normal integration
    pw = gausslegendre.(rules.npgs)
    for i in eachindex(rules.npgs)
        csis, omegas = calc_csis2d(pw[i][1], pw[i][2])
        push!(rules.gp,gauss_points_type(csis,omegas))
    end

    # Near integration (Differente number of points)
    pnear = gausslegendre(rules.npg_near)
    csis, omegas  = calc_csis2d(pnear[1], pnear[2])
    rules.gp_near = gauss_points_type(csis,omegas)

    # Singular integration
    psing = gausslegendre(rules.npg_sing)
    csis, omegas = calc_csis2d(psing[1],psing[2])

    N_lin = calc_N_gen([-1,1],csis)
    dNc_lin = calc_N_gen([-1,1],csis;dg=:dNdc)
    dNe_lin = calc_N_gen([-1,1],csis;dg=:dNde)

    if mesh.eltype == 2
        nperm = 3
    elseif mesh.eltype < 2
        nperm = 1
    end

    csi_sing, Jb_sing = csis_sing(mesh.offset, N_lin, dNc_lin,dNe_lin,mesh.eltype)
    
    nregs = Int(length(Jb_sing[:,1])/length(omegas))
    omega_sing = repeat(omegas,nregs)

    for i in 1:nperm
        push!(rules.csis_sing,gauss_points_type(csi_sing[:,:,i],Jb_sing[:,i].*omega_sing))
    end

    return rules
end

"""
    rules, csis_cont, csis_disc = define_rules(mesh::Mesh)
"""
function define_rules(mesh::Mesh)

    npgs = Int[4,5,6,8]
    dists = [4,2,0.5,0.2]
    npg_near = 12
    npg_sing = 8
    near_strat = :telles

    eltype_cont = mesh.eltype
    offset = mesh.offset
    if mesh.eltype == 0
        eltype_cont = 1
        offset = 1.0
    end

    csis_points = collect(range(-1,1,length = eltype_cont+1))
    csis_nodes = collect(range(-1+offset,1 - offset,length=mesh.eltype+1))


    csis_element_size = [0.0 -1.0
                    1.0 0.0
                    0.0 1.0
                    -1.0 0.0]

    N_aux = calc_N_gen(csis_points,csis_element_size)

    
    rules = Rules(csis_points,csis_nodes,npgs,dists,npg_near,npg_sing,near_strat,N_aux)
    calc_integration_rules!(rules,mesh)

    return rules

end