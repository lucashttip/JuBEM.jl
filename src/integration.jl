## Auxiliary

function calc_nJp(r, rules, points)

    nintegpoints = length(rules.gp[r].omegas)
    normals = zeros(nintegpoints,3)
    weights = zeros(nintegpoints)

    N = calc_N_gen(rules.csis_points,rules.gp[r].csis)
    dNc = calc_N_gen(rules.csis_points,rules.gp[r].csis;dg=:dNdc)
    dNe = calc_N_gen(rules.csis_points,rules.gp[r].csis;dg=:dNde)

    integ_points = N*points
    dpdcsi = dNc*points
    dpdeta = dNe*points

    for i in 1:nintegpoints
        v = cross(dpdcsi[i,:],dpdeta[i,:])
        J = norm(v)
        normals[i,:] = v./J
        # Verificar se na linha de baixo estão sendo multiplicadas as mesmas coisas
        weights[i] = rules.gp[r].omegas[i]*J
    end

    return integ_points, normals, weights

end

function calc_nJp_sing(rules,points,idx_cont, np)

    nintegpoints = length(rules.csis_sing[np].omegas)

    normals = zeros(nintegpoints,3)
    weights = zeros(nintegpoints)

    N = calc_N_gen(rules.csis_points,rules.csis_sing[np].csis)
    dNc = calc_N_gen(rules.csis_points,rules.csis_sing[np].csis;dg=:dNdc)
    dNe = calc_N_gen(rules.csis_points,rules.csis_sing[np].csis;dg=:dNde)


    integ_points = N[:,idx_cont]*points
    
    dpdcsi = dNc[:,idx_cont]*points
    dpdeta = dNe[:,idx_cont]*points

    for i in 1:nintegpoints
        v = cross(dpdcsi[i,:],dpdeta[i,:])
        J = norm(v)
        normals[i,:] = v./J
        # Verificar se na linha de baixo estão sendo multiplicadas as mesmas coisas
        weights[i] = rules.csis_sing[np].omegas[i]*J
    end

    return integ_points, normals, weights



end

function calc_nearpoints(csis, omegas,c, d, csis_cont, csis_descont,field_points)

    gp_local_near, weights_near = pontos_pesos_local_telles(csis, omegas, c[1:2],d)
    # gp_local_near, weights_near = points_weights_local_near_combined(csis, omegas, c,minimum(d))

    # gp_local_near, weights_near = pontos_pesos_local_subelem(c[1], c[2], Nc_lin, dNcdcsi_lin, dNcdeta_lin, omegas)
    N_near = calc_N_gen(csis_cont,gp_local_near)
    dNc_near = calc_N_gen(csis_cont,gp_local_near;dg=:dNdc)
    dNe_near = calc_N_gen(csis_cont,gp_local_near;dg=:dNde)
    Nd_near = calc_N_gen(csis_descont,gp_local_near)
    gauss_points_near = N_near*field_points
    normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)

    weights_near = J_near.*weights_near

    return gauss_points_near, Nd_near, normal_near, weights_near
end

## Statics

function integrate_nonsing(source, points, rules, material,normalsign)
    r, c, d = calc_dist(source, points, rules)

    if r > 0
        # NORMAL INTEGRATION
        integ_points, normals, weights = calc_nJp(r, rules, 
        points)
        Nd = calc_N_gen(rules.csis_nodes,rules.gp[r].csis)
    else
        # NEAR INTEGRATION
        integ_points, Nd, normals, weights = calc_nearpoints(rules.gp_near.csis, rules.gp_near.omegas,c, d, rules.csis_points, rules.csis_nodes, points)
    end

    normals = normals.*normalsign
    nnel = size(Nd,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)

    u,t = calc_funsol_static_vec(source,integ_points,normals,material.C_stat)

    for k in axes(Nd,2)
        for i in 1:3
            for j in 1:3
                HELEM[i,3*(k-1)+j] = sum(t[i,j,:].*Nd[:,k].*weights)
                GELEM[i,3*(k-1)+j] = sum(u[i,j,:].*Nd[:,k].*weights)
            end
        end
    end

    return HELEM, GELEM

end

function integrate_nonsing2(source, points, rules, material,normalsign)
    r, c, d = calc_dist(source, points, rules)

    if r > 0
        # NORMAL INTEGRATION
        integ_points, normals, weights = calc_nJp(r, rules, 
        points)
        Nd = calc_N_gen(rules.csis_nodes,rules.gp[r].csis)
    else
        # NEAR INTEGRATION
        integ_points, Nd, normals, weights = calc_nearpoints(rules.gp_near.csis, rules.gp_near.omegas,c, d, rules.csis_points, rules.csis_nodes, points)
    end

    normals = normals.*normalsign
    nGP = length(weights)
    nnel = size(Nd,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)

    for i in 1:nGP
        u, t = calc_funsol_static(source,integ_points[i,:],normals[i,:], material.C_stat)
            
        # @infiltrate
        for k in 1:nnel
            P = Nd[i,k]*weights[i]
            HELEM[:,3*(k-1)+1:3*k] += t.*P
            GELEM[:,3*(k-1)+1:3*k] += u.*P
        end
        
    end
    return HELEM, GELEM

end

function integrate_sing(source, points, rules::Rules, material,normalsign, idx)


    nnel = Int(length(rules.csis_nodes)^2)

    idx_cont, idx_descont, np = calc_idx_permutation(nnel,idx)

    integ_points, normals, weights = calc_nJp_sing(rules,points,idx_cont, np)

    normals = normals.*normalsign

    Nd = calc_N_gen(rules.csis_nodes,rules.csis_sing[np].csis)[:,idx_descont]

    nnel = size(Nd,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)

    u,t = calc_funsol_static_vec(source,integ_points,normals,material.C_stat)

    for k in axes(Nd,2)
        for i in 1:3
            for j in 1:3
                pesos = Nd[:,k].*weights
                if k != idx
                    HELEM[i,3*(k-1)+j] = sum(t[i,j,:].*pesos)
                end
                GELEM[i,3*(k-1)+j] = sum(u[i,j,:].*pesos)
            end
        end
    end

    return HELEM, GELEM
end

## Dynamics

function integrate_nonsing_dyn(source, points, rules, zconsts,normalsign)
    r, c, d = calc_dist(source, points, rules)

    if r > 0
        # NORMAL INTEGRATION
        integ_points, normals, weights = calc_nJp(r, rules, 
        points)
        Nd = calc_N_gen(rules.csis_nodes,rules.gp[r].csis)
    else
        # NEAR INTEGRATION
        integ_points, Nd, normals, weights = calc_nearpoints(rules.gp_near.csis, rules.gp_near.omegas,c, d, rules.csis_points, rules.csis_nodes, points)
    end

    normals = normals.*normalsign

    nGP = length(weights)
    nnel = size(Nd,2)
    zHELEM = zeros(ComplexF64,3,3*nnel)
    zGELEM = zeros(ComplexF64,3,3*nnel)

    for i in 1:nGP

        zu, zt = calc_funsol_dynamic(source,integ_points[i,:],normals[i,:], zconsts)
            
        # @infiltrate
        for k in 1:nnel
            P = Nd[i,k]*weights[i]
            zHELEM[:,3*(k-1)+1:3*k] += zt.*P
            zGELEM[:,3*(k-1)+1:3*k] += zu.*P
        end
        
    end
    return zHELEM, zGELEM

end

function integrate_sing_dyn(source, points, rules::Rules, material,zconsts,normalsign, idx)


    nnel = Int(length(rules.csis_nodes)^2)

    idx_cont, idx_descont, np = calc_idx_permutation(nnel,idx)

    integ_points, normals, weights = calc_nJp_sing(rules,points,idx_cont, np)

    normals = normals.*normalsign

    Nd = calc_N_gen(rules.csis_nodes,rules.csis_sing[np].csis)[:,idx_descont]

    nGP = length(weights)
    nnel = size(Nd,2)
    zHELEM = zeros(ComplexF64,3,3*nnel)
    zGELEM = zeros(ComplexF64,3,3*nnel)


    for i in 1:nGP

        _, t = calc_funsol_static(source,integ_points[i,:],normals[i,:], material.C_stat)
        zu, zt = calc_funsol_dynamic(source,integ_points[i,:],normals[i,:], zconsts)

        for k in 1:nnel
            P = Nd[i,k]*weights[i]
            if k != idx
                zHELEM[:,3*(k-1)+1:3*k] += zt.*P
            else
                zHELEM[:,3*(idx-1)+1:3*idx] += (zt - complex(t)).*P
            end
            zGELEM[:,3*(k-1)+1:3*k] += zu.*P
        end

    end
    return zHELEM, zGELEM
end

function integrate_rigid_body!(H,mesh)

    nnel = size(mesh.IEN,1)
    for e in 1:mesh.nelem

        for n in 1:nnel

            i = mesh.LM[3*(n-1)+1:3*n,e]
        
            H[i,i[1]] = - sum(H[i,1:3:end],dims=2)
            H[i,i[2]] = - sum(H[i,2:3:end],dims=2)
            H[i,i[3]] = - sum(H[i,3:3:end],dims=2)
        end

    end
    
    return H
end