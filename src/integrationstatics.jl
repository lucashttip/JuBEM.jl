function integrate_nonsing_static(source_node,gauss_points,N,normal,J, omegas, delta, C_stat)
    nGP = length(omegas)
    nnel = size(N,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
            
        P1 = J[i]*omegas[i]
        # @infiltrate
        for k in 1:nnel
            P = N[i,k]*P1
            HELEM[:,3*(k-1)+1:3*k] += t.*P
            GELEM[:,3*(k-1)+1:3*k] += u.*P
        end
        
    end
    return HELEM, GELEM
end

function integrate_nonsing_static2(source_node,gauss_points,N,normal,J, omegas, delta, C_stat,n)
    nGP = length(omegas)
    nnel = size(N,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    for i in 1:nGP
        u, _ = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
            
        P1 = J[i]*omegas[i]
        # @infiltrate
        for k in 1:nnel
            P = N[i,k]*P1
            if k != n
                GELEM[:,3*(k-1)+1:3*k] += u.*P
            end
        end
        
    end
    return GELEM
end

function integrate_sing_static(source_node,gauss_points_sing,N_sing,normal_sing,J_sing, omega_sing, delta, C_stat,n)

    nGP = length(omega_sing)
    nnel = size(N_sing,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, C_stat)
            
        P1 = J_sing[i]*omega_sing[i]
        for k in 1:nnel
            P = N_sing[i,k]*P1
            # @infiltrate
            if k!=n
                HELEM[:,3*(k-1)+1:3*k] += t.*P
            end
            if k == n
                GELEM[:,3*(k-1)+1:3*k] += u.*P
            end
        end
    end
    return HELEM, GELEM

end

function integrate_rigid_body!(H,mesh)

    for n in 1:mesh.nnodes

        i = mesh.ID[:,n]
       
        H[i,i[1]] = - sum(H[i,1:3:end],dims=2)
        H[i,i[2]] = - sum(H[i,2:3:end],dims=2)
        H[i,i[3]] = - sum(H[i,3:3:end],dims=2)

    end
    
    return H
end
