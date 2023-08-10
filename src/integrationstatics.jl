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

function integrate_sing_static(source_node,gauss_points_sing,normal_sing, delta, C_stat,pesos,n)

    nGP = size(pesos,1)
    nnel = size(pesos,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)

    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, C_stat)
            
        for k in 1:nnel
            # @infiltrate
            if k!=n
                HELEM[:,3*(k-1)+1:3*k] += t.*pesos[i,k]
            end
            # if k == n
                GELEM[:,3*(k-1)+1:3*k] += u.*pesos[i,k]
            # end
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