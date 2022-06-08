function integrate_nonsing_static(source_node,gauss_points,N,normal,J, omegas, delta, C_stat)
    nGP = length(omegas)
    nnel = size(N,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
            
        P1 = J[i]*omegas[i]
        for k in 1:nnel
            P = N[i,k]*P1
            HELEM[:,3*(k-1)+1:3*k] += t*P
            GELEM[:,3*(k-1)+1:3*k] += u*P
        end
        
    end
    return HELEM, GELEM
end

function integrate_sing_static_1(source_node,gauss_points,N,normal,J, omegas, delta, C_stat,n)

    nGP = length(omegas)
    nnel = size(N,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
            
        P1 = J[i]*omegas[i]
        for k in 1:nnel
            P = N[i,k]*P1
            if k!=n
                HELEM[:,3*(k-1)+1:3*k] += t*P
            end
            GELEM[:,3*(k-1)+1:3*k] += u*P
        end
    end
    return HELEM, GELEM

end

function integrate_rigid_body!(H,nnel)

    ngdl = size(H,1)

    for i in 1:ngdl
        j = i% (3*nnel)
        j == 0 ? j = 3*nnel : nothing

        idx = collect(j:3*nnel:ngdl)

        filter!(a->a!=i,idx)

        H[i,i] = - sum(H[i,idx])
    end
    return H
end