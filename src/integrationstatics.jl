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

function integrate_sing_static_1(source_node,gauss_points_sing,N_sing,normal_sing,J_sing, omega_sing, delta, C_stat,n)

    nGP = length(omega_sing)
    nnel = size(N_sing,2)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, C_stat)
            
        P1 = J_sing[i]*omega_sing[i]
        for k in 1:nnel
            P = N_sing[i,k]*P1
            if k!=n
                HELEM[:,3*(k-1)+1:3*k] += t.*P
            end
            GELEM[:,3*(k-1)+1:3*k] += u.*P
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

function integrate_rigid_body2!(H,nnel)

    nnodes = size(H,1)/3

    for n in 1:nnodes

        i1 = Int(3*(n-1)+1)
        i2 = Int(3*(n-1)+2)
        i3 = Int(3*n)

        j1 = i1% (3*nnel)
        j2 = i2% (3*nnel)
        j3 = i3% (3*nnel)

        j1 == 0 ? j1 = 3*nnel : nothing
        j2 == 0 ? j2 = 3*nnel : nothing
        j3 == 0 ? j3 = 3*nnel : nothing


        idx1 = Int.(collect(j1:3*nnel:3*nnodes))
        idx2 = Int.(collect(j2:3*nnel:3*nnodes))
        idx3 = Int.(collect(j3:3*nnel:3*nnodes))


        filter!(a->a!=i1,idx1)
        filter!(a->a!=i2,idx2)
        filter!(a->a!=i3,idx3)
        

        H[i1,i1] = - sum(H[i1,idx1])
        H[i2,i1] = - sum(H[i2,idx1])
        H[i3,i1] = - sum(H[i3,idx1])
        H[i1,i2] = - sum(H[i1,idx2])
        H[i2,i2] = - sum(H[i2,idx2])
        H[i3,i2] = - sum(H[i3,idx2])
        H[i1,i3] = - sum(H[i1,idx3])
        H[i2,i3] = - sum(H[i2,idx3])
        H[i3,i3] = - sum(H[i3,idx3])

    end
    
    return H
end