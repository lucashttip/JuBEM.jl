function calc_nonsing_static(source_node,gauss_points,N,normal,J, omega, delta, C_stat)
    nGP = length(omega)
    nnel = size(N,3)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)
    for i in 1:nGP
        for j in 1:nGP

            u, t = calc_funsol_static(source_node,gauss_points[i,j,:], normal[i,j,:], delta, C_stat)
            
            P1 = J[i,j]*omega[i]*omega[j]
            for k in 1:nnel
                P = N[i,j,k]*P1
                HELEM[:,3*(k-1)+1:3*k] += t*P
                GELEM[:,3*(k-1)+1:3*k] += u*P
            end
        end
    end
    return HELEM, GELEM
end

function calc_sing_static_1(source_node,gauss_points,Nd,normal,J, omega, delta, C_stat)



end