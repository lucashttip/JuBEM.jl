
function integrate_nonsing_dynamic(source_node,gauss_points,N,normal,J, omegas, delta, zconsts)
    nGP = length(omegas)
    nnel = size(N,2)
    zHELEM = zeros(ComplexF64,3,3*nnel)
    zGELEM = zeros(ComplexF64,3,3*nnel)
    for i in 1:nGP
        zu, zt = calc_funsol_dynamic(source_node,gauss_points[i,:], normal[i,:], delta, zconsts)
        
        P1 = J[i]*omegas[i]
        for k in 1:nnel
            P = N[i,k]*P1
            zHELEM[:,3*(k-1)+1:3*k] += zt*P
            zGELEM[:,3*(k-1)+1:3*k] += zu*P
        end
    end
    return zHELEM, zGELEM
end

function integrate_nonsing_dynamic2(source_node,gauss_points,N,normal,J, omegas, delta, zconsts,n)
    nGP = length(omegas)
    nnel = size(N,2)
    zGELEM = zeros(ComplexF64,3,3*nnel)
    for i in 1:nGP
        zu, zt = calc_funsol_dynamic(source_node,gauss_points[i,:], normal[i,:], delta, zconsts)
        
        P1 = J[i]*omegas[i]
        for k in 1:nnel
            P = N[i,k]*P1
            if k != n
                zGELEM[:,3*(k-1)+1:3*k] += zu.*P
            end
        end
        
    end
    return zGELEM
end

function integrate_sing_dynamic(source_node,gauss_points_sing,normal_sing, delta, zconsts,pesos,n)

    nGP = size(pesos,1)
    nnel = size(pesos,2)
    zHELEM = zeros(ComplexF64,3,3*nnel)
    zGELEM = zeros(ComplexF64,3,3*nnel)

    for i in 1:nGP
        zu, zt = calc_funsol_dynamic(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, zconsts)
            
        for k in 1:nnel
            # @infiltrate
            if k!=n
                zHELEM[:,3*(k-1)+1:3*k] += zt.*pesos[i,k]
            end
            if k == n
                zGELEM[:,3*(k-1)+1:3*k] += zu.*pesos[i,k]
            end
        end
    end
    return zHELEM, zGELEM

end


function integrate_sing_dynamic2(source_node,gauss_points_sing,normal_sing, delta, zconsts,C_stat,pesos,n)

    nGP = size(pesos,1)
    nnel = size(pesos,2)
    zHELEM = zeros(ComplexF64,3,3*nnel)

    for i in 1:nGP
        _, zt = calc_funsol_dynamic(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, zconsts)
        _, t = calc_funsol_static(source_node,gauss_points_sing[i,:], normal_sing[i,:], delta, C_stat)
   
        zHELEM[:,3*(n-1)+1:3*n] += (zt - complex(t)).*pesos[i,n]
    end
    return zHELEM

end

function integrate_const_dynamic(source_node, gauss_points, normal, J, omegas, delta, zconsts)
    nGP = length(omegas)
    zHELEM = zeros(ComplexF64,3,3)
    zGELEM = zeros(ComplexF64,3,3)
    for i in 1:nGP
        zu, zt = calc_funsol_dynamic(source_node,gauss_points[i,:], normal[i,:], delta, zconsts)
        
        P = J[i]*omegas[i]

        zHELEM += zt.*P
        zGELEM += zu.*P
        
    end
    return zHELEM, zGELEM
end

function integrate_const_sing_dynamic(source_node,gauss_points,normal,J,omegas, delta, zconsts,C_stat)

    nGP = size(J,1)
    zHELEM = zeros(ComplexF64,3,3)

    for i in 1:nGP
        _, zt = calc_funsol_dynamic(source_node,gauss_points[i,:], normal[i,:], delta, zconsts)
        _, t = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
   
        zHELEM += (zt - complex(t)).*J[i]*omegas[i]
    end
    return zHELEM

end