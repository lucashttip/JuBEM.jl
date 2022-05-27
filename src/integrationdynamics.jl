
function calc_nonsing(source_node,gauss_points,field_normal,N,normal,J, omega, zGe, zSwv, zPwv, delta, fr, zconsts)
    nGP = length(omega)
    nnel = size(N,3)
    ZHELEM = zeros(ComplexF64,3,3*nnel)
    ZGELEM = zeros(ComplexF64,3,3*nnel)
    for i in 1:nGP
        for j in 1:nGP
    
            zu, zt = calc_funsol_dynamic(source_node,gauss_points[i,j,:], normal[i,j,:], zGe, zSwv, zPwv, delta, fr, zconsts)
    
            for k in 1:nnel
                P = N[i,j,k]*J[i,j]*omega[i]*omega[j]
                ZHELEM[:,3*(k-1)+1:3*k] = zt*P
                ZGELEM[:,3*(k-1)+1:3*k] = zu*P
            end
        end
    end
    return ZHELEM, ZGELEM
end



    