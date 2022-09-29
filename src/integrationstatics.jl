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

function integrate_const_static(source_node, gauss_points, normal, J, omegas, delta, C_stat)
    nGP = length(omegas)
    HELEM = zeros(3,3)
    GELEM = zeros(3,3)
    for i in 1:nGP
        u, t = calc_funsol_static(source_node,gauss_points[i,:], normal[i,:], delta, C_stat)
            
        P = J[i]*omegas[i]

        HELEM += t.*P
        GELEM += u.*P
        
    end
    return HELEM, GELEM
end

function integration_const_raw(source_node, points, normal, csis, omegas, delta, C_stat)

    ngp = length(csis)

    HELEM = zeros(3,3)
    GELEM = zeros(3,3)

    P = zeros(2,4)
    F = zeros(4)
    XJ = zeros(2,3)
    gp = zeros(3)

        for JG in 1:ngp

            G2=csis[JG]
            P2=omegas[JG]
            SP=1.0+G2
            SM=1.0-G2
            P[1,1] = -0.25*SM
            P[1,2] =  0.25*SM
            P[1,3] =  0.25*SP
            P[1,4] = -0.25*SP

            for ig in 1:ngp

                G1=csis[ig]
                P1=omegas[ig]
                RP=1.0+G1
                RM=1.0-G1
                F[1]=0.25*RM*SM
                F[2]=0.25*RP*SM
                F[3]=0.25*RP*SP
                F[4]=0.25*RM*SP
                P[2,1]=-0.25*RM
                P[2,2]=-0.25*RP
                P[2,3]= 0.25*RP
                P[2,4]= 0.25*RM


                # ! *
                # ! * CALCULA A RELAÇÃO ENTRE AS COORDENADAS CARTESIANAS E HOMOGสNEAS
                # ! *
                for i in 1:2
                    for j in 1:3
                    TEMP=0.0
                        for k in 1:4
                            TEMP=TEMP+P[i,k]*points[k,j]
                        end
                    XJ[i,j]=TEMP
                    end
                end
                # ! *
                # ! * CALCULA O JACOBIANO
                # ! *
                DET=sqrt((XJ[1,2]*XJ[2,3]-XJ[2,2]*XJ[1,3])^2 + (XJ[2,1]*XJ[1,3]-XJ[1,1]*XJ[2,3])^2 + (XJ[1,1]*XJ[2,2]-XJ[2,1]*XJ[1,2])^2  )

                # if DET - 1.0e-7 < 0
                #     error("NONSING : ERRO, JACOBIANO NULO OU NEGATIVO = ")
                # end

                # ! *
                # ! * CALCULA AS COORDENADAS DO PONTO DE INTEGRAÇÃO
                # ! *
                gp .= 0.0
                for i in 1:4
                    gp[1]=gp[1]+points[i,1]*F[i]
                    gp[2]=gp[2]+points[i,2]*F[i]
                    gp[3]=gp[3]+points[i,3]*F[i]
                end
                # ! *
                # ! * ACIONA ROTINA QUE CALCULA A SOLUวรO FUNDAMENTAL DINÂMICA 3D
                # ! *
                u, t = calc_funsol_static(source_node,gp, normal, delta, C_stat)

                P12=P1*P2*DET
                HELEM = HELEM + t*P12
                GELEM = GELEM + u*P12


            end
        end
        return HELEM, GELEM

end