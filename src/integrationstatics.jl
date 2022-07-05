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
            # @infiltrate
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

        j1 = i1% (3)
        j2 = i2% (3)
        j3 = i3% (3)

        j1 == 0 ? j1 = 3 : nothing
        j2 == 0 ? j2 = 3 : nothing
        j3 == 0 ? j3 = 3 : nothing


        idx1 = Int.(collect(j1:3:3*nnodes))
        idx2 = Int.(collect(j2:3:3*nnodes))
        idx3 = Int.(collect(j3:3:3*nnodes))


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

function integrate_non_static2(source_node,points, csi, omega, delta,normal, C_stat,nnel,ex,ey)
    nGP = length(omega)
    HELEM = zeros(3,3*nnel)
    GELEM = zeros(3,3*nnel)

    P = zeros(2,4)
    F =zeros(4)
    XJ = zeros(2,3)
    N = zeros(4)
    gp = zeros(3)

    for JG in 1:nGP

        G2=csi[JG]
        P2=omega[JG]
        SP=1.0+G2
        SM=1.0-G2
        P[1,1]=-0.25*SM
        P[1,2]= 0.25*SM
        P[1,3]= 0.25*SP
        P[1,4]=-0.25*SP

        

        for ig in 1:nGP

            G1=csi[ig]
            P1=omega[ig]
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


            N1xi = (1-ex-G1)/(2*(1-ex))
            N2xi = (1-ex+G1)/(2*(1-ex))
            N1eta = (1-ey-G2)/(2*(1-ey))
            N2eta = (1-ey+G2)/(2*(1-ey))

            N[1] = N1xi*N1eta
            N[2] = N2xi*N1eta
            N[3] = N2xi*N2eta
            N[4] = N1xi*N2eta


            for i in 1:2
                for j in 1:3
                TEMP=0.0
                    for k in 1:4
                        TEMP=TEMP+P[i,k]*points[j,k]
                    end 
                XJ[i,j]=TEMP
                end 
            end 
            # ! *
            # ! * CALCULA O JACOBIANO
            # ! *
            DET=sqrt((XJ[1,2]*XJ[2,3]-XJ[2,2]*XJ[1,3])^2 + (XJ[2,1]*XJ[1,3]-XJ[1,1]*XJ[2,3])^2 + (XJ[1,1]*XJ[2,2]-XJ[2,1]*XJ[1,2])^2 )

            if (DET < 10^(-7)) 
                    error("Jacobiano muito pequeno")
            end

            gp .= 0.0
            for  i in 1:4
                gp[1]=gp[1]+points[1,i]*F[i]
                gp[2]=gp[2]+points[2,i]*F[i]
                gp[3]=gp[3]+points[3,i]*F[i]
            end

            u, t = calc_funsol_static(source_node,gp, normal, delta, C_stat)

            P12=P1*P2*DET

            for n in 1:4
                # HELEM = HELEM + t*P12
                # GELEM = GELEM + u*P12
                HELEM[:,3*(n-1)+1:3*n] += t.*P12*N[n]
                GELEM[:,3*(n-1)+1:3*n] += u.*P12*N[n]
            end
        end
    end

    return HELEM, GELEM    
end