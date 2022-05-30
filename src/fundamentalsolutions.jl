function calc_funsol_static(source_node,field_node, normal, delta, C_stat)
    u = zeros(3,3)
    t = zeros(3,3)

    disp = field_node-source_node
    R = norm(disp)
    R2 = R^2
    # ! Derivatives with respect to coordinates
    Rd = disp./R
    # ! Derivative with respect to normal
    Rdn = dot(Rd,normal)
    C1 = C_stat[1]
    C2 = C_stat[2]
    C3 = C_stat[3]
    C4 = C_stat[4]

    for j in 1:3
        for i in 1:3 
            u[i,j] = (C1/R)*(C2*delta[i,j]+Rd[i]*Rd[j])
            t[i,j] = (C3/R2) * (Rdn * (C4*delta[i,j] + 3.0*Rd[i]*Rd[j]) + C4*(Rd[j]*normal[i] - Rd[i]*normal[j]))
        end
    end

    return u, t

end


function calc_funsol_dynamic(source_node,field_node, normal, delta, zconsts)
    
    zu = zeros(ComplexF64,3,3)
    zt = zeros(ComplexF64,3,3)    
    
    disp = field_node - source_node

    R = norm(disp)
    R2 = R^2
    Rd = disp./R
    Rdn = dot(Rd,normal)
    zWi = zconsts.zWi
    zC0 = zconsts.zC0
    zC1 = zconsts.zC1
    zC2 = zconsts.zC2
    zKP = zconsts.zKP
    zKS = zconsts.zKS

    ZZP=zKP*R
    ZZS=zKS*R
    ZEZP=exp(ZZP)
    ZEZS=exp(ZZS)
    ZP2=ZZP*ZZP
    ZS2=ZZS*ZZS

    ZFHI=(1.0+1.0/ZS2-1.0/ZZS)*ZEZS/R - zC2*(1.0/ZP2-1.0/ZZP)*ZEZP/R
    ZCAPPA=(1.0+3.0/ZS2-3.0/ZZS)*ZEZS/R - zC2*(1.0+3.0/ZP2-3.0/ZZP)*ZEZP/R
    ZFHIDR=(-2.0+ZZS+3.0/ZZS-3.0/ZS2)*ZEZS/R^2 - zC2*(-1.0+3.0/ZZP-3.0/ZP2)*ZEZP/R^2
    ZCAPPADR=(ZZS-4.0+9.0/ZZS-9.0/ZS2)*ZEZS/R^2 - zC2*(ZZP-4.0+9.0/ZZP-9.0/ZP2)*ZEZP/R^2

    ZAA=ZFHIDR-ZCAPPA/R
    ZBB=4.0*ZCAPPA/R-2.0*ZCAPPADR
    ZCC=(zC1-2.0)*(ZAA+5.0-1*ZBB-3.0*ZCAPPA/R)-2.0*ZCAPPA/R      
    
    for j in 1:3
        for i in 1:3
            zu[i,j]=zC0*(ZFHI*delta[i,j]-ZCAPPA*Rd[j]*Rd[i])
            zt[i,j]=(1.0/(4.0*pi))*((ZAA*(Rdn*delta[i,j]+Rd[j]*normal[i])) + Rd[i]*Rd[j]*Rdn*ZBB+Rd[i]*normal[j]*ZCC)
        end
    end

    return zu, zt
end




