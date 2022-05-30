
function calc_n_J_matrix(dNdcsi, dNdeta, points)

    lx = size(dNdcsi,1)
    ly = size(dNdcsi,2)

    normal = zeros(lx,ly,3)
    J = zeros(lx,ly)

    for i in 1:lx
        for j in 1:ly
            # @infiltrate
            normal[i,j,:], J[i,j] = calc_n_J(dNdcsi[i,j,:],dNdeta[i,j,:],points)
        end
    end

    return normal, J

end


function calc_n_J(dNdcsi, dNdeta,points)
    dpdcsi = dNdcsi'*points
    dpdeta = dNdeta'*points
    # @infiltrate
    v = cross(dpdcsi', dpdeta')
    J = norm(v)
    n = v./J

    return n, J

end


function calc_static_constants(material)

    C_stat = zeros(4)
    
    C_stat[1]=1.0/(16.0*pi*material.Ge*(1.0-material.Nu))
    C_stat[2]=3.0-(4.0*material.Nu)
    C_stat[3]=-1.0/(8.0*pi*(1.0-material.Nu))
    C_stat[4]=1.0-(2.0*material.Nu)

    return C_stat

end