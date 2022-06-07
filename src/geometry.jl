
function calc_n_J_matrix(dNdcsi, dNdeta, points)

    l = size(dNdcsi,1)


    normal = zeros(l,3)
    J = zeros(l)

    for i in 1:l
            # @infiltrate
            normal[i,:], J[i] = calc_n_J(dNdcsi[i,:],dNdeta[i,:],points)
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

function calc_area(points)

    if size(points,2) < 3
        points = [points zeros(size(points,1))]
    end


    p1 = points[1,:]
    p2 = points[2,:]
    p3 = points[3,:]
    p4 = points[4,:]

    l1 = p2-p1
    l2 = p3-p2
    l3 = p4-p3
    l4 = p1-p4

    a1 = norm(cross(l1,l2))/2
    a2 = norm(cross(l3,l4))/2
    
    a = a1+a2

    return a

end

function calc_static_constants(material)

    C_stat = zeros(4)
    
    C_stat[1]=1.0/(16.0*pi*material.Ge*(1.0-material.Nu))
    C_stat[2]=3.0-(4.0*material.Nu)
    C_stat[3]=-1.0/(8.0*pi*(1.0-material.Nu))
    C_stat[4]=1.0-(2.0*material.Nu)

    return C_stat

end