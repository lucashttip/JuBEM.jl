
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

function calc_n_J_matrix_sing(dNcdcsi, dNcdeta, points,node_sing)

    l1 = size(dNcdcsi,1)
    l2 = 4*l1


    normal = zeros(l2,3)
    J = zeros(l2)
    k = 1
    for i in 1:4
        if i==1
            i1 = 1
            i2 = 2
        elseif i ==2
            i1 = 2
            i2 = 3
        elseif i ==3
            i1 = 3
            i2 = 4
        elseif i ==4
            i1 = 4
            i2 = 1
        end
        p1 = points[i1,:]
        p2 = points[i2,:]
        new_points = [p1';p2';node_sing';node_sing']

        for j in 1:l1
                # @infiltrate
                normal[k,:], J[k] = calc_n_J(dNcdcsi[j,:],dNcdeta[j,:],new_points)
                k = k+1
        end
    end
    return normal, J

end


function calc_n_J(dNdcsi, dNdeta,points)
    dpdcsi = dNdcsi'*points
    dpdeta = dNdeta'*points
    v = cross(dpdcsi', dpdeta')
    J = norm(v,2)
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


function calc_dist(source, points,dists,csis_cont)

    dist = Inf

    for i in eachindex(points[:,1])
        d = norm((source - points[i,:]))
        if d < dist
            dist = d
        end
    end

    l = maximum([norm(points[2,:] - points[1,:]),norm(points[3,:] - points[2,:])])

    d = dist/l
    c = []
    if d < 1
        p1, c, dist = find_closest_dist(points,source,csis_cont)
        d = dist/l
    end
    r = 0
    for i in eachindex(dists)
        if d > dists[i]
            r = i
            break
        end
    end

    return r, c
end

function find_closest_dist(points,source_point,csis_cont)

    tol = 1e-3
    itermax = 30

    # Primeiro ponto: ponto central do elemento.
    l = [0, 0, 0]

    N = calc_N_matrix(csis_cont,[l[1] l[2]])
    p1 = vec(N*points)

    for iter = 1:itermax

        dNcsi = calc_dNdcsi_matrix(csis_cont,[l[1] l[2]])
        dNeta = calc_dNdeta_matrix(csis_cont,[l[1] l[2]])
        # Encontrar vetores do plano tangente ao elemento que passa pelo ponto p1:
        ncsi = vec(dNcsi*points)
        neta = vec(dNeta*points)
        nzeta = cross(ncsi,neta)

        # Vetor entre ponto fonte e p1
        d = source_point - p1

        l = l + [ncsi neta nzeta]\d

        N = calc_N_matrix(csis_cont,[l[1] l[2]])

        p2 = vec(N*points)

        if norm(p2 - p1) < tol
            p1 = p2
            break
        end
        if iter == itermax
            error("Erro encontrando ponto local mais próximo do fonte.")
        end
        p1 = p2
    end

    for i in 1:2
    if l[i] > 1
        l[i] = 1
    end
    if l[i]<-1
        l[i] = -1
    end
    end

    N = calc_N_matrix(csis_cont,[l[1] l[2]])

    p1 = vec(N*points)
    d = source_point - p1
    dist = norm(d)

    return p1, l, dist
end