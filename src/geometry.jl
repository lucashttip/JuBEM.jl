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

"""
    r, c, d = calc_dist(source, points, rules)
    
    Computes the minimum distance d from the source point to the elemnent defined by points. r is the index of the rule to be used and c is the local coordinates of the closes point 
"""
function calc_dist(source, points, rules)

    dist = Inf

    for i in axes(points,1)
        d = norm((source - points[i,:]))
        if d < dist
            dist = d
        end
    end

    aux_points = rules.N_aux*points

    l = [norm(aux_points[2,:] - aux_points[4,:]),norm(aux_points[3,:] - aux_points[1,:])]

    d = dist./l
    c = []
    if any(d .< 2.0)
        c,dist = findmind_optim(points,source,rules.csis_points)              # Apenas optim
        d = dist./l
    end
    
    r = 0
    for i in eachindex(rules.dists)
        if all(d .> rules.dists[i])
            r = i
            break
        end
    end

    return r, c, d
end

function calc_dist2(source, points,dists,csis_cont)

    dist = Inf

    for i in eachindex(points[:,1])
        d = norm((source - points[i,:]))
        if d < dist
            dist = d
        end
    end

    l = [norm(points[2,:] - points[1,:]),norm(points[3,:] - points[2,:])]

    d = dist./l
    c = []
    if any(d .< 1.0)
<<<<<<< HEAD
        # c, dist = findmind_projection(points,source,csis_cont)   # Apenas algo projeção 
        c,dist = findmind_optim(points,source,csis_cont)              # Apenas optim
        # c,dist = findmind_combined(points,source,csis_cont)    # Combinado
        # c, dist = bialecki(points,source_point,csis_cont)       # Algoritmo baseado no paper do Bialecki
=======
        # dist = Inf

        for i in 1:4
            if i < 4
                p1 = i
                p2 = i+1
            else
                p1 = 4
                p2 = 1
            end
            xv = points[p2,:] - points[p1,:]
            pv = source - points[p1,:]
            alfa = dot(xv,pv)
            if alfa>0.0 && alfa <=1.0
                xp = points[p1,:] + alfa*xv
            elseif alfa<=0.0
                xp = points[p1,:]
            else
                xp = points[p2,:]
            end

            d2 = norm(xp-source)
            if d2<=dist
                dist=d2
                if i == 1
                    xi = 2*norm((xp - points[p1,:]))/norm((points[p2,:]-points[p1,:])) - 1.0
                    c = [xi,-1.]
                elseif i == 2
                    eta = 2*norm((xp - points[p1,:]))/norm((points[p2,:]-points[p1,:])) - 1.0
                    c = [1.,eta]
                elseif i == 3
                    xi = 2*norm((xp - points[p2,:]))/norm((points[p1,:]-points[p2,:])) - 1.0
                    c = [xi,1.]
                elseif i == 4
                    eta = 2*norm((xp - points[p2,:]))/norm((points[p1,:]-points[p2,:])) - 1.0
                    c = [-1.,eta]
                end
            end

        end

>>>>>>> clean-main
        d = dist./l
    end
    
    r = 0
    for i in eachindex(dists)
        if all(d .> dists[i])
            r = i
            break
        end
    end

    return r, c, d
end