function calc_FJ(x,points,source,csis_cont)

    d, G, H = d∇d∇2d(x[1], x[2], points, source, csis_cont)

    F = zeros(2)
    J = zeros(2,2)

    F[1] = 2*d*G[1]
    F[2] = 2*d*G[2]

    J[1,1] = 2*G[1]*G[1]+2*d*H[1,1] 
    J[1,2] = 2*(G[1]*G[2] + d*H[1,2])
    J[2,1] = J[1,2] 
    J[2,2] = 2*G[2]*G[2]+2*d*H[2,2]

    return F,J
end

function bialecki(points,source_point,csis_cont)
    # Only works for bilinear elements as of now

    x, converged = bialecki_general(points,source_point,csis_cont)

    if !converged
        # Fazer procura nas linhas
        x = [0,0]
        d = Inf
        converged_online = false
        for i in 1:2
            for j in 1:2
                x2,converged = bialecki_online(points,source_point,csis_cont,[i,csis_cont[j]])
                if converged
                    N = calc_N_matrix(csis_cont,[x2[1] x2[2]])
                    p1 = vec(N*points)
                    dist = norm(p1-source_point)
                    if dist < d
                        d = dist
                        x = x2
                    end
                    converged_online = true
                end
            end
        end

        if !converged_online
            x = [0,0]
            d = Inf
            for i in 1:2
                for j in 1:2
                    x2 = [csis_cont[i],csis_cont[j]]
                    N = calc_N_matrix(csis_cont,[x2[1] x2[2]])
                    p1 = vec(N*points)
                    dist = norm(p1-source_point)
                    if dist < d
                        d = dist
                        x = x2
                    end
                end
            end
        end

    end

    N = calc_N_matrix(csis_cont,[x[1] x[2]])
    p1 = vec(N*points)

    dist = norm(p1-source_point)


    return x, dist

end

function bialecki_general(points,source_point,csis_cont)

    x = [0,0]

    l1 = norm(points[2,:] - points[1,:])
    l2 = norm(points[4,:] - points[1,:])
    l = minimum([l1,l2])

    itermax = 8
    iter = 1
    tol = 1e-2*l
    t = Inf
    converged = true

    while t>tol && iter <= itermax

        F, J = calc_FJ(x,points,source_point,csis_cont)

        dx = inv(J)*F
        t = norm(dx)

        x = x - dx

        iter = iter+1
    end

    if t > tol || iter > itermax
        converged = false
    end
    if any(abs.(x) .> 1.0)
        converged = false
    end

    return x, converged

end

function bialecki_online(points,source_point,csis_cont,constraint)

    x = [0.0,0.0]
    x[constraint[1]] = constraint[2]

    if constraint[1] == 1
        idxfree = 2
    else
        idxfree = 1
    end

    l1 = norm(points[2,:] - points[1,:])
    l2 = norm(points[4,:] - points[1,:])
    l = minimum([l1,l2])

    itermax = 10
    iter = 1
    tol = 1e-2*l
    t = Inf
    converged = true

    while t>tol && iter <= itermax

        F, J = calc_FJ(x,points,source_point,csis_cont)

        dx = F[idxfree]/J[idxfree,idxfree]
        t = abs(dx)

        x[idxfree] = x[idxfree] - dx

        iter = iter+1
    end

    if t > tol || iter > itermax
        converged = false
    end
    if any(abs.(x) .> 1.0)
        converged = false
    end

    return x, converged


end