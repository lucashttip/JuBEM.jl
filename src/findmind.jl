using Optim


function ∇d!(G,csi, eta; points, source, csis_cont)
    
    N = JuBEM.calc_N_matrix(csis_cont, [csi eta])
    Ncsi = calc_dNdcsi_matrix(csis_cont, [csi eta])
    Neta = calc_dNdeta_matrix(csis_cont, [csi eta])
    p = vec(N*points)    
    d1 = sqrt(sum((p - source).^2))

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)

    G[1] = 1/(2*d1) * (2*pcsi'*p - 2*pcsi'*source)
    G[2] = 1/(2*d1) * (2*peta'*p - 2*peta'*source)
    return G
end



function d(csi, eta; points, source, csis_cont)
    
    N = calc_N_matrix(csis_cont, [csi eta])
    p = vec(N*points)    
    d = sqrt(sum((p - source).^2))

    return d
end


function findmind(points, source_point, csis_cont)

    lower = [-1.0,-1.0]

    upper = [1.0,1.0]
    # res = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])

    res = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), (G,x) -> ∇d!(G,x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])



    return Optim.minimizer(res), minimum(res)
end