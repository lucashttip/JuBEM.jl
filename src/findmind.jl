using Optim

function ∇d!(G,csi, eta; points, source, csis_cont)
    
    N = calc_N_gen(csis_cont, [csi eta])
    Ncsi = calc_N_gen(csis_cont, [csi eta];dg=:dNdc)
    Neta = calc_N_gen(csis_cont, [csi eta];dg=:dNde)
    p = vec(N*points)    
    d1 = sqrt(sum((p - source).^2))

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)

    G[1] = 1/(2*d1) * (2*pcsi'*p - 2*pcsi'*source)
    G[2] = 1/(2*d1) * (2*peta'*p - 2*peta'*source)
    return G
end

function d(csi, eta; points, source, csis_cont)
    
    N = calc_N_gen(csis_cont, [csi eta])
    p = vec(N*points)    
    d = sqrt(sum((p - source).^2))

    return d
end

function fg!(F,G,x; points, source, csis_cont)
    # do common computations here
    # ...
    csi = x[1]
    eta = x[2]

    N = JuBEM.calc_N_gen(csis_cont, [csi eta])
    Ncsi = calc_N_gen(csis_cont, [csi eta];dg=:dNdc)
    Neta = calc_N_gen(csis_cont, [csi eta];dg=:dNde)

    p = vec(N*points)    
    d = norm(p - source)

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)

    if G !== nothing
      # code to compute gradient here
      # writing the result to the vector G
        G[1] = 1/(2*d) * (2*pcsi'*p - 2*pcsi'*source)
        G[2] = 1/(2*d) * (2*peta'*p - 2*peta'*source)
    end
    if F !== nothing
      # value = ... code to compute objective function
      return d
    end
end

function d∇d∇2d(csi, eta, points, source, csis_cont)

    N = calc_N_gen(csis_cont, [csi eta])
    Ncsi = calc_N_gen(csis_cont, [csi eta];dg=:dNdc)
    Neta = calc_N_gen(csis_cont, [csi eta];dg=:dNde)
    Ncsi2 = calc_N_gen(csis_cont, [csi eta];dg = :d2Ndc2)
    Neta2 = calc_N_gen(csis_cont, [eta csi];dg = :d2Ndc2)
    Ncsieta = calc_N_gen(csis_cont, [csi eta];dg = :d2Ndce)
    p = vec(N*points)    
    d = norm(p - source)

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)
    pcsi2 = vec(Ncsi2*points)
    peta2 = vec(Neta2*points)
    pcsieta = vec(Ncsieta*points)

    f = 1/(2*d)
    fl = -1/(4*d^3)
    gcsi = (2*pcsi'*p - 2*pcsi'*source)
    geta = (2*peta'*p - 2*peta'*source)
    gcsi2 = 2*(dot(pcsi,pcsi) + dot(p,pcsi2)) - 2*dot(source,pcsi2)
    geta2 = 2*(dot(peta,peta) + dot(p,peta2)) - 2*dot(source,peta2)
    gcsieta = 2*(dot(pcsi,peta) + dot(p,pcsieta)) - 2*dot(source,pcsieta)

    H = zeros(2,2)
    G = zeros(2)

    # code to compute gradient here
    # writing the result to the vector G
    H[1,1] = fl*gcsi*gcsi + gcsi2*f
    H[1,2] = fl*gcsi*geta + gcsieta*f
    H[2,1] = H[1,2]
    H[2,2] = fl*geta*geta + geta2*f

    # code to compute gradient here
    # writing the result to the vector G
    G[1] = f * gcsi
    G[2] = f * geta

    return d, G, H

end

function findmind_optim(points, source_point, csis_cont)

    lower = [-1.0,-1.0]
    upper = [1.0,1.0]
    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0], Fminbox(BFGS()))

    return Optim.minimizer(res), minimum(res)
end