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

function find_closest_dist(points,source_point,csis_cont)
    tol = 1e-3
    itermax = 30

    # Primeiro ponto: ponto central do elemento.
    l = [0, 0, 0]

    N = calc_N_matrix(csis_cont,[l[1] l[2]])
    p1 = vec(N*points)
    t = Inf

    iterout = 0

    for iter = 1:itermax

        dNcsi = calc_dNdcsi_matrix(csis_cont,[l[1] l[2]])
        dNeta = calc_dNdeta_matrix(csis_cont,[l[1] l[2]])
        # Encontrar vetores do plano tangente ao elemento que passa pelo ponto p1:
        ncsi = vec(dNcsi*points)
        neta = vec(dNeta*points)
        nzeta = cross(ncsi,neta)

        # Vetor entre ponto fonte e p1
        d = source_point - p1

        try
            l = l + [ncsi neta nzeta]\d
        catch
            error("Erro no cálculo do ponto mais próximo.")
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

        p2 = vec(N*points)

        t = norm(p2 - p1)

        if t < tol

            if (abs(l[1]) == 1.0 && abs(l[2]) < 1.0)

                if l[1] == 1
                    a = points[2,:]
                    b = points[3,:]
                else
                    a = points[1,:]
                    b = points[4,:]
                end
                p1, c = d_pointline(a,b,source_point)
                if c > 1.0
                    c = 1.0
                elseif c < 1.0
                    c = -1.0
                end
                l[2] = c

            elseif (abs(l[2]) == 1.0 && abs(l[1]) < 1.0)

                if l[2] == 1
                    a = points[4,:]
                    b = points[3,:]
                else
                    a = points[1,:]
                    b = points[2,:]
                end
                p1, c = d_pointline(a,b,source_point)

                if c > 1.0
                    c = 1.0
                elseif c < 1.0
                    c = -1.0
                end
                l[1] = c

            else
                p1 = p2
            end
            break
        end
        if iter == itermax
            error("Erro encontrando ponto local mais próximo do fonte.")
        end
        p1 = p2
    end

    N = calc_N_matrix(csis_cont,[l[1] l[2]])

    p1 = vec(N*points)
    d = source_point - p1
    dist = norm(d)

    # return l, p1, d
    return l,dist
end

function d_pointline(a,b,s)

    ab = b - a

    alfa = (dot(ab,s) - dot(ab,a))/dot(ab,ab)

    c = a + alfa*ab


    return c, 2*alfa-1
end

function fg!(F,G,x; points, source, csis_cont)
    # do common computations here
    # ...
    csi = x[1]
    eta = x[2]

    N = JuBEM.calc_N_matrix(csis_cont, [csi eta])
    Ncsi = calc_dNdcsi_matrix(csis_cont, [csi eta])
    Neta = calc_dNdeta_matrix(csis_cont, [csi eta])
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

function optimgrad_x0(points,source_point,x0, csis_cont)
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,x0, Fminbox(BFGS()))

    return Optim.minimizer(res), minimum(res)
end

function findmind_combined(points,source, csis_cont)

    l,dist = find_closest_dist(points,source,csis_cont)

    c =l[1:2]

    for i in 1:2
        if c[i] == 1
            c[i] += -1e-3
        end
        if c[i] == -1
            c[i] += 1e-3
        end
    end

    l,dist = optimgrad_x0(points,source,c,csis_cont)

    return l,dist

end