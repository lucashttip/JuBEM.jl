using Revise
using JuBEM
using Plots
using Test
using LinearAlgebra
using Optim
using Infiltrator

using BenchmarkTools


## Functions

function teste_in(points,source_point)
    csis_cont = -1.0:2.0:1.0
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

            if (norm(l[1]) == 1.0 && norm(l[2]) < 1.0)

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
                end
                l[2] = c

            elseif (norm(l[2]) == 1.0 && norm(l[1]) < 1.0)

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

    # for i in 1:2
    #     if l[i] > 1
    #         l[i] = 1
    #     end
    #     if l[i]<-1
    #         l[i] = -1
    #     end
    # end

    N = calc_N_matrix(csis_cont,[l[1] l[2]])

    p1 = vec(N*points)
    d = source_point - p1
    dist = norm(d)

    # return l, p1, d
    return p1,dist, l
end

function d_pointline(a,b,s)

    ab = b - a

    alfa = (dot(ab,s) - dot(ab,a))/dot(ab,ab)

    c = a + alfa*ab


    return c, 2*alfa-1
end

function teste_out(points,source_point)
    csis_cont = -1.0:2.0:1.0
    tol = 1e-3
    itermax = 30

    # Primeiro ponto: ponto central do elemento.
    l = [0, 0, 0]

    N = calc_N_matrix(csis_cont,[l[1] l[2]])
    p1 = vec(N*points)
    p0 = p1

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

    # return l, p1, d
    return p1,dist
end

function ∇d!(G,csi, eta; points, source, csis_cont)
    
    N = JuBEM.calc_N_matrix(csis_cont, [csi eta])
    Ncsi = calc_dNdcsi_matrix(csis_cont, [csi eta])
    Neta = calc_dNdeta_matrix(csis_cont, [csi eta])
    p = vec(N*points)    
    d = norm(p - source)

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)

    G[1] = 1/(2*d) * (2*pcsi'*p - 2*pcsi'*source)
    G[2] = 1/(2*d) * (2*peta'*p - 2*peta'*source)
    return G
end

function ∇2d!(H,csi, eta; points, source, csis_cont)
    
    N = JuBEM.calc_N_matrix(csis_cont, [csi eta])
    Ncsi = calc_dNdcsi_matrix(csis_cont, [csi eta])
    Neta = calc_dNdeta_matrix(csis_cont, [csi eta])
    Ncsi2 = JuBEM.calc_N_gen(csis_cont, [csi eta];dg = :d2Ndc2)
    Neta2 = JuBEM.calc_N_gen(csis_cont, [eta csi];dg = :d2Ndc2)
    Ncsieta = JuBEM.calc_N_gen(csis_cont, [csi eta];dg = :d2Ndce)
    p = vec(N*points)    
    d1 = sqrt(sum((p - source).^2))

    pcsi = vec(Ncsi*points)
    peta = vec(Neta*points)
    pcsi2 = vec(Ncsi2*points)
    peta2 = vec(Neta2*points)
    pcsieta = vec(Ncsieta*points)

    f = 1/(2*d1)
    fl = -1/(4*d1^3)
    gcsi = (2*pcsi'*p - 2*pcsi'*source)
    geta = (2*peta'*p - 2*peta'*source)
    gcsi2 = 2*(dot(pcsi,pcsi) + dot(p,pcsi2)) - 2*dot(source,pcsi2)
    geta2 = 2*(dot(peta,peta) + dot(p,peta2)) - 2*dot(source,peta2)
    gcsieta = 2*(dot(pcsi,peta) + dot(p,pcsieta)) - 2*dot(source,pcsieta)

    H[1,1] = fl*gcsi*gcsi + gcsi2*f
    H[1,2] = fl*gcsi*geta + gcsieta*f
    H[2,1] = H[1,2]
    H[2,2] = fl*geta*geta + geta2*f
    return H
end

function d(csi, eta; points, source, csis_cont)
    
    N = calc_N_matrix(csis_cont, [csi eta])
    p = vec(N*points)    
    d = norm(p - source)

    return d
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

function fgh!(F,G,H,x; points, source, csis_cont)
    # do common computations here
    # ...
    csi = x[1]
    eta = x[2]

    N = JuBEM.calc_N_matrix(csis_cont, [csi eta])
    Ncsi = calc_dNdcsi_matrix(csis_cont, [csi eta])
    Neta = calc_dNdeta_matrix(csis_cont, [csi eta])
    Ncsi2 = JuBEM.calc_N_gen(csis_cont, [csi eta];dg = :d2Ndc2)
    Neta2 = JuBEM.calc_N_gen(csis_cont, [eta csi];dg = :d2Ndc2)
    Ncsieta = JuBEM.calc_N_gen(csis_cont, [csi eta];dg = :d2Ndce)
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

    if H !== nothing
        # code to compute gradient here
        # writing the result to the vector G
        H[1,1] = fl*gcsi*gcsi + gcsi2*f
        H[1,2] = fl*gcsi*geta + gcsieta*f
        H[2,1] = H[1,2]
        H[2,2] = fl*geta*geta + geta2*f
    end
    if G !== nothing
      # code to compute gradient here
      # writing the result to the vector G
        G[1] = f * gcsi
        G[2] = f * geta
    end
    if F !== nothing
      # value = ... code to compute objective function
      return d
    end
end

function teste_optim_simple(points,source_point)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])
    csi,eta = Optim.minimizer(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist

end

function teste_optim_grad(points,source_point)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0], Fminbox(BFGS()))
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist, csi, eta

end

function teste_optim_grad2(points,source_point;initial_csis)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    N0 = calc_N_matrix(csis_cont,initial_csis)
    p0 = N0*points

    k = 0
    d = Inf
    for j in axes(p0,1)
        dist = norm(p0[j,:]-source_point)
        if dist < d
            k = j
            d = dist
        end
    end

    x0 = initial_csis[k,:]

    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,x0, Fminbox(BFGS()))
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist

end

function teste_optim2(points,source_point;ino = LBFGS())

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0], Fminbox(ino))
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist, c

end

function teste_optimgh(points,source_point)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(Optim.only_fgh!((F,G,H,x) -> fgh!(F,G,H,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0])
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist, c

end

function teste_optimipnewton(points,source_point,x0)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(Optim.only_fgh!((F,G,H,x) -> fgh!(F,G,H,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,x0,IPNewton())
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist, c

end

function teste_new(points,source_point)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    d = 100
    itermax = 300
    iter = 1
    tol = 1e-4
    x = [0.0 0.0]
    G = [0.0 0.0]
    s = 0.1

    N = calc_N_matrix(csis_cont, x)
    p = vec(N*points) 
    d = norm(p-source_point)


    while iter <= itermax

        ∇d!(G,x[1],x[2]; points=points, source=source_point, csis_cont=csis_cont)

        dx = G*s

        x = x - dx

        for i in 1:2
            if norm(x[i]) > 1.0
                x[i] = sign(x[i])*1.0
            end
        end

        N = calc_N_matrix(csis_cont, x)
        p = vec(N*points) 
        d2 = norm(p-source_point)


        if norm(d2-d) < tol
            d = d2
            break
        end
        
        d = d2
        iter = iter+1
    end


    return p, d
end

function teste_optimgrad_x0(points,source_point,x0)
    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

   

    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,x0, Fminbox(BFGS()))
    csi,eta = Optim.minimizer(res)
    c = Optim.converged(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist
end

function teste_combined(points,source)

    p1,dist, l = teste_in(points,source)

    c =l[1:2]

    for i in 1:2
        if c[i] == 1
            c[i] += -1e-3
        end
        if c[i] == -1
            c[i] += 1e-3
        end
    end

    p1,dist = teste_optimgrad_x0(points,source,c)
    # p1,dist = teste_optimipnewton(points,source,c)


    return p1, dist, l

end

function plot_d(points,source_point,p1)
    Plots.scatter(points[:,1], points[:,2],label = "points")
    Plots.plot!([points[:,1]; points[1,1]], [points[:,2]; points[1,2]],label = "element")
    Plots.scatter!([source_point[1]], [source_point[2]], label = "source")
    plt = Plots.scatter!([p1[1]], [p1[2]], label = "closest found", marker = :cross)
    plt = Plots.plot!([p1[1],source_point[1]], [p1[2],source_point[2]], label = string("d = ",norm(p1-source_point)), marker = :cross, aspect_ratio=:equal)
    plt
end

function solve_plot(f,points,source_point)

    p1,d = f(points,source_point)

    Plots.scatter(points[:,1], points[:,2],label = "points")
    Plots.plot!([points[:,1]; points[1,1]], [points[:,2]; points[1,2]],label = "element")
    Plots.scatter!([source_point[1]], [source_point[2]], label = "source")
    plt = Plots.scatter!([p1[1]], [p1[2]], label = "closest found", marker = :cross)
    plt = Plots.plot!([p1[1],source_point[1]], [p1[2],source_point[2]], label = string("d = ",d), marker = :cross, aspect_ratio=:equal)

    return plt

end

## Tests

begin #1
    points = [
        -0.383226     -0.6  0.0
    -0.151206     -0.709665  0.0
    -2.75258e-12  -1.0       0.0
    -0.25     -1.5       0.0
    ]

    source_point = [
        0.21677192715792648
        -0.853764880489005
        0.0
    ]
end

begin #1.2
    points = [
        -0.383226     -0.709665  0.0
    -0.151206     -0.709665  0.0
    -2.75258e-12  -1.0       0.0
    -0.4     -1.0       0.0
    ]

    source_point = [
        0.21677192715792648
        -0.853764880489005
        0.0
    ]
end

begin #2
    points = [
        -0.00782734  -0.416598   0.0
        -0.211126    -0.0495134  0.0
        0.0486605    0.0586859  0.0
        0.148951    -0.13939    0.0
    ]

    source_point = [
        0.5323000645344033
        -0.2504690686492321
        0.0
    ]
end

begin #3
    points = [
        -2  -2  0.0
        2   -2  0.0
        2   2   0.0
        -2  2   0.0
    ]

    source_point = [
        0.5323000645344033
        -0.2504690686492321
        0.3
    ]
end

begin #4
    points = [
        -0.383226     -0.709281  0.0
    -0.151206     -0.709665  0.0
    -2.75258e-12  -1.0       0.0
    -0.333333     -1.0       0.0
    ]

    source_point = [
        -0.2
        -0.9
        0.3
    ]
end

begin
    initial_csis = [
        -0.95 -0.95
        0.95 -0.95
        0.95 0.95
        -0.95 0.95
        0 -0.95
        0.95 0
        0 0.95
        -0.95 0
        0.0 0.0
    ]
end

solve_plot(teste_in,points,source_point)
# solve_plot(teste_out,points,source_point)
solve_plot(teste_optim_simple,points,source_point)
solve_plot(teste_optim_grad,points,source_point)
solve_plot(teste_optimgh,points,source_point)
solve_plot(teste_optim2,points,source_point)
solve_plot(teste_new,points,source_point)
solve_plot(teste_optimipnewton,points,source_point)
solve_plot((p,s) -> teste_optim_grad2(p,s;initial_csis = initial_csis),points,source_point)
solve_plot(teste_combined,points,source_point)




@btime teste_in(points,source_point)
@btime teste_optim_grad(points,source_point)
@btime teste_combined(points,source_point)
@btime teste_optim_simple(points,source_point)
@btime teste_new(points,source_point)
@btime teste_optim2(points,source_point;ino = BFGS())
@btime teste_optimgh(points,source_point)
@btime teste_optimipnewton(points,source_point,[0.0,0.0])
@btime teste_optim_grad2(points,source_point;initial_csis=initial_csis)

inos = [NelderMead(), SimulatedAnnealing(), BFGS(), LBFGS(), ConjugateGradient(), GradientDescent(), MomentumGradientDescent(), AcceleratedGradientDescent()]






res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0])





## Testsets

points1 = [
        -0.383226     -0.709281  0.0
        -0.151206     -0.709665  0.0
        -2.75258e-12  -1.0       0.0
        -0.333333     -1.0       0.0
]

source_point1 = [
    0.21677192715792648
    -0.853764880489005
    0.0
]

d1 = 0.265

points2 = [
    -0.00782734  -0.416598   0.0
    -0.211126    -0.0495134  0.0
    0.0486605    0.0586859  0.0
    0.148951    -0.13939    0.0
]

source_point2 = [
    0.5323000645344033
    -0.2504690686492321
    0.0
]
    
d2 = 0.4

points3 = [
    -2  -2  0.0
    2   -2  0.0
    2   2   0.0
    -2  2   0.0
]

source_point3 = [
    0.5323000645344033
    -0.2504690686492321
    0.3
]
d3 = 0.3


@testset "teste_in" begin
    _,d1c = teste_in(points1,source_point1)
    _,d2c = teste_in(points2,source_point2)
    _,d3c = teste_in(points3,source_point3)


    @test d1c <= d1
    @test d2c <= d2
    @test d3c <= d3
end


@testset "teste_optim" begin
    _,d1c = teste_optim(points1,source_point1)
    _,d2c = teste_optim(points2,source_point2)
    _,d3c = teste_optim(points3,source_point3)


    @test d1c <= d1
    @test d2c <= d2
    @test d3c <= d3
end