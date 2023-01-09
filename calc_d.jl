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


        if norm(p2 - p1) < tol
            p1 = p2
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
    return p1,dist
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
    d = norm(p - source)

    return d
end

function teste_optim(points,source_point)

    csis_cont = -1.0:2.0:1.0
    lower = [-1.0,-1.0]
    upper = [1.0,1.0]

    res = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), (G,x) -> ∇d!(G,x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])
    csi,eta = Optim.minimizer(res)
    N = calc_N_matrix(csis_cont,[csi eta])
    p1 = vec(N*points)

    dist = norm(p1-source_point)

    return p1, dist

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
        plot_d(points,source_point,p)
        @infiltrate


        if norm(d2-d) < tol
            d = d2
            break
        end
        
        d = d2
        iter = iter+1
    end


    return p, d
end

function plot_d(points,source_point,p1)
    Plots.scatter(points[:,1], points[:,2],label = "points")
    Plots.plot!([points[:,1]; points[1,1]], [points[:,2]; points[1,2]],label = "element")
    Plots.scatter!([source_point[1]], [source_point[2]], label = "source")
    plt = Plots.scatter!([p1[1]], [p1[2]], label = "closest found", marker = :cross)
    plt = Plots.plot!([p1[1],source_point[1]], [p1[2],source_point[2]], label = string("d = ",norm(p1-source_point)), marker = :cross)
    plt
end

function solve_plot(f,points,source_point)

    p1,d = f(points,source_point)

    Plots.scatter(points[:,1], points[:,2],label = "points")
    Plots.plot!([points[:,1]; points[1,1]], [points[:,2]; points[1,2]],label = "element")
    Plots.scatter!([source_point[1]], [source_point[2]], label = "source")
    plt = Plots.scatter!([p1[1]], [p1[2]], label = "closest found", marker = :cross)
    plt = Plots.plot!([p1[1],source_point[1]], [p1[2],source_point[2]], label = string("d = ",d), marker = :cross)

    return plt

end

## Tests

begin #1
    points = [
        -0.383226     -0.709281  0.0
    -0.151206     -0.709665  0.0
    -2.75258e-12  -1.0       0.0
    -0.333333     -1.0       0.0
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

solve_plot(teste_in,points,source_point)
# solve_plot(teste_out,points,source_point)
solve_plot(teste_optim,points,source_point)

solve_plot(teste_new,points,source_point)


@btime teste_in(points,source_point)
@btime teste_optim(points,source_point)
@btime teste_new(points,source_point)






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