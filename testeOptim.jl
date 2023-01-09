using JuBEM
using Optim
import Plots
using BenchmarkTools

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

csis_cont = -1.0:2.0:1.0

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

lower = [-1.0,-1.0]
upper = [1.0,1.0]


@btime optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])

@btime optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), (G,x) -> ∇d!(G,x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])


res1 = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])

res2 = optimize(x -> d(x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), (G,x) -> ∇d!(G,x[1],x[2]; points = points, source = source_point, csis_cont = csis_cont), lower, upper,[0.0,0.0])




csi,eta = Optim.minimizer(res)
N = calc_N_matrix(csis_cont,[csi eta])
p_closest = vec(N*points)

l, p1, dist = teste(points,source_point,csis_cont)

@btime teste(points,source_point,csis_cont)


Plots.scatter(points[:,1], points[:,2],label = "points")
Plots.plot!([points[:,1]; points[1,1]], [points[:,2]; points[1,2]],label = "element")
Plots.scatter!([source_point[1]], [source_point[2]], label = "source")
plt = Plots.scatter!([p_closest[1]], [p_closest[2]], label = "Optim", marker = :star5)
plt = Plots.scatter!([p1[1]], [p1[2]], label = "teste", marker = :cross)

display(plt)