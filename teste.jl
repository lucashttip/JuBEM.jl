
using Infiltrator
using JuBEM
using LinearAlgebra
import Plots

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




function teste(points,source_point)
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
            @infiltrate
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
            @infiltrate
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

    return l, p1, d
end


l, p1, d = teste(points,source_point,csis_cont)


function calc_dNdcsicsi(csisij,ordem,kij,csi,eta)

    nnel = ordem^2

    dNdcsi = zeros(nnel)

    for k in 1:nnel
        dNdcsi[k] = 0
        i = kij[k][1]
        j = kij[k][2]
        


        for n in 1:ordem
            dN1 = 0
            for m in 1:ordem
                dN2 = 1
                for l in 1:ordem
                    if l != i && l != m && l != n
                        dN2 = dN2*(csi - csisij[l])/(csisij[i] - csisij[l])
                    end
                end
                if m != i && m != n
                    dN1 = dN1 + (dN2/(csisij[i] - csisij[m]))
                end
            end

            if n != i
                dNdcsi[k] = dNdcsi[k] + (dN1/(csisij[i] - csisij[n]))
            end
        end

        
        for l in 1:ordem
            if l != j
                dNdcsi[k] = dNdcsi[k]*(eta - csisij[l])/(csisij[j] - csisij[l])
            end
        end
    end
    return dNdcsi
end

function calc_dNdcsieta(csisij,ordem,kij,csi,eta)

    nnel = ordem^2

    dNdcsi = zeros(nnel)

    for k in 1:nnel
        dNdcsi[k] = 0
        i = kij[k][1]
        j = kij[k][2]
        


        for m in 1:ordem
            dN = 1
            for l in 1:ordem
                if l != i && l != m
                    dN = dN*(csi - csisij[l])/(csisij[i] - csisij[l])
                end
            end
            if m != i
                dNdcsi[k] = dNdcsi[k] + (dN/(csisij[i] - csisij[m]))
            end
        end


        dN1 = 0
        for m in 1:ordem
            dN = 1
            for l in 1:ordem
                if l != j && l != m
                    dN = dN*(csi - csisij[l])/(csisij[j] - csisij[l])
                end
            end
            if m != j
                dN1 = dN1 + (dN/(csisij[j] - csisij[m]))
            end
        end
        dNdcsi[k] = dNdcsi[k]*dN1

    end
    
    return dNdcsi
end


using Plots
function plottest(c,IEN)

    scatter(c[:,1],c[:,2],color=:black,label=:none)
    plt = scatter!([c[end,1]],[c[end,2]],color=:red,label=:none)

    for e in axes(IEN,2)
        plt = plot!(c[[IEN[:,e];IEN[1,e]],1],c[[IEN[:,e];IEN[1,e]],2],color=:black,label=:none)
    end

    return plt
end