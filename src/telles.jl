
const CRB = (0.85 + 0.24*log(0.05))/0.05
# calc_r(D) = D > 0.05 ? 0.85 + 0.24*log(D) : CRB*D

function calc_r(D)
    if D < 0.0
        error("Valor de D incorreto, deve ser maior que zero e foi: ", D) 
    elseif D >=0.0 && D < 0.05
        return CRB*D
    elseif D>=0.05 && D < 1.3
        return 0.85 + 0.24*log(D)
    elseif D >= 1.3 && D < 3.618
        return 0.893 + 0.0832*log(D)
    else
        return 1
    end
end

function pontos_pesos_local_telles(csis, omegas, source,D)

    npg = size(csis,1)
    pontos_gauss = zeros(npg,2)
    pesos = zeros(npg)

    rb = calc_r.(D)

    q = (1 ./(2 .*(1 .+2 .*rb))) .* ((source.*(3 .-2 .*rb) .- (2 .*source.^3)./(1 .+ 2 .*rb)).*(1 ./(1 .+2 .*rb)) .- source)

    p = (1 ./(3 .*(1 .+ 2 .*rb).^2)).*(4 .*rb.*(1 .-rb) .+ 3 .*(1 .- source.^2))

    aux = sqrt.(q.^2 + p.^3)

    gb = cbrt.(-q .+ aux) +cbrt.(-q .-aux) .+ (source./(1 .+ 2 .*rb)) 

    Q = 1 .+ 3 .*(gb.^2) 

    a = (1 .-rb) ./Q
    b = -3 .*(1 .- rb) .*gb./Q
    c =(rb .+ 3 .*(gb.^2))./Q
    d = -b


    for i in 1:npg
            csi = csis[i,1]
            eta = csis[i,2]

            gamma = [csi, eta]
            gp = a.*gamma.^3 + b.*gamma.^2 + c.*gamma .+ d
            
            J = 3 .*a.*gamma.^2 .+ 2 .*b.*gamma .+c

            J = J[1]*J[2]

            pontos_gauss[i,:] = gp
            pesos[i] = omegas[i]*J
    end
    return pontos_gauss, pesos
end

function telles1(csis, omegas, d)

    npg = size(csis,1)
    pontos_gauss = zeros(npg,2)
    pesos = zeros(npg)

    rb = calc_r(d)


    q = (1/(2*(1+2*rb)))* (((3 -2 *rb)-2/(1+ 2*rb))*(1/(1+2*rb)) - 1)

    q = (1/(2*(1+2*rb)))*(((3 -2 *rb) - 2/(1 + 2*rb))*(1 /(1 +2 *rb)) - 1)

    p = (1 /(3 *(1 + 2 *rb)^2))*(4*rb*(1-rb))

    aux = sqrt(q^2 + p^3)

    gb = cbrt(-q + aux) +cbrt(-q -aux) + (1/(1 + 2*rb)) 

    Q = 1 + 3 *(gb^2) 

    a = (1 -rb) /Q
    b = -3 *(1 - rb) *gb/Q
    c =(rb + 3 *(gb^2))/Q
    d = -b


    for i in 1:npg
            csi = csis[i,1]
            eta = csis[i,2]

            gamma = a*eta^3 + b*eta^2 + c*eta + d
            
            J = 3*a*eta^2 + 2*b*eta + c

            pontos_gauss[i,:] = [csi,gamma]
            pesos[i] = omegas[i]*J
    end
    return pontos_gauss, pesos
end

# TODO: Função interessante para realizar experimentos
function points_weights_local_near_combined(csis, omegas, source,d)

    pontos_telles, pesos_telles = telles1(csis,omegas,d)

    csis_lin = [-1.0,1.0]

    c,IEN = divide_elem(source)

    nt = size(IEN,2)
    np = length(pesos_telles)
    nc = np*nt

    pontos = zeros(nc,2)
    pesos = zeros(nc)

    N = calc_N_gen(csis_lin,pontos_telles)
    dNc = calc_N_gen(csis_lin,pontos_telles;dg=:dNdc)
    dNe = calc_N_gen(csis_lin,pontos_telles;dg=:dNde)


    for t in 1:nt
        pontos[np*(t-1)+1:np*t,:] = N*c[IEN[:,t],:]
        _,J = calc_n_J_matrix(dNc, dNe, [c[IEN[:,t],:] zeros(4)])

        pesos[np*(t-1)+1:np*t] = J.*pesos_telles
    end

    return pontos, pesos
end