
CRB = (0.85 + 0.24*log(0.05))/0.05
calc_r(D) = D > 0.05 ? 0.85 + 0.24*log(D) : CRB*D
function pontos_pesos_local_telles(csis, omegas, source,D)

    npg = length(csis)
    npontos = npg^2
    pontos_gauss = zeros(npontos,2)
    pesos = zeros(npontos)

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
        for j in 1:npg
            csi = csis[i]
            eta = csis[j]

            gamma = [csi, eta]
            gp = a.*gamma.^3 + b.*gamma.^2 + c.*gamma .+ d
            
            J = 3 .*a.*gamma.^2 .+ 2 .*b.*gamma .+c

            J = J[1]*J[2]

            idx = npg*(i-1) + j
            pontos_gauss[idx,:] = gp
            pesos[idx] = omegas[i]*omegas[j]*J
        end
    end
    return pontos_gauss, pesos
end