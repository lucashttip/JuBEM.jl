

function generate_mesh!(mesh)

    tipo = mesh.eltype
    offset = mesh.offset

    # Elemento não constante
    if tipo != 0 
        # Elemento descontínuo
        if offset != 0.0
            generate_desc_mesh!(mesh)
        # Elemento contínuo
        else
            generate_cont_mesh!(mesh)
        end
    # Elemento constante
    else
        generate_const_mesh!(mesh)
    end

    return mesh
end

function generate_desc_mesh!(mesh)
    eltype = mesh.eltype
    offset = mesh.offset
    nel = mesh.nelem

    csis_cont = range(-1.0,1.0,length = eltype+1)
    # Implementar detecção de erros para csis_descont (limite do offset)
    csis_descont = range(-1.0+offset, 1.0-offset, length = eltype+1)

    nnel = (eltype+1)^2

    k = calc_k(nnel)
    N = calc_N_matriz(csis_cont,csis_descont)

    mesh.IEN = reshape(1:nnel*nel,nnel,nel)
    mesh.nodes = zeros(nel*nnel,4)
    mesh.nodes[:,1] = 1:nel*nnel

    for i in 1:nel
        points_in_elem = mesh.IEN_geo[:,i]
        mesh.nodes[mesh.IEN[:,i],2:end] = generate_nodes_in_elem(N,mesh.points[points_in_elem,2:end],k)
    end

    return mesh
end

function generate_cont_mesh!(mesh)

    return mesh
end

function generate_const_mesh!(mesh)
    
    return mesh
end

function generate_nodes_in_elem(N,p,k)

    pd = zeros(size(p,1),3)
    for i in 1:length(k)
        pd[i,:] = N[k[i][1], k[i][2],:]'*p
    end
    return pd
end
