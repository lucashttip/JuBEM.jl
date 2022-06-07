

function generate_mesh!(mesh::mesh_type)

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
    # Generate nodes, IEN, ID and LM

    eltype = mesh.eltype
    offset = mesh.offset
    nel = mesh.nelem
    
    csis_cont = range(-1.0,1.0,length = eltype+1)
    # Implementar detecção de erros para csis_descont (limite do offset)
    csis_descont = range(-1.0+offset, 1.0-offset, length = eltype+1)
    
    nnel = (eltype+1)^2
    nnodes = nel*nnel

    k = calc_k(nnel)
    N = calc_N_matrix(csis_cont,csis_descont)

    mesh.nnodes = nnodes
    mesh.ID = reshape(1:3*nnodes,3,nnodes)
    mesh.LM = reshape(1:3*nnodes,3*nnel,nel)
    mesh.IEN = reshape(1:nnodes,nnel,nel)
    mesh.nodes = zeros(nnodes,4)
    mesh.nodes[:,1] = 1:nnodes

    for i in 1:nel
        points_in_elem = mesh.IEN_geo[:,i]
        mesh.nodes[mesh.IEN[:,i],2:end] = generate_nodes_in_elem(N,mesh.points[points_in_elem,2:end],k)
    end

    return mesh
end

function generate_cont_mesh!(mesh)
    # Generate nodes, IEN, ID and LM

    return mesh
end

function generate_const_mesh!(mesh)
    # Generate nodes, IEN, ID and LM
    nel = mesh.nelem

    mesh.nodes = zeros(nel,4)
    mesh.IEN = zeros(1,nel)
    mesh.ID = zeros(3,nel)
    mesh.LM = zeros(3,nel)
    mesh.nnodes = nel

    mesh.nodes[:,1] = 1:nel
    mesh.IEN[:] = 1:nel
    mesh.ID = reshape(1:3*nel,size(mesh.ID))
    mesh.LM = mesh.ID

    for e in 1:nel
        points = mesh.points[mesh.IEN_geo[:,e],2:end]
        mesh.nodes[e,2:end] = sum(points,dims=1)./4
    end


    return mesh
end


function generate_nodes_in_elem(N,p,k)

    pd = zeros(size(p,1),3)
    for i in 1:length(k)
        pd[i,:] = N[k[i][1], k[i][2],:]'*p
    end
    return pd
end

function generate_points_in_elem(N,p)

    pd = zeros(size(N,1),size(N,2),size(p,2))
    k = 1
    for i in 1:size(N,1)
        for j in 1:size(N,2)
            pd[i,j,:] = N[i, j,:]'*p
            k = k+1
        end
    end
    
    return pd
end
