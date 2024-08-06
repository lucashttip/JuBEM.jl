function generate_mesh!(mesh::Mesh)

    order = mesh.eltype
    offset = mesh.offset

    # Elemento não constante
    if order != 0 
        # Elemento descontínuo
        if offset != 0.0
            generate_disc_mesh!(mesh)
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

function generate_disc_mesh!(mesh)
    # Generate nodes, IEN, ID and LM

    order = mesh.eltype
    offset = mesh.offset
    nel_geo = mesh.nelem_geo
    nel = mesh.nelem
    
    csis_cont = range(-1.0,1.0,length = order+1)
    # Implementar detecção de erros para csis_descont (limite do offset)
    if offset < 0.0 || offset >=1.0
        error("Offset should be >= 0.0 and < 1.0")
    end

    csis_disc = range(-1.0+offset, 1.0-offset, length = order+1)
    
    nnel = (order+1)^2
    nnodes = nel_geo*nnel

    mesh.nnodes = nnodes
    mesh.ID = reshape(1:3*nnodes,3,nnodes)
    mesh.LM = zeros(1:3*nnodes,3*nnel,nel)
    # mesh.IEN = reshape(1:nnodes,nnel,nel)
    mesh.nodes = zeros(nnodes,4)
    mesh.nodes[:,1] = 1:nnodes

    k = map_k(nnel)

    nodalcsis = zeros(nnel,2)
    for i in eachindex(k)
        nodalcsis[i,:] = [csis_disc[k[i][1]], csis_disc[k[i][2]]]
    end

    N = calc_N_gen(csis_cont, nodalcsis;dg=:N)

    gdl = 1
    for i in 1:nel_geo
        ## TODO: PAREI AQUI

        points_in_elem = mesh.IEN_geo[:,i]
        mesh.nodes[mesh.IEN[:,i],2:end] = N*mesh.points[points_in_elem,2:end]
    end

    return mesh
end

function generate_cont_mesh!(mesh)
    
    order = mesh.eltype
    # Generate nodes, IEN, ID and LM
    nel = mesh.nelem
    nnel = (order+1)^2
    nnodes = mesh.npoints

    mesh.nnodes = nnodes
    mesh.nodes = mesh.points
    mesh.IEN = mesh.IEN_geo
    mesh.ID = zeros(3,nnodes)
    mesh.LM = zeros(3*nnel,nel)
    
    mesh.ID = reshape(1:3*nnodes,size(mesh.ID))

    for e in 1:nel
            mesh.LM[:,e] = mesh.ID[:,mesh.IEN[:,e]][:]
    end

    return mesh

end

function generate_const_mesh!(mesh)
    # Generate nodes, IEN, ID and LM
    nel = mesh.nelem
    mesh.offset = 1.0
    
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

function generate_points_in_elem(N,p)
   
    pd = N*p

    return pd
    
end