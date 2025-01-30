function remove_EE!(mesh::Mesh,assembly::Assembly,problem::Problem)

    bctypeidxee = findall(x->x == 0 , problem.bctype[:,2])
    tagidxee = findall(x->x ∈ bctypeidxee , problem.taginfo[:,2])
    elemidxee = findall(x->x ∈ tagidxee, mesh.tag)
    elemidxnonee = setdiff(1:mesh.nelem, elemidxee)
    elemgeoidxnonee = unique(abs.(mesh.EEN[elemidxnonee]))
    # @infiltrate

    nnel = size(mesh.IEN,1)
    nelem_new = length(elemidxnonee)
    nelemgeo_new = length(elemgeoidxnonee)
    nnodes_new = mesh.nnodes - nnel*length(elemidxee)

    ID_new = zeros(Int32,size(mesh.ID))
    LM_new = zeros(Int32,3*nnel,nelem_new)
    IEN_new = zeros(Int32,nnel,nelem_new)
    IEN_geo_new = zeros(Int32,nnel,nelemgeo_new)
    tag_new = zeros(Int32,nelem_new)
    EEN_new = zeros(Int32,nelem_new)

    
    # Remover as colunas de LM, IEN, IEN_geo
    mesh.nelem = nelem_new

    posID = 1
    posLM = 1
    nodesidx = []

    for e in 1:nelemgeo_new
        IEN_geo_new[:,e] = mesh.IEN_geo[:,elemgeoidxnonee[e]]
        idxsEEN = findall(elemgeoidxnonee[e] .== abs.(mesh.EEN))
        # @infiltrate
        mesh.EEN[idxsEEN] = sign.(mesh.EEN[idxsEEN]).*e
        # @infiltrate
    end

    for e in 1:nelem_new
        IEN_new[:,e] = mesh.IEN[:,elemidxnonee[e]]
        # elemgeo = abs(mesh.EEN[elemidxnonee[e]])
        # IEN_geo_new[:,e] = mesh.IEN_geo[:,elemgeo]
        EEN_new[e] = mesh.EEN[elemidxnonee[e]]
        for i in 1:nnel
            nodeidx = IEN_new[i,e]
            if nodeidx ∉ nodesidx
                ID_new[:,nodeidx] = [posID, posID+1, posID+2]
                push!(nodesidx,nodeidx)
                posID = posID+3
            end
            idxs = 3*(i-1)+1:3*i
            LM_new[idxs,e] = [posLM, posLM+1, posLM+2]
            posLM = posLM + 3
        end
    end

    # EEN_new = mesh.EEN[elemidxnonee]
    mesh.EEN = EEN_new
    tag_new = mesh.tag[elemidxnonee]

    nDofsnew = maximum(LM_new)

    G = zeros(nDofsnew, nDofsnew)
    H = zeros(nDofsnew, nDofsnew)
    # Remover as colunas de G e de H

    for e in 1:nelem_new
        for se in 1:nelem_new
            G[LM_new[:,se][:],LM_new[:,e][:]] = assembly.G[mesh.LM[:,elemidxnonee[se]][:],mesh.LM[:,elemidxnonee[e]][:]]
            H[LM_new[:,se][:],LM_new[:,e][:]] = assembly.H[mesh.LM[:,elemidxnonee[se]][:],mesh.LM[:,elemidxnonee[e]][:]]
        end
    end
 

    assembly.G = G
    assembly.H = H

    # Update mesh
    mesh.nelem = nelem_new
    mesh.nnodes = nnodes_new
    mesh.ID = ID_new
    mesh.LM = LM_new
    mesh.IEN = IEN_new
    mesh.IEN_geo = IEN_geo_new
    mesh.tag = tag_new
    mesh.EEN = EEN_new

    # Definicao de bodies
    nbodies = maximum(problem.taginfo[:,4])
    bodies_new = []
    for b in 1:nbodies
        tagidxs = findall(problem.taginfo[:,4].==b)
        elems = findall([x ∈ tagidxs for x in mesh.tag])
        push!(bodies_new,elems)
    end

    mesh.bodies = bodies_new
    
    return mesh, assembly
end



function remove_EE_mesh!(mesh)

    nonee = findall(x->x != 0 , mesh.bc[:,1])
    
    # Remover as colunas de LM, IEN, IEN_geo
    mesh.nelem = length(nonee)
    mesh.LM = mesh.LM[:,nonee]
    mesh.IEN = mesh.IEN[:,nonee]
    mesh.IEN_geo = mesh.IEN_geo[:,nonee]
    mesh.bc = mesh.bc[nonee,:]
    mesh.bcvalue = mesh.bcvalue[nonee,:]

    # Remover as colunas de ID, nodes, bc e bcvalue
    mesh.nnodes = length(unique(mesh.IEN[:]))  #Apenas para descontínuo
    mesh.ID = mesh.ID[:,unique(sort(mesh.IEN[:]))]
    mesh.nodes = mesh.nodes[unique(sort(mesh.IEN[:])),:]

    return mesh
end