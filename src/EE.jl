function remove_EE!(mesh,solver_var)

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

    # Remover as colunas de G e de H
    solver_var.G = solver_var.G[mesh.LM[:],mesh.LM[:]]
    solver_var.H = solver_var.H[mesh.LM[:],mesh.LM[:]]


    return mesh, solver_var
end