function remove_EE!(mesh,solver_var)

    ee = findall(x->x==0mesh.bc[:,1])

    return mesh, solver_var
end