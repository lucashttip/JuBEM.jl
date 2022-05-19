# incomplete
# Inspirado no pacote MeshIO.jl: https://github.com/JuliaIO/MeshIO.jl/blob/master/src/io/msh.jl
@enum MSHBlockType MSHFormatBlock MSHPhysicalNamesBlock MSHNodesBlock MSHElementsBlock MSHUnknownBlock MSHMaterialBlock MSHFrequenciesBlock MSHMeshTypeBlock MSHForcesBlock MSHEntitiesBlock

function readdatafromfile(inp_file)
    fid = open(inp_file, "r")

    # Reading information about the mesh
    data = split(read(fid,String),'\n')
    close(fid)
    return data
end

function locatesubstring(sbstr,strarr)
    i = 1
    for line in strarr
        if line ==  sbstr
            return i
        end 
        i = i+1
    end
end

function readmsh(inp_file)
    data = readdatafromfile(inp_file)
    l = 1
    nmat = parse(Int,data[1])
    
    mesh = mesh_type()
    material = []
    problem = problem_type()
    solver_var = solver_var_type()


    for i in 2:1+nmat
        line = split(data[i],' ')
        Ge = parse(Float64, line[1])
        Nu = parse(Float64, line[2])
        Dam = parse(Float64, line[3])
        Rho = parse(Float64, line[4])
        material = [material; material_table_type(Ge,Nu,Dam,Rho)]
    end
    
    l = 2+nmat

    nfrs = parse(Int,data[l])

    for i in 1:nfrs
        values = parse.(Float64,split(data[l+i], ' '))
        if i == 1
            problem.fr_range = [problem.fr_range;values[1];values[2]]
        else
            problem.fr_range = [problem.fr_range;values[2]]
        end
        problem.nFr = [problem.nFr; Int32(values[3])]
    end

    l = l + nfrs +1

    solver_var.nGP = parse(Int,data[l])

    l = l+1
    mesh.offset = parse(Float64,data[l])

    l = l+1
    mesh.eltype = parse(Int,data[l])

    l = l+2
    F = parse.(Float64,data[l:l+5])

    l = locatesubstring("\$PhysicalNames", data) +1

    n_phys = parse.(Int,data[l])

    bc = zeros(n_phys,2)
    bcvalue = zeros(n_phys,3)
    mat = zeros(n_phys,2)

    bcvalue.=0.0

    for i in 1:n_phys

    end

    # # Reading which are the rigid elements
    # data = ""
    # for i in 1:nr_rigid_elements
    #     data = string(data,readline(io_in;keep=true))
    # end
    # rigid_elements = readdlm(IOBuffer(data), '\t', Int, '\n')

    # # Skipping lines
    # for i in 1:(nr_be - nr_rigid_elements+2)
    #     readline(io_in)
    # end

    # # Reading Points
    # data = ""
    # for i in 1:nr_external_point
    #     data = string(data,join(split(readline(io_in)),'\t'),'\n')
    # end
    # points = readdlm(IOBuffer(data),'\t', Float64, '\n'; skipblanks=true)

    # # Reading elements
    # data = ""
    # for i in 1:nr_be+2
    #     data = string(data,join(split(readline(io_in)),'\t'),'\n')
    # end
    # elem = Int.(readdlm(IOBuffer(data),'\t', Float64, '\n'; skipblanks=true))

    # data = ""
    # for line in readlines(io_in;keep=true)
    #     data = string(data,line)
    # end
    # boundary_conditions =readdlm(IOBuffer(data),'\t', Float64, '\n'; skipblanks=true)

    # close(io_in)
    return mesh,material,problem,solver_var

    # return mesh, material, problem, solver_var
    
end


function read_msh(inp_file)
    
    io = open(inp_file,"r")

    material = material_table_type[]
    mesh = mesh_type()
    problem = problem_type()
    solver_var = solver_var_type()

    F = zeros(6)

    while !eof(io)
        BlockType = parse_blocktype!(io)
        if BlockType == MSHMaterialBlock
            parse_materials!(io, material)
        elseif BlockType == MSHFrequenciesBlock
            parse_frequencies!(io, problem)
        elseif BlockType == MSHMeshTypeBlock
            parse_meshtype!(io, mesh, solver_var)
        elseif BlockType == MSHForcesBlock
            parse_forces!(io,F)
        elseif BlockType == MSHPhysicalNamesBlock
            bc, bcvalue, mat = parse_physicalnames(io)
        # elseif BlockType == MSHEntitiesBlock
        #     parse_entities!(io)
        # elseif BlockType == MSHNodesBlock
        #     parse_nodes!(io)
        # elseif BlockType == MSHElementsBlock
        #     parse_elements!(io)
        else
            skip_block!(io)
        end
    end

    return mesh,material,problem,solver_var
end

function parse_blocktype!(io)
    header = readline(io)
    if header == "\$Material"
        return MSHMaterialBlock
    elseif header == "\$Frequencies"
        return MSHFrequenciesBlock
    elseif header == "\$MeshType"
        return MSHMeshTypeBlock
    elseif header == "\$Forces"
        return MSHForcesBlock
    elseif header == "\$MeshFormat"
        return MSHFormatBlock
    elseif header == "\$PhysicalNames"
        return MSHPhysicalNamesBlock
    elseif header == "\$Entities"
        return MSHEntitiesBlock
    elseif header == "\$Nodes"
        return MSHNodesBlock
    elseif header == "\$Elements"
        return MSHElementsBlock
    else
        return MSHUnknownBlock
    end
end

function parse_materials!(io, material)
    nmat = parse(Int,readline(io))
    for i in 1:nmat
        Ge, Nu, Dam, Rho = parse.(Float64, split(readline(io)))
        push!(material,material_table_type(Ge,Nu,Dam,Rho))
    end
    endblock = readline(io)
    if endblock != "\$EndMaterial"
        error("expected end block tag, got $endblock")
    end
    return material
end

function parse_frequencies!(io, problem)
    nfr_range = parse(Int,readline(io))
    for i in 1:nfr_range
        fi, ff, fr_range = parse.(Float64, split(readline(io)))
        if i == 1
            problem.fr_range = [fi; ff]
        else
            problem.fr_range = [problem.fr_range; ff]
        end
        problem.nFr = [problem.nFr; fr_range]
    end
    endblock = readline(io)
    if endblock != "\$EndFrequencies"
        error("expected end block tag, got $endblock")
    end
    return problem
end

function parse_meshtype!(io, mesh, solver_var)
    solver_var.nGP = parse(Int16, readline(io))
    mesh.offset = parse(Float64, readline(io))
    mesh.eltype = parse(Int8, readline(io))
    endblock = readline(io)
    if endblock != "\$EndMeshType"
        error("expected end block tag, got $endblock")
    end
    return mesh, solver_var
end

function parse_forces!(io, F)
    for i in 1:6
        F[i] = parse(Float64, readline(io))
    end
    if endblock != "\$EndForces"
        error("expected end block tag, got $endblock")
    end
    return mesh, solver_var
    return F
end

function parse_physicalnames(io)
    
    nphys = parse(Int,readline(io))

    bc = zeros(nphys,2)
    bcvalue = zeros(nphys,3)
    mat = zeros(nphys,2)

    for i in 1:nphys
        data = split(replace(readline(io), "\"" =>""))
        bc[i,1] = parse(Float64, data[2])
        mat[i,:] = [parse(Float64, data[2]); parse(Float64, data[end])]
        bcvalue[i,:] = parse.(Float64,data[2:4])
        if data[3] == "u"
            bc[i,2] = 1
        elseif data[3] == "t"
            bc[i,2] = 2
        elseif data[3] == "r"
            bc[i,3] = 3
        end
    end
    if endblock != "\$EndPhysicalNames"
        error("expected end block tag, got $endblock")
    end
    return mesh, solver_var
    return bc, bcvalue, mat
end

function parse_entities!(io)
    
end

function parse_nodes!(io)
    
end

function parse_elements!(io)
    
end

function skip_block!(io)
    while true
        line = readline(io)
        if length(line) < 4
            continue
        end
        if line[1:4] == "\$End"
            break
        end
    end
    return nothing
end