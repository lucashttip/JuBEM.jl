# incomplete

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