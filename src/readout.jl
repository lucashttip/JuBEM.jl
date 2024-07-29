function readvars_out(filename)
    filename = string(filename,".h5")
    fid = h5open(filename, "r")

    material = Material[]
    mesh = Mesh()
    problem = Problem()
    solver_var = Assembly()

    groupnames = keys(fid)

    vars = [mesh, problem, solver_var]

    for var in vars
        varname = replace(string(typeof(var)),"_type"=>"")
        if varname in groupnames
            parse_var!(fid[varname], var)
        end
    end

    nm, names = get_nums(groupnames,"material")

    nm = Int.(nm)

    for n in nm
        push!(material,Material(0,0,0,0))
        parse_var!(fid[names[n]], material[n])
    end


    close(fid)

    return mesh, material, problem, solver_var
end

function getnoderes_out(filename,node)
    fid = h5open(filename, "r")

    groupnames = keys(fid)
    freqs, names = get_nums(groupnames,"freq")
    p = sortperm(freqs)
    freqs = freqs[p]
    names = names[p]

    nf = length(freqs)
    u = zeros(ComplexF64,nf,3)
    t = zeros(ComplexF64,nf,3)

    for n in 1:nf
        u[n,:] = read(fid[names[n]],"u")[node,:]
        t[n,:] = read(fid[names[n]],"t")[node,:]
    end


    close(fid)

    return u,t,freqs
end

function geturb_out(filename,dim)
    fid = h5open(filename, "r")

    groupnames = keys(fid)
    freqs, names = get_nums(groupnames,"freq")
    p = sortperm(freqs)
    freqs = freqs[p]
    names = names[p]

    nf = length(freqs)
    u = zeros(ComplexF64,nf)

    for n in 1:nf
        u[n] = read(fid[names[n]],"urb")[dim]
    end


    close(fid)

    return u,freqs
end

function getflex_out(filename)
    fid = h5open(filename, "r")

    groupnames = keys(fid)
    freqs, names = get_nums(groupnames,"freq")
    p = sortperm(freqs)
    freqs = freqs[p]
    names = names[p]

    nf = length(freqs)
    N = zeros(ComplexF64,nf,6,6)

    for n in 1:nf
        N[n,:,:] = read(fid[names[n]],"N")
    end


    close(fid)

    return N,freqs
end

function getfreqres_out(filename, freq)
    filename = string(filename,".h5")
    fid = h5open(filename, "r")

    groupname = string("freq_",freq)
    
    u = read(fid[groupname],"u")
    t = read(fid[groupname],"t")

    close(fid)

    return u,t
end

function parse_var!(group,var)
    dsnames = keys(group)

    for dsname in dsnames
        setfield!(var,Symbol(dsname),read(group,dsname))
    end

    return var
end

function get_nums(groupnames, str)

    nums = Float64[]
    names = []

    for name in groupnames
        if occursin(str,name)
            num = replace(name,string(str,"_")=>"")
            push!(nums,parse(Float64,num))
            push!(names, name)
        end

    end

    return nums, names
end

function get_value_out(filename,groupname,fieldname)
    filename = string(filename,".h5")
    fid = h5open(filename, "r")

    x = read(fid[groupname],fieldname)

    close(fid)

    return x
end