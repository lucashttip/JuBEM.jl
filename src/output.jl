function output_vars_h5(filename, mesh::Mesh, problem::Problem, solver_var::Assembly, materials::Array{Material,1})
    filename = string(filename,".h5")

    vars = [mesh, problem, solver_var]

    h5open(filename, "w") do file
        for var in vars
            groupname = replace(string(typeof(var)),"_type"=>"")
            g = create_group(file, groupname) # create a group
            fields = fieldnames(typeof(var))
            for field in fields
                field_str = string(field)
                g[field_str] = getfield(var,field)     
            end   
        end
        i = 1
        for material in materials
            groupname = string("material_",i)
            g = create_group(file, groupname) # create a group
            fields = fieldnames(typeof(material))
            for field in fields
                field_str = string(field)
                g[field_str] = getfield(material,field)     
            end   
            i +=1
        end
        create_group(file,"times")
    end
end

function output_time(filename,time,label)
    groupname = "times"
    filename = string(filename,".h5")

    h5open(filename, "cw") do file
        g = file[groupname] # access group
        g[label] = time
    end
end

function output_freq_h5(filename, freq, u, t, urb = [])

    groupname = string("freq_",freq)
    filename = string(filename,".h5")

    h5open(filename, "cw") do file
        g = create_group(file, groupname) # create a group
        g["u"] = u
        g["t"] = t
        if !isempty(urb)
            g["urb"] = urb
        end
    end
end

function output_freqflex_h5(filename, freq, N)

    groupname = string("freq_",freq)
    filename = string(filename,".h5")

    h5open(filename, "cw") do file
        g = create_group(file, groupname) # create a group
        g["N"] = N
    end
end

function output_frb_h5(filename,u,freq)
    groupname = string("freq_",freq)
    filename = string(filename,".h5")

    h5open(filename, "cw") do file
    g = create_group(file, groupname) # create a group
    g["u"] = u                                  
    end
end