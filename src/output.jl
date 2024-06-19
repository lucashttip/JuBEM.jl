function output_vars(filename, mesh::Mesh, problem::Problem, materials::Array{Material,1}, solver_var::Assembly = Assembly())
    
    vars = [mesh, problem, solver_var]

    h5open(filename, "w") do file
        for var in vars
            groupname = string(typeof(var))
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
    filename = string(filename)

    h5open(filename, "cw") do file
        g = file[groupname] # access group
        g[label] = time
    end
end

function output_solution(filename, sol::Solution)
    
    freq = sol.freq
    u = sol.u
    t = sol.t
    urb = sol.urb
    time = sol.time

    groupname = string("freq_",freq)

    h5open(filename, "cw") do file
        g = create_group(file, groupname) # create a group
        g["u"] = u
        g["t"] = t
        g["freq"] = freq
        g["time"] = time
        if !isempty(urb)
            g["urb"] = urb
        end
    end
end

function output_freqflex(filename, freq, N)

    groupname = string("freq_",freq)

    h5open(filename, "cw") do file
        g = create_group(file, groupname) # create a group
        g["N"] = N
    end
end
