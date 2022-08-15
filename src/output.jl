function output_vars_h5(filename, mesh::mesh_type, problem::problem_type, solver_var::solver_var_type, materials::Array{material_table_type,1})
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
            fields = fieldnames(typeof(var))
            for field in fields
                field_str = string(field)
                g[field_str] = getfield(var,field)     
            end   
            i +=1
        end
    end
end



function output_freq_h5(filename, u, t, freq)

    groupname = string("freq_",freq)
    filename = string(filename,".h5")

    h5open(filename, "cw") do file
    g = create_group(file, groupname) # create a group
    g["u"] = u                  
    g["t"] = t                  
    end
end

