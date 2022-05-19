module JuBEM

# Write your package code here.

using Revise, DelimitedFiles



include("typedefinitions.jl")
export material_table_type, mesh_type, problem_type, solver_var_type, cmplx_consts

include("input.jl")
export readmsh, read_msh

end
