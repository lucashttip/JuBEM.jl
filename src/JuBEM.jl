module JuBEM

# Write your package code here.

using Revise, DelimitedFiles



includet("typedefinitions.jl")
export material_table_type, mesh_type, problem_type, solver_var_type, cmplx_consts

includet("Input.jl")

end
