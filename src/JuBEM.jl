module JuBEM

# Write your package code here.

using Revise, DelimitedFiles, FastGaussQuadrature



include("typedefinitions.jl")
export material_table_type, mesh_type, problem_type, solver_var_type, cmplx_consts

include("input.jl")
export read_msh

include("formfunctions.jl")

end
