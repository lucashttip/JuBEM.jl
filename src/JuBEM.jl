module JuBEM

# Write your package code here.

using Revise, DelimitedFiles, FastGaussQuadrature



include("typedefinitions.jl")
export material_table_type, mesh_type, problem_type, solver_var_type, cmplx_consts

include("input.jl")
export read_msh

include("formfunctions.jl")
export calc_N_matriz, calc_dNdcsi_matriz, calc_k,  calc_G, remap_N

include("derivedata.jl")
export derive_data!

include("generatemesh.jl")
export generate_nodes_in_elem, generate_mesh!

include("geometry.jl")
export calc_normal

end
