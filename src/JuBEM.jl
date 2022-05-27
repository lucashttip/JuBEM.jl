module JuBEM

# Write your package code here.

using Revise, DelimitedFiles, FastGaussQuadrature, LinearAlgebra, Infiltrator


include("typedefinitions.jl")
export material_table_type, mesh_type, problem_type, solver_var_type, cmplx_consts

include("input.jl")
export read_msh

include("formfunctions.jl")
export calc_N_matrix, calc_dNdcsi_matrix, calc_k, calc_dNdeta_matrix, calc_G, remap_N, calc_N_matrix_descont, calc_u_in_N

include("derivedata.jl")
export derive_data!

include("generatemesh.jl")
export generate_nodes_in_elem, generate_mesh!, generate_points_in_elem

include("geometry.jl")
export calc_n_J_matrix

end
