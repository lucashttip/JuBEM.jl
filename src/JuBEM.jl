module JuBEM

# Write your package code here.

using Revise, DelimitedFiles, FastGaussQuadrature, LinearAlgebra, Infiltrator
using WriteVTK, HDF5
using ProgressMeter
import Base.:(==)

include("typedefinitions.jl")
export Material, Mesh, Problem, Svar, cmplx_consts

include("input.jl")
export read_msh

include("shapefunctions.jl")
export calc_N_matrix, calc_N_gen
export calc_csis_grid

include("elementsubdivision.jl")

include("telles.jl")

include("derivedata.jl")
export derive_data!

include("generatemesh.jl")
export generate_nodes_in_elem, generate_mesh!, generate_points_in_elem

include("geometry.jl")

include("fundamentalsolutions.jl")
export calc_funsol_static, calc_funsol_dynamic

include("integrationdynamics.jl")
include("integrationstatics.jl")
include("integrationconstants.jl")

include("solverstatics.jl")
include("solverdynamics.jl")

include("solver.jl")
export calc_GH!, solvestatic, solvedynamic, solve, solve_flex_dyn

include("applyBC.jl")
export applyBC!, calc_utpoints

include("EE.jl")
include("rbmotion.jl")

include("posprocessor.jl")
export calc_interior_static, calc_interior_static_const

include("plotting.jl")
export view_mesh, animate_res_freq

include("writevtk.jl")
export writevtk

include("output.jl")
export output_freq_h5, output_vars_h5

include("readout.jl")
export readvars_out, getnoderes_out, getfreqres_out, geturb_out, get_value_out,getflex_out

include("findmind.jl")
include("distfuncs.jl")

end
