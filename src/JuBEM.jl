module JuBEM

# Write your package code here.

using Revise, DelimitedFiles, FastGaussQuadrature, LinearAlgebra, Infiltrator
using WriteVTK, HDF5
using ProgressMeter
import Base.:(==)

include("typedefinitions.jl")
export Material, Mesh, Problem, Assembly, Solution, cmplx_consts

include("input.jl")
export read_msh, read_problem

include("shapefunctions.jl")
export calc_N_gen

include("elementsubdivision.jl")
include("telles.jl")

include("derivedata.jl")
export derive_data!

include("generatemesh.jl")
export generate_mesh!

include("geometry.jl")

include("fundamentalsolutions.jl")
export calc_funsol_static, calc_funsol_dynamic

include("integration.jl")
include("integrationdynamics.jl")
include("integrationstatics.jl")
include("integrationrules.jl")
include("integrationconstants.jl")

include("assembly.jl")
export statics_assembly, dynamics_assembly!

include("solverstatics.jl")
include("solverdynamics.jl")

include("solver.jl")
export calc_GH!, solvestatic, solvedynamic, solve, solve_flex_dyn

include("applyBC.jl")
export applyBC!, calc_utpoints

include("EE.jl")
# include("rbmotion.jl")

include("posprocessor.jl")
export calc_interior_static, calc_interior_static_const

include("plotting.jl")
export plot_disp
# export view_mesh, animate_res_freq

include("writevtk.jl")
export writevtk

include("output.jl")
export output_vars, output_solution

include("readout.jl")
export readvars_out, getnoderes_out, getfreqres_out, geturb_out, get_value_out,getflex_out

include("findmind.jl")

end
