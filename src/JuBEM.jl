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
include("integrationrules.jl")

include("assembly.jl")
export statics_assembly, dynamics_assembly!

# include("solver.jl") # Only for future references

include("applyBC.jl")
export applyBC_rb, applyBC_simple, returnut_rb, returnut_simple

include("EE.jl")
export remove_EE!
# include("rbmotion.jl")

# include("posprocessor.jl") # TODO: NEEDS CHANGES FOR IT TO WORK AGAIN

include("plotting.jl")
export plot_disp, plot_meshtags

include("writevtk.jl")
export writevtk

include("output.jl")
export output_vars, output_solution

include("readout.jl")
export readvars_out, getnoderes_out, getfreqres_out, geturb_out, get_value_out,getflex_out, getfreqs_out

include("findmind.jl")

end
