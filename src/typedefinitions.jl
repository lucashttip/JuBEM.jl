
"""
Abstract type that is parent to all other JuBEm types.
Serves the purpose of facilitating functions such as copy
"""
abstract type JuBEMtypes end

"""
    Material(Ge,Nu,Dam,Rho)

It is the data structure that contains properties of a particular material.
It's a mutable struct.

## Arguments

  - `Ge`: Shear modulus.
  - `Nu`: Poisson coefficient.
  - `Dam`: Damping coefficient.
  - `Rho`: Density.

## Other fields

  - `zGe`: Complex property of Shear Modulus for viscoelasticity.
  - `zSwv`: Complex shear wave velocity of viscoelastic material.
  - `zPwv`: Complex pressure wave velocity of viscoelastic material.
"""
mutable struct Material <: JuBEMtypes
    Ge :: Float64
    Nu :: Float64
    Dam :: Float64
    Rho :: Float64
    C_stat :: Vector{Float64}
    zGe :: ComplexF64
    zSwv :: ComplexF64
    zPwv :: ComplexF64
    Material(Ge,Nu,Dam,Rho) = new(Ge,Nu,Dam,Rho, Float64[], 0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im)
end

"""
    Mesh

It is the data structure that contains all the mesh data. Information need to be read from a .msh file generated from Gmsh. Can be created empty by calling `Mesh()`.
It's a mutable struct.

## Fields

  ### Block 1: Variables that hold integer information of whole mesh
  - `nelem` :: Int64
  - `npoints` :: Int64
  - `nnodes` :: Int64
  - `eltype` :: Int8
  
  ### Block 2: Variables that define the mesh
  - `offset` :: Float64
  - `ID` :: Array{Int32,2}    
  - `IEN_geo` :: Array{Int32,2}   
  - `IEN` :: Array{Int32,2}   
  - `LM` :: Array{Int32,2}    
  - `points` :: Array{Float64,2}  
  - `nodes` :: Array{Float64,2}
  - `tag` :: Array{Int16,1} This associates elements with material and bc
"""
mutable struct Mesh <: JuBEMtypes
    # Block 1: Variables that hold integer information of whole mesh
    nelem :: Int64
    npoints :: Int64
    nnodes :: Int64
    eltype :: Int8

    ## Block 2: Variables that define the mesh
    offset :: Float64
    ID :: Array{Int32,2}    
    IEN_geo :: Array{Int32,2}   
    IEN :: Array{Int32,2}   
    LM :: Array{Int32,2}    
    points :: Array{Float64,2}  
    nodes :: Array{Float64,2}
    
    ## Block 3: Variables that define material and boundary conditions
    tag :: Array{Int16,1}  
    tagnames :: Array{String,1}

    Mesh() = new(
        0,0,0,0,    # 
        0,  # This is the start of Block 2
        Array{Int32,2}(undef,0,0),  # ID
        Array{Int32,2}(undef,0,0),  # IEN_geo
        Array{Int32,2}(undef,0,0),  # IEN
        Array{Int32,2}(undef,0,0),  # LM
        Array{Float64,2}(undef,0,0),# points
        Array{Float64,2}(undef,0,0),# nodes
        Int16[],    # tag This is the start of Block 3
        String[]                          # tagnames
    )
end

"""
    Problem
Type that holds the information of the problem (frequencies). Can be created empty by calling `Problem()`.
It's a mutable struct.

## Fields

  - `bctype` :: Array{Int16,2} - Contains information of boundary condition of each element. Since each node has 3 degrees of freedom, this array is size nelem x 3. 0 correspond to enclosing element, 1 is u-type bc, 2 is t-type bc and 3 is rigid-body
  - `bcvalue` :: Array{Float64,2} - Contains information of the boundary condition value. Has no meaning for enclosing elements, is prescribed value of u and t for these and for rigid body is the position of geometric center of rigid body.
  - `taginfo` :: Array{Int16,2}
  - `forces` :: Array{Float64,2}
  - `nFr` :: Array{Int32,1}
  - `fr_range` :: Array{Float64,1}
  - `frequencies` :: Array{Float64,1}
 
"""
mutable struct Problem <: JuBEMtypes
    name :: String
    bctype :: Array{Int16,2}
    bcvalue :: Array{Float64,2}
    forces :: Array{Float64,2}
    taginfo :: Array{Int16,2}
    nFr :: Array{Int32,1}
    fr_range :: Array{Float64,1}
    frequencies :: Array{Float64,1}
    Problem() = new("",Array{Int16,2}(undef,0,0),
    Array{Float64,2}(undef,0,0),
    Array{Float64,2}(undef,0,0),
    Array{Int16,2}(undef,0,0),
    [],[],[])
end


"""
    Assembly
Type that holds the information of the solver (gauss points and matrices).
Can be created empty by calling `Assembly()`.
It's a mutable struct.

## Fields

  - `nGP` :: Int16
  - `csi` :: Array{Float64,1}
  - `omega` :: Array{Float64,1}
  - `H` :: Array{Float64,2}
  - `G` :: Array{Float64,2}
  - `zH` :: Array{ComplexF64,2}
  - `zG` :: Array{ComplexF64,2}
  - `ma`
  - `zvetsol`
 
"""
mutable struct Assembly <: JuBEMtypes
    nGP :: Int16
    csi :: Array{Float64,1}
    omega :: Array{Float64,1}
    H :: Array{Float64,2}
    G :: Array{Float64,2}
    zH :: Array{ComplexF64,2}
    zG :: Array{ComplexF64,2}
    ma
    zvetsol
    Assembly() = new(0,[],[],Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{ComplexF64,2}(undef,0,0),Array{ComplexF64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,1}(undef,0))
end

"""
    Solution{T} <: JuBEMtypes
Type that holds the the solution of the simulation.
Can be created empty by calling `Assembly()`.
It's a mutable struct.

## Fields

  - `u` :: Array{T,2}
  - `t` :: Array{T,2}
  - `urb` :: Array{T,2}

"""
mutable struct Solution{T} <: JuBEMtypes
    u :: Array{T,2}
    t :: Array{T,2}
    urb :: Array{T,2}
end

"""
    cmplx_consts(material :: Material, fr :: Float64)
Type that holds complex constants for the dynamics fundamental solution.
Can be created empty by calling `Assembly()`.
It's a mutable struct.


## Arguments

  - `material`: Data of type Material
  - `fr`: frequency

## Fields

  - `zWi` :: ComplexF64
  - `zC0` :: ComplexF64
  - `zC1` :: ComplexF64
  - `zC2` :: ComplexF64
  - `zKP` :: ComplexF64
  - `zKS` :: ComplexF64
"""
mutable struct cmplx_consts <: JuBEMtypes
    zWi :: ComplexF64
    zC0 :: ComplexF64
    zC1 :: ComplexF64
    zC2 :: ComplexF64
    zKP :: ComplexF64
    zKS :: ComplexF64

    function cmplx_consts(material :: Material, fr :: Float64)
        zWi= 0.0 + fr*im
        zC0=(1.0 + 0.0im)/((4.0*pi + 0.0im)*material.zGe)
        zC1=(material.zPwv/material.zSwv)^2.0
        zC2=(material.zSwv/material.zPwv)^2.0
        zKP=-zWi/material.zPwv
        zKS=-zWi/material.zSwv
        new(zWi, zC0, zC1, zC2, zKP, zKS)
    end
end

mutable struct gauss_points_type <: JuBEMtypes
    csis :: Array{Float64,2}
    omegas :: Array{Float64,1}
end

"""
    Rules()

## Fields
- `csis_points` :: Vector{Float64}
- `csis_nodes` :: Vector{Float64}
- `npgs` :: Array{Int64,1}
- `dists` :: Array{Float64,1}
- `gp` :: Array{gauss_points_type,1}
- `npg_near` :: Int64
- `near_strat` :: Symbol
- `N_aux` :: Array{Float64,2}
- `gp_near` :: gauss_points_type
- `npg_sing` :: Int64
- `csis_sing` :: Array{gauss_points_type,1}
"""
mutable struct Rules <: JuBEMtypes

    csis_points :: Vector{Float64}
    csis_nodes :: Vector{Float64}
    npgs :: Array{Int64,1}
    dists :: Array{Float64,1}
    gp :: Array{gauss_points_type,1}
    npg_near :: Int64
    npg_sing :: Int64
    csis_sing :: Array{gauss_points_type,1}
    near_strat :: Symbol
    N_aux :: Array{Float64,2}
    gp_near :: gauss_points_type
    Rules(csis_points:: Vector{Float64}, csis_nodes:: Vector{Float64}, npgs::Array{Int64,1},dists::Array{Float64,1}, npg_near::Int64, npg_sing::Int64,near_strat::Symbol, N_aux :: Array{Float64,2}) = new(csis_points,csis_nodes,npgs,dists,gauss_points_type[], npg_near,npg_sing,gauss_points_type[],near_strat,N_aux)
end



==(a::JuBEMtypes,b::JuBEMtypes) = begin
    if typeof(a) != typeof(b)
        return false
    end
    fields = fieldnames(typeof(a))

    for field in fields
        f1 = getfield(a,field)
        f2 = getfield(b,field)
        if f1 != f2
            return false
        end
    end
    return true
end

function copy(mesh::Mesh)
    mesh2 = Mesh()


    fields = fieldnames(Mesh)

    for field in fields
        f1 = getfield(mesh,field)
        setfield!(mesh2,field,f1)
    end
    return mesh2
end