# type :: material_table_type
# real(kind=dp) :: Ge         !> Shear modulus
# real(kind=dp) :: Nu         !> Poisson coefficient
# real(kind=dp) :: Dam        !> Damping coefficient
# real(kind=dp) :: Rho        !> Specific mass
# complex(kind=dp) :: zGe     !> Shear modulus (complex)
# complex(kind=dp) :: zSwv    !> S-wave velocity (complex)
# complex(kind=dp) :: zPwv    !> P-wave velocity (complex)
# end type material_table_type

# "Type that holds the information of material tables"

abstract type JuBEMtypes end

mutable struct material_table_type <: JuBEMtypes
    Ge :: Float64
    Nu :: Float64
    Dam :: Float64
    Rho :: Float64
    zGe :: ComplexF64
    zSwv :: ComplexF64
    zPwv :: ComplexF64
    material_table_type(Ge,Nu,Dam,Rho) = new(Ge,Nu,Dam,Rho, 0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im)
end

# type :: mesh_type
# integer :: nelem       !> Number of boundary elements
# integer :: npoints      !> Number of geometrical points
# integer :: eltype       !> Type of elements (constant, linear, quadratic)
# integer, dimension (:,:), allocatable :: ID             !> Matrix containing degrees of freedom
# integer, dimension (:,:), allocatable :: IEN_geo        !> Geometrical incidence matrix
# integer, dimension(:), allocatable :: IEN               !> Incidence matrix
# integer, dimension(:,:), allocatable :: LM              !> LM matrix
# integer, dimension(:), allocatable :: bc                !> Array of the boundary conditions tyoe (1 = u known; 2 = t known, 3 = rigid body)
# integer, dimension(:), allocatable :: mat                !> material of the element
# real(kind=dp) :: offset                                !> mesh nodes offset (if 0, mesh is continuous)
# real(kind=dp), dimension(:,:), allocatable :: points    !> Array of the geometrical points
# real(kind=dp), dimension(:,:), allocatable :: nodes     !> Array of the physical nodes
# real(kind=dp), dimension(:), allocatable :: bcvalue     !> Array of the boundary conditions values
# complex(kind=dp), dimension(:), allocatable :: zbcvalue     !> Array of the boundary conditions values (complex)
# end type mesh_type
# "Type that holds the information of  3D general mesh"
mutable struct mesh_type <: JuBEMtypes
    nelem :: Int64
    npoints :: Int64
    nnodes :: Int64
    eltype :: Int8
    offset :: Float64
    ID :: Array{Int32,2}    
    IEN_geo :: Array{Int32,2}   
    IEN :: Array{Int32,2}   
    LM :: Array{Int32,2}    
    bc :: Array{Int16,2}    
    material :: Array{Int16,1}  
    points :: Array{Float64,2}  
    nodes :: Array{Float64,2}   
    bcvalue :: Array{Float64,2} 
    zbcvalue
    forces :: Array{Float64,2}
    mesh_type() = new(0,0,0,0,0,Array{Int32,2}(undef,0,0), Array{Int32,2}(undef,0,0), Array{Int32,2}(undef,0,0), Array{Int32,2}(undef,0,0),
    Array{Int16,2}(undef,0,0), Array{Int16,1}(undef,0), Array{Float64,2}(undef,0,0), Array{Float64,2}(undef,0,0), 
    Array{Float64,2}(undef,0,0), Array{Float64,1}(undef,0), Array{Float64,2}(undef,0,0))
end

# !> @brief Type that holds the information of 
# !> the problem (frequencies)
# type :: problem_type
# integer, dimension(:), allocatable :: nFr       !> Number of frequencies
# real(kind=dp) :: frequency
# real(kind=dp), dimension(:), allocatable :: fr_range    !> Frequency ranges
# real(kind=dp), dimension(:), allocatable :: frequencies !> Frequencies array
# end type problem_type
mutable struct problem_type <: JuBEMtypes
    nFr :: Array{Int32,1}
    fr_range :: Array{Float64,1}
    frequencies :: Array{Float64,1}
    problem_type() = new([],[],[])
end

# !> @brief Type that holds the information of 
# !> the solver (gauss points and matrices)
# type :: solver_var_type
# integer :: nGP  !> Number of Gauss points for integration
# real(kind=dp), dimension(:), allocatable :: csi, omega !> Gauss points and weights
# complex(kind=dp), dimension(:,:), allocatable :: H  !> Complex H matrix
# complex(kind=dp), dimension(:,:), allocatable :: G  !> Complex G matrix
# complex(kind=dp), dimension(:), allocatable :: zvetsol
# complex(kind=dp), dimension(:,:), allocatable :: zma
# end type solver_var_type
mutable struct solver_var_type <: JuBEMtypes
    nGP :: Int16
    csi :: Array{Float64,1}
    omega :: Array{Float64,1}
    H :: Array{Float64,2}
    G :: Array{Float64,2}
    zH :: Array{ComplexF64,2}
    zG :: Array{ComplexF64,2}
    ma
    zvetsol
    solver_var_type() = new(0,[],[],Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{ComplexF64,2}(undef,0,0),Array{ComplexF64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,1}(undef,0))
end


# !> @brief Type that holds complex constants for the dynamics fuyndamental solution
# type :: cmplx_consts
# complex(kind=dp) :: zWi
# complex(kind=dp) :: zC0
# complex(kind=dp) :: zC1
# complex(kind=dp) :: zC2
# complex(kind=dp) :: zKP
# complex(kind=dp) :: zKS
# end type cmplx_consts
mutable struct cmplx_consts <: JuBEMtypes
    zWi :: ComplexF64
    zC0 :: ComplexF64
    zC1 :: ComplexF64
    zC2 :: ComplexF64
    zKP :: ComplexF64
    zKS :: ComplexF64

    function cmplx_consts(material :: material_table_type, fr :: Float64)
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
    csi :: Array{Float64,2}
    omega :: Array{Float64,1}
end

mutable struct integration_rules_type <: JuBEMtypes

    npgs :: Array{Int64,1}
    dists :: Array{Float64,1}
    gp :: Array{gauss_points_type,1}
    npg_near :: Int64
    npg_sing :: Int64
    gp_near :: gauss_points_type
    integration_rules_type(npgs::Array{Int64,1},dists::Array{Float64,1}, npg_near::Int64, npg_sing::Int64) = new(npgs,dists,gauss_points_type[], npg_near,npg_sing)
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

function copy(mesh::mesh_type)
    mesh2 = mesh_type()


    fields = fieldnames(mesh_type)

    for field in fields
        f1 = getfield(mesh,field)
        setfield!(mesh2,field,f1)
    end
    return mesh2
end