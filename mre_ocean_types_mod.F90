!> Shared derived type definitions for MREs targeting nvfortran 26.3 compiler bugs.
!!
!! These types are simplified analogues of the MOM6 types used in
!! MOM_lateral_mixing_coeffs.F90, specifically preserving the structural features
!! that are suspected to be relevant to each compiler bug:
!!
!!   1. Nested derived types (hor_index_type nested inside grid_type)
!!   2. Pointer members in grid_type and mixing_cs_type
!!   3. Allocatable 1D array members (vert_grid_type%g_prime)
!!   4. Allocatable 2D array members with non-unit lower bounds (staggered grids)
!!
!! Analogues of the original MOM6 types:
!!   grid_type       <- ocean_grid_type   (MOM_grid.F90)
!!   vert_grid_type  <- verticalGrid_type (MOM_verticalGrid.F90)
!!   mixing_cs_type  <- VarMix_CS        (MOM_lateral_mixing_coeffs.F90)
module mre_ocean_types_mod
  implicit none

  !> Horizontal index bounds nested inside grid_type.
  !! Analogue of MOM6's hor_index_type from MOM_hor_index.F90.
  type :: hor_index_type
    integer :: isc, iec   !< Start/end tracer i-indices, computational domain
    integer :: jsc, jec   !< Start/end tracer j-indices, computational domain
    integer :: isd, ied   !< Start/end tracer i-indices, data domain (includes halos)
    integer :: jsd, jed   !< Start/end tracer j-indices, data domain
    integer :: IscB, IecB !< Start/end staggered i-indices, computational domain
    integer :: JscB, JecB !< Start/end staggered j-indices, computational domain
    integer :: IsdB, IedB !< Start/end staggered i-indices, data domain
    integer :: JsdB, JedB !< Start/end staggered j-indices, data domain
  end type hor_index_type

  !> Simplified ocean grid type.
  !! Analogue of MOM6's ocean_grid_type from MOM_grid.F90. Preserves:
  !!   - Integer index bounds as scalar members
  !!   - A nested derived-type member (HI)
  !!   - A pointer member (domain_id, simulating the Domain pointer)
  !!   - Allocatable 2D array members at both tracer and staggered points
  type :: grid_type
    integer :: isc, iec   !< Computational domain tracer i-bounds
    integer :: jsc, jec   !< Computational domain tracer j-bounds
    integer :: isd, ied   !< Data domain tracer i-bounds (includes halos)
    integer :: jsd, jed   !< Data domain tracer j-bounds
    integer :: IscB, IecB !< Computational domain staggered (u-point) i-bounds
    integer :: JscB, JecB !< Computational domain staggered (v-point) j-bounds
    integer :: IsdB, IedB !< Data domain staggered i-bounds
    integer :: JsdB, JedB !< Data domain staggered j-bounds

    !> Nested horizontal index structure.
    !! Matches the HI member of ocean_grid_type. The nested type structure
    !! is relevant to the LLC compilation bug.
    type(hor_index_type) :: HI

    !> Pointer member simulating ocean_grid_type's Domain and US pointers.
    !! The presence of pointer members in the struct is relevant to how the
    !! compiler generates LLVM IR for GPU kernels that implicitly capture this type.
    integer, pointer :: domain_id => NULL()

    ! Allocatable 2D arrays at h-points (tracer/cell-center points)
    real, allocatable :: bathyT(:,:)    !< Ocean bottom depth [Z ~> m]
    real, allocatable :: meanSL(:,:)    !< Mean sea level [Z ~> m]

    ! Allocatable 2D arrays at u-points (staggered in i, IsdB:IedB x jsd:jed)
    real, allocatable :: mask2dCu(:,:)  !< Ocean mask at u-points [nondim]
    real, allocatable :: OBCmaskCu(:,:) !< OBC-aware mask at u-points [nondim]
    real, allocatable :: IdxCu(:,:)     !< 1/dx at u-points [L-1 ~> m-1]

    ! Allocatable 2D arrays at v-points (staggered in j, isd:ied x JsdB:JedB)
    real, allocatable :: mask2dCv(:,:)  !< Ocean mask at v-points [nondim]
    real, allocatable :: OBCmaskCv(:,:) !< OBC-aware mask at v-points [nondim]
    real, allocatable :: IdyCv(:,:)     !< 1/dy at v-points [L-1 ~> m-1]
  end type grid_type

  !> Simplified vertical grid type.
  !! Analogue of MOM6's verticalGrid_type from MOM_verticalGrid.F90.
  type :: vert_grid_type
    integer :: ke                    !< Number of vertical layers
    real :: H_subroundoff            !< Sub-roundoff thickness [H ~> m]
    real :: Angstrom_H               !< One-Angstrom thickness in H units [H ~> m]
    real :: Angstrom_Z               !< One-Angstrom thickness in Z units [Z ~> m]
    real :: dZ_subroundoff           !< Sub-roundoff depth [Z ~> m]
    real, allocatable :: g_prime(:)  !< Reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2]
  end type vert_grid_type

  !> Simplified variable mixing control structure.
  !! Analogue of MOM6's VarMix_CS from MOM_lateral_mixing_coeffs.F90. Preserves:
  !!   - Logical flag members
  !!   - Integer and real scalar members
  !!   - A pointer member (diag_id, simulating the diag_ctrl pointer)
  !!   - Allocatable 2D array members that are explicitly mapped to GPU via
  !!     !$omp target enter/exit data directives
  type :: mixing_cs_type
    logical :: initialized                  !< True if initialized
    logical :: calculate_Eady_growth_rate   !< If true, calculate Eady growth rate
    logical :: full_depth_Eady_growth_rate  !< If true, use full column depth as denominator
    logical :: OBC_friendly                 !< If true, handle open boundary condition points
    integer :: VarMix_Ktop                  !< Shallowest layer for downward depth integrals
    real :: h_min_N2                        !< Minimum thickness for N2 denominator [H ~> m]
    real :: Visbeck_S_max                   !< Upper bound on slope magnitude [Z L-1 ~> nondim]

    !> Pointer member simulating VarMix_CS's diag_ctrl pointer.
    integer, pointer :: diag_id => NULL()

    real, allocatable :: SN_u(:,:)  !< S*N at u-points [T-1 ~> s-1]
    real, allocatable :: SN_v(:,:)  !< S*N at v-points [T-1 ~> s-1]
  end type mixing_cs_type

end module mre_ocean_types_mod
