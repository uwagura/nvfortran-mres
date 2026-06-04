!> External module type used by MRE_external_ptr_silent_noinline.F90
!!
!! This is intentionally trivial — a handful of integers — to demonstrate that the
!! content of the pointed-to type is irrelevant. The mere presence of a POINTER to
!! any type from a separately compiled module in a derived type argument is sufficient
!! to cause nvfortran's inliner to silently skip inlining.
!!
!! In the real MOM6 case, the analog is:
!!   type(MOM_domain_type), pointer :: Domain => NULL()
!! inside ocean_grid_type (MOM_grid.F90), where MOM_domain_type is from MOM_domains.
module external_type_mod
  implicit none

  !> A trivial external type. Its content does not matter — any type from a
  !! separately compiled module triggers the silent inlining failure.
  type, public :: external_type
    integer :: id       !< An identifier
    integer :: nprocs   !< Number of processors
    logical :: active   !< Whether this domain is active
  end type external_type

end module external_type_mod
