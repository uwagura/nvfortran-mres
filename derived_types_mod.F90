module derived_types_mod
  implicit none
#ifdef USE_DERIVED_TYPES
  type :: my_type1
    real, allocatable :: arr(:)
    integer :: n
    real :: scalar
  end type my_type1

  type :: my_type2
    real, allocatable :: matrix(:,:)
    integer :: m, n
    real :: value
  end type my_type2
#endif
end module derived_types_mod
