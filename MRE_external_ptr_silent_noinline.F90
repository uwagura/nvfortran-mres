!> MRE: Silent non-inlining when derived type contains a pointer to an external module type
!!
!! Root cause: When a subroutine argument's derived type contains a POINTER to a type
!! defined in a SEPARATELY COMPILED MODULE, nvfortran silently skips inlining that
!! subroutine even when -Minline=name:...,reshape is passed. No diagnostic message is
!! produced at the call site — not "inlined", not "subprogram not inlined", nothing.
!!
!! Demonstrated here with two variants in a single file:
!!   - Variant A (WORKS):   grid_no_ptr — derived type with no pointer to external type
!!                          → gradKE_no_ptr IS inlined; message appears in -Minfo=inline output
!!   - Variant B (SILENT):  grid_with_ptr — derived type with a pointer to external_type
!!                          → gradKE_with_ptr is NOT inlined; NO message in -Minfo=inline output
!!
!! The external_type is intentionally trivial (a few integers) to show that the issue is
!! structural (the pointer member's type origin) not the content of the pointed-to type.
!!
!! Compile command that reproduces the bug (requires nvfortran 26.3 or later):
!!   source <path>/envs/linux-nvidia_2026_3.env
!!   # Step 1: compile external module
!!   nvfortran -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -r8 \
!!     -c external_type_mod.F90 -o external_type_mod.o
!!   # Step 2: compile main MRE
!!   nvfortran -O2 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -r8 \
!!     -Minfo=inline -Minline=name:gradKE_no_ptr,reshape \
!!                   -Minline=name:gradKE_with_ptr,reshape \
!!     -I. external_type_mod.o MRE_external_ptr_silent_noinline.F90 -o mre_extptr.exe
!!
!! Expected output (from -Minfo=inline):
!!   gradke_no_ptr:      <- compiler processes gradke_no_ptr
!!   gradke_with_ptr:    <- compiler processes gradke_with_ptr
!!   caller_no_ptr:
!!     N, gradke_no_ptr inlined, size=..., file ...    <- Variant A: INLINED (message present)
!!   caller_with_ptr:
!!                                                     <- Variant B: SILENT (no message at all)
!!
!! The Makefile targets `mre_extptr` and `mre_extptr_baseline` compile this MRE.
!! This bug relates to gradKE in MOM6's MOM_CoriolisAdv.F90, where ocean_grid_type
!! contains `type(MOM_domain_type), pointer :: Domain => NULL()` from MOM_domains module.

! ============================================================
! VARIANT A: grid type with NO pointer to external module type
! → gradKE_no_ptr DOES inline
! ============================================================
module variant_a_mod
  implicit none

  !> Grid type containing only scalar integers and allocatable arrays (all local).
  !! No pointer to any externally-defined derived type.
  type :: grid_no_ptr
    integer :: isc, iec, jsc, jec  !< Computational domain bounds
    integer :: isd, ied, jsd, jed  !< Data domain bounds
    integer :: IscB, IecB, JscB, JecB  !< Staggered (B-grid) computational bounds
    integer :: IsdB, IedB, JsdB, JedB  !< Staggered (B-grid) data bounds
    real, allocatable :: areaCu(:,:)         !< U-cell areas [m2]
    real, allocatable :: areaCv(:,:)         !< V-cell areas [m2]
    real, allocatable :: IareaT(:,:)         !< Inverse T-cell areas [m-2]
    real, allocatable :: IdxCu_OBCmask(:,:)  !< 1/dx at u-points [m-1]
  end type grid_no_ptr

contains

  !> Compute kinetic energy and its gradient (simplified gradKE analog)
  subroutine gradKE_no_ptr(u, v, KE, KEx, G)
    type(grid_no_ptr), intent(in) :: G
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed),  intent(in)  :: u
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB), intent(in)  :: v
    real, dimension(G%isd:G%ied,   G%jsd:G%jed),   intent(out) :: KE
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed),   intent(out) :: KEx
    integer :: i, j, Isq, Ieq, Jsq, Jeq
    Isq=G%IscB; Ieq=G%IecB; Jsq=G%JscB; Jeq=G%JecB
    do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
      KE(i,j) = (G%areaCu(i,j)*u(i,j)**2 + G%areaCv(i,j)*v(i,j)**2)*0.25*G%IareaT(i,j)
    end do
    do concurrent (j=G%jsc:G%jec, i=Isq:Ieq)
      KEx(i,j) = (KE(i+1,j)-KE(i,j))*G%IdxCu_OBCmask(i,j)
    end do
  end subroutine gradKE_no_ptr

  !> Outer loop that calls gradKE_no_ptr (caller analog)
  subroutine caller_no_ptr(u3, v3, KEx_out, G, nz)
    type(grid_no_ptr), intent(in) :: G
    integer,           intent(in) :: nz
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   nz), intent(in)  :: u3
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB, nz), intent(in)  :: v3
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   nz), intent(out) :: KEx_out
    real, dimension(G%isd:G%ied,   G%jsd:G%jed)   :: KE
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed)   :: KEx
    integer :: k, i, j
    !$omp target enter data map(alloc: KE, KEx)
    do k = 1, nz
      call gradKE_no_ptr(u3(:,:,k), v3(:,:,k), KE, KEx, G) ! <-- INLINED
      do concurrent (j=G%jsc:G%jec, i=G%IscB:G%IecB)
        KEx_out(i,j,k) = KEx(i,j)
      end do
    end do
    !$omp target exit data map(delete: KE, KEx)
  end subroutine caller_no_ptr

end module variant_a_mod


! ============================================================
! VARIANT B: grid type WITH a pointer to external module type
! → gradKE_with_ptr does NOT inline; compiler is SILENT
! ============================================================
module variant_b_mod
  use external_type_mod, only : external_type  ! Type defined in separately compiled module
  implicit none

  !> Grid type identical to grid_no_ptr except it also contains a pointer to an
  !! externally-defined derived type. This single addition silences the inliner.
  type :: grid_with_ptr
    type(external_type), pointer :: Domain => NULL()  ! <-- the trigger
    integer :: isc, iec, jsc, jec  !< Computational domain bounds
    integer :: isd, ied, jsd, jed  !< Data domain bounds
    integer :: IscB, IecB, JscB, JecB  !< Staggered (B-grid) computational bounds
    integer :: IsdB, IedB, JsdB, JedB  !< Staggered (B-grid) data bounds
    real, allocatable :: areaCu(:,:)         !< U-cell areas [m2]
    real, allocatable :: areaCv(:,:)         !< V-cell areas [m2]
    real, allocatable :: IareaT(:,:)         !< Inverse T-cell areas [m-2]
    real, allocatable :: IdxCu_OBCmask(:,:)  !< 1/dx at u-points [m-1]
  end type grid_with_ptr

contains

  !> Compute kinetic energy and its gradient (same logic as gradKE_no_ptr)
  subroutine gradKE_with_ptr(u, v, KE, KEx, G)
    type(grid_with_ptr), intent(in) :: G
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed),  intent(in)  :: u
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB), intent(in)  :: v
    real, dimension(G%isd:G%ied,   G%jsd:G%jed),   intent(out) :: KE
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed),   intent(out) :: KEx
    integer :: i, j, Isq, Ieq, Jsq, Jeq
    Isq=G%IscB; Ieq=G%IecB; Jsq=G%JscB; Jeq=G%JecB
    do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
      KE(i,j) = (G%areaCu(i,j)*u(i,j)**2 + G%areaCv(i,j)*v(i,j)**2)*0.25*G%IareaT(i,j)
    end do
    do concurrent (j=G%jsc:G%jec, i=Isq:Ieq)
      KEx(i,j) = (KE(i+1,j)-KE(i,j))*G%IdxCu_OBCmask(i,j)
    end do
  end subroutine gradKE_with_ptr

  !> Outer loop that calls gradKE_with_ptr (caller analog)
  subroutine caller_with_ptr(u3, v3, KEx_out, G, nz)
    type(grid_with_ptr), intent(in) :: G
    integer,             intent(in) :: nz
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   nz), intent(in)  :: u3
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB, nz), intent(in)  :: v3
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   nz), intent(out) :: KEx_out
    real, dimension(G%isd:G%ied,   G%jsd:G%jed)   :: KE
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed)   :: KEx
    integer :: k, i, j
    !$omp target enter data map(alloc: KE, KEx)
    do k = 1, nz
      call gradKE_with_ptr(u3(:,:,k), v3(:,:,k), KE, KEx, G) ! <-- SILENT: NOT inlined
      do concurrent (j=G%jsc:G%jec, i=G%IscB:G%IecB)
        KEx_out(i,j,k) = KEx(i,j)
      end do
    end do
    !$omp target exit data map(delete: KE, KEx)
  end subroutine caller_with_ptr

end module variant_b_mod


! ============================================================
! Driver program
! ============================================================
program test_external_ptr_inline
  use variant_a_mod
  use variant_b_mod
  implicit none

  integer, parameter :: nx=20, ny=20, nz=3, halo=2

  type(grid_no_ptr)   :: G_a
  type(grid_with_ptr) :: G_b

  real, allocatable :: u3(:,:,:), v3(:,:,:), KEx_out(:,:,:)

  ! Variant A setup
  G_a%isc=1+halo; G_a%iec=nx+halo; G_a%jsc=1+halo; G_a%jec=ny+halo
  G_a%isd=1; G_a%ied=nx+2*halo; G_a%jsd=1; G_a%jed=ny+2*halo
  G_a%IscB=G_a%isc-1; G_a%IecB=G_a%iec; G_a%JscB=G_a%jsc-1; G_a%JecB=G_a%jec
  G_a%IsdB=G_a%isd-1; G_a%IedB=G_a%ied; G_a%JsdB=G_a%jsd-1; G_a%JedB=G_a%jed
  allocate(G_a%areaCu(G_a%IsdB:G_a%IedB, G_a%jsd:G_a%jed))
  allocate(G_a%areaCv(G_a%isd:G_a%ied,   G_a%JsdB:G_a%JedB))
  allocate(G_a%IareaT(G_a%isd:G_a%ied,   G_a%jsd:G_a%jed))
  allocate(G_a%IdxCu_OBCmask(G_a%IsdB:G_a%IedB, G_a%jsd:G_a%jed))
  G_a%areaCu=1.0; G_a%areaCv=1.0; G_a%IareaT=0.25; G_a%IdxCu_OBCmask=1.0

  ! Variant B setup (same bounds as A)
  G_b%isc=G_a%isc; G_b%iec=G_a%iec; G_b%jsc=G_a%jsc; G_b%jec=G_a%jec
  G_b%isd=G_a%isd; G_b%ied=G_a%ied; G_b%jsd=G_a%jsd; G_b%jed=G_a%jed
  G_b%IscB=G_a%IscB; G_b%IecB=G_a%IecB; G_b%JscB=G_a%JscB; G_b%JecB=G_a%JecB
  G_b%IsdB=G_a%IsdB; G_b%IedB=G_a%IedB; G_b%JsdB=G_a%JsdB; G_b%JedB=G_a%JedB
  allocate(G_b%areaCu(G_b%IsdB:G_b%IedB, G_b%jsd:G_b%jed))
  allocate(G_b%areaCv(G_b%isd:G_b%ied,   G_b%JsdB:G_b%JedB))
  allocate(G_b%IareaT(G_b%isd:G_b%ied,   G_b%jsd:G_b%jed))
  allocate(G_b%IdxCu_OBCmask(G_b%IsdB:G_b%IedB, G_b%jsd:G_b%jed))
  G_b%areaCu=1.0; G_b%areaCv=1.0; G_b%IareaT=0.25; G_b%IdxCu_OBCmask=1.0

  allocate(u3(G_a%IsdB:G_a%IedB, G_a%jsd:G_a%jed, nz))
  allocate(v3(G_a%isd:G_a%ied,   G_a%JsdB:G_a%JedB, nz))
  allocate(KEx_out(G_a%IsdB:G_a%IedB, G_a%jsd:G_a%jed, nz))
  call random_number(u3); call random_number(v3)

  print *, "=== Variant A (grid_no_ptr: no external pointer) ==="
  call caller_no_ptr(u3, v3, KEx_out, G_a, nz)
  print *, "  KEx sum:", sum(KEx_out)

  print *, "=== Variant B (grid_with_ptr: pointer to external module type) ==="
  call caller_with_ptr(u3, v3, KEx_out, G_b, nz)
  print *, "  KEx sum:", sum(KEx_out)

end program test_external_ptr_inline
