!> Minimal Reproducible Example (MRE) for nvfortran 26.3 LLC code-generation error:
!!   LLVM ERROR: use of undefined value '%g'
!!
!! This MRE reproduces a compile-time crash (in the LLVM backend / llc) when
!! compiling calc_slope_functions_using_just_e from MOM_lateral_mixing_coeffs.F90
!! with:
!!   nvfortran -O0 -mp=gpu -stdpar=gpu -gpu=mem:separate,cc80,sm_80 -Mnofma -Minfo=all -r8
!!
!! === The error ===
!!   /tmp/nvfortranNNNNNN-N.ll:NNN:NN: error: use of undefined value '%g'
!!         %NNN = getelementptr i8, ptr %g, i64 NNNN, !dbg !NNNNN
!! (followed by a similar line for each allocatable member of G accessed in the kernel)
!!
!! === Root cause (hypothesis) ===
!! The LLC error is triggered by the combination of:
!!   1. The derived type G (grid_type) is NOT explicitly mapped via
!!      !$omp target enter data. G has a nested derived type (HI) and a
!!      pointer member (domain_id), matching the structure of MOM6's ocean_grid_type.
!!   2. G's allocatable member arrays (G%bathyT, G%meanSL, G%OBCmaskCu) ARE
!!      accessed inside an outer do concurrent(j=js:je) loop that is offloaded
!!      to the GPU via -stdpar=gpu. The compiler therefore generates implicit
!!      copyin/copy directives for individual struct members.
!!   3. Separately, G%IdxCu and G%IdyCv are accessed in an !$OMP parallel do
!!      region with inner do concurrent loops, which also causes the compiler to
!!      generate "Generating implicit copy(g) [if not already present]" diagnostics.
!!   4. The conflict between implicitly copying the whole struct G (from step 3)
!!      and individually copying its allocatable members (from steps 2 and 3)
!!      causes the LLC to fail to define %g in the generated GPU kernel function,
!!      producing the "use of undefined value '%g'" error.
!!   5. A runtime logical (use_dztot = CS%full_depth_Eady_growth_rate) controls
!!      which branch of an if/else is taken, but the compiler statically emits
!!      code for BOTH branches and therefore maps G%bathyT, G%meanSL, AND
!!      G%OBCmaskCu regardless of the runtime value of use_dztot.
!!
!! === To compile (expected to FAIL with LLC error) ===
!!   cd /scratch/cimes/uw4770/nvfortran-mres
!!   make mre1
!!
!! === Workaround (A/B) ===
!!   Explicitly mapping G before the do concurrent avoids the LLC error:
!!     !$omp target enter data map(to: G)   ! <-- uncomment in calc_slope_mre
!!   (Tested: compilation succeeds with this addition.)
!!
!! === Original code ===
!!   MOM6 subroutine calc_slope_functions_using_just_e
!!   File: MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
!!   Compiler: nvfortran 26.3 (NVHPC 2026.3)
program mre1_llc_undefined_value
  use mre_ocean_types_mod
  implicit none

  ! Grid dimensions (small domain for MRE; halos are 2 cells wide)
  integer, parameter :: NI = 20   ! Data domain i-size (1:NI)
  integer, parameter :: NJ = 20   ! Data domain j-size (1:NJ)
  integer, parameter :: NK = 10   ! Number of vertical layers

  type(grid_type)     :: G
  type(vert_grid_type) :: GV
  type(mixing_cs_type) :: CS

  real, allocatable :: h(:,:,:)  ! Layer thickness [H ~> m],          h(isd:ied, jsd:jed, ke)
  real, allocatable :: e(:,:,:)  ! Interface position (positive up) [Z ~> m], e(isd:ied, jsd:jed, ke+1)

  integer :: i, j, k

  ! --- Set up index bounds (symmetric memory layout, 2-cell halos) ---
  G%isd = 1  ; G%ied = NI    ; G%jsd = 1  ; G%jed = NJ
  G%isc = 3  ; G%iec = NI-2  ; G%jsc = 3  ; G%jec = NJ-2

  ! u-staggered: IsdB = isd-1 (extra column on low side in symmetric memory)
  G%IsdB = 0  ; G%IedB = NI   ; G%JsdB = 0  ; G%JedB = NJ
  G%IscB = G%isc-1 ; G%IecB = G%iec
  G%JscB = G%jsc-1 ; G%JecB = G%jec

  ! Copy bounds into nested HI structure
  G%HI%isc = G%isc ; G%HI%iec = G%iec ; G%HI%jsc = G%jsc ; G%HI%jec = G%jec
  G%HI%isd = G%isd ; G%HI%ied = G%ied ; G%HI%jsd = G%jsd ; G%HI%jed = G%jed
  G%HI%IscB = G%IscB ; G%HI%IecB = G%IecB
  G%HI%JscB = G%JscB ; G%HI%JecB = G%JecB
  G%HI%IsdB = G%IsdB ; G%HI%IedB = G%IedB
  G%HI%JsdB = G%JsdB ; G%HI%JedB = G%JedB

  ! Allocate and initialize grid arrays
  allocate(G%IdxCu(G%IsdB:G%IedB, G%jsd:G%jed))
  allocate(G%IdyCv(G%isd:G%ied,   G%JsdB:G%JedB))
  allocate(G%OBCmaskCu(G%IsdB:G%IedB, G%jsd:G%jed))
  allocate(G%OBCmaskCv(G%isd:G%ied,   G%JsdB:G%JedB))
  allocate(G%mask2dCu(G%IsdB:G%IedB, G%jsd:G%jed))
  allocate(G%mask2dCv(G%isd:G%ied,   G%JsdB:G%JedB))
  allocate(G%bathyT(G%isd:G%ied, G%jsd:G%jed))
  allocate(G%meanSL(G%isd:G%ied, G%jsd:G%jed))
  G%IdxCu(:,:)    = 1.0e-5   ! ~100 km spacing: 1/100000 m-1
  G%IdyCv(:,:)    = 1.0e-5
  G%OBCmaskCu(:,:) = 1.0 ; G%OBCmaskCv(:,:) = 1.0
  G%mask2dCu(:,:)  = 1.0 ; G%mask2dCv(:,:)  = 1.0
  G%bathyT(:,:)    = 1000.0  ! 1000 m deep
  G%meanSL(:,:)    = 0.0

  ! --- Set up vertical grid ---
  GV%ke             = NK
  GV%H_subroundoff  = 1.0e-17
  GV%Angstrom_H     = 1.0e-10
  GV%Angstrom_Z     = 1.0e-10
  GV%dZ_subroundoff = 1.0e-17
  allocate(GV%g_prime(NK+1))
  GV%g_prime(:) = 0.02   ! reduced gravity ~ 0.02 m s-2

  ! --- Set up mixing control structure ---
  CS%initialized                = .true.
  CS%calculate_Eady_growth_rate = .true.
  CS%full_depth_Eady_growth_rate = .true.   ! controls the if/else branch in the kernel
  CS%OBC_friendly               = .false.
  CS%VarMix_Ktop                = 2
  CS%h_min_N2                   = 1.0e-3
  CS%Visbeck_S_max              = 0.01
  allocate(CS%SN_u(G%IsdB:G%IedB, G%jsd:G%jed))
  allocate(CS%SN_v(G%isd:G%ied,   G%JsdB:G%JedB))
  CS%SN_u(:,:) = 0.0 ; CS%SN_v(:,:) = 0.0

  ! --- Allocate layer and interface arrays ---
  allocate(h(G%isd:G%ied, G%jsd:G%jed, NK))
  allocate(e(G%isd:G%ied, G%jsd:G%jed, NK+1))
  do k = 1, NK
    h(:,:,k) = 100.0   ! uniform 100 m layers
  enddo
  do k = 1, NK+1
    e(:,:,k) = real(NK+1-k) * 100.0  ! interfaces at 1000, 900, ..., 0 m
  enddo

  ! Call the MRE subroutine — expected to trigger the LLC error at *compile* time.
  call calc_slope_mre(h, G, GV, CS, e)

  print *, "MRE1: If this line is reached, the LLC bug was not triggered."
  print *, "SN_u(isc,jsc) =", CS%SN_u(G%isc, G%jsc)

contains

  !> Reproduces the structure of calc_slope_functions_using_just_e.
  !!
  !! Key structural features preserved:
  !!   1. G (grid_type) is NOT mapped via !$omp target enter data; only local
  !!      arrays and CS%SN_u/SN_v are explicitly mapped.
  !!   2. An !$OMP parallel do region over k with inner do concurrent loops
  !!      accesses G%IdxCu and G%IdyCv (allocatable members), causing the
  !!      compiler to emit "Generating implicit copy(g)" diagnostics.
  !!   3. An outer do concurrent(j=js:je) — offloaded via -stdpar=gpu — contains:
  !!        a. A nested sequential do k loop with an inner do concurrent
  !!        b. A runtime if/else block (use_dztot) whose two arms access
  !!           G%OBCmaskCu, G%bathyT, G%meanSL, and GV%dZ_subroundoff
  !!   4. G has a nested derived type member (HI) and a pointer member (domain_id)
  !!
  !! To work around the LLC error, uncomment the explicit map(to: G) line below.
  subroutine calc_slope_mre(h, G, GV, CS, e)
    type(grid_type),      intent(inout) :: G   !< Horizontal grid; intentionally NOT mapped
    type(vert_grid_type), intent(in)    :: GV  !< Vertical grid
    real, dimension(G%isd:G%ied, G%jsd:G%jed, GV%ke),   intent(inout) :: h  !< Layer thickness [H]
    type(mixing_cs_type), intent(inout) :: CS  !< Mixing control structure
    real, dimension(G%isd:G%ied, G%jsd:G%jed, GV%ke+1), intent(in)    :: e  !< Interface position [Z]

    ! Local arrays with explicit non-unit lower bounds matching the original subroutine.
    ! These are mapped to the device via !$omp target enter data below.
    real :: E_x(G%IsdB:G%IedB, G%jsd:G%jed)                      ! i-slope at u-points [Z L-1]
    real :: E_y(G%isd:G%ied,   G%JsdB:G%JedB)                    ! j-slope at v-points [Z L-1]
    real :: dz_tot(G%isd:G%ied, G%jsd:G%jed)                     ! Full column thickness [Z]
    real :: S2N2_u_local(G%IsdB:G%IedB, G%jsd:G%jed, GV%ke)      ! S^2*N^2*H at u-pts [T-2 H]
    real :: S2N2_v_local(G%isd:G%ied,   G%JsdB:G%JedB, GV%ke)    ! S^2*N^2*H at v-pts [T-2 H]

    ! Scalar locals
    real :: H_cutoff   ! Minimum layer thickness for masking [H]
    real :: dZ_cutoff  ! Minimum water column depth for masking [Z]
    real :: h_neglect  ! Negligibly small thickness [H]
    real :: S2         ! Squared interface slope [Z2 L-2]
    real :: Hup, Hdn   ! Layer thicknesses above and below interface [H]
    real :: H_geom     ! Geometric mean of Hup and Hdn [H]
    real :: h1, h2     ! Temporary full-column depths [Z]
    logical :: use_dztot

    integer :: is, ie, js, je, nz, i, j, k

    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
    h_neglect = GV%H_subroundoff
    H_cutoff  = real(2*nz) * (GV%Angstrom_H + h_neglect)
    dZ_cutoff = real(2*nz) * (GV%Angstrom_Z + GV%dZ_subroundoff)

    ! Runtime logical: controls which branch of the if/else in the outer
    ! do concurrent is taken. The compiler maps data for BOTH branches
    ! statically regardless of this value, which is part of the bug trigger.
    use_dztot = CS%full_depth_Eady_growth_rate

    ! Map local output arrays and CS%SN_u/v explicitly.
    ! G is intentionally NOT mapped — this is the key condition for the LLC bug.
    ! A/B workaround: uncommenting the line below avoids the LLC error.
    ! !$omp target enter data map(to: G)
    !$omp target enter data map(alloc: E_x, E_y, S2N2_u_local, S2N2_v_local)
    !$omp target enter data map(to: CS%SN_u, CS%SN_v)
    !$omp target enter data map(alloc: dz_tot) if (use_dztot)

    if (use_dztot) then
      do concurrent( j=js-1:je+1, i=is-1:ie+1 )
        dz_tot(i,j) = e(i,j,1) - e(i,j,nz+1)
      enddo
    endif

    ! !$OMP parallel do region with inner do concurrent loops.
    ! The inner loops access G%IdxCu and G%IdyCv (allocatable members of G),
    ! causing the compiler to emit:
    !   "Generating implicit copyin(g%idxcu(:,:),g) [if not already present]"
    ! This is one of the two contexts that reference G without an explicit map.
    !$OMP parallel do default(shared) private(E_x, E_y, S2, Hdn, Hup, H_geom)
    do k = nz, CS%VarMix_Ktop, -1

      ! Interface slope at u-points: accesses G%IdxCu (allocatable member of G)
      do concurrent( j=js-1:je+1, i=is-1:ie )
        E_x(i,j) = (e(i+1,j,k) - e(i,j,k)) * G%IdxCu(i,j)
        if (min(h(i,j,k), h(i+1,j,k)) < H_cutoff) E_x(i,j) = 0.0
      enddo

      ! Interface slope at v-points: accesses G%IdyCv (allocatable member of G)
      do concurrent( j=js-1:je, i=is-1:ie+1 )
        E_y(i,j) = (e(i,j+1,k) - e(i,j,k)) * G%IdyCv(i,j)
        if (min(h(i,j,k), h(i,j+1,k)) < H_cutoff) E_y(i,j) = 0.0
      enddo

      ! S^2*N^2*H at u-points. LOCAL clause privatizes the scalar temporaries.
      ! Accesses GV%g_prime(k) (allocatable 1D member of GV) and CS%h_min_N2.
      do concurrent( j=js:je, i=is-1:ie ) local( S2, Hdn, Hup, H_geom )
        S2 = E_x(i,j)**2 + 0.25*( &
              ( E_y(i,  j  )**2 + E_y(i+1,j-1)**2 ) + &
              ( E_y(i+1,j  )**2 + E_y(i,  j-1)**2 ) )
        if (min(h(i,j,k-1), h(i+1,j,k-1), h(i,j,k), h(i+1,j,k)) < H_cutoff) S2 = 0.0
        Hdn    = 2.0*h(i,  j,k)*h(i,  j,k-1) / (h(i,  j,k) + h(i,  j,k-1) + h_neglect)
        Hup    = 2.0*h(i+1,j,k)*h(i+1,j,k-1) / (h(i+1,j,k) + h(i+1,j,k-1) + h_neglect)
        H_geom = sqrt(Hdn * Hup)
        S2N2_u_local(i,j,k) = (H_geom * S2) * (GV%g_prime(k) / max(Hdn, Hup, CS%h_min_N2))
      enddo

      ! S^2*N^2*H at v-points. Same structure.
      do concurrent( j=js-1:je, i=is:ie ) local( S2, Hdn, Hup, H_geom )
        S2 = E_y(i,j)**2 + 0.25*( &
              ( E_x(i,  j  )**2 + E_x(i-1,j+1)**2 ) + &
              ( E_x(i,  j+1)**2 + E_x(i-1,j  )**2 ) )
        if (min(h(i,j,k-1), h(i,j+1,k-1), h(i,j,k), h(i,j+1,k)) < H_cutoff) S2 = 0.0
        Hdn    = 2.0*h(i,j,  k)*h(i,j,  k-1) / (h(i,j,  k) + h(i,j,  k-1) + h_neglect)
        Hup    = 2.0*h(i,j+1,k)*h(i,j+1,k-1) / (h(i,j+1,k) + h(i,j+1,k-1) + h_neglect)
        H_geom = sqrt(Hdn * Hup)
        S2N2_v_local(i,j,k) = (H_geom * S2) * (GV%g_prime(k) / max(Hdn, Hup, CS%h_min_N2))
      enddo

    enddo  ! k

    ! Outer do concurrent(j=js:je) — offloaded to GPU via -stdpar=gpu.
    ! The inner do concurrent(i) loops access CS%SN_u (explicitly mapped) and
    ! S2N2_u_local (explicitly mapped), but the subsequent if/else block accesses
    ! G%OBCmaskCu, G%bathyT, G%meanSL, and GV%dZ_subroundoff. Since G is NOT
    ! explicitly mapped, the compiler generates implicit copy directives for
    ! individual allocatable members of G. Combined with the implicit copy(g)
    ! generated by the !$OMP parallel do above, this causes the LLC error:
    !   "use of undefined value '%g'"
    do concurrent( j=js:je )

      do concurrent( i=is-1:ie )
        CS%SN_u(i,j) = 0.0
      enddo

      do k=nz, CS%VarMix_Ktop, -1 ; do concurrent( i=is-1:ie )
        CS%SN_u(i,j) = CS%SN_u(i,j) + S2N2_u_local(i,j,k)
      enddo ; enddo

      ! Runtime branch: both arms access G%... members. The compiler statically
      ! maps data for both arms regardless of use_dztot's runtime value.
      if (use_dztot) then
        ! Accesses G%OBCmaskCu (allocatable member) and GV%dZ_subroundoff (scalar)
        do concurrent( i=is-1:ie )
          CS%SN_u(i,j) = G%OBCmaskCu(i,j) * sqrt( CS%SN_u(i,j) / &
                          max(dz_tot(i,j), dz_tot(i+1,j), GV%dZ_subroundoff) )
        enddo
      else
        ! Accesses G%bathyT, G%meanSL, G%OBCmaskCu (all allocatable members of G)
        do concurrent( i=is-1:ie ) local( h1, h2 )
          h1 = max(G%meanSL(i,  j) + G%bathyT(i,  j), 0.0)
          h2 = max(G%meanSL(i+1,j) + G%bathyT(i+1,j), 0.0)
          if ( min(h1, h2) > dZ_cutoff ) then
            CS%SN_u(i,j) = G%OBCmaskCu(i,j) * sqrt( CS%SN_u(i,j) / max(h1, h2) )
          else
            CS%SN_u(i,j) = 0.0
          endif
        enddo
      endif

    enddo  ! outer do concurrent(j=js:je) for u

    ! Analogous outer do concurrent(j=js-1:je) for v-points
    do concurrent( j=js-1:je )

      do concurrent( i=is:ie )
        CS%SN_v(i,j) = 0.0
      enddo

      do k=nz, CS%VarMix_Ktop, -1 ; do concurrent( i=is:ie )
        CS%SN_v(i,j) = CS%SN_v(i,j) + S2N2_v_local(i,j,k)
      enddo ; enddo

      if (use_dztot) then
        do concurrent( i=is:ie )
          CS%SN_v(i,j) = G%OBCmaskCv(i,j) * sqrt( CS%SN_v(i,j) / &
                          max(dz_tot(i,j), dz_tot(i,j+1), GV%dZ_subroundoff) )
        enddo
      else
        do concurrent( i=is:ie ) local( h1, h2 )
          h1 = max(G%meanSL(i,j  ) + G%bathyT(i,j  ), 0.0)
          h2 = max(G%meanSL(i,j+1) + G%bathyT(i,j+1), 0.0)
          if ( min(h1, h2) > dZ_cutoff ) then
            CS%SN_v(i,j) = G%OBCmaskCv(i,j) * sqrt( CS%SN_v(i,j) / max(h1, h2) )
          else
            CS%SN_v(i,j) = 0.0
          endif
        enddo
      endif

    enddo  ! outer do concurrent(j=js-1:je) for v

    !$omp target exit data map(delete: E_x, E_y, S2N2_u_local, S2N2_v_local)
    !$omp target exit data map(from: CS%SN_u, CS%SN_v)
    if (use_dztot) then
      !$omp target exit data map(delete: dz_tot)
    endif

  end subroutine calc_slope_mre

end program mre1_llc_undefined_value
