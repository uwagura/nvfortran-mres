!> Minimal Reproducible Example (MRE) for nvfortran 26.3 runtime segfault when using
!! the LOCAL clause with a 1-D array in an outer do concurrent that is offloaded
!! to the GPU via -stdpar=gpu.
!!
!! === The error ===
!! Segmentation fault (signal 11) at runtime when executing the GPU kernel generated
!! for the outer do concurrent(j=js:je) local(H_u) loop.
!! The program runs silently when local(H_u) is removed (with incorrect results,
!! since H_u is shared across j-iterations without privatization).
!! The compiler also emits this diagnostic (not an error, just a warning):
!!   "Parallelization would require privatization of array h_u(:)"
!!
!! === Root cause (hypothesis) ===
!! The segfault is triggered by the combination of:
!!   1. H_u is a 1-D local array declared in the subroutine (not a module variable),
!!      with bounds G%IsdB:G%IedB (e.g., 0:20 for a 20-element domain with halos).
!!   2. H_u is explicitly mapped to the GPU device via:
!!        !$omp target enter data map(alloc: ..., H_u, H_v)
!!      This places H_u in GPU global memory.
!!   3. Inside an outer do concurrent(j=js:je) — offloaded via -stdpar=gpu —
!!      H_u is declared as LOCAL:
!!        do concurrent(j=js:je) local(H_u)
!!      The LOCAL clause is semantically correct: each j-iteration needs its own
!!      independent copy of H_u to accumulate column-average quantities over k.
!!   4. The combination of H_u being in GPU global memory (from the target map)
!!      AND being privatized per j-iteration (from the LOCAL clause) appears to
!!      cause a memory access violation at runtime. The compiler may be generating
!!      incorrect GPU code for the LOCAL privatization when the array is already
!!      present on the device.
!!
!! === Compiler diagnostics (from -Minfo=all) ===
!!   calc_visbeck_mre:
!!      NNN, Parallelization would require privatization of array h_u(:)
!!      NNN, Parallelization would require privatization of array h_v(:)
!!
!! === To compile (should succeed) ===
!!   cd /scratch/cimes/uw4770/nvfortran-mres
!!   make mre2
!!
!! === To run (requires GPU node, expected to SEGFAULT) ===
!!   ./mre2.exe
!!
!! === Workaround (A/B) ===
!! Removing local(H_u) from the do concurrent avoids the segfault:
!!   do concurrent(j=js:je)            ! <-- no local(H_u); wrong results but no crash
!! Alternatively, NOT mapping H_u via target enter data and using only local(H_u)
!! may also avoid the crash, but that changes the data management model.
!!
!! === Original code ===
!!   MOM6 subroutine calc_Visbeck_coeffs_old
!!   File: MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
!!   Compiler: nvfortran 26.3 (NVHPC 2026.3)
program mre2_do_concurrent_local
  use mre_ocean_types_mod
  implicit none

  ! Grid dimensions (small domain; halos are 2 cells wide)
  integer, parameter :: NI = 20   ! Data domain i-size (1:NI)
  integer, parameter :: NJ = 20   ! Data domain j-size (1:NJ)
  integer, parameter :: NK = 10   ! Number of vertical layers

  type(grid_type)     :: G
  type(vert_grid_type) :: GV
  type(mixing_cs_type) :: CS

  real, allocatable :: h(:,:,:)        ! Layer thickness [H ~> m]
  real, allocatable :: slope_x(:,:,:)  ! i-direction interface slope [Z L-1] at u-pts (IsdB:IedB, jsd:jed, ke+1)
  real, allocatable :: slope_y(:,:,:)  ! j-direction interface slope [Z L-1] at v-pts (isd:ied, JsdB:JedB, ke+1)
  real, allocatable :: N2_u(:,:,:)     ! Buoyancy frequency squared at u-pts [T-2] (IsdB:IedB, jsd:jed, ke+1)
  real, allocatable :: N2_v(:,:,:)     ! Buoyancy frequency squared at v-pts [T-2] (isd:ied, JsdB:JedB, ke+1)

  integer :: i, j, k, n

  ! --- Set up index bounds (symmetric memory layout, 2-cell halos) ---
  G%isd = 1  ; G%ied = NI    ; G%jsd = 1  ; G%jed = NJ
  G%isc = 3  ; G%iec = NI-2  ; G%jsc = 3  ; G%jec = NJ-2

  G%IsdB = 0  ; G%IedB = NI   ; G%JsdB = 0  ; G%JedB = NJ
  G%IscB = G%isc-1 ; G%IecB = G%iec
  G%JscB = G%jsc-1 ; G%JecB = G%jec

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
  G%IdxCu(:,:)     = 1.0e-5
  G%IdyCv(:,:)     = 1.0e-5
  G%OBCmaskCu(:,:) = 1.0 ; G%OBCmaskCv(:,:) = 1.0
  G%mask2dCu(:,:)  = 1.0 ; G%mask2dCv(:,:)  = 1.0
  G%bathyT(:,:)    = 1000.0
  G%meanSL(:,:)    = 0.0

  ! --- Set up vertical grid ---
  GV%ke             = NK
  GV%H_subroundoff  = 1.0e-17
  GV%Angstrom_H     = 1.0e-10
  GV%Angstrom_Z     = 1.0e-10
  GV%dZ_subroundoff = 1.0e-17
  allocate(GV%g_prime(NK+1))
  GV%g_prime(:) = 0.02

  ! --- Set up mixing control structure ---
  CS%initialized                = .true.
  CS%calculate_Eady_growth_rate = .true.
  CS%full_depth_Eady_growth_rate = .false.
  CS%OBC_friendly               = .false.
  CS%VarMix_Ktop                = 2
  CS%h_min_N2                   = 1.0e-3
  CS%Visbeck_S_max              = 0.01
  allocate(CS%SN_u(G%IsdB:G%IedB, G%jsd:G%jed))
  allocate(CS%SN_v(G%isd:G%ied,   G%JsdB:G%JedB))
  CS%SN_u(:,:) = 0.0 ; CS%SN_v(:,:) = 0.0

  ! --- Allocate input arrays ---
  allocate(h(G%isd:G%ied, G%jsd:G%jed, NK))
  allocate(slope_x(G%IsdB:G%IedB, G%jsd:G%jed, NK+1))
  allocate(slope_y(G%isd:G%ied,   G%JsdB:G%JedB, NK+1))
  allocate(N2_u(G%IsdB:G%IedB, G%jsd:G%jed, NK+1))
  allocate(N2_v(G%isd:G%ied,   G%JsdB:G%JedB, NK+1))
  do k = 1, NK
    h(:,:,k) = 100.0
  enddo
  do k = 1, NK+1
    slope_x(:,:,k) = 1.0e-4 * sin(real(k))   ! small slopes
    slope_y(:,:,k) = 1.0e-4 * cos(real(k))
    N2_u(:,:,k)    = 1.0e-4                   ! N^2 ~ 10^-4 s^-2 (stable)
    N2_v(:,:,k)    = 1.0e-4
  enddo

  call calc_visbeck_mre(h, slope_x, slope_y, N2_u, N2_v, G, GV, CS)

  print *, "MRE2: If this line is reached, the segfault was not triggered."
  print *, "SN_u(isc,jsc) =", CS%SN_u(G%isc, G%jsc)

contains

  !> Reproduces the structure of calc_Visbeck_coeffs_old.
  !!
  !! Key structural features preserved:
  !!   1. Local 1-D arrays H_u(G%IsdB:G%IedB) and H_v(G%isd:G%ied) are declared
  !!      as subroutine-local arrays (not module-level allocatables).
  !!   2. H_u and H_v are explicitly mapped to GPU device memory via:
  !!        !$omp target enter data map(alloc: h4_u, h4_v, S2_u, S2_v, H_u, H_v)
  !!   3. An outer do concurrent(j=js:je) local(H_u) loop accumulates H_u(I)
  !!      over a sequential inner k-loop, then normalizes by H_u(I).
  !!      H_u must be private to each j-iteration for correct results.
  !!   4. Analogous do concurrent(J=js-1:je) local(H_v) for v-points.
  !!   5. h4_u / h4_v 3-D arrays are initialized via nested do concurrent(K=2:nz)
  !!      loops accessing G%mask2dCu and G%mask2dCv (allocatable members).
  !!   6. OBC_dir_u / OBC_dir_v integer 2-D arrays are declared locally and
  !!      accessed in the inner k-loop (they remain zero since OBC_friendly=.false.).
  !!
  !! Note: H_v also crashes with local(H_v); it is commented out in the original
  !! code as a workaround. Both are reproduced here with local() to demonstrate
  !! both failure cases.
  subroutine calc_visbeck_mre(h, slope_x, slope_y, N2_u, N2_v, G, GV, CS)
    type(grid_type),      intent(inout) :: G   !< Horizontal grid
    type(vert_grid_type), intent(in)    :: GV  !< Vertical grid
    real, dimension(G%isd:G%ied,   G%jsd:G%jed,   GV%ke),   intent(inout) :: h  !< Layer thickness [H]
    type(mixing_cs_type), intent(inout) :: CS  !< Mixing control structure
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   GV%ke+1), intent(in) :: slope_x  !< i-slope [Z L-1]
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB, GV%ke+1), intent(in) :: slope_y  !< j-slope [Z L-1]
    real, dimension(G%IsdB:G%IedB, G%jsd:G%jed,   GV%ke+1), intent(in) :: N2_u     !< N^2 at u-pts [T-2]
    real, dimension(G%isd:G%ied,   G%JsdB:G%JedB, GV%ke+1), intent(in) :: N2_v     !< N^2 at v-pts [T-2]

    ! Local arrays mapped to the device via !$omp target enter data below.
    real :: h4_u(G%IsdB:G%IedB, G%jsd:G%jed,   GV%ke+1)  ! Thickness weight at u-pts
    real :: h4_v(G%isd:G%ied,   G%JsdB:G%JedB, GV%ke+1)  ! Thickness weight at v-pts
    real :: S2_u(G%IsdB:G%IedB, G%jsd:G%jed)              ! Mean-square slope at u-pts [Z2 L-2]
    real :: S2_v(G%isd:G%ied,   G%JsdB:G%JedB)            ! Mean-square slope at v-pts [Z2 L-2]

    ! 1-D local arrays that must be PRIVATE to each j-iteration.
    ! These are mapped to the device below (H_u also appears in the target map in
    ! the original code), and then privatized with LOCAL in the outer do concurrent.
    ! The combination of the device map and the LOCAL clause triggers a segfault.
    real :: H_u(G%IsdB:G%IedB)  ! Column-integrated thickness weight at u-pts [H]
    real :: H_v(G%isd:G%ied)    ! Column-integrated thickness weight at v-pts [H]

    ! Integer 2-D arrays for OBC direction flags (0 = no OBC, stays zero here
    ! since CS%OBC_friendly = .false.). Accessed in the inner k-loop; NOT mapped.
    ! The compiler generates implicit copyin for these inside the do concurrent.
    integer :: OBC_dir_u(G%IsdB:G%IedB, G%jsd:G%jed)   ! OBC direction flag at u-pts
    integer :: OBC_dir_v(G%isd:G%ied,   G%JsdB:G%JedB) ! OBC direction flag at v-pts

    ! Scalar locals
    real :: S2max       ! Max slope squared for slope clipping [Z2 L-2]
    real :: Hdn, Hup    ! Geometric mean thicknesses [H]
    real :: H_geom      ! sqrt(Hdn * Hup) [H]
    real :: wSE, wNW    ! Thickness weights (SE and NW neighbors)
    real :: wNE, wSW    ! Thickness weights (NE and SW neighbors)
    real :: S2, N2      ! Local slope squared and buoyancy frequency squared

    integer :: is, ie, js, je, nz, i, j, k

    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
    S2max = CS%Visbeck_S_max**2

    ! Initialize OBC direction arrays to zero (CS%OBC_friendly = .false.)
    OBC_dir_u(:,:) = 0
    OBC_dir_v(:,:) = 0

    ! Map the large 3-D arrays, the 2-D output accumulator arrays, and the
    ! 1-D H_u/H_v arrays to the device. Note: H_u and H_v are mapped here
    ! exactly as in the original calc_Visbeck_coeffs_old. This map, combined
    ! with the LOCAL clause in the outer do concurrent below, is the suspected
    ! trigger for the runtime segfault.
    !$omp target enter data map(alloc: h4_u, h4_v, S2_u, S2_v, H_u, H_v)
    !$omp target enter data map(to: CS%SN_u, CS%SN_v)

    ! Initialize h4_u and h4_v: weighted product of the 4 thicknesses surrounding
    ! each interface velocity-point. Uses G%mask2dCu/mask2dCv (allocatable members
    ! of G) without an explicit target data map on G — the compiler generates
    ! implicit copyin for these inside the do concurrent.
    do concurrent(K=2:nz)
      do concurrent( j=js-1:je+1, I=is-1:ie )
        h4_u(I,j,K) = G%mask2dCu(I,j) * ( (h(I,j,K)*h(I+1,j,K)) * (h(I,j,K-1)*h(I+1,j,K-1)) )
      enddo
      do concurrent( J=js-1:je,   i=is-1:ie+1 )
        h4_v(i,J,K) = G%mask2dCv(i,J) * ( (h(i,J,K)*h(i,J+1,K)) * (h(i,J,K-1)*h(i,J+1,K-1)) )
      enddo
    enddo

    ! Outer do concurrent over j for u-points.
    ! UMW NOTE (from original): H_u should be local to each j-iteration since it only
    ! contains i-indices and accumulates column values for j=js:je independently.
    ! However, using local(H_u) here causes a runtime segfault in nvfortran 26.3.
    ! Removing local(H_u) avoids the crash but produces incorrect results.
    do concurrent( j=js:je ) local(H_u)

      ! Initialize accumulators for this j-row
      do concurrent( I=is-1:ie )
        CS%SN_u(I,j) = 0.0
        S2_u(I,j)    = 0.0
        H_u(I)       = 0.0
      enddo

      ! Accumulate over vertical layers
      do K=2,nz ; do concurrent( I=is-1:ie ) &
        & local(Hdn, Hup, H_geom, wSE, wNW, wNE, wSW, S2, N2)

        Hdn    = sqrt( h(I,  j,K) * h(I+1,j,K)   )
        Hup    = sqrt( h(I,  j,K-1) * h(I+1,j,K-1) )
        H_geom = sqrt( Hdn * Hup )

        wSE = h4_v(I+1,j-1,K) ; wNW = h4_v(I,  j,  K)
        wNE = h4_v(I+1,j,  K) ; wSW = h4_v(I,  j-1,K)

        ! Apply OBC direction corrections (all zero when OBC_friendly=.false.)
        if (OBC_dir_u(I,j) == 1) then    ! OBC_DIRECTION_E
          wSE = 0.0 ; wNE = 0.0
          H_geom = sqrt( h(I,j,K) * h(I,j,K-1) )
        elseif (OBC_dir_u(I,j) == -1) then  ! OBC_DIRECTION_W
          wSW = 0.0 ; wNW = 0.0
          H_geom = sqrt( h(I+1,j,K) * h(I+1,j,K-1) )
        endif

        S2 = slope_x(I,j,K)**2 + &
             ( ( (wNW * slope_y(I,  j,  K)**2) + (wSE * slope_y(I+1,j-1,K)**2) ) + &
               ( (wNE * slope_y(I+1,j,  K)**2) + (wSW * slope_y(I,  j-1,K)**2) ) ) / &
             ( ( (wSE + wNW) + (wNE + wSW) ) + GV%H_subroundoff**4 )

        if (S2max > 0.0) S2 = S2 * S2max / (S2 + S2max)
        N2 = max(0.0, N2_u(I,j,K))

        CS%SN_u(I,j) = CS%SN_u(I,j) + sqrt(S2 * N2) * H_geom
        S2_u(I,j)    = S2_u(I,j)    + S2 * H_geom
        H_u(I)       = H_u(I)       + H_geom

      enddo ; enddo  ! K, I

      ! Normalize by the thickness-weighted column sum stored in H_u(I)
      do concurrent( I=is-1:ie )
        if (H_u(I) > 0.0) then
          CS%SN_u(I,j) = G%OBCmaskCu(I,j) * CS%SN_u(I,j) / H_u(I)
          S2_u(I,j)    = G%OBCmaskCu(I,j) * S2_u(I,j)    / H_u(I)
        else
          CS%SN_u(I,j) = 0.0
          S2_u(I,j)    = 0.0
        endif
      enddo

    enddo  ! outer do concurrent(j=js:je) local(H_u)

    ! Outer do concurrent over J for v-points.
    ! UMW NOTE: Same issue as H_u above — local(H_v) causes a segfault.
    ! Commented out in the original code; included here to reproduce both cases.
    do concurrent( J=js-1:je ) local(H_v)

      do concurrent( i=is:ie )
        CS%SN_v(i,J) = 0.0
        S2_v(i,J)    = 0.0
        H_v(i)       = 0.0
      enddo

      do K=2,nz ; do concurrent( i=is:ie ) &
        & local(Hdn, Hup, H_geom, wSE, wNW, wNE, wSW, S2, N2)

        Hdn    = sqrt( h(i,J,  K) * h(i,J+1,K)   )
        Hup    = sqrt( h(i,J,  K-1) * h(i,J+1,K-1) )
        H_geom = sqrt( Hdn * Hup )

        wSE = h4_u(i,  J,  K) ; wNW = h4_u(i-1,J+1,K)
        wNE = h4_u(i,  J+1,K) ; wSW = h4_u(i-1,J,  K)

        if (OBC_dir_v(i,J) == 2) then    ! OBC_DIRECTION_N
          wNW = 0.0 ; wNE = 0.0
          H_geom = sqrt( h(i,J,K) * h(i,J,K-1) )
        elseif (OBC_dir_v(i,J) == -2) then  ! OBC_DIRECTION_S
          wSW = 0.0 ; wSE = 0.0
          H_geom = sqrt( h(i,J+1,K) * h(i,J+1,K-1) )
        endif

        S2 = slope_y(i,J,K)**2 + &
             ( ( (wNW * slope_x(i-1,J+1,K)**2) + (wSE * slope_x(i,  J,  K)**2) ) + &
               ( (wNE * slope_x(i,  J+1,K)**2) + (wSW * slope_x(i-1,J,  K)**2) ) ) / &
             ( ( (wSE + wNW) + (wNE + wSW) ) + GV%H_subroundoff**4 )

        if (S2max > 0.0) S2 = S2 * S2max / (S2 + S2max)
        N2 = max(0.0, N2_v(i,J,K))

        CS%SN_v(i,J) = CS%SN_v(i,J) + sqrt(S2 * N2) * H_geom
        S2_v(i,J)    = S2_v(i,J)    + S2 * H_geom
        H_v(i)       = H_v(i)       + H_geom

      enddo ; enddo  ! K, i

      do concurrent( i=is:ie )
        if (H_v(i) > 0.0) then
          CS%SN_v(i,J) = G%OBCmaskCv(i,J) * CS%SN_v(i,J) / H_v(i)
          S2_v(i,J)    = G%OBCmaskCv(i,J) * S2_v(i,J)    / H_v(i)
        else
          CS%SN_v(i,J) = 0.0
          S2_v(i,J)    = 0.0
        endif
      enddo

    enddo  ! outer do concurrent(J=js-1:je) local(H_v)

    !$omp target exit data map(delete: h4_u, h4_v, S2_u, S2_v, H_u, H_v)
    !$omp target exit data map(from: CS%SN_u, CS%SN_v)

  end subroutine calc_visbeck_mre

end program mre2_do_concurrent_local
