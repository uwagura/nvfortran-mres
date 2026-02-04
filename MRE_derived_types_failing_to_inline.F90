!> Minimal Reproducible Example (MRE) for testing gradKE-like subroutine calls
!! This program mimics the structure of CorAdCalc with a call to gradKE-like
!! subroutine within its k loop. The intent is to test performance impacts
!! of subroutine calls within nested do concurrent structures.

module constants_mod
  #ifdef USE_DERIVED_TYPES
  use derived_types_mod
  #endif
  implicit none

  ! Simulation parameters
  integer, parameter :: iterations = 100

  contains

  !> Subroutine that mimics the structure of CorAdCalc
  !! It contains a k-loop with do concurrents over i,j followed by
  !! a call to gradKE_mre
  #ifdef USE_DERIVED_TYPES
  subroutine coradcalc_mre(u, v, h, CAu, CAv, grid_metric_x, grid_metric_y, nx, ny, nz, dt1, dt2)
  #else
  subroutine coradcalc_mre(u, v, h, CAu, CAv, grid_metric_x, grid_metric_y, nx, ny, nz)
  #endif

    implicit none
    integer, intent(in) :: nx, ny, nz
    ! Double gyre uses dynamic memory, so, the SZ?_ statements expand to this
    real, intent(in) :: u(:,:,:), v(:,:,:)
    real, intent(in) :: h(:,:,:)
    real, intent(out) :: CAu(:,:,:), CAv(:,:,:)
    real, intent(in) :: grid_metric_x(:,:), grid_metric_y(:,:)
  #ifdef USE_DERIVED_TYPES
    type(my_type1), intent(in) :: dt1
    type(my_type2), intent(in) :: dt2
  #endif
    
    real :: u_iminus, u_iplus, v_jminus, v_jplus
    real :: area_factor

    real :: KE(nx,ny), KEx(nx,ny), KEy(nx,ny)
    real :: rel_vort(nx,ny), coriolis_term(nx,ny)
    integer :: i, j, k
    real :: temp1, temp2

    CAu = 0.0
    CAv = 0.0

    !$omp target enter data map(alloc: CAu, CAv)
    !$omp target enter data map(alloc: rel_vort, coriolis_term, temp1, temp2)
    !$omp target enter data map(alloc: KEy, KEx, KE)

    !$omp target enter data map(alloc: u_iminus, u_iplus, v_jminus, v_jplus)
    !$omp target enter data map(alloc: area_factor)

    ! Main k-loop structure similar to CorAdCalc
    !  !$omp target teams distribute parallel do private(rel_vort, coriolis_term, KE, KEx, KEy)
    !$omp target teams loop private(rel_vort, coriolis_term, KE, KEx, KEy)
    do k = 1, nz

      ! First do concurrent block - calculate relative vorticity
      do concurrent (j=2:ny-1, i=2:nx-1)
        rel_vort(i, j) = (v(i+1, j, k) - v(i-1, j, k)) * grid_metric_x(i, j) - &
                         (u(i, j+1, k) - u(i, j-1, k)) * grid_metric_y(i, j)
      enddo

      ! Second do concurrent block - initialize working arrays
      do concurrent (j=2:ny-1, i=2:nx-1)
        coriolis_term(i, j) = rel_vort(i, j) * h(i, j, k)
      enddo
      
      ! Call to subroutine similar to gradKE
      ! This is the problematic call that may hurt GPU performance
      #ifdef USE_DERIVED_TYPES
      call gradKE_mre(u, v, h, KE, KEx, KEy, k, nx, ny, nz, dt1, dt2)
      #else
      call gradKE_mre(u, v, h, KE, KEx, KEy, k, nx, ny, nz)
      #endif


      ! Third do concurrent block - use results from gradKE call
      do concurrent (j=2:ny-1, i=2:nx-1)
        temp1 = coriolis_term(i, j) * v(i, j, k)
        temp2 = KEx(i, j) + grid_metric_x(i, j) * h(i, j, k)
        CAu(i, j, k) = temp1 - temp2
      enddo

      ! Fourth do concurrent block - similar computation for v component
      do concurrent (j=2:ny-1, i=2:nx-1)
        temp1 = coriolis_term(i, j) * u(i, j, k)
        temp2 = KEy(i, j) + grid_metric_y(i, j) * h(i, j, k)
        CAv(i, j, k) = -temp1 - temp2
      enddo

    enddo

    !$omp target exit data map(delete: u_iminus, u_iplus, v_jminus, v_jplus)
    !$omp target exit data map(delete: area_factor)

    !$omp target exit data map(delete: KEy, KEx, KE)
    !$omp target exit data map(delete: rel_vort, coriolis_term, temp1, temp2)
    !$omp target exit data map(delete: CAu, CAv)

  end subroutine coradcalc_mre

  !> Subroutine that mimics the structure of gradKE
  !! It contains do concurrents over i,j that access multiple indices
  !! (like i-1, i+1, j-1, j+1) similar to the actual gradKE implementation
  #ifdef USE_DERIVED_TYPES
  subroutine gradKE_mre(u, v, h, KE, KEx, KEy, k, nx, ny, nz, dt1, dt2)
  #else
  subroutine gradKE_mre(u, v, h, KE, KEx, KEy, k, nx, ny, nz)
  #endif
    !$omp declare target

    implicit none

    integer, intent(in) :: k, nx, ny, nz
    real, intent(in) :: u(:,:,:), v(:,:,:)
    real, intent(in) :: h(:,:,:)
    #ifdef USE_DERIVED_TYPES
    type(my_type1), intent(in) :: dt1
    type(my_type2), intent(in) :: dt2
    #endif
    real, intent(out) :: KE(:,:), KEx(:,:), KEy(:,:)

    integer :: i, j
    real :: u_iminus, u_iplus, v_jminus, v_jplus
    real :: area_factor

    ! Calculate KE (Kinetic energy)
    ! This loop mimics the KE_ARAKAWA scheme in gradKE
    ! It accesses neighboring indices: i-1, i, j-1, j
    do concurrent (j=2:ny-1, i=2:nx-1)
      ! Access multiple indices similar to actual gradKE
      u_iminus = u(i-1, j, k)
      u_iplus = u(i, j, k)
      v_jminus = v(i, j-1, k)
      v_jplus = v(i, j, k)

      ! Compute KE with multi-index access
      area_factor = 0.25  ! Mimics grid metric factors
      KE(i, j) = area_factor * ((u_iminus * u_iminus + u_iplus * u_iplus) + &
                                (v_jminus * v_jminus + v_jplus * v_jplus)) * &
                 (1.0 / max(h(i-1, j, k), 0.001)) + &
                 (1.0 / max(h(i, j, k), 0.001))
    enddo

    ! Calculate KE gradient in x-direction
    ! This mimics: KEx(I,j) = (KE(i+1,j) - KE(i,j)) * G%IdxCu(I,j)
    ! Note: bounds are 2:nx-2 to avoid accessing KE(nx, j) which would be out of bounds
    do concurrent (j=2:ny-1, i=2:nx-2)
      KEx(i, j) = (KE(i+1, j) - KE(i, j)) * 1.0  ! 1.0 mimics grid metric
    enddo

    ! Calculate KE gradient in y-direction
    ! This mimics: KEy(i,J) = (KE(i,j+1) - KE(i,j)) * G%IdyCv(i,J)
    ! Note: bounds are 2:ny-2 to avoid accessing KE(i, ny) which would be out of bounds
    do concurrent (j=2:ny-2, i=2:nx-1)
      KEy(i, j) = (KE(i, j+1) - KE(i, j)) * 1.0  ! 1.0 mimics grid metric
    enddo

  end subroutine gradKE_mre

end module constants_mod

program MRE_gradKE
  use constants_mod
  #ifdef USE_DERIVED_TYPES
  use derived_types_mod
  #endif

  implicit none

  ! Declare parameters at program start
  integer , parameter :: nx = 288, ny = 288, nz = 100
  
  ! Array declarations - main arrays
  real, allocatable :: u(:,:,:), v(:,:,:)
  real, allocatable :: h(:,:,:)
  real, allocatable :: result_u(:,:,:), result_v(:,:,:)

  ! Intermediate arrays
  real, allocatable :: grid_metric_x(:,:), grid_metric_y(:,:)
  integer :: iter
  real :: start_time, end_time

  ! Derived type variables
  #ifdef USE_DERIVED_TYPES
  type(my_type1) :: dt1
  type(my_type2) :: dt2
  #endif

  ! Allocate arrays
  allocate(u(nx, ny, nz))
  allocate(v(nx, ny, nz))
  allocate(h(nx, ny, nz))
  allocate(result_u(nx, ny, nz))
  allocate(result_v(nx, ny, nz))
  allocate(grid_metric_x(nx, ny))
  allocate(grid_metric_y(nx, ny))

  ! Allocate and initialize derived type fields
  #ifdef USE_DERIVED_TYPES
  dt1%n = 10
  dt1%scalar = 1.23
  allocate(dt1%arr(dt1%n))
  dt1%arr = 42.0

  dt2%m = 5
  dt2%n = 8
  dt2%value = 3.14
  allocate(dt2%matrix(dt2%m, dt2%n))
  dt2%matrix = 7.0
  #endif

  ! Initialize arrays with random values
  call random_seed()
  call random_number(u)
  call random_number(v)
  call random_number(h)
  call random_number(grid_metric_x)
  call random_number(grid_metric_y)

  u = u * 10.0
  v = v * 10.0
  h = h * 100.0 + 1.0  ! Ensure h > 1
  grid_metric_x = grid_metric_x + 0.5  ! Ensure > 0
  grid_metric_y = grid_metric_y + 0.5  ! Ensure > 0

  result_u = 0.0
  result_v = 0.0

  ! Map variables to device
  !$omp target enter data map(to: u,v,h)
  !$omp target enter data map(to: result_u, result_v)
  ! !$omp target enter data map(to: nx, ny, nz)
  !$omp target enter data map(to: grid_metric_x, grid_metric_y)

  ! Timing wrapper
  call cpu_time(start_time)

  ! Main loop - iterations for performance testing
  do iter = 1, iterations
    #ifdef USE_DERIVED_TYPES
    call coradcalc_mre(u, v, h, result_u, result_v, grid_metric_x, grid_metric_y, nx, ny, nz, dt1, dt2)
    #else
    call coradcalc_mre(u, v, h, result_u, result_v, grid_metric_x, grid_metric_y, nx, ny, nz)
    #endif
  enddo

  call cpu_time(end_time)

  !$omp target exit data map(from: u,v,h)
  !$omp target exit data map(from: result_u, result_v)
  ! !$omp target exit data map(from: nx, ny, nz)
  !$omp target exit data map(from: grid_metric_x, grid_metric_y)

  print *, "Total time for ", iterations, " iterations: ", (end_time - start_time), " seconds"
  print *, "Average time per iteration: ", (end_time - start_time) / real(iterations), " seconds"
  print *, "Sample result_u value: ", result_u(50, 50, 25)
  print *, "Sample result_v value: ", result_v(50, 50, 25)

  deallocate(u, v, h, result_u, result_v, grid_metric_x, grid_metric_y)
  #ifdef USE_DERIVED_TYPES
  deallocate(dt1%arr)
  deallocate(dt2%matrix)
  #endif

end program MRE_gradKE
