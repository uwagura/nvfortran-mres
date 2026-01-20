!> Minimal Reproducible Example (MRE) for testing gradKE-like subroutine calls
!! This program mimics the structure of CorAdCalc with a call to gradKE-like
!! subroutine within its k loop. The intent is to test performance impacts
!! of subroutine calls within nested do concurrent structures.

module constants_mod
  implicit none

  ! Simulation parameters
  integer, parameter :: nx = 288, ny = 288, nz = 100
  integer, parameter :: iterations = 100

  contains

  !> Subroutine that mimics the structure of CorAdCalc
  !! It contains a k-loop with do concurrents over i,j followed by
  !! a call to gradKE_mre
  subroutine coradcalc_mre(u, v, h, CAu, CAv, grid_metric_x, grid_metric_y)

    implicit none
    real, intent(in) :: u(nx, ny, nz), v(nx, ny, nz)
    real, intent(in) :: h(nx, ny, nz)
    real, intent(out) :: CAu(nx, ny, nz), CAv(nx, ny, nz)
    real, intent(in) :: grid_metric_x(nx, ny), grid_metric_y(nx, ny)
    
    real :: u_iminus, u_iplus, v_jminus, v_jplus
    real :: area_factor

    real :: KE(nx, ny), KEx(nx, ny), KEy(nx, ny)
    real :: rel_vort(nx, ny), coriolis_term(nx, ny)
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
    ! !$omp target teams distribute parallel do private(rel_vort, coriolis_term)
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

    !   ! Call to subroutine similar to gradKE
    !   ! This is the problematic call that may hurt GPU performance
    !   call gradKE_mre(u, v, h, KE, KEx, KEy, k)

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
  subroutine gradKE_mre(u, v, h, KE, KEx, KEy, k)
    !$omp declare target

    implicit none

    integer, intent(in) :: k
    real, intent(in) :: u(nx, ny, nz), v(nx, ny, nz)
    real, intent(in) :: h(nx, ny, nz)
    real, intent(out) :: KE(nx, ny), KEx(nx, ny), KEy(nx, ny)

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

  implicit none

  ! Array declarations - main arrays
  real, allocatable :: u(:,:,:), v(:,:,:)
  real, allocatable :: h(:,:,:)
  real, allocatable :: result_u(:,:,:), result_v(:,:,:)

  ! Intermediate arrays
  real, allocatable :: grid_metric_x(:,:), grid_metric_y(:,:)
  integer :: iter
  real :: start_time, end_time

  ! Allocate arrays
  allocate(u(nx, ny, nz))
  allocate(v(nx, ny, nz))
  allocate(h(nx, ny, nz))
  allocate(result_u(nx, ny, nz))
  allocate(result_v(nx, ny, nz))
  allocate(grid_metric_x(nx, ny))
  allocate(grid_metric_y(nx, ny))

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
  !$omp target enter data map(to: nx, ny, nz)
  !$omp target enter data map(to: grid_metric_x, grid_metric_y)

  ! Timing wrapper
  call cpu_time(start_time)

  ! Main loop - iterations for performance testing
  do iter = 1, iterations
    call coradcalc_mre(u, v, h, result_u, result_v, grid_metric_x, grid_metric_y)
  enddo

  call cpu_time(end_time)

  !$omp target exit data map(from: u,v,h)
  !$omp target exit data map(from: result_u, result_v)
  !$omp target exit data map(from: nx, ny, nz)
  !$omp target exit data map(from: grid_metric_x, grid_metric_y)

  print *, "Total time for ", iterations, " iterations: ", (end_time - start_time), " seconds"
  print *, "Average time per iteration: ", (end_time - start_time) / real(iterations), " seconds"
  print *, "Sample result_u value: ", result_u(50, 50, 25)
  print *, "Sample result_v value: ", result_v(50, 50, 25)

  deallocate(u, v, h, result_u, result_v, grid_metric_x, grid_metric_y)

end program MRE_gradKE