!> Minimal Reproducible Example (MRE) for testing calls to subroutines that
!! do work over layers that are called within k-loops inside of target teams
!! regions

module constants_mod
  implicit none

  ! Simulation parameters
  integer, parameter :: nx = 288, ny = 288, nz = 100
  integer, parameter :: iterations = 100

  contains

  !> Subroutine that contains a k-loop with do concurrents over i,j followed by
  !! a call to compute_grad_mre and more do concurrents
  subroutine advection_calc_mre(u, v, h, adv_u, adv_v, metric_x, metric_y, ke, kex, key, curl, coriolis)
    implicit none
    real, intent(in) :: u(nx, ny, nz), v(nx, ny, nz), h(nx, ny, nz)
    real, intent(out) :: adv_u(nx, ny, nz), adv_v(nx, ny, nz)
    real, intent(in) :: metric_x(nx, ny), metric_y(nx, ny)
    real, intent(inout) :: ke(nx, ny), kex(nx, ny), key(nx, ny)
    real, intent(inout) :: curl(nx, ny), coriolis(nx, ny)

    ! Local variables
    integer :: i, j, k
    real :: temp1, temp2

    adv_u = 0.0
    adv_v = 0.0
   
    !$omp target enter data map(alloc: adv_u, adv_v)

    ! Main k-loop structure
    !$omp target teams loop private(curl, coriolis, ke, kex, key)
    do k = 1, nz

      ! Calculate needed quantities 
      do concurrent (j=2:ny-1, i=2:nx-1)
        curl(i, j) = (v(i+1, j, k) - v(i-1, j, k)) * metric_x(i, j) - &
                         (u(i, j+1, k) - u(i, j-1, k)) * metric_y(i, j)
        coriolis(i, j) = curl(i, j) * h(i, j, k)
      enddo
      
      ! Call to subroutine similar to compute_grad
      call compute_grad_mre(u, v, h, ke, kex, key, k)

      ! Do work after function call
      do concurrent (j=2:ny-1, i=2:nx-1)
        temp1 = coriolis(i, j) * v(i, j, k)
        temp2 = kex(i, j) + metric_x(i, j) * h(i, j, k)
        adv_u(i, j, k) = temp1 - temp2
    
        temp1 = coriolis(i, j) * u(i, j, k)
        temp2 = key(i, j) + metric_y(i, j) * h(i, j, k)
        adv_v(i, j, k) = -temp1 - temp2
      enddo

    enddo
    
    !$omp target exit data map(delete: adv_u, adv_v)

  end subroutine advection_calc_mre

  !> Subroutine containing do concurrents over i,j that access multiple indices
  subroutine compute_grad_mre(u, v, h, ke, kex, key, k)
    !$omp declare target
    implicit none

    integer, intent(in) :: k
    real, intent(in) :: u(nx, ny, nz), v(nx, ny, nz), h(nx, ny, nz)
    real, intent(out) :: ke(nx, ny), kex(nx, ny), key(nx, ny)

    ! Local variables
    integer :: i, j
    real :: area_factor

    ! Calculate ke while accessing neighboring indices: i-1, i, j-1, j
    area_factor = 0.25 
    do concurrent (j=2:ny-1, i=2:nx-1)
      ke(i, j) = area_factor * ((u(i-1, j, k) * u(i-1, j, k) + u(i, j, k) * u(i, j, k)) + &
                                (v(i, j-1, k) * v(i, j-1, k) + v(i, j, k) * v(i, j, k))) * &
                 (1.0 / max(h(i-1, j, k), 0.001)) + &
                 (1.0 / max(h(i, j, k), 0.001))
    enddo

    ! Calculate ke gradient in x-direction
    do concurrent (j=2:ny-1, i=2:nx-2)
      kex(i, j) = (ke(i+1, j) - ke(i, j)) * 1.0
    enddo

    ! Calculate ke gradient in y-direction
    do concurrent (j=2:ny-2, i=2:nx-1)
      key(i, j) = (ke(i, j+1) - ke(i, j)) * 1.0
    enddo

  end subroutine compute_grad_mre

end module constants_mod

program MRE_advection
  use constants_mod

  implicit none

  ! Array declarations - main arrays
  real, allocatable :: u(:,:,:), v(:,:,:), h(:,:,:)
  real, allocatable :: result_u(:,:,:), result_v(:,:,:)
  real, allocatable :: metric_x(:,:), metric_y(:,:)
  real, allocatable :: ke(:,:), kex(:,:), key(:,:)
  real, allocatable :: curl(:,:), coriolis(:,:)

  integer :: iter
  real :: start_time, end_time

  ! Allocate arrays
  allocate(u(nx, ny, nz))
  allocate(v(nx, ny, nz))
  allocate(h(nx, ny, nz))
  allocate(result_u(nx, ny, nz))
  allocate(result_v(nx, ny, nz))
  allocate(metric_x(nx, ny))
  allocate(metric_y(nx, ny))
  allocate( ke(nx, ny) )
  allocate( kex(nx, ny) )
  allocate( key(nx, ny) )
  allocate( curl(nx, ny) )
  allocate( coriolis(nx, ny) )


  ! Initialize arrays with random values
  call random_seed()
  call random_number(u)
  call random_number(v)
  call random_number(h)
  call random_number(metric_x)
  call random_number(metric_y)

  u = u * 10.0
  v = v * 10.0
  h = h * 100.0 + 1.0  ! Ensure h > 1
  metric_x = metric_x + 0.5  ! Ensure > 0
  metric_y = metric_y + 0.5  ! Ensure > 0

  result_u = 0.0
  result_v = 0.0

  ! Map variables to device
  !$omp target enter data map(to: u,v,h)
  !$omp target enter data map(to: result_u, result_v)
  !$omp target enter data map(to: metric_x, metric_y)
  !$omp target enter data map(alloc: curl, coriolis )
  !$omp target enter data map(alloc: key, kex, ke)


  ! Timing wrapper
  call cpu_time(start_time)

  ! Main loop - iterations for performance testing
  do iter = 1, iterations
    call advection_calc_mre(u, v, h, result_u, result_v, metric_x, metric_y, ke, kex, key, curl, coriolis)  
  enddo

  call cpu_time(end_time)

  !$omp target exit data map(from: u,v,h)
  !$omp target exit data map(from: result_u, result_v)
  !$omp target exit data map(from: metric_x, metric_y)
  !$omp target exit data map(delete: key, kex, ke)
  !$omp target exit data map(delete: curl, coriolis )

  print *, "Total time for ", iterations, " iterations: ", (end_time - start_time), " seconds"
  print *, "Average time per iteration: ", (end_time - start_time) / real(iterations), " seconds"

  deallocate(u, v, h, result_u, result_v, metric_x, metric_y)
  deallocate( ke, kex, key, curl, coriolis )

end program MRE_advection
