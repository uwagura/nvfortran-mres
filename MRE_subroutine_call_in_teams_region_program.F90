program MRE_advection
  use constants_mod

  implicit none

  ! Array declarations - main arrays
  real, allocatable :: u(:,:,:), v(:,:,:), h(:,:,:)
  real, allocatable :: result_u(:,:,:), result_v(:,:,:)
  real, allocatable :: metric_x(:,:), metric_y(:,:)

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

  ! Timing wrapper
  call cpu_time(start_time)

  ! Main loop - iterations for performance testing
  do iter = 1, iterations
    call advection_calc_mre(u, v, h, result_u, result_v, metric_x, metric_y)
  enddo

  call cpu_time(end_time)

  !$omp target exit data map(from: u,v,h)
  !$omp target exit data map(from: result_u, result_v)
  !$omp target exit data map(from: metric_x, metric_y)

  print *, "Total time for ", iterations, " iterations: ", (end_time - start_time), " seconds"
  print *, "Average time per iteration: ", (end_time - start_time) / real(iterations), " seconds"

  deallocate(u, v, h, result_u, result_v, metric_x, metric_y)

end program MRE_advection
