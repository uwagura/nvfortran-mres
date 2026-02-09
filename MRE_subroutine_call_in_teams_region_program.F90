program MRE_advection
  use chain_mod1
  #ifdef USE_DERIVED_TYPES
  use derived_types_mod
  #endif

  implicit none

  ! Array dimensions
   integer :: nx = 288, ny = 288, nz = 100

  ! Array declarations - main arrays
  real, allocatable :: u(:,:,:), v(:,:,:), h(:,:,:)
  real, allocatable :: result_u(:,:,:), result_v(:,:,:)
  real, allocatable :: metric_x(:,:), metric_y(:,:)

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
  allocate(metric_x(nx, ny))
  allocate(metric_y(nx, ny))

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
    #ifdef USE_DERIVED_TYPES
    call call_chain1(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz, dt1, dt2)
    #else
    call call_chain1(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz)
    #endif
  enddo

  call cpu_time(end_time)

  !$omp target exit data map(from: u,v,h)
  !$omp target exit data map(from: result_u, result_v)
  !$omp target exit data map(from: metric_x, metric_y)

  print *, "Total time for ", iterations, " iterations: ", (end_time - start_time), " seconds"
  print *, "Average time per iteration: ", (end_time - start_time) / real(iterations), " seconds"

  deallocate(u, v, h, result_u, result_v, metric_x, metric_y)
  #ifdef USE_DERIVED_TYPES
  deallocate(dt1%arr)
  deallocate(dt2%matrix)
  #endif

end program MRE_advection
