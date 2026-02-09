module layer_routine_mod
  implicit none
  contains

  subroutine layer_routine(full, nx, ny, nz)
    real, intent(inout) :: full(:,:,:)
    integer, intent(in) :: nx, ny, nz
    
    ! local variables
    real, dimension(nx,ny) :: layer
    integer :: i, j, k


    ! Work on layers
    #ifdef USE_LOCAL
    do concurrent ( k = 1:nz ) local( layer )
    #else
    do concurrent ( k = 1:nz )
    #endif 
      ! Give layer unique values for this layer
      do concurrent ( i = 1:nx, j = 1:ny)
          layer(i, j) = real(k)
      enddo

      ! Set total to the value of the layer
      do concurrent ( i = 1:nx, j = 1:ny)
          full(i, j, k) = layer(i, j)
      enddo

    enddo

  end subroutine layer_routine

end module layer_routine_mod


program MRE_do_concurrent_local
  use layer_routine_mod
  implicit none
  
  ! Parameters
  integer, parameter :: nz = 100, ny = 100, nx = 100
  
  ! Variables
  integer :: k,j,i
  real , allocatable :: total(:,:,:)
  ! real , allocatable :: layer(:,:)

  ! Allocate arrays
  allocate( total(nx, ny ,nz ) )
  ! allocate( layer(nx, ny ) )

  ! Initialize array values
  do concurrent( i = 1:nx, j = 1:ny, k = 1:nz)
    total(i, j, k) = 0.0
  enddo

  ! ! Initialize layer to zero
  ! do concurrent( i = 1:nx, j = 1:ny )
  !   layer(i, j) = 0.0
  ! enddo

  call layer_routine( total, nx, ny, nz )

  ! Print sum of total for verification
  PRINT *, "Sum of total: ", sum(total)

  ! Deallocate arrays
  deallocate( total )
  ! deallocate( layer )

end program MRE_do_concurrent_local
  
  
