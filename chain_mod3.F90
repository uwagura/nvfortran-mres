module chain_mod3
  use chain_mod4
  use chain_mod4
  #ifdef USE_DERIVED_TYPES
  use derived_types_mod
  #endif
  implicit none

contains

  #ifdef USE_DERIVED_TYPES
  subroutine call_chain3(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz, dt1, dt2)
  #else
  subroutine call_chain3(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz )
  #endif
    implicit none
    integer, intent(in) :: nx, ny, nz
    real, intent(in) :: u(nx, ny, nz), v(nx, ny, nz), h(nx, ny, nz)
    real, intent(out) :: result_u(nx, ny, nz), result_v(nx, ny, nz)
    real, intent(in) :: metric_x(nx, ny), metric_y(nx, ny)
    #ifdef USE_DERIVED_TYPES
    type(my_type1), intent(in) :: dt1
    type(my_type2), intent(in) :: dt2
    #endif

    #ifdef USE_DERIVED_TYPES
    call call_chain4(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz, dt1, dt2)
    #else
    call call_chain4(u, v, h, result_u, result_v, metric_x, metric_y, nx, ny, nz)
    #endif

  end subroutine call_chain3

end module chain_mod3
