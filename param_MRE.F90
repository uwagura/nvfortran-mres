! Minimal example to reproduce issue reading arrays
! nested inside of derived types. Compile with: 
!
! nvfortran -Minfo=inline -o param_MRE.exe param_MRE.F90
! 
! to see the expected output, and then compile with 
! 
! nvfortran -Minfo=inline -Minline -o param_MRE.exe param_MRE.F90
!
! to see the error
module param_mre
  implicit none

  ! Derived type for a line of text
  type :: file_line_type
    character(len=:), allocatable :: line
  end type file_line_type

  ! Derived type for file data (array of file lines)
  type :: file_data_type
    type(file_line_type), allocatable :: fln(:)
    integer :: num_lines = 0
  end type file_data_type

  ! Derived type for parameter file (array of file_data_type)
  type :: param_file_type
    type(file_data_type), allocatable :: param_data(:)
    integer :: nfiles = 0
  end type param_file_type

contains


  ! Subroutine to read lines from a file and store in param_file_type
  subroutine populate_param_data(filename, pf)
    character(len=*), intent(in) :: filename
    type(param_file_type), intent(inout) :: pf
    integer :: iunit, n, i
    character(len=1024) :: buffer
    integer :: ios

    ! Count the number of lines in the file
    open(newunit=iunit, file=filename, status='old', action='read')
    n = 0
    do
      read(iunit, '(A)', iostat=ios) buffer
      if (ios /= 0) exit
      n = n + 1
    end do
    
    ! Store information about the number of files, and the number
    ! of lines in the file
    rewind(iunit)
    pf%nfiles = 1
    allocate(pf%param_data(1))
    pf%param_data(1)%num_lines = n
    allocate(pf%param_data(1)%fln(n))

    ! Actually store the lines of the file into to the data structure
    do i = 1, n
      read(iunit, '(A)', iostat=ios) buffer
      pf%param_data(1)%fln(i)%line = trim(buffer)
    end do
    close(iunit)
  end subroutine populate_param_data

  ! Function to find the longest line in the parameter file
  function max_input_line_length(pf) result(max_len)
    type(param_file_type), intent(in) :: pf
    integer :: max_len, i
    max_len = 0
    
    ! Iterate throught store lines to find the longest one. 
    ! Print statements to debug issue when inlining this function
    print *, "Line values inside of max_input_line_length:"
    do i = 1, pf%param_data(1)%num_lines
      print *, "\tLine ", i, ": '", pf%param_data(1)%fln(i)%line, "'"
      max_len = max(max_len, len_trim(pf%param_data(1)%fln(i)%line))
    end do
  end function max_input_line_length

end module param_mre

program test_mre
  use param_mre
  implicit none

  type(param_file_type) :: pf
  integer :: i, max_len

  ! Read parameter file
  call populate_param_data("param_MRE_input.txt", pf)

  ! Print lines after populate_param_data
  print *, "After populate_param_data:"
  do i = 1, pf%param_data(1)%num_lines
    print *, "\tLine ", i, ": '", pf%param_data(1)%fln(i)%line, "'"
  end do

  ! Call max_input_line_length and print lines again
  max_len = max_input_line_length(pf)
  print *, "Max line length: ", max_len

end program test_mre
