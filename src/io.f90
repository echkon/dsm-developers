!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the the DSM code	     !
!------------------------------------------------------------!

module dsm_io
  !! Module to handle operations related to file input and output.

  use dsm_constants, only: dp
  implicit none

  private

!#ifdef MPI
!  include 'mpif.h'
!#endif

  integer, public, save           :: stdout
  character(len=50), public, save :: seedname
  integer, parameter, public      :: maxlen = 255
  logical, public, save           :: post_proc_flag
  character(len=10), public, parameter:: dsm_version = '1.0.0'

  public :: io_file_unit
  public :: io_date
  !public :: io_warning

contains

  !===========================================
  function io_file_unit()
    !===========================================
    !
    !! Returns an unused unit number
    !! so we can later open a file on that unit.
    !
    !===========================================
    implicit none

    integer :: io_file_unit, unit
    logical :: file_open

    unit = 9
    file_open = .true.
    do while (file_open)
      unit = unit + 1
      inquire (unit, OPENED=file_open)
    end do

    io_file_unit = unit

    return
  end function io_file_unit

  !=======================================================
  subroutine io_date(cdate, ctime)
  !=======================================================
  !
  !! Returns two strings containing the date and the time
  !
  !=======================================================
    implicit none
    character(len=11), intent(out) :: cdate
    !! The date
    character(len=11), intent(out) :: ctime
    !! The time

    character(len=3), dimension(12) :: months
    data months/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
    integer date_time(8)
    !
    call date_and_time(values=date_time)
    !
    write (cdate, '(i2,"-",a3,"-",i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime, '(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)

  end subroutine io_date

  !=======================================================
  !subroutine io_warning
  !=======================================================
  !
  !! To warn something
  !
  !=======================================================
  !use dsm_parameters, only: qu_status, qu_noise

  !implicit none
  !character(len=100) :: warn
 
  !if ((qu_status == 'mixed').and.(qu_noise .ne. 0.0)) then 
  !   print *, ' '
  !   print *, '=========================================='
  !   print *, ' Warning: You are runing mixed state with noise.'
  !   print *, ' In this case, qu_noise will be ignored.'
  !   print *, '=========================================='
  !   print *, ' '
  !endif
  
  !end subroutine io_warning

end module dsm_io
