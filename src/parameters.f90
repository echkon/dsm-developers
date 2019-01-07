!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the dsm code 	     !
!------------------------------------------------------------!

module dsm_parameters
  !! This module contains parameters to control the actions of dsm.
  !! Also routines to read the parameters and write them out again.

  use dsm_constants, only: dp
  use dsm_io, only: stdout

  implicit none

  private

  character(len=40), public, save :: qu_state
  character(len=40), public, save :: qu_status
  real(kind=dp), public, save :: qu_noise
  integer, public, save :: qu_dims
  integer, public, save :: qu_nums
  integer, public, save :: qu_exci

  character(len=40), public, save :: job
  real(kind=dp), public, save :: theta
  integer, public, save :: num_copies
  integer, public, save :: num_iter

  logical, public, save :: plot_his
  logical, public, save :: check_cdf
  integer, public, save :: plot_N
  real(kind=dp), public, save :: plot_min
  real(kind=dp), public, save :: plot_max
  
  public :: param_read
  public :: param_write
  public :: param_gen_random
  
contains

   !==================================================================!
   subroutine param_read
      !==================================================================!
      !! Read parameters from the root node only			       !
      !===================================================================
      use dsm_constants, only: dp
      use dsm_io, only: io_file_unit

      implicit none

      call param_in_file

   end subroutine param_read

   !===================================================================
      subroutine param_write
      !==================================================================!
      !! write dsm parameters to stdout
      !===================================================================

      implicit none

      integer :: i

      write (stdout, *)
      write (stdout, '(1x,a)') 'read quantum state'
      write (stdout, *) '------------------'
      write (stdout, '(3a12)') '#qu_state',' #qu_dims/nums','#excitation'
      write (stdout, '(2x,a12,i4,a1,i4,i4)')qu_state,qu_dims,'/',qu_nums,qu_exci

      write (stdout, *)
      write (stdout, '(1x,a)') 'read control parematers'
      write (stdout, *) '-----------------------'
      write (stdout, *) '#job_type #theta'
      write (stdout, *) job, theta 
      write (stdout, *) '#num_copies  #num_iter'
      write (stdout, '(2x,2i8)') num_copies,num_iter

      write (stdout, *)
      write (stdout, '(1x,a)') 'Plot_histogram'
      write (stdout, *) '-----------------------'
      write (stdout, *) '#plot_his  #check_cdf #plot_N  #plot_min  #plot_max'
      write (stdout, '(2x,2l4,i10,2f14.5)') plot_his,check_cdf,plot_N,plot_min,plot_max

   end subroutine param_write

   !==================================================================!
   subroutine param_gen_random
      !==================================================================!
      !! To generate a random number			      
      !===================================================================
      use dsm_constants, only: dp
      
      implicit none

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
      deallocate(seed)

   end subroutine param_gen_random


   !=======================================!
      subroutine param_in_file
      !=======================================!
      !! Load the input.in file
      !=======================================!

      use dsm_io, only: stdout,io_file_unit

      implicit none

      integer           :: in_unit, ios, line
      character(len=40) :: dummy,keyword, temp_qu_state,key_nums
      integer           :: pos

      ios = 0
      line = 0

      in_unit = io_file_unit()
      open (in_unit, file='input.in', form='formatted', status='old')

      do while (ios == 0)
         read(in_unit,'(A)', iostat = ios) dummy
	 if (ios == 0) then
	    line = line + 1
	    pos = scan(dummy, ' ')
	    keyword = dummy(1:pos)
	    dummy = dummy(pos+1:)

	    select case (keyword)

	    ! Read the true quantum state
	    case ('qu_state') 
		read(dummy, *, iostat = ios)temp_qu_state
	    case ('qu_status')
		 read(dummy, *, iostat = ios) qu_status
	    case ('qu_noise')
                 read(dummy, *, iostat = ios) qu_noise

	    ! Read control parameteres
	    case ('job')
		 read(dummy,*, iostat = ios) job
	    case ('theta')
                 read(dummy,*, iostat = ios) theta
	    case ('num_copies')
                 read(dummy, *, iostat = ios) num_copies
	    case ('num_iter')
                 read(dummy, *, iostat = ios) num_iter

	    ! Read histogram data
	    case ('plot_his')
                 read(dummy, *, iostat = ios) plot_his
	    case ('check_cdf')
                 read(dummy, *, iostat = ios) check_cdf
	    case ('plot_N')
                 read(dummy, *, iostat = ios) plot_N
	    case ('plot_min')
                 read(dummy, *, iostat = ios) plot_min
	    case ('plot_max')
                 read(dummy, *, iostat = ios) plot_max

	    case default
	    end select
	 endif
      enddo
      
      pos = scan(temp_qu_state, '_')
      qu_state = temp_qu_state(1:pos-1) !to get state's name
      key_nums = temp_qu_state(pos+1:len(temp_qu_state)) !to get N_k (k excitations)

      pos = scan(key_nums, '_')
      if (pos .eq. 0) then 
	 read(key_nums, *, iostat = ios) qu_nums
      else
         dummy = key_nums(1:pos-1)
         read(dummy, *, iostat = ios) qu_nums

         dummy = key_nums(pos+1:len(key_nums))
         read(dummy, *, iostat = ios) qu_exci
      endif

      !!! qu_nums becomes qu_dims if the state is random
      if (qu_state .eq. 'random') then
	qu_dims = qu_nums
	qu_nums= 0
      endif

      !!! Warming if qu_exci .ne. 0 and qu_state is not Dicke
      if ((qu_state.ne.'Dicke').and.(qu_exci.ne.0)) then
	 print *, ' '
	 print *, '***Warming: your input state include excide. It will be ignored.'
         print *, ' '
	 endif
 
      close (in_unit)

   end subroutine param_in_file
end module dsm_parameters
