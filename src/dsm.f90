!-*- program: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       DSM	                             !
!                                                            !
!          The code for direct state measurememnt            !
!                 written by: Le Bin Ho       	             !
!                                                            !
! Please cite                                                !
!                                                            !
! [ref] L. B. Ho, 					     !
!       "Improving direct state measurements                 !
!       by using rebits in real enlarged Hilbert spaces,"    !
!       Physics Letter A    				     !
!       https://doi.org/10.1016/j.physleta.2018.10.047 	     !
!                                                            !
! The webpage of the DSM code is www.		             !
!                                                            !
! The DSM code is hosted on GitHub:               	     !
!                                                            !
! https://github.com/				             !
!------------------------------------------------------------!

program dsmpr
  !! The main dsm program

  use dsm_constants
  use dsm_parameters
  use dsm_io
  use dsm_utilities
  use dsm_qustate
  use dsm_wenst

  implicit none

  character(len=11)  :: cdate, ctime

  integer	:: i
  real(kind=dp) :: r,ave
  real(kind=dp),allocatable :: f(:)

  call utili_gen_random

  stdout = io_file_unit()
  open (unit=stdout, file='output.out')
  write(stdout,*)
  write(stdout,*)'================================================='

  call io_date(cdate,ctime)
  write (stdout, *) 'dsm: Execution started on ',cdate, ' at ', ctime

  call param_read
  call param_write
  call qustate_get_true_state
  call qustate_print_true_state
  call wenst_run_job

  allocate(f(1:plot_N))
  if (plot_his) then
     do i = 1, plot_N
        call random_number(r)
        f(i) = r
     enddo
  endif
  if (check_cdf) then
     do i = 1, plot_N
        call random_number(r)
        f(i) = -log(1-r) 
     enddo
  endif
     call utili_get_his_save(f,ave)
     write (stdout, *)'Average of histogram',ave


  close (stdout)

end program dsmpr

