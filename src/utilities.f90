!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the dsm code 	     !
!------------------------------------------------------------!

module dsm_utilities
  !! This module provices some utilities:
  !! Tensor product 
  !! Generate a random number and plot histogram

  !use dsm_parameters, only: plot_his,plot_data,plot_min,plot_max
  !use dsm_constants, only: dp
  !use dsm_io, only: stdout

  implicit none

  private

  public :: utili_fide
  public :: utili_kron  
  public :: utili_tensor_prod
  public :: utili_gen_random
  public :: utili_get_his_save
  public :: utili_get_his

contains

   !=====================================================!
   function utili_fide(a,b)
     ! *** to provide reference fidelity of a and b (b pure) ***
     use dsm_constants, only: dp

     implicit none

     real(kind=dp)	 :: utili_fide
     complex(kind=dp)	 :: a(:,:),b(:)
     real(kind=dp),allocatable :: tb(:,:)
     integer		 :: i,j, get_dims

     get_dims = size(b)
     utili_fide = 0.0
     do i = 1, get_dims !when in function, matrix a,b run from 1
	do j = 1, get_dims
     	   utili_fide = utili_fide + conjg(b(i))*a(i,j)*b(j)
	enddo
     enddo

   end function

   !====================================!
   function utili_kron(i,j)
     ! *** to provide kronecker delta ***
     use dsm_constants, only: dp

     implicit none
     
     integer 	:: utili_kron, i,j

     if (i==j) then
        utili_kron = 1
     else
	utili_kron = 0
     endif
 
   end function


  !===========================================!
   subroutine utili_tensor_prod(vecP,vecA,vecB)
     ! *** to provide tensor product ***
     use dsm_constants, only: dp
 
     implicit none

     integer :: i,j,r,ra,rb,c,ca,cb
     integer :: vecA(:),vecB(:)
     integer,allocatable :: vecP(:)

     ra = size(vecA)
     rb = size(vecB)

     if (allocated(vecP)) deallocate(vecP)
     allocate (vecP(ra*rb))
     r = 0
     do i = 1, ra
	vecP(r+1:r+rb) = vecA(i)*vecB
	r = r + rb
     enddo
     end subroutine utili_tensor_prod

   !======================================!
      subroutine utili_gen_random
      ! *** to generate a random number ***		      
      
      implicit none

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
      deallocate(seed)

   end subroutine utili_gen_random

   
   !=======================================!
      subroutine utili_get_his_save(f,ave)
      !=======================================!
      !! To get histogram of f and average 
      !  and save to a file, using for check his
      !=======================================!

      use dsm_constants, only: dp
      use dsm_parameters, only: plot_his,check_cdf,plot_min,plot_max
      use dsm_io, only: stdout,io_file_unit

      implicit none

      integer           :: in_unit
      integer		:: i,temp,get_dims
      real(kind=dp)	:: step,lower_bound, upper_bound
      real(kind=dp)	:: ave
      real(kind=dp)     :: f(:)

      get_dims = size(f)
      step = 0.005
      i = 1

      in_unit = io_file_unit()
      open (in_unit, file='hist.out')
      write(in_unit,*) '#data for plot histogram'
      ave = 0.0

1026  lower_bound = plot_min + (i-1)*step
      upper_bound = plot_min + i*step
      temp = count((f > lower_bound).and.(f <= upper_bound))
      ave  = ave + temp
      
      if (plot_his.or.check_cdf) write(in_unit,*)(lower_bound+step/2.0),temp/(step*get_dims)

      if (upper_bound < plot_max) then
	 i = i + 1
	 goto 1026
      endif
      ave = ave/(step*get_dims)/float(i)

      close (in_unit)

   end subroutine utili_get_his_save

  !=======================================!
      subroutine utili_get_his(f,dims,step,ave)
      !=======================================!
      !! To get average of f
      !=======================================!

      use dsm_constants, only: dp
      use dsm_parameters, only: plot_his,check_cdf,plot_min,plot_max
      use dsm_io, only: stdout,io_file_unit

      implicit none

      integer           :: in_unit
      integer           :: i,temp,dims
      real(kind=dp)     :: step,lower_bound, upper_bound
      real(kind=dp)     :: ave
      real(kind=dp)     :: f(:)

      i = 1
      ave = 0.0

107  lower_bound = plot_min + (i-1)*step
      upper_bound = plot_min + i*step
      temp = count((f > lower_bound).and.(f <= upper_bound))
      ave  = ave + temp

      if (upper_bound < plot_max) then
         i = i + 1
         goto 107
      endif
      ave = ave/(step*dims)/float(i)

   end subroutine utili_get_his

end module dsm_utilities
