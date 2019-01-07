!-*- mode: F90 -*-!

!---------------------------------------------------!
! This file is distributed as part of the DSM code  !
!---------------------------------------------------!

module dsm_constants

   !! This module contains some constants
  
   implicit none
   private

   integer, parameter, public		:: dp = kind(1.0d0)
   real(kind=dp), parameter, public 	:: pi = 3.141592653589793238462643383279_dp
   
   complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp, 1.0_dp)
   complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp, 0.0_dp)
   complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp, 0.0_dp)
 
   real(kind=dp), allocatable, public  :: Ide(:,:)

   !public :: Ide()

contains
  
   !==============================================
   function Ide(n)
   ! To create an identify matrix dims n
   
   implicit none

   integer :: i,n
   allocate (Ide(0:n-1,0:n-1))
 
   Ide = 0.0
   do i = 0, n-1
      Ide(i,i) = 1.0
   enddo
   end function Ide

end module dsm_constants
