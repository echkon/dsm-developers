!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the dsm code 	     !
!------------------------------------------------------------!

module dsm_wenst
  !! This module calculates scheme for weak and strong dsm.

  use dsm_constants, only: dp
  !use dsm_parameters, only: qu_state, qu_status, qu_noise, qu_dims

  implicit none

  private

  integer			 :: get_dims
  real(kind=dp),allocatable  	 :: prob(:,:,:)
  real(kind=dp),allocatable	 :: rs(:,:,:)
  real(kind=dp),allocatable   	 :: rho10(:,:)
  real(kind=dp),allocatable   	 :: rho11(:,:)
  real(kind=dp),allocatable	 :: qurec(:,:)
  real(kind=dp)			 :: trace,fide

  public :: wenst_run_job
  
contains

  !==================================================================!
  subroutine wenst_run_job
   !==================================================================!
   !! This sub. runs as main program beforehand(note for myself)
   !! This sub. runs some DSM scheme, e.g., weak, strong,...			       
   !===================================================================
   use dsm_constants, only: dp
   use dsm_parameters, only: job,num_iter
   use dsm_qustate, only: true_state

   implicit none

   integer		:: i
   real(kind=dp)        :: itr(1:num_iter),ifide(1:num_iter)
   real(kind=dp)	:: ave_tr,ave_fide
   real(kind=dp)	:: err_tr,err_fide

   get_dims = sqrt(1.0*size(true_state))

   call get_prob
   do i = 1, num_iter
      call get_rho10_rho11
      call get_qurec
      call get_trace_fidelity

      itr(i) = trace
      ifide(i) = fide
   enddo
   ave_tr = sum(itr)/float(num_iter)
   ave_fide = sum(ifide)/float(num_iter)
   
   err_tr = 0.0
   err_fide = 0.0
   do i = 1, num_iter
      err_tr = err_tr + abs(itr(i)-ave_tr)
      err_fide = err_fide + abs(ifide(i)-ave_fide)
   enddo
   err_tr = err_tr/float(num_iter)
   err_fide = err_fide/float(num_iter)

   print *, ave_tr,err_tr,'trace and errors'
   print *, ave_fide,err_fide,'fidelity and errors'
  end subroutine wenst_run_job

!=================================================
!some common rubroutines
!  
  subroutine get_prob
  ! *** to get prob. ***
    use dsm_constants, only: dp,pi,cmplx_i
    use dsm_parameters, only: theta
    use dsm_qustate, only: true_state

    implicit none

    integer             :: x,tx,y,p
    real(kind=dp)       :: sin_tt,sin_2tt
    complex(kind=dp)    :: temp_prob(0:3,0:get_dims-1,0:get_dims-1)
 
    sin_tt  = sin(pi*theta)
    sin_2tt = sin(pi*theta/2.0)**2

    allocate (prob(0:4,0:get_dims-1,0:get_dims-1))

    do x = 0, get_dims-1
       do p = 0, get_dims-1

	  ! *** temp_prob_0 ***
          temp_prob(0,x,p) = 0.0
          do tx = 0, get_dims-1
             do y = 0, get_dims-1
                temp_prob(0,x,p) = temp_prob(0,x,p) + &
                        true_state(tx,y)*exp(2.0*pi*cmplx_i*(y-tx)*p/get_dims)
             enddo
	enddo
          do y = 0, get_dims-1
             temp_prob(0,x,p) = temp_prob(0,x,p) - &
                  2.0*sin_2tt*( &
                  true_state(x,y)*exp(2.0*pi*cmplx_i*(y-x)*p/get_dims)+&
                  true_state(y,x)*exp(-2.0*pi*cmplx_i*(y-x)*p/get_dims))
          enddo
          temp_prob(0,x,p) = temp_prob(0,x,p) + &
                             4*sin_2tt**2*true_state(x,x)
          temp_prob(0,x,p) = temp_prob(0,x,p)/get_dims
          
          ! *** temp_prob_1 ***
          temp_prob(1,x,p) = 0.0
          do y = 0, get_dims-1
             temp_prob(1,x,p) = temp_prob(1,x,p) + &
                  true_state(x,y)*exp(2.0*pi*cmplx_i*(y-x)*p/get_dims)
          enddo
          temp_prob(1,x,p) = temp_prob(1,x,p) - &
                             2*sin_2tt*true_state(x,x)
          temp_prob(1,x,p) = temp_prob(1,x,p)*sin_tt/get_dims
          
          ! *** temp_prob_2 ***
          temp_prob(2,x,p) = conjg(temp_prob(1,x,p))
          
          ! *** temp_prob_3 ***
          temp_prob(3,x,p)= sin_tt**2*true_state(x,x)/get_dims
        enddo
    enddo
    prob(0,:,:) = real((temp_prob(0,:,:)+temp_prob(1,:,:)&
                 +temp_prob(2,:,:)+temp_prob(3,:,:))/2.0) !Prob.+
    prob(1,:,:) = (temp_prob(0,:,:)-temp_prob(1,:,:)&
                  -temp_prob(2,:,:)+temp_prob(3,:,:))/2.0 !Prob.-
    prob(2,:,:) = (temp_prob(0,:,:)+cmplx_i*temp_prob(1,:,:)&
                  -cmplx_i*temp_prob(2,:,:)+temp_prob(3,:,:))/2.0 !Prob.L
    prob(3,:,:) = (temp_prob(0,:,:)-cmplx_i*temp_prob(1,:,:)&
                  +cmplx_i*temp_prob(2,:,:)+temp_prob(3,:,:))/2.0 !Prob.R
    prob(4,:,:) = temp_prob(3,:,:) !Prob.1

   end subroutine get_prob

!===============================
   subroutine get_bisection
   ! *** calcuate bisection ***
    use dsm_constants, only: dp

    implicit none

    integer 		:: i,j,k
    real(kind=dp)	:: r,a,b,c	
    real(kind=dp)       :: fm(0:4,0:get_dims-1,0:get_dims-1)

    if (allocated(rs)) deallocate(rs)
    allocate (rs(0:4,0:get_dims-1,0:get_dims-1))

    do i = 0, 4
       do j = 0, get_dims-1
	  do k = 0, get_dims-1
	     call random_number(r)
!	     a = 0.0
!	     b = xmax(i,j,k)
!	     rs(i,j,k) = 0.0
!11	     c = (a+b)/2.0
!	     fm(i,j,k) = c/xmax(i,j,k)-r
!	     if (abs(fm(i,j,k)) > 0.0001) then
!	        if (fm(i,j,k) > 0.0) then
!		   b = c
!	        else
!		   a = c
!	        endif
!	        goto 11
!	     endif
!	     rs(i,j,k) = r/c
	     rs(i,j,k) = r/prob(i,j,k)
	  enddo
       enddo
    enddo

  end subroutine get_bisection 
!===================================
  subroutine get_rho10_rho11
    ! *** to get rho10 and rho11 ***
    use dsm_utilities, only : utili_get_his
    use dsm_constants, only: dp,cmplx_i
    use dsm_parameters, only: num_copies,job

    implicit none

    integer             :: i,x,p,j,k
    integer             :: N
    real(kind=dp)       :: temp_rs(0:4,0:get_dims-1,0:get_dims-1)
    real(kind=dp)       :: ave_rs(0:4,0:get_dims-1,0:get_dims-1)
    real(kind=dp),allocatable :: Nrs(:,:,:,:)

    if (job == 'strong') N = num_copies/3
    if (job == 'weak')  N = num_copies/2

    if (allocated(rho10)) deallocate(rho10)
    if (allocated(rho11)) deallocate(rho11)
    if (allocated(Nrs)) deallocate(Nrs)

    allocate (rho10(0:get_dims-1,0:get_dims-1))
    allocate (rho11(0:get_dims-1,0:get_dims-1))
    allocate (NRs(1:N,0:4,0:get_dims-1,0:get_dims-1))

    do i = 1, N
       call get_bisection
       Nrs(i,:,:,:) = rs(:,:,:)
    enddo
    ! 
    do i = 0, 4
       do j = 0, get_dims-1
          do k = 0, get_dims-1
             call utili_get_his(Nrs(:,i,j,k),N,0.05d0,ave_rs(i,j,k))
	  enddo
       enddo
    enddo

    do x = 0, get_dims-1
       do p = 0, get_dims-1
          rho10(x,p) = 1/2.0*(ave_rs(0,x,p)-ave_rs(1,x,p)-&
                        cmplx_i*(ave_rs(2,x,p)-ave_rs(3,x,p)))
          rho11(x,p) = ave_rs(4,x,p)
       enddo
    enddo
  end subroutine get_rho10_rho11

 !=============================
 subroutine get_qurec
   ! *** to get quantum state reconstructed ***
   use dsm_constants, only: dp,pi,cmplx_i
   use dsm_parameters, only: theta,job
   use dsm_utilities, only: utili_kron
 
   implicit none
   
   integer 		:: x,y,p
   real(kind=dp)	:: rnorm    

   if (allocated(qurec)) deallocate(qurec)
   allocate (qurec(0:get_dims-1,0:get_dims-1))

   ! *** for strong ***
   if (job == 'strong') then
      do x = 0, get_dims-1
         do y = 0, get_dims-1
	    qurec(x,y) = get_dims*tan(pi*theta/2.0)*&
	 		utili_kron(x,y)*rho11(x,0)
	    do p = 0, get_dims-1
	       qurec(x,y) = qurec(x,y)+exp(2*pi*cmplx_i*(x-y)*p)*&
			 rho10(x,p)
	    enddo
         enddo
      enddo
      qurec = qurec/sin(pi*theta)
   endif
   ! *** for weak ***
   if (job == 'weak') then
    qurec = 0.0
    do x = 0, get_dims-1
         do y = 0, get_dims-1
            do p = 0, get_dims-1
               qurec(x,y) = qurec(x,y)+exp(2*pi*cmplx_i*(x-y)*p)*&
                         rho10(x,p)
            enddo
         enddo
      enddo
      qurec = qurec!/(pi*theta)
   endif
   rnorm = 0.0
   do x = 0, get_dims-1
      rnorm = rnorm + qurec(x,x)
   enddo 
   qurec = qurec/rnorm
 end subroutine get_qurec
 !===================================
 subroutine get_trace_fidelity
  ! *** to calculate the trace distance
  !	and fidelity ***
  use dsm_constants, only: dp,cmplx_1
  use dsm_parameters, only: qu_state
  use dsm_utilities, only: utili_fide
  use dsm_qustate, only: true_state,GHZ,W,Dicke

   implicit none

   integer              :: i
   
   trace = 0.0 
   do i = 0, get_dims-1
      trace = trace+abs(true_state(i,i)-qurec(i,i))
   enddo
   trace = trace/2.0

   ! *** compute fidelity ***
   if (qu_state == 'GHZ') fide = utili_fide(cmplx_1*qurec,GHZ)
   if (qu_state == 'W') fide = utili_fide(cmplx_1*qurec,W)
   if (qu_state == 'Dicke') fide = utili_fide(cmplx_1*qurec,Dicke)

 end subroutine get_trace_fidelity

end module dsm_wenst

