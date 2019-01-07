!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the dsm code 	     !
!------------------------------------------------------------!

module dsm_qustate
  !! This module contains the true quantum state for dsm.

  use dsm_constants, only: dp
  !use dsm_parameters, only: qu_state, qu_status, qu_noise, qu_dims

  implicit none

  private

  complex(kind=dp),allocatable,public,save :: true_state(:,:)
  complex(kind=dp),allocatable,public,save :: GHZ(:)
  complex(kind=dp),allocatable,public,save :: W(:)
  complex(kind=dp),allocatable,public,save :: Dicke(:)

  public :: qustate_get_true_state
  public :: qustate_print_true_state
  
contains

  !========================================!
  subroutine qustate_get_true_state
    ! *** to get the true quantum state	***		       
    use dsm_parameters, only: qu_state

    implicit none

    if (qu_state == 'random') call get_random_state
    if (qu_state == 'GHZ') call get_GHZ_state
    if (qu_state == 'W') call get_W_state
    if (qu_state == 'Dicke') call get_D_state

  end subroutine qustate_get_true_state

!==================================!
  subroutine get_random_state
    ! *** to get a random state ***
    use dsm_constants, only: dp,cmplx_i, Ide 
    use dsm_io, only: stdout
    use dsm_parameters, only: qu_state, qu_status, qu_noise, qu_dims

    implicit none

    integer 		:: i,j
    real(kind=dp)	:: r, norm
    real(kind=dp)       :: re_psi(0:qu_dims-1)
    real(kind=dp)       :: im_psi(0:qu_dims-1)
    real(kind=dp) 	:: re_true_state(0:qu_dims-1,0:qu_dims-1)
    real(kind=dp)       :: im_true_state(0:qu_dims-1,0:qu_dims-1)

    allocate (true_state(0:qu_dims-1,0:qu_dims-1))

    do i = 0,qu_dims-1
       call random_number(r)
       re_psi(i) = r
       call random_number(r)
       im_psi(i) = r
    enddo
    norm = 0.0
    do i = 0,qu_dims-1
       norm = norm + re_psi(i)**2+im_psi(i)**2
    enddo
    re_psi = re_psi/sqrt(norm)
    im_psi = im_psi/sqrt(norm)
   
    do i = 0,qu_dims-1
       do j = 0, qu_dims-1
   	 true_state(i,j) = complex(re_psi(i),im_psi(i))*complex(re_psi(j),-im_psi(j))
       enddo
    enddo
    
    true_state = (1-qu_noise)*true_state + qu_noise*Ide(qu_dims)/qu_dims

    if (qu_status == 'mixed') then
       re_true_state = 0.0
       im_true_state = 0.0

       do i = 0,qu_dims-1
          do j = i,qu_dims-1
             call random_number(r)
             re_true_state(i,j) = r
   	     if (i .ne. j) then 
   	       !only i ≠ j are non_zero
               call random_number(r)
               im_true_state(i,j) = r
   	     endif
          enddo
       enddo

       norm = 0.0
       do i = 0,qu_dims-1
          norm = norm + re_true_state(i,i)
       enddo

       do i = 0, qu_dims-1
	  do j = 0, qu_dims-1
       	     true_state(i,j) = (re_true_state(i,j)+cmplx_i*im_true_state(i,j))/norm
	     true_state(j,i) = conjg(true_state(i,j))
	  enddo
       enddo
    
    endif       

  end subroutine get_random_state
  
  !===================================================================
  subroutine get_GHZ_state
    !==================================================================!
    !! To get a GHZ state
    !===================================================================

    use dsm_constants, only: dp,cmplx_i, Ide
    use dsm_io, only: stdout
    use dsm_parameters, only: qu_state, qu_noise, qu_nums

    implicit none

    integer             	    :: i,j
    real(kind=dp)     		    :: r, norm
    integer		     	    :: temp_qu_dims

    temp_qu_dims = 2**qu_nums
    allocate(true_state(0:temp_qu_dims-1,0:temp_qu_dims-1))
    allocate(GHZ(0:temp_qu_dims-1))

    GHZ = 0.0
    GHZ(0) = 1.0/sqrt(2.0)
    GHZ(temp_qu_dims-1) = 1.0/sqrt(2.0)

    true_state = 0.0
    true_state(0,0) = 0.5
    true_state(0,temp_qu_dims-1) = 0.5
    true_state(temp_qu_dims-1,0) = 0.5
    true_state(temp_qu_dims-1,temp_qu_dims-1) = 0.5

    true_state = (1-qu_noise)*true_state + qu_noise*Ide(temp_qu_dims)/temp_qu_dims 
 
  end subroutine get_GHZ_state

!===================================================================
  subroutine get_W_state
    !==================================================================!
    !! To get a W state
    !===================================================================

    use dsm_constants, only: dp,cmplx_i, Ide
    use dsm_io, only: stdout
    use dsm_parameters, only: qu_state, qu_noise, qu_nums

    implicit none

    integer               :: i,j
    real(kind=dp)         :: r, norm
    integer               :: temp_qu_dims

    temp_qu_dims = 2**qu_nums
    allocate (true_state(0:temp_qu_dims-1,0:temp_qu_dims-1))
    allocate(W(0:temp_qu_dims-1))

    true_state = 0.0
    W = 0

    W(0) = 0
    W(1) = 1
    W(2) = 1
    W(3) = 0
    do i = 3, qu_nums
       W(2**(i-1)) = 1
    enddo
    
    do i = 0, temp_qu_dims - 1
       do j = 0, temp_qu_dims - 1
          true_state(i,j) = W(i)*W(j)
       enddo
    enddo
    true_state = (1-qu_noise)*true_state + qu_noise*Ide(temp_qu_dims)/temp_qu_dims 
    
  end subroutine get_W_state 

!===================================================================
  subroutine get_D_state
    !==================================================================!
    !! To get a Dicke state. 
    !! D = |000..0111...1> ||(qu_nums-qu_exci)/2 and (qu_nums+qu_exci)/2
    !===================================================================

    use dsm_constants, only: dp,cmplx_i, Ide
    use dsm_io, only: stdout
    use dsm_parameters, only: qu_state, qu_noise, qu_nums,qu_exci
    use dsm_utilities

    implicit none

    integer               :: i,j
    real(kind=dp)         :: r, norm
    integer               :: temp_qu_dims
    integer		  :: vecu(2)
    integer		  :: vecd(2)
    integer		  :: vecI(1) = 1
    integer,allocatable   :: temp_vec(:),temp_vec_0(:)

    temp_qu_dims = 2**qu_nums
    allocate (true_state(0:temp_qu_dims-1,0:temp_qu_dims-1))

    vecu(1) = 1
    vecu(2) = 0
    vecd(1) = 0
    vecd(2) = 1

    true_state = 0.0
 
    if (qu_exci == 0) then
	allocate (temp_vec(0:temp_qu_dims-1))
        temp_vec = 0
        temp_vec(temp_qu_dims-1) = 1
    endif
   
    if (qu_nums == qu_exci) then
        allocate (temp_vec(0:temp_qu_dims-1))
        temp_vec = 0
        temp_vec(0) = 1
    endif

    if ((qu_exci.ne.0).and.(qu_nums.ne.qu_exci)) then   
       temp_vec = vecu    
       do i = 1,qu_nums-qu_exci-1
          call utili_tensor_prod(temp_vec_0,temp_vec,vecu)
          if (allocated(temp_vec)) deallocate(temp_vec)
          allocate (temp_vec(size(temp_vec_0)))
          temp_vec = temp_vec_0 
       enddo
       do i = 1,qu_exci
          call utili_tensor_prod(temp_vec_0,temp_vec,vecd)
          if (allocated(temp_vec)) deallocate(temp_vec)
          allocate (temp_vec(size(temp_vec_0)))
          temp_vec = temp_vec_0
       enddo
    endif

    Dicke = temp_vec

    do i = 0, temp_qu_dims - 1
       do j = 0, temp_qu_dims - 1
          true_state(i,j) = temp_vec(i)*temp_vec(j)
       enddo
    enddo
    true_state = (1-qu_noise)*true_state + qu_noise*Ide(temp_qu_dims)/temp_qu_dims
    
  end subroutine get_D_state 
  
  !==================================================================!
  subroutine qustate_print_true_state
    !==================================================================!
    !! To print the true quantum state
    !===================================================================
    use dsm_parameters, only: qu_state,qu_dims,qu_nums
    use dsm_io, only: stdout

    implicit none

    integer 	:: i,j
    
    if (qu_dims == 0) qu_dims = 2**qu_nums
    write(stdout,*)''
    write(stdout,'(3a10)')'Print the ',qu_state,' state.'
    do i = 0, qu_dims-1
       do j = 0, qu_dims-1
	  write(stdout,'(a11,i2,a1,i2,a3,2F10.5)') 'true_state[',i,',',j,']=',true_state(i,j)
        enddo
    enddo
  !101 format(a10,i2,a1,i2,a3,F10)
  end subroutine qustate_print_true_state

end module dsm_qustate

