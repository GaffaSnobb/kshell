module ss_method
  use constant, only: ci, cz, c1, kwf, kdim, max_int4, maxchar
  use model_space
  use class_stopwatch, only: start_stopwatch, stop_stopwatch, &
       time_orth, time_qr_decomp, time_zpares, time_shift_ss, time_restart, &
       time_tmp2
  use wavefunction, only: dot_product_global, type_vec_p, wf_alloc_vec
  use bridge_partitions, only: dealloc_shift_block
  use block_lanczos, only: block_dot_product_real_global, gram_schmidt_qr, &
       dot_product_real_global
#ifdef MPI
  use mpi
#endif
  !$ use omp_lib  
  implicit none

  ! stored for shifted COCG cocg_solve, cocg_shift
  real(8) :: s_pos, s_bb, s_tol
  integer :: s_itr
  integer(kdim) :: s_dim
  real(8), allocatable :: s_alpha(:), s_beta(:), s_rr(:)
  real(kwf), allocatable, target :: s_r(:,:)

  
  private
  public :: zpares_level_dens, zpares_diag

contains

  subroutine zpares_level_dens(matv_blck, &
       emin, emax, nrow_local, n_rhs, estimations, Vmat, maxiter)
    use zpares
    use eigdensity_estimator_mod
    interface 
       subroutine matv_blck(nb, v1, v2)
         ! v2 = mat * v1
         use constant, only: kwf
         integer, intent(in) :: nb
         real(kwf), intent(in)  :: v1(:,:)
         real(kwf), intent(out) :: v2(:,:)
       end subroutine matv_blck
    end interface
    real(8), intent(in) :: emin, emax
    integer, intent(in) :: nrow_local, n_rhs
    real(kwf), intent(inout), allocatable :: Vmat(:,:) 
    real(8), intent(out) :: estimations(:)
    integer, intent(in), optional :: maxiter
    integer :: comm, miter, n_quad, nb, reso, i, n_block
    logical :: is_proj
    real(8) :: dif
    
    miter = 100
    if ( present(maxiter) ) miter = maxiter

    n_block = n_rhs
    n_quad = 16 ! 8
    reso = size(estimations)
    ! seed_val = emin + 30.d0 
    seed_val = min(emax, emin + 30.d0) ! ad hoc
    dif = (emax - emin) / reso
    if (myrank==0) then 
       write(*,'(/,a)') 'eigenvalue count estimation start'
       write(*,'(1x,a,f10.3  )') 'emin   = ', emin
       write(*,'(1x,a,f10.3  )') 'emax   = ', emax
       write(*,'(1x,a,f10.3,/)') 'skip   = ', dif
    end if

    comm = 0
#ifdef MPI 
    comm  = mpi_comm_world
#endif 

    if ( allocated( Vmat ) ) then
       call eigdensity_estimator(comm, matv_blck, n_quad, n_rhs, &
            emin, emax, reso, nrow_local, miter, estimations, Vmat)
    else
       call eigdensity_estimator(comm, matv_blck, n_quad, n_rhs, &
            emin, emax, reso, nrow_local, miter, estimations)     
    end if


    if ( myrank == 0 ) then
       dif = (emax - emin) / reso
       write(*,*)
       write(*,*) 'Parameters ----------------------------------------'

       write(*,'(1x,a,f10.3)') 'emin   = ', emin
       write(*,'(1x,a,f10.3)') 'emax   = ', emax
       write(*,'(1x,a,f10.3)') 'skip   = ', dif
       write(*,'(1x,a,f10.3)') 'E_seed = ', seed_val
       write(*,'(1x,a,I5)') 'Number of sub-intervals =     ', reso
       write(*,'(1x,a,I5)') 'Number of quadrature points = ', n_quad
       write(*,'(1x,a,I5)') 'Number of sample vectors =    ', n_rhs
       write(*,'(1x,a,I5)') 'Max iteration =               ', miter

       write(*,*) 'Results    ----------------------------------------'

       write(*,'(1x,a)') 'Estimation of eigenvalue counts'
       write(*,*) '    i      Energy       Ex.     eigenvalue count'
       do i = 1, reso
          write(*,'(I6, f12.3, f12.3, f18.3)') &
               i, emin+dif*(i-0.5d0), dif*(i-0.5d0), estimations(i)
       end do

       write(*,'( /,a,f18.3,/ )') &
            ' total eigenvalue counts ', sum(estimations)
    end if


    
  end subroutine zpares_level_dens


  subroutine zpares_diag(matv_blck, emin, emax, &
       nrow_local, n_block, num_ev, eval, evec, maxiter, is_proj)
    ! solve eigenvalue problem at  "emin < E < emax"
    use zpares
    interface 
       subroutine matv_blck(nb, v1, v2)
         ! v2 = mat * v1
         use constant, only: kwf
         integer, intent(in) :: nb
         real(kwf), intent(in)  :: v1(:,:)
         real(kwf), intent(out) :: v2(:,:)
       end subroutine matv_blck
    end interface
    real(8), intent(in) :: emin, emax
    integer, intent(in) :: nrow_local, n_block
    integer, intent(out) :: num_ev
    real(8), intent(inout), allocatable :: eval(:)
    real(kwf), intent(inout), allocatable :: evec(:,:)
    integer, intent(in) :: maxiter
    logical, intent(in), optional :: is_proj
    type(zpares_prm) :: prm
    integer :: ncv, left, right, info, i, j, k, n=0
    complex(8) :: z
    real(8), allocatable :: rwork(:,:), eigvec(:,:)
    complex(8), allocatable :: cwork(:,:), b(:), z_array(:), sx(:,:,:)
    real(8), allocatable :: eigval(:),res(:)
    real(kwf), allocatable :: bv(:,:)
    integer :: counter_solve, counter
    real(8) :: norm, x

    call zpares_init(prm)

#ifdef MPI 
    ! parallel in matvec
    prm%high_comm = mpi_comm_self
    prm%low_comm  = mpi_comm_world
#endif 

    prm%L = n_block ! block size
    prm%Lmax = n_block
    !prm%M   ! num of moment degree, def: 16
    ! prm%M = 4 ! 8
    ! prm%extract = ZPARES_EXTRACT_EM
    ! prm%calc_res = .false.
    ! prm%trim_out = .false.
    ! prm%trim_spu = .false.
    ! prm%N = 16
    if (present(is_proj)) prm%get_projection = is_proj
    if (prm%get_projection) then
       prm%calc_res = .false.
!       prm%calc_res = .true.
       prm%M = 1
       prm%trim_out = .true.
!       prm%trim_out = .false.
    end if

    ncv = zpares_get_ncv(prm)
    if (myrank==0) write(*,'(a,i5)') 'SS-method, ncv=', ncv


    allocate( eigval(ncv), res(ncv) )
    allocate( eigvec(nrow_local, ncv), rwork(nrow_local, prm%Lmax), &
         cwork(nrow_local, prm%Lmax) )

    num_ev = 0
    if ( allocated( evec ) ) then 
       num_ev = n_block
       eigvec(:, :num_ev) = evec(:, :num_ev)
       prm%user_source = .true.
    end if

    counter_solve = 1

    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call start_stopwatch(time_zpares)
       call zpares_drcisyev(prm, nrow_local, z, rwork, &
            cwork, emin, emax, num_ev, eigval, eigvec, res, info)
       call stop_stopwatch(time_zpares)

       select case (prm%itask)
       case(ZPARES_TASK_GET_ALL_Z)

          ! get contour integral points

          if ( .not. allocated(z_array) ) then
             counter = 1
             allocate( z_array(prm%num_z_local) )
             allocate( sx(nrow_local, prm%num_z_local, prm%Lmax) )
             call print_memory_usage()
          end if
          z_array(counter) = z
          counter = counter + 1
          
       case(ZPARES_TASK_FACTO)

          ! Here the user factorizes (z*B - A)
          ! At the next return from zpares_zrcihegv,
          ! prm%itask==ZPARES_TASK_SOLVE with the same z

       case(ZPARES_TASK_SOLVE)

          ! i = prm%ws; j = prm%ws+prm%nc-1
          ! Here the user solves (z*B - A) X = cwork(:,i:j) 
          ! The solution X should be stored in cwork(:,i:j)
          ! At the next return from zpares_zrcihegv,
          ! prm%itask==ZPARES_TASK_FACTO_H with the same z is returned

          ! do i = prm%ws, prm%ws+prm%nc-1
          !    b = cwork(:,i)
          !    call cocg_complex_solve(z, b, cwork(:,i))
          ! end do

          if (counter_solve == 1) then
             ! Store right hand side vector
             ! rhs(:) = cwork(:,prm%ws+i-1) 
             ! Solve shifted linear systems
             ! (z_j I - A) x_j = rhs
             ! with shifted_cocg method
             ! z_j (j=1,...,prm%num_z_local) is stored in z_array(j)
             !   call shifted_cocg(rhs, z_array, sx(:,:,i))
             ! solution x_j is stored in sx(:,j,i)
             !
             allocate( bv(nrow_local, prm%nc) )

             bv(:, :) = cwork(:, prm%ws : prm%ws+prm%nc-1) 

             call block_cocg_solve_shift(matv_blck, &
                  z_array, emax, bv, sx, maxiter=maxiter)

             deallocate( bv )
          end if
          ! Store solution in cwork
          cwork(:,prm%ws:prm%ws+prm%nc-1) = sx(:,counter_solve,1:prm%nc)

          counter_solve = counter_solve + 1         

       case(ZPARES_TASK_MULT_A)
          ! iw = prm%ws; jw = prm%ws+prm%nc-1
          ! ix = prm%xs; jx = prm%xs+prm%nc-1
          ! Here the user performs matrix-vector multiplications:
          ! rwork(:,iw:jw) = A*X(:,ix:jx)

          call matv_blck( prm%nc, &
               eigvec( : , prm%xs:prm%xs+prm%nc-1 ), &
               rwork( : , prm%ws : prm%ws+prm%nc-1 ) )

       end select
    end do

    if (prm%get_projection .and. num_ev /= n_block) &
         write(*,*) 'ERROR in projection', num_ev, n_block


    call dealloc_shift_block()
    if (num_ev /= size(evec, 2)) then
       if ( allocated( evec ) ) deallocate( evec )
       allocate( evec(nrow_local, num_ev) ) 
    end if

    evec(:, :num_ev) = eigvec(:, :num_ev)

    if (.not. prm%get_projection) then
       if (myrank==0) then 
          write(*,*)
          write(*,'(a,i5,a,f8.3,a,f8.3)') '# of eigenvectors', num_ev, &
               ' between ',emin,' and ',emax
          do i = 1, num_ev
             write(*,*) 'eigval,res=', eigval(i), res(i)
          end do
          write(*,*)
       end if
    end if

    if (allocated(eval)) deallocate(eval)
    allocate( eval(num_ev) )
    eval = eigval(:num_ev)

    call zpares_finalize( prm )

    deallocate( eigval, res )
    deallocate( eigvec, rwork, cwork )


  contains

    subroutine print_memory_usage()
      if (myrank /= 0) return
      write(*,'(/,a,f10.3,a,f10.3)') &
           'SS eigenvalue range  ', emin,'  < E < ', emax
      write(*,'(a,i3,a,i3,a,i6)') 'Block size = ', n_block, &
           '    Mesh points = ', prm%N, '  ncv = ', ncv
      write(*,'(a,i3)') ' num_z_local',prm%num_z_local
      x = dble(nrow_local) * (ncv*8.d0 + prm%Lmax*8.d0 + prm%Lmax*16.d0 &
           + prm%num_z_local*prm%Lmax*16.d0) * 1.d-9
      write(*,'(a,f8.3,a,/)') 'Memory size in SS: ', x, 'GB'
    end subroutine print_memory_usage

  end subroutine zpares_diag




  subroutine block_cocg_solve_shift(matv_blck, &
       z_array, rz, b, sx, maxiter, tol)
    ! TODO : b => qk and destructive to save memory
    !
    ! Shifted Block COCG method
    !
    !        solve (rz*I - H) * x_i = b_i
    !    and return sx of (z_j*I - H) * sx_(j,i) = b_i
    !
    ! rz, b in real   
    ! z, x  in complex
    ! 
    ! Ref. Y. Futamura et al., VECPAR 2012, Vol.7851 226-235, 2013.
    !    http://dx.doi.org/10.1007/978-3-642-38718-0_23
    ! 
    use lib_matrix, only: gen_eye, inverse
    interface 
       subroutine matv_blck(nb, v1, v2)
         ! v2 = mat * v1
         use constant, only: kwf
         integer, intent(in) :: nb
         real(kwf), intent(in)  :: v1(:,:)
         real(kwf), intent(out) :: v2(:,:)
       end subroutine matv_blck
    end interface
    real(8), intent(in) :: rz
    complex(8), intent(in) :: z_array(:)
    real(kwf), intent(in) :: b(:,:)
    complex(kwf), intent(out) :: sx(:,:,:)
    integer, intent(in), optional :: maxiter
    real(8), intent(in), optional :: tol
    integer :: nb, i, j, k, miter, nz
!    integer(kdim) :: ld
    integer :: ld
    real(8), allocatable :: xk(:,:), pk(:,:), qk(:,:), apk(:,:)
    real(8), allocatable :: alpha(:,:,:), beta(:,:,:), &
         rho(:,:,:), delta(:,:,:), t(:,:)
    real(8), allocatable :: eye(:,:)
    complex(8), allocatable :: x_sj(:,:,:), p_sj(:,:,:), t_sj(:,:,:)
    complex(8), allocatable :: a_sj(:,:,:), b_sj(:,:,:), xi_sj(:,:,:,:)
    real(8) :: s_bb
    logical :: is_conv 

    nb = size(b, 2) ! block size
    ld = size(b, 1)

    miter = 1000
    if (present(maxiter)) miter = maxiter
    s_tol = 1.d-7
    if (present(tol)) s_tol = tol
    nz = size(z_array)

    allocate( xk(ld, nb),  pk(ld, nb), qk(ld, nb), apk(ld, nb) )
    allocate( alpha(nb, nb, -1:miter), beta(nb, nb, -1:miter), &
         rho(nb, nb, 0:miter), delta(nb, nb, 0:miter), t(nb,nb) )
    allocate( x_sj(ld, nb, nz), p_sj(ld, nb, nz), t_sj(ld, nb, nz) )
    allocate( a_sj(nb, nb, nz), b_sj(nb, nb, nz), xi_sj(nb, nb, nz, -1:miter) )
    allocate( eye(nb,nb) )

    eye = gen_eye(nb)
    xk = 0.d0
    alpha(:,:,-1) = eye
    !$omp parallel do
    do i = 1, nb
       qk(:,i) = b(:,i)
    end do
    call gram_schmidt_qr(qk, rho(:,:,0))
    delta(:,:,0) = rho(:,:,0)
    s_bb = sqrt(sum(delta(:,:,0)**2))

    pk = qk
    !$omp parallel do
    do j = 1, nz
       xi_sj(:,:,j,-1) = eye
       xi_sj(:,:,j, 0) = rho(:,:,0)
       p_sj(:,:,j) = qk
    end do


    do k = 0, miter-1

       call matvec_block_real(pk, apk)
       apk = rz * pk - apk

!       alpha(:,:,k) = inverse( matmul( transpose(pk), apk ) )
       call block_dot_product_real_global(pk, apk, alpha(:,:,k))
       alpha(:,:,k) = inverse( alpha(:,:,k) )

       ! xk = xk + matmul(pk, matmul(alpha(:,:,k), delta(:,:,k)))
       call dgemm('n', 'n', nb, nb, nb, 1.d0, &
            alpha(:,:,k), nb, delta(:,:,k), nb, 0.d0, t, nb)
       call dgemm('n', 'n', ld, nb, nb, 1.d0, &
            pk, ld, t, nb, 1.d0, xk, ld)
       
       ! qk = qk - matmul(apk, alpha(:,:,k))
       call dgemm('n', 'n', ld, nb, nb, -1.d0, &
            apk, ld, alpha(:,:,k), nb, 1.d0, qk, ld)

       call gram_schmidt_qr( qk, rho(:,:,k+1) )

       ! delta(:,:,k+1) = matmul(rho(:,:,k+1), delta(:,:,k) )
       call dgemm('n', 'n', nb, nb, nb, 1.d0, &
            rho(:,:,k+1), nb, delta(:,:,k), nb, 0.d0, &
            delta(:,:,k+1), nb)

       ! pk = qk + matmul( pk, transpose(rho(:,:,k+1)) )
       call dcopy(ld*nb, qk, 1, apk, 1)
       call dgemm('n', 't', ld, nb, nb, 1.d0, &
            pk, ld, rho(:,:,k+1), nb, 1.d0, apk, ld)
       call dcopy(ld*nb, apk, 1, pk, 1)

       call shift()
       ! call shift_pure()

       call check_conv( delta(:,:,k+1), is_conv )
       if (is_conv) exit

    end do


    if (.not. is_conv) then
       s_itr = miter - 1 
       if (myrank==0)  write(*,*) 'COCG NOT converged', s_itr
    end if

    do i = 1, nb 
       sx(:,:,i) = x_sj(:,i,:)
    end do

    deallocate( xk, pk, qk, apk)
    deallocate( alpha, beta, rho, delta, t )
    deallocate( x_sj,  p_sj, t_sj)
    deallocate( a_sj,  b_sj, xi_sj )


  contains


    subroutine check_conv(dlt, is_conv)
      real(8), intent(in) :: dlt(:,:) 
      logical, intent(out) :: is_conv
      real(8) :: res
      res = sqrt(sum(dlt(:,:)**2)) / s_bb
      if (myrank==0) write(*,'(a,i5,100f18.10)') 'cocg res',k,res       
      if ( res < s_tol ) then
         if (myrank==0) write(*,'(a,i5,f15.8)') 'COCG converged', k, res
         is_conv = .true.
         return
      end if
      is_conv = .false.
    end subroutine check_conv


    subroutine shift()
      !$ use omp_lib
      integer :: j
      complex(8) :: shft
      real(8) :: at(nb, nb), bt(nb, nb), inv_a_kp(nb, nb)
      integer :: nt, nn, cnk, n

      call start_stopwatch(time_shift_ss)

      at = matmul( alpha(:,:,k), inverse( rho(:,:,k+1) ) )
      bt = matmul( inverse( alpha(:,:,k) ), transpose(rho(:,:,k+1)) )
      inv_a_kp =  inverse( alpha(:,:,k-1) )



      !$omp parallel do private(j, shft)
      do j = 1, nz
         shft = z_array(j) - rz

         xi_sj(:,:,j,k+1) =  &
              inverse( eye + shft * alpha(:,:,k) &
              +  matmul( rho(:,:,k) &
              &     - matmul( xi_sj(:,:,j,k), inverse( xi_sj(:,:,j,k-1) ) ), &
              &  matmul( inv_a_kp, &
              &  matmul( transpose( rho(:,:,k) ), &
              &          alpha(:,:,k) ) ) ) )
         xi_sj(:,:,j,k+1) = matmul( rho(:,:,k+1),     xi_sj(:,:,j,k+1) )
         xi_sj(:,:,j,k+1) = matmul( xi_sj(:,:,j,k+1), xi_sj(:,:,j,k) )

         a_sj(:,:,j)  = matmul( at, xi_sj(:,:,j,k+1) )
         b_sj(:,:,j)  = matmul( a_sj(:,:,j), &
              &         matmul( inverse( xi_sj(:,:,j,k) ), &
              &                 bt ) )
      end do


      nt = 1 
      !$ nt = omp_get_max_threads()
 
      ! if (.false.) then
      ! if (.true.) then
      if ( mod(nz, nt) == 0 ) then

         !$omp parallel do 
         do j = 1, nz
            ! x_sj(:,:,j) = x_sj(:,:,j) + matmul( p_sj(:,:,j), a_sj(:,:,j) )
            call zgemm('n', 'n', ld, nb, nb, c1, &
                 p_sj(:,:,j), ld, a_sj(:,:,j), nb, c1, x_sj(:,:,j), ld)
            ! p_sj(:,:,j) = qk          + matmul( p_sj(:,:,j), b_sj(:,:,j) )
            call d_z_copy(ld*nb, qk, t_sj(:,:,j))  ! 2.7sec
            call zgemm('n', 'n', ld, nb, nb, c1, &
                 p_sj(:,:,j), ld, b_sj(:,:,j), nb, c1, t_sj(:,:,j), ld)
            call zcopy(ld*nb, t_sj(:,:,j), 1, p_sj(:,:,j), 1)
         end do

      else

         ! tuned for any number of threads
         cnk = (ld - 1)/nt + 1
         do j = 1, nz
            !$omp parallel do private(i, n)
            do i = 0, nt-1
               n = min(cnk,  ld - i*cnk)
               if (n <= 0) cycle
               ! x_sj(:,:,j) = x_sj(:,:,j) + matmul( p_sj(:,:,j), a_sj(:,:,j) )
               call zgemm('n', 'n', n, nb, nb, c1, &
                    p_sj(i*cnk+1,1,j), ld, a_sj(:,:,j), nb, c1, x_sj(i*cnk+1,1,j), ld)
               ! p_sj(:,:,j) = qk          + matmul( p_sj(:,:,j), b_sj(:,:,j) )
               t_sj(i*cnk+1:i*cnk+n,:,j) = qk(i*cnk+1:i*cnk+n,:)
               call zgemm('n', 'n', n, nb, nb, c1, &
                    p_sj(i*cnk+1,1,j), ld, b_sj(:,:,j), nb, c1, t_sj(i*cnk+1,1,j), ld)
               p_sj(i*cnk+1:i*cnk+n,:,j) = t_sj(i*cnk+1:i*cnk+n,:,j)
            end do
         end do

      endif

      call stop_stopwatch(time_shift_ss)

    end subroutine shift

    subroutine inside_shift()


    end subroutine inside_shift


    subroutine shift_pure()
      ! shift before tuning
      integer :: j
      complex(8) :: shft

      call start_stopwatch(time_shift_ss)

      do j = 1, nz
         shft = z_array(j) - rz

         xi_sj(:,:,j,k+1) =  &
              inverse( eye + shft * alpha(:,:,k) &
              +  matmul( rho(:,:,k) &
              &     - matmul( xi_sj(:,:,j,k), inverse( xi_sj(:,:,j,k-1) ) ), &
              &  matmul( inverse( alpha(:,:,k-1) ), &
              &  matmul( transpose( rho(:,:,k) ), &
              &          alpha(:,:,k) ) ) ) )
         xi_sj(:,:,j,k+1) = matmul( rho(:,:,k+1),     xi_sj(:,:,j,k+1) )
         xi_sj(:,:,j,k+1) = matmul( xi_sj(:,:,j,k+1), xi_sj(:,:,j,k) )

         a_sj(:,:,j) &
              = matmul( alpha(:,:,k), &
              & matmul( inverse( rho(:,:,k+1) ), &
              &         xi_sj(:,:,j,k+1) ) )
         b_sj(:,:,j) &
              = matmul( a_sj(:,:,j), &
              & matmul( inverse( xi_sj(:,:,j,k) ), &
              & matmul( inverse( alpha(:,:,k) ), &
              &         transpose( rho(:,:,k+1) ) ) ) )

         x_sj(:,:,j) = x_sj(:,:,j) + matmul( p_sj(:,:,j), a_sj(:,:,j) )
         p_sj(:,:,j) = qk          + matmul( p_sj(:,:,j), b_sj(:,:,j) )
      end do

      call stop_stopwatch(time_shift_ss)

    end subroutine shift_pure


    subroutine matvec_block_real(p, q)
      ! x = H*p
      real(8), intent(in) :: p(:,:)
      real(8), intent(out) :: q(:,:)
      integer :: i

      call matv_blck(size(p,2), p, q)

    end subroutine matvec_block_real

    
    subroutine d_z_copy(n, d, z)
      integer, intent(in) :: n
      real(8), intent(in) :: d(*)
      complex(8), intent(out) :: z(*)
      integer :: i
      !omp parallel do 
      do i = 1, n
         z(i) = d(i)
      end do
    end subroutine d_z_copy


  end subroutine block_cocg_solve_shift



end module ss_method

