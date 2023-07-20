module eigdensity_estimator_mod
  use class_stopwatch

  double precision :: seed_val = 0d0

contains  

  subroutine eigdensity_estimator(comm, matvec, n_quad, n_rhs, emin, emax &
       , resolution, nrow_local, imax, estimations, source_matrix)
    !
    ! written by Y. Futamura, Tsukuba 2015/04/20
    !
    ! NOTE: source_matrix is destructive to suppress memory usage
    !     major memory capacity ... nrow_local*n_rhs*4*8 bytes
    ! w/o MPI, realized by N. Shimizu, CNS-Tokyo 2015/04/29
    ! 
    implicit none
    interface
       subroutine matvec(n_vec, source, destination)
         integer, intent(in) :: n_vec
         double precision, intent(in) :: source(:,:)
         double precision, intent(out) :: destination(:,:)
       end subroutine matvec
    end interface
    integer, intent(in) :: comm, n_quad, n_rhs, imax, resolution, nrow_local
    double precision, intent(in) :: emin, emax
    double precision, intent(out) :: estimations(resolution)
    double precision, optional, intent(inout) :: source_matrix(:,:)

    double precision, allocatable :: centers(:), radii(:), residual(:,:)
    integer :: i, j, k, rank, iter, ierr
    complex(kind(0d0)) :: trace
    complex(kind(0d0)), allocatable :: omegas(:), thetas(:), VX(:,:,:)
    double precision, allocatable :: c_source_matrix(:,:)
    
    rank = 0
#ifdef MPI
    call MPI_COMM_RANK(comm, rank, ierr)
#endif
    allocate(centers(resolution), radii(resolution) )
    allocate(omegas(n_quad*resolution), thetas(n_quad*resolution))
    allocate(VX(n_rhs,n_rhs,n_quad*resolution), residual(n_rhs,n_quad*resolution))
    
    radii(1) = (emax - emin)/(resolution*2d0)
    centers(1) = emin + radii(1)
    do i = 2, resolution
       centers(i) = centers(i-1) + 2*radii(1)
       radii(i) = radii(i-1)
    end do
    
    do i = 1, resolution
       call calc_omega(centers(i),radii(i),n_quad,thetas((i-1)*n_quad+1:i*n_quad),omegas((i-1)*n_quad+1:i*n_quad)) 
    end do

    if (present(source_matrix)) then
       call shifted_block_CG_rQ_bilinear &
            (comm, matvec, source_matrix, omegas, 1d-4, 1d-4, imax, &
            VX, iter, residual, resolution, n_quad, thetas, radii)
    else
       allocate( c_source_matrix(nrow_local,n_rhs) )
       call create_hutch_samples(c_source_matrix, nrow_local, n_rhs, rank)
       call shifted_block_CG_rQ_bilinear &
            (comm, matvec, c_source_matrix, omegas, 1d-4, 1d-4, imax, &
            VX, iter, residual, resolution, n_quad, thetas, radii)
       deallocate( c_source_matrix )
    end if

    !$omp parallel do private(i, j, k, trace)
    do i = 1,resolution
       estimations(i) = 0d0
       do j = 1, n_quad
          trace = (0d0, 0d0)
          do k = 1, n_rhs
            trace = trace + VX(k,k,(i-1)*n_quad+j)
          end do
          estimations(i) = estimations(i) &
               + thetas((i-1)*n_quad+j)*trace/dcmplx(dble(n_rhs),0d0)
       end do
       estimations(i) = radii(i)*estimations(i)/n_quad
    end do

    return
    
  end subroutine eigdensity_estimator
  
  subroutine calc_omega(gamma,rho,N,theta,omega)
    implicit none
    real(8),intent(in) :: gamma 
    real(8),intent(in) :: rho
    integer,intent(in) :: N

    complex(8),intent(out) :: theta(:),omega(:)
    
    integer :: j
    complex(8) :: i
    real(8) :: PI
    intrinsic dATAN,dcos,dsin
    PI = 4.0D0*dATAN(1.0D0)

    i = (0.0D0,1.0D0)
    do j = 1,N
       !theta(j) = exp(2.0D0 * PI * i * ((j - 1) + 0.5D0) / N)
       theta(j) = dcmplx(dcos(2.0D0 * PI * ((j - 1) + 0.5D0) / N),dsin(2.0D0 * PI * ((j - 1) + 0.5D0) / N))
       omega(j) = gamma + rho * theta(j)
    end do
  end subroutine calc_omega

  subroutine projection_zpares(comm, nrow_local, n_rhs, matvec, emin, emax, imax, X, n_rhs_part_in)
    use zpares
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    interface
       subroutine matvec(n_vec, source, destination)
         integer, intent(in) :: n_vec
         double precision, intent(in) :: source(:,:)
         double precision, intent(out) :: destination(:,:)
       end subroutine matvec
    end interface
    integer, intent(in) :: comm, nrow_local, n_rhs, imax
    double precision, intent(in) :: emin, emax
    double precision, intent(inout) :: X(:,:)
    integer, optional, intent(in) :: n_rhs_part_in
    
    type(zpares_prm) :: prm
    integer :: info, num_ev, ncv, counter, counter_solve, n_rhs_part, linsol_n_rhs, i, j
    double precision, allocatable :: res(:), rwork(:,:), eigval(:)
    complex(kind(0d0)) :: z
    complex(kind(0d0)), allocatable :: cwork(:,:), z_array(:), sx(:,:,:), rhs(:,:)    

    call start_stopwatch(time_proj_ss)
    
    if ( present(n_rhs_part_in) ) then
       n_rhs_part = n_rhs_part_in
    else
       n_rhs_part = n_rhs
    end if

    if (size(X,1) /= nrow_local) stop 'ERROR projection_zpares'

    call zpares_init(prm)
    prm%L = n_rhs
    prm%Lmax = n_rhs
    prm%N = 32
    ! prm%N = 16
    prm%M = 1

#ifdef MPI
    prm%high_comm = MPI_COMM_SELF
    prm%low_comm = comm
#endif 

    ncv = zpares_get_ncv(prm) ! ncv must be prm%L since prm%M==1
    allocate(eigval(ncv), res(ncv))
    allocate(rwork(nrow_local, prm%Lmax), cwork(nrow_local, prm%Lmax))
    num_ev = n_rhs
    prm%user_source = .true.
    prm%get_projection = .true.
    
    counter_solve = 1
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call zpares_drcisyev &
            (prm, nrow_local, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info)
       select case (prm%itask)
       case(ZPARES_TASK_GET_ALL_Z)
          if ( .not. allocated(z_array) ) then
             counter = 1
             allocate(z_array(prm%num_z_local))
             allocate(sx(nrow_local, prm%Lmax, prm%num_z_local))
          end if
          z_array(counter) = z
          counter = counter + 1
          
       case(ZPARES_TASK_FACTO)

          ! Do nothing
          
       case(ZPARES_TASK_SOLVE)
          
          if (counter_solve == 1) then
             allocate(rhs(nrow_local, prm%nc))
             rhs(:,1:prm%nc) = cwork(:,prm%ws:prm%ws+prm%nc-1)
             z_array(:) = -z_array(:)

             do i = 1, n_rhs, n_rhs_part

                call stop_stopwatch(time_proj_ss)

                j = min(i + n_rhs_part - 1, n_rhs)
                call shifted_block_CG_rQ_naive(comm, matvec, rhs(:,i:j), z_array, 1d-4, 1d-4, imax, sx(:,i:j,:))
!                call shifted_block_CG_rQ_unroll_iter(comm, matvec, rhs(:,i:j), z_array, 1d-4, 1d-4, imax, sx(:,i:j,:))

                call start_stopwatch(time_proj_ss)

             end do

             sx(:,:,:) = -sx(:,:,:)
          end if
          cwork(:,prm%ws:prm%ws+prm%nc-1) = sx(:,1:prm%nc,counter_solve)

          counter_solve = counter_solve + 1         

       case(ZPARES_TASK_MULT_A)
                    
          call matvec(prm%nc, X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))
          ! but this will not be called

       end select
    end do

    call stop_stopwatch(time_proj_ss)
    
  end subroutine projection_zpares

  subroutine shifted_block_CG_rQ_unroll_iter(comm,matvec,B,sigma,epsmax,epsmin,imax,X)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    interface
       subroutine matvec(n_vec, source, destination)
         integer, intent(in) :: n_vec
         double precision, intent(in) :: source(:,:)
         double precision, intent(out) :: destination(:,:)
       end subroutine matvec
    end interface

    integer, parameter :: unroll_step = 8

    integer, intent(in) :: comm, imax
    real(8), intent(in) :: epsmax,epsmin
    complex(8), intent(in) :: B(:,:), sigma(:)
    complex(8), intent(out) :: X(:,:,:)

    integer :: i, j, k, n, L, m, iter, ierr, maxidx, conv_count, myrank, nprocs
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp
    real(8), allocatable :: res_seed(:), Bnorm_inv(:), residual(:,:), r_P_seed(:,:), r_AP(:,:) &
         ,r_tmp(:,:), r_tmp2(:,:)
    
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    complex(8), allocatable :: &
         X_seed(:,:), Q(:,:,:), P(:,:,:), P_seed(:,:), ts_tmp(:,:) &
         , AP(:,:), PAP(:,:), PAP_old(:,:), alpha(:,:) &
         , delta(:,:), rho(:,:), rho_old(:,:), eye(:,:), blk_tmp1(:,:), blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:) &
         , XP_coef(:,:,:,:), ts_tmp_dble(:,:), blk_tmp_dble1(:,:), blk_tmp_dble2(:,:)
    logical, allocatable :: conv_flag(:)

    myrank = 0
    nprocs = 1
#ifdef MPI
    call MPI_COMM_RANK(comm, myrank, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)
#endif

    n = size(X,1)
    L = size(X,2)
    m = size(X,3)

    allocate(X_seed(n,L),Q(n,L,unroll_step+1),P_seed(n,L),P(n,L,m),AP(n,L), &
         ts_tmp(n,L),PAP(L,L),PAP_old(L,L),alpha(L,L) &
         ,delta(L,L),rho(L,L),rho_old(L,L),eye(L,L),blk_tmp1(L,L),blk_tmp2(L,L) &
         ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
         ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L), iter_array(m) &
         ,XP_coef(L,2*L,unroll_step+1,m),ts_tmp_dble(n,2*L) &
         ,blk_tmp_dble1(L,2*L),blk_tmp_dble2(L,2*L),conv_flag(m),residual(L,m) &
         ,r_P_seed(n,L), r_AP(n,L))

    X_seed = (0d0,0d0)
    X = (0d0,0d0)
    eye = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do
    Q(:,:,unroll_step+1) = B
    call z_MGS_QR(comm,Q(:,:,unroll_step+1),delta)
    
    rho = delta
    P_seed = Q(:,:,unroll_step+1)
    
    do k=1,m
       P(:,:,k) = Q(:,:,unroll_step+1)
       xi1(:,:,k) = rho
       xi2(:,:,k) = eye
    end do
    xitld_hat = (0d0,0d0)
    call norm2_as_block(comm,B,Bnorm_inv,n,L)
    Bnorm_inv = 1d0 / Bnorm_inv
    PAP = (0d0,0d0)
    residual = 1d0
    conv_flag = .false.
    iter_array = unroll_step + 1

    do iter=1,imax
       r_P_seed = P_seed
       call matvec(L,r_P_seed,r_AP)
       AP = r_AP

       ! !$OMP parallel do
       AP = AP + seed_val*P_seed
       ! !$OMP end parallel do
       PAP_old = PAP
       call z_matmat_CN(comm,P_seed,AP,PAP)
       call inv(PAP,alpha,ipiv,work)
       call z_matmat_NX('N',zero,alpha,delta,blk_tmp1)
       call z_matmat_NX('N',one,P_seed,blk_tmp1,X_seed)

       rho_old = rho
       
       ! !$OMP parallel do
       Q(:,:,modulo(iter-1,unroll_step+1)+1) = Q(:,:,modulo(iter-2,unroll_step+1)+1)
       ! !$OMP end parallel do

       call z_matmat_NX('N',one,AP,-alpha,Q(:,:,modulo(iter-1,unroll_step+1)+1))

       call z_MGS_QR(comm,Q(:,:,mod(iter-1,unroll_step+1)+1),rho)

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta = blk_tmp1

       call z_matmat_NX('C',zero,P_seed,rho,ts_tmp)
       ! !$OMP parallel do
       P_seed = Q(:,:,mod(iter-1,unroll_step+1)+1) + ts_tmp
       ! !$OMP end parallel do

       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed = res_seed*Bnorm_inv

       do k=1,m
          if ( .not. conv_flag(k) ) then
             xi2(:,:,k) = xi1(:,:,k)
             blk_tmp1 = eye - xitld_hat(:,:,k)
             call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
             call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)
             call inv(eye - (seed_val - sigma(k))*alpha + blk_tmp1,xitld_hat(:,:,k),ipiv,work)
             call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
             call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
             call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)
             XP_coef(:,1:L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1
             
             call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
             call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
             call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)
             XP_coef(:,L+1:2*L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1
             
             call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
             residual(:,k) = residual(:,k)*Bnorm_inv
          end if
       end do
       
       

       if ( mod(iter,unroll_step+1) == 0 ) then
          do k=1,m
             if ( .not. conv_flag(k) ) then
                do j=unroll_step,1,-1
                   blk_tmp_dble2 = XP_coef(:,:,j+1,k)
                   call z_matmat_NX('N',zero,XP_coef(:,L+1:2*L,j,k),blk_tmp_dble2,blk_tmp_dble1)
                   XP_coef(:,1:L,j,k) = XP_coef(:,1:L,j,k) + blk_tmp_dble1(:,1:L)
                   XP_coef(:,L+1:2*L,j,k) = blk_tmp_dble1(:,L+1:2*L)
                end do
             end if
          end do

          do k=1,m
             if ( .not. conv_flag(k) ) then
                ! !$OMP parallel do
                ts_tmp_dble(:,1:L) = X(:,:,k)
                ! !$OMP end parallel do
!                 call zcopy(L*n,X(:,:,k),1,ts_tmp_dble(:,1:L),1)
                ! !$OMP parallel do
                ts_tmp_dble(:,L+1:2*L) = Q(:,:,unroll_step+1)
                ! !$OMP end parallel do
!                 call zcopy(L*n,Q(:,:,unroll_step+1),1,ts_tmp_dble(:,L+1:2*L),1)
                
                do j=unroll_step,1,-1
                   call z_matmat_NX('N',one,Q(:,:,j),XP_coef(:,:,j+1,k),ts_tmp_dble)
                end do
                call z_matmat_NX('N',one,P(:,:,k),XP_coef(:,:,1,k),ts_tmp_dble)
                   
                ! !$OMP parallel do
                X(:,:,k) = ts_tmp_dble(:,1:L)
                ! !$OMP end parallel do
!                 call zcopy(L*n,ts_tmp_dble(:,1:L),1,X(:,:,k),1)
                ! !$OMP parallel do
                P(:,:,k) = ts_tmp_dble(:,L+1:2*L)
                ! !$OMP end parallel do
!                 call zcopy(L*n,ts_tmp_dble(:,L+1:2*L),1,P(:,:,k),1)
             end if
          end do
             
          ! conv_count = 0
          ! do k=1,m
          !    conv_flag(k) = maxval(residual(:,k)) < epsmin
          !    if ( maxval(residual(:,k)) > epsmin ) then
          !       iter_array(k) = iter + unroll_step + 1
          !    end if
             
          !    if ( iter_array(k) < iter ) then
          !       conv_count = conv_count + 1
          !    end if
          ! end do

          ! if ( maxval(residual) < epsmax .or. conv_count > m/2 ) then
          if ( maxval(residual) < epsmax ) then
             exit
          end if          
       end if
    end do
    
    ! compute true residual
    ! allocate(r_tmp(n,2), r_tmp2(n,2))
    ! do i=1,m       
    !    maxidx = maxloc(residual(:,i),1)
    !    r_tmp(:,1) = real(X(:,maxidx,i),kind(0d0))
    !    r_tmp(:,2) = aimag(X(:,maxidx,i))
    !    call matvec(2,r_tmp(:,1:2),r_tmp2(:,1:2))
    !    AP(:,maxidx) = cmplx(r_tmp2(:,1),r_tmp2(:,2),kind(0d0))
    !    AP(:,maxidx) = B(:,maxidx) - (AP(:,maxidx) + sigma(i)*X(:,maxidx,i))       

    !    res_true = dsqrt(dble(z_dot(comm,AP(:,maxidx),AP(:,maxidx))))
    !    res_true = res_true*Bnorm_inv(maxidx)
    !    write(*,*) i, residual(maxidx,i), res_true
    ! end do

  end subroutine shifted_block_CG_rQ_unroll_iter

  subroutine shifted_block_CG_rQ_naive(comm,matvec,B,sigma,epsmax,epsmin,imax,X)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    interface
       subroutine matvec(n_vec, source, destination)
         integer, intent(in) :: n_vec
         double precision, intent(in) :: source(:,:)
         double precision, intent(out) :: destination(:,:)
       end subroutine matvec
    end interface

    integer, intent(in) :: imax,comm
    real(8), intent(in) :: epsmax,epsmin
    complex(8), intent(inout) :: B(:,:) ! destructive
    complex(8) :: sigma(:)
    complex(8), intent(out) :: X(:,:,:)

    integer :: i, j, k, n, L, m, ierr, maxidx, conv_count, iter
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp, resm
    real(8), allocatable :: res_seed(:), Bnorm_inv(:), residual(:,:), r_P_seed(:,:), r_AP(:,:) &
         , r_tmp(:,:), r_tmp2(:,:)
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    real(8), allocatable :: P_seed(:,:), AP(:,:), Q(:,:), ts_tmp(:,:)
    complex(8), allocatable :: &
         X_seed(:,:), P(:,:,:) &
         , PAP(:,:), PAP_old(:,:), alpha(:,:) &
         , delta(:,:), rho(:,:), rho_old(:,:), eye(:,:), blk_tmp1(:,:), blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:)
    real(8), allocatable :: r_PAP(:,:), r_rho(:,:), r_alpha(:,:), r_delta(:,:)
    integer :: myrank, nprocs

    myrank = 0
    nprocs = 1
#ifdef MPI
    call MPI_COMM_RANK(comm, myrank, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)
#endif

    
    n = size(X,1)
    L = size(X,2)
    m = size(X,3)

    allocate(X_seed(n,L),P_seed(n,L),P(n,L,m),AP(n,L) &
         ,Q(n,L),ts_tmp(n,L),PAP(L,L),PAP_old(L,L),alpha(L,L) &
         ,delta(L,L),rho(L,L),rho_old(L,L),eye(L,L),blk_tmp1(L,L),blk_tmp2(L,L) &
         ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
         ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L),iter_array(m),residual(L,m) &
         ,r_P_seed(n,L), r_AP(n,L))
    allocate( r_PAP(L,L), r_rho(L,L), r_alpha(L,L), r_delta(L,L) )

    X_seed = (0d0,0d0)
    X = (0d0,0d0)
    eye = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do
    Q = real(B,kind(0d0))
       
    call d_norm2_as_block(comm,Q,Bnorm_inv,n,L)
    Bnorm_inv = 1d0 / Bnorm_inv
    
    call d_MGS_QR(comm,Q,r_delta)
    delta = r_delta
    
    rho = delta
    P_seed = Q
    
    do k=1,m
       P(:,:,k) = Q
       xi1(:,:,k) = rho
       xi2(:,:,k) = eye
    end do
    xitld_hat = (0d0,0d0)
    PAP = (0d0,0d0)
    residual = 1d0
    
    do iter=1,imax

       call matvec(L,P_seed,AP)

!       AP = seed_val*P_seed - AP
       AP = seed_val*P_seed + AP

       PAP_old = PAP
       call d_matmat_CN(comm,P_seed,AP,r_PAP)
       PAP = r_PAP
       call inv(PAP,alpha,ipiv,work)

       rho_old = rho
       r_alpha = alpha
              
       call d_matmat_NX('N',1d0,AP,-r_alpha,Q)

       call d_MGS_QR(comm,Q,r_rho)
       rho = r_rho

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta = blk_tmp1

       call d_matmat_NX('T',0d0,P_seed,r_rho,ts_tmp)
       P_seed = Q + ts_tmp
       
       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed = res_seed*Bnorm_inv
       
       call start_stopwatch(time_shift_ss)
       do k=1,m
          if ( maxval(residual(:,k)) > epsmin ) then             
             xi2(:,:,k) = xi1(:,:,k)
             blk_tmp1 = eye - xitld_hat(:,:,k)
             call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
             call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)
             call inv(eye - (seed_val - sigma(k))*alpha + blk_tmp1, &
                  xitld_hat(:,:,k), ipiv, work)
             call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
             call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
             call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)

             call z_matmat_NX('N',one,P(:,:,k),blk_tmp1,X(:,:,k))

             call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
             call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
             call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)

             call z_matmat_NX('N',zero,P(:,:,k),blk_tmp1,B)
             P(:,:,k) = Q + B

             call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
             residual(:,k) = residual(:,k)*Bnorm_inv
          end if
          if ( maxval(residual(:,k)) > epsmin ) then
             iter_array(k) = iter
          end if
       end do
       call stop_stopwatch(time_shift_ss)

       resm = maxval(residual)
       if (myrank==0) write(*,'(a,i5,f15.8)') 'Proj. COCG res', &
            iter,resm
       if ( resm < epsmax ) exit

    end do

    ! compute true residual
    ! allocate(r_tmp(n,2), r_tmp2(n,2))
    ! do i=1,m       
    !    maxidx = maxloc(residual(:,i),1)
    !    r_tmp(:,1) = real(X(:,maxidx,i),kind(0d0))
    !    r_tmp(:,2) = aimag(X(:,maxidx,i))
    !    call matvec(2,r_tmp(:,1:2),r_tmp2(:,1:2))
    !    AP(:,maxidx) = cmplx(r_tmp2(:,1),r_tmp2(:,2),kind(0d0))
    !    AP(:,maxidx) = B(:,maxidx) - (AP(:,maxidx) + sigma(i)*X(:,maxidx,i))       

    !    res_true = dsqrt(dble(z_dot(comm,AP(:,maxidx),AP(:,maxidx))))
    !    res_true = res_true*Bnorm_inv(maxidx)
    !    write(*,*) i, residual(maxidx,i), res_true
    ! end do

  end subroutine shifted_block_CG_rQ_naive


  subroutine shifted_block_CG_rQ_bilinear(comm,matvec,Q,sigma,epsmax,epsmin,imax,X,iter,residual,nc,nq,thetas,rhos)
    implicit  none
#ifdef MPI
    include 'mpif.h'
#endif
    interface
       subroutine matvec(n_vec, source, destination)
         integer, intent(in) :: n_vec
         double precision, intent(in) :: source(:,:)
         double precision, intent(out) :: destination(:,:)
       end subroutine matvec
    end interface
    integer, intent(in) :: comm, imax
    real(8), intent(in) :: epsmax,epsmin
    real(8), intent(inout), target :: Q(:,:) ! matrix B, destructive
    complex(8), intent(in) :: sigma(:)
    integer, intent(out) :: iter
    real(8), intent(out) :: residual(:,:)
    complex(8), intent(out) :: X(:,:,:)
    integer, intent(in) :: nc,nq
    complex(8), intent(in) :: thetas(:)
    real(8), intent(in) :: rhos(:)
    

    integer :: i, j, k, n, L, m, ierr, maxidx, conv_count, shift_count, nprocs, myrank
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp
    real(8), allocatable :: res_seed(:), Bnorm_inv(:)
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    real(8), allocatable :: P_seed(:,:), AP(:,:), ts_tmp(:,:)
    complex(8), allocatable :: &
         P(:,:,:), PAP(:,:), PAP_old(:,:), alpha(:,:) &
         , delta(:,:), rho(:,:), rho_old(:,:), eye(:,:) &
         , blk_tmp1(:,:), blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:) &
         , ts_tmp_blk(:,:)
    real(8), allocatable :: r_PAP(:,:), r_rho(:,:), r_alpha(:,:), r_delta(:,:)

    complex(8) :: trace
    complex(8), allocatable :: ests(:)
    real(8), allocatable :: est_v(:,:)
    integer :: showfreq
    real(8) :: s_bb, avg, err, t

    myrank = 0
    nprocs = 1
#ifdef MPI
    call MPI_COMM_RANK(comm, myrank, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)
#endif

    call start_stopwatch(time_ld)
    showfreq = 50
    allocate(ests(nc))

    shift_count = 0

    n = size(Q,1)
    L = size(Q,2)
    m = size(X,3)
    allocate( est_v(L, nc) )

    allocate( P_seed(n,L), AP(n,L), ts_tmp(n,L) )
    allocate( P(L,L,m), ts_tmp_blk(L,L) &
         ,PAP(L,L),PAP_old(L,L),alpha(L,L) &
         ,delta(L,L),rho(L,L),rho_old(L,L),eye(L,L) &
         ,blk_tmp1(L,L),blk_tmp2(L,L) &
         ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
         ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L), iter_array(m) )
    allocate( r_PAP(L,L), r_rho(L,L), r_alpha(L,L), r_delta(L,L) )

    X = (0d0,0d0)
    eye = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do
!    Q = B

    call d_norm2_as_block(comm,Q,Bnorm_inv,n,L)
    Bnorm_inv = 1d0 / Bnorm_inv

    call d_MGS_QR(comm, Q, r_delta)
    delta = r_delta
    s_bb = sqrt(sum(r_delta**2))
    
    rho = delta
    P_seed = Q

    
    do k=1,m       
       P(:,:,k) = rho
       xi1(:,:,k) = rho
       xi2(:,:,k) = eye
    end do
    xitld_hat = (0d0,0d0)
    PAP = (0d0,0d0)
    residual(:,:) = 1d0
    iter_array(:) = 0
    
    do iter=1,imax

       call stop_stopwatch(time_ld)
       call matvec(L, P_seed, AP)
       call start_stopwatch(time_ld)


       AP = seed_val*P_seed - AP
       PAP_old = PAP
       call d_matmat_CN(comm, P_seed, AP, r_PAP)
       PAP = r_PAP
       call inv(PAP,alpha,ipiv,work)

       rho_old = rho
       
       r_alpha = alpha
       call d_matmat_NX('N', 1d0, AP, -r_alpha, Q)
       

       call d_MGS_QR(comm, Q, r_rho)
       rho = r_rho 

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta = blk_tmp1
       if (myrank==0) write(*,'(a,i6,2f18.10)') &
            ' COCG iter, res ', iter, sqrt(sum(dble(delta)**2))/s_bb, &
            sqrt(sum(dble(r_rho)**2))/s_bb

       call d_matmat_NX('T', 0d0, P_seed, r_rho, ts_tmp)
       P_seed = Q + ts_tmp

       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed = res_seed*Bnorm_inv

       ! call start_stopwatch(time_tmp)
       do k=1,m
          if ( mod((k-1),nprocs) == myrank ) then
             if ( maxval(residual(:,k)) > epsmin ) then             
                xi2(:,:,k) = xi1(:,:,k)
                blk_tmp1 = eye - xitld_hat(:,:,k)
                call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
                call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
                call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
                call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)
                call inv(eye - (seed_val - sigma(k))*alpha + blk_tmp1,xitld_hat(:,:,k),ipiv,work)
                call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
                call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
                call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)
                call z_matmat_NX('N',one,P(:,:,k),blk_tmp1,X(:,:,k))
                call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
                call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
                call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)
                call z_matmat_NX('N',zero,P(:,:,k),blk_tmp1,ts_tmp_blk)
                P(:,:,k) = ts_tmp_blk
                call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
                residual(:,k) = residual(:,k)*Bnorm_inv
                shift_count = shift_count + 1
             end if
             if ( maxval(residual(:,k)) > epsmax ) then
                iter_array(k) = iter
             end if
          else
             residual(:,k) = 0d0
          end if
       end do
       ! call stop_stopwatch(time_tmp)
#ifdef MPI
       call mpi_allreduce(MPI_IN_PLACE,residual,L*m,mpi_real8,mpi_sum,comm,ierr)
#endif
       ! if ( myrank == 0 .and. mod(iter,showfreq) == 0) then
       !    ! write residual norms
       !    write(365,*) 'iter = ', iter
       !    write(365,'(a5,1x,a10,1x,a10,1x,a15)') &
       !         'id','real','imag','rec res'
       !    do i=1,m
       !       maxidx = maxloc(residual(:,i),1)
       !       write(365,'(I5,1x,1pe10.2,1x,1pe10.2,1x,1pe15.7)') &
       !            i,real(sigma(i)),aimag(sigma(i)),maxval(residual(:,i))
       !    end do
       !    write(365,*)
       ! end if
          
       if ( mod(iter,showfreq) == 0) then
          !$OMP parallel do private(i, k, trace, j)
          do i = 1, nc
             do k = 1, L
                trace = (0d0, 0d0)
                do j = 1, nq
                   trace = trace + X(k,k,(i-1)*nq+j) * thetas((i-1)*nq+j)
                end do
                est_v(k,i) = trace * rhos(i) / nq
             end do
          end do
          !$OMP end parallel do

#ifdef MPI
          call mpi_allreduce(MPI_IN_PLACE, est_v, size(est_v), mpi_real8, &
               mpi_sum, comm, ierr)
#endif

          if (myrank==0) then 
             write(*,'(a,2i8)') "Estimation, iter = ", iter, nc
             write(*,*) '    i         eigenvalue count     error'
             t = 1.d0 / sqrt(dble(max(L-1, 1)))
             do i = 1, nc
                avg = sum(est_v(:,i)) / dble(L)
                err = sqrt( sum(est_v(:,i)**2)/dble(L) - avg**2 ) * t
                write(*,'(i6, 2f18.3)') i, avg, err                     
             end do
          end if
       end if

!        conv_count = 0
!        do k=1,m
!           if ( iter_array(k) < iter ) then
!              conv_count = conv_count + 1
!           end if
!        end do
       if ( maxval(residual) < epsmax ) then !.or. conv_count > m/2 ) then
          exit
       end if
    end do

    call stop_stopwatch(time_ld)
    
#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,iter_array,m,mpi_integer,mpi_sum,comm,ierr)
    call mpi_allreduce(MPI_IN_PLACE,X,L*L*m,mpi_complex16,mpi_sum,comm,ierr)
#endif

    ! if ( myrank == 0 ) then
    !    write(*,*) 'Result of shifted block CG rQ bilinear method'
    !    write(*,'(a5,1x,a10,1x,a10,1x,a6,1x,a15,1x,a4)') &
    !         'id','real','imag','iter','rec res','rhs'
    ! end if
    ! do i=1,m
    !    maxidx = maxloc(residual(:,i),1)       
    !    if ( myrank == 0 ) then
    !       write(*,'(I5,1x,1pe10.2,1x,1pe10.2,1x,I6,1x,1pe15.7,1x,I4)') &
    !            i,real(sigma(i)),aimag(sigma(i)),iter_array(i),residual(maxidx,i),maxidx
    !    end if
    ! end do
  end subroutine shifted_block_CG_rQ_bilinear

  subroutine z_MGS_QR(comm,Q,R)
    implicit none
    ! QR factorization via modified Gram-Schmidt method
    integer, intent(in) :: comm
    complex(8),intent(out) :: Q(:,:),R(:,:)
    integer i,j,k,s1,s2
    s1 = size(Q,1)
    s2 = size(Q,2)
    
    R = (0d0,0d0)
    do i=1,s2
       R(i,i) = dsqrt(dble(z_dot(comm,Q(:,i),Q(:,i))))
       !$OMP parallel do
       do k=1,s1
          Q(k,i) = Q(k,i) / R(i,i)
       end do
       !$OMP end parallel do
       do j=i+1,s2
          R(i,j) = z_dot(comm,Q(:,i),Q(:,j))
          !$OMP parallel do
          do k=1,s1
             Q(k,j) = Q(k,j) - R(i,j)*Q(k,i)
          end do
          !$OMP end parallel do
       end do
    end do
  end subroutine z_MGS_QR

  subroutine d_MGS_QR(comm,Q,R)
    implicit none
    ! QR factorization via modified Gram-Schmidt method
    integer, intent(in) :: comm
    real(8), intent(out) :: Q(:,:),R(:,:)
    integer i,j,k,s1,s2
    s1 = size(Q,1)
    s2 = size(Q,2)
    
    R = 0.d0
    do i = 1, s2
       R(i,i) = sqrt( d_dot(comm, Q(:,i), Q(:,i)) )
       !$OMP parallel do
       do k=1,s1
          Q(k,i) = Q(k,i) / R(i,i)
       end do
       !$OMP end parallel do
       do j=i+1,s2
          R(i,j) = d_dot(comm, Q(:,i), Q(:,j))
          !$OMP parallel do
          do k=1,s1
             Q(k,j) = Q(k,j) - R(i,j)*Q(k,i)
          end do
          !$OMP end parallel do
       end do
    end do
  end subroutine d_MGS_QR

  subroutine norm2_as_block(comm,V,norm,n,L)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    complex(8), intent(in) :: V(:,:)
    integer, intent(in) :: comm,n,L
    real(8), intent(out) :: norm(:)
    
    integer :: i,ierr
    real(8) :: temp(L)

    norm = 0D0
    do i=1,n
       temp = abs(V(i,:))
       norm = norm + temp*temp
    end do
#ifdef MPI
    temp = norm
    call mpi_allreduce(temp,norm,L,mpi_real8,mpi_sum,comm,ierr)
#endif
    norm = dsqrt(dble(norm))
  end subroutine norm2_as_block

  subroutine d_norm2_as_block(comm,V,norm,n,L)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    real(8), intent(in) :: V(:,:)
    integer, intent(in) :: comm,n,L
    real(8), intent(out) :: norm(:)
    
    integer :: i,ierr
    real(8) :: temp(L)

    norm = 0D0
    do i=1,n
       temp = abs(V(i,:))
       norm = norm + temp*temp
    end do
#ifdef MPI
    temp = norm
    call mpi_allreduce(temp,norm,L,mpi_real8,mpi_sum,comm,ierr)
#endif
    norm = dsqrt(dble(norm))
  end subroutine d_norm2_as_block


  subroutine norm2_as_block_serial(V,norm,n,L)
    implicit none
    complex(8), intent(in) :: V(:,:)
    integer, intent(in) :: n,L
    real(8), intent(out) :: norm(:)
    
    integer :: i,ierr
    real(8) :: temp(L)

    norm = 0D0
    do i=1,n
       temp = abs(V(i,:))
       norm = norm + temp*temp
    end do
    norm = dsqrt(dble(norm))
  end subroutine norm2_as_block_serial

  subroutine z_matmat_NX(trans,beta,X,Y,mat)
    implicit none
    character, intent(in) :: trans
    complex(8), intent(in) :: beta,X(:,:),Y(:,:)
    complex(8), intent(out) :: mat(:,:)

    integer :: ierr,m,n,k
    complex(8), parameter :: one = (1d0,0d0)
    m = size(X,1)
    n = size(mat,2)
    k = size(X,2)
    call ZGEMM('N',trans,m,n,k,one,X,m,Y,k,beta,mat,m)
  end subroutine z_matmat_NX

  subroutine d_matmat_NX(trans,beta,X,Y,mat)
    implicit none
    character, intent(in) :: trans
    real(8), intent(in) :: beta, X(:,:), Y(:,:)
    real(8), intent(out) :: mat(:,:)
    integer :: ierr,m,n,k
    character :: t

    m = size(X,1)
    n = size(mat,2)
    k = size(X,2)
    t = trans
    if (trans=='C' .or. trans=='c') t = 'T'
    call DGEMM('N', t, m, n, k, 1d0, X, m, Y, k, beta, mat, m)
  end subroutine d_matmat_NX

  subroutine inv(X,mat,ipiv,work)
    implicit none
    complex(8), intent(in) :: X(:,:)
    complex(8), intent(out) :: mat(:,:)
    integer, intent(out) :: ipiv(:)
    complex(8), intent(out) :: work(:)
    integer :: bs,info
    
    bs = size(X,1)
    mat = X
    call ZGETRF(bs,bs,mat,bs,ipiv,info)
    call ZGETRI(bs,mat,bs,ipiv,work,bs,info)
  end subroutine inv

  subroutine z_matmat_CN(comm,X,Y,mat)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: comm
    complex(8), intent(in) :: X(:,:),Y(:,:)
    complex(8), intent(out) :: mat(:,:)
    
    integer :: ierr,ms,bsX,bsY
    complex(8) :: zero,one
    complex(8),allocatable :: tmp(:,:)
    zero = (0d0,0d0); one = (1d0,0d0)
    ms = size(X,1)
    bsX = size(X,2)
    bsY = size(Y,2)
    allocate(tmp(bsX,bsY))
    call ZGEMM('C','N',bsX,bsY,ms,one,X,ms,Y,ms,zero,tmp,bsX)
#ifdef MPI
    call mpi_allreduce(tmp,mat,bsX*bsY,mpi_complex16,mpi_sum,comm,ierr)
#else
    mat = tmp
#endif
    deallocate(tmp)
  end subroutine z_matmat_CN


  subroutine d_matmat_CN(comm,X,Y,mat)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: comm
    real(8), intent(in) :: X(:,:),Y(:,:)
    real(8), intent(out) :: mat(:,:)
    
    integer :: ierr,ms,bsX,bsY
    real(8),allocatable :: tmp(:,:)

    ms = size(X,1)
    bsX = size(X,2)
    bsY = size(Y,2)
    call DGEMM('T','N',bsX,bsY,ms,1d0,X,ms,Y,ms,0d0,mat,bsX)
#ifdef MPI
    allocate( tmp(bsX, bsY) )
    tmp = mat
    call mpi_allreduce(tmp,mat,bsX*bsY,mpi_real8,mpi_sum,comm,ierr)
    deallocate(tmp)
#endif
  end subroutine d_matmat_CN


  complex(8) function z_dot(comm,x,y)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: comm
    complex(8), intent(in) :: x(:),y(:)
    
    integer :: i,ierr
    complex(8) :: tmp
    
    complex(8) :: zdotc

    tmp = (0d0,0d0)
    !$OMP parallel do reduction(+:tmp)
    do i=1,size(x)
       tmp = tmp + dconjg(x(i))*y(i)
    end do
    !$OMP end parallel do
#ifdef MPI
    call mpi_allreduce(tmp,z_dot,1,mpi_complex16,mpi_sum,comm,ierr)
#else
    z_dot = tmp
#endif
  end function z_dot


  function d_dot(comm,x,y) result (r)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: comm
    real(8), intent(in) :: x(:), y(:)
    real(8) :: r
    
    integer :: i,ierr
    real(8) :: tmp

    tmp = 0.d0
    !$OMP parallel do reduction(+:tmp)
    do i = 1, size(x)
       tmp = tmp + x(i)*y(i)
    end do
    !$OMP end parallel do
#ifdef MPI
    call mpi_allreduce(tmp, r, 1,mpi_real8, mpi_sum, comm, ierr)
#else
    r = tmp
#endif
  end function d_dot

  subroutine create_hutch_samples(V, nrow, ncol, rank)
    implicit none
    integer, intent(in) :: nrow, ncol, rank
    double precision, intent(out) :: V(nrow,*)
    integer :: i, j
    call create_rand_matrix(V, nrow, ncol, rank)
    !$OMP parallel do private(i, j) 
    do j = 1, ncol
       do i = 1, nrow
          V(i,j) = sign(1d0, real(V(i,j), kind(0d0)))
       end do
    end do    
    !$OMP end parallel do
  end subroutine

  subroutine create_rand_matrix(V, nrow, ncol, rank)
    implicit none
    integer, intent(in) :: nrow, ncol, rank
    double precision, intent(out) :: V(nrow,*)
    integer :: iseed(4)
    iseed(1) = modulo(rank-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(rank-4096, 4096) ! must be between 0 and 4095
    iseed(3) = modulo(rank, 4096) ! must be between 0 and 4095
    iseed(4) = 1 ! must be between 0 and 4095 and odd
    call DLARNV(2, iseed, nrow*ncol, V)    
  end subroutine

end module eigdensity_estimator_mod
