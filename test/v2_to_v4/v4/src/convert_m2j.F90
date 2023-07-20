module mod_rank_mat
  implicit none

  type rank_mat 
     integer :: rank
     integer, allocatable :: n(:) ! dimension of each rank
     real(8), allocatable :: p(:) ! value in 1-dim vector
  end type rank_mat

contains

  subroutine init_rank_mat(self, ndim)
    type(rank_mat), intent(inout) :: self
    integer, intent(in) :: ndim(:)

    self%rank = size(ndim)
    allocate( self%n( self%rank ) )
    self%n(:) = ndim(:)
    allocate( self%p( product( ndim ) ) )
    self%p = 0.d0

  end subroutine init_rank_mat


  subroutine set_rank_mat(self, loc, v) 
    type(rank_mat), intent(inout) :: self
    integer, intent(in) :: loc(:)
    real(8), intent(in) :: v
    integer :: i, n, m, ii( self%rank )
    
    n = 0
    m = 1
    do i = 1, self%rank
       n = n + m * loc(i)
       m = m * self%n(i)
    end do
    self%p(n) = v
  end subroutine set_rank_mat


  subroutine get_rank_mat(self, loc, v) 
    type(rank_mat), intent(inout) :: self
    integer, intent(in) :: loc(:)
    real(8), intent(out) :: v
    integer :: i, n, m, ii( self%rank )
    
    n = 0
    m = 1
    do i = 1, self%rank
       n = n + m * loc(i)
       m = m * self%n(i)
    end do
    v = self%p(n)
  end subroutine get_rank_mat


  subroutine prod_rankmat_mat(self, nni, other, out)
    ! out = sum_ni  self_i,j,..., ni, ... * other_{ni,nj}
    ! sum of 'nni'-th index
    ! other ... dense rank-2 matrix
    type(rank_mat), intent(in) :: self
    integer, intent(in) :: nni
    real(8), intent(in) :: other(:,:)
    type(rank_mat), intent(inout) :: out
    integer :: n1, n2, nk, nl, ndim(size(self%n)), i

    nk = size(other, 1)
    nl = size(other, 2)
    ndim(:) = self%n(:)
    ndim(nni) = nl
    if (nk /= self%n(nni)) stop 'inconsistent prod_rankmat_mat'

    if ( allocated( out%p )) then 
       if (any(self%n /= ndim)) call close_rank_mat(out)
    end if
    if (.not. allocated( out%p )) call init_rank_mat(out, ndim)

    n1 = 1
    if (nni > 1        ) n1 = product( self%n(:nni-1) )
    n2 = 1
    if (nni < self%rank) n2 = product( self%n(nni+1:) )


    if (.false.) then
       ! if ( n1 == 1 ) then 
       ! TODO : performance tuning for small ni case
       call dgemm( 'T', 'N', nl, n2, nk, 1.d0, &
            other, nk, self%p, nk, 0.d0, &
            out%p, nk)
    else
       do i = 1, n2
          call dgemm( 'N', 'N', n1, nl, nk, 1.d0, &
               self%p( n1*nk*i+1 ), n1, other, nk, 0.d0, &
               out%p(  n1*nk*i+1 ), n1)
       end do
    end if

  end subroutine prod_rankmat_mat


  subroutine close_rank_mat(self)
    type(rank_mat), intent(inout) :: self
    deallocate( self%n, self%p )
  end subroutine close_rank_mat
  

end module mod_rank_mat






module jscheme
  use partition, only: type_ptn_pn
  use constant, only: kwf, kdim, kmbit, maxchar
  use model_space
  use wavefunction, only: type_vec_p
  use mod_rank_mat
  implicit none


  type sb_m
     integer :: n ! dimension
     integer, allocatable :: mbit(:), jj_v_a(:,:)
     real(8), allocatable :: p(:,:)
  end type sb_m

  type sb_n
     type(sb_m), allocatable :: m(:)
  end type sb_n

  type type_seniority_basis
     type (sb_n), allocatable :: n(:)
  end type type_seniority_basis
  
  type (type_seniority_basis), allocatable :: singlejv(:)


  type mbit_singlej_m
     integer :: n
     integer(kmbit), allocatable :: mbit(:)
     integer, allocatable :: jva(:,:) ! J, v, alpha
     real(8), allocatable :: c_jva_mm(:,:)
  end type mbit_singlej_m

  type type_mbit_singlej_n
     type(mbit_singlej_m), allocatable :: mm(:)
     integer :: maxm
  end type type_mbit_singlej_n

  type(type_mbit_singlej_n), allocatable :: mbit_sj(:,:)
  integer, allocatable :: j_list(:)
  
  character(maxchar) :: fn_jbasis = 'seniority_basis.dat'



  
contains

  subroutine init_jscheme(ptn)
    type(type_ptn_pn), intent(in) :: ptn
    integer :: i, j, k, n, jt(n_jorb_pn), jv(n_jorb_pn)
    integer :: mj, ni, nj, sj, mm, nf, nnf, loop, mb, m, jj, vs, alph, na
    integer, allocatable :: mm_orb(:)

    mj = maxval(jorb) + 1
    jt(:) = jorb(:)
    n = 0
    do i = 1, n_jorb_pn
       j = minloc( jt, 1 )
       if (n /= 0) then
          if (jt(j) == jv(n)) then 
             jt(j) = mj
             cycle
          end if
       end if
       n = n + 1
       jv(n) = jt(j)
       jt(j) = mj
    end do
    allocate( j_list(n) )
    j_list(:) = jv(:n)  ! jorb without duplicacy

    if (myrank==0) write(*,*) "j_list", j_list

    allocate( mbit_sj( maxval(jorb), 0:min( maxval(ptn%n_ferm), maxval(j_list)+1 ) ) ) 

    ! prepare M-scheme bits in single j
    allocate( mm_orb( maxval(j_list)+1 ) )

    do i = 1, size(j_list)
       sj = j_list(i)
       do j = 1, sj+1
          mm_orb(j) = -sj + 2*(j-1)
       end do
       
       nnf = min( maxval(ptn%n_ferm), sj + 1)
       do nf = 0, nnf
          mm = - sum( (/( mm_orb(k), k=1, nf )/) )
          allocate( mbit_sj( sj, nf )%mm( -mm:mm ) )
          mbit_sj( sj, nf )%maxm = mm
          mbit_sj(sj, nf)%mm(:)%n = 0
       end do

       do loop = 1, 2

          do nf = 0, nnf
             mbit_sj(sj, nf)%mm(:)%n = 0
          end do

          do mb = 0, 2**(sj+1)-1
             nf = popcnt(mb)
             if (nf > nnf) cycle

             m = 0
             do k = 1, sj+1
                if (btest(mb, k-1)) m = m + mm_orb(k)
             end do

             mbit_sj(sj, nf)%mm(m)%n = mbit_sj(sj, nf)%mm(m)%n  + 1
             if (loop==2) mbit_sj(sj, nf)%mm(m)%mbit( mbit_sj(sj, nf)%mm(m)%n ) = mb
          end do

          if (loop==2) cycle
          do nf = 0, nnf
             do m = -mbit_sj(sj, nf)%maxm, mbit_sj(sj, nf)%maxm, 2
                allocate( mbit_sj(sj, nf)%mm(m)%mbit( mbit_sj(sj, nf)%mm(m)%n ) )
             end do
          end do

       end do
    end do

    deallocate( mm_orb )

    ! prepare J-scheme basis in single j
    do i = 1, size(j_list)
       sj = j_list(i)
       nnf = min( maxval(ptn%n_ferm), sj + 1)
       do nf = 0, nnf
          do m = -mbit_sj(sj, nf)%maxm, mbit_sj(sj, nf)%maxm, 2
             n = mbit_sj(sj, nf)%mm(m)%n 
             allocate( mbit_sj(sj, nf)%mm(m)%jva(3, n) )
             allocate( mbit_sj(sj, nf)%mm(m)%c_jva_mm(n, n) )
             mbit_sj(sj, nf)%mm(m)%c_jva_mm(:,:) = 0.d0
             k = 0
             do jj = abs(m),  mbit_sj(sj, nf)%maxm, 2
                do vs = mod(nf, 2), min(nf, sj/2+2), 2
                   do alph = 1, calc_dim_vj_s_j(sj, jj, vs)
                      k = k + 1
                      mbit_sj(sj, nf)%mm(m)%jva(1, k) = jj
                      mbit_sj(sj, nf)%mm(m)%jva(2, k) = vs
                      mbit_sj(sj, nf)%mm(m)%jva(3, k) = alph
                   end do
                end do
             end do
             if (k /= mbit_sj(sj, nf)%mm(m)%n) stop "single-j jscheme dim error"

          end do
       end do
    end do

    call read_seniority_basis( fn_jbasis )


    

  end subroutine init_jscheme



  subroutine finalize_jscheme(ptn)
    type(type_ptn_pn), intent(in) :: ptn

    deallocate( j_list )

  end subroutine finalize_jscheme



  subroutine read_seniority_basis(fname)
    character(len=*),intent(in) :: fname
    integer, parameter :: lun=20, nsize=300
    integer :: sj, ssj
    integer :: nf, v, alpha, jj, mm, ns, i, j, ierr, nj, k, l, m, kjva, km
    integer :: mbit(nsize)
    real(8) :: vc(nsize)

    open(lun, file=fname, form='unformatted', status='old', access='stream')

    ! singlejv(j)%n(n)%(m)%c-matrix, dim, m-scheme bit list     

    do while (.true.) 
       read(lun, iostat=ierr) sj
       if (ierr /= 0) exit
       read(lun) nj
       do j = 1, nj
          read(lun) nf, v, alpha, jj, mm, ns
          if (ns > nsize) stop 'increase nsize in read_seniority_basis'
          do i = 1, ns
             read(lun) mbit(i), vc(i)
          end do

          if (.not.( any(sj == j_list) )) cycle
          if ( nf > ubound(mbit_sj, 2)) cycle

          do kjva = 1, mbit_sj(sj, nf)%mm(mm)%n
             if ( all( mbit_sj(sj, nf)%mm(mm)%jva(:, kjva) &
                  &    == (/ jj, v, alpha+1 /) ) ) goto 100
          end do
          stop 'NOT found JJ, v, alpha'
100       continue

          do i = 1, ns
             do km = 1, mbit_sj(sj, nf)%mm(mm)%n
                if ( mbit_sj(sj, nf)%mm(mm)%mbit(km) == mbit(i) ) goto 200
             end do
             mbit_sj(sj, nf)%mm(mm)%c_jva_mm(kjva, km) = vc(i)
          end do
          stop 'NOT found M-scheme bit'
200       continue

          
       end do
    end do

    close(lun)
  end subroutine read_seniority_basis




  recursive function young_t(h, w, n) result (r)
    ! Young tableau
    ! number of states of n fermions in j orbit with Jz=M is 
    ! young(2j+1-n, n, (2j+1-n)*n/2)-M)
    integer, intent(in) :: h, w, n
    integer :: r

    r = 0
    if ( h<0 .or. w<0 .or. n<0 ) return
    if ( n == 0 ) then 
       r = 1
       return
    end if

    r = young_t(h, w-1, n-h) + young_t(h-1, w, n)
  end function young_t


  function calc_dim_vj_s_j(sj, jj, v) result (r)
    integer, intent(in) :: sj, jj, v
    integer :: r
    ! J-scheme dimension of J,v  for single-j(sj) orbit
    r  =   young_t(sj+1-v, v,   ((sj+1-v)*v-jj)/2 ) &
         - young_t(sj+1-v, v,   ((sj+1-v)*v-jj)/2-1 ) &
         - young_t(sj+3-v, v-2, ((sj+3-v)*(v-2)-jj)/2 ) &
         + young_t(sj+3-v, v-2, ((sj+3-v)*(v-2)-jj)/2-1 )
  end function calc_dim_vj_s_j



  subroutine main_convert_m2j(ptn, evec)
    type(type_ptn_pn), intent(in) :: ptn
    type(type_vec_p), intent(in) :: evec(:)




  end subroutine main_convert_m2j




  subroutine convert_m2_sj_ptn(ptn, id, evec, eout)
    type(type_ptn_pn), intent(in) :: ptn
    integer, intent(in) :: id
    real(kwf), intent(in) :: evec(:)
    real(8), intent(out) :: eout(:)
    eout = 0.d0
    ! olp(:) = 0.d0
    !   !$omp parallel do private(id, mq, idpn, n) reduction (+: olp)
    !   do id = ptn%idl_start, ptn%idl_end
    !      idpn(:) = ptn%pidpnM_pid_srt(1:2, ptn%pid_dpl2srt(id))
    !      n = n_ptn(idpn(1), 1) + n_ptn(idpn(2), 2) 
    !      do mq = ptn%local_dim_acc_start(id),  ptn%local_dim_acc(id)
    !         olp(n) = olp(n) + self%p(mq) * self%p(mq)
    !      end do
    !   end do
    

  end subroutine convert_m2_sj_ptn




end module jscheme







program convert_m2j
  ! convert_m2j.exe foo.snt bar.ptn M-scheme.wav J-scheme.jwv
  !$ use omp_lib
#ifdef MPI
  use mpi
#endif
  use constant, only: kwf, kdim, kmbit, nmbit, maxchar, c_no_init, pi
  use jscheme
  use partition, only: init_partition, finalize_partition, deploy_partition
  use wavefunction, only: type_vec_p, load_wf, dot_product_global
  use bp_io, only: bp_load_wf
  implicit none
  character(len=maxchar) :: fn_snt, fn_ptn, fn_wav, fn_out
  integer, parameter :: lunint=11, lunptn=12, &
       lunwv=13, lunout=14
  type(type_ptn_pn) :: ptn
  integer :: n_eig, mtotal
  type(type_vec_p), allocatable :: evec(:), vec_sj(:)
  

#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank/=0) is_debug = .false.
#endif
  !$ if(myrank==0) write(*,'(1a,1i3)') "OpenMP  # of threads=", omp_get_max_threads()

  call getarg(1, fn_snt)
  call getarg(2, fn_ptn)
  call getarg(3, fn_wav)
  call getarg(4, fn_out)
#ifdef MPI
  call mpi_bcast(fn_snt,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_ptn,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_wav,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_out,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
#endif

#ifdef MPI
  call init_mpi_shift_reduce()
#endif

  if (nv_shift == 0) nv_shift = 1

  if (myrank==0) then 
     write(*,'(a,4i3)') &
          "compile conf. kwf, kdim, kmbit, nmbit =", kwf, kdim, kmbit, nmbit
     write(*,'(/,2a)') ' M-scheme wave function : ', trim(fn_wav)
     write(*,'(/,a)')  ' is converted into '
     write(*,'(/,2a)') ' J-scheme wave function : ', trim(fn_out)
  end if



  open(lunint, file=fn_snt, status='old')
  call read_sps(lunint)


  ! read header of wave functions
  open(lunwv, file=fn_wav, form='unformatted', &
       status='old', access='stream')
  open(lunwv, file=fn_wav, form='unformatted', &
       status='old', access='stream')
  read(lunwv) n_eig
  read(lunwv) mtotal
  close(lunwv)


  open(lunptn, file=fn_ptn, status='old')
  if (myrank==0) write(*,'("set partition_file=",1a)') trim(fn_ptn)
  call init_partition(ptn, lunptn, mtotal)
  close(lunptn)

  call deploy_partition(ptn)

  allocate( evec(n_eig) )
  
  call bp_load_wf(fn_wav, evec, ptn, fn_ptn, mtotal)

  
  call init_jscheme(ptn)







  call finalize_partition(ptn)

end program convert_m2j
