module operator_mscheme
  !
  ! m-scheme operator
  !
  !$ use omp_lib
  use model_space, only: korb, itorb, n_jorb, n_jorb_pn, jorb, iporb, n_morb, &
       & n_morb_pn, morbn, korbn, myrank, jorbn
  use operator_jscheme, only : opr_j, jcouple, jcouplemax
  use rotation_group, only :dcg
  implicit none

  private
  ! definition of structures
  public ::  opr_m, opr_m_p, idx_nocc2b, jz_2b_v, v_2b, jz_2b_idx, add_opr_m, &
       opr_m_eff_charge, print_operator_mscheme, idx_gt, idx_2b, print_op_1b
  ! main procedures
  public :: operator_j2m, opr_m_one_crt, opr_m_two_crt, opr_m_one_anh, &
       finalize_opr_m
  public :: type_type_id_connect_ptn
  public :: n_id_idx ! # of idx_nocc2b(ipn, n1, n2)%md,p(m)%id, p,n,pp,nn

  ! TBTD
  public :: init_tbtd_op, clear_op, get_cpld_tbtd, get_cpld_obtd
  public :: init_idx_tbtd, finalize_idx_tbtd
  public :: idx_pph, idx_p, tbtd_container, init_copy_opr_m, print_tbtd_op

  ! index for opr_m, initialized in init_operator_mscheme
  ! (i,j) = idx_nocc2b(ipn, n1, n2)%md,p(m)%idx(:,ij)
  type idx_2b
     integer, allocatable :: idx(:,:)  ! m1,m2=idx(:,ij)
     integer :: id                     ! for p-n int. jump store
  end type idx_2b
  integer :: n_id_idx(4)

  type jz_2b_idx
     type(idx_2b), allocatable :: mp(:)  ! index of pairing combination
     type(idx_2b), allocatable :: md(:)  ! index of density combination
  end type jz_2b_idx

  type(jz_2b_idx), target, allocatable :: idx_nocc2b(:,:,:)

  ! md : 1: c_p+ c_n,  2: c_n+ c_p
  type(jz_2b_idx), target, allocatable :: idx_gt(:,:,:)


  ! 2-particle 1-hole index for ppp, nnn  Note : md is doulbled, half integer
  type(jz_2b_idx), target, allocatable :: idx_pph(:,:,:,:)
  type(jz_2b_idx), target, allocatable :: idx_p(:,:)



!!! definition for opr_m, operator in m-scheme
  type v_1b
     real(8), allocatable :: v(:)     ! one body term
  end type v_1b

  type jz_1b_v
     type(v_1b), allocatable :: m(:)  ! Jz
  end type jz_1b_v

  type v_2b
     real(8), allocatable :: v(:,:)   ! v((m1,m2),(m3,m4))
  end type v_2b

  type jz_2b_v
     type(v_2b), allocatable :: m(:)  ! Jz
  end type jz_2b_v

  type type_id_connect_ptn
     integer :: n
     integer, allocatable :: id(:)
  end type type_id_connect_ptn

  type type_type_id_connect_ptn
     type(type_id_connect_ptn), allocatable :: idl(:)
  end type type_type_id_connect_ptn

  type opr_m
     integer :: irank
     integer :: nbody
     integer :: ipr1_type
     ! one-body diagonal, rank 0 if nbody==2
     type(v_1b) :: spe(2)
     ! op(pn,j1,j2,j3,j4)%m(m)%v(j12,j34) j1<=j2, j3<=j4  TBME
     type(jz_2b_v), allocatable :: nocc2b(:,:,:,:,:)
     logical :: is_j_square = .false.  ! flag if JJ
     !
     ! general rank one-body operator in density decomposition
     type(jz_1b_v), allocatable :: nocc1b(:,:,:)
     !
     ! creation operator / annihilation operator
     integer :: crt_orb, crt_idx
     real(8) :: crt_v
     !
     ! partition information
     type(type_type_id_connect_ptn), allocatable :: mlmr(:,:)
     !
     ! zero-body term
     real(8) :: e0 = 0.d0
     !
     ! 3-body monopole term
     integer :: n_three_body_mp = 0
     integer, allocatable :: idx_three_body_mp(:,:,:)
     real(8), allocatable :: v_three_body_mp(:)
     !
     ! Jz for two-body transitiion density, not doubled
     integer :: mm = 0
  end type opr_m

  type opr_m_p
     type(opr_m), pointer :: p
  end type opr_m_p



  type tbtd_container
     integer :: n12, n1234
     integer, allocatable :: k1234jj(:,:), k12jj(:,:)
     real(8), allocatable :: v1234jj(:),   v12jj(:)
  end type tbtd_container


contains


  subroutine operator_j2m(oj, om)
    !
    ! get an m-scheme operator "om" from j-scheme operator "oj"
    !
    type(opr_j), intent(in) :: oj
    type(opr_m), intent(out) :: om
    integer :: ipn, iprty
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, ip12, m1, m2, m3, m4, mm1, mm2
    integer :: n, n1, n2, n3, n4, n12, k13
    integer :: jcpl, jcplmin, jcplmax, ij12, ij34
    integer :: nj, i, j, k, l, ij, kl, ik, jl, mm, md, maxm, m12
    real(8) :: v, c12, c34
    type(jz_2b_idx), pointer :: noc1, noc2
    logical :: is
    real(8) :: t(n_jorb_pn, n_jorb_pn)

    om%irank = oj%irank
    om%mm    = 0
    om%nbody = oj%nbody
    om%ipr1_type = oj%ipr1_type
    om%is_j_square = .false.
    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    nj = maxval(n_jorb)

    ! --- two-body beta ---
    if (om%nbody == -12 .or. om%nbody == -13) then
       call j2m_tb_beta()
       return
    end if

    ! -------------- 1-body rank-k operator -----------------------
    if (om%nbody == 1) then
       allocate( om%nocc1b(2, nj, nj) )
       iprty = 1
       if (om%ipr1_type /= 1) iprty = -1

       do k1 = 1, n_jorb_pn
          n1 = k1
          if ( n1 > n_jorb(1) ) n1 = n1 - n_jorb(1)
          j1 = jorb(k1)
          ipn = 1
          if (itorb(k1) == 1) ipn = 2
          do k2 = 1, n_jorb_pn
             if (itorb(k1) /= itorb(k2)) cycle
             if (iporb(k1)*iporb(k2) /= iprty) cycle
!             if (abs(oj%p1%v(k1,k2)) < 1.d-8) cycle
             n2 = k2
             if ( n2 > n_jorb(1) ) n2 = n2 - n_jorb(1)
             j2 = jorb(k2)
             noc1 => idx_nocc2b(ipn, n1, n2)
             if ( .not. allocated(noc1%md) ) cycle
             m1 = max( -om%irank, lbound(noc1%md, 1) )
             m2 = min(  om%irank, ubound(noc1%md, 1) )
             allocate( om%nocc1b(ipn,n1,n2)%m(m1:m2) )
             do md = m1, m2
                n = size(noc1%md(md)%idx, 2)
                allocate(om%nocc1b(ipn,n1,n2)%m(md)%v(n))
                om%nocc1b(ipn, n1, n2)%m(md)%v = 0.d0
                do ij = 1, n
                   i = noc1%md(md)%idx(1, ij)
                   j = noc1%md(md)%idx(2, ij)
                   mm1 = morbn(i, ipn)
                   mm2 = morbn(j, ipn)
                   om%nocc1b(ipn, n1, n2)%m(md)%v(ij) &
                        = oj%p1%v(k1, k2) &
                        * dcg(j2, mm2, 2*om%irank, 2*md, j1, mm1) &
                        / sqrt(dble(j1+1))
                end do
                if ( maxval(abs(om%nocc1b(ipn, n1, n2)%m(md)%v)) < 1.d-8 ) &
                     deallocate( om%nocc1b(ipn, n1, n2)%m(md)%v )
             end do
          end do
       end do
       return
    end if

    ! --- beta decay 1-body rank-k operator ---
    if (om%nbody == -10 .or. om%nbody == -11) then
       call j2m_gt()
       return
    end if


    ! -------------- 2-body rank-0 operator -----------------------
    if (oj%irank/=0) stop "not implimented yet 2-body rank-0"
    allocate( om%nocc2b(3, nj, nj, nj, nj) )


    ! one-body rank=0 operator
    do ipn = 1, 2
       allocate( om%spe(ipn)%v(n_jorb(ipn)) )
       do k1 = 1, n_jorb(ipn)
          n1 = k1
          if (ipn==2) n1 = k1 + n_jorb(1)
          om%spe(ipn)%v(k1) = oj%p1%v(n1,n1)/sqrt(dble(jorb(n1)+1))
       end do
    end do

    ! not implemented non-diagonal OBME yet
    t(:,:) = oj%p1%v(:,:)
    forall(k1=1:n_jorb_pn) t(k1,k1) = 0.d0
    if (maxval(abs(t)) > 1.d-5) stop 'non-diag OBME error'

    if (om%nbody == 1) return

    !  p-p TBME in pairing combination
    ! do n1 = 1, n_jorb(1)
    !    do n2 = n1, n_jorb(1)
    ipn = 1
    !$omp parallel do private( n12, n1, n2, n3, n4, ip12, maxm, mm, ij, &
    !$omp &                    i, j, k, l, j1, j2, j3, j4, m1, m2, m3, m4, v, &
    !$omp &                    c12, c34, jcplmin, jcplmax, jcpl, ij12, ij34 ) &
    !$omp &           schedule( dynamic )
    do n12 = 1, n_jorb(1) * n_jorb(1)
       n1 = (n12 - 1) / n_jorb(1) + 1
       n2 = mod(n12 - 1, n_jorb(1)) + 1
       if (n1>n2) cycle

       if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle
       ! noc1 => idx_nocc2b(ipn,n1,n2)
       ip12 = iporb(n1)*iporb(n2)
       if (ip12==-1) ip12 = 2
       do n3 = 1, n_jorb(1)
          do n4 = n3, n_jorb(1)
             if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
             if (iporb(n1)*iporb(n2) /= iporb(n3)*iporb(n4)) cycle
             ! noc2 => idx_nocc2b(ipn,n3,n4)
             maxm = min( ubound(idx_nocc2b(ipn,n1,n2)%mp,1), &
                  &      ubound(idx_nocc2b(ipn,n3,n4)%mp,1) )
             allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
             do mm = -maxm, maxm
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                     size(idx_nocc2b(ipn,n1,n2)%mp(mm)%idx, 2), &
                     size(idx_nocc2b(ipn,n3,n4)%mp(mm)%idx, 2) ))
                do ij = 1,  size(idx_nocc2b(ipn,n1,n2)%mp(mm)%idx, 2)
                   i = idx_nocc2b(ipn,n1,n2)%mp(mm)%idx(1,ij)
                   j = idx_nocc2b(ipn,n1,n2)%mp(mm)%idx(2,ij)
                   j1 = jorb(n1)
                   j2 = jorb(n2)
                   m1 = morbn(i,ipn)
                   m2 = morbn(j,ipn)
                   do kl = 1,  size(idx_nocc2b(ipn,n3,n4)%mp(mm)%idx, 2)
                      k = idx_nocc2b(ipn,n3,n4)%mp(mm)%idx(1,kl)
                      l = idx_nocc2b(ipn,n3,n4)%mp(mm)%idx(2,kl)
                      j3 = jorb(n3)
                      j4 = jorb(n4)
                      m3 = morbn(k,ipn)
                      m4 = morbn(l,ipn)
                      v = 0.0d0
                      c12 = 1.0d0
                      if (n1 == n2) c12 = sqrt(2.0d0)
                      c34 = 1.0d0
                      if (n3 == n4) c34 = sqrt(2.0d0)
                      jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                      jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                      do jcpl = jcplmin, jcplmax
                         ij12 = jcouple(jcpl,ip12,ipn)%idxrev(n1,n2)
                         ij34 = jcouple(jcpl,ip12,ipn)%idxrev(n3,n4)
                         if (ij12*ij34 == 0) cycle
                         ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                         v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                              & * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                              & * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                      end do
                      om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                           & = c12 * c34 * v
                   end do
                end do
             end do
          end do
       end do
    end do

    !  n-n TBME in pairing combination
!    do k1 = n_jorb(1)+1, n_jorb_pn
!       n1 = k1 - n_jorb(1)
!       do k2 = k1, n_jorb_pn
!          n2 = k2 - n_jorb(1)
    ipn = 2
    !$omp parallel do private( n12, n1, n2, n3, n4, k1, k2, k3, k4, &
    !$omp &                    ip12, maxm, mm, ij, &
    !$omp &                    i, j, k, l, j1, j2, j3, j4, m1, m2, m3, m4, v, &
    !$omp &                    c12, c34, jcplmin, jcplmax, jcpl, ij12, ij34 ) &
    !$omp &           schedule( dynamic )
    do n12 = 1, n_jorb(2) * n_jorb(2)
       n1 = (n12 - 1) / n_jorb(2) + 1
       n2 = mod(n12 - 1, n_jorb(2)) + 1
       if (n1 > n2) cycle
       k1 = n1 + n_jorb(1)
       k2 = n2 + n_jorb(1)
       if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle
       ! noc1 => idx_nocc2b(ipn, n1, n2)
       ip12 = iporb(k1)*iporb(k2)
       if (ip12==-1) ip12 = 2
       do k3 = n_jorb(1)+1, n_jorb_pn
          n3 = k3 - n_jorb(1)
          do k4 = k3, n_jorb_pn
             n4 = k4 - n_jorb(1)
             if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
             if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
             ! noc2 => idx_nocc2b(ipn, n3, n4)
             maxm = min( ubound(idx_nocc2b(ipn, n1, n2)%mp,1), &
                  &      ubound(idx_nocc2b(ipn, n3, n4)%mp,1))
             allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
             do mm = -maxm, maxm
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                     size(idx_nocc2b(ipn, n1, n2)%mp(mm)%idx, 2), &
                     size(idx_nocc2b(ipn, n3, n4)%mp(mm)%idx, 2) ))
                do ij = 1,  size(idx_nocc2b(ipn, n1, n2)%mp(mm)%idx, 2)
                   i = idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(1,ij)
                   j = idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(2,ij)
                   j1 = jorb(k1)
                   j2 = jorb(k2)
                   m1 = morbn(i,ipn)
                   m2 = morbn(j,ipn)
                   do kl = 1,  size(idx_nocc2b(ipn, n3, n4)%mp(mm)%idx, 2)
                      k = idx_nocc2b(ipn, n3, n4)%mp(mm)%idx(1,kl)
                      l = idx_nocc2b(ipn, n3, n4)%mp(mm)%idx(2,kl)
                      j3 = jorb(k3)
                      j4 = jorb(k4)
                      m3 = morbn(k,ipn)
                      m4 = morbn(l,ipn)
                      v = 0.0d0
                      c12 = 1.0d0
                      if (n1 == n2) c12 = sqrt(2.0d0)
                      c34 = 1.0d0
                      if (n3 == n4) c34 = sqrt(2.0d0)
                      jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                      jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                      do jcpl = jcplmin, jcplmax
                         ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                         ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k3,k4)
                         if (ij12*ij34 == 0) cycle
                         ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                         v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                              & * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                              & * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                      end do
                      om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                           & = c12 * c34 * v
                   end do
                end do
             end do
          end do
       end do
    end do


!    end do

    !  p-n TBME in density combination
    ipn = 3
!    do k1 = 1, n_jorb(1)
!       n1 = k1
!       do k3 = 1, n_jorb(1)
!          n3 = k3
    !$omp parallel do private( k13, k1, k2, k3, k4, n1, n2, n3, n4, ip12, &
    !$omp &                    maxm, mm, ik, jl, i, j, k, l, j1, j2, j3, j4, &
    !$omp &                    m1, m2, m3, m4, m12, v, jcplmin, jcplmax, jcpl, ij12, ij34) &
    !$omp &           schedule (dynamic)
    do k13 = 1, n_jorb(1)*n_jorb(1)
       k1 = (k13 - 1) / n_jorb(1) + 1
       k3 = mod(k13 - 1, n_jorb(1)) + 1
       n1 = k1
       n3 = k3
       if (.not. allocated(idx_nocc2b(1,n1,n3)%md)) cycle
       ! noc1 => idx_nocc2b(1, n1, n3)
       do k2 = n_jorb(1)+1, n_jorb_pn
          n2 = k2 - n_jorb(1)
          ip12 = iporb(k1)*iporb(k2)
          if (ip12==-1) ip12 = 2
          do k4 = n_jorb(1)+1, n_jorb_pn
             if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
             n4 = k4 - n_jorb(1)
             if (.not. allocated(idx_nocc2b(2,n2,n4)%md)) cycle
             ! noc2 => idx_nocc2b(2, n2, n4)
             maxm = min(ubound(idx_nocc2b(1, n1, n3)%md,1), ubound(idx_nocc2b(2, n2, n4)%md,1))
             allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
             do mm = -maxm, maxm
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                     size(idx_nocc2b(1, n1, n3)%md( mm)%idx, 2), &
                     size(idx_nocc2b(2, n2, n4)%md(-mm)%idx, 2) ))
                do ik = 1, size(idx_nocc2b(1, n1, n3)%md(mm)%idx, 2)
                   i = idx_nocc2b(1, n1, n3)%md(mm)%idx(1,ik)
                   k = idx_nocc2b(1, n1, n3)%md(mm)%idx(2,ik)
                   j1 = jorb(k1)
                   j3 = jorb(k3)
                   m1 = morbn(i,1)
                   m3 = morbn(k,1)
                   do jl = 1,  size(idx_nocc2b(2, n2, n4)%md(-mm)%idx, 2)
                      j = idx_nocc2b(2, n2, n4)%md(-mm)%idx(1,jl)
                      l = idx_nocc2b(2, n2, n4)%md(-mm)%idx(2,jl)
                      j2 = jorb(k2)
                      j4 = jorb(k4)
                      m2 = morbn(j,2)
                      m4 = morbn(l,2)
                      m12 = (m1+m2)/2
                      v = 0.0d0
                      jcplmin = max(abs(m12), abs(j1-j2)/2, abs(j3-j4)/2)
                      jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                      do jcpl = jcplmin, jcplmax
                         ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                         ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k3,k4)
                         if (ij12*ij34 == 0) cycle
                         ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                         v = v +  oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                              & * dcg(j1,m1,j2,m2,2*jcpl,2*m12) &
                              & * dcg(j3,m3,j4,m4,2*jcpl,2*m12)
                      end do
                      om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ik, jl) = v
                   end do
                end do
             end do
          end do
       end do
    end do

    !return
    ! deallocate zero partition
    do ipn = 1, 3
       do n1 = 1, maxval(n_jorb)
          do n2 = 1, maxval(n_jorb)
             do n3 = 1, maxval(n_jorb)
                do n4 = 1, maxval(n_jorb)
                   if (.not. allocated(om%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                   is = .true.
                   do mm = lbound(om%nocc2b(ipn,n1,n2,n3,n4)%m, 1), &
                        ubound(om%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
                      if ( maxval(abs( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )) > 1.d-8) then
                         is = .false.
                         cycle
                      end if
                      deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )
                   end do
                   if (is) deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m )
                end do
             end do
          end do
       end do
    end do

  contains

    subroutine j2m_gt()
      call init_op_m_gt()
      allocate( om%nocc1b(2, nj, nj) )
      iprty = 1
      if (om%ipr1_type /= 1) iprty = -1

      do k1 = 1, n_jorb_pn
         n1 = k1
         if ( n1 > n_jorb(1) ) n1 = n1 - n_jorb(1)
         j1 = jorb(k1)
         ipn = 1
         if (itorb(k1) == 1) ipn = 2
         do k2 = 1, n_jorb_pn
            if (itorb(k1)==itorb(k2)) cycle
            n2 = k2
            if ( n2 > n_jorb(1) ) n2 = n2 - n_jorb(1)
            j2 = jorb(k2)
            if (iporb(k1)*iporb(k2) /= iprty) cycle
            noc1 => idx_gt(ipn, n1, n2)
            !             if (abs(oj%p1%v(k1,k2)) < 1.d-8) cycle

            if ( .not. allocated(noc1%md) ) cycle
            m1 = max( -om%irank, lbound(noc1%md, 1) )
            m2 = min(  om%irank, ubound(noc1%md, 1) )
            allocate( om%nocc1b(ipn,n1,n2)%m(m1:m2) )
            do md = m1, m2
               n = size( noc1%md(md)%idx, 2 )
               allocate( om%nocc1b(ipn, n1, n2)%m(md)%v(n) )
               om%nocc1b(ipn, n1, n2)%m(md)%v = 0.d0
               do ij = 1, n
                  i = noc1%md(md)%idx(1, ij)
                  j = noc1%md(md)%idx(2, ij)
                  mm1 = morbn(i, ipn)
                  mm2 = morbn(j, 3-ipn)
                  om%nocc1b(ipn, n1, n2)%m(md)%v(ij) &
                       = oj%p1%v(k1, k2) &
                       * dcg(j2, mm2, 2*om%irank, 2*md, j1, mm1) &
                       / sqrt(dble(j1+1))
               end do
               if ( maxval(abs(om%nocc1b(ipn, n1, n2)%m(md)%v)) < 1.d-8 ) &
                    deallocate( om%nocc1b(ipn, n1, n2)%m(md)%v )
            end do
            is = .true.
            do md = m1, m2
               if ( allocated( om%nocc1b(ipn, n1, n2)%m(md)%v ) ) is = .false.
            end do
            if (is) deallocate( om%nocc1b(ipn, n1, n2)%m )
         end do
      end do
    end subroutine j2m_gt


    subroutine j2m_tb_beta()
      integer :: ipn, inp, nj1, nj2
      if (om%nbody==-12) then
         ipn = 1
         inp = 2
         nj1 = 0
         nj2 = n_jorb(1)
      else
         ipn = 2
         inp = 1
         nj1 = n_jorb(1)
         nj2 = 0
      end if
      allocate( om%nocc2b(1, nj, nj, nj, nj) )

      !  pp-nn TBME in pairing combination
      do n1 = 1, n_jorb(ipn)
         k1 = n1 + nj1
         do n2 = n1, n_jorb(ipn)
            k2 = n2 + nj1
            if (.not. allocated(idx_nocc2b(ipn, n1, n2)%mp)) cycle
            noc1 => idx_nocc2b(ipn, n1, n2)
            ip12 = iporb(k1)*iporb(k2)
            if (ip12==-1) ip12 = 2
            do n3 = 1, n_jorb(inp)
               k3 = n3 + nj2
               do n4 = n3, n_jorb(inp)
                  k4 = n4 + nj2
                  if (.not. allocated(idx_nocc2b(inp, n3, n4)%mp)) cycle
                  if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
                  noc2 => idx_nocc2b(inp, n3, n4)
                  maxm = min(ubound(noc1%mp,1), ubound(noc2%mp,1))
                  allocate(om%nocc2b(1,n1,n2,n3,n4)%m(-maxm:maxm))
                  do mm = -maxm, maxm
                     allocate(om%nocc2b(1,n1,n2,n3,n4)%m(mm)%v( &
                          size(noc1%mp(mm)%idx, 2), &
                          size(noc2%mp(mm)%idx, 2) ))
                     do ij = 1, size(noc1%mp(mm)%idx, 2)
                        i = noc1%mp(mm)%idx(1,ij)
                        j = noc1%mp(mm)%idx(2,ij)
                        j1 = jorb(k1)
                        j2 = jorb(k2)
                        m1 = morbn(i,ipn)
                        m2 = morbn(j,ipn)
                        do kl = 1, size(noc2%mp(mm)%idx, 2)
                           k = noc2%mp(mm)%idx(1,kl)
                           l = noc2%mp(mm)%idx(2,kl)
                           j3 = jorb(k3)
                           j4 = jorb(k4)
                           m3 = morbn(k,inp)
                           m4 = morbn(l,inp)
                           v = 0.0d0
                           c12 = 1.0d0
                           if (n1 == n2) c12 = sqrt(2.0d0)
                           c34 = 1.0d0
                           if (n3 == n4) c34 = sqrt(2.0d0)
                           jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                           jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                           do jcpl = jcplmin, jcplmax
                              ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                              ij34 = jcouple(jcpl,ip12,inp)%idxrev(k3,k4)
                              if (ij12*ij34 == 0) cycle
                              ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                              v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                                   * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                                   * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                           end do
                           om%nocc2b(1,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                                = c12 * c34 * v
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do


    end subroutine j2m_tb_beta


  end subroutine operator_j2m



  subroutine finalize_opr_m(op)
    type(opr_m), intent(inout) :: op
    integer :: ipn, n1, n2, n3, n4, mm, nj

    nj = maxval(n_jorb)

    if (allocated(op%nocc1b)) then
       do ipn = 1, 2
          do n1 = 1, size(op%nocc1b, 2)
             do n2 = 1, size(op%nocc1b, 3)
                if ( .not. allocated( op%nocc1b(ipn, n1, n2)%m) ) cycle
                do mm = lbound( op%nocc1b(ipn,n1,n2)%m, 1 ), &
                     ubound(op%nocc1b(ipn,n1,n2)%m, 1)
                   if ( .not. allocated( op%nocc1b(ipn, n1, n2)%m(mm)%v ) ) cycle
                   deallocate( op%nocc1b(ipn, n1, n2)%m(mm)%v )
                end do
                deallocate( op%nocc1b(ipn, n1, n2)%m )
             end do
          end do
       end do
       deallocate(op%nocc1b)
    end if


    if (allocated( op%nocc2b )) then
       do ipn = lbound(op%nocc2b, 1), ubound(op%nocc2b, 1)
          do n1 = 1, nj
             do n2 = 1, nj
                do n3 = 1, nj
                   do n4 = 1, nj
                      if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                      do mm = lbound(op%nocc2b(ipn,n1,n2,n3,n4)%m, 1), &
                           ubound(op%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
                         if ( allocated( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v ) ) &
                              deallocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )
                      end do
                      deallocate( op%nocc2b(ipn,n1,n2,n3,n4)%m )
                   end do
                end do
             end do
          end do
       end do
       deallocate( op%nocc2b )
    end if

    op%n_three_body_mp = 0
    if (allocated( op%idx_three_body_mp )) &
         deallocate( op%idx_three_body_mp )
    if (allocated( op%v_three_body_mp )) &
         deallocate( op%v_three_body_mp )

  end subroutine finalize_opr_m



  subroutine init_copy_opr_m(fr, op)
    type(opr_m), intent(in) :: fr
    type(opr_m), intent(out) :: op
    integer :: ipn, n1, n2, n3, n4, mm, nj, npn, nml, nmu


    op%irank       = fr%irank
    op%mm          = fr%mm
    op%nbody       = fr%nbody
    op%ipr1_type   = fr%ipr1_type
    op%is_j_square = fr%is_j_square

    nj = maxval(n_jorb)

    if (allocated(fr%nocc1b)) then
       npn = size(fr%nocc1b, 1)
       allocate( op%nocc1b(npn, nj, nj) )
       do n2 = 1, nj
          do n1 = 1, nj
             do ipn = 1, npn
                if ( .not. allocated( fr%nocc1b(ipn, n1, n2)%m) ) cycle
                nml = lbound( fr%nocc1b(ipn,n1,n2)%m, 1 )
                nmu = ubound( fr%nocc1b(ipn,n1,n2)%m, 1 )
                allocate( op%nocc1b(ipn, n1, n2)%m(nml:nmu) )
                do mm = nml, nmu
                   if ( .not. allocated( fr%nocc1b(ipn, n1, n2)%m(mm)%v ) ) cycle
                   allocate( op%nocc1b(ipn, n1, n2)%m(mm)%v( &
                        size(fr%nocc1b(ipn, n1, n2)%m(mm)%v ) ) )
                   op%nocc1b(ipn, n1, n2)%m(mm)%v(:) &
                        = fr%nocc1b(ipn, n1, n2)%m(mm)%v(:)
                end do
             end do
          end do
       end do
    end if


    if (allocated( fr%nocc2b )) then
       npn = size(fr%nocc2b, 1)
       allocate( op%nocc2b(npn, nj, nj, nj, nj) )
       do n4 = 1, nj
          do n3 = 1, nj
             do n2 = 1, nj
                do n1 = 1, nj
                   do ipn = 1, npn
                      if (.not. allocated(fr%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                      nml = lbound( fr%nocc2b(ipn,n1,n2,n3,n4)%m, 1 )
                      nmu = ubound( fr%nocc2b(ipn,n1,n2,n3,n4)%m, 1 )
                      allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(nml:nmu) )
                      do mm = nml, nmu
                         if ( .not. allocated( fr%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v ) ) cycle
                         allocate( op%nocc2b(ipn, n1, n2, n3, n4)%m(mm)%v( &
                              size(fr%nocc2b(ipn, n1, n2, n3, n4)%m(mm)%v, 1), &
                              size(fr%nocc2b(ipn, n1, n2, n3, n4)%m(mm)%v, 2) ) )
                         op%nocc2b(ipn, n1, n2, n3, n4)%m(mm)%v(:,:) &
                              = fr%nocc2b(ipn, n1, n2, n3, n4)%m(mm)%v(:,:)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if


    op%n_three_body_mp = fr%n_three_body_mp

    if ( allocated( fr%idx_three_body_mp ) ) then
       allocate( op%idx_three_body_mp( &
            size( fr%idx_three_body_mp, 1 ), &
            size( fr%idx_three_body_mp, 2 ), &
            size( fr%idx_three_body_mp, 3 ) ) )
       op%idx_three_body_mp(:,:,:) = fr%idx_three_body_mp(:,:,:)
    end if

    if (allocated( fr%v_three_body_mp )) then
       allocate( op%v_three_body_mp( size(fr%v_three_body_mp, 1) ) )
       op%v_three_body_mp(:) =  fr%v_three_body_mp(:)
    end if

  end subroutine init_copy_opr_m


  subroutine init_operator_mscheme()
    integer :: k1, k2, n1, n2, i, j, m1, m2, n, maxm, maxmd, nj, mm, ns, ipn
    integer, allocatable :: m1inv(:), m2inv(:)

    nj = maxval(n_jorb)
    allocate( idx_nocc2b(2, nj, nj) )
    allocate(m1inv(-maxval(jorb):maxval(jorb)), &
         m2inv(-maxval(jorb):maxval(jorb)))
    m1inv = 0
    m2inv = 0
    n_id_idx(:) = 0
    ! pairing combination, mp = m1+m2  for p-p, n-n int.
    ! density combination, md = m1-m2  for p-n int.
    do k1 = 1, n_jorb_pn
       if (itorb(k1)==-1) then
          ipn = 1
          n1 = k1
       else
          ipn = 2
          n1 = k1 - n_jorb(1)
       end if
       forall(i=1:n_morb(ipn), korbn(i,ipn)==k1) m1inv(morbn(i,ipn)) = i
       do k2 = 1, n_jorb_pn
          if (itorb(k1)/=itorb(k2)) cycle
          if (itorb(k1)==-1) then
             n2 = k2
          else
             n2 = k2 - n_jorb(1)
          end if
          forall(i=1:n_morb(ipn), korbn(i,ipn)==k2) m2inv(morbn(i,ipn)) = i
          maxm = (jorb(k1)+jorb(k2))/2
          if (k1==k2) maxm = maxm - 1
          maxmd = (jorb(k1)+jorb(k2))/2
          allocate( idx_nocc2b(ipn, n1, n2)%mp(-maxm:maxm), &
               idx_nocc2b(ipn, n1, n2)%md(-maxmd:maxmd) )
          do mm = -maxm, maxm
             ns = min(maxm-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             if (k1==k2) ns = (maxm-abs(mm))/2+1
             allocate( idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(2,ns) )
             n_id_idx(ipn+2) = n_id_idx(ipn+2) + 1
             idx_nocc2b(ipn, n1, n2)%mp(mm)%id = n_id_idx(ipn+2)
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = mm*2 - m1
                if (abs(m2)>jorb(k2)) cycle
                if (k1==k2 .and. m1inv(m1)>=m2inv(m2)) cycle
                n = n + 1
                idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(:,n) &
                     = (/ m1inv(m1), m2inv(m2) /)
             end do
             !  do i = 1, ns
             !    write(*,'(a,8i3)') "PP,NN  pairing index",i, ns, ipn, n1, n2, mm, &
             !          & idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(:,i)
             ! end do
             if (n/=ns) stop "error"
          end do
          do mm = -maxmd, maxmd
             ns = min(maxmd-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             allocate( idx_nocc2b(ipn, n1, n2)%md(mm)%idx(2,ns) )
             n_id_idx(ipn) = n_id_idx(ipn) + 1
             idx_nocc2b(ipn, n1, n2)%md(mm)%id = n_id_idx(ipn)
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = m1 - mm*2
                if (abs(m2)>jorb(k2)) cycle
                n = n + 1
                idx_nocc2b(ipn, n1, n2)%md(mm)%idx(:,n) &
                     & = (/ m1inv(m1), m2inv(m2) /)
             end do
!             do i = 1, ns
!               write(*,*) "PP,NN density index",i, ns, ipn, n1, n2, mm, &
!                     & idx_nocc2b(ipn, n1, n2)%md(mm)%idx(:,i)
!             end do
             if (n/=ns) stop "error"
          end do
       end do
    end do

    deallocate( m1inv, m2inv )

  end subroutine init_operator_mscheme


  subroutine opr_m_one_crt(om, k, mm)
    ! one-particle creation operator  orbit: k, jz:mm
    type(opr_m), intent(inout) :: om
    integer, intent(in) :: k, mm
    integer :: i, n, ipn

    ipn = 1
    om%nbody = -1
    n = k
    if (k > n_jorb(1)) then
       ipn = 2
       om%nbody = -2
       n = k - n_jorb(1)
    end if
    om%crt_idx = 0
    do i = 1, n_morb(ipn)
       if (korbn(i,ipn) /= k) cycle
       if (morbn(i,ipn) /= mm) cycle
       om%crt_idx = i
       exit
    end do
    if (om%crt_idx==0) stop "ERROR: opr_m_create"
    om%crt_orb = n
    om%crt_v = 1.d0
    om%ipr1_type = iporb(k)

  end subroutine opr_m_one_crt



  subroutine opr_m_one_anh(om, k, mm)
    ! one-particle annihilation operator  orbit: k, jz:-mm
    ! N.B. c_(j,-m) act as J,+M
    type(opr_m), intent(inout) :: om
    integer, intent(in) :: k, mm
    integer :: i, n, ipn

    ipn = 1
    om%nbody = -6
    n = k
    if (k > n_jorb(1)) then
       ipn = 2
       om%nbody = -7
       n = k - n_jorb(1)
    end if
    om%crt_idx = 0
    do i = 1, n_morb(ipn)
       if (korbn(i,ipn) /= k) cycle
       if (morbn(i,ipn) /= -mm) cycle
       om%crt_idx = i
       exit
    end do
    if (om%crt_idx==0) stop "ERROR: opr_m_anh"
    om%crt_orb = n
    om%crt_v = 1.d0
    om%ipr1_type = iporb(k)
  end subroutine opr_m_one_anh



  subroutine opr_m_two_crt(om, k1, k2, irank)
    ! two-particle creation operator  orbit: k1, k2, irank
    !  1/sqrt(1+delta12) [c+ c+]^(irank)
    type(opr_m), intent(inout) :: om
    integer, intent(in) :: k1, k2, irank
    integer :: i, j, ij, n1, n2, ipn, j1, j2, m1, m2, mp, n, mm1, mm2
    integer :: jpn, kpn
    real(8) :: f
    type(jz_2b_idx), pointer :: noc1, noc2

    if (      k1 <= n_jorb(1) .and. k2 <= n_jorb(1) ) then
       ipn = 1
       om%nbody = -3
       n1 = k1
       n2 = k2
       jpn = 1
       kpn = 1
    else if ( k1 > n_jorb(1) .and. k2 > n_jorb(1) ) then
       ipn = 2
       om%nbody = -4
       n1 = k1 - n_jorb(1)
       n2 = k2 - n_jorb(1)
       jpn = 2
       kpn = 2
    else if ( k1 <= n_jorb(1) .and. k2 > n_jorb(1) ) then
       ipn = 3
       om%nbody = -5
       n1 = k1
       n2 = k2 - n_jorb(1)
       jpn = 1
       kpn = 2
       call init_op_m_gt()
    else
       stop "error in opr_m_two_crt"
    end if

    f = 1.d0
    if (k1 == k2) f = sqrt( 0.5d0 )
    om%irank = irank
    om%ipr1_type = iporb(k1) * iporb(k2)

    allocate( om%nocc1b(3, n_jorb(jpn), n_jorb(kpn)) )
    j1 = jorb(k1)
    j2 = jorb(k2)
    if ( j1+j2 < irank*2 .or. abs(j1-j2) > irank*2 ) stop "ERROR two_crt irank"

    if (ipn == 3) then
       noc1 => idx_gt(jpn, n1, n2)
    else
       noc1 => idx_nocc2b(ipn, n1, n2)
    end if

    if ( .not. allocated(noc1%mp) ) stop "error in two_crt"
    m1 = max( lbound(noc1%mp, 1), -om%irank )
    m2 = min( ubound(noc1%mp, 1),  om%irank )
    allocate( om%nocc1b(ipn, n1, n2)%m( m1:m2 ) )
    do mp = m1, m2
       n = size( noc1%mp(mp)%idx, 2 )
       allocate( om%nocc1b(ipn, n1, n2)%m(mp)%v(n) )
       om%nocc1b(ipn, n1, n2)%m(mp)%v = 0.d0
       do ij = 1, n
          i = noc1%mp(mp)%idx(1, ij)
          j = noc1%mp(mp)%idx(2, ij)
          if (korbn(i, jpn) /= k1 .or. korbn(j, kpn) /= k2) cycle
          mm1 = morbn(i, jpn)
          mm2 = morbn(j, kpn)
          if (k1 == k2) then
             om%nocc1b(ipn, n1, n2)%m(mp)%v(ij) &
                  = f * (  dcg(j1, mm1, j2, mm2, 2*om%irank, 2*mp)  &
                  &      - dcg(j2, mm2, j1, mm1, 2*om%irank, 2*mp) )
          else
             om%nocc1b(ipn, n1, n2)%m(mp)%v(ij) &
                  = f *    dcg(j1, mm1, j2, mm2, 2*om%irank, 2*mp)
          end if

       end do
       if ( maxval(abs(om%nocc1b(ipn, n1, n2)%m(mp)%v)) < 1.d-8 ) &
            deallocate(om%nocc1b(ipn, n1, n2)%m(mp)%v)
    end do

  end subroutine opr_m_two_crt


  subroutine init_op_m_gt()
    ! idx_gt : index of c_p+ c_n, c_n+ c_p
    integer :: k1, k2, n1, n2, i, j, m1, m2, n, maxm, maxmd
    integer :: ni, nj, mm, ns, ipn, inp
    integer, allocatable :: m1inv(:), m2inv(:)

    if ( allocated(idx_gt) ) return

    allocate( idx_gt( 2, maxval(n_jorb), maxval(n_jorb) ) )
    allocate( m1inv(-maxval(jorb):maxval(jorb)), &
         m2inv(-maxval(jorb):maxval(jorb)) )
    m1inv = 0
    m2inv = 0
    ! p-n density combination md = m1-m2
    !     pairing combination mp = m1+m2
    do k1 = 1, n_jorb_pn
       if (itorb(k1)==-1) then
          ipn = 1
          n1 = k1
       else
          ipn = 2
          n1 = k1 - n_jorb(1)
       end if
       inp = 3 - ipn
       do i = 1, n_morb(ipn)
          if (korbn(i,ipn)==k1) m1inv(morbn(i,ipn)) = i
       end do
       do k2 = 1, n_jorb_pn
          if (itorb(k1)==itorb(k2)) cycle
          if (inp == 1) then
             n2 = k2
          else
             n2 = k2 - n_jorb(1)
          end if
          do i = 1, n_morb(inp)
             if (korbn(i,inp)==k2) m2inv(morbn(i,inp)) = i
          end do
          maxm  = (jorb(k1)+jorb(k2))/2
          maxmd = (jorb(k1)+jorb(k2))/2
          allocate( idx_gt(ipn, n1, n2)%mp(-maxm:maxm), &
               idx_gt(ipn, n1, n2)%md(-maxmd:maxmd) )
          do mm = -maxm, maxm
             ns = min(maxm-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             allocate( idx_gt(ipn, n1, n2)%mp(mm)%idx(2, ns) )
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = mm*2 - m1
                if (abs(m2)>jorb(k2)) cycle
                n = n + 1
!                idx_gt(ipn, n1, n2)%mp(mm)%idx(:,n) &
!                     = (/ m1inv(m1), m2inv(m2) /)
                idx_gt(ipn, n1, n2)%mp(mm)%idx(1,n) = m1inv(m1)
                idx_gt(ipn, n1, n2)%mp(mm)%idx(2,n) = m2inv(m2)
             end do
!             do i = 1, ns
!               write(*,*) "cP-dN   index",i, ns, ipn, n1, n2, mm, &
!                     & idx_gt(ipn, n1, n2)%mp(mm)%idx(:,i)
!             end do
             if (n/=ns) stop "error idx_gt mm"
          end do
          do mm = -maxmd, maxmd
             ns = min(maxmd-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             allocate( idx_gt(ipn, n1, n2)%md(mm)%idx(2,ns) )
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = m1 - mm*2
                if (abs(m2)>jorb(k2)) cycle
                n = n + 1
                idx_gt(ipn, n1, n2)%md(mm)%idx(1,n) = m1inv(m1)
                idx_gt(ipn, n1, n2)%md(mm)%idx(2,n) = m2inv(m2)
             end do
             ! do i = 1, ns
             !    write(*,*) "Gamow Teller density index",i, ns, ipn, n1, n2, mm, &
             !         idx_gt(ipn, n1, n2)%md(mm)%idx(:,i)
             ! end do
             if (n/=ns) stop "error idx_gt md"
          end do
       end do
    end do

    deallocate( m1inv, m2inv )

  end subroutine init_op_m_gt



  subroutine add_opr_m(op, v, oa)
    ! op = op + v * oa
    type(opr_m), intent(inout) :: op
    type(opr_m), intent(in) :: oa
    real(8), intent(in) :: v
    integer :: ipn, k1, n1, n2, n3, n4, mm, n, it, nj, nn, npn, ml, mu

    if (op%nbody /= oa%nbody) stop "error add_opr_m rank nbody"
    if (op%irank /= oa%irank) stop "error add_opr_m rank"
    if (op%ipr1_type /= oa%ipr1_type) stop "error add_opr_m ipr1_type"

    nj = maxval(n_jorb)

    if ( allocated( oa%nocc1b ) ) then
       npn = size(oa%nocc1b, 1)
       if (.not. allocated(op%nocc1b)) allocate( op%nocc1b(npn, nj, nj) )
       !$omp parallel do private(nn, n1, n2, ipn, n, mm, ml, mu) &
       !$omp& schedule(dynamic)
       do nn = 0, nj**2-1
          n2 = nn / nj + 1
          n1 = mod(nn, nj) + 1
          do ipn = 1, npn
             if (.not. allocated(oa%nocc1b(ipn,n1,n2)%m)) cycle
             ml = lbound(op%nocc1b(ipn,n1,n2)%m, 1)
             mu = ubound(op%nocc1b(ipn,n1,n2)%m, 1)
             if (.not. allocated(op%nocc1b(ipn,n1,n2)%m)) &
                  allocate( op%nocc1b(ipn,n1,n2)%m(ml:mu) )
             do mm = ml, mu
                if (.not. allocated(oa%nocc1b(ipn,n1,n2)%m(mm)%v)) cycle
                n = size( oa%nocc1b(ipn,n1,n2)%m(mm)%v )
                if (n == 0) cycle
                if (.not. allocated(op%nocc1b(ipn,n1,n2)%m(mm)%v)) then
                   allocate( op%nocc1b(ipn,n1,n2)%m(mm)%v( &
                        size(oa%nocc1b(ipn,n1,n2)%m(mm)%v) ) )
                   op%nocc1b(ipn,n1,n2)%m(mm)%v(:) = 0.d0
                end if
                op%nocc1b(ipn,n1,n2)%m(mm)%v &
                     = op%nocc1b(ipn,n1,n2)%m(mm)%v &
                     + v * oa%nocc1b(ipn,n1,n2)%m(mm)%v
             end do
          end do
       end do
    end if

    do ipn = 1, 2
       if (.not. allocated(oa%spe(ipn)%v)) cycle
       if (.not. allocated(op%spe(ipn)%v)) then
          allocate( op%spe(ipn)%v( size(oa%spe(ipn)%v) ) )
          op%spe(ipn)%v(:) = 0d0
       end if
       op%spe(ipn)%v(:) = op%spe(ipn)%v(:) + v * oa%spe(ipn)%v(:)
    end do

    if ( allocated(oa%nocc2b) ) then
       npn = size(oa%nocc2b, 1)
       if (.not. allocated(op%nocc2b)) allocate( op%nocc2b(npn, nj, nj, nj, nj) )
       !$omp parallel do private(nn, n1, n2, n3, n4, ipn, n, mm, ml, mu) &
       !$omp& schedule(dynamic)
       do nn = 0, nj**4-1
          n4 = nn / nj**3 + 1
          n3 = mod(nn, nj**3) / nj**2 + 1
          n2 = mod(nn, nj**2) / nj + 1
          n1 = mod(nn, nj) + 1
          do ipn = 1, npn
             if (.not. allocated(oa%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
             ml = lbound(oa%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
             mu = ubound(oa%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
             if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m)) then
                if(myrank==0) write(*,*) "WARNING: allocate M",ipn,n1,n2,n3,n4
                allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(ml:mu) )
             end if

             do mm = ml, mu

                if (.not. allocated(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v)) cycle
                if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v) ) then
                   if (myrank==0) write(*,*) &
                        "WARNING: allocate",ipn,n1,n2,n3,n4,mm
                   allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                        size(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v, 1), &
                        size(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v, 2) ))
                   op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
                end if

                op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) &
                     = op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) &
                     + v * oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:)

             end do
          end do
       end do
    end if

  end subroutine add_opr_m



  subroutine opr_m_eff_charge(op, oa, eff_charge)
    ! op = oa(p)*eff_charge(p) + oa(n)*eff_charge(n)
    type(opr_m), intent(inout) :: op
    type(opr_m), intent(in) :: oa
    real(8), intent(in) :: eff_charge(2)
    integer :: ipn, n1, n2, mm, nj, m1, m2

    if (oa%nbody /= 1 .and. oa%nbody/=-10) &
         stop "not implemented opr_m_eff_charge"
    op%nbody = oa%nbody
    op%irank = oa%irank
    op%ipr1_type = oa%ipr1_type

    nj = maxval(n_jorb)
    allocate( op%nocc1b(2, nj, nj) )
    do n2 = 1, nj
       do n1 = 1, nj
          do ipn = 1, 2
             if (.not. allocated(oa%nocc1b(ipn, n1, n2)%m)) cycle
             m1 = lbound(oa%nocc1b(ipn, n1, n2)%m, 1)
             m2 = ubound(oa%nocc1b(ipn, n1, n2)%m, 1)
             allocate( op%nocc1b(ipn, n1, n2)%m(m1:m2) )
             do mm = m1, m2
                if (.not. allocated(oa%nocc1b(ipn, n1, n2)%m(mm)%v)) cycle
                allocate( op%nocc1b(ipn, n1, n2)%m(mm)%v( &
                     size(oa%nocc1b(ipn, n1, n2)%m(mm)%v, 1) ))
                op%nocc1b(ipn, n1, n2)%m(mm)%v &
                     = oa%nocc1b(ipn, n1, n2)%m(mm)%v * eff_charge(ipn)
             end do
          end do
       end do
    end do
  end subroutine opr_m_eff_charge





  subroutine init_tbtd_op(op, mmm, iprty, n_ferm_l, n_ferm_r)
    !
    ! initialize "op" to store <|c+i c+j cl ck|>
    !    for TBTD, Jz=mmm/2, parity=iprty
    !
    !   NOTE : mmm is doubled integer
    !
    type(opr_m), intent(out) :: op
    integer, intent(in) :: mmm, iprty
    integer, intent(in) :: n_ferm_l(2), n_ferm_r(2)
    integer :: ipn
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, mm1, mm2, mu
    integer :: n, n1, n2, n3, n4, n12, k13, mm, nj
    type(jz_2b_idx), pointer :: noc1

    if ( sum(n_ferm_l) == sum(n_ferm_r) .and. &
         abs(n_ferm_l(1) - n_ferm_r(1)) == 1 ) then
       call init_tbtd_pn(op, mmm, iprty, (n_ferm_l(1)-n_ferm_r(1) == 1) )
       return
    else if ( n_ferm_l(1) == n_ferm_r(1) + 2 .and. &
         n_ferm_l(2) == n_ferm_r(2) - 2 ) then
       call init_tbtd_ppnn(op, mmm, iprty)
       return
    else if ( any(n_ferm_l /= n_ferm_r) ) then
       stop "ERROR TBTD not implemented"
    end if

    op%irank = -1  ! full rank

    if (mod(mmm, 2) /= 0) stop "error in init_tbtd_op"
    mu = mmm / 2
    op%mm    = mu
    op%nbody =  5
    op%ipr1_type = iprty
    op%is_j_square = .false.
    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    nj = maxval(n_jorb)


    ! -------------- 1-body operator -----------------------
    allocate( op%nocc1b(2, nj, nj) )

    do k1 = 1, n_jorb_pn
       n1 = k1
       if ( n1 > n_jorb(1) ) n1 = n1 - n_jorb(1)
       ipn = 1
       if (itorb(k1) == 1) ipn = 2
       do k2 = 1, n_jorb_pn
          if (itorb(k1) /= itorb(k2)) cycle
          if (iporb(k1)*iporb(k2) /= iprty) cycle

          n2 = k2
          if ( n2 > n_jorb(1) ) n2 = n2 - n_jorb(1)
          j2 = jorb(k2)
          noc1 => idx_nocc2b(ipn, n1, n2)
          if ( .not. allocated(noc1%md) ) cycle
          mm1 = lbound(noc1%md, 1)
          mm2 = ubound(noc1%md, 1)
          ! only for delta_M = op%mm
          if (mm1>op%mm .or. mm2<op%mm) cycle
          allocate( op%nocc1b(ipn,n1,n2)%m(op%mm:op%mm) )
          n = size(noc1%md(op%mm)%idx, 2)
          allocate(op%nocc1b(ipn,n1,n2)%m(op%mm)%v(n))
          op%nocc1b(ipn, n1, n2)%m(op%mm)%v = 0.d0
       end do
    end do



    allocate( op%nocc2b(3, nj, nj, nj, nj) )
    ! p-p TBME in pairing combination
    ipn = 1
    do n1 = 1, n_jorb(ipn)
       do n2 = n1, n_jorb(ipn)
          if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle

          do n3 = 1, n_jorb(ipn)
             do n4 = n3, n_jorb(ipn)
                if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
                if (iporb(n1)*iporb(n2)*iporb(n3)*iporb(n4) /= iprty) cycle
                mm1 = max( lbound(idx_nocc2b(ipn,n1,n2)%mp,1), &
                     mu  + lbound(idx_nocc2b(ipn,n3,n4)%mp,1) )
                mm2 = min( ubound(idx_nocc2b(ipn,n1,n2)%mp,1), &
                     mu  + ubound(idx_nocc2b(ipn,n3,n4)%mp,1) )

                allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm1:mm2))
                do mm = mm1, mm2
                   allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                        size(idx_nocc2b(ipn,n1,n2)%mp(mm)%idx, 2), &
                        size(idx_nocc2b(ipn,n3,n4)%mp(mm-mu)%idx, 2) ))
                   op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
                end do

             end do
          end do
       end do
    end do

    !  n-n TBME in pairing combination
    ipn = 2
    do n12 = 1, n_jorb(ipn) * n_jorb(ipn)
       n1 = (n12 - 1) / n_jorb(2) + 1
       n2 = mod(n12 - 1, n_jorb(2)) + 1
       if (n1 > n2) cycle
       k1 = n1 + n_jorb(1)
       k2 = n2 + n_jorb(1)
       if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle
       ! do k3 = n_jorb(1)+1, n_jorb_pn  ! fail in Fugaku
       !    n3 = k3 - n_jorb(1)
       !    do k4 = k3, n_jorb_pn
       !       n4 = k4 - n_jorb(1)
       do n3 = 1, n_jorb(2)
          k3 = n3 + n_jorb(1)
          do n4 = n3, n_jorb(2)
             k4 = n4 + n_jorb(1)
             if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
             if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) /= iprty) cycle
             mm1 = max( lbound(idx_nocc2b(ipn,n1,n2)%mp,1), &
                  mu  + lbound(idx_nocc2b(ipn,n3,n4)%mp,1) )
             mm2 = min( ubound(idx_nocc2b(ipn,n1,n2)%mp,1), &
                  mu  + ubound(idx_nocc2b(ipn,n3,n4)%mp,1) )

             allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm1:mm2))
             do mm = mm1, mm2
                allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                     size(idx_nocc2b(ipn, n1, n2)%mp(mm   )%idx, 2), &
                     size(idx_nocc2b(ipn, n3, n4)%mp(mm-mu)%idx, 2) ))
                op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
             end do
          end do
       end do
    end do

    !  p-n TBME in density combination
    ipn = 3
    do k1 = 1, n_jorb(1)
       n1 = k1
       do k3 = 1, n_jorb(1)
          n3 = k3
          if (.not. allocated(idx_nocc2b( 1,n1,n3)%md)) cycle

          ! do k2 = n_jorb(1)+1, n_jorb_pn  ! fail in Fugaku
          !    n2 = k2 - n_jorb(1)
          !    do k4 = n_jorb(1)+1, n_jorb_pn
          !       n4 = k4 - n_jorb(1)
          do n2 = 1, n_jorb(2)
             k2 = n2 + n_jorb(1)
             do n4 = 1, n_jorb(2)
                k4 = n4 + n_jorb(1)

                if (.not. allocated(idx_nocc2b( 2,n2,n4)%md)) cycle
                if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) /= iprty) cycle

                mm1 = max( lbound(idx_nocc2b( 1,n1,n3)%md,1), &
                     mu  - ubound(idx_nocc2b( 2,n2,n4)%md,1) )
                mm2 = min( ubound(idx_nocc2b( 1,n1,n3)%md,1), &
                     mu  - lbound(idx_nocc2b( 2,n2,n4)%md,1) )

                allocate( op%nocc2b( 3,n1,n2,n3,n4)%m(mm1:mm2) )
                do mm = mm1, mm2
                   allocate( op%nocc2b( 3,n1,n2,n3,n4)%m(mm)%v( &
                        size(idx_nocc2b(1, n1, n3)%md(   mm)%idx, 2), &
                        size(idx_nocc2b(2, n2, n4)%md(mu-mm)%idx, 2) ))
                   op%nocc2b( 3,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
                end do

             end do
          end do
       end do
    end do

  end subroutine init_tbtd_op



  subroutine clear_op(op)
    !
    ! clear  "op" to store for the TBTD and OBTD
    !
    type(opr_m), intent(inout) :: op
    integer :: ipn
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, mm1, mm2
    integer :: n, n1, n2, n3, n4, n12, k13, mm, nj, nn

    nj = maxval(n_jorb)

    if (allocated( op%nocc1b )) then
       !$omp parallel do private(n1, n2, ipn, mm) schedule(dynamic)
       do nn = 0, nj**2-1
          n2 = nn / nj + 1
          n1 = mod(nn, nj) + 1
          do ipn = 1, 2
             if (.not. allocated( op%nocc1b(ipn,n1,n2)%m )) cycle
             do mm = lbound( op%nocc1b(ipn,n1,n2)%m, 1 ), &
                  ubound( op%nocc1b(ipn,n1,n2)%m, 1 )
                if (.not. allocated(op%nocc1b(ipn,n1,n2)%m(mm)%v) ) cycle
                op%nocc1b(ipn,n1,n2)%m(mm)%v(:) = 0.d0
             end do
          end do
       end do
    end if

    if (allocated( op%nocc2b )) then
       !$omp parallel do private(n1, n2, n3, n4, ipn, mm) schedule(dynamic)
       do nn = 0, nj**4-1
          n4 = nn / nj**3 + 1
          n3 = mod(nn, nj**3) / nj**2 + 1
          n2 = mod(nn, nj**2) / nj + 1
          n1 = mod(nn, nj) + 1
          do ipn = lbound(op%nocc2b, 1), ubound(op%nocc2b, 1)
             if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
             do mm = lbound( op%nocc2b(ipn,n1,n2,n3,n4)%m, 1 ), &
                  ubound( op%nocc2b(ipn,n1,n2,n3,n4)%m, 1 )
                if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v) ) cycle
                op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
             end do
          end do
       end do
    end if

    if (op%n_three_body_mp /= 0) stop 'error not implemented in clear_op'

  end subroutine clear_op



  subroutine init_idx_tbtd(jl, jr, iprty, nbody, tidx)
    ! index list for OBTD and TBTD
    integer, intent(in) :: jl, jr, iprty, nbody
    type(tbtd_container), intent(inout) :: tidx
    ! integer, parameter :: max_k1234jj=1000000, max_k12jj=1000
    integer, parameter :: max_k1234jj=10000000, max_k12jj=1000
    integer, allocatable :: t12jj(:,:), t1234jj(:,:)
    integer :: ii, jj, npn, k1, k2, k3, k4

    if (jl < 0 .or. jr < 0) return

    ! OBTD index
    allocate( t12jj(  3, max_k12jj  ) )

    ii = 0
    do jj  = abs(jl-jr)/2, (jl+jr)/2

       if (nbody == 5 .or. nbody == 10) then
          call klist1b( 1,           n_jorb(1), &
               &        1,           n_jorb(1))  ! c+p cp
          call klist1b( n_jorb(1)+1, n_jorb_pn, &
               &        n_jorb(1)+1, n_jorb_pn)  ! c+n cn
       elseif (nbody == 11) then
          call klist1b( 1,           n_jorb(1), &
               &        n_jorb(1)+1, n_jorb_pn)  ! c+p cn
       elseif (nbody == 12) then
          call klist1b( n_jorb(1)+1, n_jorb_pn, &
               &        1,           n_jorb(1) ) ! c+n cp
       end if

    end do

    tidx%n12 = ii
    allocate( tidx%k12jj(3, ii),  tidx%v12jj(ii) )
    tidx%k12jj(:, :) = t12jj(:, :ii)
    tidx%v12jj(:)    = 0.d0
    deallocate( t12jj )


    ! TBTD index
    allocate( t1234jj(7, max_k1234jj) )
    ii = 0

    do jj  = abs(jl-jr)/2, (jl+jr)/2
       if (     nbody ==  5) then
          call klist2b( 1, n_jorb(1), 1, n_jorb(1), &
               &        1, n_jorb(1), 1, n_jorb(1) )  ! pppp
          call klist2b( n_jorb(1)+1, n_jorb_pn, n_jorb(1)+1, n_jorb_pn, &
               &        n_jorb(1)+1, n_jorb_pn, n_jorb(1)+1, n_jorb_pn )  ! nnnn
          call klist2b( 1, n_jorb(1), n_jorb(1)+1, n_jorb_pn, &
               &        1, n_jorb(1), n_jorb(1)+1, n_jorb_pn )  ! pnpn
       else if (nbody == 11) then
          call klist2b( 1, n_jorb(1), 1, n_jorb(1), &
               &        1, n_jorb(1), n_jorb(1)+1, n_jorb_pn )   ! pppn
          call klist2b( 1, n_jorb(1),           n_jorb(1)+1, n_jorb_pn, &
               &        n_jorb(1)+1, n_jorb_pn, n_jorb(1)+1, n_jorb_pn )  ! pnnn
       else if (nbody == 12) then
          call klist2b( n_jorb(1)+1, n_jorb_pn, n_jorb(1)+1, n_jorb_pn, &
               &        1, n_jorb(1),           n_jorb(1)+1, n_jorb_pn )  ! nnpn
          call klist2b( 1, n_jorb(1),           n_jorb(1)+1, n_jorb_pn, &
               &        1, n_jorb(1),           1, n_jorb(1) )  ! pnpp
       else if (nbody == 13) then
          call klist2b( 1, n_jorb(1),           1, n_jorb(1), &
               &        n_jorb(1)+1, n_jorb_pn, n_jorb(1)+1, n_jorb_pn )  ! ppnn
       end if
    end do
    tidx%n1234 = ii

    allocate( tidx%k1234jj(7, ii), tidx%v1234jj(ii) )
    tidx%k1234jj(:,:) = t1234jj(:,:ii)
    tidx%v1234jj(:)   = 0.d0
    deallocate(t1234jj)

  contains

    subroutine klist1b(ki1, ki2, kj1, kj2)
      ! output : ii, t12jj
      integer, intent(in) :: ki1, ki2, kj1, kj2
      integer :: k1, k2

      do k1 = ki1, ki2
         do k2 = kj1, kj2

            if (iporb(k1)*iporb(k2) /= iprty) cycle
            if (2*jj > jorb(k1)+jorb(k2)) cycle
            if (2*jj < abs(jorb(k1)-jorb(k2))) cycle

            ii = ii + 1
            if (ii > max_k12jj) stop 'increase max_k12jj'
            t12jj(:, ii) = (/ k1, k2, jj /)

         end do
      end do

    end subroutine klist1b


    subroutine klist2b(ki1, ki2, kj1, kj2, kl1, kl2, km1, km2)
      ! output : ii, t1234jj
      integer, intent(in) :: ki1, ki2, kj1, kj2, kl1, kl2, km1, km2
      integer :: j12, j34, k1, k2, k3, k4

      do k1 = ki1, ki2
         do k2 = kj1, kj2
            if (k1 > k2) cycle
            do j12 = abs(jorb(k1) - jorb(k2))/2, (jorb(k1) + jorb(k2))/2
               if (k1==k2 .and. mod(j12,2)==1) cycle

               do k3 = kl1, kl2
                  do k4 = km1, km2
                     if (k3 > k4) cycle
                     if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) &
                          /=  iprty) cycle

                     do j34 = abs(jorb(k3)-jorb(k4))/2, &
                          (jorb(k3)+jorb(k4))/2

                        if (k3==k4 .and. mod(j34,2)==1) cycle
                        if ( jj < abs(j12-j34) ) cycle
                        if ( jj > j12+j34 ) cycle

                        ii = ii + 1
                        if (ii > max_k1234jj) stop 'increase max_k1234jj'
                        t1234jj(:, ii) = (/ k1, k2, k3, k4, j12, j34, jj /)

                     end do
                  end do
               end do
            end do
         end do
      end do

    end subroutine klist2b

  end subroutine init_idx_tbtd


  subroutine finalize_idx_tbtd(tidx)
    type(tbtd_container), intent(inout) :: tidx

    if ( allocated( tidx%k12jj ) ) deallocate( tidx%k12jj )
    if ( allocated( tidx%v12jj ) ) deallocate( tidx%v12jj )
    if ( allocated( tidx%k1234jj ) ) deallocate( tidx%k1234jj )
    if ( allocated( tidx%v1234jj ) ) deallocate( tidx%v1234jj )

  end subroutine finalize_idx_tbtd


  subroutine get_cpld_tbtd(op, k1, k2, k3, k4, j12, j34, jj, v)
    !
    ! output TBTD from stored M-scheme operator "op"
    ! without factor (-1)^{Ji-Mi} / <Ji,Mi,Jf,-Mf|jj,Mi-Mf>
    !
    !  k1, k2, k3, k4 : orbit number
    !  j12, j34       : two-body coupled-j, not doubled
    !  jj             : rank, not doubled
    !
    !  Ref. transit-2b.pdf
    !
    type(opr_m), intent(in) :: op
    integer, intent(in) :: k1, k2, k3, k4, j12, j34, jj
    real(8), intent(out) :: v
    integer :: nj, mu, n1, n2, n3, n4, j1, j2, j3, j4


    if ( all( op%nbody /= (/ 5, 11, 12, 13 /) ) ) stop 'ERROR in get_cpld_tbtd op%nbody'

    v = 0.d0
    if ( abs(op%mm) > jj ) return

    nj = n_jorb(1)
    mu = op%mm ! not doubled

    n1 = k1
    n2 = k2
    n3 = k3
    n4 = k4
    if (k1 > nj) n1 = k1 - nj
    if (k2 > nj) n2 = k2 - nj
    if (k3 > nj) n3 = k3 - nj
    if (k4 > nj) n4 = k4 - nj
    j1 = jorb(k1)
    j2 = jorb(k2)
    j3 = jorb(k3)
    j4 = jorb(k4)

    if      (k1 <= nj .and. k2 <= nj .and. k3 <= nj .and. k4 <= nj) then
       call pair_sum(1, 1, 1, v)      ! pppp
    else if (k1 >  nj .and. k2 >  nj .and. k3 >  nj .and. k4 >  nj) then
       call pair_sum(2, 2, 2, v)      ! nnnn
    else if (k1 <= nj .and. k2 >  nj .and. k3 <= nj .and. k4 >  nj) then
       call density_sum(3, 1, 2, v)   ! pnpn
    else if (k1 <= nj .and. k2 <= nj .and. k3 >  nj .and. k4 >  nj) then
       call pair_sum(1, 1, 2, v)      ! ppnn
    else if (k1 <= nj .and. k2 <= nj .and. k3 <= nj .and. k4 >  nj) then
       call pppn_sum(v)               ! pppn
    else if (k1 <= nj .and. k2 >  nj .and. k3 >  nj .and. k4 >  nj) then
       call pnnn_sum(v)               ! pppn
    else
       stop 'ERROR in get_cpld_tbtd pn type'
    end if

    if (k1==k2) v = v * sqrt(2.d0)
    if (k3==k4) v = v * sqrt(2.d0)


  contains

    subroutine pair_sum(ipn_v, ipn, jpn, v)
      integer, intent(in) :: ipn_v, ipn, jpn
      real(8), intent(out) :: v
      integer :: mm1, mm2, mm, ij, kl, i, j, k, l, m1, m2, m3, m4

      v = 0.d0
      mm1 = max( lbound(idx_nocc2b(ipn, n1, n2)%mp,1), &
           mu  + lbound(idx_nocc2b(jpn, n3, n4)%mp,1) )
      mm2 = min( ubound(idx_nocc2b(ipn, n1, n2)%mp,1), &
           mu  + ubound(idx_nocc2b(jpn ,n3, n4)%mp,1) )

      do mm = mm1, mm2
         do ij = 1,  size(idx_nocc2b(ipn, n1, n2)%mp(mm)%idx, 2)
            i = idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(1,ij)
            j = idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(2,ij)
            m1 = morbn(i, ipn)
            m2 = morbn(j, ipn)
            if (abs(m1+m2) > 2*j12) cycle

            do kl = 1,  size(idx_nocc2b(jpn, n3, n4)%mp(mm-mu)%idx, 2)
               k = idx_nocc2b(jpn, n3, n4)%mp(mm-mu)%idx(1,kl)
               l = idx_nocc2b(jpn, n3, n4)%mp(mm-mu)%idx(2,kl)
               m3 = morbn(k, jpn)
               m4 = morbn(l, jpn)
               if (abs(m3+m4) > 2*j34) cycle

               v = v + op%nocc2b(ipn_v,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                    * dcg(j1, m1, j2, m2, 2*j12, m1+m2) &
                    * dcg(j3, m3, j4, m4, 2*j34, m3+m4) &
                    * dcg(2*j12, m1+m2, 2*j34, -m3-m4, 2*jj, 2*mu) &
                    * (-1)**(j34 - (m3 + m4)/2)

            end do
         end do
      end do

    end subroutine pair_sum


     subroutine density_sum(ipn_v, ipn, jpn, v)
       integer, intent(in) :: ipn_v, ipn, jpn
       real(8), intent(out) :: v
       integer :: mm1, mm2, mm, ik, jl, i, j, k, l, m1, m2, m3, m4

       v = 0.d0
       mm1 = max( lbound(idx_nocc2b( ipn, n1, n3)%md,1), &
            mu -  ubound(idx_nocc2b( jpn, n2, n4)%md,1) )
       mm2 = min( ubound(idx_nocc2b( ipn, n1, n3)%md,1), &
            mu -  lbound(idx_nocc2b( jpn, n2, n4)%md,1) )

       do mm = mm1, mm2
          do ik = 1, size(idx_nocc2b(ipn, n1, n3)%md(mm)%idx, 2)
             i = idx_nocc2b(ipn, n1, n3)%md(mm)%idx(1, ik)
             k = idx_nocc2b(ipn, n1, n3)%md(mm)%idx(2, ik)
             m1 = morbn(i, ipn)
             m3 = morbn(k, ipn)
             do jl = 1,  size(idx_nocc2b(jpn, n2, n4)%md(mu-mm)%idx, 2)
                j = idx_nocc2b(jpn, n2, n4)%md(mu-mm)%idx(1, jl)
                l = idx_nocc2b(jpn, n2, n4)%md(mu-mm)%idx(2, jl)
                m2 = morbn(j, jpn)
                m4 = morbn(l, jpn)
                if (abs(m1+m2) > 2*j12) cycle
                if (abs(m3+m4) > 2*j34) cycle

                v = v + op%nocc2b(ipn_v, n1, n2, n3, n4)%m(mm)%v(ik, jl) &
                     * dcg(j1, m1, j2, m2, 2*j12, m1+m2) &
                     * dcg(j3, m3, j4, m4, 2*j34, m3+m4) &
                     * dcg(2*j12, m1+m2, 2*j34, -m3-m4, 2*jj, 2*mu) &
                     * (-1)**(j34 - (m3 + m4)/2)

             end do
          end do
       end do
     end subroutine density_sum


     subroutine pppn_sum(v)
       real(8), intent(out) :: v
       integer :: mm1, mm2, mm, ijk, kl, i, j, k, l, m1, m2, m3, m4, md
       real(8) :: x

       v = 0.d0
       mm1 = lbound(op%nocc2b(1,n1,n2,n3,n4)%m, 1)
       mm2 = ubound(op%nocc2b(1,n1,n2,n3,n4)%m, 1)

       do mm = mm1, mm2, 2
          do ijk = 1, size(idx_pph(1, n1, n2, n3)%md(mm)%idx, 2)
             i = idx_pph(1, n1, n2, n3)%md(mm)%idx(1,ijk)
             j = idx_pph(1, n1, n2, n3)%md(mm)%idx(2,ijk)
             k = idx_pph(1, n1, n2, n3)%md(mm)%idx(3,ijk)
             m1 = morbn(i, 1)
             m2 = morbn(j, 1)
             m3 = morbn(k, 1)
             if (abs(m1+m2) > 2*j12) cycle
             md = mm - mu*2

             do kl = 1, size(idx_p(2, n4)%md(md)%idx, 2)
                l = idx_p(2, n4)%md(md)%idx(1, kl)
                m4 = morbn(l, 2)
                if (abs(m3+m4) > 2*j34) cycle

                v = v + op%nocc2b(1,n1,n2,n3,n4)%m(mm)%v(ijk, kl) &
                     * dcg(j1, m1, j2, m2, 2*j12, m1+m2) &
                     * dcg(j3, m3, j4, m4, 2*j34, m3+m4) &
                     * dcg(2*j12, m1+m2, 2*j34, -m3-m4, 2*jj, 2*mu) &
                     * (-1)**(j34 - (m3 + m4)/2)
             end do

          end do
       end do

     end subroutine pppn_sum


     subroutine pnnn_sum(v)
       real(8), intent(out) :: v
       integer :: mm1, mm2, mm, ij, kl, i, j, k, l, m1, m2, m3, m4, md

       v = 0.d0
       mm1 = lbound(op%nocc2b(2,n1,n2,n3,n4)%m, 1)
       mm2 = ubound(op%nocc2b(2,n1,n2,n3,n4)%m, 1)

       do mm = mm1, mm2, 2
          do ij = 1,  size(idx_p(1, n1)%md(mm)%idx, 2)
             i = idx_p(1, n1)%md(mm)%idx(1, ij)
             m1 = morbn(i, 1)
             md = mm - mu*2

             do kl = 1,  size(idx_pph(2, n3, n4, n2)%md(md)%idx, 2)
                k = idx_pph(2, n3, n4, n2)%md(md)%idx(1, kl)
                l = idx_pph(2, n3, n4, n2)%md(md)%idx(2, kl)
                j = idx_pph(2, n3, n4, n2)%md(md)%idx(3, kl)

                m2 = morbn(j, 2)
                m3 = morbn(k, 2)
                m4 = morbn(l, 2)

                if (abs(m1+m2) > 2*j12) cycle
                if (abs(m3+m4) > 2*j34) cycle

                v = v + op%nocc2b(2,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                     * dcg(j1, m1, j2, m2, 2*j12, m1+m2) &
                     * dcg(j3, m3, j4, m4, 2*j34, m3+m4) &
                     * dcg(2*j12, m1+m2, 2*j34, -m3-m4, 2*jj, mu*2) &
                     * (-1)**(j34 - (m3 + m4)/2)
             end do
          end do
       end do

     end subroutine pnnn_sum


  end subroutine get_cpld_tbtd






  subroutine get_cpld_obtd(op, k1, k2, jj, v)
    !
    ! output OBTD from stored M-scheme operator "op"
    ! without factor (-1)^{Ji-Mi} / <Ji,Mi,Jf,-Mf|jj,Mi-Mf>
    !
    !  k1, k2 : orbit number
    !  jj     : rank, not doubled
    !
    !  see transit-2b.pdf
    !
    type(opr_m), intent(in) :: op
    integer, intent(in) :: k1, k2, jj
    real(8), intent(out) :: v
    integer :: n1, n2, nj, mu

    v = 0.d0
    if (op%nbody /= 5 .and. op%nbody /= 11 ) stop 'ERROR in get_cpld_obtd'
    if (abs(op%mm) > 2*jj) return

    nj = n_jorb(1)
    mu = op%mm

    n1 = k1
    n2 = k2
    if (k1 > nj) n1 = k1 - nj
    if (k2 > nj) n2 = k2 - nj

    if      (k1 <= nj .and. k2 <= nj) then  ! pp
       call pp_sum(1, v)
    else if (k1 >  nj .and. k2 >  nj) then  ! nn
       call pp_sum(2, v)
    else if (k1 <= nj .and. k2 >  nj) then  ! pn
       call pn_sum(   v)
    else
       stop 'ERROR in get_cpld_obtd ipn'
    end if

  contains

    subroutine pp_sum(ipn, v)
      integer, intent(in) :: ipn
      real(8), intent(out) :: v
      integer :: ij, i, j, j1, j2, m1, m2

      v = 0.d0
      do ij = 1,  size(idx_nocc2b(ipn, n1, n2)%md(mu)%idx, 2)
         i = idx_nocc2b(ipn, n1, n2)%md(mu)%idx(1, ij)
         j = idx_nocc2b(ipn, n1, n2)%md(mu)%idx(2, ij)
         j1 = jorb(k1)
         j2 = jorb(k2)
         m1 = morbn(i, ipn)
         m2 = morbn(j, ipn)

         v = v + op%nocc1b(ipn,n1,n2)%m(mu)%v(ij) &
              * dcg(j1, m1, j2, -m2, 2*jj, m1-m2) &
              * (-1)**((j2-m2)/2)
      end do
    end subroutine pp_sum


    subroutine pn_sum( v )
      real(8), intent(out) :: v
      integer :: md, i, j, j1, j2, m1, m2, mm1, mm2

      v = 0.d0
      mm1 = max( lbound( idx_p(1, n1)%md, 1), mu*2 + lbound( idx_p(2, n2)%md, 1) )
      mm2 = min( ubound( idx_p(1, n1)%md, 1), mu*2 + ubound( idx_p(2, n2)%md, 1) )

      do md = mm1, mm2, 2
         i = idx_p(1, n1)%md(md     )%idx(1, 1)
         j = idx_p(2, n2)%md(md-mu*2)%idx(1, 1)
         j1 = jorb(k1)
         j2 = jorb(k2)
         m1 = morbn(i, 1)
         m2 = morbn(j, 2)
         if (abs(m1-m2) > 2*jj) cycle

         v = v + op%nocc1b(1, n1, n2)%m(md)%v(1) &
              * dcg(j1, m1, j2, -m2, 2*jj, m1-m2) &
              * (-1)**((j2-m2)/2)

      end do
    end subroutine pn_sum

  end subroutine get_cpld_obtd







  subroutine init_tbtd_pn(op, mmm, iprty, is_pn)
    !
    ! initialize "op" for Gamow-Teller like
    !   to store <|c+i(p) cj(n)|>,
    !            <|c+i(p) c+j(p) cl(p) ck(n)|> and <|c+i(p) c+j(n) cl(n) ck(n)|>
    !    for TBTD, Jz=mm, parity=iprty
    !
    type(opr_m), intent(out) :: op
    integer, intent(in) :: mmm, iprty
    logical, intent(in) :: is_pn  !  cp+ cn-
    integer :: mu, ipn, i, j, k, l, n
    integer :: k1, k2, k3, k4, mm1, mm2, nn1, nn2
    integer :: n1, n2, n3, n4, n12, k13, mm, nj
    type(jz_2b_idx), pointer :: noc1
    integer :: md, maxm, jpn


    op%irank = -1  ! full rank
    if (mod(mmm, 2) /= 0) stop "error in init_tbtd_pn"
    op%mm  = mmm / 2
    mu = op%mm
    if (is_pn) then
       op%nbody = 11
    else
       stop " not implemented init_tbtd_pn"
       op%nbody = 12
    end if
    op%ipr1_type = iprty
    op%is_j_square = .false.


    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    nj = maxval(n_jorb)

    if (.not. allocated(idx_pph)) call init_tbtd_pph_idx()

    allocate( op%nocc2b(2, nj, nj, nj, nj) )

    if (is_pn) then
       ipn = 1  ! p-p-p x n , p x n-n-n
    else
       ipn = 2  ! n-n-n x p , n x p-p-p
    end if
    jpn = 3 - ipn


    ! -------------- 1-body operator -----------------------
    allocate( op%nocc1b(2, nj, nj) )

    do k1 = 1, n_jorb_pn
       if ( k1 <= n_jorb(1) ) then
          if (ipn /= 1) cycle
          n1 = k1
       else
          if (ipn /= 2) cycle
          n1 = k1 - n_jorb(1)
       end if
       do k2 = 1, n_jorb_pn
          if ( k2 <= n_jorb(1) ) then
             if (jpn /= 1) cycle
             n2 = k2
          else
             if (jpn /= 2) cycle
             n2 = k2 - n_jorb(1)
          end if
          if (iporb(k1)*iporb(k2) /= iprty) cycle

          mm1 = max(  lbound(idx_p(ipn,n1)%md,1), &
               mu*2 + lbound(idx_p(jpn,n2)%md,1) )
          mm2 = min(  ubound(idx_p(ipn,n1)%md,1), &
               mu*2 + ubound(idx_p(jpn,n2)%md,1) )

          allocate(op%nocc1b(ipn,n1,n2)%m(mm1:mm2))
          do mm = mm1, mm2, 2
             nn1 = size( idx_p(ipn,n1)%md(mm     )%idx, 2 )
             nn2 = size( idx_p(jpn,n2)%md(mm-mu*2)%idx, 2 )
             if (nn1 == 0 .or. nn2 == 0) cycle
             if (nn1 /= 1 .or. nn2 /= 1) stop 'error occ1b'
             allocate( op%nocc1b(ipn,n1,n2)%m(mm)%v(nn1))
             op%nocc1b(ipn,n1,n2)%m(mm)%v(:) = 0.d0
          end do

       end do
    end do

    ! two-body operator
    ! p-p-p x n
    do n1 = 1, n_jorb(ipn)
      do n2 = n1, n_jorb(ipn)
        do n3 = 1, n_jorb(ipn)
          k1 = n1
          k2 = n2
          k3 = n3
          if (.not. allocated(idx_pph(ipn,n1,n2,n3)%md)) cycle

          do n4 = 1, n_jorb(jpn)
            k4 = n4 + n_jorb(1)
            if (.not. allocated(idx_p(jpn, n4)%md)) cycle
            if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) /= iprty) cycle

            mm1 = max( lbound(idx_pph(ipn,n1,n2,n3)%md,1), &
                mu*2  + lbound(idx_p(  jpn,n4)%md,1) )
            mm2 = min( ubound(idx_pph(ipn,n1,n2,n3)%md,1), &
                mu*2  + ubound(idx_p(  jpn,n4)%md,1) )

            allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm1:mm2))
            do mm = mm1, mm2, 2
              nn1 = size(idx_pph(ipn,n1,n2,n3)%md(mm     )%idx, 2)
              nn2 = size(idx_p(  jpn,n4      )%md(mm-mu*2)%idx, 2)
              if ( nn1 < 1 .or. nn2 < 1) cycle
              allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( nn1, nn2 ) )
              op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
            end do
          end do
        end do
      end do
    end do

    ! p x n-n-n
    do n1 = 1, n_jorb(ipn)
      k1 = n1
      if (.not. allocated(idx_p(ipn, n1)%md)) cycle
      do n2 = 1, n_jorb(jpn)
        do n3 = 1, n_jorb(jpn)
          do n4 = n3, n_jorb(jpn)
            k2 = n2 + n_jorb(1)
            k3 = n3 + n_jorb(1)
            k4 = n4 + n_jorb(1)
            if (.not. allocated(idx_pph(jpn,n3,n4,n2)%md)) cycle

            if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) /= iprty) cycle

            mm1 = max(  lbound(idx_p(  ipn,n1)%md,1), &
                mu*2 + lbound(idx_pph(jpn,n3,n4,n2)%md,1) )
            mm2 = min(  ubound(idx_p(  ipn,n1)%md,1), &
                mu*2 + ubound(idx_pph(jpn,n3,n4,n2)%md,1) )

            allocate(op%nocc2b(jpn,n1,n2,n3,n4)%m(mm1:mm2))
            do mm = mm1, mm2, 2
              nn1 = size(idx_p(  ipn,n1      )%md(mm)%idx, 2)
              nn2 = size(idx_pph(jpn,n3,n4,n2)%md(mm-mu*2)%idx, 2)
              if ( nn1 < 1 .or. nn2 < 1) cycle
              allocate( op%nocc2b(jpn,n1,n2,n3,n4)%m(mm)%v( nn1, nn2 ) )
              op%nocc2b(jpn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
            end do

          end do
        end do
      end do
    end do

    ! two-body operator
    ! p-p-p x n
    !do n1 = 1, n_jorb(ipn)
    !   do n2 = n1, n_jorb(ipn)
    !      do n3 = 1, n_jorb(ipn)
    !         if (.not. allocated(idx_pph(ipn,n1,n2,n3)%md)) cycle

    !         do n4 = 1, n_jorb(jpn)
    !            if (.not. allocated(idx_p(jpn, n4)%md)) cycle
    !            if (iporb(n1)*iporb(n2)*iporb(n3)*iporb(n4) /= iprty) cycle

    !            mm1 = max( lbound(idx_pph(ipn,n1,n2,n3)%md,1), &
    !                 mu*2  + lbound(idx_p(  jpn,n4)%md,1) )
    !            mm2 = min( ubound(idx_pph(ipn,n1,n2,n3)%md,1), &
    !                 mu*2  + ubound(idx_p(  jpn,n4)%md,1) )

    !            allocate(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm1:mm2))
    !            do mm = mm1, mm2, 2
    !               nn1 = size(idx_pph(ipn,n1,n2,n3)%md(mm     )%idx, 2)
    !               nn2 = size(idx_p(  jpn,n4      )%md(mm-mu*2)%idx, 2)
    !               if ( nn1 < 1 .or. nn2 < 1) cycle
    !               allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( nn1, nn2 ) )
    !               op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
    !            end do
    !         end do
    !      end do
    !   end do
    !end do

    !! p x n-n-n
    !do n1 = 1, n_jorb(ipn)
    !   if (.not. allocated(idx_p(ipn, n1)%md)) cycle
    !   do n2 = 1, n_jorb(jpn)
    !      do n3 = 1, n_jorb(jpn)
    !         do n4 = n3, n_jorb(jpn)
    !            if (.not. allocated(idx_pph(jpn,n3,n4,n2)%md)) cycle

    !            if (iporb(n1)*iporb(n2)*iporb(n3)*iporb(n4) /= iprty) cycle

    !            mm1 = max(  lbound(idx_p(  ipn,n1)%md,1), &
    !                 mu*2 + lbound(idx_pph(jpn,n3,n4,n2)%md,1) )
    !            mm2 = min(  ubound(idx_p(  ipn,n1)%md,1), &
    !                 mu*2 + ubound(idx_pph(jpn,n3,n4,n2)%md,1) )

    !            allocate(op%nocc2b(jpn,n1,n2,n3,n4)%m(mm1:mm2))
    !            do mm = mm1, mm2, 2
    !               nn1 = size(idx_p(  ipn,n1      )%md(mm)%idx, 2)
    !               nn2 = size(idx_pph(jpn,n3,n4,n2)%md(mm-mu*2)%idx, 2)
    !               if ( nn1 < 1 .or. nn2 < 1) cycle
    !               allocate( op%nocc2b(jpn,n1,n2,n3,n4)%m(mm)%v( nn1, nn2 ) )
    !               op%nocc2b(jpn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
    !            end do

    !         end do
    !      end do
    !   end do
    !end do

  end subroutine init_tbtd_pn


  subroutine init_tbtd_pph_idx()
    ! construct 2p1h operator index idx_pph, and 1p operator idx_p
    integer :: maxm, nj, ipn, i, j, k, n1, n2, n3, mm, mm1, mm2, n, nk
    integer, allocatable :: mt(:,:,:,:,:,:,:), nmt(:,:,:,:,:)
    integer, parameter :: maxi=2000

    nj = maxval(jorb)

    allocate( idx_pph(2, nj, nj, nj) )
    allocate( idx_p(2, nj) )

    maxm = maxval(jorb)*3
    allocate( nmt(   2, nj, nj, nj, -maxm:maxm), &
         mt(3, maxi, 2, nj, nj, nj, -maxm:maxm) )

    maxm = 0
    nmt(:,:,:,:,:) = 0
    do ipn = 1, 2
       nk = 0
       if (ipn == 2) nk = n_jorb(1)
       do i = 1, n_morb(ipn)
          n1 = korbn(i, ipn) - nk
          do j = i+1, n_morb(ipn)  ! i<j
             n2 = korbn(j, ipn) - nk
             do k = 1, n_morb(ipn)
                n3 = korbn(k, ipn) - nk

                mm = morbn(i, ipn) + morbn(j, ipn) - morbn(k, ipn)
                n = nmt(ipn, n1, n2, n3, mm) + 1
                if (n > maxi) stop 'increase maxi'
                nmt(ipn, n1, n2, n3, mm) = n
                mt(:, n, ipn, n1, n2, n3, mm) = (/ i, j, k /)
                if (abs(mm) > maxm) maxm = abs(mm)

             end do
          end do
       end do
    end do

    do ipn = 1, 2
       do n1 = 1, n_jorb(ipn)
          do n2 = n1, n_jorb(ipn)
             do n3 = 1, n_jorb(ipn)
                mm1 =  1
                mm2 = -1
                do mm = -maxm, maxm, 2
                   n = nmt(ipn, n1, n2, n3, mm)
                   if (n == 0) cycle
                   if (mm1 > mm) mm1 = mm
                   if (mm2 < mm) mm2 = mm
                end do
                if (mm1 > mm2) cycle

                allocate( idx_pph(ipn, n1, n2, n3)%md(mm1:mm2)  )

                do mm = mm1, mm2, 2
                   n = nmt(ipn, n1, n2, n3, mm)
                   if (n == 0) cycle
                   allocate( idx_pph(ipn, n1, n2, n3)%md(mm)%idx(3, n) )
                   idx_pph(ipn,n1,n2,n3)%md(mm)%idx(:, :) &
                        = mt(:, :n, ipn, n1, n2, n3, mm)
                end do

             end do
          end do
       end do
    end do

    deallocate(nmt, mt)


    do ipn = 1, 2
       do n1 = 1, n_jorb(ipn)
          j = jorbn(n1, ipn)
          allocate( idx_p(ipn, n1)%md(-j:j) )
       end do

       do i = 1, n_morb(ipn)
          n1 = korbn(i,ipn)
          if (ipn == 2) n1 = n1 - n_jorb(1)
          mm = morbn(i, ipn)
          allocate( idx_p(ipn, n1)%md(mm)%idx(1,1) )
          idx_p(ipn, n1)%md(mm)%idx(1,1) = i
       end do
    end do

  end subroutine init_tbtd_pph_idx






  subroutine init_tbtd_ppnn(op, mmm, iprty)
    !
    ! initialize "op" for double-beta decay
    !   to store <|c+i(p) c+j(p) cl(p) ck(n)|>
    !       and  <|c+i(p) c+j(n) cl(n) ck(n)|>
    !    for TBTD, Jz=mm, parity=iprty
    !
    !
    type(opr_m), intent(out) :: op
    integer, intent(in) :: mmm, iprty
    integer :: mu, ipn, i, j, k, l, n
    integer :: k1, k2, k3, k4, mm1, mm2, nn1, nn2
    integer :: n1, n2, n3, n4, n12, k13, mm, nj
    type(jz_2b_idx), pointer :: noc1
    integer :: md, maxm, jpn


    op%irank = -1  ! full rank
    if (mod(mmm, 2) /= 0) stop "error in init_tbtd_pn"
    op%mm  = mmm / 2
    mu = op%mm
    op%nbody = 13
    op%ipr1_type = iprty
    op%is_j_square = .false.

    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    nj = maxval(n_jorb)

    allocate( op%nocc2b(1, nj, nj, nj, nj) )

    ipn = 1
    jpn = 3 - ipn

    !  p-p, n-n pairing combination
    do k1 = 1, n_jorb(1)
       n1 = k1
       do k2 = n1, n_jorb(1)
          n2 = k2
          do k3 = n_jorb(1)+1, n_jorb_pn
             n3 = k3 - n_jorb(1)
             do k4 = n_jorb(1)+1, n_jorb_pn
                n4 = k4 - n_jorb(1)

                if (iporb(k1)*iporb(k2)*iporb(k3)*iporb(k4) /= iprty) cycle

                mm1 = max( lbound(idx_nocc2b( 1,n1,n2)%mp,1), &
                     mu  - ubound(idx_nocc2b( 2,n3,n4)%mp,1) )
                mm2 = min( ubound(idx_nocc2b( 1,n1,n2)%mp,1), &
                     mu  - lbound(idx_nocc2b( 2,n3,n4)%mp,1) )

                allocate( op%nocc2b( 1,n1,n2,n3,n4)%m(mm1:mm2) )
                do mm = mm1, mm2
                   allocate( op%nocc2b( 1,n1,n2,n3,n4)%m(mm)%v( &
                        size(idx_nocc2b(1, n1, n2)%mp(   mm)%idx, 2), &
                        size(idx_nocc2b(2, n3, n4)%mp(mu-mm)%idx, 2) ))
                   op%nocc2b( 1,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
                end do

             end do
          end do
       end do
    end do

  end subroutine init_tbtd_ppnn


  subroutine print_operator_mscheme(om)
    !
    ! print m-scheme operator opr_m for debugg
    !
    type(opr_m), intent(in) :: om
    integer :: ipn, m, n1, n2, n3, n4, nj

    if (myrank /= 0) return
    nj = maxval(n_jorb)
    write(*,'(/,a,/)') '### m-scheme matrix elements ###'
    write(*,*) '### one-body matrix elements ###'
    do ipn = 1, 2
       do m = lbound(om%spe(ipn)%v, 1), ubound(om%spe(ipn)%v, 1)
          write(*,'(2i5,1f10.5)') m, ipn, om%spe(ipn)%v(m)
       end do
    end do
    write(*,*) '### index of pairing two-body term ###'
    do ipn = 1, 2
       do n2 = 1, nj
          do n1 = 1, nj
             if (.not. allocated(idx_nocc2b(ipn, n1, n2)%mp)) cycle
             do m = lbound(idx_nocc2b(ipn, n1, n2)%mp, 1), &
                  ubound(idx_nocc2b(ipn, n1, n2)%mp, 1)
                write(*,'(4i3, 1000i4)') ipn, n1,n2, m, &
                     idx_nocc2b(ipn, n1, n2)%mp(m)%idx

             end do
          end do
       end do
    end do

    write(*,*) '### two-body matrix elements ###'
    do ipn = 1, 3
       do n4 = 1, nj
          do n3 = 1, nj
             do n2 = 1, nj
                do n1 = 1, nj
                   if (.not. allocated(om%nocc2b(ipn, n1, n2, n3, n4)%m)) cycle
                   do m = lbound(om%nocc2b(ipn, n1, n2, n3, n4)%m, 1), &
                        ubound(om%nocc2b(ipn, n1, n2, n3, n4)%m, 1)
                      write(*,'(6i3,1000f10.5)') ipn, n1,n2,n3,n4,m, &
                           om%nocc2b(ipn, n1, n2, n3, n4)%m(m)%v
                   end do
                end do
             end do
          end do
       end do
    end do
    write(*,*)

  end subroutine print_operator_mscheme


  subroutine print_tbtd_op(op)
    !
    ! initialize "op" to store for the TBTD, Jz=mm, parity=iprty
    !   NOTE : mm is integer, NOT doubled
    !
    type(opr_m), intent(in) :: op
    integer :: ipn
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, mm1, mm2
    integer :: n, n1, n2, n3, n4, n12, k13, mm, nj

    if (myrank/=0) return

    nj = maxval(n_jorb)

    write(*,*) "OBTD m-scheme m.e."
    do ipn = lbound(op%nocc1b, 1), ubound(op%nocc1b, 1)
       do n2 = 1, nj
          do n1 = 1, nj
             if (.not. allocated(op%nocc1b(ipn,n1,n2)%m)) cycle
             do mm = lbound( op%nocc1b(ipn,n1,n2)%m, 1 ), &
                  ubound( op%nocc1b(ipn,n1,n2)%m, 1 )
                write(*,'(3i3,1000f8.4)') &
                     n1, n2, ipn, op%nocc1b(ipn,n1,n2)%m(mm)%v(:)
             end do
          end do
       end do
    end do


    write(*,*) "TBTD m-scheme m.e."
    do ipn = lbound(op%nocc2b, 1), ubound(op%nocc2b, 1)
       do n4 = 1, nj
          do n3 = 1, nj
             do n2 = 1, nj
                do n1 = 1, nj
                   if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                   do mm = lbound( op%nocc2b(ipn,n1,n2,n3,n4)%m, 1 ), &
                        ubound( op%nocc2b(ipn,n1,n2,n3,n4)%m, 1 )
                      write(*,'(5i3,1000f8.4)') n1, n2, n3, n4, ipn, &
                           op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine print_tbtd_op




  subroutine print_op_1b(op)
    ! 
    type(opr_m), intent(in) :: op
    integer :: ipn, n1, n2, mm, nj, m1, m2

    if (op%nbody /= 1) stop "ERROR print_op_1b"
    if (myrank/=0) return
    write(*,'(/,a,/)') "********** one-body operator ***************"
    write(*,*) "nbody, irank, ipr1_type",  op%nbody, op%irank, op%ipr1_type
    write(*,'(a)') "ipn,n1,n2,m1,m2, mm, vals"
    nj = maxval(n_jorb)
    do ipn = 1, 2
       do n2 = 1, nj
          do n1 = 1, nj
             write(*,'(a,3i3)')"---------",ipn,n1,n2
             if (.not. allocated(op%nocc1b(ipn, n1, n2)%m)) then
                write(*,'(a,3i3)')"skip", ipn, n1, n2
                cycle
             end if
             m1 = lbound(op%nocc1b(ipn, n1, n2)%m, 1)
             m2 = ubound(op%nocc1b(ipn, n1, n2)%m, 1)
             
             if (m1>m2) write(*,'(a,5i3)') "skipm", ipn, n1, n2, m1, m2
             do mm = m1, m2
                if (.not. allocated(op%nocc1b(ipn, n1, n2)%m(mm)%v)) write(*,'(a,4i3)') "not alloc",ipn,n1,n2,mm
                if (.not. allocated(op%nocc1b(ipn, n1, n2)%m(mm)%v)) cycle
                write(*,'(6i3,100f10.5)') ipn, n1, n2, m1, m2, mm, op%nocc1b(ipn, n1, n2)%m(mm)%v(:)
             end do
          end do
       end do
    end do
    write(*,'(/,a,/)') "******************************************"
  end subroutine print_op_1b

end module operator_mscheme
