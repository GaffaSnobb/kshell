module operator_jscheme
  !
  ! j-scheme operator
  ! 
  use model_space
  use harmonic_oscillator, only: init_ho_by_mass, init_ho_by_hw, hbar_omega
  implicit none

  private
! definition of structures
  public :: index_jcouple, opr_j  
! global variables
  public :: jcouple, jcouplemax                    
! main procedures
  public :: read_intfile, set_opr_j, get_opr_j_2b, non_diag_ob_2_tbme
! useful for debugging
  public :: print_operator_jscheme 
! for compatibility of old compilers, NOT use from outside
  public :: opr_1j, opr_2j0

  type index_jcouple
     ! 
     ! index for the coupling of [a+_i*a+_j]^J for given J, iprty (parity) 
     ! and ipn (p-p, n-n or p-n)
     !
     ! n : the number of (i,j)
     ! idx(2,n) : index to i and j
     ! idxrev : reverse index (i,j) to n
     ! 
     integer :: n
     integer, allocatable :: idx(:,:)
     integer, allocatable :: idxrev(:,:)     
  end type index_jcouple

  type(index_jcouple), target, allocatable :: jcouple(:,:,:) ! jcouple(J,iprty,ipn)
  integer :: jcouplemax = -1

   type opr_1j
      !
      ! one-body reduced matrix element : v(i,j) = <i||v||j>
      ! 
      ! i and j: j-scheme orbits
      !
      real(8), allocatable :: reduced_matrix_element(:, :)
   end type opr_1j

   type opr_2j0
      !
      ! rank-zero two-body part for given J, iprty (parity) 
      ! and ipn (p-p, n-n or p-n) : v(ij, kl) = <ij|V|kl> 
      !   
      ! ij and kl : index of (i,j)
      !
      real(8), allocatable :: reduced_matrix_element(:,:)
   end type opr_2j0

   type opr_j
      !
      ! onebody : one-body part
      ! p2(J,iprty,ipn) : with J, iprty and ipn (pp, nn, or pn) designated
      ! 
      type(opr_1j) :: onebody
      type(opr_2j0), allocatable :: p2(:,:,:)
      integer :: nbody ! = 1 for 1b only, or = 2 for 2b
      integer :: irank ! rank
      integer :: ipr1_type  ! parity of the operator: = 1 for even, or 2 for odd
   end type opr_j

  ! kini(:) = 1, n_j_orbitals(1)+1    kfin(:) = n_j_orbitals(1), n_j_orbitals_pn
  integer :: kini(2), kfin(2)

contains

   subroutine read_intfile(lun_interaction, operator_j)
      !
      ! read hamiltomian from (LUN = lun_interaction)
      !
      ! traditional shell model calc.
      ! method1 = 0  ... no mass dependence for one-body m.e.
      ! method2 = 0  ... no mass dependence for TBME
      ! method2 = 1  ... (mass/im0)**pwr mass dependence for TBME
      !
      ! method1 =  1  ... no mass dependence, read hw 
      !
      ! no core calc.
      ! method1 = 10 ... one-body dependence hw*(A-1)/A for kinetic energy
      ! method2 = 10 ... V + Trel * two-body dependence hw/A for kinetic energy
      !
      integer, intent(in) :: lun_interaction
      type(opr_j), intent(out) :: operator_j
      integer :: num1, num2, method1, method2
      integer :: i, k1, k2, k3, k4, jj, ipn, irank, nbody
      integer :: isign12, isign34, ij12, ij34, iprty12, iprty34, iprty
      integer :: ipr1_type
      real(8) :: dep_mass1, dep_mass2, v, x, factor
      real(8), parameter :: eps = 1.0d-14
      write(*, *) "INSIDE read_intfile"
      if (jcouplemax < 0) call init_operator_jscheme()
      irank = 0
      nbody = 2
      ipr1_type = 1  ! positive parity
      call init_opj(irank, nbody, ipr1_type, operator_j) ! allocate 

   ! one-body
      call read_firstline_1b(num1, dep_mass1, method1) ! determine mass dependence
      ! write(*, *) "AFTER read_firstline_1b"
      ! write(*, *) "num1: ", num1, "dep_mass1: ", dep_mass1, "method1:", method1
      ! write(*, *) "AIDS!!!"
      do i = 1, num1
         read(lun_interaction,*) k1, k2, v
         ! write(*, *) "k1: ", k1, "k2: ", k2, "v: ", v
         ! detect invalid input
         ! if (k1 > n_j_orbitals_pn .or. k2 > n_j_orbitals_pn .or. j_orbitals(k1) /= j_orbitals(k2)) &
         !    & then
         !    write(*,'(1a, 1i3, 1a, 1i3)') &
         !       &  'error [read_intfile]: invalid 1-body ME  k1 = ', &
         !       & k1, '  k2 = ', k2
         !    stop
         ! end if
         ! if (abs(operator_j%onebody%reduced_matrix_element(k1,k2)) > eps .or. abs(operator_j%onebody%reduced_matrix_element(k2,k1)) > eps) &
         !    & then
         !    write(*,'(1a, 1i3, 1a, 1i3, 1a)') &
         !       & 'error [read_intfile]: 1-body ME  k1 = ', &
         !       & k1, '  k2 = ', k2, ' appeared multiple times'
         !    stop          
         ! end if
         factor = sqrt(dble(j_orbitals(k1)+1))
         ! write(*, *) "j_orbitals(k1): ", j_orbitals(k1)
         ! write(*, *) "factor: ", factor
         operator_j%onebody%reduced_matrix_element(k1, k2) = v * dep_mass1 * factor
         operator_j%onebody%reduced_matrix_element(k2, k1) = v * dep_mass1 * factor
         ! write(*, *) "matrix_element(k1, k2): ", operator_j%onebody%reduced_matrix_element(k1, k2)
      end do
      ! two-body
      call read_firstline_2b(num2, dep_mass2, method2) ! fix mass dependence
      ! write(*, *) "num2: ", num2, "dep_mass2: ", dep_mass2, "method2: ", method2
      do i = 1, num2 ! num2 is the number of two-body matrix elements.
         if (method2==10) then
            read(lun_interaction, *) k1, k2, k3, k4, jj, v, x
            v = v + x * dep_mass2
         else
            ! k1, k2, k3, k4 are the four orbital indices, jj is the ?
            ! and v is the red. mat. elem.
            read(lun_interaction, *) k1, k2, k3, k4, jj, v
            v = v * dep_mass2
         end if
         if (parity_orbitals(k1) == parity_orbitals(k2)) then
            iprty12 = 1
         else
            iprty12 = 2
         end if
         if (parity_orbitals(k3) == parity_orbitals(k4)) then
            iprty34 = 1
         else
            iprty34 = 2
         end if
         if (iprty12 /= iprty34) then
            stop 'error [read_intfile]: an invalid 2b ME (parity) exists'
         end if
         iprty = iprty12
         isign12 = 1
         isign34 = 1
         if (k1 > k2) call swap(k1, k2, jj, isign12)
         if (k3 > k4) call swap(k3, k4, jj, isign34)
         v = v * isign12 * isign34
         if (k2 <= kfin(1) .and. k4 <= kfin(1)) then
            ! Only proton orbitals.
            ipn = 1
         else if (k1 >= kini(2) .and. k3 >= kini(2)) then
            ! Only neutron orbitals.
            ipn = 2
         else if (k1 <= kfin(1) .and. k3 <= kfin(1) .and. k2 >= kini(2) .and. &
            & k4 >= kini(2)) then
            ! Both proton and neutron orbitals.
            ipn = 3
         else
            stop 'error [read_intfile]: an invalid 2b ME (pn) exists'
         end if

         ij12 = jcouple(jj,iprty,ipn)%idxrev(k1,k2)
         ij34 = jcouple(jj,iprty,ipn)%idxrev(k3,k4)
         if (ij12 == 0) then  ! Just a sanity test!
            ! If the orbital indices k1 and k2 are defined, then
            ! idxrev(k1, k2) must have an entry
            write(*,'(2a, 1i3, 1a, 1i3, 1a, 1i3, 1a, 1i3, 1a, 1i3)') &
               & 'error [read_intfile]: no index found for', &
               & ' J = ', jj, ' iprty = ', iprty, ' ipn = ', ipn, &
               & ' k1 = ', k1, ' k2 = ', k2 
            stop
         end if
         if (ij34 == 0) then  ! Just a sanity test!
            write(*,'(2a, 1i3, 1a, 1i3, 1a, 1i3, 1a, 1i3, 1a, 1i3)') &
               & 'error [read_intfile]: no index found for', &
               & ' J = ', jj, ' iprty = ', iprty, ' ipn = ', ipn, &
               & ' k3 = ', k3, ' k4 = ', k4 
            stop
         end if
         operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij12,ij34) = v
         operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij34,ij12) = v
      end do
      ! write(*, *) "parity orbitals: ", parity_orbitals
      ! write(*, *) "kfin: ", kfin
      ! write(*, *) "kfin(1): ", kfin(1)
      ! write(*, *) "kini: ", kini
      ! write(*, *) "kini(1): ", kini(1)
      ! write(*, *) "LOOOOOLLZZZLZLZLZLZLZL"
      ! write(*, *) jcouple(4, 1, 3)%idxrev
      ! do jj = 0, 5
      !    write(*, *) operator_j%p2(jj, 1, 3)%reduced_matrix_element
      ! end do
      ! write(*, *) "jj = 0, iprty = 1, ipn = 3"
      ! write(*, *) jcouple(0, 1, 3)%idxrev
      ! write(*, *) "jj = 1, iprty = 1, ipn = 3"
      ! write(*, *) jcouple(1, 1, 3)%idxrev
      ! write(*, *) "jj = 2, iprty = 1, ipn = 3"
      ! write(*, *) jcouple(2, 1, 3)%idxrev
      ! write(*, *) "jj = 0, iprty = 1, ipn = 1"
      ! write(*, *) jcouple(0, 1, 1)%idxrev
      ! write(*, *) "jj = 1, iprty = 1, ipn = 1"
      ! write(*, *) jcouple(1, 1, 1)%idxrev
      ! write(*, *) "jj = 2, iprty = 1, ipn = 1"
      ! write(*, *) jcouple(2, 1, 1)%idxrev
         
      ! stop
      contains
         ! private subrouitnes used in read_intfile
         subroutine read_firstline_1b(num1, dep_mass1, method1)
            integer, intent(out) :: num1,  method1
            real(8), intent(out) :: dep_mass1
            real(8) :: hw
            write(*, *) "INSIDE read_firstline_1b"
            call skip_comment(lun_interaction)
            read(lun_interaction, *) num1, method1
            ! write(*, *) "INSIDE read_firstline_1b"
            ! write(*, *) "num1: ", num1, "method1: ", method1
            
            if (method1 == 0) then 
               dep_mass1 = 1.d0
            
            else if (method1 == 1) then 
               backspace(lun_interaction)
               read(lun_interaction, *) num1, method1, hw
               dep_mass1 = 1.d0
            
               if (myrank == 0) write(*,'(1a, 1f7.3)') &
                     & 'No mass dependence,  hw from intfile  hw = ', hw
            
                  if (hbar_omega < eps) call init_ho_by_hw(hw)
            
            else if (method1 == 10) then 
               backspace(lun_interaction)
               read(lun_interaction, *) num1, method1, hw
               dep_mass1 = hw*dble(mass-1)/dble(mass)
            
               if (myrank == 0) write(*,'(1a, 1i4, 1a, 1f7.3)') &
                     & 'K.E.: hw*(A-1)/A dependence in 1-body m.e.    A=', &
                     & mass,'  hw=', hw
            
                  if (hbar_omega < eps) call init_ho_by_hw(hw)
            else
               stop 'error [read_intfile]: 1-body parameter not supported'
            end if
            
            call skip_comment(lun_interaction)
         end subroutine read_firstline_1b

         subroutine read_firstline_2b(n, dep, method2)
            integer, intent(out) :: n, method2
            real(8), intent(out) :: dep
            integer :: im0
            real(8) :: pwr, hw
            write(*, *) "INSIDE read_firstline_2b"

            call skip_comment(lun_interaction)
            read(lun_interaction, *) n, method2
            if (method2 == 0) then
               dep = 1.0d0
            else if (method2 == 1) then 
               backspace(lun_interaction)
               read(lun_interaction, *) n, method2, im0, pwr
               dep = (dble(mass)/dble(im0))**pwr
               if (myrank==0) write(*,'(1a, 1i4, 1a, 1f12.8)') &
                  & 'TBME mass dependence (mass/', im0, ')^', pwr
            else if (method2 == 10) then 
               backspace(lun_interaction)
               read(lun_interaction, *) n, method2, hw
               dep = hw/dble(mass)
               if (myrank==0) write(*,'(1a, 1i4, 1a, 1f7.3)') &
                  & 'K.E.: hw/A dependence in 7th column of TBME.  A=' &
                  & ,mass,'  hw=', hw
               if (hbar_omega < eps) call init_ho_by_hw(hw)
            else
               stop 'error [read_intfile]: 2-body parameter not supported'
            end if
            call skip_comment(lun_interaction)
            ! if (hbar_omega < eps) call init_ho_by_mass(1, mass)
         end subroutine read_firstline_2b
         ! end of private subrouitnes used in read_intfile
   end subroutine read_intfile


  subroutine non_diag_ob_2_tbme(operator_j)
    ! one-body nondiagonal matrix elements to two-body int.
    type(opr_j), intent(inout) :: operator_j
    integer :: jj, ipn, iprty, n, ij12, ij34, k1, k2, k3, k4
    real(8) :: v, x, t(n_j_orbitals_pn, n_j_orbitals_pn)
    write(*, *) "INSIDE non_diag_ob_2_tbme"

    if (operator_j%irank/=0 .or. operator_j%nbody/=2) stop 'not in non_diag_ob_2_tbme'
    t(:,:) = operator_j%onebody%reduced_matrix_element(:,:)
    forall(k1 = 1 : n_j_orbitals_pn) t(k1, k1) = 0.d0
    if (maxval(abs(t)) < 1.d-8) return

    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             if (n == 0) cycle
             do ij12 = 1, n
                do ij34 = ij12, n
                   k1 = jcouple(jj,iprty,ipn)%idx(1,ij12)
                   k2 = jcouple(jj,iprty,ipn)%idx(2,ij12)
                   k3 = jcouple(jj,iprty,ipn)%idx(1,ij34)
                   k4 = jcouple(jj,iprty,ipn)%idx(2,ij34)
                   v = 0.d0
                   x = 1.d0 / dble(n_ferm_pn - 1)
                   if (k1==k2) x = x / sqrt(2.d0)
                   if (k3==k4) x = x / sqrt(2.d0)
                   if (k1/=k3 .and. k2==k4) &
                        v = v + x * operator_j%onebody%reduced_matrix_element(k1,k3) / sqrt(j_orbitals(k1)+1.d0)
                   if (k1==k3 .and. k2/=k4) &
                        v = v + x * operator_j%onebody%reduced_matrix_element(k2,k4) / sqrt(j_orbitals(k2)+1.d0)
                   if (mod(j_orbitals(k3)+j_orbitals(k4)-2*jj, 4)==0) x = -x
                   if (k1/=k4 .and. k2==k3) &
                        v = v + x * operator_j%onebody%reduced_matrix_element(k1,k4) / sqrt(j_orbitals(k1)+1.d0)
                   if (k1==k4 .and. k2/=k3) &
                        v = v + x * operator_j%onebody%reduced_matrix_element(k2,k3) / sqrt(j_orbitals(k2)+1.d0)
                   ! if (abs(v)>1.d-8) write(*,*) "non-diag. 1-body to TBME", &
                   !    k1,k2,k3,k4,jj, &
                   !    operator_j%p2(jj,iprty,ipn)%v(ij12,ij34)," + ", &
                   !    v," = ", operator_j%p2(jj,iprty,ipn)%v(ij12,ij34)+v
                   operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij12,ij34) &
                        = operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij12,ij34) + v
                   if (ij12/=ij34) &
                        operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij34,ij12) &
                        = operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij34,ij12) + v
                end do
             end do
          end do
       end do
    end do
    forall(k1=1:n_j_orbitals_pn, k2=1:n_j_orbitals_pn, k1/=k2) operator_j%onebody%reduced_matrix_element(k1,k2) = 0.d0
  end subroutine non_diag_ob_2_tbme


  subroutine set_opr_j(operator_j, func_1, func_2, irank, ipr1_type, nbody)
    !
    ! set operator_j by one-body function func_1, two-body function func_2
    ! input:
    !  func_1(n1,l1,j1,t1, n2,l2,j2,t2) = <nlj_1|| op^irank ||nlj_2> 
    !                                       reduced matrix element
    !  func_2(n1,l1,j1,t1, n2,l2,j2,t2, n3,l3,j3,t3, n4,l4,j4,t4, JJ)
    !     = <nlj_1 nlj_2| op^irank | nlj_3 nlj_4>_JJ   TBME
    !
    !  nbody =  0   copy
    !           1   one-body int. cp+ cp,   cn+ cn
    !           2   two-body int. c+c+cc
    !           5   two-body transition density (init_tbtd_op, container)
    !          10   one-body transition density (init_obtd_beta, container)
    !          11   two-body transition density for cp+ cn type
    !          12   two-body transition density for cn+ cp type
    !          -1   cp+     for s-factor
    !          -2   cn+     for s-factor
    !          -3   cp+ cp+ for 2p s-factor 
    !          -4   cn+ cn+ for 2n s-factor 
    !          -5   cp+ cn+ reserved, NOT yet available
    !          -6   cp      not available 
    !          -7   cn      not available 
    !          -10  cp+ cn  for beta decay
    !          -11  cn+ cp  for beta decay (only in set_ob_channel)
    !          -12  cp+ cp+ cn cn  for 0v-bb decay
    !          -13  cn+ cn+ cp cp  for 0v-bb decay (not yet used)
    !
    real(8), external :: func_1
    real(8), external, optional :: func_2
    integer, intent(in), optional :: irank, ipr1_type, nbody
    type(opr_j), intent(out) :: operator_j
    integer :: iirank, iipr1_type
    integer :: k1, k2, k3, k4, jj, ipn, inbody
    integer :: ij12, ij34, iprty, n
    write(*, *) "INSIDE set_opr_j"
    
    if (jcouplemax < 0) call init_operator_jscheme()
    inbody = 1
    if (present(func_2)) inbody = 2
    if (present(nbody)) inbody = nbody
    iirank = 0      ! default: rank 0
    if (present(irank)) iirank = irank
    iipr1_type = 1  ! default: positive parity
    if (present(ipr1_type)) iipr1_type = ipr1_type
    if (inbody == -12 .or. inbody == -13) then
       call set_opr_j_cpcp_dndn()
       return
    end if
    call init_opj(iirank, inbody, iipr1_type, operator_j) ! allocate 
    ! one-body
    do k1 = 1, n_j_orbitals_pn
       do k2 = 1, n_j_orbitals_pn
          operator_j%onebody%reduced_matrix_element(k1,k2) = func_1(n_orbitals(k1), &
               & l_orbitals(k1), j_orbitals(k1), isospin_orbitals(k1), &
               & n_orbitals(k2), l_orbitals(k2), j_orbitals(k2), isospin_orbitals(k2))
       end do
    end do
    ! two-body
    if (.not. present(func_2)) return
    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             if (n == 0) cycle
             do ij12 = 1, n
                do ij34 = ij12, n
                   k1 = jcouple(jj,iprty,ipn)%idx(1,ij12)
                   k2 = jcouple(jj,iprty,ipn)%idx(2,ij12)
                   k3 = jcouple(jj,iprty,ipn)%idx(1,ij34)
                   k4 = jcouple(jj,iprty,ipn)%idx(2,ij34)
                   operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij12,ij34)  = func_2( &
                        & n_orbitals(k1), l_orbitals(k1), j_orbitals(k1), isospin_orbitals(k1), &
                        & n_orbitals(k2), l_orbitals(k2), j_orbitals(k2), isospin_orbitals(k2), &
                        & n_orbitals(k3), l_orbitals(k3), j_orbitals(k3), isospin_orbitals(k3), &
                        & n_orbitals(k4), l_orbitals(k4), j_orbitals(k4), isospin_orbitals(k4), jj)
                   if (ij12/=ij34) &
                        & operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij34,ij12)  &
                        & = operator_j%p2(jj,iprty,ipn)%reduced_matrix_element(ij12,ij34) 
                end do
             end do
          end do
       end do
    end do

  contains

    subroutine set_opr_j_cpcp_dndn()
      integer :: jcpl, iprty, ipn, ij12, ij34, n, m
      operator_j%irank = iirank
      operator_j%nbody = inbody
      operator_j%ipr1_type = iipr1_type
      write(*, *) "INSIDE set_opr_j_cpcp_dndn"
      ! init_opj 
      allocate(operator_j%p2(0:jcouplemax, 2, 2))
      do iprty = 1, 2
         do jcpl = 0, jcouplemax
            do ipn = 1, 2
               n = jcouple(jcpl, iprty, ipn  )%n
               m = jcouple(jcpl, iprty, 3-ipn)%n
               allocate( operator_j%p2(jcpl, iprty, ipn)%reduced_matrix_element(n, m) )
               operator_j%p2(jcpl, iprty, ipn)%reduced_matrix_element(:,:) = 0.d0
            end do
         end do
      end do

      do jj = 0, jcouplemax
         do iprty = 1, 2
            do ipn = 1, 2
               do ij12 = 1, jcouple(jj, iprty, ipn)%n
                  do ij34 = 1, jcouple(jj, iprty, 3-ipn)%n
                     k1 = jcouple(jj, iprty, ipn  )%idx(1, ij12)
                     k2 = jcouple(jj, iprty, ipn  )%idx(2, ij12)
                     k3 = jcouple(jj, iprty, 3-ipn)%idx(1, ij34)
                     k4 = jcouple(jj, iprty, 3-ipn)%idx(2, ij34)
                     operator_j%p2(jj, iprty, ipn)%reduced_matrix_element(ij12, ij34)  = func_2( &
                          n_orbitals(k1), l_orbitals(k1), j_orbitals(k1), isospin_orbitals(k1), &
                          n_orbitals(k2), l_orbitals(k2), j_orbitals(k2), isospin_orbitals(k2), &
                          n_orbitals(k3), l_orbitals(k3), j_orbitals(k3), isospin_orbitals(k3), &
                          n_orbitals(k4), l_orbitals(k4), j_orbitals(k4), isospin_orbitals(k4), jj)
                  end do
               end do
            end do
         end do
      end do
    end subroutine set_opr_j_cpcp_dndn

  end subroutine set_opr_j

  subroutine print_operator_jscheme(operator_j)
    !
    ! useful for debugging
    !
    type(opr_j) :: operator_j
    integer :: ipn, iprty, jc, num, n, ij12, ij34, numtot
    integer :: k1, k2, k3, k4
    real(8), parameter :: eps = 1.0d-14

    if (myrank==0) write(*,*)
    if (myrank==0) write(*,'(1a)') '### one-body reduced matrix elements ###'
    do k1 = 1, n_j_orbitals_pn
       do k2 = 1, n_j_orbitals_pn
          if (abs(operator_j%onebody%reduced_matrix_element(k1,k2)) > eps) then
             if (myrank==0) write(*,'(3x, 2i4, 1f12.5)') &
                  & k1, k2, operator_j%onebody%reduced_matrix_element(k1,k2)
          end if
       end do
    end do

    if (operator_j%nbody /= 2) return

    numtot = 0
    do ipn = 1, 3
       if (ipn == 1) then
          if (myrank==0) write(*, '(1a)') '### two-body p-p matrix elements ###'
       else if (ipn == 2) then
          if (myrank==0) write(*, '(1a)') '### two-body n-n matrix elements ###'
       else if (ipn == 3) then
          if (myrank==0) write(*, '(1a)') '### two-body p-n matrix elements ###'
       end if

       do jc = 0, jcouplemax
       do iprty = 1, 2
          n = jcouple(jc,iprty,ipn)%n
          if (n > 0) then             
             num = (n*(n+1))/2
             numtot = numtot + num
             if (myrank==0) write(*,'(1a, 1i3, 1a, 1i1, 1a, 1i6, 1a)') &
                  & '* J = ', jc, ' iprty = ', iprty, ' : ', num, &
                  & ' MEs'
             do ij12 = 1, n
                do ij34 = ij12, n
                   k1 = jcouple(jc,iprty,ipn)%idx(1,ij12)
                   k2 = jcouple(jc,iprty,ipn)%idx(2,ij12)
                   k3 = jcouple(jc,iprty,ipn)%idx(1,ij34)
                   k4 = jcouple(jc,iprty,ipn)%idx(2,ij34)
                   if (myrank==0) write(*,'(3x, 4i4, 1f12.5)') &
                        & k1, k2, k3, k4, operator_j%p2(jc,iprty,ipn)%reduced_matrix_element(ij12,ij34)
                end do
             end do
          end if
       end do
       end do
    end do
    if (myrank==0) write(*,'(1a, 1i6)') '# of total 2 body ME =', numtot
  end subroutine print_operator_jscheme


  function get_opr_j_2b(op, k1, k2, k3, k4, jcpl) result (r)
    type(opr_j), intent(in) :: op
    integer, intent(in) :: k1, k2, k3, k4, jcpl
    real(8) :: r
    integer :: isgn12, isgn34, ip12, ij12, ij34, ipn
    write(*, *) "INSIDE get_opr_j_2b"
    r = 0.d0
    if (parity_orbitals(k1)*parity_orbitals(k2) /= parity_orbitals(k3)*parity_orbitals(k4)) return
    ip12 = 1
    if (parity_orbitals(k1) /= parity_orbitals(k2)) ip12 = 2
    if (isospin_orbitals(k1)==-1 .and. isospin_orbitals(k2)==-1 .and. isospin_orbitals(k3)==-1 .and. &
         & isospin_orbitals(k4)==-1) then
       ipn = 1
    else if (isospin_orbitals(k1)==1 .and. isospin_orbitals(k2)==1 .and. isospin_orbitals(k3)==1 &
         & .and. isospin_orbitals(k4)==1) then
       ipn = 2
    else if (isospin_orbitals(k1)==-1 .and. isospin_orbitals(k2)== 1 &
         .and. isospin_orbitals(k3)==-1 .and. isospin_orbitals(k4)== 1) then
       ipn = 3
    else if (isospin_orbitals(k1)== 1 .and. isospin_orbitals(k2)==-1 &
         .and. isospin_orbitals(k3)== 1 .and. isospin_orbitals(k4)==-1) then
       ipn = 3 !sign?
    else
       return
    end if    
    if (.not. allocated(jcouple(jcpl,ip12,ipn)%idxrev)) return
    if (k1 <= k2) then
       ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
       isgn12 = 1
    else
       ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k2,k1)
       isgn12 = (-1)**((j_orbitals(k1)+j_orbitals(k2))/2-jcpl+1)
    end if
    if (k3 <= k4) then
       ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k3,k4)
       isgn34 = 1
    else
       ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k4,k3)
       isgn34 = (-1)**((j_orbitals(k3)+j_orbitals(k4))/2-jcpl+1)
    end if
    if (ij12==0 .or. ij34==0) return
    r = op%p2(jcpl,ip12,ipn)%reduced_matrix_element(ij12,ij34) * isgn12 * isgn34
  end function get_opr_j_2b


!!!!!!!!!!!! private routines !!!!!!!!!!!!!!!

   subroutine init_operator_jscheme()
      !
      ! initialize indices
      !
      integer :: jc_i, jc_f, jinterval, jcpl, loop
      integer :: n
      integer :: k1, k2, k1ini, k1fin, k2ini, k2fin
      integer :: ipn, iprty
      type(index_jcouple), pointer :: jc_p=>null()
      integer :: tmorb(size(morb))
      write(*, *) "INSIDE init_operator_jscheme"

      kini(:) = (/         1, n_j_orbitals(1)+1 /)
      kfin(:) = (/ n_j_orbitals(1), n_j_orbitals_pn   /)
      ! kini(:) = (/         1, n_proton_orbitals + 1 /)
      ! kfin(:) = (/ n_proton_orbitals, n_orbitals   /)
      


      ! determine maximum J 
      tmorb = morb  ! morb is an array of possible jz values.
      n = maxloc(tmorb, 1)   ! Find index of the element of max value along dim 1. Like np.argmax
      jcouplemax = tmorb(n)  ! (usda example) At this point, jcouplemax is 5
      tmorb(n) = minval(tmorb) - 1 ! Why??
      jcouplemax = jcouplemax + maxval(tmorb)  ! (usda example) At this point, jcouplemax is 5 + 5 = 10.
      allocate(jcouple(0:jcouplemax, 2, 3))

    do ipn = 1, 3
       do loop = 1, 2
          if (loop == 2) then
             do iprty = 1, 2
             do jcpl = 0, jcouplemax
                jc_p => jcouple(jcpl,iprty,ipn)
                n = jc_p%n
               !  write(*, *) "n: ", n
                allocate(jc_p%idx(2,n))
                if (ipn == 1 .or. ipn == 2) then
                  ! for ipn = 1, idxrev is n_proton_orbs x n_proton_orbs
                  ! and indices are 1 to n_proton_orbitals (both dims).
                  ! for ipn = 2, idxrev is n_neutron_orbs x n_neutron_orbs
                  ! and indices are n_proton_orbitals + 1 to n_orbitals (both dims).
                   allocate(jc_p%idxrev( kini(ipn):kfin(ipn), kini(ipn):kfin(ipn) ))
                else if (ipn == 3) then
                  ! for ipn = 3, idxrev is n_proton_orbs x n_neutron_orbs
                  ! and indices are 1 to n_proton_orbitals and
                  ! n_proton_orbs + 1 to n_orbitals.
                   allocate(jc_p%idxrev(kini(1):kfin(1),kini(2):kfin(2)))
                end if
                jc_p%idx(:,:) = 0
                jc_p%idxrev(:,:) = 0
             end do
             end do
          end if
          jcouple(:,:,ipn)%n = 0

          if (ipn == 1 .or. ipn == 2) then
             k1ini = kini(ipn)
             k1fin = kfin(ipn)
             k2ini = kini(ipn)
             k2fin = kfin(ipn)
          else if (ipn == 3) then
             k1ini = kini(1)
             k1fin = kfin(1)
             k2ini = kini(2)
             k2fin = kfin(2)                 
          end if
          do k1 = k1ini, k1fin
             do k2 = max(k1, k2ini), k2fin
                if (parity_orbitals(k1) == parity_orbitals(k2)) then
                   iprty = 1
                else
                   iprty = 2
                end if
                if (k1 /= k2) then
                   jc_i = abs(j_orbitals(k1)-j_orbitals(k2))/2
                   jc_f = (j_orbitals(k1)+j_orbitals(k2))/2
                   jinterval = 1
                else
                   jc_i = 0
                   jc_f = j_orbitals(k1) - 1
                   jinterval = 2
                end if
                if (jc_f > jcouplemax) then
                   write(*,'(2a, 1i3, 1a, 1i3)') &
                        & 'error [init_operator_jscheme]: jc_f exceeded jcouplemax', &
                        & '  jc_f = ', jc_f, ' jcouplemax = ', jcouplemax
                   stop
                end if
                do jcpl = jc_i, jc_f, jinterval
                   jc_p => jcouple(jcpl,iprty,ipn)
                   jc_p%n = jc_p%n + 1
                   n = jc_p%n
                   if (loop == 2) then
                      jc_p%idx(:,n) = (/ k1, k2 /)
                      jc_p%idxrev(k1,k2) = n
                   end if
                end do
             end do
          end do

       end do
    end do
  end subroutine init_operator_jscheme

  subroutine swap(k1, k2, jj, isign)
    integer, intent(inout) :: k1, k2
    integer, intent(in) :: jj
    integer, intent(out) :: isign
    integer :: ktmp

    ktmp = k1
    k1 = k2
    k2 = ktmp
    isign = (-1) ** ((j_orbitals(k1)+j_orbitals(k2))/2 - jj + 1)
  end subroutine swap

  subroutine init_opj(irank, nbody, ipr1_type, operator_j)
    ! Allocate space for all operator_j arrays.
    integer, intent(in) :: irank, nbody, ipr1_type
    type(opr_j), intent(inout) :: operator_j
    integer :: ipn, iprty, jcpl, n
    write(*, *) "INSIDE init_opj"
    operator_j%irank = irank
    operator_j%nbody = nbody
    operator_j%ipr1_type = ipr1_type
       
    if (.not. allocated(operator_j%onebody%reduced_matrix_element)) &
         allocate(operator_j%onebody%reduced_matrix_element(n_j_orbitals_pn, n_j_orbitals_pn))
    operator_j%onebody%reduced_matrix_element(:,:) = 0.0d0

    if (nbody == 2) then
       if (.not. allocated(operator_j%p2)) &
            allocate(operator_j%p2(0:jcouplemax,2,3))
       do ipn = 1, 3
          do iprty = 1, 2
             do jcpl = 0, jcouplemax
                n = jcouple(jcpl,iprty,ipn)%n
                if (.not. allocated(operator_j%p2(jcpl,iprty,ipn)%reduced_matrix_element)) &
                     allocate(operator_j%p2(jcpl,iprty,ipn)%reduced_matrix_element(n,n))
                operator_j%p2(jcpl,iprty,ipn)%reduced_matrix_element(:,:) = 0.0d0
             end do
          end do
       end do
    end if
  end subroutine init_opj
  
end module operator_jscheme
