module model_space
  use constant, only : kwf, kdim

  implicit none
  integer :: n_morb(2), n_morb_pn=0
  integer :: n_j_orbitals(2), n_j_orbitals_pn
  integer :: n_core(2), n_core_pn=-1
  integer :: ki(2), kf(2)
  integer, allocatable :: n_orbitals(:), l_orbitals(:), j_orbitals(:), isospin_orbitals(:), parity_orbitals(:)
  integer, allocatable :: morb(:), korb(:)
  integer, allocatable :: morbn(:,:), korbn(:,:), j_orbitalsn(:,:), parity_orbitalsn(:,:)
  integer :: n_ferm(2), n_ferm_pn=0, mass
  logical :: is_debug=.false.

  integer :: myrank=0, nprocs=1, ierr=0 ! MPI environment
  integer :: nprocs_shift=1, myrank_shift=0, myrank_reduce=0
  integer :: nprocs_reduce=1, nv_shift=0
  integer :: mycomm_shift, mycomm_reduce
  logical :: is_mpi = .false.

  private :: check_orb, set_core_orbit

  integer, parameter :: n_orb_char=8
  character(n_orb_char), allocatable  :: corb(:)

  real(8) :: sum_rule = 0.d0
  integer, private :: n_alloc_vec = 0, max_n_alloc_vec = 0

  ! core orbits
  integer :: n_nljt_core=0
  integer, allocatable :: nljt_core(:,:)

contains
  
   subroutine read_sps(lunsps, verbose)
      !
      !  read definition of the model space from sps file of LUN = lunsps
      !
      integer, intent(in) :: lunsps
      logical, intent(in), optional :: verbose
      character(1) :: c1, c2
      character(1), parameter :: com1 = '!', com2 = '#'
      integer :: ipn, k, m, i, index

      ! skip comment
      read(lunsps,'(2a1)') c1, c2
      do while (c1 == com1 .or. c1 == com2 .or. c2 == com1 .or. c2 == com2) 
         read(lunsps,'(2a1)') c1, c2
      end do
      backspace(lunsps)
      
      ! read data
      read(lunsps,*) (n_j_orbitals(ipn), ipn = 1, 2), (n_core(ipn), ipn = 1, 2)
      n_j_orbitals_pn = n_j_orbitals(1) + n_j_orbitals(2)
      n_core_pn = abs(n_core(1)) + abs(n_core(2))
      write(*, *) "INSIDE read_sps"
      ! write(*, *) "n_j_orbitals(1): ", n_j_orbitals(1), "n_j_orbitals(2): ", n_j_orbitals(2)
      ! write(*, *) "n_core(1): ", n_core(1), "n_core(2): ", n_core(2)

      allocate(n_orbitals(1:n_j_orbitals_pn))
      allocate(l_orbitals(1:n_j_orbitals_pn))
      allocate(j_orbitals(1:n_j_orbitals_pn))
      allocate(isospin_orbitals(1:n_j_orbitals_pn))
      allocate(parity_orbitals(1:n_j_orbitals_pn))

      ki(1) = 1
      kf(1) = n_j_orbitals(1)
      ki(2) = n_j_orbitals(1) + 1
      kf(2) = n_j_orbitals_pn
      do ipn = 1, 2
         do k = ki(ipn), kf(ipn)
            read(lunsps,*) index, n_orbitals(k), l_orbitals(k), j_orbitals(k), isospin_orbitals(k)
            if (index /= k) then
               write(*,'(1a, 1i3, 1a)') 'error [read_sps]: index of', k, &
                  & '-th orbit'
               stop
            end if
            if (isospin_orbitals(k) /= 2*ipn-3) then
               write(*, '(1a, 1i3, 1a, 1i3)') 'error [read_sps]: tz of the', k, &
                  & '-th orbit must be', 2*ipn-3
               stop
            end if
            call check_orb(k, n_orbitals(k), l_orbitals(k), j_orbitals(k))
            parity_orbitals(k) = (-1)**l_orbitals(k)
         end do
      end do

      ! set m-orbit
      n_morb(1:2) = 0
      do k = 1, n_j_orbitals_pn
         ! From 1 to the total number of model space orbitals.
         ! Isospin can be -1 or 1. Thus, ipn can be (-1+3)/2 = 1 or
         ! (1+3)/2 = 2.
         ipn = (isospin_orbitals(k)+3)/2
         ! write(*, *) "ipn: ", ipn
         ! write(*, *) "n_morb(ipn): ", n_morb(ipn)
         n_morb(ipn) = n_morb(ipn) + j_orbitals(k) + 1
         ! write(*, *) "n_morb(ipn): ", n_morb(ipn)
      end do
      n_morb_pn = n_morb(1) + n_morb(2)

      allocate(morb(1:n_morb_pn), korb(1:n_morb_pn), corb(1:n_morb_pn))
      allocate(morbn(maxval(n_morb),2), korbn(maxval(n_morb),2))
      allocate(j_orbitalsn(maxval(n_j_orbitals),2), parity_orbitalsn(maxval(n_j_orbitals),2) )

      i = 0
      do k = 1, n_j_orbitals_pn
         do m = -j_orbitals(k), j_orbitals(k), 2
            i = i + 1
            morb(i) = m
            korb(i) = k
         end do
         call char_orbit(n_orbitals(k), l_orbitals(k), j_orbitals(k), isospin_orbitals(k), corb(k))
      end do

      i = 0
      ipn = 1
      do k = 1, n_j_orbitals_pn
         if (k==n_j_orbitals(1)+1) then
            ipn = 2
            i = n_j_orbitals(1)
         end if
         j_orbitalsn(k-i,ipn) = j_orbitals(k)
         parity_orbitalsn(k-i, ipn) = parity_orbitals(k)
      end do


      i = 0
      ipn = 1
      do k = 1, n_j_orbitals_pn
         if (k==n_j_orbitals(1)+1) then 
            i = 0
            ipn = 2
         end if
         do m = -j_orbitals(k), j_orbitals(k), 2
            i = i + 1
            morbn(i,ipn) = m
            korbn(i,ipn) = k
         end do
      end do

      call set_core_orbit()

      if (present(verbose)) then
         if (.not. verbose) return
      end if
      if (myrank==0) then 
         write(*,*)
         write(*,*)'model space'
         write(*,*)'  k,  n,  l,  j, tz,  p, 2n+l'
         do k = 1, n_j_orbitals_pn
            write(*,'(7i4, 3a, 8a)') k,n_orbitals(k),l_orbitals(k),j_orbitals(k),isospin_orbitals(k),parity_orbitals(k), & 
               & 2*n_orbitals(k)+l_orbitals(k),"   ", corb(k)
         end do
         write(*,*)
      end if
   end subroutine read_sps


  subroutine char_orbit(n, l, j, it, corb)
    integer, intent(in) :: n, l, j, it
    character(n_orb_char), intent(out) :: corb
    character(1), parameter :: cl(0:16) = &
         & (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', &
         'l', 'm', 'n', 'o', 'p', 'q', 'r', 's'/)
    character(2), parameter :: cj(1:40) = &
         & (/'_1', '_2', '_3', '_4', '_5', '_6', '_7', '_8', '_9', '10', &
         &   '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', &
         &   '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', &
         &   '31', '32', '33', '34', '35', '36', '37', '38', '39', '40' /)
    character(1), parameter :: cn(0:16) = (/'0', '1', '2', '3', '4', '5', &
         '6', '7', '8', '9', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)
    character(1), parameter :: cit(-1:1) = (/'p',' ','n'/)

    corb = cit(it)// ' '//cn(n)//cl(l)//cj(j)//'/2'

  end subroutine char_orbit

  subroutine check_orb(k, n, l, j)

    implicit none
    integer, intent(in) :: k, n, l, j

    if (n < 0 .or. l < 0 .or. j < 0) then
       write(*,'(1a,1i3,1a)') &
            & 'error [check_orb]: (n l j) must be non negative for', &
            & k, '-th orbit' 
       stop
    end if

    if (mod(j, 2) /= 1) then
       write(*,'(1a,1i3,1a)') 'error [check_orb]: invalid j for', &
            & k, '-th orbit'
       stop
    end if

    if (abs(j-2*l) /= 1) then
       write(*,'(1a,1i3,1a)')  'error [check_orb]: |2*l-j| must be 1 for', &
            & k, '-th orbit'
       stop
    end if

  end subroutine check_orb


  subroutine set_n_ferm(np, nn, imass)
    integer, intent(in) :: np, nn, imass
    if (n_core_pn==-1) stop "read_sps should be called before set_n_ferm"
    write(*, *) "INSIDE set_n_ferm"
    n_ferm(1) = np
    n_ferm(2) = nn
    n_ferm_pn = sum(n_ferm)
    if (imass==0) then
       ! mass = n_core_pn + np + nn
       mass = n_core_pn 
       if (n_core(1) >= 0) then
          mass = mass + np
       else
          mass = mass - np
       end if
       if (n_core(2) >= 0) then
          mass = mass + nn
       else
          mass = mass - nn
       end if
          
    else
       mass = imass
    end if
    if (myrank==0) write(*,'(a, 2i3, a, i3, a, 2i5)') &
         'N. of valence protons and neutrons = ', np,nn, '   mass=',mass,&
         "   n,z-core ", n_core(1), n_core(2)

   write(*, *) "END set_n_ferm"
  end subroutine set_n_ferm


  subroutine skip_comment(lun)
    integer, intent(in) :: lun
    character(1) :: c1, c2
    character(1), parameter :: com1 = '!', com2 = '#'
    ! skip comment
    read(lun,'(2a1)') c1, c2
    do while (c1 == com1 .or. c1 == com2 .or. c2 == com1 .or. c2 == com2) 
       read(lun,'(2a1)') c1, c2
    end do
    backspace(lun)
  end subroutine skip_comment
  
   subroutine allocate_l_vec(v, n)
      ! Allocate v and add 1 to the n_alloc_vec counter.
      ! Raises error if v is already associated.
      ! Updates the max_n_alloc_vec counter.
      real(kwf), pointer, intent(inout) :: v(:)
      integer(kdim) :: n
      if (associated(v)) stop "error in allocate_l_vec"
      allocate( v(n) )
      if (.not. associated(v)) stop "*** ERROR: memory insufficient ***"
      n_alloc_vec = n_alloc_vec + 1
      ! if (myrank==0) write(*,*)"allocate lanc_vec",n_alloc_vec
      max_n_alloc_vec = max(max_n_alloc_vec, n_alloc_vec)
   end subroutine allocate_l_vec

   subroutine deallocate_l_vec(v)
      ! De-allocate v and subtract 1 from the n_alloc_vec counter.
      ! Raises error if v is not associated.
      real(kwf), pointer, intent(inout) :: v(:)
      if (.not. associated(v)) stop "error in deallocate_l_vec"
      deallocate( v )
      n_alloc_vec = n_alloc_vec - 1
      ! if (myrank==0) write(*,*)"deallocate lanc_vec",n_alloc_vec
   end subroutine deallocate_l_vec

  subroutine print_max_l_vec()
    if (myrank /= 0) return
    write(*,*)
    write(*,*) "maximum num of allocated lanczos vec.", &
         max_n_alloc_vec
    write(*,*) "present num of allocated lanczos vec.", &
         n_alloc_vec
    write(*,*)
  end subroutine print_max_l_vec


  subroutine set_core_orbit()
    integer ::  it, itc, maxosc, nsum, nosc, &
         nc, lc, jc, k, i
    character(n_orb_char) :: corb

    n_nljt_core = 0
    if (n_core(1)<0 .or. n_core(2)<0) return

    allocate( nljt_core(4, sum(n_core)) )
    do it = 1, 2
       if (n_core(it) == 0) cycle
       itc = it*2-3
       maxosc = int(sqrt( dble(n_core(it)) ))+2
       nsum = 0
       outer: do nosc = 0, maxosc
          do nc = 0, nosc/2+1
             lc = nosc - 2*nc
             core: do jc = 2*lc+1, max(2*lc-1, 1), -2
                do k = 1, n_j_orbitals_pn
                   if ( nc == n_orbitals(k) .and. lc  == l_orbitals(k) .and. &
                        jc == j_orbitals(k) .and. itc == isospin_orbitals(k) ) cycle core
                end do
                n_nljt_core = n_nljt_core + 1
                nljt_core(:, n_nljt_core) = (/ nc, lc, jc, itc /)
!                if (myrank==0) write(*,'(a, 4i3)') 'core orbit ', nc,lc,jc,itc
                nsum = nsum + jc + 1
                if (nsum == n_core(it)) exit outer
                if (nsum > n_core(it)) stop "n_core error"
             end do core
          end do
       end do outer
    end do


    if (myrank/=0) return

    write(*,*)
    do it = 1, 2
       if (it==1) write(*,'(a)', advance='no') ' proton  core'
       if (it==2) write(*,'(a)', advance='no') ' neutron core'
       write(*,'(i3,a)', advance='no') n_core(it),', orbit:'
       do i = 1, n_nljt_core
          if (it*2-3 /= nljt_core(4,i)) cycle
          call char_orbit( &
               nljt_core(1,i), &
               nljt_core(2,i), &
               nljt_core(3,i), &
               nljt_core(4,i), corb)
          write(*,'(2a)', advance='no') corb(2:)
       end do
       write(*,*)
    end do

  end subroutine set_core_orbit

end module model_space
