module jscheme
  use constant, only: maxchar, c_no_init, kdim
  use model_space, only: jorb, norb, lorb, n_jorb, n_jorb_pn, &
       myrank, skip_comment, n_ferm_pn
  use partition, only: compare_nocc
  implicit none

  !
  logical :: is_sum_alpha = .true. ! sum
  ! logical :: is_sum_alpha = .false. ! each j-scheme basis

  ! dimension of single-j orbit with v(seniority), J(total ang.)
  ! dim_vj(sj, v, 2J) j^n, J, v 
  integer, allocatable :: dim_jvjj(:,:,:)

  integer, allocatable :: idx_j_n(:,:,:)

  integer, allocatable :: order_c(:,:)

  type type_ptn_j  ! partition for proton (or neutron)
     integer :: n_id   ! number of partition id's
     integer, allocatable :: nocc(:,:)  
  end type type_ptn_j


  ! proton-neutron combined partitions
  type type_ptn_pn_j ! self 
     integer :: n_ferm(2), iprty, mtotal, max_jj
     type(type_ptn_j) :: pn(2)
     integer :: n_pidpn
     ! pidpnM_pid(:,id) = idp, idn
     integer, allocatable :: pidpn_pid(:,:)
     integer(kdim), allocatable :: jdim(:), jdim_acc(:)
  end type type_ptn_pn_j

  character(len=maxchar) :: nl2char=c_no_init

  ! store j-scheme index
  type type_store_top_j
     integer :: n_top = 0, n_store = 0
     integer, allocatable :: occ(:,:), jlist(:,:), jclist(:,:), &
          vlist(:,:), alist(:,:)
     real(8), allocatable :: val_list(:)
  end type type_store_top_j


  integer, parameter :: max_nj_idx = 1000000, max_jc_idx = 1000000

contains

  subroutine init_j_partition(self, lun)
    type(type_ptn_pn_j), intent(out) :: self
    integer, intent(in) :: lun
    integer :: ipn, i, j

    call skip_comment(lun)
    read(lun, *) self%n_ferm(1), self%n_ferm(2), self%iprty

    ! read partition information of proton, neutron sectors
    call skip_comment(lun)
    read(lun, *) self%pn(1)%n_id, self%pn(2)%n_id

    do ipn = 1, 2
       allocate( self%pn(ipn)%nocc( n_jorb(ipn), self%pn(ipn)%n_id ) )
       call skip_comment(lun)
       do i = 1, self%pn(ipn)%n_id
          read(lun, *) j, self%pn(ipn)%nocc(:, i)
          if (i>2) then
             if (compare_nocc( self%pn(ipn)%nocc(:,i-1), &
                  self%pn(ipn)%nocc(:,i) ) /= 1) &
                  stop "error order pp (or nn) partition"
          end if
          if (i/=j) stop "error in partition file" 
       end do
    end do


    ! read partition of p-n combination
    call skip_comment(lun)
    read(lun, *) self%n_pidpn
    allocate( self%pidpn_pid(2, self%n_pidpn) )
    do i = 1, self%n_pidpn
       read(lun, *) self%pidpn_pid(1,i), self%pidpn_pid(2,i)
       if (i>2) then
          if (compare_nocc(self%pidpn_pid(1:2,i-1),  &
               self%pidpn_pid(1:2,i)) /= 1) &
               stop "error order p-n partition"
       end if
    end do

  end subroutine init_j_partition
  


  subroutine set_dim_jvjj()
    integer :: j, v, jj, i, n, maxsj

    maxsj = maxval(jorb)
    allocate( dim_jvjj( 1:maxsj, 0:(maxsj+1)/2, 0:(maxsj/2+1)**2) )  
    dim_jvjj(:,:,:) = 0

    !$omp parallel do private(j, v, jj) schedule(dynamic)
    do j = 1, maxsj, 2
       do v = 0, (j+1)/2
          do jj = (j+1-v)*v, 0, -2
             dim_jvjj(j, v, jj) &
                  = young_t( j+1-v, v,   ((j+1-v)*v-jj)/2 ) &
                  - young_t( j+1-v, v,   ((j+1-v)*v-jj)/2-1 ) &
                  - young_t( j+3-v, v-2, ((j+3-v)*(v-2)-jj)/2 ) &
                  + young_t( j+3-v, v-2, ((j+3-v)*(v-2)-jj)/2-1 )
!             if (dim_jvjj(j, v, jj) /=0) &
!                  write(*,*) 'dim_jvjj',j,v,jj,dim_jvjj(j, v, jj) 

          end do
       end do
    end do
  end subroutine set_dim_jvjj


  subroutine jlist_from_occ(occ, nj_idx, j_idx)
    integer, intent(in) :: occ(:)
    integer, intent(out) :: nj_idx, j_idx(:,:)
    integer :: i, j, k, n, jj, nj_idx0, nc, ic, c1, c2, c3
    integer :: j1, j2, j3
    integer, allocatable :: njlist(:), jlist(:,:)
    integer, allocatable :: j_idx0(:,:)
    integer, parameter :: nnj = 300

    allocate( njlist(n_jorb_pn), jlist(nnj, n_jorb_pn) )

    ! list of possible J in each orbits in a partition
    njlist(:) = 0
    do i = 1, n_jorb_pn
       n = occ(i)
       j = jorb(i)
       do jj = mod(n,2), (j+1-n)*n, 2
          njlist(i) = njlist(i) + 1
          if ( njlist(i) > size(jlist,1) ) stop 'increase size of jlist'
          jlist( njlist(i), i ) = jj
       end do
    end do

    nj_idx = product( njlist )
    if ( nj_idx > size(j_idx, 2) ) stop 'increase size of j_idx '
    allocate( j_idx0( size(j_idx, 1), size(j_idx, 2) ) )
    ! (J1, J2, ... Jn)
    do k = 1, nj_idx
       n = k - 1
       do i = n_jorb_pn, 1, -1
          j_idx(i, k) = jlist( mod(n, njlist(i))+1, i )
          n = n / njlist(i)
       end do
    end do

    deallocate( njlist, jlist )
    
  end subroutine jlist_from_occ


  subroutine jclist_from_jlist(jlist, jtot, &
       njc_idx, jc_idx )
    integer, intent(in)  :: jtot, jlist(:)
    integer, intent(out) :: njc_idx, jc_idx(:,:)
    integer :: i, j, k, n, jj, njc, nc, ic, c1, c2, c3
    integer :: j1, j2, j3
    integer, allocatable :: jc_idx0(:,:)

    allocate( jc_idx0( size(jc_idx,1), size(jc_idx,2) ) )
    njc_idx = 1
    
    ! (J1, J2, ... Jn, J12, J(12)3, ... , Jtot)
    do nc = 1, n_jorb_pn-1
       njc = njc_idx
       if (nc/=1) jc_idx0(:nc-1, :njc) = jc_idx(:nc-1, :njc)
       c1 = order_c(1, nc)
       c2 = order_c(2, nc)
       c3 = order_c(3, nc) - n_jorb_pn
       if (c1 <= n_jorb_pn) j1 = jlist(c1)
       if (c2 <= n_jorb_pn) j2 = jlist(c2)
       if (nc /= c3) stop 'error coupling order'

       njc_idx = 0
       do i = 1, njc
          if (c1 > n_jorb_pn) j1 = jc_idx0(c1-n_jorb_pn, i)
          if (c2 > n_jorb_pn) j2 = jc_idx0(c2-n_jorb_pn, i)
          do j3 = abs(j1-j2), j1+j2, 2
             if ( nc == n_jorb_pn-1 .and. j3 /= jtot ) cycle
             njc_idx = njc_idx + 1
             if (njc_idx > size(jc_idx,2)) stop "increase size of jc_idx"
             if (nc>1) jc_idx(:nc-1, njc_idx) = jc_idx0(:nc-1, i)
             jc_idx(nc, njc_idx) = j3
          end do

       end do
    end do

    deallocate( jc_idx0 )

    if (njc_idx > 1) call qsort_jclist( jc_idx, 1, njc_idx )

!    do i = 1, njc_idx
!       write(*,'(a,100i3)') "jclist ", jc_idx(:,i)
!    end do

  contains
   
    recursive subroutine qsort_jclist(jclist, left, right)
      ! quick sort of jclist
      integer, intent(inout) :: jclist(:,:)
      integer, intent(in) :: left, right
      integer :: i, j, pvt(size(jclist,1)), t(size(jclist,1))

      pvt = jclist(:, (left+right)/2 )
      i = left
      j = right
      do while (.true.)
         do while ( compare_nocc( jclist(:,i), pvt ) == 1)
            i = i + 1
         end do
         do while ( compare_nocc( pvt, jclist(:,j) ) == 1)
            j = j - 1
         end do
         if (i >= j) exit
         t = jclist(:,i)
         jclist(:,i) = jclist(:,j)
         jclist(:,j) = t
         i = i + 1
         j = j - 1
      end do
      if (left < i-1  ) call qsort_jclist(jclist, left, i-1)
      if (j+1  < right) call qsort_jclist(jclist, j+1,  right)

    end subroutine qsort_jclist
    
  end subroutine jclist_from_jlist
  

  subroutine v_a_from_occ_jlist(occ, jlist, nvlist, vlist, dalist)
    ! seniority and dimension of alpha 
    integer, intent(in) :: occ(:), jlist(:)
    integer, intent(out) :: nvlist(:), vlist(:,:), dalist(:,:)
    integer :: i, n, j, jj, v, d

    nvlist(:) = 0
    do i = 1, n_jorb_pn
       n = occ(i)
       j = jorb(i)
       jj = jlist(i)
       do v = mod(n, 2), min(n,j+1-n), 2
          if ( dim_jvjj(j, v, jj) == 0 ) cycle
          nvlist(i) = nvlist(i) + 1
          if ( nvlist(i) > size(vlist, 1) ) stop "increas size of nvlist, dalist"
          vlist(nvlist(i), i) = v
          dalist(nvlist(i), i) = dim_jvjj(j, v, jj)
       end do
    end do
    
  end subroutine v_a_from_occ_jlist
  


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



  subroutine set_order_coupling()
    integer :: i, n, orbp, orbn, c

    allocate( order_c(3, n_jorb_pn-1) )
    c = 0    
    
    orbp = 1
    do i = 2, n_jorb(1)
       c = c + 1
       order_c(1, c) = orbp
       order_c(2, c) = i
       order_c(3, c) = c + n_jorb_pn
       orbp = n_jorb_pn + c
    end do

    orbn = n_jorb(1)+1
    do i = n_jorb(1)+2, n_jorb_pn
       c = c + 1
       order_c(1, c) = orbn
       order_c(2, c) = i
       order_c(3, c) = c + n_jorb_pn
       orbn = n_jorb_pn + c
    end do

    c = c + 1
    order_c(1, c) = orbp
    order_c(2, c) = orbn
    order_c(3, c) = c + n_jorb_pn

    if (myrank/=0) return
    write(*,*)
    write(*,*) ' *** coupling order ***'
    do i = 1, c
       write(*,'(i3,a,i3,a,i3)') &
            order_c(1,i),' x ', order_c(2,i), ' = ', order_c(3,i)
    end do
    write(*,*)

  end subroutine set_order_coupling


  subroutine print_jscheme_index(occ, jlist, jclist, vlist, alist, val)
    integer, intent(in) :: occ(:), jlist(:), jclist(:), &
         vlist(:), alist(:)
    real(8), intent(in) :: val
    integer :: i
    character :: c
    character(len=9) :: lc = 'spdfghijk'

    if (nl2char == c_no_init) then
       nl2char = ''
       do i = 1, n_jorb_pn
          nl2char( 3*i:3*i ) = lc(lorb(i)+1:lorb(i)+1)
          write(nl2char( 3*i-1:3*i-1 ), '(i1)') norb(i)
       end do
    end if
    
    write(*,'(a,100i3)') "index: ", (/( i, i=1, n_jorb_pn )/)
    write(*,'(a,a)')     "nl   : ", trim(nl2char)
    write(*,'(a,100i3)') "2j   : ", jorb
    write(*,'(a,100i3)') "n    : ", occ
    write(*,'(a,100i3)') "2J   : ", jlist
    write(*,'(a,100i3)', advance='no') "Jcoup:    ", &
         jclist(:n_jorb(1)-1)
    write(*,'(a,100i3)') "   ", jclist(n_jorb(1):)
    write(*,'(a,100i3)') "v    : ", vlist
    if (is_sum_alpha) then
       write(*,'(a,f10.7)') "val^2: ", val
    else
       write(*,'(a,100i3)') "alpha: ", alist
       write(*,'(a,f10.7)') "value: ", val
       write(*,'(a,f10.7)') "val^2 : ", val**2
    end if
    write(*,*)

  end subroutine print_jscheme_index



  subroutine init_store_top_j(self, ntop)
    type(type_store_top_j), intent(inout) :: self
    integer, intent(in) :: ntop
    integer :: i
    
    self%n_store = 0
    if (self%n_top == ntop) return
    self%n_top = ntop
    
    if (allocated( self%jlist)) deallocate( self%val_list, &
         self%occ, self%jlist, self%vlist, self%alist )

    allocate( &
         self%val_list(self%n_top), &
         self%occ(   n_jorb_pn,   self%n_top), &
         self%jlist( n_jorb_pn,   self%n_top), &
         self%jclist(n_jorb_pn-1, self%n_top), &
         self%vlist( n_jorb_pn,   self%n_top ), &
         self%alist( n_jorb_pn,   self%n_top ) )
    
  end subroutine init_store_top_j

  
  subroutine store_top_j(self, val, occ, jlist, jclist, vlist, alist)
    type(type_store_top_j), intent(inout) :: self
    real(8), intent(in) :: val
    integer, intent(in) :: occ(:), jlist(:), jclist(:), vlist(:)
    integer ,intent(in), optional :: alist(:)
    integer :: i, n, ns1, ns2

    if (self%n_store == self%n_top) then
       if (present(alist)) then
          if (val**2 < self%val_list(self%n_store)**2) return
       else
          if (val < self%val_list(self%n_store)) return
       end if
    end if

    n = 0 
    do i = self%n_store, 1, -1
       n = i
       if (present(alist)) then
          if ( val**2 < self%val_list(i)**2 ) exit
       else
          if ( val < self%val_list(i) ) exit
       end if
       n = 0
    end do

    n = n + 1
    if (self%n_store == self%n_top) then
       ns1 = self%n_store
       ns2 = self%n_store - 1
    else
       ns1 = self%n_store + 1
       ns2 = self%n_store
    end if

    if ( n /= self%n_store + 1) then 
       self%val_list(n+1:ns1) = self%val_list(n:ns2)
       self%occ(  :, n+1:ns1) = self%occ(  :, n:ns2)
       self%jlist(:, n+1:ns1) = self%jlist(:, n:ns2)
       self%jclist(:,n+1:ns1) = self%jclist(:,n:ns2)
       self%vlist(:, n+1:ns1) = self%vlist(:, n:ns2)
       if (present(alist)) self%alist(:, n+1:ns1) = self%alist(:, n:ns2)
    end if

    self%val_list(n) = val
    self%occ(  :,n) = occ
    self%jlist(:,n) = jlist
    self%jclist(:,n)= jclist
    self%vlist(:,n) = vlist
    if (present(alist)) self%alist(:,n) = alist

    self%n_store = ns1
    
  end subroutine store_top_j


  subroutine print_top_j(self)
    type(type_store_top_j), intent(in) :: self
    integer :: i

    do i = 1, self%n_store
       call print_jscheme_index(self%occ(:,i), &
            self%jlist(:,i), self%jclist(:,i), &
            self%vlist(:,i), &
            self%alist(:,i), self%val_list(i) &
            )
    end do
    
  end subroutine print_top_j


  subroutine jscheme_dim_ptn(self, jtotal)
    type(type_ptn_pn_j), intent(inout) :: self
    integer, intent(in) :: jtotal
    integer :: id
    integer(kdim) :: n
    
    if (.not. allocated( self%jdim )) then
       allocate(self%jdim(self%n_pidpn))
       allocate(self%jdim_acc(0 : self%n_pidpn))
    end if

    !$omp parallel do
    do id = 1, self%n_pidpn
       call inner(id, self%jdim(id))
    end do
    
    self%jdim_acc(0) = 0
    do id = 1, self%n_pidpn
       self%jdim_acc(id) = self%jdim_acc(id-1) + self%jdim(id)
    end do

    write(*,'(a,i3,i15,/)') 'j-scheme dim', jtotal, self%jdim_acc(self%n_pidpn)

  contains
    
    subroutine inner(id, dim)
      integer, intent(in) :: id
      integer(kdim), intent(out) :: dim
      integer :: ij_idx, nj_idx, njc_idx
      integer :: i, n, iiv, ijc_idx
      integer :: occs(n_jorb_pn), dan(n_jorb_pn), nvlist(n_jorb_pn)
      integer, allocatable :: j_idx(:, :), jc_idx(:,:), &
           vlist(:,:), dalist(:,:)

      allocate( j_idx( n_jorb_pn,   max_nj_idx) )
      allocate( jc_idx(n_jorb_pn-1, max_jc_idx), &
           vlist(max_jc_idx, n_jorb_pn), &
           dalist(max_jc_idx, n_jorb_pn) )

      dim = 0
      occs(:n_jorb(1))   = self%pn(1)%nocc( :, self%pidpn_pid(1, id) )
      occs(n_jorb(1)+1:) = self%pn(2)%nocc( :, self%pidpn_pid(2, id) )
      
      call jlist_from_occ( occs, nj_idx, j_idx )

      do ij_idx = 1, nj_idx

         call jclist_from_jlist( j_idx(:, ij_idx),  &
              jtotal, njc_idx, jc_idx)
         if (njc_idx == 0) cycle

         call v_a_from_occ_jlist(occs, &
              j_idx(:, ij_idx), nvlist, vlist, dalist)

         do ijc_idx = 1, njc_idx 
            do iiv = 1, product(nvlist)
               n = iiv - 1
               do i = n_jorb_pn, 1, -1
                  dan(i) = dalist( mod(n, nvlist(i))+1, i)
                  n = n / nvlist(i)
               end do
               dim = dim + product(dan)
            end do
         end do
      end do
    end subroutine inner

  end subroutine jscheme_dim_ptn



  subroutine check_vec_ptn(self, jtotal, id, vec, storage_top_j, rto_v, rto_JpJn)
    type(type_ptn_pn_j), intent(in) :: self
    integer, intent(in) :: jtotal, id
    real(8), intent(in) :: vec(:)
    type(type_store_top_j), intent(inout) :: storage_top_j
    real(8), intent(out) :: rto_v(:), rto_JpJn(:,:)
    integer :: ij_idx, nj_idx, njc_idx, iida
    integer :: i, n, iiv, ijc_idx, nv
    integer :: vj(n_jorb_pn), dan(n_jorb_pn), &
         nvlist(n_jorb_pn), aj(n_jorb_pn), occs(n_jorb_pn)
    integer, allocatable :: j_idx(:, :), jc_idx(:,:), &
         vlist(:,:), dalist(:,:)
    integer(kdim) :: ndim, Jp, Jn
    real(8) :: val

    rto_v(:) = 0.d0
    rto_JpJn(:,:) = 0d0
    
    occs(:n_jorb(1))   = self%pn(1)%nocc( :, self%pidpn_pid(1, id) )
    occs(n_jorb(1)+1:) = self%pn(2)%nocc( :, self%pidpn_pid(2, id) )
    ! write(*,'(a,100i3)') "occ    ", occs
    
    allocate( j_idx( n_jorb_pn,   max_nj_idx) )
    allocate( jc_idx(n_jorb_pn-1, max_jc_idx), &
         vlist(max_jc_idx, n_jorb_pn), &
         dalist(max_jc_idx, n_jorb_pn) )

    call jlist_from_occ( occs, nj_idx, j_idx )
!    ndim = self%jdim_acc(id-1) 
    ndim = self%jdim_acc(id)  - self%jdim(id)

    do ij_idx = 1, nj_idx

       call jclist_from_jlist( j_idx(:, ij_idx),  &
            jtotal, njc_idx, jc_idx)
       if (njc_idx == 0) cycle

       call v_a_from_occ_jlist(occs, &
            j_idx(:, ij_idx), nvlist, vlist, dalist)
           
       do ijc_idx = 1, njc_idx 

          do iiv = 1, product(nvlist)
             n = iiv - 1
             do i = n_jorb_pn, 1, -1
                vj(i) = vlist( mod(n, nvlist(i))+1, i)
                dan(i) = dalist( mod(n, nvlist(i))+1, i)
                n = n / nvlist(i)
             end do

             nv = sum(vj)
             if (is_sum_alpha) then
                n = product(dan)
                val = sum( vec(ndim+1:ndim+n) ** 2 )
                ndim = ndim + n 
                rto_v(nv+1) = rto_v(nv+1) + val
                
                !$omp critical
                call store_top_j(storage_top_j, val, &
                     occs, j_idx(:, ij_idx), jc_idx(:,ijc_idx), &
                     vj )
                !$omp end critical
                Jp = jc_idx( n_jorb(1) - 1, ijc_idx ) + 1
                Jn = jc_idx( n_jorb_pn - 2, ijc_idx ) + 1
                rto_JpJn(Jp, Jn) = rto_JpJn(Jp, Jn) + val

                cycle
             end if

             
             do iida = 1, product(dan)
                n = iida - 1
                do i = n_jorb_pn, 1, -1
                   aj(i) = mod(n, dan(i)) + 1
                   n = n / dan(i)
                end do

                ndim = ndim + 1
                val = vec(ndim)
                rto_v(nv+1) = rto_v(nv+1) + val**2
                
                !$omp critical
                ! if (val**2> 0.0001) then
                !    call print_jscheme_index(occs, j_idx(:, ij_idx), &
                !         vj, aj, val)
                ! endif
                call store_top_j(storage_top_j, val, &
                     occs, j_idx(:, ij_idx), jc_idx(:,ijc_idx), &
                     vj, aj )
                !$omp end critical

             end do
          end do
       end do
    end do

  end subroutine check_vec_ptn
  

end module jscheme





program print_jwav
! #ifdef MPI
!   use mpi
! #endif
  !$ use omp_lib, only : omp_get_max_threads
  use constant, only: kwf, kdim, maxchar, c_no_init, max_n_jorb
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_jorb_pn, n_jorb, n_ferm, n_core, is_debug, &
       jorb, lorb, korb, norb, itorb, iporb, is_mpi, n_ferm_pn
  use model_space, only: m_mass=>mass, print_max_l_vec, &
       nprocs_reduce, nprocs_shift, nv_shift, is_mpi, allocate_l_vec, deallocate_l_vec
  use jscheme
  implicit none
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn_int, fn_nml, fn_ptn, fn_wav
  integer :: n_eig, iargc, mtot, i, id, istate, n
  type(type_ptn_pn_j) :: ptn
  real(8), allocatable :: e_val(:)
  integer, allocatable :: jj_val(:)
  type(type_store_top_j) :: storage_top_j
  real(8), allocatable :: vec(:)

  
! #ifdef MPI
!   write(*,*) ' not yet implemented MPI'
!   is_mpi = .true.
!   call mpi_init(ierr)
!   call mpi_comm_size(mpi_comm_world, nprocs, ierr)
!   call mpi_comm_rank(mpi_comm_world, myrank, ierr)
!   if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') &
!        "nprocs",nprocs,"    myrank", myrank
!   if (myrank/=0) is_debug = .false.
! #endif
  !$ if (myrank==0) write(*,'(1a,1i3,/)') "OpenMP  # of threads=", omp_get_max_threads()

  if (iargc() < 3) &
       stop 'usage: print_jwav.exe foo.snt bar.ptn fn.jwv'

  call getarg(1, fn_int)
  call getarg(2, fn_ptn)
  call getarg(3, fn_wav)
  
! #ifdef MPI
!   call mpi_bcast(fn_int,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
!   call mpi_bcast(fn_ptn,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
!   call mpi_bcast(fn_wav,  maxchar, mpi_character, 0, mpi_comm_world, ierr)

!   call init_mpi_shift_reduce()
! #endif

  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)
  close(lunint)

  
  open(lunptn, file=fn_ptn)
  call init_j_partition(ptn, lunptn)
  close(lunptn)

  call set_n_ferm(ptn%n_ferm(1), ptn%n_ferm(2), 0)

  call set_dim_jvjj()
  call set_order_coupling()
  
  
  ! -------------------------------
  if (myrank==0) then
     write(*,*)
     write(*,*) 'snt file : ', trim(fn_int)
     write(*,*) 'ptn file : ', trim(fn_ptn)
     write(*,*) 'jwv file : ', trim(fn_wav)
     write(*,'( /, "Z= ", i3, "   N= ", i3 )')  &
          ptn%n_ferm(1)+n_core(1), ptn%n_ferm(2)+n_core(2)
     write(*,'(a,i3,/)') ' parity= ', ptn%iprty
     write(*,*)
  end if

  open(lunwv, file=fn_wav, form='unformatted', &
       status='old', access='stream')
  read(lunwv) n_eig, mtot
  allocate( e_val(n_eig), jj_val(n_eig) )
  read(lunwv) e_val
  read(lunwv) jj_val

  if (myrank==0) then
     write(*,'(a,i3,a,i3,a)') 'n_eig= ', n_eig, '   M= ', mtot, '/2'
     write(*,*) '  i     JP       E(MeV)'
     do istate = 1, n_eig
        write(*,'(i4,i5,3a,f12.5)') istate, jj_val(istate),'/2', &
             prty2char(ptn%iprty), ' ', e_val(istate)
     end do
     write(*,*)
  end if


  call main()


  close(lunwv)
  
  ! -------------------------------

! #ifdef MPI
!   call mpi_finalize(ierr)
! #endif


contains

  subroutine main()
    integer :: istate, n, id, i, j, maxi, maxj, mini, minj
    real(8) :: rto_v(n_ferm_pn+1), r, t(n_ferm_pn+1)
    integer, parameter :: maxJp=100
    real(8) :: rto_JpJn(maxJp, maxJp), tJpJn(maxJp, maxJp)

    do istate = 1, n_eig
       if (myrank==0) then
          write(*,*) '---------------------------------'
          write(*,*) '  i     JP       E(MeV)'
          write(*,'(i4,i5,3a,f12.5)') istate, jj_val(istate),'/2', &
               prty2char(ptn%iprty), ' ', e_val(istate)
          write(*,*)
       end if

       if (istate == 1) then
          call jscheme_dim_ptn(ptn, jj_val(istate))
       else
          if (jj_val(istate-1) /= jj_val(istate)) &
               call jscheme_dim_ptn(ptn, jj_val(istate))
       end if

       allocate( vec( ptn%jdim_acc( ptn%n_pidpn ) ) )
       read(lunwv) vec(:)

       r = dot_product(vec, vec)
       if (abs(r-1.d0)>1.d-5) write(*,*) '*** WARNING *** norm', r
       

       call init_store_top_j(storage_top_j, 10)

       rto_v(:) = 0.d0
       rto_JpJn(:,:) = 0d0
       n = ptn%n_pidpn
       !$omp parallel do private(t, tJpJn) reduction(+: rto_v, rto_JpJn) schedule(dynamic)
       do id = 1, n
          call check_vec_ptn(ptn, jj_val(istate), id, vec, &
               storage_top_j, t, tJpJn)
          rto_v = rto_v + t
          rto_JpJn = rto_JpJn + tJpJn
       end do

       deallocate(vec)
       call print_top_j(storage_top_j)

       if (myrank==0) then
          write(*,*) "*** seniority ratio *** "
          write(*,*) "  v :  ratio "
          do i = 0, n_ferm_pn
             if (rto_v(i+1) == 0.d0) cycle
             write(*,'(i3,a,f10.7)') i, ' : ',rto_v(i+1)
!             write(*,*) i, ' : ',rto_v(i+1)
          end do
          write(*,'(a,f10.7)') " sum: ", sum(rto_v)

          maxi = 0
          mini = maxJp
          do i = 1, maxJp
             if (sum(rto_JpJn(i,:))==0d0) cycle
             if ( i > maxi ) maxi = i
             if ( i < mini ) mini = i
          end do

          maxj = 0
          minj = maxJp
          do j = 1, maxJp
             if (sum(rto_JpJn(:,j))==0d0) cycle
             if ( j > maxj ) maxj = j
             if ( j < minj ) minj = j
          end do
          
          write(*,*)
          write(*,*) "Jp Jn: ratio"
          write(*,'(a)') '2*Jp    2*Jn->'
          write(*, '(100i8)', advance='no') (/( j-1, j=minj, maxj, 2 )/)
          write(*,*) '     sum'
          do i = mini, maxi, 2
             write(*,'(i3,100f8.4)') i-1, rto_JpJn(i, minj:maxj:2), sum(rto_JpJn(i, minj:maxj:2))
          end do
          write(*,'(a,100f8.4)') 'sum', sum(rto_JpJn(:, minj:maxj:2), dim=1)
          write(*,*)
          
       end if
       

    end do

  end subroutine main

  

  function prty2char(i) result (r)
    integer, intent(in) :: i
    character(1) :: r
    select case (i)
    case(-1)
       r = '-'
    case(1)
       r = '+'
    case default
       stop 'error prty2char'
    end select
  end function prty2char





end program print_jwav

  
