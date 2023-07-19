!
! print M-scheme bit and wave functions
!   usage:  ./print_mbit_wf.exe foo.snt bar.ptn foobar.wav
!

module mod_print_mbit_wf
  use constant
  use model_space
  use partition
  use wavefunction, only : type_vec_p
  implicit none

contains

  subroutine print_mbitwf(ptn, evec)
    type(type_ptn_pn), intent(in) :: ptn
    type(type_vec_p), intent(in) :: evec(:)
    integer :: id, idp, idn, mmp, mmn, ip, in, idimp, idimn, neig
    integer(kdim) :: ist, ien, ndim
    integer(kmbit) :: mbp, mbn
    integer :: ipn, i, j, k, l, m, n
    real(kwf) :: v(size(evec))


    write(*,*) "order of single-particle states"
    write(*,*) " idx:  orbit          Jz"
    n = 0
    do k = 1, n_jorb_pn
       do m = -jorb(k), jorb(k), 2
          n = n + 1
          write(*,'(i5,3a,i5,a)') n,': ', corb(k),'  Jz=', m,'/2'
       end do
       if (k==n_jorb(1).or.k==n_jorb_pn) write(*,*) "---------------------"
    end do
    write(*,*) 

    write(*,*) "    dim:  proton occ. neutron occ.  v(dim,1), v(dim,2), ..."

    neig = size(evec)

    do id = ptn%idl_start, ptn%idl_end
       idp = ptn%pidpnM_pid(1,id)
       idn = ptn%pidpnM_pid(2,id)
       mmp = ptn%pidpnM_pid(3,id)
       mmn = ptn%mtotal - mmp
       ist = ptn%local_dim_acc_start(id)
       ien = ptn%local_dim_acc(id)
       idimp = ptn%pn(1)%id(idp)%mz(mmp)%n
       idimn = ptn%pn(2)%id(idn)%mz(mmn)%n

       ndim = ist 
       do in = 1, idimn
          mbn = ptn%pn(2)%id(idn)%mz(mmn)%mbit(in)
          do ip = 1, idimp
             mbp = ptn%pn(1)%id(idp)%mz(mmp)%mbit(ip)
             v(:) = (/( evec(i)%p(ndim), i=1, neig )/)
             call print_mbitv(ndim, mbp, mbn, v)
             ndim = ndim + 1
          end do
       end do
    end do

  end subroutine print_mbitwf

  subroutine print_mbitv(ndim, mbp, mbn, v)
    integer(kdim) :: ndim
    integer(kmbit), intent(in) :: mbp, mbn
    real(kwf), intent(in) :: v(:)
    integer :: i
    write(*,'(i8,a)', advance='no') ndim, ": "
    do i = 1, n_morb(1)
       if ( btest(mbp,i) ) then 
          write(*,'(a)', advance='no') "x"
       else
          write(*,'(a)', advance='no') "-"
       end if
    end do
    write(*,'(a)', advance='no') " "
    do i = 1, n_morb(2)
       if ( btest(mbn,i) ) then 
          write(*,'(a)', advance='no') "x"
       else
          write(*,'(a)', advance='no') "-"
       end if
    end do
    write(*,'(1000f12.8)') v(:)

  end subroutine print_mbitv

end module mod_print_mbit_wf


program print_mbit_wf
  use constant
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_ferm, mass
  use partition
  use wavefunction, only: type_vec_p, wf_alloc_vec
  use bp_io, only: bp_load_wf
  use mod_print_mbit_wf
  implicit none
  integer, parameter :: lunint=11, lunptn=12, lunwv=13
  character(maxchar) :: fn_int, fn_ptn, fn_load_wave
  integer :: mtot, n_eigen
  type(type_ptn_pn) :: ptn
  type(type_vec_p), allocatable :: evec(:)
  integer :: iargc 

  if (myrank==0) then 
     if (iargc() < 3) stop "usage:  ./print_mbit_wf.exe foo.snt bar.ptn foobar.wav "
     call getarg(1, fn_int)
     call getarg(2, fn_ptn)
     call getarg(3, fn_load_wave)
  end if

  if (nv_shift == 0) nv_shift = 1
  if (nprocs /= 1) stop "NOT compatible with MPI"

  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)
  close(lunint)


  ! read header of wave functions
  open(lunwv, file=fn_load_wave, form='unformatted', &
       status='old', access='stream')
  read(lunwv) n_eigen
  read(lunwv) mtot
  close(lunwv)
  

  open(lunptn, file=fn_ptn, status='old')
  if (myrank==0) write(*,'(1a,1i3,2a)') "set partition Mtotal=", mtot, &
       "  partition_file= ", trim(fn_ptn)
  call init_partition(ptn, lunptn, mtot)
  close(lunptn)

  if (myrank==0) write(*,'(3/,a,i5,a)') "*** total M subspace ***, M= ", mtot,"/2"
  if (myrank==0) write(*,'(a,i5,3/)') "# of eigenvectors : ", n_eigen

  call set_n_ferm(ptn%n_ferm(1), ptn%n_ferm(2), mass)

  call deploy_partition(ptn, verbose=.true.)

  allocate( evec( n_eigen ) )

  call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot)

  call print_mbitwf(ptn, evec)

    
end program print_mbit_wf


