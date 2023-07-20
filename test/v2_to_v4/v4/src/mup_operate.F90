!
!   ./mup_operate.exe foo.snt bar.ptn input.wav output.wav
!
!       |M+1> = J^+ |M>  by  Ladder operator
!
!

program mup_oeprate
#ifdef MPI
  use mpi
#endif
  use constant, only: kwf, kdim, kmbit, maxchar, c_no_init
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_ferm
  use model_space, only: m_mass=>mass, print_max_l_vec, &
       nprocs_reduce, nprocs_shift, allocate_l_vec, deallocate_l_vec, nv_shift, is_mpi
  use class_stopwatch
  use interaction, only: read_interaction, jtensor
  use partition, only: type_ptn_pn, init_partition, deploy_partition, &
       cost_from_localdim, finalize_partition
  use wavefunction, only: type_vec_p, load_wf, dot_product_global
  use bridge_partitions, only: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, &
       init_bp_operator, bp_operate, finalize_bp_operator, &
       init_mpi_shift_reduce
  use bp_io, only: bp_save_wf, bp_load_wf
  !$ use omp_lib, only : omp_get_max_threads
  implicit none
  integer, parameter :: maxj=50
  type(type_ptn_pn), target :: ptnlr(-maxj:maxj)
  type(type_ptn_pn), pointer :: ptnl, ptnr
  type(type_bridge_partitions) :: bp
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn
  !
  type(type_vec_p), allocatable :: evec_l(:), evec_r(:)
  real(8), allocatable :: cost(:)
  real(8) :: x
  integer :: i, n, mtotl, mtotr, n_eig_l, n_eig_r
  !
  character(len=maxchar) :: fn_int, fn_ptn, fn_load_wave, fn_save_wave

#ifdef MPI
  is_mpi = .true.
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') &
       "nprocs", nprocs, "    myrank", myrank
#endif
  !$ if(myrank==0) write(*,'(1a,1i3)') "OpenMP  # of threads=", omp_get_max_threads()

  call start_stopwatch(time_total, is_reset=.true.)

  ! default parameters
  fn_ptn = c_no_init        ! file name of partition
  fn_load_wave = c_no_init  ! file name for load
  fn_save_wave = c_no_init  ! file name for save, M+1
  mtotl = -1                 ! 2*Jz for save w.f.
  nv_shift = 0
  !

  if (myrank == 0) then
     n = 0 
     do i = 1, iargc()
        call getarg(i, fn)
        if (fn(:1) == "-") cycle
        n = n + 1
        if (n==1) fn_int = fn 
        if (n==2) fn_ptn = fn 
        if (n==3) fn_load_wave = fn 
        if (n==4) fn_save_wave = fn 
     end do
  end if
  
  if (myrank == 0 .and. n < 4) then 
     write(*,'(/,a,/)') &
          'usage: mup_operate.exe foo.snt bar.ptn input.wav output.wav'
     stop 
  end if

#ifdef MPI
  call mpi_bcast(fn_int,        maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_ptn,        maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_load_wave,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
  call mpi_bcast(fn_save_wave,  maxchar, mpi_character, 0, mpi_comm_world, ierr)

  call init_mpi_shift_reduce()
#endif
  if (nv_shift == 0) nv_shift = 1


  if (myrank==0) write(*,'(a,3i3)') "compile conf. kwf, kdim, kmbit =", &
       kwf, kdim, kmbit

  ! read header of wave functions
  open(lunwv, file=fn_load_wave, form='unformatted', &
       status='old', access='stream')
  read(lunwv) n_eig_r
  read(lunwv) mtotr
  close(lunwv)

  allocate(evec_r(n_eig_r))

  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)

  if (myrank==0) write(*,'("set partition_file=",1a)') trim(fn_ptn)

  if (mtotr > maxj) stop 'ERROR: increase maxj'
  ptnr => ptnlr(mtotr)

  open(lunptn, file=fn_ptn)
  call init_partition(ptnr, lunptn, mtotr, verbose=.false.)
  close(lunptn)

  allocate( cost(ptnr%n_pidpnM) )
  call cost_from_localdim(ptnr, ptnr%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptnr, cost, verbose=.false.)
  deallocate(cost)


  call set_n_ferm(ptnr%n_ferm(1), ptnr%n_ferm(2), 0)
  call read_interaction(lunint)
  close(lunint)


  call bp_load_wf(fn_load_wave, evec_r, ptnr, fn_ptn, mtotr)

  call m_up_1()

  ! call print_max_l_vec()
  call finalize()

contains

  subroutine m_up_1()
    integer :: jj, i
    ! |M+1> = J+ |M>
    mtotl = mtotr + 2
    if (mtotr > maxj) stop 'ERROR: increase maxj'
    ptnl => ptnlr(mtotl)

    open(lunptn, file=fn_ptn)
    call init_partition(ptnl, lunptn, mtotl, verbose=.false.)
    close(lunptn)

    allocate( cost(ptnl%n_pidpnM) )
    call cost_from_localdim(ptnl, ptnl%pidpnM_pid_srt, nprocs, cost)
    call deploy_partition(ptnl, cost)
    deallocate(cost)


    call init_bridge_partitions(bp, ptnl, ptnr)
    call init_bp_operator(bp, jtensor, verbose=.false.)

    if (myrank==0) write(*,'(/5x,a35,a,a40/)') fn_load_wave, &
         '=> ',fn_save_wave

    allocate( evec_l(n_eig_r) )
    n_eig_l = 0
    do i = 1, n_eig_r
       jj = evec_r(i)%jj

       if (jj < mtotl) then
          if (myrank==0) write(*,'(a,i3,a,i3,a,i3,a,f8.3)') &
               'skip |J=',jj, '/2 M=',mtotr, '/2> ',i,'-th ',&
               evec_r(i)%eval
          cycle
       end if

       n_eig_l = n_eig_l + 1
       if (myrank==0) then 
           write(*,'(a,i3,a,i3,a,i3,a,f8.3,a,i3,a,i3,a,i3,a)') &
            'move |J=', jj, '/2,M=',mtotr, '/2> ', &
            i,'-th ', evec_r(i)%eval, &
            '  => |J=', jj, '/2,M=',mtotl, '/2> ', &
            n_eig_l,'-th w.f.'
        end if
       call allocate_l_vec( evec_l(n_eig_l)%p, ptnl%max_local_dim )
       evec_l(n_eig_l)%jj   = evec_r(i)%jj
       evec_l(n_eig_l)%eval = evec_r(i)%eval
       call bp_operate(bp, evec_l(n_eig_l), jtensor, evec_r(i))

    end do

    if (n_eig_l==0) then 
       if (myrank==0) write(*,*) "No w.f. to be saved"
       call finalize()
    end if

    if (myrank==0) write(*,*)

    do i = 1, n_eig_l
       x = dot_product_global(evec_l(i), evec_l(i))
       !     if (myrank==0) write(*,*) "normalization ",i,x
       evec_l(i)%p = 1.d0/sqrt(x) * evec_l(i)%p
    end do

    call finalize_bp_operator(bp, jtensor)
    call finalize_bridge_partitions(bp)

    call bp_save_wf(fn_save_wave, evec_l(:n_eig_l), ptnl)
    ! , fn_ptn, mtotl)
    
    do i = 1, n_eig_l
       if (associated(evec_l(i)%p)) call deallocate_l_vec(evec_l(i)%p)
    end do
    do i = 1, n_eig_r
       if (associated(evec_r(i)%p)) call deallocate_l_vec(evec_r(i)%p)
    end do

  end subroutine m_up_1

  
  subroutine finalize()
    call stop_stopwatch(time_total)
    if (myrank==0) print "(A, F10.3,/)", &
         "total time it took was:", time_total%time
#ifdef MPI
    call mpi_finalize(ierr)
#endif
    stop
  end subroutine finalize
  


end program mup_oeprate



