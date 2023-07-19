!
! calculate B(E{rank}),B(M{rank-1}) of two wave functions 
!
! ./transit_em.exe foobar.input
!  Thanks to Yusuke Tsunoda (CNS Tokyo)
!

program transit_em
#ifdef MPI
  use mpi
#endif
  use constant, only: kwf, kdim, kmbit, maxchar, c_no_init, pi
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_jorb_pn, n_jorb, n_ferm, n_core, is_debug, &
       jorbn, jorb, lorb, korb, norb, itorb, iporb
  use model_space, only: m_mass=>mass, print_max_l_vec, nv_shift
  use class_stopwatch
  use interaction, only: read_interaction, hamltn, ham_cm, j_square, &
       set_em, metensor, mmstensor, mmltensor
  use operator_jscheme, only: opr_j, set_opr_j
  use operator_mscheme, only: opr_m, add_opr_m, opr_m_p, opr_m_one_crt
  use partition, only: type_ptn_pn, init_partition, deploy_partition, &
       cost_from_localdim
  use wavefunction, only: type_vec_p, load_wf, dot_product_global
  use bridge_partitions, only: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, &
       init_bp_operator, bp_operate, finalize_bp_operator, &
       ex_val, init_mpi_shift_reduce
  use bp_expc_val, only: bp_ex_vals_pn, bp_ex_vals_ij
  use bp_io, only: bp_load_wf
  use rotation_group, only: dcg
  !$ use omp_lib, only : omp_get_max_threads
  implicit none
  type(type_ptn_pn), target :: ptnl, ptnr ! partition information
  type(type_bridge_partitions) :: bp
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn_int, fn_nml, &
       fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r, ctmp
  integer :: hw_type, mass
  real(8) :: eff_charge(2), gl(2), gs(2), e1_charge(2)
  !
  type(type_vec_p), allocatable :: evec_l(:), evec_r(:)
  real(8), allocatable :: cost(:), evs(:,:,:,:), &
       evv(:,:), evm(:,:), evs_p_ij(:,:,:,:,:), evs_n_ij(:,:,:,:,:)
  real(8) :: x, y
  integer :: i, j, n, mtotl, mtotr, n_eig_l, n_eig_r, nop, minm, maxm, jl, jr
  type(opr_m_p), allocatable :: ops(:)
  logical :: is_print_reduced_me = .false.
  ! logical :: is_print_reduced_me = .true.
  !
  integer :: rank
  namelist /input/ fn_int, hw_type,  &
       fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r, &
       eff_charge, gl, gs, e1_charge, mass, rank, nv_shift


#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
!  write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank/=0) is_debug = .false.
#endif
  !$ if(myrank==0) write(*,'(1a,1i3)') "OpenMP  # of threads=", omp_get_max_threads()

! default parameters
  mass = 0             ! mass number, optional
  hw_type = 1          ! harmonic oscillator formula
  fn_ptn_l = c_no_init        ! file name of left partition
  fn_load_wave_l = c_no_init  ! file name of left wave function
  fn_ptn_r = c_no_init        ! file name of right partition
  fn_load_wave_r = c_no_init  ! file name of right wave function
  eff_charge = (/ 1.d0,  0.d0 /) ! effective charges for E2 operator
  e1_charge = (/ 0.d0, 0.d0 /)  ! effective charges for E1 operator
  gl(:) = (/1.d0, 0.d0/) ! gyromagnetic ratios for orbital angular momentum
  gs(:) = (/5.586d0, -3.826d0/) ! gyromagnetic ratios for spin
  rank = 2 ! rank of electromagnetic transitions E{rank}, M{rank-1}
  nv_shift = 0           ! # of vectors for shift 
!
  call getarg(1, fn_nml)
#ifdef MPI
  call mpi_bcast(fn_nml,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
#endif
  open(lunnml, file=fn_nml, status='old')
  read(lunnml, nml=input)  
  close(lunnml)


#ifdef MPI
  call init_mpi_shift_reduce()
#endif

  if (nv_shift == 0) nv_shift = 1
  

  if (myrank==0) write(*,nml=input)
  if (myrank==0) write(*,'(a,3i3)') "compile conf. kwf, kdim, kmbit =",kwf,kdim,kmbit


  ! read header of wave functions
  open(lunwv, file=fn_load_wave_l, form='unformatted', status='old', access='stream')
  read(lunwv) n_eig_l
  read(lunwv) mtotl
  close(lunwv)

  open(lunwv, file=fn_load_wave_r, form='unformatted', status='old', access='stream')
  read(lunwv) n_eig_r
  read(lunwv) mtotr
  close(lunwv)

  
  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)

  open(lunptn, file=fn_ptn_l)
  if (myrank==0) write(*,'("set left partition_file=",1a)') trim(fn_ptn_l)
  call init_partition(ptnl, lunptn, mtotl)
  close(lunptn)

  open(lunptn, file=fn_ptn_r)
  if (myrank==0) write(*,'("set right partition_file=",1a)') trim(fn_ptn_r)
  call init_partition(ptnr, lunptn, mtotr)
  close(lunptn)

  call set_n_ferm(ptnl%n_ferm(1), ptnl%n_ferm(2), mass)
  call read_interaction(lunint, hw_type=hw_type)
  close(lunint)

  allocate( evec_l(n_eig_l), evec_r(n_eig_r))

  allocate( cost(ptnl%n_pidpnM) )
  call cost_from_localdim(ptnl, ptnl%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptnl, cost)
  deallocate(cost)
!  call deploy_partition(ptnl)

  allocate( cost(ptnr%n_pidpnM) )
  call cost_from_localdim(ptnr, ptnr%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptnr, cost)
  deallocate(cost)
!  call deploy_partition(ptnr)


  if (myrank==0) then
     x = ptnl%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for left global Lanczos vector:", x, " GB"
     x = ptnl%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory / process is:", x, " GB "
     y = x * n_eig_l
     write(*,*)
     x = ptnr%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for right global Lanczos vector:", x, " GB"
     x = ptnr%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory / process is:", x, " GB "
     write(*,*)
     y = y + x * n_eig_r
     write(*,'(1a,1f10.3,1a)') "Total Memory / process is:", y, " GB "
     write(*,*)
  end if

  call start_stopwatch(time_total, is_reset=.true.)

  call bp_load_wf(fn_load_wave_l, evec_l, ptnl, fn_ptn_l, mtotl)
  call bp_load_wf(fn_load_wave_r, evec_r, ptnr, fn_ptn_r, mtotr)
!  call load_wf(evec_l, ptnl, fn_load_wave_l, is_sorted=.true.)
!  call load_wf(evec_r, ptnr, fn_load_wave_r, is_sorted=.true.)

  call init_bridge_partitions(bp, ptnl, ptnr)
  

  do i = 1, n_eig_l
     x = dot_product_global(evec_l(i), evec_l(i))
     if (abs(x-1.d0) > 1.d-4) then
        if (myrank==0) write(*,*) "Warning ... normalization left",i,x
        evec_l(i)%p = 1.d0/sqrt(x) * evec_l(i)%p
     end if
  end do
  do i = 1, n_eig_r
     x = dot_product_global(evec_r(i), evec_r(i))
     if (abs(x-1.d0) > 1.d-4) then
        if (myrank==0) write(*,*) "Warning ... normalization right",i,x
        evec_r(i)%p = 1.d0/sqrt(x) * evec_r(i)%p
     end if
  end do

  minm = abs( abs(mtotl) - abs(mtotr)) / 2
  maxm = ( abs(mtotl) + abs(mtotr) ) / 2
  allocate( evv(n_eig_l, n_eig_r), evm(n_eig_l, n_eig_r) )

  call set_em(rank)

  if (minm <= rank .and. ptnl%iprty == ptnr%iprty*(-1)**rank &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_el()
     
  if (minm <= rank-1 .and. ptnl%iprty == ptnr%iprty*(-1)**rank &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_ml()
  
  

  call finalize_bridge_partitions(bp)
  call stop_stopwatch(time_total)
  if (myrank==0) print "(A, F10.3,/)", &
         "total time it took was:", time_total%time

  call print_max_l_vec()

#ifdef MPI
  call mpi_finalize(ierr)
#endif
  

contains


  subroutine calc_el()
    ! El transition  
    integer :: i, j, jl, jr
    real(8) :: x
    call init_bp_operator(bp, metensor)
    nop = 1
    allocate( ops(nop), evs(2, nop, n_eig_l, n_eig_r))
    evs = 0.d0
    ops(nop)%p => metensor
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > rank*2) cycle
          x = dcg(jr, mtotr, rank*2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call bp_ex_vals_pn(bp, evec_l(i), ops, evs(:,:,i,j), evec_r(j))
          evs(:,:,i,j) = evs(:,:,i,j) * sqrt(dble(jl+1)) / x
       end do
    end do
    evv(:,:) = eff_charge(1)*evs(1,1,:,:) + eff_charge(2)*evs(2,1,:,:)
    evm = 0.d0
    if (myrank==0) then
       write(*,*)
       write(*,'(a,i1,a,i2,a,2f8.4,a,2i3)') " E", rank, " transition  e^2 fm^", rank*2, "  eff_charge=", &
            eff_charge, " parity",ptnl%iprty, ptnr%iprty
       call print_trans(evv, evm, rank*2)
    end if

    call finalize_bp_operator(bp, metensor)
  end subroutine calc_el


  
  subroutine calc_ml()
    ! Ml transition
    integer :: i, j, jl, jr
    real(8) :: x
    call init_bp_operator(bp, mmltensor)
    call init_bp_operator(bp, mmstensor)
    nop = 2
    if (allocated(ops)) deallocate(ops, evs)
    if (allocated(evs_p_ij)) deallocate(evs_p_ij, evs_n_ij)
    allocate( ops(nop), evs(2, nop, n_eig_l, n_eig_r), &
         evs_p_ij(n_jorb(1), n_jorb(1), nop, n_eig_l, n_eig_r ), &
         evs_n_ij(n_jorb(2), n_jorb(2), nop, n_eig_l, n_eig_r ) )

    evs_p_ij = 0.d0
    evs_n_ij = 0.d0
    evs = 0.d0
    ops(1)%p => mmltensor
    ops(2)%p => mmstensor
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl-jr) > rank*2-2) cycle
          x = dcg(jr, mtotr, rank*2-2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call bp_ex_vals_ij(bp, evec_l(i), ops, &
               evs_p_ij(:,:,:,i,j), evs_n_ij(:,:,:,i,j), evec_r(j))
          evs_p_ij(:,:,:,i,j) = evs_p_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          evs_n_ij(:,:,:,i,j) = evs_n_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          do n = 1, nop
             evs(1,n,i,j) = sum(evs_p_ij(:,:,n,i,j))
             evs(2,n,i,j) = sum(evs_n_ij(:,:,n,i,j))
          end do
       end do
    end do

    evv(:,:) = gl(1)*evs(1,1,:,:) + gl(2)*evs(2,1,:,:) &
         + gs(1)*evs(1,2,:,:) + gs(2)*evs(2,2,:,:)
    evm = 0.d0
    if (myrank==0) then
       write(*,*)
       write(*,'(a,i1,a,i2,a,4f8.4,a,2i3)') " M", rank-1, &
            " transition  mu_N^2 fm^", rank*2-4, "  gl,gs=", gl, gs, &
            " parity",ptnl%iprty, ptnr%iprty
       call print_trans(evv, evm, rank*2-2)
    end if


    call finalize_bp_operator(bp, mmltensor)
    call finalize_bp_operator(bp, mmstensor)
  end subroutine calc_ml


  subroutine print_trans(evv, evm, mple)
    real(8), intent(in) :: evv(:,:), evm(:,:)
    integer, intent(in) :: mple
    integer :: i, j, jl, jr
    real(8) :: x, y
    character(len=maxchar) :: fmt

    write(*,*) '2xJi      Ei      2xJf     Ef       Ex       Mred.      B(EM )->     B(EM)<-     Mom.'
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0) cycle
          if (abs(jl-jr) > mple) cycle
          if (abs(mtotl-mtotr) > mple) cycle
          x = dcg(jr, mtotr, mple, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle

          x = evv(i, j)
          y = evm(i, j)
          if (mple <= 2) then 
             fmt = '(2(i2,"(",i4,")",f9.3), f8.3,5f12.5)'
          else
             fmt = '(2(i2,"(",i4,")",f9.3), f8.3,5f12.2)'
          end if
          write(*, fmt) &
               evec_l(i)%jj, i, evec_l(i)%eval, &
               evec_r(j)%jj, j, evec_r(j)%eval, &
               evec_r(j)%eval - evec_l(i)%eval, &
               x, x**2/dble(evec_l(i)%jj+1), x**2/dble(evec_r(j)%jj+1), y
       end do
    end do
    write(*,*)
  end subroutine print_trans



end program transit_em



