!
! calculate transition of two wave functions 
!
! ./transit.exe foobar.input  
!

module mod_transit
#ifdef MPI
  use mpi
#endif
  use constant, only: kwf, kdim, kmbit, nmbit, maxchar, c_no_init, pi
  use class_stopwatch
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_jorb_pn, n_jorb, n_ferm, n_core, is_debug, &
       jorb, lorb, korb, norb, itorb, iporb
  use model_space, only: m_mass=>mass, print_max_l_vec, &
       nprocs_reduce, nprocs_shift, nv_shift, is_mpi, allocate_l_vec, deallocate_l_vec
  use interaction, only: read_interaction, j_square, jtensor, &
       r2y2, r1y1, r3y3, ltensor, stensor, set_gt, gt_m, set_fm, fm_m, &
       set_ff, ff_0rs, ff_0sp, ff_1r, ff_1rs, ff_1p, ff_2rs, &
       ff_0rs_d, ff_1r_d, ff_1rs_d, r3y1
  use operator_jscheme, only: opr_j, set_opr_j
  use operator_mscheme, only: opr_m, opr_m_p, opr_m_one_crt, &
       opr_m_eff_charge, add_opr_m, print_op_1b
  use partition, only: type_ptn_pn, init_partition, copy_partition, &
       deploy_partition, cost_from_localdim, finalize_partition
  use wavefunction, only: type_vec_p, load_wf, dot_product_global, wf_alloc_vec, &
       hw_proj
  use bridge_partitions, only: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, &
       init_bp_operator, bp_operate, finalize_bp_operator, &
       ex_val, init_mpi_shift_reduce
  use bp_expc_val, only: bp_ex_vals_pn, bp_ex_vals_ij
  use bp_io, only: bp_load_wf
  use rotation_group, only: dcg
  !$ use omp_lib, only : omp_get_max_threads
  implicit none
  type(type_ptn_pn), target :: ptnl, ptnr_target ! partition information
  type(type_ptn_pn), pointer :: ptnr
  type(type_bridge_partitions) :: bp
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, &
       lunwv=13, lun_rme=14
  character(len=maxchar) :: fn_int, fn_nml, &
       fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r
  integer :: hw_type, mass
  real(8) :: eff_charge(2), gl(2), gs(2), e1_charge(2), e3_charge(2)
  logical, allocatable :: is_calc_lr(:,:)
  !
  type(type_vec_p), allocatable :: evec_l(:), evec_r(:)
  real(8), allocatable :: cost(:), evs(:,:,:,:), &
       evv(:,:), evm(:,:), evs_p_ij(:,:,:,:,:), evs_n_ij(:,:,:,:,:)
  integer ::  mtotl, mtotr, n_eig_l, n_eig_r
  type(opr_m_p), allocatable :: ops(:)
  logical :: is_fermi_trn=.false.
  logical :: is_obtd = .false.  ! for one-body transition density
  logical :: is_tbtd = .false.  ! for two-body transition density
  integer :: ihw_proj
  integer :: n_eig_lr_pair(2,1000) ! pair to be computed for TBTD

  ! m-up
  real(8), parameter :: thd_mup = 1.d-8
  integer :: mupr = -1
  type(type_ptn_pn) :: ptnm
  type(type_vec_p), allocatable :: evec_m(:)


contains

  subroutine transit_main()
    integer :: i, j, k, n, minm, maxm
    real(8) :: x, xl, xr
    integer :: access

    namelist /input/ fn_int, hw_type,  &
         fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r, &
         eff_charge, gl, gs, e1_charge, e3_charge, mass, is_fermi_trn, &
         is_obtd, is_tbtd, nprocs_reduce, nv_shift, ihw_proj, &
         n_eig_lr_pair

#ifdef MPI
    is_mpi = .true.
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
    !  write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
    if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
    if (myrank/=0) is_debug = .false.
#endif

    call start_stopwatch(time_total, is_reset=.true.)
    
    !$ if(myrank==0) write(*,'(1a,1i3)') "OpenMP  # of threads=", omp_get_max_threads()

    ! default parameters
    mass = 0             ! mass number, optional
    hw_type = 1          ! harmonic oscillator formula
    fn_ptn_l = c_no_init        ! file name of left partition
    fn_load_wave_l = c_no_init  ! file name of left wave function
    fn_ptn_r = c_no_init        ! file name of right partition
    fn_load_wave_r = c_no_init  ! file name of right wave function
    eff_charge = (/ 1.d0, 0.d0 /)  ! effective charges for E2 operator
    e1_charge  = (/ 0.d0, 0.d0 /)  ! effective charges for E1 operator
    e3_charge  = (/ 0.d0, 0.d0 /)  ! effective charges for E3 operator
    gl(:) = (/1.d0, 0.d0/) ! gyromagnetic ratios for orbital angular momentum
    gs(:) = (/5.586d0, -3.826d0/) ! gyromagnetic ratios for spin
    nprocs_reduce = 1      ! # of processes for MPI reduction
    nv_shift = 0           ! # of vectors for shift 
    ihw_proj = -1          ! projection to hw excited subspace, temporal
    n_eig_lr_pair(:,:) = -1  ! pair states to be computed for TBTD
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

    if (myrank==0) then 
       write(*,nml=input)
       write(*,'(a,4i3)') "compile conf. kwf, kdim, kmbit, nmbit =", &
            kwf, kdim, kmbit, nmbit
       write(*,'(/,2a)') 'fn_load_wave_l = ', trim(fn_load_wave_l)
       write(*,'(  2a)') 'fn_load_wave_r = ', trim(fn_load_wave_r)
    end if
    if (access(fn_load_wave_l, 'r') /= 0) then
       if (myrank==0) write(*,*) "NOT found ", trim(fn_load_wave_l), " in transit.exe"
       stop "file not found"
    end if
    if (access(fn_load_wave_r, 'r') /= 0) then
       if (myrank==0) write(*,*) "NOT found ", trim(fn_load_wave_r), " in transit.exe"
       stop "file not found"
    end if
       


    ! read header of wave functions
    open(lunwv, file=fn_load_wave_l, form='unformatted', &
         status='old', access='stream')
    read(lunwv) n_eig_l
    read(lunwv) mtotl
    close(lunwv)

    open(lunwv, file=fn_load_wave_r, form='unformatted', &
         status='old', access='stream')
    read(lunwv) n_eig_r
    read(lunwv) mtotr
    close(lunwv)


    open(lunint, file=fn_int, status='old')
    call read_sps(lunint)

    open(lunptn, file=fn_ptn_l)
    if (myrank==0) write(*,'("set left partition_file=",1a)') trim(fn_ptn_l)
    call init_partition(ptnl, lunptn, mtotl)
    close(lunptn)

    allocate( cost(ptnl%n_pidpnM) )
    call cost_from_localdim(ptnl, ptnl%pidpnM_pid_srt, nprocs, cost)
    call deploy_partition(  ptnl, cost )
    deallocate( cost )

    if (fn_ptn_r == fn_ptn_l .and. mtotr == mtotl) then
       if (myrank==0) write(*,'(a,/)') 'set right partition => left partition'
       ptnr => ptnl
    else
       ptnr => ptnr_target

       open(lunptn, file=fn_ptn_r)
       if (myrank==0) write(*,'("set right partition_file=",1a)') trim(fn_ptn_r)
       call init_partition(ptnr, lunptn, mtotr)
       close(lunptn)

       allocate( cost(ptnr%n_pidpnM) )
       call cost_from_localdim(ptnr, ptnr%pidpnM_pid_srt, nprocs, cost)
       call deploy_partition(ptnr, cost)
       deallocate(cost)
    end if


    call set_n_ferm(ptnl%n_ferm(1), ptnl%n_ferm(2), mass)
    call read_interaction(lunint, hw_type=hw_type)
    close(lunint)

    if ( all(e1_charge == 0.d0) ) then
       e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(m_mass)
       e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(m_mass)
    end if

    if ( all(e3_charge == 0.d0) ) e3_charge = eff_charge

    allocate( evec_l(n_eig_l), evec_r(n_eig_r) )

    allocate( is_calc_lr(n_eig_l, n_eig_r) )
    if (n_eig_lr_pair(1,1) == -1) then 
       is_calc_lr(:,:) = .true.
       if (myrank==0) write(*,*) "all combination of states to be computed."
    else
       is_calc_lr(:,:) = .false.
       do n = 1, size(n_eig_lr_pair, 2)
          i = n_eig_lr_pair(1, n)
          j = n_eig_lr_pair(2, n)
          if (i<0 .or. j<0) exit
          if (i > n_eig_l .or. j > n_eig_r) cycle
          if (i==0 .and. j==0) then  
             do k = 1, min(n_eig_l, n_eig_r)
                is_calc_lr(k, k) = .true. ! n_eig_lr_pair = 0,0 : diagonal only
             end do
          else if (i==0) then
             is_calc_lr(:, j) = .true.
          else if (j==0) then 
             is_calc_lr(i, :) = .true.
          else
             is_calc_lr(i, j) = .true.
          end if
       end do
       if (myrank==0) then
          write(*,'(/,a)') ' -------------------------'
          do j = 1, n_eig_l
             do i = 1, n_eig_l
                if (.not. is_calc_lr(i,j)) cycle
                write(*,'(a,3i5)') "pair states to be computed : ", i, j
             end do
          end do
          write(*,'(a,/)') ' -------------------------'
       end if
    end if

    

    if (myrank==0) then
       xl = ptnl%ndim*kwf/1024.d0/1024.d0/1024.d0
       write(*,'(1a,1f10.3,1a)') &
            "Memory for left global Lanczos vector:", xl, " GB"
       xl = ptnl%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
       write(*,'(1a,1f10.3,1a)') "Memory / process is:", xl, " GB "
       write(*,*)
       xr = ptnr%ndim*kwf/1024.d0/1024.d0/1024.d0
       write(*,'(1a,1f10.3,1a)') &
            "Memory for right global Lanczos vector:", xr, " GB"
       xr = ptnr%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
       write(*,'(1a,1f10.3,1a)') "Memory / process is:", xr, " GB "
       write(*,*)
       i = n_eig_l 
       j = n_eig_r
#ifdef MPI
       i = i + nprocs_reduce - 1
       j = j + nprocs_shift - 1
       i = i + 1 ! vltmp
#endif
       x = xl * i + xr * j
       write(*,'(1a,1f10.3,1a)') "Total Memory / process is:", x, " GB "
       write(*,*)
    end if


    call bp_load_wf(fn_load_wave_l, evec_l, ptnl, fn_ptn_l, mtotl)
    call bp_load_wf(fn_load_wave_r, evec_r, ptnr, fn_ptn_r, mtotr)
    ! call load_wf(evec_l, ptnl, fn_load_wave_l, is_sorted=.true.)
    ! call load_wf(evec_r, ptnr, fn_load_wave_r, is_sorted=.true.)

    call init_bridge_partitions(bp, ptnl, ptnr, verbose=.true.)

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
       if (ihw_proj /= -1) call hw_proj(evec_r(i), ihw_proj)
    end do

    minm = abs( abs(mtotl) - abs(mtotr)) / 2
    maxm = ( abs(mtotl) + abs(mtotr) ) / 2
    allocate( evv(n_eig_l, n_eig_r), evm(n_eig_l, n_eig_r) )

    if (myrank==0) then
       write(*,'(a,5i4)') 'left  Z,N,A,M,prty: ', &
            ptnl%n_ferm(1)+n_core(1), &
            ptnl%n_ferm(2)+n_core(2), &
            sum(ptnl%n_ferm + n_core), &
            ptnl%mtotal, &
            ptnl%iprty          
       write(*,'(a,5i4)') 'right Z,N,A,M,prty: ', &
            ptnr%n_ferm(1)+n_core(1), &
            ptnr%n_ferm(2)+n_core(2), &
            sum(ptnr%n_ferm + n_core), &
            ptnr%mtotal, &
            ptnr%iprty          
    end if


    if (is_tbtd .and. mtotl == mtotr .and. ptnl%iprty == ptnr%iprty &
         ! if (.false. .and. mtotl == mtotr .and. ptnl%iprty == ptnr%iprty &
         !  if (mtotl == mtotr .and. ptnl%iprty == ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) ) then
       call calc_olp()
    end if

    if (is_tbtd) then
       call calc_tbtd()
       goto 999
    end if


    

    if ( .false. .and. &
         ! if ( .true. .and. &
         minm == 0 .and. mtotl==0 .and. ptnl%iprty == ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) ) then
       call calc_e0()
       goto 999
    end if


    if (minm <= 2 .and. ptnl%iprty == ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_e2()

    if (minm <= 1 .and. ptnl%iprty == ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_m1()

    if (minm <= 1 .and. ptnl%iprty /= ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_e1()

!    if (minm <= 1 .and. ptnl%iprty /= ptnr%iprty &
!         .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_isd()

    if (minm <= 3 .and. ptnl%iprty /= ptnr%iprty &
         .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_e3()


    if ( minm<= 1 .and. abs(ptnl%n_ferm(1)-ptnr%n_ferm(1))==1 &
         .and. ptnl%iprty == ptnr%iprty &
         .and. sum(ptnl%n_ferm) == sum(ptnr%n_ferm) ) call calc_gt()


    !  call calc_gttype_olp(3, 8)

    ! if ( .false. .and. &
    if ( .true. .and. &
         abs(ptnl%n_ferm(1)-ptnr%n_ferm(1))==1 &
         .and. sum(ptnl%n_ferm) == sum(ptnr%n_ferm) ) &
         call calc_gttype_obtd()


    if ( is_fermi_trn .and. &
         minm== 0 .and. abs(ptnl%n_ferm(1)-ptnr%n_ferm(1))==1 &
         .and. ptnl%iprty == ptnr%iprty &
         .and. sum(ptnl%n_ferm) == sum(ptnr%n_ferm) ) call calc_fm()

    if ( minm<= 2 .and. abs(ptnl%n_ferm(1)-ptnr%n_ferm(1))==1 &
         .and. ptnl%iprty /= ptnr%iprty &
         .and. sum(ptnl%n_ferm) == sum(ptnr%n_ferm) ) call calc_ff()


    if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)+1) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) )   call calc_sf(1)
    if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+1) ) call calc_sf(2)

    if ( (ptnl%n_ferm(1)+1 == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) .and. myrank==0)  &
         write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for s-factor"

    if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2)+1 == ptnr%n_ferm(2)) .and. myrank==0) &
         write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for s-factor"


    if (       (ptnl%n_ferm(1) == ptnr%n_ferm(1)+2) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) )   call calc_tna(1)
    if (       (ptnl%n_ferm(1) == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+2) ) call calc_tna(2)
    if (       (ptnl%n_ferm(1) == ptnr%n_ferm(1)+1) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+1) ) call calc_tna(3)

    if ( (ptnl%n_ferm(1)+2 == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) .and. myrank==0)  &
         write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for TNA"

    if ( (ptnl%n_ferm(1)+1 == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2)+1 == ptnr%n_ferm(2)) .and. myrank==0) &
         write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for TNA"

    ! if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)+2) &
    !      .and. (ptnl%n_ferm(2)+2 == ptnr%n_ferm(2)) &
    !      .and. mtotl == mtotr ) &
    !      call calc_nme_0v(1)
    ! if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)+2) &
    !     .and. (ptnl%n_ferm(2)+2 == ptnr%n_ferm(2)) ) &
    !     call calc_tbtd_ppnn()

    if ( (ptnl%n_ferm(1)+2 == ptnr%n_ferm(1)) &
         .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+2) &
         .and. mtotl == mtotr ) &
         write(*,'(/,a,/)') &
         "ERROR: swap fn_read_wave_l and _r for double beta decay"

999 continue  

    call finalize_bridge_partitions(bp)

    call finalize_mup()

    call stop_stopwatch(time_total)
    if (myrank==0) write(*, '(/, a, f10.3,/)') &
         "total elapsed time:", time_total%time

    call print_max_l_vec()

    call print_summary_stopwatch()
    
#ifdef MPI
    call mpi_finalize(ierr)
#endif

  end subroutine transit_main

  
  
  subroutine calc_olp()
    ! overlap just for check
    integer :: i, j
    real(8) :: x
    real(8), allocatable :: ev(:,:)

    allocate( ev(n_eig_l, n_eig_r) )
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          ev(i,j) = dot_product_global(evec_l(i), evec_r(j))
       end do
    end do
    if (myrank==0) then
       write(*,'(/,a)') " overlap "
       write(*,'(3x,1000i8)') (/( j, j=1, n_eig_r )/)
       do i = 1, n_eig_l 
          write(*,'(i3,1000f8.4)') i,ev(i,:)
       end do
       ev = ev ** 2
       write(*,'(/,a)') " overlap probability "
       write(*,'(3x,1000i8)', advance='no') (/( j, j=1, n_eig_r )/)
       write(*,'(a)') '      sum '
       do i = 1, n_eig_l 
          write(*,'(i3,1001f8.4)') i,ev(i,:), sum(ev(i,:))
       end do
       write(*,'(a3,1000f8.4)') 'sum', &
            (/( sum(ev(:,j)), j = 1, n_eig_r )/)
       write(*,*)
    end if

    do i = 1, min(n_eig_l, n_eig_r)
       ev(i,i) = ev(i,i) - 1.d0
    end do
    if (myrank==0) then
       write(*,'(a,e10.2)') "largest diff", maxval(abs(ev))
    end if

    deallocate( ev ) 

  end subroutine calc_olp



  subroutine calc_e2()
    ! E2 transition  
    use interaction, only : r2y2_func1
    integer :: i, j, jl, jr
    real(8) :: x
    character(len=maxchar) :: comment
    type(opr_m) :: op_e2
    
    if (is_obtd) call calc_obtd('E2', &
         r2y2_func1, 2, 1, 1, eff_charge)

    call opr_m_eff_charge(op_e2, r2y2, eff_charge)

    write(comment,'(a,2f8.4,a,2i3)') &
         " E2 transition  e^2 fm^4  eff_charge=", &
         eff_charge, " parity",ptnl%iprty, ptnr%iprty
    call calc_transit_op(op_e2, sqrt(16.d0*pi/5.d0), comment)


  end subroutine calc_e2



  subroutine calc_m1()
    ! M1 transition
    use interaction, only : ltensor_func1, stensor_func1
    integer :: i, j, jl, jr
    real(8) :: x
    logical :: is_mup(n_eig_l, n_eig_r)
    character(len=maxchar) :: comment
    type(opr_m) :: op_m1, op_s

    if (is_obtd) then
       call calc_obtd('L ', &
            ltensor_func1, 1, 1, 1, gl/sqrt(4.d0*pi/3.d0))
       call calc_obtd('S ', &
            stensor_func1, 1, 1, 1, gs/sqrt(4.d0*pi/3.d0))
    end if

    call opr_m_eff_charge(op_m1, ltensor, gl/sqrt(4.d0*pi/3.d0))
    call opr_m_eff_charge(op_s,  stensor, gs/sqrt(4.d0*pi/3.d0))
    ! call print_op_1b(op_s)
    call add_opr_m(op_m1, 1.d0, op_s)


    write(comment,'(a,4f8.4,a,2i3)') &
         " M1 transition  mu_N^2  gl,gs=", gl, gs, &
         " parity", ptnl%iprty, ptnr%iprty
    call calc_transit_op(op_m1, sqrt(4.d0*pi/3.d0), comment)

  end subroutine calc_m1



  subroutine calc_isd()
    !
    ! Isoscalar dipole
    ! Note : only r3Y1, not r3Y1 - <r2>r1Y1
    !
    use interaction, only : r3y1_func1
    integer :: i, j, jl, jr
    real(8) :: x
    character(len=maxchar) :: comment
    type(opr_m) :: op
    real(8) :: ec(2) = (/ 1.d0, 1.d0 /)

    if (is_obtd) call calc_obtd( 'ISD', &
         r3y1_func1, 1, 1, -1, ec )

    call opr_m_eff_charge(op, r3y1, ec)

    write(comment,'(a,2f8.4,a,2i3)') &
         " Isoscalar Dipole   fm^6  eff_charge=", &
         eff_charge, " parity", ptnl%iprty, ptnr%iprty
    call calc_transit_op(op, 0.d0, comment)


  end subroutine calc_isd



  subroutine calc_mup_ex_vals(is_mup, ops)
    !
    ! |ini J, M+1> = 1/N * J+ |ini J, M>
    !   <fin |op|ini J, M+1> at evv
    !
    logical, intent(in) :: is_mup(:,:)
    type(opr_m_p), intent(inout) :: ops(:)
    integer :: i, j, jl, jr, n, nop
    real(8) :: x

    call finalize_bridge_partitions(bp)

    call set_mup(is_mup)

    call init_bridge_partitions(bp, ptnl, ptnm)
    do i = 1, size(ops)
       call init_bp_operator(bp, ops(i)%p)
    end do

    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          if (.not. is_mup(i,j)) cycle
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > ops(1)%p%irank*2) cycle
          x = dcg(jr, mupr, ops(1)%p%irank*2, mtotl-mupr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          !
          call bp_ex_vals_ij(bp, evec_l(i), ops, &
               evs_p_ij(:,:,:,i,j), evs_n_ij(:,:,:,i,j), evec_m(j))
          evs_p_ij(:,:,:,i,j) = evs_p_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          evs_n_ij(:,:,:,i,j) = evs_n_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          do n = 1, size(ops)
             evs(1,n,i,j) = sum(evs_p_ij(:,:,n,i,j))
             evs(2,n,i,j) = sum(evs_n_ij(:,:,n,i,j))
          end do
          !
       end do
    end do

    do i = 1, size(ops)
       call finalize_bp_operator(bp, ops(i)%p)
    end do
    call finalize_bridge_partitions(bp)

    call init_bridge_partitions(bp, ptnl, ptnr)
    
  end subroutine calc_mup_ex_vals



  subroutine calc_e1()
    ! E1 transition  
    use interaction, only : r1y1_func1
    integer :: i, j, jl, jr
    real(8) :: x
    character(len=maxchar) :: comment
    type(opr_m) :: op_e1

    if (is_obtd) call calc_obtd('E1', &
         r1y1_func1, 1, 1, -1, e1_charge)

    call opr_m_eff_charge(op_e1, r1y1, e1_charge)
    write(comment, '(a,2f8.4,a,2i3)') &
         " E1 transition  e^2 fm^2  e1_charge=", &
         e1_charge, " parity", ptnl%iprty, ptnr%iprty
    call calc_transit_op(op_e1, 0.d0, comment)
  end subroutine calc_e1





  subroutine calc_e3()
    ! E3 transition  
    use interaction, only : r3y3_func1
    integer :: i, j, jl, jr
    real(8) :: x
    character(len=maxchar) :: comment
    type(opr_m) :: op_e3

    if (is_obtd) call calc_obtd('E3', &
         r3y3_func1, 3, 1, -1, e3_charge)

    call opr_m_eff_charge(op_e3, r3y3, e3_charge)
    write(comment, '(a,2f8.4,a,2i3)') &
         " E3 transition  e^2 fm^6  e3_charge=", &
         e3_charge, " parity", ptnl%iprty, ptnr%iprty
    call calc_transit_op(op_e3, 0.d0, comment)
  end subroutine calc_e3


  subroutine calc_gt()
    use interaction, only : gt_func1
    ! Gamow-Teller transition
    character(len=maxchar) :: comment
    integer :: irank=1, iprty=1, nbody=-10

    if (is_obtd) then
       nbody = -10
       if ( ptnl%n_ferm(1)-ptnr%n_ferm(1) == -1)  nbody = -11
       call calc_obtd('GT', gt_func1, 1, nbody, 1)
    end if

    call set_gt()
    write(comment,'(a,2i3)') &
         " GT transition  <f||sigma t^(+-)||i>    parity", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(gt_m, 0.d0, comment)
  end subroutine calc_gt



  subroutine calc_e0()
    use interaction, only : set_e0, r2y0, r2y0_func1
    ! E0 transition
    character(len=maxchar) :: comment
    real(8) :: e_charge(2)
    type(opr_m) :: op_e0

    ! e_charge = (/ 1.d0, 0.d0 /)
    e_charge = eff_charge
    if (myrank==0) write(*,'(a,2f10.5)') 'E0 charge ', e_charge

    if (is_obtd) then
       call calc_obtd('E0', r2y0_func1, 0, 1, 1, e_charge )
    end if

    call set_e0()

    call opr_m_eff_charge(op_e0, r2y0, e_charge)
    write(comment,'(a,2i3)') &
         " E0 transition  <f||E0||i>    parity", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(op_e0, 1.d0, comment)

  end subroutine calc_e0



  subroutine calc_gttype_obtd()
    ! Gamow-Teller-type (p^+ n)  OBTD for any rank
    use constant, only : max_int4
    use interaction, only : dummy_func1
    character(len=maxchar) :: comment
    integer :: nbody=-10
    character(4) :: op_type
    integer :: il, jjl, ir, jjr, minrl, maxrl, minrk, maxrk, irank, iprty


    nbody = -10
    if ( ptnl%n_ferm(1) - ptnr%n_ferm(1) == -1)  nbody = -11
    iprty = ptnl%iprty * ptnr%iprty

    minrl = max_int4
    maxrl = -1
    do il = 1, n_jorb(1)
       jjl = jorb(il)
       do ir = n_jorb(1)+1, n_jorb_pn
          jjr = jorb(ir)
          if ( iprty /= iporb(il)*iporb(ir) ) cycle
          minrl = min( abs(jjl - jjr)/2, minrl )
          maxrl = max(    (jjl + jjr)/2, maxrl )
       end do
    end do

    minrk = max_int4
    maxrk = -1
    do il = 1, n_eig_l
       jjl = evec_l(il)%jj
       do ir = 1, n_eig_r
          jjr = evec_r(ir)%jj
          minrk = min( abs(jjl - jjr)/2, minrk )
          maxrk = max(    (jjl + jjr)/2, maxrk )
       end do
    end do

    minrk = max(minrl, minrk)
    maxrk = min(maxrl, maxrk)

    do irank = minrk, maxrk
       write(op_type,'(i0)') irank
       op_type = 'R' // trim(op_type)
       call calc_obtd(op_type, dummy_func1, irank, nbody, iprty)
    end do

  end subroutine calc_gttype_obtd


  subroutine calc_gttype_olp(iorbl, iorbr)
    ! <Jl | * 1/N ( (il^+ ir)^Jc |Jr> )  
    ! N ... normalization factor
    ! TODO .... j-projection
    use interaction, only : set_ob_channel
    use interaction, only : dummy_func1
    integer, intent(in) :: iorbl, iorbr ! orbit number 
    character(len=maxchar) :: comment
    integer :: nbody=-10
    character(4) :: op_type
    integer :: il, jjl, ir, jjr, minrl, maxrl, minrk, maxrk, irank, iprty
    type(opr_m), allocatable :: ops(:)
    real(8) :: x(2), r, olp
    real(8), allocatable :: evv(:,:)
    integer, allocatable :: ij_orb(:,:)
    type(type_vec_p) :: vt

    nbody = -10
    if ( ptnl%n_ferm(1) - ptnr%n_ferm(1) == -1)  nbody = -11

    iprty = ptnl%iprty * ptnr%iprty
    if ( iprty /= iporb(iorbl)*iporb(iorbr) ) &
         stop "parity mismatch in calc_gttype_olp"


    jjl = jorb(iorbl)
    jjr = jorb(iorbr)
    minrl = abs(jjl - jjr)/2
    maxrl = (jjl + jjr)/2

    minrk = 1000000
    maxrk = -1
    do il = 1, n_eig_l
       jjl = evec_l(il)%jj
       do ir = 1, n_eig_r
          jjr = evec_r(ir)%jj
          minrk = min( abs(jjl - jjr)/2, minrk )
          maxrk = max(    (jjl + jjr)/2, maxrk )
       end do
    end do

    minrk = max(minrl, minrk)
    maxrk = min(maxrl, maxrk)

    allocate( evv(n_eig_l, n_eig_r) )
    if (myrank==0) write(*,'(a,2i3)') '  *** calc_gttpe_olp *** ', minrl, maxrk
    if (myrank==0) write(*,'(a,i2,a,i2,a)') &
         ' # <il |  1/sqrt(N) [c+_',iorbl,', c_',iorbr,']^rank | ir > '
    if (myrank==0) write(*,*) '#        Jl(il),  Jr(ir),rank,  1/sqrt(N)*<il|ir>,  N '

    do irank = minrk, maxrk
       call set_ob_channel(irank, iprty, nbody, ops, ij_orb, iorbl, iorbr)
       if (size(ops) /= 1) stop 'error in calc_gtype_olp'

       call init_bp_operator(bp, ops(1))
       call wf_alloc_vec(vt, ptnl)

       do ir = 1, n_eig_r

          call bp_operate(bp, vt, ops(1), evec_r(ir))
          olp = dot_product_global( vt, vt )

          do il = 1, n_eig_l
             r = dot_product_global(evec_l(il), vt)
             evv(il, ir) = r / sqrt(olp)
          end do

       end do
       call deallocate_l_vec(vt%p)
       call finalize_bp_operator(bp, ops(1))

       do il = 1, n_eig_l
          do ir = 1, n_eig_r 
             if (myrank==0) write(*,'(a, 3i5, 2f12.6)') 'NormGTtype', &
                  il, ir, irank, evv(il, ir), olp
          end do
       end do

    end do

    deallocate(evv)

  end subroutine calc_gttype_olp




  subroutine calc_obtd(op_type, func_1, irank, nbody, iprty, e_charge)
    use interaction, only : set_ob_channel
    use model_space, only : mass
    ! one-body transition density 
    character(len=*), intent(in) :: op_type
    real(8), external :: func_1
    integer, intent(in) :: irank, nbody, iprty
    real(8), intent(in), optional :: e_charge(2)
    !
    character(len=maxchar) :: fn_rme, ctmp1, ctmp2
    integer :: i, n, il, ir, jjl, jjr, ipn, i1, i2
    type(opr_m), allocatable :: ops(:)
    integer, allocatable :: ij_orb(:,:)
    real(8) :: x(2), r
    real(8), allocatable :: evvs(:,:,:), rm(:)
    !
    if (.false.) then
       !    if (nbody==1) then
       ! TODO ex_vals_ij 
    else
       call set_ob_channel(irank, iprty, nbody, ops, ij_orb)
       n = size(ops)
       allocate( evvs(n, n_eig_l, n_eig_r), rm(n) )
       do i = 1, n
          call calc_transit_op(ops(i), 0.d0)
          evvs(i, :, :) = evv
       end do
    end if

    if (myrank == 0) then

       do i = 1, n
          i1 = ij_orb(1,i)
          i2 = ij_orb(2,i)
          rm(i) = func_1( &
               norb(i1), lorb(i1), jorb(i1), itorb(i1), &
               norb(i2), lorb(i2), jorb(i2), itorb(i2) )
       end do

       ctmp1 =  fn_load_wave_l(:len_trim(fn_load_wave_l)-4)
       do i = 1, len(ctmp1)
          if ( ctmp1(i:i) == '/' ) ctmp1(i:i) = '_'
       end do
       ctmp2 =  fn_load_wave_r(:len_trim(fn_load_wave_r)-4)
       do i = 1, len(ctmp2)
          if ( ctmp2(i:i) == '/' ) ctmp2(i:i) = '_'
       end do

       fn_rme = 'OBTD_' // trim(op_type) // '_' // &
            trim(ctmp1) // '_' // &
            trim(ctmp2) // '.dat'

       open(lun_rme, file=fn_rme)
       write(lun_rme,*) '#  --- orbit numbers ---'
       write(lun_rme,*) '#   idx      n,   l,  2j,  2tz'
       do i = 1, n_jorb_pn
          write(lun_rme,'(a,i5,a,4i5)') ' # ', i, '  ', &
               norb(i), lorb(i), jorb(i), itorb(i)
       end do
       write(lun_rme, '(/,a,i3,a,i3)') &
            ' # rank ', irank, ' parity ', iprty
       write(lun_rme, '(/,a)') ' # ' // trim(op_type) &
            // ' reduced metrix element' 
       write(lun_rme, '(a,i5)') ' # of elements = ',n
       if (present(e_charge)) write(lun_rme, '(a,2f11.5)') &
            ' # eff. charge = ', e_charge
       write(lun_rme,*)

       do il = 1, n_eig_l
          jjl = evec_l(il)%jj
          do ir = 1, n_eig_r
             jjr = evec_r(ir)%jj
             if (     jjl + jjr  < 2*irank .or. &
                  abs(jjl - jjr) > 2*irank ) cycle
             write(lun_rme, '(a,i3,a,i5,a,i3,a,i5,a)') &
                  'w.f.  J1=', jjl, '/2(', il,&
                  ')     J2=', jjr, '/2(', ir, ')'
             x = 0.d0
             do i = 1, n
                ipn = (itorb(ij_orb(1,i)) + 3)/2
                x(ipn) = x(ipn) + evvs(i,il,ir) * rm(i)
             end do
             if (present(e_charge)) then
                r = sum(x * e_charge)**2
             else
                r = sum(x)**2
             end if
             write(lun_rme, '(a,2f12.5)') &
                  ' B(' // trim(op_type) &
                  // ';=>), B(' // trim(op_type) // ' ;<=) ', &
                  r/dble(jjl + 1), r/dble(jjr + 1)
             write(lun_rme, '(a,2i5,2f12.5)') &
                  ' <||'//trim(op_type)//'||>  ', il, ir, x

             if (op_type=='E0') write(lun_rme, '(a,f12.5,a,2f12.5/)') &
                  ' rho^2*1000  = ', r / (1.2d0 * mass**(1.d0/3.d0))**4 * 1000.d0, &
                  '    rho_p, rho_n = ', x(:) / (1.2d0 * mass**(1.d0/3.d0))**2

             write(lun_rme, '(/,a)') '  i  j      OBTD    <i||' &
                  // trim(op_type) // '||j>  OBTD*<||>'
             do i = 1, n
                write(lun_rme,'(2i3,3f12.5)') ij_orb(:,i), &
                     evvs(i,il,ir), rm(i), evvs(i,il,ir)*rm(i)
             end do
             write(lun_rme, '(/,/)')
          end do
       end do
       close(lun_rme)

    end if

    deallocate(ij_orb, ops)

  end subroutine calc_obtd


  subroutine calc_fm()
    ! Fermi transition
    character(len=maxchar) :: comment
    call set_fm()
    write(comment,'(a,2i3)') &
         " Fermi transition  <f||1 t^(+-)||i>    parity", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(fm_m, 0.d0, comment)
  end subroutine calc_fm



  subroutine calc_transit_op(op, c_mom, comment)
    ! transition for general operator
    type(opr_m), intent(inout) :: op
    character(len=*), optional :: comment
    real(8), intent(in) :: c_mom ! coeff. for moment
    integer :: i, j, jl, jr
    real(8) :: x
    logical,allocatable  :: is_mup(:,:)
    real(8), allocatable :: fc(:,:)
    type(type_vec_p) :: vt

    evv(:,:) = 0.d0

    allocate( is_mup(n_eig_l, n_eig_r), fc(n_eig_l, n_eig_r))
    is_mup(:,:) = .false.
    fc(:,:) = 0.d0

    call init_bp_operator(bp, op)
    call wf_alloc_vec(vt, bp%ptnl)

    do j = 1, n_eig_r
       jr = evec_r(j)%jj
       if (jr < 0) cycle
       do i = 1, n_eig_l
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          if (abs(jl-jr) > op%irank*2 .or. op%irank*2 > jl+jr ) cycle
          
          x = dcg(jr, mtotr, op%irank*2, mtotl-mtotr, jl, mtotl)
          
          if (abs(x) < 1.d-8) then
             if (jr >= mtotr+2) is_mup(i,j) = .true.
             cycle
          end if

          fc(i,j) = sqrt(dble(jl+1)) / x
       end do
    end do


    do j = 1, n_eig_r
       if ( all( fc(:,j) == 0.d0 ) ) cycle
       call bp_operate(bp, vt, op, evec_r(j))
       do i = 1, n_eig_l
          if ( fc(i,j) == 0.d0 ) cycle
          evv(i,j) = fc(i,j) * dot_product_global(evec_l(i), vt)
       end do
    end do

    call deallocate_l_vec(vt%p)
    deallocate(fc)

    call finalize_bp_operator(bp, op)

    ! avoid <J, 0, lambda, 0 | J 0 > = 0
    if (any(is_mup)) call calc_mup(is_mup, op)
    
    deallocate( is_mup )

    evm = 0.d0
    if (mtotl == mtotr .and. c_mom /= 0.d0) then
       do i = 1, min(n_eig_l, n_eig_r)
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          evm(i,i) = evv(i,i) * c_mom &
               * dcg(jl, jl, op%irank*2, 0, jl, jl) / sqrt(jl+1.d0)
       end do
    end if

    if (present(comment) .and. myrank==0) then 
       write(*,'(/,a)') trim(comment)
       call print_trans(evv, evm, op%irank)
    end if

  end subroutine calc_transit_op


  subroutine calc_mup(is_mup, op)
    !
    ! |ini J, M+1> = 1/N * J+ |ini J, M>
    !   <fin |op|ini J, M+1> at evv
    !
    logical, intent(in) :: is_mup(:,:)
    type(opr_m), intent(inout) :: op
    integer :: i, j, jl, jr
    real(8) :: x
    type(type_vec_p) :: vt

    call finalize_bridge_partitions(bp)

    call set_mup(is_mup)

    call init_bridge_partitions(bp, ptnl, ptnm)

    call init_bp_operator(bp, op)
    ! 
    call wf_alloc_vec(vt, bp%ptnl)

    do j = 1, n_eig_r
       jr = evec_r(j)%jj
       if ( .not. any(is_mup(:,j)) ) cycle
       call bp_operate( bp, vt, op, evec_m(j) )

       do i = 1, n_eig_l
          jl = evec_l(i)%jj
          if ( .not. is_mup(i,j) ) cycle

          x = dcg(jr, mupr, op%irank*2, mtotl-mupr, jl, mtotl)
          if ( abs(x) < 1.d-8 ) then 
             if ( myrank == 0 ) write(*,'(a,4i3)') &
                  'WARNING: J+^2|> required, skip ',i,jl,j,jr
             cycle
          end if
          evv(i,j) = dot_product_global(evec_l(i), vt) * sqrt(dble(jl+1)) / x 
       end do
    end do

    call deallocate_l_vec(vt%p)
    
    call finalize_bp_operator(bp, op)
    call finalize_bridge_partitions(bp)

    call init_bridge_partitions(bp, ptnl, ptnr)

  end subroutine calc_mup



  
  subroutine set_mup(is_mup)
    !
    ! |ini J, M+1> = 1/N * J+ |ini J, M> at evec_m, ptnm
    !
    logical, intent(in) :: is_mup(:,:)
    integer :: i, j, jl, jr, nj, jlist(n_eig_r)
    real(8) :: x
    type(type_vec_p) :: vt

    if ( .not. allocated(evec_m) ) then
       mupr = mtotr + 2
       allocate( evec_m(n_eig_r) )

       open(lunptn, file=fn_ptn_r)
       call init_partition(ptnm, lunptn, mupr, verbose=.false.)
       close(lunptn)

       allocate( cost(ptnm%n_pidpnM) )
       call cost_from_localdim(ptnm, ptnm%pidpnM_pid_srt, nprocs, cost)
       call deploy_partition(ptnm, cost, verbose=.false.)
       deallocate(cost)
    end if

    nj = 0
    do j = 1, n_eig_r
       if ( associated(evec_m(j)%p) ) cycle
       if ( .not. any(is_mup(:,j)) ) cycle
       if ( mupr > evec_r(j)%jj ) cycle
       nj = nj + 1
       jlist(nj) = j
    end do

    if (nj == 0) return

    call init_bridge_partitions(bp, ptnm, ptnr)
    call init_bp_operator(bp, jtensor, verbose=.false.)

    if (myrank == 0) write(*,'(a)', advance='no') "mup called   w.f.: "

    do i = 1, nj
       j = jlist(i)
       call allocate_l_vec( evec_m(j)%p, ptnm%max_local_dim )
       call bp_operate( bp, evec_m(j), jtensor, evec_r(j) )
       x = dot_product_global( evec_m(j), evec_m(j) )
       evec_m(j)%p = - 1.d0/sqrt(x) * evec_m(j)%p ! phase : tensor->ladder
       if (myrank == 0) write(*,'(i4)', advance='no') j
    end do
    
    if (myrank == 0) write(*,*)

    call finalize_bp_operator(bp, jtensor)
    call finalize_bridge_partitions(bp)

  end subroutine set_mup


  
  subroutine finalize_mup()
    integer :: i
    if ( .not. allocated(evec_m) ) return
    call finalize_partition(ptnm)
    do i = 1, size(evec_m)
       if (associated(evec_m(i)%p)) call deallocate_l_vec(evec_m(i)%p)
    end do
  end subroutine finalize_mup



  subroutine print_mup_called(is_mup, comment)
    ! print J+|> needed
    logical, intent(in) :: is_mup(:,:)
    character(len=*), intent(in) :: comment
    integer :: i, j

    if (myrank /= 0) return
    write(*,'(/,a)',advance='no') trim(comment)
!    do i = 1, min(n_eig_l, 10)
!       do j = 1, min(n_eig_r, 10)
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          if (.not. is_mup(i,j)) cycle
          write(*,'(a,i3,a,i3,a)',advance='no') '  (',i,',',j,')  '
       end do
    end do
    write(*,*)
    write(*,*)
  end subroutine print_mup_called



  subroutine calc_ff()
    ! first forbidden transition
    use interaction, only : dummy_func1
    integer :: i, j, jl, jr
    character(len=maxchar) :: comment
    integer :: nbody


    ! ff_0rs, ff_0sp, ff_1r, ff_1rs, ff_1p, ff_2rs
    call set_ff()
    
    
    if (is_obtd) then
       nbody = -10
       if ( ptnl%n_ferm(1)-ptnr%n_ferm(1) == -1)  nbody = -11
       call calc_obtd('FF', dummy_func1, 0, nbody, -1)
       call calc_obtd('FF', dummy_func1, 1, nbody, -1)
       call calc_obtd('FF', dummy_func1, 2, nbody, -1)
    end if


    ! Rank 0
    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 0, M0rs=<f||[rC1 x sigma]^(0)t^-||i>   ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_0rs, 0.d0, comment)


    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 0, M0sp=<f||[sigma x nabla]^(0)t^-||i> ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_0sp, 0.d0, comment)

    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 0, M0rsd=<f||2/3I[rC1 x sigma]^(0)t^-||i> ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_0rs_d, 0.d0, comment)

    ! Rank 1
    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 1, M1r =<f||[rC1 t^-||i>               ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_1r,  0.d0, comment)

    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 1, M1rs=<f||[rC1 x sigma]^(1)t^-||i>   ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_1rs, 0.d0, comment)

    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 1, M1p =<f||nabla t^-||i>              ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_1p, 0.d0, comment)

    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 1, M1rd =<f||2/3I[rC1 t^-||i>          ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_1r_d,  0.d0, comment)

    write(comment,'(a,2i3)') &
         "First Forbidden, Rank 1, M1rsd=<f||2/3I[rC1 x sigma]^(1)t^-||i>   ", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_1rs_d, 0.d0, comment)


    ! Rank 2
    write(comment,'(a,2i3)') &
         "First forbidden, Rank 2, M2rs=<f|| [rC1 x sigma]^(2) t^- ||i>", &
         ptnl%iprty, ptnr%iprty
    call calc_transit_op(ff_2rs, 0.d0, comment)




  end subroutine calc_ff



  subroutine calc_sf(ipn)
    ! one-particle Spectroscopic factor
    use model_space, only : char_orbit, n_orb_char
    integer,  intent(in) :: ipn
    integer :: i, j, mm, n, k, jl, jr, iprty, no, nj
    type(opr_m) :: crt
    character(n_orb_char) :: corb
    real(8) :: x
    logical, allocatable :: is_mup(:,:,:)
    real(8), allocatable :: cg(:,:,:)

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty

    if (myrank==0) then 
       write(*,'(/,/,/,a)') ' ------------------------------------ '
       write(*,'(a)')       ' --- Spectroscopic factor , C^2*S --- '
       write(*,'(a,/)')     ' ------------------------------------ '
       if (ipn==1) then
          write(*,*) trim(fn_load_wave_l), "=", trim(fn_load_wave_r), " + p"
       else
          write(*,*) trim(fn_load_wave_l), "=", trim(fn_load_wave_r), " + n"
       end if
       write(*,*)
    end if


    ! 
    allocate( is_mup( n_eig_l, n_eig_r, n_jorb(ipn) ) )
    allocate( cg(     n_eig_l, n_eig_r, n_jorb(ipn) ) )
    is_mup(:,:,:) = .false.
    cg(    :,:,:) = 0.d0

    nj = 0
    if (ipn==2) nj = n_jorb(1)

    do n = 1, n_jorb(ipn)
       k = n + nj
       if ( iporb(k) /= iprty ) cycle
       do i = 1, n_eig_l
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          do j = 1, n_eig_r
             jr = evec_r(j)%jj
             if (jr < 0 .or. abs(jl - jr) > jorb(k) .or. jl+jr < jorb(k)) cycle
             
             if (abs(mm) > jorb(k)) then
                is_mup(i, j, n) = .true.
                cycle
             end if
             
             x = dcg(jr, mtotr, jorb(k), mtotl-mtotr, jl, mtotl)
             if (abs(x) < 1.d-8) then
                is_mup(i, j, n) = .true.
                cycle
             end if
             
             cg(i, j, n) = x
          end do
       end do
    end do


    do n = 1, n_jorb(ipn)
       k = n + nj
       if ( iporb(k) /= iprty ) cycle

       call char_orbit(norb(k), lorb(k), jorb(k), itorb(k), corb)
       if (myrank==0) then
          write(*,'(/,/,a,/)') '         n  l 2j 2tz'
          write(*,'(1a, 4i3, 2a)') "orbit : ", &
               norb(k), lorb(k), jorb(k), itorb(k), '    ', corb
          write(*,*)
       end if
       if (abs(mm) <= jorb(k)) then 
          call opr_m_one_crt(crt, k, mm)
          call init_bp_operator(bp, crt)
       end if

       if (myrank==0) then 
          write(*,'(a,2i3)') " S-factor parity", ptnl%iprty, ptnr%iprty
          write(*,*) '2xJf      Ef      2xJi     Ei       Ex       C^2*S '
       end if

       evv = 0.d0
       do i = 1, n_eig_l
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          do j = 1, n_eig_r
             jr = evec_r(j)%jj
             if (jr < 0 .or. abs(jl - jr) > jorb(k) .or. jl + jr < jorb(k)) cycle

             ! x = dcg(jr, mtotr, jorb(k), mtotl-mtotr, jl, mtotl)
             ! if (abs(x) < 1.d-8) cycle
             if (is_mup(i, j, n)) then
                if (myrank==0) write(*, '(2(i2,"(",i4,")",f9.3), f8.3, a)') &
                     evec_l(i)%jj, i, evec_l(i)%eval, &
                     evec_r(j)%jj, j, evec_r(j)%eval, &
                     evec_r(j)%eval - evec_l(i)%eval, &
                     "    ------"
                cycle
             end if
             
             call ex_val(bp, evec_l(i), crt, evv(i,j), evec_r(j))
             ! evv(i,j) = (evv(i,j) * sqrt(dble(jl+1)) / x)**2  / dble(jl+1)
             evv(i,j) = (evv(i,j) / cg(i, j, n))** 2 

             if (myrank==0) write(*, '(2(i2,"(",i4,")",f9.3), f8.3, f10.4)') &
                  evec_l(i)%jj, i, evec_l(i)%eval, &
                  evec_r(j)%jj, j, evec_r(j)%eval, &
                  evec_r(j)%eval - evec_l(i)%eval, &
                  evv(i,j)

          end do
       end do

       ! do i = 1, n_eig_l
       !    write(*,*)"sum_r", i, sum(evv(i,:))
       ! end do
       ! do j = 1, n_eig_r
       !    write(*,*)"sum_l", j, sum(evv(:,j))
       ! end do

       if (abs(mm) <= jorb(k)) call finalize_bp_operator(bp, crt)
    end do
    if (myrank==0) then
       if (any(is_mup)) write(*,'(/,a,/)') " *** use mup_operate.exe to calculate '------' *** "
    end if
    if (myrank==0) write(*,'(/,/)')

          
  end subroutine calc_sf


  subroutine calc_tna(ipn)
    ! two-nucleon transfer ampletude
    use model_space, only : char_orbit, n_orb_char
    use operator_mscheme, only: opr_m_two_crt, finalize_opr_m
    integer,  intent(in) :: ipn
    integer :: i, j, mm, n, k1, k2, j1, j2, irank, jl, jr, iprty, no, kk(4)
    type(opr_m) :: op
    character(n_orb_char) :: corb1, corb2
    real(8) :: x

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty

    if (myrank==0) then 
       write(*,'(/,/,/,a)') ' ------------------------------------ '
       write(*,'(a)')       ' --- TNA for 2-particle s-factor  --- '
       write(*,'(a,/)')     ' ------------------------------------ '
    end if

    select case (ipn)
    case (1) 
       kk = (/           1, n_jorb(1),             1, n_jorb(1) /)
       if (myrank==0) write(*,*) trim(fn_load_wave_l)," = ",trim(fn_load_wave_r)," + 2p "
    case (2) 
       kk = (/ n_jorb(1)+1, n_jorb_pn,   n_jorb(1)+1, n_jorb_pn /) 
       if (myrank==0) write(*,*) trim(fn_load_wave_l)," = ",trim(fn_load_wave_r)," + 2n "
    case (3)
       kk = (/           1, n_jorb(1),   n_jorb(1)+1, n_jorb_pn /) 
       if (myrank==0) write(*,*) trim(fn_load_wave_l)," = ",trim(fn_load_wave_r)," + pn "
    case default
       stop "ERROR in calc_tna"
    end select

    if (myrank==0) then
       write(*,*)
       write(*,'(a,2i3)') " TNA,  parity", ptnl%iprty, ptnr%iprty
       write(*,*) '2xJf      Ef      2xJi     Ei       Ex       TNA '
    end if

    do k1 = kk(1), kk(2)
       j1 = jorb(k1)
       do k2 = kk(3), kk(4)
          if (k1 > k2) cycle
          j2 = jorb(k2)
          if ( iporb(k1) * iporb(k2) /= iprty ) cycle

          do irank = abs(j1-j2)/2, (j1+j2)/2
             if (k1 == k2 .and. mod(irank, 2)==1) cycle
             if (abs(mm) > irank*2) cycle

             no = 0
             do i = 1, n_eig_l
                jl = evec_l(i)%jj
                if (jl < 0) cycle
                do j = 1, n_eig_r
                   jr = evec_r(j)%jj
                   if (jr < 0 .or. abs(jl - jr) > irank*2) cycle
                   x = dcg(jr, mtotr, irank*2, mtotl-mtotr, jl, mtotl)
                   if (abs(x) < 1.d-8) cycle 
                   no = no + 1
                end do
             end do
             if (no == 0) cycle

             call char_orbit(norb(k1), lorb(k1), jorb(k1), itorb(k1), corb1)
             call char_orbit(norb(k2), lorb(k2), jorb(k2), itorb(k2), corb2)
             if (myrank==0) then
                write(*,'(/,a,i3,x,2a,i3,x,2a,i3,a)') &
                     'TNA : [ ', k1, corb1, ' x ', k2, corb2, ' ]^(', irank, ')'
                write(*,*)
             end if

             call opr_m_two_crt( op, k1, k2, irank )
             call init_bp_operator(bp, op) 

             evv = 0.d0
             do i = 1, n_eig_l
                jl = evec_l(i)%jj
                if (jl < 0) cycle
                do j = 1, n_eig_r
                   jr = evec_r(j)%jj
                   if (jr < 0 .or. abs(jl - jr) > irank*2) cycle
                   x = dcg(jr, mtotr, irank*2, mtotl-mtotr, jl, mtotl)

                   if (abs(x) < 1.d-8) cycle

                   call ex_val(bp, evec_l(i), op, evv(i,j), evec_r(j)) 
                   evv(i,j) = evv(i,j) / x

                   if (myrank==0) write(*, &
                                ! '(2(i2,"(",i4,")",f9.3), f8.3, 2f10.4)' ) &
                        '(2(i2,"(",i4,")",f9.3), f8.3, f11.5)' ) &
                        evec_l(i)%jj, i, evec_l(i)%eval, &
                        evec_r(j)%jj, j, evec_r(j)%eval, &
                        evec_r(j)%eval - evec_l(i)%eval, &
                        evv(i,j)


                end do
             end do

             call finalize_bp_operator(bp, op)
             call finalize_opr_m( op ) 

          end do

       end do
    end do
    if (myrank==0) write(*,'(/,/)')

  end subroutine calc_tna




  subroutine calc_nme_0v(ipn)
    ! Nuclear matrix element for neutrinoless double-beta decay
    use model_space, only : char_orbit
    use interaction, only : set_nme_0v, nme_0v, set_stst0
    integer,  intent(in) :: ipn
    integer :: i, j, mm, n, k, jl, jr, iprty
    real(8) :: x

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty
    if (mm /= 0 .or. iprty /= 1) return

    if (ipn/=1) stop 'not implemented calc_nme_0v'

    call set_nme_0v(ipn)
    !    call set_stst0(nme_0v)
    call init_bp_operator(bp, nme_0v, verbose=.true.)

    if (myrank==0) then 
       write(*,*)
       write(*,*) ' Nuclear matrix elements '
       write(*,*) '    for neutrinoless double-beta decay'
       write(*,*) trim(fn_load_wave_l), " | cp+ cp+ cn cn | ", trim(fn_load_wave_r)
       write(*,*)
    end if

    evv = 0.d0
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > nme_0v%irank*2) cycle
          x = dcg(jr, mtotr, nme_0v%irank*2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call ex_val(bp, evec_l(i), nme_0v, evv(i,j), evec_r(j))
          !          evv(i,j) = evv(i,j) * sqrt(dble(jl+1)) / x
       end do
    end do

    evm = 0.d0
    if (myrank==0) then
       write(*,*)
       write(*,'(a)') " Nuclear matrix element shown at Mred. (NOT reduced m.e.)"
       call print_trans(evv, evm, nme_0v%irank)
    end if
    if (myrank==0) write(*,*)
    call finalize_bp_operator(bp, nme_0v)

  end subroutine calc_nme_0v




  subroutine print_trans(evv, evm, irank)
    real(8), intent(in) :: evv(:,:), evm(:,:)
    integer, intent(in) :: irank
    integer :: i, j, jl, jr
    real(8) :: x, y

    write(*,*) '2xJf      Ef      2xJi     Ei       Ex       Mred.    B(EM )->   B(EM)<-   Mom.'
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl-jr) > irank*2 &
               .or. irank*2 > jl+jr ) cycle
          x = evv(i, j)
          y = evm(i, j)
          !          if (x == 0.d0 .and. y == 0.d0) cycle
          write(*, '(2(i2,"(",i4,")",f9.3), f8.3,5f10.4)') &
               evec_l(i)%jj, i, evec_l(i)%eval, &
               evec_r(j)%jj, j, evec_r(j)%eval, &
               evec_r(j)%eval - evec_l(i)%eval, &
               x, x**2/dble(evec_l(i)%jj+1), x**2/dble(evec_r(j)%jj+1), y
       end do
    end do
    write(*,*)
  end subroutine print_trans


  subroutine print_reduced_me(c1, c2, c3, c4, ev1, ev2, ev3, ev4)
    character(len=*), intent(in), optional, target :: c1, c2, c3, c4
    real(8), intent(in), optional, target :: &
         ev1(:,:,:,:), ev2(:,:,:,:), ev3(:,:,:,:), ev4(:,:,:,:)
    integer :: i, j, iop, nop=4, n
    character(len=:), pointer :: cp
    real(8), pointer :: ev(:,:)

    write(*,'(/,a,/)') "reduced matrix elements"
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          do iop = 1, nop
             select case (iop)
             case (1)
                if (.not. present(c1 )) cycle
                if (.not. present(ev1)) cycle
                cp => c1
                ev => ev1(:,:,i,j)
             case (2)
                if (.not. present(c2 )) cycle
                if (.not. present(ev2)) cycle
                cp => c2
                ev => ev2(:,:,i,j)
             case (3)
                if (.not. present(c3 )) cycle
                if (.not. present(ev3)) cycle
                cp => c3
                ev => ev3(:,:,i,j)
             case (4)
                if (.not. present(c4 )) cycle
                if (.not. present(ev4)) cycle
                cp => c4
                ev => ev4(:,:,i,j)
             case default
                stop "not implemented print_reduced_me"
             end select
             write(*,'("<",i5,f9.3,"||",a,"||",i5,f9.3,"> (Ex.", &
                  & f9.3,")     ",f12.5)') &
                  i, evec_l(i)%eval, trim(cp), j, evec_r(j)%eval, &
                  evec_r(j)%eval - evec_l(i)%eval, sum(ev(:,:))
             write(*,'(1000i8)', advance='no') &
                  (/( n, n=1, size(ev,1) )/)
             write(*,'(a)') '  sum '
             do n = 1, size(ev,1)
                write(*,'(i3,1000f8.4)') n, ev(n,:), sum(ev(n,:))
             end do
             write(*,'(a,1000f8.4)') 'sum', &
                  (/( sum(ev(:,n)), n=1, size(ev,2) )/), sum(ev)
          end do
       end do
    end do
  end subroutine print_reduced_me


  subroutine check_mup(jl, ml, jr, mr, is_mup )
    integer, intent(in) :: jl, ml, jr, mr
    logical, intent(out) :: is_mup
    integer :: jj, kk

    is_mup = .false.
    if (jl < 0 .or. (jr < 0)) return

    do jj = abs(jl-jr), jl+jr, 2

       if ( abs(dcg(jr, mr, jl, -ml, jj, mr-ml)) < thd_mup  &
            .and. jr >= mr+2) then 

          do kk = abs(jl-jr), jl+jr, 2
             if ( abs(dcg(jr, mr+2, jl, -ml, jj, mr+2-ml)) < thd_mup) then 
                write(*,*)'increase M twice', jl, mr, jl, -ml, jj, mr-ml
                stop 
             end if
          end do

          is_mup = .true. 
          return

       end if
    end do

  end subroutine check_mup


  subroutine calc_tbtd()
    ! two-body transition density and one-body transition density
    !  in any rank
    use operator_mscheme, only : init_tbtd_op, clear_op, &
         get_cpld_tbtd, get_cpld_obtd, finalize_opr_m, &
         init_idx_tbtd, finalize_idx_tbtd, tbtd_container
    use bp_expc_val, only : bp_ex_val_tbtd 
    integer :: i, j, jl, ml, jr, mr, mm, iprty
    real(8) :: x, c
    character(len=maxchar) :: comment
    type(opr_m) :: op
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, j12, j34, jj, mu, nb
    integer :: ii, ipn, kst(4)

    type(tbtd_container), allocatable :: tbtdlr(:,:)

    ! m-up
    type(type_ptn_pn) :: ptnm
    type(type_vec_p), allocatable :: evec_m(:)
    logical :: is_mup(n_eig_l, n_eig_r)

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty

    call init_tbtd_op(op, mm, iprty, ptnl%n_ferm, ptnr%n_ferm)

    if (myrank==0) then 
       write(*,*) 
       write(*,*) '--------------------------------------'
       write(*,*) '      Two-body transition density   '
       write(*,*) trim(fn_load_wave_l), " | TBTD | ", trim(fn_load_wave_r)
       write(*,'(a,i3,a,i2)') ' operator MM= ', mm, ' , parity =', iprty
       write(*,'(a,i3)') '   op. type= ', op%nbody
       write(*,*) '--------------------------------------'
       write(*,*)
    end if

    allocate( tbtdlr(n_eig_l, n_eig_r) )
    nb = op%nbody ! ad hoc
    !$omp parallel do private(i, j)
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          if (.not. is_calc_lr(i,j)) cycle
          
          call init_idx_tbtd( evec_l(i)%jj, evec_r(j)%jj, iprty,  nb, &
               tbtdlr(i, j) )
          call check_mup( evec_l(i)%jj, ptnl%mtotal, &
               evec_r(j)%jj, ptnr%mtotal, is_mup(i, j) )
          
       end do
    end do

    call init_bp_operator( bp, op, verbose=.false. )


    do i = 1, n_eig_l
       do j = 1, n_eig_r

          if (.not. is_calc_lr(i,j)) cycle
          call inside(ptnr%mtotal, evec_r)

       end do
    end do

    call finalize_bp_operator(bp, op)
    call finalize_opr_m(op)


    if (any(is_mup)) then 
       if (myrank==0) write(*,*) ' *** called mup ***'

       call finalize_bridge_partitions(bp)

       mu = ptnr%mtotal + 2 
       open(lunptn, file=fn_ptn_r)
       call init_partition(ptnm, lunptn, mu)
       close(lunptn)

       allocate( cost(ptnm%n_pidpnM) )
       call cost_from_localdim( ptnm, ptnm%pidpnM_pid_srt, nprocs, cost )
       call deploy_partition(   ptnm, cost, verbose=.false.  )
       deallocate(cost)


       mm = ptnl%mtotal - ptnm%mtotal
       iprty = ptnl%iprty * ptnm%iprty

       call init_tbtd_op(op, mm, iprty, ptnl%n_ferm, ptnr%n_ferm)

       call init_bridge_partitions(bp, ptnm, ptnr)
       call init_bp_operator(bp, jtensor, verbose=.true.)
       allocate( evec_m(n_eig_r) )

       do j = 1, n_eig_r
          if (.not. any(is_mup(:,j))) cycle
          call allocate_l_vec( evec_m(j)%p, ptnm%max_local_dim )
          call bp_operate(bp, evec_m(j), jtensor, evec_r(j))
          x = dot_product_global( evec_m(j), evec_m(j))
          evec_m(j)%p = - 1.d0/sqrt(x) * evec_m(j)%p ! phase : tensor->ladder
          evec_m(j)%jj = evec_r(j)%jj
       end do

       call finalize_bp_operator(bp, jtensor)
       call finalize_bridge_partitions(bp)

       call init_bridge_partitions(bp, ptnl, ptnm)
       call init_bp_operator(bp, op, verbose=.true.)

       do i = 1, n_eig_l
          do j = 1, n_eig_r
             if (.not. is_calc_lr(i,j)) cycle
             if (.not. is_mup(i,j)) cycle

             call inside(ptnm%mtotal, evec_m)

          end do
       end do

       call finalize_bp_operator(bp, op)
       call finalize_opr_m(op)
       call finalize_bridge_partitions(bp)
       call finalize_partition(ptnm)
       call init_bridge_partitions(bp, ptnl, ptnr)
       do i = 1, size(evec_m)
          if (associated(evec_m(i)%p)) call deallocate_l_vec(evec_m(i)%p)
       end do
    end if

    call start_stopwatch(time_print_tbtd)
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          if (.not. is_calc_lr(i,j)) cycle
          call print_res(i, j, tbtdlr(i, j))
       end do
    end do
    call stop_stopwatch(time_print_tbtd)

    !$omp parallel do private(i, j)
    do j = 1, n_eig_r
       do i = 1, n_eig_l
          if (.not. is_calc_lr(i,j)) cycle
          call finalize_idx_tbtd( tbtdlr(i, j) )
       end do
    end do

    deallocate( tbtdlr )

  contains


    subroutine inside(mr, evec_r)
      integer, intent(in) :: mr
      type(type_vec_p), intent(inout) :: evec_r(:)
      real(8) :: tt

      jl = evec_l(i)%jj
      if (jl < 0) return
      ml = ptnl%mtotal

      jr = evec_r(j)%jj
      if (jr < 0) return

      call clear_op( op )

      call bp_ex_val_tbtd( bp, evec_l(i), evec_r(j), op )

      ! call print_tbtd_op( op )

      call start_stopwatch(time_cpld_tbtd)

      !$omp parallel do private(ii, k1, k2, jj, c, x, tt) schedule(dynamic)
      do ii = 1, tbtdlr(i,j)%n12
         if ( tbtdlr(i,j)%v12jj(ii) /= 0.d0 ) cycle
         ! tt =  tbtdlr(i,j)%v12jj(ii) 
         
         k1  = tbtdlr(i,j)%k12jj(1, ii)
         k2  = tbtdlr(i,j)%k12jj(2, ii)
         jj  = tbtdlr(i,j)%k12jj(3, ii)

         if ( abs(jr-jl) > 2*jj .or. jr+jl < 2*jj ) cycle
         if ( abs(mr-ml) > 2*jj ) cycle

         c = dcg( jr, mr, jl, -ml, jj*2, mr-ml )
         if (abs(c) < thd_mup) cycle

         c = (-1d0)**( (jr-mr)/2 ) / c

         call get_cpld_obtd( op, k1, k2, jj, x )
         tbtdlr(i,j)%v12jj(ii) = c * x
         
         !if (tt /= 0.d0) then
         !   if (abs(tbtdlr(i,j)%v12jj(ii) - tt)>1.d-3) &
         !        write(*,*) "warning diff mup ob", tt, tbtdlr(i,j)%v12jj(ii) 
         !end if

      end do


      !$omp parallel do private(ii, k1, k2, k3, k4, j12, j34, jj, c, x, tt) &
      !$omp & schedule(dynamic)
      do ii = 1, tbtdlr(i,j)%n1234
         ! if ( tbtdlr(i,j)%v1234jj(ii) /= 0.d0 ) cycle
         tt =  tbtdlr(i,j)%v1234jj(ii)
         
         k1  = tbtdlr(i,j)%k1234jj(1, ii)
         k2  = tbtdlr(i,j)%k1234jj(2, ii)
         k3  = tbtdlr(i,j)%k1234jj(3, ii)
         k4  = tbtdlr(i,j)%k1234jj(4, ii)
         j12 = tbtdlr(i,j)%k1234jj(5, ii)
         j34 = tbtdlr(i,j)%k1234jj(6, ii)
         jj  = tbtdlr(i,j)%k1234jj(7, ii)

         if ( abs(mr-ml) > 2*jj ) cycle

         c = dcg( jr, mr, jl, -ml, jj*2, mr-ml )
         if (abs(c) < thd_mup) cycle

         c = (-1d0)**( (jr-mr)/2 ) / c
         call get_cpld_tbtd( op, k1, k2, k3, k4, j12, j34, jj, x )

         tbtdlr(i,j)%v1234jj(ii) = c * x

         if (tt /= 0.d0) then
           if (abs(tbtdlr(i,j)%v1234jj(ii) - tt) > 1.d-3) then
              if (myrank==0) write(*,'(a,2f10.5,4i3,a,3i3,a,2i3)') &
                   "warning  mup tb diff.", tt, tbtdlr(i,j)%v1234jj(ii),&
                   k1,k2,k3,k4,' |',j12,j34,jj, ' |',i,j
              ! try "tol=1.d-8" for better precision w.f.
           end if
           tbtdlr(i,j)%v1234jj(ii) = tt ! better result 
         end if
         
      end do

      call stop_stopwatch(time_cpld_tbtd)

    end subroutine inside



    subroutine print_res(i, j, ct)
      integer, intent(in) :: i, j
      integer :: jl, jr
      type(tbtd_container), intent(in) :: ct
      if (myrank /= 0) return

      jl = evec_l(i)%jj
      jr = evec_r(j)%jj
      if (ct%n12 > 0) then
         write(*, '(/,a,i3,a,i5,a,i3,a,i5,a)') &
              'w.f.  J1=', jl, '/2(', i,&
              ')     J2=', jr, '/2(', j, ')'
         write(*,'(a,i5)') &
              '        i  j :Rank:    w.f.  :      OBTD           ',  ct%n12
      end if
      do ii = 1, ct%n12
         if (ct%v12jj(ii) == 0.d0) cycle
         write(*,'(a,2i3,a,i3,a,2i4,a,f14.7)') 'OBTD: ', &
              ct%k12jj(1:2,ii), ' :', ct%k12jj(3,ii), ' : ', &
              i, j, ' : ', ct%v12jj(ii)
      end do

      if (ct%n1234 > 0) then
         write(*, '(/,a,i3,a,i5,a,i3,a,i8,a)') &
              'w.f.  J1=', jl, '/2(', i,&
              ')     J2=', jr, '/2(', j, ')'
         write(*,'(a,i5)') &
              '        i  j  k  l :J_ij J_kl Rank: w.f. :     TBTD    ', ct%n1234
      end if
      do ii = 1, ct%n1234
         if (ct%v1234jj(ii) == 0.d0) cycle
         write(*,'(a,4i3,a,3i3,a,2i4,a,f14.7)') &
              'TBTD: ', ct%k1234jj(1:4,ii), ' :', &
              ct%k1234jj(5:7,ii), ' : ', i, j, ' : ', ct%v1234jj(ii)
      end do
    end subroutine print_res
    

  end subroutine calc_tbtd


  
end module mod_transit



program transit
  use mod_transit

  call transit_main()

end program transit



