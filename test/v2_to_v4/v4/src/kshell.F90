!
! m-scheme shell model code with partition truncation
! ./kshell.exe foobar.input  
!   or 
! mpirun ./kshell.exe foobar.input
!

module kshell_func
  ! module to use subroutine as dummy object
  use constant, only: kwf
  use operator_mscheme, only: opr_m
  use bridge_partitions, only: type_bridge_partitions, bp_operate
  use bp_block, only:  bp_operate_block
  use wavefunction, only: dot_product_global, type_vec_p
  implicit none
  private
  public :: set_kshell_func, matvec, dotprod, matvec_jj, &
       matvec_block, matvec_block_jj
  type(type_bridge_partitions), pointer :: wf_save
  type(opr_m), pointer :: hamltn_save, j_square_save
  
contains

  subroutine set_kshell_func(wf, hamltn, j_square)
    type(type_bridge_partitions), target :: wf
    type(opr_m), target :: hamltn, j_square
    wf_save => wf
    hamltn_save => hamltn
    j_square_save => j_square
  end subroutine set_kshell_func

  subroutine matvec(v1, v2)
    ! in MPI v1 might be broken by shift communication
    type(type_vec_p), intent(inout) :: v1
    type(type_vec_p), intent(inout) :: v2
    call bp_operate(wf_save, v2, hamltn_save, v1)
  end subroutine matvec

  subroutine dotprod(v1, v2, r)
    type(type_vec_p), intent(in) :: v1
    type(type_vec_p), intent(in) :: v2
    real(8), intent(out) :: r
    r = dot_product_global(v1, v2)
  end subroutine dotprod

  subroutine matvec_jj(v1, v2)
    type(type_vec_p), intent(inout) :: v1
    type(type_vec_p), intent(inout) :: v2
    call bp_operate(wf_save, v2, j_square_save, v1)
  end subroutine matvec_jj


  subroutine matvec_block(nb, v1, v2)
    integer, intent(in) :: nb
    real(kwf), intent(in)  :: v1(:,:)
    real(kwf), intent(out) :: v2(:,:)
    call bp_operate_block(wf_save, v2(:,:nb), &
         hamltn_save, v1(:,:nb))
  end subroutine matvec_block

  subroutine matvec_block_jj(nb, v1, v2)
    integer, intent(in) :: nb
    real(kwf), intent(in)  :: v1(:,:)
    real(kwf), intent(out) :: v2(:,:)
    call bp_operate_block(wf_save, v2(:,:nb), &
         j_square_save, v1(:,:nb))
  end subroutine matvec_block_jj

end module kshell_func




program kshell
#ifdef MPI
  use mpi
#endif
  use constant, only: pi, kwf, kdim, kmbit, nmbit, maxchar, c_no_init, max_int4, max_n_jorb
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, is_mpi, ierr, n_jorb_pn, n_jorb, sum_rule, &
       nprocs_reduce, nprocs_shift, n_ferm, n_core, is_debug
  use model_space, only: m_mass=>mass, allocate_l_vec, deallocate_l_vec, &
       print_max_l_vec, nv_shift, corb, jorb
  use class_stopwatch
  use interaction, only: read_interaction, hamltn, ham_cm, j_square, t_square, &
       r2y2, jtensor, ltensor, stensor, r1y1, set_gt, gt_m, &
       set_r1y1t, set_sd0t, set_sd1t, set_sd2t, &
       r1y1t_m, sd0t_m, sd1t_m, sd2t_m, &
       set_ry_sum, sum_rank, rkyk_square, r1y1_f_square, set_stst0, &
       set_three_body_monopole, set_ob_ij, r3y3, r3y1
  use operator_mscheme, only: opr_m, opr_m_p, add_opr_m, opr_m_eff_charge, &
       opr_m_one_crt, opr_m_one_anh
  use partition, only: type_ptn_pn, init_partition, &
       finalize_partition, deploy_partition, &
       cost_from_localdim
  use wavefunction, only: type_vec_p, ex_occ_orb, wf_alloc_vec, &
       wf_random_vec, ratio_nocc, hw_ratio_nocc, ph_ratio_nocc, &
       inf_entropy, inf_entropy_pn, ratio_nocc_orbs
  use bridge_partitions, only: type_bridge_partitions, &
       init_bridge_partitions, finalize_bridge_partitions, &
       init_bp_operator, finalize_bp_operator, bp_operate, &
       ex_val, eig_residual, init_mpi_shift_reduce
  use bp_expc_val, only: bp_ex_vals_pn 
  use bp_io, only: bp_save_wf, bp_load_wf
  use lanczos, only: lanczos_main, max_lanc_vec_doublej, set_lanczos_tmp_fn
  use rotation_group, only: dcg
  use lib_matrix, only: set_seed, set_rank_seed, gaussian_random_mat
  use kshell_func, only: set_kshell_func, matvec, dotprod, matvec_jj, &
       matvec_block, matvec_block_jj
  use block_lanczos, only: tr_block_lanczos, maxiter_jjrefine
#ifdef ZPARES
  use ss_method, only: zpares_diag, zpares_level_dens
  use eigdensity_estimator_mod, only : projection_zpares
#endif
  !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num
  implicit none
  type(type_ptn_pn), target :: ptn, ptn_init, ptn_srt ! partition information
  type(type_bridge_partitions) :: wf
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn_int, fn_nml, fn_ptn, fn_ptn_init, &
       fn_save_wave, fn_load_wave, op_type_init, ctmp, fn_three_body_mp
  character :: cp = '?'
  integer :: mtot, hw_type, n_eigen, &
       n_restart_vec, max_lanc_vec, maxiter, Jguess, mtot_init, n_eig_init, &
       mass, mode_lv_hdd, neig_load_wave, tt_proj
  real(8) :: beta_cm, eff_charge(2), gl(2), gs(2), e1_charge(2), tol, add_randinit=0.d0
  logical :: is_double_j, is_load_snapshot, is_calc_tbdm
  !
  type(type_vec_p), allocatable :: evec(:)
  real(kwf), allocatable :: bl_evec(:,:)
  real(8), allocatable :: eval(:), cost(:), occ(:), evs(:,:)
  real(8), allocatable :: estimations(:)
  type(opr_m) :: op_init_wf, op_init_s
  type(opr_m_p), allocatable :: ops(:)
  real(8) :: x, c, hcm
  integer :: i, j, n, provided, required, irank
  integer :: access
  !
  real(8), parameter :: ss_e_range_init = -1.d8
  real(8) :: ss_e_range(2) = ss_e_range_init, skip_ld = 0.d0
  integer :: n_block = -1, n_ld_reso = -1, nn, comm=0
  integer(kdim) :: mq
  logical :: is_leveldens = .false.,  is_ss_diag = .false., &
       is_bl_lan = .false., is_lanczos = .true.
  integer :: orbs_ratio(max_n_jorb)
  logical :: is_h2 = .false.
  logical :: is_matvecone_dump = .false. ! save v' = H*v  v=(1,1,1,1,...,1)
  type(type_vec_p) :: vone
  !
  namelist /input/ fn_int, fn_ptn, fn_ptn_init, &
       mtot, hw_type, n_eigen, &
       n_restart_vec, max_lanc_vec, maxiter, is_double_j, &
       fn_save_wave, fn_load_wave, is_load_snapshot, &
       beta_cm, eff_charge, gl, gs, e1_charge, op_type_init, mass, &
       mode_lv_hdd, is_calc_tbdm, tol, neig_load_wave, &
       ss_e_range, n_block, n_ld_reso, skip_ld, &
       nv_shift, nprocs_reduce, tt_proj, add_randinit, orbs_ratio, &
       is_matvecone_dump, fn_three_body_mp, is_h2
  

#ifdef MPI
  is_mpi = .true.
#if defined (_OPENMP) && defined(SPARC) 
  required = mpi_thread_serialized
  call mpi_init_thread(required, provided, ierr)
  if (provided < required) write(*,*) "***** warning in mpi_init_thread *****"
#else  
  call mpi_init(ierr)
#endif

  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  ! write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "MPI : # of processes=",nprocs,"    myrank", myrank
  if (myrank/=0) is_debug = .false.
#endif
  !$ if (myrank==0) write(*,'(1a,1i3,/)') "OpenMP :  # of threads=", omp_get_max_threads()

  call start_stopwatch(time_total, is_reset=.true.)
  call start_stopwatch(time_preproc, is_reset=.true.)
  
  call set_seed()
!  call set_seed(is_clock=.true.)
  call set_rank_seed(myrank)


  ! default parameters -------------
  mass = 0             ! mass number, optional
  mtot = 0             ! Jz * 2
  n_eigen = 1          ! # of eigenvalues to be otained
  n_restart_vec = 10   ! # of Lanczos vectors for thick-restart Lanczos
  max_lanc_vec = 100   ! max. # of vectors for thick-restart Lanczos
  maxiter = 300        ! max. # of iteration for thick-restart Lanczos
  hw_type = 1          ! harmonic oscillator formula
  is_double_j = .false. ! double Lanczos for J-projection
  is_load_snapshot = .false. ! snapshot restart at Thick-restart dump files
  fn_save_wave = c_no_init  ! file name of save wave functions 
  fn_ptn_init  = c_no_init  ! partion file for loading w.f. (def. fn_ptn)
  fn_load_wave = c_no_init  ! file name of load wave functions 
  neig_load_wave = 1   ! n-th wave function at "fn_load_wave" for initial (-1 : all states)
  beta_cm = 0.d0        ! Lawson parameter beta_cm (= beta*hw/A like OXBASH)
  eff_charge = (/ 1.d0, 0.d0 /) ! effective charges for E2, Q-moment
  e1_charge  = (/ 0.d0, 0.d0 /) ! effective charges for E1 
  gl = (/1.d0, 0.d0/)  ! gyromagnetic ratios for orbital angular momentum
  gs = (/5.586d0, -3.826d0/) ! gyromagnetic ratios for spin
  op_type_init = c_no_init ! operate init w.f. E2, E1, M1, GT for strength function
  tol = 1.d-6      ! convergence condition for Lanczos method
  ! if (kwf==8) tol = 1.d-7
  mode_lv_hdd = 0  ! see lanczos_main for the description
  is_calc_tbdm = .false. ! two-body density matrix
  nprocs_reduce = 1      ! # of process for MPI reduction
  nv_shift = 0           ! # of vectors for shift 
  tt_proj = -1           ! isospin projection for LSF 2*T
  orbs_ratio(:) = 0      ! orbit number for showing ratio of w.f.
  fn_three_body_mp = c_no_init  ! file name for three-body monopole intearction
  ! -----------------------------------

  call getarg(1, fn_nml)
#ifdef MPI
  call mpi_bcast(fn_nml,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
#endif
  open(lunnml, file=fn_nml, status='old')
  read(lunnml, nml=input)  
  close(lunnml)
#ifdef MPI
  call init_mpi_shift_reduce()
  comm = mpi_comm_world
#endif

  call print_mem_status('01-init')

  if (n_eigen > n_restart_vec) n_restart_vec = n_eigen
  if (n_eigen >= max_lanc_vec) max_lanc_vec = n_eigen + 1 
  if (max_lanc_vec <= n_restart_vec) stop "max_lanc_vec should be larger than n_restart_vec"
  if (n_block > 0) mode_lv_hdd = 0
  if (nv_shift == 0) nv_shift = 1

  if (myrank==0) write(*,nml=input)
  if (myrank==0) write(*,'(a,4i3)') &
       "compile conf. kwf, kdim, kmbit, nmbit =", kwf, kdim, kmbit, nmbit

  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)

  call print_mem_status('015-readsps')
  
  open(lunptn, file=fn_ptn, status='old')
  if (myrank==0) write(*,'(1a,1i3,2a)') "set partition Mtotal=", mtot, &
       "  partition_file= ", trim(fn_ptn)
  call init_partition(ptn, lunptn, mtot)
  close(lunptn)

  call print_mem_status('02-part')

  if (ptn%iprty ==  1) cp = '+'
  if (ptn%iprty == -1) cp = '-'
  if (myrank==0) write(*,'(a,i3,2a,/)') &
       'M = ', ptn%mtotal, '/2  :  parity = ', cp

  call set_n_ferm(ptn%n_ferm(1), ptn%n_ferm(2), mass)
  call read_interaction(lunint, hw_type=hw_type)
  close(lunint)

  if (fn_three_body_mp /= c_no_init) &
       call set_three_body_monopole(fn_three_body_mp, hamltn)


  call print_mem_status('03-hamil')

  if (e1_charge(1) == 0.d0 .and. e1_charge(2) == 0.d0) then
     e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(m_mass)
     e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(m_mass)
  end if

  if (beta_cm/=0.d0 .and. (n_core(1)<0 .or. n_core(2)<0)) &
       stop 'NOT implemented Hcm with hole state'
  if (beta_cm /= 0.d0) call add_opr_m( hamltn, beta_cm, ham_cm)

  call print_mem_status('04-add_opr_m')

  allocate( cost(ptn%n_pidpnM) )
  call cost_from_localdim(ptn, ptn%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptn, cost, verbose=.true.)
  deallocate(cost)

  if (ptn%ndim<=0) stop "ERROR: 0 dimension"

  call print_mem_status('05-init_ptn')

  n = n_eigen
  if (n_block > 0) n = n_block
  allocate( eval(n), evec(n))
  eval = 0.d0
  do i = 1, size(evec)
     evec(i)%ptn => ptn
  end do


  if (myrank==0) then
     x = ptn%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for one global Lanczos vector:", x, " GB"
     x = ptn%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     n = max_lanc_vec
     if (mode_lv_hdd == 2) n = 2
     if (mode_lv_hdd == 1) n = max(n_eigen, 2)
     if (n_block > 0) then
        n = n + n_block * 2
        if (is_double_j) n = n + n_block * maxiter_jjrefine
     else
        if (is_double_j) n = n + max_lanc_vec_doublej
     end if
#ifdef MPI
     if (n_block <= 0) then 
        n = n + nprocs_reduce - 1 + nv_shift - 1
        n = n + 1 ! vltmp
     else ! adhoc
        n = n + (nprocs_reduce + nv_shift - 2)*n_block
     end if
#endif
     write(*,'(a,f10.3,a,i6,a,f10.3,a)') &
          "Memory / process is:", x, " GB x ", n, " = ", x*n, " GB"
     write(*,'(a,f10.3,a,/)') &
          "Total Memory for Lanczos vectors:", x*n*nprocs, " GB"
  end if


  if (n_block > 0) then
     is_lanczos = .false.
     if ( abs(ss_e_range(1) - ss_e_range_init) > 1.d-8 ) then
        if ( n_ld_reso <= 0 ) then
           is_ss_diag = .true.
           if (myrank==0) write(*,*) " SS-method solver "
        else
           is_leveldens = .true.
           if (myrank==0) write(*,*) " Stochastic estimation of Level density "
        end if
     else
        is_bl_lan = .true.
        if (myrank==0) write(*,*) " Block Lanczos method "
     end if

     if ( ptn%local_dim > max_int4 ) then
        write(*,*) "ERROR: local_dim > max_int4"
        goto 999
     end if

     allocate( bl_evec( ptn%local_dim, n_block ) )
  end if



  call print_mem_status('06-bl')


  ! sorted partition for load and save
  if ((.not. is_load_snapshot) .and. fn_load_wave /= c_no_init) then

     if ( access(fn_load_wave, 'r') /= 0 ) stop "NOT found fn_load_wave"
    
     ! read header of wave functions
     open(lunwv, file=fn_load_wave, form='unformatted', &
          status='old', access='stream')
     read(lunwv) n_eig_init
     read(lunwv) mtot_init
     close(lunwv)

     if (neig_load_wave > n_eig_init) stop 'ERROR: neig_load_wave'

     if ( op_type_init == c_no_init .or. op_type_init == "copy" ) then
        if (mtot_init /= mtot) stop 'ERROR mtot_init /= mtot, use mup_operate.exe '
        op_init_wf%nbody = 0
        op_init_wf%irank = 0

     else if (op_type_init == "E1" .or. op_type_init == "e1") then

        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(E1)|M=", mtot_init, "/2>"
           write(*,'(1a,2f9.5,/)') "effective charge for E1 ", e1_charge
        end if
        call opr_m_eff_charge(op_init_wf, r1y1, e1_charge)

     else if (op_type_init == "E2" .or. op_type_init == "e2") then

        if (myrank==0) then 
           write(*,'(/,a,i3,a)') "initial vec = T(E2)|M=", &
                mtot_init, "/2>"
           write(*,'(1a,2f9.5,/)') "effective charge", eff_charge
        end if
        call opr_m_eff_charge(op_init_wf, r2y2, eff_charge)

     else if (op_type_init == "E3" .or. op_type_init == "e3") then

        if (myrank==0) then 
           write(*,'(/,a,i3,a)') "initial vec = T(E3)|M=", &
                mtot_init, "/2>"
           write(*,'(1a,2f9.5,/)') "effective charge", eff_charge
        end if
        call opr_m_eff_charge(op_init_wf, r3y3, eff_charge)

     else if (op_type_init == "M1" .or. op_type_init == "m1") then

        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(M1)|M=",mtot_init,"/2>"
           write(*,'(1a,2f9.5)') "gl = ", gl
           write(*,'(1a,2f9.5,/)') "gs = ", gs
        end if
        call opr_m_eff_charge(op_init_wf, ltensor, gl / sqrt(4.d0*pi/3.d0))
        call opr_m_eff_charge(op_init_s,  stensor, gs / sqrt(4.d0*pi/3.d0))
        call add_opr_m(op_init_wf, 1.d0, op_init_s)

     else if (op_type_init == "ISD" .or. op_type_init == "isd") then
        ! isoscalar dipole r3Y1
        if (myrank==0) then 
           write(*,'(/,a,i3,a)') "initial vec = T(ISD)|M=", &
                mtot_init, "/2>"
           write(*,'(1a,2f9.5,/)') "effective charge", 1.d0, 1.d0 
        end if
        call opr_m_eff_charge(op_init_wf, r3y1, (/ 1.d0, 1.d0 /) )

     else if (op_type_init == "GT" .or. op_type_init == "gt") then

        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(GT)|M=",mtot_init,"/2>"
        end if
        call set_gt()
        call opr_m_eff_charge(op_init_wf, gt_m, (/1.d0, 1.d0/))

     else if (op_type_init == "r1Y1t" .or. op_type_init == "r1y1t") then

        if (myrank==0) write(*,'(/,1a,1i3,1a)') &
             "initial vec = T(r1Y1t)|M=",mtot_init,"/2>"
        call set_r1y1t()
        call opr_m_eff_charge(op_init_wf, r1y1t_m, (/1.d0, 1.d0/))

     else if (op_type_init == "SD0t" .or. op_type_init == "sd0t") then

        if (myrank==0) write(*,'(/,1a,1i3,1a)') &
             "initial vec = T(SD0t)|M=",mtot_init,"/2>"
        call set_sd0t()
        call opr_m_eff_charge(op_init_wf, sd0t_m, (/1.d0, 1.d0/))

     else if (op_type_init == "SD1t" .or. op_type_init == "sd1t") then

        if (myrank==0) write(*,'(/,1a,1i3,1a)') &
             "initial vec = T(SD1t)|M=",mtot_init,"/2>"
        call set_sd1t()
        call opr_m_eff_charge(op_init_wf, sd1t_m, (/1.d0, 1.d0/))

     else if (op_type_init == "SD2t" .or. op_type_init == "sd2t") then

        if (myrank==0) write(*,'(/,1a,1i3,1a)') &
             "initial vec = T(SD2t)|M=",mtot_init,"/2>"
        call set_sd2t()
        call opr_m_eff_charge(op_init_wf, sd2t_m, (/1.d0, 1.d0/))

     else if (op_type_init == "DGT0" .or. op_type_init == "dgt0") then
        if (myrank==0) write(*,'(/,1a,1i3,1a)') &
             "initial vec = T( [st*st]^(0) )|M=",mtot_init,"/2>"
        call set_stst0(op_init_wf)

     else if (op_type_init(:5) == 'gtorb') then
        !
        ! dump w.f.  1/N * [ c^i c_j ]^(irank)| Jinit>
        ! e.g. gtorb_03_08_02  for SNV [ c^3 c_8 ]^(2)
        !
        read(op_type_init( 7:8 ), *) i
        read(op_type_init(10:11), *) j
        read(op_type_init(13:14), *) irank
        if (myrank==0) write(*,'(/,a,i3,a,i3,a,i3,a,i3,a)') &
             'initial vec = T( [',i, '*', j, ']^(',irank,') )|M=', &
             mtot_init,'/2>'
        if (myrank==0) write(*,'(/,a,i3)') 'neig_load_wave = ', neig_load_wave 

        call set_ob_ij(i, j, irank, op_init_wf)
        call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
             fn_ptn_init, mtot_init, op_init_wf, op_type_init, &
             neig_load_wave, tt_proj)

        call bp_save_wf(fn_save_wave, evec(:n_eigen), ptn)

#ifdef MPI
        call mpi_finalize(ierr)
#endif
        stop 

     else if (op_type_init(:3) == 'sp_') then
        !
        ! dump w.f.  1/N *  c^+_i | Jinit>
        ! e.g. sp_03   c^+_3 |Jinit>
        !
        read(op_type_init( 4:5 ), *) i
        if (myrank==0) write(*,'(/,a,i3,a,i3,a,i3,a,i3,a)') &
             'initial vec = c+_',i,'  ( j =',jorb(i),'/2) )|M=', &
             mtot_init,'/2>'
        if (myrank==0) write(*,'(/,a,i3)') 'neig_load_wave = ', neig_load_wave 

        call opr_m_one_crt(op_init_wf, i, mtot - mtot_init)

        call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
             fn_ptn_init, mtot_init, op_init_wf, op_type_init, &
             neig_load_wave, tt_proj)

        if (maxiter==0) then
           call bp_save_wf(fn_save_wave, evec(:n_eigen), ptn)
#ifdef MPI
           call mpi_finalize(ierr)
#endif
           stop
        end if

     else if (op_type_init(:3) == 'sd_') then
        !
        ! one-particle annihilation
        ! dump w.f.  1/N *  c_i | Jinit>  
        ! e.g. sd_03   c_3 |Jinit>
        !
        read(op_type_init( 4:5 ), *) i
        if (myrank==0) write(*,'(/,a,i3,a,i3,a,i3,a,i3,a)') &
             'initial vec = c_',i,'  ( j =',jorb(i),'/2) )|M=', &
             mtot_init,'/2>'
        if (myrank==0) write(*,'(/,a,i3)') 'neig_load_wave = ', neig_load_wave 

        call opr_m_one_anh(op_init_wf, i, mtot - mtot_init)
        ! call opr_m_one_anh(op_init_wf, i,  mtot_init - mtot)

        call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
             fn_ptn_init, mtot_init, op_init_wf, op_type_init, &
             neig_load_wave, tt_proj)

        if (maxiter==0) then
           call bp_save_wf(fn_save_wave, evec(:n_eigen), ptn)
#ifdef MPI
           call mpi_finalize(ierr)
#endif
           stop
        end if
        
     else

        
        stop "not implemented op_type_init"

     end if

     if (fn_ptn_init == c_no_init) fn_ptn_init = fn_ptn

     if ( neig_load_wave <= -1 .and. n_eig_init > size(evec) ) then
        deallocate( eval, evec ) 
        allocate( eval(n_eig_init), evec(n_eig_init) )
        eval = 0.d0
        do i = 1, size(evec)
           evec(i)%ptn => ptn
        end do
     end if

     call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
          fn_ptn_init, mtot_init, op_init_wf, op_type_init, &
          neig_load_wave, tt_proj )
     
     if (neig_load_wave == -1) call compress_init_vecs(evec, max(1, n_block))

     if (myrank==0) write(*,*) ' load w.f. finished'
     
     if (add_randinit /= 0.d0) &
          call add_random_init_vecs(evec, max(1, n_block), add_randinit)
     
  end if
  
  call print_mem_status('07-load')

  
  call init_bridge_partitions(wf, ptn, verbose=.true.)
  call print_mem_status('071-init-bp')
  call init_bp_operator(wf, j_square)
  call print_mem_status('072-j-sq')
  call init_bp_operator(wf, hamltn, verbose=.true.)
  call print_mem_status('073-hamltn')
  call set_kshell_func(wf, hamltn, j_square)
  call print_mem_status('074-set-j-sq')

  call stop_stopwatch(time_preproc)


  if (is_matvecone_dump) then
     if (n_eigen /= 1) stop 'set n_eigen = 1'
     if (myrank==0) write(*,*) 'matvecone_dump called'
     vone%ptn => ptn
     call wf_alloc_vec( vone, ptn )
     call wf_alloc_vec( evec(1), ptn )
     vone%p(:) = 0._kwf
     vone%p(:ptn%local_dim) = 1._kwf

     call matvec(vone, evec(1))

     if (myrank==0) write(*,*) 'matvec passed'

     goto 999
  end if



!  if (.true.) then
  if (.false.) then
     call ex_val_snt()
     goto 999
  end if


  if (is_h2) then  ! energy variance 
     call ex_val_h2()
     goto 999
  end if
     
  
  call print_mem_status('08-init-bp')


  if (is_leveldens) then

     do j = 1, size(evec)
        if (associated(evec(j)%p)) deallocate( evec(j)%p )
     end do
     do j = 1, n_block
        call gaussian_random_mat( ptn%local_dim, bl_evec(:,j) )
        !$omp parallel do
        do mq = 1, ptn%local_dim
           bl_evec(mq, j) = sign(1._kwf, bl_evec(mq, j))
        end do
     end do

  elseif (is_ss_diag .or. is_bl_lan) then

     do j = 1, n_block
        if ( associated(evec(j)%p) ) then
           bl_evec(:,j) = evec(j)%p
        else
           call gaussian_random_mat( ptn%local_dim, bl_evec(:,j) )
        end if
     end do
     do j = 1, size(evec)
        if (associated(evec(j)%p)) deallocate( evec(j)%p )
     end do

  end if
     

  call print_mem_status('09-initwf-bl')


  ! J-projection for SS-method and level density estimation
  if ( ( is_ss_diag .or. is_leveldens) .and. is_double_j) then

     x = mtot*(mtot+2)*0.25d0
     if (myrank==0) write(*,'(/,a,f8.4,i3/)') &
          "J projection for initial state ", x, n_block

     nn = 4 ! max block size to compute at once for memory saving
     ! nn = 1
     ! nn = 6
     if (nn > n_block) nn = n_block

     do i = 1, n_block, nn
        n = min(nn, n_block-i+1)

        if (myrank==0) write(*,'(/,a,i3,/)') &
             ' J-projection with block size ', n

#ifdef ZPARES
        call projection_zpares(comm, int(ptn%local_dim), &
             n, matvec_block_jj, &
             x-0.5d0, x+0.5d0, maxiter, bl_evec(:,i:i+n-1))
#else
        stop 'enable ZPARES'
#endif /* ZPARES */

        ! call zpares_diag(matvec_block_jj, &
        !      x-0.5d0, x+0.5d0, int(ptn%local_dim), &
        !      n, n_eigen, teval, &
        !      tevec, maxiter, is_proj=.true.)

     end do
     ! write(*,'(a,100f10.4)')"projected norm",i,sum(abs(bl_evec(:,i)))

  end if


  if (is_leveldens) then
     ! stochastic estimation of level density

     if ( ss_e_range(2) <= ss_e_range(1)+1.d-8 .and. skip_ld==0.d0 ) &
          skip_ld = 1.d0
     if ( skip_ld > 0.d0 ) &
          ss_e_range(2) = ss_e_range(1) + dble(n_ld_reso) * skip_ld
     allocate( estimations(n_ld_reso) )
#ifdef ZPARES
     call zpares_level_dens( matvec_block, &
          ss_e_range(1), ss_e_range(2), int(ptn%local_dim), n_block, &
          estimations, bl_evec, maxiter )
#else
     stop 'enable ZPARES'
#endif /* ZPARES */
     goto 999

  end if


  write(ctmp,'(a,"_",i0)') trim(fn_ptn),mtot
  call set_lanczos_tmp_fn(ctmp, fn_save_wave/=c_no_init)

  if ( is_ss_diag ) then
#ifdef ZPARES
     call zpares_diag(matvec_block, &
          ss_e_range(1), ss_e_range(2), int(ptn%local_dim), &
          n_block, n_eigen, eval, bl_evec, maxiter)

     if (n_eigen <= 0) then
        if (myrank == 0) write(*,*) 'eigenvalue NOT found'
        goto 999
     end if
#else
     stop 'enable ZPARES'
#endif /* ZPARES */
  end if

  if ( is_bl_lan ) then

     if (kwf==4) stop 'Not implemented kwf=4 and block Lanczos'

     if (is_double_j) then
        call tr_block_lanczos( matvec_block, matvec_block_jj, &
             ptn%local_dim, n_block, n_eigen, eval, bl_evec, &
             max_lanc_vec=max_lanc_vec, maxiter=maxiter, &
             tol=tol, n_res_vec=n_restart_vec, &
             is_load_snapshot=is_load_snapshot, &
             eval_jj=mtot*(mtot+2)*0.25d0 )
     else
        call tr_block_lanczos( matvec_block, matvec_block_jj, &
             ptn%local_dim, n_block, n_eigen, eval, bl_evec, &
             max_lanc_vec=max_lanc_vec, maxiter=maxiter, &
             tol=tol, n_res_vec=n_restart_vec, &
             is_load_snapshot=is_load_snapshot )
     end if

  end if

  if ( is_ss_diag .or. is_bl_lan ) then
     if ( allocated(evec) ) then 
        do i = 1, size(evec) 
           if (associated(evec(i)%p)) deallocate(evec(i)%p)
        end do
        deallocate(evec)
     end if
     allocate( evec(n_eigen)  )
     do i = 1, size(evec)
        evec(i)%ptn => ptn
        call wf_alloc_vec( evec(i), ptn )
        !$omp parallel do
        do mq = 1, ptn%local_dim
           evec(i)%p(mq) = bl_evec(mq, i)
        end do
        do mq = ptn%local_dim+1, ptn%max_local_dim
           evec(i)%p(mq) = 0._kwf
        end do
     end do
     
     deallocate( bl_evec )
  end if


  if (is_lanczos) then
     if (.not. associated(evec(1)%p) .and. .not. is_load_snapshot) &
          call wf_random(ptn, evec(1), mtot, is_double_j)

     call print_mem_status('10-initwf')

     if (maxiter > 0) then

        ! if (.not. is_load_snapshot) call compress_init_vecs(evec, 1)

        if (myrank==0) write(*,*)
        if (myrank==0) write(*,*) '*** Lanczos start ***'
        if (myrank==0) write(*,*)


        if (is_double_j) then ! double lanczos
           call lanczos_main(matvec, dotprod, ptn%max_local_dim, eval, evec, &
                matvec_jj, eval_jj=dble(mtot*(mtot+2)*0.25d0), n_eig=n_eigen, &
                is_load_snapshot=is_load_snapshot, &
                n_restart_vec=n_restart_vec, max_lanc_vec=max_lanc_vec, &
                maxiter=maxiter, mode_lv_hdd=mode_lv_hdd, tol=tol)
        else ! diag w/o jj
           call lanczos_main(matvec, dotprod, ptn%max_local_dim, eval, evec, &
                matvec_jj, is_load_snapshot=is_load_snapshot, n_eig=n_eigen, &
                n_restart_vec=n_restart_vec, max_lanc_vec=max_lanc_vec, &
                maxiter=maxiter, mode_lv_hdd=mode_lv_hdd, tol=tol)
        end if

        do i = 1, n_eigen
           if (.not. associated(evec(i)%p)) then 
              n_eigen = i-1
              exit
           end if
           if (myrank==0) write(*,'(a,i6,f12.6)') "lanczos eigenvalues",i, eval(i)
        end do

     end if

  end if
     
  call print_mem_status('11-lanczos')

  if (myrank==0) write(*, '(/, a, f10.3, a, f10.3, a, /)') &
       "total time it took was:", &
       get_ctime_stopwatch(time_total) ," sec. ", &
       get_ctime_stopwatch(time_total) / 3600.d0, " hours"
  call stop_stopwatch(time_total)
  call print_summary_stopwatch()
  call start_stopwatch(time_total)

  if ( fn_save_wave == c_no_init .and. &
       op_type_init /= c_no_init .and. op_type_init /= "copy" ) &
       goto 999 


  do i = 1, n_eigen
     if ( associated(evec(i)%p) ) cycle
     n_eigen = i - 1
     exit
  end do
  

  if (neig_load_wave == 0 .and. maxiter == 0) then
     do i = 1, n_eigen
        call eig_residual( wf, eval(i), evec(i), hamltn, x )
        if (myrank==0) write(*,'(a,i8,f12.6,ES12.3)') "eigval,res=", i, eval(i), x
     end do
  end if


  do i = 1, n_eigen
     if (maxiter > 0) then
        evec(i)%eval = eval(i)
     else ! just read wave function
        call ex_val(wf, evec(i), hamltn, evec(i)%eval)
     end if
  end do

  call print_mem_status('12-postlan')

  call finalize_bp_operator(wf, hamltn)

  call print_mem_status('13-fin-hamltn')


  allocate( occ(n_jorb_pn) )

  if (myrank==0) then
     write(*,*)
     write(*,'(a,2f7.3)') ' effective charges ',eff_charge
     write(*,'(a,4f8.4)')  ' gl,gs = ', gl, gs
     if (maxval(orbs_ratio) > 0) then
        write(*,'(a)', advance='no') ' orbits for ratio : '
        do i = 1, size(orbs_ratio)
           if ( orbs_ratio(i) == 0 ) cycle
           write(*,'(i3,4a)',advance='no') orbs_ratio(i), ':', corb(orbs_ratio(i)), ','
        end do
        write(*,*)
     end if
     write(*,'(" p orbit ", 50a)')  (/( corb(i)(3:) // " ", i= 1, n_jorb(1) )/)
     write(*,'(" n orbit ", 50a)')  (/( corb(i)(3:) // " ", i= n_jorb(1)+1, n_jorb_pn )/)
  end if

  call init_bp_operator(wf, t_square)
  call print_mem_status('14-init-bp-tt')
  call init_bp_operator(wf, r2y2)
  call print_mem_status('142-init-bp-r2y2')
  call init_bp_operator(wf, ltensor)
  call print_mem_status('15-init-bp-ltensor')
!  call init_bp_operator(wf, jtensor)
  call init_bp_operator(wf, stensor)
  call print_mem_status('16-init-bp-st')
  if (beta_cm/=0.d0) call init_bp_operator(wf, ham_cm)
  call print_mem_status('17-init-bp-cm')
  n = 3
  allocate( ops(n), evs(2,n) )

  ! E1 sum rule
!  call set_ry_sum(1, e1_charge)
!  call init_bp_operator(wf, rkyk_square)
!  call init_bp_operator(wf, r1y1_f_square)
  ! E2 sum rule
  ! call set_ry_sum(2, eff_charge) 
  ! call init_bp_operator(wf, rkyk_square)


  if (myrank==0) write(*,'(a)') "-------------------------------------------------"

  do i = 1, n_eigen
     call ex_val(wf, evec(i), j_square, x)
     Jguess = guess_J_from_JJ(ptn, x)
     evec(i)%jj = Jguess
     if (beta_cm/=0.d0) then 
        call ex_val(wf, evec(i), ham_cm, hcm)
        evec(i)%eval = evec(i)%eval - beta_cm * hcm
     end if
     if (myrank==0) then 
        write(*,'(i4,a,1f12.5,a,f12.5,a,i3, a, i2)') &
             i, '  <H>:',evec(i)%eval, '  <JJ>:',x, '  J:',Jguess, &
             "/2  prty ", ptn%iprty
     end if

     call ex_val(wf, evec(i), t_square, x)
     x = x + t_square%e0 
     evec(i)%tt = guess_J_from_JJ(ptn, x)
     
     if (myrank==0) then
        if (beta_cm/=0.d0) then 
           write(*,'(4x,a,1f12.5,a,1f12.5,a,i3,a)') '<Hcm>:', hcm, &
                '  <TT>:', x, '  T:', evec(i)%tt, "/2"
        else
           write(*,'(23x,a,1f12.5,a,i3,a)') ' <TT>:', x, &
                '  T:', evec(i)%tt, "/2"
        end if
     end if

     call ex_occ_orb(evec(i), occ)
     if (myrank==0) write(*,'(" <p Nj>", 50f7.3)') occ(:n_jorb(1))
     if (myrank==0) write(*,'(" <n Nj>", 50f7.3)') occ(n_jorb(1)+1:)

     if (Jguess > 0) then
! call ex_val(wf, evec(i), r2y2, x)
! if (myrank==0) write(*,'(1a, 100f7.3)') " < Q > ", x*dsqrt(16.0d0*pi/5.0d0)

        ops(1)%p => r2y2
        ops(2)%p => ltensor
        ops(3)%p => stensor
        call bp_ex_vals_pn(wf, evec(i), ops, evs)
        evs(:,1) = evs(:,1) * dsqrt(16.0d0*pi/5.0d0)
        if (myrank==0) then 
           x = dot_product(evs(:,1), eff_charge)
           c = dcg(Jguess, mtot, 4, 0, Jguess, mtot)
           if (abs(c)> 1.d-8) then
              c = dcg(Jguess, Jguess, 4, 0, Jguess, Jguess) / c
              write(*,'(1a, 1f9.3, 1a, 1f9.3, 1a, 1f9.3)') &
                   "   <Qp> ", evs(1,1)*c, "   <Qn> ", evs(2,1)*c, &
                   "   <eQ> ", x*c
           end if
           c = dcg(Jguess, mtot, 2, 0, Jguess, mtot)
           if (abs(c)> 1.d-8) then 
              c = dcg(Jguess, Jguess, 2, 0, Jguess, Jguess) / c
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <Lp> ", evs(1,2)*c, "   <Ln> ", evs(2,2)*c
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <sp> ", evs(1,3)*c, "   <sn> ", evs(2,3)*c
              x = (dot_product(evs(:,2),gl) + dot_product(evs(:,3),gs))
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <gm> ", x*c, "   <Jz> ", sum(evs(:,2)) + sum(evs(:,3))
           end if
        end if
     end if

     ! ratio of occupation
     if (maxval(orbs_ratio)==0) then
        call hw_ratio_nocc(evec(i))
     else
        call ratio_nocc_orbs(evec(i), orbs_ratio)
     end if
     ! call ratio_nocc(evec(i))
     ! call ph_ratio_nocc(evec(i), 1, min(n_jorb(1)+1, n_jorb_pn))

     ! Information entropy 
     ! x = inf_entropy( evec(i) )
     ! if (myrank==0) write(*,'(a,f12.5)')  'Information entropy   ', x
     ! call inf_entropy_pn( evec(i), r2 )
     ! if (myrank==0) write(*,'(a,2f12.5)') 'Information entropy pn', r2
     
    
     if ( sum_rank == 1 ) then
        if ( evec(i)%jj == 0 ) then
           if (myrank==0) write(*,'(a,2f8.3)') &
                'B(E1) sum rule : Eff. charge ',e1_charge
           call ex_val(wf, evec(i), rkyk_square, x)
           if ( myrank == 0) write(*,'(a, f12.5)')  "<rY*P*rY>    ", x
           call ex_val(wf, evec(i), r1y1_f_square, x)
           x = x + r1y1_f_square%e0
           if (myrank==0) write(*,'(a, f12.5)')  "<rY * rY>_fs ", x
        end if
     elseif ( sum_rank == 2 ) then
        if (myrank == 0) write(*,'(a,2f8.3)') &
             'B(E2) sum rule : Eff. charge ', eff_charge
        call ex_val(wf, evec(i), rkyk_square, x)
        if (myrank==0) write(*,'(a, f12.5)')  "<r2Y2*P*r2Y2>", x
     end if

        
     if (myrank==0) write(*,'(1a)') &
          "-------------------------------------------------"
  end do

  call print_mem_status('18-summary')

  call finalize_bp_operator(wf, rkyk_square)
  call finalize_bp_operator(wf, j_square)
  call finalize_bp_operator(wf, t_square)
  call finalize_bp_operator(wf, r2y2)
  call finalize_bp_operator(wf, ltensor)
!  call finalize_bp_operator(wf, jtensor)
  call finalize_bp_operator(wf, stensor)

  call print_mem_status('19-fin-bp-op')

999 continue

  call finalize_bp_operator(wf, hamltn)

  if (beta_cm /= 0.d0) call finalize_bp_operator(wf, ham_cm)

  call finalize_bridge_partitions(wf)

  call print_mem_status('20-fin-bps')

  if ( associated(evec(1)%p) .and. fn_save_wave /= c_no_init) &
       call bp_save_wf(fn_save_wave, evec(:n_eigen), ptn)

  call print_mem_status('205-save-wf')

  
  if (is_calc_tbdm) then 
     call init_bridge_partitions(wf, ptn)
     call calc_tbme_from_tbtd(wf, evec)
     call finalize_bridge_partitions(wf)
  end if


  call print_max_l_vec()

  call stop_stopwatch(time_total)
  call print_summary_stopwatch()

  call print_mem_status('21-final')

#ifdef MPI
  call mpi_finalize(ierr)
#endif
  
contains


  function guess_J_from_JJ(ptn, x) result(jj)
    type(type_ptn_pn), intent(in) :: ptn
    real(8), intent(in) :: x 
    integer :: jj, i, imin
    real(8), parameter :: eps=0.1
    real(8) :: y
    jj = -1 
    imin = 0
    if ( mod(sum(ptn%n_ferm), 2) == 1 ) imin = 1
    
    do i = imin, 1000, 2
       y = x-i*(i+2)/4.d0
       if ( abs(y) < eps) then
          jj = i
          return
       end if
       if ( y < 0.d0) return
    end do
  end function guess_J_from_JJ


  subroutine wf_random(ptn, v, mtot, is_double_j)
    use bp_io, only: v_remove_j
    type(type_ptn_pn), intent(in) :: ptn
    type(type_vec_p), intent(inout) :: v
    integer, intent(in) :: mtot
    logical, intent(in) :: is_double_j
    integer :: jj
    type(type_vec_p) :: vt
    real(8) :: x

    call wf_random_vec(v, ptn)

    !  J-projection, not work well at kwf==4
    ! if (is_double_j .and. kwf==8) then  
    if (.false. .and. is_double_j .and. kwf==8) then  
       call wf_alloc_vec(vt, ptn)
       do jj = ptn%max_jj, abs(mtot)+2, -2
          if (myrank==0) write(*,'(a,1i3,a)') "J-projection remove ",jj,"/2"
          call matvec_jj(v, vt)
          call v_remove_j(v, jj, mtot, vt)
       end do
       call deallocate_l_vec( vt%p )
       call dotprod(v, v, x)
       v%p = 1.d0 / sqrt(x) * v%p
    end if
  end subroutine wf_random




  subroutine compress_init_vecs(evec, nb)
    type(type_vec_p), intent(inout) :: evec(:)
    integer, intent(in) :: nb
    integer :: i, n, iv, ii
    integer(kdim) :: mq
    real(8) :: x

    if (.not. associated(evec(1)%p)) stop 'ERROR compress_init_vecs'
    
    n = 0
    do iv = nb+1, size(evec)
       if ( associated(evec(iv)%p) ) then
          i = mod(iv-1, nb) + 1
          if (.not. associated(evec(i)%p)) stop 'ERROR compress_init_vecs'
          evec(i)%p = evec(i)%p + evec(iv)%p 
          n = n + 1
       else
          exit
       end if
    end do

    if (n == 0) return

    if (myrank==0) write(*,'(a,a,i5,a,i5,a/)') 'load initial ', &
         trim(fn_load_wave), nb+n, ' vecs => ', nb, ' vecs.'

    do iv = 1, nb
       evec(iv)%eval = 0.d0
       call dotprod(evec(iv), evec(iv), x)
       x = 1.d0 / sqrt(x)

       !$omp parallel do private (mq)
       do mq = 1, size( evec(iv)%p, kind=kwf )
          evec(iv)%p(mq) = evec(iv)%p(mq) * x
       end do
    end do

  end subroutine compress_init_vecs


  subroutine add_random_init_vecs(evec, nb, eps)
    use lib_matrix, only: gaussian_random_mat
    type(type_vec_p), intent(inout) :: evec(:)
    integer, intent(in) :: nb
    real(8), intent(in) :: eps
    integer :: i, n, iv, ii
    integer(kdim) :: mq
    real(8) :: x, ep
    type(type_vec_p) :: vt
    

    if (.not. associated(evec(1)%p)) stop 'ERROR compress_init_vecs'

    call allocate_l_vec( vt%p, ptn%max_local_dim )
    vt%p(ptn%local_dim:) = 0._kwf

    
    ep = eps / sqrt( dble(ptn%ndim) )
    if (myrank==0) write(*,'(/,a,i5,2f16.12/)') &
         'add random num. to initiall vec ', nb, eps, ep

    do iv = 1, nb
       call gaussian_random_mat( ptn%local_dim, vt%p ) 
    
       evec(iv)%p = evec(iv)%p + ep * vt%p
       call dotprod(evec(iv), evec(iv), x)
       x = 1.d0 / sqrt(x)

       !$omp parallel do private (mq)
       do mq = 1, size( evec(iv)%p, kind=kwf )
          evec(iv)%p(mq) = evec(iv)%p(mq) * x
       end do
    end do

    call deallocate_l_vec(vt%p)

  end subroutine add_random_init_vecs



  subroutine calc_tbme_from_tbtd(bp, evec)
    !   < T_ii >,  < V_ijklJ >  
    !   < T_ij > is not implemented
    use model_space
    use operator_jscheme, only : jcouple, jcouplemax
    use interaction, only : hamltn, hamltn_j
    use operator_mscheme, only: opr_m, operator_j2m, opr_m_p, &
         print_operator_mscheme, finalize_opr_m
    use operator_mscheme, only : init_tbtd_op, clear_op, get_cpld_tbtd
    use bp_expc_val, only : bp_ex_val_tbtd
    type(type_vec_p), intent(inout) :: evec(:)
    type(type_bridge_partitions), intent(inout) :: bp
    integer :: jj, ipn, iprty, n, ij12, ij34, k1, k2, k3, k4
    integer :: i, nop, iop
    real(8) :: v, x, e, occ(n_jorb_pn, size(evec))
    type(opr_m) :: op
    integer, allocatable :: idxs(:,:)
    real(8), allocatable :: vs(:), evs(:,:)

    nop = 0
    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             nop = nop + (n*(n+1)) / 2
          end do
       end do
    end do

    if (myrank==0) write(*,'(/,a,i8,i5)') 'TBME Num ', nop

    allocate( idxs(7,nop), vs(nop), evs(nop, size(evec)) )

    iop = 0 
    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             if (n == 0) cycle
             do ij12 = 1, n
                do ij34 = ij12, n
                   iop = iop + 1
                   if (iop > nop) stop 'error in calc_tbme_from_tbtd'
                   k1 = jcouple(jj, iprty, ipn)%idx(1, ij12)
                   k2 = jcouple(jj, iprty, ipn)%idx(2, ij12)
                   k3 = jcouple(jj, iprty, ipn)%idx(1, ij34)
                   k4 = jcouple(jj, iprty, ipn)%idx(2, ij34)
                   v = hamltn_j%p2(jj, iprty, ipn)%v(ij12, ij34) 
                   idxs(:,iop) = (/ k1, k2, k3, k4, jj, iprty, ipn /)
                   vs(iop) = v 
                end do
             end do
          end do
       end do
    end do


    do i = 1, size(evec) 
       if (.not. associated(evec(i)%p)) cycle
       call ex_occ_orb( evec(i), occ(:,i) )
    end do


    call init_tbtd_op(op, 0, 1, ptn%n_ferm, ptn%n_ferm)
    call init_bp_operator(bp, op, verbose=.true.)

    
    do i = 1, size(evec) 
       if (.not. associated(evec(i)%p)) cycle

       call clear_op( op )
       call bp_ex_val_tbtd( bp, evec(i), evec(i), op )

       call start_stopwatch(time_cpld_tbtd)

       !$omp parallel do private(iop, k1, k2, k3, k4, jj, x)
       do iop = 1, nop
          k1 = idxs(1, iop)
          k2 = idxs(2, iop)
          k3 = idxs(3, iop)
          k4 = idxs(4, iop)
          jj = idxs(5, iop)
          call get_cpld_tbtd( op, k1, k2, k3, k4, jj, jj, 0, x )
          if (k1 /= k3 .or. k2 /= k4) x = x * 2.d0
          evs(iop, i) = x * sqrt(2*jj + 1d0) 
       end do

       call stop_stopwatch(time_cpld_tbtd)
       
    end do

    call finalize_bp_operator(bp, op)


    if (myrank/=0) return ! print results

    call start_stopwatch(time_print_tbtd)

    write(*,'(a,f10.3,a)') &
         "time computation calc_tbdm", time_tmp%time, ' sec.'
    
    do i = 1, size(evec) 
       if (.not. associated(evec(i)%p)) cycle
       write(*,'(/,a,i7,/)') " ***** TBDM information ***** state = ", i
       write(*,'(a, i7,i3,i7,f12.5/)' )  "JJ, prty, TT,  E = ", &
            evec(i)%jj, ptn%iprty, evec(i)%tt, evec(i)%eval

       write(*,'(a, i5)') "OBDM,  i   j      Eij  &
            &           <E>  ",  n_jorb_pn
       e = 0.d0
       do k1 = 1, n_jorb_pn
          v = hamltn_j%p1%v(k1, k1) / sqrt(dble(jorb(k1)+1))
          write(*,'(a, 2i4, 2f14.7)') "OBDM ", k1, k1, v, occ(k1, i)
          e = e + v * occ(k1, i)
       end do

       write(*,'(a, i10)') &
            "TBDM,  i   j   k   l   J prty pn V_ijkl     <V>    ", nop
       do iop = 1, nop
          if (myrank==0) write(*,'(a, 7i4, 2f14.7)') &
               "TBDM ", idxs(:, iop), vs(iop), evs(iop, i)
          e = e + vs(iop) * evs(iop, i)
       end do
       write(*,'(a,2f12.5)') "confirm energy", e, evec(i)%eval
    end do
    
    call stop_stopwatch(time_print_tbtd)
    
  end subroutine calc_tbme_from_tbtd
  



  subroutine ex_val_snt()
    integer :: i, j
    type(type_vec_p) :: vt
    real(8) :: x, y
    if (myrank==0) write(*,'(/,a,/)') 'compute <vl| snt-file | vr> '
    call allocate_l_vec( vt%p, ptn%max_local_dim )

    do j = 1, n_eigen
       call matvec(evec(j), vt)
       do i = 1, n_eigen
          call dotprod(evec(i), vt, x)
          call dotprod(evec(i), evec(j), y)
          if (myrank==0) write(*,'(a,i3,f8.3,a,i3,f8.3,a,f10.5,a,f10.5)') &
               'i=',i, evec(i)%eval, ' j=',j, evec(j)%eval, &
               '    <i|snt|j>= ', x, '  <i|j> ',y
       end do
    end do
    call deallocate_l_vec( vt%p )
  end subroutine ex_val_snt



  subroutine ex_val_h2()
    ! energy variance
    integer :: i
    type(type_vec_p) :: vt
    real(8) :: x, y

    if (myrank==0) write(*,'(/,a,/)') '***  compute <|H^2|> - <|H|>^2  *** '
    call allocate_l_vec( vt%p, ptn%max_local_dim )

    do i = 1, n_eigen
       call matvec(evec(i), vt)
       call dotprod(vt, vt, x)
       call dotprod(evec(i), vt, y)
       
       ! if ( myrank==0 .and. abs( evec(i)%eval - y ) > 1.d-3 ) &
       !      write(*,'(a,i3,2f8.3)') "WARNING difference energy", i, evec(i)%eval, y 
       if ( myrank==0 ) &
            write(*,'(a,i3,3f12.4)') " energy variance  i, <H^2>-<H>^2,  <H^2>, <H> : ", &
            i, x-y**2, x, y
    end do

    call deallocate_l_vec( vt%p )

  end subroutine ex_val_h2

  
end program kshell

