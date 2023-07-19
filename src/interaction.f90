module interaction
  !$ use omp_lib
  use constant, only: pi
  use model_space
  use harmonic_oscillator, only: radius_j, nabla_j, init_ho_by_mass, &
       r_i_dash_nl, radius_power, radius_j
  use sp_matrix_element, only: r2_redmat, l_redmat, s_redmat, &
       ry_redmat, y_redmat, r3y1_redmat
  use rotation_group, only : init_rotation_group, dcg, d6j
  use operator_jscheme, only: opr_j, read_intfile, jcouple, set_opr_j, &
       non_diag_ob_2_tbme, print_operator_jscheme
  use operator_mscheme, only: opr_m, operator_j2m, print_operator_mscheme
  !
  implicit none
  private
  public :: read_interaction  ! read interaction file and set the following operators
  ! share operators in opr_m form
  public :: hamltn, jtensor, ltensor, stensor, r2y2, r1y1, r3y3, &
       ham_cm, num_orb, r3y1
  public :: hamltn_j, j_square, t_square, fm_m, nme_0v, &
       set_gt, gt_m, set_stst0

  ! first fobidden
  public :: set_fm, set_nme_0v, set_ff
  public :: set_r1y1t, set_sd0t, set_sd1t, set_sd2t, &
       r1y1t_m, sd0t_m, sd1t_m, sd2t_m

  type(opr_m), save :: hamltn, ham_cm, jtensor, ltensor, stensor, &
       r2y2, r1y1, r3y3, j_square, t_square, gt_m, fm_m, nme_0v, r3y1
  type(opr_m), save :: r1y1t_m, sd0t_m, sd1t_m, sd2t_m
  target :: jtensor, ltensor, stensor, r2y2, r1y1, r3y3, r3y1
  type(opr_m), allocatable :: num_orb(:)
  type(opr_j) :: hamltn_j, gt_j, fm_j, nme_0v_j
  type(opr_j), allocatable :: num_orb_j(:)
  integer :: n_num_orb  ! internal status for generate N_occ very private
  ! first forbidden beta decay
  public :: ff_0rs, ff_0sp, ff_1r, ff_1rs, ff_1p, ff_2rs, &
       ff_0rs_d, ff_1r_d, ff_1rs_d
  type(opr_m), save :: ff_0rs, ff_0sp, ff_1r, ff_1rs, ff_1p, ff_2rs, &
       ff_0rs_d, ff_1r_d, ff_1rs_d
  type(opr_j) :: ff_0rs_j, ff_0sp_j, ff_1r_j, ff_1rs_j, ff_1p_j, ff_2rs_j, &
       ff_0rs_d_j, ff_1r_d_j, ff_1rs_d_j
  type(opr_j) :: r1y1t_j, sd0t_j, sd1t_j, sd2t_j

  type(opr_j) :: r2y0_j
  type(opr_m), save :: r2y0
  target :: r2y0
  public :: r2y0_j, r2y0, r2y0_func1, set_e0

  public :: set_em, metensor, mmstensor, mmltensor
  type(opr_m), save :: metensor, mmstensor, mmltensor
  type(opr_j) :: metensor_j, mmstensor_j, mmltensor_j
  target :: metensor, mmstensor, mmltensor
  integer :: n_rank  ! internal status for generate M(EM) very private


  ! E1, E2 sum rule
  public :: set_ry_sum, rkyk_square, r1y1_f_square, sum_rank
  type(opr_m), save :: rkyk_square, r1y1_f_square
  integer :: sum_rank = -1  ! internal 
  real(8) :: sum_charge(2) ! internal 

  public :: gt_func1, r2y2_func1, r1y1_func1, r3y3_func1, &
       ltensor_func1, stensor_func1, &
       set_ob_channel, dummy_func1, set_ob_ij, r3y1_func1
  ! for  single-channel operator, private
  integer :: nljt_sc(4,2)

  ! for nme_0v, private
  integer :: n_nme_0v = 0
  integer, allocatable :: nme_0v_ijklJ(:,:)
  real(8), allocatable :: nme_0v_v(:)

  public :: set_three_body_monopole

contains

  subroutine read_interaction(lun, hw_type)
    !
    ! read interaction in lun and
    ! and set only hamltn, ham_cm, jtensor, ltensor, stensor, r2y2  
    !
    integer, intent(in) :: lun
    integer, intent(in), optional :: hw_type
    type(opr_j) :: jtensor_j, ltensor_j, stensor_j, r2y2_j, r1y1_j, &
         r3y3_j, hcm_j, j_square_j, t_square_j, r3y1_j
    !
    call init_rotation_group(maxval(jorb))
    !
    call read_intfile(lun, hamltn_j)
    call non_diag_ob_2_tbme(hamltn_j)
    ! call print_operator_jscheme(hamltn_j)
    call operator_j2m(hamltn_j, hamltn)
    ! call print_operator_mscheme(hamltn)
    if (present(hw_type)) call init_ho_by_mass(hw_type, mass)

    ! JJ
    call set_opr_j(j_square_j, j_square_func1, j_square_func2)
    call operator_j2m(j_square_j, j_square)
    j_square%is_j_square = .true.

    ! TT
    call set_opr_j(t_square_j, t_square_func1, t_square_func2)
    call operator_j2m(t_square_j, t_square)
    t_square%e0 = t_square_func0()
    
    ! Center of Mass hamiltonian
    call set_opr_j(hcm_j, hcm_func1, hcm_func2)
    call operator_j2m(hcm_j, ham_cm)

    ! J^(1)
    call set_opr_j(jtensor_j, jtensor_func1, irank=1)
    call operator_j2m(jtensor_j, jtensor)

    ! L^(1)
    call set_opr_j(ltensor_j, ltensor_func1, irank=1)
    call operator_j2m(ltensor_j, ltensor)

    ! S^(1)
    call set_opr_j(stensor_j, stensor_func1, irank=1)
    call operator_j2m(stensor_j, stensor)

    ! r2Y2^(2)
    call set_opr_j(r2y2_j, r2y2_func1, irank=2)
    call operator_j2m(r2y2_j, r2y2)

    ! r1Y1^(1)
    call set_opr_j(r1y1_j, r1y1_func1, irank=1, ipr1_type=-1)
    call operator_j2m(r1y1_j, r1y1)

    ! r3Y3^(3)
    call set_opr_j(r3y3_j, r3y3_func1, irank=3, ipr1_type=-1)
    call operator_j2m(r3y3_j, r3y3)

    ! r3Y1^(1)
    call set_opr_j(r3y1_j, r3y1_func1, irank=1, ipr1_type=-1)
    call operator_j2m(r3y1_j, r3y1)
    

  
  end subroutine read_interaction


 
  subroutine set_fm()
    ! set Fermi transition operator
    call set_opr_j(fm_j, fm_func1, irank=0, ipr1_type=1, nbody=-10)
    call operator_j2m(fm_j, fm_m)
  end subroutine set_fm

  subroutine set_gt()
    ! set Gamow-Teller operator
    call set_opr_j(gt_j, gt_func1, irank=1, ipr1_type=1, nbody=-10)
    ! call print_operator_jscheme(gt_j)
    call operator_j2m(gt_j, gt_m)
  end subroutine set_gt


  subroutine set_r1y1t()
    ! set r1Y1*tau 
    call set_opr_j(r1y1t_j, r1y1t_func1, irank=1, ipr1_type=-1, nbody=-10)
    ! call print_operator_jscheme(r1y1t_j)
    call operator_j2m(r1y1t_j, r1y1t_m)
  end subroutine set_r1y1t

  subroutine set_sd0t()
    ! set spin-dipole rank=0 * tau 
    call set_opr_j(sd0t_j, sd0t_func1, irank=0, ipr1_type=-1, nbody=-10)
    ! call print_operator_jscheme(sd0t_j)
    call operator_j2m(sd0t_j, sd0t_m)
  end subroutine set_sd0t

  subroutine set_sd1t()
    ! set spin-dipole rank=1 * tau 
    call set_opr_j(sd1t_j, sd1t_func1, irank=1, ipr1_type=-1, nbody=-10)
    ! call print_operator_jscheme(sd1t_j)
    call operator_j2m(sd1t_j, sd1t_m)
  end subroutine set_sd1t

  subroutine set_sd2t()
    ! set spin-dipole rank=2 * tau 
    call set_opr_j(sd2t_j, sd2t_func1, irank=2, ipr1_type=-1, nbody=-10)
    ! call print_operator_jscheme(sd2t_j)
    call operator_j2m(sd2t_j, sd2t_m)
  end subroutine set_sd2t


  subroutine set_e0
    ! E0 operator
    call set_opr_j(r2y0_j, r2y0_func1, irank=0)
    call operator_j2m(r2y0_j, r2y0)
  end subroutine set_e0




  function jtensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| J || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    real(8) :: rj
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    rj = dble(j1)*0.5d0
    r = sqrt(rj*(rj+1.d0)*(2.d0*rj+1.d0))
  end function jtensor_func1

  function ltensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| L || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2) return
    r = l_redmat(n1, l1, j1, n2, l2, j2)
  end function ltensor_func1

  function stensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| S || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2) return
    r = s_redmat(n1, l1, j1, n2, l2, j2)
  end function stensor_func1

  function num_orb_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| N_orb (n_num_orb) || nljt2>
    ! Note: n_num_orb state dependent
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    if (n1/=norb(n_num_orb) .or. l1/=lorb(n_num_orb) &
         & .or. j1/=jorb(n_num_orb) .or. t1/=itorb(n_num_orb)) return
    r = sqrt(dble(j1 + 1))
  end function num_orb_func1

  function r2y2_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r2Y2 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)/=mod(l2,2)) return
    r = ry_redmat(2, n1, l1, j1, n2, l2, j2)
  end function r2y2_func1

  function r1y1_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r1Y1 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)==mod(l2,2)) return
    r = ry_redmat(1, n1, l1, j1, n2, l2, j2)
  end function r1y1_func1

  function r3y3_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r3Y3 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)==mod(l2,2)) return
    r = ry_redmat(3, n1, l1, j1, n2, l2, j2)
  end function r3y3_func1


  function r3y1_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r3Y1 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)==mod(l2,2)) return
    r = r3y1_redmat(n1, l1, j1, n2, l2, j2)
  end function r3y1_func1
  
  
  function hcm_func1(n1, l1, j1, t1,  n2, l2, j2, t2) result (r)
    ! <nljt1|| H_cm ||nljt2>  one-body part
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: i, nc, lc, jc, tc, jj
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    r = dble(2*n1 + l1) * sqrt(dble(j1+1))

    do i = 1, n_nljt_core
       nc = nljt_core(1, i)
       lc = nljt_core(2, i)
       jc = nljt_core(3, i)
       tc = nljt_core(4, i)

       do jj = abs(jc-j1)/2, abs(jc+j1)/2
          r = r + dble(2*jj+1) &
               * hcm_func2(nc, lc, jc, tc, n1, l1, j1, t1, &
               &           nc, lc, jc, tc, n2, l2, j2, t2, jj) &
               / sqrt(dble(j1+1))
       end do
    end do
  end function hcm_func1


  
  function hcm_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <j1, j2 | H_cm (2-body) |j3 j4>_JJ 
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2/=t3+t4) return
    if ( mod(l1+l2, 2)/=mod(l3+l4, 2)) return    
    r = (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
         &     * ( - nabla_j(n1,l1,j1,n3,l3,j3)*nabla_j(n2,l2,j2,n4,l4,j4) &
         &         + radius_j(n1,l1,j1,n3,l3,j3)*radius_j(n2,l2,j2,n4,l4,j4))
    if (t3==t4) &  ! exchange term
         r = r - (-1.d0)**((j3+j4)/2-JJ) &
         &     * (-1.d0)**((j2+j4)/2-JJ) * d6j(j1, j2, 2*JJ, j3, j4, 2) &
         &     * ( - nabla_j(n1,l1,j1,n4,l4,j4)*nabla_j(n2,l2,j2,n3,l3,j3) &
         &         + radius_j(n1,l1,j1,n4,l4,j4)*radius_j(n2,l2,j2,n3,l3,j3) )
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  end function hcm_func2


  function sqr_radius_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! return zero, point-particle radius, one-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0 
  end function sqr_radius_func1

  function sqr_radius_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! point-particle matter radius, 1/A^2 \sum_i r_i^2  
    !    in intrinsic wavefunction 
    !  = 1/A^2 <j1, j2 | (r_i - r_j)^2 |j3 j4>_JJ 
    use constant, only: pi
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "sqr_radius_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    r = r /(dble(n_ferm_pn**2))
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)

  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | (r_1 - r_2)^2 | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = r2_redmat(n1,l1,j1,n3,l3,j3)/sqrt(dble(j1+1)) &
           &   * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
           & + r2_redmat(n2,l2,j2,n4,l4,j4)/sqrt(dble(j2+1)) &
           &   * delta(n1,l1,j1,t1, n3,l3,j3,t3) &
           & - 2.d0 * (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           & * ry_redmat(1,n1,l1,j1,n3,l3,j3)*ry_redmat(1,n2,l2,j2,n4,l4,j4) &
           & * 4.d0*pi/3.d0
    end function direct_term
  end function sqr_radius_func2


  function delta(n1,l1,j1,t1, n2,l2,j2,t2) result (r)
    ! Kronecker's delta
    integer, intent(in) :: n1,l1,j1,t1, n2,l2,j2,t2
    real(8) :: r
    r = 0.d0 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r=1.d0
  end function delta


  function chg_radius_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! point-proton radius, one-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0 
    if (mod(l1, 2) /= mod(l2, 2)) return
    if (t1==-1 .and. t2==-1) then
       r = r2_redmat(n1,l1,j1,n2,l2,j2) &
            * dble(n_ferm_pn**2 - 2*n_ferm_pn + n_ferm(1)) &
            / dble(n_ferm(1)*n_ferm_pn**2)
    else if (t1==1 .and. t2==1) then
       r = r2_redmat(n1,l1,j1,n2,l2,j2) / dble(n_ferm_pn**2)
    end if
  end function chg_radius_func1

  function chg_radius_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! point-proton matter radius, two-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "sqr_charge_radius_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (t1==-1 .and. t2==-1 .and. t3==-1 .and. t4==-1) then
       r = r * dble(2*n_ferm(1)-4*n_ferm_pn)/dble(n_ferm(1)*n_ferm_pn**2)
    else if (t1==1 .and. t2==1 .and. t3==1 .and. t4==1) then
       r = r * 2.d0/dble(n_ferm_pn**2)
    else 
       r = r * 2.d0*dble((n_ferm(1)-n_ferm_pn))/dble(n_ferm(1)*n_ferm_pn**2)
    end if
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | r_1 * r_2 | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           * ry_redmat(1,n1,l1,j1,n3,l3,j3)*ry_redmat(1,n2,l2,j2,n4,l4,j4) &
           * 4.d0*pi/3.d0
    end function direct_term
  end function chg_radius_func2


  function j_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nlj1 || JJ || njl2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r, rj
    rj = dble(j1)*0.5d0
    r = delta(n1,l1,j1,t1, n2,l2,j2,t2) * rj*(rj + 1.d0) * sqrt(dble(j1+1))
  end function j_square_func1

  function j_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <nlj1 njl2 | JJ | nlj3 nlj4>_J
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2 /= t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | J * J | j3j4>_JJ  direct 
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r, rj1, rj2
      rj1 = dble(j1)*0.5d0
      rj2 = dble(j2)*0.5d0
      r = 2.d0 * (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           * delta(n1,l1,j1,t1, n3,l3,j3,t3) &
           * sqrt( rj1*(rj1 + 1.d0)*(2.d0*rj1 + 1.d0) ) &
           * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
           * sqrt( rj2*(rj2 + 1.d0)*(2.d0*rj2 + 1.d0) ) 
    end function direct_term
  end function j_square_func2


  function t_square_func0() result (r)
    real(8) :: r
    r = abs(n_core(2) - n_core(1) ) * 0.5d0
    r = r * (r + 1d0)    
  end function t_square_func0
  

  function t_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nlj1 || TT || njl2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: ic, nc, lc, jc, itc, jj
    
    r = delta(n1,l1,j1,t1, n2,l2,j2,t2) * 0.75d0 * sqrt(dble(j1+1))

    if ( n_core(1) == n_core(2) ) return
    
    do ic = 1, n_nljt_core
       nc  = nljt_core(1, ic)
       lc  = nljt_core(2, ic)
       jc  = nljt_core(3, ic)
       itc = nljt_core(4, ic)

       do jj = abs(jc-j1)/2, abs(jc+j1)/2
          r = r + dble(2*jj+1) * t_square_func2(nc, lc, jc, itc, n1, l1, j1, t1, &
               & nc, lc, jc, itc, n2, l2, j2, t2, jj) / sqrt(dble(j1+1))
       end do
    end do
    
  end function t_square_func1
  

  function t_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <nlj1 njl2 | TT | nlj3 nlj4>_J
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2 /= t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ)
    r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ)
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | T * T | j3j4>_JJ  direct
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      if (t1==t2) then
         r = 0.5d0 * delta(n1,l1,j1,t1, n3,l3,j3,t3) * delta(n2,l2,j2,t2, n4,l4,j4,t4)
      else
         r = - 0.5d0 * delta(n1,l1,j1,t1, n3,l3,j3,t3) * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
              + delta(n1,l1,j1,t1, n3,l3,j3,t4) * delta(n2,l2,j2,t2, n4,l4,j4,t3)
      end if
    end function direct_term
  end function t_square_func2


  function fm_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body Fermi transition op.  <nljt1|| t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    if (n1/=n2 .or. l1/=l2 .or. j1 /= j2) return
    r = sqrt(j1 + 1.d0)
  end function fm_func1

  function gt_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body Gamow-Teller op.  <nljt1|| simga*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = 2.d0 * s_redmat(n1, l1, j1, n2, l2, j2)
  end function gt_func1



  function r2y0_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! E0 transition matrix element
    ! one-body proton <nljt1|| r2 || nljt2> ! / R^2  R=1.2A^(1/3)
    ! NOT r2Y0
    ! 
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    ! if (t1 /= -1 .or. t2 /= -1) return
    r = r2_redmat(n1,l1,j1, n2,l2,j2)  ! / (1.2d0 * mass**(1.d0/3.d0))**2
  end function r2y0_func1
  



! ---------------------------------------------------------------



  function r1y1t_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body r1Y1*tau   <nljt1|| r1Y1*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    if (mod(l1,2)==mod(l2,2)) return
    r = ry_redmat(1, n1, l1, j1, n2, l2, j2)
  end function r1y1t_func1


  function sd0t_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body spin-dipole rank=0 tau operator 
    ! <nljt1|| [rY1*sigma]^0*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(0, n1, l1, j1, n2, l2, j2) * sqrt( 3.d0 / 4.d0 / pi)
  end function sd0t_func1

  function sd1t_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body spin-dipole rank=1 tau operator 
    ! <nljt1|| [rY1*sigma]^0*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(1, n1, l1, j1, n2, l2, j2) * sqrt( 3.d0 / 4.d0 / pi)
  end function sd1t_func1

  function sd2t_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body spin-dipole rank=0 tau operator 
    ! <nljt1|| [rY1*sigma]^0*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(2, n1, l1, j1, n2, l2, j2) * sqrt( 3.d0 / 4.d0 / pi)
  end function sd2t_func1


! ---------------------------------------------------------------

  subroutine set_ry_sum(irank, charge_input)
    integer, intent(in) :: irank
    real(8), intent(in) :: charge_input(:)
    type(opr_j) :: rkyk_square_j, r1y1_f_square_j
    
    sum_charge(:) = charge_input(:)
    sum_rank = irank

    ! er^kY^k*er^kY^k, intermediate state in model space
    call set_opr_j(rkyk_square_j, rkyk_square_func1, rkyk_square_func2)
    call non_diag_ob_2_tbme(rkyk_square_j)
    call operator_j2m(rkyk_square_j, rkyk_square)

    if (irank==1) then
       ! er1Y1*er1Y1, intermediate state in full space
       call set_opr_j(r1y1_f_square_j, r1y1_f_square_func1, rkyk_square_func2)
       call non_diag_ob_2_tbme(r1y1_f_square_j)
       call operator_j2m(r1y1_f_square_j, r1y1_f_square)
       r1y1_f_square%e0 = r1y1_f_square_func0() 
    end if

  end subroutine set_ry_sum



  function rkyk_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! E(sum_irank) sum rule : erkyk*P*erkyk one-body term 
    !   P ... Projection opertor to model space
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: i, nk, lk, jk
    r = 0.d0 
    if (mod(l1, 2) /= mod(l2, 2)) return
    if (t1 /= t2) return
    if (j1 /= j2) return
    do i = 1, n_jorb_pn
       if (itorb(i) /= t1) cycle
       nk = norb(i)
       lk = lorb(i)
       jk = jorb(i)
       r = r + ry_redmat(sum_rank, n1, l1, j1, nk, lk, jk) &
            *  ry_redmat(sum_rank, n2, l2, j2, nk, lk, jk)
    end do
    r = r / sqrt(dble(j2+1.d0)) * sum_charge( (t1+1)/2+1 )**2
    
  end function rkyk_square_func1



  function rkyk_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! E(k) sum rule : er1y1*er1y1 two-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "rkyk_square_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,   JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | erY(k) * erY(k) | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = 2.d0 * (-1.d0)**((j2+j3)/2-JJ) &
           * d6j( j1, j2, 2*JJ, j4, j3, 2*sum_rank ) &
           * ry_redmat(sum_rank, n1,l1,j1, n3,l3,j3) & 
           * ry_redmat(sum_rank, n2,l2,j2, n4,l4,j4) &
           * sum_charge((t1+1)/2+1) * sum_charge((t2+1)/2+1)
    end function direct_term
  end function rkyk_square_func2



  function r1y1_f_square_func0() result (r)
    ! er1y1*er1y1 zero-body term (constant)
    ! intermediate state is in full space (beyond model space)
    ! <core | rY1 * rY1 | core>
    use model_space, only: n_core
    real(8) :: r
    integer :: i, j, k, nsum, maxosc, nosc, nnljtc
    integer :: jj, nc, lc, jc, it, itc, n1, l1, j1, t1, n2, l2, j2, t2, sk

    if (sum_rank /= 1) stop 'ERROR in ry1y_f_square_func0'
    if (n_core(1)<0 .or. n_core(2)<0) stop 'Not implmented r1y1_f_square_func0'

    r = 0.d0
    do i = 1, n_nljt_core
       n1 = nljt_core(1, i)
       l1 = nljt_core(2, i)
       j1 = nljt_core(3, i)
       t1 = nljt_core(4, i)
       r = r + 3.d0/4.d0/pi * r2_redmat( n1,l1,j1, n1,l1,j1 ) &
            * sum_charge( (t1+1)/2+1 )**2 * sqrt( j1 + 1d0 )

       do j = i, n_nljt_core
          n2 = nljt_core(1, j)
          l2 = nljt_core(2, j)
          j2 = nljt_core(3, j)
          t2 = nljt_core(4, j)
          
          sk = 1
          if (i==j) sk = 2
          do jj = abs(j1-j2)/2, abs(j1+j2)/2, sk
             r = r + dble(2*jj+1) &
                  * rkyk_square_func2( &
                  &   n1, l1, j1, t1, n2, l2, j2, t2, &
                  &   n1, l1, j1, t1, n2, l2, j2, t2, jj)
          end do
       end do
    end do

  end function r1y1_f_square_func0



  function r1y1_f_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    !
    ! er1y1*er1y1 one-body term 
    ! intermediate state is in full space (beyond model space)
    !
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: i, jj, nc, lc, jc, tc
    r = 0.d0 
    if (mod(l1, 2) /= mod(l2, 2)) return
    if (t1 /= t2) return
    if (j1 /= j2) return

    if (sum_rank /= 1) stop 'ERROR in ry1y_f_square_func1'

    r = 3.d0/4.d0/pi * r2_redmat( n1,l1,j1, n2,l2,j2 ) &
         * sum_charge( (t1+1)/2+1 )**2 

    do i = 1, n_nljt_core
       nc = nljt_core(1, i)
       lc = nljt_core(2, i)
       jc = nljt_core(3, i)
       tc = nljt_core(4, i)

       do jj = abs(jc-j1)/2, abs(jc+j1)/2
          r = r + dble(2*jj+1) &
               * rkyk_square_func2( &
               &  nc, lc, jc, tc, n1, l1, j1, t1, &
               &  nc, lc, jc, tc, n2, l2, j2, t2, jj) &
               / sqrt( j1 + 1d0 )
       end do
    end do

  end function r1y1_f_square_func1




  subroutine set_nme_0v(ipn)
    ! Nuclear matrix element for 0-nu double-beta decay, closure approx.
    integer, intent(in) :: ipn
    call set_opr_j(nme_0v_j, dummy_func1, nme_0v_func2, nbody=-11-ipn)
    call operator_j2m(nme_0v_j, nme_0v)
  end subroutine set_nme_0v


  function dummy_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! dummy for one-body, nme_0v
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
  end function dummy_func1

  function nme_0v_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! Nuclear matrix element for 0-nu double-beta decay, closure approx.
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    integer :: i
    logical :: is_found 
    
    r = 0.d0
    if (t1/=t2 .or. t3/=t4 .or. t1==t3) return

    if (n_nme_0v == 0) call read_nme_file()

    is_found = .false.
    
    do i = 1, n_nme_0v
       if ( all( nme_0v_ijklJ(:,i) == (/ &
            n1, l1, j1, n2, l2, j2, &
            n3, l3, j3, n4, l4, j4, JJ /) )) then 
          r = nme_0v_v(i)
          is_found = .true. 
          exit
       end if
    end do

    if (myrank==0 .and. .not. is_found) &
         write(*,'(a,17i3)') 'WARNING: not found nme_0v', &
         n1, l1, j1, t1,  n2, l2, j2, t2, &
         n3, l3, j3, t3,  n4, l4, j4, t4,  JJ

!    if (myrank==0) write(*,'(a,17i3,f10.5)') 'read mele ', &
!         n1, l1, j1, t1,  n2, l2, j2, t2, n3, l3, j3, t3,  n4, l4, j4, t4,  JJ, r
    

  contains

    subroutine read_nme_file()
      integer, parameter :: lun=14
      integer :: nii, lii, jii, njj, ljj, jjj, &
              nkk, lkk, jkk, nll, lll, jll, i_jj
      real(8) :: v
      integer :: i

      ! adhoc file name
      open(lun, file='nme.dat', status='old')

      n_nme_0v = 0
      do while (.true.)
         read(lun, *, end=10) &
              nii, lii, jii, njj, ljj, jjj, &
              nkk, lkk, jkk, nll, lll, jll, i_jj, v
         n_nme_0v = n_nme_0v + 1
      end do
10    continue

      rewind( lun )

      allocate( nme_0v_ijklJ(13, n_nme_0v), nme_0v_v(n_nme_0v) )

      do i = 1, n_nme_0v
         read(lun, *) &
              nii, lii, jii, njj, ljj, jjj, &
              nkk, lkk, jkk, nll, lll, jll, i_jj, v
         nme_0v_ijklJ(:, i) = (/ &
              nii, lii, jii, njj, ljj, jjj, &
              nkk, lkk, jkk, nll, lll, jll, i_jj /)
         nme_0v_v(i) = v
      end do

      close(lun)

    end subroutine read_nme_file
    
  end function nme_0v_func2



  subroutine set_ff()
    ! set first forbidden operators
    ! Rank 0,  rs, sp
    call set_opr_j(ff_0rs_j, ff_0rs_func1, irank=0, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_0rs_j, ff_0rs)
    call set_opr_j(ff_0sp_j, ff_0sp_func1, irank=0, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_0sp_j, ff_0sp)
    ! Ranks 0, rs' 
    call set_opr_j(ff_0rs_d_j, ff_0rs_d_func1, irank=0, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_0rs_d_j, ff_0rs_d)
    ! Rank 1,  r, rs,  p
    call set_opr_j(ff_1r_j,  ff_1r_func1,  irank=1, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_1r_j,  ff_1r)
    call set_opr_j(ff_1rs_j, ff_1rs_func1, irank=1, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_1rs_j, ff_1rs)
    call set_opr_j(ff_1p_j,  ff_1p_func1,  irank=1, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_1p_j,  ff_1p)
    ! Rank 1 r', rs'
    call set_opr_j(ff_1r_d_j,  ff_1r_d_func1,  irank=1, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_1r_d_j,  ff_1r_d)
    call set_opr_j(ff_1rs_d_j, ff_1rs_d_func1, irank=1, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_1rs_d_j, ff_1rs_d)
    ! Rank 2 unique first forbidden
    call set_opr_j(ff_2rs_j, ff_2rs_func1, irank=2, ipr1_type=-1, nbody=-10)
    call operator_j2m(ff_2rs_j, ff_2rs)
  end subroutine set_ff

  ! Rank 0
  function ff_0rs_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| (rC1*simga)^0 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(0, n1, l1, j1, n2, l2, j2)
  end function ff_0rs_func1

  function ff_0sp_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| (sigma*nabla)^0 * t(+,-) || nljt2>
    use harmonic_oscillator, only: bpar
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: jt
    r = 0.d0
    if (t1==t2) return
    do jt = max(1, 2*l1-1), 2*l1+1, 2
       r = r + d6j(2, 2, 0, j2, j1, jt) &
            * 2.d0 * s_redmat(n1, l1, j1, n1, l1, jt) &
            * nabla_j(n1, l1, jt, n2, l2, j2) / bpar
    enddo
    r = r * (-1)**((j1+j2)/2) 

    ! if (abs(r)<1.d-8) return
    ! x = 0.d0
    ! do jt = max(1, 2*l2-1), 2*l2+1, 2
    !    x = x + d6j(2, 2, 0, j2, j1, jt) &
    !         * nabla_j(n1, l1, j1, n2, l2, jt) / bpar &
    !         * 2.0d0 * s_redmat(n2, l2, jt, n2, l2, j2)
    ! enddo
    ! x = x * (-1)**((j1+j2)/2) 
    ! write(*,*)"check",2*n1+l1,2*n2+l2,r,x
  end function ff_0sp_func1

  function ff_0rs_d_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| 2/3I(r)*r(C1*simga)^0 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(0, n1, l1, j1, n2, l2, j2, is_i_dash=.true.)
  end function ff_0rs_d_func1


  ! Rank 1
  function ff_1r_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| rC1 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = ry_redmat(1, n1, l1, j1, n2, l2, j2) * sqrt( 4.d0*pi/3.d0 )
  end function ff_1r_func1

  function ff_1rs_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| (rC1*simga)^1 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(1, n1, l1, j1, n2, l2, j2)
  end function ff_1rs_func1


  function ff_1p_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    use harmonic_oscillator, only: bpar
    ! <nljt1|| nabla * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = nabla_j(n1, l1, j1, n2, l2, j2) / bpar
  end function ff_1p_func1


  function ff_1r_d_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| 2/3I*r*C1 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = r_i_dash_nl(n1, l1, n2, l2) &
         * y_redmat(1, l1, j1, l2, j2) * sqrt( 4.d0*pi/3.d0 )
  end function ff_1r_d_func1

  function ff_1rs_d_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| 2/3I*r*(C1*simga)^1 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(1, n1, l1, j1, n2, l2, j2, is_i_dash=.true.)
  end function ff_1rs_d_func1

  !Rank 2
  function ff_2rs_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nljt1|| (rC1*simga)^2 * t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = c1sigk(2, n1, l1, j1, n2, l2, j2)
  end function ff_2rs_func1


  function c1sigk( k, n1, l1, j1, n2, l2, j2, is_i_dash ) result (r)
    !
    !   <n1 l1 j1|| r[C^(1)*sigma]^k ||n2 l2 j2>
    !   or, if (is_i_dash) 
    !   <n1 l1 j1|| 2/3I(1,1,1,1,;r)*r[C^(1)*sigma]^k ||n2 l2 j2>
    !
    integer, intent(in) :: k, n1, l1, j1, n2, l2, j2
    logical, intent(in), optional :: is_i_dash
    real(8) :: r
    integer :: jt
    logical :: is
    real(8) :: x

    if (k<0 .or. k>2) stop 'error [c1sigk]: k should be 0 to 2'

    r = 0.d0
    if ( mod(l1,2) == mod(l2,2) ) return

    is = .false. 
    if (present(is_i_dash)) is = is_i_dash
    if (is) then 
       x = r_i_dash_nl(n1, l1, n2, l2)
    else
       if (abs(l1-l2)>1 .or. l1+l2<1) return
       x = radius_power(1, n1, l1, n2, l2)
    end if
       
    do jt = max(1, 2*l2-1), 2*l2+1, 2
       r = r + d6j(2, 2, 2*k, j2, j1, jt) &
            * x * y_redmat(1, l1, j1, l2, jt) &
            * 2.0d0 * s_redmat(n2, l2, jt, n2, l2, j2)
    enddo
    r = r * (-1)**(k+(j1+j2)/2) * sqrt( 2.d0*k + 1.d0 ) &
         * sqrt( 4.d0*pi/3.d0 )
  end function c1sigk



!------------------------------------------------
  subroutine set_stst0(op)
    ! set [sigma*tau+ sigma*tau+]^(0) for double Gamow-Teller 
    type(opr_m), intent(out) :: op
    type(opr_j) :: op_j
    ! set 
    call set_opr_j(op_j, dummy_func1, stst0_func2, nbody=-12)
    ! call print_operator_jscheme(sd2t_j)
    call operator_j2m(op_j, op)
  end subroutine set_stst0


  function stst0_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <nlj1 njl2 | [sigma*tau sigmat*tau]^(0) | nlj3 nlj4>_J
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r

    r = 0.d0

    if ( any( (/ t1, t2, t3, t4 /) /= (/ -1, -1, 1, 1 /) ) ) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return

    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
    ! write(*,*) '<j1j2||DGT0||j3j4>_J',n1, l1, j1, t1,  n2, l2, j2, t2, &
    !      & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ, r

  contains

    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | [sigma * sigma]^(0) | j3j4>_JJ  direct 
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r

      r = 2.d0 * (-1.d0)**((j2+j3)/2-JJ) &
           * d6j( j1, j2, 2*JJ, j4, j3, 2*1 ) &
           * s_redmat(n1, l1, j1, n3, l3, j3) &
           * s_redmat(n2, l2, j2, n4, l4, j4) &
           * 4.d0 
      r = - r / sqrt(3.d0)
    end function direct_term
  end function stst0_func2




!  ---------------------------------------
  subroutine set_em(rank)
    integer, intent(in) :: rank
    n_rank = rank
    call set_opr_j(metensor_j, metensor_func1, irank=rank, &
         ipr1_type=(-1)**n_rank)
    call operator_j2m(metensor_j, metensor)
    call set_opr_j(mmstensor_j, mmstensor_func1, irank=rank-1, &
         ipr1_type=(-1)**n_rank)
    call operator_j2m(mmstensor_j, mmstensor)
    call set_opr_j(mmltensor_j, mmltensor_func1, irank=rank-1, &
         ipr1_type=(-1)**n_rank)
    call operator_j2m(mmltensor_j, mmltensor)
  end subroutine set_em


  function metensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r^{n_rank}Y_{n_rank} || nljt2>
    use harmonic_oscillator, only: radius_power
    use constant, only: pi
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1+l2+n_rank,2)==1) return
    if (n_rank+l1<l2 .or. n_rank+l2<l1) return
    r = (-1)**((j1-1)/2+l1+l2) * sqrt((j1+1)*(j2+1)/4.d0/pi) * dcg(j1,1,j2,-1,n_rank*2,0) &
        * radius_power(n_rank, n1, l1, n2, l2)
  end function metensor_func1

  function mmstensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| M_s(M{n_rank-1}) || nljt2>
    ! M = g_s M_s + g_l M_l
    use harmonic_oscillator, only: radius_power
    use constant, only: pi
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1+l2+n_rank,2)==1) return
    if (n_rank-1+l1<l2 .or. n_rank-1+l2<l1) return
    r = (-1)**((j1-1)/2+l1+l2) * sqrt((j1+1)*(j2+1)/4.d0/pi) * dcg(j1,1,j2,-1,n_rank*2-2,0) &
        * (n_rank-1 + (-1)**((1-j1)/2+l1) * (j1+1)/2 + (-1)**((1-j2)/2+l2) * (j2+1)/2) &
        / 2.d0 * radius_power(n_rank-2, n1, l1, n2, l2)
  end function mmstensor_func1

  function mmltensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1 || M_l(M{n_rank-1}) || nljt2>
    ! M = g_s M_s + g_l M_l
    use harmonic_oscillator, only: radius_power
    use constant, only: pi
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1+l2+n_rank,2)==1) return
    if (n_rank-1+l1<l2 .or. n_rank-1+l2<l1) return
    r = (-1)**((j1-1)/2+l1+l2) * sqrt((j1+1)*(j2+1)/4.d0/pi) * dcg(j1,1,j2,-1,n_rank*2-2,0) &
        * (n_rank-1 + (-1)**((1-j1)/2+l1) * (j1+1)/2 + (-1)**((1-j2)/2+l2) * (j2+1)/2) &
        * (-1 + 1.0/n_rank * ((-1)**((1-j1)/2+l1) * (j1+1)/2 + (-1)**((1-j2)/2+l2) * (j2+1)/2)) &
        * radius_power(n_rank-2, n1, l1, n2, l2)
  end function mmltensor_func1

!  ---------------------------------------



  subroutine set_ob_ij(ii, ij, irank, op)
    ! operator  [c+_ii x  c_ij]^irank 
    integer, intent(in) :: ii, ij, irank
    type(opr_m), intent(out) :: op
    integer :: iprty, nbody
    type(opr_j) :: oj
    
    iprty = iporb(ii)*iporb(ij)
    if (     jorb(ii)+jorb(ij)  < 2*irank  .or. &
         abs(jorb(ii)-jorb(ij)) > 2*irank ) stop 'rank failed in set_ob_ij'

    if (     itorb(ii) == itorb(ij) ) then 
       nbody =   1
    else if (itorb(ii) == -1 .and. itorb(ij) ==  1 ) then 
       nbody = -10
    else if (itorb(ii) ==  1 .and. itorb(ij) == -1 ) then 
       nbody = -11
    else
       stop 'ERROR: set_ob_ij'
    end if

    nljt_sc(:,1) =  (/ norb(ii), lorb(ii), jorb(ii), itorb(ii) /)
    nljt_sc(:,2) =  (/ norb(ij), lorb(ij), jorb(ij), itorb(ij) /)

    call set_opr_j(oj, single_channel_func1, &
         irank=irank, ipr1_type=iprty, nbody=nbody)
    call operator_j2m(oj, op)

  end subroutine set_ob_ij



  subroutine set_ob_channel(irank, iprty, nbody, ops, ij_orb, iorbl, iorbr)
    integer, intent(in) :: irank, iprty, nbody
    type(opr_m), intent(out), allocatable :: ops(:)
    integer, intent(out), allocatable :: ij_orb(:,:)
    integer, optional, intent(in) :: iorbl, iorbr
    integer :: i, j, it, jt, iloop, iop
    type(opr_j) :: oj
    
    do iloop = 1, 2
       iop = 0
       do i = 1, n_jorb_pn
          if (present(iorbl)) then
             if (i /= iorbl) cycle
          end if
          nljt_sc(:,1) =  (/ norb(i), lorb(i), jorb(i), itorb(i) /)
          do j = 1, n_jorb_pn
             if (present(iorbr)) then
                if ( j /= iorbr ) cycle
             end if
             if ( iporb(i)*iporb(j) /= iprty ) cycle
             if (     jorb(i)+jorb(j)  < 2*irank  .or. &
                  abs(jorb(i)-jorb(j)) > 2*irank ) cycle


             if ( nbody == 1 ) then
                if ( itorb(i) /= itorb(j) ) cycle
             else if ( nbody == -10 ) then
                if ( itorb(i) /= -1 .or. itorb(j) /=  1 ) cycle
             else if (  nbody == -11) then
                if ( itorb(i) /=  1 .or. itorb(j) /= -1 ) cycle
             else
                stop 'ERROR: set_ob_channel'
             end if

             iop = iop + 1

             if (iloop==2) then
                nljt_sc(:,2) = (/ norb(j), lorb(j), jorb(j), itorb(j) /)
                ij_orb(:,iop) = (/ i, j /)
                call set_opr_j(oj, single_channel_func1, &
                     irank=irank, ipr1_type=iprty, nbody=nbody)
                call operator_j2m(oj, ops(iop))
             end if

          end do
       end do
       if (iloop == 1) allocate( ops(iop), ij_orb(2,iop) )
    end do
  end subroutine set_ob_channel



  function single_channel_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| J || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if ( any( (/ n1,l1,j1,t1 /) /= nljt_sc(:,1) ) ) return
    if ( any( (/ n2,l2,j2,t2 /) /= nljt_sc(:,2) ) ) return
    r = 1.d0
  end function single_channel_func1


  !--------------- set 3-body monopole interaction ---

  subroutine set_three_body_monopole(fname, op)
    character(len=*), intent(in) :: fname
    type(opr_m), intent(inout) :: op
    integer, parameter :: lun=14
    integer :: i, j, n, ijk(3), ipn
    real(8) :: v

    open(lun, file=fname, status='old')
    call skip_comment(lun)
    read(lun,*) n
    op%n_three_body_mp = n

    allocate( &
         op%idx_three_body_mp(2,3,n), &
         op%v_three_body_mp(n) )
    
    do i = 1, n

       call skip_comment(lun)
       read(lun,*) ijk(1), ijk(2), ijk(3), v
       
       if ( ijk(1) > ijk(2) .or. ijk(2) > ijk(3)) &
            stop "ERROR in fn_three_body_mp"

       op%v_three_body_mp(i) = v
       do j = 1, 3
          if (ijk(j) <= n_jorb(1)) then
             op%idx_three_body_mp(1,j,i) = ijk(j)
             op%idx_three_body_mp(2,j,i) = 1
          else
             op%idx_three_body_mp(1,j,i) = ijk(j) - n_jorb(1)
             op%idx_three_body_mp(2,j,i) = 2
          end if
       end do

       if (myrank==0) write(*,'(a,6i3,f12.5)') 'three-body-monpole', &
            op%idx_three_body_mp(1,1,i), op%idx_three_body_mp(2,1,i), &
            op%idx_three_body_mp(1,2,i), op%idx_three_body_mp(2,2,i), &
            op%idx_three_body_mp(1,3,i), op%idx_three_body_mp(2,3,i), &
            op%v_three_body_mp(i)
       
    end do

    close(lun)

  end subroutine set_three_body_monopole
  
end module interaction
