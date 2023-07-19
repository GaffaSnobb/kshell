! define block size at compile
#define NB02
#define NB03
#define NB04
! #define NB05
! #define NB06
! #define NB07
#define NB08
#define NB16
! #define NB24
#define NB32
! #define NB48
#define NB64
! ---



module bp_block
  use constant, only: kwf, kdim, kmbit, nmbit, nmlen, max_int4
  use model_space
  use partition, only: type_ptn_pn
  use partition, only: init_partition, type_ptn_pn, type_mbit, &
       bin_srch_nocc
  use wavefunction, only: dot_product_global, type_vec_p, wf_alloc_vec
  use operator_mscheme, only: opr_m, v_2b, idx_nocc2b, idx_gt, idx_2b
  use class_stopwatch
  use bridge_partitions
#ifdef MPI
  use mpi
#endif
  !$ use omp_lib  
  implicit none

  private
  public :: bp_operate_block

contains

  subroutine bp_operate_block(self, vl, op, vr)
    ! block version : vl = op * vr 
    type(type_bridge_partitions), intent(inout) :: self
    real(kwf), intent(out) :: vl(:,:)
    type(opr_m), intent(inout) :: op
    real(kwf), intent(in) :: vr(:,:)
    integer :: i, n, nb, k
    
    if (.not. allocated(op%mlmr)) stop "Error: call init_bp_operator"
    if (size(vl, 2) /= size(vr,2)) stop "error in bp_operate_block"

    n = size(vl, 2)
    i = 1
    
    do while ( n-i+1 > 0 )

       k = n-i+1

       if (.false.) then 
#ifdef NB64
       elseif (64 <= k) then 
          nb = 64
          call bp_operate_block64(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB48
       elseif (48 <= k) then 
          nb = 48
          call bp_operate_block48(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB32
       elseif (32 <= k) then 
          nb = 32
          call bp_operate_block32(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB24
       elseif (24 <= k) then 
          nb = 24
          call bp_operate_block24(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB16
       elseif (16 <= k) then 
          nb = 16
          call bp_operate_block16(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB08
       elseif (8 <= k) then 
          nb = 8
          call bp_operate_block8(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB07
       elseif (7 <= k) then 
          nb = 7
          call bp_operate_block7(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB06
       elseif (6 <= k) then 
          nb = 6
          call bp_operate_block6(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB05
       elseif (5 <= k) then 
          nb = 5
          call bp_operate_block5(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB04
       elseif (4 <= k) then 
          nb = 4
          call bp_operate_block4(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB03
       elseif (3 <= k) then 
          nb = 3
          call bp_operate_block3(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
#ifdef NB02
       elseif (2 <= k) then 
          nb = 2
          call bp_operate_block2(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
#endif
       else
          nb = 1
          call bp_operate_block1(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       end if

       ! if (myrank==0) write(*,'(a,i3,a)') 'operate block ',nb,' called'
       i = i + nb

    end do

  end subroutine bp_operate_block

#ifdef NB64
#define NBLOCK 64
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB48
#define NBLOCK 48
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB32
#define NBLOCK 32
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB24
#define NBLOCK 24
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB16
#define NBLOCK 16
#include "bp_block_inc.F90"
#undef NBLOCK
#endif 

#ifdef NB08
#define NBLOCK 8
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB07
#define NBLOCK 7
#include "bp_block_inc.F90"
#undef NBLOCK
#endif


#ifdef NB06
#define NBLOCK 6
#include "bp_block_inc.F90"
#undef NBLOCK
#endif

#ifdef NB05
#define NBLOCK 5
#include "bp_block_inc.F90"
#undef NBLOCK
#endif 

#ifdef NB04
#define NBLOCK 4
#include "bp_block_inc.F90"
#undef NBLOCK
#endif 

#ifdef NB03
#define NBLOCK 3
#include "bp_block_inc.F90"
#undef NBLOCK
#endif 


#ifdef NB02
#define NBLOCK 2
#include "bp_block_inc.F90"
#undef NBLOCK
#endif 

#define NBLOCK 1
#include "bp_block_inc.F90"
#undef NBLOCK


end module bp_block
