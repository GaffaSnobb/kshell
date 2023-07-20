module tests
    implicit none
    public :: dot_product_global

    type type_mbit 
        integer :: n                     ! number of states in m-scheme of id
        integer(8), allocatable :: mbit(:) ! bit repr. of m-scheme states
    end type type_mbit

    type type_ptn_id
        type(type_mbit), allocatable :: mz(:) ! bit of m-state for Jz=mm
        integer :: min_m, max_m, iprty ! 2*max of Jz, J &  parity
    end type type_ptn_id

    type type_ptn  ! partition for proton (or neutron)
        integer :: n_id   ! number of partition id's
        type(type_ptn_id), allocatable :: id(:) 
        !  nocc(korb, id)  number of occupation in orbits
        integer, allocatable :: nocc(:,:)  
        integer(8) :: n_mbit ! total number of M-scheme bits
        integer :: n_targeted = 0
    end type type_ptn

    type type_ptn_pn ! self 
        integer :: n_ferm(2), iprty, mtotal, max_jj
        type(type_ptn), pointer :: pn(:) => null()  ! pn(2)
        integer :: n_pidpnM
        ! pidpnM_pid(:,id) = idp, idn, Mp*2 (deployed order)
        integer, allocatable :: pidpnM_pid(:,:)     
        ! sorted partion id for binary search
        integer, allocatable :: pidpnM_pid_srt(:,:) 
        ! index trans. of sorted to that of deployed, and reversed
        integer, allocatable :: pid_srt2dpl(:) 
        integer, allocatable :: pid_dpl2srt(:) ! index reversed trans.
        integer(8) :: ndim, max_ndim_pid, max_ndim_pid_pn(2)
        integer(8), allocatable :: ndim_pid(:), ndim_srt(:), ndim_srt_acc(:)
        integer :: n_nocc
        integer, allocatable :: srt2nocc(:) ! # partition ID => nocc ID in ptn file
        ! for MPI
        integer :: idl_start, idl_end ! , ntask
        integer, allocatable :: rank2ntask(:), pid2rank(:)
        integer(8) :: local_dim, max_local_dim
        integer(8), allocatable :: local_dim_acc(:), local_dim_acc_start(:)
    end type type_ptn_pn
    
    type type_vec_p
        integer :: jj = -1, tt = -1
        real(8) :: eval = 0.d0
        real(8), pointer :: p(:) => null()
        type(type_ptn_pn), pointer :: ptn => null()
    end type type_vec_p
    
    contains
    function dot_product_global(self, rwf) result (r)
        ! <self|rwf>  self and rwf are global vectors
        !  "r" should be real(8), not real(kwf)
        type(type_vec_p), intent(in) :: self, rwf
        real(8) :: r
        integer(8) :: mq

        r = 0.d0
        !$omp parallel do private(mq) reduction (+: r)
        do mq = 1, size(self%p, kind=8)
            r = r + self%p(mq) * rwf%p(mq)
        end do
end function dot_product_global

end module tests

program test_dot_product_global
    use tests, only: type_vec_p, dot_product_global
    implicit none
  
    integer :: ierr, my_rank, num_procs, test_case
    real(8), parameter :: eps = 1.0e-6
    type(type_vec_p) :: vec1, vec2
    real(8) :: result, expected_result
    integer :: i
  
    ! Test cases
    real(8), dimension(6, 4) :: test_vectors = reshape([ &
         1.0, 2.0, 3.0, 4.0, &
         5.0, 6.0, 7.0, 8.0, &
         1.0, 0.0, 0.0, 0.0, &
         0.0, 1.0, 0.0, 0.0, &
         2.0, 3.0, 1.0, 1.0, &
         4.0, 2.0, 2.0, 3.0], [6, 4])
  
    ! Initialize MPI
    ! call mpi_init(ierr)
    ! call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    ! call mpi_comm_size(mpi_comm_world, num_procs, ierr)
  
    ! Allocate memory for vec1 and vec2
    allocate(vec1%p(4), vec2%p(4))
  
    do test_case = 1, 5
      ! Initialize vec1 and vec2 with test case values
      vec1%p = test_vectors(2 * test_case - 1, 1:4)
      vec2%p = test_vectors(2 * test_case, 1:4)
  
      ! Calculate the dot product of vec1 and vec2 using dot_product_global
      result = dot_product_global(vec1, vec2)
  
      ! Calculate the expected result
      expected_result = sum(vec1%p * vec2%p)
  
      ! Check if the result is close to the expected result
      if (abs(result - expected_result) > eps) then
        write(*,*) "Test", test_case, "FAILED on process", my_rank, ": Expected", expected_result, "but got", result
      else
        write(*,*) "Test", test_case, "PASSED on process", my_rank
      end if
    end do
  
    ! Finalize MPI
    ! call mpi_finalize(ierr)
  
  end program test_dot_product_global

!   program benchmark_dot_product_global
!     use omp_lib
!     use tests, only: type_vec_p, dot_product_global
!     implicit none
  
!     integer :: ierr, my_rank, num_procs, num_threads
!     type(type_vec_p) :: vec1, vec2
!     real(8) :: result
!     integer :: vec_size, i, max_threads
!     real(8) :: start_time, end_time
  
!     ! Set the size of the vectors
!     vec_size = 100000000
  
!     ! Initialize MPI
!     ! call mpi_init(ierr)
!     ! call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
!     ! call mpi_comm_size(mpi_comm_world, num_procs, ierr)
  
!     ! Allocate memory for vec1 and vec2
!     start_time = omp_get_wtime()
!     allocate(vec1%p(vec_size), vec2%p(vec_size))
!     end_time = omp_get_wtime()
!     if (my_rank == 0) write(*,*) "Execution time for allocate: ", end_time - start_time, "seconds"
  
!     ! Initialize vec1 and vec2 with random values
!     start_time = omp_get_wtime()
!     call random_number(vec1%p)
!     call random_number(vec2%p)
!     end_time = omp_get_wtime()
!     if (my_rank == 0) write(*,*) "Execution time for random_number: ", end_time - start_time, "seconds"
  
!     ! Determine the maximum number of threads
!     max_threads = omp_get_max_threads()
  
!     ! Loop over the number of threads
!     do num_threads = 1, max_threads, 2
!       ! Set the number of threads for this test
!       call omp_set_num_threads(num_threads)
  
!       ! Perform the dot product calculation and measure the execution time
!       if (my_rank == 0) write(*,*) "Running test with", num_threads, "threads"
!       start_time = omp_get_wtime()
!       result = dot_product_global(vec1, vec2)
!       end_time = omp_get_wtime()
  
!       ! Output the execution time for this test
!       if (my_rank == 0) write(*,*) "Execution time for", num_threads, "threads:", end_time - start_time, "seconds"
!     end do
  
!     ! Finalize MPI
!     ! call mpi_finalize(ierr)
  
!   end program benchmark_dot_product_global
  