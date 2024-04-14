program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO, RHUGE
  use driver
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=2
  integer            :: itest
  integer            :: n_dim
!
  call u%init('test python_driver')
!
  call n_dims(n_dim)
  print*,n_dim
!
  do itest=1,NTEST
    call test1(n_dim)
  enddo
!
  call test2(n_dim)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(n_dim)
    integer(IK), intent(in)  :: n_dim
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: threads = 2
    integer(IK), parameter   :: n = 52, M = 8
    real(RK)                 :: X(n_dim, n, m), Y(n_dim, n, m, threads)
    integer(IK)              :: n_header, n_int, n_float
    integer(IK), allocatable :: h(:), si(:, :)
    real(RK), allocatable    :: sr(:, :)
    integer(IK)              :: i, j, k, l
!
    call add_molecule(n, M, s, [[5, 6, 7, 8, 1, 2, 3, 4], [(i, i=9, n)]])
!
    do k = -1, 1, 2
      do j = -1, 1, 2
        do i = -1, 1, 2
          X(:, :, (i + 1) * 2 + (j + 1) + (k + 1) / 2 + 1) = sample(n_dim, n, [i, j, k])
        end do
      end do
    end do
!
    do l = 1, threads
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            Y(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l) = sample(n_dim, n, [i, j, k])
          end do
        end do
      end do
      Y(:, :, :, l) = 0.2 * X + 0.8 * Y(:, :, :, l)
    end do
!
    call state_vector_lengthes(n_header, n_int, n_float)
    print*, n_header, n_int, n_float
!
    allocate(h(n_header))
    allocate(si(n_int, threads))
    allocate(sr(n_float, threads))
!
    call batch_run(n_dim, n * M, threads, n_header, n_int, n_float, X, Y, &
                   999.0_RK, 0.0_RK, -1, .FALSE., h, si, sr)
    print'(10(I6))',h
    do i = 1, threads
      print'(*(I4))', si(:, i)
    end do
    do i = 1, threads
      print'(5(F12.3))',sr(:, i)
    end do
    call clear_molecule()
!
  end subroutine test1
!
  subroutine test2(n_dim)
    integer(IK), intent(in)  :: n_dim
    integer(IK), parameter   :: n_apm = 8, n_mol = 8, n_target = 1000
    integer(IK), parameter   :: n_atom = n_apm * n_mol
    real(RK)                 :: X(n_dim, n_apm, n_mol), Y(n_dim, n_apm, n_mol, n_target)
    integer(IK)              :: n_job, n_mem
    integer(IK)              :: nn_index
    real(RK)                 :: nn_value
    real(RK), allocatable    :: W(:)
!
    call RANDOM_NUMBER(X)
    call RANDOM_NUMBER(Y)
!   Y(:, :, :, 11) = 0.9 * Y(:, :, :, 11) + 0.1 * X
!   Y(:, :, :, 22) = 0.7 * Y(:, :, :, 22) + 0.3 * X
!   Y(:, :, :, 56) = 0.6 * Y(:, :, :, 56) + 0.4 * X
    Y(:, :, :, 83) = 0.4 * Y(:, :, :, 83) + 0.6 * X
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
!
    call n_jobs(n_job)
    call workmemory_length(n_mem)
    print*, n_job, n_mem
    allocate (w(n_job * n_mem))
!
    call nearest_neighbor(n_dim, n_atom, n_target, n_mem, n_job,&
 &                        X, Y, W, RHUGE, ZERO, -1, &
 &                        nn_index, nn_value)
!
    print*,nn_index, nn_value
!
    call clear_molecule()
!
  end subroutine test2
!
  function sample(d, n, com) result(res)
    integer, intent(in)  :: d, n, com(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) + 8 * com
    enddo
  end function sample
!
end program main
