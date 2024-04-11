program main
  use blas_lapack_interface, only: D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use driver
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=2
  integer            :: itest
!
  call u%init('test python_driver')
!
  do itest=1,NTEST
    call test1()
  enddo
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: threads = 2
    integer(IK), parameter   :: n = 52, M = 8, f = 12, g = 8
    real(RK)                 :: X(d, n, m), Y(d, n, m, threads)
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
          X(:, :, (i + 1) * 2 + (j + 1) + (k + 1) / 2 + 1) = sample(d, n, [i, j, k])
        end do
      end do
    end do
!
    do l = 1, threads
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            Y(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l) = sample(d, n, [i, j, k])
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
    call batch_run(D, n * M, threads, n_header, n_int, n_float, X, Y, &
                   999.0_RK, 0.0_RK, -1, .FALSE., h, si, sr)
    print'(10(I6))',h
    do i = 1, threads
      print'(*(I4))', si(:, i)
    end do
    do i = 1, threads
      print'(*(F12.3))',sr(:, i)
    end do
    call clear_molecule()
!
  end subroutine test1
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
