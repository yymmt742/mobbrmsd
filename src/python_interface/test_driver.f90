program main
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
  do itest = 1, NTEST
    call test1()
  end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer(IK), parameter :: d = 3
    integer(IK), parameter :: s = 1
    integer(IK), parameter :: threads = 2
    integer(IK), parameter :: m = 52, n = 8, f = 12, g = 8
    real(RK)               :: X(d, m, n), Y(d, m, n, threads)
    real(RK)               :: res(threads), ratio(threads)
    integer(IK)            :: nsearch(threads)
    integer(IK)            :: i, j, k, l
!
    call add_molecule(m, n, s, [[5, 6, 7, 8, 1, 2, 3, 4], [(i, i=9, m)]])
    call setup()
!
    do k = -1, 1, 2
      do j = -1, 1, 2
        do i = -1, 1, 2
          X(:, :, (i + 1) * 2 + (j + 1) + (k + 1) / 2 + 1) = sample(d, m, [i, j, k])
        end do
      end do
    end do
!
    do l = 1, threads
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            Y(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l) = sample(d, m, [i, j, k])
          end do
        end do
      end do
      Y(:, :, :, l) = 0.8 * X + 0.2 * Y(:, :, :, l)
    end do
!
    call run(X, Y, threads, res, ratio, nsearch)
    print'(5f9.3)', res
    print'(5f9.3)', ratio
    call clear()
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
