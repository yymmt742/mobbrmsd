program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use driver
  use mod_unittest
  implicit none
  type(unittest) :: u
! integer, parameter :: NTEST=25
! integer            :: itest
!
  call u%init('test symRMSD')
  call test1()
!
! do itest = 1, NTEST
!   call test1()
! end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter         :: d = 3
    integer, parameter         :: s = 1
    integer, parameter         :: threads = 10
    integer, parameter         :: m = 50, n = 18, f = 12, g = 18
    real(RK)                   :: X(d, m, n), Y(d, m, n, threads)
    real(RK)                   :: res(threads)
    integer                    :: i, j, k, l
!
    call add_molecule(m, n, s, [(i, i=1, 0)])
    call setup()
!
    do k = -1, 1, 2
      do j = -1, 1
        do i = -1, 1
          X(:, :, (i + 1) * 6 + (j + 1) * 2 + (k + 1) / 2 + 1) = sample(d, m, [i, j, k])
        end do
      end do
    end do
!
    do l = 1, threads
      do k = -1, 1, 2
        do j = -1, 1
          do i = -1, 1
            Y(:, :, (k + 1) / 2 * 9 + (j + 1) * 3 + (i + 1) + 1, l) = sample(d, m, [i, j, k])
          end do
        end do
      end do
    end do
!
    call run(X, Y, threads, res)
    print*,res
!
  end subroutine test1
!
  function sample(d, n, com) result(res)
    integer, intent(in)  :: d, n, com(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) + 5 * com
    enddo
  end function sample
!
end program main
