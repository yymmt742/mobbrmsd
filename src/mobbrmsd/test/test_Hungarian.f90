program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_Hungarian
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer(IK), parameter :: NTEST = 25
  integer(IK)            :: itest
#ifdef USE_REAL32
  integer(IK), parameter :: place = 3
#else
  integer(IK), parameter :: place = 7
#endif
!
  call u%init('test Hungarian square')
  do itest = 1, NTEST
    call test1()
  end do
  call u%init('test Hungarian general')
  do itest = 1, NTEST
    call test2()
  end do
  call u%init('test Hungarian general transpose')
  do itest = 1, NTEST
    call test3()
  end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer(IK), parameter :: N_TEST = 10
    integer(IK), parameter :: d = 3
    integer(IK), parameter :: n = 3
    real(RK)               :: C(n, n), dum(1)
    real(RK)               :: minsp
!
    call RANDOM_NUMBER(C)
    C = C - 5
    minsp = 999
    minsp = MIN(minsp, SP(n, n, [1_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, SP(n, n, [1_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, SP(n, n, [2_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, SP(n, n, [2_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, SP(n, n, [3_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, SP(n, n, [3_IK, 2_IK, 1_IK], C))
    call Hungarian(-n, -n, C, dum)
    block
      real(RK) :: W(NINT(dum(1)))
      call Hungarian(n, n, C, W)
      call u%assert_almost_equal(W(1), minsp, 'W1 = HV   ', place=place)
    end block
!
  end subroutine test1
!
  subroutine test2()
    integer(IK), parameter :: N_TEST = 10
    integer(IK), parameter :: d = 3
    integer(IK), parameter :: m = 3
    integer(IK), parameter :: n = 5
    real(RK)               :: C(m, n), dum(1)
    real(RK)               :: minsp
!
    call RANDOM_NUMBER(C)
    minsp = 999
    minsp = MIN(minsp, SP(m, n, [1_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 3_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 3_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 4_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 3_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 3_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 4_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 4_IK, 5_IK], C))

    minsp = MIN(minsp, SP(m, n, [1_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 4_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 5_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [1_IK, 5_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 4_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 5_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 5_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 5_IK, 4_IK], C))

    minsp = MIN(minsp, SP(m, n, [2_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 3_IK, 5_IK], C))

    minsp = MIN(minsp, SP(m, n, [2_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [2_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [3_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 5_IK, 3_IK], C))

    minsp = MIN(minsp, SP(m, n, [3_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 3_IK, 4_IK], C))

    minsp = MIN(minsp, SP(m, n, [3_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, SP(m, n, [4_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, SP(m, n, [5_IK, 4_IK, 3_IK], C))
!
    call Hungarian(m, -n, C, dum)
    block
      real(RK) :: W(NINT(dum(1)))
      call Hungarian(m, n, C, W)
      call u%assert_almost_equal(W(1), minsp, 'W1 = HV   ', place=place)
    end block
!
!
  end subroutine test2
!
  subroutine test3()
    integer(IK), parameter :: N_TEST = 10
    integer(IK), parameter :: d = 3
    integer(IK), parameter :: m = 5
    integer(IK), parameter :: n = 3
    real(RK)               :: C(m, n), dum(1)
    real(RK)               :: minsp
!
    call RANDOM_NUMBER(C)
    minsp = 999
    minsp = MIN(minsp, PS(m, n, [1_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 3_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 3_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 4_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 3_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 3_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 4_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 4_IK, 5_IK], C))
!
    minsp = MIN(minsp, PS(m, n, [1_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 4_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 5_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [1_IK, 5_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 4_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 5_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 5_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 5_IK, 4_IK], C))
!
    minsp = MIN(minsp, PS(m, n, [2_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 1_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 2_IK, 5_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 3_IK, 5_IK], C))
!
    minsp = MIN(minsp, PS(m, n, [2_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [2_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 5_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [3_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 5_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 5_IK, 3_IK], C))
!
    minsp = MIN(minsp, PS(m, n, [3_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 1_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 1_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 1_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 2_IK, 3_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 2_IK, 4_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 3_IK, 4_IK], C))
!
    minsp = MIN(minsp, PS(m, n, [3_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 2_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 3_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 4_IK, 1_IK], C))
    minsp = MIN(minsp, PS(m, n, [4_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 3_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 4_IK, 2_IK], C))
    minsp = MIN(minsp, PS(m, n, [5_IK, 4_IK, 3_IK], C))
!
    call Hungarian(m, -n, C, dum)
    block
      real(RK) :: W(NINT(dum(1)))
      call Hungarian(m, n, C, W)
      call u%assert_almost_equal(W(1), minsp, 'W1 = HV   ', place=place)
    end block
!
  end subroutine test3
!
  pure function SP(m, n, ix, C) result(res)
    integer, intent(in)  :: m, n, ix(m)
    real(RK), intent(in) :: C(m, n)
    integer              :: i
    real(RK)             :: res
    res = 0D0
    do i = 1, m
      res = res + C(i, ix(i))
    end do
  end function SP
!
  pure function PS(m, n, ix, C) result(res)
    integer, intent(in)  :: m, n, ix(n)
    real(RK), intent(in) :: C(m, n)
    integer              :: i
    real(RK)             :: res
    res = 0D0
    do i = 1, n
      res = res + C(ix(i), i)
    end do
  end function PS
!
end program main
