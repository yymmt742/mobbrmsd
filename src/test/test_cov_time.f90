program main
  use, intrinsic :: ISO_FORTRAN_ENV
  implicit none
!
  integer, parameter :: N = 10000000
  integer, parameter :: NRUN = 10
  integer, parameter :: P = 6
  integer, parameter :: NI(P) = [10, 100, 1000, 10000, 100000, 100000]
  integer, parameter :: MI(P) = [10000, 1000, 100, 10, 1, 1]
  integer, parameter :: DI(P) = [10, 10, 10, 10, 10, 10]
  real(REAL64)       :: RET, T(P), DX(N), DY(N)
  integer            :: J, K
!
  interface
    pure subroutine cov(d, n, x, y, res)
      use, intrinsic :: ISO_FORTRAN_ENV
      integer, intent(in)         :: d, n
      real(REAL64), intent(in)    :: x(d, n), y(d, n)
      real(REAL64), intent(inout) :: res(d, d)
    end subroutine cov
  end interface
!
  T   = 0.0_REAL64
  RET = 0.0_REAL64
!
  do J = 1, NRUN
!
    call RANDOM_NUMBER(DX)
    call RANDOM_NUMBER(DY)
!
    do K = 1, P
!
      call wallclock_time(NI(K), DI(K), MI(K), DX, DY, RET, T(K))
!
    end do
  end do
!
  print'(F24.12,*(F12.6))', RET, T
!
contains
!
  subroutine wallclock_time(N, D, M, DX, DY, RET, time)
    integer, intent(in)         :: N, D, M
    real(REAL64), intent(in)    :: DX(N, D, M), DY(N, D, M)
    real(REAL64), intent(inout) :: RET, time
    real(REAL64)                :: DTMP(D,D)
    integer                     :: t1, t2, t_rate, t_max, diff
    integer                     :: I
!
    call SYSTEM_CLOCK(t1)
!
    RET = 0.0_REAL64
!
    do I = 1, M
      DTMP = 0.0_REAL64
      call cov(d, n, DX(:,:,I), DY(:,:,I), DTMP)
      RET = RET + SUM(DTMP)
    end do
!
    RET = RET / (M * N)
!
    call SYSTEM_CLOCK(t2, t_rate, t_max)
    if (t2 < t1) then
      diff = (t_max - t1) + t2 + 1
    else
      diff = t2 - t1
    end if
!
    time = time + diff / real(t_rate, REAL64)
!
  end subroutine wallclock_time
!
end program main
