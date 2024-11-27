!| mobbrmsd_SSWAP interchanges two vectors.
!  uses unrolled loops for increments equal to 1.
!
!  Level 1 Blas routine.
!  Reference SSWAP is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
pure subroutine mobbrmsd_SSWAP(N, SX, INCX, SY, INCY)
  implicit none
  integer, intent(in)     :: N
!!  Number of elements in input vector(s)
!!
  real(RK), intent(inout) :: SX(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)     :: INCX
!!  Storage spacing between elements of SX
!!
  real(RK), intent(inout) :: SY(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!!
  integer, intent(in)     :: INCY
!!  INCY is INTEGER
!!
  real(RK) :: STEMP
  integer :: I, IX, IY, M, MP1
  intrinsic :: MOD
!
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
! code for both increments equal to 1
!
! clean - up loop
!
    M = MOD(N, 3)
    if (M /= 0) then
      do I = 1, M
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
      end do
      if (N < 3) return
    end if
    MP1 = M + 1
    do I = MP1, N, 3
      STEMP = SX(I)
      SX(I) = SY(I)
      SY(I) = STEMP
      STEMP = SX(I + 1)
      SX(I + 1) = SY(I + 1)
      SY(I + 1) = STEMP
      STEMP = SX(I + 2)
      SX(I + 2) = SY(I + 2)
      SY(I + 2) = STEMP
    end do
  else
!
! code for unequal increments or equal increments not equal to 1
!
    IX = 1
    IY = 1
    if (INCX < 0) IX = (-N + 1) * INCX + 1
    if (INCY < 0) IY = (-N + 1) * INCY + 1
    do I = 1, N
      STEMP = SX(IX)
      SX(IX) = SY(IY)
      SY(IY) = STEMP
      IX = IX + INCX
      IY = IY + INCY
    end do
  end if
  return
end
