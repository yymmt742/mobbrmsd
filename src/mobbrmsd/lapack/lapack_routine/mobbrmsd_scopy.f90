!| copies \( x \) to \( y \)
!
!  mobbrmsd_DCOPY copies a vector, \( x \), to a vector, \( y \),
!  and uses unrolled loops for increments equal to 1.
!
!  reference DCOPY is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!   November 2017
!
!  jack dongarra, linpack, 3/11/78.
!  modified 12/3/93, array(1) declarations changed to array(*)
!
pure subroutine mobbrmsd_SCOPY(N, SX, INCX, SY, INCY)
  implicit none
  integer, intent(in) :: N
!!  number of elements in input vector(s)
!!
  real(RK), intent(in)  :: SX(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in) :: INCX
!!  storage spacing between elements of SX
!!
  real(RK), intent(out) :: SY(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!!
  integer, intent(in) :: INCY
!!  storage spacing between elements of SY
!!
  integer :: I, IX, IY, M, MP1
  intrinsic :: MOD
!
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
!      code for both increments equal to 1
!      clean-up loop
!
    M = MOD(N, 7)
    if (M /= 0) then
      do I = 1, M
        SY(I) = SX(I)
      end do
      if (N < 7) return
    end if
    MP1 = M + 1
    do I = MP1, N, 7
      SY(I) = SX(I)
      SY(I + 1) = SX(I + 1)
      SY(I + 2) = SX(I + 2)
      SY(I + 3) = SX(I + 3)
      SY(I + 4) = SX(I + 4)
      SY(I + 5) = SX(I + 5)
      SY(I + 6) = SX(I + 6)
    end do
  else
!
!      code for unequal increments or equal increments
!        not equal to 1
!
    IX = 1
    IY = 1
    if (INCX < 0) IX = (-N + 1) * INCX + 1
    if (INCY < 0) IY = (-N + 1) * INCY + 1
    do I = 1, N
      SY(IY) = SX(IX)
      IX = IX + INCX
      IY = IY + INCY
    end do
  end if
  return
end

