!| applies a plane rotation.
!
!  Reference SROT is provided by [netlib.org](http://www.netlib.org/lapack/).
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
pure subroutine mobbrmsd_DROT(N, DX, INCX, DY, INCY, C, S)
  implicit none
  integer, intent(in)     :: N
!!         number of elements in input vector(s)
!!
  real(RK), intent(inout) :: DX(*)
!!          DX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)     :: INCX
!!         storage spacing between elements of SX
!!
  real(RK), intent(inout) :: DY(*)
!!          DY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!!
  integer, intent(in)     :: INCY
!!         storage spacing between elements of SY
!!
  real(RK), intent(in) :: C
!!          C is REAL
!!
  real(RK), intent(in) :: S
!!          S is REAL
!!
  real(RK)               :: DTEMP
  integer                :: I, IX, IY
!
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
!       code for both increments equal to 1
!
    do I = 1, N
      DTEMP = C * DX(I) + S * DY(I)
      DY(I) = C * DY(I) - S * DX(I)
      DX(I) = DTEMP
    end do
  else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
    IX = 1
    IY = 1
    if (INCX < 0) IX = (-N + 1) * INCX + 1
    if (INCY < 0) IY = (-N + 1) * INCY + 1
    do I = 1, N
      DTEMP = C * DX(IX) + S * DY(IY)
      DY(IY) = C * DY(IY) - S * DX(IX)
      DX(IX) = DTEMP
      IX = IX + INCX
      IY = IY + INCY
    end do
  end if
  return
!
!     End of mobbrmsd_DROT
!
end subroutine mobbrmsd_DROT

