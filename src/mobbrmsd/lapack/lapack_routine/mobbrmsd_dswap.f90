!| mobbrmsd_DSWAP interchanges two vectors.
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
pure subroutine mobbrmsd_DSWAP(N, DX, INCX, DY, INCY)
  implicit none
  integer, intent(in)     :: N
!!  Number of elements in input vector(s)
!!
  real(RK), intent(inout) :: DX(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)     :: INCX
!!  Storage spacing between elements of SX
!!
  real(RK), intent(inout) :: DY(*)
!!  REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!!
  integer, intent(in)     :: INCY
!!  INCY is INTEGER
!!
  integer                :: I, IX, IY, M, MP1
  intrinsic              :: MOD
!
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
    M = MOD(N, 3)
    if (M /= 0) then
      do concurrent(I=1:M)
        block
          real(RK) :: DTEMP
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
        end block
      end do
      if (N < 3) return
    end if
    MP1 = M + 1
    do concurrent(I=MP1:N:3)
      block
        real(RK) :: DTEMP
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
      end block
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
!
    do concurrent(I=0:N - 1)
      block
        real(RK) :: DTEMP
        integer  :: IIX, IIY
        IIX = IX + I * INCX
        IIY = IY + I * INCY
        DTEMP = DX(IIX)
        DX(IIX) = DY(IIY)
        DY(IIX) = DTEMP
      end block
    end do
  end if
  return
!
!     End of mobbrmsd_DSWAP
!
end subroutine mobbrmsd_DSWAP

