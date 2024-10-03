!> \brief \b DCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
pure subroutine DCOPY(N, DX, INCX, DY, INCY)
! use LA_CONSTANTS, only: RK => DP
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)   :: INCX, INCY, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)  :: DX(*)
  real(RK), intent(out) :: DY(*)
!     ..
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
!        code for both increments equal to 1
!
!        clean-up loop
!
    block
      integer   :: I, M, MP1
      intrinsic :: MOD
!
      M = MOD(N, 7)
      if (M /= 0) then
        do concurrent(I=1:M)
          DY(I) = DX(I)
        end do
        if (N < 7) return
      end if
!
      MP1 = M + 1
!
      do concurrent(I=MP1:N:7)
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
      end do
!
    end block
!
  elseif ((INCX == 0) .and. (INCY == 0)) then
!
    DY(1) = DX(1)
!
  elseif (INCX == 0) then
!
    block
      integer   :: I, LY, UY
!
      if (INCY < 0) then
        LY = (-N + 1) * INCY + 1
        UY = 1
      else
        LY = 1
        UY = (N - 1) * INCY + 1
      end if
!
      do concurrent(I=LY:UY:INCY)
        DY(I) = DX(1)
      end do
!
    end block
!
  elseif (INCY == 0) then
!
    if (INCX < 0) then
      DY(1) = DX(1)
    else
      DY(1) = DX((N - 1) * INCY + 1)
    end if
!
  else
!
!   code for unequal increments or equal increments
!     not equal to 1
!
    block
      integer   :: I, IX, IY
!
      IX = 1
      IY = 1
      if (INCX < 0) IX = (-N + 1) * INCX + 1
      if (INCY < 0) IY = (-N + 1) * INCY + 1
!
      do concurrent(I=0:N - 1)
        block
          integer :: IIX, IIY
          IIX = IX + I * INCX
          IIY = IY + I * INCY
          DY(IIY) = DX(IIX)
        end block
      end do
!
    end block
!
  end if
!
!     End of DCOPY
!
end subroutine DCOPY
