!> \brief \b DROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,S
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
!>    DROT applies a plane rotation.
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
!> \param[in,out] DX
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
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
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
pure subroutine DROT(N, DX, INCX, DY, INCY, C, S)
! use LA_CONSTANTS, only: RK => dp
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in)    :: C, S
  integer, intent(in)     :: INCX, INCY, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: DX(*), DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  real(RK)               :: DTEMP
  integer                :: I, IX, IY
!     ..
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
!     End of DROT
!
end subroutine DROT

