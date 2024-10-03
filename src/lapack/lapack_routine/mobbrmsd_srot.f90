!> \brief \b mobbrmsd_SROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SROT(N,SX,INCX,SY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       REAL C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    applies a plane rotation.
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
!> \param[in,out] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in,out] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL
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
!> \date November 2017
!
!> \ingroup single_blas_level1
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
pure subroutine mobbrmsd_SROT(N, SX, INCX, SY, INCY, C, S)
  implicit none
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
  real(RK), intent(in) :: C, S
  integer, intent(in) :: INCX, INCY, N
!..
!..Array Arguments..
  real(RK), intent(inout) :: SX(*), SY(*)
!..
!
!  =====================================================================
!
!..Local Scalars..
  real(RK) :: STEMP
  integer :: I, IX, IY
!..
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
! code for both increments equal to 1
!
    do I = 1, N
      STEMP = C * SX(I) + S * SY(I)
      SY(I) = C * SY(I) - S * SX(I)
      SX(I) = STEMP
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
      STEMP = C * SX(IX) + S * SY(IY)
      SY(IY) = C * SY(IY) - S * SX(IX)
      SX(IX) = STEMP
      IX = IX + INCX
      IY = IY + INCY
    end do
  end if
  return
end
