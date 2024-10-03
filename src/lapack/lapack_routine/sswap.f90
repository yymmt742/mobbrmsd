!> \brief \b SSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
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
!>    SSWAP interchanges two vectors.
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
pure subroutine SSWAP(N, SX, INCX, SY, INCY)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
  integer, intent(in) :: INCX, INCY, N
!..
!..Array Arguments..
  real, intent(inout) ::  SX(*), SY(*)
!..
!
!  =====================================================================
!
!..Local Scalars..
  real    :: STEMP
  integer :: I, IX, IY, M, MP1
!..
!..intrinsic Functions..
  intrinsic :: MOD
!..
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
    !
    !code for both increments equal to 1
    !
    !
    !clean - up loop
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
    !code for unequal increments or equal increments not equal
    !to 1
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
