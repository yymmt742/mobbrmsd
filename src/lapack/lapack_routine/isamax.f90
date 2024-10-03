!> \brief\b ISAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!        http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!   INTEGER FUNCTION ISAMAX(N,SX,INCX)
!
!   .. Scalar Arguments ..
!   INTEGER INCX,N
!   ..
!   .. Array Arguments ..
!   REAL SX(*)
!   ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ISAMAX finds the index of the first element having maximum absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!> N is integer
!> number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!> SX is real array, dimension(1 + (N - 1)! ABS(INCX))
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!> INCX is integer
!> storage spacing between elements of SX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ.of Tennessee
!> \author Univ.of California Berkeley
!> \author Univ.of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2017
!
!> \ingroup aux_blas
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!> jack dongarra, linpack, 3 / 11 / 78.
!> modified 3 / 93 to return if incx <= 0.
!> modified 12 / 3 / 93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
pure function ISAMAX(N, SX, INCX)
  implicit none
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! November 2017
!
! .. Scalar Arguments ..
  integer, intent(in) :: INCX, N
! ..
! .. Array Arguments ..
  real, intent(in)    :: SX(*)
! ..
  integer             :: ISAMAX
!
!  =====================================================================
!
! .. Local Scalars ..
  real                :: SMAX
  integer             :: I, IX
! ..
! .. Intrinsic Functions ..
  intrinsic ABS
! ..
  ISAMAX = 0
  if (N < 1 .or. INCX <= 0) return
  ISAMAX = 1
  if (N == 1) return
  if (INCX == 1) then
!
!    code for increment equal to 1
!
    SMAX = ABS(SX(1))
    do I = 2, N
      if (ABS(SX(I)) > SMAX) then
        ISAMAX = I
        SMAX = ABS(SX(I))
      end if
    end do
  else
!
!    code for increment not equal to 1
!
    IX = 1
    SMAX = ABS(SX(1))
    IX = IX + INCX
    do I = 2, N
      if (ABS(SX(IX)) > SMAX) then
        ISAMAX = I
        SMAX = ABS(SX(IX))
      end if
      IX = IX + INCX
    end do
  end if
  return
end
