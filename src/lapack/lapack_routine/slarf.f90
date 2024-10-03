!> \brief \b SLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       REAL               TAU
!       ..
!       .. Array Arguments ..
!       REAL               C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARF applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
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
!> \date December 2016
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
pure subroutine SLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  character, intent(in) :: SIDE
  integer, intent(in)   :: INCV, LDC, M, N
  real, intent(in)      :: TAU
!..
!..Array Arguments..
  real, intent(in)    :: V(*)
  real, intent(inout) :: C(LDC, *)
  real, intent(out)   :: WORK(*)
!..
!
!  =====================================================================
!
!..Parameters..
  real, parameter :: ONE = 1.0E+0, ZERO = 0.0E+0
!..
!..Local Scalars..
  logical :: APPLYLEFT
  integer :: I, LASTV, LASTC
!..
  interface
! .. External Functions ..
    include 'lsame.h'
    include 'ilaslr.h'
    include 'ilaslc.h'
! .. External Subroutines ..
    include 'sgemv.h'
    include 'sger.h'
  end interface
!..
!..Executable Statements..
!
  APPLYLEFT = LSAME(SIDE, 'L')
  LASTV = 0
  LASTC = 0
  if (TAU /= ZERO) then
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
    if (APPLYLEFT) then
      LASTV = M
    else
      LASTV = N
    end if
    if (INCV > 0) then
      I = 1 + (LASTV - 1) * INCV
    else
      I = 1
    end if
!     Look for the last non-zero row in V.
    do while (LASTV > 0 .and. V(I) == ZERO)
      LASTV = LASTV - 1
      I = I - INCV
    end do
    if (APPLYLEFT) then
!     Scan for the last non-zero column in C(1:lastv,:).
      LASTC = ILASLC(LASTV, N, C, LDC)
    else
!     Scan for the last non-zero row in C(:,1:lastv).
      LASTC = ILASLR(M, LASTV, C, LDC)
    end if
  end if
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
  if (APPLYLEFT) then
    !
    !Form H * C
    !
    if (LASTV > 0) then
      !
      !w(1:lastc, 1): = C(1:lastv, 1:lastc)**T * v(1:lastv, 1)
      !
      call SGEMV('Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
          &      ZERO, WORK, 1)
      !
      !C(1:lastv, 1:lastc): = C(...) - v(1:lastv, 1) * w(1:lastc, 1)**T
      !
      call SGER(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC)
    end if
  else
    !
    !Form C * H
    !
    if (LASTV > 0) then
      !
      !w(1:lastc, 1): = C(1:lastc, 1:lastv) * v(1:lastv, 1)
      !
      call SGEMV('No transpose', LASTC, LASTV, ONE, C, LDC, &
          &      V, INCV, ZERO, WORK, 1)
      !
      !C(1:lastc, 1:lastv): = C(...) - w(1:lastc, 1) * v(1:lastv, 1)**T
      !
      call SGER(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC)
    end if
  end if
  return
!
!end of SLARF
!
end
