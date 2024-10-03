!> \brief \b DLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dlarf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dlarf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dlarf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       real(RK)           ::   TAU
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARF applies a real elementary reflector H to a real m by n matrix
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
!>          V is real(RK)           :: array, dimension
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
!>          TAU is real(RK)           ::
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is real(RK)           :: array, dimension (LDC,N)
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
!>          WORK is real(RK)           :: array, dimension
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
pure subroutine DLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: SIDE
  integer, intent(in)      :: INCV, LDC, M, N
  real(RK), intent(in)     :: TAU
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: V(*)
  real(RK), intent(inout)  :: C(LDC, *)
  real(RK), intent(out)    :: WORK(LDC, *)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  logical                 :: APPLYLEFT
  integer                 :: I, LASTV, LASTC
!     ..
!     .. Parameters ..
! real(RK), parameter      :: ONE = 1.0_RK
! real(RK), parameter      :: ZERO = 0.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dgemv.h'
!   include 'dger.h'
!     .. External Functions ..
!   include 'lsame.h'
!   include 'iladlr.h'
!   include 'iladlc.h'
! end interface
!     ..
!     .. Executable Statements ..
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
      LASTC = ILADLC(LASTV, N, C, LDC)
    else
!     Scan for the last non-zero row in C(:,1:lastv).
      LASTC = ILADLR(M, LASTV, C, LDC)
    end if
  end if
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
  if (APPLYLEFT) then
!
!        Form  H * C
!
    if (LASTV > 0) then
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
      call DGEMV('Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,&
     &     ZERO, WORK, 1)
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
!
      call DGER(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC)
    end if
  else
!
!        Form  C * H
!
    if (LASTV > 0) then
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
      call DGEMV('No transpose', LASTC, LASTV, ONE, C, LDC,&
     &     V, INCV, ZERO, WORK, 1)
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
!
      call DGER(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC)
    end if
  end if
  return
!
!     End of DLARF
!
end subroutine DLARF
