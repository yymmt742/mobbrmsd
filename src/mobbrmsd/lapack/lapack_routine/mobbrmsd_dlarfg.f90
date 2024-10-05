!| mobbrmsd_DLARFG generates an elementary reflector (Householder matrix).
!  mobbrmsd_DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H**T * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v**T ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLARFG(N, ALPHA, X, INCX, TAU)
  implicit none
  integer, intent(in)     :: INCX
!!          The increment between elements of X. INCX > 0.
  integer, intent(in)     :: N
!!          The order of the elementary reflector.
  real(RK), intent(inout) :: ALPHA
!!          On entry, the value alpha. <br>
!!          On exit, it is overwritten with the value beta. <br>
  real(RK), intent(out)   :: TAU
!!          The value tau.
  real(RK), intent(inout) :: X(*)
!!          X is real(RK)           :: array, dimension (1+(N-2)*abs(INCX))
!!          On entry, the vector x.
!!          On exit, it is overwritten with the vector v.
  real(RK)             :: BETA, RSAFMN, SAFMIN, XNORM
  integer              :: J, KNT
  intrinsic            :: ABS, SIGN
!     ..
!     .. Parameters ..
! real(RK), parameter   :: ZERO = 0.0_RK
! real(RK), parameter   :: ONE = 1.0_RK
!     ..
!     ..
! interface
!     .. External Functions ..
!   include 'dlamch.h'
!   include 'dlapy2.h'
!   include 'dnrm2.h'
!     .. External Subroutines ..
!   include 'dscal.h'
! end interface
!     ..
!     .. Executable Statements ..
!
  if (N <= 1) then
    TAU = ZERO
    return
  end if
!
  XNORM = mobbrmsd_DNRM2(N - 1, X, INCX)
!
  if (XNORM == ZERO) then
!
! H  =  I
!
    TAU = ZERO
  else
!
!        general case
!
    BETA = -SIGN(mobbrmsd_DLAPY2(ALPHA, XNORM), ALPHA)
    SAFMIN = mobbrmsd_DLAMCH('S') / mobbrmsd_DLAMCH('E')
    KNT = 0
    if (ABS(BETA) < SAFMIN) then
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
      RSAFMN = ONE / SAFMIN
10    continue
      KNT = KNT + 1
      call mobbrmsd_DSCAL(N - 1, RSAFMN, X, INCX)
      BETA = BETA * RSAFMN
      ALPHA = ALPHA * RSAFMN
      if ((ABS(BETA) < SAFMIN) .and. (KNT < 20)) GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
      XNORM = mobbrmsd_DNRM2(N - 1, X, INCX)
      BETA = -SIGN(mobbrmsd_DLAPY2(ALPHA, XNORM), ALPHA)
    end if
    TAU = (BETA - ALPHA) / BETA
    call mobbrmsd_DSCAL(N - 1, ONE / (ALPHA - BETA), X, INCX)
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
    do J = 1, KNT
      BETA = BETA * SAFMIN
    end do
    ALPHA = BETA
  end if
!
  return
!
!     End of mobbrmsd_DLARFG
!
end subroutine mobbrmsd_DLARFG

