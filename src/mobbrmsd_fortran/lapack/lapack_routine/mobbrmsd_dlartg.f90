!| mobbrmsd_DLARTG generates a plane rotation with real cosine and real sine.
!
! mobbrmsd_DLARTG generates a plane rotation so that
!
!    [  C  S  ]  .  [ F ]  =  [ R ]
!    [ -S  C  ]     [ G ]     [ 0 ]
!
! where C**2 + S**2 = 1.
!
! The mathematical formulas used for C and S are
!    R = sign(F) * sqrt(F**2 + G**2)
!    C = F / R
!    S = G / R
! Hence C >= 0. The algorithm used to compute these quantities
! incorporates scaling to avoid overflow or underflow in computing the
! square root of the sum of squares.
!
! This version is discontinuous in R at F = 0 but it returns the same
! C and S as ZLARTG for complex inputs (F,0) and (G,0).
!
! This is a more accurate version of the BLAS1 routine mobbrmsd_DROTG,
! with the following other differences:
!    F and G are unchanged on return.
!    If G=0, then C=1 and S=0.
!    If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any
!       floating point operations (saves work in mobbrmsd_DBDSQR when
!       there are zeros on the diagonal).
!
! If F exceeds G in magnitude, C will be positive.
!
! Below, RK=>dp stands for double precision from LA_CONSTANTS module.
!
! Further Details:
!
!  Anderson E. (2017)
!  Algorithm 978: Safe Scaling in the Level 1 BLAS
!  [ACM Trans Math Softw 44:1--28](https://doi.org/10.1145/3061665)
!
!  Reference DLARTG is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
! Weslley Pereira, University of Colorado Denver, USA
!
pure elemental subroutine mobbrmsd_DLARTG(f, g, c, s, r)
  implicit none
  real(RK), intent(in)  :: f
!!  The first component of vector to be rotated.
!!
  real(RK), intent(in)  :: g
!!  The second component of vector to be rotated.
!!
  real(RK), intent(out) :: c
!!  The cosine of the rotation.
!!
  real(RK), intent(out) :: s
!!  The sine of the rotation.
!!
  real(RK), intent(out) :: r
!!  The nonzero component of the rotated vector.
!!
  real(RK)             :: d, f1, fs, g1, gs, p, u, uu
  intrinsic            :: abs, sign, sqrt
!
  f1 = ABS(f)
  g1 = ABS(g)
  if (g == zero) then
    c = one
    s = zero
    r = f
  else if (f == zero) then
    c = zero
    s = SIGN(one, g)
    r = g1
  else if (f1 > RTMIN .and. f1 < RTMAX .and. &
           g1 > RTMIN .and. g1 < RTMAX) then
    d = SQRT(f * f + g * g)
    p = one / d
    c = f1 * p
    s = g * SIGN(p, f)
    r = SIGN(d, f)
  else
    u = MIN(SAFMAX, MAX(SAFMIN, f1, g1))
    uu = one / u
    fs = f * uu
    gs = g * uu
    d = SQRT(fs * fs + gs * gs)
    p = one / d
    c = ABS(fs) * p
    s = gs * SIGN(p, f)
    r = SIGN(d, f) * u
  end if
  return
end subroutine mobbrmsd_DLARTG

