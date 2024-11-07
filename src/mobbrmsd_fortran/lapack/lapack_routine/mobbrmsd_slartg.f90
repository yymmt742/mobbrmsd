!| mobbrmsd_SLARTG generates a plane rotation with real cosine and real sine.
!
! mobbrmsd_SLARTG generates a plane rotation so that
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
!  Reference SLARTG is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
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
pure elemental subroutine mobbrmsd_SLARTG(f, g, c, s, r)
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
  real(RK)  :: D, F1, FS, G1, GS, U, RTMIN, RTMAX
  intrinsic :: ABS, SIGN, SQRT
!  ..
  RTMIN = SQRT(SAFMIN)
  RTMAX = SQRT(SAFMAX / 2)
!
  F1 = ABS(F)
  G1 = ABS(G)
  if (G == ZERO) then
    C = ONE
    S = ZERO
    R = F
  else if (F == ZERO) then
    C = ZERO
    S = SIGN(ONE, G)
    R = G1
  else if (F1 > RTMIN .and. F1 < RTMAX .and. &
           G1 > RTMIN .and. G1 < RTMAX) then
    D = SQRT(F * F + G * G)
    C = F1 / D
    R = SIGN(D, F)
    S = G / R
  else
    U = MIN(SAFMAX, MAX(SAFMIN, F1, G1))
    FS = F / U
    GS = G / U
    D = SQRT(FS * FS + GS * GS)
    C = ABS(FS) / D
    R = SIGN(D, F)
    S = GS / R
    R = R * U
  end if
  return
end

