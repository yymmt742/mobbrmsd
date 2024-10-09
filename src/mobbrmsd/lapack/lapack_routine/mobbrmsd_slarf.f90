!| mobbrmsd_SLARF applies an elementary reflector to a general rectangular matrix.
!
!  mobbrmsd_SLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v**T
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Reference SLARF is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
  implicit none
  character, intent(in) :: SIDE
!!  = 'L': form  H * C
!!
!!  = 'R': form  C * H
!!
  integer, intent(in)     :: M
!!  The number of rows of the matrix C.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix C.
!!
  real(RK), intent(in)    :: V(*)
!!  REAL array, dimension
!!             (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!!          or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!!
!!  The vector v in the representation of H. V is not used if
!!  TAU = 0.
!!
  integer, intent(in)     :: INCV
!!  The increment between elements of v. INCV <> 0.
!!
  real(RK), intent(in)    :: TAU
!!  The value tau in the representation of H.
!!
  integer, intent(in)     :: LDC
!!  The leading dimension of the array C. LDC >= max(1,M).
!!
  real(RK), intent(inout) :: C(LDC, *)
!!  REAL array, dimension (LDC,N)
!!
!!  On entry, the m by n matrix C.
!!
!!  On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!!  or C * H if SIDE = 'R'.
!!
  real(RK), intent(out)   :: WORK(LDC, *)
!!  REAL array, dimension
!!                 (N) if SIDE = 'L'
!!              or (M) if SIDE = 'R'
  logical :: APPLYLEFT
  integer :: I, LASTV, LASTC
! interface
!   include 'lsame.h'
!   include 'ilaslr.h'
!   include 'ilaslc.h'
!   include 'sgemv.h'
!   include 'sger.h'
! end interface
!
  APPLYLEFT = mobbrmsd_LSAME(SIDE, 'L')
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
      LASTC = mobbrmsd_ILASLC(LASTV, N, C, LDC)
    else
!     Scan for the last non-zero row in C(:,1:lastv).
      LASTC = mobbrmsd_ILASLR(M, LASTV, C, LDC)
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
      call mobbrmsd_SGEMV('Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
          &      ZERO, WORK, 1)
      !
      !C(1:lastv, 1:lastc): = C(...) - v(1:lastv, 1) * w(1:lastc, 1)**T
      !
      call mobbrmsd_SGER(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC)
    end if
  else
    !
    !Form C * H
    !
    if (LASTV > 0) then
      !
      !w(1:lastc, 1): = C(1:lastc, 1:lastv) * v(1:lastv, 1)
      !
      call mobbrmsd_SGEMV('No transpose', LASTC, LASTV, ONE, C, LDC, &
          &      V, INCV, ZERO, WORK, 1)
      !
      !C(1:lastc, 1:lastv): = C(...) - w(1:lastc, 1) * v(1:lastv, 1)**T
      !
      call mobbrmsd_SGER(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC)
    end if
  end if
  return
!
!end of mobbrmsd_SLARF
!
end

