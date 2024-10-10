!| mobbrmsd_DLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr.
!
!  mobbrmsd_DLASQ1 computes the singular values of a real N-by-N bidiagonal
!  matrix with diagonal D and off-diagonal E. The singular values
!  are computed to high relative accuracy, in the absence of
!  denormalization, underflow and overflow. The algorithm was first
!  presented in
!
!  "Accurate singular values and differential qd algorithms" by K. V.
!  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
!  1994,
!
!  and the present implementation is described in "An implementation of
!  the dqds Algorithm (Positive Case)", LAPACK Working Note.
!
!  Reference DLASQ1 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLASQ1(N, D, E, WORK, INFO)
  implicit none
  integer, intent(in)     :: N
!!  The number of rows and columns in the matrix. N >= 0.
!!
  real(RK), intent(inout) :: D(*)
!!  DOUBLE PRECISION array, dimension (N)
!!
!!  On entry, D contains the diagonal elements of the
!!  bidiagonal matrix whose SVD is desired. On normal exit,
!!  D contains the singular values in decreasing order.
!!
  real(RK), intent(inout) :: E(*)
!!  DOUBLE PRECISION array, dimension (N)
!!
!!  On entry, elements E(1:N-1) contain the off-diagonal elements
!!  of the bidiagonal matrix whose SVD is desired.
!!  On exit, E is overwritten.
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (4*N)
!!
  integer, intent(out)    :: INFO
!!  = 0: successful exit
!!
!!  < 0: if INFO = -i, the i-th argument had an illegal value
!!
!!  \> 0: the algorithm failed
!!       = 1, a split was marked by a positive value in E
!!       = 2, current block of Z not diagonalized after 100*N
!!            iterations (in inner while loop)  On exit D and E
!!            represent a matrix with the same singular values
!!            which the calling subroutine could use to finish the
!!            computation, or even feed back into mobbrmsd_DLASQ1
!!       = 3, termination criterion of outer while loop not met
!!            (program created more than N unreduced blocks)
!!
  integer            :: I, IINFO
  real(RK)           :: EPS, SCALE, SAFMIN, SIGMN, SIGMX
  intrinsic          :: ABS, MAX, SQRT
! interface
!   include 'dcopy.h'
!   include 'dlas2.h'
!   include 'dlascl.h'
!   include 'dlasq2.h'
!   include 'dlasrt.h'
!   include 'dlamch.h'
! end interface
!
  INFO = 0
  if (N < 0) then
    INFO = -1
    !CALL XERBLA( 'DLASQ1', -INFO )
    return
  else if (N == 0) then
    return
  else if (N == 1) then
    D(1) = ABS(D(1))
    return
  else if (N == 2) then
    call mobbrmsd_DLAS2(D(1), E(1), D(2), SIGMN, SIGMX)
    D(1) = SIGMX
    D(2) = SIGMN
    return
  end if
!
!     Estimate the largest singular value.
!
  SIGMX = ZERO
  do I = 1, N - 1
    D(I) = ABS(D(I))
    SIGMX = MAX(SIGMX, ABS(E(I)))
  end do
  D(N) = ABS(D(N))
!
!     Early return if SIGMX is zero (matrix is already diagonal).
!
  if (SIGMX == ZERO) then
    call mobbrmsd_DLASRT('D', N, D, IINFO)
    return
  end if
!
  do I = 1, N
    SIGMX = MAX(SIGMX, D(I))
  end do
!
!     Copy D and E into WORK (in the Z format) and scale (squaring the
!     input data makes scaling by a power of the radix pointless).
!
  EPS = mobbrmsd_DLAMCH('Precision')
  SAFMIN = mobbrmsd_DLAMCH('Safe minimum')
  SCALE = SQRT(EPS / SAFMIN)
  call mobbrmsd_DCOPY(N, D, 1, WORK(1), 2)
  call mobbrmsd_DCOPY(N - 1, E, 1, WORK(2), 2)
  call mobbrmsd_DLASCL('G', 0, 0, SIGMX, SCALE, 2 * N - 1, 1, WORK, 2 * N - 1, IINFO)
!
!     Compute the q's and e's.
!
  do concurrent(I=1:2 * N - 1)
    WORK(I) = WORK(I)**2
  end do
!       do 30 I = 1, 2 * N - 1
!         WORK(I) = WORK(I)**2
!30     continue
  WORK(2 * N) = ZERO
!
  call mobbrmsd_DLASQ2(N, WORK, INFO)
!
  if (INFO == 0) then
    do concurrent(I=1:N)
      D(I) = SQRT(WORK(I))
    end do
!         do 40 I = 1, N
!           D(I) = SQRT(WORK(I))
!40       continue
    call mobbrmsd_DLASCL('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO)
  else if (INFO == 2) then
!
!   Maximum number of iterations exceeded.  Move data from WORK
!   into D and E so the calling subroutine can try to finish
!
    do I = 1, N
      D(I) = SQRT(WORK(2 * I - 1))
      E(I) = SQRT(WORK(2 * I))
    end do
    call mobbrmsd_DLASCL('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO)
    call mobbrmsd_DLASCL('G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO)
  end if
!
  return
!
! End of mobbrmsd_DLASQ1
!
end subroutine mobbrmsd_DLASQ1

