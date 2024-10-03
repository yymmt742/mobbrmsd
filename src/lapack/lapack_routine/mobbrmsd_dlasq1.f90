!> \brief \b mobbrmsd_DLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLASQ1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DLASQ1( N, D, E, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLASQ1 computes the singular values of a real N-by-N bidiagonal
!> matrix with diagonal D and off-diagonal E. The singular values
!> are computed to high relative accuracy, in the absence of
!> denormalization, underflow and overflow. The algorithm was first
!> presented in
!>
!> "Accurate singular values and differential qd algorithms" by K. V.
!> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
!> 1994,
!>
!> and the present implementation is described in "An implementation of
!> the dqds Algorithm (Positive Case)", LAPACK Working Note.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>        The number of rows and columns in the matrix. N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>        On entry, D contains the diagonal elements of the
!>        bidiagonal matrix whose SVD is desired. On normal exit,
!>        D contains the singular values in decreasing order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>        On entry, elements E(1:N-1) contain the off-diagonal elements
!>        of the bidiagonal matrix whose SVD is desired.
!>        On exit, E is overwritten.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>        = 0: successful exit
!>        < 0: if INFO = -i, the i-th argument had an illegal value
!>        > 0: the algorithm failed
!>             = 1, a split was marked by a positive value in E
!>             = 2, current block of Z not diagonalized after 100*N
!>                  iterations (in inner while loop)  On exit D and E
!>                  represent a matrix with the same singular values
!>                  which the calling subroutine could use to finish the
!>                  computation, or even feed back into mobbrmsd_DLASQ1
!>             = 3, termination criterion of outer while loop not met
!>                  (program created more than N unreduced blocks)
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
pure subroutine mobbrmsd_DLASQ1(N, D, E, WORK, INFO)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: N
  integer, intent(out)    :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: D(*), E(*)
  real(RK), intent(out)   :: WORK(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  integer            :: I, IINFO
  real(RK)           :: EPS, SCALE, SAFMIN, SIGMN, SIGMX
!     ..
!     .. Parameters ..
! real(RK), parameter :: ZERO = 0.0_RK
!     ..
! interface
!   include 'dcopy.h'
!   include 'dlas2.h'
!   include 'dlascl.h'
!   include 'dlasq2.h'
!   include 'dlasrt.h'
!   !include 'xerbla.h'
!   include 'dlamch.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic          :: ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
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
!     Maximum number of iterations exceeded.  Move data from WORK
!     into D and E so the calling subroutine can try to finish
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
!     End of mobbrmsd_DLASQ1
!
end subroutine mobbrmsd_DLASQ1
