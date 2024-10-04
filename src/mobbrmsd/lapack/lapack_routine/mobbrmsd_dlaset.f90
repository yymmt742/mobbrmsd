!> \brief \b mobbrmsd_DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLASET + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dlaset.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dlaset.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dlaset.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       pure subroutine mobbrmsd_DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       real(RK)           ::   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set; the strictly lower
!>                      triangular part of A is not changed.
!>          = 'L':      Lower triangular part is set; the strictly upper
!>                      triangular part of A is not changed.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is real(RK)           ::
!>          The constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is real(RK)           ::
!>          The constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          On exit, the leading m-by-n submatrix of A is set as follows:
!>
!>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_DLASET(UPLO, M, N, ALPHA, BETA, A, LDA)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: UPLO
  integer, intent(in)      :: LDA, M, N
  real(RK), intent(in)     :: ALPHA, BETA
!     ..
!     .. Array Arguments ..
  real(RK), intent(out)    :: A(LDA, *)
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
  integer                 :: I, J
!     .. Intrinsic Functions ..
  intrinsic               :: MIN
!     ..
! interface
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!     ..
!     .. Executable Statements ..
!
  if (mobbrmsd_LSAME(UPLO, 'U')) then
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
    do J = 2, N
      do I = 1, MIN(J - 1, M)
        A(I, J) = ALPHA
      end do
    end do
!
  else if (mobbrmsd_LSAME(UPLO, 'L')) then
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
    do J = 1, MIN(M, N)
      do I = J + 1, M
        A(I, J) = ALPHA
      end do
    end do
!
  else
!
!        Set the leading m-by-n submatrix to ALPHA.
!
    do concurrent(J=1:N, I=1:M)
      A(I, J) = ALPHA
    end do
!         do 60 J = 1, N
!           do 50 I = 1, M
!             A(I, J) = ALPHA
!50         continue
!60       continue
  end if
!
!     Set the first min(M,N) diagonal elements to BETA.
!
  do concurrent(I=1:MIN(M, N))
    A(I, I) = BETA
  end do
!       do 70 I = 1, MIN(M, N)
!         A(I, I) = BETA
!70     continue
!
  return
!
!     End of mobbrmsd_DLASET
!
end subroutine mobbrmsd_DLASET

