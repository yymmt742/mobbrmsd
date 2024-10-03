!> \brief \b mobbrmsd_DLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dlacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dlacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dlacpy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
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
!> \param[in] A
!> \verbatim
!>          A is real(RK)           :: array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is real(RK)           :: array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
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
pure subroutine mobbrmsd_DLACPY(UPLO, M, N, A, LDA, B, LDB)
! use LA_CONSTANTS, only: RK => DP
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character, intent(in) :: UPLO
  integer, intent(in)   :: LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)  :: A(LDA, *)
  real(RK), intent(out) :: B(LDB, *)
!     ..
!
!  =====================================================================
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
    do concurrent(J=1:N)
      do concurrent(I=1:MIN(J, M))
        B(I, J) = A(I, J)
      end do
    end do
  else if (mobbrmsd_LSAME(UPLO, 'L')) then
    do concurrent(J=1:N)
      do concurrent(I=J:M)
        B(I, J) = A(I, J)
      end do
    end do
  else
    do concurrent(J=1:N, I=1:M)
      B(I, J) = A(I, J)
    end do
  end if
!
!     End of mobbrmsd_DLACPY
!
end subroutine mobbrmsd_DLACPY
