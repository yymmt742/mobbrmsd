!> \brief \b mobbrmsd_SLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SLACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacpy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SLACPY copies all or part of a two-dimensional matrix A to another
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
!>          A is REAL array, dimension (LDA,N)
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
!>          B is REAL array, dimension (LDB,N)
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
!> \date December 2016
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_SLACPY(UPLO, M, N, A, LDA, B, LDB)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  character, intent(in) :: UPLO
  integer, intent(in)   :: LDA, LDB, M, N
!..
!..Array Arguments..
  real(RK), intent(in)  :: A(LDA, *)
  real(RK), intent(out) :: B(LDB, *)
!..
!
!  =====================================================================
!
!..Local Scalars..
  integer :: I, J
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: MIN
!..
!..Executable Statements..
!
  if (mobbrmsd_LSAME(UPLO, 'U')) then
    do J = 1, N
      do I = 1, MIN(J, M)
        B(I, J) = A(I, J)
      end do
    end do
  else if (mobbrmsd_LSAME(UPLO, 'L')) then
    do J = 1, N
      do I = J, M
        B(I, J) = A(I, J)
      end do
    end do
  else
    do J = 1, N
      do I = 1, M
        B(I, J) = A(I, J)
      end do
    end do
  end if
  return
  !
  !end of mobbrmsd_SLACPY
  !
end
