!> \brief \b SGETRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGETRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGETRF computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 3 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
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
!> \ingroup realGEcomputational
!
!  =====================================================================
pure subroutine SGETRF(M, N, A, LDA, IPIV, INFO)
  implicit none
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!.. Scalar Arguments ..
  integer, intent(in)  :: LDA, M, N
  integer, intent(out) :: INFO
!..
!.. Array Arguments ..
  integer, intent(out) :: IPIV(*)
  real(RK), intent(inout) :: A(LDA, *)
!..
!
!  =====================================================================
!..
!..Local Scalars..
  integer :: I, IINFO, J, JB, NB
!
!..Parameters..
! real(RK), parameter :: ONE = 1.0E+0
!..
! interface
!..external Functions..
!   include 'ilaenv.h'
!..external Subroutines..
!   include 'sgemm.h'
!   include 'sgetrf2.h'
!   include 'slaswp.h'
!   include 'strsm.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN
!..
!..Executable Statements..
!
!Test the input parameters.
!
  INFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  end if
  if (INFO /= 0) then
!   call XERBLA('SGETRF', -INFO)
    return
  end if
!
!Quick return if possible
!
  if (M == 0 .or. N == 0) return
!
!Determine the block size for this environment.
!
  NB = ILAENV(1, 'SGETRF', ' ', M, N, -1, -1)
  if (NB <= 1 .or. NB >= MIN(M, N)) then
    !
    !use unblocked code.
    !
    call SGETRF2(M, N, A, LDA, IPIV, INFO)
  else
    !
    !use blocked code.
    !
    do J = 1, MIN(M, N), NB
      JB = MIN(MIN(M, N) - J + 1, NB)
      !
      !Factor diagonal and subdiagonal blocks and test for exact
      !singularity.
      !
      call SGETRF2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO)
      !
      !Adjust INFO and the pivot indices.
      !
      if (INFO == 0 .and. IINFO > 0) INFO = IINFO + J - 1
      do I = J, MIN(M, J + JB - 1)
        IPIV(I) = J - 1 + IPIV(I)
      end do
      !
      !Apply interchanges to columns 1:J - 1.
      !
      call SLASWP(J - 1, A, LDA, J, J + JB - 1, IPIV, 1)
      !
      if (J + JB <= N) then
        !
        !Apply interchanges to columns J + JB:N.
        !
        call SLASWP(N - J - JB + 1, A(1, J + JB), LDA, J, J + JB - 1, IPIV, 1)
        !
        !Compute block row of U.
        !
        call STRSM('Left', 'Lower', 'No transpose', 'Unit', JB, &
            & N - J - JB + 1, ONE, A(J, J), LDA, A(J, J + JB), LDA)
        if (J + JB <= M) then
          !
          !Update trailing submatrix.
          !
          call SGEMM('No transpose', 'No transpose', M - J - JB + 1, &
              &      N - J - JB + 1, JB, -ONE, A(J + JB, J), LDA, &
              &      A(J, J + JB), LDA, ONE, A(J + JB, J + JB), LDA)
        end if
      end if
    end do
  end if
  return
  !
  !end of SGETRF
  !
end
