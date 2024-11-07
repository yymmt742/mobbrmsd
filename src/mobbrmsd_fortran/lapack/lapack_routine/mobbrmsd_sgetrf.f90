!| computes an LU factorization using partial pivoting.
!
!  mobbrmsd_SGETRF computes an LU factorization of a general
!  \( M \) -by- \( N \) matrix  \( \mathbf{A} \)
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!
!  \[ \mathbf{A} = \mathbf{P L U} \]
!
!  where \( \mathbf{P} \) is a permutation matrix,
!  \( \mathbf{L} \) is lower triangular with unit
!  diagonal elements (lower trapezoidal if \( m > n \) ),
!  and \( \mathbf{U} \) is upper
!  triangular (upper trapezoidal if \( m < n \) ).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!  Reference SGETRF is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK computational routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SGETRF(M, N, A, LDA, IPIV, INFO)
  implicit none
  integer, intent(in)      :: M
!!  The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)      :: N
!!  The number of columns of the matrix A.  N >= 0.
!!
  integer, intent(in)      :: LDA
!!  The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout)  :: A(LDA, *)
!!  A is REAL array, dimension (LDA,N)
!!
!!  On entry, the M-by-N matrix to be factored.
!!
!!  On exit, the factors L and U from the factorization
!!  A = P*L*U; the unit diagonal elements of L are not stored.
!!
  integer, intent(out)     :: IPIV(*)
!!  INTEGER array, dimension (min(M,N))
!!
!!  The pivot indices; for 1 <= i <= min(M,N), row i of the
!!  matrix was interchanged with row IPIV(i).
!!
  integer, intent(out)     :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
!!  \> 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!!        has been completed, but the factor U is exactly
!!        singular, and division by zero will occur if it is used
!!        to solve a system of equations.
!!
  integer   :: I, IINFO, J, JB, NB
  intrinsic :: MAX, MIN
! interface
!   include 'ilaenv.h'
!   include 'sgemm.h'
!   include 'sgetrf2.h'
!   include 'slaswp.h'
!   include 'strsm.h'
! end interface
!
! Test the input parameters.
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
    return
  end if
!
! Quick return if possible
!
  if (M == 0 .or. N == 0) return
!
! Determine the block size for this environment.
!
  NB = mobbrmsd_ILAENV(1, 'SGETRF', ' ', M, N, -1, -1)
  if (NB <= 1 .or. NB >= MIN(M, N)) then
!
! use unblocked code.
!
    call mobbrmsd_SGETRF2(M, N, A, LDA, IPIV, INFO)
  else
!
! use blocked code.
!
    do J = 1, MIN(M, N), NB
      JB = MIN(MIN(M, N) - J + 1, NB)
!
! Factor diagonal and subdiagonal blocks and test for exact
! singularity.
!
      call mobbrmsd_SGETRF2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO)
!
! Adjust INFO and the pivot indices.
!
      if (INFO == 0 .and. IINFO > 0) INFO = IINFO + J - 1
      do I = J, MIN(M, J + JB - 1)
        IPIV(I) = J - 1 + IPIV(I)
      end do
!
! Apply interchanges to columns 1:J - 1.
!
      call mobbrmsd_SLASWP(J - 1, A, LDA, J, J + JB - 1, IPIV, 1)
!
      if (J + JB <= N) then
!
! Apply interchanges to columns J + JB:N.
!
        call mobbrmsd_SLASWP(N - J - JB + 1, A(1, J + JB), LDA, J, J + JB - 1, IPIV, 1)
!
! Compute block row of U.
!
        call mobbrmsd_STRSM('Left', 'Lower', 'No transpose', 'Unit', JB, &
            & N - J - JB + 1, ONE, A(J, J), LDA, A(J, J + JB), LDA)
        if (J + JB <= M) then
!
! Update trailing submatrix.
!
          call mobbrmsd_SGEMM('No transpose', 'No transpose', M - J - JB + 1, &
              &      N - J - JB + 1, JB, -ONE, A(J + JB, J), LDA, &
              &      A(J, J + JB), LDA, ONE, A(J + JB, J + JB), LDA)
        end if
      end if
    end do
  end if
  return
!
!end of mobbrmsd_SGETRF
!
end

