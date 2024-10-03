!> \brief \b mobbrmsd_DGETRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE pure subroutine mobbrmsd_DGETRF2( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(RK)           ::   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DGETRF2 computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = min(m,n)/2
!>        [  A21 | A22  ]       n2 = n-n1
!>
!>                                       [ A11 ]
!> The subroutine calls itself to factor [ --- ],
!>                                       [ A12 ]
!>                 [ A12 ]
!> do the swaps on [ --- ], solve A12, update A22,
!>                 [ A22 ]
!>
!> then calls itself to factor A22 and do the swaps on A21.
!>
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
!>          A is real(RK)           :: array, dimension (LDA,N)
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
!> \ingroup doubleGEcomputational
!
!  =====================================================================
pure recursive subroutine mobbrmsd_DGETRF2(M, N, A, LDA, IPIV, INFO)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)      :: LDA, M, N
  integer, intent(out)     :: INFO
!     ..
!     .. Array Arguments ..
  integer, intent(out)     :: IPIV(*)
  real(RK), intent(inout)  :: A(LDA, *)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  real(RK), parameter      :: ZERO = 0.0_RK
  real(RK), parameter      :: ONE = 1.0_RK
!     ..
!     .. Local Scalars ..
  real(RK)                :: SFMIN, TEMP
  integer                 :: I, IINFO, N1, N2
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: MAX, MIN
!     ..
! interface
!     .. External Subroutines ..
!   include 'dgemm.h'
!   include 'dscal.h'
!   include 'dlaswp.h'
!   include 'dtrsm.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'dlamch.h'
!   include 'idamax.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
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
    !CALL XERBLA( 'DGETRF2', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (M == 0 .or. N == 0) return

  if (M == 1) then
!
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO
!
    IPIV(1) = 1
    if (A(1, 1) == ZERO) INFO = 1
!
  else if (N == 1) then
!
!        Use unblocked code for one column case
!
!
!        Compute machine safe minimum
!
    SFMIN = mobbrmsd_DLAMCH('S')
!
!        Find pivot and test for singularity
!
    I = mobbrmsd_IDAMAX(M, A(1, 1), 1)
    IPIV(1) = I
    if (A(I, 1) /= ZERO) then
!
!           Apply the interchange
!
      if (I /= 1) then
        TEMP = A(1, 1)
        A(1, 1) = A(I, 1)
        A(I, 1) = TEMP
      end if
!
!           Compute elements 2:M of the column
!
      if (ABS(A(1, 1)) >= SFMIN) then
        call mobbrmsd_DSCAL(M - 1, ONE / A(1, 1), A(2, 1), 1)
      else
        do concurrent(I=1:M - 1)
          A(1 + I, 1) = A(1 + I, 1) / A(1, 1)
        end do
      end if
!
    else
      INFO = 1
    end if
!
  else
!
!      Use recursive code
!
    N1 = MIN(M, N) / 2
    N2 = N - N1
!
!           [ A11 ]
!    Factor [ --- ]
!           [ A21 ]
!
    call mobbrmsd_DGETRF2(M, N1, A, LDA, IPIV, IINFO)

    if (INFO == 0 .and. IINFO > 0) INFO = IINFO
!
!                          [ A12 ]
!    Apply interchanges to [ --- ]
!                          [ A22 ]
!
    call mobbrmsd_DLASWP(N2, A(1, N1 + 1), LDA, 1, N1, IPIV, 1)
!
!    Solve A12
!
    call mobbrmsd_DTRSM('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, &
   &            A(1, N1 + 1), LDA)
!
!    Update A22
!
    call mobbrmsd_DGEMM('N', 'N', M - N1, N2, N1, -ONE, A(N1 + 1, 1), LDA, &
   &            A(1, N1 + 1), LDA, ONE, A(N1 + 1, N1 + 1), LDA)
!
!    Factor A22
!
    call mobbrmsd_DGETRF2(M - N1, N2, A(N1 + 1, N1 + 1), LDA, IPIV(N1 + 1), &
   &              IINFO)
!
!    Adjust INFO and the pivot indices
!
    if (INFO == 0 .and. IINFO > 0) INFO = IINFO + N1
    do concurrent(I=N1 + 1:MIN(M, N))
      IPIV(I) = IPIV(I) + N1
    end do
!
!    Apply interchanges to A21
!
    call mobbrmsd_DLASWP(N1, A(1, 1), LDA, N1 + 1, MIN(M, N), IPIV, 1)
!
  end if
  return
!
!     End of mobbrmsd_DGETRF2
!
end subroutine mobbrmsd_DGETRF2

