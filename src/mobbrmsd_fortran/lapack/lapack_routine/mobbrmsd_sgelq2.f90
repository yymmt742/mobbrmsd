!| mobbrmsd_SGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.
!
!  mobbrmsd_SGELQ2 computes an LQ factorization of a real m-by-n matrix A:
!
!     A = ( L 0 ) *  Q
!
!  where:
!
!     Q is a n-by-n orthogonal matrix;
!     L is a lower-triangular m-by-m matrix;
!     0 is a m-by-(n-m) zero matrix, if m < n.
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!  and tau in TAU(i).
!
!  Reference SGELQ2 is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK computational routine (version 3.9.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
pure subroutine mobbrmsd_SGELQ2(M, N, A, LDA, TAU, WORK, INFO)
  implicit none
  integer, intent(in)  :: M
!!          The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)  :: N
!!          The number of columns of the matrix A.  N >= 0.
!!
  integer, intent(in)  :: LDA
!!          LDA is INTEGER
!!          The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!          A is REAL array, dimension (LDA,N)
!!
!!          On entry, the m by n matrix A.
!!
!!          On exit, the elements on and below the diagonal of the array
!!          contain the m by min(m,n) lower trapezoidal matrix L (L is
!!          lower triangular if m <= n); the elements above the diagonal,
!!          with the array TAU, represent the orthogonal matrix Q as a
!!          product of elementary reflectors (see Further Details).
!!
  real(RK), intent(out)   :: TAU(*)
!!          TAU is REAL array, dimension (min(M,N))
!!          The scalar factors of the elementary reflectors (see Further
!!          Details).
!!
  real(RK), intent(out)   :: WORK(*)
!!          WORK is REAL array, dimension (M)
!!
  integer, intent(out)    :: INFO
!!          = 0: successful exit
!!
!!          < 0: if INFO = -i, the i-th argument had an illegal value
!!
  intrinsic :: MAX, MIN
  real(RK)  :: AII
  integer   :: I, K
! real(RK), parameter :: ONE = 1.0E+0
! interface
!   include 'slarf.h'
!   include 'slarfg.h'
! end interface
!
! Test the input arguments
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
! call XERBLA('SGELQ2', -INFO)
    return
  end if
!
  K = MIN(M, N)
!
  do I = 1, K
!
! Generate elementary reflector H(i) to annihilate A(i, i + 1:n)
!
    call mobbrmsd_SLARFG(N - I + 1, A(I, I), A(I, MIN(I + 1, N)), LDA, TAU(I))
    if (I < M) then
!
! Apply H(i) to A(i + 1:m, i:n) from the right
!
      AII = A(I, I)
      A(I, I) = ONE
      call mobbrmsd_SLARF('Right', M - I, N - I + 1, A(I, I), LDA, TAU(I), A(I + 1, I), LDA, WORK)
      A(I, I) = AII
    end if
  end do
  return
!
! end of mobbrmsd_SGELQ2
!
end

