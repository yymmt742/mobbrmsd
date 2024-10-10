!| mobbrmsd_DGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.
!
!  mobbrmsd_DGEQR2 computes a QR factorization of a real m-by-n matrix A:
!
!     A = Q * ( R ),
!             ( 0 )
!
!  where:
!
!     Q is a m-by-m orthogonal matrix;
!     R is an upper-triangular n-by-n matrix;
!     0 is a (m-n)-by-n zero matrix, if m > n.
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  Reference DGEQR2 is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DGEQR2(M, N, A, LDA, TAU, WORK, INFO)
  implicit none
  integer, intent(in)     :: M
!!  The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix A.  N >= 0.
!!
  integer, intent(in)     :: LDA
!!  The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,N)
!!
!!  On entry, the m by n matrix A.
!!
!!  On exit, the elements on and above the diagonal of the array
!!  contain the min(m,n) by n upper trapezoidal matrix R (R is
!!  upper triangular if m >= n); the elements below the diagonal,
!!  with the array TAU, represent the orthogonal matrix Q as a
!!  product of elementary reflectors (see Further Details).
!!
  real(RK), intent(out)   :: TAU(*)
!!  DOUBLE PRECISION array, dimension (min(M,N))
!!  The scalar factors of the elementary reflectors (see Further
!!  Details).
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (N)
!!
  integer, intent(out)    :: INFO
!!  = 0: successful exit
!!
!!  < 0: if INFO = -i, the i-th argument had an illegal value
!!
  integer                 :: I, K
  real(RK)                :: AII
  intrinsic          :: MAX, MIN
! interface
!   include 'dlarf.h'
!   include 'dlarfg.h'
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
    return
  end if
!
  K = MIN(M, N)
!
  do I = 1, K
!
!   Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
    call mobbrmsd_DLARFG(M - I + 1, A(I, I), A(MIN(I + 1, M), I), 1, TAU(I))
!
    if (I < N) then
!
!     Apply H(i) to A(i:m,i+1:n) from the left
!
      AII = A(I, I)
      A(I, I) = ONE
      call mobbrmsd_DLARF('Left', M - I + 1, N - I, A(I, I), 1, TAU(I), &
          &      A(I, I + 1), LDA, WORK)
      A(I, I) = AII
    end if
  end do
  return
!
!     End of mobbrmsd_DGEQR2
!
end subroutine mobbrmsd_DGEQR2

