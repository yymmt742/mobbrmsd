!| multiplies a general matrix by the orthogonal matrix (unblocked algorithm).
!
!  mobbrmsd_DORM2R multiplies a general matrix by the orthogonal matrix
!  from a QR factorization determined by mobbrmsd_DGEQRF.
!  mobbrmsd_SORM2R overwrites the general real \( m \) by \( n \) matrix \( C \) with
!
!  \( Q C \) if SIDE = 'L' and TRANS = 'N', or
!
!  \( Q ^ \top C \) if SIDE = 'L' and TRANS = 'T', or
!
!  \( C Q \) if SIDE = 'R' and TRANS = 'N', or
!
!  \( C Q ^ \top \) if SIDE = 'R' and TRANS = 'T'.
!
!  where \( Q \) is a real orthogonal matrix defined
!  as the product of \( k \) elementary reflectors
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ k
!  \]
!
! as returned by mobbrmsd_DGEQRF.
! \( Q \) is of order \( m \) if SIDE = 'L' and of order n if SIDE = 'R'.
!
!  Reference DORM2R is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DORM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, &
     &                          C, LDC, WORK, INFO)
  implicit none
  character, intent(in)   :: SIDE
!!  = 'L': apply Q or Q**T from the Left
!!
!!  = 'R': apply Q or Q**T from the Right
!!
  character, intent(in)   :: TRANS
!!  = 'N': apply Q  (No transpose)
!!
!!  = 'T': apply Q**T (Transpose)
!!
  integer, intent(in)     :: M
!!  The number of rows of the matrix C. M >= 0.
!!
  integer, intent(in)     :: N
!!  The number of columns of the matrix C. N >= 0.
!!
  integer, intent(in)     :: K
!!  The number of elementary reflectors whose product defines
!!  the matrix Q.
!!
!!  If SIDE = 'L', M >= K >= 0;
!!  if SIDE = 'R', N >= K >= 0.
!!
  integer, intent(in)     :: LDA
!!   The leading dimension of the array A.
!!
!!   If SIDE = 'L', LDA >= max(1,M);
!!
!!   if SIDE = 'R', LDA >= max(1,N).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  REAL array, dimension (LDA,K)
!!
!!  The i-th column must contain the vector which defines the
!!  elementary reflector H(i), for i = 1,2,...,k, as returned by
!!  mobbrmsd_SGEQRF in the first k columns of its array argument A.
!!  A is modified by the routine but restored on exit.
!!
  real(RK), intent(in)    :: TAU(*)
!!  REAL array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_SGEQRF.
!!
  integer, intent(in)     :: LDC
!!  The leading dimension of the array C. LDC >= max(1,M).
!!
  real(RK), intent(inout) :: C(LDC, *)
!!  REAL array, dimension (LDC,N)
!!
!!  On entry, the m by n matrix C.
!!
!!  On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!!
  real(RK), intent(out)   :: WORK(*)
!!  REAL array, dimension
!!
!!  (N) if SIDE = 'L',
!!
!!  (M) if SIDE = 'R'
!!
  integer, intent(out)    :: INFO
!!  = 0: successful exit
!!
!!  < 0: if INFO = -i, the i-th argument had an illegal value
!!
  logical                 :: LEFT, NOTRAN
  integer                 :: I, I1, I2, I3, IC, JC, MI, NI, NQ
  real(RK)                :: AII
  intrinsic               :: MAX
!
! interface
!     .. External Subroutines ..
!   include 'dlarf.h'
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!
! Test the input arguments
!
  INFO = 0
  LEFT = mobbrmsd_LSAME(SIDE, 'L')
  NOTRAN = mobbrmsd_LSAME(TRANS, 'N')
!
!     NQ is the order of Q
!
  if (LEFT) then
    NQ = M
  else
    NQ = N
  end if
  if (.not. LEFT .and. .not. mobbrmsd_LSAME(SIDE, 'R')) then
    INFO = -1
  else if (.not. NOTRAN .and. .not. mobbrmsd_LSAME(TRANS, 'T')) then
    INFO = -2
  else if (M < 0) then
    INFO = -3
  else if (N < 0) then
    INFO = -4
  else if (K < 0 .or. K > NQ) then
    INFO = -5
  else if (LDA < MAX(1, NQ)) then
    INFO = -7
  else if (LDC < MAX(1, M)) then
    INFO = -10
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DORM2R', -INFO )
    return
  end if
!
!     Quick return if possible
!
  if (M == 0 .or. N == 0 .or. K == 0) return
!
  if ((LEFT .and. .not. NOTRAN) .or. (.not. LEFT .and. NOTRAN)) then
    I1 = 1
    I2 = K
    I3 = 1
  else
    I1 = K
    I2 = 1
    I3 = -1
  end if
!
  if (LEFT) then
    NI = N
    JC = 1
  else
    MI = M
    IC = 1
  end if
!
  do I = I1, I2, I3
    if (LEFT) then
!
!           H(i) is applied to C(i:m,1:n)
!
      MI = M - I + 1
      IC = I
    else
!
!           H(i) is applied to C(1:m,i:n)
!
      NI = N - I + 1
      JC = I
    end if
!
!        Apply H(i)
!
    AII = A(I, I)
    A(I, I) = ONE
    call mobbrmsd_DLARF(SIDE, MI, NI, A(I, I), 1, TAU(I), C(IC, JC), &
   &            LDC, WORK)
    A(I, I) = AII
  end do
  return
!
!     End of mobbrmsd_DORM2R
!
end subroutine mobbrmsd_DORM2R

