!| generates an matrix \( Q \in \mathbb{R} ^ {m \times n} \) with orthonormal rows
!
!  mobbrmsd_DORGL2 generates an \( m \) by \( n \) real matrix \( Q \)
!  with orthonormal rows, which is defined as the first \( m \) rows
!  of a product of \( k \) elementary reflectors of order \( n \).
!
!  \[
!    Q  =  H _ k \cdots H _ 2 H _ 1
!  \]
!
!  as returned by mobbrmsd_DGELQF.
!
!  Reference DORGL2 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DORGL2(M, N, K, A, LDA, TAU, WORK, INFO)
  implicit none
  integer, intent(in)     :: M
!! The number of rows of the matrix Q. M >= 0.
!!
  integer, intent(in)     :: N
!! The number of columns of the matrix Q. N >= M.
!!
  integer, intent(in)     :: K
!! The number of elementary reflectors whose product defines the
!! matrix Q. M >= K >= 0.
!!
  integer, intent(in)     :: LDA
!!  The first dimension of the array A. LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,N)
!!
!!  On entry, the i-th row must contain the vector which defines
!!  the elementary reflector H(i), for i = 1,2,...,k, as returned
!!  by mobbrmsd_DGELQF in the first k rows of its array argument A.
!!
!!  On exit, the m-by-n matrix Q.
!!
  real(RK), intent(in)    :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_DGELQF.
!!
  real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (M)
!!
  integer, intent(out)    :: INFO
!!  = 0: successful exit
!!
!!  < 0: if INFO = -i, the i-th argument has an illegal value
!!
  integer                 :: I, J, L
  intrinsic               :: MAX
! interface
!   include 'dlarf.h'
!   include 'dscal.h'
! end interface
!
! Test the input arguments
!
  inFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < M) then
    INFO = -2
  else if (K < 0 .or. K > M) then
    INFO = -3
  else if (LDA < MAX(1, M)) then
    INFO = -5
  end if
  if (INFO /= 0) then
    return
  end if
!
! Quick return if possible
!
  if (M <= 0) return
!
  if (K < M) then
!
!   Initialise rows k+1:m to rows of the unit matrix
!
    do J = 1, N
      do L = K + 1, M
        A(L, J) = ZERO
      end do
      if (J > K .and. J <= M) A(J, J) = ONE
    end do
  end if
!
  do I = K, 1, -1
!
!   Apply H(i) to A(i:m,i:n) from the right
!
    if (I < N) then
      if (I < M) then
        A(I, I) = ONE
        call mobbrmsd_DLARF('Right', M - I, N - I + 1, A(I, I), LDA, &
       &            TAU(I), A(I + 1, I), LDA, WORK)
      end if
      call mobbrmsd_DSCAL(N - I, -TAU(I), A(I, I + 1), LDA)
    end if
    A(I, I) = ONE - TAU(I)
!
!   Set A(i,1:i-1) to zero
!
    do L = 1, I - 1
      A(I, L) = ZERO
    end do
  end do
  return
!
! End of mobbrmsd_DORGL2
!
end subroutine mobbrmsd_DORGL2

