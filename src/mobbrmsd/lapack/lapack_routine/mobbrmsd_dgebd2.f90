!| mobbrmsd_DGEBD2 reduces a real general m by n matrix A to upper or lower
!  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!  If m >= n,
!
!     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n,
!
!     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The contents of A on exit are illustrated by the following examples:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!    (  v1  v2  v3  v4  v5 )
!
!   where d and e denote diagonal and off-diagonal elements of B, vi
!   denotes an element of the vector defining H(i), and ui an element of
!   the vector defining G(i).
!
!     Reference IDAMAX is provided by [](http://www.netlib.org/lapack/)
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DGEBD2(M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO)
  implicit none
  integer, intent(in)     :: M
!!          The number of rows in the matrix A.  M >= 0.
  integer, intent(in)     :: N
!!          The number of columns in the matrix A.  N >= 0.
  integer, intent(in)     :: LDA
!!          The leading dimension of the array A.  LDA >= max(1,M).
  real(RK), intent(inout) :: A(LDA, *)
!!          A is real(RK)           :: array, dimension (LDA,N) <br>
!!          On entry, the m by n general matrix to be reduced. <br>
!!          On exit, <br>
!!          if m >= n, the diagonal and the first superdiagonal are
!!            overwritten with the upper bidiagonal matrix B; the
!!            elements below the diagonal, with the array TAUQ, represent
!!            the orthogonal matrix Q as a product of elementary
!!            reflectors, and the elements above the first superdiagonal,
!!            with the array TAUP, represent the orthogonal matrix P as
!!            a product of elementary reflectors; <br>
!!          if m < n, the diagonal and the first subdiagonal are
!!            overwritten with the lower bidiagonal matrix B; the
!!            elements below the first subdiagonal, with the array TAUQ,
!!            represent the orthogonal matrix Q as a product of
!!            elementary reflectors, and the elements above the diagonal,
!!            with the array TAUP, represent the orthogonal matrix P as
!!            a product of elementary reflectors. <br>
  real(RK), intent(out)   :: D(*)
!!          D is real(RK)           :: array, dimension (min(M,N)) <br>
!!          The diagonal elements of the bidiagonal matrix B: <br>
!!          D(i) = A(i,i).
  real(RK), intent(out)   :: E(*)
!!          E is real(RK)           :: array, dimension (min(M,N)-1) <br>
!!          The off-diagonal elements of the bidiagonal matrix B: <br>
!!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; <br>
!!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  real(RK), intent(out)   :: TAUQ(*)
!!          TAUQ is real(RK)           :: array, dimension (min(M,N)) <br>
!!          The scalar factors of the elementary reflectors which
!!          represent the orthogonal matrix Q. See Further Details.
  real(RK), intent(out)   :: TAUP(*)
!!          TAUP is real(RK)           :: array, dimension (min(M,N)) <br>
!!          The scalar factors of the elementary reflectors which
!!          represent the orthogonal matrix P. See Further Details.
  real(RK), intent(out)   :: WORK(*)
!!          WORK is real(RK)           :: array, dimension (max(M,N))
  integer, intent(out)    :: INFO
!!          INFO is INTEGER <br>
!!          = 0: successful exit. <br>
!!          < 0: if INFO = -i, the i-th argument had an illegal value. <br>
  integer :: I
  intrinsic :: MAX, MIN
!     ..
!     .. Parameters ..
! real(RK), parameter     :: ZERO = 0.0_RK
! real(RK), parameter     :: ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dlarf.h'
!   include 'dlarfg.h'
!   !include 'xerbla.h'
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
  if (INFO < 0) then
!   !CALL XERBLA( 'DGEBD2', -INFO )
    return
  end if
!
  if (M >= N) then
!
!        Reduce to upper bidiagonal form
!
    do I = 1, N
!
!           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
      call mobbrmsd_DLARFG(M - I + 1, A(I, I), A(MIN(I + 1, M), I), 1, TAUQ(I))
      D(I) = A(I, I)
      A(I, I) = ONE
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
      if (I < N)&
&         call mobbrmsd_DLARF('Left', M - I + 1, N - I, A(I, I), 1, TAUQ(I),&
&                     A(I, I + 1), LDA, WORK)
      A(I, I) = D(I)
!
      if (I < N) then
!
!              Generate elementary reflector G(i) to annihilate
!              A(i,i+2:n)
!
        call mobbrmsd_DLARFG(N - I, A(I, I + 1), A(I, MIN(I + 2, N)), LDA, TAUP(I))
        E(I) = A(I, I + 1)
        A(I, I + 1) = ONE
!
!              Apply G(i) to A(i+1:m,i+1:n) from the right
!
        call mobbrmsd_DLARF('Right', M - I, N - I, A(I, I + 1), LDA,&
&                     TAUP(I), A(I + 1, I + 1), LDA, WORK)
        A(I, I + 1) = E(I)
      else
        TAUP(I) = ZERO
      end if
    end do
  else
!
!        Reduce to lower bidiagonal form
!
    do I = 1, M
!
!           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
!
      call mobbrmsd_DLARFG(N - I + 1, A(I, I), A(I, MIN(I + 1, N)), LDA, TAUP(I))
      D(I) = A(I, I)
      A(I, I) = ONE
!
!           Apply G(i) to A(i+1:m,i:n) from the right
!
      if (I < M)&
 &      call mobbrmsd_DLARF('Right', M - I, N - I + 1, A(I, I), LDA,&
 &                 TAUP(I), A(I + 1, I), LDA, WORK)
      A(I, I) = D(I)
!
      if (I < M) then
!
!            Generate elementary reflector H(i) to annihilate
!            A(i+2:m,i)
!
        call mobbrmsd_DLARFG(M - I, A(I + 1, I), A(MIN(I + 2, M), I), 1, TAUQ(I))
        E(I) = A(I + 1, I)
        A(I + 1, I) = ONE
!
!              Apply H(i) to A(i+1:m,i+1:n) from the left
!
        call mobbrmsd_DLARF('Left', M - I, N - I, A(I + 1, I), 1, TAUQ(I),&
&                  A(I + 1, I + 1), LDA, WORK)
        A(I + 1, I) = E(I)
      else
        TAUQ(I) = ZERO
      end if
    end do
  end if
  return
!
!     End of mobbrmsd_DGEBD2
!
end subroutine mobbrmsd_DGEBD2

