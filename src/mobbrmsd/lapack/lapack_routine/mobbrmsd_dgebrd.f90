!| mobbrmsd_DGEBRD reduces a general real M-by-N matrix A to upper or lower
!  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!
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
!  reference DBDSQR is provided by http://www.netlib.org/lapack/explore-html/
!  \author Univ. of Tennessee
!  \author Univ. of California Berkeley
!  \author Univ. of Colorado Denver
!  \author NAG Ltd.
pure subroutine mobbrmsd_DGEBRD(M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO)
  implicit none
  integer, intent(in)     :: M
!!          The number of rows in the matrix A.  M >= 0.
  integer, intent(in)     :: N
!!          The number of columns in the matrix A.  N >= 0.
  integer, intent(in)     :: LDA
!!          The leading dimension of the array A.  LDA >= max(1,M).
  real(RK), intent(inout) :: A(LDA, *)
!!          On entry, the M-by-N general matrix to be reduced.
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
!!            a product of elementary reflectors.
  real(RK), intent(out)   :: D(*)
!!          The diagonal elements of the bidiagonal matrix B:
!!          D(i) = A(i,i).
  real(RK), intent(out)   :: E(*)
!!          The off-diagonal elements of the bidiagonal matrix B:
!!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  real(RK), intent(out)   :: TAUP(*)
!!          TAUP is real(RK)           :: array, dimension (min(M,N))
!!          The scalar factors of the elementary reflectors which
!!          represent the orthogonal matrix P. See Further Details.
  real(RK), intent(out)   :: TAUQ(*)
!!          The scalar factors of the elementary reflectors which
!!          represent the orthogonal matrix Q. See Further Details.
  integer, intent(in)     :: LWORK
!!          The length of the array WORK.  LWORK >= max(1,M,N).
!!          For optimum performance LWORK >= (M+N)*NB, where NB
!!          is the optimal blocksize.
!!
!!          If LWORK = -1, then a workspace query is assumed; the routine
!!          only calculates the optimal size of the WORK array, returns
!!          this value as the first entry of the WORK array, and no error
!!          message related to LWORK is issued by XERBLA.
  real(RK), intent(out)   :: WORK(*)
!!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  integer, intent(out)    :: INFO
!!          = 0:  successful exit
!
  logical                :: LQUERY
  integer                :: I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, &
 &                          MINMN, NB, NBMIN, NX, WS
  intrinsic              :: DBLE, MAX, MIN
!     ..
!     .. Parameters ..
! real(RK), parameter     ::   ONE = 1.0_RK
!     ..
! interface
!     .. External Subroutines ..
!   include 'dgebd2.h'
!   include 'dgemm.h'
!   include 'dlabrd.h'
!     .. External Functions ..
!   include 'ilaenv.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  NB = MAX(1, mobbrmsd_ILAENV(1, 'DGEBRD', ' ', M, N, -1, -1))
  LWKOPT = (M + N) * NB
  WORK(1) = DBLE(LWKOPT)
  LQUERY = (LWORK == -1)
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  else if (LWORK < MAX(1, M, N) .and. .not. LQUERY) then
    INFO = -10
  end if
  if (INFO < 0) then
!   CALL XERBLA( 'DGEBRD', -INFO )
    return
  else if (LQUERY) then
    return
  end if
!
!     Quick return if possible
!
  MINMN = MIN(M, N)
  if (MINMN == 0) then
    WORK(1) = 1
    return
  end if
!
  WS = MAX(M, N)
  LDWRKX = M
  LDWRKY = N
!
  if (NB > 1 .and. NB < MINMN) then
!
!        Set the crossover point NX.
!
    NX = MAX(NB, mobbrmsd_ILAENV(3, 'DGEBRD', ' ', M, N, -1, -1))
!
!        Determine when to switch from blocked to unblocked code.
!
    if (NX < MINMN) then
      WS = (M + N) * NB
      if (LWORK < WS) then
!
!              Not enough work space for the optimal NB, consider using
!              a smaller block size.
!
        NBMIN = mobbrmsd_ILAENV(2, 'DGEBRD', ' ', M, N, -1, -1)
        if (LWORK >= (M + N) * NBMIN) then
          NB = LWORK / (M + N)
        else
          NB = 1
          NX = MINMN
        end if
      end if
    end if
  else
    NX = MINMN
  end if
!
  do I = 1, MINMN - NX, NB
!
!        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
!        the matrices X and Y which are needed to update the unreduced
!        part of the matrix
!
    call mobbrmsd_DLABRD(M - I + 1, N - I + 1, NB, A(I, I), LDA, D(I), E(I), &
        &                TAUQ(I), TAUP(I), WORK, LDWRKX, &
        &                WORK(LDWRKX * NB + 1), LDWRKY)
!
!        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
!        of the form  A := A - V*Y**T - X*U**T
!
    call mobbrmsd_DGEMM('No transpose', 'Transpose', M - I - NB + 1, N - I - NB + 1,&
        &               NB, -ONE, A(I + NB, I), LDA,&
        &               WORK(LDWRKX * NB + NB + 1), LDWRKY, ONE,&
        &               A(I + NB, I + NB), LDA)
    call mobbrmsd_DGEMM('No transpose', 'No transpose', M - I - NB + 1, N - I - NB + 1,&
        &               NB, -ONE, WORK(NB + 1), LDWRKX, A(I, I + NB), LDA,&
        &               ONE, A(I + NB, I + NB), LDA)
!
!        Copy diagonal and off-diagonal elements of B back into A
!
    if (M >= N) then
      do J = I, I + NB - 1
        A(J, J) = D(J)
        A(J, J + 1) = E(J)
      end do
    else
      do J = I, I + NB - 1
        A(J, J) = D(J)
        A(J + 1, J) = E(J)
      end do
    end if
  end do
!
!     Use unblocked code to reduce the remainder of the matrix
!
  call mobbrmsd_DGEBD2(M - I + 1, N - I + 1, A(I, I), LDA, D(I), E(I),&
      &                TAUQ(I), TAUP(I), WORK, IINFO)
  WORK(1) = WS
  return
!
!     End of mobbrmsd_DGEBRD
!
end subroutine mobbrmsd_DGEBRD

