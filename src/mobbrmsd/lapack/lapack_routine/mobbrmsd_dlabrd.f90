!| mobbrmsd_DLABRD reduces the first \( nb \)  rows and columns of
!  a general matrix to a bidiagonal form.
!
!  mobbrmsd_DLABRD reduces the first \( NB \) rows and columns of
!  a real general \( m \) by \( n \) matrix \( A \) to upper
!  or lower bidiagonal form by an orthogonal transformation
!  \( Q ^ \top A P \), and returns the matrices \( X \) and \( Y \)
!  which are needed to apply the transformation to
!  the unreduced part of \( A \).
!
!  If \( m >= n \), \( A \) is reduced to upper bidiagonal form;
!  if \( m < n \), to lower  bidiagonal form.
!
!  This is an auxiliary routine called by mobbrmsd_DGEBRD
!
!  The matrices \( Q \) and \( P \) are represented
!  as products of elementary reflectors:
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ {n_b}
!  \]
!
!  and
!
!  \[
!     P = G _ 1 G _ 2 \cdots G _ {n_b}
!  \]
!
!  Each \( H _ i \) and \( G _ i \) has the form:
!
!  \[ H _ i = I - \tau _ q * v * v^\top \]
!
!  and
!
!  \[ G _ i = I - \tau _ p * u * u^\top \]
!
!  where \( \tau _ q \) and \( \tau _ p \) are real scalars,
!  and \( v \) and \( u \) are real vectors.
!
!  If \( m \ge n \) , v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored
!  on exit in A(i:m,i); u(1:i) = 0, u(i+1) = 1,
!  and u(i+1:n) is stored on exit in A(i,i+1:n);
!  \( \tau _ q \) is stored in TAUQ(i) and \( \tau _ p \) in TAUP(i).
!
!  If \( m < n \), v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored
!  on exit in A(i+2:m,i); u(1:i-1) = 0, u(i) = 1,
!  and u(i:n) is stored on exit in A(i,i+1:n);
!  \( \tau _ q \) is stored in TAUQ(i) and \( \tau _ p \) in TAUP(i).
!
!  The elements of the vectors \( v \) and \( u \) together form
!  the \( m \)-by-\( nb \), matrix \( V \),
!  and the nb-by-n matrix \( nb = 2 U ^ \top \) which are needed,
!  with X and Y, to apply
!  the transformation to the unreduced part of the matrix,
!  using a block update of the form:
!  \( A := A - V Y ^ {\top} - X U ^ {\top} \).
!
!  The contents of \( A \) on exit are illustrated
!  by the following examples with \( nb = 2 \):
!
!  \( m = 6 \) and \( n = 5 \), \( (m > n) \):
!
!  \[
!     \left (
!       \begin{array}{}
!         1   & 1   & u1 & u1  &  u1  \\
!         v1  & 1   & 1  & u2  &  u2  \\
!         v1  & v2  & a  & a   &  a   \\
!         v1  & v2  & a  & a   &  a   \\
!         v1  & v2  & a  & a   &  a   \\
!         v1  & v2  & a  & a   &  a   \\
!       \end{array}{}
!     \right )
!  \]
!
!  \( m = 5 \) and \( n = 6 \), \( (m < n) \):
!
!  \[
!     \left (
!       \begin{array}{}
!         1  &  u1 &  u1 &  u1 &  u1 &  &  u1 \\
!         1  &  1  &  u2 &  u2 &  u2 &  &  u2 \\
!         v1 &  1  &  a  &  a  &  a  &  &  a  \\
!         v1 &  v2 &  a  &  a  &  a  &  &  a  \\
!         v1 &  v2 &  a  &  a  &  a  &  &  a  \\
!       \end{array}{}
!     \right )
!  \]
!
!  Reference DLABRD is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  where a denotes an element of the original matrix which is unchanged,
!  vi denotes an element of the vector defining H(i), and ui an element
!  of the vector defining G(i).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLABRD(M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY)
  implicit none
  integer, intent(in)     :: M
!!  The number of rows in the matrix A.
!!
  integer, intent(in)     :: N
!!  The number of columns in the matrix A.
!!
  integer, intent(in)     :: NB
!!  The number of leading rows and columns of A to be reduced.
!!
  integer, intent(in)     :: LDA
!!  The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!  REAL array, dimension (LDA,N)
!!
!!  On entry, the m by n general matrix to be reduced.
!!
!!  On exit, the first NB rows and columns of the matrix are
!!  overwritten; the rest of the array is unchanged.
!!
!!  If m >= n, elements on and below the diagonal in the first NB
!!    columns, with the array TAUQ, represent the orthogonal
!!    matrix Q as a product of elementary reflectors; and
!!    elements above the diagonal in the first NB rows, with the
!!    array TAUP, represent the orthogonal matrix P as a product
!!    of elementary reflectors.
!!
!!  If m < n, elements below the diagonal in the first NB
!!    columns, with the array TAUQ, represent the orthogonal
!!    matrix Q as a product of elementary reflectors, and
!!    elements on and above the diagonal in the first NB rows,
!!    with the array TAUP, represent the orthogonal matrix P as
!!    a product of elementary reflectors.
!!
  real(RK), intent(out)   :: D(*)
!!  REAL array, dimension (NB)
!!  The diagonal elements of the first NB rows and columns of
!!  the reduced matrix.  D(i) = A(i,i).
!!
  real(RK), intent(out)   :: E(*)
!!  REAL array, dimension (NB)
!!  The off-diagonal elements of the first NB rows and columns of
!!  the reduced matrix.
!!
  real(RK), intent(out)   :: TAUQ(*)
!!  REAL array, dimension (NB)
!!  The scalar factors of the elementary reflectors which
!!  represent the orthogonal matrix Q. See Further Details.
!!
  real(RK), intent(out)   :: TAUP(*)
!!  TAUP is REAL array, dimension (NB)
!!  The scalar factors of the elementary reflectors which
!!  represent the orthogonal matrix P. See Further Details.
!!
  integer, intent(in)     :: LDX
!!  The leading dimension of the array X. LDX >= max(1,M).
!!
  real(RK), intent(out)   :: X(LDX, *)
!!  REAL array, dimension (LDX,NB)
!!  The m-by-nb matrix X required to update the unreduced part
!!  of A.
!!
  integer, intent(in)     :: LDY
!!  The leading dimension of the array Y. LDY >= max(1,N).
!!
  real(RK), intent(out)   :: Y(LDY, *)
  integer                :: I
  intrinsic              :: MIN
! interface
!   include 'dgemv.h'
!   include 'dlarfg.h'
!   include 'dscal.h'
! end interface
!
! Quick return if possible
!
  if (M <= 0 .or. N <= 0) return
!
  if (M >= N) then
!
!        Reduce to upper bidiagonal form
!
    do I = 1, NB
!
!           Update A(i:m,i)
!
      call mobbrmsd_DGEMV('No transpose', M - I + 1, I - 1, -ONE, A(I, 1),&
     &            LDA, Y(I, 1), LDY, ONE, A(I, I), 1)
      call mobbrmsd_DGEMV('No transpose', M - I + 1, I - 1, -ONE, X(I, 1),&
     &            LDX, A(1, I), 1, ONE, A(I, I), 1)
!
!           Generate reflection Q(i) to annihilate A(i+1:m,i)
!
      call mobbrmsd_DLARFG(M - I + 1, A(I, I), A(MIN(I + 1, M), I), 1,&
     &             TAUQ(I))
      D(I) = A(I, I)
      if (I < N) then
        A(I, I) = ONE
!
!              Compute Y(i+1:n,i)
!
        call mobbrmsd_DGEMV('Transpose', M - I + 1, N - I, ONE, A(I, I + 1),&
       &            LDA, A(I, I), 1, ZERO, Y(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', M - I + 1, I - 1, ONE, A(I, 1), LDA,&
       &            A(I, I), 1, ZERO, Y(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', N - I, I - 1, -ONE, Y(I + 1, 1),&
       &            LDY, Y(1, I), 1, ONE, Y(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', M - I + 1, I - 1, ONE, X(I, 1), LDX,&
       &            A(I, I), 1, ZERO, Y(1, I), 1)
        call mobbrmsd_DGEMV('Transpose', I - 1, N - I, -ONE, A(1, I + 1),&
       &            LDA, Y(1, I), 1, ONE, Y(I + 1, I), 1)
        call mobbrmsd_DSCAL(N - I, TAUQ(I), Y(I + 1, I), 1)
!
!              Update A(i,i+1:n)
!
        call mobbrmsd_DGEMV('No transpose', N - I, I, -ONE, Y(I + 1, 1),&
       &            LDY, A(I, 1), LDA, ONE, A(I, I + 1), LDA)
        call mobbrmsd_DGEMV('Transpose', I - 1, N - I, -ONE, A(1, I + 1),&
       &            LDA, X(I, 1), LDX, ONE, A(I, I + 1), LDA)
!
!              Generate reflection P(i) to annihilate A(i,i+2:n)
!
        call mobbrmsd_DLARFG(N - I, A(I, I + 1), A(I, MIN(I + 2, N)),&
       &             LDA, TAUP(I))
        E(I) = A(I, I + 1)
        A(I, I + 1) = ONE
!
!              Compute X(i+1:m,i)
!
        call mobbrmsd_DGEMV('No transpose', M - I, N - I, ONE, A(I + 1, I + 1),&
      &             LDA, A(I, I + 1), LDA, ZERO, X(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', N - I, I, ONE, Y(I + 1, 1), LDY,&
       &            A(I, I + 1), LDA, ZERO, X(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', M - I, I, -ONE, A(I + 1, 1),&
       &            LDA, X(1, I), 1, ONE, X(I + 1, I), 1)
        call mobbrmsd_DGEMV('No transpose', I - 1, N - I, ONE, A(1, I + 1),&
       &            LDA, A(I, I + 1), LDA, ZERO, X(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', M - I, I - 1, -ONE, X(I + 1, 1),&
       &            LDX, X(1, I), 1, ONE, X(I + 1, I), 1)
        call mobbrmsd_DSCAL(M - I, TAUP(I), X(I + 1, I), 1)
      end if
    end do
  else
!
!      Reduce to lower bidiagonal form
!
    do I = 1, NB
!
!         Update A(i,i:n)
!
      call mobbrmsd_DGEMV('No transpose', N - I + 1, I - 1, -ONE, Y(I, 1),&
     &            LDY, A(I, 1), LDA, ONE, A(I, I), LDA)
      call mobbrmsd_DGEMV('Transpose', I - 1, N - I + 1, -ONE, A(1, I), LDA,&
     &            X(I, 1), LDX, ONE, A(I, I), LDA)
!
!         Generate reflection P(i) to annihilate A(i,i+1:n)
!
      call mobbrmsd_DLARFG(N - I + 1, A(I, I), A(I, MIN(I + 1, N)), LDA,&
     &            TAUP(I))
      D(I) = A(I, I)
      if (I < M) then
        A(I, I) = ONE
!
!            Compute X(i+1:m,i)
!
        call mobbrmsd_DGEMV('No transpose', M - I, N - I + 1, ONE, A(I + 1, I),&
       &            LDA, A(I, I), LDA, ZERO, X(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', N - I + 1, I - 1, ONE, Y(I, 1), LDY,&
       &            A(I, I), LDA, ZERO, X(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', M - I, I - 1, -ONE, A(I + 1, 1),&
       &            LDA, X(1, I), 1, ONE, X(I + 1, I), 1)
        call mobbrmsd_DGEMV('No transpose', I - 1, N - I + 1, ONE, A(1, I),&
       &            LDA, A(I, I), LDA, ZERO, X(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', M - I, I - 1, -ONE, X(I + 1, 1),&
       &            LDX, X(1, I), 1, ONE, X(I + 1, I), 1)
        call mobbrmsd_DSCAL(M - I, TAUP(I), X(I + 1, I), 1)
!
!            Update A(i+1:m,i)
!
        call mobbrmsd_DGEMV('No transpose', M - I, I - 1, -ONE, A(I + 1, 1),&
       &            LDA, Y(I, 1), LDY, ONE, A(I + 1, I), 1)
        call mobbrmsd_DGEMV('No transpose', M - I, I, -ONE, X(I + 1, 1),&
       &            LDX, A(1, I), 1, ONE, A(I + 1, I), 1)
!
!            Generate reflection Q(i) to annihilate A(i+2:m,i)
!
        call mobbrmsd_DLARFG(M - I, A(I + 1, I), A(MIN(I + 2, M), I), 1,&
       &            TAUQ(I))
        E(I) = A(I + 1, I)
        A(I + 1, I) = ONE
!
!            Compute Y(i+1:n,i)
!
        call mobbrmsd_DGEMV('Transpose', M - I, N - I, ONE, A(I + 1, I + 1),&
       &            LDA, A(I + 1, I), 1, ZERO, Y(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', M - I, I - 1, ONE, A(I + 1, 1), LDA,&
       &            A(I + 1, I), 1, ZERO, Y(1, I), 1)
        call mobbrmsd_DGEMV('No transpose', N - I, I - 1, -ONE, Y(I + 1, 1),&
       &            LDY, Y(1, I), 1, ONE, Y(I + 1, I), 1)
        call mobbrmsd_DGEMV('Transpose', M - I, I, ONE, X(I + 1, 1), LDX,&
       &            A(I + 1, I), 1, ZERO, Y(1, I), 1)
        call mobbrmsd_DGEMV('Transpose', I, N - I, -ONE, A(1, I + 1), LDA,&
       &            Y(1, I), 1, ONE, Y(I + 1, I), 1)
        call mobbrmsd_DSCAL(N - I, TAUQ(I), Y(I + 1, I), 1)
      end if
    end do
  end if
!
!     End of mobbrmsd_DLABRD
!
end subroutine mobbrmsd_DLABRD
