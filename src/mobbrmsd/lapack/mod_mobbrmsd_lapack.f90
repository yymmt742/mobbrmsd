!| mod_mobbrmsd_lapack is lapack interface.
!
!  Only routines used for mob calculations are supported.
!  Here is a list of the routines currently required:
!
!  - xGEMM
!
!  - xGESVD
!
!  - xORMQR
!
!  - xGEQRF
!
!  - xGETRF
!
!  If the EXTERNAL_LAPACK macro is enabled,
!  the double precision ones will be replaced
!  with the corresponding routines from the external library.
!
!  The interfaces follow the standard [lapack api](https://www.netlib.org/lapack/).
!
module mod_mobbrmsd_lapack
#ifdef USE_REAL32
  use mod_mobbrmsd_lapack_routines_sp, only: &
    &   SGEMM => mobbrmsd_SGEMM, &
    &   SGESVD => mobbrmsd_SGESVD, &
    &   SORMQR => mobbrmsd_SORMQR, &
    &   SGEQRF => mobbrmsd_SGEQRF, &
    &   SGETRF => mobbrmsd_SGETRF
  implicit none
  private
  public :: SGEMM
  public :: SGESVD
  public :: SORMQR
  public :: SGEQRF
  public :: SGETRF
#else
#ifdef USE_EXTERNAL_LAPACK
!
  integer, parameter   :: RK = SELECTED_REAL_KIND(15)
!
  interface
!| performs \( C \gets \alpha \text{op}( A ) \text{op}( B ) + \beta C \)
!
!  DGEMM performs one of the matrix-matrix operations
!
!  \[ C \gets \alpha \text{op}( A ) \text{op}( B ) + \beta C, \]
!
!  where \( \text{op}( X ) \) is one of
!  \( \text{op}( X ) = X \) or \( \text{op}( X ) = X ^ {\top} \),
!  \( \alpha \) and \( \beta \) are scalars,
!  and \( A \), \( B \) and \( C \) are matrices,
!  with \( \text{op}( A ) \in \mathbb{R} ^ { m \times k } \),
!  \( \text{op}( B ) \in \mathbb{R} ^ { k \times n } \),
!  and \( C \in \mathbb{R} ^ { m \times n } \).
!
    pure subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      import RK
      character, intent(in) :: TRANSA
!!  On entry, TRANSA specifies the form of op( A ) to be used in
!!  the matrix multiplication as follows:
!!
!!     TRANSA = 'N' or 'n',  op( A ) = A.
!!
!!     TRANSA = 'T' or 't',  op( A ) = A**T.
!!
!!     TRANSA = 'C' or 'c',  op( A ) = A**T.
!!
      character, intent(in) :: TRANSB
!!  On entry, TRANSB specifies the form of op( B ) to be used in
!!  the matrix multiplication as follows:
!!
!!     TRANSB = 'N' or 'n',  op( B ) = B.
!!
!!     TRANSB = 'T' or 't',  op( B ) = B**T.
!!
!!     TRANSB = 'C' or 'c',  op( B ) = B**T.
!!
      integer, intent(in)      :: M
!!  On entry,  M  specifies  the number  of rows  of the  matrix
!!  op( A )  and of the  matrix  C.  M  must  be at least  zero.
!!
      integer, intent(in)      :: N
!!  On entry,  N  specifies the number  of columns of the matrix
!!  op( B ) and the number of columns of the matrix C. N must be
!!  at least zero.
!!
      integer, intent(in)      :: K
!!  On entry,  K  specifies  the number of columns of the matrix
!!  op( A ) and the number of rows of the matrix op( B ). K must
!!  be at least  zero.
!!
      real(RK), intent(in)     :: ALPHA
!!  On entry, ALPHA specifies the scalar alpha.
!!
      integer, intent(in)      :: LDA
!!  On entry, LDA specifies the first dimension of A as declared
!!  in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!!  LDA must be at least  max( 1, m ), otherwise  LDA must be at
!!  least  max( 1, k ).
!!
      real(RK), intent(in)     :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
!!  k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!!  Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!!  part of the array  A  must contain the matrix  A,  otherwise
!!  the leading  k by m  part of the array  A  must contain  the
!!  matrix A.
!!
      integer, intent(in)      :: LDB
!!  On entry, LDB specifies the first dimension of B as declared
!!  in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!!  LDB must be at least  max( 1, k ), otherwise  LDB must be at
!!  least  max( 1, n ).
!!
      real(RK), intent(in)     :: B(LDB, *)
!!  DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
!!  n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!!  Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!!  part of the array  B  must contain the matrix  B,  otherwise
!!  the leading  n by k  part of the array  B  must contain  the
!!  matrix B.
!!
      real(RK), intent(in)     :: BETA
!! BETA is DOUBLE PRECISION.
!!  On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!!  supplied as zero then C need not be set on input.
!!
      integer, intent(in)      :: LDC
!!  On entry, LDC specifies the first dimension of C as declared
!!  in  the  calling  (sub)  program.   LDC  must  be  at  least
!!  max( 1, m ).
!!
      real(RK), intent(inout)  :: C(LDC, *)
!!  DOUBLE PRECISION array, dimension ( LDC, N )
!!  Before entry, the leading  m by n  part of the array  C must
!!  contain the matrix  C,  except when  beta  is zero, in which
!!  case C need not be set on entry.
!!  On exit, the array  C  is overwritten by the  m by n  matrix
!!  ( alpha*op( A )*op( B ) + beta*C ).
!!
    end subroutine DGEMM
!
!| mobbrmsd_DGESVD computes the singular value decomposition (SVD) of a real
!  M-by-N matrix A, optionally computing the left and/or right singular
!  vectors. The SVD is written
!
!  \begin{equation}
!     \mathbf{A} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{\top}
!  \end{equation}
!
!  where \(\mathbf{\Sigma}\) is an M-by-N matrix which is zero except for its
!  \(\min(m,n)\) diagonal elements, \(\mathbf{U}\) is an M-by-M orthogonal matrix, and
!  \(\mathbf{V}\) is an N-by-N orthogonal matrix.  The diagonal elements of \(\mathbf{\Sigma}\)
!  are the singular values of \(\mathbf{A}\); they are real and non-negative, and
!  are returned in descending order.  The first \(\min(m,n)\) columns of
!  \(\mathbf{U}\)  and \(\mathbf{V}\) are the left and right singular vectors of \(\mathbf{A}\).
!
!  Note that the routine returns \(\mathbf{V}^{\top}\), not \(\mathbf{V}\).
!
    pure subroutine DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
        &                  VT, LDVT, WORK, LWORK, INFO)
      import RK
      character, intent(in) :: JOBU
!!          Specifies options for computing all or part of the matrix \(\mathbf{U}\)
!!
!!          = 'A':  all M columns of U are returned in array U
!!
!!          = 'S':  the first min(m,n) columns of U (the left singular
!!                  vectors) are returned in the array U.
!!
!!          = 'O':  the first min(m,n) columns of U (the left singular
!!                  vectors) are overwritten on the array A.
!!
!!          = 'N':  no columns of U (no left singular vectors) are
!!                  computed.
!!
      character, intent(in) :: JOBVT
!!          Specifies options for computing all or part of the matrix
!!          \(\mathbf{V}^{\top}\)
!!
!!          = 'A':  all N rows of \(\mathbf{V}^\top\) are returned in the array VT
!!
!!          = 'S':  the first min(m,n) rows of \(\mathbf{V}^\top\) (the right singular
!!                  vectors) are returned in the array VT
!!
!!          = 'O':  the first min(m,n) rows of \(\mathbf{V}^\top\) (the right singular
!!                  vectors) are overwritten on the array A
!!
!!          = 'N':  no rows of \(\mathbf{V}^\top\) (no right singular vectors) are
!!                  computed
!!
!!          JOBVT and JOBU cannot both be 'O'.
!!
      integer, intent(in)   :: LDA
!!          The leading dimension of the array A.  LDA >= max(1,M).
!!
      integer, intent(in)   :: LDU
!!          The leading dimension of the array U.  LDU >= 1; if
!!          JOBU = 'S' or 'A', LDU >= M.
!!
      integer, intent(in)   :: LDVT
!!          The leading dimension of the array VT.  LDVT >= 1; if
!!          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
!!
      integer, intent(in)   :: LWORK
!!          The dimension of the array WORK.
!!
!!          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
!!             - PATH 1  (M much larger than N, JOBU='N')
!!             - PATH 1t (N much larger than M, JOBVT='N')
!!
!!          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths
!!
!!          For good performance, LWORK should generally be larger.
!!
!!          If LWORK = -1, then a workspace query is assumed; the routine
!!          only calculates the optimal size of the WORK array, returns
!!          this value as the first entry of the WORK array.
!!
      integer, intent(in)   :: M
!!          The number of rows of the input matrix A.  M >= 0.
!!
      integer, intent(in)   :: N
!!          The number of columns of the input matrix A.  N >= 0.
!!
      real(RK), intent(inout)   :: A(LDA, *)
!!          A is REAL array, dimension (LDA,N)
!!
!!          On entry, the M-by-N matrix \(\mathbf{A}\).
!!          On exit,
!!
!!          if JOBU = 'O',  \(\mathbf{A}\) is overwritten with the first min(m,n)
!!                          columns of U (the left singular vectors,
!!                          stored columnwise).
!!
!!          if JOBVT = 'O', \(\mathbf{A}\) is overwritten with the first min(m,n)
!!                          rows of V**T (the right singular vectors,
!!                          stored rowwise).
!!
!!          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of \(\mathbf{A}\)
!!                          are destroyed.
!!
      real(RK), intent(out)     :: S(*)
!!          S is REAL array, dimension (min(M,N))
!!
!!          The singular values of \(\mathbf{A}\), sorted so that S(i) >= S(i+1).
!!
      real(RK), intent(out)     :: U(LDU, *)
!!          U is REAL array, dimension (LDU,UCOL)
!!
!!          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!!
!!          If JOBU = 'A', U contains the M-by-M orthogonal matrix \(\mathbf{U}\).
!!
!!          if JOBU = 'S', U contains the first min(m,n) columns of \(\mathbf{U}\)
!!          (the left singular vectors, stored columnwise).
!!
!!          if JOBU = 'N' or 'O', U is not referenced.
!!
      real(RK), intent(out)     :: VT(LDVT, *)
!!          VT is REAL array, dimension (LDVT,N)
!!
!!          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
!!           \(\mathbf{V}^\top\);
!!
!!          if JOBVT = 'S', VT contains the first min(m,n) rows of
!!           \(\mathbf{V}^\top\)
!!           (the right singular vectors, stored rowwise);
!!
!!          if JOBVT = 'N' or 'O', VT is not referenced.
!!
      real(RK), intent(out)     :: WORK(*)
!!          WORK is REAL array, dimension (MAX(1,LWORK))
!!
!!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!!
!!          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!!          superdiagonal elements of an upper bidiagonal matrix \(\mathbf{B}\)
!!          whose diagonal is in S (not necessarily sorted). \(\mathbf{B}\)
!!          satisfies \(\mathbf{A} = \mathbf{U} \mathbf{B} \mathbf{V}^{\top}\),
!!          so it has the same singular values as \(\mathbf{A}\),
!!          and singular vectors related by \(\mathbf{U}\) and \(\mathbf{V}\).
!!
      integer, intent(out)  :: INFO
!!          =0 :  successful exit.
!!
!!          <0 :  if INFO = -i, the i-th argument had an illegal value.
!!
!!          \>0 : if mobbrmsd_SBDSQR did not converge, INFO specifies how many
!!                superdiagonals of an intermediate bidiagonal form \(\mathbf{B}\)
!!                did not converge to zero. See the description of WORK
!!                above for details.
!!
    end subroutine DGESVD
!
!| multiply an orthogonal matrix \( Q = H _ 1 H _ 2 \cdots H _ k \).
!
!  DORMQR overwrites the general real
!  \( m \)-by-\( n \) matrix \( C \) with
!
!  | SIDE  | TRANS |                  |
!  | :---: | :---: |   :---:          |
!  |  'L'  |  'N'  | \( Q C \)        |
!  |  'R'  |  'N'  | \( C Q \)        |
!  |  'L'  |  'T'  | \( Q ^ \top C \) |
!  |  'R'  |  'T'  | \( C Q ^ \top \) |
!
!  where \( Q \) is a real orthogonal matrix defined
!  as the product of \( k \) elementary reflectors
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ k
!  \]
!
!  as returned by DORMQR.
!  \( Q \) is of order \( m \) if SIDE = 'L'
!  and of order \( n \) if SIDE = 'R'.
!
    pure subroutine DORMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, &
                   &                C, LDC, WORK, LWORK, INFO)
      import RK
      character, intent(in) :: SIDE
!!  = 'L': apply Q or Q**T from the Left;
!!
!!  = 'R': apply Q or Q**T from the Right.
!!
      character, intent(in) :: TRANS
!!  = 'N':  No transpose, apply Q;
!!
!!  = 'T':  Transpose, apply Q**T.
!!
      integer, intent(in)   :: M
!!  The number of rows of the matrix C. M >= 0.
!!
      integer, intent(in)   :: N
!!  The number of columns of the matrix C. N >= 0.
!!
      integer, intent(in)   :: K
!!  The number of elementary reflectors whose product defines
!!  the matrix Q.
!!
!!  If SIDE = 'L', M >= K >= 0;
!!
!!  if SIDE = 'R', N >= K >= 0.
!!
      integer, intent(in)   :: LDA
!!  The leading dimension of the array A.
!!
!!  If SIDE = 'L', LDA >= max(1,M);
!!
!!  if SIDE = 'R', LDA >= max(1,N).
!!
      real(RK), intent(inout)  :: A(LDA, *)
!!  DOUBLE PRECISION array, dimension (LDA,K)
!!
!!  The i-th column must contain the vector which defines the
!!  elementary reflector H(i), for i = 1,2,...,k, as returned by
!!  mobbrmsd_DGEQRF in the first k columns of its array argument A.
!!
      integer, intent(in)      :: LDC
!!  The leading dimension of the array C. LDC >= max(1,M).
!!
      real(RK), intent(in)     :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i), as returned by mobbrmsd_DGEQRF.
!!
      real(RK), intent(inout)  :: C(LDC, *)
!!  DOUBLE PRECISION array, dimension (LDC,N)
!!
!!  On entry, the M-by-N matrix C.
!!
!!  On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!!
      real(RK), intent(out)    :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
      integer, intent(in)      :: LWORK
!!  The dimension of the array WORK.
!!
!!  If SIDE = 'L', LWORK >= max(1,N);
!!
!!  if SIDE = 'R', LWORK >= max(1,M).
!!
!!  For good performance, LWORK should generally be larger.
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array, and no error
!!  message related to LWORK is issued by XERBLA.
!!
      integer, intent(out)     :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
    end subroutine DORMQR
!
!| computes a QR factorization of a real \( m \)-by-\( n \) matrix \( A \):
!
!  DGEQRF computes a QR factorization
!  of a real \( m \)-by-\( n \) matrix \( A \):
!
!  \[
!    A = Q
!    \left(
!      \begin{array}{}
!        R \\
!        0 \\
!      \end{array}
!    \right)
!  \]
!  where:
!
!  \( Q \) is a \( m \)-by-\( m \) orthogonal matrix;
!  \( R \) is an upper-triangular \( n \)-by-\( n \) matrix;
!  \( 0 \) is a \( (m-n) \)-by-\( n \) zero matrix, if \( m > n \).
!
!  The matrix \( Q \) is represented as a product of elementary reflectors
!
!  \[
!     Q = H _ 1 H _ 2 \cdots H _ k,
!  \]
!
!   where \( k = \min(m,n) \).
!   Each  \( H _ i \) has the form
!
!  \[
!     H _ i = I - \tau v v ^ \top,
!  \]
!
!   where \( \tau \) is a real scalar,
!   and \( v \) is a real vector with
!   \( v _ j = 0, 1 \le j \le i-1) \) and \( v _ i = 1 \);
!   v(i+1:m) is stored on exit in A(i+1:m,i), and \( \tau \) in TAU(i).
!
    pure subroutine DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
      import RK
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
!!  On entry, the M-by-N matrix A.
!!
!!  On exit, the elements on and above the diagonal of the array
!!  contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!!  upper triangular if m >= n); the elements below the diagonal,
!!  with the array TAU, represent the orthogonal matrix Q as a
!!  product of min(m,n) elementary reflectors (see Further
!!  Details).
!!
      real(RK), intent(out)   :: TAU(*)
!!  DOUBLE PRECISION array, dimension (min(M,N))
!!  The scalar factors of the elementary reflectors (see Further
!!  Details).
!!
      real(RK), intent(out)   :: WORK(*)
!!  DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!
      integer, intent(in)     :: LWORK
!!  The dimension of the array WORK.
!!
!!  LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
!!  For optimum performance LWORK >= N*NB, where NB is
!!  the optimal blocksize.
!!
!!  If LWORK = -1, then a workspace query is assumed; the routine
!!  only calculates the optimal size of the WORK array, returns
!!  this value as the first entry of the WORK array.
!!
      integer, intent(out)    :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
    end subroutine DGEQRF
!
!| computes an LU factorization using partial pivoting.
!
!  mobbrmsd_DGETRF computes an LU factorization of a general
!  \( m \)-by-\( n \) matrix \( A \)
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!
!  \[ A = P L U \]
!
!  where \( P \) is a permutation matrix,
!  \( L \) is lower triangular with unit diagonal elements
!  (lower trapezoidal if \( m > n \) ),
!  and \( U \) is upper triangular
!  (upper trapezoidal if \( m < n \) ).
!
    pure subroutine DGETRF(M, N, A, LDA, IPIV, INFO)
      import RK
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
!!  A is DOUBLE PRECISION array, dimension (LDA,N)
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
    end subroutine DGETRF
  end interface
#else
  use mod_mobbrmsd_lapack_routines_dp, only: &
    &   DGEMM => mobbrmsd_DGEMM, &
    &   DGESVD => mobbrmsd_DGESVD, &
    &   DORMQR => mobbrmsd_DORMQR, &
    &   DGEQRF => mobbrmsd_DGEQRF, &
    &   DGETRF => mobbrmsd_DGETRF
#endif
  public :: DGEMM
  public :: DGESVD
  public :: DORMQR
  public :: DGEQRF
  public :: DGETRF
#endif
end module mod_mobbrmsd_lapack

