!| mobbrmsd_SGEMM performs one of the matrix-matrix operations
!
!  \[ \mathbf{C} \gets \alpha \text{op}( \mathbf{A} ) \text{op}( \mathbf{B} ) + \beta \mathbf{C}, \]
!
!  where \( \text{op}( \mathbf{X} ) \) is one of
!  \( \text{op}( \mathbf{X} ) = \mathbf{X} \) or \( \text{op}( \mathbf{X} ) = \mathbf{X} ^ {\top} \),
!  \( \alpha \) and \( \beta \) are scalars, and \( \mathbf{A} \), \( \mathbf{B} \) and \( \mathbf{C} \) are matrices,
!  with \( \text{op}( \mathbf{A} ) \) an \( m \) by \( k \) matrix,
!  \( \text{op}( \mathbf{B} ) \) a \( k \) by \( n \) matrix,
!  and \( \mathbf{B} \) an \( m \) by \( n \) matrix.
!
!  Reference SGEMM is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
pure subroutine mobbrmsd_SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  implicit none
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
  intrinsic :: MAX
  real(RK) :: TEMP
  integer  :: I, INFO, J, L, NCOLA, NROWA, NROWB
  logical  :: NOTA, NOTB
!
! Set NOTA and NOTB as true if A and B respectively are not
! transposed and set NROWA, NCOLA and NROWB as the number of rows
! and columns of A and the number of rows of B respectively.
!
  NOTA = mobbrmsd_LSAME(TRANSA, 'N')
  NOTB = mobbrmsd_LSAME(TRANSB, 'N')
  if (NOTA) then
    NROWA = M
    NCOLA = K
  else
    NROWA = K
    NCOLA = M
  end if
  if (NOTB) then
    NROWB = K
  else
    NROWB = N
  end if
!
! Test the input parameters.
!
  INFO = 0
  if ((.not. NOTA) .and. (.not. mobbrmsd_LSAME(TRANSA, 'C')) .and. (.not. mobbrmsd_LSAME(TRANSA, 'T'))) then
    INFO = 1
  else if ((.not. NOTB) .and. (.not. mobbrmsd_LSAME(TRANSB, 'C')) .and. (.not. mobbrmsd_LSAME(TRANSB, 'T'))) then
    INFO = 2
  else if (M < 0) then
    INFO = 3
  else if (N < 0) then
    INFO = 4
  else if (K < 0) then
    INFO = 5
  else if (LDA < MAX(1, NROWA)) then
    INFO = 8
  else if (LDB < MAX(1, NROWB)) then
    INFO = 10
  else if (LDC < MAX(1, M)) then
    INFO = 13
  end if
  if (INFO /= 0) then
!call XERBLA('SGEMM ', INFO)
    return
  end if
!
! Quick return if possible.
!
  if ((M == 0) .or. (N == 0) .or. (((ALPHA == ZERO) .or. (K == 0)) .and. (BETA == ONE))) return
!
! And if alpha == zero.
!
  if (ALPHA == ZERO) then
    if (BETA == ZERO) then
      do J = 1, N
        do I = 1, M
          C(I, J) = ZERO
        end do
      end do
    else
      do J = 1, N
        do I = 1, M
          C(I, J) = BETA * C(I, J)
        end do
      end do
    end if
    return
  end if
!
!Start the operations.
!
  if (NOTB) then
    if (NOTA) then
!
!Form C: = alpha * A * B + beta * C.
!
      do J = 1, N
        if (BETA == ZERO) then
          do I = 1, M
            C(I, J) = ZERO
          end do
        else if (BETA /= ONE) then
          do I = 1, M
            C(I, J) = BETA * C(I, J)
          end do
        end if
        do L = 1, K
          TEMP = ALPHA * B(L, J)
          do I = 1, M
            C(I, J) = C(I, J) + TEMP * A(I, L)
          end do
        end do
      end do
    else
!
!Form C: = alpha * A**T * B + beta * C
!
      do J = 1, N
        do I = 1, M
          TEMP = ZERO
          do L = 1, K
            TEMP = TEMP + A(L, I) * B(L, J)
          end do
          if (BETA == ZERO) then
            C(I, J) = ALPHA * TEMP
          else
            C(I, J) = ALPHA * TEMP + BETA * C(I, J)
          end if
        end do
      end do
    end if
  else
    if (NOTA) then
!
! Form C: = alpha * A * B**T + beta * C
!
      do J = 1, N
        if (BETA == ZERO) then
          do I = 1, M
            C(I, J) = ZERO
          end do
        else if (BETA /= ONE) then
          do I = 1, M
            C(I, J) = BETA * C(I, J)
          end do
        end if
        do L = 1, K
          TEMP = ALPHA * B(J, L)
          do I = 1, M
            C(I, J) = C(I, J) + TEMP * A(I, L)
          end do
        end do
      end do
    else
!
! Form C: = alpha * A**T * B**T + beta * C
!
      do J = 1, N
        do I = 1, M
          TEMP = ZERO
          do L = 1, K
            TEMP = TEMP + A(L, I) * B(J, L)
          end do
          if (BETA == ZERO) then
            C(I, J) = ALPHA * TEMP
          else
            C(I, J) = ALPHA * TEMP + BETA * C(I, J)
          end if
        end do
      end do
    end if
  end if
!
  return
!
! end of mobbrmsd_SGEMM.
!
end

