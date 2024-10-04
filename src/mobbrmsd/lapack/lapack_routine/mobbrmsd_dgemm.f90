!> \brief \b mobbrmsd_DGEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension ( LDC, N )
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, m ).
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
!> \ingroup double_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  implicit none
! use LA_CONSTANTS, only: RK => dp
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in)     :: ALPHA, BETA
  integer, intent(in)      :: K, LDA, LDB, LDC, M, N
  character(*), intent(in) :: TRANSA, TRANSB
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: A(LDA, *), B(LDB, *)
  real(RK), intent(inout)  :: C(LDC, *)
!     ..
!
!  =====================================================================
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: MAX
!     ..
!     .. Local Scalars ..
  real(RK) :: TEMP
  integer  :: I, INFO, J, L, NROWA, NROWB
  logical  :: NOTA, NOTB
!     ..
!     .. Parameters ..
! real(RK), parameter      :: ONE = 1.0_RK
! real(RK), parameter      :: ZERO = 0.0_RK
!
! interface
!   include 'lsame.h'
!   !include 'xerbla.h'
! end interface
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA and NROWB  as the number of rows of  A
!     and  B  respectively.
!
  NOTA = mobbrmsd_LSAME(TRANSA, 'N')
  NOTB = mobbrmsd_LSAME(TRANSB, 'N')
  if (NOTA) then
    NROWA = M
  else
    NROWA = K
  end if
  if (NOTB) then
    NROWB = K
  else
    NROWB = N
  end if
!
!     Test the input parameters.
!
  INFO = 0
  if ((.not. NOTA) .and. (.not. mobbrmsd_LSAME(TRANSA, 'C')) .and. &
 &    (.not. mobbrmsd_LSAME(TRANSA, 'T'))) then
    INFO = 1
  else if ((.not. NOTB) .and. (.not. mobbrmsd_LSAME(TRANSB, 'C')) .and. &
 &         (.not. mobbrmsd_LSAME(TRANSB, 'T'))) then
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
    !CALL XERBLA('DGEMM ',INFO)
    return
  end if
!
!     Quick return if possible.
!
  if ((M == 0) .or. (N == 0) .or. &
 &    (((ALPHA == ZERO) .or. (K == 0)) .and. (BETA == ONE))) return
!
!     And if  alpha.eq.zero.
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
!     Start the operations.
!
  if (NOTB) then
    if (NOTA) then
!
!         Form  C := alpha*A*B + beta*C.
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
!           Form  C := alpha*A**T*B + beta*C
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
!           Form  C := alpha*A*B**T + beta*C
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
!           Form  C := alpha*A**T*B**T + beta*C
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
!     End of mobbrmsd_DGEMM
!
end subroutine mobbrmsd_DGEMM

