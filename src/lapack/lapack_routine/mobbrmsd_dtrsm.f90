!> \brief \b mobbrmsd_DTRSM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       pure subroutine mobbrmsd_DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       real(RK)           :: ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       real(RK)           :: A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DTRSM  solves one of the matrix equations
!>
!>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!>
!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!>
!> The matrix X is overwritten on B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry, SIDE specifies whether op( A ) appears on the left
!>           or right of X as follows:
!>
!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!>
!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is real(RK)           ::.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is real(RK)           :: array, dimension ( LDA, k ),
!>           where k is m when SIDE = 'L' or 'l'
!>             and k is n when SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is real(RK)           :: array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!>           overwritten by the solution matrix  X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DTRSM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: DIAG, SIDE, TRANSA, UPLO
  real(RK), intent(in)     :: ALPHA
  integer, intent(in)      :: LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)     :: A(LDA, *)
  real(RK), intent(inout)  :: B(LDB, *)
!     ..
!  =====================================================================
!     .. Intrinsic Functions ..
  intrinsic MAX
!     ..
!     .. Local Scalars ..
  real(RK)                :: TEMP
  integer                 :: I, INFO, J, K, NROWA
  logical                 :: LSIDE, NOUNIT, UPPER
!     ..
!     .. Parameters ..
! real(RK), parameter      :: ONE = 1.0_RK
! real(RK), parameter      :: ZERO = 0.0_RK
!
! interface
!     .. External Subroutines ..
!   !include 'xerbla.h'
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!
!     Test the input parameters.
!
  LSIDE = mobbrmsd_LSAME(SIDE, 'L')
  if (LSIDE) then
    NROWA = M
  else
    NROWA = N
  end if
  NOUNIT = mobbrmsd_LSAME(DIAG, 'N')
  UPPER = mobbrmsd_LSAME(UPLO, 'U')
!
  INFO = 0
  if ((.not. LSIDE) .and. (.not. mobbrmsd_LSAME(SIDE, 'R'))) then
    INFO = 1
  else if ((.not. UPPER) .and. (.not. mobbrmsd_LSAME(UPLO, 'L'))) then
    INFO = 2
  else if ((.not. mobbrmsd_LSAME(TRANSA, 'N')) .and.&
 &         (.not. mobbrmsd_LSAME(TRANSA, 'T')) .and.&
 &         (.not. mobbrmsd_LSAME(TRANSA, 'C'))) then
    INFO = 3
  else if ((.not. mobbrmsd_LSAME(DIAG, 'U')) .and. (.not. mobbrmsd_LSAME(DIAG, 'N'))) then
    INFO = 4
  else if (M < 0) then
    INFO = 5
  else if (N < 0) then
    INFO = 6
  else if (LDA < MAX(1, NROWA)) then
    INFO = 9
  else if (LDB < MAX(1, M)) then
    INFO = 11
  end if
  if (INFO /= 0) then
    !CALL XERBLA('mobbrmsd_DTRSM ',INFO)
    return
  end if
!
!     Quick return if possible.
!
  if (M == 0 .or. N == 0) return
!
!     And when  alpha.eq.zero.
!
  if (ALPHA == ZERO) then
    do J = 1, N
      do I = 1, M
        B(I, J) = ZERO
      end do
    end do
    return
  end if
!
!     Start the operations.
!
  if (LSIDE) then
    if (mobbrmsd_LSAME(TRANSA, 'N')) then
!
!           Form  B := alpha*inv( A )*B.
!
      if (UPPER) then
        do J = 1, N
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, J) = ALPHA * B(I, J)
            end do
          end if
          do K = M, 1, -1
            if (B(K, J) /= ZERO) then
              if (NOUNIT) B(K, J) = B(K, J) / A(K, K)
              do I = 1, K - 1
                B(I, J) = B(I, J) - B(K, J) * A(I, K)
              end do
            end if
          end do
        end do
      else
        do J = 1, N
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, J) = ALPHA * B(I, J)
            end do
          end if
          do K = 1, M
            if (B(K, J) /= ZERO) then
              if (NOUNIT) B(K, J) = B(K, J) / A(K, K)
              do I = K + 1, M
                B(I, J) = B(I, J) - B(K, J) * A(I, K)
              end do
            end if
          end do
        end do
      end if
    else
!
!           Form  B := alpha*inv( A**T )*B.
!
      if (UPPER) then
        do J = 1, N
          do I = 1, M
            TEMP = ALPHA * B(I, J)
            do K = 1, I - 1
              TEMP = TEMP - A(K, I) * B(K, J)
            end do
            if (NOUNIT) TEMP = TEMP / A(I, I)
            B(I, J) = TEMP
          end do
        end do
      else
        do J = 1, N
          do I = M, 1, -1
            TEMP = ALPHA * B(I, J)
            do K = I + 1, M
              TEMP = TEMP - A(K, I) * B(K, J)
            end do
            if (NOUNIT) TEMP = TEMP / A(I, I)
            B(I, J) = TEMP
          end do
        end do
      end if
    end if
  else
    if (mobbrmsd_LSAME(TRANSA, 'N')) then
!
!           Form  B := alpha*B*inv( A ).
!
      if (UPPER) then
        do J = 1, N
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, J) = ALPHA * B(I, J)
            end do
          end if
          do K = 1, J - 1
            if (A(K, J) /= ZERO) then
              do I = 1, M
                B(I, J) = B(I, J) - A(K, J) * B(I, K)
              end do
            end if
          end do
          if (NOUNIT) then
            TEMP = ONE / A(J, J)
            do I = 1, M
              B(I, J) = TEMP * B(I, J)
            end do
          end if
        end do
      else
        do J = N, 1, -1
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, J) = ALPHA * B(I, J)
            end do
          end if
          do K = J + 1, N
            if (A(K, J) /= ZERO) then
              do I = 1, M
                B(I, J) = B(I, J) - A(K, J) * B(I, K)
              end do
            end if
          end do
          if (NOUNIT) then
            TEMP = ONE / A(J, J)
            do I = 1, M
              B(I, J) = TEMP * B(I, J)
            end do
          end if
        end do
      end if
    else
!
!           Form  B := alpha*B*inv( A**T ).
!
      if (UPPER) then
        do K = N, 1, -1
          if (NOUNIT) then
            TEMP = ONE / A(K, K)
            do I = 1, M
              B(I, K) = TEMP * B(I, K)
            end do
          end if
          do J = 1, K - 1
            if (A(J, K) /= ZERO) then
              TEMP = A(J, K)
              do I = 1, M
                B(I, J) = B(I, J) - TEMP * B(I, K)
              end do
            end if
          end do
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, K) = ALPHA * B(I, K)
            end do
          end if
        end do
      else
        do K = 1, N
          if (NOUNIT) then
            TEMP = ONE / A(K, K)
            do I = 1, M
              B(I, K) = TEMP * B(I, K)
            end do
          end if
          do J = K + 1, N
            if (A(J, K) /= ZERO) then
              TEMP = A(J, K)
              do I = 1, M
                B(I, J) = B(I, J) - TEMP * B(I, K)
              end do
            end if
          end do
          if (ALPHA /= ONE) then
            do I = 1, M
              B(I, K) = ALPHA * B(I, K)
            end do
          end if
        end do
      end if
    end if
  end if
!
  return
!
!     End of mobbrmsd_DTRSM
!
end subroutine mobbrmsd_DTRSM

