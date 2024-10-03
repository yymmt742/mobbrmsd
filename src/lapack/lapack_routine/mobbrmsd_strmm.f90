!> \brief \b mobbrmsd_STRMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_STRMM  performs one of the matrix-matrix operations
!>
!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!>
!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
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
!>          ALPHA is REAL
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
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
!>          B is REAL array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B,  and  on exit  is overwritten  by the
!>           transformed matrix.
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
!> \date December 2016
!
!> \ingroup single_blas_level3
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
pure subroutine mobbrmsd_STRMM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
  implicit none
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  integer, intent(in)   ::  LDA, LDB, M, N
  character, intent(in) ::  DIAG, SIDE, TRANSA, UPLO
  real, intent(in)      :: ALPHA
!..
!..Array Arguments..
  real(RK), intent(in)    :: A(LDA, *)
  real(RK), intent(inout) :: B(LDB, *)
!..
!
!  =====================================================================
!
!..intrinsic Functions..
  intrinsic :: MAX
!..
!..Local Scalars..
  real(RK) :: TEMP
  integer :: I, INFO, J, K, NROWA
  logical :: LSIDE, NOUNIT, UPPER
!..
!..Parameters..
! real, parameter :: ONE = 1.0E+0, ZERO = 0.0E+0
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
! end interface
!..
!
!Test the input parameters.
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
  &(.not. mobbrmsd_LSAME(TRANSA, 'T')) .and.&
  &(.not. mobbrmsd_LSAME(TRANSA, 'C'))) then
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
!   call XERBLA('STRMM ', INFO)
    return
  end if
!
! Quick return if possible.
!
  if (M == 0 .or. N == 0) return
!
! And when alpha == zero.
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
! Start the operations.
!
  if (LSIDE) then
    if (mobbrmsd_LSAME(TRANSA, 'N')) then
!
! Form B: = alpha * A * B.
!
      if (UPPER) then
        do J = 1, N
          do K = 1, M
            if (B(K, J) /= ZERO) then
              TEMP = ALPHA * B(K, J)
              do I = 1, K - 1
                B(I, J) = B(I, J) + TEMP * A(I, K)
              end do
              if (NOUNIT) TEMP = TEMP * A(K, K)
              B(K, J) = TEMP
            end if
          end do
        end do
      else
        do J = 1, N
          do K = M, 1, -1
            if (B(K, J) /= ZERO) then
              TEMP = ALPHA * B(K, J)
              B(K, J) = TEMP
              if (NOUNIT) B(K, J) = B(K, J) * A(K, K)
              do I = K + 1, M
                B(I, J) = B(I, J) + TEMP * A(I, K)
              end do
            end if
          end do
        end do
      end if
    else
!
! Form B: = alpha * A**T * B.
!
      if (UPPER) then
        do J = 1, N
          do I = M, 1, -1
            TEMP = B(I, J)
            if (NOUNIT) TEMP = TEMP * A(I, I)
            do K = 1, I - 1
              TEMP = TEMP + A(K, I) * B(K, J)
            end do
            B(I, J) = ALPHA * TEMP
          end do
        end do
      else
        do J = 1, N
          do I = 1, M
            TEMP = B(I, J)
            if (NOUNIT) TEMP = TEMP * A(I, I)
            do K = I + 1, M
              TEMP = TEMP + A(K, I) * B(K, J)
            end do
            B(I, J) = ALPHA * TEMP
          end do
        end do
      end if
    end if
  else
    if (mobbrmsd_LSAME(TRANSA, 'N')) then
!
! Form B: = alpha * B * A.
!
      if (UPPER) then
        do J = N, 1, -1
          TEMP = ALPHA
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = 1, M
            B(I, J) = TEMP * B(I, J)
          end do
          do K = 1, J - 1
            if (A(K, J) /= ZERO) then
              TEMP = ALPHA * A(K, J)
              do I = 1, M
                B(I, J) = B(I, J) + TEMP * B(I, K)
              end do
            end if
          end do
        end do
      else
        do J = 1, N
          TEMP = ALPHA
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = 1, M
            B(I, J) = TEMP * B(I, J)
          end do
          do K = J + 1, N
            if (A(K, J) /= ZERO) then
              TEMP = ALPHA * A(K, J)
              do I = 1, M
                B(I, J) = B(I, J) + TEMP * B(I, K)
              end do
            end if
          end do
        end do
      end if
    else
!
! Form B: = alpha * B * A**T.
!
      if (UPPER) then
        do K = 1, N
          do J = 1, K - 1
            if (A(J, K) /= ZERO) then
              TEMP = ALPHA * A(J, K)
              do I = 1, M
                B(I, J) = B(I, J) + TEMP * B(I, K)
              end do
            end if
          end do
          TEMP = ALPHA
          if (NOUNIT) TEMP = TEMP * A(K, K)
          if (TEMP /= ONE) then
            do I = 1, M
              B(I, K) = TEMP * B(I, K)
            end do
          end if
        end do
      else
        do K = N, 1, -1
          do J = K + 1, N
            if (A(J, K) /= ZERO) then
              TEMP = ALPHA * A(J, K)
              do I = 1, M
                B(I, J) = B(I, J) + TEMP * B(I, K)
              end do
            end if
          end do
          TEMP = ALPHA
          if (NOUNIT) TEMP = TEMP * A(K, K)
          if (TEMP /= ONE) then
            do I = 1, M
              B(I, K) = TEMP * B(I, K)
            end do
          end if
        end do
      end if
    end if
  end if
!
  return
!
! end of mobbrmsd_STRMM.
!
end
