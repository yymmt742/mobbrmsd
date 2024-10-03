!> \brief \b mobbrmsd_DGER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       real(RK)           :: ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       real(RK)           :: A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DGER   performs the rank 1 operation
!>
!>    A := alpha*x*y**T + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is real(RK)           ::.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is real(RK)           :: array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is real(RK)           :: array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is real(RK)           :: array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
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
!> \ingroup double_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
! use LA_CONSTANTS, only: RK => dp
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in)    :: ALPHA
  integer, intent(in)     :: INCX, INCY, LDA, M, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(in)    :: X(*), Y(*)
  real(RK), intent(inout) :: A(LDA, *)
!     ..
!
!  =====================================================================
!..
!.. Local Scalars ..
  real(RK)             :: TEMP
  integer              :: I, INFO, IX, J, JY, KX
!..
  intrinsic            :: MAX
!..
!
!.. Parameters ..
! real(RK), parameter   :: ZERO = 0.0_RK
!..
! Test the input parameters.
!
  INFO = 0
  if (M < 0) then
    INFO = 1
  else if (N < 0) then
    INFO = 2
  else if (INCX == 0) then
    INFO = 5
  else if (INCY == 0) then
    INFO = 7
  else if (LDA < MAX(1, M)) then
    INFO = 9
  end if
  if (INFO /= 0) then
!    CALL XERBLA('mobbrmsd_DGER  ',INFO)
    return
  end if
!
!     Quick return if possible.
!
  if ((M == 0) .or. (N == 0) .or. (ALPHA == ZERO)) return
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if (INCY > 0) then
    JY = 1
  else
    JY = 1 - (N - 1) * INCY
  end if
  if (INCX == 1) then
    do J = 1, N
      if (Y(JY) /= ZERO) then
        TEMP = ALPHA * Y(JY)
        do I = 1, M
          A(I, J) = A(I, J) + X(I) * TEMP
        end do
      end if
      JY = JY + INCY
    end do
  else
    if (INCX > 0) then
      KX = 1
    else
      KX = 1 - (M - 1) * INCX
    end if
    do J = 1, N
      if (Y(JY) /= ZERO) then
        TEMP = ALPHA * Y(JY)
        IX = KX
        do I = 1, M
          A(I, J) = A(I, J) + X(IX) * TEMP
          IX = IX + INCX
        end do
      end if
      JY = JY + INCY
    end do
  end if
!
  return
!
!     End of mobbrmsd_DGER
!
end subroutine mobbrmsd_DGER

