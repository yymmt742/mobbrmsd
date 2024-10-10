!| mobbrmsd_DGER performs the rank 1 operation
!
!  \[ A \gets \alpha x y ^ {\top} + A \],
!
!  where \( \alpha \) is a scalar,  \( x \) is an  \( m \) element vector,
!  \( y \) is an \( n \) element vector,
!  and \( A \) is an \( m \) by \( n \) matrix.
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
pure subroutine mobbrmsd_DGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
  implicit none
  integer, intent(in)     :: M
!!           On entry, M specifies the number of rows of the matrix A.
!!           M must be at least zero.
!!
  integer, intent(in)     :: N
!!           On entry, N specifies the number of columns of the matrix A.
!!           N must be at least zero.
!!
  real(RK), intent(in)    :: ALPHA
!!           On entry, ALPHA specifies the scalar alpha.
!!
  real(RK), intent(in)    :: X(*)
!!          DOUBLE PRECISION array, dimension at least
!!           ( 1 + ( m - 1 )*abs( INCX ) ).
!!           Before entry, the incremented array X must contain the m
!!           element vector x.
!!
  integer, intent(in)     :: INCX
!!           On entry, INCX specifies the increment for the elements of
!!           X. INCX must not be zero.
!!
  real(RK), intent(in)    :: Y(*)
!!          DOUBLE PRECISION array, dimension at least
!!           ( 1 + ( n - 1 )*abs( INCY ) ).
!!           Before entry, the incremented array Y must contain the n
!!           element vector y.
!!
  integer, intent(in)     :: INCY
!!           On entry, INCY specifies the increment for the elements of
!!           Y. INCY must not be zero.
!!
  integer, intent(in)     :: LDA
!!           On entry, LDA specifies the first dimension of A as declared
!!           in the calling (sub) program. LDA must be at least
!!           max( 1, m ).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!           DOUBLE PRECISION array, dimension ( LDA, N )
!!           Before entry, the leading m by n part of the array A must
!!           contain the matrix of coefficients. On exit, A is
!!           overwritten by the updated matrix.
!!
  real(RK)             :: TEMP
  integer              :: I, INFO, IX, J, JY, KX
  intrinsic            :: MAX
!
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
!    CALL XERBLA('DGER  ',INFO)
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

