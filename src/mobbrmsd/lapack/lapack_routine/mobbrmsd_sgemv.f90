!|  mobbrmsd_SGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!  The vector and matrix arguments are not referenced when N = 0, or M = 0
!
!  Reference SGEMV is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!  -- Written on 22-October-1986.
!      Jack Dongarra, Argonne National Lab.
!      Jeremy Du Croz, Nag Central Office.
!      Sven Hammarling, Nag Central Office.
!      Richard Hanson, Sandia National Labs.
!
pure subroutine mobbrmsd_SGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  implicit none
  character(*), intent(in) :: TRANS
!!           On entry, TRANS specifies the operation to be performed as
!!           follows:
!!
!!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!!
!!              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!!
!!              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!!
  integer, intent(in)      :: M
!!           On entry, M specifies the number of rows of the matrix A.
!!           M must be at least zero.
!!
  integer, intent(in)      :: N
!!           On entry, N specifies the number of columns of the matrix A.
!!           N must be at least zero.
!!
  real(RK), intent(in)     :: ALPHA
!!           On entry, ALPHA specifies the scalar alpha.
!!
  integer, intent(in)      :: LDA
!!           On entry, LDA specifies the first dimension of A as declared
!!           in the calling (sub) program. LDA must be at least
!!           max( 1, m ).
!!
  real(RK), intent(in)     :: A(LDA, *)
!!          A is real(RK)           :: array, dimension ( LDA, N )
!!
!!           Before entry, the leading m by n part of the array A must
!!           contain the matrix of coefficients.
!!
  integer, intent(in)      :: INCX
!!          INCX is INTEGER
!!           On entry, INCX specifies the increment for the elements of
!!           X. INCX must not be zero.
  real(RK), intent(in)     :: X(*)
!!          X is real(RK)           :: array, dimension at least
!!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!!           and at least
!!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!!           Before entry, the incremented array X must contain the
!!           vector x.
!!
  real(RK), intent(in)     :: BETA
!!           On entry, BETA specifies the scalar beta. When BETA is
!!           supplied as zero then Y need not be set on input.
!!
  integer, intent(in)      :: INCY
!!           On entry, INCY specifies the increment for the elements of
!!           Y. INCY must not be zero.
!!
  real(RK), intent(inout)  :: Y(*)
!!          Y is real(RK)           :: array, dimension at least
!!
!!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!!           and at least
!!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!!
!!           Before entry with BETA non-zero, the incremented array Y
!!           must contain the vector y. On exit, Y is overwritten by the
!!           updated vector y.
!!
  intrinsic :: MAX
  real(RK)  :: TEMP
  integer   :: I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
! real(RK), parameter :: ONE = 1.0E+0, ZERO = 0.0E+0
! interface
! .. External Functions ..
!   include 'lsame.h'
! end interface
!
! Test the input parameters.
!
  INFO = 0
  if (.not. mobbrmsd_LSAME(TRANS, 'N') .and. .not. mobbrmsd_LSAME(TRANS, 'T') .and. .not. mobbrmsd_LSAME(TRANS, 'C')) then
    INFO = 1
  else if (M < 0) then
    INFO = 2
  else if (N < 0) then
    INFO = 3
  else if (LDA < MAX(1, M)) then
    INFO = 6
  else if (INCX == 0) then
    INFO = 8
  else if (INCY == 0) then
    INFO = 11
  end if
  if (INFO /= 0) then
!   call XERBLA('SGEMV ', INFO)
    return
  end if
!
! Quick return if possible.
!
  if ((M == 0) .or. (N == 0) .or. ((ALPHA == ZERO) .and. (BETA == ONE))) return
!
! Set LENX and LENY, the lengths of the vectors x and y, and set
! up the start points in X and Y.
!
  if (mobbrmsd_LSAME(TRANS, 'N')) then
    LENX = N
    LENY = M
  else
    LENX = M
    LENY = N
  end if
  if (INCX > 0) then
    KX = 1
  else
    KX = 1 - (LENX - 1) * INCX
  end if
  if (INCY > 0) then
    KY = 1
  else
    KY = 1 - (LENY - 1) * INCY
  end if
!
! Start the operations.In this version the elements of A are
! accessed sequentially with one pass through A.
!
! First form y: = beta * y.
!
  if (BETA /= ONE) then
    if (INCY == 1) then
      if (BETA == ZERO) then
        do I = 1, LENY
          Y(I) = ZERO
        end do
      else
        do I = 1, LENY
          Y(I) = BETA * Y(I)
        end do
      end if
    else
      IY = KY
      if (BETA == ZERO) then
        do I = 1, LENY
          Y(IY) = ZERO
          IY = IY + INCY
        end do
      else
        do I = 1, LENY
          Y(IY) = BETA * Y(IY)
          IY = IY + INCY
        end do
      end if
    end if
  end if
  if (ALPHA == ZERO) return
  if (mobbrmsd_LSAME(TRANS, 'N')) then
!
! Form y: = alpha * A * x + y.
!
    JX = KX
    if (INCY == 1) then
      do J = 1, N
        TEMP = ALPHA * X(JX)
        do I = 1, M
          Y(I) = Y(I) + TEMP * A(I, J)
        end do
        JX = JX + INCX
      end do
    else
      do J = 1, N
        TEMP = ALPHA * X(JX)
        IY = KY
        do I = 1, M
          Y(IY) = Y(IY) + TEMP * A(I, J)
          IY = IY + INCY
        end do
        JX = JX + INCX
      end do
    end if
  else
!
! Form y: = alpha * A**T * x + y.
!
    JY = KY
    if (INCX == 1) then
      do J = 1, N
        TEMP = ZERO
        do I = 1, M
          TEMP = TEMP + A(I, J) * X(I)
        end do
        Y(JY) = Y(JY) + ALPHA * TEMP
        JY = JY + INCY
      end do
    else
      do J = 1, N
        TEMP = ZERO
        IX = KX
        do I = 1, M
          TEMP = TEMP + A(I, J) * X(IX)
          IX = IX + INCX
        end do
        Y(JY) = Y(JY) + ALPHA * TEMP
        JY = JY + INCY
      end do
    end if
  end if
!
  return
!
! end of mobbrmsd_SGEMV.
!
end

