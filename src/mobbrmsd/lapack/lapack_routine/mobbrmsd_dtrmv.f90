!| mobbrmsd_DTRMV performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A**T*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  reference DTRMV is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
pure subroutine mobbrmsd_DTRMV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
  implicit none
  character, intent(in) :: UPLO
!!           On entry, UPLO specifies whether the matrix is an upper or
!!           lower triangular matrix as follows:
!!
!!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!!
!!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!!
  character, intent(in) :: TRANS
!!           On entry, TRANS specifies the operation to be performed as
!!           follows:
!!
!!              TRANS = 'N' or 'n'   x := A*x.
!!
!!              TRANS = 'T' or 't'   x := A**T*x.
!!
!!              TRANS = 'C' or 'c'   x := A**T*x.
!!
  character, intent(in) :: DIAG
!!           On entry, DIAG specifies whether or not A is unit
!!           triangular as follows:
!!
!!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!!
!!              DIAG = 'N' or 'n'   A is not assumed to be unit
!!                                  triangular.
!!
  integer, intent(in)   :: N
!!           On entry, N specifies the order of the matrix A.
!!           N must be at least zero.
  integer, intent(in)   :: LDA
!!           On entry, LDA specifies the first dimension of A as declared
!!           in the calling (sub) program. LDA must be at least
!!           max( 1, n ).
  real(RK), intent(in)    :: A(LDA, *)
!!          A is REAL array, dimension ( LDA, N )
!!
!!           Before entry with  UPLO = "U" or "u", the leading n by n
!!           upper triangular part of the array A must contain the upper
!!           triangular matrix and the strictly lower triangular part of
!!           A is not referenced.
!!           Before entry with UPLO = "L" or "l", the leading n by n
!!           lower triangular part of the array A must contain the lower
!!           triangular matrix and the strictly upper triangular part of
!!           A is not referenced.
!!           Note that when  DIAG = "U" or "u", the diagonal elements of
!!           A are not referenced either, but are assumed to be unity.
!!
  integer, intent(in)   :: INCX
!!           On entry, INCX specifies the increment for the elements of
!!           X. INCX must not be zero.
  real(RK), intent(inout) :: X(*)
!!           X is REAL array, dimension at least
!!           ( 1 + ( n - 1 )*abs( INCX ) ).
!!           Before entry, the incremented array X must contain the n
!!           element vector x. On exit, X is overwritten with the
!!           transformed vector x.
!!
  real(RK)  :: TEMP
  integer   :: I, INFO, IX, J, JX, KX
  logical   :: NOUNIT
  intrinsic :: MAX
! interface
!   include 'lsame.h'
! end interface
!
! Test the input parameters.
!
  INFO = 0
  if (.not. mobbrmsd_LSAME(UPLO, 'U') .and. .not. mobbrmsd_LSAME(UPLO, 'L')) then
    INFO = 1
  else if (.not. mobbrmsd_LSAME(TRANS, 'N') .and. .not. mobbrmsd_LSAME(TRANS, 'T') .and.&
 &         .not. mobbrmsd_LSAME(TRANS, 'C')) then
    INFO = 2
  else if (.not. mobbrmsd_LSAME(DIAG, 'U') .and. .not. mobbrmsd_LSAME(DIAG, 'N')) then
    INFO = 3
  else if (N < 0) then
    INFO = 4
  else if (LDA < MAX(1, N)) then
    INFO = 6
  else if (INCX == 0) then
    INFO = 8
  end if
  if (INFO /= 0) then
    return
  end if
!
!     Quick return if possible.
!
  if (N == 0) return
!
  NOUNIT = mobbrmsd_LSAME(DIAG, 'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  if (INCX <= 0) then
    KX = 1 - (N - 1) * INCX
  else if (INCX /= 1) then
    KX = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if (mobbrmsd_LSAME(TRANS, 'N')) then
!
!        Form  x := A*x.
!
    if (mobbrmsd_LSAME(UPLO, 'U')) then
      if (INCX == 1) then
        do J = 1, N
          if (X(J) /= ZERO) then
            TEMP = X(J)
            do I = 1, J - 1
              X(I) = X(I) + TEMP * A(I, J)
            end do
            if (NOUNIT) X(J) = X(J) * A(J, J)
          end if
        end do
      else
        JX = KX
        do J = 1, N
          if (X(JX) /= ZERO) then
            TEMP = X(JX)
            IX = KX
            do I = 1, J - 1
              X(IX) = X(IX) + TEMP * A(I, J)
              IX = IX + INCX
            end do
            if (NOUNIT) X(JX) = X(JX) * A(J, J)
          end if
          JX = JX + INCX
        end do
      end if
    else
      if (INCX == 1) then
        do J = N, 1, -1
          if (X(J) /= ZERO) then
            TEMP = X(J)
            do I = N, J + 1, -1
              X(I) = X(I) + TEMP * A(I, J)
            end do
            if (NOUNIT) X(J) = X(J) * A(J, J)
          end if
        end do
      else
        KX = KX + (N - 1) * INCX
        JX = KX
        do J = N, 1, -1
          if (X(JX) /= ZERO) then
            TEMP = X(JX)
            IX = KX
            do I = N, J + 1, -1
              X(IX) = X(IX) + TEMP * A(I, J)
              IX = IX - INCX
            end do
            if (NOUNIT) X(JX) = X(JX) * A(J, J)
          end if
          JX = JX - INCX
        end do
      end if
    end if
  else
!
!        Form  x := A**T*x.
!
    if (mobbrmsd_LSAME(UPLO, 'U')) then
      if (INCX == 1) then
        do J = N, 1, -1
          TEMP = X(J)
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = J - 1, 1, -1
            TEMP = TEMP + A(I, J) * X(I)
          end do
          X(J) = TEMP
        end do
      else
        JX = KX + (N - 1) * INCX
        do J = N, 1, -1
          TEMP = X(JX)
          IX = JX
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = J - 1, 1, -1
            IX = IX - INCX
            TEMP = TEMP + A(I, J) * X(IX)
          end do
          X(JX) = TEMP
          JX = JX - INCX
        end do
      end if
    else
      if (INCX == 1) then
        do J = 1, N
          TEMP = X(J)
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = J + 1, N
            TEMP = TEMP + A(I, J) * X(I)
          end do
          X(J) = TEMP
        end do
      else
        JX = KX
        do J = 1, N
          TEMP = X(JX)
          IX = JX
          if (NOUNIT) TEMP = TEMP * A(J, J)
          do I = J + 1, N
            IX = IX + INCX
            TEMP = TEMP + A(I, J) * X(IX)
          end do
          X(JX) = TEMP
          JX = JX + INCX
        end do
      end if
    end if
  end if
!
  return
!
!     End of mobbrmsd_DTRMV
!
end subroutine mobbrmsd_DTRMV

