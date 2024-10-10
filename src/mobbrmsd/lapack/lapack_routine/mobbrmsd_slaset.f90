!| mobbrmsd_SLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Reference SLASET is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK driver routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SLASET(UPLO, M, N, ALPHA, BETA, A, LDA)
  character, intent(in)  :: UPLO
!!          Specifies the part of the matrix A to be set.
!!
!!          = 'U':      Upper triangular part is set; the strictly lower
!!                      triangular part of A is not changed.
!!
!!          = 'L':      Lower triangular part is set; the strictly upper
!!                      triangular part of A is not changed.
!!
!!          Otherwise:  All of the matrix A is set.
!!
  integer, intent(in)    :: M
!!          The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)    :: N
!!          The number of columns of the matrix A.  N >= 0.
!!
  real(RK), intent(in)   :: ALPHA
!!          The constant to which the offdiagonal elements are to be set.
!!
  real(RK), intent(in)   :: BETA
!!          The constant to which the diagonal elements are to be set.
!!
  integer, intent(in)    :: LDA
!!          The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!          A is REAL array, dimension (LDA,N)
!!          On exit, the leading m-by-n submatrix of A is set as follows:
!!
!!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!!
!!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!!
!!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!!
  integer :: I, J
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: MIN
!..
!..Executable Statements..
!
  if (mobbrmsd_LSAME(UPLO, 'U')) then
!
! Set the strictly upper triangular or trapezoidal part of the
! array to ALPHA.
!
    do J = 2, N
      do I = 1, MIN(J - 1, M)
        A(I, J) = ALPHA
      end do
    end do
!
  else if (mobbrmsd_LSAME(UPLO, 'L')) then
!
! Set the strictly lower triangular or trapezoidal part of the
! array to ALPHA.
!
    do J = 1, MIN(M, N)
      do I = J + 1, M
        A(I, J) = ALPHA
      end do
    end do
!
  else
!
! Set the leading m - by - n submatrix to ALPHA.
!
    do J = 1, N
      do I = 1, M
        A(I, J) = ALPHA
      end do
    end do
  end if
!
! Set the first MIN(M, N) diagonal elements to BETA.
!
  do I = 1, MIN(M, N)
    A(I, I) = BETA
  end do
!
  return
!
!end of mobbrmsd_SLASET
!
end

