!| mobbrmsd_SLACPY copies all or part of one two-dimensional array to another.
!
!  mobbrmsd_SLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Reference DBDSQR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SLACPY(UPLO, M, N, A, LDA, B, LDB)
  character, intent(in) :: UPLO
!!  Specifies the part of the matrix A to be copied to B.
!!
!!  = 'U':      Upper triangular part
!!
!!  = 'L':      Lower triangular part
!!
!!  Otherwise:  All of the matrix A
!!
  integer, intent(in)   :: M
!!  The number of rows of the matrix A.  M >= 0.
!!
  integer, intent(in)   :: N
!!  The number of columns of the matrix A.  N >= 0.
!!
  integer, intent(in)   :: LDA
!!  The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(in)  :: A(LDA, *)
!!  REAL array, dimension (LDA,N)
!!
!!  The m by n matrix A.  If UPLO = 'U', only the upper triangle
!!  or trapezoid is accessed; if UPLO = 'L', only the lower
!!  triangle or trapezoid is accessed.
!!
  integer, intent(in)   :: LDB
!!  The leading dimension of the array B.  LDB >= max(1,M).
!!
  real(RK), intent(out) :: B(LDB, *)
!!  REAL array, dimension (LDB,N)
!!
!!  On exit, B = A in the locations specified by UPLO.
!!
  integer :: I, J
  intrinsic :: MIN
! interface
! .. External Functions ..
!   include 'lsame.h'
! end interface
!
  if (mobbrmsd_LSAME(UPLO, 'U')) then
    do J = 1, N
      do I = 1, MIN(J, M)
        B(I, J) = A(I, J)
      end do
    end do
  else if (mobbrmsd_LSAME(UPLO, 'L')) then
    do J = 1, N
      do I = J, M
        B(I, J) = A(I, J)
      end do
    end do
  else
    do J = 1, N
      do I = 1, M
        B(I, J) = A(I, J)
      end do
    end do
  end if
  return
!
! end of mobbrmsd_SLACPY
!
end

