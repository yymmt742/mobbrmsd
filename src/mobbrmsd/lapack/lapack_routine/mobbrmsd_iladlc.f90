!| mobbrmsd_ILADLC scans a matrix A for its last non-zero column.
!
!  Reference ILADLC is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure function mobbrmsd_ILADLC(M, N, A, LDA)
  implicit none
  integer, intent(in)  :: M
!! The number of rows of the matrix A.
!!
  integer, intent(in)  :: N
!! The number of columns of the matrix A.
!!
  integer, intent(in)  :: LDA
!! The leading dimension of the array A. LDA >= max(1,M).
!!
  real(RK), intent(in) :: A(LDA, *)
!! DOUBLE PRECISION array, dimension (LDA,N)
!! The m by n matrix A.
!!
  integer              :: mobbrmsd_ILADLC
!! The last non-zero column of A.
!!
  integer              :: I
!
! Quick test for the common case where one corner is non-zero.
  if (N == 0) then
    mobbrmsd_ILADLC = N
  else if (A(1, N) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILADLC = N
  else
! Now scan each column from the end, returning with the first non-zero.
    do mobbrmsd_ILADLC = N, 1, -1
      do I = 1, M
        if (A(I, mobbrmsd_ILADLC) /= ZERO) return
      end do
    end do
  end if
  return
end function mobbrmsd_ILADLC

