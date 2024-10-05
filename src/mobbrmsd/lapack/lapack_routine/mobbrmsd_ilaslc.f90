!| mobbrmsd_ILASLC scans a matrix A for its last non-zero column.
!
!  Reference ILASLC is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! June 2017
!
pure function mobbrmsd_ILASLC(M, N, A, LDA)
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
!! A is DOUBLE PRECISION array, dimension (LDA,N)
!! The m by n matrix A.
!!
  integer              :: mobbrmsd_ILASLC
!! The last non-zero column of A.
!!
  integer :: I
!
! Quick test for the common case where one corner is non-zero.
  if (N == 0) then
    mobbrmsd_ILASLC = N
  else if (A(1, N) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILASLC = N
  else
! Now scan each column from the end, returning with the first non-zero.
    do mobbrmsd_ILASLC = N, 1, -1
      do I = 1, M
        if (A(I, mobbrmsd_ILASLC) /= ZERO) return
      end do
    end do
  end if
  return
end

