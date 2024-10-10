!| mobbrmsd_ILASLR scans a matrix for its last non-zero row.
!
!  Reference ILASLR is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! December 2016
!
pure function mobbrmsd_ILASLR(M, N, A, LDA)
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
  real(RK), intent(in) ::  A(LDA, *)
!! A is REAL array, dimension(LDA, N)
!!
!! The m by n matrix A.
!!
  integer :: mobbrmsd_ILASLR
!! Index of the last non-zero row of A.
!!
  integer :: I, J
!
! Quick test for the common case where one corner is non-zero.
!
  if (M == 0) then
    mobbrmsd_ILASLR = M
  elseif (A(M, 1) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILASLR = M
  else
!
! Scan up each column tracking the last zero row seen.
!
    mobbrmsd_ILASLR = 0
    do J = 1, N
      I = M
      do while ((A(MAX(I, 1), J) == ZERO) .and. (I >= 1))
        I = I - 1
      end do
      mobbrmsd_ILASLR = MAX(mobbrmsd_ILASLR, I)
    end do
  end if
  return
end

