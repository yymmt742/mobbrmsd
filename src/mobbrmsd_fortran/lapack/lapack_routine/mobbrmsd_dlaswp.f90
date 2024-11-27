!| mobbrmsd_DLASWP performs a series of row interchanges on a general rectangular matrix.
!
!  mobbrmsd_DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Reference DLASWP is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK driver routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!    R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
pure subroutine mobbrmsd_DLASWP(N, A, LDA, K1, K2, IPIV, INCX)
  implicit none
  integer, intent(in)     :: N
!!          The number of columns of the matrix A.
!!
  integer, intent(in)     :: LDA
!!          The leading dimension of the array A.
!!
  real(RK), intent(inout) :: A(LDA, *)
!!          A is DOUBLE PRECISION array, dimension (LDA,N)
!!
!!          On entry, the matrix of column dimension N to which the row
!!          interchanges will be applied.
!!
!!          On exit, the permuted matrix.
!!
  integer, intent(in)     :: K1
!!          The first element of IPIV for which a row interchange will
!!          be done.
!!
  integer, intent(in)     :: K2
!!          (K2-K1+1) is the number of elements of IPIV for which a row
!!          interchange will be done.
!!
  integer, intent(in)     :: IPIV(*)
!!          IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX))
!!
!!          The vector of pivot indices. Only the elements in positions
!!          K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed.
!!          IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be
!!          interchanged.
!!
  integer, intent(in)     :: INCX
!!          The increment between successive values of IPIV. If INCX
!!          is negative, the pivots are applied in reverse order.
!!
  integer  :: I, I1, I2, INC, IP, IX, IX0, J, K, N32
  real(RK) :: TEMP
!
!     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
!     K1 through K2.
!
  if (INCX > 0) then
    IX0 = K1
    I1 = K1
    I2 = K2
    INC = 1
  else if (INCX < 0) then
    IX0 = K1 + (K1 - K2) * INCX
    I1 = K2
    I2 = K1
    INC = -1
  else
    return
  end if
!
  N32 = (N / 32) * 32
  if (N32 /= 0) then
    do J = 1, N32, 32
      IX = IX0
      do I = I1, I2, INC
        IP = IPIV(IX)
        if (IP /= I) then
          do K = J, J + 31
            TEMP = A(I, K)
            A(I, K) = A(IP, K)
            A(IP, K) = TEMP
          end do
        end if
        IX = IX + INCX
      end do
    end do
  end if
  if (N32 /= N) then
    N32 = N32 + 1
    IX = IX0
    do I = I1, I2, INC
      IP = IPIV(IX)
      if (IP /= I) then
        do K = N32, N
          TEMP = A(I, K)
          A(I, K) = A(IP, K)
          A(IP, K) = TEMP
        end do
      end if
      IX = IX + INCX
    end do
  end if
!
  return
!
!     End of mobbrmsd_DLASWP
!
end subroutine mobbrmsd_DLASWP

