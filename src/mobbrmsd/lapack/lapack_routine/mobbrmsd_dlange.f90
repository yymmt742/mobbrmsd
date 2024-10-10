!| mobbrmsd_DLANGE returns the value of the 1-norm, Frobenius norm,
!  infinity-norm, or the largest absolute value of
!  any element of a general rectangular matrix.
!
!  mobbrmsd_DLANGE returns the value of the one norm,
!  or the Frobenius norm, or the  infinity norm,
!  or the  element of  largest absolute value  of a
!  real matrix \( A \).
!
!  \[
!    RES =
!    \left \{
!      \begin{array}{}
!        \max( \lvert A _ {i,j} \rvert ) & \text{if NORM} = \text{'M' or 'm'} \\
!        \text{norm1}( A )               & \text{if NORM} = \text{'1', 'O' or 'o'} \\
!        \text{normI}( A )               & \text{if NORM} = \text{'I' or 'i'} \\
!        \text{normF}( A )               & \text{if NORM} = \text{'F', 'f', 'E' or 'e'} \\
!      \end{array}
!    \right .
!  \]
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLANGE(NORM, M, N, A, LDA, RES, WORK)
  character, intent(in)    :: NORM
!!   Specifies the value to be returned in mobbrmsd_DLANGE as described
!!   above.
!!
  integer, intent(in)      :: M
!!   The number of rows of the matrix A.  M >= 0.  When M = 0,
!!   mobbrmsd_DLANGE is set to zero.
!!
  integer, intent(in)      :: N
!!   The number of columns of the matrix A.  N >= 0.  When N = 0,
!!   mobbrmsd_DLANGE is set to zero.
!!
  integer, intent(in)      :: LDA
!!   The leading dimension of the array A.  LDA >= max(M,1).
!!
  real(RK), intent(inout)  :: A(LDA, *)
!!   DOUBLE PRECISION array, dimension (LDA,N)
!!   The m by n matrix A.
!!
  real(RK), intent(out)    :: RES
!!   Return value
!!
  real(RK), intent(out)    :: WORK(*)
!!   WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!!   where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!!   referenced.
!!
  integer                  :: I, J
  real(RK)                 :: SCL, SSM, TEMP
  intrinsic                :: ABS, MIN, SQRT
! interface
!     .. External Subroutines ..
!   include 'dlassq.h'
!     .. External Functions ..
!   include 'lsame.h'
!   include 'disnan.h'
! end interface
!
  if (MIN(M, N) == 0) then
    RES = ZERO
  else if (mobbrmsd_LSAME(NORM, 'M')) then
!
!   Find max(abs(A(i,j))).
!
    RES = ZERO
    do J = 1, N
      do I = 1, M
        TEMP = ABS(A(I, J))
        if (RES < TEMP .or. IEEE_IS_NAN(TEMP)) RES = TEMP
      end do
    end do
  else if ((mobbrmsd_LSAME(NORM, 'O')) .or. (NORM == '1')) then
!
!   Find norm1(A).
!
    RES = ZERO
    do J = 1, N
      SSM = ZERO
      do I = 1, M
        SSM = SSM + ABS(A(I, J))
      end do
      if (RES < SSM .or. IEEE_IS_NAN(SSM)) RES = SSM
    end do
  else if (mobbrmsd_LSAME(NORM, 'I')) then
!
!   Find normI(A).
!
    do I = 1, M
      WORK(I) = ZERO
    end do
    do J = 1, N
      do I = 1, M
        WORK(I) = WORK(I) + ABS(A(I, J))
      end do
    end do
    RES = ZERO
    do I = 1, M
      TEMP = WORK(I)
      if (RES < TEMP .or. IEEE_IS_NAN(TEMP)) RES = TEMP
    end do
  else if ((mobbrmsd_LSAME(NORM, 'F')) .or. (mobbrmsd_LSAME(NORM, 'E'))) then
!
!   Find normF(A).
!
    SCL = ZERO
    SSM = ONE
    do J = 1, N
      call mobbrmsd_DLASSQ(M, A(1, J), 1, SCL, SSM)
    end do
    RES = SCL * SQRT(SSM)
  end if
!
  return
!
! End of mobbrmsd_DLANGE
!
end subroutine mobbrmsd_DLANGE

