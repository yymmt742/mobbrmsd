!| mobbrmsd_DLASR applies a sequence of plane rotations to a general rectangular matrix.
!
!  mobbrmsd_DLASR applies a sequence of plane rotations to
!  a real matrix \( A \), from either the left or the right.
!
!  When SIDE = 'L', the transformation takes the form
!
!  \[
!     A \gets PA
!  \]
!
!  and when SIDE = 'R', the transformation takes the form
!
!  \[
!     A \gets AP ^ \top
!  \]
!
!  where \( P \) is an orthogonal matrix consisting of
!  a sequence of z plane rotations,
!  with \( z = M \) when SIDE = 'L' and \( z = N \)
!  when SIDE = 'R', and \( P ^ \top \) is the transpose of \( P \).
!
!  When DIRECT = 'F' (Forward sequence), then
!
!  \[
!     P = P _ {z-1} \cdots P _ 2 P _ 1
!  \]
!
!  and when DIRECT = 'B' (Backward sequence), then
!
!  \[
!     P = P _ 1 P _ 2 \cdots  P _ {z-1}
!  \]
!
!  where \( P _ k \) is a plane rotation matrix defined by the 2-by-2 rotation
!
!  \[
!     R_k =
!     \left (
!     \begin{array}{}
!       c _ k &  s _ k \\
!      -s _ k &  c _ k \\
!     \end{array}
!     \right )
!  \]
!
!  When PIVOT = 'V' (Variable pivot), the rotation is performed
!  for the plane \( (k,k+1) \), so \( P _ k \) has the form
!
!  \[
!     P_k =
!     \left (
!     \begin{array}{}
!        1   &        &     &        &        &     &        &    \\
!            & \ddots &     &        &        &     &        &    \\
!            &        &  1  &        &        &     &        &    \\
!            &        &     &  c _ k &  s _ k &     &        &    \\
!            &        &     & -s _ k &  c _ k &     &        &    \\
!            &        &     &        &        &  1  &        &    \\
!            &        &     &        &        &     & \ddots &    \\
!            &        &     &        &        &     &        &  1 \\
!     \end{array}
!     \right )
!  \]
!
!  where \( R _ k  \) appears as a rank-2 modification to
!  the identity matrix in rows and columns \( k \) and \( k + 1 \).
!
!  When PIVOT = 'T' (Top pivot), the rotation is performed
!  for the plane \( (1,k+1) \), so \( P _ k \) has the form
!
!  \[
!     P_k =
!     \left (
!     \begin{array}{}
!      c _ k &     &        &     & s _ k &     &        &    \\
!            &  1  &        &     &       &     &        &    \\
!            &     & \ddots &     &       &     &        &    \\
!            &     &        &  1  &       &     &        &    \\
!     -s _ k &     &        &     & c _ k &     &        &    \\
!            &     &        &     &       &  1  &        &    \\
!            &     &        &     &       &     & \ddots &    \\
!            &     &        &     &       &     &        &  1 \\
!     \end{array}
!     \right )
!  \]
!
!  where \( R_k \) appears in rows and columns 1 and k+1.
!
!  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!  performed for the plane \( (k,z) \), giving \( P _ k \) the form
!
!  \[
!     P_k =
!     \left (
!     \begin{array}{}
!       1  &        &     &        &     &        &     &       \\
!          & \ddots &     &        &     &        &     &       \\
!          &        &  1  &        &     &        &     &       \\
!          &        &     &  c _ k &     &        &     & s _ k \\
!          &        &     &        &  1  &        &     &       \\
!          &        &     &        &     & \ddots &     &       \\
!          &        &     &        &     &        &  1  &       \\
!          &        &     & -s _ k &     &        &     & c _ k \\
!     \end{array}
!     \right )
!  \]
!
!  where \( R _ k \) appears in rows and columns \( k \) and \( z \).
!  The rotations are performed without ever forming \( P _ k \) explicitly.
!
!  Reference DLASR is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLASR(SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA)
  implicit none
  character, intent(in) :: SIDE
!!   Specifies whether the plane rotation matrix P is applied to
!!   A on the left or the right.
!!
!!   = 'L':  Left, compute A := P*A
!!
!!   = 'R':  Right, compute A:= A*P**T
!!
  character, intent(in) :: PIVOT
!!   Specifies the plane for which P(k) is a plane rotation
!!   matrix.
!!
!!   = 'V':  Variable pivot, the plane (k,k+1)
!!
!!   = 'T':  Top pivot, the plane (1,k+1)
!!
!!   = 'B':  Bottom pivot, the plane (k,z)
!!
  character, intent(in) :: DIRECT
!!   Specifies whether P is a forward or backward sequence of
!!   plane rotations.
!!
!!   = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!!
!!   = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!!
  integer, intent(in)     :: M
!!  The number of rows of the matrix A.  If m <= 1, an immediate
!!  return is effected.
!!
  integer, intent(in)     :: N
!!    The number of columns of the matrix A.  If n <= 1, an
!!    immediate return is effected.
!!
  real(RK), intent(in)    :: C(*)
!!   DOUBLE PRECISION array, dimension
!!           (M-1) if SIDE = 'L'
!!           (N-1) if SIDE = 'R'
!!
!!   The cosines c(k) of the plane rotations.
!!
  real(RK), intent(in)    :: S(*)
!!   DOUBLE PRECISION array, dimension
!!           (M-1) if SIDE = 'L'
!!           (N-1) if SIDE = 'R'
!!   The sines s(k) of the plane rotations.  The 2-by-2 plane
!!   rotation part of the matrix P(k), R(k), has the form
!!   R(k) = (  c(k)  s(k) )
!!          ( -s(k)  c(k) ).
!!
  integer, intent(in)     :: LDA
!!   The leading dimension of the array A.  LDA >= max(1,M).
!!
  real(RK), intent(inout) :: A(LDA, *)
!!   DOUBLE PRECISION array, dimension (LDA,N)
!!
!!   The M-by-N matrix A.  On exit, A is overwritten by P*A if
!!   SIDE = 'L' or by A*P**T if SIDE = 'R'.
!!
  integer                 :: I, INFO, J
  real(RK)                :: CTEMP, STEMP, TEMP
  intrinsic               :: MAX
! interface
!   include "lsame.h"
! end interface
!
! Test the input parameters
!
  INFO = 0
  if (.not. (mobbrmsd_LSAME(SIDE, 'L') .or. mobbrmsd_LSAME(SIDE, 'R'))) then
    INFO = 1
  else if (.not. (mobbrmsd_LSAME(PIVOT, 'V') .or. mobbrmsd_LSAME(PIVOT, 'T') .or. mobbrmsd_LSAME(PIVOT, 'B'))) then
    INFO = 2
  else if (.not. (mobbrmsd_LSAME(DIRECT, 'F') .or. mobbrmsd_LSAME(DIRECT, 'B'))) then
    INFO = 3
  else if (M < 0) then
    INFO = 4
  else if (N < 0) then
    INFO = 5
  else if (LDA < MAX(1, M)) then
    INFO = 9
  end if
  if (INFO /= 0) then
    !CALL XERBLA( 'DLASR ', INFO )
    return
  end if
!
!     Quick return if possible
!
  if ((M == 0) .or. (N == 0)) return
!
  if (mobbrmsd_LSAME(SIDE, 'L')) then
!
!        Form  P * A
!
    if (mobbrmsd_LSAME(PIVOT, 'V')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 1, M - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J + 1, I)
              A(J + 1, I) = CTEMP * TEMP - STEMP * A(J, I)
              A(J, I) = STEMP * TEMP + CTEMP * A(J, I)
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = M - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J + 1, I)
              A(J + 1, I) = CTEMP * TEMP - STEMP * A(J, I)
              A(J, I) = STEMP * TEMP + CTEMP * A(J, I)
            end do
          end if
        end do
      end if
    else if (mobbrmsd_LSAME(PIVOT, 'T')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 2, M
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = CTEMP * TEMP - STEMP * A(1, I)
              A(1, I) = STEMP * TEMP + CTEMP * A(1, I)
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = M, 2, -1
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = CTEMP * TEMP - STEMP * A(1, I)
              A(1, I) = STEMP * TEMP + CTEMP * A(1, I)
            end do
          end if
        end do
      end if
    else if (mobbrmsd_LSAME(PIVOT, 'B')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 1, M - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = STEMP * A(M, I) + CTEMP * TEMP
              A(M, I) = CTEMP * A(M, I) - STEMP * TEMP
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = M - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = STEMP * A(M, I) + CTEMP * TEMP
              A(M, I) = CTEMP * A(M, I) - STEMP * TEMP
            end do
          end if
        end do
      end if
    end if
  else if (mobbrmsd_LSAME(SIDE, 'R')) then
!
!        Form P**T
!
    if (mobbrmsd_LSAME(PIVOT, 'V')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 1, N - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J + 1)
              A(I, J + 1) = CTEMP * TEMP - STEMP * A(I, J)
              A(I, J) = STEMP * TEMP + CTEMP * A(I, J)
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = N - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J + 1)
              A(I, J + 1) = CTEMP * TEMP - STEMP * A(I, J)
              A(I, J) = STEMP * TEMP + CTEMP * A(I, J)
            end do
          end if
        end do
      end if
    else if (mobbrmsd_LSAME(PIVOT, 'T')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 2, N
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = CTEMP * TEMP - STEMP * A(I, 1)
              A(I, 1) = STEMP * TEMP + CTEMP * A(I, 1)
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = N, 2, -1
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = CTEMP * TEMP - STEMP * A(I, 1)
              A(I, 1) = STEMP * TEMP + CTEMP * A(I, 1)
            end do
          end if
        end do
      end if
    else if (mobbrmsd_LSAME(PIVOT, 'B')) then
      if (mobbrmsd_LSAME(DIRECT, 'F')) then
        do J = 1, N - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = STEMP * A(I, N) + CTEMP * TEMP
              A(I, N) = CTEMP * A(I, N) - STEMP * TEMP
            end do
          end if
        end do
      else if (mobbrmsd_LSAME(DIRECT, 'B')) then
        do J = N - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = STEMP * A(I, N) + CTEMP * TEMP
              A(I, N) = CTEMP * A(I, N) - STEMP * TEMP
            end do
          end if
        end do
      end if
    end if
  end if
!
  return
!
!     End of mobbrmsd_DLASR
!
end subroutine mobbrmsd_DLASR

