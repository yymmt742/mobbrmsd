!| mobbrmsd_SLARFT forms the triangular factor \( T \)
!  of a block reflector \( H = I - v ^ \top v H \).
!
!  mobbrmsd_SLARFT forms the triangular factor \( T \)
!  of a real block reflector \( H \) of order \( n \),
!  which is defined as a product of \( k \) elementary reflectors.
!
!  If DIRECT = 'F', \( H = H _ 1 H _ 2 \cdots  H _ k \) and
!  \( T \) is upper triangular;
!
!  If DIRECT = 'B', \( H = H _ k \cdots H _ 2 H _ 1 \) and
!  \( T \) is lower triangular;
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  \( H _ i \) is stored in the \( i \)-th column of the array \( V \), and
!
!  \[
!     H  =  I - V T V ^ {\top}
!  \]
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the \( i \)-th row of the array \( V \), and
!
!  \[
!     H  =  I - V ^ {\top} T V
!  \]
!
!  The shape of the matrix \( V \) and the storage of the vectors which
!  define the \( H _ i \) is best illustrated by the following example
!  with \( n = 5 \) and \( k = 3 \).
!  The elements equal to 1 are not stored.
!
!   DIRECT = 'F' and STOREV = 'C':
!
!  \[
!     V  =
!     \left (
!       \begin{array}{}
!          1     &       &       \\
!          v _ 1 & 1     &       \\
!          v _ 1 & v _ 2 & 1     \\
!          v _ 1 & v _ 2 & v _ 3 \\
!          v _ 1 & v _ 2 & v _ 3 \\
!       \end{array}
!     \right)
!  \]
!
!   DIRECT = 'F' and STOREV = 'R':
!
!  \[
!     V  =
!     \left (
!       \begin{array}{}
!          1     & v _ 1 & v _ 1 & v _ 1 & v _ 1 \\
!                & 1     & v _ 2 & v _ 2 & v _ 2 \\
!                &       & 1     & v _ 3 & v _ 3 \\
!       \end{array}
!     \right)
!  \]
!
!   DIRECT = 'B' and STOREV = 'C':
!
!  \[
!     V  =
!     \left (
!       \begin{array}{}
!          v _ 1 & v _ 2 & v _ 3 \\
!          v _ 1 & v _ 2 & v _ 3 \\
!          1     & v _ 2 & v _ 3 \\
!                & 1     & v _ 3 \\
!                &       & 1     \\
!       \end{array}
!     \right)
!  \]
!
!   DIRECT = 'B' and STOREV = 'R':
!
!  \[
!     V  =
!     \left (
!       \begin{array}{}
!          v _ 1 & v _ 1 & 1     &       &       \\
!          v _ 2 & v _ 2 & v _ 2 & 1     &       \\
!          v _ 3 & v _ 3 & v _ 3 & v _ 3 & 1     \\
!       \end{array}
!     \right)
!  \]
!
!  Reference SLARFT is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!
pure subroutine mobbrmsd_SLARFT(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
  implicit none
  character, intent(in)    :: DIRECT
!!  Specifies the order in which the elementary reflectors are
!!  multiplied to form the block reflector:
!!
!!  = 'F': H = H(1) H(2) . . . H(k) (Forward)
!!
!!  = 'B': H = H(k) . . . H(2) H(1) (Backward)
!!
  character, intent(in)    :: STOREV
!!  Specifies how the vectors which define the elementary
!!  reflectors are stored:
!!
!!  = 'C': columnwise
!!
!!  = 'R': rowwise
!!
  integer, intent(in)      :: N
!!  The order of the block reflector H. N >= 0.
!!
  integer, intent(in)      :: K
!!  The order of the triangular factor T (= the number of
!!  elementary reflectors). K >= 1.
!!
  integer, intent(in)      :: LDV
!!  The leading dimension of the array V.
!!
!!  If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!!
  real(RK), intent(in)     :: V(LDV, *)
!!  DOUBLE PRECISION array, dimension
!!
!!  (LDV,K) if STOREV = 'C'
!!
!!  (LDV,N) if STOREV = 'R'
!!
!!  The matrix V.
!!
  real(RK), intent(in)     :: TAU(*)
!!  DOUBLE PRECISION array, dimension (K)
!!
!!  TAU(i) must contain the scalar factor of the elementary
!!  reflector H(i).
!!
  integer, intent(in)      :: LDT
!!  The leading dimension of the array T. LDT >= K.
!!
  real(RK), intent(out)    :: T(LDT, *)
!!  DOUBLE PRECISION array, dimension (LDT,K)
!!
!!  The k by k triangular factor T of the block reflector.
!!  If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!!  lower triangular. The rest of the array is not used.
!!
  integer :: I, J, PREVLASTV, LASTV
! interface
!   include 'lsame.h'
!   include 'sgemv.h'
!   include 'strmv.h'
! end interface
!
! Quick return if possible
!
  if (N == 0) return
!
  if (mobbrmsd_LSAME(DIRECT, 'F')) then
    PREVLASTV = N
    do I = 1, K
      PREVLASTV = MAX(I, PREVLASTV)
      if (TAU(I) == ZERO) then
!
! H(i) = I
!
        do J = 1, I
          T(J, I) = ZERO
        end do
      else
!
! general case
!
        if (mobbrmsd_LSAME(STOREV, 'C')) then
! Skip any trailing zeros.
          do LASTV = N, I + 1, -1
            if (V(LASTV, I) /= ZERO) exit
          end do
          do J = 1, I - 1
            T(J, I) = -TAU(I) * V(I, J)
          end do
          J = MIN(LASTV, PREVLASTV)
!
! T(1:i - 1, i): = -tau(i) * V(i:j, 1:i - 1)**T * V(i:j, i)
!
          call mobbrmsd_SGEMV('Transpose', J - I, I - 1, -TAU(I), &
          & V(I + 1, 1), LDV, V(I + 1, I), 1, ONE, &
          & T(1, I), 1)
        else
! Skip any trailing zeros.
          do LASTV = N, I + 1, -1
            if (V(I, LASTV) /= ZERO) exit
          end do
          do J = 1, I - 1
            T(J, I) = -TAU(I) * V(J, I)
          end do
          J = MIN(LASTV, PREVLASTV)
!
! T(1:i - 1, i): = -tau(i) * V(1:i - 1, i:j) * V(i, i:j)**T
!
          call mobbrmsd_SGEMV('No transpose', I - 1, J - I, -TAU(I), &
              &      V(1, I + 1), LDV, V(I, I + 1), LDV, ONE, T(1, I), 1)
        end if
!
! T(1:i - 1, i): = T(1:i - 1, 1:i - 1) * T(1:i - 1, i)
!
        call mobbrmsd_STRMV('Upper', 'No transpose', 'Non-unit', I - 1, T, &
            &      LDT, T(1, I), 1)
        T(I, I) = TAU(I)
        if (I > 1) then
          PREVLASTV = MAX(PREVLASTV, LASTV)
        else
          PREVLASTV = LASTV
        end if
      end if
    end do
  else
    PREVLASTV = 1
    do I = K, 1, -1
      if (TAU(I) == ZERO) then
!
! H(i) = I
!
        do J = I, K
          T(J, I) = ZERO
        end do
      else
!
! general case
!
        if (I < K) then
          if (mobbrmsd_LSAME(STOREV, 'C')) then
! Skip any leading zeros.
            do LASTV = 1, I - 1
              if (V(LASTV, I) /= ZERO) exit
            end do
            do J = I + 1, K
              T(J, I) = -TAU(I) * V(N - K + I, J)
            end do
            J = MAX(LASTV, PREVLASTV)
!
! T(i + 1:k, i) = -tau(i) * V(j:n - k + i, i + 1:k)**T * V(j:n - k + i, i)
!
            call mobbrmsd_SGEMV('Transpose', N - K + I - J, K - I, -TAU(I), &
                &      V(J, I + 1), LDV, V(J, I), 1, ONE, &
                &      T(I + 1, I), 1)
          else
! Skip any leading zeros.
            do LASTV = 1, I - 1
              if (V(I, LASTV) /= ZERO) exit
            end do
            do J = I + 1, K
              T(J, I) = -TAU(I) * V(J, N - K + I)
            end do
            J = MAX(LASTV, PREVLASTV)
!
! T(i + 1:k, i) = -tau(i) * V(i + 1:k, j:n - k + i) * V(i, j:n - k + i)**T
!
            call mobbrmsd_SGEMV('No transpose', K - I, N - K + I - J, &
                &      -TAU(I), V(I + 1, J), LDV, V(I, J), LDV, &
                &      ONE, T(I + 1, I), 1)
          end if
!
! T(i + 1:k, i): = T(i + 1:k, i + 1:k) * T(i + 1:k, i)
!
          call mobbrmsd_STRMV('Lower', 'No transpose', 'Non-unit', K - I, &
              &      T(I + 1, I + 1), LDT, T(I + 1, I), 1)
          if (I > 1) then
            PREVLASTV = MIN(PREVLASTV, LASTV)
          else
            PREVLASTV = LASTV
          end if
        end if
        T(I, I) = TAU(I)
      end if
    end do
  end if
  return
!
! end of mobbrmsd_SLARFT
!
end

