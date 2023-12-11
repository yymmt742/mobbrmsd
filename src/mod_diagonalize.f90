module mod_diagonalize
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: simultaneous_diagonalize
!
  interface
    include 'dgemm.h'
    include 'dsyev.h'
  end interface
!
contains
!
!| Calculates the optimal size of the WORK array
  pure subroutine simultaneous_diagonalize(n, m, X, U, L, rwrk)
    integer(IK), intent(in) :: n, m
    real(RK), intent(in)    :: X(n, n, *)
    real(RK), intent(inout) :: U(n, *), L(n, *)
    real(RK), intent(inout) :: rwrk(*)
    integer(IK)             :: i, j, info, lw, nn
!
    nn = n * n
!
!   if n<1, returns the memory size corresponding to |n|.
!
    call DSYEV('V', 'L', ABS(n), rwrk, ABS(n), rwrk, rwrk, -1, info)
    lw = rwrk(1) + 2 * nn
    if (n < 1) then
      rwrk(1) = lw
      return
    end if
!
!   copy -X to P(n,n) in lower matrix
    call copy_to_lower(n, n, 1, X, U)
!
!   decompose P**T (-X) P = (-DX'), P is stored to U.
    call DSYEV('V', 'L', n, U, n, L, rwrk, lw, info)
!
    do j = 2, m
!   get P**T (-Y) P = -R ( R is symmetric, because X is symmetric ), (-R) is stored to rwrk(1:)
      call diagonalize(n, U, X(1, 1, j), rwrk)
!
!   decompose Q**T (-R) Q = (-DY), Q is stored to rwrk(1:)
      call decompose_R(n, lw, L, 1D-4, rwrk, L(1, j), rwrk(nn + 1))
!     call DSYEV( 'V', 'L', n, rwrk, n, DY, rwrk(nn+1:), SIZE(rwrk(nn+1:)), info )
!
!   get U = PQ
!
      U(:,:n) = MATMUL(U(:,:n), RESHAPE(rwrk(:nn), [n, n]))
!
!   get (-DY) = Q**T (-DY') Q = U**T (-Y) U
!
      call proc_vector(n, rwrk, L(1, j))
!
!   change signs.
    end do
!
    do concurrent(i=1:n, j=1:m)
      L(i, j) = -L(i, j)
    end do
!
  end subroutine simultaneous_diagonalize
!
  pure subroutine proc_vector(n, U, V)
    integer(IK), intent(in)    :: n
    real(RK), intent(in)    :: U(n, n)
    real(RK), intent(inout) :: V(n)
    integer(IK)                 :: i
!
    V(:) = [(DOT_PRODUCT(V(:) * U(:, i), U(:, i)), i=1, n)]
!
  end subroutine proc_vector
!
  pure subroutine copy_to_lower(n, m, l, S, D)
    integer(IK), intent(in) :: n, m, l
    real(RK), intent(in)    :: S(m, n)
    real(RK), intent(inout) :: D(n, n)
    integer(IK)             :: i, j, k
!
    j = l - 1; k = l + n - 1
    do i = 1, n
      D(i:n, i) = -S(i + j:k, i)
    end do
!
  end subroutine copy_to_lower
!
  pure subroutine decompose_R(n, lw, LY, deg, R, LX, W)
    integer(IK), intent(in) :: n, lw
    real(RK), intent(in)    :: LY(n), deg
    real(RK), intent(inout) :: R(n, n), LX(n), W(*)
    integer(IK)             :: i, l, u
!
    u = 0
!
    do while (u < n)
      l = u + 1
      do while (u < n)
        if (ABS(LY(u + 1) - LY(l)) > deg) exit
        u = u + 1
      end do
!
      if (u > l) then
        call decompose_block(l, u, u - l + 1, n, lw, R(1, l), LX(l), W)
      else
        LX(l) = R(l, l)
        do concurrent(i=1:n)
          R(i, l) = MERGE(ONE, ZERO, i==l)
        enddo
      end if
    end do
!
  end subroutine decompose_R
!
  pure subroutine decompose_block(l, u, d, n, lw, R, LX, W)
    integer(IK), intent(in) :: l, u, d, n, lw
    real(RK), intent(inout) :: R(n, d), LX(d), W(*)
    integer(IK)             :: i, dd, info
!
    dd = d * d
!
    call copy_to_lower(d, n, l, R, W)
!
    W(:dd) = -W(:dd)
!
    call DSYEV('V', 'L', d, W(1), d, W(dd + 1), W(dd + d + 1), lw, info)
!
!   R(l:u,:) = RESHAPE( W( 1:dd ), [d,d] )
    do concurrent(i = 1:d)
      R(l:u, i) = W((i - 1) * d + 1:i * d)
    end do
!
    do concurrent(i=1:d)
      LX(i) = -W(dd + i)
    end do
!
  end subroutine decompose_block
!
  pure subroutine diagonalize(n, U, M, D)
    integer(IK), intent(in)    :: n
    real(RK), intent(in)    :: U(n, n), M(n, n)
    real(RK), intent(inout) :: D(n, n)
    D = MATMUL(MATMUL(TRANSPOSE(U), M), U)
  end subroutine diagonalize
!
end module mod_diagonalize
