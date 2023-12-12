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
   subroutine simultaneous_diagonalize(n, m, X, U, L, rwrk)
    integer(IK), intent(in) :: n
    !!! matrix dimension
    integer(IK), intent(in) :: m
    !!! number of matrices
    real(RK), intent(in)    :: X(n * n, *)
    !!! matcies X(n, n, m)
    real(RK), intent(inout) :: U(*)
    !!! simultaneous right eigenvector U(n, n)
    real(RK), intent(inout) :: L(n, *)
    !!! eigenvalues L(n, m), L(:, i) corresponding to X(:,:,i).
    real(RK), intent(inout) :: rwrk(*)
    !!! work array
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
    call copy_to_lower(n, n, X, U)
!
!   decompose P**T (-X) P = (-L'), P is stored to U.
    call DSYEV('V', 'L', n, U, n, L, rwrk, lw, info)
!
    do j = 2, m
!     get P**T (-Y) P = -R ( R is symmetric, because X is symmetric ), (-R) is stored to rwrk(1:)
      call diagonalize(n, n, n, U, X(1, j), rwrk(1), rwrk(nn + 1))
      !call diagonalize(n, U, X(1, 1, j), rwrk)
!
!     decompose Q**T (-R) Q = (-DY), Q is stored to rwrk(1:)
      call decompose_R(n, lw, L, 1D-4, rwrk, U, L(1, j), rwrk(nn + 1))
!
print'(10f6.1)',rwrk(:n*n)
print*
return
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
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: U(n, n)
    real(RK), intent(inout) :: V(n)
    integer(IK)             :: i
!
    V(:) = [(DOT_PRODUCT(V(:) * U(:, i), U(:, i)), i=1, n)]
!
  end subroutine proc_vector
!
  pure subroutine copy_to_lower(n, lx, X, Z)
    integer(IK), intent(in) :: n, lx
    real(RK), intent(in)    :: X(lx, *)
    real(RK), intent(inout) :: Z(n, *)
    integer(IK)             :: j
    do concurrent(j=1:n)
      block
        integer(IK) :: i
        do concurrent(i=j:n)
          Z(i, j) = -X(i, j)
        end do
      end block
    end do
  end subroutine copy_to_lower
!
  pure subroutine diagonalize(n, lu, lx, U, X, D, W)
    integer(IK), intent(in) :: n, lu, lx
    real(RK), intent(in)    :: U(*), X(*)
    real(RK), intent(inout) :: D(*), W(*)
    call DGEMM('T', 'N', n, n, n, ONE, U, lu, X, lx, ZERO, W, n)
    call DGEMM('N', 'N', n, n, n, ONE, W, n, U, lu, ZERO, D, n)
  end subroutine diagonalize
!
  pure subroutine decompose_R(n, lw, LY, deg, R, U, LX, W)
    integer(IK), intent(in) :: n, lw
    real(RK), intent(in)    :: LY(n), deg
    real(RK), intent(inout) :: R(n, n), U(n, *), LX(n), W(*)
    integer(IK)             :: i, lb, ub, d, dd
!
    ub = 0
    lb = 1
!
    do while (ub < n)
      ub = ub + 1
      do while (ub < n)
        if (ABS(LY(ub + 1) - LY(lb)) > deg) exit
        ub = ub + 1
      end do
!
      if (ub > lb) then
        d = ub - lb + 1
        dd = d * d
        call decompose_block(d, n, lw - dd - d, W, W(dd + 1), W(dd + d + 1))
        do concurrent(i=1:d)
          LX(lb + i - 1) = W(dd + i)
        end do
        call DGEMM('N', 'N', n, d, d, ONE, U(1, lb), n, W, d, ZERO, W(dd + 1), n)
!       get U = PQ
        call copy(n * d, W(dd + 1), U(1, lb))
        lb = ub + 1
      else
        LX(lb) = R(lb, lb)
      end if
    end do
!
  end subroutine decompose_R
!
  pure subroutine decompose_block(d, n, lw, R, LX, W)
    integer(IK), intent(in) :: d, n, lw
    real(RK), intent(inout) :: R(*), LX(d), W(*)
    integer(IK)             :: i, dd, info
!
    dd = d * d
    call copy_to_lower(d, d, R, W)
    call DSYEV('V', 'L', d, W(1), d, W(dd + 1), W(dd + d + 1), lw, info)
    call copy(dd, W, R)
!
    do concurrent(i=1:d)
      LX(i) = -W(dd + i)
    end do
!
  end subroutine decompose_block
!
  pure subroutine copy(n, X, Z)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Z(*)
    integer(IK)             :: i
    do concurrent(i=1:n)
      Z(i) = X(i)
    end do
  end subroutine copy
!
end module mod_diagonalize
