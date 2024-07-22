program main
  use mod_params, only: RK, IK
  use mod_dimspec_functions, only: D, DD, setup_dimension
  use mod_unittest
  use mod_testutil
  implicit none
  type(unittest) :: u
#ifdef USE_REAL32
  integer, parameter :: place = 3
#else
  integer, parameter :: place = 6
#endif
!
  call u%init('test testutils')
  call setup_dimension(4)
  call test1()
!
  call u%finish_and_terminate()
!
contains
  subroutine test1()
    real(RK) :: X(D, 10), X_(D, 10), Y(D, 3, 4)
    real(RK) :: CX(D), CY(D)
    real(RK) :: C(D, D), GC(DD + 1)
    real(RK) :: R(D, D), S(D, D), E(D, D)
    real(RK) :: G, T
    integer  :: i, j
    X = sample(SIZE(X, 2))
    Y = sample(SIZE(Y, 2), SIZE(Y, 3))
!   do i = 1, SIZE(X, 2)
!     print'(*(f9.3))', X(:, i)
!   end do
!   print *
!   do j = 1, SIZE(Y, 3)
!     do i = 1, SIZE(Y, 2)
!       print'(*(f9.3))', Y(:, i, j)
!     end do
!   end do
!   print *
!
    call centering(SIZE(X, 2), X)
    call centering(SIZE(Y, 2), SIZE(Y, 3), Y)
    CX = 0.0_RK
    CY = 0.0_RK
    do i = 1, SIZE(X, 2)
      CX = CX + X(:, i)
    end do
    do j = 1, SIZE(Y, 3)
      do i = 1, SIZE(Y, 2)
        CY = CY + Y(:, i, j)
      end do
    end do
    call u%assert_is_zero(CX, 'CX', place=place)
    call u%assert_is_zero(CY, 'CY', place=place)
!
    C = covmat(1)
    GC = gcov(1)
!   do i = 1, SIZE(C, 2)
!     print'(*(f9.3))', C(:, i)
!   end do
!   do i = 1, D
!     print'(*(f9.3))', GC(D * (i - 1) + 1:D * i)
!   end do
!
    R = SO()
!
    call u%assert_is_eye(MATMUL(R, TRANSPOSE(R)), 'R@RT', place=place)
    call u%assert_is_eye(MATMUL(TRANSPOSE(R), R), 'RT@R', place=place)
!
    E = eye()
    call u%assert_is_eye(E, 'Eye', place=place)
!
    X = sample(SIZE(X, 2))
    X_ = sample(SIZE(X_, 2))
    call centering(SIZE(X, 2), X)
    call centering(SIZE(X_, 2), X_)
    G = autovar(SIZE(X, 2), X, X_)
    call u%assert_almost_equal(G, SUM(X * X) + SUM(X_ * X_), 'autovar', place=place)
!
    do i = 1, 100
      C = covmat(10)
      if (D == 2) then
        T = C(1, 1) * C(2, 2) - C(1, 2) * C(2, 1)
      elseif (D == 3) then
        T = C(1, 1) * (C(2, 2) * C(3, 3) - C(2, 3) * C(3, 2)) &
       &  - C(1, 2) * (C(2, 1) * C(3, 3) - C(2, 3) * C(3, 1)) &
       &  + C(1, 3) * (C(2, 1) * C(3, 2) - C(2, 2) * C(3, 1))
      else
        T = C(1, 1) * (C(2, 2) * (C(3, 3) * C(4, 4) - C(3, 4) * C(4, 3)) &
       &             - C(2, 3) * (C(3, 2) * C(4, 4) - C(3, 4) * C(4, 2)) &
       &             + C(2, 4) * (C(3, 2) * C(4, 3) - C(3, 3) * C(4, 2))) &
       &  - C(2, 1) * (C(1, 2) * (C(3, 3) * C(4, 4) - C(3, 4) * C(4, 3)) &
       &             - C(1, 3) * (C(3, 2) * C(4, 4) - C(3, 4) * C(4, 2)) &
       &             + C(1, 4) * (C(3, 2) * C(4, 3) - C(3, 3) * C(4, 2))) &
       &  + C(3, 1) * (C(1, 2) * (C(2, 3) * C(4, 4) - C(2, 4) * C(4, 3)) &
       &             - C(1, 3) * (C(2, 2) * C(4, 4) - C(2, 4) * C(4, 2)) &
       &             + C(1, 4) * (C(2, 2) * C(4, 3) - C(2, 3) * C(4, 2))) &
       &  - C(4, 1) * (C(1, 2) * (C(2, 3) * C(3, 4) - C(2, 4) * C(3, 3)) &
       &             - C(1, 3) * (C(2, 2) * C(3, 4) - C(2, 4) * C(3, 2)) &
       &             + C(1, 4) * (C(2, 2) * C(3, 3) - C(2, 3) * C(3, 2)))
      end if
      G = det_sign(C)
      call u%assert_almost_equal(SIGN(1.0_RK, T), G, 'det', place=place)
    end do
!
    do i = 1, 100
      X = sample(SIZE(X, 2))
      R = SO()
      X_ = MATMUL(R, X)
      call kabsch(MATMUL(X, TRANSPOSE(X_)), S)
      call u%assert_almost_equal([R], [TRANSPOSE(S)], 'kabsch', place=place)
    end do
!
  end subroutine test1
end program main

