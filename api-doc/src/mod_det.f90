module mod_det
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: det, det_sign
!
  interface
    include 'dgetrf.h'
  end interface
!
  interface det
    module procedure :: det_full, det_part
  end interface det
!
  interface det_sign
    module procedure :: det_sign_copy, det_sign_full, det_sign_part
  end interface det_sign
!
contains
!
!| calculate determinant of square matrix x.
   pure subroutine det_full(d, x)
     integer(IK), intent(in) :: d
    !! matrix dimension
     real(RK), intent(inout) :: x(*)
    !! square matrix, on exit, x(1) is assigned the determinant of x,
    !! and the other elements are undefined.
!
     if (d <= 1) then
!
       return
!
     elseif (d == 2) then
!
       x(1)  = x(1) * x(4) - x(2) * x(3)
!
     elseif (d == 3) then
!
       x(1) = x(1) * (x(5) * x(9) - x(8) * x(6)) +&
            & x(4) * (x(8) * x(3) - x(2) * x(9)) +&
            & x(7) * (x(2) * x(6) - x(5) * x(3))
!
     else
!
       call det_part(d, x, d)
!
     end if
!
   end subroutine det_full
!
!| calculate determinant of square matrix x, with leading dimension.
   pure subroutine det_part(d, x, ld)
     integer(IK), intent(in) :: d
    !! matrix dimension
     real(RK), intent(inout) :: x(*)
    !! square matrix, on exit, x(1) is assigned the determinant of x,
    !! and the other elements are undefined.
     integer(IK), intent(in) :: ld
    !! leading dimension
!
     if (d <= 1) then
!
       return
!
     elseif (d == 2) then
       block
         integer(IK) :: l2
         l2 = MAX(d, ld) + 1
         x(1) = x(1) * x(l2 + 1) - x(2) * x(l2)
       end block
     elseif (d == 3) then
       block
         integer(IK) :: l2, l3
         l2 = MAX(d, ld)
         l3 = l2 + l2
         l2 = l2 + 1
         l3 = l3 + 1
!&<
         x(1) = x( 1) * (x(1 +l2) * x(2 +l3) - x(1 +l3) * x(2 +l2)) +&
              & x(l2) * (x(1 +l3) * x(2 + 1) - x(1 + 1) * x(2 +l3)) +&
              & x(l3) * (x(1 + 1) * x(2 +l2) - x(1 +l2) * x(2 + 1))
!>&
       end block
     else
       block
         integer(IK) :: i, j, k, ipiv(d)
         k = MAX(d, ld)
         call DGETRF(d, d, x, k, ipiv(1:d), j)
         if (MODULO(COUNT([(ipiv(i) == i, i=1, d)]), 2) == 1) x(1) = -x(1)
         ipiv(1) = k + 1
         k = k * d
         do i = 2, d
           j = k - ipiv(1)
           x(j) = x(j) * x(k)
           k = j
         end do
       end block
     end if
!
   end subroutine det_part
!
!| calculate determinant sign of square matrix x, with leading dimension.
!
   pure subroutine det_sign_copy(d, x, w)
     integer(IK), intent(in) :: d
    !! matrix dimension
     real(RK), intent(in)    :: x(*)
    !! d * d square matrix.
     real(RK), intent(inout) :: w(*)
    !! work array, on exit, w(1) is assigned the determinant sign of x.
     w(:d * d) = x(:d * d)
     call det_sign_full(d, w)
   end subroutine det_sign_copy
!
   pure subroutine det_sign_full(d, x)
     integer(IK), intent(in) :: d
    !! matrix dimension
     real(RK), intent(inout) :: x(*)
    !! square matrix, on exit, x(1) is assigned the determinant sign of x, <br>
    !! and the other elements are undefined.
!
     if (d < 1) then
!
       return
!
     elseif (d == 1) then
!
       x(1) = SIGN(ONE, x(1))
!
     elseif (d == 2) then
!
       x(1) = SIGN(ONE, x(1) * x(4) - x(2) * x(3))
!
     elseif (d == 3) then
!
       x(1) = SIGN(ONE, x(1) * (x(5) * x(9) - x(8) * x(6)) +&
         &              x(4) * (x(8) * x(3) - x(2) * x(9)) +&
         &              x(7) * (x(2) * x(6) - x(5) * x(3))  )
!
     else
!
       call det_sign_part(d, x, d)

     end if
!
   end subroutine det_sign_full
!
!| calculate determinant sign of square matrix x, with leading dimension.
!
   pure subroutine det_sign_part(d, x, ld)
     integer(IK), intent(in) :: d
    !! matrix dimension
     real(RK), intent(inout) :: x(*)
    !! square matrix, on exit, x(1) is assigned the determinant sign of x, <br>
    !! and the other elements are undefined.
     integer(IK), intent(in) :: ld
    !! leading dimension
!
     if (d < 1) then
!
       return
!
     elseif (d == 1) then
!
       x(1) = SIGN(ONE, x(1))
!
     elseif (d == 2) then
!
       block
         integer(IK) :: l2
         l2 = MAX(d, ld) + 1
         x(1) = SIGN(ONE, x(1) * x(l2 + 1) - x(2) * x(l2))
       end block
!
     elseif (d == 3) then
!
       block
         integer(IK) :: l2, l3
         l2 = MAX(d, ld)
         l3 = l2 + l2
         l2 = l2 + 1
         l3 = l3 + 1
!&<
         x(1) = SIGN(ONE, x( 1) * (x(1 +l2) * x(2 +l3) - x(1 +l3) * x(2 +l2)) +&
           &              x(l2) * (x(1 +l3) * x(2 + 1) - x(1 + 1) * x(2 +l3)) +&
           &              x(l3) * (x(1 + 1) * x(2 +l2) - x(1 +l2) * x(2 + 1)))
!>&
       end block
!
     else
!
       block
         integer(IK) :: i, j, k, ipiv(d)
         k = MAX(d, ld)
         call DGETRF(d, d, x, k, ipiv, j)
         ipiv(1) = COUNT([(ipiv(i) == i, i=1, d)])
         j = 1
         k = k + 1
         do i = 1, d
           if (x(j) <= ZERO) ipiv(1) = ipiv(1) + 1
           j = j + k
         end do
         if( MODULO(ipiv(1), 2)==0)then
           x(1) = ONE
         else
           x(1) = -ONE
         endif
       end block
!
     end if
!
   end subroutine det_sign_part
!
 end module mod_det
