module mod_svd
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO
  implicit none
  private
  public :: svd_worksize, svd
!
  interface
    include 'dgesvd.h'
  end interface
!
contains
!
!| Calculates the optimal size of the WORK array
   pure function svd_worksize(d) result(res)
    integer(IK), intent(in) :: d
    !! matrix dimension
    real(RK)                :: dum(1)
    integer(IK)             :: res, info
!
      if(d<1)then
        res = 0
      else
        call DGESVD('A', 'A', d, d, dum, d, dum, dum, d, dum, d, dum, -1, info)
        res = NINT(dum(1))
      endif
!
  end function svd_worksize
!
!| singular value decomposition of square matrix x.
   pure subroutine svd(d, x, s, u, vt, w)
    integer(IK), intent(in) :: d
    !! matrix dimension
    real(RK), intent(inout) :: x(*)
    !! d*d square matrix, on exit, x is destroyed.
    real(RK), intent(inout) :: s(*)
    !! d vector, returns singular value.
    real(RK), intent(inout) :: u(*)
    !! d*d square matrix, returns left eigen vector.
    real(RK), intent(inout) :: vt(*)
    !! d*d square matrix, returns transpose right eigen vector.
    real(RK), intent(inout) :: w(*)
    !! work array, must be greater than svd_worksize(d).
    integer(IK)             :: lw, info
!
      lw = svd_worksize(d)
      if (lw <= 0 .or. d <= 0) return
      call DGESVD('A', 'A', d, d, x, d, s, u, d, vt, d, w, lw, info)
!
  end subroutine svd
!
end module mod_svd
