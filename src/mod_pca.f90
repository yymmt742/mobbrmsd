module mod_pca
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO
  use mod_optarg
  implicit none
  private
  public :: pca_worksize, pca
!
  interface
    include 'dgemm.h'
    include 'dsyev.h'
  end interface
!
contains
!
!| Calculates the optimal size of the WORK array
   pure function pca_worksize(d) result(res)
    integer(IK), intent(in) :: d
    !! matrix dimension
    real(RK)                :: dum(1)
    integer(IK)             :: res, info
!
      if (d < 1) then
        res = 0
      else
        call DSYEV('V', 'L', d, dum, d, dum, dum, -1, info)
        res = NINT(dum(1))
      endif
!
  end function pca_worksize
!
!| singular value decomposition of square matrix x.
  pure subroutine pca(T, d, n, x, u, v, w)
    logical, intent(in)     :: T
    !! transpose flag
    integer(IK), intent(in) :: d
    !! matrix dimension
    integer(IK), intent(in) :: n
    !! data dimension
    real(RK), intent(in)    :: x(d, *)
    !! d*n square matrix, on exit, x is destroyed.
    real(RK), intent(inout) :: u(*)
    !! d*d square matrix, returns eigen vector.
    real(RK), intent(inout) :: v(*)
    !! eigen value
    real(RK), intent(inout) :: w(*)
    !! work array, must be greater than svd_worksize(d).
    integer(IK)             :: lw, info
!
      if (T) then
        lw = pca_worksize(n)
        if (lw <= 0 .or. n <= 0) return
        if (d <= 0) then
          call zfill(n * n, u)
          call zfill(n, v)
          return
        end if
        call DGEMM('T', 'N', n, n, d,-ONE, x, d, x, d, ZERO, u, n)
        call DSYEV('V', 'L', n, u, n, v, w, lw, info)
        v(:n) = -v(:n)
      else
        lw = pca_worksize(d)
        if (lw <= 0 .or. d <= 0) return
        if (n <= 0) then
          call zfill(d * d, u)
          call zfill(d, v)
          return
        end if
        call DGEMM('N', 'T', d, d, n,-ONE, x, d, x, d, ZERO, u, d)
        call DSYEV('V', 'L', d, u, d, v, w, lw, info)
        v(:d) = -v(:d)
      end if
!
  end subroutine pca
!
  pure subroutine zfill(d, x)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      x(i) = ZERO
    end do
  end subroutine zfill
!
end module mod_pca
