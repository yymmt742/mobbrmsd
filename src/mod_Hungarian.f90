module mod_Hungarian
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: Hungarian
!
contains
!
!| Hungarian method
  pure subroutine Hungarian(n, C, piv, W)
    integer(IK), intent(in)    :: n
    !! matrix dimension
    real(RK), intent(in)       :: C(*)
    !! n*n score matrix.
    integer(IK), intent(inout) :: piv(*)
    !! pivot index
    real(RK), intent(inout)    :: W(*)
    !! work array
!
    if (n < 0)then
      ! query work array size
      W(1) = ABS(n + n) + 1
    elseif (n == 0)then
      return
    elseif (n == 1) then
      piv(1) = 1
    elseif (n == 2) then
      if (C(1) + C(4) >= C(2) + C(3))then
        piv(1) = 1
        piv(2) = 2
      else
        piv(1) = 2
        piv(2) = 1
      endif
    else
      block
        integer(IK) :: iw(3 * (n + 1)), i
!
        call get_piv(n, n + 1, C, iw(1), iw(n + 2), iw(n + n + 3), W(1), W(n+1))
!
        do concurrent(i=1:n)
          piv(i) = iw(i)
        end do
!
      end block
    end if
!
  end subroutine Hungarian
!
  pure subroutine get_piv(n, n1, C, piv, vis, prv, y, cij)
    integer(IK), intent(in)    :: n, n1
    real(RK), intent(in)       :: C(n, n)
    integer(IK), intent(inout) :: piv(n1), vis(n1), prv(n1)
    real(RK), intent(inout)    :: y(*), cij(*)
    integer(IK)                :: i, ic, ix, j
!
    piv = -1
!
    do j = 1, n
!
      do concurrent(i=1:n1)
        vis(i) = 0
      end do
!
      do concurrent(i=1:n1)
        prv(i) = -1
      end do
!
      do concurrent(i=1:n)
        cij(i) = RHUGE
      end do
      cij(n1) = ZERO
!
      ic = n1
      piv(ic) = j
!
      do while (piv(ic) /= -1)
        block
          real(RK) :: minc, edge
          minc = RHUGE
          vis(ic) = 1
          ix = -1
          do i = 1, n
            if (vis(i)==0) then
              edge = -C(i, piv(ic)) - y(i)
              if (ic /= n1) edge = edge + C(ic, piv(ic)) + y(ic)
              if (cij(i) > cij(ic) + edge) then
                prv(i) = ic
                cij(i) = cij(ic) + edge
              end if
              if (minc > cij(i)) then
                ix = i
                minc = cij(i)
              end if
            end if
          end do
        end block
        ic = ix
      end do
!
      do concurrent(i = 1:n)
        if(i/=ic) cij(i) = MIN(cij(i), cij(ic))
      end do
      do concurrent(i = 1:n)
        y(i) = y(i) + cij(i)
      end do
!
      do while(ic /= n1)
        i = prv(ic)
        piv(ic) = piv(i)
        ic = i
      end do
!
    end do
!
  end subroutine get_piv
!
end module mod_Hungarian
