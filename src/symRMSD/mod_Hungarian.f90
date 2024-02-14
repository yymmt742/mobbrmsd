module mod_Hungarian
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: Hungarian, Hungarian_value
!
contains
!
!| Calculate linear assignment minimum cost using Hungarian method.
  pure function Hungarian_value(m, n, C) result(res)
    integer(IK), intent(in) :: m
    !! matrix dimension 1, must be m > 0.
    integer(IK), intent(in) :: n
    !! matrix dimension 2, must be n > 0 and m <= n.
    real(RK), intent(in)    :: C(*)
    !! score matrix C(m, n).
    real(RK)                :: res
    if (m < 1) then
      res = ZERO
      return
    elseif (m > n) then
      res = ZERO
      return
    elseif (m == 1) then
      res = MINVAL(C(:n))
    else
      block
        real(RK)    :: W(n + n + 1)
        integer(IK) :: iw(n + n + n + 3)
        call get_piv(m, n, C, iw(1), iw(n + 2), iw(n + n + 3), W(1), W(n + 1), res)
      end block
    end if
  end function Hungarian_value
!
!| Calculate linear assignment minimum cost using Hungarian method with pivot.
  pure subroutine Hungarian(m, n, C, W)
    integer(IK), intent(in)    :: m
    !! matrix dimension 1, must be m > 0.
    integer(IK), intent(in)    :: n
    !! matrix dimension 2, must be n > 0 and m <= n.
    real(RK), intent(in)       :: C(*)
    !! score matrix C(m, n).
    real(RK), intent(inout)    :: W(*)
    !! work array.
!
    if (n < 0)then
      ! query work array size
      W(1) = ABS(n + n) + 2
    elseif (m < 1.or.m > n)then
      return
    elseif (m == 1) then
      W(1) = MINVAL(C(:n))
    else
      block
        integer(IK) :: iw(n + n + n + 3)
        call get_piv(m, n, C, iw(1), iw(n + 2), iw(n + n + 3), W(2), W(n + 2), W(1))
      end block
    end if
!
  end subroutine Hungarian
!
  pure subroutine get_piv(m, n, C, piv, is_visited, prv, y, cij, res)
    integer(IK), intent(in)    :: m, n
    real(RK), intent(in)       :: C(m, n)
    integer(IK), intent(inout) :: piv(*), is_visited(*), prv(*)
    real(RK), intent(inout)    :: y(*), cij(*), res
    real(RK)                   :: minc, edge, cedg
    integer(IK)                :: n1, i, j, ic, ix
!
    n1 = n + 1
    do concurrent(i=1:n1)
      piv(i) = -1
    end do
!
    do concurrent(i=1:n1)
      y(i) = ZERO
    end do
!
    res = ZERO
!
    do j = 1, m
!
      do concurrent(i=1:n1)
        is_visited(i) = 0
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
        minc = RHUGE
        is_visited(ic) = 1
        ix = -1
        do i = 1, n
          if (is_visited(i) == 0) then
            edge = C(piv(ic), i) - y(i)
            if (ic < n1) edge = edge - C(piv(ic), ic) + y(ic)
            cedg = cij(ic) + edge
            if (cij(i) > cedg) then
              prv(i) = ic
              cij(i) = cedg
            end if
            if (minc > cij(i)) then
              ix = i
              minc = cij(i)
            end if
          end if
        end do
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
      res = res + y(ic)
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
