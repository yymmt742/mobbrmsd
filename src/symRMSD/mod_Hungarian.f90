!
!| Calculate the minimum linear assignment cost using Hungarian method.
module mod_Hungarian
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: Hungarian_worksize
  public :: Hungarian
!
contains
!
!| query work array size
  pure elemental function Hungarian_worksize(m, n) result(res)
    !| matrix dimension 1.
    integer(IK), intent(in) :: m
    !| matrix dimension 2.
    integer(IK), intent(in) :: n
    integer(IK)             :: res
    if (n > 0 .and. m > 0) then
      res = MAX(n, m) * 2 + 5
    else
      res = 0
    end if
  end function Hungarian_worksize
!
!| Calculate the minimum linear assignment cost using Hungarian method.
!  If m and n are different, the sum of the linear assignments of the smaller is returned.
!  If m>0 and n>0, W(1) stores the minimum linear assignment cost.
!  If m==0 or n==0, do nothing.
!  If (m<0 or n<0) and (|m|>0 and |n|>0), W(1) stores the required memory size for C(|m|,|n|).
  pure subroutine Hungarian(m, n, C, W)
    !| matrix dimension 1.
    integer(IK), intent(in)    :: m
    !| matrix dimension 2.
    integer(IK), intent(in)    :: n
    !| score matrix C(m, n).
    real(RK), intent(in)       :: C(*)
    !| work array.
    real(RK), intent(inout)    :: W(*)
!
    if (n == 0 .or. m == 0) then
      return
    elseif (n < 0 .or. m < 0) then
      ! query work array size
      W(1) = MAX(ABS(n), ABS(m)) * 2 + 5
    elseif (m <= n)then
      block
        integer(IK) :: iw(n + n + n + 3)
        call get_piv(m, n, C, iw(1), iw(n + 2), iw(n + n + 3), &
       &             W(5), W(n + 5), W(2), W(3), W(4), W(1))
      end block
    elseif (m > n)then
      block
        integer(IK) :: iw(m + m + m + 3)
        call get_piv_T(m, n, C, iw(1), iw(m + 2), iw(m + m + 3), &
       &               W(5), W(m + 5), W(2), W(3), W(4), W(1))
      end block
    end if
!
  end subroutine Hungarian
!
  pure subroutine get_piv(m, n, C, piv, is_visited, prv, y, cij, minc, edge, cedg, res)
    integer(IK), intent(in)    :: m, n
    real(RK), intent(in)       :: C(m, n)
    integer(IK), intent(inout) :: piv(*), is_visited(*), prv(*)
    real(RK), intent(inout)    :: y(*), cij(*), minc, edge, cedg, res
    integer(IK)                :: l, i, j, ic, ix
!
    l = n + 1
    do concurrent(i=1:l)
      piv(i) = -1
    end do
!
    do concurrent(i=1:l)
      y(i) = ZERO
    end do
!
    res = ZERO
!
    do j = 1, m
!
      do concurrent(i=1:l)
        is_visited(i) = 0
      end do
!
      do concurrent(i=1:l)
        prv(i) = -1
      end do
!
      do concurrent(i=1:n)
        cij(i) = RHUGE
      end do
      cij(l) = ZERO
!
      ic = l
      piv(ic) = j
!
      do while (piv(ic) /= -1)
        minc = RHUGE
        is_visited(ic) = 1
        ix = -1
        do i = 1, n
          if (is_visited(i) == 0) then
            edge = C(piv(ic), i) - y(i)
            if (ic < l) edge = edge - C(piv(ic), ic) + y(ic)
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
      do while(ic /= l)
        i = prv(ic)
        piv(ic) = piv(i)
        ic = i
      end do
!
    end do
!
  end subroutine get_piv
!
  pure subroutine get_piv_T(m, n, C, piv, is_visited, prv, y, cij, minc, edge, cedg, res)
    integer(IK), intent(in)    :: m, n
    real(RK), intent(in)       :: C(m, n)
    integer(IK), intent(inout) :: piv(*), is_visited(*), prv(*)
    real(RK), intent(inout)    :: y(*), cij(*), minc, edge, cedg, res
    integer(IK)                :: l, i, j, ic, ix
!
    l = m + 1
    do concurrent(i=1:l)
      piv(i) = -1
    end do
!
    do concurrent(i=1:l)
      y(i) = ZERO
    end do
!
    res = ZERO
!
    do j = 1, n
!
      do concurrent(i=1:l)
        is_visited(i) = 0
      end do
!
      do concurrent(i=1:l)
        prv(i) = -1
      end do
!
      do concurrent(i=1:m)
        cij(i) = RHUGE
      end do
      cij(l) = ZERO
!
      ic = l
      piv(ic) = j
!
      do while (piv(ic) /= -1)
        minc = RHUGE
        is_visited(ic) = 1
        ix = -1
        do i = 1, m
          if (is_visited(i) == 0) then
            edge = C(i, piv(ic)) - y(i)
            if (ic < l) edge = edge - C(ic, piv(ic)) + y(ic)
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
      do concurrent(i = 1:m)
        if(i/=ic) cij(i) = MIN(cij(i), cij(ic))
      end do
      do concurrent(i = 1:m)
        y(i) = y(i) + cij(i)
      end do
!
      res = res + y(ic)
!
      do while(ic /= l)
        i = prv(ic)
        piv(ic) = piv(i)
        ic = i
      end do
!
    end do
!
  end subroutine get_piv_T
!
end module mod_Hungarian
