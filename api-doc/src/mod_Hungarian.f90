module mod_Hungarian
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: Hungarian, Hungarian_value
!
contains
!
  pure function Hungarian_value(n, C) result(res)
    integer(IK), intent(in) :: n
    !! matrix dimension
    real(RK), intent(in)    :: C(*)
    !! n*n score matrix.
    real(RK)                :: res
    if (n < 1) then
      res = ZERO
      return
    elseif (n == 1) then
      res = C(1)
    elseif (n == 2) then
      if (C(1) + C(4) >= C(2) + C(3)) then
        res = C(1) + C(4)
      else
        res = C(2) + C(3)
      end if
    else
      block
        real(RK)    :: W(n + n + 1)
        integer(IK) :: iw(3 * (n + 1)), i, j
!
        call get_piv(n, n + 1, C, iw(1), iw(n + 2), iw(n + n + 3), W(1), W(n + 1))
        res = ZERO
        do i = 1, n
          j = i + n * (iw(i) - 1)
          res = res + C(j)
        end do
      end block
    end if
  end function Hungarian_value
!
!| Hungarian method
  pure  subroutine Hungarian(n, C, W, piv)
    integer(IK), intent(in)              :: n
    !! matrix dimension
    real(RK), intent(in)                 :: C(*)
    !! pivot index
    real(RK), intent(inout)              :: W(*)
    !! work array
    integer(IK), intent(inout), optional :: piv(*)
    !! n*n score matrix.
!
    if (n < 0)then
      ! query work array size
      W(1) = ABS(n + n) + 1
    elseif (n == 0)then
      return
    elseif (n == 1) then
      W(1) = C(1)
      if (PRESENT(piv)) piv(1) = 1
    elseif (n == 2) then
      W(1) = C(1) + C(4)
      W(2) = C(2) + C(3)
      if (W(1) >= W(2))then
        if (PRESENT(piv)) then
          piv(1) = 1
          piv(2) = 2
        end if
      else
        W(2) = W(1)
        if (PRESENT(piv)) then
          piv(1) = 2
          piv(2) = 1
        end if
      endif
    else
      block
        integer(IK) :: iw(n + n + n + 3), n1, i, j
!
        n1 = n + 1
        call get_piv(n, n1, C, iw(1), iw(n + 2), iw(n + n + 3), W(1), W(n + 1))
!
        W(1) = ZERO
        do i = 1, n
          j = i + n * (iw(i) - 1)
          W(1) = W(1) + C(j)
        end do
!
        if (PRESENT(piv)) then
          do concurrent(i=1:n)
            piv(i) = iw(i)
          end do
        end if
!
      end block
    end if
!
  end subroutine Hungarian
!
  pure subroutine get_piv(n, n1, C, piv, is_visited, prv, y, cij)
    integer(IK), intent(in)    :: n, n1
    real(RK), intent(in)       :: C(n, n)
    integer(IK), intent(inout) :: piv(n1), is_visited(n1), prv(n1)
    real(RK), intent(inout)    :: y(*), cij(*)
    real(RK)                   :: minc, edge, cedg
    integer(IK)                :: i, ic, ix, j
!
    do concurrent(i=1:n1)
      piv(i) = -1
    end do
!
    do concurrent(i=1:n1)
      y(i) = ZERO
    end do
!
    do j = 1, n
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
            edge = C(i, piv(ic)) - y(i)
            if (ic < n1) edge = edge - C(ic, piv(ic)) + y(ic)
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
