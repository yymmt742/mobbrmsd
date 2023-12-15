module mod_Hungarian
  use mod_params, only : IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: Hungarian
!
contains
!
!| Hungarian
  subroutine Hungarian(n, C, P)
    integer(IK), intent(in) :: n
    !! matrix dimension
    real(RK), intent(in)    :: C(n, *)
    !! n*n score matrix.
    real(RK), intent(inout) :: P(n, *)
    !! n*n permutation matrix
!
    if (n < 1)then
      return
    elseif (n == 1) then
      P(1, 1) = ONE
    elseif (n == 2) then
      if (C(1, 1) + C(2, 2) >= C(1, 2) + C(2, 1))then
        P(1, 1) = ONE
        P(2, 1) = ZERO
        P(1, 2) = ZERO
        P(2, 2) = ONE
      else
        P(1, 1) = ZERO
        P(2, 1) = ONE
        P(1, 2) = ONE
        P(2, 2) = ZERO
      endif
    else
      block
        integer(IK) :: iw(3 * (n + 1)), i, j
!
        call get_piv(n, n + 1, C, iw(1), iw(n + 2), iw(n + n + 3), P(1, 1), P(1, 2))
!
        do concurrent(i=1:n, j=1:n)
          P(i, j) = MERGE(ONE, ZERO, iw(j) == i)
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
