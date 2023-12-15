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
    real(RK), intent(in)    :: C(n, n)
    !! d*n square matrix, on exit, x is destroyed.
    real(RK), intent(inout) :: P(n, n)
    !! n*m permutation matrix
    real(RK)                :: h(n), answers(n), ans_cur
    integer(IK)             :: j_cur, w_cur, job(n + 1)
!
    P = MAXVAL(C) - C
!   print'(10f9.3)',P
!   print*
!   call sweep(n, P)
    print'(10f9.3)',P
    print*
    job = -1
    ans_cur = ZERO
    h = ZERO
!
    do j_cur = 1, n
      block
        logical     :: vis(n + 1)
        integer(IK) :: prv(n + 1)
        real(RK)    :: dist(n), min_dist, edge
        integer(IK) :: w, w_next
        w_cur = n
        job(w_cur) = j_cur
        dist = RHUGE
        dist(n) = ZERO
        prv = -1
print*,j_cur, w_cur, job
        do
        !do while (job(w_cur)/=-1)
          min_dist = RHUGE
          vis(w_cur) = .TRUE.
          w_next = -1
print*, vis
          do w = 1, n
print*, w, w_cur, vis(w)
            if (.not. vis(w)) then
              edge = P(w, job(w_cur)) - h(w)
              if (w_cur/=n) then
                edge = edge - P(w_cur, job(w_cur)) - h(w_cur)
              end if
              if (dist(w) > dist(w_cur) + edge) then
                dist(w) = dist(w_cur) + edge
                prv(w) = w_cur
              end if
              if (min_dist > dist(w)) then
                min_dist = dist(w)
                w_next = w
              end if
              print'(2i4,3f9.3)',w, job(w_cur), P(w, job(w_cur)), h(w), edge
            end if
          end do
print*, w_next
          w_cur = w_next
          if(job(w_cur)==-1) exit
        end do
print*,prv
print*,w_cur
        do w = 1, n
          if (dist(w) > dist(w_cur)) then
            dist(w) = dist(w_cur)
            h(w) = h(w) + dist(w)
          end if
        end do
        ans_cur = ans_cur + dist(w_cur)
        w = prv(w_cur)
        do
          w = prv(w_cur)
          job(w_cur) = job(w)
          w_cur = w
          if(w_cur == n) exit
        end do
      end block
      answers(j_cur) = ans_cur
    end do
!
  end subroutine Hungarian

!
  pure subroutine sweep(n, P)
    integer(IK), intent(in) :: n
    real(RK), intent(inout) :: P(n, n)
    integer(IK)             :: i
    do concurrent(i=1:n)
      block
        real(RK) :: mv
        mv = MINVAL(P(:, i))
        P(:, i) = P(:, i) - mv
      end block
    end do
    do concurrent(i=1:n)
      block
        real(RK) :: mv
        mv = MINVAL(P(i, :))
        P(i, :) = P(i, :) - mv
      end block
    end do
  end subroutine sweep
!
end module mod_Hungarian
