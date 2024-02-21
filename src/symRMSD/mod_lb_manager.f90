!| Module for manage lowerbound matrices.
module mod_lb_manager
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: gemm=>DGEMM, dot=>DDOT, copy=>DCOPY
  use mod_mol_block
  use mod_c_matrix
  use mod_s_matrix
  use mod_lowerbound
  implicit none
  private
  public :: worksize_lb_manager
  public :: lb_manager
!
!| lb_manager<br>
!    By default, memory is allocated as follows.<br>
!    |-C1-|----W1----|<br>
!    |----|--C2--|--------W2-------|<br>
!    |-----------|--C3--|-W3-|<br>
!    Therefore, the maximum memory allocation size is MAX( SUM_i^I |Ci| + |W_I| ).
  type lb_manager
    private
    sequence
    type(c_matrix) :: c
    type(s_matrix) :: s
  end type lb_manager
!
  interface lb_manager
    module procedure lb_manager_new
  end interface lb_manager
!
contains
!
!| Constructer
  pure elemental function lb_manager_new(b) result(res)
    !| b :: mol_block, must be initialized.
    type(mol_block), intent(in) :: b
    type(c_matrix)              :: res
!
    res%c  = c_matrix(b)
    res%s  = s_matrix(b)
!
  end function lb_manager_new
!
!| Inquire worksize of s_matrix.
  pure elemental function worksize_lowerbound(p, b) result(res)
    !| p :: level
    integer(IK), intent(in)    :: p
    !| b :: mol_block
    type(mol_block),intent(in) :: b
    integer(IK)                :: res
    integer(IK)                :: n1, n2
    n1 = MAX(b%x%n, b%y%n) - p
    n2 = MIN(b%x%n, b%y%n) - p
    res = MAX(worksize_Hungarian(n1, n2), 1 + worksize_sdmin())
  end function worksize_lowerbound
!
! pure subroutine d_matrix_partial_eval(a, p, iprm, isym, ires, W, LT, H, C, LF, LB)
!   type(d_matrix), intent(in)        :: a
!   integer(IK), intent(in)           :: p, iprm, isym, ires(*)
!   real(RK), intent(in)              :: W(*)
!   real(RK), intent(inout)           :: LT, H, C(*)
!   real(RK), intent(inout), optional :: LF, LB
!
!   LT = ZERO
!   if (PRESENT(LF)) LF = ZERO
!   if (PRESENT(LB)) LB = ZERO
!
!   if (p < 0 .or. a%g < p) return
!   if (p < a%g) call setminus_eval(a, p, ires, W, LT)
!   if (PRESENT(LB)) LB = LT
!
!   if (p == 0) return
!   if (iprm < 1 .or. a%g < iprm) return
!   if (isym < 0 .or. a%s <= isym) return
!
!   block
!     integer(IK) :: ih, ic
!     ih = a%c + (iprm - 1) * a%cb + (p - 1) * a%cl
!     ic = ih + 1 + DD * isym
!     if (PRESENT(LF)) then
!       call partial_eval(a%s, a%g, a%nw2, p, W(ih), W(ic), LF, H, C)
!       LT = LT + LF
!     else
!       call partial_eval(a%s, a%g, a%nw2, p, W(ih), W(ic), LT, H, C)
!     end if
!   end block
!
! contains
!
!   pure subroutine setminus_eval(a, p, ires, W, res)
!     type(d_matrix), intent(in) :: a
!     integer(IK), intent(in)    :: p, ires(*)
!     real(RK), intent(in)       :: W(*)
!     real(RK), intent(inout)    :: res
!     integer(IK)                :: nw, nr, rr
!
!     if (p < 0 .or. a%g < p) then
!       res = ZERO
!     elseif (p == 0) then
!       nw = nwork(a%g)
!       block
!         real(RK) :: T(nw)
!         call Hungarian(a%g, a%g, W(a%z), T(1))
!         res = T(1)
!       end block
!     else
!       nr = a%g - p
!       rr = nr * nr
!       nw = rr + nwork(nr)
!       block
!         integer(IK) :: iw, iz
!         real(RK)    :: T(nw)
!         iw = 1 + rr
!         iz = a%z + a%g * (p - 1)
!         call pack_Z(a%g, nr, ires, W(a%z), T(1))
!         call Hungarian(nr, nr, T(1), T(iw))
!         res = T(iw)
!       end block
!     end if
!
!   end subroutine setminus_eval
!
!   pure subroutine partial_eval(s, n, nw, p, H, C, LF, HP, CP)
!     integer(IK), intent(in) :: s, n, nw
!     integer(IK), intent(in) :: p
!     real(RK), intent(in)    :: H, C(*)
!     real(RK), intent(inout) :: LF, HP, CP(*)
!     real(RK)                :: W(nw)
!     integer(IK), parameter  :: it = 1
!     integer(IK), parameter  :: ic = 2
!     integer(IK)             :: iw
!
!     iw = ic + DD
!
!!! update H and C
!     w(it) = HP + H
!     HP = w(it)
!     call add(DD, CP, C, W(ic))
!     call DCOPY(DD, W(ic), 1, CP, 1)
!
!!! get squared displacement
!     call estimate_sdmin(w(it), w(ic), W(iw))
!     LF = LF + w(iw)
!
!   end subroutine partial_eval
!
!   pure elemental function nwork(np) result(res)
!     integer(IK), intent(in) :: np
!     real(RK)                :: W(1)
!     integer(IK)             :: res
!     call Hungarian(-np, -np, W, W)
!     res = NINT(W(1), IK)
!   end function nwork
!
!   pure subroutine pack_Z(n, nr, ires, Z, T)
!     integer(IK), intent(in) :: n, nr, ires(nr)
!     real(RK), intent(in)    :: Z(n, *)
!     real(RK), intent(inout) :: T(nr, nr)
!     integer(IK)             :: i, j
!     do concurrent(i=1:nr, j=1:nr)
!       T(i, j) = Z(ires(i), n - nr + j)
!     end do
!   end subroutine pack_Z
!
! end subroutine d_matrix_partial_eval
!
!!! d_matrix_list
!
!| Constructor of d_matrix_list
! pure function d_matrix_list_new(b, p) result(res)
!   type(mol_block_list), intent(in)     :: b
!   integer(IK), intent(in)              :: p
!   type(d_matrix_list)                  :: res
!   integer(IK)                          :: i, ip
!
!   res%l = b%nspecies()
!   res%h = p              ! H(1)
!   res%v = res%h + 1      ! V(1)
!   res%c = res%v + 1      ! C(d*d)
!   res%o = res%c + DD     ! O(L+1)
!   ip = res%o + res%l + 1
!
!   res%nk = worksize_sdmin()
!
!   allocate (res%m(res%l))
!
!   do i = 1, res%l
!     res%m(i) = d_matrix(ip, res%nk, b%b(i))
!     ip = ip + d_matrix_memsize(res%m(i))
!   end do
!
! end function d_matrix_list_new
!
! pure function d_matrix_list_memsize(this) result(res)
!   class(d_matrix_list), intent(in) :: this
!   integer(IK)                      :: res
!   if (ALLOCATED(this%m)) then
!     res = SUM(d_matrix_memsize(this%m)) + (this%l + 1) + DD + 2
!   else
!     res = 0
!   end if
! end function d_matrix_list_memsize
!
! pure elemental function d_matrix_list_n_depth(this) result(res)
!   class(d_matrix_list), intent(in) :: this
!   integer(IK)                      :: res
!   if (ALLOCATED(this%m)) then
!     res = SUM(this%m%g)
!   else
!     res = 0
!   end if
! end function d_matrix_list_n_depth
!
! pure subroutine d_matrix_list_eval(this, rot, X, Y, W)
!   class(d_matrix_list), intent(in) :: this
!   type(mol_symmetry), intent(in)   :: rot(*)
!   real(RK), intent(in)             :: X(*), Y(*)
!   real(RK), intent(inout)          :: W(*)
!   integer(IK)                      :: i
!
!   if (.not. ALLOCATED(this%m)) return
!
!!! estimate H_fix, V_fix, C_fix. (mol_blocks for g<n or (g=1, s=1))
!!! if set is empty, H_fix=0, V_fix=0, C_fix=0.
!
!   call fixpoints_eval(this%nk, this%l, this%m, X, Y,  &
!  &                    W(this%h), W(this%v), W(this%c))
!
!!! estimate floating mol blocks.
!!! O is filled by lower bounds for each blocks.
!
!   do concurrent(i=1:this%l)
!     call d_matrix_eval(this%m(i), rot(i), X, Y, W)
!     block
!       integer(IK) :: j, ires(this%m(i)%g)
!       do concurrent(j=1:this%m(i)%g)
!         ires(j) = j
!       end do
!       j = this%o + i - 1
!       call d_matrix_partial_eval(this%m(i), 0, 0, 0, ires, W, W(j), W(j), W(j))
!     end block
!   end do
!
!!! transform O to accumulated form.
!
!   do concurrent(i=this%l - 1:1:-1)
!     W(this%o + i - 1) = W(this%o + i - 1) + W(this%o + i)
!   end do
!
!   W(this%o + this%l) = ZERO
!
! contains
!
!   pure subroutine fixpoints_eval(nk, l, m, X, Y, H, V, C)
!     integer(IK), intent(in)    :: nk, l
!     type(d_matrix), intent(in) :: m(l)
!     real(RK), intent(in)       :: X(*), Y(*)
!     real(RK), intent(inout)    :: H, V, C(*)
!     integer(IK)                :: t(l), p(l), q(l)
!     integer(IK), parameter     :: ig = 1
!     integer(IK), parameter     :: iv = 2
!     integer(IK), parameter     :: ic = 3
!     integer(IK)                :: i, iw, ix, iy, mn, dmn, dmn2, nw
!
!     do concurrent(i=1:l)
!       t(i) = m(i)%m * (m(i)%n - m(i)%g)
!     end do
!
!     mn = SUM(t)
!     dmn = D * mn
!     dmn2 = dmn + dmn
!
!     iw = ic + DD
!     ix = ic + DD
!     iy = ix + dmn
!     nw = 2 + DD + MAX(dmn + dmn, DD + nk)
!
!     do concurrent(i=1:l)
!       t(i) = t(i) * D
!     end do
!     p(1) = 0
!     do i = 2, l
!       p(i) = p(i - 1) + t(i - 1)
!     end do
!     do concurrent(i=1:l)
!       q(i) = m(i)%x + m(i)%dm * m(i)%g
!     end do
!
!     block
!       real(RK) :: W(nw)
!
!       do concurrent(i=1:l)
!         block
!           integer(IK) :: px
!           px = p(i) + ix
!           call DCOPY(t(i), X(q(i)), 1, W(px), 1)
!         end block
!       end do
!       do concurrent(i=1:l)
!         block
!           integer(IK) :: py
!           py = p(i) + iy
!           call DCOPY(t(i), Y(q(i)), 1, W(py), 1)
!         end block
!       end do
!
!       W(ig) = ddot(dmn2, W(ix), 1, W(ix), 1)
!
!       if (mn > 0) then
!         call DGEMM('N', 'T', D, D, mn, ONE, W(iy), D, W(ix), D, ZERO, W(ic), D)
!         call estimate_sdmin(W(ig), W(ic), W(iw))
!         w(iv) = w(iw)
!       else
!         call zfill(DD, W(ic))
!         w(iv) = ZERO
!       end if
!
!       H = W(ig)
!       V = W(ig) - W(iv) - W(iv)
!       call DCOPY(DD, W(ic), 1, C, 1)
!
!     end block
!
!   end subroutine fixpoints_eval
!
! end subroutine d_matrix_list_eval
!
! pure subroutine d_matrix_list_partial_eval(this, p, perm, iprm, isym, W, LT, H, C, LF, LB)
!   class(d_matrix_list), intent(in)  :: this
!   integer(IK), intent(in)           :: p, perm(*), iprm, isym
!   real(RK), intent(in)              :: W(*)
!   real(RK), intent(inout)           :: LT, H, C(*)
!   real(RK), intent(inout), optional :: LF, LB
!   integer(IK)                       :: ispc, iofs, nres
!
!   call p_index(this, p, ispc, iofs)
!   if (ispc < 1) return
!
!   nres = this%m(ispc)%g - iofs
!   block
!     integer(IK) :: i, jper, ires(nres)
!     jper = perm(p - 1 + iprm)
!     do concurrent(i=1:iprm - 1)
!       ires(i) = perm(p + i - 1)
!     end do
!     do concurrent(i=iprm:nres)
!       ires(i) = perm(p + i)
!     end do
!     call d_matrix_partial_eval(this%m(ispc), iofs, jper, isym, ires, W, LT, H, C, LF, LB)
!   end block
!
! contains
!
!   pure elemental subroutine p_index(this, p, ispc, iofs)
!     type(d_matrix_list), intent(in) :: this
!     integer(IK), intent(in)         :: p
!     integer(IK), intent(inout)      :: ispc, iofs
!     integer(IK)                     :: i, q, r
!
!     ispc = 0
!     iofs = 0
!     r = p
!     q = 1
!
!     do i = 1, this%l
!       if (q > p) return
!       ispc = i
!       iofs = r
!       r = r - this%m(i)%g
!       q = q + this%m(i)%g
!     end do
!
!     if (q <= p) ispc = 0
!
!   end subroutine p_index
!
! end subroutine d_matrix_list_partial_eval
!
! pure elemental subroutine d_matrix_list_clear(this)
!   class(d_matrix_list), intent(inout) :: this
!   this%l = 0
!   if (ALLOCATED(this%m)) deallocate (this%m)
! end subroutine d_matrix_list_clear
!
! pure elemental subroutine d_matrix_list_destroy(this)
!   type(d_matrix_list), intent(inout) :: this
!   call d_matrix_list_clear(this)
! end subroutine d_matrix_list_destroy
!
!!! util
!
! pure subroutine add(d, A, B, C)
!   integer(IK), intent(in) :: d
!   real(RK), intent(in)    :: A(*), B(*)
!   real(RK), intent(inout) :: C(*)
!   integer(IK)             :: i
!   do concurrent(i=1:d)
!     C(i) = A(i) + B(i)
!   end do
! end subroutine add
!
! pure subroutine zfill(d, x)
!   integer(IK), intent(in) :: d
!   real(RK), intent(inout) :: x(*)
!   integer(IK)             :: i
!   do concurrent(i=1:d)
!     x(i) = ZERO
!   end do
! end subroutine zfill
!
end module mod_lb_manager

