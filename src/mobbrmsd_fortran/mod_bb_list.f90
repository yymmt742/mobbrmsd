!| Manage multiple bb_blocks and execute multi-component Branch-and-Bound.
module mod_bb_list
  use mod_params, only: IK, RK, ONE => RONE, TWO => RTWO, HALF => RHALF, ZERO => RZERO, RHUGE
  use mod_dimspec_functions, only: D, DD, compute_com
  use mod_bb_block
  use mod_rotation
  implicit none
  private
  public :: bb_list
  public :: bb_list_init
  public :: bb_list_memsize
  public :: bb_list_n_block
  public :: bb_list_n_atoms
  public :: bb_list_log_n_nodes
  public :: bb_list_setup
  public :: bb_list_run
  public :: bb_list_swap_y
  public :: bb_list_swap_indices
  public :: bb_list_rotation_matrix
  public :: bb_list_is_finished
  public :: bb_list_destroy
  public :: bb_list_INDEX_TO_AUTOCORR
  public :: bb_list_INDEX_TO_UPPERBOUND
  public :: bb_list_INDEX_TO_LOWERBOUND
  public :: bb_list_INDEX_TO_N_EVAL
  public :: bb_list_INDEX_TO_LOG_N_COMB
!
  integer(IK), parameter :: header_size = 1
  integer(IK), parameter :: bb_list_NUMBER_OF_SPEACIES = 1 ! number of block
!
  integer(IK), parameter :: header_sttsize = 1
  integer(IK), parameter :: bb_list_INDEX_TO_SPEACIES = 1
  integer(IK), parameter :: bb_list_INDEX_TO_BESTSTATE = 2 ! pointer to best state vector
!
  integer(IK), parameter :: bb_list_INDEX_TO_AUTOCORR = 1
  integer(IK), parameter :: bb_list_INDEX_TO_UPPERBOUND = 2
  integer(IK), parameter :: bb_list_INDEX_TO_LOWERBOUND = 3
  integer(IK), parameter :: bb_list_INDEX_TO_N_EVAL = 4
  integer(IK), parameter :: bb_list_INDEX_TO_LOG_N_COMB = 5
  integer(IK), parameter :: header_memsize = 5
!| This derived type is mainly used for passing during initialization.
  type bb_list
    sequence
    integer(IK), allocatable :: q(:)
    !! integer array
    integer(IK), allocatable :: s(:)
    !! work integer array
  end type bb_list
!
!| Constructer
  interface bb_list
    module procedure bb_list_new
  end interface bb_list
!
contains
!| Constructer
  pure function bb_list_new(n, blk) result(res)
    integer, intent(in)        :: n
    type(bb_block), intent(in) :: blk(n)
    type(bb_list)              :: res
    call bb_list_init(res, n, blk)
  end function bb_list_new
!
!| Constructer
  pure subroutine bb_list_init(this, n, blk)
    type(bb_list), intent(inout) :: this
    integer, intent(in)          :: n
    type(bb_block), intent(in)   :: blk(n)
    integer(IK)                  :: q(header_size), s(header_sttsize)
    integer(IK)                  :: pq(n), pa(n), px(n), ps(n), pw(n)
    integer(IK)                  :: i, j, nstat
    associate ( &
   &  nb => q(bb_list_NUMBER_OF_SPEACIES), &
   &  sb => s(bb_list_INDEX_TO_SPEACIES) &
   &  )
      nb = SIZE(blk)
      sb = 0
!
      j = header_size + 5 * nb + 1
      do i = 1, nb
        pq(i) = j
        j = j + SIZE(blk(i)%q)
      end do
      nstat = SUM([(bb_block_statesize(blk(i)%q), i=1, nb)])
      j = header_sttsize + nstat + 1
      do i = 1, nb
        ps(i) = j
        j = j + SIZE(blk(i)%s)
      end do
      j = header_memsize + 1
      do i = 1, nb
        pw(i) = j
        j = j + bb_block_memsize(blk(i)%q) + bb_block_worksize(blk(i)%q)
      end do
      j = 1
      do i = 1, nb
        px(i) = j
        j = j + bb_block_molsize(blk(i)%q)
      end do
      j = 1
      do i = 1, nb
        pa(i) = j
        j = j + bb_block_natm(blk(i)%q)
      end do
!
      allocate (this%q, source=[q, pq, ps, pw, px, pa, [(blk(i)%q, i=1, SIZE(blk))]])
      allocate (this%s, source=[s, [(-1, i=1, nstat)], [(blk(i)%s, i=1, SIZE(blk))]])
    end associate
  end subroutine bb_list_init
!
!| Inquire worksize of f_matrix.
  pure function bb_list_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res, i, j
    associate (n => q(bb_list_NUMBER_OF_SPEACIES))
      res = header_memsize
      j = q_pointer(q)
      do i = 1, n
        res = res + bb_block_memsize(q(q(j))) + bb_block_worksize(q(q(j)))
        j = j + 1
      end do
    end associate
  end function bb_list_memsize
!
!| Setup
  pure subroutine bb_list_setup(q, s, X, Y, W, remove_com, sort_by_g)
    integer(IK), intent(in)       :: q(*)
    !! header
    integer(IK), intent(inout)    :: s(*)
    !! state
    real(RK), intent(in)          :: X(*)
    !! reference coordinate
    real(RK), intent(in)          :: Y(*)
    !! target coordinate
    real(RK), intent(inout)       :: W(*)
    !! work array
    logical, intent(in), optional :: remove_com
    !! if true, remove centroids.
    !  default [.true.]
    logical, intent(in), optional :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate.
    !  default [.true.]
    real(RK)                      :: CX(D), CY(D)
    !! centroids
    integer(IK)                   :: i, ps, pq, px, pw, n_atoms
    associate ( &
       n_block => q(bb_list_NUMBER_OF_SPEACIES), &
   &   sb => s(bb_list_INDEX_TO_SPEACIES), &
   &   ac => W(bb_list_INDEX_TO_AUTOCORR), &
   &   ub => W(bb_list_INDEX_TO_UPPERBOUND), &
   &   lb => W(bb_list_INDEX_TO_LOWERBOUND), &
   &   nv => W(bb_list_INDEX_TO_N_EVAL), &
   &   cm => W(bb_list_INDEX_TO_LOG_N_COMB) &
   &  )
      sb = 0
      ub = ZERO
      !lb = -RHUGE
      nv = ZERO
      ps = s_pointer(q)
      px = x_pointer(q)
      pw = w_pointer(q)
      pq = q_pointer(q)
!
      n_atoms = bb_list_n_atoms(q)
!
      if (PRESENT(remove_com)) then
        if (remove_com) then
          call compute_com(D, n_atoms, X, CX)
          call compute_com(D, n_atoms, Y, CY)
        else
          CX = ZERO; CY = ZERO
        end if
      else
        call compute_com(D, n_atoms, X, CX)
        call compute_com(D, n_atoms, Y, CY)
      end if
!
      do concurrent(i=0:n_block - 1)
        call bb_block_setup( &
       &  q(q(pq + i)), &
       &  X(q(px + i)), &
       &  Y(q(px + i)), &
       &  CX, CY, &
       &  s(q(ps + i)), W(q(pw + i)), &
       &  zfill=(i == 0),&
       &  sort_by_g=sort_by_g &
       & )
      end do
!
      ac = ZERO
      do i = n_block - 1, 0, -1
        call bb_block_set_ub_offset(W(q(pw + i)), ac)
        ac = ac + bb_block_lowerbound(W(q(pw + i)))
      end do
      ac = ZERO
      do i = 0, n_block - 1
        ac = ac + bb_block_autocorr(q(q(pq + i)), W(q(pw + i)))
      end do
!
      lb = -HALF * ac
!
      cm = ZERO
      do i = 0, n_block - 1
        cm = cm + bb_block_log_ncomb(q(q(pq + i)))
      end do
      call save_state(q(pq), q(ps), q, s)
    end associate
  end subroutine bb_list_setup
!
!| run branch and bound
  pure subroutine bb_list_run(q, s, W, cutoff, ub_cutoff, difflim, maxeval, difflim_absolute)
    integer(IK), intent(in)           :: q(*)
    !! header
    integer(IK), intent(inout)        :: s(*)
    !! state
    real(RK), intent(inout)           :: W(*)
    !! work array
    real(RK), intent(in), optional    :: cutoff
    !! The search ends when lowerbound is determined
    !  to be greater than to cutoff (in RMSD).
    real(RK), intent(in), optional    :: ub_cutoff
    !! The search ends when lowerbound is determined
    !  to be greater than to cutoff (in RMSD).
    real(RK), intent(in), optional    :: difflim
    !! The search ends when the difference
    !  between the lower and upper gap is less than
    !  \(\text{threshold} = (G / 2) * \text{difflim}^2\).
    !  default = 0.0.
    integer(IK), intent(in), optional :: maxeval
    !! The search ends when ncount exceeds maxiter.
    !  If maxeval=0, run only once, and early return.
    logical, intent(in), optional     :: difflim_absolute
    !! If True, use absolute difflim value for diff threshold.
    !  Since the gap between the upperbound and lowerbound corresponds to SD,
    !  difflim is converted to the RMSD scale, such that
    !  \(\text{threshold} = \text{difflim}^2 * n_atoms\).
    !  default False.
    real(RK)                          :: coff, ucut, diff, nlim
    integer(IK)                       :: pq, ps, pw
    associate ( &
   &   b => s(bb_list_INDEX_TO_SPEACIES), &
   &   ac => W(bb_list_INDEX_TO_AUTOCORR), &
   &   ub => W(bb_list_INDEX_TO_UPPERBOUND), &
   &   lb => W(bb_list_INDEX_TO_LOWERBOUND), &
   &   nv => W(bb_list_INDEX_TO_N_EVAL) &
   &  )
      b = MAX(b, 1)
      pq = q_pointer(q)
      ps = s_pointer(q)
      pw = w_pointer(q)
!
      if (PRESENT(cutoff)) then
        coff = HALF * (cutoff**2 * bb_list_n_atoms(q) - ac)
      else
        coff = RHUGE
      end if
!
      if (PRESENT(ub_cutoff)) then
        ucut = HALF * (ub_cutoff**2 * bb_list_n_atoms(q) - ac)
      else
        ucut = RHUGE
      end if
!
      if (PRESENT(difflim)) then
        if (PRESENT(difflim_absolute)) then
          if (difflim_absolute) then
            diff = MAX(ZERO, difflim**2 * bb_list_n_atoms(q)) ! use absolute delta (as RMSD).
          else
            diff = MAX(ZERO, HALF * ac * difflim) ! use relative delta.
          end if
        else
          diff = MAX(ZERO, HALF * ac * difflim) ! use relative delta.
        end if
      else
        diff = ZERO
      end if
!
      if (PRESENT(maxeval)) then
        if (maxeval == 0) then
          call run_bb(q(pq), q(ps), q(pw), q, s, W) ! run only once, and early return.
          return
        elseif (maxeval > 0) then
          nlim = real(maxeval, RK)  ! finite run
        else
          nlim = RHUGE ! unlimited run if maxeval < 0
        end if
      else
        nlim = RHUGE ! unlimited run
      end if
!
      call update_lowerbound(b, q(pq), q(ps), q(pw), q, s, W)
!
      do while (nv < nlim &
          .and. lb < coff &
     &    .and. lb + diff <= ub)
        call run_bb(q(pq), q(ps), q(pw), q, s, W)
        if (bb_list_is_finished(q, s) .or. ub >= ucut) exit
      end do
    end associate
  end subroutine bb_list_run
!
  pure subroutine run_bb(pq, ps, pw, q, s, W)
    integer(IK), intent(in)    :: pq(*), ps(*), pw(*), q(*)
    integer(IK), intent(inout) :: s(*)
    real(RK), intent(inout)    :: W(*)
    associate ( &
   &   n => q(bb_list_NUMBER_OF_SPEACIES), &
   &   b => s(bb_list_INDEX_TO_SPEACIES), &
   &   ub => W(bb_list_INDEX_TO_UPPERBOUND), &
   &   lb => W(bb_list_INDEX_TO_LOWERBOUND), &
   &   nv => W(bb_list_INDEX_TO_N_EVAL) &
   &  )
!
!     Expansion process
!
      do
        call bb_block_descend(ub, q(pq(b)), s(ps(b)), W(pw(b)))
        if (b == n .or. bb_block_tree_is_empty(q(pq(b)), s(ps(b)))) exit
        b = b + 1
        call bb_block_inheritance(q(pq(b)), s(ps(b)), W(pw(b)), &
       &                          q(pq(b - 1)), s(ps(b - 1)), W(pw(b - 1)))
      end do
!
!     Update upperbound and state
!
      if (bb_block_is_bottom(q(pq(b)), s(ps(b))) .and. b == n) then
        block
          real(RK) :: cv
          cv = bb_block_current_value(q(pq(b)), s(ps(b)), W(pw(b)))
          if (ub > cv) then
            ub = cv
            call save_state(pq, ps, q, s)
          end if
        end block
      end if
!
!     Update lowerbound
!
      call update_lowerbound(b, pq, ps, pw, q, s, W)
!
!     ascend process
!
      do
        call bb_block_ascend(ub, q(pq(b)), s(ps(b)), W(pw(b)))
        if (b == 1 .or. bb_block_is_left(ub, q(pq(b)), s(ps(b)), W(pw(b)))) exit
        b = b - 1
      end do
!
      block
        integer(IK) :: i
        nv = ZERO
        do i = 1, n
          nv = nv + bb_block_evaluation_count(W(pw(i)))
        end do
      end block
    end associate
  end subroutine run_bb
!
  pure subroutine update_lowerbound(n, pq, ps, pw, q, s, W)
    integer(IK), intent(in)    :: n, pq(n), ps(n), pw(n)
    integer(IK), intent(in)    :: q(*), s(*)
    real(RK), intent(inout)    :: W(*)
    real(RK)                   :: lv
    integer(IK)                :: b
    associate ( &
   &  ub => W(bb_list_INDEX_TO_UPPERBOUND), &
   &  lb => W(bb_list_INDEX_TO_LOWERBOUND) &
   &  )
      lv = RHUGE
      do b = 1, n
        lv = MIN(lv, bb_block_lowest_value(q(pq(b)), s(ps(b)), W(pw(b))))
      end do
      lb = MIN(MAX(lb, lv), ub)
    end associate
  end subroutine update_lowerbound
!
  pure subroutine save_state(pq, ps, q, s)
    integer(IK), intent(in)    :: pq(*), ps(*)
    integer(IK), intent(in)    :: q(*)
    integer(IK), intent(inout) :: s(*)
    integer(IK)                :: i, j
    associate (n => q(bb_list_NUMBER_OF_SPEACIES))
      j = bb_list_INDEX_TO_BESTSTATE
      do i = 1, n
        call bb_block_save_state(q(pq(i)), s(ps(i)), s(j))
        j = j + bb_block_statesize(q(pq(i)))
      end do
    end associate
  end subroutine save_state
!
!| Swap target coordinate.
  pure subroutine bb_list_swap_y(q, s, Y)
    integer(IK), intent(in) :: q(*)
    !! header
    integer(IK), intent(in) :: s(*)
    !! state
    real(RK), intent(inout) :: Y(*)
    !! target coordinate
    integer(IK)             :: i, pb, pq, ps, px
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      px = x_pointer(q)
      ps = s_pointer(q)
      pq = q_pointer(q)
      pb = bb_list_INDEX_TO_BESTSTATE
!
      do i = 0, n_block - 1
        call bb_block_swap_y(q(q(pq + i)), s(q(ps + i)), s(pb), Y(q(px + i)))
        pb = pb + bb_block_statesize(q(q(pq + i)))
      end do
    end associate
  end subroutine bb_list_swap_y
!
!| Swap target coordinate.
  pure subroutine bb_list_swap_indices(q, s, IX)
    integer(IK), intent(in)    :: q(*)
    !! header
    integer(IK), intent(in)    :: s(*)
    !! state
    integer(IK), intent(inout) :: IX(*)
    !! swap indices
    integer(IK)             :: i, pb, pq, ps, pa
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      pa = a_pointer(q)
      ps = s_pointer(q)
      pq = q_pointer(q)
      pb = bb_list_INDEX_TO_BESTSTATE
!
      do concurrent(i=1:bb_list_n_atoms(q))
        IX(i) = i
      end do
      do i = 0, n_block - 1
        call bb_block_swap_indices(q(q(pq + i)), s(q(ps + i)), s(pb), IX(q(pa + i)))
        pb = pb + bb_block_statesize(q(q(pq + i)))
      end do
    end associate
  end subroutine bb_list_swap_indices
!
!| Sum covariance matrix by saved state z.
  pure subroutine bb_list_rotation_matrix(q, s, W, R)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! state
    real(RK), intent(in)    :: W(*)
    !! main memory
    real(RK), intent(inout) :: R(*)
    !! rotation matrix
    real(RK)                :: G, C(DD), V(rotation_worksize())
    integer(IK)             :: i, pb, pq, pw
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      pb = bb_list_INDEX_TO_BESTSTATE
      pw = w_pointer(q)
      pq = q_pointer(q)
!
      G = ZERO
      C = ZERO
!
      do i = 0, n_block - 1
        call bb_block_covmat_add(q(q(pq + i)), s(pb), W(q(pw + i)), G, C)
        pb = pb + bb_block_nmol(q(q(pq + i)))
      end do
!
      call estimate_rotation(G, C, R, V)
    end associate
  end subroutine bb_list_rotation_matrix
!
!| Returns number of molecular blocks.
  pure function bb_list_n_block(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      res = n_block
    end associate
  end function bb_list_n_block
!
!| Returns number of total atoms.
  pure function bb_list_n_atoms(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: res, i, j
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      res = 0
      j = q_pointer(q)
      do i = 1, n_block
        res = res + bb_block_natm(q(q(j)))
        j = j + 1
      end do
    end associate
  end function bb_list_n_atoms
!
!| Returns the logarithm of the total number of nodes.
  pure function bb_list_log_n_nodes(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block.
    integer(IK)             :: i, j
    real(RK)                :: res
    associate (n_block => q(bb_list_NUMBER_OF_SPEACIES))
      res = ZERO
      j = q_pointer(q)
      do i = 1, n_block
        res = res + bb_block_log_ncomb(q(q(j)))
        j = j + 1
      end do
    end associate
  end function bb_list_log_n_nodes
!
!| Returns bb is finished.
  pure function bb_list_is_finished(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! header
    integer(IK), intent(in) :: s(*)
    !! state
    logical                 :: res
    integer(IK)             :: bq, bs
    associate (sb => s(bb_list_INDEX_TO_SPEACIES))
!     early return
      res = sb == 1; if (.not. res) return
      bq = q(q_pointer(q) + sb - 1)
      bs = q(s_pointer(q) + sb - 1)
      res = bb_block_tree_is_finished(q(bq), s(bs))
    end associate
  end function bb_list_is_finished
!
!| destractor
  pure elemental subroutine bb_list_destroy(this)
    type(bb_list), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine bb_list_destroy
!
!  ---
!
  pure function s_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + q(bb_list_NUMBER_OF_SPEACIES) + 1
  end function s_pointer
!
  pure function w_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 2 * q(bb_list_NUMBER_OF_SPEACIES) + 1
  end function w_pointer
!
  pure function x_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 3 * q(bb_list_NUMBER_OF_SPEACIES) + 1
  end function x_pointer
!
  pure function a_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 4 * q(bb_list_NUMBER_OF_SPEACIES) + 1
  end function a_pointer
!
  pure function q_pointer(q) result(res)
    integer(IK), intent(in) :: q(*)
    integer(IK)             :: res
    res = header_size + 1
  end function q_pointer
!
end module mod_bb_list

