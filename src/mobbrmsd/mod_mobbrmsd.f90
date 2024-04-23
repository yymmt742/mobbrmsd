!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
!$ use omp_lib
  use blas_lapack_interface, only: D, setup_dimension
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_bb_list
  use mod_bb_block
  use mod_mobbrmsd_header
  use mod_mobbrmsd_state
  use mod_forbar
  use mod_forbar_collections
  implicit none
  public :: setup_dimension
  public :: mobbrmsd_input
  public :: mol_block_input
  public :: mol_block_input_add
  public :: mobbrmsd
  public :: mobbrmsd_num_threads
  public :: mobbrmsd_header
  public :: mobbrmsd_state
  public :: mobbrmsd_run
  public :: mobbrmsd_restart
  public :: mobbrmsd_batch_run
  public :: mobbrmsd_nearest_neighbor
  public :: mobbrmsd_min_span_tree
  public :: mobbrmsd_is_finished
!
!| mobbrmsd_input
  type mobbrmsd_input
    private
    type(mol_block_input), allocatable :: blk(:)
  contains
    procedure :: add => mobbrmsd_input_add
    final     :: mobbrmsd_input_destroy
  end type mobbrmsd_input
!
!| mol_block_input (for python interface)
  type mol_block_input
    private
    integer(IK) :: m
    !! number of atoms per molecule
    integer(IK) :: n
    !! number of molecule
    integer(IK), allocatable :: sym(:, :)
    !! molecular symmetry
  contains
    final :: mol_block_input_destroy
  end type mol_block_input
!
!| mobbrmsd
  type mobbrmsd
    type(mobbrmsd_header) :: h
    !! mobbrmsd_header
    type(mobbrmsd_state)  :: s
    !! mobbrmsd_state
  contains
    final     :: mobbrmsd_destroy
  end type mobbrmsd
!
  interface mobbrmsd
    module procedure mobbrmsd_new, mobbrmsd_new_from_block
  end interface mobbrmsd
!
contains
!
!| add molecule
  pure subroutine mobbrmsd_input_add(this, m, n, sym)
    class(mobbrmsd_input), intent(inout) :: this
    !! mol_block_input array
    integer(IK), intent(in)              :: m
    !! number of atoms per molecule
    integer(IK), intent(in)              :: n
    !! number of molecule
    integer(IK), intent(in), optional    :: sym(:, :)
    !! molecular symmetry, sym(m, s-1)
    call mol_block_input_add(this%blk, m, n, sym)
  end subroutine mobbrmsd_input_add
!
!| destractor
  pure elemental subroutine mobbrmsd_input_destroy(this)
    type(mobbrmsd_input), intent(inout) :: this
    if (ALLOCATED(this%blk)) deallocate (this%blk)
  end subroutine mobbrmsd_input_destroy
!
! ------
!
!| add molecule
  pure subroutine mol_block_input_add(this, m, n, sym)
    type(mol_block_input), allocatable, intent(inout) :: this(:)
    !! mol_block_input array
    integer(IK), intent(in)            :: m
    !! number of atoms per molecule
    integer(IK), intent(in)            :: n
    !! number of molecule
    integer(IK), intent(in), optional  :: sym(:, :)
    !! molecular symmetry
    type(mol_block_input), allocatable :: blocks(:)
    integer(IK)                        :: nblock, i
!
    nblock = 1
    if (ALLOCATED(this)) nblock = nblock + SIZE(this)
!
    allocate (blocks(nblock))
    do concurrent(i=1:nblock - 1)
      blocks(i)%m = this(i)%m
      blocks(i)%n = this(i)%n
      call MOVE_ALLOC(from=this(i)%sym, to=blocks(i)%sym)
    end do
!
    blocks(nblock)%m = m
    blocks(nblock)%n = n
    if (PRESENT(sym)) then
      if (SIZE(sym) > 0) blocks(nblock)%sym = sym
    end if
!
    call MOVE_ALLOC(from=blocks, to=this)
!
  end subroutine mol_block_input_add
!
!| destractor
  pure elemental subroutine mol_block_input_destroy(this)
    type(mol_block_input), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
  end subroutine mol_block_input_destroy
!
! ------
!
!| constructor
  pure elemental function mobbrmsd_new(inp) result(res)
    type(mobbrmsd_input), intent(in) :: inp
    !! mobbrmsd_input
    type(mobbrmsd)                   :: res
!
    res = mobbrmsd_new_from_block(inp%blk)
!
  end function mobbrmsd_new
!
!| constructor, from mol_block_input. (for python interface)
  pure function mobbrmsd_new_from_block(blocks) result(res)
    type(mol_block_input), intent(in) :: blocks(:)
    !! mol_block_input array
    type(mobbrmsd)                    :: res
    integer(IK)                       :: nblock
!
    nblock = SIZE(blocks)
!
    block
      type(bb_block) :: bbblk(nblock)
      type(bb_list)  :: bblst
      integer(IK)    :: i
!
      do concurrent(i=1:nblock)
        if (ALLOCATED(blocks(i)%sym)) then
          bbblk(i) = bb_block(blocks(i)%m, blocks(i)%n, sym=blocks(i)%sym)
        else
          bbblk(i) = bb_block(blocks(i)%m, blocks(i)%n)
        end if
      end do
!
      bblst = bb_list(bbblk)
      res%h = mobbrmsd_header(bblst%q, bblst%s)
      res%s = mobbrmsd_state(res%h)

    end block
!
  end function mobbrmsd_new_from_block
!
!| returns omp_get_num_threads
  function mobbrmsd_num_threads() result(res)
    integer(IK) :: res
    !$omp parallel
    res = MAX(omp_get_num_threads(), 1)
    !$omp end parallel
  end function mobbrmsd_num_threads
!
!| run mobbrmsd
  subroutine mobbrmsd_run( &
 &             header, state, &
 &             X, Y, W, &
 &             cutoff, difflim, maxeval)
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(in)                 :: Y(*)
    !! target coordinate
    real(RK), intent(inout), optional    :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
!
    if (PRESENT(W)) then
!
      call bb_list_setup(header%q, &
     &                   state%s,  &
     &                   X,  &
     &                   Y,  &
     &                   W)
      call mobbrmsd_restart(header, &
     &                      state,  &
     &                      W, &
     &                      cutoff=cutoff, &
     &                      difflim=difflim, &
     &                      maxeval=maxeval)
!
    else
!
      block
        real(RK), allocatable :: T(:)
        allocate (T(header%memsize()))
        call bb_list_setup(header%q, &
       &                   state%s,  &
       &                   X,  &
       &                   Y,  &
       &                   T)
        call mobbrmsd_restart(header, &
       &                      state,  &
       &                      T, &
       &                      cutoff=cutoff, &
       &                      difflim=difflim, &
       &                      maxeval=maxeval)
      end block
!
    end if
!
  end subroutine mobbrmsd_run
!
!| run mobbrmsd
  subroutine mobbrmsd_restart(header, state, W, &
 &                                 cutoff, difflim, maxeval)
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(inout)              :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
!
    call bb_list_run(header%q, state%s, W, cutoff=cutoff, difflim=difflim, maxeval=maxeval)
    call state%update(header, W)
!
  end subroutine mobbrmsd_restart
!
  !| batch parallel run
  subroutine mobbrmsd_batch_run(n_target, header, state, &
  &                             X, Y, W, &
  &                             cutoff, difflim, maxeval, rotate_y)
    integer(IK), intent(in)              :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(inout)              :: Y(*)
    !! target coordinate
    real(RK), intent(inout)              :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional        :: rotate_y
    !! The search ends when ncount exceeds maxiter.
    integer(kind=IK)                     :: i, itgt, ijob, ypnt, wpnt, ldy, ldw
!
    i = 0
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    !$omp parallel private(itgt, ijob, ypnt, wpnt)
    do
      !$omp critical
      i = i + 1
      itgt = i
      !$omp end critical
      if (itgt > n_target) exit
      wpnt = ldw * omp_get_thread_num() + 1
      ypnt = (itgt - 1) * ldy + 1
      call mobbrmsd_run(header, state(itgt), X, Y(ypnt), W(wpnt), &
     &                  cutoff=cutoff, difflim=difflim, maxeval=maxeval)
      if (rotate_y) call state(itgt)%rotation(header, Y(ypnt))
    end do
    !$omp end parallel
!
  end subroutine mobbrmsd_batch_run
!
  !| batch parallel nearest_neighbor
  subroutine mobbrmsd_nearest_neighbor(n_target, header, state, X, Y, W, &
 &                                     cutoff, difflim, maxeval, &
 &                                     mask, nnval, nnidx)
    integer(IK), intent(in)             :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)           :: X(*)
   !! reference coordinate
    real(kind=RK), intent(in)           :: y(*)
   !! target coordinate
    real(kind=RK), intent(inout)        :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: mask(n_target)
    !! If .false., skip the calculation.
    real(RK), intent(out), optional     :: nnval
    !! nearest neighbor index
    integer(IK), intent(out), optional  :: nnidx
    !! nearest neighbor index
    real(kind=RK)                       :: cutoff_global, ub
    integer(kind=IK)                    :: i, itgt, ypnt, wpnt, ldy, ldw
!
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    if (PRESENT(cutoff)) then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    end if
!
    if (PRESENT(mask)) then
      cutoff_global = MIN(MINVAL(state%upperbound(), mask), cutoff_global)
    else
      cutoff_global = MIN(MINVAL(state%upperbound()), cutoff_global)
    end if
!
    i = 0
!
    !$omp parallel private(itgt, ub)
    do
      !$omp critical
      i = i + 1
      itgt = i
      ub = cutoff_global
      !$omp end critical
!
      if (itgt > n_target) exit
!
      if (PRESENT(mask)) then
        if (.not. mask(itgt)) cycle
      end if
      if (bb_list_is_finished(header%q, state(itgt)%s)) cycle
      if (ub < state(itgt)%lowerbound()) cycle
!
      wpnt = ldw * omp_get_thread_num() + 1
      ypnt = (itgt - 1) * ldy + 1
!
      call mobbrmsd_run(header, state(itgt), &
     &                  X, Y(ypnt), W(wpnt), &
     &                  cutoff=ub, &
     &                  difflim=difflim, &
     &                  maxeval=maxeval)
      ub = state(itgt)%upperbound()
!
      !$omp critical
      cutoff_global = MIN(ub, cutoff_global)
      !$omp end critical
    end do
    !$omp end parallel
!
    if (PRESENT(mask)) then
      if (PRESENT(nnval)) nnval = MINVAL(state%upperbound(), mask)
      if (PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1, mask)
    else
      if (PRESENT(nnval)) nnval = MINVAL(state%upperbound())
      if (PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1)
    end if
!
  end subroutine mobbrmsd_nearest_neighbor
!
!| min_span_tree construction
  subroutine mobbrmsd_min_span_tree(n_target, header, state, X, W, &
 &                                  cutoff, difflim, maxeval, &
 &                                  edges, weights, show_progress, &
 &                                  verbose &
 )
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target, n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)           :: X(*)
    !! coordinate sequence
    real(kind=RK), intent(inout)        :: W(*)
    !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    logical, intent(in), optional       :: show_progress
    !! if true, show progress bar
    logical, intent(in), optional       :: verbose
    !! show progress bar
    type(forbar)                        :: fbar
    logical                             :: verbose_
    logical                             :: mask(n_target)
    integer(kind=IK)                    :: list(2, n_target)
    real(RK)                            :: vval(n_target - 1), nnval, cutoff_
    integer(kind=IK)                    :: i, j, xpnt, ldx, nnidx
!
    ldx = header%n_dims() * header%n_atoms()
!
! Initialize header
    do concurrent(i=1:n_target, j=1:n_target)
      state(i, j) = mobbrmsd_state(header)
    end do
!
    mask(:) = .true.
    mask(1) = .false.
    list(1, 1) = 0
    list(2, 1) = 1
!
    if (PRESENT(verbose)) then
      verbose_ = verbose
    else
      verbose_ = .false.
    end if
!
    if (verbose_) then
      fbar = progress_bar(limit=n_target - 1)
    end if
!
    do j = 1, n_target - 1
!
      if (verbose_) call fbar%update_and_display()
!
      vval(j) = RHUGE
      if (PRESENT(cutoff)) then
        cutoff_ = MIN(RHUGE, MAX(cutoff, ZERO))
      else
        cutoff_ = RHUGE
      end if
!
      do i = 1, j
!
        xpnt = (list(2, i) - 1) * ldx + 1
!
        call mobbrmsd_nearest_neighbor( &
       &  n_target, header, &
       &  state(:, list(2, i)), &
       &  X(xpnt), X, W, &
       &  cutoff=cutoff_, &
       &  difflim=difflim, &
       &  maxeval=maxeval, &
       &  mask=mask, &
       &  nnval=nnval,&
       &  nnidx=nnidx)
!
        if (nnval < vval(j)) then
          vval(j) = nnval
          cutoff_ = MIN(cutoff_, vval(j))
          list(1, j + 1) = list(2, i)
          list(2, j + 1) = nnidx
        end if
!
      end do
!
      mask(list(2, j + 1)) = .false.
!
    end do
!
    if (verbose_) call fbar%clean_up()
!
    do j = 1, n_target
      do i = 1, j - 1
        if (state(i, j)%upperbound() < state(j, i)%upperbound()) then
          state(j, i) = state(i, j)
        else
          state(i, j) = state(j, i)
        end if
      end do
    end do
!
    if (PRESENT(edges)) then
      do concurrent(i=1:n_target - 1)
        edges(1, i) = list(1, i + 1)
        edges(2, i) = list(2, i + 1)
      end do
    end if
!
    if (PRESENT(weights)) then
      do concurrent(i=1:n_target - 1)
        weights(i) = vval(i)
      end do
    end if
!
  end subroutine mobbrmsd_min_span_tree
!
!| Returns bb process is finished.
  pure function mobbrmsd_is_finished(header, state) result(res)
    class(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd_header
    class(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    logical                            :: res
!
    res = bb_list_is_finished(header%q, state%s)
!
  end function mobbrmsd_is_finished
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

