module driver
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mobbrmsd_header, only: &
    mobbrmsd_header
  use mod_mobbrmsd_state, only: &
 &  mobbrmsd_state, &
 &  RN => mobbrmsd_state_INDEX_TO_RCP_N_ATOMS, &
 &  UB => mobbrmsd_state_INDEX_TO_UPPERBOUND, &
 &  LB => mobbrmsd_state_INDEX_TO_LOWERBOUND, &
 &  LR => mobbrmsd_state_INDEX_TO_LOG_RATIO, &
 &  NE => mobbrmsd_state_INDEX_TO_N_EVAL, &
 &  RT => mobbrmsd_state_INDEX_TO_ROTMAT
  use mod_mobbrmsd, only: &
 &  mobbrmsd, &
 &  mobbrmsd_num_threads, &
 &  mobbrmsd_run, &
 &  mobbrmsd_restart, &
 &  mobbrmsd_batch_run, &
 &  mobbrmsd_min_span_tree, &
 &  mobbrmsd_is_finished, &
 &  mol_block_input, &
 &  mol_block_input_add, &
 &  setup_dimension_ => setup_dimension

  implicit none
  private
  public setup_dimension
  public add_molecule
  public clear_molecule
  public n_atoms
  public workmemory_lengthes
  public state_vector_lengthes
  public rmsd
  public bounds
  public n_eval
  public log_eval_ratio
  public is_finished
  public run
  public restart
  public rotate_y
  public batch_run
  public min_span_tree

  type(mol_block_input), allocatable :: blocks(:)

contains

  !| setup dimension
  subroutine setup_dimension(d)
    integer(kind=IK), intent(in) :: d

    call setup_dimension_(d)

  end subroutine setup_dimension

  !| add molecule
  subroutine add_molecule(n, M, s, sym)
    integer(kind=IK), intent(in) :: n
    integer(kind=IK), intent(in) :: M
    integer(kind=IK), intent(in) :: s
    integer(kind=IK), intent(in), optional :: sym(n * (s - 1))

    if (.not. ALLOCATED(blocks)) allocate (blocks(0))
    call mol_block_input_add(blocks, n, M, RESHAPE(sym, [n, s - 1]))

  end subroutine add_molecule

  !| clear molecule
  subroutine clear_molecule()

    if (ALLOCATED(blocks)) deallocate (blocks)

  end subroutine clear_molecule

  !| Returns spatial dimension and total number of atoms
  pure subroutine n_atoms(n_dim, n_atom)
    integer(kind=IK), intent(out) :: n_dim
    integer(kind=IK), intent(out) :: n_atom
    type(mobbrmsd)                :: mob

    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      n_dim = mob%h%n_dims()
      n_atom = mob%h%n_atoms()
    else
      n_dim = 0
      n_atom = 0
    end if

  end subroutine n_atoms

  !| Returns total number of atoms
  subroutine workmemory_lengthes(n_mem, n_job)
    integer(kind=IK), intent(out) :: n_mem
    integer(kind=IK), intent(out) :: n_job
    type(mobbrmsd)                :: mob

    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      n_mem = mob%h%memsize()
    else
      n_mem = 0
    end if

    n_job = mobbrmsd_num_threads()

  end subroutine workmemory_lengthes

  !| return header and state vector lengthes
  pure subroutine state_vector_lengthes(n_head, n_int, n_float)
    integer(kind=IK), intent(out) :: n_head
    integer(kind=IK), intent(out) :: n_int
    integer(kind=IK), intent(out) :: n_float
    type(mobbrmsd)                :: mob

    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      n_head = SIZE(mob%h%dump())
      n_int = SIZE(mob%s%dump())
      n_float = SIZE(mob%s%dump_real())
    else
      n_head = 0
      n_int = 0
      n_float = 0
    end if

  end subroutine state_vector_lengthes

  !| return rmsd
  subroutine rmsd(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res
    res = SQRT(float_states(RN) * MAX(float_states(UB), ZERO))
  end subroutine rmsd

  !| return bounds
  subroutine bounds(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res(2)
    res(1) = float_states(UB)
    res(2) = float_states(LB)
  end subroutine bounds

  !| return n_eval
  subroutine n_eval(n_float, float_states, res)
    integer(kind=IK), intent(in)  :: n_float
    real(kind=RK), intent(in)     :: float_states(n_float)
    integer(kind=IK), intent(out) :: res
    res = NINT(float_states(NE), IK)
  end subroutine n_eval

  !| return log_eval_ratio
  subroutine log_eval_ratio(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res
    res = float_states(LR)
  end subroutine log_eval_ratio

  !| inquire bb is finished
  subroutine is_finished(n_head, n_int, n_float, header, int_states, float_states, res)
    integer(kind=IK), intent(in)  :: n_head, n_int, n_float
    integer(kind=IK), intent(in)  :: header(n_head)
    integer(kind=IK), intent(in)  :: int_states(n_int)
    real(kind=RK), intent(in)     :: float_states(n_float)
    logical, intent(out)          :: res
    type(mobbrmsd_header)         :: h
    type(mobbrmsd_state)          :: s

    call h%load(header)
    call s%load(int_states, float_states)
    res = mobbrmsd_is_finished(h, s)

  end subroutine is_finished

  !| restart with working memory
  subroutine restart(n_head, n_int, n_float, n_mem, &
 &                   header, int_states, float_states, w, &
 &                   cutoff, difflim, maxeval)
    integer(kind=IK), intent(in)    :: n_head
    integer(kind=IK), intent(in)    :: n_int
    integer(kind=IK), intent(in)    :: n_float
    integer(kind=IK), intent(in)    :: n_mem
    integer(kind=IK), intent(in)    :: header(n_head)
    integer(kind=IK), intent(inout) :: int_states(n_int)
    real(kind=RK), intent(inout)    :: float_states(n_float)
    real(kind=RK), intent(inout)    :: W(n_mem)
   !! work memory
    integer(kind=IK), intent(in)    :: maxeval
    real(kind=RK), intent(in)       :: cutoff
    real(kind=RK), intent(in)       :: difflim
    type(mobbrmsd_header)           :: h
    type(mobbrmsd_state)            :: s

    call h%load(header)
    call s%load(int_states, float_states)
    call mobbrmsd_restart(h, s, W, cutoff, difflim, maxeval)

    int_states = s%dump()
    float_states = s%dump_real()

  end subroutine restart

  !| single run with working memory
  subroutine rotate_y(n_dim, n_atom, n_head, n_int, n_float,&
 &                    header, int_states, float_states, Y)
    integer(kind=IK), intent(in) :: n_dim
    integer(kind=IK), intent(in) :: n_atom
    integer(kind=IK), intent(in) :: n_head
    integer(kind=IK), intent(in) :: n_int
    integer(kind=IK), intent(in) :: n_float
    integer(kind=IK), intent(in) :: header(n_head)
    integer(kind=IK), intent(in) :: int_states(n_int)
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(inout) :: Y(n_dim, n_atom)
    type(mobbrmsd_header)        :: h
    type(mobbrmsd_state)         :: s

    call h%load(header)
    call s%load(int_states, float_states)
    call s%rotation(h, Y)

  end subroutine rotate_y

  !| single run with working memory
  subroutine run(n_dim, n_atom, n_head, n_int, n_float, n_mem, X, Y, W, &
 &               cutoff, difflim, maxeval, rotate_y, &
 &               header, int_states, float_states)
    integer(kind=IK), intent(in)  :: n_dim
    integer(kind=IK), intent(in)  :: n_atom
    integer(kind=IK), intent(in)  :: n_head
    integer(kind=IK), intent(in)  :: n_int
    integer(kind=IK), intent(in)  :: n_float
    integer(kind=IK), intent(in)  :: n_mem
    real(kind=RK), intent(in)     :: X(n_dim, n_atom)
   !! reference coordinate
    real(kind=RK), intent(inout)  :: Y(n_dim, n_atom)
   !! target coordinate
    real(kind=RK), intent(inout)  :: W(n_mem)
   !! work memory
    integer(kind=IK), intent(in)  :: maxeval
    real(kind=RK), intent(in)     :: cutoff
    real(kind=RK), intent(in)     :: difflim
    logical, intent(in)           :: rotate_y
    integer(kind=IK), intent(out) :: header(n_head)
    integer(kind=IK), intent(out) :: int_states(n_int)
    real(kind=RK), intent(out)    :: float_states(n_float)
    type(mobbrmsd)                :: mob

    mob = mobbrmsd(blocks)
    call mobbrmsd_run(mob%h, mob%s, X, Y, w, cutoff, difflim, maxeval)
    if (rotate_y) call mob%s%rotation(mob%h, Y)

    header = mob%h%dump()
    int_states = mob%s%dump()
    float_states = mob%s%dump_real()

  end subroutine run

  !| batch parallel run
  subroutine batch_run( &
               n_dim, n_atom, n_target, &
 &             n_head, n_int, n_float, n_mem, n_job,&
 &             X, Y, W, &
 &             cutoff, difflim, maxeval, rotate_y, &
 &             header, int_states, float_states)
    integer(kind=IK), intent(in)  :: n_dim
    integer(kind=IK), intent(in)  :: n_atom
    integer(kind=IK), intent(in)  :: n_target
    integer(kind=IK), intent(in)  :: n_head
    integer(kind=IK), intent(in)  :: n_int
    integer(kind=IK), intent(in)  :: n_float
    integer(kind=IK), intent(in)  :: n_mem
    integer(kind=IK), intent(in)  :: n_job
    real(kind=RK), intent(in)     :: X(n_dim, n_atom)
   !! reference coordinate
    real(kind=RK), intent(inout)  :: Y(n_dim, n_atom, n_target)
   !! target coordinate
    real(kind=RK), intent(inout)  :: W(n_mem, n_job)
   !! work array
    integer(kind=IK), intent(in)  :: maxeval
    real(kind=RK), intent(in)     :: cutoff
    real(kind=RK), intent(in)     :: difflim
    logical, intent(in)           :: rotate_y
    integer(kind=IK), intent(out) :: header(n_head)
    integer(kind=IK), intent(out) :: int_states(n_int, n_target)
    real(kind=RK), intent(out)    :: float_states(n_float, n_target)
    type(mobbrmsd)                :: mob
    type(mobbrmsd_state)          :: s(n_target)
    integer(kind=IK)              :: i

    mob = mobbrmsd(blocks)
    do concurrent(i=1:n_target)
      s(i) = mob%s
    end do

    call mobbrmsd_batch_run( &
   &  n_target, mob%h, s, X, Y, W, &
   &  cutoff, difflim, maxeval, rotate_y &
   &  )

    header = mob%h%dump()
    do concurrent(i=1:n_target)
      int_states(:, i) = s(i)%dump()
      float_states(:, i) = s(i)%dump_real()
    end do

  end subroutine batch_run

  subroutine min_span_tree(n_dim, n_atom, n_target, n_head, &
 &                         n_int, n_float, n_mem, n_job,&
 &                         x, w, cutoff, difflim,  &
 &                         maxeval, verbose, &
 &                         edges, weights, header, &
 &                         int_states, float_states)
    integer(kind=IK), intent(in)      :: n_dim
    integer(kind=IK), intent(in)      :: n_atom
    integer(kind=IK), intent(in)      :: n_target
    integer(kind=IK), intent(in)      :: n_head
    integer(kind=IK), intent(in)      :: n_int
    integer(kind=IK), intent(in)      :: n_float
    integer(kind=IK), intent(in)      :: n_mem
    integer(kind=IK), intent(in)      :: n_job
    real(kind=RK), intent(in)         :: x(n_dim, n_atom, n_target)
   !! reference coordinate
    real(kind=RK), intent(inout)      :: W(n_mem, n_job)
   !! work memory
    real(kind=RK), intent(in)         :: cutoff
    real(kind=RK), intent(in)         :: difflim
    integer(kind=IK), intent(in)      :: maxeval
    logical, intent(in)               :: verbose
    integer(kind=IK), intent(out)     :: edges(2, n_target - 1)
    real(kind=RK), intent(out)        :: weights(n_target - 1)
    integer(kind=IK), intent(out)     :: header(n_head)
    integer(kind=IK), intent(out)     :: int_states(n_int, n_target, n_target)
    real(kind=RK), intent(out)        :: float_states(n_float, n_target, n_target)
    type(mobbrmsd)                    :: mob
    type(mobbrmsd_state), allocatable :: s(:, :)
    integer(kind=IK)                  :: i, j

    mob = mobbrmsd(blocks)
    allocate (s(n_target, n_target))

    call mobbrmsd_min_span_tree( &
   &  n_target, mob%h, s, X, W, &
   &  cutoff=cutoff, &
   &  difflim=difflim, &
   &  maxeval=maxeval, &
   &  edges=edges, &
   &  weights=weights, &
   &  verbose=verbose)

    header = mob%h%dump()
    do concurrent(i=1:n_target, j=1:n_target)
      int_states(:, i, j) = s(i, j)%dump()
      float_states(:, i, j) = s(i, j)%dump_real()
    end do

  end subroutine min_span_tree

end module driver

