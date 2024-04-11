module driver
  !$ use omp_lib
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd, only : mobbrmsd, &
 &                         mobbrmsd_n_dims, &
 &                         mobbrmsd_run, &
 &                         mobbrmsd_swap_y, &
 &                         mol_block_input, &
 &                         mol_block_input_add, &
 &                         mobbrmsd_state, &
 &                         mobbrmsd_header, &
 &                         setup_dimension_ => setup_dimension
!
  implicit none
  private
  public setup_dimension
  public add_molecule
  public clear_molecule
  public n_atoms
  public n_dims
  public workmemory_length
  public state_vector_lengthes
  public rmsd
  public bounds
  public n_eval
  public log_eval_ratio
  public run
  public batch_run
!
  type(mol_block_input), allocatable :: blocks(:)
!
contains
!
  !| setup dimension
  subroutine setup_dimension(d)
    integer(kind=ik), intent(in) :: d
!
    call setup_dimension_(d)
!
  end subroutine setup_dimension
!
  !| add molecule
  subroutine add_molecule(n, M, s, sym)
    integer(kind=ik), intent(in) :: n
    integer(kind=ik), intent(in) :: M
    integer(kind=ik), intent(in) :: s
    integer(kind=ik), intent(in), optional :: sym(n * (s - 1))
!
    if (.not. ALLOCATED(blocks)) allocate (blocks(0))
    call mol_block_input_add(blocks, n, M, reshape(sym, [n, s - 1]))
!
  end subroutine add_molecule
!
  !| clear molecule
  subroutine clear_molecule()
!
    if (ALLOCATED(blocks)) deallocate (blocks)
!
  end subroutine clear_molecule
!
  !| Returns total number of atoms
  pure subroutine n_dims(n_dim)
    integer(kind=ik), intent(out) :: n_dim
!
    n_dim = mobbrmsd_n_dims()
!
  end subroutine n_dims
!
  !| Returns total number of atoms
  pure subroutine n_atoms(n_atom)
    integer(kind=ik), intent(out) :: n_atom
    type(mobbrmsd)                :: mob
!
    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      n_atom = mob%h%n_atoms()
    else
      n_atom = 0
    end if
!
  end subroutine n_atoms
!
  !| Returns total number of atoms
  pure subroutine workmemory_length(nmem)
    integer(kind=ik), intent(out) :: nmem
    type(mobbrmsd)                :: mob
    mob = mobbrmsd(blocks)
    nmem = mob%h%memsize()
  end subroutine workmemory_length
!
  pure subroutine state_vector_lengthes(n_header, n_int, n_float)
    integer(kind=ik), intent(out) :: n_header
    integer(kind=ik), intent(out) :: n_int
    integer(kind=ik), intent(out) :: n_float
    type(mobbrmsd)                :: mob
!
    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      n_header = SIZE(mob%h%dump())
      n_int = SIZE(mob%s%dump())
      n_float = SIZE(mob%s%dump_real())
    else
      n_header = 0
      n_int = 0
      n_float = 0
    endif
!
  end subroutine state_vector_lengthes
!
!| return rmsd
  subroutine rmsd(int_states, float_states, res, n_int, n_float)
    integer(kind=ik), intent(in) :: n_int
    integer(kind=ik), intent(in) :: n_float
    integer(kind=ik), intent(in) :: int_states(n_int)
    real(kind=rk), intent(in)    :: float_states(n_float)
    real(kind=rk), intent(out)   :: res
    type(mobbrmsd_state)         :: s
!
    call s%load(n_int, int_states, float_states)
    res = s%rmsd()
!
  end subroutine rmsd
!
!| return bounds
  subroutine bounds(int_states, float_states, res, n_int, n_float)
    integer(kind=ik), intent(in) :: n_int
    integer(kind=ik), intent(in) :: n_float
    integer(kind=ik), intent(in) :: int_states(n_int)
    real(kind=rk), intent(in)    :: float_states(n_float)
    real(kind=rk), intent(out)   :: res(2)
    type(mobbrmsd_state)         :: s
!
    call s%load(n_int, int_states, float_states)
    res(1) = s%lowerbound()
    res(2) = s%upperbound()
!
  end subroutine bounds
!
!| return n_eval
  subroutine n_eval(int_states, float_states, res, n_int, n_float)
    integer(kind=ik), intent(in)  :: n_int
    integer(kind=ik), intent(in)  :: n_float
    integer(kind=ik), intent(in)  :: int_states(n_int)
    real(kind=rk), intent(in)     :: float_states(n_float)
    integer(kind=ik), intent(out) :: res
    type(mobbrmsd_state)          :: s
!
    call s%load(n_int, int_states, float_states)
    res = s%n_eval()
!
  end subroutine n_eval
!
!| return log_eval_ratio
  subroutine log_eval_ratio(int_states, float_states, res, n_int, n_float)
    integer(kind=ik), intent(in) :: n_int
    integer(kind=ik), intent(in) :: n_float
    integer(kind=ik), intent(in) :: int_states(n_int)
    real(kind=rk), intent(in)    :: float_states(n_float)
    real(kind=rk), intent(out)   :: res
    type(mobbrmsd_state)         :: s
!
    call s%load(n_int, int_states, float_states)
    res = s%log_eval_ratio()
!
  end subroutine log_eval_ratio
!
  !| single run with working memory
  subroutine run(n_dim, n_atom, n_header, n_int, n_float, n_mem, x, y, w,&
 &               cutoff, difflim, maxeval, header, int_states, float_states)
    integer(kind=ik), intent(in)  :: n_dim
    integer(kind=ik), intent(in)  :: n_atom
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: n_int
    integer(kind=ik), intent(in)  :: n_float
    integer(kind=ik), intent(in)  :: n_mem
    real(kind=rk), intent(in)     :: x(n_dim, n_atom)
    !! reference coordinate
    real(kind=rk), intent(in)     :: y(n_dim, n_atom)
    !! target coordinate
    real(kind=rk), intent(inout)  :: W(n_mem)
    !! work memory
    integer(kind=ik), intent(in)  :: maxeval
    real(kind=rk), intent(in)     :: cutoff
    real(kind=rk), intent(in)     :: difflim
    integer(kind=ik), intent(out) :: header(n_header)
    integer(kind=ik), intent(out) :: int_states(n_int)
    real(kind=rk), intent(out)    :: float_states(n_float)
    type(mobbrmsd)                :: mob
!
    mob = mobbrmsd(blocks)
    call mobbrmsd_run(mob%h, mob%s, X, Y, w, cutoff, difflim, maxeval)
!
    header = mob%h%dump()
    int_states = mob%s%dump()
    float_states = mob%s%dump_real()

  end subroutine run
!
  !| batch parallel run
  subroutine batch_run(n_dim, n_atom, ntarget, n_header, n_int, n_float, x, y, &
 &                     cutoff, difflim, maxeval, rotate_y, header, int_states, float_states)
    integer(kind=ik), intent(in)  :: n_dim
    integer(kind=ik), intent(in)  :: n_atom
    integer(kind=ik), intent(in)  :: ntarget
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: n_int
    integer(kind=ik), intent(in)  :: n_float
    real(kind=rk), intent(in)     :: x(n_dim, n_atom)
    !! reference coordinate
    real(kind=rk), intent(inout)  :: y(n_dim, n_atom, ntarget)
    !! target coordinate
    integer(kind=ik), intent(in)  :: maxeval
    real(kind=rk), intent(in)     :: cutoff
    real(kind=rk), intent(in)     :: difflim
    logical, intent(in)           :: rotate_y
    integer(kind=ik), intent(out) :: header(n_header)
    integer(kind=ik), intent(out) :: int_states(n_int, ntarget)
    real(kind=rk), intent(out)    :: float_states(n_float, ntarget)
    real(kind=rk), allocatable    :: W(:, :)
    integer(kind=ik)              :: nmem, njob
    type(mobbrmsd)                :: mob
    integer(kind=ik)              :: i
!
    mob = mobbrmsd(blocks)
    nmem = mob%h%memsize()
    header = mob%h%dump()
!
    !$omp parallel
    njob = MIN(ntarget, MAX(omp_get_num_threads(), 1))
    !$omp end parallel
!
    allocate (W(nmem, njob))
!
    !$omp parallel do
    do i = 1, ntarget
      block
        integer(kind=ik)     :: ijob
        type(mobbrmsd_state) :: s
!
        ijob = omp_get_thread_num() + 1
        s = mob%s
        call mobbrmsd_run(mob%h, s, X, Y(1, 1, i), w(1, ijob), cutoff, difflim, maxeval)
        if(rotate_y) call mobbrmsd_swap_y(mob%h, s, Y(1, 1, i))
!
        int_states(:, i) = s%dump()
        float_states(:, i) = s%dump_real()
!
      end block
    end do
    !$omp end parallel do

  end subroutine batch_run
!
end module driver
