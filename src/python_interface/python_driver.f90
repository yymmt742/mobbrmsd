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
  public state_vector_lengthes
  public run
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
  subroutine add_molecule(m, n, s, sym)
    integer(kind=ik), intent(in) :: m
    integer(kind=ik), intent(in) :: n
    integer(kind=ik), intent(in) :: s
    integer(kind=ik), intent(in), optional :: sym(m * (s - 1))
    type(mobbrmsd)               :: mob
!
    if (.not. ALLOCATED(blocks)) allocate (blocks(0))
    call mol_block_input_add(blocks, m, n, reshape(sym, [m, s - 1]))
!
  end subroutine add_molecule
!
  !| clear molecule
  subroutine clear_molecule(m, n, s, sym)
    integer(kind=ik), intent(in) :: m
    integer(kind=ik), intent(in) :: n
    integer(kind=ik), intent(in) :: s
    integer(kind=ik), intent(in), optional :: sym(m * (s - 1))
    type(mobbrmsd)               :: mob
!
    if (ALLOCATED(blocks)) deallocate (blocks)
!
  end subroutine clear_molecule
!
  !| Returns total number of atoms
  pure subroutine n_dims(ndim)
    integer(kind=ik), intent(out) :: ndim
!
    ndim = mobbrmsd_n_dims()
!
  end subroutine n_dims
!
  !| Returns total number of atoms
  pure subroutine n_atoms(natom)
    integer(kind=ik), intent(out) :: natom
    type(mobbrmsd)                :: mob
!
    if (ALLOCATED(blocks)) then
      mob = mobbrmsd(blocks)
      natom = mob%h%n_atoms()
    else
      natom = 0
    end if
!
  end subroutine n_atoms
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
  !| run mobbrmsd
  subroutine run(ndim, natom, ntarget, n_header, n_int, n_float, x, y, &
 &               cutoff, difflim, maxeval, rotate_y, header, int_states, float_states)
    integer(kind=ik), intent(in)  :: ndim
    integer(kind=ik), intent(in)  :: natom
    integer(kind=ik), intent(in)  :: ntarget
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: n_int
    integer(kind=ik), intent(in)  :: n_float
    real(kind=rk), intent(in)     :: x(ndim, natom)
    !! reference coordinate
    real(kind=rk), intent(inout)  :: y(ndim, natom, ntarget)
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
        integer(kind=ik)     :: ijob, pnt
        type(mobbrmsd_state) :: s
        real(kind=rk)        :: rat(3)
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

  end subroutine run
!
end module driver
