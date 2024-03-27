module driver
  !$ use omp_lib
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  implicit none
  private
  public add_molecule
  public run
!
  type(mobbrmsd_input), save :: inp
!
contains
!
  subroutine add_molecule(m, n, s, sym)
    integer(kind=ik), intent(in) :: m
    integer(kind=ik), intent(in) :: n
    integer(kind=ik), intent(in) :: s
    integer(kind=ik), intent(in) :: sym(*)
!
    call inp%add_molecule(m, n, sym=RESHAPE(sym(:m * (s - 1)), [m, s - 1]))
!
  end subroutine add_molecule
!
  !| run mobbrmsd
  subroutine run(x, y, cutoff, difflim, maxeval, n, rmsd, upper, lower, log_ratio, neval)
    real(kind=rk), intent(in)     :: x(*)
    !! reference coordinate
    real(kind=rk), intent(inout)  :: y(*)
    !! target coordinate
    integer(kind=ik), intent(in)  :: maxeval
    real(kind=rk), intent(in)     :: cutoff
    real(kind=rk), intent(in)     :: difflim
    integer(kind=ik), intent(in)  :: n
    real(kind=rk), intent(out)    :: rmsd(n)
    real(kind=rk), intent(out)    :: upper(n)
    real(kind=rk), intent(out)    :: lower(n)
    real(kind=rk), intent(out)    :: log_ratio(n)
    integer(kind=ik), intent(out) :: neval(n)
    real(kind=rk), allocatable    :: W(:, :)
    integer(kind=ik)              :: nmem
    integer(kind=ik)              :: njob
    type(mobbrmsd)                :: mob
    integer(kind=ik)              :: i, dmn
!
    mob = mobbrmsd(inp)
    nmem = mob%h%memsize()
    dmn = mob%h%n_atoms()
!
    !$omp parallel
    njob = MAX(omp_get_num_threads(), 1)
    !$omp end parallel
!
    allocate (W(nmem, njob))
!
    !$omp parallel do
    do i = 1, n
      block
        integer(kind=ik)     :: ijob, pnt
        type(mobbrmsd_state) :: s
        real(kind=rk)        :: rat(3)
!
        pnt = (i - 1) * dmn + 1
        ijob = omp_get_thread_num() + 1
        s = mob%s
        call mobrmsd_run(mob%h, s, X, Y(pnt), w(1, ijob), cutoff, difflim, maxeval)
!
        rmsd(i) = s%rmsd()
        upper(i) = s%upperbound()
        lower(i) = s%lowerbound()
        log_ratio(i) = s%log_eval_ratio()
!
      end block
    end do
    !$omp end parallel do

  end subroutine run
!
end module driver
