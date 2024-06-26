!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd_header
  use mod_dimspec_functions, only: D
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_bb_list
  use mod_bb_block
  implicit none
  public :: mobbrmsd_header
  public :: mobbrmsd_header_init
!| mobbrmsd_header
  type mobbrmsd_header
    private
    integer(IK)              :: d
    !! spatial dimension
    integer(IK), allocatable, public :: q(:)
    !! header array
    integer(IK), allocatable :: s(:)
    !! state template
  contains
    procedure :: n_dims => mobbrmsd_header_n_dims
    procedure :: n_block => mobbrmsd_header_n_block
    procedure :: n_atoms => mobbrmsd_header_n_atoms
    procedure :: log_n_nodes => mobbrmsd_header_log_n_nodes
    procedure :: frac_n_nodes => mobbrmsd_header_frac_n_nodes
    procedure :: exp_n_nodes => mobbrmsd_header_exp_n_nodes
    procedure :: memsize => mobbrmsd_header_memsize
    procedure :: state_template => mobbrmsd_header_state_template
    procedure :: dump => mobbrmsd_header_dump
    procedure :: load => mobbrmsd_header_load
    final     :: mobbrmsd_header_destroy
  end type mobbrmsd_header
!
  interface mobbrmsd_header
    module procedure :: mobbrmsd_header_new
  end interface mobbrmsd_header
!
contains
!| Returns spatial dimension
  pure function mobbrmsd_header_new(q, s) result(res)
    integer(IK), intent(in) :: q(:)
    !! mobbrmsd_header sequence
    integer(IK), intent(in) :: s(:)
    !! mobbrmsd_state template sequence
    type(mobbrmsd_header)   :: res
    call mobbrmsd_header_init(res, SIZE(q), q, SIZE(s), s)
  end function mobbrmsd_header_new
!
!| Initializer, for NVHPC.
  pure subroutine mobbrmsd_header_init(this, nq, q, ns, s)
    type(mobbrmsd_header), intent(inout) :: this
    integer(IK), intent(in)              :: nq
    !! dimension of q
    integer(IK), intent(in)              :: q(nq)
    !! mobbrmsd_header sequence
    integer(IK), intent(in)              :: ns
    !! dimension of s
    integer(IK), intent(in)              :: s(ns)
    !! mobbrmsd_state template sequence
    this%d = D
    this%q = q
    this%s = s
  end subroutine mobbrmsd_header_init
!
!| Returns spatial dimension
  pure elemental function mobbrmsd_header_n_dims(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK)                        :: res
    res = this%d
  end function mobbrmsd_header_n_dims
!
!| Returns number of molecular blocks
  pure elemental function mobbrmsd_header_n_block(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK)                        :: res
    res = bb_list_n_block(this%q)
  end function mobbrmsd_header_n_block
!
!| Returns header_memsize
  pure elemental function mobbrmsd_header_memsize(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK)                        :: res
    res = bb_list_memsize(this%q)
  end function mobbrmsd_header_memsize
!
!| Returns n_atoms
  pure elemental function mobbrmsd_header_n_atoms(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK)                        :: res
    res = bb_list_n_atoms(this%q)
  end function mobbrmsd_header_n_atoms
!
!| Returns log_n_nodes
  pure elemental function mobbrmsd_header_log_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    real(RK)                           :: res
    res = bb_list_log_n_nodes(this%q)
  end function mobbrmsd_header_log_n_nodes
!
!| returns number of nodes in fraction.
  pure function mobbrmsd_header_frac_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    real(RK)                           :: tmp, res
    tmp = LN_TO_L10 * bb_list_log_n_nodes(this%q)
    res = TEN**(tmp - real(INT(tmp), RK))
  end function mobbrmsd_header_frac_n_nodes
!
!| returns number of nodes in exp.
  pure function mobbrmsd_header_exp_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK)                        :: res
    res = INT(LN_TO_L10 * bb_list_log_n_nodes(this%q), IK)
  end function mobbrmsd_header_exp_n_nodes
!
!| dump state template as integer array
  pure function mobbrmsd_header_state_template(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK), allocatable           :: res(:)
!
    allocate (res, source=this%s)
!
  end function mobbrmsd_header_state_template
!
!| dump header as integer array
  pure function mobbrmsd_header_dump(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! this
    integer(IK), allocatable           :: res(:)
!
    allocate (res, source=[this%d, SIZE(this%q), this%q, this%s])
!
  end function mobbrmsd_header_dump
!
!| load integer array as header
  pure subroutine mobbrmsd_header_load(this, q)
    class(mobbrmsd_header), intent(inout) :: this
    !! this
    integer(IK), intent(in)               :: q(:)
    !! header array
    integer(IK)                           :: sq, i
    this%d = q(1)
    sq = q(2)
    if (sq < 1) then
      this%q = [(0, i=1, 0)]
      this%s = [(0, i=1, 0)]
      return
    else
      this%q = q(3:2 + sq)
      this%s = q(3 + sq:)
    end if
  end subroutine mobbrmsd_header_load
!
!| destructer
  pure elemental subroutine mobbrmsd_header_destroy(this)
    type(mobbrmsd_header), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine mobbrmsd_header_destroy
!
end module mod_mobbrmsd_header

