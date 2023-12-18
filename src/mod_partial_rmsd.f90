!1e| Calculate min_{R,P,Q1,...,Q_n} |X-R@Y@f(P)@g(Q1,Q2,...,Qn)|^2 for R\in\R^{d,d}, P\in \R^{n,n} and Qi\in \R^{m,m}.
!  f(P) = I x P ( I is the identity and \in R^{m,m} ).
!  g(Q1,Q2,...,Qn) = Diag(Q1 Q2 ... Qn).
!  Here, R@R^T=I, det(R)=1, P@P^T=I, and Q@Q^T=I are satisfied.
module mod_partial_rmsd
  use mod_params, only: IK, RK, ONE => RONE, FOUR => RFOUR, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_Kabsch
  implicit none
  private
  public :: partial_rmsd
!
  type partial_rmsd
    private
    integer(IK)           :: d   = 0
    integer(IK)           :: n   = 0
    integer(IK)           :: fix = 0
    integer(IK)           :: yxt = 0
    integer(IK)           :: rot = 0
    real(RK), allocatable :: w(:)
  contains
!
    procedure         :: append  => partial_rmsd_append
    procedure         :: trace   => partial_rmsd_trace
    procedure         :: sd      => partial_rmsd_sd
    procedure         :: msd     => partial_rmsd_msd
    procedure         :: rmsd    => partial_rmsd_rmsd
    final             :: partial_rmsd_destroy
!
  end type partial_rmsd
!
  interface partial_rmsd
    module procedure partial_rmsd_new
  end interface partial_rmsd
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!| Calculate partial_rmsd.
  pure function partial_rmsd_new(d, n, X, Y) result(res)
    integer(IK), intent(in)           :: d, n
    !! mol_block
    real(RK), intent(in)              :: X(*)
    !! reference array.
    real(RK), intent(in)              :: Y(*)
    !! target array.
    type(partial_rmsd)                :: res
    real(RK)                          :: w(nwork(d))
    integer(IK)                       :: dd, dn
!
    if (d < 1 .or. n < 1) then
      allocate (res%w(0)); return
    end if
!
    dd = d * d
    dn = d * n
!
    res%d = d
    res%n = n
!
    res%fix = 1
    res%yxt = res%fix + 1
    res%rot = res%yxt + dd
!
    allocate (res%w(dd + dd + 1))
!
    res%w(res%fix) = ddot(dn, X, X) + ddot(dn, Y, Y)
!
!!! get optimal R
    call DGEMM('N', 'T', res%d, res%d, n, ONE, Y, res%d, X, res%d, ZERO, w, res%d)
    call Kabsch(res%d, w, w(dd + 1), w(dd + dd + 1))
    call copy(dd, w, res%w(res%yxt))
    call copy(dd, w(dd + 1), res%w(res%rot))
!
  end function partial_rmsd_new
!
  pure function partial_rmsd_append(this, n, X, Y) result(res)
    class(partial_rmsd), intent(in) :: this
    integer(IK), intent(in)         :: n
    !! mol_block
    real(RK), intent(in)            :: X(*)
    !! reference array.
    real(RK), intent(in)            :: Y(*)
    !! target array.
    type(partial_rmsd)              :: res
    real(RK)                        :: w(nwork(this%d))
    integer(IK)                     :: dd, dn
!
    if(this%d<1)then
      res = partial_rmsd_new(this%d, n, X, Y)
      return
    endif
!
    res%d = this%d
    res%n = this%n + n
!
    dd = res%d * res%d
    dn = res%d * n
!
    res%fix = this%fix
    res%yxt = this%yxt
    res%rot = this%rot
!
    allocate (res%w(dd + dd + 1))
!
    res%w(res%fix) = this%w(this%fix) + ddot(dn, X, X) + ddot(dn, Y, Y)
!
!!! get optimal R
    call copy(dd, this%w(this%yxt), w)
    call DGEMM('N', 'T', res%d, res%d, n, ONE, Y, res%d, X, res%d, ONE, w, res%d)
    call Kabsch(res%d, w, w(dd + 1), w(dd + dd + 1))
    call copy(dd, w, res%w(res%yxt))
    call copy(dd, w(dd + 1), res%w(res%rot))
!
  end function partial_rmsd_append
!
  pure elemental function partial_rmsd_trace(this) result(res)
    class(partial_rmsd), intent(in) :: this
    real(RK)                        :: res
    if (this%fix < 1) then
      res = ZERO
    else
      res = ddot(this%d**2, this%w(this%yxt), this%w(this%rot))
    endif
  end function partial_rmsd_trace
!
  pure elemental function partial_rmsd_sd(this) result(res)
    class(partial_rmsd), intent(in) :: this
    real(RK)                        :: res
    if (this%fix < 1) then
      res = ZERO
    else
      res = -this%trace()
      res = ABS(this%w(this%fix) + res + res)
    endif
  end function partial_rmsd_sd
!
  pure elemental function partial_rmsd_msd(this) result(res)
    class(partial_rmsd), intent(in) :: this
    real(RK)                        :: res
    if(this%n<1)then
      res = ZERO
    else
      res = this%sd() / this%n
    endif
  end function partial_rmsd_msd
!
  pure elemental function partial_rmsd_rmsd(this) result(res)
    class(partial_rmsd), intent(in) :: this
    real(RK)                        :: res
    res = SQRT(this%msd())
  end function partial_rmsd_rmsd
!
  pure elemental subroutine partial_rmsd_destroy(this)
    type(partial_rmsd), intent(inout) :: this
    this%d = 0
    this%n = 0
    this%fix = 0
    this%yxt = 0
    this%rot = 0
    if (ALLOCATED(this%w)) deallocate (this%w)
  end subroutine partial_rmsd_destroy
!
!!! utils
!
  pure elemental function nwork(d) result(res)
    integer(IK), intent(in) :: d
    integer(IK)             :: res
    res = 2 * d * d + Kabsch_worksize(d)
  end function nwork
!
  pure subroutine copy(d, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = source(i)
    end do
  end subroutine copy
!
  pure function ddot(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    if(d==1)then
      res = X(1) * Y(1)
    elseif(d==4)then
      res = X(1) * Y(1) + X(2) * Y(2) + &
        &   X(3) * Y(3) + X(4) * Y(4)
    elseif(d==9)then
      res = X(1) * Y(1) + X(2) * Y(2) + X(3) * Y(3) +&
        &   X(4) * Y(4) + X(5) * Y(5) + X(6) * Y(6) +&
        &   X(7) * Y(7) + X(8) * Y(8) + X(9) * Y(9)
    else
      block
        integer(IK) :: i
        res = ZERO
        do i = 1, d
          res = res + X(i) * Y(i)
        end do
      end block
    endif
  end function ddot
!
end module mod_partial_rmsd
