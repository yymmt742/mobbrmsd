module mod_permutation
use mod_params, only : IK, RK
use,intrinsic :: ISO_FORTRAN_ENV, only : I8 => INT64
implicit none
private
public :: permutation
!
  type :: permutation
!
    private
    integer( IK )                    :: n, r
    integer( IK ),allocatable        :: cy(:), ix(:), ip
    integer( IK ),allocatable,public :: id(:)
!
  contains
!
    procedure         :: next    => permutation_next
    procedure         :: reset   => permutation_reset
    procedure         :: endl    => permutation_endl
    final             :: permutation_destroy
!
  end type permutation
!
  interface permutation
    module procedure permutation_new
  end interface permutation
!
  real(RK),parameter              :: ULIM = LOG( REAL( HUGE( 0_I8 ), RK ) )
!
contains
!
  pure elemental function permutation_new( n, r ) result( res )
  integer( IK ),intent( in )          :: n
  integer( IK ),intent( in ),optional :: r
  type( permutation )                 :: res
!
    res%n = n
!
    if( PRESENT(r) )then ;   res%r = r
    else                 ;   res%r = n
    endif
!
    if( res%n<1     ) res%n = 0
    if( res%n<res%r ) res%r = 0
!
    ALLOCATE( res%id( res%r ), res%ix( 0 ), res%cy( 0 ) )
!
    call res%reset()
!
  end function permutation_new
!
  pure elemental subroutine permutation_next( this )
  class( permutation ),intent( inout ) :: this
  integer( IK )                        :: i, j, swp
!
    do while( .not.this%endl() )
!
      i = this%ip
!
      this%cy( i ) = this%cy( i ) - 1
!
      if( this%cy( i ) < 1 )then
!
        this%ix( i: ) = [ this%ix( i+1: ), this%ix( i:i ) ]
        this%cy( i )  = this%n + 1 - i
        this%ip       = this%ip - 1
!
      else
!
        j             = this%n + 1 - this%cy( i )
        swp           = this%ix( i )
        this%ix( i )  = this%ix( j )
        this%ix( j )  = swp
        this%ip       = this%r
!
        this%id( : )  = this%ix( :this%r )
!
        RETURN
!
       endif
!
    enddo
!
  end subroutine permutation_next
!
  pure elemental subroutine permutation_reset(this)
  class( permutation ),intent( inout ) :: this
  integer( IK )                        :: i
!
    this%ix = [( i, i = 1, this%n )]
    this%cy = [( i, i = this%n, this%n - this%r + 1, -1 )]
    this%id = [( i, i = 1, this%r )]
!
    this%ip = this%r
!
  end subroutine permutation_reset
!
  pure elemental function permutation_endl(this) result(res)
  class( permutation ),intent( in ) :: this
  logical                           :: res
!
    res = this%ip < 1
!
  end function permutation_endl
!
  pure elemental subroutine permutation_destroy(this)
  type( permutation ),intent( inout ) :: this
!
    this%n = 0
    this%r = 0
!
    if( ALLOCATED(this%id) ) DEALLOCATE( this%id )
    if( ALLOCATED(this%ix) ) DEALLOCATE( this%ix )
    if( ALLOCATED(this%cy) ) DEALLOCATE( this%cy )
!
  end subroutine permutation_destroy
!
end module mod_permutation
!
!| Utility functions for testing.
module mod_testutil
  use blas_lapack_interface, only : D, DD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, PI => RPI, RHUGE
  use mod_permutation
  use mod_rotation
  implicit none
  private
  public :: sample
  public :: covmat
  public :: gcov
  public :: SO
  public :: eye
  public :: sd
  public :: swp
  public :: brute_sd
  public :: brute_sd_double
!
  interface sample
    module procedure :: sample_2, sample_3
  end interface sample
!
contains
!
  function sample_2(n) result(res)
    integer(IK), intent(in) :: n
    real(RK)                :: cnt(D), res(D, n)
    integer(IK)             :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / real(n, RK)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    end do
  end function sample_2
!
  function sample_3(m, n) result(res)
    integer(IK), intent(in) :: m, n
    real(RK)                :: cnt(D), res(D, m, n)
    integer(IK)             :: i, j
    call RANDOM_NUMBER(res)
    cnt = SUM(RESHAPE(res, [D, m * n]), 2) / real(m * n, RK)
    do concurrent(i=1:m, j=1:n)
      res(:, i, j) = res(:, i, j) - cnt
    end do
  end function sample_3
!
  function covmat(n) result(res)
    integer(IK), intent(in) :: n
    real(RK)                :: res(D, D)
    res(:, :) = MATMUL(sample(n), TRANSPOSE(sample(n)))
  end function covmat
!
  function gcov(n) result(res)
    integer(IK), intent(in) :: n
    real(RK)                :: res(DD + 1)
    real(RK)                :: x(D, n), y(D, n)
    x = sample(n)
    y = sample(n)
    res(1) = SUM(x * x) + SUM(y * y)
    res(2:) = [MATMUL(y, TRANSPOSE(x))]
  end function gcov
!
  function SO() result(res)
    real(RK) :: res(D, D)
    select case (D)
    case (1)
      res(1, 1) = ONE
    case (2)
      call SO2(res)
    case (3)
      call SO3(res)
    case default
      res = eye()
    end select
  end function SO
!
  subroutine SO2(res)
    real(RK), intent(inout) :: res(2, *)
    real(RK)                :: a
    call RANDOM_NUMBER(a)
    a = a * PI
    res(:, 1) = [COS(a), SIN(a)]
    res(:, 2) = [-res(2, 1), res(1, 1)]
  end subroutine SO2
!
  subroutine SO3(res)
    real(RK), intent(inout) :: res(3, *)
    real(RK)                :: a(4), c, t, s
    call RANDOM_NUMBER(a)
    a(:3) = a(:3) / SQRT(DOT_PRODUCT(a(:3), a(:3)))
    a(4) = (a(4) + a(4)) * PI
    c = COS(a(4))
    t = ONE - c
    s = SIN(a(4))
    res(:, 1) = [c + t * a(1) * a(1), t * a(1) * a(2) - s * a(3), t * a(1) * a(3) + s * a(2)]
    res(:, 2) = [t * a(1) * a(2) + s * a(3), c + t * a(2) * a(2), t * a(2) * a(3) - s * a(1)]
    res(:, 3) = [t * a(1) * a(3) - s * a(2), t * a(2) * a(3) + s * a(1), c + t * a(3) * a(3)]
  end subroutine SO3
!
  pure function eye() result(res)
    real(RK)           :: res(D, D)
    integer            :: i, j
    do concurrent(j=1:D, i=1:D)
      res(i, j) = MERGE(ONE, ZERO, i == j)
    enddo
  end function eye
!
  pure function sd(n, X, Y) result(res)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: X(D, n), Y(D, n)
    real(RK)                :: G, C(D, D), res
    integer(IK)             :: nw
    G = SUM(X * X) + SUM(Y * Y)
    C = MATMUL(Y, TRANSPOSE(X))
    nw = sdmin_worksize()
    block
      real(rk) :: w(nw)
      call estimate_sdmin(G, C, w)
      res = w(1)
    end block
  end function sd
!
  pure function brute_sd(m, n, s, sym, X, Y) result(res)
    integer(IK), intent(in) :: m, n, s, sym(m * (s - 1))
    real(RK), intent(in)    :: X(D, m, n), Y(D, m, n)
    real(RK)                :: res
    type(permutation)       :: per
    integer(IK)             :: map(n)
    per = permutation(n, n)
    res = RHUGE
    map = 1
    do while (.not. per%endl())
      do
        res = MIN(res, sd(m * n, X, pws(m, n, s, per%id, map, sym, Y)))
        call map_next(n, s, map)
        if (ALL(map == 1)) exit
      enddo
      call per%next()
    end do
  end function brute_sd
!
  pure function brute_sd_double(m1, n1, s1, sym1, m2, n2, s2, sym2, X1, Y1, X2, Y2) result(res)
    integer(IK), intent(in) :: m1, n1, s1, sym1(m1 * (s1 - 1))
    integer(IK), intent(in) :: m2, n2, s2, sym2(m2 * (s2 - 1))
    real(RK), intent(in)    :: X1(D, m1, n1), Y1(D, m1, n1)
    real(RK), intent(in)    :: X2(D, m2, n2), Y2(D, m2, n2)
    real(RK)                :: res
    real(RK)                :: Z1(D, m1, n1), Z2(D, m2, n2)
    type(permutation)       :: per1, per2
    integer(IK)             :: map1(n1), map2(n2), nz
    res = RHUGE
    nz = m1 * n1 + m2 * n2

    per2 = permutation(n2, n2)
    map2 = 1
    do while (.not. per2%endl())
      do
        Z2 = pws(m2, n2, s2, per2%id, map2, sym2, Y2)
        per1 = permutation(n1, n1)
        map1 = 1
        do while (.not. per1%endl())
          do
            Z1 = pws(m1, n1, s1, per1%id, map1, sym1, Y1)
            res = MIN(res, sd(nz, RESHAPE([X1, X2], [D, nz]), RESHAPE([Z1, Z2], [D, nz])))
            call map_next(n1, s1, map1)
            if (ALL(map1 == 1)) exit
          end do
          call per1%next()
        end do
        call map_next(n2, s2, map2)
        if (ALL(map2 == 1)) exit
      end do
      call per2%next()
    end do
  end function brute_sd_double
!
  pure subroutine map_next(n, s, map)
    integer(IK), intent(in)    :: n, s
    integer(IK), intent(inout) :: map(n)
    integer(IK)                :: i
    do i = 1, n
      if (map(i) < s)then
        map(i) = map(i) + 1
        return
      endif
      map(i) = 1
    end do
  end subroutine map_next
!
  pure function swp(m, n, s, per, map, sym, X) result(res)
    integer(IK), intent(in) :: m, n, s, per(n), map(n), sym(m * (s - 1))
    real(RK), intent(in)    :: X(D, m, n)
    real(RK)                :: res(D, m, n)
    integer(IK)             :: i, sym1(m, s)
    sym1 = RESHAPE([[(i, i=1, m)], sym], SHAPE(sym1))
    do i = 1, n
      res(:, sym1(:, map(i)), per(i)) = X(:, :, i)
    end do
  end function swp
!
  pure function pws(m, n, s, per, map, sym, X) result(res)
    integer(IK), intent(in) :: m, n, s, per(n), map(n), sym(m * (s - 1))
    real(RK), intent(in)    :: X(D, m, n)
    real(RK)                :: res(D, m, n)
    integer(IK)             :: i, sym1(m, s)
    sym1 = RESHAPE([[(i, i=1, m)], sym], SHAPE(sym1))
    do i = 1, n
      res(:, :, i) = X(:, sym1(:, map(i)), per(i))
    end do
  end function pws
!
end module mod_testutil

