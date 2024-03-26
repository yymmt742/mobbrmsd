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
  use mod_params, only: D, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_permutation
  use mod_rotation_matrix
  implicit none
  private
  public :: sample
  public :: covmat
  public :: gcov
  public :: SO3
  public :: eye
  public :: swp
  public :: brute_sd
!
contains
!
  function sample(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: cnt(d)
    real(RK)                :: res(d, n)
    integer(IK)             :: i
    call RANDOM_NUMBER(res)
    cnt = SUM(res, 2) / n
    do concurrent(i=1:n)
      res(:, i) = res(:, i) - cnt
    end do
  end function sample
!
  function covmat(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: res(d, d)
    res(:, :) = MATMUL(sample(d, n), TRANSPOSE(sample(d, n)))
  end function covmat
!
  function gcov(d, n) result(res)
    integer(IK), intent(in) :: d, n
    real(RK)                :: res(d * d + 1)
    real(RK)                :: x(d, n), y(d, n)
    x = sample(d, n)
    y = sample(d, n)
    res(1) = SUM(x * x) + SUM(y * y)
    res(2:) = [MATMUL(y, TRANSPOSE(x))]
  end function gcov
!
  function SO3() result(res)
    real(RK) :: a(4), c, t, s, res(3, 3)
    call RANDOM_NUMBER(a)
    a(:3) = a(:3) / SQRT(DOT_PRODUCT(a(:3), a(:3)))
    a(4) = (a(4) + a(4)) * ACOS(0.0_RK)
    c = COS(a(4))
    t = ONE - c
    s = SIN(a(4))
    res(:, 1) = [c + t * a(1) * a(1), t * a(1) * a(2) - s * a(3), t * a(1) * a(3) + s * a(2)]
    res(:, 2) = [t * a(1) * a(2) + s * a(3), c + t * a(2) * a(2), t * a(2) * a(3) - s * a(1)]
    res(:, 3) = [t * a(1) * a(3) - s * a(2), t * a(2) * a(3) + s * a(1), c + t * a(3) * a(3)]
  end function SO3
!
  pure function eye(d) result(res)
    integer,intent(in) :: d
    real(RK)           :: res(d, d)
    integer            :: i, j
    do concurrent(j=1:d, i=1:d)
      res(i, j) = MERGE(ONE, ZERO, i == j)
    enddo
  end function eye
!
  pure function sd(m, n, X, Y) result(res)
    integer(IK), intent(in) :: m, n
    real(RK), intent(in)    :: X(D, m, n), Y(D, m, n)
    real(RK)                :: G, C(D, D), W(100), res
    G = SUM(X * X) + SUM(Y * Y)
    C = MATMUL(RESHAPE(Y, [D, m * n]), TRANSPOSE(RESHAPE(X, [D, m * n])))
    call estimate_sdmin(G, C, w)
    res = w(1)
  end function sd
!
  pure function brute_sd(m, n, s, sym, X, Y) result(res)
    integer(IK), intent(in) :: m, n, s, sym(m, s)
    real(RK), intent(in)    :: X(D, m, n), Y(D, m, n)
    real(RK)                :: res
    type(permutation)       :: per
    integer(IK)             :: map(n)
    per = permutation(n, n)
    res = RHUGE
    map = 1
    do while (.not. per%endl())
      do
        res = MIN(res, sd(m, n, X, swp(m, n, s, per%id, map, sym, Y)))
        call map_next(n, s, map)
        if (ALL(map == 1)) exit
      enddo
      call per%next()
    end do
  end function brute_sd
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
    integer(IK), intent(in) :: m, n, s, per(n), map(n), sym(m, s)
    real(RK), intent(in)    :: X(D, m, n)
    real(RK)                :: res(D, m, n)
    integer(IK)             :: i
    do i = 1, n
      res(:, sym(:, map(i)), per(i)) = X(:, :, i)
    end do
  end function swp
!
end module mod_testutil

