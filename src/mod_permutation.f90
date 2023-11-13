module mod_permutation
use mod_params, only : IK, RK
use,intrinsic :: ISO_FORTRAN_ENV, only : I8 => INT64
implicit none
private
public :: permutation
public :: npr, log_npr, factorial, log_factorial
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
    procedure         :: npr     => permutation_npr
    procedure         :: log_npr => permutation_log_npr
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
  pure function factorial( n ) result(res)
  integer( IK ),intent( in )          :: n
  real( RK )                          :: res
    res = EXP( log_npr( n ) )
  end function factorial
!
  pure function npr( n, r ) result(res)
  integer( IK ),intent( in )          :: n
  integer( IK ),intent( in ),optional :: r
  real( RK )                          :: res
    res = EXP( log_npr( n, r ) )
  end function npr
!
  pure function log_factorial( n ) result(res)
  integer( IK ),intent( in )          :: n
  real( RK )                          :: res
    res = log_npr( n )
  end function log_factorial
!
  pure function log_npr( n, r ) result(res)
  integer( IK ),intent( in )          :: n
  integer( IK ),intent( in ),optional :: r
  real( RK )                          :: res
  integer( IK )                       :: l,i
!
    res = 0D0
!
    if( PRESENT( r ) )then ; l = MAX( n - r + 1, 2 )
    else                   ; l = 2
    endif
!
    if( n<2 .or. n<l ) RETURN
!
    res = SUM( LOG( REAL( [( i, i=l,n )], RK ) ) )
!
  end function log_npr
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
  pure elemental function permutation_npr( this ) result(res)
  class( permutation ),intent( in ) :: this
  real( RK )                        :: res
    res = npr( this%n, this%r )
  end function permutation_npr
!
  pure elemental function permutation_log_npr( this ) result(res)
  class( permutation ),intent( in ) :: this
  real( RK )                        :: res
    res = log_npr( this%n, this%r )
  end function permutation_log_npr
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
