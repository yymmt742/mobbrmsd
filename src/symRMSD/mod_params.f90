!| Module with constant collection.
module mod_params
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                R4 => REAL32,  &
    &                R8 => REAL64,  &
    &                RQ => REAL128, &
    &                I1 => INT8,    &
    &                I2 => INT16,   &
    &                I4 => INT32,   &
    &                I8 => INT64,   &
    &                STDIN => INPUT_UNIT,   &
    &                STDOUT => OUTPUT_UNIT, &
    &                STDERR => ERROR_UNIT
  implicit none
  private
  public  :: I1, I2, I4, I8
  public  :: R4, R8, RQ
  public  :: RK, IK, LK
  public  :: STDIN, STDOUT, STDERR
  public  :: RZERO, RONE, RHALF, RFOUR, RHUGE
  public  :: D, DD
! public  :: NEWLINE, CNULL, CARRET, ESCSEQ
! public  :: FS_BOLD, FS_WEAK, FS_UNDER_LINE, FS_INVERT, FS_CROSSED_OUT, FS_RESET
! public  :: FC_BLACK, FC_RED, FC_GREEN, FC_YELLOW
! public  :: FC_MAGENTA, FC_CYAN, FC_WHITE
!
!&<
!
  integer, parameter          :: IK = KIND(0)
  !! Selected integer kind.
  integer, parameter          :: RK = KIND(0.0_R8)
  !! Selected real kind.
  integer, parameter          :: LK = KIND(.true.)
  !! Selected logical kind.
!
  real(RK), parameter         :: RZERO = 0.0_RK
  !! Real zero.
  real(RK), parameter         :: RONE  = 1.0_RK
  !! Real one.
  real(RK), parameter         :: RHALF = 0.5_RK
  !! Real 1/2.
  real(RK), parameter         :: RFOUR = 4.0_RK
  !! Real four.
!
  real(RK), parameter         :: RHUGE = HUGE(RZERO)
  !! Real large number.
!
! character(*), parameter     :: NEWLINE        = NEW_LINE(' ')
! !! Line break char.
! character(*), parameter     :: CNULL          = CHAR(0)
! !! Null char.
! character(*), parameter     :: CARRET         = CHAR(13)
! !! Carriage return.
! character(*), parameter     :: ESCSEQ         = CHAR(27)
! !! Escape sequence.
!
! character(*), parameter     :: FS_BOLD        = ESCSEQ//'[1m'
! character(*), parameter     :: FS_WEAK        = ESCSEQ//'[2m'
! character(*), parameter     :: FS_UNDER_LINE  = ESCSEQ//'[4m'
! character(*), parameter     :: FS_INVERT      = ESCSEQ//'[7m'
! character(*), parameter     :: FS_CROSSED_OUT = ESCSEQ//'[9m'
! character(*), parameter     :: FS_RESET       = ESCSEQ//'[0m'
!
! character(*), parameter     :: FC_BLACK       = ESCSEQ//'[30m'
! character(*), parameter     :: FC_RED         = ESCSEQ//'[31m'
! character(*), parameter     :: FC_GREEN       = ESCSEQ//'[32m'
! character(*), parameter     :: FC_YELLOW      = ESCSEQ//'[33m'
! character(*), parameter     :: FC_BLUE        = ESCSEQ//'[34m'
! character(*), parameter     :: FC_MAGENTA     = ESCSEQ//'[35m'
! character(*), parameter     :: FC_CYAN        = ESCSEQ//'[36m'
! character(*), parameter     :: FC_WHITE       = ESCSEQ//'[37m'
!
  integer, save               :: D  = 3
  !! Spatial dimension
  integer, save               :: DD = 9
  !! Square spatial dimension
!
!&>
!
end module mod_params
