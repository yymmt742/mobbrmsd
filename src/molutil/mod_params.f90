module mod_params
  use, intrinsic :: ISO_FORTRAN_ENV, only:  &
    &                STDIN => INPUT_UNIT,  &
    &                STDOUT => OUTPUT_UNIT, &
    &                STDERR => ERROR_UNIT,  &
    &                R4 => REAL32,  &
    &                R8 => REAL64,  &
    &                RQ => REAL128, &
    &                I1 => INT8, &
    &                I2 => INT16, &
    &                I4 => INT32, &
    &                I8 => INT64
  implicit none
  private
  public  :: R4, I8, RQ
  public  :: I1, I2, I4, I8
  public  :: RK, IK, LK
  public  :: STDIN, STDOUT, STDERR
  public  :: RZERO, RONE, RHALF, RFOUR, RHUGE
  public  :: NEWLINE, CNULL, CARRET, ESCSEQ
  public  :: FS_BOLD, FS_WEAK, FS_UNDER_LINE, FS_INVERT, FS_CROSSED_OUT, FS_RESET
  public  :: FC_BLACK, FC_RED, FC_GREEN, FC_YELLOW
  public  :: FC_MAGENTA, FC_CYAN, FC_WHITE
!
!&<
!
  integer, parameter          :: IK = I4
  integer, parameter          :: RK = R8
  integer, parameter          :: LK = KIND(.true.)
!
  real(RK), parameter         :: RZERO = 0.0_RK
  real(RK), parameter         :: RONE  = 1.0_RK
  real(RK), parameter         :: RHALF = 0.5_RK
  real(RK), parameter         :: RFOUR = 4.0_RK
!
  real(RK), parameter         :: RHUGE = HUGE(RZERO)
!
  character(*), parameter     :: NEWLINE        = NEW_LINE(' ')
  character(*), parameter     :: CNULL          = CHAR(0)
  character(*), parameter     :: CARRET         = CHAR(13)
  character(*), parameter     :: ESCSEQ         = CHAR(27)
!
  character(*), parameter     :: FS_BOLD        = ESCSEQ//'[1m'
  character(*), parameter     :: FS_WEAK        = ESCSEQ//'[2m'
  character(*), parameter     :: FS_UNDER_LINE  = ESCSEQ//'[4m'
  character(*), parameter     :: FS_INVERT      = ESCSEQ//'[7m'
  character(*), parameter     :: FS_CROSSED_OUT = ESCSEQ//'[9m'
  character(*), parameter     :: FS_RESET       = ESCSEQ//'[0m'
!
  character(*), parameter     :: FC_BLACK       = ESCSEQ//'[30m'
  character(*), parameter     :: FC_RED         = ESCSEQ//'[31m'
  character(*), parameter     :: FC_GREEN       = ESCSEQ//'[32m'
  character(*), parameter     :: FC_YELLOW      = ESCSEQ//'[33m'
  character(*), parameter     :: FC_BLUE        = ESCSEQ//'[34m'
  character(*), parameter     :: FC_MAGENTA     = ESCSEQ//'[35m'
  character(*), parameter     :: FC_CYAN        = ESCSEQ//'[36m'
  character(*), parameter     :: FC_WHITE       = ESCSEQ//'[37m'
!
!&>
!
end module mod_params
