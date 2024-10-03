!> \brief \b LA_CONSTANTS is a module for the scaling constants for the compiled Fortran single and double precisions
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date May 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!> Nick Papior, Technical University of Denmark, DK
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
module LA_CONSTANTS
!  -- LAPACK auxiliary module --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  Standard constants for
  use, intrinsic :: ISO_FORTRAN_ENV, only: SP => REAL32, DP => REAL64
  implicit none
  public
!
!&<
!
! integer, parameter      :: SP      = KIND(1.E0)
!
  real(SP), parameter     :: SZERO   = 0.0_SP
  real(SP), parameter     :: SHALF   = 0.5_SP
  real(SP), parameter     :: SONE    = 1.0_SP
  real(SP), parameter     :: STWO    = 2.0_SP
  real(SP), parameter     :: STHREE  = 3.0_SP
  real(SP), parameter     :: SFOUR   = 4.0_SP
  real(SP), parameter     :: SEIGHT  = 8.0_SP
  real(SP), parameter     :: STEN    = 10.0_SP
  complex(SP), parameter  :: CZERO   = (0.0_SP, 0.0_SP)
  complex(SP), parameter  :: CHALF   = (0.5_SP, 0.0_SP)
  complex(SP), parameter  :: CONE    = (1.0_SP, 0.0_SP)
  character(1), parameter :: SPREFIX = 'S'
  character(1), parameter :: CPREFIX = 'C'

!  Scaling constants
  real(SP), parameter     :: SULP    = EPSILON(0._SP)
  real(SP), parameter     :: SEPS    = sulp * 0.5_SP
  real(SP), parameter     :: SSAFMIN = real(RADIX(0._SP), SP)**MAX( MINEXPONENT(0._SP) - 1, &
                                            1 - MAXEXPONENT(0._SP) &
                                           )
  real(SP), parameter     :: SSAFMAX = SONE / SSAFMIN
  real(SP), parameter     :: SSMLNUM = SSAFMIN / SULP
  real(SP), parameter     :: SBIGNUM = SSAFMAX * SULP
  real(SP), parameter     :: SRTMIN  = SQRT(SSMLNUM)
  real(SP), parameter     :: SRTMAX  = SQRT(SBIGNUM)
!
!  Blue's scaling constants
  real(SP), parameter     :: STSML   = real(RADIX(0._SP), SP)**CEILING( &
                                       (MINEXPONENT(0._SP) - 1) * 0.5_SP)
  real(SP), parameter     :: STBIG   = real(RADIX(0._SP), SP)**FLOOR( &
                             (MAXEXPONENT(0._SP) - DIGITS(0._SP) + 1) * 0.5_SP)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
  real(SP), parameter     :: SSSML   = real(RADIX(0._SP), SP)**(-FLOOR( &
                                            (MINEXPONENT(0._SP) - DIGITS(0._SP)) * 0.5_SP))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
  real(SP), parameter     :: ssbig   = real(RADIX(0._SP), SP)**(-CEILING( &
                                            (MAXEXPONENT(0._SP) + DIGITS(0._SP) - 1) * 0.5_SP))

!  Standard constants for
! integer, parameter      :: DP = KIND(1.D0)

  real(DP), parameter     :: DZERO   = 0.0_DP
  real(DP), parameter     :: DHALF   = 0.5_DP
  real(DP), parameter     :: DONE    = 1.0_DP
  real(DP), parameter     :: DTWO    = 2.0_DP
  real(DP), parameter     :: DTHREE  = 3.0_DP
  real(DP), parameter     :: DFOUR   = 4.0_DP
  real(DP), parameter     :: DEIGHT  = 8.0_DP
  real(DP), parameter     :: DTEN    = 10.0_DP
  complex(DP), parameter  :: ZZERO   = (0.0_DP, 0.0_DP)
  complex(DP), parameter  :: ZHALF   = (0.5_DP, 0.0_DP)
  complex(DP), parameter  :: ZONE    = (1.0_DP, 0.0_DP)
  character(1), parameter :: DPREFIX = 'D'
  character(1), parameter :: ZPREFIX = 'Z'

!  Scaling constants
  real(DP), parameter     :: DULP = EPSILON(0._DP)
  real(DP), parameter     :: DEPS = DULP * 0.5_DP
  real(DP), parameter     :: DSAFMIN = real(RADIX(0._DP), dp)**MAX( &
                                            MINEXPONENT(0._DP) - 1, &
                                            1 - MAXEXPONENT(0._DP) &
                                           )
  real(DP), parameter     :: DSAFMAX = DONE / DSAFMIN
  real(DP), parameter     :: DSMLNUM = DSAFMIN / DULP
  real(DP), parameter     :: DBIGNUM = DSAFMAX * DULP
  real(DP), parameter     :: DRTMIN  = SQRT(DSMLNUM)
  real(DP), parameter     :: DRTMAX  = SQRT(DBIGNUM)

!  Blue's scaling constants
  real(DP), parameter     :: DTSML = real(RADIX(0._DP), dp)**CEILING( &
                                          (MINEXPONENT(0._DP) - 1) * 0.5_DP)
  real(DP), parameter     :: DTBIG = real(RADIX(0._DP), dp)**FLOOR( &
                                          (MAXEXPONENT(0._DP) - DIGITS(0._DP) + 1) * 0.5_DP)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
  real(DP), parameter     :: DSSML = real(RADIX(0._DP), dp)**(-FLOOR( &
                                          (MINEXPONENT(0._DP) - DIGITS(0._DP)) * 0.5_DP))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
  real(DP), parameter     :: DSBIG = real(RADIX(0._DP), dp)**(-CEILING( &
                                          (MAXEXPONENT(0._DP) + DIGITS(0._DP) - 1) * 0.5_DP))

!&>
end module LA_CONSTANTS
