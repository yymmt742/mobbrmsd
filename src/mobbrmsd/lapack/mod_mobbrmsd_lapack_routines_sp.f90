!| mod_mobbrmsd_lapack_routines is an internal implementation of the lapack routines.
!  It only supports the routines necessary to compute mobbrmsd and its dependencies.
!
!  Reference lapack is provided by [netlib](http://www.netlib.org/lapack/).
!
!  The interfaces follow the standard [lapack api](http://www.netlib.org/lapack/explore-html/).
!
module mod_mobbrmsd_lapack_routines_sp
  use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_IS_NAN
  implicit none
  private
  public :: mobbrmsd_ieeeck
  public :: mobbrmsd_ilaenv
  public :: mobbrmsd_iparmq
  public :: mobbrmsd_lsame
  public :: mobbrmsd_slamch
  public :: mobbrmsd_isamax
  public :: mobbrmsd_ilaslc
  public :: mobbrmsd_ilaslr
  public :: mobbrmsd_sbdsqr
  public :: mobbrmsd_scombssq
  public :: mobbrmsd_scopy
  public :: mobbrmsd_sgebd2
  public :: mobbrmsd_sgebrd
  public :: mobbrmsd_sgelq2
  public :: mobbrmsd_sgelqf
  public :: mobbrmsd_sgemm
  public :: mobbrmsd_sgemv
  public :: mobbrmsd_sgeqr2
  public :: mobbrmsd_sgeqrf
  public :: mobbrmsd_sger
  public :: mobbrmsd_sgesvd
  public :: mobbrmsd_sgetrf
  public :: mobbrmsd_sgetrf2
  public :: mobbrmsd_slabrd
  public :: mobbrmsd_slacpy
  public :: mobbrmsd_slange
  public :: mobbrmsd_slapy2
  public :: mobbrmsd_slarf
  public :: mobbrmsd_slarfb
  public :: mobbrmsd_slarfg
  public :: mobbrmsd_slarft
  public :: mobbrmsd_slartg
  public :: mobbrmsd_slas2
  public :: mobbrmsd_slascl
  public :: mobbrmsd_slaset
  public :: mobbrmsd_slasq1
  public :: mobbrmsd_slasq2
  public :: mobbrmsd_slasq3
  public :: mobbrmsd_slasq4
  public :: mobbrmsd_slasq5
  public :: mobbrmsd_slasq6
  public :: mobbrmsd_slasr
  public :: mobbrmsd_slasrt
  public :: mobbrmsd_slassq
  public :: mobbrmsd_slasv2
  public :: mobbrmsd_slaswp
  public :: mobbrmsd_snrm2
  public :: mobbrmsd_sorg2r
  public :: mobbrmsd_sorgbr
  public :: mobbrmsd_sorgl2
  public :: mobbrmsd_sorglq
  public :: mobbrmsd_sorgqr
  public :: mobbrmsd_sorm2r
  public :: mobbrmsd_sormbr
  public :: mobbrmsd_sorml2
  public :: mobbrmsd_sormlq
  public :: mobbrmsd_sormqr
  public :: mobbrmsd_srot
  public :: mobbrmsd_sscal
  public :: mobbrmsd_sswap
  public :: mobbrmsd_strmm
  public :: mobbrmsd_strmv
  public :: mobbrmsd_strsm
!&<
!| Single precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(6)
  character, parameter :: PREFIX = 'S'
!
  real(RK), parameter  :: ZERO    = 0.0_RK
  real(RK), parameter  :: QURTR   = 0.250_RK
  real(RK), parameter  :: HALF    = 0.5_RK
  real(RK), parameter  :: ONE     = 1.0_RK
  real(RK), parameter  :: TWO     = 2.0_RK
  real(RK), parameter  :: THREE   = 3.0_RK
  real(RK), parameter  :: THIRD   = 0.3330_RK
  real(RK), parameter  :: FOUR    = 4.0_RK
  real(RK), parameter  :: EIGHT   = 8.0_RK
  real(RK), parameter  :: TEN     = 10.0_RK
  real(RK), parameter  :: HUNDRD  = 100.0_RK
!
!  Scaling constants
  real(RK), parameter  :: ULP    = EPSILON(ZERO)
  real(RK), parameter  :: EPS    = ULP * HALF
  real(RK), parameter  :: SAFMIN = real(RADIX(ZERO), RK)**MAX( MINEXPONENT(ZERO) - 1, &
                                       & 1 - MAXEXPONENT(ZERO) )
  real(RK), parameter  :: SAFMAX = ONE / SAFMIN
  real(RK), parameter  :: SMLNUM = SAFMIN / ULP
  real(RK), parameter  :: BIGNUM = SAFMAX * ULP
  real(RK), parameter  :: RTMIN  = SQRT(SMLNUM)
  real(RK), parameter  :: RTMAX  = SQRT(BIGNUM)
!  Blue's scaling constants
  real(RK), parameter  :: TSML   = real(RADIX(ZERO), RK)**CEILING( &
                                 & (MINEXPONENT(ZERO) - 1) * HALF)
  real(RK), parameter  :: TBIG   = real(RADIX(ZERO), RK)**FLOOR( &
                                 & (MAXEXPONENT(ZERO) - DIGITS(ZERO) + 1) * HALF)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
  real(RK), parameter  :: SSML   = real(RADIX(ZERO), RK)**(-FLOOR( &
                                 & (MINEXPONENT(ZERO) - DIGITS(ZERO)) * HALF))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
  real(RK), parameter  :: SBIG   = real(RADIX(ZERO), RK)**(-CEILING( &
                                 & (MAXEXPONENT(ZERO) + DIGITS(ZERO) - 1) * HALF))
!&>
contains
!
  include "lapack_routine/mobbrmsd_ieeeck.f90"
  include "lapack_routine/mobbrmsd_ilaenv.f90"
  include "lapack_routine/mobbrmsd_iparmq.f90"
  include "lapack_routine/mobbrmsd_lsame.f90"
  include "lapack_routine/mobbrmsd_slamch.f90"
  include "lapack_routine/mobbrmsd_isamax.f90"
  include "lapack_routine/mobbrmsd_ilaslc.f90"
  include "lapack_routine/mobbrmsd_ilaslr.f90"
  include "lapack_routine/mobbrmsd_sbdsqr.f90"
  include "lapack_routine/mobbrmsd_scombssq.f90"
  include "lapack_routine/mobbrmsd_scopy.f90"
  include "lapack_routine/mobbrmsd_sgebd2.f90"
  include "lapack_routine/mobbrmsd_sgebrd.f90"
  include "lapack_routine/mobbrmsd_sgelq2.f90"
  include "lapack_routine/mobbrmsd_sgelqf.f90"
  include "lapack_routine/mobbrmsd_sgemm.f90"
  include "lapack_routine/mobbrmsd_sgemv.f90"
  include "lapack_routine/mobbrmsd_sgeqr2.f90"
  include "lapack_routine/mobbrmsd_sgeqrf.f90"
  include "lapack_routine/mobbrmsd_sger.f90"
  include "lapack_routine/mobbrmsd_sgesvd.f90"
  include "lapack_routine/mobbrmsd_sgetrf.f90"
  include "lapack_routine/mobbrmsd_sgetrf2.f90"
  include "lapack_routine/mobbrmsd_slabrd.f90"
  include "lapack_routine/mobbrmsd_slacpy.f90"
  include "lapack_routine/mobbrmsd_slange.f90"
  include "lapack_routine/mobbrmsd_slapy2.f90"
  include "lapack_routine/mobbrmsd_slarf.f90"
  include "lapack_routine/mobbrmsd_slarfb.f90"
  include "lapack_routine/mobbrmsd_slarfg.f90"
  include "lapack_routine/mobbrmsd_slarft.f90"
  include "lapack_routine/mobbrmsd_slartg.f90"
  include "lapack_routine/mobbrmsd_slas2.f90"
  include "lapack_routine/mobbrmsd_slascl.f90"
  include "lapack_routine/mobbrmsd_slaset.f90"
  include "lapack_routine/mobbrmsd_slasq1.f90"
  include "lapack_routine/mobbrmsd_slasq2.f90"
  include "lapack_routine/mobbrmsd_slasq3.f90"
  include "lapack_routine/mobbrmsd_slasq4.f90"
  include "lapack_routine/mobbrmsd_slasq5.f90"
  include "lapack_routine/mobbrmsd_slasq6.f90"
  include "lapack_routine/mobbrmsd_slasr.f90"
  include "lapack_routine/mobbrmsd_slasrt.f90"
  include "lapack_routine/mobbrmsd_slassq.f90"
  include "lapack_routine/mobbrmsd_slasv2.f90"
  include "lapack_routine/mobbrmsd_slaswp.f90"
  include "lapack_routine/mobbrmsd_snrm2.f90"
  include "lapack_routine/mobbrmsd_sorg2r.f90"
  include "lapack_routine/mobbrmsd_sorgbr.f90"
  include "lapack_routine/mobbrmsd_sorgl2.f90"
  include "lapack_routine/mobbrmsd_sorglq.f90"
  include "lapack_routine/mobbrmsd_sorgqr.f90"
  include "lapack_routine/mobbrmsd_sorm2r.f90"
  include "lapack_routine/mobbrmsd_sormbr.f90"
  include "lapack_routine/mobbrmsd_sorml2.f90"
  include "lapack_routine/mobbrmsd_sormlq.f90"
  include "lapack_routine/mobbrmsd_sormqr.f90"
  include "lapack_routine/mobbrmsd_srot.f90"
  include "lapack_routine/mobbrmsd_sscal.f90"
  include "lapack_routine/mobbrmsd_sswap.f90"
  include "lapack_routine/mobbrmsd_strmm.f90"
  include "lapack_routine/mobbrmsd_strmv.f90"
  include "lapack_routine/mobbrmsd_strsm.f90"
end module mod_mobbrmsd_lapack_routines_sp

