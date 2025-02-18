!| mod_mobbrmsd_lapack_routines is an internal implementation of the lapack routines.
!  It only supports the routines necessary to compute mobbrmsd and its dependencies.
!
!  Reference lapack is provided by [netlib](http://www.netlib.org/lapack/).
!
!  The interfaces follow the standard [lapack api](http://www.netlib.org/lapack/explore-html/).
!
module mod_mobbrmsd_lapack_routines_dp
  use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_IS_NAN
  implicit none
  private
  public :: mobbrmsd_ieeeck
  public :: mobbrmsd_ilaenv
  public :: mobbrmsd_iparmq
  public :: mobbrmsd_lsame
  public :: mobbrmsd_dlamch
  public :: mobbrmsd_idamax
  public :: mobbrmsd_iladlc
  public :: mobbrmsd_iladlr
  public :: mobbrmsd_dbdsqr
  public :: mobbrmsd_dcopy
  public :: mobbrmsd_dgebd2
  public :: mobbrmsd_dgebrd
  public :: mobbrmsd_dgelq2
  public :: mobbrmsd_dgelqf
  public :: mobbrmsd_dgemm
  public :: mobbrmsd_dgemv
  public :: mobbrmsd_dgeqr2
  public :: mobbrmsd_dgeqrf
  public :: mobbrmsd_dger
  public :: mobbrmsd_dgesvd
  public :: mobbrmsd_dgetrf
  public :: mobbrmsd_dgetrf2
  public :: mobbrmsd_dlabrd
  public :: mobbrmsd_dlacpy
  public :: mobbrmsd_dlange
  public :: mobbrmsd_dlapy2
  public :: mobbrmsd_dlarf
  public :: mobbrmsd_dlarfb
  public :: mobbrmsd_dlarfg
  public :: mobbrmsd_dlarft
  public :: mobbrmsd_dlartg
  public :: mobbrmsd_dlas2
  public :: mobbrmsd_dlascl
  public :: mobbrmsd_dlaset
  public :: mobbrmsd_dlasq1
  public :: mobbrmsd_dlasq2
  public :: mobbrmsd_dlasq3
  public :: mobbrmsd_dlasq4
  public :: mobbrmsd_dlasq5
  public :: mobbrmsd_dlasq6
  public :: mobbrmsd_dlasr
  public :: mobbrmsd_dlasrt
  public :: mobbrmsd_dlassq
  public :: mobbrmsd_dlasv2
  public :: mobbrmsd_dlaswp
  public :: mobbrmsd_dnrm2
  public :: mobbrmsd_dorg2r
  public :: mobbrmsd_dorgbr
  public :: mobbrmsd_dorgl2
  public :: mobbrmsd_dorglq
  public :: mobbrmsd_dorgqr
  public :: mobbrmsd_dorm2r
  public :: mobbrmsd_dormbr
  public :: mobbrmsd_dorml2
  public :: mobbrmsd_dormlq
  public :: mobbrmsd_dormqr
  public :: mobbrmsd_drot
  public :: mobbrmsd_dscal
  public :: mobbrmsd_dswap
  public :: mobbrmsd_dtrmm
  public :: mobbrmsd_dtrmv
  public :: mobbrmsd_dtrsm
!&<
!| Double precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(15)
  character, parameter :: PREFIX = 'D'
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
!
  real(RK), parameter  :: ULP    = EPSILON(ZERO)
  real(RK), parameter  :: EPS    = ULP * HALF
  real(RK), parameter  :: SAFMIN = real(RADIX(ZERO), RK)**MAX( MINEXPONENT(ZERO) - 1, &
                                       & 1 - MAXEXPONENT(ZERO) )
  real(RK), parameter  :: SAFMAX = ONE / SAFMIN
  real(RK), parameter  :: SMLNUM = SAFMIN / ULP
  real(RK), parameter  :: BIGNUM = SAFMAX * ULP
  real(RK), parameter  :: RTMIN  = SQRT(SMLNUM)
  real(RK), parameter  :: RTMAX  = SQRT(BIGNUM)
!
!  Blues scaling constants
!
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
  include "lapack_routine/mobbrmsd_dlamch.f90"
  include "lapack_routine/mobbrmsd_idamax.f90"
  include "lapack_routine/mobbrmsd_iladlc.f90"
  include "lapack_routine/mobbrmsd_iladlr.f90"
  include "lapack_routine/mobbrmsd_dbdsqr.f90"
  include "lapack_routine/mobbrmsd_dcopy.f90"
  include "lapack_routine/mobbrmsd_dgebd2.f90"
  include "lapack_routine/mobbrmsd_dgebrd.f90"
  include "lapack_routine/mobbrmsd_dgelq2.f90"
  include "lapack_routine/mobbrmsd_dgelqf.f90"
  include "lapack_routine/mobbrmsd_dgemm.f90"
  include "lapack_routine/mobbrmsd_dgemv.f90"
  include "lapack_routine/mobbrmsd_dgeqr2.f90"
  include "lapack_routine/mobbrmsd_dgeqrf.f90"
  include "lapack_routine/mobbrmsd_dger.f90"
  include "lapack_routine/mobbrmsd_dgesvd.f90"
  include "lapack_routine/mobbrmsd_dgetrf.f90"
  include "lapack_routine/mobbrmsd_dgetrf2.f90"
  include "lapack_routine/mobbrmsd_dlabrd.f90"
  include "lapack_routine/mobbrmsd_dlacpy.f90"
  include "lapack_routine/mobbrmsd_dlange.f90"
  include "lapack_routine/mobbrmsd_dlapy2.f90"
  include "lapack_routine/mobbrmsd_dlarf.f90"
  include "lapack_routine/mobbrmsd_dlarfb.f90"
  include "lapack_routine/mobbrmsd_dlarfg.f90"
  include "lapack_routine/mobbrmsd_dlarft.f90"
  include "lapack_routine/mobbrmsd_dlartg.f90"
  include "lapack_routine/mobbrmsd_dlas2.f90"
  include "lapack_routine/mobbrmsd_dlascl.f90"
  include "lapack_routine/mobbrmsd_dlaset.f90"
  include "lapack_routine/mobbrmsd_dlasq1.f90"
  include "lapack_routine/mobbrmsd_dlasq2.f90"
  include "lapack_routine/mobbrmsd_dlasq3.f90"
  include "lapack_routine/mobbrmsd_dlasq4.f90"
  include "lapack_routine/mobbrmsd_dlasq5.f90"
  include "lapack_routine/mobbrmsd_dlasq6.f90"
  include "lapack_routine/mobbrmsd_dlasr.f90"
  include "lapack_routine/mobbrmsd_dlasrt.f90"
  include "lapack_routine/mobbrmsd_dlassq.f90"
  include "lapack_routine/mobbrmsd_dlasv2.f90"
  include "lapack_routine/mobbrmsd_dlaswp.f90"
  include "lapack_routine/mobbrmsd_dnrm2.f90"
  include "lapack_routine/mobbrmsd_dorg2r.f90"
  include "lapack_routine/mobbrmsd_dorgbr.f90"
  include "lapack_routine/mobbrmsd_dorgl2.f90"
  include "lapack_routine/mobbrmsd_dorglq.f90"
  include "lapack_routine/mobbrmsd_dorgqr.f90"
  include "lapack_routine/mobbrmsd_dorm2r.f90"
  include "lapack_routine/mobbrmsd_dormbr.f90"
  include "lapack_routine/mobbrmsd_dorml2.f90"
  include "lapack_routine/mobbrmsd_dormlq.f90"
  include "lapack_routine/mobbrmsd_dormqr.f90"
  include "lapack_routine/mobbrmsd_drot.f90"
  include "lapack_routine/mobbrmsd_dscal.f90"
  include "lapack_routine/mobbrmsd_dswap.f90"
  include "lapack_routine/mobbrmsd_dtrmm.f90"
  include "lapack_routine/mobbrmsd_dtrmv.f90"
  include "lapack_routine/mobbrmsd_dtrsm.f90"
end module mod_mobbrmsd_lapack_routines_dp

