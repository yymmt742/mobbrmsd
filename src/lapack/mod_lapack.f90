module mod_lapack
  implicit none
  private
#ifdef USE_REAL32
  public :: SGEMM, &
          & SGESVD, &
          & SORMQR, &
          & SGEQRF, &
          & SGETRF
#else
  public :: DGEMM, &
          & DGESVD, &
          & DORMQR, &
          & DGEQRF, &
          & DGETRF
#endif
!
!&<
#ifdef USE_REAL32
!| Single precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(6)
#else
!| Double precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(15)
#endif
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
#ifdef USE_REAL32
  character, parameter :: PREFIX = 'S'
#else
  character, parameter :: PREFIX = 'D'
#endif
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
!
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
  pure elemental function LA_ISNAN(x)
    real(RK), intent(in) :: x
    logical              :: LA_ISNAN
    LA_ISNAN = (x /= x)
  end function LA_ISNAN
!
  include "lapack_routine/ieeeck.f90"
  include "lapack_routine/ilaenv.f90"
  include "lapack_routine/iparmq.f90"
  include "lapack_routine/lsame.f90"
#ifdef USE_REAL32
  include "lapack_routine/slamch.f90"
  include "lapack_routine/isamax.f90"
  include "lapack_routine/ilaslc.f90"
  include "lapack_routine/ilaslr.f90"
  include "lapack_routine/sbdsqr.f90"
  include "lapack_routine/scombssq.f90"
  include "lapack_routine/scopy.f90"
  include "lapack_routine/sgebd2.f90"
  include "lapack_routine/sgebrd.f90"
  include "lapack_routine/sgelq2.f90"
  include "lapack_routine/sgelqf.f90"
  include "lapack_routine/sgemm.f90"
  include "lapack_routine/sgemv.f90"
  include "lapack_routine/sgeqr2.f90"
  include "lapack_routine/sgeqrf.f90"
  include "lapack_routine/sger.f90"
  include "lapack_routine/sgesvd.f90"
  include "lapack_routine/sgetrf.f90"
  include "lapack_routine/sgetrf2.f90"
  include "lapack_routine/sisnan.f90"
  include "lapack_routine/slabrd.f90"
  include "lapack_routine/slacpy.f90"
  include "lapack_routine/slaisnan.f90"
  include "lapack_routine/slange.f90"
  include "lapack_routine/slapy2.f90"
  include "lapack_routine/slarf.f90"
  include "lapack_routine/slarfb.f90"
  include "lapack_routine/slarfg.f90"
  include "lapack_routine/slarft.f90"
  include "lapack_routine/slartg.f90"
  include "lapack_routine/slas2.f90"
  include "lapack_routine/slascl.f90"
  include "lapack_routine/slaset.f90"
  include "lapack_routine/slasq1.f90"
  include "lapack_routine/slasq2.f90"
  include "lapack_routine/slasq3.f90"
  include "lapack_routine/slasq4.f90"
  include "lapack_routine/slasq5.f90"
  include "lapack_routine/slasq6.f90"
  include "lapack_routine/slasr.f90"
  include "lapack_routine/slasrt.f90"
  include "lapack_routine/slassq.f90"
  include "lapack_routine/slasv2.f90"
  include "lapack_routine/slaswp.f90"
  include "lapack_routine/snrm2.f90"
  include "lapack_routine/sorg2r.f90"
  include "lapack_routine/sorgbr.f90"
  include "lapack_routine/sorgl2.f90"
  include "lapack_routine/sorglq.f90"
  include "lapack_routine/sorgqr.f90"
  include "lapack_routine/sorm2r.f90"
  include "lapack_routine/sormbr.f90"
  include "lapack_routine/sorml2.f90"
  include "lapack_routine/sormlq.f90"
  include "lapack_routine/sormqr.f90"
  include "lapack_routine/srot.f90"
  include "lapack_routine/sscal.f90"
  include "lapack_routine/sswap.f90"
  include "lapack_routine/strmm.f90"
  include "lapack_routine/strmv.f90"
  include "lapack_routine/strsm.f90"
#else
  include "lapack_routine/dlamch.f90"
  include "lapack_routine/idamax.f90"
  include "lapack_routine/iladlc.f90"
  include "lapack_routine/iladlr.f90"
  include "lapack_routine/dbdsqr.f90"
  include "lapack_routine/dcopy.f90"
  include "lapack_routine/dgebd2.f90"
  include "lapack_routine/dgebrd.f90"
  include "lapack_routine/dgelq2.f90"
  include "lapack_routine/dgelqf.f90"
  include "lapack_routine/dgemm.f90"
  include "lapack_routine/dgemv.f90"
  include "lapack_routine/dgeqr2.f90"
  include "lapack_routine/dgeqrf.f90"
  include "lapack_routine/dger.f90"
  include "lapack_routine/dgesvd.f90"
  include "lapack_routine/dgetrf.f90"
  include "lapack_routine/dgetrf2.f90"
  include "lapack_routine/disnan.f90"
  include "lapack_routine/dlabrd.f90"
  include "lapack_routine/dlacpy.f90"
  include "lapack_routine/dlaisnan.f90"
  include "lapack_routine/dlange.f90"
  include "lapack_routine/dlapy2.f90"
  include "lapack_routine/dlarf.f90"
  include "lapack_routine/dlarfb.f90"
  include "lapack_routine/dlarfg.f90"
  include "lapack_routine/dlarft.f90"
  include "lapack_routine/dlartg.f90"
  include "lapack_routine/dlas2.f90"
  include "lapack_routine/dlascl.f90"
  include "lapack_routine/dlaset.f90"
  include "lapack_routine/dlasq1.f90"
  include "lapack_routine/dlasq2.f90"
  include "lapack_routine/dlasq3.f90"
  include "lapack_routine/dlasq4.f90"
  include "lapack_routine/dlasq5.f90"
  include "lapack_routine/dlasq6.f90"
  include "lapack_routine/dlasr.f90"
  include "lapack_routine/dlasrt.f90"
  include "lapack_routine/dlassq.f90"
  include "lapack_routine/dlasv2.f90"
  include "lapack_routine/dlaswp.f90"
  include "lapack_routine/dnrm2.f90"
  include "lapack_routine/dorg2r.f90"
  include "lapack_routine/dorgbr.f90"
  include "lapack_routine/dorgl2.f90"
  include "lapack_routine/dorglq.f90"
  include "lapack_routine/dorgqr.f90"
  include "lapack_routine/dorm2r.f90"
  include "lapack_routine/dormbr.f90"
  include "lapack_routine/dorml2.f90"
  include "lapack_routine/dormlq.f90"
  include "lapack_routine/dormqr.f90"
  include "lapack_routine/drot.f90"
  include "lapack_routine/dscal.f90"
  include "lapack_routine/dswap.f90"
  include "lapack_routine/dtrmm.f90"
  include "lapack_routine/dtrmv.f90"
  include "lapack_routine/dtrsm.f90"
#endif
!
end module mod_lapack

