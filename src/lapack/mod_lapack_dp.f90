module mod_lapack_sp
  implicit none
  private
  public :: DGEMM, &
          & DGESVD, &
          & DORMQR, &
          & DGEQRF, &
          & DGETRF
!&<
  !| Single precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(15)
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
  character, parameter :: DPREFIX = 'D'
!  Scaling constants
  real(RK), parameter  :: DULP    = EPSILON(ZERO)
  real(RK), parameter  :: DEPS    = DULP * HALF
  real(RK), parameter  :: DSAFMIN = real(RADIX(ZERO), RK)**MAX( MINEXPONENT(ZERO) - 1, &
                                        & 1 - MAXEXPONENT(ZERO) )
  real(RK), parameter  :: DSAFMAX = ONE / DSAFMIN
  real(RK), parameter  :: DSMLNUM = DSAFMIN / DULP
  real(RK), parameter  :: DBIGNUM = DSAFMAX * DULP
  real(RK), parameter  :: DRTMIN  = SQRT(DSMLNUM)
  real(RK), parameter  :: DRTMAX  = SQRT(DBIGNUM)
!
!  Blue's scaling constants
  real(RK), parameter  :: DTSML   = real(RADIX(ZERO), RK)**CEILING( &
                                      & (MINEXPONENT(ZERO) - 1) * HALF)
  real(RK), parameter  :: DTBIG   = real(RADIX(ZERO), RK)**FLOOR( &
                                      & (MAXEXPONENT(ZERO) - DIGITS(ZERO) + 1) * HALF)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
  real(RK), parameter  :: DSSML   = real(RADIX(ZERO), RK)**(-FLOOR( &
                                   &     (MINEXPONENT(ZERO) - DIGITS(ZERO)) * HALF))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
  real(RK), parameter  :: DSBIG   = real(RADIX(ZERO), RK)**(-CEILING( &
                                   &    (MAXEXPONENT(ZERO) + DIGITS(ZERO) - 1) * HALF))
!&>
contains
  include "util/dlamch.f90"
  include "lapack_routine/ieeeck.f90"
  include "lapack_routine/ilaenv.f90"
  include "lapack_routine/iladlc.f90"
  include "lapack_routine/iladlr.f90"
  include "lapack_routine/iparmq.f90"
  include "lapack_routine/idamax.f90"
  include "lapack_routine/lsame.f90"
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
!
  pure elemental function LA_ISNAN(x)
    real(RK), intent(in) :: x
    logical              :: LA_ISNAN
    LA_ISNAN = (x /= x)
  end function LA_ISNAN
!
end module mod_lapack_sp

