module lapack_sp
  implicit none
  private
  public :: SGEMM, &
          & SGESVD, &
          & SORMQR, &
          & SGEQRF, &
          & SGETRF
!&<
  !| Single precision kind.
  integer, parameter   :: RK = SELECTED_REAL_KIND(6)
!
  real(RK), parameter  :: ZERO    = 0.0_RK
  real(RK), parameter  :: QURTR   = 0.250E0
  real(RK), parameter  :: THIRD   = 0.3330E0
  real(RK), parameter  :: HALF    = 0.5_RK
  real(RK), parameter  :: ONE     = 1.0_RK
  real(RK), parameter  :: TWO     = 2.0_RK
  real(RK), parameter  :: THREE   = 3.0_RK
  real(RK), parameter  :: FOUR    = 4.0_RK
  real(RK), parameter  :: EIGHT   = 8.0_RK
  real(RK), parameter  :: TEN     = 10.0_RK
  real(RK), parameter  :: HUNDRD  = 100.0E0
  character, parameter :: SPREFIX = 'S'

!  Scaling constants
  real(RK), parameter  :: SULP    = EPSILON(ZERO)
  real(RK), parameter  :: SEPS    = SULP * HALF
  real(RK), parameter  :: SSAFMIN = real(RADIX(ZERO), RK)**MAX( MINEXPONENT(ZERO) - 1, &
                                        & 1 - MAXEXPONENT(ZERO) )
  real(RK), parameter  :: SSAFMAX = ONE / SSAFMIN
  real(RK), parameter  :: SSMLNUM = SSAFMIN / SULP
  real(RK), parameter  :: SBIGNUM = SSAFMAX * SULP
  real(RK), parameter  :: SRTMIN  = SQRT(SSMLNUM)
  real(RK), parameter  :: SRTMAX  = SQRT(SBIGNUM)
!
!  Blue's scaling constants
  real(RK), parameter  :: STSML   = real(RADIX(ZERO), RK)**CEILING( &
                                      & (MINEXPONENT(ZERO) - 1) * HALF)
  real(RK), parameter  :: STBIG   = real(RADIX(ZERO), RK)**FLOOR( &
                                      & (MAXEXPONENT(ZERO) - DIGITS(ZERO) + 1) * HALF)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
  real(RK), parameter  :: SSSML   = real(RADIX(ZERO), RK)**(-FLOOR( &
                                   &     (MINEXPONENT(ZERO) - DIGITS(ZERO)) * HALF))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
  real(RK), parameter  :: SSBIG   = real(RADIX(ZERO), RK)**(-CEILING( &
                                   &    (MAXEXPONENT(ZERO) + DIGITS(ZERO) - 1) * HALF))
!&>
contains
  include "util/slamch.f90"
  include "lapack_routine/ieeeck.f90"
  include "lapack_routine/ilaenv.f90"
  include "lapack_routine/ilaslc.f90"
  include "lapack_routine/ilaslr.f90"
  include "lapack_routine/iparmq.f90"
  include "lapack_routine/isamax.f90"
  include "lapack_routine/lsame.f90"
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
end module lapack_sp

