!| mod_mobbrmsd_lapack is lapack interface.
!  Only routines used for mob calculations are supported.
!  Here is a list of the routines currently required:
!
!  - xGEMM
!  - xGESVD
!  - xORMQR
!  - xGEQRF
!  - xGETRF
!
!  If the EXTERNAL_LAPACK macro is enabled,
!  the double precision ones will be replaced
!  with the corresponding routines from the external library.
!  The interfaces follow the standard [lapack api](https://www.netlib.org/lapack/).
module mod_mobbrmsd_lapack
#ifdef USE_REAL32
  use mod_mobbrmsd_lapack_routines_s, only: &
    &   SGEMM => mobbrmsd_SGEMM, &
    &   SGESVD => mobbrmsd_SGESVD, &
    &   SORMQR => mobbrmsd_SORMQR, &
    &   SGEQRF => mobbrmsd_SGEQRF, &
    &   SGETRF => mobbrmsd_SGETRF
#else
  use mod_mobbrmsd_lapack_routines_d, only: &
    &   DGEMM => mobbrmsd_DGEMM, &
    &   DGESVD => mobbrmsd_DGESVD, &
    &   DORMQR => mobbrmsd_DORMQR, &
    &   DGEQRF => mobbrmsd_DGEQRF, &
    &   DGETRF => mobbrmsd_DGETRF
#endif
  implicit none
  private
#ifdef USE_REAL32
  public :: SGEMM
  public :: SGESVD
  public :: SORMQR
  public :: SGEQRF
  public :: SGETRF
#else
  public :: DGEMM
  public :: DGESVD
  public :: DORMQR
  public :: DGEQRF
  public :: DGETRF
#endif
end module mod_mobbrmsd_lapack

