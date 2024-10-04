module mod_mobbrmsd_lapack
  use mod_mobbrmsd_lapack_routines, only: RK
#ifdef USE_REAL32
  use mod_mobbrmsd_lapack_routines, only: &
    &   SGEMM => mobbrmsd_SGEMM, &
    &   SGESVD => mobbrmsd_SGESVD, &
    &   SORMQR => mobbrmsd_SORMQR, &
    &   SGEQRF => mobbrmsd_SGEQRF, &
    &   SGETRF => mobbrmsd_SGETRF
#else
  use mod_mobbrmsd_lapack_routines, only: &
    &   DGEMM => mobbrmsd_DGEMM, &
    &   DGESVD => mobbrmsd_DGESVD, &
    &   DORMQR => mobbrmsd_DORMQR, &
    &   DGEQRF => mobbrmsd_DGEQRF, &
    &   DGETRF => mobbrmsd_DGETRF
#endif
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
end module mod_mobbrmsd_lapack

