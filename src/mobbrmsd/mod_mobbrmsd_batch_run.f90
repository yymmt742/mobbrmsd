!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd_batch_run
!$ use omp_lib
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_mobbrmsd_header
  use mod_mobbrmsd_state
  use mod_mobbrmsd
  implicit none
  public :: mobbrmsd_batch_run
  public :: mobbrmsd_batch_tri_run
!
contains
  !| batch parallel run
  subroutine mobbrmsd_batch_run(n_reference, n_target, header, state, &
  &                             X, Y, W, &
  &                             cutoff, difflim, maxeval, &
  &                             remove_com, sort_by_g, &
  &                             rotate_y)
    integer(IK), intent(in)              :: n_reference
    !! number of reference coordinates
    integer(IK), intent(in)              :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(*)
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(inout)              :: Y(*)
    !! target coordinate
    real(RK), intent(inout)              :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional        :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional        :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    logical, intent(in), optional        :: rotate_y
    !! The search ends when ncount exceeds maxiter.
    integer(kind=IK)                     :: i, ipnt, ijob, xpnt, ypnt, wpnt, ldx, ldw, nlim
!
    if (n_reference < 1 .or. n_target < 1) return
    ldx = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
    nlim = n_reference * n_target
    i = 0
!
    !$omp parallel private(ipnt, ijob, xpnt, ypnt, wpnt)
    do
      !$omp critical
      i = i + 1
      ipnt = i
      !$omp end critical
      if (ipnt > nlim) exit
      xpnt = MODULO(ipnt, n_reference) * ldx + 1
      ypnt = ipnt / n_reference * ldx + 1
      wpnt = ldw * omp_get_thread_num() + 1
      call mobbrmsd_run(header, state(ipnt), X(xpnt), Y(ypnt), W(wpnt), &
     &                  cutoff=cutoff, difflim=difflim, maxeval=maxeval, &
     &                  remove_com=remove_com, sort_by_g=sort_by_g &
     &      )
    end do
    !$omp end parallel
!
    if (rotate_y) then
      do concurrent(i=0:n_target - 1)
        ipnt = i * n_reference + 1
        ypnt = i * ldx + 1
        call state(ipnt)%rotation(header, Y(ypnt))
      end do
    end if
  end subroutine mobbrmsd_batch_run
!
  !| batch parallel tri run
  subroutine mobbrmsd_batch_tri_run( &
  &            n_target, header, state, &
  &            X, W, &
  &            cutoff, difflim, maxeval, &
  &            remove_com, sort_by_g, &
  &            n_lower, n_upper)
    integer(IK), intent(in)              :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(*)
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(inout)              :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional        :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional        :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    integer(IK), intent(in), optional    :: n_lower
    !! Specify the lower limit of the range to be calculated. Default [1].
    integer(IK), intent(in), optional    :: n_upper
    !! Specify the upper limit of the range to be calculated. Default [(n_target - 1) * n_target / 2].
    integer(kind=IK)                     :: i, ipnt, ijob, spnt, xpnt, ypnt, wpnt, ldx, ldw, nmin, nlim
    if (n_target < 2) return
    ldx = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
    nmin = 0
    if (PRESENT(n_lower)) nmin = MAX(nmin, n_lower - 1)
    nlim = (n_target - 1) * n_target / 2
    if (PRESENT(n_upper)) nlim = MIN(nlim, MAX(i + 1, n_upper))
    i = nmin
    !$omp parallel private(ipnt, ijob, spnt, xpnt, ypnt, wpnt)
    do
      !$omp critical
      i = i + 1
      ipnt = i
      !$omp end critical
      if (ipnt > nlim) exit
      spnt = ipnt - nmin
      xpnt = INT(0.5_RK * (SQRT(real(8 * ipnt - 7, RK)) - ONE), IK) + 1 ! cantor_pair inverse
      ypnt = ipnt - xpnt * (xpnt - 1) / 2 - 1
      xpnt = xpnt * ldx + 1
      ypnt = ypnt * ldx + 1
      wpnt = ldw * omp_get_thread_num() + 1
      call mobbrmsd_run(header, state(spnt), X(xpnt), X(ypnt), W(wpnt), &
     &                  cutoff=cutoff, difflim=difflim, maxeval=maxeval, &
     &                  remove_com=remove_com, sort_by_g=sort_by_g &
     &      )
    end do
    !$omp end parallel
  end subroutine mobbrmsd_batch_tri_run
end module mod_mobbrmsd_batch_run

