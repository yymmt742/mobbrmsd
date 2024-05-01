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
  use mod_forbar
  use mod_forbar_collections
  implicit none
  public :: mobbrmsd_batch_run
!
contains
!
  !| batch parallel run
  subroutine mobbrmsd_batch_run(n_target, header, state, &
  &                             X, Y, W, &
  &                             cutoff, difflim, maxeval, &
  &                             remove_com, sort_by_g, &
  &                             rotate_y)
    integer(IK), intent(in)              :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(n_target)
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
    integer(kind=IK)                     :: i, itgt, ijob, ypnt, wpnt, ldy, ldw
!
    i = 0
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    !$omp parallel private(itgt, ijob, ypnt, wpnt)
    do
      !$omp critical
      i = i + 1
      itgt = i
      !$omp end critical
      if (itgt > n_target) exit
      wpnt = ldw * omp_get_thread_num() + 1
      ypnt = (itgt - 1) * ldy + 1
      call mobbrmsd_run(header, state(itgt), X, Y(ypnt), W(wpnt), &
     &                  cutoff=cutoff, difflim=difflim, maxeval=maxeval, &
     &                  remove_com=remove_com, sort_by_g=sort_by_g &
     &      )
      if (rotate_y) call state(itgt)%rotation(header, Y(ypnt))
    end do
    !$omp end parallel
!
  end subroutine mobbrmsd_batch_run
!
end module mod_mobbrmsd_batch_run

