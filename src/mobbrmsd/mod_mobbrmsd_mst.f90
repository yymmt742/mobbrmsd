!| Configure a minimum global tree with mobbrmsd. <br>
!  MST construction is based on the Prims algorithm,
!  but the useless calculations are reduced using the cutoff possibilities of mobbrmsd.
module mod_mobbrmsd_mst
!$ use omp_lib
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_mobbrmsd_state
  use mod_mobbrmsd
  implicit none
  private
  public :: mobbrmsd_min_span_tree
!
contains
!| minimum spanning tree construction.
  subroutine mobbrmsd_min_span_tree( &
 &             n_target, &
 &             header, &
 &             state, &
 &             X, &
 &             W, &
 &             cutoff, &
 &             ub_cutoff, &
 &             difflim, &
 &             maxeval, &
 &             remove_com, &
 &             sort_by_g, &
 &             difflim_absolute, &
 &             edges, &
 &             weights &
 &            )
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd), intent(in)          :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target * (n_target - 1) / 2)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)           :: X(*)
    !! coordinate sequence
    real(kind=RK), intent(inout)        :: W(*)
!   !! work memory, must be larger than header%memsize() * n_target * (n_target-1) / 2
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: ub_cutoff
    !! The search ends when upperbound is determined to be greater than to ub_cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    logical, intent(in), optional       :: difflim_absolute
    !! if true, use absolute difflim. default [.false.]
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    integer(kind=IK)                    :: f(n_target), pi, pj
    real(RK)                            :: ub, rmsd, cutoff_global
    integer(kind=IK)                    :: i, j, k, ilcl, jlcl, nlim, spnt, xpnt, ypnt, wpnt, l0, ldx, ldw
!
    l0 = n_target * (n_target - 1) / 2
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    ldw = mobbrmsd_memsize(header)
!
    !$omp parallel do private(spnt, xpnt, ypnt, wpnt)
    do j = 1, n_target
      do i = 1, j - 1
        call cantor_pair(i, j, spnt)
        xpnt = (i - 1) * ldx + 1
        ypnt = (j - 1) * ldx + 1
        wpnt = (spnt - 1) * ldw + 1
        call mobbrmsd_run( &
       &       header, &
       &       state(spnt), &
       &       X(xpnt), &
       &       X(ypnt), &
       &       W(wpnt), &
       &       cutoff=cutoff, &
       &       ub_cutoff=ub_cutoff, &
       &       difflim=difflim, &
       &       maxeval=1, &
       &       remove_com=remove_com, &
       &       sort_by_g=sort_by_g, &
       &       difflim_absolute=difflim_absolute  &
       &      )
      end do
    end do
    !$omp end parallel do
!
    if (PRESENT(cutoff)) then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    end if
!
    do concurrent(i=1:n_target)
      f(i) = i
    end do
    k = MINLOC(mobbrmsd_state_lowerbound_as_rmsd(state(:l0)), 1)
    call cantor_pair_inverse(k, i, j)
    pi = 0
    pj = 0
    f(1) = f(i)
    f(i) = 1
!
    do k = 1, n_target - 1
      ub = cutoff_global
      do i = 1, k
        do j = k + 1, n_target
          call cantor_pair(f(i), f(j), spnt)
          ub = MIN(mobbrmsd_state_rmsd(state(spnt)), ub)
        end do
      end do
      nlim = k * (n_target - k)
      i = 0
      !$omp parallel private(ilcl, jlcl, spnt, xpnt, ypnt, wpnt)
      do
        !$omp critical
        ilcl = i
        i = i + 1
        !$omp end critical
        if (ilcl >= nlim) exit
        jlcl = 1 + MODULO(ilcl, k)
        ilcl = k + 1 + ilcl / k
        call cantor_pair(f(ilcl), f(jlcl), spnt)
        if (.not. mobbrmsd_state_is_finished(state(spnt))) then
          wpnt = (spnt - 1) * ldw + 1
          call mobbrmsd_restart( &
         &       header, &
         &       state(spnt), &
         &       W(wpnt), &
         &       cutoff=ub, &
         &       ub_cutoff=ub_cutoff, &
         &       difflim=difflim, &
         &       maxeval=maxeval &
         &      )
        end if
        !$omp critical
        ub = MIN(ub, mobbrmsd_state_rmsd(state(spnt)))
        !$omp end critical
      end do
      !$omp end parallel
!
      ub = RHUGE
      do i = 1, k
        do j = k + 1, n_target
          call cantor_pair(f(i), f(j), spnt)
          rmsd = mobbrmsd_state_rmsd(state(spnt))
          if (rmsd < ub) then
            ub = rmsd
            pi = i
            pj = j
          end if
        end do
      end do
!
      if (PRESENT(edges)) then
        if (f(pi) < f(pj)) then
          edges(1, k) = f(pi)
          edges(2, k) = f(pj)
        else
          edges(1, k) = f(pj)
          edges(2, k) = f(pi)
        end if
      end if
      if (PRESENT(weights)) then
        weights(k) = ub
      end if
      j = f(k + 1)
      f(k + 1) = f(pj)
      f(pj) = j
    end do
!
  end subroutine mobbrmsd_min_span_tree
!
! ---
!
  pure elemental subroutine cantor_pair(i, j, k)
    integer(IK), intent(in)    :: i, j
    integer(IK), intent(inout) :: k
    if (i > j) then
      k = (i - 2) * (i - 1) / 2 + j
    else
      k = (j - 2) * (j - 1) / 2 + i
    end if
  end subroutine cantor_pair
!
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = INT((SQRT(real(8 * k - 7, RK)) - ONE), IK) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
end module mod_mobbrmsd_mst

