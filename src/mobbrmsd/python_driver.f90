module driver
!$ use omp_lib
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mobbrmsd
  use mod_mobbrmsd_state
  use mod_mobbrmsd_batch_run
  use mod_mobbrmsd_mst

  implicit none
  private
  public decode_attributes
  public decode_header
  public setup_dimension_
  public run
  public restart
  public permutation_indices
  public rotate_y
  public batch_run
  public batch_run_tri
  public min_span_tree

contains
!   return attributes
!   1 : n_dim
!   2 : n_atom
!   3 : n_header
!   4 : n_int
!   5 : n_float
!   6 : n_rot
!   7 : n_mem
!   8 : n_job
  subroutine decode_attributes(l, seq, att)
    integer(kind=ik), intent(in)  :: l
    integer(kind=ik), intent(in)  :: seq(l)
    integer(kind=ik), intent(out) :: att(8)
    type(mobbrmsd_input)          :: inp
    type(mobbrmsd)                :: mobb
    integer                       :: i, n, m, s, ns

    i = 1
    do while (i < l)
      n = MAX(seq(i), 1)
      m = MAX(seq(i + 1), 1)
      s = MAX(seq(i + 2), 1)
      ns = n * (s - 1)
      call mobbrmsd_input_add_molecule(inp, n, m, RESHAPE(seq(i + 3:i + 2 + ns), [n, s - 1]))
      i = i + 3 + ns
    end do
    call mobbrmsd_init(mobb, inp)

    call mobbrmsd_attributes(mobb &
                          &, n_dim=att(1) &
                          &, n_atom=att(2) &
                          &, n_mem=att(7) &
                          &, n_header=att(3) &
                          &, n_int=att(4) &
                          &, n_float=att(5) &
                          & )
    att(6) = att(1) * att(1)
    att(8) = 1
    !$omp parallel
    if (omp_get_thread_num() == 0) att(8) = omp_get_num_threads()
    !$omp end parallel
  end subroutine decode_attributes

  pure subroutine decode_header(l, seq, n_header, header)
    integer(kind=ik), intent(in)  :: l
    integer(kind=ik), intent(in)  :: seq(l)
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(out) :: header(n_header)
    type(mobbrmsd)                :: mobb
    type(mobbrmsd_input)          :: inp
    integer                       :: i, n, m, s, ns

    i = 1
    do while (i < l)
      n = MAX(seq(i), 1)
      m = MAX(seq(i + 1), 1)
      s = MAX(seq(i + 2), 1)
      ns = n * (s - 1)
      call mobbrmsd_input_add_molecule(inp, n, m, RESHAPE(seq(i + 3:i + 2 + ns), [n, s - 1]))
      i = i + 3 + ns
    end do
    call mobbrmsd_init(mobb, inp)

    header = mobbrmsd_dump(mobb)

  end subroutine decode_header

  subroutine setup_dimension_(d)
    integer(kind=ik), intent(in) :: d

    call setup_dimension(d)

  end subroutine setup_dimension_

  pure subroutine run(n_header &
                   &, n_int &
                   &, n_float &
                   &, n_rot &
                   &, header &
                   &, X &
                   &, Y &
                   &, W &
                   &, ropts &
                   &, iopts &
                   &, remove_com &
                   &, sort_by_g &
                   &, difflim_absolute &
                   &, rotate_y &
                   &, get_rotation &
                   &, int_states &
                   &, float_states &
                   &, rotation &
                   & )
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: n_int
    integer(kind=ik), intent(in)  :: n_float
    integer(kind=ik), intent(in)  :: n_rot
    integer(kind=ik), intent(in)  :: header(n_header)
    real(kind=rk), intent(in)     :: X(*)
    real(kind=rk), intent(inout)  :: Y(*)
    real(kind=rk), intent(inout)  :: W(*)
    real(kind=rk), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=ik), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    logical, intent(in)           :: rotate_y
    logical, intent(in)           :: get_rotation
    integer(kind=ik), intent(out) :: int_states(n_int)
    real(kind=rk), intent(out)    :: float_states(n_float)
    real(kind=rk), intent(out)    :: rotation(n_rot)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s
    call mobbrmsd_load(h, header)
    call mobbrmsd_run(h &
                   &, s &
                   &, X &
                   &, Y &
                   &, w &
                   &, cutoff=ropts(1) &
                   &, ub_cutoff=ropts(2) &
                   &, difflim=ropts(3) &
                   &, maxeval=iopts(1) &
                   &, remove_com=remove_com &
                   &, sort_by_g=sort_by_g &
                   &, difflim_absolute=difflim_absolute &
                   &, get_rotation=(get_rotation .or. rotate_y) &
                   & )
    if (rotate_y) call mobbrmsd_swap_and_rotation(h, s, Y)

    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)
    if (get_rotation) then
      rotation = mobbrmsd_state_dump_rotation(s)
    else
      rotation = ZERO
    end if
  end subroutine run

  pure subroutine restart(n_header &
                       &, n_int &
                       &, n_float &
                       &, n_rot &
                       &, header &
                       &, int_states &
                       &, float_states &
                       &, W &
                       &, rotation &
                       &, ropts &
                       &, iopts &
                       &, difflim_absolute &
                       &, get_rotation &
                       & )
    integer(kind=ik), intent(in)    :: n_header
    integer(kind=ik), intent(in)    :: n_int
    integer(kind=ik), intent(in)    :: n_float
    integer(kind=ik), intent(in)    :: n_rot
    integer(kind=ik), intent(in)    :: header(n_header)
    integer(kind=ik), intent(inout) :: int_states(n_int)
    real(kind=rk), intent(inout)    :: float_states(n_float)
    real(kind=rk), intent(inout)    :: W(*)
    real(kind=rk), intent(inout)    :: rotation(n_rot)
    real(kind=rk), intent(in)       :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=ik), intent(in)    :: iopts(*) ! 1 maxeval
    logical, intent(in)             :: difflim_absolute
    logical, intent(in)             :: get_rotation
    type(mobbrmsd)                  :: h
    type(mobbrmsd_state)            :: s
    call mobbrmsd_load(h, header)
    if (get_rotation) then
      call mobbrmsd_state_load(s, int_states, float_states, rotation)
    else
      call mobbrmsd_state_load(s, int_states, float_states)
    end if
    call mobbrmsd_restart(h &
                       &, s &
                       &, W &
                       &, cutoff=ropts(1) &
                       &, ub_cutoff=ropts(2) &
                       &, difflim=ropts(3) &
                       &, maxeval=iopts(1) &
                       &, difflim_absolute=difflim_absolute &
                       & )
    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)
    if (get_rotation) rotation = mobbrmsd_state_dump_rotation(s)
  end subroutine restart

  pure subroutine rotate_y(n_header &
                        &, n_int &
                        &, n_float &
                        &, n_rot &
                        &, header &
                        &, int_states &
                        &, float_states &
                        &, rotation &
                        &, Y &
                        & )
    integer(kind=ik), intent(in) :: n_header
    integer(kind=ik), intent(in) :: n_int
    integer(kind=ik), intent(in) :: n_float
    integer(kind=ik), intent(in) :: n_rot
    integer(kind=ik), intent(in) :: header(n_header)
    integer(kind=ik), intent(in) :: int_states(n_int)
    real(kind=rk), intent(in)    :: float_states(n_float)
    real(kind=rk), intent(in)    :: rotation(n_rot)
    real(kind=rk), intent(inout) :: Y(*)
    type(mobbrmsd)               :: h
    type(mobbrmsd_state)         :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, z=float_states, r=rotation)
    call mobbrmsd_swap_and_rotation(h, s, Y)

  end subroutine rotate_y

  pure subroutine permutation_indices(n_header &
                                   &, n_int &
                                   &, n_float &
                                   &, header &
                                   &, int_states &
                                   &, float_states &
                                   &, IX &
                                   & )
    integer(kind=ik), intent(in)    :: n_header
    integer(kind=ik), intent(in)    :: n_int
    integer(kind=ik), intent(in)    :: n_float
    integer(kind=ik), intent(in)    :: header(n_header)
    integer(kind=ik), intent(in)    :: int_states(n_int)
    real(kind=rk), intent(in)       :: float_states(n_float)
    integer(kind=ik), intent(inout) :: IX(*)
    type(mobbrmsd)                  :: h
    type(mobbrmsd_state)            :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, z=float_states)
    call mobbrmsd_swap_indices(h, s, IX, base=0_IK)

  end subroutine permutation_indices

  subroutine batch_run(n_reference &
                    &, n_target &
                    &, n_chunk &
                    &, n_lower &
                    &, n_header &
                    &, header &
                    &, X &
                    &, Y &
                    &, W &
                    &, ropts &
                    &, iopts &
                    &, remove_com &
                    &, sort_by_g &
                    &, difflim_absolute &
                    &, rmsd &
                    &, log_n_eval &
                    & )
    integer(kind=ik), intent(in)  :: n_reference
    integer(kind=ik), intent(in)  :: n_target
    integer(kind=ik), intent(in)  :: n_chunk
    integer(kind=ik), intent(in)  :: n_lower
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: header(n_header)
    real(kind=rk), intent(in)     :: X(*)
    real(kind=rk), intent(inout)  :: Y(*)
    real(kind=rk), intent(inout)  :: W(*)
    real(kind=rk), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=ik), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    real(kind=rk), intent(out)    :: rmsd(n_chunk)
    real(kind=rk), intent(out)    :: log_n_eval
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=ik)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_run(n_reference &
                         &, n_target &
                         &, h &
                         &, s &
                         &, X &
                         &, Y &
                         &, W &
                         &, cutoff=ropts(1) &
                         &, ub_cutoff=ropts(2) &
                         &, difflim=ropts(3)&
                         &, maxeval=iopts(1) &
                         &, remove_com=remove_com &
                         &, sort_by_g=sort_by_g &
                         &, difflim_absolute=difflim_absolute &
                         &, n_lower=n_lower &
                         &, n_upper=n_lower + n_chunk - 1 &
                         &     )

    do concurrent(i=1:n_chunk)
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
    log_n_eval = mobbrmsd_state_log_sum_n_eval(n_chunk, s)
  end subroutine batch_run

  subroutine batch_run_tri(n_target &
                        &, n_chunk &
                        &, n_lower &
                        &, n_header &
                        &, header &
                        &, X &
                        &, W &
                        &, ropts &
                        &, iopts &
                        &, remove_com &
                        &, sort_by_g &
                        &, difflim_absolute &
                        &, rmsd &
                        &, log_n_eval &
                        & )
    integer(kind=ik), intent(in)  :: n_target
    integer(kind=ik), intent(in)  :: n_chunk
    integer(kind=ik), intent(in)  :: n_lower
    integer(kind=ik), intent(in)  :: n_header
    integer(kind=ik), intent(in)  :: header(n_header)
    real(kind=rk), intent(in)     :: X(*)
    real(kind=rk), intent(inout)  :: W(*)
    real(kind=rk), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=ik), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    real(kind=rk), intent(out)    :: rmsd(n_chunk)
    real(kind=rk), intent(out)    :: log_n_eval
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=ik)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_tri_run(n_target &
                             &, h &
                             &, s &
                             &, X &
                             &, W &
                             &, cutoff=ropts(1) &
                             &, ub_cutoff=ropts(2) &
                             &, difflim=ropts(3) &
                             &, maxeval=iopts(1) &
                             &, remove_com=remove_com &
                             &, sort_by_g=sort_by_g &
                             &, difflim_absolute=difflim_absolute &
                             &, n_lower=n_lower &
                             &, n_upper=n_lower + n_chunk - 1 &
                             & )
    do concurrent(i=1:SIZE(s))
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
    log_n_eval = mobbrmsd_state_log_sum_n_eval(n_chunk, s)
  end subroutine batch_run_tri

  subroutine min_span_tree(n_target &
                        &, n_header &
                        &, header &
                        &, X &
                        &, ropts &
                        &, iopts &
                        &, remove_com &
                        &, sort_by_g &
                        &, verbose &
                        &, edges &
                        &, weights &
                        &, log_n_eval &
                        & )
    integer(kind=ik), intent(in)      :: n_target
    integer(kind=ik), intent(in)      :: n_header
    integer(kind=ik), intent(in)      :: header(n_header)
    real(kind=rk), intent(in)         :: X(*)
    real(kind=rk), intent(in)         :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=ik), intent(in)      :: iopts(*) ! 1 n_work
    logical, intent(in)               :: remove_com
    logical, intent(in)               :: sort_by_g
    logical, intent(in)               :: verbose
    integer(kind=ik), intent(out)     :: edges(2, n_target - 1)
    real(kind=rk), intent(out)        :: weights(n_target - 1)
    real(kind=rk), intent(out)        :: log_n_eval
    type(mobbrmsd)                    :: h

    call mobbrmsd_load(h, header)
    call mobbrmsd_min_span_tree(n_target &
                             &, h &
                             &, X &
                             &, n_work=iopts(1) &
                             &, remove_com=remove_com &
                             &, sort_by_g=sort_by_g &
                             &, verbose=verbose &
                             &, edges=edges &
                             &, weights=weights &
                             &, log_n_eval=log_n_eval &
                             &)

  end subroutine min_span_tree
end module driver

