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
  public rmsd
  public autocorr
  public bounds
  public n_eval
  public log_eval_ratio
  public is_finished
  public run
  public restart
  public rotate_y
  public batch_run
  public batch_run_tri
  public min_span_tree

contains

  pure subroutine decode_input(l, seq, mobb)
    integer(kind=IK), intent(in)  :: l
    integer(kind=IK), intent(in)  :: seq(l)
    type(mobbrmsd), intent(inout) :: mobb
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
  end subroutine decode_input

!|  return attributes
!   1 : n_dim
!   2 : n_atom
!   3 : n_header
!   4 : n_int
!   5 : n_float
!   6 : n_rot
!   7 : n_mem
!   8 : n_job
  subroutine decode_attributes(l, seq, att)
    integer(kind=IK), intent(in)  :: l
    integer(kind=IK), intent(in)  :: seq(l)
    integer(kind=IK), intent(out) :: att(8)
    type(mobbrmsd)                :: mobb
    call decode_input(l, seq, mobb)
    call mobbrmsd_attributes( &
   &       mobb, &
   &       n_dim=att(1), &
   &       n_atom=att(2), &
   &       n_mem=att(7), &
   &       n_header=att(3), &
   &       n_int=att(4), &
   &       n_float=att(5) &
   &      )
    att(6) = att(1) * att(1)
    att(8) = 1
    !$omp parallel
    if (omp_get_thread_num() == 0) att(8) = omp_get_num_threads()
    !$omp end parallel
  end subroutine decode_attributes

  pure subroutine decode_header(l, seq, n_header, header)
    integer(kind=IK), intent(in)  :: l
    integer(kind=IK), intent(in)  :: seq(l)
    integer(kind=IK), intent(in)  :: n_header
    integer(kind=IK), intent(out) :: header(n_header)
    type(mobbrmsd)                :: mobb

    call decode_input(l, seq, mobb)
    header = mobbrmsd_dump(mobb)

  end subroutine decode_header

  !| setup dimension
  subroutine setup_dimension_(d)
    integer(kind=IK), intent(in) :: d

    call setup_dimension(d)

  end subroutine setup_dimension_

  !| return rmsd
  pure subroutine rmsd(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res
    res = SQRT(float_states(mobbrmsd_state_RECIPROCAL_OF_N) * &
   &      MAX(ZERO, &
   &          (float_states(mobbrmsd_state_INDEX_TO_AUTOCORR) &
   & + 2._RK * float_states(mobbrmsd_state_INDEX_TO_UPPERBOUND))))
  end subroutine rmsd

  !| return autocorr
  pure subroutine autocorr(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res
    res = float_states(mobbrmsd_state_INDEX_TO_AUTOCORR)
  end subroutine autocorr

  !| return bounds
  pure subroutine bounds(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res(2)
    res(1) = float_states(mobbrmsd_state_INDEX_TO_UPPERBOUND)
    res(2) = float_states(mobbrmsd_state_INDEX_TO_LOWERBOUND)
  end subroutine bounds

  !| return n_eval
  pure subroutine n_eval(n_float, float_states, res)
    integer(kind=IK), intent(in)  :: n_float
    real(kind=RK), intent(in)     :: float_states(n_float)
    integer(kind=IK), intent(out) :: res
    res = NINT(float_states(mobbrmsd_state_INDEX_TO_N_EVAL), IK)
  end subroutine n_eval

  !| return log_eval_ratio
  pure subroutine log_eval_ratio(n_float, float_states, res)
    integer(kind=IK), intent(in) :: n_float
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(out)   :: res
    res = float_states(mobbrmsd_state_INDEX_TO_LOG_RATIO)
  end subroutine log_eval_ratio

  !| inquire bb is finished
  pure subroutine is_finished(n_head, n_int, n_float, header, int_states, float_states, res)
    integer(kind=IK), intent(in)  :: n_head, n_int, n_float
    integer(kind=IK), intent(in)  :: header(n_head)
    integer(kind=IK), intent(in)  :: int_states(n_int)
    real(kind=RK), intent(in)     :: float_states(n_float)
    logical, intent(out)          :: res
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s
    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, float_states)
    res = mobbrmsd_is_finished(h, s)
  end subroutine is_finished

  !| single run with working memory
  pure subroutine run( &
 &                  n_header, &
 &                  n_int, &
 &                  n_float, &
 &                  n_rot, &
 &                  header, &
 &                  X, &
 &                  Y, &
 &                  W, &
 &                  ropts, &
 &                  iopts, &
 &                  remove_com, &
 &                  sort_by_g, &
 &                  difflim_absolute, &
 &                  rotate_y, &
 &                  get_rotation, &
 &                  int_states, &
 &                  float_states, &
 &                  rotation &
 &                 )
    integer(kind=IK), intent(in)  :: n_header
    !! header length
    integer(kind=IK), intent(in)  :: n_int
    integer(kind=IK), intent(in)  :: n_float
    integer(kind=IK), intent(in)  :: n_rot
    integer(kind=IK), intent(in)  :: header(n_header)
    real(kind=RK), intent(in)     :: X(*)
    !! reference coordinate
    real(kind=RK), intent(inout)  :: Y(*)
    !! target coordinate
    real(kind=RK), intent(inout)  :: W(*)
    !! work memory
    real(kind=RK), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=IK), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    logical, intent(in)           :: rotate_y
    logical, intent(in)           :: get_rotation
    integer(kind=IK), intent(out) :: int_states(n_int)
    real(kind=RK), intent(out)    :: float_states(n_float)
    real(kind=RK), intent(out)    :: rotation(n_rot)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s
    call mobbrmsd_load(h, header)
    call mobbrmsd_run( &
      &    h, &
      &    s, &
      &    X, &
      &    Y, &
      &    w, &
      &    cutoff=ropts(1), &
      &    ub_cutoff=ropts(2), &
      &    difflim=ropts(3), &
      &    maxeval=iopts(1),&
      &    remove_com=remove_com,&
      &    sort_by_g=sort_by_g,&
      &    difflim_absolute=difflim_absolute, &
      &    get_rotation=(get_rotation .or. rotate_y) &
      &   )
    if (rotate_y) call mobbrmsd_swap_and_rotation(h, s, Y)

    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)
    if (get_rotation) then
      rotation = mobbrmsd_state_dump_rotation(s)
    else
      rotation = ZERO
    end if
  end subroutine run

  !| restart with working memory
  !pure subroutine restart( &
  subroutine restart( &
 &             n_header, &
 &             n_int, &
 &             n_float, &
 &             n_rot, &
 &             header, &
 &             int_states, &
 &             float_states, &
 &             W, &
 &             rotation, &
 &             ropts, &
 &             iopts, &
 &             difflim_absolute, &
 &             get_rotation &
 &           )
    integer(kind=IK), intent(in)    :: n_header
    integer(kind=IK), intent(in)    :: n_int
    integer(kind=IK), intent(in)    :: n_float
    integer(kind=IK), intent(in)    :: n_rot
    integer(kind=IK), intent(in)    :: header(n_header)
    integer(kind=IK), intent(inout) :: int_states(n_int)
    real(kind=RK), intent(inout)    :: float_states(n_float)
    !! work memory
    real(kind=RK), intent(inout)    :: W(*)
    !! work memory
    real(kind=RK), intent(inout)    :: rotation(n_rot)
    !! rotation matrix
    real(kind=RK), intent(in)       :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=IK), intent(in)    :: iopts(*) ! 1 maxeval
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
    call mobbrmsd_restart( &
      &    h, &
      &    s, &
      &    W, &
      &    cutoff=ropts(1), &
      &    ub_cutoff=ropts(2), &
      &    difflim=ropts(3), &
      &    maxeval=iopts(1), &
      &    difflim_absolute=difflim_absolute &
          )
    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)
    if (get_rotation) rotation = mobbrmsd_state_dump_rotation(s)
  end subroutine restart

  !| compute rotration respect to state.
  pure subroutine rotate_y( &
 &                  n_header, &
 &                  n_int, &
 &                  n_float, &
 &                  n_rot, &
 &                  header, &
 &                  int_states, &
 &                  float_states, &
 &                  rotation, &
 &                  Y &
 &                )
    integer(kind=IK), intent(in) :: n_header
    integer(kind=IK), intent(in) :: n_int
    integer(kind=IK), intent(in) :: n_float
    integer(kind=IK), intent(in) :: n_rot
    integer(kind=IK), intent(in) :: header(n_header)
    integer(kind=IK), intent(in) :: int_states(n_int)
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(in)    :: rotation(n_rot)
    real(kind=RK), intent(inout) :: Y(*)
    type(mobbrmsd)               :: h
    type(mobbrmsd_state)         :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, float_states, rotation)
    call mobbrmsd_swap_and_rotation(h, s, Y)

  end subroutine rotate_y

  !| batch parallel run
  subroutine batch_run( &
 &             n_reference, &
 &             n_target, &
 &             n_chunk, &
 &             n_lower, &
 &             n_header, &
 &             header, &
 &             X, &
 &             Y, &
 &             W, &
 &             ropts, &
 &             iopts, &
 &             remove_com, &
 &             sort_by_g, &
 &             difflim_absolute, &
 &             rmsd &
 &           )
    integer(kind=IK), intent(in)  :: n_reference
    integer(kind=IK), intent(in)  :: n_target
    integer(kind=IK), intent(in)  :: n_chunk
    integer(kind=IK), intent(in)  :: n_lower
    integer(kind=IK), intent(in)  :: n_header
    integer(kind=IK), intent(in)  :: header(n_header)
    real(kind=RK), intent(in)     :: X(*)
   !! reference coordinate
    real(kind=RK), intent(inout)  :: Y(*)
   !! target coordinate
    real(kind=RK), intent(inout)  :: W(*)
   !! work array
    real(kind=RK), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=IK), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    real(kind=RK), intent(out)    :: rmsd(n_chunk)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=IK)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_run( &
   &       n_reference, n_target, h, s, &
   &       X, Y, W, &
   &       cutoff=ropts(1), &
   &       ub_cutoff=ropts(2), &
   &       difflim=ropts(3),&
   &       maxeval=iopts(1), &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       difflim_absolute=difflim_absolute, &
   &       n_lower=n_lower, &
   &       n_upper=n_lower + n_chunk - 1 &
   &     )

    do concurrent(i=1:n_chunk)
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
  end subroutine batch_run

  !| batch parallel tri run
  subroutine batch_run_tri( &
 &    n_target, &
 &    n_chunk, &
 &    n_lower, &
 &    n_header, &
 &    header, &
 &    X, &
 &    W, &
 &    ropts, &
 &    iopts, &
 &    remove_com, &
 &    sort_by_g, &
 &    difflim_absolute, &
 &    rmsd &
 &  )
    integer(kind=IK), intent(in)  :: n_target
    integer(kind=IK), intent(in)  :: n_chunk
    integer(kind=IK), intent(in)  :: n_lower
    integer(kind=IK), intent(in)  :: n_header
    integer(kind=IK), intent(in)  :: header(n_header)
    real(kind=RK), intent(in)     :: X(*)
   !! reference and target coordinate
    real(kind=RK), intent(inout)  :: W(*)
   !! work array
    real(kind=RK), intent(in)     :: ropts(*) ! 1 cutoff 2 ub_cutoff 3 difflim
    integer(kind=IK), intent(in)  :: iopts(*) ! 1 maxeval
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: difflim_absolute
    real(kind=RK), intent(out)    :: rmsd(n_chunk)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=IK)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_tri_run( &
   &       n_target, h, s, X, W, &
   &       cutoff=ropts(1), &
   &       ub_cutoff=ropts(2), &
   &       difflim=ropts(3), &
   &       maxeval=iopts(1), &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       difflim_absolute=difflim_absolute, &
   &       n_lower=n_lower, &
   &       n_upper=n_lower + n_chunk - 1 &
   &     )
    do concurrent(i=1:SIZE(s))
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
  end subroutine batch_run_tri

  subroutine min_span_tree( &
 &    n_target, &
 &    n_header, &
 &    n_int, &
 &    n_float, &
 &    header, &
 &    X, &
 &    W, &
 &    ropts, &
 &    iopts, &
 &    remove_com, &
 &    sort_by_g, &
 &    difflim_absolute, &
 &    edges, &
 &    weights &
 &  )
    integer(kind=IK), intent(in)      :: n_target
    integer(kind=IK), intent(in)      :: n_header
    integer(kind=IK), intent(in)      :: n_int
    integer(kind=IK), intent(in)      :: n_float
    integer(kind=IK), intent(in)      :: header(n_header)
    real(kind=RK), intent(in)         :: X(*)
   !! reference coordinate
    real(kind=RK), intent(inout)      :: W(*)
   !! work memory
    real(kind=RK), intent(in)         :: ropts(*) ! 1 cutoff 2 ub_cutoff 2 difflim
    integer(kind=IK), intent(in)      :: iopts(*) ! 1 maxeval
    logical, intent(in)               :: remove_com
    logical, intent(in)               :: sort_by_g
    logical, intent(in)               :: difflim_absolute
    integer(kind=IK), intent(out)     :: edges(2, n_target - 1)
    real(kind=RK), intent(out)        :: weights(n_target - 1)
    type(mobbrmsd)                    :: h
    type(mobbrmsd_state), allocatable :: s(:, :)

    call mobbrmsd_load(h, header)
    allocate (s(n_target, n_target))

    call mobbrmsd_min_span_tree( &
   &       n_target, h, s, X, W, &
   &       cutoff=ropts(1), &
   &       ub_cutoff=ropts(2), &
   &       difflim=ropts(3), &
   &       maxeval=iopts(1), &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       difflim_absolute=difflim_absolute, &
   &       edges=edges, &
   &       weights=weights &
   &    )

  end subroutine min_span_tree
end module driver

