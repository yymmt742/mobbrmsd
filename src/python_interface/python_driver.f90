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

  subroutine decode_attributes(l, seq, &
 &                  n_dim, n_atom, &
 &                  n_mem, n_job, &
 &                  n_header, n_int, n_float)
    integer(kind=IK), intent(in)  :: l
    integer(kind=IK), intent(in)  :: seq(l)
    integer(kind=IK), intent(out) :: n_dim
    integer(kind=IK), intent(out) :: n_atom
    integer(kind=IK), intent(out) :: n_mem
    integer(kind=IK), intent(out) :: n_job
    integer(kind=IK), intent(out) :: n_header
    integer(kind=IK), intent(out) :: n_int
    integer(kind=IK), intent(out) :: n_float
    type(mobbrmsd)                :: mobb
    call decode_input(l, seq, mobb)
    call mobbrmsd_attributes( &
   &       mobb, &
   &       n_dim=n_dim, &
   &       n_atom=n_atom, &
   &       n_mem=n_mem, &
   &       n_header=n_header, &
   &       n_int=n_int, &
   &       n_float=n_float &
   &      )
    !$omp parallel
    if (omp_get_thread_num() == 0) n_job = omp_get_num_threads()
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
 &                  header, &
 &                  X, &
 &                  Y, &
 &                  W, &
 &                  cutoff, &
 &                  difflim, &
 &                  maxeval, &
 &                  remove_com, &
 &                  sort_by_g, &
 &                  rotate_y, &
 &                  int_states, &
 &                  float_states &
 &                 )
    integer(kind=IK), intent(in)  :: n_header
    !! header length
    integer(kind=IK), intent(in)  :: n_int
    integer(kind=IK), intent(in)  :: n_float
    integer(kind=IK), intent(in)  :: header(n_header)
    real(kind=RK), intent(in)     :: X(*)
    !! reference coordinate
    real(kind=RK), intent(inout)  :: Y(*)
    !! target coordinate
    real(kind=RK), intent(inout)  :: W(*)
    !! work memory
    integer(kind=IK), intent(in)  :: maxeval
    real(kind=RK), intent(in)     :: cutoff
    real(kind=RK), intent(in)     :: difflim
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    logical, intent(in)           :: rotate_y
    integer(kind=IK), intent(out) :: int_states(n_int)
    real(kind=RK), intent(out)    :: float_states(n_float)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_run(h, s, X, Y, w, &
      &               cutoff=cutoff, &
      &               difflim=difflim, &
      &               maxeval=maxeval,&
      &               remove_com=remove_com,&
      &               sort_by_g=sort_by_g&
      &   )
    if (rotate_y) call mobbrmsd_swap_and_rotation(h, s, Y)

    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)

  end subroutine run

  !| restart with working memory
  pure subroutine restart( &
 &    n_header, n_int, n_float, &
 &    header, &
 &    int_states, float_states, &
 &    W, &
 &    cutoff, difflim, maxeval &
 &  )
    integer(kind=IK), intent(in)    :: n_header
    integer(kind=IK), intent(in)    :: n_int
    integer(kind=IK), intent(in)    :: n_float
    integer(kind=IK), intent(in)    :: header(n_header)
    integer(kind=IK), intent(inout) :: int_states(n_int)
    real(kind=RK), intent(inout)    :: float_states(n_float)
    !! work memory
    real(kind=RK), intent(inout)    :: W(*)
    !! work memory
    integer(kind=IK), intent(in)    :: maxeval
    real(kind=RK), intent(in)       :: cutoff
    real(kind=RK), intent(in)       :: difflim
    type(mobbrmsd)                  :: h
    type(mobbrmsd_state)            :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, float_states)
    call mobbrmsd_restart(h, s, W, cutoff, difflim, maxeval)

    int_states = mobbrmsd_state_dump(s)
    float_states = mobbrmsd_state_dump_real(s)

  end subroutine restart

  !| compute rotration respect to state.
  pure subroutine rotate_y( &
 &    n_header, n_int, n_float, &
 &    header, int_states, float_states, Y &
 &  )
    integer(kind=IK), intent(in) :: n_header
    integer(kind=IK), intent(in) :: n_int
    integer(kind=IK), intent(in) :: n_float
    integer(kind=IK), intent(in) :: header(n_header)
    integer(kind=IK), intent(in) :: int_states(n_int)
    real(kind=RK), intent(in)    :: float_states(n_float)
    real(kind=RK), intent(inout) :: Y(*)
    type(mobbrmsd)               :: h
    type(mobbrmsd_state)         :: s

    call mobbrmsd_load(h, header)
    call mobbrmsd_state_load(s, int_states, float_states)
    call mobbrmsd_swap_and_rotation(h, s, Y)

  end subroutine rotate_y

  !| batch parallel run
  subroutine batch_run( &
 &    n_reference, n_target, &
 &    n_chunk, n_lower, &
 &    n_header, header, &
 &    X, Y, W, &
 &    cutoff, difflim, maxeval, &
 &    remove_com, sort_by_g, &
 &    rmsd &
 &  )
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
    integer(kind=IK), intent(in)  :: maxeval
    real(kind=RK), intent(in)     :: cutoff
    real(kind=RK), intent(in)     :: difflim
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    real(kind=RK), intent(out)    :: rmsd(n_chunk)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=IK)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_run( &
   &       n_reference, n_target, h, s, &
   &       X, Y, W, &
   &       cutoff=cutoff, &
   &       difflim=difflim,&
   &       maxeval=maxeval, &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       n_lower=n_lower, &
   &       n_upper=n_lower + n_chunk - 1 &
   &     )

    do concurrent(i=1:n_chunk)
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
  end subroutine batch_run

  !| batch parallel tri run
  subroutine batch_run_tri( &
 &    n_target, n_chunk, n_lower, &
 &    n_header, header, &
 &    X, W, &
 &    cutoff, difflim, maxeval, &
 &    remove_com, sort_by_g, &
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
    integer(kind=IK), intent(in)  :: maxeval
    real(kind=RK), intent(in)     :: cutoff
    real(kind=RK), intent(in)     :: difflim
    logical, intent(in)           :: remove_com
    logical, intent(in)           :: sort_by_g
    real(kind=RK), intent(out)    :: rmsd(n_chunk)
    type(mobbrmsd)                :: h
    type(mobbrmsd_state)          :: s(n_chunk)
    integer(kind=IK)              :: i
    call mobbrmsd_load(h, header)
    call mobbrmsd_batch_tri_run( &
   &       n_target, h, s, X, W, &
   &       cutoff=cutoff, &
   &       difflim=difflim, &
   &       maxeval=maxeval, &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       n_lower=n_lower, &
   &       n_upper=n_lower + n_chunk - 1 &
   &     )
    do concurrent(i=1:SIZE(s))
      rmsd(i) = mobbrmsd_state_rmsd(s(i))
    end do
  end subroutine batch_run_tri

  subroutine min_span_tree( &
 &    n_target, n_header, &
 &    n_int, n_float, &
 &    header, &
 &    X, W, cutoff, difflim, maxeval,  &
 &    remove_com, sort_by_g, &
 &    edges, weights &
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
    real(kind=RK), intent(in)         :: cutoff
    real(kind=RK), intent(in)         :: difflim
    integer(kind=IK), intent(in)      :: maxeval
    logical, intent(in)               :: remove_com
    logical, intent(in)               :: sort_by_g
    integer(kind=IK), intent(out)     :: edges(2, n_target - 1)
    real(kind=RK), intent(out)        :: weights(n_target - 1)
    type(mobbrmsd)                    :: h
    type(mobbrmsd_state), allocatable :: s(:, :)

    call mobbrmsd_load(h, header)
    allocate (s(n_target, n_target))

    call mobbrmsd_min_span_tree( &
   &       n_target, h, s, X, W, &
   &       cutoff=cutoff, &
   &       difflim=difflim, &
   &       maxeval=maxeval, &
   &       remove_com=remove_com, &
   &       sort_by_g=sort_by_g, &
   &       edges=edges, &
   &       weights=weights &
   &    )

  end subroutine min_span_tree
end module driver

