program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO, RHUGE
  use driver
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test python_driver')
  call test0()
  call test1()
!
  call u%init('test python_driver, min_span_tree')
  call test2()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: n_target = 20
    integer(IK), parameter   :: n_apm = 8, n_mol = 8
    integer(IK), parameter   :: seq(3 + n_apm) = [n_apm, n_mol, s, 5, 6, 7, 8, 1, 2, 3, 4]
    real(RK), allocatable    :: X(:, :, :, :)
    integer(IK)              :: att(8)
    integer(IK)              :: n_dim, n_atm, n_header, n_int, n_float, n_rot, n_mem, n_job
    integer(IK), allocatable :: h(:)
    real(RK), allocatable    :: rmsd(:), W(:, :)
    real(RK)                 :: log_n_eval
    integer(IK)              :: i, j, k, l
!
    call decode_attributes(SIZE(seq), seq, att)
    print'(*(I6))', att
    n_dim = att(1)
    n_atm = att(2)
    n_header = att(3)
    n_int = att(4)
    n_float = att(5)
    n_rot = att(6)
    n_mem = att(7)
    n_job = att(8)

    allocate (h(n_header))
    allocate (rmsd((n_target - 1) * n_target / 2))
    allocate (W(n_mem, n_job))
    allocate (X(n_dim, n_apm, n_mol, n_target))
!
    do k = -1, 1, 2
      do j = -1, 1, 2
        do i = -1, 1, 2
          X(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, 1) = sample(n_dim, n_apm, [i, j, k])
        end do
      end do
    end do
    do l = 2, n_target
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            X(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l) = 0.9_RK * X(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l - 1)&
           &                                                    + 0.1_RK * sample(n_dim, n_apm, [i, j, k])
          end do
        end do
      end do
    end do
!
    call decode_header(SIZE(seq), seq, n_header, h)
    call batch_run_tri( &
   &  n_target, n_target * (n_target - 1) / 2, 1, &
   &  n_header, h, &
   &  X, W, &
   &  [999.0_RK, 999.0_RK, 0.0_RK], [-1], .true., .true., .false., &
   &  rmsd, log_n_eval)
!
    k = 0
    do j = 1, n_target
      do i = 1, j - 1
        k = k + 1
        print'(2I4,F16.9)', i, j, rmsd(k)
        !print'(F9.6,F9.1,F16.9)', sr(1, k) * (sr(2, k) + 2 * sr(3, k)), sr(5, k), EXP(sr(6, k))
      end do
      print *
    end do
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test0
!
  subroutine test1()
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: n_reference = 2
    integer(IK), parameter   :: n_target = 300
    integer(IK), parameter   :: n_apm = 8, n_mol = 8
    integer(IK), parameter   :: seq(3 + n_apm) = [n_apm, n_mol, s, 5, 6, 7, 8, 1, 2, 3, 4]
    integer(IK)              :: att(8)
    integer(IK)              :: n_dim, n_atm, n_header, n_int, n_float, n_rot, n_mem, n_job
    real(RK), allocatable    :: X(:, :, :, :), Y(:, :, :, :)
    integer(IK), allocatable :: h(:)
    real(RK), allocatable    :: W(:, :), rmsd(:, :)
    real(RK)                 :: log_n_eval
    integer(IK)              :: i, j, k, l
!
    call decode_attributes(SIZE(seq), seq, att)
    print'(*(I6))', att
    n_dim = att(1)
    n_atm = att(2)
    n_header = att(3)
    n_int = att(4)
    n_float = att(5)
    n_rot = att(6)
    n_mem = att(7)
    n_job = att(8)

    allocate (h(n_header))
    allocate (rmsd(n_reference, n_target))
    allocate (W(n_mem, n_job))
    allocate (X(n_dim, n_apm, n_mol, n_reference))
    allocate (Y(n_dim, n_apm, n_mol, n_target))
!
    do l = 1, n_reference
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            X(:, :, (i + 1) * 2 + (j + 1) + (k + 1) / 2 + 1, l) = sample(n_dim, n_apm, [i, j, k])
          end do
        end do
      end do
    end do
!
    do l = 1, n_target
      do k = -1, 1, 2
        do j = -1, 1, 2
          do i = -1, 1, 2
            Y(:, :, (k + 1) * 2 + (i + 1) + (j + 1) / 2 + 1, l) = sample(n_dim, n_apm, [i, j, k])
          end do
        end do
      end do
      Y(:, :, :, l) = 0.2 * X(:, :, :, 1) + 0.8 * Y(:, :, :, l)
    end do
!
    call decode_header(SIZE(seq), seq, n_header, h)
!
    call batch_run( &
   &  n_reference, n_target, &
   &  n_reference * n_target, &
   &  1, n_header, h, &
   &  X, Y, W, &
   &  [999.0_RK, 999.0_RK, 0.0_RK], [-1], .true., .true., .false., &
   &  rmsd, log_n_eval)
!
    do j = 1, n_target
      do i = 1, n_reference
        print'(2I4,F16.9)', i, j, rmsd(i, j)
      end do
      print *
    end do
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test1
!
  subroutine test2()
    integer(IK), parameter   :: n_apm = 8, n_mol = 5, n_sym = 2, n_target = 30
    integer(IK), parameter   :: seq(3 + n_apm) = [n_apm, n_mol, n_sym, 5, 6, 7, 8, 1, 2, 3, 4]
    real(RK), allocatable    :: X(:, :, :, :)
    integer(IK)              :: att(8)
    integer(IK)              :: n_dim, n_atm, n_header, n_int, n_float, n_rot, n_mem, n_job
    integer(IK)              :: edges(2, n_target - 1)
    real(RK)                 :: weights(n_target - 1), log_n_eval
    integer(IK), allocatable :: header(:), int_states(:, :, :)
    real(RK), allocatable    :: float_states(:, :, :)
    integer(IK)              :: i
!
    call decode_attributes(SIZE(seq), seq, att)
    n_dim = att(1)
    n_atm = att(2)
    n_header = att(3)
    n_int = att(4)
    n_float = att(5)
    n_rot = att(6)
    n_mem = att(7)
    n_job = att(8)
!
    print'(*(I6))', n_dim, n_atm, n_header, n_int, n_float, n_rot, n_mem, n_job

    allocate (X(n_dim, n_apm, n_mol, n_target))
    allocate (header(n_header))
    allocate (int_states(n_int, n_target, n_target))
    allocate (float_states(n_float, n_target, n_target))
!
    call RANDOM_NUMBER(X)
    call decode_header(SIZE(seq), seq, n_header, header)
!
    call min_span_tree(n_target &
                    &, n_header &
                    &, header &
                    &, X &
                    &, ropts=[0.0_RK] &
                    &, iopts=[0_IK] &
                    &, remove_com=.true. &
                    &, sort_by_g=.true. &
                    &, verbose=.true. &
                    &, edges=edges &
                    &, weights=weights &
                    &, log_n_eval=log_n_eval &
                    & )
!
    do i = 1, n_target - 1
      print *, "# ", edges(:, i), weights(i)
    end do
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
!
  end subroutine test2
!
  function sample(d, n, com) result(res)
    integer, intent(in)  :: d, n, com(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) + 8 * com
    end do
  end function sample
!
end program main

