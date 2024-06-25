program main
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
    integer(IK)              :: n_atm
    integer(IK)              :: n_dim
    real(RK), allocatable    :: X(:, :, :, :)
    integer(IK)              :: n_header, n_int, n_float, n_job, n_mem
    integer(IK), allocatable :: h(:), si(:, :)
    real(RK), allocatable    :: sr(:, :), W(:, :)
    integer(IK)              :: i, j, k, l
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
    call n_atoms(n_dim, n_atm)
!
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
    call workmemory_lengthes(n_mem, n_job)
    call state_vector_lengthes(n_header, n_int, n_float)
    print *, n_dim, n_atm
    print *, n_mem, n_job, n_header, n_int, n_float
!
    allocate (h(n_header))
    allocate (si(n_int, (n_target - 1) * n_target / 2))
    allocate (sr(n_float, (n_target - 1) * n_target / 2))
    allocate (W(n_mem, n_job))
!
    call batch_run_tri( &
   &  n_target, n_header, n_int, n_float, &
   &  n_target * (n_target - 1) / 2, 1, &
   &  X, W, &
   &  999.0_RK, 0.0_RK, -1, .true., .true., &
   &  h, si, sr)
!
    k = 0
    do j = 1, n_target
      do i = 1, j - 1
        k = k + 1
        print'(F9.2,F9.1,F16.9)', sr(2, k) + 2 * sr(3, k), sr(5, k), EXP(sr(6, k))
        !print'(F9.6,F9.1,F16.9)', sr(1, k) * (sr(2, k) + 2 * sr(3, k)), sr(5, k), EXP(sr(6, k))
      end do
      print *
    end do
    call clear_molecule()
!
  end subroutine test0
!
  subroutine test1()
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: n_reference = 2
    integer(IK), parameter   :: n_target = 1000
    integer(IK), parameter   :: n_apm = 8, n_mol = 8
    integer(IK)              :: n_atm
    integer(IK)              :: n_dim
    real(RK), allocatable    :: X(:, :, :, :), Y(:, :, :, :)
    integer(IK)              :: n_header, n_int, n_float, n_job, n_mem
    integer(IK), allocatable :: h(:), si(:, :, :)
    real(RK), allocatable    :: sr(:, :, :), W(:, :)
    integer(IK)              :: i, j, k, l
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
    call n_atoms(n_dim, n_atm)
!
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
    call workmemory_lengthes(n_mem, n_job)
    call state_vector_lengthes(n_header, n_int, n_float)
    print *, n_dim, n_atm
    print *, n_mem, n_job, n_header, n_int, n_float
!
    allocate (h(n_header))
    allocate (si(n_int, n_reference, n_target))
    allocate (sr(n_float, n_reference, n_target))
    allocate (W(n_mem, n_job))
!
    call batch_run( &
   &  n_reference, n_target, n_header, n_int, n_float, &
   &  n_reference * n_target, 1, &
   &  X, Y, W, &
   &  999.0_RK, 0.0_RK, -1, .true., .true., &
   &  h, si, sr)
!
    do j = 1, n_target
      do i = 1, n_reference
        print'(F9.6,4(F9.1),F16.9)', sr(:5, i, j), EXP(sr(6, i, j))
      end do
      print *
    end do
    call clear_molecule()
!
  end subroutine test1
!
  subroutine test2()
    integer(IK), parameter   :: n_apm = 8, n_mol = 5, n_target = 50
    real(RK), allocatable    :: X(:, :, :, :)
    integer(IK)              :: n_dim, n_atm
    integer(IK)              :: n_header, n_int, n_float, n_job, n_mem
    integer(IK)              :: edges(2, n_target - 1)
    real(RK)                 :: weights(n_target - 1)
    real(RK), allocatable    :: W(:)
    integer(IK), allocatable :: header(:), int_states(:, :, :)
    real(RK), allocatable    :: float_states(:, :, :)
    integer(IK)              :: i
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
    call n_atoms(n_dim, n_atm)
    call workmemory_lengthes(n_mem, n_job)
    call state_vector_lengthes(n_header, n_int, n_float)
!
    print'(*(I4))', n_dim, n_atm, n_header, n_int, n_float, n_job, n_mem
!
    allocate (X(n_dim, n_apm, n_mol, n_target))
    allocate (w(n_job * n_mem))
    allocate (header(n_header))
    allocate (int_states(n_int, n_target, n_target))
    allocate (float_states(n_float, n_target, n_target))
!
    call RANDOM_NUMBER(X)
!
    call min_span_tree( &
 &    n_target, n_header, n_int, n_float,&
 &    X, W, RHUGE, ZERO, -1, .true., .true., &
 &    edges, weights, header, int_states, float_states)
!
    do i = 1, n_target - 1
      print *, edges(:, i), weights(i)
    end do
    print *
!
    call clear_molecule()
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
