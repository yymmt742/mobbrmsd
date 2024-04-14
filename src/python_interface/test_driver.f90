program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO, RHUGE
  use driver
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer        :: n_dim, n_atom
!
  call u%init('test python_driver')
!
  call n_atoms(n_dim, n_atom)
  print*,n_dim
!
  call test1(n_dim)
!
  call u%init('test python_driver, nearest_neighbor')
  call test2(n_dim)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(n_dim)
    integer(IK), intent(in)  :: n_dim
    integer(IK), parameter   :: s = 2
    integer(IK), parameter   :: n_target = 1000
    integer(IK), parameter   :: n_apm = 8, n_mol = 8
    integer(IK), parameter   :: n_atm = n_apm * n_mol
    real(RK)                 :: X(n_dim, n_apm, n_mol), Y(n_dim, n_apm, n_mol, n_target)
    integer(IK)              :: n_header, n_int, n_float
    integer(IK)              :: n_job, n_mem
    integer(IK), allocatable :: h(:), si(:, :)
    real(RK), allocatable    :: sr(:, :), W(:, :)
    integer(IK)              :: i, j, k, l
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
!
    do k = -1, 1, 2
      do j = -1, 1, 2
        do i = -1, 1, 2
          X(:, :, (i + 1) * 2 + (j + 1) + (k + 1) / 2 + 1) = sample(n_dim, n_apm, [i, j, k])
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
      Y(:, :, :, l) = 0.2 * X + 0.8 * Y(:, :, :, l)
    end do
!
    call workmemory_lengthes(n_mem, n_job)
    call state_vector_lengthes(n_header, n_int, n_float)
    print*, n_header, n_int, n_float
!
    allocate(h(n_header))
    allocate(si(n_int, n_target))
    allocate(sr(n_float, n_target))
    allocate(W(n_mem, n_job))
!
    call batch_run(n_dim, n_atm, n_target, n_header, n_int, n_float, n_mem, n_job, &
   &               X, Y, W, &
   &               999.0_RK, 0.0_RK, -1, .FALSE., h, si, sr)
    print'(10(I6))',h
   !do i = 1, n_target
   !  print'(*(I4))', si(:, i)
   !end do
    do i = 1, n_target
      print'(5(F12.3))', sr(:4, i), EXP(sr(5, i))
    end do
    call clear_molecule()
!
  end subroutine test1
!
  subroutine test2(n_dim)
    integer(IK), intent(in)  :: n_dim
    integer(IK), parameter   :: n_apm = 8, n_mol = 8, n_target = 1000
    integer(IK), parameter   :: n_atom = n_apm * n_mol
    real(RK)                 :: X(n_dim, n_apm, n_mol), Y(n_dim, n_apm, n_mol, n_target)
    integer(IK)              :: n_job, n_mem
    integer(IK)              :: nn_index
    real(RK)                 :: nn_value
    real(RK), allocatable    :: W(:)
!
    call RANDOM_NUMBER(X)
    call RANDOM_NUMBER(Y)
    Y(:, :, :, 83) = 0.3 * Y(:, :, :, 83) + 0.7 * X
!
    call add_molecule(n_apm, n_mol, 2, [5, 6, 7, 8, 1, 2, 3, 4])
!
    call workmemory_lengthes(n_mem, n_job)
    print*, n_job, n_mem
    allocate (w(n_job * n_mem))
!
    call nearest_neighbor(n_dim, n_atom, n_target, n_mem, n_job,&
 &                        X, Y, W, RHUGE, ZERO, -1, &
 &                        nn_index, nn_value)
!
    print*,nn_index, nn_value
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
    enddo
  end function sample
!
end program main
