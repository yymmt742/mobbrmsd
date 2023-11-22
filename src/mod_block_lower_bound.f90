!| Calculate min_{R,P,Q1,...,Q_n} |X-R@Y@f(P)@g(Q1,Q2,...,Qn)|^2 for R\in\R^{d,d}, P\in \R^{n,n} and Qi\in \R^{m,m}.
!  f(P) = I x P ( I is the identity and \in R^{m,m} ).
!  g(Q1,Q2,...,Qn) = Diag(Q1 Q2 ... Qn).
!  Here, R@R^T=I, det(R)=1, P@P^T=I, and Q@Q^T=I are satisfied.
module mod_block_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, FOUR => RFOUR, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_rmsd
  use mod_det
  use mod_cov
  use mod_Kabsch
  use mod_Procrustes
  implicit none
  private
  public :: block_lower_bound_worksize, block_lower_bound
!
  interface
    include 'dgemm.h'
  end interface
!
  integer(IK), parameter :: DEF_maxiter   = 1000_IK
  real(RK), parameter    :: DEF_threshold = 1.0E-12_RK
!
contains
!
!| Calculate work array size for d*d matrix.
  pure function block_lower_bound_worksize(d, m, n, free_m_indices, free_n_indices) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)           :: m
    !! matrix dimension 2.
    integer(IK), intent(in)           :: n
    !! matrix dimension 3.
    integer(IK), intent(in)           :: free_m_indices(:)
    !! rotation index, must be SIZE(free_indices) <= m
    integer(IK), intent(in)           :: free_n_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    integer(IK)                       :: res
    integer(IK)                       :: f, g
!
    if (d < 1 .or. m < 1 .or. n < 1) then
      res = 0
    else
      f = SIZE(free_m_indices)
      g = SIZE(free_n_indices)
      res = 2 + 3 * d * m * n + MAX(d * d * 2 + Kabsch_worksize(d), (m * n)**2 + g * (2 * f * f + procrustes_worksize(f)))
    end if
!
  end function block_lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2, Q is block rotation matrix
  subroutine block_lower_bound(d, m, n, free_m_indices, free_n_indices, X, Y, w, maxiter, threshold)
    integer(IK), intent(in)           :: d
    !! matrix dimension 1.
    integer(IK), intent(in)           :: m
    !! matrix dimension 2.
    integer(IK), intent(in)           :: n
    !! matrix dimension 3.
    integer(IK), intent(in)           :: free_m_indices(:)
    !! rotation index, must be SIZE(free_indices) <= m
    integer(IK), intent(in)           :: free_n_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    real(RK), intent(in)              :: X(*)
    !! reference d*m*n array
    real(RK), intent(in)              :: Y(*)
    !! target d*m*n array
    real(RK), intent(inout)           :: w(*)
    !! work array, must be larger than lower_bound_worksize(d, n, free_indices)
    integer(IK), intent(in), optional :: maxiter
    !! iteration limit, default = 1000
    real(RK), intent(in), optional    :: threshold
    !! iteration limit, default = 1E-12
    integer(IK)                       :: maxiter_
    real(RK)                          :: threshold_
    integer(IK)                       :: n_indices(n), m_indices(m)
    integer(IK)                       :: f, g, mg, mn, dd, dm, df, ff, dfg, dmg, dmn, mmn, mmgg
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, rcov, rrot, cov, mrot, nrot
    integer(IK)                       :: i, j, iw1, iw2
!
    real(RK)                          :: R(d, d), P(SIZE(free_n_indices), SIZE(free_n_indices)), LAM(4)
    real(RK)                          :: Q(SIZE(free_m_indices), SIZE(free_m_indices), SIZE(free_n_indices))
    real(RK)                          :: GR(d, d), GP(SIZE(free_n_indices), SIZE(free_n_indices)), GLAM(4)
    real(RK)                          :: GQ(SIZE(free_m_indices), SIZE(free_m_indices), SIZE(free_n_indices))
    integer(IK)                       :: k
!
    f = SIZE(free_m_indices)
    g = SIZE(free_n_indices)
print*,d, m, n, f, g
!
    if (f < 1 .or. g < 1) then
      w(1) = rmsd(d, m * n, X, Y)
      return
    end if
!
    mg = m * g
    mn = m * n
    dm = d * m
    dd = d * d
    df = d * f
    ff = f * f
    dmn = d * m * n
    mmn = m * m * n
    mmgg = mg * mg
    dfg = d * f * g
    dmg = d * m * g
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1      ! copy of X
    iy   = ix   + dmn    ! copy of Y
    iz   = iy   + dmn    ! Z = P@Y
    rcov = iz   + dmn    ! covariance matrix X^TYPQ
    rrot = rcov + dd     ! rotation matrix R
    iw1  = rrot + dd
    cov  = iz   + dmn    ! covariance matrix XRY^T
    mrot = cov  + mmgg   ! rotation matrix P
    nrot = mrot + ff * g ! rotation matrix Q
    iw2 = nrot + g * g
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) * dmn
!
    call calc_swap_indices(m, f, free_m_indices, m_indices)
    call calc_swap_indices(n, g, free_n_indices, n_indices)
!
    call pack_matrix(d, m, n, m_indices, n_indices, X, w(ix))
    call pack_matrix(d, m, n, m_indices, n_indices, Y, w(iy))
!
    do concurrent(j=0:dmn-1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = -RHUGE
!
do concurrent(i=1:g,j=1:g)
  P(i, j) = MERGE(ONE, ZERO, i == j)
enddo
do concurrent(i=1:f,j=1:f,k=1:g)
  Q(i, j, k) = MERGE(ONE, ZERO, i == j)
enddo
do concurrent(i=1:d,j=1:d)
  R(i, j) = MERGE(ONE, ZERO, i == j)
enddo
LAM = ONE
!
print*, objective(d, m, n, f, g, P, Q, R, LAM, X, Y)
!
call random_number(P)! ; P = ZERO
call random_number(Q)! ; Q = ZERO
call random_number(R)
call random_number(LAM)
print'(15f6.1)',PxQ(m, n, f, g, P, Q)
call gradient(d, m, n, f, g, P, Q, R, LAM, X, Y, GP, GQ, GR, GLAM)
print'(2f9.3)', GP
print*
print'(3f9.3)', GQ
print*
print'(3f9.3)', GR
print*
print'(4f9.3)', GLAM
print*
call numgrad(d, m, n, f, g, P, Q, R, LAM, X, Y, GP, GQ, GR, GLAM)
print'(2f9.3)', GP
print*
print'(3f9.3)', GQ
print*
print'(3f9.3)', GR
print*
print'(4f9.3)', GLAM
print*
return
!
do i = 1, 10
  call gradient(d, m, n, f, g, P, Q, R, LAM, X, Y, GP, GQ, GR, GLAM)
  R = R - GR * 0.001
  P = P - GP * 0.001
  Q = Q - GQ * 0.001
  LAM = LAM - GLAM * 0.001
  print *, objective(d, m, n, f, g, P, Q, R, LAM, X, Y)
end do
print'(2f9.3)', P
print*
print'(3f9.3)', Q
print*
print'(3f9.3)', R
print*
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
!!    M = X^T@Z
      call DGEMM('T', 'N', mg, mg, d, ONE, w(ix), d, w(iz), d, ZERO, w(cov), mg)
      call block_rotation(f, g, m, mg, w(cov), w(mrot), w(nrot), maxiter_, threshold_)
!
!!    R rotation
!!    M = X(d,mn)@Z^T(mn,d)@P^T@Q^T
      w(iz:iz + dmg - 1) = w(iy:iy + dmg - 1)
      call P_rotation(d, m, n, f, g, w(mrot), w(iz))
      call Q_rotation(dm, n, g, w(nrot), w(iz))
      call DGEMM('N', 'T', d, d, mn, ONE, w(ix), d, w(iz), d, ZERO, w(rcov), d)
!!    get optimal R
      call Kabsch(d, w(rcov), w(rrot), w(iw1))
!
!!    trace(MR) = trace(MR^T)
!
      w(conv) = ZERO
      do j = 0, dd - 1
        w(conv) = w(conv) + w(rcov + j) * w(rrot + j)
      end do
!
      if (ABS(w(prev) - w(conv)) < threshold_) exit
!
!!    Z = R@Y
      call DGEMM('N', 'N', d, mn, d, ONE, w(rrot), d, w(iy), d, ZERO, w(iz), d)
!
    enddo
!
    w(iz:iz + dmg - 1) = w(iy:iy + dmg - 1)
    call P_rotation(d, m, n, f, g, w(mrot), w(iz))
    call Q_rotation(dm, n, g, w(nrot), w(iz))
    w(conv) = rmsd(d, mn, w(ix), w(iz))
!
  contains
!
    pure subroutine pack_matrix(d, m, n, m_indices, n_indices, source, dest)
      integer(IK), intent(in) :: d, m, n, m_indices(m), n_indices(n)
      real(RK), intent(in)    :: source(d, m, n)
      real(RK), intent(inout) :: dest(d, m, n)
      integer(IK)             :: i, j, k
!
      do concurrent(k=1:n, j=1:m, i=1:d)
        dest(i, j, k) = source(i, m_indices(j), n_indices(k))
      end do
!
    end subroutine pack_matrix
!
  end subroutine block_lower_bound
!
  pure subroutine block_rotation(f, g, m, mg, cov, mrot, nrot, maxiter, threshold)
    integer(IK), intent(in) :: f, g, m, mg, maxiter
    real(RK), intent(in)    :: cov(mg, mg)
    real(RK), intent(in)    :: threshold
    real(RK), intent(inout) :: mrot(f, f, g), nrot(g, g)
    real(RK)                :: mcov(f, f), ncov(g, g), ncov0(g, g)
    real(RK)                :: w(procrustes_worksize(f))
    real(RK)                :: prev, conv
    integer(IK)             :: i, j, k
!
    do concurrent(j=1:g, k=1:g)
      block
        integer(IK) :: jm, km
        jm = (j - 1) * m
        km = (k - 1) * m
        ncov0(j, k) = SUM(cov(jm + f + 1:jm + m, km + 1:km + f))&
       &            + SUM(cov(jm + 1:jm + m, km + f + 1:km + m))
      end block
    end do
!
    prev = -RHUGE
    conv = -RHUGE
!
    do concurrent(i=1:g, j=1:g)
      nrot(i, j) = MERGE(1, 0, i == j)
    end do
!
    do concurrent(i=1:f, j=1:f, k=1:g)
      mrot(i, j, k) = MERGE(1, 0, i == j)
    end do
!
    do i = 1, maxiter
!
      if (prev < conv) prev = conv
!
!!    P rotation
!
      do concurrent(j=1:g)
        mcov = ZERO
        do k = 1, g
          mcov(:, :) = mcov(:, :) &
         &           + nrot(k, j) * cov((k - 1) * m + 1:(k - 1) * m + f, (k - 1) * m + 1:(k - 1) * m + f)
        end do
        call Procrustes(f, mcov, mrot(:, :, j), w)
      end do
!
!!    Q rotation
      do concurrent(j=1:g, k=1:g)
        block
          integer(IK) :: jm, km
          jm = (j - 1) * m
          km = (k - 1) * m
          ncov(j, k) = SUM(MATMUL(cov(jm + 1:jm + f, km + 1:km + f), TRANSPOSE(mrot(:, :, k))))
        end block
      end do
!
      call Procrustes(g, ncov, nrot, w)
!
      conv = SUM(ncov * nrot)
      if (ABS(prev - conv) < threshold) exit
!
    end do
!
  end subroutine block_rotation
!
  pure subroutine P_rotation(d, m, n, f, g, r, z)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: r(f, f, g)
    real(RK), intent(inout) :: z(d, m, n)
    real(RK)                :: t(d, f)
    integer(IK)             :: i
!
    do i = 1, g
      t = MATMUL(z(:, :f, i), TRANSPOSE(r(:, :, i)))
      z(:, :f, i) = t
      !z(:, :f, i) = MATMUL(z(:, :f, i), TRANSPOSE(r(:, :, i)))
    end do
!
  end subroutine P_rotation
!
  pure subroutine Q_rotation(dm, n, g, r, z)
    integer(IK), intent(in) :: dm, n, g
    real(RK), intent(in)    :: r(g, g)
    real(RK), intent(inout) :: z(dm, n)
!
    z(:, :g) = MATMUL(z(:, :g), TRANSPOSE(r))
!
  end subroutine Q_rotation
!
  pure subroutine calc_swap_indices(n, f, free_indices, res)
    integer(IK), intent(in)    :: n, f, free_indices(f)
    integer(IK), intent(inout) :: res(n)
    integer(IK)                :: i, j, g
!
    do concurrent(i=1:f)
      res(i) = free_indices(i)
    end do
!
    g = n - f
    j = f + 1
!
    do i = 1, n
      if (ANY(i == res(:f))) cycle
      res(j) = i
      if (j == n) exit
      j = j + 1
    end do
!
  end subroutine calc_swap_indices
!
   function objective(d, m, n, f, g, P, Q, R, L, X, Y) result(res)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: P(g, g), Q(f, f, g), R(d, d), L(4)
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n)
    real(RK)                :: Cji(d, d)
    real(RK)                :: res
    integer(IK)             :: i, j
!
    Cji = ZERO
!
    do j = 1, g
      do i = 1, g
        Cji = Cji + P(i, j) * MATMUL(MATMUL(X(:, :f, j), Q(:, :, i)), TRANSPOSE(Y(:, :f, j)))
      end do
      Cji = Cji + SUM(P(:, j)) * MATMUL(X(:, f + 1:, j), TRANSPOSE(Y(:, f + 1:, j)))
    end do
    Cji = Cji + MATMUL(RESHAPE(X(:, :, g + 1:), [d, m * (n - g)]), TRANSPOSE(RESHAPE(Y(:, :, g + 1:), [d, m * (n - g)])))
!
    res = SUM(R * Cji)
!
    return
Cji = MATMUL(MATMUL(RESHAPE(X, [d, m * n]), PxQ(m, n, f, g, P, Q)), TRANSPOSE(RESHAPE(Y, [d, m * n])))
res = SUM(R * Cji)
!   res = res - L(1) * constraint_1(g, P)
!   do i = 1, g
!     res = res - L(2) * constraint_1(f, Q(:, :, i))
!   end do
!   res = res - L(3) * constraint_1(d, R)
!   res = res - L(4) * constraint_2(d, R)
!
  end function objective
!
   subroutine gradient(d, m, n, f, g, P, Q, R, L, X, Y, GP, GQ, GR, GL)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: P(g, g), Q(f, f, g), R(d, d), L(4)
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n)
    real(RK), intent(inout) :: GP(g, g), GQ(f, f, g), GR(d, d), GL(4)
    real(RK)                :: Cij(f, f), TR(d, d)
    integer(IK)             :: i, j, k
!
!!! GP = - lambda_P * d(g(P))/dP
    call constraint_1_grad(g, P, GP)
    do concurrent(i=1:g, j=1:g)
      GP(i, j) = - L(1) * GP(i, j)
    end do
!
!!! GQ_I = - lambda_Q * d(g(Q_I))/dQ_I
    do concurrent(k=1:g)
      call constraint_1_grad(f, Q(:, :, k), GQ(:, :, k))
      do concurrent(i=1:f, j=1:f)
        GQ(i, j, k) = - L(2) * GQ(i, j, k)
      end do
    enddo
!
!!! GR = - lambda_R * d(g(R))/dR - lambda_|R| * d(h(R))/dR
    call constraint_1_grad(d, R, GR)
    call constraint_2_grad(d, R, TR)
    do concurrent(i=1:d, j=1:d)
      GR(i, j) = - L(3) * GR(i, j) - L(4) * TR(i, j)
    enddo
!
GP = ZERO
GQ = ZERO
GR = ZERO
GL = ZERO
!
    do j = 1, g
      do i = 1, g
        Cij = MATMUL(MATMUL(TRANSPOSE(X(:, :f, i)), R), Y(:, :f, j))
        GP(j, i) = GP(j, i) + SUM(Cij * Q(:, :, j))
        GQ(:, :, j) = GQ(:, :, j) + P(j, i) * Cij
        TR = MATMUL(Y(:, :f, j), MATMUL(Q(:, :, i), TRANSPOSE(X(:, :f, j)))) &
           + MATMUL(Y(:, f + 1:, j), TRANSPOSE(X(:, f + 1:, j)))
        GR(:, :) = GR(:, :) + P(j, i) * TR(:, :)
      end do
    end do
!
    do i = g + 1, n
      GR(:, :) = GR(:, :) + MATMUL(Y(:, :, i), TRANSPOSE(X(:, :, i)))
    end do
!
GR = MATMUL(MATMUL(RESHAPE(X, [d, m * n]), PxQ(m, n, f, g, P, Q)), TRANSPOSE(RESHAPE(Y, [d, m * n])))
return
!

    GL(1) = constraint_1(g, P)
    GL(2) = ZERO
    do i = 1, g
      GL(2) = GL(2) + constraint_1(f, Q(:, :, i))
    end do
    GL(3) = constraint_1(d, R)
    GL(4) = constraint_2(d, R)
!
  end subroutine gradient
!
  pure function PxQ(m, n, f, g, P, Q) result(res)
    integer(IK), intent(in) :: m, n, f, g
    real(RK), intent(in)    :: P(g, g), Q(f, f, g)
    real(RK)                :: res(m * n, m * n), eye(m - f, m - f)
    integer(IK)             :: i, j
    do concurrent(i=1:m - f, j=1:m - f)
      eye(i, j) = MERGE(ONE, ZERO, i==j)
    end do
    do concurrent(i=1:m * n, j=1:m * n)
      res(i, j) = MERGE(ONE, ZERO, i==j)
    end do
    do concurrent(i=1:g, j=1:g)
      res(1 + m * (i - 1):f + m * (i - 1), 1 + m * (j - 1):f + m * (j - 1)) = P(i, j) * Q(:, :, j)
      res(f + 1 + m * (i - 1):m * i, f + 1 + m * (j - 1):m * j) = P(i, j) * eye
    end do
  end function PxQ
!
  pure function constraint_1(n, A) result(res)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: A(n, n)
    real(RK)                :: AA(n, n)
    real(RK)                :: res
    AA = MATMUL(A, TRANSPOSE(A))
    res = SUM(AA * AA) - 2 * SUM(A * A) + n
  end function constraint_1
!
  pure subroutine constraint_1_grad(n, A, G)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: A(n, n)
    real(RK), intent(inout) :: G(n, n)
    real(RK)                :: B(n, n)
    integer(IK)             :: i
!
!!! B = A@A^T - I
    call DGEMM('N', 'T', n, n, n, ONE, A, n, A, n, ZERO, B, n)
    do concurrent(i=1:n)
      B(i, i) = B(i, i) - ONE
    enddo
!!! G = B@A = 4 *(A@A^T - I)@A
    call DGEMM('N', 'N', n, n, n, FOUR, B, n, A, n, ZERO, G, n)
!
  end subroutine constraint_1_grad
!
  pure function constraint_2(n, A) result(res)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: A(n, n)
    real(RK)                :: res
    res = (det_(n, A) - ONE)**2
  end function constraint_2
!
  pure subroutine constraint_2_grad(n, A, G)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: A(n, n)
    real(RK), intent(inout) :: G(n, n)
    real(RK)                :: tmp
    integer(IK)             :: i, j
!
    tmp = det_(n, A) - ONE
!
    do concurrent(i=1:n, j=1:n)
      block
        real(RK)    :: P(n - 1, n - 1)
        integer(IK) :: k, l
        do concurrent(k=1:n - 1, l=1:n - 1)
          block
            integer(IK) :: ik, jl
            ik = MERGE(k, k + 1, k < i)
            jl = MERGE(l, l + 1, l < j)
            P(k, l) = A(ik, jl)
          end block
        end do
        G(i, j) = tmp * det_(n - 1, P)
        G(i, j) = G(i, j) + G(i, j)
      end block
    end do
!
  end subroutine constraint_2_grad
!
   subroutine numgrad(d, m, n, f, g, P, Q, R, L, X, Y, GP, GQ, GR, GL)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: P(g, g), Q(f, f, g), R(d, d), L(4)
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n)
    real(RK), intent(inout) :: GP(g, g), GQ(f, f, g), GR(d,d), GL(4)
    real(RK), parameter     :: dif = 1.0E-8_RK
    integer(IK)             :: i, j, k
    do concurrent(i=1:g, j=1:g)
      block
        real(RK) :: t(g, g), pz, mz
        t = P
        t(i, j) = t(i, j) + dif
        pz = objective(d, m, n, f, g, t, Q, R, L, X, Y)
        t(i, j) = t(i, j) - dif - dif
        mz = objective(d, m, n, f, g, t, Q, R, L, X, Y)
        GP(i, j) = (pz - mz) / (dif + dif)
      end block
    end do
    do concurrent(i=1:f, j=1:f, k=1:g)
      block
        real(RK) :: t(f, f, g), pz, mz
        t = Q
        t(i, j, k) = t(i, j, k) + dif
        pz = objective(d, m, n, f, g, P, t, R, L, X, Y)
        t(i, j, k) = t(i, j, k) - dif - dif
        mz = objective(d, m, n, f, g, P, t, R, L, X, Y)
        GQ(i, j, k) = (pz - mz) / (dif + dif)
      end block
    end do
    do concurrent(i=1:d, j=1:d)
      block
        real(RK) :: t(d, d), pz, mz
        t = R
        t(i, j) = t(i, j) + dif
        pz = objective(d, m, n, f, g, P, Q, t, L, X, Y)
        t(i, j) = t(i, j) - dif - dif
        mz = objective(d, m, n, f, g, P, Q, t, L, X, Y)
        GR(i, j) = (pz - mz) / (dif + dif)
      end block
    end do
    do concurrent(i=1:4)
      block
        real(RK) :: t(4), pz, mz
        t = L
        t(i) = t(i) + dif
        pz = objective(d, m, n, f, g, P, Q, R, t, X, Y)
        t(i) = t(i) - dif - dif
        mz = objective(d, m, n, f, g, P, Q, R, t, X, Y)
        GL(i) = (pz - mz) / (dif + dif)
      end block
    end do
  end subroutine numgrad
end module mod_block_lower_bound
