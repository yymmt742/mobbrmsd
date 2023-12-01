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
  public :: mol_block
  public :: block_lower_bound_worksize, block_lower_bound
!
  type mol_block
    sequence
    integer(IK) :: m = 1
    integer(IK) :: n = 1
    integer(IK) :: f = 1
    integer(IK) :: g = 1
  end type mol_block
!
  type mol_block_
    private
    sequence
    integer(IK) :: m
    integer(IK) :: n
    integer(IK) :: f
    integer(IK) :: g
    integer(IK) :: mg
    integer(IK) :: mn
    integer(IK) :: df
    integer(IK) :: dm
    integer(IK) :: ff
    integer(IK) :: gg
    integer(IK) :: dmn
    integer(IK) :: dmg
    integer(IK) :: m_f
    integer(IK) :: ffg
    integer(IK) :: ffgg
    integer(IK) :: dm_df
    integer(IK) :: proc_g
    integer(IK) :: proc_f
    integer(IK) :: ix
    integer(IK) :: iw
    integer(IK) :: ip
    integer(IK) :: iq
    integer(IK) :: rest
    integer(IK) :: xtry
    integer(IK) :: cost
    integer(IK) :: cp
    integer(IK) :: wp
    integer(IK) :: cq
    integer(IK) :: wq
    integer(IK) :: bp
    integer(IK) :: bq
    integer(IK) :: bs
  end type mol_block_
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
  pure elemental function invalid(b) result(res)
    type(mol_block), intent(in) :: b
    logical                     :: res
    res = (b%m < 1) .or. (b%n < 1) .or. (b%m < b%f) .or. (b%n < b%g)
  end function invalid
!
  pure elemental function b2b_(d, x, w, b) result(res)
    integer(IK), intent(in)     :: d
!!! spatial dimension
    integer(IK), intent(in)     :: x, w
!!! pointer to x, y, p, q, w
    type(mol_block), intent(in) :: b
    type(mol_block_)            :: res
    res%m = b%m
    res%n = b%n
    res%g = b%g
    res%f = b%f
    res%mg = res%m * res%g
    res%mn = res%m * res%n
    res%df = d * res%f
    res%dm = d * res%m
    res%ff = res%f * res%f
    res%gg = res%g * res%g
    res%ffg = res%ff * res%g
    res%dmn = d * res%mn
    res%dmg = d * res%mg
    res%m_f = res%m - res%f
    res%ffgg = res%ff * res%gg
    res%dm_df = d * res%m_f
    res%proc_g = Procrustes_worksize(res%g)
    res%proc_f = Procrustes_worksize(res%f)
!!! relative address to X
    res%ix = x
!!! memory block, with absolute index.
!!! covariance matrices and residue matrix, REST(g, g), XTRY(f, f, g, g)
    res%iw   = w
    res%rest = w
    res%xtry = res%rest + res%gg
    res%cost = res%xtry + res%ffgg
!!! rotation matrix P(g,g)
    res%ip = res%cost + 1
!!! rotation matrix Q(f,f,g)
    res%iq = res%ip + res%gg
!!! covariance matrix for P estimation, CP(g, g), w(*)
    res%cp = res%iq + res%ffg
    res%wp = res%cp + res%gg
!!! covariance matrix for Q estimation, CQ(f, f, g), w(*)
    res%cq = res%cp
    res%wq = res%cq + res%ff
    res%bp = res%gg + res%proc_g
    res%bq = res%ff + res%proc_f
    res%bs = res%gg + res%ffgg + 1 + res%gg + res%ffg + MAX(res%bp, res%bq * res%g)
  end function b2b_
!
!| Calculate work array size for d*d matrix.
  pure function block_lower_bound_worksize(d, s, b) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)           :: s
    !! block size
    type(mol_block), intent(in)       :: b(s)
    !! mol_block
    type(mol_block_)                  :: b_(s)
    integer(IK)                       :: dd, dmn, res
!
    if (d < 1 .or. s < 1 .or. ANY(invalid(b))) then
      res = 0
    else
      b_ = b2b_(d, 1, 1, b)
      dd = d * d
      dmn = d * SUM(b_%mn)
      res = SUM(b_%bs) + MAX(dmn, dd + dd + Kabsch_worksize(d)) + dmn + 1
    end if
!
  end function block_lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2, Q is block rotation matrix
  pure subroutine block_lower_bound(d, s, b, X, Y, w, maxiter, threshold, R0)
    integer(IK), intent(in)           :: d
    !! matrix dimension 1.
    integer(IK), intent(in)           :: s
    !! matrix dimension 1.
    type(mol_block), intent(in)       :: b(s)
    !! mol_block
    real(RK), intent(in)              :: X(*)
    !! reference d*sum_i (mi*ni) array
    real(RK), intent(in)              :: Y(*)
    !! target d*sum_i (mi*ni) array
    real(RK), intent(inout)           :: w(*)
    !! work array, must be larger than lower_bound_worksize(d, n, free_indices)
    integer(IK), intent(in), optional :: maxiter
    !! iteration limit, default = 1000
    real(RK), intent(in), optional    :: threshold
    !! iteration limit, default = 1E-12
    real(RK), intent(in), optional    :: R0(*)
    !! initial Rotation matrix
    integer(IK)                       :: maxiter_
    type(mol_block_)                  :: b_(s)
    integer(IK)                       :: dd, mn, mg, dmn
    integer(IK)                       :: cost, prev, thre, cf, cr, rw, r
    integer(IK)                       :: iy, iz
    integer(IK)                       :: i, j
!
    if (d < 1 .or. s < 1 .or. ANY(invalid(b))) then
      w(1) = rmsd(d, SUM(b%m * b%n), X, Y)
      return
    end if
!
    dd = d * d
    mn = SUM(b%m * b%n)
    mg = SUM(b%m * b%g)
    dmn = d * mn
!
    cost = 1
    prev = 2
    thre = 3
    iy = thre + 1                ! copy of Y
    iz = iy + dmn                ! copy of Y (caution, interference with others)
    cf = iy + dmn                ! covariance matrix X^TYPQ, independent of P and Q.
    cr = cf + dd                 ! covariance matrix X^TYPQ
    rw = cr + dd                 ! covariance matrix X^TYPQ
    r  = rw + Kabsch_worksize(d) ! rotation matrix R
!
    maxiter_ = optarg(maxiter,   DEF_maxiter)
    w(thre) = optarg(threshold, DEF_threshold) * dmn
!
    block
      integer(IK) :: iw, iix
      iix = 1
      iw = MAX(r + dd, cf + dmn)
      do i = 1, s
        b_(i) = b2b_(d, iix, iw, b(i))
        iix = iix + b_(i)%dmn
        iw  = iw + b_(i)%bs
      end do
    end block
!
    call get_C_fix(d, s, b_, X, Y, w(cf))
!
    w(cost) = -RHUGE
    call copy(dmn, Y, w(iy))
!
!!! R = R0 or R = I
    if (PRESENT(R0)) then
      call copy(dd, R0, w(r))
    else
      call eye(d, w(r))
    end if
!
    do i = 1, maxiter_
      w(prev) = w(cost)
      do concurrent(j=1:s)
        call update_PQ(d, b_(j), X(b_(j)%ix), Y(b_(j)%ix), w(r), maxiter_, w(thre), &
       &               w(b_(j)%ip), w(b_(j)%iq), W)
      end do
      call update_R(d, dd, s, b_, X, Y, w(cf), w(r), W(cr), w(cost), W(iy), W(rw), W)
      if (ABS(w(cost) - w(prev)) < w(thre)) exit
    end do
!
    call R_rotation(d, mn, w(r), Y, w(iy))
    do concurrent(i=1:s)
      block
        integer(IK) :: iiy, iiz
        iiy = iy + b_(i)%ix - 1
        iiz = iz + b_(i)%ix - 1
        call PQ_rotation(d, b_(i), w(b_(i)%ip), w(b_(i)%iq), w(iiy), w(iiz))
        call copy(b_(i)%dmg, w(iiz), w(iiy))
      end block
    end do
!
    w(1) = rmsd(d, mn, X, w(iy))
!
  contains
!
    pure subroutine calc_swap_indices(n, f, free_indices, res)
      integer(IK), intent(in)    :: n, f, free_indices(f)
      integer(IK), intent(inout) :: res(n)
      integer(IK)                :: i, fi, n_f
      do concurrent(i=1:f)
        res(i) = free_indices(i)
      end do
      n_f = n - f
      fi  = f + 1
      do i = 1, n
        if (ANY(i == res(:f))) cycle
        res(fi) = i
        if (fi == n) exit
        fi = fi + 1
      end do
    end subroutine calc_swap_indices
!
  end subroutine block_lower_bound
!
  pure subroutine get_C_fix(d, s, b, X, Y, CF)
    integer(IK), intent(in)      :: d, s
    type(mol_block_), intent(in) :: b(s)
    real(RK), intent(in)         :: X(*), Y(*)
    real(RK), intent(inout)      :: CF(d, d)
    integer(IK)                  :: i, ix, mn_mg
      call zfill(d * d, CF)
      do i = 1, s
        ix = b(i)%ix + b(i)%dmg
        mn_mg = b(i)%mn - b(i)%mg
        call DGEMM('N', 'T', d, d, mn_mg, ONE, X(ix), d, Y(ix), d, ONE, CF, d)
      end do
  end subroutine get_C_fix
!
  pure subroutine update_PQ(d, b, X, Y, R, maxiter, threshold, P, Q, W)
    integer(IK), intent(in)      :: d, maxiter
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: X(*), Y(*), R(d, d), threshold
    real(RK), intent(inout)      :: P(b%g, b%g), Q(b%f, b%f, b%g)
    real(RK), intent(inout)      :: W(*)
    real(RK)                     :: prev
    integer(IK)                  :: i, j, l
!
    w(b%cost) = -RHUGE
!
    call update_XTRY(d, b, X, Y, R, w(b%xtry), w(b%rest))
!
    call eye(b%g, P)
    do concurrent(i=1:b%g)
      call eye(b%f, Q(1, 1, i))
    enddo
!
    do l = 1, maxiter
!
      prev = w(b%cost)
!
!!!   CP = RESTR
      call copy(b%gg, w(b%rest), w(b%cp))
      do concurrent(i=1:b%g, j=1:b%g)
!!!     CP(i,j) = trace(Qi@Cij) = ddot(Qi^T, Cij)
        block
          integer(IK) :: icp, ixtry
          icp = b%cp + i + (j - 1) * b%g - 1
          ixtry = b%xtry + b%ff * ((i - 1) + b%g * (j - 1))
          w(icp) = w(icp) + ddot(b%ff, Q(1, 1, i), w(ixtry))
        end block
      end do
!!!   Get P = UV^T
      call Procrustes(b%g, w(b%cp), P, w(b%wp))
!!!   tr(P^T, CP) = ddot(P, CP)
      w(b%cost) = ddot(b%gg, P, w(b%cp))
      if (ABS(w(b%cost) - prev) < threshold) exit
!
      do concurrent(j=1:b%g)
        block
          integer(IK) :: icq, iwq, ixtry
          icq = b%cq + b%bq * (j - 1)
          iwq = b%wq + b%bq * (j - 1)
          ixtry = b%xtry + b%ff * (j - 1)
!!        CQ^T = sum_i Pi * Cji
          call zfill(b%ff, w(icq))
          do i = 1, b%g
            call add(b%ff, P(j, i), w(ixtry), w(icq))
            ixtry = ixtry + b%ffg
          end do
!!!       Get Qi^T = UV^T
          call Procrustes(b%f, w(icq), Q(1, 1, j), w(iwq))
        end block
      end do
    end do
!
  end subroutine update_PQ
!
  pure subroutine update_R(d, dd, s, b, X, Y, CF, R, CR, trace, Z, RW, W)
    integer(IK), intent(in)      :: d, dd, s
    type(mol_block_), intent(in) :: b(s)
    real(RK), intent(in)         :: X(*), Y(*), CF(*)
    real(RK), intent(inout)      :: R(*), CR(*), trace, Z(*), RW(*), W(*)
    integer(IK)                  :: i
!!    W = CF + YPQXT
      do concurrent(i = 1:s)
!!      PQ rotation
        call PQ_rotation(d, b(i), w(b(i)%ip), w(b(i)%iq), Y(b(i)%ix), Z(b(i)%ix))
      enddo
!
      call copy(dd, CF, CR)
      do i = 1, s
!!      CR = X(d,mn)@Z^T(mn,d)@P^T@Q^T
        call DGEMM('N', 'T', d, d, b(i)%mg, ONE, X(b(i)%ix), d, Z(b(i)%ix), d, ONE, CR, d)
      enddo
!!    get optimal R
      call Kabsch(d, CR, R, RW)
      trace = ddot(dd, CR, R)
  end subroutine update_R
!
  pure subroutine R_rotation(d, n, R, X, W)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: R(*), X(*)
    real(RK), intent(inout) :: W(*)
    call DGEMM('N', 'N', d, n, d, ONE, R, d, X, d, ZERO, W, d)
  end subroutine R_rotation
!
  pure subroutine PQ_rotation(d, b, P, Q, Y, Z)
    integer(IK), intent(in)      :: d
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: P(*), Q(*), Y(*)
    real(RK), intent(inout)      :: Z(*)
    integer(IK)                  :: i
!
    call DGEMM('N', 'N', b%dm, b%g, b%g, ONE, Y, b%dm, P, b%g, ZERO, Z, b%dm)
!
    do concurrent(i=1:b%g)
      block
        integer(IK) :: iq, ix
        real(RK)    :: Ti(b%df)
        ix = (i - 1) * b%dm + 1
        iq = (i - 1) * b%ff + 1
        call DGEMM('N', 'T', d, b%f, b%f, ONE, Z(ix), d, Q(iq), b%f, ZERO, Ti, d)
        call copy(b%df, Ti, Z(ix))
      end block
    end do
!
  end subroutine PQ_rotation
!
  pure subroutine update_XTRY(d, b, X, Y, R, XTRY, RESTR)
    integer(IK), intent(in)      :: d
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: X(d, b%m, b%n), Y(d, b%m, b%n), R(d, d)
    real(RK), intent(inout)      :: XTRY(b%f, b%f, b%g, b%g), RESTR(b%g, b%g)
    integer(IK)                  :: i, j
!
    if (b%f > 0) then
      do concurrent(i=1:b%g, j=1:b%g)
        block
          real(RK) :: temp(b%df)
          call DGEMM('T', 'N', b%f, d, d, ONE, X(1, 1, i), d, R, d, ZERO, temp, b%f)
          call DGEMM('N', 'N', b%f, b%f, d, ONE, temp, b%f, Y(1, 1, j), d, ZERO, XTRY(1, 1, i, j), b%f)
        end block
      end do
    end if
!
    if (b%m_f < 1) return
!
    do concurrent(i=1:b%g, j=1:b%g)
      block
        real(RK) :: RY(b%dm_df)
        call DGEMM('N', 'N', d, b%m_f, d, ONE, R, d, Y(1, b%f + 1, j), d, ZERO, RY, d)
        RESTR(i, j) = ddot(b%dm_df, X(1, b%f + 1, i), RY)
      end block
    end do
!
  end subroutine update_XTRY
!
  pure subroutine copy(d, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = source(i)
    end do
  end subroutine copy
!
  pure subroutine add(d, c, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: c, source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = dest(i) + c * source(i)
    end do
  end subroutine add
!
  pure function ddot(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    integer(IK)             :: i
    res = ZERO
    do i = 1, d
      res = res + X(i) * Y(i)
    end do
  end function ddot
!
  pure subroutine eye(d, X)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: X(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:d)
      X(i, j) = MERGE(ONE, ZERO, i == j)
    end do
  end subroutine eye
!
  pure subroutine zfill(d, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = ZERO
    end do
  end subroutine zfill
!
end module mod_block_lower_bound
