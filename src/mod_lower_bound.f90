!| Calculate min_{R,P} |X-RYP|^2
!  Here, RR^T=PP^T=I and det(R)=1 are satisfied.
module mod_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_rmsd
  use mod_cov
  use mod_Kabsch
  use mod_Procrustes
  implicit none
  private
  public :: lower_bound_worksize, lower_bound
!         & block_lower_bound_worksize, block_lower_bound
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
  pure function lower_bound_worksize(d, n, free_indices) result(res)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    integer(IK), intent(in)       :: free_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    integer(IK)                   :: res
    integer(IK)                   :: f
!
    if (d < 1 .or. n < 1) then
      res = 0
    else
      f = SIZE(free_indices)
      res = 2 + 3 * d * n + MAX(d * d * 2 + Kabsch_worksize(d), f * f * 2 + procrustes_worksize(f))
    end if
!
  end function lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2
  pure subroutine lower_bound(d, n, free_indices, X, Y, w, maxiter, threshold)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)           :: n
    !! matrix row dimension.
    integer(IK), intent(in)           :: free_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    real(RK), intent(in)              :: X(*)
    !! reference d*n array
    real(RK), intent(in)              :: Y(*)
    !! target d*n array
    real(RK), intent(inout)           :: w(*)
    !! work array, must be larger than lower_bound_worksize(d, n, free_indices)
    integer(IK), intent(in), optional :: maxiter
    !! iteration limit, default = 1000
    real(RK), intent(in), optional    :: threshold
    !! iteration limit, default = 1E-12
    integer(IK)                       :: i, j, maxiter_
    real(RK)                          :: threshold_
    integer(IK)                       :: f, dn, dd, df, ff
    integer(IK)                       :: n_indices(n)
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, ncov, nrot
    integer(IK)                       :: iw1, iw2
!
    f = SIZE(free_indices)
    if (f < 1) then
      w(1) = rmsd(d, n, X, Y)
      return
    end if
!
    dn = d * n
    dd = d * d
    df = d * f
    ff = f * f
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1  ! copy of X
    iy   = ix   + dn ! copy of Y
    iz   = iy   + dn ! Z = P@Y
    iw1  = iz   + dn
    nrot = iz   + dn ! rotation matrix Q
    ncov = nrot + ff ! covariance matrix, X^TY
    iw2  = ncov + ff
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) ** 2 * n
!
    call calc_swap_indices(n, f, free_indices, n_indices)
!
    call pack_matrix(d, n, n_indices, X, w(ix))
    call pack_matrix(d, n, n_indices, Y, w(iy))
!
    do concurrent(j=0:dn - 1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = -RHUGE
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
      call R_rotation(d, n, w(ix), w(iy), w(iz), w(iw1))
!
!!    M = X^T@Z
      call DGEMM('T', 'N', f, f, d, ONE, w(ix), d, w(iz), d, ZERO, w(ncov), f)
!!    calc P from M
      call Procrustes(f, w(ncov), w(nrot), w(iw2))
!!    trace(XRYP)
      w(conv) = SUM(w(ncov:ncov + ff - 1) * w(nrot:nrot + ff - 1)) &
     &        + SUM(w(ix + df:ix + dn - 1) * w(iz + df:iz + dn - 1))
!
      if (ABS(w(prev) - w(conv)) < threshold_)then
        call DGEMM('N', 'T', d, f, f, ONE, w(iz), d, w(nrot), f, ZERO, w(iy), d)
        do concurrent(j=df:dn - 1)
          w(iy + j) = w(iz + j)
        end do
        exit
      endif
!
!!    Z = Y@P
      call DGEMM('N', 'T', d, f, f, ONE, w(iy), d, w(nrot), f, ZERO, w(iz), d)
      do concurrent(j=df:dn - 1)
        w(iz + j) = w(iy + j)
      end do
!
    enddo
!
    w(conv) = rmsd(d, n, w(ix), w(iy))
!
  contains
!
    pure subroutine pack_matrix(d, n, n_indices, source, dest)
      integer(IK), intent(in) :: d, n, n_indices(n)
      real(RK), intent(in)    :: source(d, n)
      real(RK), intent(inout) :: dest(d, n)
      integer(IK)             :: i, j
!
      do concurrent(j=1:n, i=1:d)
        dest(i, j) = source(i, n_indices(j))
      end do
!
    end subroutine pack_matrix
!
  end subroutine lower_bound
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
    integer(IK)                       :: f, g, mg, mn, dd, df, ff, dfg, dmn, mmn, mmgg
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, ncov, nrot, bs
    integer(IK)                       :: i, j, k, iw1, iw2
!
    f = SIZE(free_m_indices)
    g = SIZE(free_n_indices)
!
    if (f < 1 .or. g < 1) then
      w(conv) = rmsd(d, m * n, X, Y)
      return
    end if
!
    mg = m * g
    mn = m * n
    dd = d * d
    df = d * f
    ff = f * f
    dmn = d * m * n
    mmn = m * m * n
    mmgg = mg * mg
    dfg = d * f * g
    bs = ff + Procrustes_worksize(f)
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1   ! copy of X
    iy   = ix   + dmn ! copy of Y
    iz   = iy   + dmn ! Z = P@Y
    iw1  = iz   + dmn
    ncov = iz   + dmn ! covariance matrix XRY^T
    nrot = ncov + mmgg ! rotation matrix Q
    iw2  = nrot + ff * g
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) ** 2 * n
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
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
      call R_rotation(d, mn, w(ix), w(iy), w(iz), w(iw1))
!
      !!  M = X^T@Z
      call DGEMM('T', 'N', mg, mg, d, ONE, w(ix), d, w(iz), d, ZERO, w(ncov), mg)
!
!     do concurrent(j=0:g-1)
!       block
!         integer(IK) :: ncov_j, nrot_j, iw_j
!         ncov_j = ncov + j * (mmn + m)
!         nrot_j = nrot + j * ff
!         iw_j = iw2 + j * bs
!         call Procrustes(f, w(ncov_j), w(nrot_j), w(iw_j), ldcov=mn)
!         do concurrent(k=0:n - 1)
!           block
!             integer(IK) :: ncov_jk, iw_jk
!             ncov_jk = ncov + k * mmn + j * m
!             iw_jk = iw_j + k * ff
!             call DGEMM('N', 'T', f, f, f, ONE, w(ncov_jk), mn, w(nrot_j), f, ZERO, w(iw_jk), d)
!           end block
!         end do
!       end block
!     end do
!
!     w(conv) = sd(d, mn, w(ix), w(iz))
!     if (ABS(w(prev) - w(conv)) < threshold_) exit
!
!     do concurrent(j=0:g-1)
!       call DGEMM('N', 'T', d, f, f, ONE, w(iy + j * df), d, w(nrot + j * ff), f, ZERO, w(iz + j * df), d)
!     enddo

!     do concurrent(j=dfg:dmn-1)
!       w(iz + j) = w(iy + j)
!     enddo
!
    enddo
!
    w(conv) = rmsd(d, mn, w(ix), w(iy))
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
!   pure subroutine block_rotation(d, f, mn, M, P, w)
!     integer(IK), intent(in)           :: d, f, mn
!     real(RK), intent(in)              :: M(*)
!     real(RK), intent(inout)           :: P(*), w(*)
!     integer(IK)                       :: i, iw, df
!
!!  Z' = Z@P
!     call DGEMM('N', 'T', d, f, f, ONE, z, d, P, f, ZERO, w(iw), d)
!!  Z = Z'
!     iw = iw - 1
!     do concurrent(i=1:df)
!       z(i) = w(iw + i)
!     end do
!
!   end subroutine block_rotation
!
  end subroutine block_lower_bound
!
  pure subroutine R_rotation(d, n, x, y, z, w)
    integer(IK), intent(in)           :: d, n
    real(RK), intent(in)              :: x(*), y(*)
    real(RK), intent(inout)           :: z(*), w(*)
    integer(IK)                       :: ir, iw
!
      ir = d * d + 1
      iw = d * d + ir
!
!!    M = X@Z^T
      call DGEMM('N', 'T', d, d, n, ONE, x, d, z, d, ZERO, w, d)
!!    calc R from M
      call Kabsch(d, w, w(ir), w(iw))
!!    Z = R@Y
      call DGEMM('N', 'N', d, n, d, ONE, w(ir), d, y, d, ZERO, z, d)
!
  end subroutine R_rotation
!
! subroutine P_rotation(d, m, n, f, g, mg, cov, w, maxiter, threshold)
!   integer(IK), intent(in) :: d, m, n, f, g, mg, maxiter
!   real(RK), intent(in)    :: cov(mg, mg), threshold
!   real(RK), intent(inout) :: w(*)
!   integer(IK)             :: ff, mm, bs, iw1, iw2, covn
!   integer(IK)             :: i, j, k
!
!!
!!        | A_11 A_12 A_13 ... A_1g |          | b_11 b_12 b_13 ... b_1m |
!!        | A_21 A_12 A_13 ... A_2g |          | b_21 b_22 b_23 ... b_2m |
!!  cov = | A_31 A_32 A_33 ... A_3g |   A_IJ = | b_31 b_32 b_33 ... b_3m |
!!        |  :    :    :        :   |          |  :    :    :        :   |
!!        | A_g1 A_g2 A_g3 ... A_gg |          | b_m1 b_m2 b_m3 ... b_mm |
!!
!!               | b_11 ... b_1f |
!!   B_IJ(f,f) = |  :        :   |   i, j = 1,2,...,g
!!               | b_f1 ... b_ff |,
!!
!!   Pj minimizes trace[ B_jj@P_j ] and satisfies Pj@Pj^T=I
!!
!!   N      = { M1@P1 + M2@P_2 + ... + Mg@P_g,               if k, l \in 1,2,...,f
!!            { A_11 + ... + A_g1 + ... + A_1g + ... + A_gg, otherwise
!!
!!   where, M_J(f,f)  = B_1J + B_2J + ... + B_nJ
!!
!!   Q minimizes trace[ N@Q ] and satisfies Q@Q^T=I
!!
!
!     ff = f * f
!     mm = m * m
!     bs = 3 * ff + Procrustes_worksize(f)
!     iw1 = 1
!     covn = iw1 + bs * g
!     iw2 = covn + ff
!
!     do concurrent(i = 0:g-1)
!       block
!         integer(IK) :: ib, im, ic
!         im = iw1 + i * bs
!         ib = im + ff
!         ic = m * i + 1
!         call pack_fcov(f, m, g, cov( 1, ic), w(im))
!         call pack_fcov(f, m, 1, cov(ic, ic), w(ib))
!       end block
!     end do
!
!     do concurrent(i=covn:covn + ff - 1)
!       w(i) = ZERO
!     enddo
!     do j = 0, g-1
!       do i = 0, g-1
!         do concurrent(k = 0:mm-1)
!           w(covn + k) = cov(i*m,j)
!         enddo
!       enddo
!     enddo
!
!     do i = 1, maxiter
!
!       do concurrent(j=0:g - 1)
!         block
!           integer(IK) :: jcov, mcov, rot, jw
!           jcov = iw1 + j * bs
!           mcov = jcov + ff
!           jrot = jcov + ff
!           jw   = jrot + ff
!           call Procrustes(f, w(jcov), p(jrot), w(jw))
!         end block
!       end do
!

!       do j=1,g
!         block
!           integer(IK) :: jcov, mcov, rot, jw
!         end block
!       enddo
!
!     end do
!
! contains
!
!   pure subroutine pack_fcov(f, m, s, cov, fcov)
!     integer(IK), intent(in) :: f, m, s
!     real(RK), intent(in)    :: cov(m * g, m)
!     real(RK), intent(inout) :: fcov(f, f)
!     integer(IK)             :: i, j, k
!     do concurrent(i=1:f, j=1:f)
!       fcov(i, j) = ZERO
!     end do
!     do k = 0, s - 1
!       do concurrent(i=1:f, j=1:f)
!         fcov(i, j) = fcov(i, j) + cov(i + k * m, j)
!       end do
!     end do
!   end subroutine pack_fcov
!
!   pure subroutine pack_ncov(m, g, cov, ncov)
!     integer(IK), intent(in) :: f, m
!     real(RK), intent(in)    :: cov(m * g, m * g)
!     real(RK), intent(inout) :: ncov(m, m)
!     integer(IK)             :: i, j, k, l
!     do l = 0, g - 1
!       do k = 0, g - 1
!         do concurrent(i=1:m, j=1:m)
!           ncov(i, j) = jcov(i, j) + cov(i + k * m, j + k * m)
!         end do
!       end do
!     end do
!   end subroutine pack_ncov
!
! end subroutine P_rotation
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
end module mod_lower_bound
