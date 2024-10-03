!> \brief \b mobbrmsd_DLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix associated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLASQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DLASQ2( N, Z, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLASQ2 computes all the eigenvalues of the symmetric positive
!> definite tridiagonal matrix associated with the qd array Z to high
!> relative accuracy are computed to high relative accuracy, in the
!> absence of denormalization, underflow and overflow.
!>
!> To see the relation of Z to the tridiagonal matrix, let L be a
!> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
!> let U be an upper bidiagonal matrix with 1's above and diagonal
!> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
!> symmetric tridiagonal to which it is similar.
!>
!> Note : mobbrmsd_DLASQ2 defines a logical variable, IEEE, which is true
!> on machines which follow ieee-754 floating-point standard in their
!> handling of infinities and NaNs, and false otherwise. This variable
!> is passed to mobbrmsd_DLASQ3.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>        The number of rows and columns in the matrix. N >= 0.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( 4*N )
!>        On entry Z holds the qd array. On exit, entries 1 to N hold
!>        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
!>        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
!>        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
!>        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
!>        shifts that failed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>        = 0: successful exit
!>        < 0: if the i-th argument is a scalar and had an illegal
!>             value, then INFO = -i, if the i-th argument is an
!>             array and the j-entry had an illegal value, then
!>             INFO = -(i*100+j)
!>        > 0: the algorithm failed
!>              = 1, a split was marked by a positive value in E
!>              = 2, current block of Z not diagonalized after 100*N
!>                   iterations (in inner while loop).  On exit Z holds
!>                   a qd array with the same eigenvalues as the given Z.
!>              = 3, termination criterion of outer while loop not met
!>                   (program created more than N unreduced blocks)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Local Variables: I0:N0 defines a current unreduced segment of Z.
!>  The shifts are accumulated in SIGMA. Iteration count is in ITER.
!>  Ping-pong is controlled by PP (alternates between 0 and 1).
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DLASQ2(N, Z, INFO)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     ::  N
  integer, intent(out)    ::  INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: Z(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  logical            :: IEEE
  integer            :: I0, I1, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, &
 &                      K, KMIN, N0, N1, NBIG, NDIV, NFAIL, PP, SPLT, &
 &                      TTYPE
  real(RK)           :: D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN, &
 &                      DN1, DN2, E, EMAX, EMIN, EPS, G, OLDEMN, QMAX, &
 &                      QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL,     &
 &                      TOL2, TRACE, ZMAX, TEMPE, TEMPQ
!     ..
!     .. Parameters ..
  real(RK), parameter :: CBIAS = 1.50_RK
! real(RK), parameter :: ZERO = 0.0_RK
! real(RK), parameter :: HALF = 0.5_RK
! real(RK), parameter :: ONE = 1.0_RK
! real(RK), parameter :: TWO = 2.0_RK
! real(RK), parameter :: FOUR = 4.0_RK
! real(RK), parameter :: HUNDRD = 100.0_RK
!     ..
! interface
!   include 'dlasq3.h'
!   include 'dlasrt.h'
!   include 'xerbla.h'
!   include 'ilaenv.h'
!   include 'dlamch.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic            :: ABS, DBLE, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!     (in case mobbrmsd_DLASQ2 is not called by mobbrmsd_DLASQ1)
!
  INFO = 0
  EPS = mobbrmsd_DLAMCH('Precision')
  SAFMIN = mobbrmsd_DLAMCH('Safe minimum')
  TOL = EPS * HUNDRD
  TOL2 = TOL**2
!
  if (N < 0) then
    INFO = -1
    !CALL XERBLA( 'mobbrmsd_DLASQ2', 1 )
    return
  else if (N == 0) then
    return
  else if (N == 1) then
!
!        1-by-1 case.
!
    if (Z(1) < ZERO) then
      INFO = -201
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
    end if
    return
  else if (N == 2) then
!
!        2-by-2 case.
!
    if (Z(1) < ZERO) then
      INFO = -201
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
      return
    else if (Z(2) < ZERO) then
      INFO = -202
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
      return
    else if (Z(3) < ZERO) then
      INFO = -203
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
      return
    else if (Z(3) > Z(1)) then
      D = Z(3)
      Z(3) = Z(1)
      Z(1) = D
    end if
    Z(5) = Z(1) + Z(2) + Z(3)
    if (Z(2) > Z(3) * TOL2) then
      T = HALF * ((Z(1) - Z(3)) + Z(2))
      S = Z(3) * (Z(2) / T)
      if (S <= T) then
        S = Z(3) * (Z(2) / (T * (ONE + SQRT(ONE + S / T))))
      else
        S = Z(3) * (Z(2) / (T + SQRT(T) * SQRT(T + S)))
      end if
      T = Z(1) + (S + Z(2))
      Z(3) = Z(3) * (Z(1) / T)
      Z(1) = T
    end if
    Z(2) = Z(3)
    Z(6) = Z(2) + Z(1)
    return
  end if
!
!     Check for negative data and compute sums of q's and e's.
!
  Z(2 * N) = ZERO
  EMIN = Z(2)
  QMAX = ZERO
  ZMAX = ZERO
  D = ZERO
  E = ZERO
!
  do K = 1, 2 * (N - 1), 2
    if (Z(K) < ZERO) then
      INFO = -(200 + K)
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
      return
    else if (Z(K + 1) < ZERO) then
      INFO = -(200 + K + 1)
      !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
      return
    end if
    D = D + Z(K)
    E = E + Z(K + 1)
    QMAX = MAX(QMAX, Z(K))
    EMIN = MIN(EMIN, Z(K + 1))
    ZMAX = MAX(QMAX, ZMAX, Z(K + 1))
  end do
  if (Z(2 * N - 1) < ZERO) then
    INFO = -(200 + 2 * N - 1)
    !CALL XERBLA( 'mobbrmsd_DLASQ2', 2 )
    return
  end if
  D = D + Z(2 * N - 1)
  QMAX = MAX(QMAX, Z(2 * N - 1))
  ZMAX = MAX(QMAX, ZMAX)
!
!     Check for diagonality.
!
  if (E == ZERO) then
    do K = 2, N
      Z(K) = Z(2 * K - 1)
    end do
    call mobbrmsd_DLASRT('D', N, Z, IINFO)
    Z(2 * N - 1) = D
    return
  end if
!
  TRACE = D + E
!
!     Check for zero data.
!
  if (TRACE == ZERO) then
    Z(2 * N - 1) = ZERO
    return
  end if
!
!     Check whether the machine is IEEE conformable.
!
  IEEE = (mobbrmsd_ILAENV(10, 'mobbrmsd_DLASQ2', 'N', 1, 2, 3, 4) == 1)
!
!     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
!
  do K = 2 * N, 2, -2
    Z(2 * K) = ZERO
    Z(2 * K - 1) = Z(K)
    Z(2 * K - 2) = ZERO
    Z(2 * K - 3) = Z(K - 1)
  end do
!
  I0 = 1
  N0 = N
!
!     Reverse the qd-array, if warranted.
!
  if (CBIAS * Z(4 * I0 - 3) < Z(4 * N0 - 3)) then
    IPN4 = 4 * (I0 + N0)
    do I4 = 4 * I0, 2 * (I0 + N0 - 1), 4
      TEMP = Z(I4 - 3)
      Z(I4 - 3) = Z(IPN4 - I4 - 3)
      Z(IPN4 - I4 - 3) = TEMP
      TEMP = Z(I4 - 1)
      Z(I4 - 1) = Z(IPN4 - I4 - 5)
      Z(IPN4 - I4 - 5) = TEMP
    end do
  end if
!
!     Initial split checking via dqd and Li's test.
!
  PP = 0
!
  do K = 1, 2
!
    D = Z(4 * N0 + PP - 3)
    do I4 = 4 * (N0 - 1) + PP, 4 * I0 + PP, -4
      if (Z(I4 - 1) <= TOL2 * D) then
        Z(I4 - 1) = -ZERO
        D = Z(I4 - 3)
      else
        D = Z(I4 - 3) * (D / (D + Z(I4 - 1)))
      end if
    end do
!
!        dqd maps Z to ZZ plus Li's test.
!
    EMIN = Z(4 * I0 + PP + 1)
    D = Z(4 * I0 + PP - 3)
    do I4 = 4 * I0 + PP, 4 * (N0 - 1) + PP, 4
      Z(I4 - 2 * PP - 2) = D + Z(I4 - 1)
      if (Z(I4 - 1) <= TOL2 * D) then
        Z(I4 - 1) = -ZERO
        Z(I4 - 2 * PP - 2) = D
        Z(I4 - 2 * PP) = ZERO
        D = Z(I4 + 1)
      else if (SAFMIN * Z(I4 + 1) < Z(I4 - 2 * PP - 2) .and. &
     &         SAFMIN * Z(I4 - 2 * PP - 2) < Z(I4 + 1)) then
        TEMP = Z(I4 + 1) / Z(I4 - 2 * PP - 2)
        Z(I4 - 2 * PP) = Z(I4 - 1) * TEMP
        D = D * TEMP
      else
        Z(I4 - 2 * PP) = Z(I4 + 1) * (Z(I4 - 1) / Z(I4 - 2 * PP - 2))
        D = Z(I4 + 1) * (D / Z(I4 - 2 * PP - 2))
      end if
      EMIN = MIN(EMIN, Z(I4 - 2 * PP))
    end do
    Z(4 * N0 - PP - 2) = D
!
!        Now find qmax.
!
    QMAX = Z(4 * I0 - PP - 2)
    do I4 = 4 * I0 - PP + 2, 4 * N0 - PP - 2, 4
      QMAX = MAX(QMAX, Z(I4))
    end do
!
!        Prepare for the next iteration on K.
!
    PP = 1 - PP
  end do
!
!     Initialise variables to pass to mobbrmsd_DLASQ3.
!
  TTYPE = 0
  DMIN1 = ZERO
  DMIN2 = ZERO
  DN = ZERO
  DN1 = ZERO
  DN2 = ZERO
  G = ZERO
  TAU = ZERO
!
  ITER = 2
  NFAIL = 0
  NDIV = 2 * (N0 - I0)
!
  do IWHILA = 1, N + 1
    if (N0 < 1) GO TO 170
!
!        While array unfinished do
!
!        E(N0) holds the value of SIGMA when submatrix in I0:N0
!        splits from the rest of the array, but is negated.
!
    DESIG = ZERO
    if (N0 == N) then
      SIGMA = ZERO
    else
      SIGMA = -Z(4 * N0 - 1)
    end if
    if (SIGMA < ZERO) then
      INFO = 1
      return
    end if
!
!        Find last unreduced submatrix's top index I0, find QMAX and
!        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
!
    EMAX = ZERO
    if (N0 > I0) then
      EMIN = ABS(Z(4 * N0 - 5))
    else
      EMIN = ZERO
    end if
    QMIN = Z(4 * N0 - 3)
    QMAX = QMIN
    do I4 = 4 * N0, 8, -4
      if (Z(I4 - 5) <= ZERO) GO TO 100
      if (QMIN >= FOUR * EMAX) then
        QMIN = MIN(QMIN, Z(I4 - 3))
        EMAX = MAX(EMAX, Z(I4 - 5))
      end if
      QMAX = MAX(QMAX, Z(I4 - 7) + Z(I4 - 5))
      EMIN = MIN(EMIN, Z(I4 - 5))
    end do
    I4 = 4
!
100 continue
    I0 = I4 / 4
    PP = 0
!
    if (N0 - I0 > 1) then
      DEE = Z(4 * I0 - 3)
      DEEMIN = DEE
      KMIN = I0
      do I4 = 4 * I0 + 1, 4 * N0 - 3, 4
        DEE = Z(I4) * (DEE / (DEE + Z(I4 - 2)))
        if (DEE <= DEEMIN) then
          DEEMIN = DEE
          KMIN = (I4 + 3) / 4
        end if
      end do
      if ((KMIN - I0) * 2 < N0 - KMIN .and. &
     &   DEEMIN <= HALF * Z(4 * N0 - 3)) then
        IPN4 = 4 * (I0 + N0)
        PP = 2
        do I4 = 4 * I0, 2 * (I0 + N0 - 1), 4
          TEMP = Z(I4 - 3)
          Z(I4 - 3) = Z(IPN4 - I4 - 3)
          Z(IPN4 - I4 - 3) = TEMP
          TEMP = Z(I4 - 2)
          Z(I4 - 2) = Z(IPN4 - I4 - 2)
          Z(IPN4 - I4 - 2) = TEMP
          TEMP = Z(I4 - 1)
          Z(I4 - 1) = Z(IPN4 - I4 - 5)
          Z(IPN4 - I4 - 5) = TEMP
          TEMP = Z(I4)
          Z(I4) = Z(IPN4 - I4 - 4)
          Z(IPN4 - I4 - 4) = TEMP
        end do
      end if
    end if
!
!        Put -(initial shift) into DMIN.
!
    DMIN = -MAX(ZERO, QMIN - TWO * SQRT(QMIN) * SQRT(EMAX))
!
!        Now I0:N0 is unreduced.
!        PP = 0 for ping, PP = 1 for pong.
!        PP = 2 indicates that flipping was applied to the Z array and
!               and that the tests for deflation upon entry in mobbrmsd_DLASQ3
!               should not be performed.
!
    NBIG = 100 * (N0 - I0 + 1)
    do IWHILB = 1, NBIG
      if (I0 > N0) GO TO 150
!
!           While submatrix unfinished take a good dqds step.
!
      call mobbrmsd_DLASQ3(I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
     &             ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
     &             DN2, G, TAU)
!
      PP = 1 - PP
!
!           When EMIN is very small check for splits.
!
      if (PP == 0 .and. N0 - I0 >= 3) then
        if (Z(4 * N0) <= TOL2 * QMAX .or. &
       &    Z(4 * N0 - 1) <= TOL2 * SIGMA) then
          SPLT = I0 - 1
          QMAX = Z(4 * I0 - 3)
          EMIN = Z(4 * I0 - 1)
          OLDEMN = Z(4 * I0)
          do I4 = 4 * I0, 4 * (N0 - 3), 4
            if (Z(I4) <= TOL2 * Z(I4 - 3) .or. &
           &    Z(I4 - 1) <= TOL2 * SIGMA) then
              Z(I4 - 1) = -SIGMA
              SPLT = I4 / 4
              QMAX = ZERO
              EMIN = Z(I4 + 3)
              OLDEMN = Z(I4 + 4)
            else
              QMAX = MAX(QMAX, Z(I4 + 1))
              EMIN = MIN(EMIN, Z(I4 - 1))
              OLDEMN = MIN(OLDEMN, Z(I4))
            end if
          end do
          Z(4 * N0 - 1) = EMIN
          Z(4 * N0) = OLDEMN
          I0 = SPLT + 1
        end if
      end if
!
    end do
!
    INFO = 2
!
!        Maximum number of iterations exceeded, restore the shift
!        SIGMA and place the new d's and e's in a qd array.
!        This might need to be done for several blocks
!
    I1 = I0
    N1 = N0
145 continue
    TEMPQ = Z(4 * I0 - 3)
    Z(4 * I0 - 3) = Z(4 * I0 - 3) + SIGMA
    do K = I0 + 1, N0
      TEMPE = Z(4 * K - 5)
      Z(4 * K - 5) = Z(4 * K - 5) * (TEMPQ / Z(4 * K - 7))
      TEMPQ = Z(4 * K - 3)
      Z(4 * K - 3) = Z(4 * K - 3) + SIGMA + TEMPE - Z(4 * K - 5)
    end do
!
!        Prepare to do this on the previous block if there is one
!
    if (I1 > 1) then
      N1 = I1 - 1
      do while ((I1 >= 2) .and. (Z(4 * I1 - 5) >= ZERO))
        I1 = I1 - 1
      end do
      SIGMA = -Z(4 * N1 - 1)
      GO TO 145
    end if

    do K = 1, N
      Z(2 * K - 1) = Z(4 * K - 3)
!
!        Only the block 1..N0 is unfinished.  The rest of the e's
!        must be essentially zero, although sometimes other data
!        has been stored in them.
!
      if (K < N0) then
        Z(2 * K) = Z(4 * K - 1)
      else
        Z(2 * K) = 0
      end if
    end do
    return
!
!        end IWHILB
!
150 continue
!
  end do
!
  INFO = 3
  return
!
!     end IWHILA
!
170 continue
!
!     Move q's to the front.
!
  do K = 2, N
    Z(K) = Z(4 * K - 3)
  end do
!
!     Sort and compute sum of eigenvalues.
!
  call mobbrmsd_DLASRT('D', N, Z, IINFO)
!
  E = ZERO
  do K = N, 1, -1
    E = E + Z(K)
  end do
!
!     Store trace, sum(eigenvalues) and information on performance.
!
  Z(2 * N + 1) = TRACE
  Z(2 * N + 2) = E
  Z(2 * N + 3) = DBLE(ITER)
  Z(2 * N + 4) = DBLE(NDIV) / DBLE(N**2)
  Z(2 * N + 5) = HUNDRD * NFAIL / DBLE(ITER)
  return
!
!     End of mobbrmsd_DLASQ2
!
end subroutine mobbrmsd_DLASQ2
