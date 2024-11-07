!| This program sets problem and machine dependent parameters
!  useful for xHSEQR and related subroutines for eigenvalue
!  problems. It is called whenever
!  mobbrmsd_IPARMQ is called with 12 <= ISPEC <= 16
!
!  Little is known about how best to choose these parameters.
!  It is possible to use different values of the parameters
!  for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!  It is probably best to choose different parameters for
!  different matrices and different parameters at different
!  times during the iteration, but this has not been
!  implemented --- yet.
!
!  The best choices of most of the parameters depend
!  in an ill-understood way on the relative execution
!  rate of xLAQR3 and xLAQR5 and on the nature of each
!  particular eigenvalue problem.  Experiment may be the
!  only practical way to determine which choices are most
!  effective.
!
!  Following is a list of default values supplied by mobbrmsd_IPARMQ.
!  These defaults may be adjusted in order to attain better
!  performance in any particular computational environment.
!
!  mobbrmsd_IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                   Default: 75. (Must be at least 11.)
!
!  mobbrmsd_IPARMQ(ISPEC=13) Recommended deflation window size.
!                   This depends on ILO, IHI and NS, the
!                   number of simultaneous shifts returned
!                   by mobbrmsd_IPARMQ(ISPEC=15).  The default for
!                   (IHI-ILO+1) <= 500 is NS.  The default
!                   for (IHI-ILO+1) > 500 is 3*NS/2.
!
!  mobbrmsd_IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!  mobbrmsd_IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                   a multi-shift QR iteration.
!
!                   If IHI-ILO+1 is ...
!
!                   greater than      ...but less    ... the
!                   or equal to ...      than        default is
!
!                           0               30       NS =   2+
!                          30               60       NS =   4+
!                          60              150       NS =  10
!                         150              590       NS =  **
!                         590             3000       NS =  64
!                        3000             6000       NS = 128
!                        6000             infinity   NS = 256
!
!               (+)  By default matrices of this order are
!                    passed to the implicit double shift routine
!                    xLAHQR.  See mobbrmsd_IPARMQ(ISPEC=12) above.   These
!                    values of NS are used only in case of a rare
!                    xLAHQR failure.
!
!               (**) The asterisks (**) indicate an ad-hoc
!                    function increasing from 10 to 64.
!
!  mobbrmsd_IPARMQ(ISPEC=16) Select structured matrix multiply.
!                   (See ISPEC=16 above for details.)
!                   Default: 3.
!
!  mobbrmsd_IPARMQ(ISPEC=17) Relative cost heuristic for blocksize selection.
!                   Expressed as a percentage.
!                   Default: 10.
!
!  Reference IPARMQ is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure elemental function mobbrmsd_IPARMQ(ISPEC, NAME, OPTS, N, ILO, IHI, LWORK)
  implicit none
  integer, intent(in)      :: ISPEC
!!  ISPEC specifies which tunable parameter mobbrmsd_IPARMQ should return.
!!
!!  - ISPEC=12: (INMIN)  Matrices of order nmin or less
!!              are sent directly to xLAHQR, the implicit
!!              double shift QR algorithm.  NMIN must be
!!              at least 11.
!!
!!  - ISPEC=13: (INWIN)  Size of the deflation window.
!!              This is best set greater than or equal to
!!              the number of simultaneous shifts NS.
!!              Larger matrices benefit from larger deflation
!!              windows.
!!
!!  - ISPEC=14: (INIBL) Determines when to stop nibbling and
!!              invest in an (expensive) multi-shift QR sweep.
!!              If the aggressive early deflation subroutine
!!              finds LD converged eigenvalues from an order
!!              NW deflation window and LD > (NW*NIBBLE)/100,
!!              then the next QR sweep is skipped and early
!!              deflation is applied immediately to the
!!              remaining active diagonal block.  Setting
!!              mobbrmsd_IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!!              multi-shift QR sweep whenever early deflation
!!              finds a converged eigenvalue.  Setting
!!              mobbrmsd_IPARMQ(ISPEC=14) greater than or equal to 100
!!              prevents TTQRE from skipping a multi-shift
!!              QR sweep.
!!
!!  - ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!!              a multi-shift QR iteration.
!!
!!  - ISPEC=16: (IACC22) mobbrmsd_IPARMQ is set to 0, 1 or 2 with the
!!              following meanings.
!!
!!    0.  During the multi-shift QR/QZ sweep,
!!        blocked eigenvalue reordering, blocked
!!        Hessenberg-triangular reduction,
!!        reflections and/or rotations are not
!!        accumulated when updating the
!!        far-from-diagonal matrix entries.
!!
!!    1.  During the multi-shift QR/QZ sweep,
!!        blocked eigenvalue reordering, blocked
!!        Hessenberg-triangular reduction,
!!        reflections and/or rotations are
!!        accumulated, and matrix-matrix
!!        multiplication is used to update the
!!        far-from-diagonal matrix entries.
!!
!!    2.  During the multi-shift QR/QZ sweep,
!!        blocked eigenvalue reordering, blocked
!!        Hessenberg-triangular reduction,
!!        reflections and/or rotations are
!!        accumulated, and 2-by-2 block structure
!!        is exploited during matrix-matrix
!!        multiplies.
!!
!!    (If xTRMM is slower than xGEMM, then
!!    mobbrmsd_IPARMQ(ISPEC=16)=1 may be more efficient than
!!    mobbrmsd_IPARMQ(ISPEC=16)=2 despite the greater level of
!!    arithmetic work implied by the latter choice.)
!!
!!  - ISPEC=17: (ICOST) An estimate of the relative cost of flops
!!              within the near-the-diagonal shift chase compared
!!              to flops within the BLAS calls of a QZ sweep.
!!
  character(*), intent(in) :: NAME
!! Name of the calling subroutine
!!
  character(*), intent(in) :: OPTS
!!  This is a concatenation of the string arguments to TTQRE.
!!
  integer, intent(in)      :: N
!!  N is the order of the Hessenberg matrix H.
!!
  integer, intent(in)      :: ILO
!!  An INTEGER
!!
  integer, intent(in)      :: IHI
!!  It is assumed that H is already upper triangular
!!  in rows and columns 1:ILO-1 and IHI+1:N.
!!
  integer, intent(in)      :: LWORK
!!  The amount of workspace available.
!!
  integer                  :: mobbrmsd_IPARMQ
!! Machine dependent parameters useful for xHSEQR.
!!
!
  integer, parameter  :: INMIN = 12
  integer, parameter  :: INWIN = 13
  integer, parameter  :: INIBL = 14
  integer, parameter  :: ISHFTS = 15
  integer, parameter  :: IACC22 = 16
  integer, parameter  :: ICOST = 17
  integer, parameter  :: NMIN = 75
  integer, parameter  :: K22MIN = 14
  integer, parameter  :: KACMIN = 14
  integer, parameter  :: NIBBLE = 14
  integer, parameter  :: KNWSWP = 500
  integer, parameter  :: RCOST = 10
!     ..
!     .. Local Scalars ..
  integer      :: NH, NS
  integer      :: I, IC, IZ
  character(6) :: SUBNAM
!     ..
!     .. Intrinsic Functions ..
  intrinsic               :: LOG, MAX, MOD, NINT, real
!     ..
!     .. Executable Statements ..
  if ((ISPEC == ISHFTS) .or. (ISPEC == INWIN) .or. (ISPEC == IACC22)) then
!
!        ==== Set the number simultaneous shifts ====
!
    NH = IHI - ILO + 1
    NS = 2
    if (NH >= 30) NS = 4
    if (NH >= 60) NS = 10
    !if (NH >= 150) NS = MAX(10, NH / NINT(LOG(real(NH)) / LOG(TWO)))
    if (NH >= 150) NS = MAX(10, NH / NINT(LOG(real(NH - 2))))
    if (NH >= 590) NS = 64
    if (NH >= 3000) NS = 128
    if (NH >= 6000) NS = 256
    NS = MAX(2, NS - MOD(NS, 2))
  else
    NS = 0
    NH = 0
  end if
!
  if (ISPEC == INMIN) then
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
    mobbrmsd_IPARMQ = NMIN
!
  else if (ISPEC == INIBL) then
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
    mobbrmsd_IPARMQ = NIBBLE
!
  else if (ISPEC == ISHFTS) then
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
    mobbrmsd_IPARMQ = NS
!
  else if (ISPEC == INWIN) then
!
!        ==== NW: deflation window size.  ====
!
    if (NH <= KNWSWP) then
      mobbrmsd_IPARMQ = NS
    else
      mobbrmsd_IPARMQ = 3 * NS / 2
    end if
!
  else if (ISPEC == IACC22) then
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
!
!        Convert NAME to upper case if the first character is lower case.
!
    mobbrmsd_IPARMQ = 0
    SUBNAM = NAME
    IC = ICHAR(SUBNAM(1:1))
    IZ = ICHAR('Z')
    if (IZ == 90 .or. IZ == 122) then
!
!           ASCII character set
!
      if (IC >= 97 .and. IC <= 122) then
        SUBNAM(1:1) = CHAR(IC - 32)
        do I = 2, 6
          IC = ICHAR(SUBNAM(I:I))
          if (IC >= 97 .and. IC <= 122) SUBNAM(I:I) = CHAR(IC - 32)
        end do
      end if
!
    else if (IZ == 233 .or. IZ == 169) then
!
!           EBCDIC character set
!
      if ((IC >= 129 .and. IC <= 137) .or. &
        & (IC >= 145 .and. IC <= 153) .or. &
        & (IC >= 162 .and. IC <= 169)) then
        SUBNAM(1:1) = CHAR(IC + 64)
        do I = 2, 6
          IC = ICHAR(SUBNAM(I:I))
          if ((IC >= 129 .and. IC <= 137) .or. &
            & (IC >= 145 .and. IC <= 153) .or. &
            & (IC >= 162 .and. IC <= 169)) SUBNAM(I:I) = CHAR(IC + 64)
        end do
      end if
!
    else if (IZ == 218 .or. IZ == 250) then
!
!           Prime machines:  ASCII+128
!
      if (IC >= 225 .and. IC <= 250) then
        SUBNAM(1:1) = CHAR(IC - 32)
        do I = 2, 6
          IC = ICHAR(SUBNAM(I:I))
          if (IC >= 225 .and. IC <= 250) SUBNAM(I:I) = CHAR(IC - 32)
        end do
      end if
    end if
!
    if (SUBNAM(2:6) == 'GGHRD' .or. SUBNAM(2:6) == 'GGHD3') then
      mobbrmsd_IPARMQ = 1
      if (NH >= K22MIN) mobbrmsd_IPARMQ = 2
    else if (SUBNAM(4:6) == 'EXC') then
      if (NH >= KACMIN) mobbrmsd_IPARMQ = 1
      if (NH >= K22MIN) mobbrmsd_IPARMQ = 2
    else if (SUBNAM(2:6) == 'HSEQR' .or. SUBNAM(2:5) == 'LAQR') then
      if (NS >= KACMIN) mobbrmsd_IPARMQ = 1
      if (NS >= K22MIN) mobbrmsd_IPARMQ = 2
    end if
!
  else if (ISPEC == ICOST) then
!
!        === Relative cost of near-the-diagonal chase vs
!            BLAS updates ===
!
    mobbrmsd_IPARMQ = RCOST
  else
!        ===== invalid value of ispec =====
    mobbrmsd_IPARMQ = -1
!
  end if
!
!     ==== End of mobbrmsd_IPARMQ ====
!
end function mobbrmsd_IPARMQ

