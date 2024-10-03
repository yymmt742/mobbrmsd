!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAENV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee infinity and NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 17:
!>               xHSEQR or related subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
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
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
pure elemental function ILAENV(ISPEC, NAME, OPTS, N1, N2, N3, N4)
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: NAME, OPTS
  integer, intent(in)      :: ISPEC, N1, N2, N3, N4
  integer                  :: ILAENV
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  integer                  :: I, IC, IZ, NB, NBMIN, NX
  logical                  :: CNAME, SNAME, TWOSTAGE
  character(1)             :: C1
  character(2)             :: C2, C4
  character(3)             :: C3
  character(16)            :: SUBNAM * 16
!     ..
!     .. Intrinsic Functions ..
  intrinsic                :: CHAR, ICHAR, INT, MIN, real
!     ..
!     .. External Functions ..
  interface
    include 'ieeeck.h'
    include 'iparmq.h'
  end interface
!     ..
!     .. Executable Statements ..
!
!     GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
!    &        130, 140, 150, 160, 160, 160, 160, 160, 160) ISPEC
!
! Invalid value for ISPEC
!
  select case (ISPEC)
  case (1, 2, 3)
!   10 CONTINUE
!
!   Convert NAME to upper case if the first character is lower case.
!
    ILAENV = 1
    SUBNAM = NAME
    IC = ICHAR(SUBNAM(1:1))
    IZ = ICHAR('Z')
    if (IZ == 90 .or. IZ == 122) then
!
!          ASCII character set
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
!          EBCDIC character set
!
      if ((IC >= 129 .and. IC <= 137) .or. &
     &    (IC >= 145 .and. IC <= 153) .or. &
     &    (IC >= 162 .and. IC <= 169)) then
        SUBNAM(1:1) = CHAR(IC + 64)
        do I = 2, 6
          IC = ICHAR(SUBNAM(I:I))
          if ((IC >= 129 .and. IC <= 137) .or. &
         &    (IC >= 145 .and. IC <= 153) .or. &
         &    (IC >= 162 .and. IC <= 169)) SUBNAM(I:I) = CHAR(IC + 64)
        end do
      end if
!
    else if (IZ == 218 .or. IZ == 250) then
!
!     Prime machines:  ASCII+128
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
    C1 = SUBNAM(1:1)
    SNAME = C1 == 'S' .or. C1 == 'D'
    CNAME = C1 == 'C' .or. C1 == 'Z'
    if (.not. (CNAME .or. SNAME)) return
    C2 = SUBNAM(2:3)
    C3 = SUBNAM(4:6)
    C4 = C3(2:3)
    TWOSTAGE = LEN(SUBNAM) >= 11 .and. SUBNAM(11:11) == '2'
!
!   GO TO ( 50, 60, 70 ) ISPEC
    select case (ISPEC)
    case (1)
!
!     50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      if (SUBNAM(2:6) == 'LAORH') then
!
!            This is for *LAORHR_GETRFNP routine
!
        if (SNAME) then
          NB = 32
        else
          NB = 32
        end if
      else if (C2 == 'GE') then
        if (C3 == 'TRF') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        else if (C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or. C3 == 'QLF') then
          if (SNAME) then
            NB = 32
          else
            NB = 32
          end if
        else if (C3 == 'QR ') then
          if (N3 == 1) then
            if (SNAME) then
!             M*N
              if ((N1 * N2 <= 131072) .or. (N1 <= 8192)) then
                NB = N1
              else
                NB = 32768 / N2
              end if
            else
              if ((N1 * N2 <= 131072) .or. (N1 <= 8192)) then
                NB = N1
              else
                NB = 32768 / N2
              end if
            end if
          else
            if (SNAME) then
              NB = 1
            else
              NB = 1
            end if
          end if
        else if (C3 == 'LQ ') then
          if (N3 == 2) then
            if (SNAME) then
!         M*N
              if ((N1 * N2 <= 131072) .or. (N1 <= 8192)) then
                NB = N1
              else
                NB = 32768 / N2
              end if
            else
              if ((N1 * N2 <= 131072) .or. (N1 <= 8192)) then
                NB = N1
              else
                NB = 32768 / N2
              end if
            end if
          else
            if (SNAME) then
              NB = 1
            else
              NB = 1
            end if
          end if
        else if (C3 == 'HRD') then
          if (SNAME) then
            NB = 32
          else
            NB = 32
          end if
        else if (C3 == 'BRD') then
          if (SNAME) then
            NB = 32
          else
            NB = 32
          end if
        else if (C3 == 'TRI') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        end if
      else if (C2 == 'PO') then
        if (C3 == 'TRF') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        end if
      else if (C2 == 'SY') then
        if (C3 == 'TRF') then
          if (SNAME) then
            if (TWOSTAGE) then
              NB = 192
            else
              NB = 64
            end if
          else
            if (TWOSTAGE) then
              NB = 192
            else
              NB = 64
            end if
          end if
        else if (SNAME .and. C3 == 'TRD') then
          NB = 32
        else if (SNAME .and. C3 == 'GST') then
          NB = 64
        end if
      else if (CNAME .and. C2 == 'HE') then
        if (C3 == 'TRF') then
          if (TWOSTAGE) then
            NB = 192
          else
            NB = 64
          end if
        else if (C3 == 'TRD') then
          NB = 32
        else if (C3 == 'GST') then
          NB = 64
        end if
      else if (SNAME .and. C2 == 'OR') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NB = 32
          end if
        else if (C3(1:1) == 'M') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NB = 32
          end if
        end if
      else if (CNAME .and. C2 == 'UN') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NB = 32
          end if
        else if (C3(1:1) == 'M') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NB = 32
          end if
        end if
      else if (C2 == 'GB') then
        if (C3 == 'TRF') then
          if (SNAME) then
            if (N4 <= 64) then
              NB = 1
            else
              NB = 32
            end if
          else
            if (N4 <= 64) then
              NB = 1
            else
              NB = 32
            end if
          end if
        end if
      else if (C2 == 'PB') then
        if (C3 == 'TRF') then
          if (SNAME) then
            if (N2 <= 64) then
              NB = 1
            else
              NB = 32
            end if
          else
            if (N2 <= 64) then
              NB = 1
            else
              NB = 32
            end if
          end if
        end if
      else if (C2 == 'TR') then
        if (C3 == 'TRI') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        else if (C3 == 'EVC') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        end if
      else if (C2 == 'LA') then
        if (C3 == 'UUM') then
          if (SNAME) then
            NB = 64
          else
            NB = 64
          end if
        end if
      else if (SNAME .and. C2 == 'ST') then
        if (C3 == 'EBZ') then
          NB = 1
        end if
      else if (C2 == 'GG') then
        NB = 32
        if (C3 == 'HD3') then
          if (SNAME) then
            NB = 32
          else
            NB = 32
          end if
        end if
      end if
      ILAENV = NB
      return
!
    case (2)
!       60 CONTINUE
!
!         ISPEC = 2:  minimum block size
!
      NBMIN = 2
      if (C2 == 'GE') then
        if (C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or. C3 == 'QLF') then
          if (SNAME) then
            NBMIN = 2
          else
            NBMIN = 2
          end if
        else if (C3 == 'HRD') then
          if (SNAME) then
            NBMIN = 2
          else
            NBMIN = 2
          end if
        else if (C3 == 'BRD') then
          if (SNAME) then
            NBMIN = 2
          else
            NBMIN = 2
          end if
        else if (C3 == 'TRI') then
          if (SNAME) then
            NBMIN = 2
          else
            NBMIN = 2
          end if
        end if
      else if (C2 == 'SY') then
        if (C3 == 'TRF') then
          if (SNAME) then
            NBMIN = 8
          else
            NBMIN = 8
          end if
        else if (SNAME .and. C3 == 'TRD') then
          NBMIN = 2
        end if
      else if (CNAME .and. C2 == 'HE') then
        if (C3 == 'TRD') then
          NBMIN = 2
        end if
      else if (SNAME .and. C2 == 'OR') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NBMIN = 2
          end if
        else if (C3(1:1) == 'M') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NBMIN = 2
          end if
        end if
      else if (CNAME .and. C2 == 'UN') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NBMIN = 2
          end if
        else if (C3(1:1) == 'M') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NBMIN = 2
          end if
        end if
      else if (C2 == 'GG') then
        NBMIN = 2
        if (C3 == 'HD3') then
          NBMIN = 2
        end if
      end if
      ILAENV = NBMIN
      return
!
!   70 CONTINUE
    case (3)
!
!         ISPEC = 3:  crossover point
!
      NX = 0
      if (C2 == 'GE') then
        if (C3 == 'QRF' .or. C3 == 'RQF' .or. C3 == 'LQF' .or. C3 == 'QLF') then
          if (SNAME) then
            NX = 128
          else
            NX = 128
          end if
        else if (C3 == 'HRD') then
          if (SNAME) then
            NX = 128
          else
            NX = 128
          end if
        else if (C3 == 'BRD') then
          if (SNAME) then
            NX = 128
          else
            NX = 128
          end if
        end if
      else if (C2 == 'SY') then
        if (SNAME .and. C3 == 'TRD') then
          NX = 32
        end if
      else if (CNAME .and. C2 == 'HE') then
        if (C3 == 'TRD') then
          NX = 32
        end if
      else if (SNAME .and. C2 == 'OR') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NX = 128
          end if
        end if
      else if (CNAME .and. C2 == 'UN') then
        if (C3(1:1) == 'G') then
          if (C4 == 'QR' .or. C4 == 'RQ' .or. C4 == 'LQ' .or. C4 == &
   &          'QL' .or. C4 == 'HR' .or. C4 == 'TR' .or. C4 == 'BR') &
   &           then
            NX = 128
          end if
        end if
      else if (C2 == 'GG') then
        NX = 128
        if (C3 == 'HD3') then
          NX = 128
        end if
      end if
      ILAENV = NX
      return
    end select
!
!   80 CONTINUE
  case (4)
!
!       ISPEC = 4:  number of shifts (used by xHSEQR)
!
    ILAENV = 6
    return
!
!   90 CONTINUE
  case (5)
!
!       ISPEC = 5:  minimum column dimension (not used)
!
    ILAENV = 2
    return
!
!  100 CONTINUE
  case (6)
!
!       ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
    ILAENV = INT(real(MIN(N1, N2)) * 1.6E0)
    return
!
! 110 CONTINUE
  case (7)
!
!       ISPEC = 7:  number of processors (not used)
!
    ILAENV = 1
    return
!
! 120 CONTINUE
  case (8)
!
!       ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
    ILAENV = 50
    return
!
! 130 CONTINUE
  case (9)
!
!       ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                   computation tree in the divide-and-conquer algorithm
!                   (used by xGELSD and xGESDD)
!
    ILAENV = 25
    return
!
! 140 CONTINUE
  case (10)
!
!       ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
!
!       ILAENV = 0
    ILAENV = 1
    if (ILAENV == 1) then
      ILAENV = IEEECK(1, 0.0, 1.0)
    end if
    return
!
! 150 CONTINUE
  case (11)
!
!       ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
!
!       ILAENV = 0
    ILAENV = 1
    if (ILAENV == 1) then
      ILAENV = IEEECK(0, 0.0, 1.0)
    end if
    return
!
! 160 CONTINUE
  case (12, 13, 14, 15, 16, 17)
!
!       12 <= ISPEC <= 17: xHSEQR or related subroutines.
!
    ILAENV = IPARMQ(ISPEC, NAME, OPTS, N1, N2, N3, N4)
    return
!
  case default
!
    ILAENV = -1
    return
!
  end select
!
! End of ILAENV
!
end function ILAENV
