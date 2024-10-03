!> \brief <b> SGESVD computes the singular value decomposition (SVD) for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGESVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU, JOBVT
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), S( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGESVD computes the singular value decomposition (SVD) of a real
!> M-by-N matrix A, optionally computing the left and/or right singular
!> vectors. The SVD is written
!>
!>      A = U * SIGMA * transpose(V)
!>
!> where SIGMA is an M-by-N matrix which is zero except for its
!> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
!> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
!> are the singular values of A; they are real and non-negative, and
!> are returned in descending order.  The first min(m,n) columns of
!> U and V are the left and right singular vectors of A.
!>
!> Note that the routine returns V**T, not V.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          Specifies options for computing all or part of the matrix U:
!>          = 'A':  all M columns of U are returned in array U:
!>          = 'S':  the first min(m,n) columns of U (the left singular
!>                  vectors) are returned in the array U;
!>          = 'O':  the first min(m,n) columns of U (the left singular
!>                  vectors) are overwritten on the array A;
!>          = 'N':  no columns of U (no left singular vectors) are
!>                  computed.
!> \endverbatim
!>
!> \param[in] JOBVT
!> \verbatim
!>          JOBVT is CHARACTER*1
!>          Specifies options for computing all or part of the matrix
!>          V**T:
!>          = 'A':  all N rows of V**T are returned in the array VT;
!>          = 'S':  the first min(m,n) rows of V**T (the right singular
!>                  vectors) are returned in the array VT;
!>          = 'O':  the first min(m,n) rows of V**T (the right singular
!>                  vectors) are overwritten on the array A;
!>          = 'N':  no rows of V**T (no right singular vectors) are
!>                  computed.
!>
!>          JOBVT and JOBU cannot both be 'O'.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if JOBU = 'O',  A is overwritten with the first min(m,n)
!>                          columns of U (the left singular vectors,
!>                          stored columnwise);
!>          if JOBVT = 'O', A is overwritten with the first min(m,n)
!>                          rows of V**T (the right singular vectors,
!>                          stored rowwise);
!>          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
!>                          are destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (min(M,N))
!>          The singular values of A, sorted so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension (LDU,UCOL)
!>          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!>          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
!>          if JOBU = 'S', U contains the first min(m,n) columns of U
!>          (the left singular vectors, stored columnwise);
!>          if JOBU = 'N' or 'O', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1; if
!>          JOBU = 'S' or 'A', LDU >= M.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,N)
!>          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
!>          V**T;
!>          if JOBVT = 'S', VT contains the first min(m,n) rows of
!>          V**T (the right singular vectors, stored rowwise);
!>          if JOBVT = 'N' or 'O', VT is not referenced.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1; if
!>          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!>          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!>          superdiagonal elements of an upper bidiagonal matrix B
!>          whose diagonal is in S (not necessarily sorted). B
!>          satisfies A = U * B * VT, so it has the same singular values
!>          as A, and singular vectors related by U and VT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
!>             - PATH 1  (M much larger than N, JOBU='N')
!>             - PATH 1t (N much larger than M, JOBVT='N')
!>          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if SBDSQR did not converge, INFO specifies how many
!>                superdiagonals of an intermediate bidiagonal form B
!>                did not converge to zero. See the description of WORK
!>                above for details.
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
!> \date April 2012
!
!> \ingroup realGEsing
!
!  =====================================================================
pure subroutine SGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
  implicit none
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
  character, intent(in) :: JOBU, JOBVT
  integer, intent(in)   :: LDA, LDU, LDVT, LWORK, M, N
  integer, intent(out)  :: INFO
!..
!..Array Arguments..
  real, intent(inout)   :: A(LDA, *)
  real, intent(out)     :: S(*), U(LDU, *), VT(LDVT, *), WORK(*)
!..
!
!  =====================================================================
!
!..Parameters..
  real, parameter :: ZERO = 0.0E0, ONE = 1.0E0
!..
!..Local Scalars..
  logical :: LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS
  logical :: WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS
  integer :: BDSPAC, BLK, CHUNK, I, IE, IERR, IR, ISCL
  integer :: ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU
  integer :: MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU
  integer :: NRVT, WRKBL
  integer :: LWORK_SGEQRF, LWORK_SORGQR_N, LWORK_SORGQR_M
  integer :: LWORK_SGEBRD, LWORK_SORGBR_P, LWORK_SORGBR_Q
  integer :: LWORK_SGELQF, LWORK_SORGLQ_N, LWORK_SORGLQ_M
  real :: ANRM, BIGNUM, EPS, SMLNUM
!..
!..Local Arrays..
  real :: DUM(1)
!..
  interface
! .. External Functions ..
    include 'lsame.h'
    include 'ilaenv.h'
    include 'slamch.h'
    include 'slange.h'
! .. External Subroutines ..
    include 'sgemm.h'
    include 'sbdsqr.h'
    include 'sgebrd.h'
    include 'sgelqf.h'
    include 'sgeqrf.h'
    include 'slacpy.h'
    include 'slascl.h'
    include 'slaset.h'
    include 'sorgbr.h'
    include 'sorglq.h'
    include 'sorgqr.h'
    include 'sormbr.h'
!   include 'xerbla.h'
  end interface
!..
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN, SQRT
!..
!..Executable Statements..
!
!Test the input arguments
!
  INFO = 0
  MINMN = MIN(M, N)
  WNTUA = LSAME(JOBU, 'A')
  WNTUS = LSAME(JOBU, 'S')
  WNTUAS = WNTUA .or. WNTUS
  WNTUO = LSAME(JOBU, 'O')
  WNTUN = LSAME(JOBU, 'N')
  WNTVA = LSAME(JOBVT, 'A')
  WNTVS = LSAME(JOBVT, 'S')
  WNTVAS = WNTVA .or. WNTVS
  WNTVO = LSAME(JOBVT, 'O')
  WNTVN = LSAME(JOBVT, 'N')
  LQUERY = (LWORK == -1)
!
  if (.not. (WNTUA .or. WNTUS .or. WNTUO .or. WNTUN)) then
    INFO = -1
  else if (.not. (WNTVA .or. WNTVS .or. WNTVO .or. WNTVN) .or. (WNTVO .and. WNTUO)) then
    INFO = -2
  else if (M < 0) then
    INFO = -3
  else if (N < 0) then
    INFO = -4
  else if (LDA < MAX(1, M)) then
    INFO = -6
  else if (LDU < 1 .or. (WNTUAS .and. LDU < M)) then
    INFO = -9
  else if (LDVT < 1 .or. (WNTVA .and. LDVT < N) .or. (WNTVS .and. LDVT < MINMN)) then
    INFO = -11
  end if
!
!Compute workspace
!(Note:Comments in the code beginning "Workspace:"describe the
!minimal amount of workspace needed at that point in the code,
!as well as the preferred amount for good performance.
!NB refers to the optimal block size for the immediately
!following subroutine, as returned by ILAENV.)
!
  if (INFO == 0) then
    MINWRK = 1
    MAXWRK = 1
    if (M >= N .and. MINMN > 0) then
!
!Compute space needed for SBDSQR
!
      MNTHR = ILAENV(6, 'SGESVD', JOBU//JOBVT, M, N, 0, 0)
      BDSPAC = 5 * N
!Compute space needed for SGEQRF
      call SGEQRF(M, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SGEQRF = INT(DUM(1))
!Compute space needed for SORGQR
      call SORGQR(M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGQR_N = INT(DUM(1))
      call SORGQR(M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGQR_M = INT(DUM(1))
!Compute space needed for SGEBRD
      call SGEBRD(N, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR)
      LWORK_SGEBRD = INT(DUM(1))
!Compute space needed for SORGBR P
      call SORGBR('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGBR_P = INT(DUM(1))
!Compute space needed for SORGBR Q
      call SORGBR('Q', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGBR_Q = INT(DUM(1))
!
      if (M >= MNTHR) then
        if (WNTUN) then
!
!Path 1(M much larger than N, JOBU='N')
!
          MAXWRK = N + LWORK_SGEQRF
          MAXWRK = MAX(MAXWRK, 3 * N + LWORK_SGEBRD)
          if (WNTVO .or. WNTVAS) MAXWRK = MAX(MAXWRK, 3 * N + LWORK_SORGBR_P)
          MAXWRK = MAX(MAXWRK, BDSPAC)
          MINWRK = MAX(4 * N, BDSPAC)
        else if (WNTUO .and. WNTVN) then
!
!Path 2(M much larger than N, JOBU='O', JOBVT='N')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_N)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = MAX(N * N + WRKBL, N * N + M * N + N)
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUO .and. WNTVAS) then
!
!Path 3(M much larger than N, JOBU='O', JOBVT='S'or
!'A')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_N)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = MAX(N * N + WRKBL, N * N + M * N + N)
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUS .and. WNTVN) then
!
!Path 4(M much larger than N, JOBU='S', JOBVT='N')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_N)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUS .and. WNTVO) then
!
!Path 5(M much larger than N, JOBU='S', JOBVT='O')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_N)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = 2 * N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUS .and. WNTVAS) then
!
!Path 6(M much larger than N, JOBU='S', JOBVT='S'or
!'A')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_N)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUA .and. WNTVN) then
!
!Path 7(M much larger than N, JOBU='A', JOBVT='N')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_M)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUA .and. WNTVO) then
!
!Path 8(M much larger than N, JOBU='A', JOBVT='O')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_M)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = 2 * N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        else if (WNTUA .and. WNTVAS) then
!
!Path 9(M much larger than N, JOBU='A', JOBVT='S'or
!'A')
!
          WRKBL = N + LWORK_SGEQRF
          WRKBL = MAX(WRKBL, N + LWORK_SORGQR_M)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, 3 * N + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = N * N + WRKBL
          MINWRK = MAX(3 * N + M, BDSPAC)
        end if
      else
!
!Path 10(M at least N, but not much larger)
!
        call SGEBRD(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR)
        LWORK_SGEBRD = INT(DUM(1))
        MAXWRK = 3 * N + LWORK_SGEBRD
        if (WNTUS .or. WNTUO) then
          call SORGBR('Q', M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR)
          LWORK_SORGBR_Q = INT(DUM(1))
          MAXWRK = MAX(MAXWRK, 3 * N + LWORK_SORGBR_Q)
        end if
        if (WNTUA) then
          call SORGBR('Q', M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR)
          LWORK_SORGBR_Q = INT(DUM(1))
          MAXWRK = MAX(MAXWRK, 3 * N + LWORK_SORGBR_Q)
        end if
        if (.not. WNTVN) then
          MAXWRK = MAX(MAXWRK, 3 * N + LWORK_SORGBR_P)
        end if
        MAXWRK = MAX(MAXWRK, BDSPAC)
        MINWRK = MAX(3 * N + M, BDSPAC)
      end if
    else if (MINMN > 0) then
!
!Compute space needed for SBDSQR
!
      MNTHR = ILAENV(6, 'SGESVD', JOBU//JOBVT, M, N, 0, 0)
      BDSPAC = 5 * M
!Compute space needed for SGELQF
      call SGELQF(M, N, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SGELQF = INT(DUM(1))
!Compute space needed for SORGLQ
      call SORGLQ(N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGLQ_N = INT(DUM(1))
      call SORGLQ(M, N, M, A, LDA, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGLQ_M = INT(DUM(1))
!Compute space needed for SGEBRD
      call SGEBRD(M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR)
      LWORK_SGEBRD = INT(DUM(1))
!Compute space needed for SORGBR P
      call SORGBR('P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGBR_P = INT(DUM(1))
!Compute space needed for SORGBR Q
      call SORGBR('Q', M, M, M, A, N, DUM(1), DUM(1), -1, IERR)
      LWORK_SORGBR_Q = INT(DUM(1))
      if (N >= MNTHR) then
        if (WNTVN) then
!
!Path 1T(N much larger than M, JOBVT='N')
!
          MAXWRK = M + LWORK_SGELQF
          MAXWRK = MAX(MAXWRK, 3 * M + LWORK_SGEBRD)
          if (WNTUO .or. WNTUAS) MAXWRK = MAX(MAXWRK, 3 * M + LWORK_SORGBR_Q)
          MAXWRK = MAX(MAXWRK, BDSPAC)
          MINWRK = MAX(4 * M, BDSPAC)
        else if (WNTVO .and. WNTUN) then
!
!Path 2T(N much larger than M, JOBU='N', JOBVT='O')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_M)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = MAX(M * M + WRKBL, M * M + M * N + M)
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVO .and. WNTUAS) then
!
!Path 3T(N much larger than M, JOBU='S'or 'A',
!JOBVT = 'O')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_M)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = MAX(M * M + WRKBL, M * M + M * N + M)
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVS .and. WNTUN) then
!
!Path 4T(N much larger than M, JOBU='N', JOBVT='S')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_M)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVS .and. WNTUO) then
!
!Path 5T(N much larger than M, JOBU='O', JOBVT='S')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_M)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = 2 * M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
          MAXWRK = MAX(MAXWRK, MINWRK)
        else if (WNTVS .and. WNTUAS) then
!
!Path 6T(N much larger than M, JOBU='S'or 'A',
!JOBVT = 'S')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_M)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVA .and. WNTUN) then
!
!Path 7T(N much larger than M, JOBU='N', JOBVT='A')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_N)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVA .and. WNTUO) then
!
!Path 8T(N much larger than M, JOBU='O', JOBVT='A')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_N)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = 2 * M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
        else if (WNTVA .and. WNTUAS) then
!
!Path 9T(N much larger than M, JOBU='S'or 'A',
!JOBVT = 'A')
!
          WRKBL = M + LWORK_SGELQF
          WRKBL = MAX(WRKBL, M + LWORK_SORGLQ_N)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SGEBRD)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_P)
          WRKBL = MAX(WRKBL, 3 * M + LWORK_SORGBR_Q)
          WRKBL = MAX(WRKBL, BDSPAC)
          MAXWRK = M * M + WRKBL
          MINWRK = MAX(3 * M + N, BDSPAC)
        end if
      else
!
!Path 10T(N greater than M, but not much larger)
!
        call SGEBRD(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR)
        LWORK_SGEBRD = INT(DUM(1))
        MAXWRK = 3 * M + LWORK_SGEBRD
        if (WNTVS .or. WNTVO) then
!Compute space needed for SORGBR P
          call SORGBR('P', M, N, M, A, N, DUM(1), DUM(1), -1, IERR)
          LWORK_SORGBR_P = INT(DUM(1))
          MAXWRK = MAX(MAXWRK, 3 * M + LWORK_SORGBR_P)
        end if
        if (WNTVA) then
          call SORGBR('P', N, N, M, A, N, DUM(1), &
          &DUM(1), -1, IERR)
          LWORK_SORGBR_P = INT(DUM(1))
          MAXWRK = MAX(MAXWRK, 3 * M + LWORK_SORGBR_P)
        end if
        if (.not. WNTUN) then
          MAXWRK = MAX(MAXWRK, 3 * M + LWORK_SORGBR_Q)
        end if
        MAXWRK = MAX(MAXWRK, BDSPAC)
        MINWRK = MAX(3 * M + N, BDSPAC)
      end if
    end if
    MAXWRK = MAX(MAXWRK, MINWRK)
    WORK(1) = MAXWRK
!
    if (LWORK < MINWRK .and. .not. LQUERY) then
      INFO = -13
    end if
  end if
!
  if (INFO /= 0) then
!   call XERBLA('SGESVD', -INFO)
    return
  else if (LQUERY) then
    return
  end if
!
!Quick return if possible
!
  if (M == 0 .or. N == 0) then
    return
  end if
!
!Get machine constants
!
  EPS = SLAMCH('P')
  SMLNUM = SQRT(SLAMCH('S')) / EPS
  BIGNUM = ONE / SMLNUM
!
!Scale A if max element outside range[SMLNUM, BIGNUM]
!
  call SLANGE('M', M, N, A, LDA, ANRM, DUM)
! ANRM = SLANGE('M', M, N, A, LDA, DUM)
  ISCL = 0
  if (ANRM > ZERO .and. ANRM < SMLNUM) then
    ISCL = 1
    call SLASCL('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR)
  else if (ANRM > BIGNUM) then
    ISCL = 1
    call SLASCL('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR)
  end if
!
  if (M >= N) then
!
!A has at least as many rows as columns.if A has sufficiently
!more rows than columns, first reduce using the QR
!decomposition(if sufficient workspace available)
!
    if (M >= MNTHR) then
!
      if (WNTUN) then
!
!Path 1(M much larger than N, JOBU='N')
!No left singular vectors to be computed
!
        ITAU = 1
        IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need 2 * N, prefer N + N * NB)
!
        call SGEQRF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), &
        &LWORK - IWORK + 1, IERR)
!
!Zero out below R
!
        if (N > 1) then
          call SLASET('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA)
        end if
        IE = 1
        ITAUQ = IE + N
        ITAUP = ITAUQ + N
        IWORK = ITAUP + N
!
!Bidiagonalize R in A
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
        call SGEBRD(N, N, A, LDA, S, WORK(IE), WORK(ITAUQ), &
        &WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
        NCVT = 0
        if (WNTVO .or. WNTVAS) then
!
!if right singular vectors desired, generate P'.
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
          call SORGBR('P', N, N, N, A, LDA, WORK(ITAUP), &
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          NCVT = N
        end if
        IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of A in A if desired
!(Workspace:need BDSPAC)
!
        call SBDSQR('U', N, NCVT, 0, 0, S, WORK(IE), A, LDA, &
        &DUM, 1, DUM, 1, WORK(IWORK), INFO)
!
!if right singular vectors desired in VT, copy them there
!
        if (WNTVAS) call SLACPY('F', N, N, A, LDA, VT, LDVT)
!
      else if (WNTUO .and. WNTVN) then
!
!Path 2(M much larger than N, JOBU='O', JOBVT='N')
!N left singular vectors to be overwritten on A and
!no right singular vectors to be computed
!
        if (LWORK >= N * N + MAX(4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
          IR = 1
          if (LWORK >= MAX(WRKBL, LDA * N + N) + LDA * N) then
!
!WORK(IU) is LDA by N, WORK(IR) is LDA by N
!
            LDWRKU = LDA
            LDWRKR = LDA
          else if (LWORK >= MAX(WRKBL, LDA * N + N) + N * N) then
!
!WORK(IU) is LDA by N, WORK(IR) is N by N
!
            LDWRKU = LDA
            LDWRKR = N
          else
!
!WORK(IU) is LDWRKU by N, WORK(IR) is N by N
!
            LDWRKU = (LWORK - N * N - N) / N
            LDWRKR = N
          end if
          ITAU = IR + LDWRKR * N
          IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
          call SGEQRF(M, N, A, LDA, WORK(ITAU), &
              &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IR) and zero out below it
!
          call SLACPY('U', N, N, A, LDA, WORK(IR), LDWRKR)
          call SLASET('L', N - 1, N - 1, ZERO, ZERO, WORK(IR + 1), LDWRKR)
!
!Generate Q in A
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
          call SORGQR(M, N, N, A, LDA, WORK(ITAU), &
              &       WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
          call SGEBRD(N, N, WORK(IR), LDWRKR, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left vectors bidiagonalizing R
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
          call SORGBR('Q', N, N, N, WORK(IR), LDWRKR,&
          &WORK(ITAUQ), WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
          IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IR)
!(Workspace:need N * N + BDSPAC)
!
          call SBDSQR('U', N, 0, N, 0, S, WORK(IE), DUM, 1,&
          &WORK(IR), LDWRKR, DUM, 1,&
          &WORK(IWORK), INFO)
          IU = IE + N
!
!Multiply Q in A by left singular vectors of R in
!WORK(IR), storing result in WORK(IU) and copying to A
!(Workspace:need N * N + 2 * N, prefer N * N + M * N + N)
!
          do I = 1, M, LDWRKU
            CHUNK = MIN(M - I + 1, LDWRKU)
            call SGEMM('N', 'N', CHUNK, N, N, ONE, A(I, 1),&
            &LDA, WORK(IR), LDWRKR, ZERO,&
            &WORK(IU), LDWRKU)
            call SLACPY('F', CHUNK, N, WORK(IU), LDWRKU, A(I, 1), LDA)
          end do
!
        else
!
!Insufficient workspace for a fast algorithm
!
          IE = 1
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          IWORK = ITAUP + N
!
!Bidiagonalize A
!(Workspace:need 3 * N + M, prefer 3 * N + (M + N) * NB)
!
          call SGEBRD(M, N, A, LDA, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left vectors bidiagonalizing A
!(Workspace:need 4 * N, prefer 3 * N + N * NB)
!
          call SORGBR('Q', M, N, N, A, LDA, WORK(ITAUQ),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in A
!(Workspace:need BDSPAC)
!
          call SBDSQR('U', N, 0, M, 0, S, WORK(IE), DUM, 1,&
          &A, LDA, DUM, 1, WORK(IWORK), INFO)
!
        end if
!
      else if (WNTUO .and. WNTVAS) then
!
!Path 3(M much larger than N, JOBU='O', JOBVT='S'or 'A')
!N left singular vectors to be overwritten on A and
!N right singular vectors to be computed in VT
!
        if (LWORK >= N * N + MAX(4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
          IR = 1
          if (LWORK >= MAX(WRKBL, LDA * N + N) + LDA * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
            LDWRKU = LDA
            LDWRKR = LDA
          else if (LWORK >= MAX(WRKBL, LDA * N + N) + N * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is N by N
!
            LDWRKU = LDA
            LDWRKR = N
          else
!
!WORK(IU) is LDWRKU by N and WORK(IR) is N by N
!
            LDWRKU = (LWORK - N * N - N) / N
            LDWRKR = N
          end if
          ITAU = IR + LDWRKR * N
          IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
          call SGEQRF(M, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to VT, zeroing out below it
!
          call SLACPY('U', N, N, A, LDA, VT, LDVT)
          if (N > 1)&
          &call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
          &VT(2, 1), LDVT)
!
!Generate Q in A
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
          call SORGQR(M, N, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          IWORK = ITAUP + N
!
!Bidiagonalize R in VT, copying result to WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
          call SGEBRD(N, N, VT, LDVT, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          call SLACPY('L', N, N, VT, LDVT, WORK(IR), LDWRKR)
!
!Generate left vectors bidiagonalizing R in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
          call SORGBR('Q', N, N, N, WORK(IR), LDWRKR,&
          &WORK(ITAUQ), WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing R in VT
!(Workspace:need N * N + 4 * N - 1, prefer N * N + 3 * N + (N - 1) * NB)
!
          call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IR) and computing right
!singular vectors of R in VT
!(Workspace:need N * N + BDSPAC)
!
          call SBDSQR('U', N, N, N, 0, S, WORK(IE), VT, LDVT,&
          &WORK(IR), LDWRKR, DUM, 1,&
          &WORK(IWORK), INFO)
          IU = IE + N
!
!Multiply Q in A by left singular vectors of R in
!WORK(IR), storing result in WORK(IU) and copying to A
!(Workspace:need N * N + 2 * N, prefer N * N + M * N + N)
!
          do I = 1, M, LDWRKU
            CHUNK = MIN(M - I + 1, LDWRKU)
            call SGEMM('N', 'N', CHUNK, N, N, ONE, A(I, 1),&
            &LDA, WORK(IR), LDWRKR, ZERO,&
            &WORK(IU), LDWRKU)
            call SLACPY('F', CHUNK, N, WORK(IU), LDWRKU,&
            &A(I, 1), LDA)
          end do
!
        else
!
!Insufficient workspace for a fast algorithm
!
          ITAU = 1
          IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need 2 * N, prefer N + N * NB)
!
          call SGEQRF(M, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to VT, zeroing out below it
!
          call SLACPY('U', N, N, A, LDA, VT, LDVT)
          if (N > 1)&
          &call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
          &VT(2, 1), LDVT)
!
!Generate Q in A
!(Workspace:need 2 * N, prefer N + N * NB)
!
          call SORGQR(M, N, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          IWORK = ITAUP + N
!
!Bidiagonalize R in VT
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
          call SGEBRD(N, N, VT, LDVT, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in A by left vectors bidiagonalizing R
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
          call SORMBR('Q', 'R', 'N', M, N, N, VT, LDVT,&
          &WORK(ITAUQ), A, LDA, WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing R in VT
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
          call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in A and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
          call SBDSQR('U', N, N, M, 0, S, WORK(IE), VT, LDVT,&
          &A, LDA, DUM, 1, WORK(IWORK), INFO)
!
        end if
!
      else if (WNTUS) then
!
        if (WNTVN) then
!
!Path 4(M much larger than N, JOBU='S', JOBVT='N')
!N left singular vectors to be computed in U and
!no right singular vectors to be computed
!
          if (LWORK >= N * N + MAX(4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IR = 1
            if (LWORK >= WRKBL + LDA * N) then
!
!WORK(IR) is LDA by N
!
              LDWRKR = LDA
            else
!
!WORK(IR) is N by N
!
              LDWRKR = N
            end if
            ITAU = IR + LDWRKR * N
            IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IR), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IR), LDWRKR)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
            &WORK(IR + 1), LDWRKR)
!
!Generate Q in A
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SORGQR(M, N, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IR), LDWRKR, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate left vectors bidiagonalizing R in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IR), LDWRKR,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IR)
!(Workspace:need N * N + BDSPAC)
!
            call SBDSQR('U', N, 0, N, 0, S, WORK(IE), DUM,&
            &1, WORK(IR), LDWRKR, DUM, 1,&
            &WORK(IWORK), INFO)
!
!Multiply Q in A by left singular vectors of R in
!WORK(IR), storing result in U
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, A, LDA,&
            &WORK(IR), LDWRKR, ZERO, U, LDU)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SORGQR(M, N, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Zero out below R in A
!
            if (N > 1) then
              call SLASET('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA)
            end if
!
!Bidiagonalize R in A
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, A, LDA, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left vectors bidiagonalizing R
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, A, LDA,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, 0, M, 0, S, WORK(IE), DUM,&
            &1, U, LDU, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        else if (WNTVO) then
!
!Path 5(M much larger than N, JOBU='S', JOBVT='O')
!N left singular vectors to be computed in U and
!N right singular vectors to be overwritten on A
!
          if (LWORK >= 2 * N * N + MAX(4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + 2 * LDA * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
              LDWRKU = LDA
              IR = IU + LDWRKU * N
              LDWRKR = LDA
            else if (LWORK >= WRKBL + (LDA + N) * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is N by N
!
              LDWRKU = LDA
              IR = IU + LDWRKU * N
              LDWRKR = N
            else
!
!WORK(IU) is N by N and WORK(IR) is N by N
!
              LDWRKU = N
              IR = IU + LDWRKU * N
              LDWRKR = N
            end if
            ITAU = IR + LDWRKR * N
            IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need 2 * N * N + 2 * N, prefer 2 * N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IU), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IU), LDWRKU)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
            &WORK(IU + 1), LDWRKU)
!
!Generate Q in A
!(Workspace:need 2 * N * N + 2 * N, prefer 2 * N * N + N + N * NB)
!
            call SORGQR(M, N, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IU), copying result to
!WORK(IR)
!(Workspace:need 2 * N * N + 4 * N,
!prefer 2 * N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IU), LDWRKU, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            call SLACPY('U', N, N, WORK(IU), LDWRKU,&
            &WORK(IR), LDWRKR)
!
!Generate left bidiagonalizing vectors in WORK(IU)
!(Workspace:need 2 * N * N + 4 * N, prefer 2 * N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IU), LDWRKU,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in WORK(IR)
!(Workspace:need 2 * N * N + 4 * N - 1,
!prefer 2 * N * N + 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, WORK(IR), LDWRKR,&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IU) and computing
!right singular vectors of R in WORK(IR)
!(Workspace:need 2 * N * N + BDSPAC)
!
            call SBDSQR('U', N, N, N, 0, S, WORK(IE),&
            &WORK(IR), LDWRKR, WORK(IU),&
            &LDWRKU, DUM, 1, WORK(IWORK), INFO)
!
!Multiply Q in A by left singular vectors of R in
!WORK(IU), storing result in U
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, A, LDA, WORK(IU), LDWRKU, ZERO, U, LDU)
!
!Copy right singular vectors of R to A
!(Workspace:need N * N)
!
            call SLACPY('F', N, N, WORK(IR), LDWRKR, A, LDA)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SORGQR(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Zero out below R in A
!
            if (N > 1) then
              call SLASET('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA)
            end if
!
!Bidiagonalize R in A
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, A, LDA, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left vectors bidiagonalizing R
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, A, LDA,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing R in A
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, A, LDA, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in A
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, N, M, 0, S, WORK(IE), A,&
            &LDA, U, LDU, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        else if (WNTVAS) then
!
!Path 6(M much larger than N, JOBU='S', JOBVT='S'
!or 'A')
!N left singular vectors to be computed in U and
!N right singular vectors to be computed in VT
!
          if (LWORK >= N * N + MAX(4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + LDA * N) then
!
!WORK(IU) is LDA by N
!
              LDWRKU = LDA
            else
!
!WORK(IU) is N by N
!
              LDWRKU = N
            end if
            ITAU = IU + LDWRKU * N
            IWORK = ITAU + N
!
!Compute A = Q * R
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IU), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IU), LDWRKU)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO, WORK(IU + 1), LDWRKU)
!
!Generate Q in A
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SORGQR(M, N, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IU), copying result to VT
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IU), LDWRKU, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            call SLACPY('U', N, N, WORK(IU), LDWRKU, VT, LDVT)
!
!Generate left bidiagonalizing vectors in WORK(IU)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IU), LDWRKU,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in VT
!(Workspace:need N * N + 4 * N - 1,
!prefer N * N + 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IU) and computing
!right singular vectors of R in VT
!(Workspace:need N * N + BDSPAC)
!
            call SBDSQR('U', N, N, N, 0, S, WORK(IE), VT,&
            &LDVT, WORK(IU), LDWRKU, DUM, 1,&
            &WORK(IWORK), INFO)
!
!Multiply Q in A by left singular vectors of R in
!WORK(IU), storing result in U
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, A, LDA, WORK(IU), LDWRKU, ZERO, U, LDU)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SORGQR(M, N, N, U, LDU, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to VT, zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, VT, LDVT)
            if (N > 1) call SLASET('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in VT
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, VT, LDVT, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left bidiagonalizing vectors
!in VT
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, VT, LDVT,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in VT
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, N, M, 0, S, WORK(IE), VT,&
            &LDVT, U, LDU, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        end if
!
      else if (WNTUA) then
!
        if (WNTVN) then
!
!Path 7(M much larger than N, JOBU='A', JOBVT='N')
!M left singular vectors to be computed in U and
!no right singular vectors to be computed
!
          if (LWORK >= N * N + MAX(N + M, 4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IR = 1
            if (LWORK >= WRKBL + LDA * N) then
!
!WORK(IR) is LDA by N
!
              LDWRKR = LDA
            else
!
!WORK(IR) is N by N
!
              LDWRKR = N
            end if
            ITAU = IR + LDWRKR * N
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Copy R to WORK(IR), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IR),&
            &LDWRKR)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
            &WORK(IR + 1), LDWRKR)
!
!Generate Q in U
!(Workspace:need N * N + N + M, prefer N * N + N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IR), LDWRKR, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in WORK(IR)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IR), LDWRKR,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IR)
!(Workspace:need N * N + BDSPAC)
!
            call SBDSQR('U', N, 0, N, 0, S, WORK(IE), DUM,&
            &1, WORK(IR), LDWRKR, DUM, 1,&
            &WORK(IWORK), INFO)
!
!Multiply Q in U by left singular vectors of R in
!WORK(IR), storing result in A
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, U, LDU,&
            &WORK(IR), LDWRKR, ZERO, A, LDA)
!
!Copy left singular vectors of A from A to U
!
            call SLACPY('F', M, N, A, LDA, U, LDU)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need N + M, prefer N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Zero out below R in A
!
            if (N > 1) then
              call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
              &A(2, 1), LDA)
            end if
!
!Bidiagonalize R in A
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, A, LDA, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left bidiagonalizing vectors
!in A
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, A, LDA,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, 0, M, 0, S, WORK(IE), DUM,&
            &1, U, LDU, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        else if (WNTVO) then
!
!Path 8(M much larger than N, JOBU='A', JOBVT='O')
!M left singular vectors to be computed in U and
!N right singular vectors to be overwritten on A
!
          if (LWORK >= 2 * N * N + MAX(N + M, 4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + 2 * LDA * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
              LDWRKU = LDA
              IR = IU + LDWRKU * N
              LDWRKR = LDA
            else if (LWORK >= WRKBL + (LDA + N) * N) then
!
!WORK(IU) is LDA by N and WORK(IR) is N by N
!
              LDWRKU = LDA
              IR = IU + LDWRKU * N
              LDWRKR = N
            else
!
!WORK(IU) is N by N and WORK(IR) is N by N
!
              LDWRKU = N
              IR = IU + LDWRKU * N
              LDWRKR = N
            end if
            ITAU = IR + LDWRKR * N
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N * N + 2 * N, prefer 2 * N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need 2 * N * N + N + M, prefer 2 * N * N + N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IU), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IU),&
            &LDWRKU)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
            &WORK(IU + 1), LDWRKU)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IU), copying result to
!WORK(IR)
!(Workspace:need 2 * N * N + 4 * N,
!prefer 2 * N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IU), LDWRKU, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            call SLACPY('U', N, N, WORK(IU), LDWRKU,&
            &WORK(IR), LDWRKR)
!
!Generate left bidiagonalizing vectors in WORK(IU)
!(Workspace:need 2 * N * N + 4 * N, prefer 2 * N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IU), LDWRKU,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in WORK(IR)
!(Workspace:need 2 * N * N + 4 * N - 1,
!prefer 2 * N * N + 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, WORK(IR), LDWRKR,&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IU) and computing
!right singular vectors of R in WORK(IR)
!(Workspace:need 2 * N * N + BDSPAC)
!
            call SBDSQR('U', N, N, N, 0, S, WORK(IE),&
            &WORK(IR), LDWRKR, WORK(IU),&
            &LDWRKU, DUM, 1, WORK(IWORK), INFO)
!
!Multiply Q in U by left singular vectors of R in
!WORK(IU), storing result in A
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, U, LDU,&
            &WORK(IU), LDWRKU, ZERO, A, LDA)
!
!Copy left singular vectors of A from A to U
!
            call SLACPY('F', M, N, A, LDA, U, LDU)
!
!Copy right singular vectors of R from WORK(IR) to A
!
            call SLACPY('F', N, N, WORK(IR), LDWRKR, A, LDA)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need N + M, prefer N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Zero out below R in A
!
            if (N > 1) then
              call SLASET('L', N - 1, N - 1, ZERO, ZERO, A(2, 1), LDA)
            end if
!
!Bidiagonalize R in A
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, A, LDA, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left bidiagonalizing vectors
!in A
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, A, LDA,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in A
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, A, LDA, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in A
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, N, M, 0, S, WORK(IE), A,&
            &LDA, U, LDU, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        else if (WNTVAS) then
!
!Path 9(M much larger than N, JOBU='A', JOBVT='S'
!or 'A')
!M left singular vectors to be computed in U and
!N right singular vectors to be computed in VT
!
          if (LWORK >= N * N + MAX(N + M, 4 * N, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + LDA * N) then
!
!WORK(IU) is LDA by N
!
              LDWRKU = LDA
            else
!
!WORK(IU) is N by N
!
              LDWRKU = N
            end if
            ITAU = IU + LDWRKU * N
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need N * N + 2 * N, prefer N * N + N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need N * N + N + M, prefer N * N + N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R to WORK(IU), zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, WORK(IU),&
            &LDWRKU)
            call SLASET('L', N - 1, N - 1, ZERO, ZERO,&
            &WORK(IU + 1), LDWRKU)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in WORK(IU), copying result to VT
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, WORK(IU), LDWRKU, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            call SLACPY('U', N, N, WORK(IU), LDWRKU, VT, LDVT)
!
!Generate left bidiagonalizing vectors in WORK(IU)
!(Workspace:need N * N + 4 * N, prefer N * N + 3 * N + N * NB)
!
            call SORGBR('Q', N, N, N, WORK(IU), LDWRKU,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in VT
!(Workspace:need N * N + 4 * N - 1,
!prefer N * N + 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of R in WORK(IU) and computing
!right singular vectors of R in VT
!(Workspace:need N * N + BDSPAC)
!
            call SBDSQR('U', N, N, N, 0, S, WORK(IE), VT,&
            &LDVT, WORK(IU), LDWRKU, DUM, 1,&
            &WORK(IWORK), INFO)
!
!Multiply Q in U by left singular vectors of R in
!WORK(IU), storing result in A
!(Workspace:need N * N)
!
            call SGEMM('N', 'N', M, N, N, ONE, U, LDU,&
            &WORK(IU), LDWRKU, ZERO, A, LDA)
!
!Copy left singular vectors of A from A to U
!
            call SLACPY('F', M, N, A, LDA, U, LDU)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + N
!
!Compute A = Q * R, copying result to U
!(Workspace:need 2 * N, prefer N + N * NB)
!
            call SGEQRF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, N, A, LDA, U, LDU)
!
!Generate Q in U
!(Workspace:need N + M, prefer N + M * NB)
!
            call SORGQR(M, M, N, U, LDU, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy R from A to VT, zeroing out below it
!
            call SLACPY('U', N, N, A, LDA, VT, LDVT)
            if (N > 1) call SLASET('L', N - 1, N - 1, ZERO, ZERO, VT(2, 1), LDVT)
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
!
!Bidiagonalize R in VT
!(Workspace:need 4 * N, prefer 3 * N + 2 * N * NB)
!
            call SGEBRD(N, N, VT, LDVT, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply Q in U by left bidiagonalizing vectors
!in VT
!(Workspace:need 3 * N + M, prefer 3 * N + M * NB)
!
            call SORMBR('Q', 'R', 'N', M, N, N, VT, LDVT,&
            &WORK(ITAUQ), U, LDU, WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in VT
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
            call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + N
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', N, N, M, 0, S, WORK(IE), VT,&
            &LDVT, U, LDU, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        end if
!
      end if
!
    else
!
!M < MNTHR
!
!Path 10(M at least N, but not much larger)
!Reduce to bidiagonal form without QR decomposition
!
      IE = 1
      ITAUQ = IE + N
      ITAUP = ITAUQ + N
      IWORK = ITAUP + N
!
!Bidiagonalize A
!(Workspace:need 3 * N + M, prefer 3 * N + (M + N) * NB)
!
      call SGEBRD(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ),&
      &WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,&
      &IERR)
      if (WNTUAS) then
!
!if left singular vectors desired in U, copy result to U
!and generate left bidiagonalizing vectors in U
!(Workspace:need 3 * N + NCU, prefer 3 * N + NCU * NB)
!
        call SLACPY('L', M, N, A, LDA, U, LDU)
        if (WNTUS) NCU = N
        if (WNTUA) NCU = M
        call SORGBR('Q', M, NCU, N, U, LDU, WORK(ITAUQ),&
        &WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTVAS) then
!
!if right singular vectors desired in VT, copy result to
!VT and generate right bidiagonalizing vectors in VT
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
        call SLACPY('U', N, N, A, LDA, VT, LDVT)
        call SORGBR('P', N, N, N, VT, LDVT, WORK(ITAUP),&
        &WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTUO) then
!
!if left singular vectors desired in A, generate left
!bidiagonalizing vectors in A
!(Workspace:need 4 * N, prefer 3 * N + N * NB)
!
        call SORGBR('Q', M, N, N, A, LDA, WORK(ITAUQ),&
        &WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTVO) then
!
!if right singular vectors desired in A, generate right
!bidiagonalizing vectors in A
!(Workspace:need 4 * N - 1, prefer 3 * N + (N - 1) * NB)
!
        call SORGBR('P', N, N, N, A, LDA, WORK(ITAUP),&
        &WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      IWORK = IE + N
      if (WNTUAS .or. WNTUO) NRU = M
      if (WNTUN) NRU = 0
      if (WNTVAS .or. WNTVO) NCVT = N
      if (WNTVN) NCVT = 0
      if ((.not. WNTUO) .and. (.not. WNTVO)) then
!
!Perform bidiagonal QR iteration, if desired, computing
!left singular vectors in U and computing right singular
!vectors in VT
!(Workspace:need BDSPAC)
!
        call SBDSQR('U', N, NCVT, NRU, 0, S, WORK(IE), VT,&
        &LDVT, U, LDU, DUM, 1, WORK(IWORK), INFO)
      else if ((.not. WNTUO) .and. WNTVO) then
!
!Perform bidiagonal QR iteration, if desired, computing
!left singular vectors in U and computing right singular
!vectors in A
!(Workspace:need BDSPAC)
!
        call SBDSQR('U', N, NCVT, NRU, 0, S, WORK(IE), A, LDA,&
        &U, LDU, DUM, 1, WORK(IWORK), INFO)
      else
!
!Perform bidiagonal QR iteration, if desired, computing
!left singular vectors in A and computing right singular
!vectors in VT
!(Workspace:need BDSPAC)
!
        call SBDSQR('U', N, NCVT, NRU, 0, S, WORK(IE), VT,&
        &LDVT, A, LDA, DUM, 1, WORK(IWORK), INFO)
      end if
!
    end if
!
  else
!
!A has more columns than rows.if A has sufficiently more
!columns than rows, first reduce using the LQ decomposition(if
!sufficient workspace available)
!
    if (N >= MNTHR) then
!
      if (WNTVN) then
!
!Path 1T(N much larger than M, JOBVT='N')
!No right singular vectors to be computed
!
        ITAU = 1
        IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need 2 * M, prefer M + M * NB)
!
        call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK),&
        &LWORK - IWORK + 1, IERR)
!
!Zero out above L
!
        call SLASET('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA)
        IE = 1
        ITAUQ = IE + M
        ITAUP = ITAUQ + M
        IWORK = ITAUP + M
!
!Bidiagonalize L in A
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
        call SGEBRD(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ),&
        &WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1,&
        &IERR)
        if (WNTUO .or. WNTUAS) then
!
!if left singular vectors desired, generate Q
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
          call SORGBR('Q', M, M, M, A, LDA, WORK(ITAUQ),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
        end if
        IWORK = IE + M
        NRU = 0
        if (WNTUO .or. WNTUAS) NRU = M
!
!Perform bidiagonal QR iteration, computing left singular
!vectors of A in A if desired
!(Workspace:need BDSPAC)
!
        call SBDSQR('U', M, 0, NRU, 0, S, WORK(IE), DUM, 1, A,&
        &LDA, DUM, 1, WORK(IWORK), INFO)
!
!if left singular vectors desired in U, copy them there
!
        if (WNTUAS) call SLACPY('F', M, M, A, LDA, U, LDU)
!
      else if (WNTVO .and. WNTUN) then
!
!Path 2T(N much larger than M, JOBU='N', JOBVT='O')
!M right singular vectors to be overwritten on A and
!no left singular vectors to be computed
!
        if (LWORK >= M * M + MAX(4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
          IR = 1
          if (LWORK >= MAX(WRKBL, LDA * N + M) + LDA * M) then
!
!WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
            LDWRKU = LDA
            CHUNK = N
            LDWRKR = LDA
          else if (LWORK >= MAX(WRKBL, LDA * N + M) + M * M) then
!
!WORK(IU) is LDA by N and WORK(IR) is M by M
!
            LDWRKU = LDA
            CHUNK = N
            LDWRKR = M
          else
!
!WORK(IU) is M by CHUNK and WORK(IR) is M by M
!
            LDWRKU = M
            CHUNK = (LWORK - M * M - M) / M
            LDWRKR = M
          end if
          ITAU = IR + LDWRKR * M
          IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
          call SGELQF(M, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IR) and zero out above it
!
          call SLACPY('L', M, M, A, LDA, WORK(IR), LDWRKR)
          call SLASET('U', M - 1, M - 1, ZERO, ZERO,&
          &WORK(IR + LDWRKR), LDWRKR)
!
!Generate Q in A
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
          call SORGLQ(M, N, M, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IR)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
          call SGEBRD(M, M, WORK(IR), LDWRKR, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing L
!(Workspace:need M * M + 4 * M - 1, prefer M * M + 3 * M + (M - 1) * NB)
!
          call SORGBR('P', M, M, M, WORK(IR), LDWRKR,&
          &WORK(ITAUP), WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
          IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of L in WORK(IR)
!(Workspace:need M * M + BDSPAC)
!
          call SBDSQR('U', M, M, 0, 0, S, WORK(IE),&
          &WORK(IR), LDWRKR, DUM, 1, DUM, 1,&
          &WORK(IWORK), INFO)
          IU = IE + M
!
!Multiply right singular vectors of L in WORK(IR) by Q
!in A, storing result in WORK(IU) and copying to A
!(Workspace:need M * M + 2 * M, prefer M * M + M * N + M)
!
          do I = 1, N, CHUNK
            BLK = MIN(N - I + 1, CHUNK)
            call SGEMM('N', 'N', M, BLK, M, ONE, WORK(IR),&
            &LDWRKR, A(1, I), LDA, ZERO,&
            &WORK(IU), LDWRKU)
            call SLACPY('F', M, BLK, WORK(IU), LDWRKU, A(1, I), LDA)
          end do
!
        else
!
!Insufficient workspace for a fast algorithm
!
          IE = 1
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          IWORK = ITAUP + M
!
!Bidiagonalize A
!(Workspace:need 3 * M + N, prefer 3 * M + (M + N) * NB)
!
          call SGEBRD(M, N, A, LDA, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing A
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
          call SORGBR('P', M, N, M, A, LDA, WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of A in A
!(Workspace:need BDSPAC)
!
          call SBDSQR('L', M, N, 0, 0, S, WORK(IE), A, LDA,&
          &DUM, 1, DUM, 1, WORK(IWORK), INFO)
!
        end if
!
      else if (WNTVO .and. WNTUAS) then
!
!Path 3T(N much larger than M, JOBU='S'or 'A', JOBVT='O')
!M right singular vectors to be overwritten on A and
!M left singular vectors to be computed in U
!
        if (LWORK >= M * M + MAX(4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
          IR = 1
          if (LWORK >= MAX(WRKBL, LDA * N + M) + LDA * M) then
!
!WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
            LDWRKU = LDA
            CHUNK = N
            LDWRKR = LDA
          else if (LWORK >= MAX(WRKBL, LDA * N + M) + M * M) then
!
!WORK(IU) is LDA by N and WORK(IR) is M by M
!
            LDWRKU = LDA
            CHUNK = N
            LDWRKR = M
          else
!
!WORK(IU) is M by CHUNK and WORK(IR) is M by M
!
            LDWRKU = M
            CHUNK = (LWORK - M * M - M) / M
            LDWRKR = M
          end if
          ITAU = IR + LDWRKR * M
          IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
          call SGELQF(M, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to U, zeroing about above it
!
          call SLACPY('L', M, M, A, LDA, U, LDU)
          call SLASET('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU)
!
!Generate Q in A
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
          call SORGLQ(M, N, M, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          IWORK = ITAUP + M
!
!Bidiagonalize L in U, copying result to WORK(IR)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
          call SGEBRD(M, M, U, LDU, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          call SLACPY('U', M, M, U, LDU, WORK(IR), LDWRKR)
!
!Generate right vectors bidiagonalizing L in WORK(IR)
!(Workspace:need M * M + 4 * M - 1, prefer M * M + 3 * M + (M - 1) * NB)
!
          call SORGBR('P', M, M, M, WORK(IR), LDWRKR,&
          &WORK(ITAUP), WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
!
!Generate left vectors bidiagonalizing L in U
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + M * NB)
!
          call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of L in U, and computing right
!singular vectors of L in WORK(IR)
!(Workspace:need M * M + BDSPAC)
!
          call SBDSQR('U', M, M, M, 0, S, WORK(IE),&
          &WORK(IR), LDWRKR, U, LDU, DUM, 1,&
          &WORK(IWORK), INFO)
          IU = IE + M
!
!Multiply right singular vectors of L in WORK(IR) by Q
!in A, storing result in WORK(IU) and copying to A
!(Workspace:need M * M + 2 * M, prefer M * M + M * N + M))
!
          do I = 1, N, CHUNK
            BLK = MIN(N - I + 1, CHUNK)
            call SGEMM('N', 'N', M, BLK, M, ONE, WORK(IR),&
            &LDWRKR, A(1, I), LDA, ZERO,&
            &WORK(IU), LDWRKU)
            call SLACPY('F', M, BLK, WORK(IU), LDWRKU,&
            &A(1, I), LDA)
          end do
!
        else
!
!Insufficient workspace for a fast algorithm
!
          ITAU = 1
          IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need 2 * M, prefer M + M * NB)
!
          call SGELQF(M, N, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to U, zeroing out above it
!
          call SLACPY('L', M, M, A, LDA, U, LDU)
          call SLASET('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU)
!
!Generate Q in A
!(Workspace:need 2 * M, prefer M + M * NB)
!
          call SORGLQ(M, N, M, A, LDA, WORK(ITAU),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IE = ITAU
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          IWORK = ITAUP + M
!
!Bidiagonalize L in U
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
          call SGEBRD(M, M, U, LDU, S, WORK(IE),&
          &WORK(ITAUQ), WORK(ITAUP),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right vectors bidiagonalizing L by Q in A
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
          call SORMBR('P', 'L', 'T', M, N, M, U, LDU,&
          &WORK(ITAUP), A, LDA, WORK(IWORK),&
          &LWORK - IWORK + 1, IERR)
!
!Generate left vectors bidiagonalizing L in U
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
          call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ),&
          &WORK(IWORK), LWORK - IWORK + 1, IERR)
          IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in A
!(Workspace:need BDSPAC)
!
          call SBDSQR('U', M, N, M, 0, S, WORK(IE), A, LDA,&
          &U, LDU, DUM, 1, WORK(IWORK), INFO)
!
        end if
!
      else if (WNTVS) then
!
        if (WNTUN) then
!
!Path 4T(N much larger than M, JOBU='N', JOBVT='S')
!M right singular vectors to be computed in VT and
!no left singular vectors to be computed
!
          if (LWORK >= M * M + MAX(4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IR = 1
            if (LWORK >= WRKBL + LDA * M) then
!
!WORK(IR) is LDA by M
!
              LDWRKR = LDA
            else
!
!WORK(IR) is M by M
!
              LDWRKR = M
            end if
            ITAU = IR + LDWRKR * M
            IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IR), zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, WORK(IR), LDWRKR)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO,&
            &WORK(IR + LDWRKR), LDWRKR)
!
!Generate Q in A
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SORGLQ(M, N, M, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IR)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IR), LDWRKR, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate right vectors bidiagonalizing L in
!WORK(IR)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IR), LDWRKR,&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of L in WORK(IR)
!(Workspace:need M * M + BDSPAC)
!
            call SBDSQR('U', M, M, 0, 0, S, WORK(IE),&
            &WORK(IR), LDWRKR, DUM, 1, DUM, 1,&
            &WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IR) by
!Q in A, storing result in VT
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IR),&
            &LDWRKR, A, LDA, ZERO, VT, LDVT)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy result to VT
!
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SORGLQ(M, N, M, VT, LDVT, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Zero out above L in A
!
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA)
!
!Bidiagonalize L in A
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, A, LDA, S, WORK(IE),&
            &WORK(ITAUQ), WORK(ITAUP),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right vectors bidiagonalizing L by Q in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, A, LDA,&
            &WORK(ITAUP), VT, LDVT,&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, 0, 0, S, WORK(IE), VT,&
            &LDVT, DUM, 1, DUM, 1, WORK(IWORK),&
            &INFO)
!
          end if
!
        else if (WNTUO) then
!
!Path 5T(N much larger than M, JOBU='O', JOBVT='S')
!M right singular vectors to be computed in VT and
!M left singular vectors to be overwritten on A
!
          if (LWORK >= 2 * M * M + MAX(4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + 2 * LDA * M) then
!
!WORK(IU) is LDA by M and WORK(IR) is LDA by M
!
              LDWRKU = LDA
              IR = IU + LDWRKU * M
              LDWRKR = LDA
            else if (LWORK >= WRKBL + (LDA + M) * M) then
!
!WORK(IU) is LDA by M and WORK(IR) is M by M
!
              LDWRKU = LDA
              IR = IU + LDWRKU * M
              LDWRKR = M
            else
!
!WORK(IU) is M by M and WORK(IR) is M by M
!
              LDWRKU = M
              IR = IU + LDWRKU * M
              LDWRKR = M
            end if
            ITAU = IR + LDWRKR * M
            IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need 2 * M * M + 2 * M, prefer 2 * M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IU), zeroing out below it
!
            call SLACPY('L', M, M, A, LDA, WORK(IU), LDWRKU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO,&
            &WORK(IU + LDWRKU), LDWRKU)
!
!Generate Q in A
!(Workspace:need 2 * M * M + 2 * M, prefer 2 * M * M + M + M * NB)
!
            call SORGLQ(M, N, M, A, LDA, WORK(ITAU),&
            &WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IU), copying result to
!WORK(IR)
!(Workspace:need 2 * M * M + 4 * M,
!prefer 2 * M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IU), LDWRKU, S,&
            &WORK(IE), WORK(ITAUQ),&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, M, WORK(IU), LDWRKU,&
            &WORK(IR), LDWRKR)
!
!Generate right bidiagonalizing vectors in WORK(IU)
!(Workspace:need 2 * M * M + 4 * M - 1,
!prefer 2 * M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IU), LDWRKU,&
            &WORK(ITAUP), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in WORK(IR)
!(Workspace:need 2 * M * M + 4 * M, prefer 2 * M * M + 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, WORK(IR), LDWRKR,&
            &WORK(ITAUQ), WORK(IWORK),&
            &LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of L in WORK(IR) and computing
!right singular vectors of L in WORK(IU)
!(Workspace:need 2 * M * M + BDSPAC)
!
            call SBDSQR('U', M, M, M, 0, S, WORK(IE),&
            &WORK(IU), LDWRKU, WORK(IR),&
            &LDWRKR, DUM, 1, WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IU) by
!Q in A, storing result in VT
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IU),&
            &LDWRKU, A, LDA, ZERO, VT, LDVT)
!
!Copy left singular vectors of L to A
!(Workspace:need M * M)
!
            call SLACPY('F', M, M, WORK(IR), LDWRKR, A, LDA)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SORGLQ(M, N, M, VT, LDVT, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Zero out above L in A
!
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA)
!
!Bidiagonalize L in A
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, A, LDA, S, WORK(IE), &
                &       WORK(ITAUQ), WORK(ITAUP), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right vectors bidiagonalizing L by Q in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, A, LDA, &
                &       WORK(ITAUP), VT, LDVT, &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors of L in A
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, A, LDA, WORK(ITAUQ), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, compute left
!singular vectors of A in A and compute right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, M, 0, S, WORK(IE), VT, &
                &       LDVT, A, LDA, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        else if (WNTUAS) then
!
!Path 6T(N much larger than M, JOBU='S'or 'A',
!JOBVT = 'S')
!M right singular vectors to be computed in VT and
!M left singular vectors to be computed in U
!
          if (LWORK >= M * M + MAX(4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + LDA * M) then
!
!WORK(IU) is LDA by N
!
              LDWRKU = LDA
            else
!
!WORK(IU) is LDA by M
!
              LDWRKU = M
            end if
            ITAU = IU + LDWRKU * M
            IWORK = ITAU + M
!
!Compute A = L * Q
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IU), zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, WORK(IU), LDWRKU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, &
                &       WORK(IU + LDWRKU), LDWRKU)
!
!Generate Q in A
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SORGLQ(M, N, M, A, LDA, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IU), copying result to U
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IU), LDWRKU, S, &
                &       WORK(IE), WORK(ITAUQ), &
                &       WORK(ITAUP), WORK(IWORK), &
                &       LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, M, WORK(IU), LDWRKU, U, LDU)
!
!Generate right bidiagonalizing vectors in WORK(IU)
!(Workspace:need M * M + 4 * M - 1,
!prefer M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IU), LDWRKU, &
                &       WORK(ITAUP), WORK(IWORK), &
                &       LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in U
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ), &
                &        WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of L in U and computing right
!singular vectors of L in WORK(IU)
!(Workspace:need M * M + BDSPAC)
!
            call SBDSQR('U', M, M, M, 0, S, WORK(IE), &
                &       WORK(IU), LDWRKU, U, LDU, DUM, 1, &
                &       WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IU) by
!Q in A, storing result in VT
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IU), &
                &      LDWRKU, A, LDA, ZERO, VT, LDVT)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SORGLQ(M, N, M, VT, LDVT, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to U, zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, U, LDU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in U
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, U, LDU, S, WORK(IE), &
                &       WORK(ITAUQ), WORK(ITAUP), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right bidiagonalizing vectors in U by Q
!in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, U, LDU, &
                &       WORK(ITAUP), VT, LDVT, &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in U
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, M, 0, S, WORK(IE), VT, &
                &       LDVT, U, LDU, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        end if
!
      else if (WNTVA) then
!
        if (WNTUN) then
!
!Path 7T(N much larger than M, JOBU='N', JOBVT='A')
!N right singular vectors to be computed in VT and
!no left singular vectors to be computed
!
          if (LWORK >= M * M + MAX(N + M, 4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IR = 1
            if (LWORK >= WRKBL + LDA * M) then
!
!WORK(IR) is LDA by M
!
              LDWRKR = LDA
            else
!
!WORK(IR) is M by M
!
              LDWRKR = M
            end if
            ITAU = IR + LDWRKR * M
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Copy L to WORK(IR), zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, WORK(IR), LDWRKR)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, WORK(IR + LDWRKR), LDWRKR)
!
!Generate Q in VT
!(Workspace:need M * M + M + N, prefer M * M + M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IR)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IR), LDWRKR, S, &
                &       WORK(IE), WORK(ITAUQ), &
                &       WORK(ITAUP), WORK(IWORK), &
                &       LWORK - IWORK + 1, IERR)
!
!Generate right bidiagonalizing vectors in WORK(IR)
!(Workspace:need M * M + 4 * M - 1,
!prefer M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IR), LDWRKR, &
                &       WORK(ITAUP), WORK(IWORK), &
                &       LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of L in WORK(IR)
!(Workspace:need M * M + BDSPAC)
!
            call SBDSQR('U', M, M, 0, 0, S, WORK(IE), &
                &       WORK(IR), LDWRKR, DUM, 1, DUM, 1, &
                &       WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IR) by
!Q in VT, storing result in A
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IR), &
                &      LDWRKR, VT, LDVT, ZERO, A, LDA)
!
!Copy right singular vectors of A from A to VT
!
            call SLACPY('F', M, N, A, LDA, VT, LDVT)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need M + N, prefer M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Zero out above L in A
!
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA)
!
!Bidiagonalize L in A
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right bidiagonalizing vectors in A by Q
!in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, A, LDA, &
                &       WORK(ITAUP), VT, LDVT, &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, 0, 0, S, WORK(IE), VT, &
                &       LDVT, DUM, 1, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        else if (WNTUO) then
!
!Path 8T(N much larger than M, JOBU='O', JOBVT='A')
!N right singular vectors to be computed in VT and
!M left singular vectors to be overwritten on A
!
          if (LWORK >= 2 * M * M + MAX(N + M, 4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + 2 * LDA * M) then
!
!WORK(IU) is LDA by M and WORK(IR) is LDA by M
!
              LDWRKU = LDA
              IR = IU + LDWRKU * M
              LDWRKR = LDA
            else if (LWORK >= WRKBL + (LDA + M) * M) then
!
!WORK(IU) is LDA by M and WORK(IR) is M by M
!
              LDWRKU = LDA
              IR = IU + LDWRKU * M
              LDWRKR = M
            else
!
!WORK(IU) is M by M and WORK(IR) is M by M
!
              LDWRKU = M
              IR = IU + LDWRKU * M
              LDWRKR = M
            end if
            ITAU = IR + LDWRKR * M
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M * M + 2 * M, prefer 2 * M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need 2 * M * M + M + N, prefer 2 * M * M + M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IU), zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, WORK(IU), LDWRKU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, WORK(IU + LDWRKU), LDWRKU)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IU), copying result to
!WORK(IR)
!(Workspace:need 2 * M * M + 4 * M,
!prefer 2 * M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IU), LDWRKU, S, WORK(IE), WORK(ITAUQ), &
                &       WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, M, WORK(IU), LDWRKU, WORK(IR), LDWRKR)
!
!Generate right bidiagonalizing vectors in WORK(IU)
!(Workspace:need 2 * M * M + 4 * M - 1,
!prefer 2 * M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IU), LDWRKU, WORK(ITAUP), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in WORK(IR)
!(Workspace:need 2 * M * M + 4 * M, prefer 2 * M * M + 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, WORK(IR), LDWRKR, WORK(ITAUQ), &
                &       WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of L in WORK(IR) and computing
!right singular vectors of L in WORK(IU)
!(Workspace:need 2 * M * M + BDSPAC)
!
            call SBDSQR('U', M, M, M, 0, S, WORK(IE), WORK(IU), LDWRKU, &
                &       WORK(IR), LDWRKR, DUM, 1, WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IU) by
!Q in VT, storing result in A
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IU), LDWRKU, VT, LDVT, ZERO, A, LDA)
!
!Copy right singular vectors of A from A to VT
!
            call SLACPY('F', M, N, A, LDA, VT, LDVT)
!
!Copy left singular vectors of A from WORK(IR) to A
!
            call SLACPY('F', M, M, WORK(IR), LDWRKR, A, LDA)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need M + N, prefer M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Zero out above L in A
!
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, A(1, 2), LDA)
!
!Bidiagonalize L in A
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right bidiagonalizing vectors in A by Q
!in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, A, LDA, WORK(ITAUP), VT, LDVT, WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in A
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, A, LDA, WORK(ITAUQ), WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in A and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, M, 0, S, WORK(IE), VT, LDVT, A, LDA, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        else if (WNTUAS) then
!
!Path 9T(N much larger than M, JOBU='S'or 'A',
!JOBVT = 'A')
!N right singular vectors to be computed in VT and
!M left singular vectors to be computed in U
!
          if (LWORK >= M * M + MAX(N + M, 4 * M, BDSPAC)) then
!
!Sufficient workspace for a fast algorithm
!
            IU = 1
            if (LWORK >= WRKBL + LDA * M) then
!
!WORK(IU) is LDA by M
!
              LDWRKU = LDA
            else
!
!WORK(IU) is M by M
!
              LDWRKU = M
            end if
            ITAU = IU + LDWRKU * M
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need M * M + 2 * M, prefer M * M + M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need M * M + M + N, prefer M * M + M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to WORK(IU), zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, WORK(IU), LDWRKU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, WORK(IU + LDWRKU), LDWRKU)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in WORK(IU), copying result to U
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, WORK(IU), LDWRKU, S, WORK(IE), WORK(ITAUQ), &
                &       WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('L', M, M, WORK(IU), LDWRKU, U, LDU)
!
!Generate right bidiagonalizing vectors in WORK(IU)
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + (M - 1) * NB)
!
            call SORGBR('P', M, M, M, WORK(IU), LDWRKU, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in U
!(Workspace:need M * M + 4 * M, prefer M * M + 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of L in U and computing right
!singular vectors of L in WORK(IU)
!(Workspace:need M * M + BDSPAC)
!
            call SBDSQR('U', M, M, M, 0, S, WORK(IE), WORK(IU), LDWRKU, U, LDU, DUM, 1, WORK(IWORK), INFO)
!
!Multiply right singular vectors of L in WORK(IU) by
!Q in VT, storing result in A
!(Workspace:need M * M)
!
            call SGEMM('N', 'N', M, N, M, ONE, WORK(IU), LDWRKU, VT, LDVT, ZERO, A, LDA)
!
!Copy right singular vectors of A from A to VT
!
            call SLACPY('F', M, N, A, LDA, VT, LDVT)
!
          else
!
!Insufficient workspace for a fast algorithm
!
            ITAU = 1
            IWORK = ITAU + M
!
!Compute A = L * Q, copying result to VT
!(Workspace:need 2 * M, prefer M + M * NB)
!
            call SGELQF(M, N, A, LDA, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
            call SLACPY('U', M, N, A, LDA, VT, LDVT)
!
!Generate Q in VT
!(Workspace:need M + N, prefer M + N * NB)
!
            call SORGLQ(N, N, M, VT, LDVT, WORK(ITAU), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Copy L to U, zeroing out above it
!
            call SLACPY('L', M, M, A, LDA, U, LDU)
            call SLASET('U', M - 1, M - 1, ZERO, ZERO, U(1, 2), LDU)
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!
!Bidiagonalize L in U
!(Workspace:need 4 * M, prefer 3 * M + 2 * M * NB)
!
            call SGEBRD(M, M, U, LDU, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Multiply right bidiagonalizing vectors in U by Q
!in VT
!(Workspace:need 3 * M + N, prefer 3 * M + N * NB)
!
            call SORMBR('P', 'L', 'T', M, N, M, U, LDU, WORK(ITAUP), VT, LDVT, WORK(IWORK), LWORK - IWORK + 1, IERR)
!
!Generate left bidiagonalizing vectors in U
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
            call SORGBR('Q', M, M, M, U, LDU, WORK(ITAUQ), WORK(IWORK), LWORK - IWORK + 1, IERR)
            IWORK = IE + M
!
!Perform bidiagonal QR iteration, computing left
!singular vectors of A in U and computing right
!singular vectors of A in VT
!(Workspace:need BDSPAC)
!
            call SBDSQR('U', M, N, M, 0, S, WORK(IE), VT, LDVT, U, LDU, DUM, 1, WORK(IWORK), INFO)
!
          end if
!
        end if
!
      end if
!
    else
!
!N < MNTHR
!
!Path 10T(N greater than M, but not much larger)
!Reduce to bidiagonal form without LQ decomposition
!
      IE = 1
      ITAUQ = IE + M
      ITAUP = ITAUQ + M
      IWORK = ITAUP + M
!
!Bidiagonalize A
!(Workspace:need 3 * M + N, prefer 3 * M + (M + N) * NB)
!
      call SGEBRD(M, N, A, LDA, S, WORK(IE), WORK(ITAUQ), WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
      if (WNTUAS) then
!
!if left singular vectors desired in U, copy result to U
!and generate left bidiagonalizing vectors in U
!(Workspace:need 4 * M - 1, prefer 3 * M + (M - 1) * NB)
!
        call SLACPY('L', M, M, A, LDA, U, LDU)
        call SORGBR('Q', M, M, N, U, LDU, WORK(ITAUQ), WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTVAS) then
!
!if right singular vectors desired in VT, copy result to
!VT and generate right bidiagonalizing vectors in VT
!(Workspace:need 3 * M + NRVT, prefer 3 * M + NRVT * NB)
!
        call SLACPY('U', M, N, A, LDA, VT, LDVT)
        if (WNTVA) NRVT = N
        if (WNTVS) NRVT = M
        call SORGBR('P', NRVT, N, M, VT, LDVT, WORK(ITAUP), &
       & WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTUO) then
!
!if left singular vectors desired in A, generate left
!bidiagonalizing vectors in A
!(Workspace:need 4 * M - 1, prefer 3 * M + (M - 1) * NB)
!
        call SORGBR('Q', M, M, N, A, LDA, WORK(ITAUQ), &
        & WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      if (WNTVO) then
!
!if right singular vectors desired in A, generate right
!bidiagonalizing vectors in A
!(Workspace:need 4 * M, prefer 3 * M + M * NB)
!
        call SORGBR('P', M, N, M, A, LDA, WORK(ITAUP), WORK(IWORK), LWORK - IWORK + 1, IERR)
      end if
      IWORK = IE + M
      if (WNTUAS .or. WNTUO) NRU = M
      if (WNTUN) NRU = 0
      if (WNTVAS .or. WNTVO) NCVT = N
      if (WNTVN) NCVT = 0
      if ((.not. WNTUO) .and. (.not. WNTVO)) then
!
! Perform bidiagonal QR iteration, if desired, computing
! left singular vectors in U and computing right singular
! vectors in VT
! (Workspace:need BDSPAC)
!
        call SBDSQR('L', M, NCVT, NRU, 0, S, WORK(IE), VT, &
        & LDVT, U, LDU, DUM, 1, WORK(IWORK), INFO)
      else if ((.not. WNTUO) .and. WNTVO) then
!
! Perform bidiagonal QR iteration, if desired, computing
! left singular vectors in U and computing right singular
! vectors in A
! (Workspace:need BDSPAC)
!
        call SBDSQR('L', M, NCVT, NRU, 0, S, WORK(IE), A, LDA, &
        & U, LDU, DUM, 1, WORK(IWORK), INFO)
      else
!
! Perform bidiagonal QR iteration, if desired, computing
! left singular vectors in A and computing right singular
! vectors in VT
! (Workspace:need BDSPAC)
!
        call SBDSQR('L', M, NCVT, NRU, 0, S, WORK(IE), VT, LDVT, A, LDA, DUM, 1, WORK(IWORK), INFO)
      end if
!
    end if
!
  end if
!
! if SBDSQR failed to converge, copy unconverged superdiagonals
! to WORK(2:MINMN)
!
  if (INFO /= 0) then
    if (IE > 2) then
      do I = 1, MINMN - 1
        WORK(I + 1) = WORK(I + IE - 1)
      end do
    end if
    if (IE < 2) then
      do I = MINMN - 1, 1, -1
        WORK(I + 1) = WORK(I + IE - 1)
      end do
    end if
  end if
!
! Undo scaling if necessary
!
  if (ISCL == 1) then
    if (ANRM > BIGNUM) call SLASCL('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR)
    if (INFO /= 0 .and. ANRM > BIGNUM) call SLASCL('G', 0, 0, BIGNUM, ANRM, MINMN - 1, 1, WORK(2), MINMN, IERR)
    if (ANRM < SMLNUM) call SLASCL('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR)
    if (INFO /= 0 .and. ANRM < SMLNUM) call SLASCL('G', 0, 0, SMLNUM, ANRM, MINMN - 1, 1, WORK(2), MINMN, IERR)
  end if
!
! return optimal workspace in WORK(1)
!
  WORK(1) = MAXWRK
!
  return
!
! end of SGESVD
!
end
