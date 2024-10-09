!| mobbrmsd_SLASRT sorts numbers in increasing or decreasing order.
!
!  Sort the numbers in \( D \) in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits \( N \) to about 2**32.
!
!  Reference SLASRT is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
pure subroutine mobbrmsd_SLASRT(ID, N, D, INFO)
  implicit none
  character, intent(in)   :: ID
!!  = 'I': sort D in increasing order;
!!
!!  = 'D': sort D in decreasing order.
!!
  integer, intent(in)     :: N
!!  The length of the array D.
!!
  real(RK), intent(inout) :: D(*)
!!  DOUBLE PRECISION array, dimension (N)
!!
!!  On entry, the array to be sorted.
!!
!!  On exit, D has been sorted into increasing order
!!  (D(1) <= ... <= D(N) ) or into decreasing order
!!  (D(1) >= ... >= D(N) ), depending on ID.
!!
  integer, intent(out)    :: INFO
!!  = 0:  successful exit
!!
!!  < 0:  if INFO = -i, the i-th argument had an illegal value
!!
  integer :: DIR, ENDD, I, J, START, STKPNT
  real(RK) :: D1, D2, D3, DMNMX, TMP
  integer :: STACK(2, 32)
  integer, parameter :: SELCT = 20
!..
! interface
!..external Functions..
!   include 'lsame.h'
! end interface
!
! Test the input parameters.
!
  INFO = 0
  DIR = -1
  if (mobbrmsd_LSAME(ID, 'D')) then
    DIR = 0
  else if (mobbrmsd_LSAME(ID, 'I')) then
    DIR = 1
  end if
  if (DIR == -1) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  end if
  if (INFO /= 0) then
!   call XERBLA('SLASRT', -INFO)
    return
  end if
!
! Quick return if possible
!
  if (N <= 1) return
!
  STKPNT = 1
  STACK(1, 1) = 1
  STACK(2, 1) = N
10 continue
  START = STACK(1, STKPNT)
  ENDD = STACK(2, STKPNT)
  STKPNT = STKPNT - 1
  if (ENDD - START <= SELCT .and. ENDD - START > 0) then
!
! do Insertion sort on D(START:ENDD)
!
    if (DIR == 0) then
!
! Sort into decreasing order
!
      LOOP30: do I = START + 1, ENDD
        do J = I, START + 1, -1
          if (D(J) > D(J - 1)) then
            DMNMX = D(J)
            D(J) = D(J - 1)
            D(J - 1) = DMNMX
          else
            exit LOOP30
          end if
        end do
      end do LOOP30
!
    else
!
! Sort into increasing order
!
      LOOP50: do I = START + 1, ENDD
        do J = I, START + 1, -1
          if (D(J) < D(J - 1)) then
            DMNMX = D(J)
            D(J) = D(J - 1)
            D(J - 1) = DMNMX
          else
            exit LOOP50
          end if
        end do
      end do LOOP50
!
    end if
!
  else if (ENDD - START > SELCT) then
!
! Partition D(START:ENDD) and stack parts, largest one first
!
! Choose partition entry as median of 3
!
    D1 = D(START)
    D2 = D(ENDD)
    I = (START + ENDD) / 2
    D3 = D(I)
    if (D1 < D2) then
      if (D3 < D1) then
        DMNMX = D1
      else if (D3 < D2) then
        DMNMX = D3
      else
        DMNMX = D2
      end if
    else
      if (D3 < D2) then
        DMNMX = D2
      else if (D3 < D1) then
        DMNMX = D3
      else
        DMNMX = D1
      end if
    end if
!
    if (DIR == 0) then
!
! Sort into decreasing order
!
      I = START - 1
      J = ENDD + 1
60    continue
70    continue
      J = J - 1
      if (D(J) < DMNMX) GO TO 70
80    continue
      I = I + 1
      if (D(I) > DMNMX) GO TO 80
      if (I < J) then
        TMP = D(I)
        D(I) = D(J)
        D(J) = TMP
        GO TO 60
      end if
      if (J - START > ENDD - J - 1) then
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = START
        STACK(2, STKPNT) = J
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = J + 1
        STACK(2, STKPNT) = ENDD
      else
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = J + 1
        STACK(2, STKPNT) = ENDD
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = START
        STACK(2, STKPNT) = J
      end if
    else
!
! Sort into increasing order
!
      I = START - 1
      J = ENDD + 1
90    continue
100   continue
      J = J - 1
      if (D(J) > DMNMX) GO TO 100
110   continue
      I = I + 1
      if (D(I) < DMNMX) GO TO 110
      if (I < J) then
        TMP = D(I)
        D(I) = D(J)
        D(J) = TMP
        GO TO 90
      end if
      if (J - START > ENDD - J - 1) then
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = START
        STACK(2, STKPNT) = J
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = J + 1
        STACK(2, STKPNT) = ENDD
      else
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = J + 1
        STACK(2, STKPNT) = ENDD
        STKPNT = STKPNT + 1
        STACK(1, STKPNT) = START
        STACK(2, STKPNT) = J
      end if
    end if
  end if
  if (STKPNT > 0) GO TO 10
  return
!
! end of mobbrmsd_SLASRT
!
end

