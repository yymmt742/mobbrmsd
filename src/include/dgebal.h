      pure subroutine DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: JOB
      integer,intent(in)      :: LDA, N
      integer,intent(out)     :: IHI, ILO, INFO
      real(wp),intent(out)    :: SCALE( * )
      real(wp),intent(inout)  :: A( LDA, * )
      end subroutine DGEBAL
