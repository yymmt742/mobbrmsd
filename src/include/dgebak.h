      pure subroutine DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,&
     &                        INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: JOB, SIDE
      integer,intent(in)      :: IHI, ILO, INFO, LDV, M, N
      real(wp),intent(in)     :: SCALE( * )
      real(wp),intent(inout)  :: V( LDV, * )
      end subroutine DGEBAK
