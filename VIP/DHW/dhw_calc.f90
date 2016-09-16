
  SUBROUTINE dhw_calc()

  IMPLICIT NONE

  INTEGER,INTENT(in) :: levstart,nx,ny,nz
  REAL(4),DIMENSION(nz,ny,nx),INTENT(in)  :: tabin
  REAL(4),DIMENSION(nz,ny,nx),INTENT(out) :: tabout

  INTEGER :: ji,jj,jk,jkstart

  tabout(:,:,:) = tabin(:,:,:)

  ! convert C to F index
  jkstart = levstart + 1

  DO ji=1,nx
    DO jj=1,ny
      DO jk=jkstart,nz



         ! Read in temp from all times 

         IF ( tabout(jk,jj,ji) > tabout(jk-1,jj,ji) ) THEN
              tabout(jk,jj,ji) = tabout(jk-1,jj,ji)
         ENDIF

      ENDDO
    ENDDO
  ENDDO
            
  END SUBROUTINE
