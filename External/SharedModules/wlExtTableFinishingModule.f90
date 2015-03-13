MODULE wlExtTableFinishingModule

  USE wlKindModule, ONLY: dp 

  implicit none

  CONTAINS

    SUBROUTINE LoneCellLocate( fail, LoneCell )

      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: fail
      LOGICAL, DIMENSION(:,:,:), INTENT(out) :: LoneCell

      INTEGER :: i, j, k
      LOGICAL, DIMENSION(3,3,3) :: CornerMask

      CornerMask = .false.
      DO k = 1, 3, 2
        DO j = 1, 3, 2
          DO i = 1, 3, 2
            CornerMask(i,j,k) = .true.
          END DO
        END DO
      END DO
      
      LoneCell = .false. 
      DO k = 2, SIZE(fail, DIM=3) - 1 
        DO j = 2, SIZE(fail, DIM=2) - 1
          DO i = 2, SIZE(fail, DIM=1) - 1

            IF ( .not.fail(i,j,k) ) CYCLE 
            IF ( ANY( fail(i-1:i+1,j-1:j+1,k-1:k+1) .and. CornerMask ) ) CYCLE
            LoneCell(i,j,k) = .true.

            WRITE(*,*) i,j,k 

          END DO
        END DO
      END DO

    END SUBROUTINE LoneCellLocate

!    SUBROUTINE LoneCellInterpolate( LoneCell, Table )
        
!    END SUBROUTINE LoneCellInterpolate

END MODULE wlExtTableFinishingModule
