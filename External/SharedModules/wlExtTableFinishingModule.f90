MODULE wlExtTableFinishingModule

  USE wlKindModule, ONLY: dp 
  USE wlInterpolationModule

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

    SUBROUTINE HoleCharacterize( fail, LinearOK, Table, iMinGradient )

      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: fail
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
      LOGICAL, DIMENSION(:,:,:,:), INTENT(out) :: LinearOK 
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: iMinGradient

      INTEGER :: i, j, k, idim 
      INTEGER :: isize, jsize, ksize 
      REAL(dp) :: Gradient
      REAL(dp) :: MinGradient

      LinearOK = .false.
      iMinGradient = 0

      isize = SIZE(fail, DIM=1)
      jsize = SIZE(fail, DIM=2)
      ksize = SIZE(fail, DIM=3)

      DO k = 1, ksize
        DO j = 1, jsize
          DO i = 1, isize 
            
            IF ( .not.fail(i,j,k) ) CYCLE
      
            MinGradient = 1.d99

            IF ( i == 1 .or. i == isize ) THEN 
            ELSE
              IF ( .not.fail(i-1,j,k) .and. .not.fail(i+1,j,k) ) LinearOK(1,i,j,k) = .true. 
            END IF
            IF ( LinearOK(1,i,j,k) ) THEN 
              Gradient = ABS(Table(i+1,j,k) - Table(i-1,j,k))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 1
              END IF  
            END IF

            IF ( j == 1 .or. j == jsize ) THEN 
            ELSE
              IF ( .not.fail(i,j-1,k) .and. .not.fail(i,j+1,k) ) LinearOK(2,i,j,k) = .true. 
            END IF
            IF ( LinearOK(2,i,j,k) ) THEN 
              Gradient = ABS(Table(i,j+1,k) - Table(i,j-1,k))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 2
              END IF  
            END IF

            IF ( k == 1 .or. k == ksize ) THEN 
            ELSE
              IF ( .not.fail(i,j,k-1) .and. .not.fail(i,j,k+1) ) LinearOK(3,i,j,k) = .true. 
            END IF
            IF ( LinearOK(3,i,j,k) ) THEN 
              Gradient = ABS(Table(i,j,k+1) - Table(i,j,k-1))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 3
              END IF  
            END IF

          END DO
        END DO
      END DO

    END SUBROUTINE HoleCharacterize

END MODULE wlExtTableFinishingModule
