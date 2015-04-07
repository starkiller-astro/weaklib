MODULE wlExtTableFinishingModule

  USE wlKindModule, ONLY: dp 
  USE wlInterpolationModule

  implicit none

  CONTAINS

    SUBROUTINE LoneCellLocate( Fail, LoneCell )

      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: Fail
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
      DO k = 2, SIZE(Fail, DIM=3) - 1 
        DO j = 2, SIZE(Fail, DIM=2) - 1
          DO i = 2, SIZE(Fail, DIM=1) - 1

            IF ( .not.Fail(i,j,k) ) CYCLE 
            IF ( ANY( Fail(i-1:i+1,j-1:j+1,k-1:k+1) .and. CornerMask ) ) CYCLE
            LoneCell(i,j,k) = .true.

            WRITE(*,*) i,j,k 

          END DO
        END DO
      END DO

    END SUBROUTINE LoneCellLocate

    SUBROUTINE HoleCharacterizeFine( Fail, Table, iMinGradient )

      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: Fail
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: iMinGradient

      INTEGER :: i, j, k, idim 
      INTEGER :: isize, jsize, ksize 
      REAL(dp) :: Gradient
      REAL(dp) :: MinGradient
      LOGICAL, DIMENSION(3) :: LinearOK 

      iMinGradient = 0

      isize = SIZE(Fail, DIM=1)
      jsize = SIZE(Fail, DIM=2)
      ksize = SIZE(Fail, DIM=3)

      DO k = 1, ksize
        DO j = 1, jsize
          DO i = 1, isize 
            
            IF ( .not.Fail(i,j,k) ) CYCLE
      
            LinearOK = .false.
            MinGradient = 1.d99

            IF ( i == 1 .or. i == isize ) THEN 
            ELSE
              IF ( .not.Fail(i-1,j,k) .and. .not.Fail(i+1,j,k) ) LinearOK(1) = .true. 
            END IF
            IF ( LinearOK(1) ) THEN 
              Gradient = ABS(Table(i+1,j,k) - Table(i-1,j,k))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 1
              END IF  
            END IF

            IF ( j == 1 .or. j == jsize ) THEN 
            ELSE
              IF ( .not.Fail(i,j-1,k) .and. .not.Fail(i,j+1,k) ) LinearOK(2)= .true. 
            END IF
            IF ( LinearOK(2) ) THEN 
              Gradient = ABS(Table(i,j+1,k) - Table(i,j-1,k))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 2
              END IF  
            END IF

            IF ( k == 1 .or. k == ksize ) THEN 
            ELSE
              IF ( .not.Fail(i,j,k-1) .and. .not.Fail(i,j,k+1) ) LinearOK(3) = .true. 
            END IF
            IF ( LinearOK(3) ) THEN 
              Gradient = ABS(Table(i,j,k+1) - Table(i,j,k-1))
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 3
              END IF  
            END IF

          END DO
        END DO
      END DO

    END SUBROUTINE HoleCharacterizeFine

    SUBROUTINE HoleCharacterizeCoarse( Fail, Repaired, Table, iMinGradient, iLimits )

      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: Fail
      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: Repaired 
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: iMinGradient
      INTEGER, DIMENSION(:,:,:,:), INTENT(out) :: iLimits

      INTEGER :: i, j, k, Low, High
      INTEGER :: isize, jsize, ksize, Test
      REAL(dp) :: Gradient
      REAL(dp) :: MinGradient
      LOGICAL, DIMENSION(3) :: LinearOK

      iMinGradient = 0

      isize = SIZE(Fail, DIM=1)
      jsize = SIZE(Fail, DIM=2)
      ksize = SIZE(Fail, DIM=3)

      DO k = 1, ksize
        DO j = 1, jsize
          DO i = 1, isize

            IF ( .not.Fail(i,j,k) .or. Repaired(i,j,k) ) CYCLE

            LinearOK = .false.
            MinGradient = 1.d99

            IF ( i == 1 .or. i == isize ) THEN
            ELSE
              Low = 0
              High = 0
              DO Test = i-1, 1, -1 
                IF( .not.Fail(Test,j,k) ) THEN
                  Low = Test
                  EXIT
                END IF
              END DO              

              DO Test = i+1, isize
                IF( .not.Fail(Test,j,k) ) THEN
                  High = Test
                  EXIT
                END IF
              END DO              
      
              IF ( Low > 0 .and. High > 0 ) LinearOK(1) = .true.

            END IF
            IF ( LinearOK(1) ) THEN
              Gradient = ABS( Table(High,j,k) - Table(Low,j,k) )
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 1
                iLimits(1,i,j,k) = Low
                iLimits(2,i,j,k) = High
              END IF
            END IF

            IF ( j == 1 .or. j == jsize ) THEN
            ELSE
              Low = 0
              High = 0
              DO Test = j-1, 1, -1
                IF( .not.Fail(i,Test,k) ) THEN
                  Low = Test
                  EXIT
                END IF
              END DO

              DO Test = j+1, jsize
                IF( .not.Fail(i,Test,k) ) THEN
                  High = Test
                  EXIT
                END IF
              END DO

              IF ( Low > 0 .and. High > 0 ) LinearOK(2) = .true.

            END IF
            IF ( LinearOK(2) ) THEN
              Gradient = ABS( Table(i,High,k) - Table(i,Low,k) )
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 2
                iLimits(1,i,j,k) = Low
                iLimits(2,i,j,k) = High
              END IF
            END IF

            IF ( k == 1 .or. k == ksize ) THEN
            ELSE
              Low = 0
              High = 0
              DO Test = k-1, 1, -1
                IF( .not.Fail(i,j,Test) ) THEN
                  Low = Test
                  EXIT
                END IF
              END DO

              DO Test = k+1, ksize
                IF( .not.Fail(i,j,Test) ) THEN
                  High = Test
                  EXIT
                END IF
              END DO

              IF ( Low > 0 .and. High > 0 ) LinearOK(3) = .true.

            END IF
            IF ( LinearOK(3) ) THEN
              Gradient = ABS( Table(i,j,High) - Table(i,j,Low) )
              IF ( Gradient < MinGradient ) THEN
                MinGradient = Gradient
                iMinGradient(i,j,k) = 3
                iLimits(1,i,j,k) = Low
                iLimits(2,i,j,k) = High
              END IF
            END IF

          END DO
        END DO
      END DO

    END SUBROUTINE HoleCharacterizeCoarse

END MODULE wlExtTableFinishingModule
