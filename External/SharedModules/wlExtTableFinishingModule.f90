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

    SUBROUTINE HoleCharacterize( fail, LinOkX, LinOkY, LinOkZ )
      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: fail
      LOGICAL, DIMENSION(:,:,:), INTENT(out) :: LinOkX 
      LOGICAL, DIMENSION(:,:,:), INTENT(out) :: LinOkY 
      LOGICAL, DIMENSION(:,:,:), INTENT(out) :: LinOkZ 

      INTEGER :: i, j, k

      LinOkX = .false.
      DO k = 1, SIZE(fail, DIM=3)
        DO j = 1, SIZE(fail, DIM=2)
          DO i = 1, SIZE(fail, DIM=1)

            IF ( .not.fail(i,j,k) ) CYCLE
            IF ( .not.fail(i-1,j,k) .and. .not.fail(i+1,j,k) .and. i/=1 .and. &
                 i/= SIZE(fail, DIM=1)  ) THEN 
            LinOkX(i,j,k) = .true.
            WRITE (*,*) "Ok to linearly interpolate in x at point =", i, j, k 
            END IF 

          END DO
        END DO
      END DO

      LinOkY = .false.
      DO k = 1, SIZE(fail, DIM=3)
        DO j = 1, SIZE(fail, DIM=2)
          DO i = 1, SIZE(fail, DIM=1)

            IF ( .not.fail(i,j,k) ) CYCLE
            IF ( .not.fail(i,j-1,k) .and. .not.fail(i,j+1,k) .and. j/=1 .and. & 
                 j/= SIZE(fail, DIM=2)  ) THEN 
            LinOkY(i,j,k) = .true.
            WRITE (*,*) "Ok to linearly interpolate in y at point =", i, j, k 
            END IF 

          END DO
        END DO
      END DO

      LinOkZ = .false.
      DO k = 1, SIZE(fail, DIM=3)
        DO j = 1, SIZE(fail, DIM=2)
          DO i = 1, SIZE(fail, DIM=1)

            IF ( .not.fail(i,j,k) ) CYCLE
            IF ( .not.fail(i,j,k-1) .and. .not.fail(i,j,k+1) .and. k/=1 .and. & 
                 k/= SIZE(fail, DIM=3)  ) THEN 
            LinOkZ(i,j,k) = .true.
            WRITE (*,*) "Ok to linearly interpolate in z at point =", i, j, k 
            END IF

          END DO
        END DO
      END DO

    END SUBROUTINE HoleCharacterize

    SUBROUTINE GradientCheck( i, j, k, LinOkX, LinOkY, LinOkZ, Table, MinGradient )

      INTEGER, INTENT(in) :: i 
      INTEGER, INTENT(in) :: j
      INTEGER, INTENT(in) :: k
      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: LinOkX
      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: LinOkY
      LOGICAL, DIMENSION(:,:,:), INTENT(in) :: LinOkZ
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table

      INTEGER, INTENT(out) :: MinGradient 
      REAL(dp) :: px0, px1, py0, py1, pz0, pz1 
      INTEGER :: il1, il2, il3

      MinGradient = 0

      IF ( LinOkX(i,j,k) .and. LinOkY(i,j,k) .and. LinOkZ(i,j,k) ) THEN

           px0 = Table( i-1, j, k ) 
           px1 = Table( i+1, j, k ) 
           py0 = Table( i, j-1, k ) 
           py1 = Table( i, j+1, k ) 
           pz0 = Table( i, j, k-1 ) 
           pz1 = Table( i, j, k+1 )

        IF ( ( LOG10(px1/px0) <= LOG10(py1/py0) ) .and. ( LOG10(px1/px0) <= LOG10(pz1/pz0) ) ) THEN  

        MinGradient = 1
             
        ELSE IF ( ( LOG10(py1/py0) <= LOG10(px1/px0) ) .and. ( LOG10(py1/py0) <= LOG10(pz1/pz0) ) ) THEN 
        MinGradient = 2

        ELSE IF ( ( LOG10(pz1/pz0) <= LOG10(py1/py0) ) .and. ( LOG10(pz1/pz0) <= LOG10(px1/px0) ) ) THEN   
        MinGradient = 3

        END IF

      ELSE IF ( LinOkX(i,j,k) .and. LinOkY(i,j,k) .and. .not.LinOkZ(i,j,k) ) THEN

           px0 = Table( i-1, j, k )
           px1 = Table( i+1, j, k )
           py0 = Table( i, j-1, k )
           py1 = Table( i, j+1, k )

        IF ( LOG10(px1/px0) <= LOG10(py1/py0) ) THEN    

        MinGradient = 1

        ELSE IF ( LOG10(py1/py0) <= LOG10(px1/px0) ) THEN  

        MinGradient = 2

        END IF

      ELSE IF ( LinOkX(i,j,k) .and. .not.LinOkY(i,j,k) .and. LinOkZ(i,j,k) ) THEN

           px0 = Table( i-1, j, k )
           px1 = Table( i+1, j, k )
           pz0 = Table( i, j, k-1 )
           pz1 = Table( i, j, k+1 )

        IF ( LOG10(px1/px0) <= LOG10(py1/py0) .and. LOG10(px1/px0) <= LOG10(pz1/pz0) ) THEN 

        MinGradient = 1

        ELSE IF ( LOG10(pz1/pz0) <= LOG10(py1/py0) .and. LOG10(pz1/pz0) <= LOG10(px1/px0) ) THEN   

        MinGradient = 3

        END IF

      ELSE IF ( .not.LinOkX(i,j,k) .and. LinOkY(i,j,k) .and. LinOkZ(i,j,k) ) THEN

           py0 = Table( i, j-1, k )
           py1 = Table( i, j+1, k )
           pz0 = Table( i, j, k-1 )
           pz1 = Table( i, j, k+1 )
             
        IF ( LOG10(py1/py0) <= LOG10(px1/px0) .and. LOG10(py1/py0) <= LOG10(pz1/pz0) ) THEN   

        MinGradient = 2

        ELSE IF ( LOG10(pz1/pz0) <= LOG10(py1/py0) .and. LOG10(pz1/pz0) <= LOG10(px1/px0) ) THEN   

        MinGradient = 3

        END IF

      ELSE IF ( LinOkX(i,j,k) .and. .not.LinOkY(i,j,k) .and. .not.LinOkZ(i,j,k) ) THEN

        MinGradient = 1

      ELSE IF ( .not.LinOkX(i,j,k) .and. LinOkY(i,j,k) .and. .not.LinOkZ(i,j,k) ) THEN

        MinGradient = 2

      ELSE IF ( .not.LinOkX(i,j,k) .and. .not.LinOkY(i,j,k) .and. LinOkZ(i,j,k) ) THEN

        MinGradient = 3

      END IF

    END SUBROUTINE GradientCheck 

END MODULE wlExtTableFinishingModule
