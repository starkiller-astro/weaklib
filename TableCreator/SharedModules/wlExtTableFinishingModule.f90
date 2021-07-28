MODULE wlExtTableFinishingModule

  USE wlKindModule, ONLY: dp 
  USE wlInterpolationModule

  implicit none

  PUBLIC LogInterpolateFine1D
  PUBLIC LogInterpolateCoarse1D
  PUBLIC LoneCellLogInterpolateSingleVariable

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

  SUBROUTINE LogInterpolateFine1D( i, j, k, iMinGradient, delta, Table, Interpolant )

    INTEGER, INTENT(in) :: i, j, k
    REAL(dp), INTENT(in) :: delta
    INTEGER, INTENT(in) :: iMinGradient
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table

    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p0, p1, epsilon, offset

    epsilon = 1.d-200

    SELECT CASE( iMinGradient )

      CASE(1)
        offset = -2.d0*MIN( 0.d0, Table( i-1, j, k ), Table( i+1, j, k ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i-1, j, k ) + offset + epsilon  )
        p1 = LOG10( Table( i+1, j, k ) + offset + epsilon )

      CASE(2)
        offset = -2.d0*MIN( 0.d0, Table( i, j-1, k ), Table( i, j+1, k ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j-1, k ) + offset + epsilon  )
        p1 = LOG10( Table( i, j+1, k ) + offset + epsilon  )

      CASE(3)
        offset = -2.d0*MIN( 0.d0, Table( i, j, k-1 ), Table( i, j, k+1 ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j, k-1 ) + offset + epsilon  )
        p1 = LOG10( Table( i, j, k+1 ) + offset + epsilon  )

    END SELECT

    Interpolant = 10.d0**( delta * p1 + ( 1.d0 - delta ) * p0 ) - offset
    WRITE (*,*) "Interpolant=", Interpolant
    WRITE (*,*) "Delta=", delta

  END SUBROUTINE LogInterpolateFine1D

  SUBROUTINE LogInterpolateCoarse1D( i, j, k, iMinGradient, iLimits, delta, Table, Interpolant )

    INTEGER, INTENT(in) :: i, j, k
    REAL(dp), INTENT(in) :: delta
    INTEGER, INTENT(in) :: iMinGradient
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    INTEGER, DIMENSION(:,:,:,:), INTENT(in) :: iLimits

    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p0, p1, offset, epsilon

    epsilon = 1.d-200

    SELECT CASE( iMinGradient )

      CASE(1)
        offset = -2.d0*MIN( 0.d0, Table(iLimits(1,i,j,k ),j,k), Table(iLimits(2,i,j,k),j,k) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( iLimits(1,i,j,k), j, k ) + offset + epsilon  )
        p1 = LOG10( Table( iLimits(2,i,j,k), j, k ) + offset + epsilon  )

      CASE(2)
        offset = -2.d0*MIN( 0.d0, Table(i,iLimits(1,i,j,k ),k), Table(i,iLimits(2,i,j,k),k) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, iLimits(1,i,j,k), k ) + offset + epsilon  )
        p1 = LOG10( Table( i, iLimits(2,i,j,k), k ) + offset + epsilon  )

      CASE(3)
        offset = -2.d0*MIN( 0.d0, Table(i,j,iLimits(1,i,j,k )), Table(i,j,iLimits(2,i,j,k)) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j, iLimits(1,i,j,k) ) + offset + epsilon  )
        p1 = LOG10( Table( i, j, iLimits(2,i,j,k) ) + offset + epsilon  )

    END SELECT
    WRITE (*,*) "p0, p1 =", p0, p1
    Interpolant = 10.d0**( delta * p1 + ( 1.d0 - delta ) * p0 ) - offset
    WRITE (*,*) "Interpolant=", Interpolant
    WRITE (*,*) "Delta=", delta

  END SUBROUTINE LogInterpolateCoarse1D

  SUBROUTINE LoneCellLogInterpolateSingleVariable( x1, x2, x3, Coordinate1, Coordinate2, &
                                           Coordinate3, LogInterp, Table, Interpolant )

    INTEGER, INTENT(in) :: x1
    INTEGER, INTENT(in) :: x2
    INTEGER, INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table


    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: il1, il2, il3
    REAL(dp) :: PreInterpolant


      il1 = x1 - 1
      il2 = x2 - 1
      il3 = x3 - 1

      p000 = LOG10( Table( il1  , il2  , il3   ) )
      p100 = LOG10( Table( il1+2, il2  , il3   ) )
      p010 = LOG10( Table( il1  , il2+2, il3   ) )
      p110 = LOG10( Table( il1+2, il2+2, il3   ) )
      p001 = LOG10( Table( il1  , il2  , il3+2 ) )
      p101 = LOG10( Table( il1+2, il2  , il3+2 ) )
      p011 = LOG10( Table( il1  , il2+2, il3+2 ) )
      p111 = LOG10( Table( il1+2, il2+2, il3+2 ) )

      WRITE (*,*) p000

      IF ( LogInterp(1) ) THEN
      delta(1) = LOG10( Coordinate1(x1) / Coordinate1(il1) ) &
                   / LOG10( Coordinate1(il1+2) / Coordinate1(il1) )
      ELSE
      delta(1) = ( Coordinate1(x1) - Coordinate1(il1) ) / &
                   ( Coordinate1(il1+2) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) ) THEN
      delta(2) = LOG10( Coordinate2(x2) / Coordinate2(il2) ) &
                   / LOG10( Coordinate2(il2+2) / Coordinate2(il2) )
      ELSE
      delta(2) = ( Coordinate2(x2) - Coordinate2(il2) ) &
                   / ( Coordinate2(il2+2) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) ) THEN
      delta(3) = LOG10( Coordinate3(x3) / Coordinate3(il3) ) &
                   / LOG10( Coordinate3(il3+2) / Coordinate3(il3) )
      ELSE
      delta(3) = ( Coordinate3(x3) - Coordinate3(il3) ) &
                   / ( Coordinate3(il3+2) - Coordinate3(il3) )
      END IF

      WRITE (*,*) delta
      Interpolant &
        = 10.d0**( &
           (  (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &     
                                   +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                   + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                   +            delta(1)  *           delta(2)  * p110 ) &
                      + delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                   +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                   +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                   +            delta(1)  *           delta(2)  * p111 ) &
                 ) )
  WRITE (*,*) Interpolant

  END SUBROUTINE LoneCellLogInterpolateSingleVariable

END MODULE wlExtTableFinishingModule
