PROGRAM wlInversionTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF,  & 
                           ReadEquationOfStateTableHDF, & 
                           WriteEquationOfStateTableHDF
  implicit none

  INTEGER  :: i, j, k
  REAL(dp) :: maxnorm, norm
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table, Temp_Table, E_Internal, &
                                         Temperature, Test_Temps
  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  ALLOCATE( Energy_Table( nPoints(2) ), Temp_Table( nPoints(2) ),          &
            E_Internal( nPoints(2) - 1 ), Temperature( nPoints(2) - 1 ),   &
            Test_Temps( nPoints(2) - 1 ) )
  
  END ASSOCIATE

  ! Internal energy for (10,:,10) in the standard res table runs from
  ! 8.22e32 to 1.1e42 
  !    nPoints = (/161,55,47/)


!    WRITE (*,*) "E_Table =", Energy_Table

  Temp_Table = ( EOSTable % TS % States(2) % Values ) 

  DO i = 1, SIZE(Temp_Table) - 1
    !E_Internal(i) = Energy_Table(i) * 1.05d0  
    Test_Temps(i) = ( Temp_Table(i+1) - Temp_Table(i) )/2 + Temp_Table(i)
  END DO

  maxnorm = 0.0d0

  DO j = 1, 47
    DO i = 1, 161
    CALL LogInterpolateSingleInputSingleVariable &
           ( i,           &
             Test_Temps,                                   &
             j,           &
             EOSTable % TS % States(1) % Values,           &
             EOSTable % TS % States(2) % Values,           &
             EOSTable % TS % States(3) % Values,           &
             LogInterp,                                    &
             EOSTable % DV % Offsets(3),                   &
             EOSTable % DV % Variables(3) % Values(:,:,:), &
             E_Internal )

      Energy_Table = 10**( EOSTable % DV % Variables(3) % Values(i,:,j) ) &
                     - EOSTable % DV % Offsets(3) 

      CALL ComputeTempFromIntEnergy( E_Internal, Energy_Table, Temp_Table, &
                          EOSTable % DV % Offsets(3), Temperature )  

      DO k = 1, SIZE(Temperature)
        norm = ABS( Temperature(k) - Test_Temps(k) ) / Test_Temps(k) 
        WRITE (*,*) norm, i, j
        IF ( norm > maxnorm) THEN 
          maxnorm = norm
        END IF
      END DO  
    END DO
  END DO
  
  WRITE (*,*) "Maxnorm=", maxnorm

!  DO i = 1, SIZE(E_Internal)

!    WRITE (*,'(i4, 4es12.5)') i, E_Internal(i), Temperature(i), Test_Temps(i),    &
!                ( Temperature(i) - Test_Temps(i) ) / Test_Temps(i)
!  END DO

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInversionTest
