PROGRAM wlRewriteEquationOfStateTest

  USE wlKindModule, ONLY: DP
  USE wlExtNumericalModule, ONLY: epsilon
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlEOSInterpolationModule, ONLY: &
      EOSTableQuery

  implicit none

  INTEGER                        :: i, i_rho, i_t, i_ye, ivar
  INTEGER, DIMENSION(3)          :: nPoints
  INTEGER                        :: nVariables
  INTEGER                        :: j
  TYPE(EquationOfStateTableType) :: EOSInTable, EOSOutTable
  INTEGER(HID_T)                 :: file_id
  INTEGER(HID_T)                 :: group_id
  REAL(DP)                       :: minvar, maxvar
  REAL(DP), ALLOCATABLE          :: Interpolants(:,:)
  REAL(DP), ALLOCATABLE          :: t(:), ye(:)
  
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF &
       ( EOSInTable, "wl-EOS-LS220-20-40-100-Lower-T.h5" )

  !! ----- USE FOR REWRITE OLD TABLE WITH ADDING MIN-MAX
!!$  nPoints = EOSInTable % DV % nPoints
!!$  nVariables = EOSInTable % DV % nVariables
!!$
!!$  CALL AllocateEquationOfStateTable( EOSOutTable, nPoints , nVariables )
!!$
!!$  EOSOutTable % TS = EOSInTable % TS
!!$  EOSOutTable % DV = EOSInTable % DV
!!$  EOSOutTable % MD = EOSInTable % MD
!!$
!!$  DO ivar = 1, EOSOutTable % DV % nVariables
!!$    minvar = MINVAL( EOSOutTable % DV % Variables(ivar) % Values )
!!$    EOSOutTable % DV % minValues(ivar) &
!!$    = 10.0d0**minvar - EOSOutTable % DV % Offsets(ivar)
!!$
!!$    maxvar = MAXVAL( EOSOutTable % DV % Variables(ivar) % Values )
!!$    EOSOutTable % DV % maxValues(ivar) &
!!$    = 10.0d0**maxvar - EOSOutTable % DV % Offsets(ivar)
!!$  END DO
!!$
!!$  CALL InitializeHDF( )
!!$  CALL WriteEquationOfStateTableHDF( EOSOutTable, &
!!$       "wl-EOS-LS220-20-40-100-Lower-T-rewrite.h5" )
!!$  CALL DeAllocateEquationOfStateTable( EOSOutTable )
  
  !! ------ USE FOR COARSING EOS TABLE
  nPoints = (/121,68,25/)
  nVariables = EOSInTable % DV % nVariables

  ALLOCATE( t(nPoints(1)) )
  ALLOCATE( ye(nPoints(1)) )
  ALLOCATE( Interpolants(nPoints(1),nVariables) )

  CALL AllocateEquationOfStateTable( EOSOutTable, nPoints , nVariables )

  EOSOutTable % TS % Names = EOSInTable % TS % Names
  EOSOutTable % TS % Units = EOSInTable % TS % Units
  EOSOutTable % TS % Indices = EOSInTable % TS % Indices
  EOSOutTable % TS % LogInterp = EOSInTable % TS % LogInterp
  EOSOutTable % TS % minValues = EOSInTable % TS % minValues
  EOSOutTable % TS % maxValues = EOSInTable % TS % maxValues

  CALL MakeLogGrid( EOSOutTable % TS % minValues(1), EOSOutTable % TS % maxValues(1),&
         EOSOutTable % TS % nPoints(1), EOSOutTable % TS % States(1) % Values)
  CALL MakeLogGrid( EOSOutTable % TS % minValues(2), EOSOutTable % TS % maxValues(2),&
         EOSOutTable % TS % nPoints(2), EOSOutTable % TS % States(2) % Values)
  CALL MakeLinearGrid( EOSOutTable % TS % minValues(3), EOSOutTable % TS % maxValues(3),&
         EOSOutTable % TS % nPoints(3), EOSOutTable % TS % States(3) % Values)

  EOSOutTable % DV % nVariables = EOSInTable % DV % nVariables
  EOSOutTable % DV % nPoints = nPoints
  EOSOutTable % DV % Names = EOSInTable % DV % Names
  EOSOutTable % DV % Units = EOSInTable % DV % Units
  EOSOutTable % DV % Indices = EOSInTable % DV % Indices

  ASSOCIATE &
  ( iRho => EOSOutTable % TS % Indices % iRho, &
    iT   => EOSOutTable % TS % Indices % iT, &
    iYe  => EOSOutTable % TS % Indices % iYe )
 
    DO i_ye = 1,nPoints(iYe)
      ye = EOSOutTable % TS % States(iYe) % Values(i_ye)
      DO i_t = 1,nPoints(iT)
        t = EOSOutTable % TS % States(iT) % Values(i_t) 
        DO i_rho = 1,nPoints(iRho)
  
          CALL EOSTableQuery &
               ( EOSOutTable % TS % States(iRho) % Values, &
                 t, ye, &
                 EOSOutTable % TS % LogInterp, &
                 EOSInTable  % TS, EOSInTable % DV, Interpolants )
  
          DO ivar = 1,EOSOutTable % DV % nVariables
             EOSOutTable % DV % Variables(ivar) % Values(:,i_t,i_ye) &
             = Interpolants(:,ivar)
          END DO
  
        END DO
      END DO
    END DO

  END ASSOCIATE

  DO ivar = 1, EOSOutTable % DV % nVariables
    minvar = MINVAL( EOSOutTable % DV % Variables(ivar) % Values )
    EOSOutTable % DV % minValues(ivar) = minvar

    maxvar = MAXVAL( EOSOutTable % DV % Variables(ivar) % Values )
    EOSOutTable % DV % maxValues(ivar) = maxvar

    EOSOutTable % DV % Offsets(ivar) = -2.d0 * MIN( 0.d0, minvar )

    EOSOutTable % DV % Variables(ivar) % Values &
      = LOG10( EOSOutTable % DV % Variables(ivar) % Values &
               + EOSOutTable % DV % Offsets(ivar) + epsilon )
  END DO

  CALL InitializeHDF( )
  CALL WriteEquationOfStateTableHDF &
       ( EOSOutTable, "wl-EOS-LS220-15-25-50-Lower-T-rewrite.h5" )
  CALL FinalizeHDF( )

  DEALLOCATE( t, ye, Interpolants )
  CALL DeAllocateEquationOfStateTable( EOSInTable )
  CALL DeAllocateEquationOfStateTable( EOSOutTable )
  
  
  WRITE (*,*) "HDF write successful"

END PROGRAM wlRewriteEquationOfStateTest
