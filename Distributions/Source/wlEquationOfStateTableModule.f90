MODULE wlEquationOfStateTableModule

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none
  PRIVATE

  TYPE, PUBLIC :: EquationOfStateTableType
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Name
    INTEGER :: nVariables
    INTEGER, DIMENSION(3) :: nPoints
    TYPE(ThermoStateType), ALLOCATABLE           :: TS 
    TYPE(DependentVariablesType), ALLOCATABLE    :: DV
  END TYPE


  PUBLIC AllocateEquationOfStateTable
  PUBLIC DeAllocateEquationOfStateTable

  SUBROUTINE AllocateEquationOfStateTable( TS, DV, nPoints, nVariables )

    TYPE(EquationOfStateTableType) :: EOSTable
    INTEGER, INTENT(in) :: nVariables
    INTEGER, DIMENSION(3), INTENT(in) :: nPoints

  END SUBROUTINE AllocateEquationOfStateTable

END MODULE wlEquationOfStateTableModule
