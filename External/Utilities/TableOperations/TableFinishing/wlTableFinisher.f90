PROGRAM wlTableFinisher

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlExtNumericalModule
  USE wlExtPhysicalConstantsModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
                  
                 
  implicit none

  INTEGER  :: i, j, k, l, count
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iMinGradient
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iLimits
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: InterpolantFine, InterpolantCoarse
  REAL(dp) :: DeltaFine, DeltaCoarse
  REAL(dp) :: minvar, x, y
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Fail 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Repaired 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: LoneCells    

  REAL(dp)                :: chem_n           ! neutron chemical potential [MeV]
  REAL(dp)                :: chem_n_save      ! neutron chemical potential [MeV]
  REAL(dp)                :: chem_p           ! proton chemical potential [MeV]
  REAL(dp)                :: chem_p_save     ! proton chemical potential [MeV]
  REAL(dp)                :: xhe              ! helium mass fraction
  REAL(dp)                :: x_neutron        ! neutron mass fraction
  REAL(dp)                :: XNUT             ! neutron mass fraction
  REAL(dp)                :: x_proton         ! proton mass fraction
  REAL(dp)                :: XPROT            ! proton mass fraction
  REAL(dp)                :: XA               ! Z/A
  REAL(dp)                :: XH               ! heavy mass fraction
  REAL(dp), PARAMETER     :: c0 = ( mb/( 2.d0 * pi * hbarc**2 ) )**1.5d0
  REAL(dp)                :: therm            ! thermal parameter in the chemical potential
  REAL(dp)                :: brydns           ! LS input baryon density
  REAL(dp)                :: tmev             ! temperature [MeV]
  REAL(dp)                :: ye               ! electron fraction
  REAL(dp), PARAMETER     :: x_min = 1.d-30   ! minimum mass fraction fraction

  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  WRITE (*,*) "Table Read", nPoints 

  ALLOCATE( Fail( nPoints(1), nPoints(2), nPoints(3) ),         &
            Repaired( nPoints(1), nPoints(2), nPoints(3) ),     & 
            LoneCells( nPoints(1), nPoints(2), nPoints(3) ),    & 
            iMinGradient( nPoints(1), nPoints(2), nPoints(3) ), &  
            iLimits( 2, nPoints(1), nPoints(2), nPoints(3) ) ) 
  
  END ASSOCIATE

  WRITE (*,*) "Logicals Allocated" 

  Repaired = .false.

  Fail(:,:,:) = EOSTable % DV % Repaired(:,:,:) < 0.0d0  
!  Fail(:,:,:) = EOSTable % DV % Variables(1) % Values(:,:,:) <= 0.0d0  &
!                .or. EOSTable % DV % Variables(2) % Values(:,:,:) <= 0.0d0 &
!                .or. EOSTable % DV % Variables(3) % Values(:,:,:) <= 0.0d0 &
!                .or. EOSTable % DV % Variables(7) % Values(:,:,:) < 0.0d0 &
!                .or. EOSTable % DV % Variables(8) % Values(:,:,:) < 0.0d0 &
!                .or. EOSTable % DV % Variables(9) % Values(:,:,:) < 0.0d0 &
!                .or. EOSTable % DV % Variables(15) % Values(:,:,:) <= 0.0d0 &
!                .or. ISNAN(EOSTable % DV % Variables(1) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(2) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(3) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(7) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(8) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(9) % Values(:,:,:))  & 
!                .or. ISNAN(EOSTable % DV % Variables(15) % Values(:,:,:))   

  !------------------------------------------------------
  ! Find dimension (iMinGradient) with smallest gradient 
  !   for interpolation across single cell hole 
  !------------------------------------------------------
  
  CALL HoleCharacterizeFine( Fail, EOSTable % DV % Variables(1) % Values, iMinGradient )

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

        IF ( .not.Fail(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE 

        DeltaFine = 0.5d0 

        DO l = 1, EOSTable % nVariables

          CALL LogInterpolateFine1D&
                 ( i, j, k, iMinGradient(i,j,k), DeltaFine, &
                 EOSTable % DV % Variables(l) % Values, InterpolantFine )

          EOSTable % DV % Variables(l) % Values(i,j,k) = InterpolantFine

          Repaired(i,j,k) = .true.

        END DO

        WRITE (*,*) InterpolantFine 

      END DO
    END DO
  END DO

  CALL HoleCharacterizeCoarse( Fail, Repaired, &
         EOSTable % DV % Variables(1) % Values, iMinGradient, iLimits ) 

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

        IF ( .not.Fail(i,j,k) .or. Repaired(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE

        WRITE (*,*) "iLimits =", iLimits(1,i,j,k), iLimits(2,i,j,k)

        SELECT CASE( iMinGradient(i,j,k) )
          CASE(1)
            DeltaCoarse = DBLE( i - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
          CASE(2)
            DeltaCoarse = DBLE( j - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
          CASE(3)
            DeltaCoarse = DBLE( k - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
        END SELECT

        DO l = 1, EOSTable % nVariables
          CALL LogInterpolateCoarse1D&
                 ( i, j, k, iMinGradient(i,j,k), iLimits, DeltaCoarse, &
                  EOSTable % DV % Variables(l) % Values, InterpolantCoarse )

          EOSTable % DV % Variables(l) % Values(i,j,k) = InterpolantCoarse
          Repaired(i,j,k) = .true.

        WRITE (*,*) InterpolantCoarse
        END DO

      END DO
    END DO
  END DO

  count = 0

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

      IF ( Repaired(i,j,k) ) THEN
        EOSTable % DV % Repaired(i,j,k) = 1
        count = count + 1
      ELSE 
        EOSTable % DV % Repaired(i,j,k) = 0
      END IF      

      END DO 
    END DO
  END DO

  !CALL MonotonicityCheck( EOSTable % DV % Variables(3) % Values(:,:,:), &
  !                        EOSTable % nPoints(1), EOSTable % nPoints(2), &
  !                        EOSTable % nPoints(3), 2, EOSTable % DV % Repaired )

!--------------------------------------------------------------------
!
!  This next bit of code is only for use with LS in CHIMERA.
!  Left to do: confirm X = Z/A
!              rewrite data
!              make sure physical constant values match/are present
!              confirm and convert internal XPROT, XNUT, etc. stuff
!
!
!
!--------------------------------------------------------------------

  DO k = 1, EOSTable % TS % nPoints(3)
    DO j = 1, EOSTable % TS % nPoints(2)
      DO i = 1, EOSTable % TS % nPoints(1)
          IF ( EOSTable % TS % States(3) % Values(k) > 0.50d0 ) THEN


        brydns            = kfm * EOSTable % TS % States(1) % Values(i)
        tmev              = kmev * EOSTable % TS % States(2) % Values(j)
        ye                = EOSTable % TS % States(3) % Values(k)

        x_neutron         = EOSTable % DV % Variables(8) % Values(i,j,k)
        XNUT              = x_neutron
        x_proton          = EOSTable % DV % Variables(7) % Values(i,j,k)
        XPROT             = x_proton

        chem_n            = EOSTable % DV % Variables(6) % Values(i,j,k)
        chem_n_save       = chem_n 
        chem_p            = EOSTable % DV % Variables(5) % Values(i,j,k)
        chem_p_save       = chem_p
        XH                = EOSTable % DV % Variables(10) % Values(i,j,k)
        XA                = ( EOSTable % DV % Variables(11) % Values(i,j,k) ) &
                            / EOSTable % DV % Variables(12) % Values(i,j,k)

          xhe           = DMAX1( one - x_neutron - x_proton - XH, zero )
          IF ( x_neutron + x_proton + 0.5d0 * xhe + XH * XA > ye ) THEN
            x_proton    = DMAX1( ye - 0.5d0 * xhe - XH * XA, x_min )
            x_neutron   = DMAX1( XNUT + XPROT - x_proton, x_min )
            therm       = brydns/( c0 * tmev * sqrt(tmev) )
            chem_n      = tmev * DLOG( half * therm * x_neutron ) - dmnp
            chem_p      = tmev * DLOG( half * therm * x_proton ) - dmnp

            EOSTable % DV % Variables(8) % Values(i,j,k) = x_neutron
            EOSTable % DV % Variables(7) % Values(i,j,k) = x_proton
            EOSTable % DV % Variables(6) % Values(i,j,k) = chem_n
            EOSTable % DV % Variables(5) % Values(i,j,k) = chem_p

            IF ( ( XNUT .ne. x_neutron ) .or. ( XPROT .ne. x_proton ) .or. &
                 ( chem_p .ne. chem_p_save ) .or. ( chem_n .ne. chem_n_save ) ) THEN
              WRITE (*,*) 'abundances and potentials changed'

            END IF 

          END IF ! x_neutron + x_proton + 0.5d0 * xhe + XH * X > yep
        END IF ! ye > 0.51d0
      END DO
    END DO
  END DO

  DO l = 1, EOSTable % nVariables
    WRITE (*,*) EOSTable % DV % Names(l)
    minvar = MINVAL( EOSTable % DV % Variables(l) % Values )
    WRITE (*,*) "minvar=", minvar
    EOSTable % DV % Offsets(l) = -2.d0 * MIN( 0.d0, minvar )
    WRITE (*,*) "Offset=", EOSTable % DV % Offsets(l)
    EOSTable % DV % Variables(l) % Values &
      = LOG10( EOSTable % DV % Variables(l) % Values &
               + EOSTable % DV % Offsets(l) + epsilon )        

  END DO

WRITE (*,*) count, " repairs " 

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinisher
