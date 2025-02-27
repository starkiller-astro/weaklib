MODULE wlOpacityFieldsModule

  USE wlKindModule, ONLY: dp
  USE wlParallelModule, ONLY: &
    myid, ierr, mpi_is_init

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: iNu_e     = 1
  INTEGER, PUBLIC, PARAMETER :: iNu_e_bar = 2
  INTEGER, PUBLIC, PARAMETER :: iNu_x     = 3
  INTEGER, PUBLIC, PARAMETER :: iNu_x_bar = 4
  INTEGER, PUBLIC, PARAMETER :: nSpecies  = 4

  ! --- NNS Scattering Kernels ( indices of Phi ( :, : ) => Kernel ( : ) ) ---

  INTEGER, PUBLIC, PARAMETER :: iNu_NNS      = 1
  INTEGER, PUBLIC, PARAMETER :: iNuBar_NNS   = 2
  INTEGER, PUBLIC, PARAMETER :: iNeutron_NNS = 1
  INTEGER, PUBLIC, PARAMETER :: iProton_NNS  = 2  

  ! --- NES Scattering Kernels ---

  INTEGER, PUBLIC, PARAMETER :: iHi0  = 1
  INTEGER, PUBLIC, PARAMETER :: iHii0 = 2
  INTEGER, PUBLIC, PARAMETER :: iHi1  = 3
  INTEGER, PUBLIC, PARAMETER :: iHii1 = 4

  ! --- Pair Kernels ---

  INTEGER, PUBLIC, PARAMETER :: iJi0  = 1
  INTEGER, PUBLIC, PARAMETER :: iJii0 = 2
  INTEGER, PUBLIC, PARAMETER :: iJi1  = 3
  INTEGER, PUBLIC, PARAMETER :: iJii1 = 4

  TYPE :: ValueType_3D
    REAL(dp), ALLOCATABLE :: Values(:,:,:)
  END type ValueType_3D

  TYPE :: ValueType_4D
    REAL(dp), ALLOCATABLE :: Values(:,:,:,:)
  END type ValueType_4D

  TYPE :: ValueType_5D
    REAL(dp), ALLOCATABLE :: Values(:,:,:,:,:)
  END type ValueType_5D

!---------------------------------------------
!
! OpacityTypeEmAb (Emission / Absorption)
!   Dependency ( E, rho, T, Ye )
!     E:   Neutrino Energy
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!
!---------------------------------------------

  TYPE, PUBLIC :: OpacityTypeEmAb
    INTEGER                         :: nOpacities
    INTEGER                         :: nPoints(4)
    CHARACTER(LEN=32),  ALLOCATABLE :: Names(:)
    CHARACTER(LEN=32),  ALLOCATABLE :: Units(:)
    REAL(dp),           ALLOCATABLE :: Offsets(:)
    TYPE(ValueType_4D), ALLOCATABLE :: Opacity(:)

    INTEGER                         :: EC_table_nOpacities
    CHARACTER(LEN=32),  ALLOCATABLE :: EC_table_Names(:)
    CHARACTER(LEN=32),  ALLOCATABLE :: EC_table_Units(:)
    REAL(dp),           ALLOCATABLE :: EC_table_spec_Offsets(:)
    REAL(dp),           ALLOCATABLE :: EC_table_rate_Offsets(:)
    TYPE(ValueType_4D), ALLOCATABLE :: EC_table_spec(:)
    TYPE(ValueType_3D), ALLOCATABLE :: EC_table_rate(:)
    REAL(dp),           ALLOCATABLE :: EC_table_rho(:)
    REAL(dp),           ALLOCATABLE :: EC_table_T(:)
    REAL(dp),           ALLOCATABLE :: EC_table_Ye(:)
    REAL(dp),           ALLOCATABLE :: EC_table_E(:)

    REAL(dp)                        :: EC_table_rho_min
    REAL(dp)                        :: EC_table_rho_max
    REAL(dp)                        :: EC_table_T_min
    REAL(dp)                        :: EC_table_T_max
    REAL(dp)                        :: EC_table_Ye_min
    REAL(dp)                        :: EC_table_Ye_max
    REAL(dp)                        :: EC_table_E_min
    REAL(dp)                        :: EC_table_E_max

    INTEGER                         :: EC_table_nRho
    INTEGER                         :: EC_table_nT
    INTEGER                         :: EC_table_nYe
    INTEGER                         :: EC_table_nE

    INTEGER                         :: np_FK    
                                    !Fischer et al 2020 full kinematics rates for
                                    !EmAb on free nucleons

    INTEGER                         :: np_FK_inv_n_decay 
                                    !Fischer et al 2020 full kinematics rates for
                                    !inverse neutron decay

    INTEGER                         :: np_isoenergetic       
                                    !EmAb on free nucleons using isoenergetic approximation
                                    !Bruenn 1985
                                    !Mezzacappa & Bruenn (1993)

    INTEGER                         :: np_non_isoenergetic             
                                    !EmAb on free nucleons taking into account recoil,
                                    !nucleon final-state blocking, and special relativity
                                    !Reddy et al 1998

    INTEGER                         :: np_weak_magnetism     
                                    !Weak magnetism corrections for EmAb on free nucleons
                                    !Horowitz 1997

    INTEGER                         :: nuclei_EC_FFN                  
                                    !EmAb on nuclei using the Fuller, Fowler,
                                    !Neuman 1982, Ap. J. 252, 715 approximation for the
                                    !heavy nucleus matrix element as given in Bruenn 1985

    INTEGER                         :: nuclei_EC_table
                                    !EmAb on nuclei using a NSE-folded LMSH table
                                    !Langanke et al. (2003), Hix et al. (2003)
  END TYPE OpacityTypeEmAb

!---------------------------------------------
!
! OpacityTypeScat (Elastic Scattering)
!   Dependency ( E, l, rho, T, Ye )
!     E:   Neutrino Energy
!     l:   Legendre Moment
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!
! OpacityTypeScat (Inelastic Neutrino-Nucleon Scattering)
!   Dependency ( E', E, l, T, MuB )
!     E':  Neutrino Energy
!     E:   Neutrino Energy
!     l:   Legendre Moment
!     T:   Temperature
!     MuB: Baryon Chemical Potential
!
! OpacityTypeScat (Inelastic Neutrino-Electron Scattering)
!   Dependency ( E', E, l, T, Eta )
!     E':  Neutrino Energy
!     E:   Neutrino Energy
!     l:   Legendre Moment
!     T:   Temperature
!     Eta: Electron Chemical Pot. / Temperature
!
!---------------------------------------------

  TYPE, PUBLIC :: OpacityTypeScat
    INTEGER                         :: nOpacities
    INTEGER                         :: nMoments
    INTEGER                         :: nPoints(5)
    !REAL(dp)                        :: rho_min
    !REAL(dp)                        :: rho_max
    CHARACTER(LEN=32),  ALLOCATABLE :: Names(:)
    CHARACTER(LEN=32),  ALLOCATABLE :: Units(:)
    REAL(dp),           ALLOCATABLE :: Offsets(:,:)
    TYPE(ValueType_5D), POINTER     :: Kernel(:) => null ( )
  END TYPE OpacityTypeScat

  TYPE, PUBLIC, EXTENDS(OpacityTypeScat) :: OpacityTypeScatIso
    INTEGER  :: weak_magnetism_corrections !Horowitz 2002
    INTEGER  :: ion_ion_corrections        !Horowitz 1997, Bruenn and Mezzacappa 1997
    INTEGER  :: many_body_corrections      !Horowith et al 2017
    INTEGER  :: np_isoenergetic       
                !Scat on free nucleons using isoenergetic approximation
                !Bruenn 1985
                !Mezzacappa & Bruenn (1993)
    INTEGER  :: np_non_isoenergetic             
                !Scat on free nucleons taking into account recoil,
                !nucleon final-state blocking, and special relativity
                !Reddy et al 1998, Bruenn et al. 2020
    REAL(DP) :: ga_strange                 !strange quark contributions
  END TYPE OpacityTypeScatIso

  TYPE, PUBLIC, EXTENDS(OpacityTypeScat) :: OpacityTypeScatNNS
    INTEGER  :: weak_magnetism_corrections !Horowitz 2002
    INTEGER  :: many_body_corrections      !Horowith et al 2017
    REAL(DP) :: ga_strange                 !strange quark contributions
    TYPE(ValueType_5D), POINTER :: Phi(:,:) => null ( )
       !-- remapping of Kernel(:) to index nu/nub and n/p
  CONTAINS
    FINAL :: Finalize_NNS
  END TYPE OpacityTypeScatNNS

  TYPE, PUBLIC, EXTENDS(OpacityTypeScat) :: OpacityTypeScatNES
    INTEGER :: NPS !include neutrino-positron scattering as well 
  END TYPE OpacityTypeScatNES

  PUBLIC :: AllocateOpacity
  PUBLIC :: DeallocateOpacity
  PUBLIC :: DescribeOpacity

  INTERFACE AllocateOpacity
    MODULE PROCEDURE AllocateOpacityTypeEmAb
    MODULE PROCEDURE AllocateOpacityTypeScat
  END INTERFACE AllocateOpacity

  INTERFACE DeallocateOpacity
    MODULE PROCEDURE DeallocateOpacityTypeEmAb
    MODULE PROCEDURE DeallocateOpacityTypeScat
  END INTERFACE DeallocateOpacity

  INTERFACE DescribeOpacity
    MODULE PROCEDURE DescribeOpacityTypeEmAb
    MODULE PROCEDURE DescribeOpacityTypeScat
  END INTERFACE DescribeOpacity

CONTAINS


  SUBROUTINE AllocateOpacityTypeEmAb( Opacity, nPoints, nOpacities )

    TYPE(OpacityTypeEmAb), INTENT(inout) :: Opacity
    INTEGER,               INTENT(in)    :: nPoints(4)
    INTEGER,               INTENT(in)    :: nOpacities

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nPoints    = nPoints

    !just electron capture on nuclei for now!
    Opacity % EC_table_nOpacities = 1
  
    ALLOCATE( Opacity % Names  (nOpacities) )
    ALLOCATE( Opacity % Units  (nOpacities) )
    ALLOCATE( Opacity % Offsets(nOpacities) )
    ALLOCATE( Opacity % Opacity(nOpacities) )

    DO i = 1, nOpacities

      ALLOCATE( Opacity % Opacity(i) % Values &
              ( nPoints(1), nPoints(2), nPoints(3), nPoints(4)) )
    END DO


    IF(Opacity % nuclei_EC_table .gt. 0) THEN

      ALLOCATE( Opacity % EC_table_Names        (Opacity % EC_table_nOpacities) )
      ALLOCATE( Opacity % EC_table_Units        (Opacity % EC_table_nOpacities) )
      ALLOCATE( Opacity % EC_table_spec_Offsets (Opacity % EC_table_nOpacities) )
      ALLOCATE( Opacity % EC_table_rate_Offsets (Opacity % EC_table_nOpacities) )
      ALLOCATE( Opacity % EC_table_spec         (Opacity % EC_table_nOpacities) )
      ALLOCATE( Opacity % EC_table_rate         (Opacity % EC_table_nOpacities) )

    ENDIF
  END SUBROUTINE AllocateOpacityTypeEmAb


  SUBROUTINE DeallocateOpacityTypeEmAb( Opacity )

    TYPE(OpacityTypeEmAb), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Opacity(i) % Values )
    END DO

    DEALLOCATE( Opacity % Opacity )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Names )

    !just electron capture on nuclei for now!
    Opacity % EC_table_nOpacities = 1

    IF(Opacity % nuclei_EC_table .gt. 0) THEN

      DEALLOCATE( Opacity % EC_table_spec_Offsets )
      DEALLOCATE( Opacity % EC_table_rate_Offsets )
      DEALLOCATE( Opacity % EC_table_Units )
      DEALLOCATE( Opacity % EC_table_Names )

      DO i = 1, Opacity % EC_table_nOpacities
        DEALLOCATE( Opacity % EC_table_spec(i) % Values )
        DEALLOCATE( Opacity % EC_table_rate(i) % Values )
      END DO
      DEALLOCATE(Opacity % EC_table_spec)
      DEALLOCATE(Opacity % EC_table_rate)

    ENDIF

  END SUBROUTINE DeallocateOpacityTypeEmAb
  

  SUBROUTINE DescribeOpacityTypeEmAb( Opacity )

    USE MPI

    TYPE(OpacityTypeEmAb), INTENT(in) :: Opacity

    INTEGER :: i

    CALL MPI_INITIALIZED(mpi_is_init, ierr)

    IF(mpi_is_init) THEN
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    ELSE
      myid = 0
    ENDIF

    IF(myid == 0) THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') ' ', 'Opacity Type EmAb on nucleons and nuclei'
      WRITE(*,'(A4,A)') ' ', '--------------'
      WRITE(*,'(A6,A13,I3.3)') &
        ' ', 'nOpacities = ', Opacity % nOpacities
      WRITE(*,'(A6,A13,4I5.4)') &
        ' ', 'nPoints    = ', Opacity % nPoints
      WRITE(*,'(A6,A13,I10.10)') &
        ' ', 'DOFs       = ', &
        Opacity % nOpacities * PRODUCT( Opacity % nPoints )

      DO i = 1, Opacity % nOpacities
        WRITE(*,*)
        WRITE(*,'(A6,A8,I3.3,A3,A)') &
          ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
        WRITE(*,'(A8,A12,A)') &
          ' ', 'Units     = ', TRIM( Opacity % Units(i) )
        WRITE(*,'(A8,A12,ES12.4E3)') &
          ' ', 'Min Value = ', MINVAL( Opacity % Opacity(i) % Values )
        WRITE(*,'(A8,A12,ES12.4E3)') &
          ' ', 'Max Value = ', MAXVAL( Opacity % Opacity(i) % Values )
        WRITE(*,'(A8,A12,ES12.4E3)') &
          ' ', 'Offset    = ', Opacity % Offsets(i)
      END DO

      IF(Opacity % nuclei_EC_table .gt. 0) THEN
        DO i = 1, Opacity % EC_table_nOpacities
          WRITE(*,*)
          WRITE(*,*) 'EC table:'
          WRITE(*,'(A6,A8,I3.3,A3,A)') &
            ' ', 'Opacity(',i,'): ', TRIM( Opacity % EC_table_Names(i) )
          WRITE(*,'(A8,A17,A)') &
            ' ', 'Units          = ', TRIM( Opacity % EC_table_Units(i) )
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Min Value spec = ', MINVAL( Opacity % EC_table_spec(i) % Values )
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Max Value spec = ', MAXVAL( Opacity % EC_table_spec(i) % Values )
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Offset spec    = ', Opacity % EC_table_spec_Offsets
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Min Value rate = ', MINVAL( Opacity % EC_table_rate(i) % Values )
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Max Value rate = ', MAXVAL( Opacity % EC_table_rate(i) % Values )
          WRITE(*,'(A8,A17,ES12.4E3)') &
            ' ', 'Offset rate    = ', Opacity % EC_table_rate_Offsets
        ENDDO
      ENDIF

      WRITE(*,*) 'EmAb on np, Full kinematics, Fischer et al 2020                        ', Opacity % np_FK
      WRITE(*,*) 'EmAb on np, Full kinematics inverse neutron decay, Fischer et al 2020  ', Opacity % np_FK_inv_n_decay
      WRITE(*,*) 'EmAb on np, isoenergetic, Bruenn 1985                                  ', Opacity % np_isoenergetic
      WRITE(*,*) 'EmAb on np, Reddy et al 1998                                           ', Opacity % np_non_isoenergetic
      WRITE(*,*) 'EmAb on np, weak magnetism corrections, Horowitz 2002                  ', Opacity % np_weak_magnetism
      WRITE(*,*) 'EmAb on nuclei,   Bruenn 1985 using FFN formalism                      ', Opacity % nuclei_EC_FFN
      WRITE(*,*) 'EmAb on nuclei,   NSE-folded LMSH EC table                             ', Opacity % nuclei_EC_table
      WRITE(*,*)

    END IF

  END SUBROUTINE DescribeOpacityTypeEmAb


  SUBROUTINE AllocateOpacityTypeScat( Opacity, nPoints, nMoments, nOpacities )

    INTEGER,               INTENT(in)    :: nMoments
    INTEGER,               INTENT(in)    :: nOpacities
    INTEGER,               INTENT(in)    :: nPoints(5)
    CLASS(OpacityTypeScat), INTENT(inout) :: Opacity

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nMoments   = nMoments
    Opacity % nPoints    = nPoints

    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Offsets(nOpacities, nMoments) )
    ALLOCATE( Opacity % Kernel(nOpacities) )

    DO i = 1, nOpacities

      ALLOCATE &
        ( Opacity % Kernel(i) % Values &
            ( nPoints(1), nPoints(2), nPoints(3), nPoints(4), nPoints(5) ) )

    END DO

    SELECT TYPE ( Opacity )
    TYPE IS ( OpacityTypeScatNNS )
      Opacity % Phi ( 1:2, 1:2 )  =>  Opacity % Kernel ( 1:4 )
    END SELECT

  END SUBROUTINE AllocateOpacityTypeScat


  SUBROUTINE DeallocateOpacityTypeScat( Opacity )

    TYPE(OpacityTypeScat), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Kernel(i) % Values )
    END DO

    DEALLOCATE( Opacity % Kernel )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeScat


  SUBROUTINE DescribeOpacityTypeScat( Opacity )

    USE MPI

    CLASS(OpacityTypeScat), INTENT(in) :: Opacity

    INTEGER :: i, l

    CALL MPI_INITIALIZED(mpi_is_init, ierr)

    IF(mpi_is_init) THEN
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    ELSE
      myid = 0
    ENDIF

    IF(myid == 0) THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') ' ', 'Opacity Type Scat'
      WRITE(*,'(A4,A)') ' ', '--------------'
      WRITE(*,'(A6,A13,I3.3)') &
        ' ', 'nOpacities = ', Opacity % nOpacities
      WRITE(*,'(A6,A13,I3.3)') &
        ' ', 'nMoments   = ', Opacity % nMoments
      WRITE(*,'(A6,A13,5I5.4)') &
        ' ', 'nPoints    = ', Opacity % nPoints
      WRITE(*,'(A6,A13,I10.10)') &
        ' ', 'DOFs       = ', &
        Opacity % nOpacities * PRODUCT( Opacity % nPoints )

      DO i = 1, Opacity % nOpacities
        WRITE(*,*)
        WRITE(*,'(A6,A8,I3.3,A3,A)') &
          ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
        WRITE(*,'(A8,A12,A)') &
          ' ', 'Units     = ', TRIM( Opacity % Units(i) )

        DO l = 1, Opacity % nMoments
          WRITE(*,*)
          IF( Opacity % nMoments .eq. 4 )THEN
            WRITE(*,'(A8,A16,I3.3)') &
              ' ', 'For Moments l = ', l
            WRITE(*,'(A8,A12,ES12.4E3)') &
              ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values(:,:,l,:,:) )
            WRITE(*,'(A8,A12,ES12.4E3)') &
              ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values(:,:,l,:,:) )
          ELSE
            WRITE(*,'(A8,A16,I3.3)') &
              ' ', 'For Moments l = ', l
            WRITE(*,'(A8,A12,ES12.4E3)') &
              ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values(:,l,:,:,:) )
            WRITE(*,'(A8,A12,ES12.4E3)') &
              ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values(:,l,:,:,:) )
          END IF
          WRITE(*,'(A8,A12,ES12.4E3)') &
            ' ', 'Offset    = ', Opacity % Offsets(i,l)
        END DO ! l = nMoment
        WRITE(*,*)
      END DO
      WRITE(*,*)

      SELECT TYPE ( Opacity )

        TYPE IS ( OpacityTypeScatIso )

          WRITE(*,*) 'Weak magnetism corrections, Horowitz 2002                                   ', &
                      Opacity % weak_magnetism_corrections
          WRITE(*,*) 'Ion-ion corrections, Horowitz 1997, Bruenn and Mezzacappa 1997              ', Opacity % ion_ion_corrections
          WRITE(*,*) 'Many body corrections, Horowitz et al 2017                                  ', Opacity % many_body_corrections
          WRITE(*,*) 'Scat on np, isoenergetic, Bruenn 1985                                       ', Opacity % np_isoenergetic
          WRITE(*,*) 'Scat on np, Reddy et al 1998, Bruenn et al. 2020                            ', Opacity % np_non_isoenergetic
          WRITE(*,*) 'Strange quark contribution to axial vector coupling constant (0.0 if no-op) ', Opacity % ga_strange

        TYPE IS ( OpacityTypeScatNNS )

          WRITE(*,*) 'Weak magnetism corrections, Horowitz 2002                                   ', &
                      Opacity % weak_magnetism_corrections
          WRITE(*,*) 'Many body corrections, Horowitz et al 2017                                  ', Opacity % many_body_corrections
          WRITE(*,*) 'Strange quark contribution to axial vector coupling constant (0.0 if no-op) ', Opacity % ga_strange

        TYPE IS ( OpacityTypeScatNES )

          WRITE(*,*) 'Neutrino-positron scattering included ', Opacity % NPS

      END SELECT

    END IF

  END SUBROUTINE DescribeOpacityTypeScat


  SUBROUTINE Finalize_NNS ( Opacity )

    type ( OpacityTypeScatNNS ), intent ( inout ) :: &
      Opacity

    NULLIFY ( Opacity % Phi )

  END SUBROUTINE Finalize_NNS


END MODULE wlOpacityFieldsModule
