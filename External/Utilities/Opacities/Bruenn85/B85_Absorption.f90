MODULE B85_Absorption

  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
    h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, cvel, ergmev,&
      cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero, one

  IMPLICIT NONE

  PUBLIC TotalNuEAbsorption, &
         TotalNuEbarAbsorption

CONTAINS

!============================== Function ======================================

  REAL(dp) FUNCTION TotalNuEAbsorption &
    ( energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )
!------------------------------------------------------------------------------
! Purpose:
!   To compute the electron-type neutrino absorptivity.
!   (1) Absorptivity = emissivity + inverse of mean path
! Ref: 
!   Mezzacappa & Bruenn(1993) AJ 405
!   Bruenn(1985) AJSS 58
!------------------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                            xheavy, xn, xp

    REAL(dp)             :: TMeV, N, qpri, nhn, npz, etapn, midFe, &
                            midE, chem_v, mpG, mnG, fexp, &
                            rop, ron, midFexpp, midFep, midEp, midCons
    REAL(dp)             :: emitnp, absornp, emitni, absorni

    TMeV   = T * kMeV 
      N    = A - Z
    qpri   = chem_n - chem_p + 3.0_dp + dmnp  ! [MeV] 3 = energy of the 1f5/2 level 
    chem_v = chem_e + chem_p - chem_n - dmnp 
    
   
     mpG   = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
     mnG   = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]

    if(n.le.34.0)               nhn = 6.0_dp
    if(n.gt.34.0.and.n.le.40.0) nhn = 40.0_dp - N
    if(n.gt.40.0)               nhn = zero

    if(z.le.20.0)               npz = zero
    if(z.gt.20.0.and.z.le.28.0) npz = z - 20.0
    if(z.gt.28.0)               npz = 8.0
    
    
!    etapn = rho * ( xn - xp ) / ( mbG * ( FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp ) )
     rop   = rho * xp / mpG                   ! Approxiation in the nondegenerate regime
     ron   = rho * xn / mnG 

!----------------------------------------------------------------------------
! Stop the function if any of etapn/rop/ron is negative
!--------------------------------------------------------------------------- 
    IF ( ( rop < zero ) .or. ( ron < zero ) ) THEN
!    IF ( (etapn < 0.0_dp ) .or. ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
      WRITE(*,*)'xn - xp is ', xn - xp
      WRITE(*,*)'fexp term is ', FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp 
      WRITE(*,*)'chem_n - chem_p is ', chem_n - chem_p 
      STOP 
    END IF

!-----------------------------------------------------------------------------
!   emissivity and absorptivity on nuclei
!-----------------------------------------------------------------------------
    IF ( xheavy * npz * nhn == zero ) THEN
       emitni  = zero
      absorni  = zero
    ELSE
       midEp   = (energy+qpri)**2 * SQRT( & 
                    MAX( 1.0_dp - ( me / (energy+qpri) )**2, 0.0_dp ) )     
      midFexpp = FEXP( (energy+qpri-chem_e) / TMeV )
       midFep  = 1.0_dp / ( midFexpp + 1.0_dp )
      midCons  = therm2 * rho * xheavy * npz * nhn * midEp / (mbG * A)

       emitni  = midCons * midFep 
      absorni  = midCons * FEXP( (chem_n + dmnp - chem_p - qpri) /TMeV ) &
               * ( 1.0_dp - midFep) 
    END IF

!-----------------------------------------------------------------------------
!   emissivity and absorptivity on nucleons
!-----------------------------------------------------------------------------
!------------------------------------------------------------------
!  Set emitnp + absornp = zero and return if both xn and xp are zero
!------------------------------------------------------------------
    IF ( xn == zero .and. xp == zero ) THEN
      TotalNuEAbsorption = emitni + absorni
      WRITE(*,*) 'xn and xp = zero'
      RETURN
    END IF

      midFe   = 1.0_dp / ( FEXP( (energy+dmnp-chem_e) / TMeV ) + 1.0_dp )
       midE   = (energy+dmnp)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+dmnp) )**2 )

      emitnp  = therm1 * rop * midE * midFe
     absornp  = therm1 * ron * midE * ( 1.0_dp - midFe )

    TotalNuEAbsorption = ( emitni + absorni ) + ( emitnp + absornp )

    IF ( ISNAN(TotalNuEAbsorption) ) THEN
      WRITE(*,*) "TotalNuEAbsorption is NAN! "
      STOP
    END IF

    RETURN

  END FUNCTION TotalNuEAbsorption


  REAL(dp) FUNCTION TotalNuEbarAbsorption &
    ( energy, rho, T, chem_e, xn, xp )
!------------------------------------------------------------------------------
! Purpose:
!   To compute the electron-type antineutrino absorptivity.
!   (1) Absorptivity = emissivity + inverse of mean path
! Ref: 
!   Bruenn(1985) AJSS 58
!------------------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in) :: energy, rho, T, chem_e, xn, xp

    REAL(dp) :: TMeV, midE, midFe, rop, ron, mpG, mnG, fexp
    REAL(dp) :: emitnp, absornp

    IF ( energy .lt. (dmnp + me) ) THEN
      
      TotalNuEbarAbsorption = zero

    ELSE

      TMeV = T * kMeV 
      mpG  = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
      mnG  = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]
      rop  = rho * xp / mpG                    ! Approxiation in the  nondegenerate regime
      ron  = rho * xn / mnG

      midE = (energy - dmnp)**2 * ( 1.0_dp - (me/(energy-dmnp))**2 )

      midFe = 1.0_dp / ( FEXP( (energy-dmnp+chem_e) / TMeV ) + 1.0_dp )

      absornp = therm1 * rop * midE * (1.0_dp - midFe)

      emitnp  = therm1 * ron * midE * midFe

      TotalNuEbarAbsorption =  emitnp + absornp 

      IF ( ISNAN(TotalNuEbarAbsorption) ) THEN
        WRITE(*,*) "TotalNuepAbsorption is NAN! "
        STOP
      END IF

      IF ( TotalNuEbarAbsorption .le. zero ) THEN
        WRITE(*,*) "TotalNuepAbsorption is Negative! "
        STOP
      END IF

    END IF
    
    RETURN

  END FUNCTION TotalNuEbarAbsorption

END MODULE B85_Absorption
