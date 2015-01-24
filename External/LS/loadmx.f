C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         ALOADMX.FOR
C    MODULE:       LOADMX
C    TYPE:         LOADMX
C
C    PURPOSE:      LOAD THE LOOK-UP TABLE FOR THE MAXWELL CONSTRUCTION
C
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/16/90
C
C    CALL LINE:    CALL LOADMX
C
C    INPUTS:       N/A
C
C    OUTPUTS       N/A
C
C    SUBROUTINE CALLS: EOS_M4C
C
C    INCLUDE FILES:  EOS_M4C.INC, MAXWEL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE LOADMX(FilePath,lscompress)
C
C
      USE eos_m4c_module
      USE maxwel_module
C
      IMPLICIT NONE
C
C
      INTEGER NTMP, NYE, NYE2, NUM_BP
      INTEGER LUN1, LUN2, KK, KMIN
      PARAMETER(LUN1=54,LUN2=55)
C
C
      INTEGER FNML1, FNML2
      CHARACTER*3   lscompress
      CHARACTER*128 FilePath
      CHARACTER*128 FNAME1, FNAME2
      CHARACTER*128 TERM_STR
      logical loaddata
C
      DOUBLE PRECISION N_SM, SYMM_M, COMP_M, BINDEM, SYM_SM, SIG_SM
      DOUBLE PRECISION N_SB, SYMM_B, COMP_B, BINDEB, SYM_SB, SIG_SB
      DOUBLE PRECISION MSCDD3, BSCDD3

      DOUBLE PRECISION :: send_buf1(12+2*NUMYE+2*NUMYE*NUMTMP)
      DOUBLE PRECISION :: send_buf2(7+2*NUMYE*NBPNTS)
      INTEGER i_extent
C
      INCLUDE 'force.inc'
      data loaddata /.false./
C
C
          if ( loaddata ) return
      FNAME1        = TRIM(FilePath)//'/max'//
     &      lscompress//'.atb'
      FNML1         = 42
      FNAME2        = TRIM(FilePath)//'/bd'//
     &      lscompress//'.atb'
      FNML2         = 41
      loaddata      = .true.
      WRITE(*,*) "Reading LS EoS data files:"
      WRITE(*,*) FNAME1
      WRITE(*,*) FNAME2

      i_extent = 12+2*NUMYE+2*NUMYE*NUMTMP

!     IF ( myid .eq. 0) THEN
C
C
C
C-----------------------------------------------------------------------
C        Read the file Maxwell construction data file
C-----------------------------------------------------------------------
C
C
C
      OPEN(UNIT=LUN1,FILE=FNAME1(1:FNML1),STATUS='OLD')
C
C
C
C
C
C
      READ(LUN1,*) N_SM, SYMM_M
      READ(LUN1,*) COMP_M,BINDEM
      READ(LUN1,*) SYM_SM, SIG_SM

      send_buf1(1) = N_SM
      send_buf1(2) = SYMM_M
      send_buf1(3) = COMP_M
      send_buf1(4) = BINDEM
      send_buf1(5) = SYM_SM
      send_buf1(6) = SIG_SM
C
C
C
C
      READ(LUN1,*) NTMP,NYE
      READ(LUN1,*) T_LOW,T_HI
      READ(LUN1,*) Y_LOW,Y_HI

      send_buf1(7) = T_LOW
      send_buf1(8) = T_HI
      send_buf1(9) = Y_LOW
      send_buf1(10) = Y_HI
C
C
C
      IF((NTMP.NE.NUMTMP).OR.(NYE.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  MXWL TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
C
      DO 101 J=1,NUMYE,1
        DO 100 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYLOW(KK,J),KK=I,KMIN,1)
 100    CONTINUE
 101  CONTINUE
C
C
      DO 103 J=1,NUMYE,1
        DO 102 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYHI(KK,J),KK=I,KMIN,1)
 102    CONTINUE
 103  CONTINUE
C
C
C
      DO 104 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (T_H(KK),KK=I,KMIN,1)
 104  CONTINUE
C
C
      DO 105 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (D_H(KK),KK=I,KMIN,1)
 105  CONTINUE
C
      READ(LUN1,*) YCUT
      READ(LUN1,*) MSCDD3
C
      send_buf1(11) = YCUT
      send_buf1(12) = MSCDD3
C
      CLOSE(UNIT=LUN1,STATUS='KEEP')

      KK = 13
      DO I = 1,NUMYE
        send_buf1(KK) = T_H(I)
        KK = KK + 1
      ENDDO
      DO I = 1,NUMYE
        send_buf1(KK) = D_H(I)
        KK = KK + 1
      ENDDO
      DO J=1,NUMYE
        DO I=1,NUMTMP
          send_buf1(KK) = BRYLOW(I,J)
          KK = KK + 1
        ENDDO
      ENDDO
      DO J=1,NUMYE
        DO I=1,NUMTMP
          send_buf1(KK) = BRYHI(I,J)
          KK = KK + 1
        ENDDO
      ENDDO
      KK = KK - 1
C
      IF (i_extent .NE. KK ) THEN
       TERM_STR = "Miss filled send_buf1 in loadmx."
       WRITE(*,*) TERM_STR
!      CALL TERMINATE( TERM_STR )
       STOP
      ENDIF

!     ENDIF !myid == 0

!     CALL MPI_BCAST( send_buf1 , i_extent, MPI_DOUBLE_PRECISION, 0,
!    &                     MPI_COMM_WORLD, ierr)


      N_SM = send_buf1(1) 
      SYMM_M = send_buf1(2) 
      COMP_M = send_buf1(3) 
      BINDEM = send_buf1(4) 
      SYM_SM = send_buf1(5) 
      SIG_SM = send_buf1(6) 
      T_LOW = send_buf1(7) 
      T_HI = send_buf1(8) 
      Y_LOW = send_buf1(9) 
      Y_HI = send_buf1(10) 
      YCUT = send_buf1(11) 
      MSCDD3 = send_buf1(12) 
      KK = 13
      DO I = 1,NUMYE
        T_H(I) = send_buf1(KK)
        KK = KK + 1
      ENDDO
      DO I = 1,NUMYE
        D_H(I) = send_buf1(KK)
        KK = KK + 1
      ENDDO
      DO J=1,NUMYE
        DO I=1,NUMTMP
          BRYLOW(I,J) = send_buf1(KK)
          KK = KK + 1
        ENDDO
      ENDDO
      DO J=1,NUMYE
        DO I=1,NUMTMP
          BRYHI(I,J) = send_buf1(KK)
          KK = KK + 1
        ENDDO
      ENDDO
      KK = KK - 1
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  MAXWELL CON. TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C        Read the file Boundary data file
C-----------------------------------------------------------------------
C
C
!     IF ( myid .eq. 0 ) THEN
C
      OPEN(UNIT=LUN2,FILE=FNAME2(1:FNML2),STATUS='OLD')
C
C
C
C
      READ(LUN2,*) N_SB,SYMM_B
      READ(LUN2,*) COMP_B,BINDEB
      READ(LUN2,*) SYM_SB,SIG_SB
C
C
C
C
C
C
      READ(LUN2,*) NUM_BP,NYE2
      READ(LUN2,*) LNL,LNH,LNC
      READ(LUN2,*) Y_LOW2,Y_HI2
C
C
      IF((NBPNTS.NE.NUM_BP).OR.(NYE2.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  BNDY TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
      IF(ABS(LNL-LNLOW).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNH-LNHI).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNC-LNCUT).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  MID CUT OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_LOW-Y_LOW2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_HI-Y_HI2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
      send_buf2(1) = N_SB
      send_buf2(2) = SYMM_B
      send_buf2(3) = COMP_B
      send_buf2(4) = BINDEB
      send_buf2(5) = SYM_SB
      send_buf2(6) = SIG_SB
C
      DO 201 J=1,NUMYE,1
        DO 200 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (LBOUND(KK,J),KK=I,KMIN,1)
 200    CONTINUE
 201  CONTINUE
C
C
      DO 203 J=1,NUMYE,1
        DO 202 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (UBOUND(KK,J),KK=I,KMIN,1)
 202    CONTINUE
 203  CONTINUE
C
      READ(LUN2,*) BSCDD3
      send_buf2(7) = BSCDD3
C
      IF(ABS(MSCDD3-BSCDD3).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  SCRDD3 VALUES ARE INCONSIST.'
        STOP
      ENDIF
C
      KK = 8
      DO J=1,NUMYE
        DO I=1,NBPNTS
          send_buf2(KK) = LBOUND(I,J)
          KK = KK + 1
        ENDDO
      ENDDO
      DO J=1,NUMYE
        DO I=1,NBPNTS
          send_buf2(KK) = UBOUND(I,J)
          KK = KK + 1
        ENDDO
      ENDDO
C
      KK = KK - 1
C
      CLOSE(UNIT=LUN2,STATUS='KEEP')
C
!     ENDIF   ! myid == 0

      i_extent = 7+2*NUMYE*NBPNTS
      IF (i_extent .NE. KK ) THEN !.AND. myid .eq. 0) THEN
       WRITE(*,*) "Miss filled send_buf2 in loadmx."
       STOP
      ENDIF

!     CALL MPI_BCAST( send_buf2 , i_extent, MPI_DOUBLE_PRECISION, 0,
!    &                     MPI_COMM_WORLD, ierr)


      N_SB = send_buf2(1)
      SYMM_B = send_buf2(2)
      COMP_B = send_buf2(3)
      BINDEB = send_buf2(4)
      SYM_SB = send_buf2(5)
      SIG_SB = send_buf2(6)
      BSCDD3 = send_buf2(7)
      KK = 8
      DO J=1,NUMYE
        DO I=1,NBPNTS
          LBOUND(I,J) = send_buf2(KK)
          KK = KK + 1
        ENDDO
      ENDDO
      DO J=1,NUMYE
        DO I=1,NBPNTS
          UBOUND(I,J) = send_buf2(KK)
          KK = KK + 1
        ENDDO
      ENDDO
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  BOUNDARY TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C                  All arrays are now loaded so return
C-----------------------------------------------------------------------
C
      SCRDD3 = BSCDD3
      N_S = N_SM
      NSUBS = N_SM
      SYMM = SYMM_M
      COMP = COMP_M
      BIND_E = BINDEM
      SYM_S = SYM_SM
      SIG_S = SIG_SM
C
c20      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
c20      DD = (COMP+2.0*SKYRMC)/(3.0*SKYRMC+9.0*BIND_E)
c20      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
c20      AA = (OVR23*SKYRMC-DD*(SKYRMC+BIND_E))/(N_S*(DD-1.0))-BB
c20      CC = (COMP+2.0*SKYRMC)/(9.0*DD*(DD-1.0)*N_S**DD)
      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
      DD = (COMP+2.0*SKYRMC+SCRDD3*(COMP-4.0*SKYRMC-18.0*BIND_E))/
     1    ((1.0-SCRDD3)*(3.0*SKYRMC+9.0*BIND_E))
      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
      CC = ((OVR3*SKYRMC+BIND_E)*(1.0+SCRDD3)**2)/((N_S**DD)*(DD-1.0))
      AA = ((OVR23*SKYRMC-DD*(SKYRMC+BIND_E)-SCRDD3*(OVR3*SKYRMC+BIND_E)
     1    )/(N_S*(DD-1.0)) )-BB
      DD3 = SCRDD3/(N_S**(DD-1.0))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  SKYRME PARAMETERS FOR THIS RUN ARE:>>'
      WRITE(*,*) 'ABCD: ',AA,BB,CC,DD,SCRDD3
      WRITE(*,*) ' Satur. density, symmetry engy, & compression
     1 mod.:'
      WRITE(*,*) N_SM, SYMM_M, COMP_M
      WRITE(*,*) N_SB, SYMM_B, COMP_B
      WRITE(*,*) ' Binding engy, surf. symm. engy, & surface
     1 tension:'
      WRITE(*,*) BINDEM,SYM_SM,SIG_SM
      WRITE(*,*) BINDEB,SYM_SB,SIG_SB
C
      WRITE(*,*)
C
C
      CALL INITFERM(FilePath)
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX: FERMI INTEGRAL TABLES ARE INITIALIZED>>'
      WRITE(*,*)
C
C
      RETURN
C
      END SUBROUTINE LOADMX
