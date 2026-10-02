C
C  #####  #    #  #####   #####    ####    ####   ######  ######   ####
C    #    #    #  #    #  #    #  #    #  #    #  #       #       #
C    #    #    #  #    #  #####   #       #    #  #####   #####    ####
C    #    #    #  #####   #    #  #       #    #  #       #            #   ###
C    #    #    #  #   #   #    #  #    #  #    #  #       #       #    #   ###
C    #     ####   #    #  #####    ####    ####   ######  #        ####    ###
c
C
c...................................................
C turbcoefs.f > hydrosealt.f Drs. Luis San Andres & Zhou Yang, TexasA&MUniv. 1993
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR I
c Technical monitor: Mr. James Walker, NASA Lewis Research Center
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Stantard Bulk-flow Model for HYDROxxx codes from Tribology Group
c at Texas A&M University. First edited on 1992.
c
C input data for turbulence parameters and surface roughness coefficients
c and subprograms for all zeroth and first order turbulent shear 
c flow coefficients.
c
c This version for friction factors based on Moody's formulae.
c
c COntains the following SUBS:
c InputRGRS, InputMoody, FNKX,FNKY, FNGAMMA, CALCKX1, CALCKY1
c                        for variable properties liquids.
C
C *****************************************************************************
C **                                                                         **
C **  Subroutine InputRGRS                                                   **
C **                                                                         **
C **  INPUTRGRS: keyboard input changes in surface roughness coefficients    **
C **                                                                         **
C *****************************************************************************
 
      SUBROUTINE INPUTRGRS

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO, DUMY
C ----------------------------------------------------------------------------
C --  INPUTRGRS code                                                       --
C ----------------------------------------------------------------------------

  903 write (6, *) '------------------------------------------------'
      write (6, *) '(3): Relative ROUGHNESS, MAX: 10% Radial CLEAR'
      write (6, *) '------------------------------------------------'
      write (6, *)
      CALL BEEPER

  140 WRITE (6, *) 'RUGR (REAL*8): Relative Rotor surface roughness' 
      WRITE (6, 141) RUGR

  141 FORMAT ('$', 'ENTER Rugr [Default (Rugr)= ',E12.5E2,']: ')

      CALL ENTERVAL(RUGR)
      RUGR=DABS(RUGR)

      IF (Rugr.GT.(0.1D0)) THEN
         write (6, *) '----------------------------------------'
         WRITE (6, *) 'ERROR: Rotor Surface roughness TOO LARGE'
         write (6, *) '       MAX ALLOWED = 10% of clearance   '
         write (6, *) ' RESULTS may be erroneous for large RUGR'
         write (6, *) '----------------------------------------'
         CALL BEEPER
         RUGR=0.0D0
         GOTO 140
      END IF

C.....................

  143 WRITE (6, *) 'RUGS (REAL*8): Relative Stator surface roughness'
      WRITE (6, 144) RUGS


  144 FORMAT ('$', 'ENTER Rugs [Default (Rugs)= ',E12.5E2,']: ')

      CALL ENTERVAL(RUGS)
      RUGS=DABS(RUGS)

      IF (RUGS.GT.(0.15D0)) THEN
         write (6, *) '-----------------------------------------'
         WRITE (6, *) 'ERROR: BearingSurface roughness TOO LARGE'
         write (6, *) '       MAX ALLOWED = 10% of clearance   '
         write (6, *) ' RESULTS may be erroneous for large RUGs'
         write (6, *) '----------------------------------------'
         CALL BEEPER
         RUGS=0.0D0
         GOTO 143
      END IF

c...............................................................

      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine InputMoody                                                  **
C **                                                                         **
C **  INPUTMOODY: keyboard input changes in Moody's friction coefficients    **
C **                                                                         **
C *****************************************************************************
 
      SUBROUTINE INPUTMOODY

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO, DUMY
C ----------------------------------------------------------------------------
C --  INPUTMOODY CODE                                                      --
C ----------------------------------------------------------------------------

  909 write (6, *) '------------------------------------------------'
      write (6, *) '(9):MODIFY Coefs. for Friction Factor formulae  '
      write (6, *) 
      write (6, *) '       f = A { 1 + [ 10**4 r/h + B/Re# ]**EXPO }'
      write (6, *) 'typical:                                        '
      write (6, *) 'A=0.001375,B=5.0E+5,EXPO=1/3=0.333 to 1/2.65=0.377'
      write (6, *) '------------------------------------------------'
      write (6, *)
      CALL BEEPER

      WRITE (6, 171) Amod
 171  FORMAT ('$', 'ENTER Coef:  A [Default (A)= ',E12.5E2,']: ')
      CALL ENTERVAL(AMOD)
      AMOD=DABS(AMOD)

      WRITE (6, 172) Bmod
 172  FORMAT ('$', 'ENTER Coef:  B  [Default (B)= ',E12.5E2,']: ')

      CALL ENTERVAL(BMOD)
      BMOD=DABS(BMOD)

      WRITE (6, 173) EXPO 
 173  FORMAT ('$', 'ENTER Coef:EXPO [Default(EXPO)=',E12.5E2,']: ')

      CALL ENTERVAL(EXPO)
      EXPO=DABS(EXPO)

      END

C *****************************************************************************
C **                                                                         **
C **  Function Fnkx                                                          **
C **                                                                         **
C **  FNKX:  Calculates shear factor Kx for turbulent flow, Sp: 2.0*Mspeed.  **
C **                                                                         **
C *****************************************************************************

      DOUBLE PRECISION FUNCTION FNKX(U,V,H,SP,REPTYP,RHO,EMU,KXJ,KXB)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /MOODY/ AMOD, BMOD, RUGRU, RUGSU, EXPO
      DOUBLE PRECISION AMOD, BMOD, RUGRU, RUGSU, EXPO

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION U,H,V,SP,REP,KXJ,KXB, REPTYP,RHO,EMU,
     +                 U2, V2, RS, RR, CR, CS, FR, FS, DUMY
C ----------------------------------------------------------------------------
C --  FNKX code                                                             --
C ----------------------------------------------------------------------------
                                 
      U2=U*U                                   ! TURBULENCE CORRECTION FACTORS 
      V2=V*V                                   ! use Moody's Formula:
      REP=REPTYP*RHO/EMU		       ! SP=OmegaxR
      RS=DSQRT(U2+V2)*H*REP                    ! Stator Reynolds number
      RR=DSQRT((U-SP)*(U-SP)+V2)*H*REP         ! Rotor Reynolds number

          CR=RUGRU*10000.0D0/H                 ! Rotor roughness coef.
          CS=RUGSU*10000.0D0/H                 ! Stator roughness coef.
          FR=AMOD*(1.0D0+(CR+BMOD/RR)**EXPO)   ! rotor friction factor
          FS=AMOD*(1.0D0+(CS+BMOD/RS)**EXPO)   ! stator friction factor

      KXJ=DABS(RR*FR)                          !
      KXB=DABS(RS*FS)                          !
      FNKX=0.5D0*(KXJ+KXB)                     !

      DUMY=12.0D0                              !
      FNKX=EMU*DMAX1(DUMY, FNKX)               ! Select laminar or turbulent
      KXJ= EMU*DMAX1(DUMY, KXJ)                !
      KXB= EMU*DMAX1(DUMY, KXB)                !

      END
c
c

C *****************************************************************************
C **                                                                         **
C **  Function Fnky                                                          **
C **                                                                         **
C **  FNKY:  Calculates shear factor Ky turbulent flow.                      **
C **                                                                         **
C *****************************************************************************

      DOUBLE PRECISION FUNCTION FNKY(U, V, H, SP, REPTYP,RHO,EMU )

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /MOODY/ AMOD, BMOD, RUGRU, RUGSU, EXPO

      DOUBLE PRECISION AMOD, BMOD, RUGRU, RUGSU, EXPO

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
     
      DOUBLE PRECISION U, V, H, SP, REP, REPTYP, RHO, EMU,
     +                 U2, V2, RS, RR, CR, CS, FR, FS, DUMY
C ----------------------------------------------------------------------------
C --  FNKY code                                                             --
C ----------------------------------------------------------------------------

      U2=U*U                                   ! Turbulent friction factors
      V2=V*V                                   ! Use Moody's formula
      REP=REPTYP*RHO/EMU                       !
      RS=DSQRT(U2+V2)*H*REP                    ! Stator Reynold's number
      RR=DSQRT((U-SP)*(U-SP)+V2)*H*REP         ! Rotor Reynold's number
          CR=RUGRU*10000.0D0/H                 ! Rotor roughness coef.
          CS=RUGSU*10000.0D0/H                 ! Stator roughness coef.
          FR=AMOD*(1.0D0+(CR+BMOD/RR)**EXPO)   ! rotor friction factor
          FS=AMOD*(1.0D0+(CS+BMOD/RS)**EXPO)   ! stator friction factor
      FNKY=0.5D0*DABS(FS*RS+FR*RR)             ! Turbulent shear coefficient

      DUMY=12.0D0                              !
      FNKY=EMU*DMAX1(DUMY, FNKY)               ! Select laminar of turbulent

      END 


C *****************************************************************************
C **                                                                         **
C **  Function Fngama                                                        **
C **                                                                         **
C **  FNGAMA:  Calculates shear stress at moving wall, S: Speed parameter/2  **
C **                                                                         **
C *****************************************************************************

      DOUBLE PRECISION FUNCTION FNGAMA(U, V, H, SP, REPTYP, RHO, EMU)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO

      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION U, V, H, REP, REPTYP, RHO, EMU,
     +       FR, FS, SP, USP, U2, V2, RS, RR, CR, CS, DUMY

C ----------------------------------------------------------------------------
C --  FNGAMA code                                                           --
C ----------------------------------------------------------------------------

      IF (H.LE.(0.0D0)) THEN     
          WRITE (6, *) 'FNGAMA: H<=0.0 STOP'
          STOP                          !
      END IF                            !

      USP=U-SP                          ! Speed parameter
      U2=U*U                            ! Turbulence correction factors
      V2=V*V                            !  Use Moody's formula
      REP=REPTYP*RHO/EMU                !
      RS=DSQRT(U2+V2)*H*REP             ! Stator Reynold's number
      RR=DSQRT(USP*USP+V2)*H*REP        ! Rotor Reynold's number
          CR=RUGR*10000.0D0/H           ! Rotor roughness coef.
          CS=RUGS*10000.0D0/H           ! Stator roughness coef.
      FR=AMOD*(1.0D0+(CR+BMOD/RR)**EXPO)! rotor friction factor
      FS=AMOD*(1.0D0+(CS+BMOD/RS)**EXPO)! stator friction factor

      FNGAMA=(U*RS*FS-USP*RR*FR)/4.0D0

      IF ((RR*FR).LT.(12.0D0)) THEN     ! for laminar flow
         FNGAMA=SP
      END IF

      FNGAMA=EMU*FNGAMA/H    

      END 

c
C ----------------------------------------------------------------------------


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calckx1                                                     **
C **                                                                         **
C **  CALCKX1:  Calculates shear factors Kx and first order coefficients.    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCKX1(U,V,H,SP,REPT,dpx,rho,emu,GUO,GUU,GUV,GUR,GUM)

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION U, V, H, SP, REP, REPT, dpx, KXO, GUO, GUU, GUV,
     +                 U2, V2, UR, UR2, RS, RR, CR, CS, BR,BS,FRO,FSO,
     +                 KXR, DUMY, GR, GS, FR1, FS1, CCR, CCS, GUR, GUM,
     +                 rho, emu, expm1
C ----------------------------------------------------------------------------
C --  CALCKX1 code                                                          --
C ----------------------------------------------------------------------------
      U2=U*U                                     ! Turbulence correction factors
      V2=V*V                                     ! Use Moody's formula:
      UR=(U-SP)                                  !
      UR2=UR*UR                                  !
      REP=REPT*rho/emu                           !
      RS=H*DSQRT(U2+V2)                          ! Stator Reynolds number/rep
      RR=H*DSQRT(UR2+V2)                         ! Rotor Reynolds number/rep
      CR=RUGR*10000.0D0/H                        ! Rotor roughness coeff.
      CS=RUGS*10000.0D0/H                        ! Stator roughness coeff.
      BR=BMOD/RR/REP                             !
      BS=BMOD/RS/REP                             !
      FRO=AMOD*(1.0D0+(CR+BR)**EXPO)             ! Zeroth rotor friction factor
      FSO=AMOD*(1.0D0+(CS+BS)**EXPO)             !        stator friction factor
      KXR=REP*RR*FRO                             !
      KXO=(REP*RS*FSO+KXR)/2.0D0                 ! Zeroth turbulent shear coeffs
      DUMY=12.0D0                                ! 
      KXO=DMAX1(DUMY, KXO)                       ! Select laminar or turbulent
      KXR=DMAX1(DUMY, KXR)                       !

      IF (KXO.GT.DUMY) THEN                      ! Turbulent flow condition
          expm1=1.0D0/expo-1.0D0                       !
          GR=-AMOD*EXPO/(DABS(FRO/AMOD-1.0D0))**expm1  !............
          GS=-AMOD*EXPO/(DABS(FSO/AMOD-1.0D0))**expm1  !
          FR1=REP*H*(FRO+GR*BR)/RR/2.0D0         !
          FS1=REP*H*(FSO+GS*BS)/RS/2.0D0         !
          CCR=GR*(REP*RR*CR+BMOD)/2.0D0          !
          CCS=GS*(REP*RS*CS+BMOD)/2.0D0          !
          GUO=U*(-KXO+CCR+CCS)-SP*(-KXR/2.0D0+CCR) ! coeffic. of cosx/sinx
          GUO=GUO*Emu/H/H                        !
          GUU=(FS1*U2+FR1*UR2+KXO/H)*Emu         ! coeffic. of U1/U2
          GUV=(FS1*U+FR1*UR)*V*Emu               ! coeffic. of V1/V2
          Dumy=Bmod*(U*(GR+GS)-SP*GR)/2.0D0/H    !
          GUR=(Dumy*emu-H*dpx)/Rho               ! coeffic of RHO1/RHO2
          GUM=-Dumy                              ! coeffic of EMU1/EMU2
          RETURN                                 !
      END IF                                     !

      GUO=(-2.0D0*U+SP)*12.0D0*emu/H/H           ! Kxo=12, Kxor=12
      GUU=12.0D0*Emu/H                           ! Laminar flow
      GUV=0.0D0                                  !
      GUR=0.0D0                                  !
      GUR=-(12.0D0*Emu*(U-SP/2.0D0)/H+H*dpx)/Rho ! corrected 2/22/95
      GUM=12.0D0*(U-SP/2.0D0)/H                  ! GUR=0 if inertialess
      RETURN                                     !

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcky1                                                     **
C **                                                                         **
C **  CALCKY1:  Calculates shear factors Ky and first order coefficients.    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCKY1(U,V,H,SP,REPT,dpy,rho,emu,GVO,GVU,GVV,GVR,GVM)

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION U,V,H,SP,REPT,dpy,REP, KYO,GVO,GVU,GVV,GVR,GVM,
     +                 U2, V2, UR, UR2, RS, RR, CR, CS, BR,BS,FRO,FSO,
     +                 DUMY,GR,GS,FR1,FS1,CCR,CCS,rho, emu, expm1
C ----------------------------------------------------------------------------
C --  CALCKY1 code                                                          --
C ----------------------------------------------------------------------------
      U2=U*U                                     ! Turbulence correction factors
      V2=V*V                                     ! Use Moodys formula:
      UR=(U-SP)                                  !
      UR2=UR*UR                                  !

      REP=REPT*rho/emu                           !
      RS=H*DSQRT(U2+V2)                          ! Stator Reynolds number/rep
      RR=H*DSQRT(UR2+V2)                         ! Rotor Reynolds number/rep
      CR=RUGR*10000.0D0/H                        ! Rotor roughness coeff.
      CS=RUGS*10000.0D0/H                        ! Stator roughness coeff.
      BR=BMOD/RR/REP                             !
      BS=BMOD/RS/REP                             !
      FRO=AMOD*(1.0D0+(CR+BR)**EXPO)             ! Zeroth rotor friction factor
      FSO=AMOD*(1.0D0+(CS+BS)**EXPO)             !        stator friction factor
      KYO=REP*(RR*FRO+RS*FSO)/2.0D0              ! Zeroth turbulent shear coefs
      DUMY=12.0D0                                ! 
      KYO=DMAX1(DUMY, KYO)                       ! Select Laminar or turbulent


      IF (KYO.GT.DUMY) THEN                      ! Turbulent flow condition
          expm1=1.0D0/expo-1.0D0                       !
          GR=-AMOD*EXPO/(DABS(FRO/AMOD-1.0D0))**expm1    !
          GS=-AMOD*EXPO/(DABS(FSO/AMOD-1.0D0))**expm1    !
          FR1=REP*H*(FRO+GR*BR)/RR/2.0D0                 !
          FS1=REP*H*(FSO+GS*BS)/RS/2.0D0                 !
          CCR=GR*(REP*RR*CR+BMOD)/2.0D0                  !
          CCS=GS*(REP*RS*CS+BMOD)/2.0D0                  !
          GVO=V*Emu*(-KYO+CCR+CCS)/H/H                   ! coeffic. of cosx/sinx
          GVU=(FS1*U+FR1*UR)*V*Emu                       ! coeffic. of U1/U2
          GVV=((FS1+FR1)*V2+KYO/H)*Emu                   ! coeffic. of V1/V2
          Dumy=Bmod*V*(GR+GS)/2.0D0/H                    !
          GVR=(Dumy*Emu-H*dpy)/Rho                       ! coeffic of Rho1/Rho2
          GVM=-Dumy                                      ! coeffic of Emu1/EMu2
          RETURN                                         !
      END IF                                             !.....................


      GVO=(-24.0D0*emu*V)/H/H                    ! Kyo=12.0*emu
      GVU=0.0D0                                  ! Laminar flow
      GVV=12.0D0*emu/H                           !
      GVR=-(12.0D0*emu*V/H+H*dpy)/Rho            ! corrected 2/22/95
      GVM=12.0D0*V/H                             ! GVR=0 if inertialess
      RETURN                                     !

      END 

C ----------------------------------------------------------------------------
C --  turbcoefs.f code                                                  --
C ----------------------------------------------------------------------------