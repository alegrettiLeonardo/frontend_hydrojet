C 7/03/95 updated for angled injection, define CJET
C 6/20/95 updated for cell depth effect LSA
C 6/28/94 updated for radial heat flow LSA
c
C NOTE: IFULL not defined until XYDATA(0) is called
C       so there is a conflict when ITPAD is first defined
C
C    #    #    #     #     #####   ####   #####    ####    #####          ######
C    #    ##   #     #       #    #    #  #    #  #          #            #
C    #    # #  #     #       #    #    #  #    #   ####      #            #####
C    #    #  # #     #       #    #    #  #####        #     #     ###    #
C    #    #   ##     #       #    #    #  #       #    #     #     ###    #
C    #    #    #     #       #     ####   #        ####      #     ###    #
C
C    #    #   #   #  #####   #####    ####        #  ######   #####
C    #    #    # #   #    #  #    #  #    #       #  #          #
C    ######     #    #    #  #    #  #    #       #  #####      #
C    #    #     #    #    #  #####   #    #       #  #          #
C    #    #     #    #    #  #   #   #    #  #    #  #          #
C    #    #     #    #####   #    #   ####    ####   ######     #
C
C  hydrojet.f Copyright Dr. Luis San Andres TexasA&MUniversity / 1995
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR III
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C
C
C =========================================================================
C Dimensionless pressures are calculated according to:
C
C  p =(P-PA)/PSA   where P: pressure [Pa]
C
C BEARING=1 (HJB) or (2) seal
C
C  PSA= PS-PA, and
C                  PA=MIN(Pleft,Pright) for hydrostatic bearings, BEARING=1
C                  PA=MIN(Pright)       for annular seals, BEARING=2
C
C
C BEARING=3 plain journal bearings or tilt-pad bearings or foil bearings
C
C IF FLUID: liquid   (IF=1,2,3,4 or 12)
C
C  PSA= EMU x OMEGA x (R/C)^2    when PS=PA
C
c  and ==> PA= MIN(Pleft, Pright)=PS
C
C IF FLUID: air  (IF=5)
C
C PSA = PS = PA   so p = P/PA
C                 and        PA = 0.0 (is then set to zero automatically)
C
C ==========================================================================
C Dimensionless fluid velocities u=U/U*, v=V/U* where
C
C     U*=PSAxC^2/(EMUxR)                ! TYP fluid velocity
C
C ==========================================================================
C For foil bearings, ac>0, PBACK = pressure on back of foil = PA
c
c
c SUBS: Cparam, Askiter, Checkprops, Locprops
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Cparam                                                      **
C **                                                                         **
C **  CPARAM:  First operations with input data, dimensionless parameters.   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CPARAM(RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /RECES/ PREC(MAXNPOCK),TREC(MAXNPOCK),QREC(MAXNPOCK),
     +               QIN,QOUT,QFACTOR
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /SPLDATA/Z(Nsl),CL(Nsl),BCL(Nsl),CCL(Nsl),DCL(Nsl),NJ
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAW, IWEAR
      COMMON /PARAM3/ EMU, RHO, PC, CD, DORIF, LOSXSI,ALPHA !RPM, PS, PA,
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr
      COMMON /LOSPAD/ LOSleadP, KLOSPad
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS

      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DEN2,VIS2,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2
      COMMON /COMPRE/ L4, PCAVI
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /PreLR/ Ple, Pri, Csel, Cser
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /FORCE0/ FFACTOR, FX, FY, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /AIR/ CP,R,AIRMU,GAMA,PR,kKS
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT
      COMMON /CONT/ IFF

      COMMON /THERMAL/ ALFT, UC, TC, EC
      COMMON /THERMID/ ISOTH
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /TISOBJD/ TBJ,TBS
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd
      COMMON /TGROOVE/ LEMDA,DELTA
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET
C..............................................................................
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
      COMMON /SWITCH/ IPROP
      COMMON /IDTHDPAD/ ITPAD
C........................................................ VARIABLE DECLARATIONS
      DOUBLE PRECISION PREC,TREC, QREC, QIN, QOUT, QFACTOR,
     +                 CP,R,AIRMU,GAMA,PR,kKS,
     +                 Z, CL, BCL, CCL, DCL, ROTPAD, TSHAFT,TSTATOR,
     +                 CLEAR, DIAM, LENGTH, LD, AR,BR,HREC,
     +                 X1,X1R,LPAD,LREC,PaDorif, LENGTHR, LENGTHL,
     +                 EXO, EYO, AXO, AYO, ZO,EWX, EWY, EWEAR, BETAW,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr,
     +                 CINLET, CEXIT,ClearO,ClearR,ClearL,YR,YL,
     +                 Pleft, Pright, Cleft, Cright,
     +                 Ple, Pri, Csel, Cser,LOSleadP, KLOSPad,
     +                 DIAORIF, CORIF, CPS,THS,BETAKS,
     +                 PRCOEFC,PRCOEFS,PLCOEFC,PLCOEFS,
     +                 ALFT, UC, TC, EC,
     +                 TSUMP,TSUMPd,RHOSU,RHOSd ,LEMDA,DELTA
      DOUBLE PRECISION RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 RHOTYP,EMUTYP,DEN2,VIS2,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2,
     +                 L4, PCAVI, HRECU, VSUP, BETA,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 AMOD, BMOD, RUGR, RUGS, EXPO,
     +                 PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 FFACTOR, FX, FY, TO, TOR, MFACTOR, MX, MY,
     +                 TEMPK, VSOUND, PRELOAD, OFFSET,ROTDEL,
     +                 AC, ETA, PBACK, LIFT, TBJ, TBS,
     +                 TBOUT, THERMALK, ROUTER, HKB,
     +                 HCELL, HCDIM, ANGLEJ, LOCJET, CJET, DPJET

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPAD, NREC, NJ, ITER, ITMAX, ITPMAX,
     +        IF, IL, IFF, IWEAR, ISYM, ICSTEP, IFULL, TILTPAD,
     +        IPRuni, IPLuni, NPRcs, NPLcs, BEARING,
     +        ISOTH ,IPROP,ITPAD

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION PI, ECC,  RADIUS, ABAR,P1, Press,Presa,
     +                 CRATIO, C2, CR2, OMEGA, DIFF,
     +                 VIS, U, Rhor, Ao, DENOM, Rhorp, Zero ,UTILT,TTILT
      INTEGER J, I, ISYMold, Zdumy ,Jstart

C ----------------------------------------------------------------------------
C --  CPARAM code                                                           --
C ----------------------------------------------------------------------------
      ISYMold=ISYM                      !
      ZERO=0.0D0                        !
      PI=DACOS(-1.0D0)                  !
      ITPAD=0                           ! =0: No update of groove temperatures
      IPROP=1                           ! YES update properties
      ICAV=1                            ! YES cavitation : IF P<Pcav -> P=Pcav
      ITURB=1                           ! YES turbulent
      PREPS=0.001D0                     ! Conv criteria for recess pressures
c    !..................................!....................................
      Ecc=DSQRT(Exo*Exo+Eyo*Eyo)        ! Magnitude of eccentricity vector
c     IF (Ecc.GE.(1.0D0)) THEN          ! Ecc > 1
c        WRITE (6, 177) EXO, EYO, Ecc   ! .......
c        CALL BEEPER                    !
c        STOP                           !
c     END IF                            !
c    !..................................!....................................

      PRATIO=DABS(PRATIO)                 !.......................
      IF (PRATIO.GE.(1.0D0)) PRATIO=0.5D0 ! RESET OFF SCALE VALUES
      ALFU=DABS(ALFU)                     !.......................
      IF (ALFU.LE.ZERO ) ALFU=0.90D0      ! Checks relaxation factors
      IF (ALFU.GT.1.0D0) ALFU=0.90D0      ! for momentum, pressure and
      ALFP=DABS(ALFP)                     ! energy equations.
      IF (ALFP.LE.ZERO   ) ALFP=0.80D0    !
      IF (ALFP.GT.1.0D0) ALFP=0.50D0      !
      ALFT=DABS(ALFT)                     !
      IF (ALFT.LE.(0.0D0)) ALFT=0.80D0    !
      IF (ALFT.GT.(1.0D0)) ALFT=0.80D0    !
      MPEPS=DABS(MPEPS)                   ! Checks if convergence criteria
      IF (MPEPS.EQ.ZERO   ) MPEPS=0.02D0  ! are low enough
      SFLOW=DABS(SFLOW)                   !
      IF (SFLOW.EQ.ZERO  )  SFLOW=0.006D0 !
      CD=DABS(CD)                         ! orifice discharge coefficient
      IF (CD.LE.Zero ) CD=1.0D0           !

c    !....................................!...............................2/6/95
      IF (BEARING.EQ.2) THEN              ! Check Pressure convergence
       IF (SFLOW.GT.5D-4) SFLOW=5D-4      !    PEPS=SFLOW for SEAL
      ELSE IF (BEARING.EQ.3) THEN         !.....
       IF (SFLOW.GT.2D-4) SFLOW=2D-4      !    PEPS=SFLOW for hydrodyn bearing
      ELSE                                !
       IF (SFLOW.LT.3D-3) SFLOW=3D-3      ! Check flow convergence for
      END IF                              !    HJB.
c    !....................................!...............................2/6/95

c    !....................................!....................................
      RUGR=DABS(RUGR)                     ! Check for maximum roughness params.
      RUGS=DABS(RUGS)                     !
      IF ((RUGR.GT.(0.10D0)).OR.(RUGS.GT.(0.15D0))) THEN
          WRITE (6, 300) RUGR, RUGS
          RETURN
      END IF
C ......................................!....................................
      IF (BEARING.EQ.1) THEN            !  FOR HJBS:
      IF (MOD(NPC, 2).EQ.1) GOTO 100    !  NPC MUST BE odd
        WRITE (6, 200) NPC              !
        CALL BEEPER                     !  JMOD gives remainder of NPC/2
        NPC=5                           !  If = 1, NPC is ODD
        CALL PAUSE                      !  If = 0, NPC is EVEN (not good)
      ELSE                              !.....
        HREC=ZERO                       ! for SEALS/JBs, no recess depth
        VSUP=ZERO                       ! or orifice supply volume.
        MODEL=2                         ! no double symmetry
      END IF                            !
c    !..................................!....................................
                                        ! CALCULATE DIMENSIONLESS PARAMETERS
 100  RADIUS=DIAM/2.0D0                 ! Journal radius
      ABAR=AR/LENGTH                    ! Dimensionless axial length of recess
      LD=LENGTH/DIAM                    ! L/D ratio
C ...........................................................................
      DIFF=DABS((LENGTHR/LENGTH)-0.50D0)

C=======
C NOTES:
c     ISYM=1 bearing is symmetric so program solves only 1/2 the bearing
c     ISYM=0 bearing is asymmetric and solution covers entire bearing.
c
c     SYMMETRY (ISYM=1) of bearing requires:
c     geometrical symmetry about bearing midplane (Z=0):
c       recess disposition, LR=LL, symmetric clearance distribution C(-Z)=C(Z)
c     symmetry in exit (discharge) pressures and end seals:
c       P(x,-LL) = P(x,LR),  CsealL = CsealR
c     symmetry in edge loss coefficients for hydrostatic bearings
c       XSI(left) = XSI(right) in axial direction
c     symmetry in journal positioning:
c       no journal misalignment AXO=AYO=0
c=======
c
c     DETERMINE TYP(ical) CLEARANCE FOR BEARING. THIS VALUE WILL BE USED
c     in ALL CALCULATIONS for definition of DIMENSIONLESS ECCENTRICITIES
c     TYP CLEARANCE = MIN VALUE OF MACHINED CLEARANCE
c
c    !..................................!
      IF (ICSTEP.eq.1) THEN             ! FOR STEP CLEARANCE BEARING
c    !..................................!
         CLEAR=DMIN1(ClearO,ClearR)     ! ........
         CINLET=ClearO
         IF ( (DIFF.LT.0.01D0) .AND.(YR.EQ.-YL).
     +        AND.(ClearL.eq.ClearR)) THEN        !==> bearing is symmetric
            ISYM=1
            CEXIT=ClearR
         ELSE
            ISYM=0
            CLEAR=DMIN1(Clear,ClearL)             !==> bearing is asymmetric
            CEXIT=0.5D0*(ClearL+ClearR)
         END IF
c    !..................................!
      ELSE                              ! CONTINUOUS CLEARANCE BEARING
c    !..................................!
         ISYM=1                         !.............................
         CALL SPLINE(NJ,Z,Cl,BCL,CCL,DCL)  !Calculate spline coefficients
         CALL SEVAL(Nj,0.0D0,Z,CL,BCL,CCL,DCL,CINLET)
         CALL SEVAL(Nj,LengthR,Z,CL,BCL,CCL,DCL,CEXIT)
         CLEAR=CL(1)                       ! Calculates MIN clearance
         Zdumy=0			   !
         DO j=1, NJ                        !
          CLEAR=DMIN1(Cl(j),Clear)         ! = TYP Clearance
          IF (Z(j).LT.0.0D0) Zdumy=-1      !
         END DO
         IF (DIFF.LT.0.01D0) THEN
            ISYM=1
            IF (Zdumy.EQ.-1) ISYM=0
         ELSE
            ISYM=0
         END IF
c    !..................................!
      END IF
c    !..................................!


c ......................................!...............................
       IF (CLEAR.LE.Zero   ) THEN       ! CLEAR= Characteristic Clearance
          WRITE (6, 301) CINLET,CEXIT   !
          STOP                          !
       END IF                           !
c ......................................!...............................
      IF (PRELOAD.GT.CLEAR) THEN        ! journal is in contact with bearing
          WRITE (6, 3011) PRELOAD,CLEAR !
          STOP                          !
      END IF                            !
c
c IF Preload <>0 then actual assembled clearance is C=CLEAR-PRELOAD.
c
c ......................................!...............................
       HCDIM=HCELL/CLEAR                ! TYP cell depth/C*    6/20/95
c ......................................!...............................
c                                       !
      CRATIO=CLEAR/RADIUS               ! Clearance ratio=c/R
      C2=CLEAR**2                       !
      CR2=CRATIO*CRATIO                 !
      OMEGA=RPM*PI/30.0D0               ! Rotating angular frequency (rad/sec)
c ......................................!..................................
c SET Characteristic discharge pressure PA [N/m2]

      IF (BEARING.EQ.2) THEN            ! TYP discharge pressure
        PA=Pright                       ! for seals
        ISYM=0                          !
      ELSE                              !........... &
        PA= DMIN1(Pleft, Pright)        ! for HJBs and hydrodynamic bearings
      END IF                            !
c ......................................!..................................
      IF (TILTPAD.EQ.1) THEN            ! for tilt-pads:
         AXO=zero                       ! SET no journal misalignment
         AYO=zero                       !
         ZO=zero                        !
      END IF                            !
c ......................................!..................................


c ......................................! Check Symmetry in Pressure BC's
      IF (ISYM.EQ.0) GOTO 11
c ......................................!..................................
      IF ( (DABS(Pright-Pleft).LE.1.0D-3).AND.       ! check symmetry in
     +     (DABS(Cleft-Cright).LE.1.0D-6)   ) THEN   ! exit pressures &
          ISYM=1                                     ! end seals
      ELSE                                           !.....................!
          ISYM=0
          GOTO 11
      END IF
c ......................................! Check Misalignment Condition
      IF ((AXO.EQ.zero).AND.(AYO.EQ.zero)) THEN
          ISYM=1
      ELSE
          ISYM=0
          GOTO 11
      END IF
c ......................................! Check SYMMETRY in Non-uniform Exit Press
      IF ((IPRuni.NE.1).AND.(IPLuni.NE.1)) THEN
        IF (NPRcs.NE.NPLcs) ISYM=0
        DO j=1, NPRcs
         IF ((DABS(PRCOEFC(j)-PLCOEFC(j)).GE.1.0D-3)) ISYM=0
         IF ((DABS(PRCOEFS(j)-PLCOEFS(j)).GE.1.0D-3)) ISYM=0
        END DO
      ELSE
        ISYM=1
      END IF
c ......................................!
 11   CONTINUE
c ......................................!

      IF ((MODEL.EQ.1).AND.(BEARING.EQ.1)) THEN
C                                       ! MODEL=1 DOUBLE ROW HJB
           ISYM=0                       ! then -> ASYMMETRIC BEARING
           Pright=Pleft                 !
           Cright=zero                  ! w/no end seals
           AXO=zero                     ! w/no misalignment
           AYO=zero                     !
           ZO=zero                      !
           LENGTHL=LENGTH-LENGTHR       !
      END IF                            !....................

c ......................................!...................................

c ......................................!...................................
c Calculate props at (S) supply, (le) left side LL, (ri) right side LR
c           and at (a) TYP Pa
c............................................................................

      Press=Ps/1.0D+06                  ! Pressure in MPa
      Presa=Pa/1.0D+06
c    !.......................!
      IF (IF.EQ.5 ) GOTO 1101           ! Air
      IF (IF.EQ.6 ) GOTO 1102           ! Water
      IF (IF.EQ.7 ) GOTO 1103           ! Oil-->Properties can be modified
      IF (IF.EQ.8 ) GOTO 1105           ! Bilinear Properties with P & T
      IF (IF.EQ.12) GOTO 1104           ! Linear Properties with pressure
c    !.......................!
      CALL CHECKPROPS        ! => for cryogenic liquids IF=1,2,3,4
c    !.......................!
      GOTO 1114

c................................................................
C --- Air as an Ideal Gas
c................................................................
 1101 CPS=CP                          ! Specific Heat (J/kg.K)
      Vsound=DSQRT(GAMA*R*TEMPK)      ! Sound of speed (m/s)
      RHOS=PS/R/TEMPK                 ! Density (supply) (Kg/m^3)
      EMUS=1.79D-5*(TEMPK/288.2D0)**.76!Viscosity ("")  (N-s/m^2)
      RHOA=PA/R/TEMPK                 ! Guessed characteristic properties
      RHOSU=RHOA                      ! Sump density
      EMUA=EMUS                       ! temperature is a constant
      RHOle=Pleft/R/TEMPK             ! Left discharge
      EMUle=EMUA                      !
      RHOri=Pright/R/TEMPK            ! Right discharge
      EMUri=EMUA                      !
      THS=(9.D0*GAMA-5.D0)/4.D0       !
     +                 /GAMA*EMUS*CPS ! Thermal conductivity (W/m.K)
      BETAKS=1.D0/TEMPK               ! Volumetric expansion coeff.(1/K)
      BETA=1.D0/PS                    ! Compressibility factor (m^2/N)
      GOTO 1114

c................................................................
C --- FOR WATER
c................................................................
 1102 RHOS=1000.D0*DEXP(-4.85D-4*     ! Density (supply) (Kg/m^3)
     + ((TEMPK-293.D0)-(Press-0.1D0)))!
      EMUS=1.005D-3*(TEMPK/293.D0)    ! Viscosity (`')  (N-s/m^2)
     +  **8.9*DEXP(4700.D0*(1.D0/TEMPK! Formulae from Sherman,F.S.
     +   -1.D0/293.D0))               ! <<Viscous Flows>>, (1990)
      RHOA=1000.D0*DEXP(-4.85D-4*     ! Guessed exit properties
     +((TEMPK-293.D0)-(Presa-0.1D0))) !
      RHOSU=RHOA                      ! Sump density
      EMUA=EMUS                       ! temperature is a constant
      RHOle=RHOA                      ! Left discharge (approx.)
      EMUle=EMUA                      !
      RHOri=RHOA                      ! Right discharge (approx.)
      EMUri=EMUA                      !
      CPS=4182.D0                     ! Specific heat (J/kg.K)
      BETAKS=2.07167D-4               ! Volumetric expansion coeff.(1/K)
      THS=0.597D0                     ! Thermal conductivity (W/m.K)
      BETA=4.6412D-10                 ! Compressibility factor (m^2/N)
      Vsound=1.0D0/DSQRT(Beta*Rhos)   ! Sound speed (m/s)

      GOTO 1114
c................................................................
C --- FOR OIL (From  Pinkus(1990),<<Thermal Effects ...>>)
c................................................................
 1103 RHOA=RHOS                       ! Guessed exit properties--
      EMUA=EMUS                       ! temperature is a constant
      RHOSU=RHOA                      ! Sump density
      RHOle=RHOA                      ! Left discharge (approx.)
      EMUle=EMUA                      !
      RHOri=RHOA                      ! Right discharge (approx.)
      EMUri=EMUA                      ! OIL is regarded as incompressible

      GOTO 1114
c................................................................
C --- FOR Bi-linear properties (with P&T) fluids--
c................................................................
C#### CARE here since TSUMP has not been defined yet ######
 1105 RHOS=134.7-2.08D-6*PS-2.714*TEMPK   ! X=a+b*P+c*T+d*P*T
     +    +0.12D-6*PS*TEMPK               ! Coeffs. a to c are based
      EMUS=0.79D-5+0.957D-12*PS-0.1297D-6 ! on  LH2 properties at
     +    *TEMPK-1.034D-14*PS*TEMPK       ! P=2.413MPa and 16.18MPa
      CPS=147218.0-8791.0D-6*PS-2747.0    ! and T=37 and 46.7 K-deg.
     +    *TEMPK+180.0D-6*PS*TEMPK        ! The bi-linear properties
      RHOA=134.7-2.08D-6*PA-2.714*TSUMP   ! model serves as a quick
     +    +0.12D-6*PA*TSUMP               ! check of THD calculation
      EMUA=0.79D-5+0.957D-12*PA-0.1297D-6 ! without involvement of
     +    *TSUMP-1.034D-14*PA*TSUMP       ! real cryogenic liquids.
      RHOSU=RHOA                          !
      BETAKS=(2.714-0.12D-6*PS)/RHOS      !
      THS=0.0917                          !
      BETA=(-2.08+0.12*TEMPK)*1.0D-6/RHOS ! m2/N

      GOTO 1114
c................................................................
c     for unknown fluid (IF=12):
c................................................................
 1104 CONTINUE                          ! calculate props for fluid
c     !.................................! with properties varying linearly with
       IF (P1prop.EQ.P2prop) THEN       ! pressure
c     !.................................!
       RHOS=RHO2
       RHOle=RHOS
       RHOri=RHOS
       RHOa=RHOS
       EMUs=EMU2
       EMUle=EMUs
       EMUri=EMUs
       EMUa=EMUs

c     !.................................!^^^ incompressible fluid
       ELSE
c     !.................................!vvv X(P) = X2 + (X1-X2) (P-P2)/(P1-P2)
       RHOs=RHO2+(RHO1-RHO2)*(Ps-P2prop)/(P1prop-P2prop)
       RHOle=RHO2+(RHO1-RHO2)*(Pleft-P2prop)/(P1prop-P2prop)
       RHOri=RHO2+(RHO1-RHO2)*(Pright-P2prop)/(P1prop-P2prop)
       RHOa=RHO2+(RHO1-RHO2)*(Pa-P2prop)/(P1prop-P2prop)
       EMUs=EMU2+(EMU1-EMU2)*(Ps-P2prop)/(P1prop-P2prop)
       EMUle=EMU2+(EMU1-EMU2)*(Pleft-P2prop)/(P1prop-P2prop)
       EMUri=EMU2+(EMU1-EMU2)*(Pright-P2prop)/(P1prop-P2prop)
       EMUa=EMU2+(EMU1-EMU2)*(Pa-P2prop)/(P1prop-P2prop)
c     !.................................!
      END IF
c     !.................................! barotropic (linear) fluid

      Vsound=0.0D0                      ! Sound  speed (m/s)
      Beta=DABS(Beta)                   ! compressibility parameter (m2/N)
      RHOSU=RHOA                        ! Sump density

      ISOTH=1                           ! Automatic isothermal fluid model


c ......................................!...................................
c     SET CHARACTERISTIC FLUID PROPERTIES:  RHOTYP AND EMUTYP
c ..........................................................................
 1114 IF (BEARING.LE.2) THEN            ! FOR HJBS & SEALS
        RHO=DABS(RHOS)                  ! FLUID properties at
        EMU=DABS(EMUS)                  ! supply conditions
      ELSE                              !
         INERP=0                        ! no inertia entrance effects
         ALPHA=0.50                     ! no swirl effects
         PS=PA                          ! suppy=sump pressure
        RHO=DABS(RHOA)                  ! fluid props at
        EMU=DABS(EMUA)                  ! exit pressure
      END IF                            !
c ......................................!...................................

c    !..................................!....................
      IF (RHO.EQ.Zero   ) THEN          ! CHECK FOR NONZERO PROPERTIES
         WRITE (6, 302) RHO             !
         STOP                           !
      END IF                            !
      IF (RHOA.LE.Zero   ) THEN         !
         WRITE (6, 303) RHOA
         STOP
      END IF                            !
      IF (EMU.EQ.Zero   ) THEN          !
         WRITE (6, 304) EMU             !
         WRITE (6, 304) EMU             !
         STOP
      END IF                            !
      IF (EMUA.LE.Zero   ) THEN         !
         WRITE (6, 305) EMUA            !
         STOP                           !
      END IF                            !
c    !..................................!....................
      IF ((MODEL.EQ.1).AND.(BEARING.EQ.1)) THEN ! MODEL=1 DOUBLE ROW HJB
          RHOle=RHOri                   !
          EMUle=EMUri                   !
      END IF                            !
                                        !
C.......................................!.....................................!
      RHOTYP=RHO                        ! Typical props used
      EMUTYP=EMU                        ! for dimensionless values
c.......................................!.....................................!

c.......................................!...................................
c                                       ! Characteristic Pressure DELTAP=PSA
      PSA=PS-PA                         !
      PATYP=PA                          !
c    !.......................!          ! dimensionless p = (P-PA)/PSA
      IF (PSA.EQ.ZERO) THEN             !
c    !.......................!          !
         PSA=OMEGA*EMU/CR2              ! EMUx OMEGA x (R/C)^2
         INERP=0                        ! no entrance edge Pdrop

         IF (IF.EQ.5) THEN  !---> for air
              PSA=PA           !     i.e. p=P/PA
              PA=0.0D0         !     use exit pressure for dim.pressure
              IF (PSA.EQ.0.0D0) THEN
                WRITE (6,*)
     +        ' $ AIR: TYP pressure PA can not be 0.0D0,  => STOP'
              STOP
            END IF           !
          END IF             !...........!......
c    !.......................!          !
      END IF                            !SET DIM Pressure
c    !.......................!          !

c ......................................!...................................
      IF ((BEARING.GE.2).OR.            ! FOR SEAL/PLAIN BEARING
     +    (NPOCKET.EQ.0)) THEN          !
         PRATIO=(PS-PA)/PSA             !
C##      IF (IF.EQ.5) PRATIO=0.0D0      !so as not to bring guespp.f into an
      END IF                            !error
c ......................................!...................................
c                                       ! dimensionless exit presures
      Ple=(Pleft-PA)/PSA                ! left side
      Pri=(Pright-PA)/PSA               ! right side
      PCAV=(PC-PA)/PSA                  ! cavitation pressure
      Pcavi=PCAV                        !
      PBACK=(PATYP-PA)/PSA              ! typ back pressure for COMPLIANCE bearing
      PATYP=PA                          ! or air bearing.

      TSUMP  =TEMPK                     ! Initial sump temperature
      IF (TSHAFT.EQ.0.0D0) TSHAFT=TEMPK   ! = TSUPPLY
      IF (TSTATOR.EQ.0.0D0) TSTATOR=TEMPK ! = TSUPPLY

c.......................................!....................................
      IF (IF.EQ.12) THEN                ! for barotropic fluid
c    !...................!              ! set coefficients for linear property .
      DEN2=DABS(RHO2)/RHO               ! eqn X(P) = X2 + (X1-X2) (P-P2)/(P1-P2)
      VIS2=DABS(EMU2)/EMU               ! DIM. density & viscosity at P2
      P1=(P1prop-PA)/PSA                !
      P2=(P2prop-PA)/PSA                !
      IF (P1prop.EQ.P2prop) THEN        !=> incompressible
         DEN12P12=0.0D0                 !
         VIS12P12=0.0D0                 !
      ELSE                              !=> barotropic fluid
         DEN12P12=(RHO1-RHO2)/(P1-P2)/RHO
         VIS12P12=(EMU1-EMU2)/(P1-P2)/EMU
      END IF
c                          .............! Calculate Sound speed [m/s]
      IF (BETA.GT.ZERO ) THEN
         Rhor=(DABS(RHOS)-RHOA)*Pratio+RHOA
         Vsound=1.0D0/DSQRT(BETA*Rhor)
      END IF
c    !.....................!            !...........................!
      ELSE
c    !.....................!            !...........................!
      P1PROP=PS                         !
      P2PROP=PA                         !
      RHO1=RHOS                         !
      RHO2=RHOA                         !
      EMU1=EMUS                         !
      EMU2=EMUA                         !
c    !...................!              !..... FLUID <> 12
      END IF
c    !...................!

c................................................................
      RHOSd=RHOSU/RHO                   ! Dimensionless Sump density
      L4=PSA*Beta                       ! compressibility ratio (dimensionless)
                                        !..................................
c     U*=PSAxC^2/(EMUxR)                ! TYP fluid velocity
      SPEED=OMEGA*EMU/(PSA*CR2)         ! Hydrodynamic Speed parameter
      MSPEED=SPEED/2.0D0                ! Mean speed parameter
      SWIRL=ALPHA*SPEED                 ! Entrance swirl speed parameter
      VIS=EMU/RHO                       ! kinematic viscosity
C                                       ! ...................................
      U=OMEGA*RADIUS                    ! rotor surface speed
C                                       ! ...................................
      RENC=U*CLEAR/VIS                  ! Circ. Reynolds number at rotor speed
      REP=C2*CRATIO*PSA/(EMU*VIS)       ! Poiseuille Reynolds number
      REY=REP*CRATIO                    ! Flow nominal Reynolds number
C                                       !
      BETU=1.0D0-ALFU                   ! Relaxation parameter
C                                       !
C#    ASPEC=BR/DIAM                     ! REC Circ. aspect ratio
      HRECD=HREC/CLEAR                  ! dim. recess depth

      QFACTOR=C2*CLEAR*PSA/VIS          ! Flow Conversion Coefficient (Kg/s)
      FFACTOR=PSA*DIAM*DIAM/4.0D0       ! Force Conversion coefficient (N)
      MFACTOR=FFACTOR*DIAM/2.0D0        ! Moment Conversion coefficient (Nm)
      TO=PSA*CLEAR*DIAM*DIAM/4.0D0      ! Torque Conversion coefficient (Nm)
C                                       !.....................................!
      UC=C2*PSA/EMU/RADIUS              ! C#--CHARACTERISTIC Velocity
      TC=TEMPK                          ! C#--CHARACTERISTIC temperature = TSUPPLY
      EC=UC*UC/TEMPK/CPS                ! Eckert No.
      TSUMPd =TEMPK/TC                  ! Dimensionless sump temperature =1
      TBJ=TSHAFT/TC                     ! dimensionless shaft & stator
      TBS=TSTATOR/TC                    ! temperatures TC=TSUPPLY
c    !.....................!            !.....................................!
      IF (ISOTH.GE.4) THEN 		!
c    !.....................!            ! For BEARING RADIAL HEAT FLOW
c##      IF (ROUTER.LT.RADIUS) ROUTER=RADIUS
         TBS=TBOUT/TC                   ! dimensionless outer temperature
         HKB=(THERMALK/RADIUS)/         ! dimensionless heat transfer coeff.
     +       (REY*CPS*EMU/CLEAR)/       ! due to solid conductivity
     +       DLOG(ROUTER/RADIUS)        !
c    !.....................!            !
      END IF				! ==: 6/28/94 LSA
c    !.....................!            !.....................................!



C                                       !.....................................!


c    !==================================!////////////////////////////////////
c     FOR Isothermal model of LH2 (IF=1) program sets always groove temperatures
c         equal to Tsupply or Tsump so no iteration is required on
c         energy transfer at grooves between pads.

      IF (ISOTH.EQ.1.OR.IF.EQ.1) THEN   ! Note: for LH2, ITPAD=0
         ITPAD=0                        ! Ggroove temperature=Tsupply
      ELSE IF (IFULL.EQ.0) THEN         !
         ITPAD=1                        ! Ggroove temperature will be updated
      END IF                            ! for pad bearings.

c    !==================================!////////////////////////////////////
c
      CJET=ZERO                         !...................................!
      DORIF=ZERO                        !
      CORIFS=ZERO                       ! SET ORIFICE PARAMETER
      IF (NPOCKET.GT.0) THEN            !
        DO J=1, NPOCKET                 !
         Ao=Cd*PI*(DIAORIF(J)*DIAORIF(J))/4.0D0 ! Equiv. orifice area (m2)
         DENom=C2*CLEAR*DSQRT(RHOS*PSA/2.0D0)
         CORIF(J)=Ao*EMUS/Denom         !    => DIM Orifice coefficient
         DORIF=DORIF+DIAORIF(J)         !
         CORIFS=CORIFS+CORIF(J)         !
        END DO                          !
        DORIF=DORIF/NPOCKET             ! average orifice diameter &
        CORIFS=CORIFS/NPOCKET           ! orifice coefficient
        CJET=PI*CD*DORIF**2*DSIN(ANGLEJ*PI/180.0D0)/
     +     (2.0D0*AR*(CLEAR+HREC))      ! Angled injection coefficient
      END IF                            !
c    !.....................!            !.....................................!

c    !.....................!            !.....................................!
C                                       ! SET ENTRANCE LOSS COEFFICIENTS
      KLOSYl=REY*(1.0D0+LOSXSIyl)/2.0D0 ! Loss factor*Rey left recess side
      KLOSYr=REY*(1.0D0+LOSXSIyr)/2.0D0 !  ''    ''       right recess ''
      KLOSXu=REY*(1.0D0+LOSXSIxu)/2.0D0 ! Loss factor Circumf. upstream recess
      KLOSXd=REY*(1.0D0+LOSXSIxd)/2.0D0 ! Loss factor Circ. downstream recess
      LOSXSI=LOSXSIyr                   !....
      KLOSPAD=0.5D0*REY*LOSleadP*MSPEED**2 ! Leading edge Pad RAM Press factor
      HRECU=HREC                        !
c ......................................! DIMensionless END Seals coeffs
      Csel=REY*Cleft/2.0D0              !
      Cser=REY*Cright/2.0D0             !
c................................................................
c                                       !Check value of Moody's coefficients
      Amod=DABS(Amod)
      Bmod=DABS(Bmod)
      Expo=DABS(Expo)
      IF (Amod.eq.Zero)  Amod=0.001375D0
      IF (Bmod.eq.Zero)  Bmod=5.0D+5
      IF (Expo.eq.Zero)  Expo=1.0D0/3.00D0

C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IF ((ISYMold.eq.1).AND.(ISYM.eq.0)) THEN   ! Set UVP mirror image to
	CALL TEMPASYM                            ! start solution for
      END IF                                     ! asymmetric HJB

C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 300     FORMAT (' ', 'PRG STOPS: Relative  Roughness coefs > 0.10',
     +            /, 'ROUGHNESS VALUES UNREALISTIC in ANY application,',
     +            /, ' RUGR: ', E12.5E2, '   RUGS: ', E12.5E2)

 301     FORMAT (' ', ' PRG STOPS: Clearance(Inlet/Exit) <= Zero',
     +           /, 'Cinlet:', E12.5E2, 'm, Cexit:', E12.5E2, 'm')

3011     FORMAT (' ', ' PRG STOPS: PAD PRELOAD',E12.5E2, 'm, is',
     +           ' larger than MIN Clearance:', E12.5E2, 'm')

 302     FORMAT (' ', ' PRG STOPS: Inlet Supply Density (Rhos):',
     +            E12.5E2, 'Kg/m3 IS LESS THAN ZERO')

 303     FORMAT (' ', ' PRG STOPS: Discharge Density (Rhoa):',
     +            E12.5E2, 'Kg/m3 IS LESS OR EQUAL to ZERO')

 304     FORMAT (' ', 'PRG STOPS: Inlet Supply Viscosity (Emus):',
     +            E12.5E2, 'Nsec/m2 is LESS or EQUAL to ZERO')

 305     FORMAT (' ', 'PRG STOPS: Discharge Viscosity (Emua):',
     +            E12.5E2, 'Nsec/m2 is LESS or EQUAL to ZERO')

 177     FORMAT (' ', 'PRG STOPS: Eccentricity >=1.0, Exo=',
     + F7.5, 2X, 'Eyo=', F7.5, 2X, 'Ecc=', F7.5)

 200  FORMAT (' ', 'Number of points on pocket circumferential MUST BE',
     +        /, ' ODD and not equal to ', I3, /, ' PROGRAM sets a',
     +           ' value NPC=5')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Askiter                                                     **
C **                                                                         **
C **  ASKITER:  Ask if changes in number of iterations are to be made.       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ASKITER

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*1 YN
      INTEGER I, IOS, ITER, ITMAX, ITPMAX

C ----------------------------------------------------------------------------
C --  ASKITER code                                                          --
C ----------------------------------------------------------------------------
C Unit 42 is used as a dumping place for translations between formats.
      OPEN (UNIT=42, STATUS='SCRATCH', IOSTAT=IOS,
     +      ERR=1234)                       !

C...........................................!....................
   10 CALL BEEPER                           !
      WRITE (6, *) 'USER may change any of these parameters.'
      WRITE (6, *) 'Press <Enter> to accept the defaults'
      WRITE (6, *)

C:::::::::::::::::::::::::::::::::::::::::::!
      DO WHILE ((YN.NE.'Y').AND.(YN.NE.'y').AND.(YN.NE.' ')) !
C:::::::::::::::::::::::::::::::::::::::::::!

          YN='Y'
c      !.....................!
          CALL INPUTEXEY
c      !.....................!

          WRITE (6, 100)                    !
  100     FORMAT ('$', 'Accept changes (Y/N) (Default: Y) ? ')
          READ (5, 110) YN                  !
  110     FORMAT (1A)                       !
C
C:::::::::::::::::::::::::::::::::::::::::::!
      END DO                                ! END DO WHILE
C:::::::::::::::::::::::::::::::::::::::::::!
                                            !
      CLOSE (UNIT=42)                       !
C...........................................!................

C##      write (6,*) '$    CALCULATIONS START'
C##      write (6,*) '$    ------------ -----'

      ITMAX=MAX0(ITMAX,99)
      ITPMAX=MAX0(ITPMAX,10)
      RETURN


C...........................................!................
                                            !
 1234 CALL DECODIOS(IOS)                    !
      WRITE (6, *) 'I/O ERROR ON INPUT DATA (unit 42)=> STOP'
      STOP                                  !

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Checkprops                                                  **
C **                                                                         **
C **  CHECKPROPS: checks that boundary pressures are within of fluid props.  **
C **              for cryogenic liquids                                      **
C *****************************************************************************

      SUBROUTINE CHECKPROPS

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM3/ EMU, RHO, PC, CD, DORIF, LOSXSI,ALPHA!, RPM, PS, PA
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /COMPRE/ L4, PCAVI
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /CONT/ IFF
C........................................................ VARIABLE DECLARATIONS
      DOUBLE PRECISION EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 L4, PCAVI, HRECU, VSUP, BETA,
     +                 PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 TEMPK, VSOUND, CPS,THS,BETAKS

      INTEGER IF, IL, IFF

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION Precess,Rhor, EMur, D, FINDD, DEltp,RHorp,
     +                 THA,CPA,THR,CPR ,DELT,DNEW

C ----------------------------------------------------------------------------
C --  CHECKPROPScode                                                   --
C ----------------------------------------------------------------------------


C ----------------------------------------------------------------------------
C IF: cryogenic liquid: IF=1 LH2, file PH2.COF
C                          2  N2,   ''  N2.COF
C                          3  O2,   ''  O2.COF
C                          4 methane,   METH.COF
C ALL subs and function called here are in program miprops.f
C
C ----------------------------------------------------------------------------

      CALL LIMITS(Ps/1.0D+06,Tempk,IL)
      IF (IL.le.0) THEN
        write (6,*) 'ERROR: Ps & T at supply OUT OF RANGE: ABORT'
        GOTO 104
      END IF

      CALL LIMITS(Pleft/1.0D+06,Tempk,IL)
      IF (IL.LE.0) THEN
        write (6,*) 'ERROR: PLeft & T at exit OUT OF RANGE: ABORT'
        GOTO 104
      END IF

      CALL LIMITS(Pright/1.0D+06,Tempk,IL)
      IF (IL.LE.0) THEN
        write (6,*) 'ERROR: PRight & T at exit OUT OF RANGE: ABORT'
        GOTO 104
      END IF

      CALL LIMITS(Pc/1.0D+06,Tempk,IL)
      IF (IL.LE.0) THEN
        write (6,*) 'ERROR: Pcavitation & T OUT OF RANGE: ABORT'
        GOTO 104
      END IF

c...............................................................
c    Calculate Density & viscosity at supply conditions
c...............................................................
      D=FINDD(Ps/1.0D+06,Tempk)
      DELT=0.01D0                       !
      DNEW=FINDD(Ps/1.0D+06,Tempk+DELT) !
      BETAKS=-(DNEW-D)/DELT/D           ! Thermal expansion coeff.(1/K)
      CALL REPRO(Ps/1.0D+06,D,Tempk,RHOS,EMUS,Vsound,THS,CPS)
c...............................................................
c    Calculate density & viscosity at discharge pressures
c...............................................................
      D=FINDD(Pleft/1.0D+06,Tempk)
      CALL REPRO(Pleft/1.0D+06,D,Tempk,RHOle,EMUle,Vsound,THA,CPA)

      D=FINDD(Pright/1.0D+06,Tempk)
      CALL REPRO(Pright/1.0D+06,D,Tempk,RHOri,EMUri,Vsound,THA,CPA)
c...............................................................
c     Calculate density & viscosity at charact. pressure Pa
      D=FINDD(Pa/1.0D+06,Tempk)
      CALL REPRO(Pa/1.0D+06,D,Tempk,RHOa,EMUa,Vsound,THA,CPA)

c...............................................................
c    Calculate density & viscosity at nominal recess pressure
c...............................................................
      Precess=(PRATIO*(Ps-Pa)+Pa)/1.0D+06
      CALL LIMITS(Precess, Tempk, IL)
      IF (IL.le.0) goto 104
      D=FINDD(Precess,Tempk)
      CALL REPRO(Precess,D,Tempk,RHOr,EMUr,Vsound,THR,CPR)
      DEltp=0.01D0*Precess
      D=FINDD(Precess+Deltp,Tempk)
      CALL REPRO(Precess+Deltp,D,Tempk,Rhorp,Emur,Vsound,THR,CPR)

c...............................................................
c     Calculate compressibility Factor at Recess pressure
c...............................................................
      Beta=(Rhorp-Rhor)/Rhor/(1.0D+06*Deltp)

      RETURN
c     ............................................................

 104  WRITE (6,109)
 109  FORMAT (' ','PROGRAM ABORTS: FLuid properties out of range',
     +        /, 'REVISE your INPUT DATA')
      STOPc     ............................................................
      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Locprops                                                    **
C **                                                                         **
C **  LOCPROPS: calculates fluid properties at values of (P,T)               **
C **                                                                         **
C *****************************************************************************

c ............................................................................

      SUBROUTINE LOCPROPS(RHO,EMU,CPP,BETAK,TH, P,TBAR)
c                    --    ^   ^   ^    ^   ^   ^   ^   --
c                    --   All Dimensionless Variables   --
c.............................................................................

      IMPLICIT NONE

      COMMON /AIR/ CP, RR, AIRMU, GAMA, PR, CKS
      COMMON /compre/ L4,PCAV
      COMMON /CONT/ IFF
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /LIQUID/ TEMPK, Vsound, IF, IL
      COMMON /PARAM3/ EMUN,RHON,PC,CD,DORIF,LOSXSI,ALPHA!,RPM,PS,PA
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DEN2,VIS2,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2
      COMMON /THERMAL/ ALFT, UC, TC ,Ec
      COMMON /THERMID/ ISOTH
      COMMON /OILCOEF/ Talpha

      DOUBLE PRECISION CP, RR, AIRMU, GAMA, PR, CKS,
     +                 EMUN,RHON,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 RHOTYP,EMUTYP,DEN2,VIS2,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 CPS,THS,BETAKS ,ALFT, UC, TC, Ec ,Talpha
      DOUBLE PRECISION L4, PCAV
      DOUBLE PRECISION RHO,EMU,P,TBAR,TH
      DOUBLE PRECISION TEMPK, Vsound ,CPP,BETAK
      INTEGER IF, IL, IFF ,ISOTH

c........................................................
c Local variable declarations
c........................................................

      DOUBLE PRECISION PMPA,D,FINDD,RDUM,EDUM,DNEW,
     +                 PRESS,TK,PBAR

C ----------------------------------------------------------------------------

C   !......................................!
      PBAR =DMAX1(P,PCAV)                  ! =(P-Pa)/Psa: dimensionless p
C   !......................................!

c    !-------------------------------------!
      IF (IF.EQ.12) THEN                   ! IF =12 barotropic fluid
c    !-------------------------------------! Liquid with linear properties
c     DEN2=RHO2/RHOTYP, DEN12P12=(RHO1-RHO2) /(p1-p2)/RHOTYP
c     VIS2=EMU2/EMUTYP, VIS12P12=(EMU1-EMU2) /(p1-p2)/EMUTYP

c(1)  Linear property liquids              ! ISOTHERMAL & Properties change
c                                          ! with pressure linearly
      RHO=DEN2+DEN12P12*(PBAR-P2)          ! Dimensionless density
      EMU=VIS2+VIS12P12*(PBAR-P2)          ! Dimensionless viscosity
      BETAK=1.D0                           ! Note: BETAK,CPP,&TH
      CPP=1.D0                             ! will not be used in
      TH =1.D0                             ! linear property liquids since ISOTH=1
      RETURN                               !
c    !-------------------------------------!
      END IF
c    !-------------------------------------!

      PRESS=PBAR*PSA+PA                    ! IN (N/M^2)
      PMPA=PRESS/1.0D+6                    ! Pressure in MegaPascals
      TK=TBAR*TC                           ! Temperature (K)


c    !-------------------------------------!
      IF (IF.eq.5) THEN  		   ! IF=5
c    !-------------------------------------!
c(2)  Air properties                       ! AIR
c					   !
      RHO=PRESS/RR/TK/RHOTYP               ! Dimensionless Dens.
      EMU=1.79D-5*(TK/288.2d0)**0.76/EMUTYP! Dimensionless vis.(Constantinescu)
      CPP=1.D0                             ! Dimensionless Specific Heat Coeff.
      BETAK=1.D0/TBAR                      ! Dimensionless Expansion Coefficient
      TH=(9.D0*GAMA-5.D0)/4.D0/GAMA*EMU*CPP! Dimensionless Conductivity
      Vsound=DSQRT(GAMA*RR*TK)             ! (From Dr. Sherman' Fluid Mechanics)
      RETURN                               !
c    !-------------------------------------!
      ELSE IF (IF.EQ.6) THEN               ! IF=6
c    !-------------------------------------!
c(3)  For water                            ! WATER
c					   !
        RHO=1000.D0*DEXP(-4.85D-4*         !
     +     ((TK-293.D0)-(PMPA-0.1D0)))     !
        RHO=RHO/RHOTYP                     ! Formulae from Sherman,F.S.(1990)
        EMU=1.005D-3*(TK/293.D0)**8.9*     ! <<Viscous Flows>>
     +  DEXP(4700.D0*(1.D0/TK-1.D0/293.D0))!
        EMU=EMU/EMUTYP                     !
	CPP=1.D0                           !
	BETAK=BETAKS*TC                    !
	TH=1.D0                            !
	RETURN                             !
c    !-------------------------------------!
      ELSE IF (IF.EQ.7) THEN               ! IF=7
c    !-------------------------------------!
c(4)  For Oil                              ! OIL
c					   ! is taken as incompressible
	RHO=1.D0                           ! Formulae from Pinkus(1991)
	EMU=DEXP(-Talpha*(TK-Tempk))       !
	CPP=1.D0                           ! MU = MUx EXP(a(T-Tx))
	BETAK=BETAKS*TC                    !
	TH=1.D0                            !
	RETURN                             !
c    !-------------------------------------!
      ELSE IF (IF.EQ.8) THEN               ! IF=8
c    !-------------------------------------! BILINEAR PROPS FLUID
c(4)  For Bilinear Properties              ! Dimensionless Properties
c                                          ! X=a+b*P+c*T+d*P*T
      RHO=(134.7-2.08*PMPA-2.714*TK        ! Coeffs. a to c are based
     +    +0.12*PMPA*TK)/RHOS              ! on  LH2 properties at
      EMU=(0.79D-5+0.957D-6*PMPA-0.1297D-6 ! P=2.413MPa and 16.18MPa
     +    *TK-1.034D-8*PMPA*TK)/EMUS       ! and T=37 and 46.7 K-deg.
      CPP=(147218.0-8791.0*PMPA-2747.0     !
     +    *TK+180.0*PMPA*TK)/CPS           ! The bi-linear properties
      BETAK=(2.714-0.12*PMPA)/RHO/RHOS*TC  ! model serves as a quick
      TH=1.0D0                             ! check of THD calculation
c                                          ! without involvement of
      RETURN                               ! real cryogenic liquids.
c    !-------------------------------------!
      ELSE
c    !-------------------------------------!
C(5)  For cryogenic liquids                !
c     CRYOGENIC liquid Use eqn. of state from prg: miprops
c     RDUM: density in Kg/m3   EDUM: viscosity in Pa-s
c
c NOTE: all subs called below are on on mipropst.f

c    !.....................................!
          CALL Limits(PMPA,TK,IL)          ! check if within range
c    !.....................................!
          IF (IL.LE.0) THEN                !
            WRITE (6, 99)                  !
            STOP                           !
          END IF                           !
c    !.....................................!

          D=FINDD(PMPA,TK)                 ! Find density a fn of P and T
	  DNEW=FINDD(PMPA,TK+0.001D0)
	  BETAK=-(DNEW-D)*TC/0.001D0/D     ! Non-D Thermal Expansion Coeff.
          CALL REPRO(PMPA,D,TK,RDUM,EDUM,Vsound,TH,CPP)
          CPP=CPP/CPS                      ! Dimensionless Specific Heat
          TH =TH/THS                       ! Dimensionless Conductivity
          RHO=RDUM/RHOTYP                  ! Dimensionless Density
          EMU=EDUM/EMUTYP                  ! Dimensionless Viscosity
	  RETURN
c    !-------------------------------!
      END IF
c    !-------------------------------!

 99   FORMAT( ' ',/,3X,'ERROR: on calculation of liquid properties',
     +           /,10X,'PROGRAM ABorts: Revise your Input Data')

      END

C===========================================================================
c last revision 6/28/94 by Dr. Luis San Andres
c updated 7/3/95for angled jet parameter.
C===========================================================================