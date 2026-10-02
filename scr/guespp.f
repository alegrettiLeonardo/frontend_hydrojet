C 7/7/95 update for angled jet injection
C ===========================================================================+
C
C  ####   #    #  ######   ####   #####   #####           ######
C #    #  #    #  #       #       #    #  #    #          #
C #       #    #  #####    ####   #    #  #    #          #####
C #  ###  #    #  #            #  #####   #####    ###    #
C #    #  #    #  #       #    #  #       #        ###    #
C  ####    ####   ######   ####   #       #        ###    #
C
C #    #   #   #  #####   #####    ####        #  ######   #####
C #    #    # #   #    #  #    #  #    #       #  #          #
C ######     #    #    #  #    #  #    #       #  #####      #
C #    #     #    #    #  #####   #    #       #  #          #
C #    #     #    #    #  #   #   #    #  #    #  #          #
C #    #     #    #####   #    #   ####    ####   ######     #
C
C  hydrojet.f Copyright Dr. Luis San Andres TexasA&MUniversity / 1995
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR III
c Technical monitor: Mr. James Walker, NASA Lewis Research Center
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Guesp                                                       **
C **                                                                         **
C **  GUESP  :  Guess of pressure field for concentric operation             **
C **            for multi-recess configuration                               **
C **                                                                         **
C *****************************************************************************
C Assumptions: Concentric Operation ex=ey=0 , Incompressible liquid,
C              uniform discharge pressures
C              uniform clearance (no taper, wear, elliptic, etc)
C              same recess pressure on all recesses = PRATIO
C              linear pressure drop axially
C              quadratic pressure drop between recesses
C              linear pressure drop between recess and lad/trail pad edges.
C ............................................................................

      SUBROUTINE GUESP(DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /PreLR/ Ple, Pri, Csel, Cser
      COMMON/LIQUID/ TEMPK, VSOUND, IFD, IL
      COMMON /JET/ ANGLEJ, XJET, CJET, DPJET

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /PADK/ KPAD

C    !................................................!
      DOUBLE PRECISION U,V,P, T, XP, XU, YP, YV,
     +                 PREC,TREC, QREC, QIN, QOUT, QFACTOR,
     +                 PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 Ple,Pri, Csel, Cser, TEMPK, VSOUND,
     +                 ANGLEJ, XJET, CJET, DPJET,RPM,PS,PA

      INTEGER NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL, IFD,IL,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP, ICASE, BEARING ,KPAD
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION KX, KY, UIN, Y2, DY12, REAX, VIN,HIN,
     +                 PEDGE, XL, XR, XLPR, XLR, POWER,
     +                 AOO, VVIN, PB, AO, UL, UR, PBDY,
     +                 FNKY, FNKX, VK, KXJ,KXB, CC, DUMY,
     +                 Pedgeold,Vnew,Rhout, RHO, EMU, ETA, AY,
     +                 Vout,Emout,Kyo, Klosy, DeltaP, Pexit,
     +                 Cseal, Pend, Sign,Zero,PGROOVE,Unew,PJ,
     +                 CKR,BTR,THR, XRML, PRup, PRdown
      INTEGER J, K, KP1, I, JL, JS, JR, DEVICE, Kstart, Kend, Kstep

C ----------------------------------------------------------------------------
C --  GUESP code: guesses fields for right and left sides of BEARING
C ----------------------------------------------------------------------------
      ZERO=0.0D0
      ICASE=1                             ! for right side of bearing
C    !....................................!
      Kstart=ISYM+(ISYM-1)*NYI            !
      DO I=1, NXT                         !
          DO J=Kstart,NYI,1               !
              T(I, J)=1.0D0               ! Initial guess for
          END DO                          ! DIMensionless Temperature field
      END DO                              ! = Tsupply ; i.e.  Ts/Ts=1.0
C    !....................................!

      Kstart=1
      Kend=NYI
      Kstep=1
c    !....................................! INLET (EDGE) film ratio
      IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN
         ETA=1.0D0/(1.0D0+HRECD)
      ELSE
         ETA=ZERO
      END IF
      HIN=1.0D0                        ! uniform clearance CLEAR/TYP CLEAR
c
C    !................! SET :..........! AT BEARING RIGHT PLANE
      KLOSY=KLOSYR
      Cseal=Cser
      Pexit=Pri       ! right end pressure
      Pend=Pexit      !
      PGROOVE=Pri
      IF (NPOCKET.GT.0) THEN           ! => for hydrostatic bearings
       DO I=1, NPOCKET                 ! set recess Pressure and temperature
         PREC(I)=PRATIO                !
         TREC(I)=1.0D0                 !
C##        DO J=-NYI,NYI               !
C##           PEDRISE(I,J)=ZERO        !
C##         END DO                     !
       END DO
      END IF
c    !....................................!

C    !....................................!
      IF ((MODEL.EQ.1).AND.(BEARING.EQ.1)) GOTO 333  ! .. ==> 2ROW HJB
c    !....................................!

 200  CONTINUE                            !
      IF (PRATIO.EQ.ZERO) THEN            !
         PEDGE=ZERO                       !
         VIN=ZERO                         !
         GOTO 103                         !
      END IF                              ! START GUESS OF AXIAL VELOCITY:
c                                         ! Iterative process
      CALL LOCPROPS(RHO,EMU,CKR,BTR,THR,PRATIO,1.0D0)
      Y2=YP(Kend)                         !
      DY12=(Y2-YP(NPA*Kstep))*Kstep       !=DeltaY between REC edge and exit plane
      KX=12.0D0*EMU                       ! Laminar shear coefficients
      KY=12.0D0*EMU                       !
      AY=KLOSY*(1.-ETA*ETA)*Inerp         ! inlet loss factor = 0.5(1+Xsi)REP*

      UIN=SWIRL                           ! Circ. mean velocity

      DeltaP=(PRATIO-Pexit)               ! Pressure difference
      IF (DeltaP.GT.ZERO) THEN            ! and direction of axial flow
          Sign=1.0D0                      !
      ELSE                                ! Pratio=(PS-PR)/PSA
          Sign=-1.0D0                     ! is an assumed value
          DeltaP=Sign*DeltaP              !
      END IF                              !..........................!
      DeltaP=DSQRT(DeltaP)

      REAX= REP*DeltaP*(RHO/EMU)/DY12     ! Axial pressure Reynolds Number

      KY=.066D0*REAX**0.75                ! Use Hirs formulae for approx
      KY=KY**(1.0D0/1.75D0)               ! Shear correction factor
      KY=DMAX1(KY,12.0D0*DUMY)            ! FIND Mean axial velocity
      VIN=Sign*DeltaP/(KY*DY12*EMU)       !
      KY=FNKY(UIN, VIN, HIN, SPEED, REP, RHO, EMU)
      I=0                                 !
      Vnew=VIN
      KYO=KY
      CALL LOCPROPS(Rhout,Emout,CKR,BTR,THR,Pexit,1.0D0)


 100  Vnew=Sign*DeltaP/DSQRT(DABS(
     +     Rho*(AY+Cseal*Rho/Rhout)+DY12*KY/VIN))
c  -- SETS recess edge axial pressure:
      PEDGE=PRATIO
      IF (Vnew.GT.ZERO) PEDGE=PRATIO-AY*Vnew**2*Rho

      CALL LOCPROPS(Rho,Emu,CKR,BTR,THR,PEDGE,1.0D0)

      KY=FNKY(Uin,Vnew, HIN, SPEED, REP, RHO, EMU)

      Vout=Vnew*Rho/Rhout

      IF ((Vout*Cseal).gt.ZERO ) THEN       ! IF end seal present
         PEND=Pexit+Cseal*Vout**2*Rhout     ! rise end pressure
      ELSE                                  !....
         PEND=Pexit
      END IF                                !

      CALL LOCPROPS(Rhout,Emout,CKR,BTR,THR,Pend,1.0D0)

      KYO=FNKY(MSPEED,Vout,HIN,SPEED, REP,Rhout,EMout)
      KY=(KYO+KY)/2.0D0

      IF(DABS(1.0D0-Vnew/Vin).le.(0.001D0)) goto 101   !=> Vaxial OK

      Vin=Vnew                            !
      I=I+1                               !in maximum 20 steps
      IF (I.LE.20) GOTO 100               !re-evaluate Vaxial

 101  CONTINUE
c    !....................................! PRINT OUT if requested
      IF (SVERB.GE.1) THEN                ! values GUESSED
          WRITE (6, 20) PRATIO, VIN, PEDGE, KY
      END IF                              !
      IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
          WRITE (1, 20) PRATIO, VIN, PEDGE, KY
      END IF                              !
C    !....................................!

 102  IF (PEDGE.LE.PCAV) THEN
          WRITE (6, *)
     +       'GUESSED Edge Recess Pressure TOO LOW, Increase Pratio'
          CALL BEEPER
          WRITE (6, 10)
   10     FORMAT ('$', 'INPUT: NEW(GUESSED) RECESS PRESSURE Pratio: ')
          READ (5, *) PRATIO
          PRATIO=DABS(PRATIO)
          PRATIO=DMIN1(PRATIO,1.0D0)
          IF (PRATIO.EQ.ZERO) PRATIO=0.50D0
          GOTO 200
      END IF

C    !....................................!
 103  IF (ICASE.EQ.2) Vin=-Vin            ! Left (BOTTOM) Side of HJB
C    !....................................!


C    !....................................!
      IF ((BEARING.GE.2).OR.              ! => SEAL AND PLAIN BEARING
     +    (NPOCKET.EQ.0)) GOTO 366        !    set guess field values
C    !....................................!


c    !====================================! FOR HYDROSTATIC BEARINGS
      DO I=1, NXT                         !
          DO J=Kstart,Kend,Kstep          !
              V(I, J)= ZERO               ! default values before
              P(I, J)=PRATIO              ! approximate calculations
              U(I, J)=SWIRL               !
          END DO                          !
      END DO                              !.................................

      JL=1                                ! .................................
      JR=NLC                              ! Calculate initial (guess) solution
      XL=XP(JL)                           ! X left of land between recess : JL=1
      XR=XP(JR)                           ! X right of land ''''          : JR=NLC
      XLPR=XL+XR                          !
      XRML=XR-XL                          !
      XLR=XL*XR                           !
      POWER=3.0D0                         ! V = f(Y^Power) on interlands
      AOO=POWER*VIN/(YV(Kend))**POWER     ! dV/dy on inter-lands
      VVIN=AOO/POWER                      ! V = VVIN * YV^3
c                                         ! NXI=NLC+NPC-2 : #nodes in X-recess+land
      JS=NXI+1                            ! at other end of recess JS=JL+NXI
      PB=PEDGE                            ! Pressure at L&R boundary
      DPJET=CJET*(1.0D0-PRATIO)           ! DeltaP due to jet injection
      PRdown=PB-DPJET*(1.0D0-XJET)        ! downstream reces pressure
      PRup  =PB+DPJET*XJET                ! upstream recess pressure
C                                         ! XJET=jet circ. loc. 0:left,.5:mid,1:0 right
C                                         ! ---------------------------------
       DO K=Kstep,Kstep*(NPA-1),Kstep     ! On inter-recess lands
          KP1=K+Kstep                     ! ---------------------------------
          VK=VVIN*(YV(KP1))**POWER        ! V velocity at YV
          AO=AOO*(YP(K))**(POWER-1.0)     ! dV/dy at k
          UL=UIN+AO*(XR-XL)/2.0D0         ! U left border, UIN=SWIRL
          UR=UIN-AO*(XR-XL)/2.0D0         ! U right border
          P(JL, K)=PRdown                 !  pressures at left & right
          P(JR, K)=PRup                   !     recesses
          P(JS, K)=PRdown                 !

          DUMY=1.0D0  !=CLEAR             ! Calculate aprox. U, V & P
          KX=FNKX(UL,VK,DUMY,SPEED,REP,RHO,EMU,KXJ,KXB) ! Calc. appx. U, V & P
          CC=AO*KX/2.0D0                  ! Kxl=Kxr

          DO J=JL+1, JR-1                 ! SWEEP on inter-recess LAND
              P(J, K)= CC*(XP(J)*(XP(J)-XLPR)+XLR)    ! Quadratic P field
     +               + ( PRdown*(XR-XP(J)) +          ! between recesses
     +                     PRup*(XP(J)-XL)    )/XRML  !
              V(J, KP1)=VK                            ! constant Vaxial
              V(J, KP1)=VK                            !
              U(J-1, K)=UIN-AO*(XU(J-1)-XLPR/2.0)/2.0
          END DO                          !
          U(JR-1, K)=UIN-AO*(XU(JR-1)-XLPR/2.0)/2.0
      END DO                              ! ----------------------------------
                                          !
                                          ! .................................
      K=Kstep*NPA                         ! On AXIAL edge of inter reccess lands
      KP1=K+Kstep                         ! ..................................
      VK=VVIN*(YV(KP1))**POWER            ! V velocity at YV
      AO=AOO*(YP(K))**(POWER-1.0)         ! dV/dy at k
      UL=UIN+AO*(XR-XL)/2.0D0             ! U left border
      UR=UIN-AO*(XR-XL)/2.0D0             ! U right border
      PB=PEDGE                            !  = Edge recess pressure
      P(JL, K)=PRdown                     ! P on left boundary
      P(JR, K)=PRup                       ! P on right boundary

      DUMY=1.0D0                          !
      KX=FNKX(UL,VK,DUMY,SPEED,REP,RHO,EMU,KXJ,KXB) ! Calc. appc U, V & P
      CC=AO*KX/2.0D0                      ! Kxl=Kxr
      DO J=JL+1, JR-1                     !
              P(J, K)= CC*(XP(J)*(XP(J)-XLPR)+XLR)    ! Quadratic P field
     +               + ( PRdown*(XR-XP(J)) +          ! between recesses
     +                     PRup*(XP(J)-XL)    )/XRML  !
          V(J, KP1)=VK                    !
          U(J-1, K)=UIN-AO*(XU(J-1)-XLPR/2.0D0)/2.0D0
      END DO                              !
      U(JR-1, K)=UIN-AO*(XU(JR-1)-XLPR/2.0D0)/2.0D0


c    !....................................!................................
      JS=NXI+1                            !=NXI+JL
      XL=XP(JR)                           ! X left of recess :
      XR=XP(JS)                           ! X right of recess:
      XRML=XR-XL                          !
C##    DO K=Kstep,Kstep*NPA               !
       DO J=JR+1, JS                      ! .................................
        P(J, K)=( PRup*(XR-XP(J)) +       ! ON AXIAL edge of recess K=NPA
     +            PRdown*(XP(J)-XL))/XRML ! linear pressure change within recess
       END DO                             ! .................................
C##    END DO


      DY12=Kstep*DY12                     ! YP(Kend)-YP(Kstep*NPA)
      Y2=Y2*Kstep                         ! YP(Kend=NYI*Kstep)
c                                         !..................................
      DO J=JL, JS                         ! On extended lands, JL=1, JS=NXI+1
          DO K=Kstep*(NPA+1),Kend,Kstep   !..................................
              DeltaP=Pend-P(J,Kstep*NPA)  !
              P(J,K)=Pend+DeltaP*(YP(K)-Y2)/DY12
              V(J,K)=VVIN*(YV(K))**POWER  ! U = SWIRL
          END DO                          !
      END DO                              !
C                                         !====================================
C.........................................!
      DO I=2, NPOCKET                     ! Copy U, V & P field to other recess
          JL=(I-1)*NXI+2                  ! ...................................
          JR=I*NXI+1                      ! IFULL=1 , PIxDIAM LANDS
          JS=1                            !
          DO J=JL, JR                     !
              JS=JS+1                     !
              DO K=Kstep,Kend,Kstep       !
                  P(J, K)=P(JS, K)        !
                  U(J, K)=U(JS, K)        !
                  V(J, K)=V(JS, K)        !
              END DO                      !
          END DO                          !
      END DO                              !

      DO J=1, NXT                         ! On outlet boundary
          P(J,Kend)=Pend                  !
      END DO                              !

      IF (IFULL.EQ.1) GOTO 444            !=> 360 DEG HJB (1pad)
C                                         ! ...................................
C.........................................! IFULL=0, a hydrostatic pad
C                                         !
      JL=1                                ! .................................
      JR=NLC                              ! ON land between Groove & RECESS edge
      XL=XP(JL)                           ! X left of inter-land
      XR=XP(JR)                           ! X right of inter-land
      XRML=XR-XL                          !

      PB=PEDGE                            ! Pressure at R boundary (recess)
      DPJET=CJET*(1.0D0-PRATIO)           ! DeltaP due to jet injection
      PRup  =PB+DPJET*Xjet                !
      DELTAP=PRup-PGROOVE
      VK=ZERO
      UL=UIN
      I=0

 1010 KX=FNKX(UL,VK,HIN,SPEED,REP,RHO,EMU,KXJ,KXB)
      Unew=MSPEED*(KXJ/KX)-DELTAP/(XRML*KX)
      IF(DABS(1.0D0-Unew/UL).le.(0.001D0)) goto 1011
      UL=Unew
      I=I+1
      IF (I.LE.20) GOTO 1010

 1011 DO K=Kstep,Kstep*NPA,Kstep          ! ON land between Groove & RECESS edge
          DO J=JL, JR                     ! use linear pressure drop
              P(J, K)=DELTAP*(XP(J)-XL)/XRML+PGROOVE
              U(J, K)=UL
              V(J, K)=VK
          END DO
      END DO                              !..................................

      JL=NXI*NPOCKET+1                    ! .................................
      JR=NXT                              ! ON land between RECESS edge & groove
      XL=XP(JL)                           ! X left of inter-land
      XR=XP(JR)                           ! X right of inter-land
      XRML=XR-XL                          !

      PB=PEDGE                            ! Pressure at L boundary
      PRdown=PB-DPJET*(1.0D0-Xjet)        !
      DELTAP=PRdown-PGROOVE
      VK=ZERO
      UR=UIN
      I=0
      UR=(UIN*KXJ+DELTAP/XLR)/KX
 1012 KX=FNKX(UR,VK,HIN,SPEED,REP,RHO,EMU,KXJ,KXB)
      Unew=MSPEED*(KXJ/KX)+DELTAP/(XRML*KX)
      IF(DABS(1.0D0-Unew/UR).le.(0.001D0)) goto 1013
      UR=Unew
      I=I+1
      IF (I.LE.20) GOTO 1012

 1013 DO K=Kstep,Kstep*NPA,Kstep          ! ON land between RECESS edge & groove
          DO J=JL, JR                     !
              P(J, K)=-DELTAP*(XP(J)-XL)/XRML+PRdown
              U(J, K)=UR
              V(J, K)=VK
          END DO
      END DO                              !..............................!

      DO K=Kstep*(NPA+1), Kend, Kstep     ! ON EXTENDED LANDS
         I=1
         DO J=JL, JR
           U(J,K)=SWIRL
           V(J,K)=V(NLC-I+1,K)
           P(J,K)=P(NLC-I+1,K)
           I=I+1
         END DO
      END DO

      GOTO 444
C.........................................!.........................
C     DOUBLE ROW HYDROSTATIC BEARING
C.........................................!.........................
 333  DO I=1, NXT                         ! Right Side
          DO J=Kstart,Kend,Kstep          ! is symmetry side for 2-row HJB
              V(I, J)=ZERO                !
              P(I, J)=PRATIO              !
              U(I, J)=SWIRL               !
          END DO                          !
      END DO                              !
      GOTO 444                            !
C.........................................!.........................
C     ANNULAR SEAL OR PLAIN BEARING       ! BEARING=2 or 3
C.........................................!.........................
 366  DO J=KSTART,KEND,KSTEP
       PJ=PEDGE*(1.0D0-YP(J)/YP(KEND))+PEND*(YP(J)/YP(Kend))
        DO I=1, NXT
            V(I,J)=VIN
            U(I,J)=SWIRL
            P(I,J)=PJ
        END DO
      END DO

C.........................................!.........................
      IF (BEARING.EQ.2) GOTO 777          ! i.e. an annular SEAL
C.........................................!.........................


c     !---------------------!
 444   IF (ISYM.EQ.0) THEN
c     !---------------------!  for asymmetric HJB & disch pressures
       IF (Icase.eq.1) THEN
        Kstep=-1
        Kstart=Kstep
        Kend=-NYI
        KLOSY=KLOSYl
        Cseal=Csel
        Pexit=Ple
        Pend=Ple
        PGROOVE=Ple
        Icase=2
        goto 200
       ELSE
	do j=1, NXT
          P(j,0)=0.5D0*(P(j,1)+P(j,-1))
          U(j,0)=0.5D0*(U(j,1)+U(j,-1))
          V(j,0)=0.5D0*(V(j,1)+V(j,-1))
        end do
       END IF
c     !---------------------!
       END IF
c     !---------------------!

C.........................................!.........................

C     !.....................!Specify Known BCs on bearing pad
 777   CALL SETBC(RPM,PS,PA)
C...........................!..............!.........................!
C##      PRINT *, 'GUESPP: Initial Pfields:'      ! ==> FOR DEBUGGING
C##        CALL PFIELDS(Device)		       	  !     PROCESS
C##      CALL PAUSE				  !     ONLY

 888  IF (KPAD.EQ.1) THEN
        write (6,*) '$ GUESP OK => CALCULATIONS START hydrojet/95'
        write (6,*) '$ ------------------------------------------'
      END IF

   20 FORMAT (' ', /, 4X, 'Concentric pressure ratio, Pratio:',
     +        E12.5E2, /, 4X, 'GUESSED: Vinlet=', E12.5E2, 3X,
     +        'Pedge=', E12.5E2, /, 4X, 'Shear Factor Ky=',
     +        E12.5E2)

      END

c ---------------------------------------------------------------- c
c revised 12/09/93 by Luis San Andres FOR SEALS/HJBs/JBs
c updated  7/7/95  for angled jet injection LSA
c ---------------------------------------------------------------- c
