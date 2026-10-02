c modifications for angled jet started on 9/13/95
c call to FORCE1 modified
c ===========================================================================

C ######     #    #####    ####    #####  #####           ######
C #          #    #    #  #          #    #    #          #
C #####      #    #    #   ####      #    #    #          #####
C #          #    #####        #     #    #####    ###    #
C #          #    #   #   #    #     #    #        ###    #
C #          #    #    #   ####      #    #        ###    #

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
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c ===================================
C NOTE: Important Modification 9/21/93
C ===================================
c
C For compressible fluids, including gases, as well as for tilting pad
c bearings and foil bearings, the force coefficients show
c a peculiar dependency on frequency. That is, the simple quadratic formulae
c to identify force coefficients as
c
c    F/X = (K - M w**2) + i C w    IS NOT VALID !
c
c Therefore, I have decided to present the force coefficients at the
c frequency of interest selected by the user.
c That is, no longer mass coefficients ARE available
c
c The force coefficients then will represent actual force impedances at
c the frequency of interest, such that
c
c F/X = Kd + i Cd w
c
c where Kd, Cd are the dynamic stiffness and damping coefficients at frequency w
c
c If you still want mass coefficients and your fluid is nearly incompressible
c then run program with opt (35) and two different frequencies well appart and
c
c determine M = (K(w<>0)-K(w=0))/w**2
c
c The only change made is set
c.............................................!..............................!
c         INERFREQ=0                          ! 9/21/93 Calculate coeffs at w always
c.............................................!..............................!
C ============================================! on COEFFIC and
C                                                  PRINTCOEF routines



C *****************************************************************************
C **                                                                         **
C **  Subroutine Coeffic                                                     **
C **                                                                         **
C **  COEFFIC:  Finds solution of first order equations for given            **
C **            frequency w and calculates dynamic coefficients for set of   **
C **            L2 parameters.                                               **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE COEFFIC(DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /AIR/ CP, RC , AIRMU, AGAMA, PR ,KS
      COMMON /PERPARM/ FXX(MAXNPOCKP1), FYY(MAXNPOCKP1),
     +                 MXX(MAXNPOCKP1), MYY(MAXNPOCKP1),
     +                 QZERO(MAXNPOCK), QRECREC(MAXNPOCK, MAXNPOCK)
      COMMON /PERPARM0/ FXXO(MAXNPOCKP1), FYYO(MAXNPOCKP1),
     +                  MXXO(MAXNPOCKP1), MYYO(MAXNPOCKP1),
     +                  QZEROO(MAXNPOCK), QRECRECO(MAXNPOCK, MAXNPOCK)
      COMMON /XPERTURB/ XFXX, XFYY, XMXX, XMYY, QXZERO(MAXNPOCK)
      COMMON /YPERTURB/ YFXX, YFYY, YMXX, YMYY, QYZERO(MAXNPOCK)
      COMMON /XPERTURB0/ XFXXO, XFYYO, XMXXO, XMYYO, QXZEROO(MAXNPOCK)
      COMMON /YPERTURB0/ YFXXO, YFYYO, YMXXO, YMYYO, QYZEROO(MAXNPOCK)

      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /RECASP/ ASPEC(MAXNPOCK)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /PARAM1/ CLEAR,DIAM,LENGTH,LD,AR, HRECC
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /THERMID/ ISOTH
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /RECPAR/ HREC, VSUP, BETA
      COMMON /COEFS/ K11, K12, C11, C12, KM11,KM12,CM11,CM12
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILTPAD

      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /FLMOM/ IMOMFLAG

      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK),
     +                QREC1(MAXNPOCK), QIN1, QOUT1
      COMMON /QRECT1/ CRECT(MAXNPOCK),CRECTR(MAXNPOCK)
c     ..................................................................
      DOUBLE COMPLEX FXX, FYY, MXX, MYY, QZERO, QRECREC,
     +               FXXO,FYYO,MXXO,MYYO,QZEROO,QRECRECO,
     +               XFXX, XFYY, XMXX, XMYY, QXZERO,
     +               XFXXO,XFYYO,XMXXO,XMYYO,QXZEROO,
     +               YFXX, YFYY, YMXX, YMYY, QYZERO,
     +               YFXXO,YFYYO,YMXXO,YMYYO,QYZEROO,
     +               prec1,trec1,qrec1,qin1,qout1,
     +               CRECT,CRECTR
      DOUBLE PRECISION XP, XU, YP, YV, COSXP, SINXP, COSXU, SINXU,
     +                 PREC,TREC, QREC, QIN, QOUT, QFACTOR,
     +                 DIAORIF, CORIF,alft,uc,tc,ec ,
     +                 RHOP,EMUP,CK,BETAK
      DOUBLE PRECISION EXO, EYO, AXO, AYO, ZO,
     +                 CP,RC,AIRMU,AGAMA,PR,KS ,CPS,THS,BETAKS,
     +                 CLEAR,DIAM,LENGTH,LD,AR,BR,HRECC,ASPEC,
     +                 FREQU, SIGMA, L1, RES,TEMPK, VSOUND,
     +                 HREC, VSUP, BETA, CINLET, CEXIT,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2, PMIN, PMAX
      DOUBLE PRECISION K11, K12, C11, C12, KM11,KM12,CM11,CM12,
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD

      INTEGER IF, IL, ICASE, NCASE, ITER, ITMAX, ITPMAX,ISOTH,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL,TILTPAD,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP, IMOMFLAG, BEARING


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE COMPLEX CREC(MAXNPOCK), CRECR(MAXNPOCK),
     +               FXcom(MAXNPOCKP1), FYcom(MAXNPOCKP1),
     +               MXcom(MAXNPOCKP1), MYcom(MAXNPOCKP1),
     +               QINV(MAXNPOCK,MAXNPOCK),QINVO(MAXNPOCK,MAXNPOCK)
      DOUBLE PRECISION FUNRY(MAXNPOCK), FUNRX(MAXNPOCK), FREQ2,
     +                 PI, AREA, FREQ, RADIUS, CRATIO, CSTIF, CDAMP,
     +                 CORIF2,XR,TR,TL,OMEGA,HRECU,L2,L2s,Betas,L4s,
     +                 bkr,dbkr ,cpr,thr,emur,CINER
      DOUBLE PRECISION KdumD, KdumC, YREC, YO,
     +                 L3, HR, VR, CRECI, L11, ZERO,Deltp, rhord,
     +                 L4, L4c, creco, rhor, whirl,
     +                 pro,Z,Fz,Gama,ALfc,Aori,Vori,Wbreak
      INTEGER I,JMI,JMA,DEVICE,J,L,INERFREQ,IFM,ISYMold,Imod
C ----------------------------------------------------------------------------
C --  COEFFIC code                                                          --
C ----------------------------------------------------------------------------

C     here check if bearing is unloaded to set coefficients = zero
c    !......................! Calculate MIN and MAX Pressures.

      CALL MINMAXP
c    !......................!
      IF (PMAX.EQ.PC) THEN
c    !......................! PAD is unloaded
      KXXD=0.0D0            ! ALL Force Coefficients are NULL
      KYYD=0.0D0
      KXYD=0.0D0
      KYXD=0.0D0
      CXXD=0.0D0
      CYYD=0.0D0
      CXYD=0.0D0
      CYXD=0.0D0
      MXXD=0.0D0
      MYYD=0.0D0
      MXYD=0.0D0
      MYXD=0.0D0
      KXXA=0.0D0
      KYYA=0.0D0
      KXYA=0.0D0
      KYXA=0.0D0
      CXXA=0.0D0
      CYYA=0.0D0
      CXYA=0.0D0
      CYXA=0.0D0
      MXXA=0.0D0
      MYYA=0.0D0
      MXYA=0.0D0
      MYXA=0.0D0
      KmXXD=0.0D0
      KmYYD=0.0D0
      KmXYD=0.0D0
      KmYXD=0.0D0
      CmXXD=0.0D0
      CmYYD=0.0D0
      CmXYD=0.0D0
      CmYXD=0.0D0
      MmXXD=0.0D0
      MmYYD=0.0D0
      MmXYD=0.0D0
      MmYXD=0.0D0
      KmXXA=0.0D0
      KmYYA=0.0D0
      KmXYA=0.0D0
      KmYXA=0.0D0
      CmXXA=0.0D0
      CmYYA=0.0D0
      CmXYA=0.0D0
      CmYXA=0.0D0
      MmXXA=0.0D0
      MmYYA=0.0D0
      MmXYA=0.0D0
      MmYXA=0.0D0
      KdXk=0.0D0
      KdYk=0.0D0
      KXdk=0.0D0
      KYdk=0.0D0
      Kddk=0.0D0
      CdXk=0.0D0
      CdYk=0.0D0
      CXdk=0.0D0
      CYdk=0.0D0
      Cddk=0.0D0
      RETURN
c    !......................!
      END IF
c    !......................!

c ----------------------------------------------------------------------------
      Zero=0.0D0
      ISYMold=ISYM
      PI=DACOS(-1.0D0)
c    !........................................!
      BR=Zero
      IF (NPOCKET.GT.0) THEN
       DO I=1, NPOCKET
        BR=BR+ASPEC(I)
       END DO
       BR=BR*DIAM/NPOCKET                     ! average circ. recess length
      END IF
c............................................ ! .............................
C                                             ! CALCULATE OPERATING PARAMETERS
C                                             !...............................
          YO=2.0D0*ZO/DIAM                    ! pivot loc. for moments
          HRECU=HREC                          !
          AREA=AR*BR                          ! Pocket area, freq in Hz
          FREQ=FREQU*2.0D0*PI                 ! w in rad/sec
          FREQ2=FREQ*FREQ                     !
          RADIUS=DIAM/2.0D0                   ! D/2
          CRATIO=CLEAR/RADIUS                 ! c/R
          SIGMA=EMU*FREQ/PSA/CRATIO**2.0      ! u w R2/(Ps-Pa)/c2
          L1=EMU*FREQ*AREA/PSA/CLEAR**2.0     ! u w Area/(Ps-Pa)/c2
          RES=RHO*FREQ*CLEAR*CLEAR/EMU        ! Squeeze Reynolds Number
          L4=PSA*BETA                         ! compressib. ratio
          OMEGA=RPM/60.0D0                    ! Rotational speed in Hz
          IF (BEARING.EQ.1) THEN
          L11=L1*RADIUS/BR                    ! Modified L1 parameter
          END IF
          CSTIF=PSA*RADIUS*RADIUS/CLEAR       ! Stiffness conversion coeffic. [N/m]
          CDAMP=Zero

          !WRITE(6,*) CSTIF, RPM, FREQ
          IF (FREQ.NE.Zero ) then
            CDAMP=CSTIF/FREQ ! Damping Conv. Coef. [Ns/m]
            CINER=CDAMP/FREQ ! Inertia Coef. [Ns˛/m]
            !WRITE(6,*) PSA, RPM, CINER
        END IF
C
c.............................................!..............................

          INERFREQ=0                          !
          IF ((INERL.EQ.1).AND.(FREQ.NE.Zero)) INERFREQ=1

C#############################################
c.............................................!..............................!
          INERFREQ=0                          ! 9/21/93 Calculate coeffs at w always
c.............................................!..............................!
C##       WRITE (6,111) FREQU
 111      FORMAT ('$ Impedances evaluated at frequency: ',
     +             E12.5E2,' Hz ==> INERFREQ=0',/,
     +            1X,78('+'))
C#############################################

          RES=RES*INERL
c.............................................!..............................
c ICASE=1: X perturbation, =2: X,Y perturbation
          IF (ICASE.GT.2) ICASE=2
          IF (ICASE.LT.1) ICASE=1

c.............................................!..............................
c for tilt-pad calculate force impedances only at given frequency
c              and perturb in X & Y directions
c              no-moment journal rotations are also allowed
c ..................................................................

          IF (TILTPAD.EQ.1) THEN
              INERFREQ=0
              ICASE=2
              GOTO 77
          END IF

c.............................................!..............................
          IF (IFULL.EQ.1) THEN                ! 360 FULL BEARING
c        !---------------------!              !......................
          IF (ICASE.EQ.1) THEN                ! Exo=Eyo=Axo=Ayo=0
             IF ((AXO.EQ.zero).AND.(AYO.EQ.zero)) THEN
                 ICASE=1
                 GOTO 33
             ELSE
                 ICASE=2
                 GOTO 77
             END IF
          ELSE                                ! ICASE=2
             GOTO 77
          END IF
c        !---------------------!              !......................
          ELSE                                ! PAD BEARING IFULL=0
c        !---------------------!              !......................
             ICASE=2
             GOTO 77
c        !---------------------!              !......................
          END IF
c        !---------------------!              !......................


c.............................................!.............................
c 	Determine Break Frequency and Reduced Damping Factor
c 	for centered HJBs with symmetric conditions only
c 	and identical orifice diameter (average value taken).
c 	Based on Dr. Luis San Andres paper:
c 	FLUID COMPRESSIBILITY EFFECTS ON THE DYNAMIC RESPONSE OF
C 	HYDROSTATIC JOURNAL BEARINGS, WEAR , Vol. 146, pp. 269-283, 1991.
c
c.............................................!.............................

 33    IF ((BETA.EQ.ZERO).OR.(RPM.EQ.ZERO).OR.
     +     (ISYM.EQ.0)   .OR.(NPOCKET.EQ.0)) GOTO 77

       CALL FINPROP(ICASE,DEVICE,RPM,PS,PA)         ! calculate props/derivatives
                                          ! ..........................
          pro=prec(1)
          Z=0.5D0*pro/(1.0D0-pro)
          Gama=(Length-Ar)*Length*Npocket*Npocket/
     /         (PI*(2.0D0*PI*Diam-Npocket*Br)*Diam)

          Fz=Z+1.0D0+2.*Gama*DSIN(PI/Npocket)**2.0
          Vr=Area*(Hrec+Cinlet)+Vsup

          call LOCPROPS(rhor,emur,CPR,BKR,THR, pro,TREC(1))

          IF (IF.eq.12.OR.IF.eq.6.OR.IF.eq.7) THEN
           Vsound=1.0D0/DSQRT(Beta*Rhor*RHO)
          ELSE
           Deltp=0.001D0
           call LOCPROPS(rhord,emur,CPR,DBKR,
     +     THR, pro+deltp,TREC(1))
           L4=(rhord-rhor)/deltp/rhor
           Beta=L4/(Psa)
          END IF
c                                              ! Reduced Damping factor

          Alfc=3.0*Npocket*Vr*pro*L4/(PI*Length*Diam*Clear)/Fz

c                                              ! Break frequency

          Aori=Cd*(PI/4.0D0)*Dorif*Dorif
          Vori=DSQRT(2.0D0*(Psa)*(1.0D0-Pro)/(Rhor*RHO))
          Wbreak=Aori*(Vsound**2)/(Vr*Vori)
          Wbreak=Wbreak*(Fz/Z)/(2.0D0*PI)
          Whirl=0.5D0/(1.0D0-ALfc)             ! Approximate whirl ratio

        IF (alfc.gt.(1.0D0)) THEN
          write (6,1190)
          if ((Device.eq.1).AND.(Dverb.ge.1)) THEN
             write(1,1190)
          end if
        ELSE
          WRITE (6,1200 ) Wbreak, Alfc, Whirl
          IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN
             WRITE (1, 1200) Wbreak, Alfc, Whirl
          END IF
        END IF
        GOTO 79

c............................................................................

 77     CALL FINPROP(ICASE,DEVICE,RPM,PS,PA)         ! calculate props/derivatives
C                                          !.............=> on findprop.f prg

 79     CONTINUE
C##     WRITE (6, 1300) FREQU, SIGMA, RES
          IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN
             WRITE (1, 1300) FREQU, SIGMA, RES
        END IF                                !
c ............................................!

      IF (NPOCKET.EQ.0) GOTO 88
C
C ------------------------------------------- ! ------------------------------
C Set vectors of recess functions             !
C ............................................!

      Betas=0.0D0
      L2s=0.0D0
      L4s=0.0D0


C    !---------------!                        !
      DO I=1, NPOCKET                         ! ::::::::::::::::::::::::::::::
c    !---------------!                        !
       Corif2=CORIF(I)*CORIF(I)               ! Orifice coefficient
       call LOCPROPS(rhor,emur,CPR,BKR,       !
     +            THR, prec(i),TREC(I))       !
c      .......................................!
c      compressibility parameter :    Beta=(dRho/DPres)*(1/Rho)

c      !......................................!.............................
         IF ((IF.EQ.12).OR.(IF.EQ.7)) THEN    ! for unknown fluid
c      !......................................! or oil
           L4=Beta*(Psa)

           IF(DABS(Rhoa/Rhos-1.0D0).lt.(0.01D0)) THEN
                 Creco=rhor*Corif2/2.0D0      ! Incompressible &
           ELSE                               ! compressible liquids
                 Creco=rhor*Corif2*(1.0D0-(1.0D0-prec(i))*L4)/2.0D0
           END IF
c      !......................................!.............................
         ELSE                                 ! for  LH2,LO2,LN2,LCH4 & air
c      !......................................!.............................
           Deltp=0.001D0                      ! barotropic cryogens
           call LOCPROPS(rhord,emur,CPR,
     +     DBKR,THR,prec(i)+deltp,TREC(I))    !
           L4=(rhord-rhor)/deltp/rhor
           Beta=L4/Psa
           Creco=rhor*Corif2*(1.0D0-(1.0D0-prec(i))*L4)/2.0D0
c      !......................................!.............................
         END IF
c      .......................................!

      CRECR(I)=DCMPLX(Creco/DABS(Qrec(i)),zero)  ! Real part of recess conductances
      CRECTR(I)=(1.D0-PREC(I))*RHOR              ! due to press and temperature
     +   *Corif2*BKR/Qrec(i)/2.D0                !................................!

C                                             ! ..............................
          JMI=(I-1)*NXI+NLC                   ! Functions required for
          JMA=I*NXI+1                         !  Funrx=INTEGRAL(RhorxCOSx) on rec
          TR=XP(JMA)                          !  Funry=INTEGRAL(RhorxSINx) on rec
          TL=XP(JMI)                          !
          XR=(TR+TL)/2.0D0                    ! Xr: coordinate of recess center

          FUNRX(I)= (DSIN(TR)-DSIN(TL))*rhor  !
          FUNRY(I)=-(DCOS(TR)-DCOS(TL))*rhor  !

          HR=CLEAR*(1.0D0+EXO*DCOS(XR)+EYO*DSIN(XR))

          BR=ASPEC(I)*DIAM
          AREA=AR*BR                          ! Recess area
          VR=AREA*(HREC+HR)+VSUP              ! Recess volume at Xr
          L2=L4*VR/CLEAR/AREA                 ! = B (Ps-Pa) Vr/c/Ar
          L3=EMU*FREQ*BETA*VR/CLEAR/CLEAR**2  ! = L1xL2 =u w Beta Vro/c3
          L2s=L2s+L2
          Betas=Betas+Beta
          L4s=L4s+L4
          CRECI=L3*Rhor                       ! Imag of Recess conductance
          CREC(I)=CRECR(I)+DCMPLX(ZERO,CRECI)
          CRECI=RHOR*BKR*SIGMA*VR/(CLEAR*RADIUS*RADIUS)
          CRECT(I)=CRECTR(I)-DCMPLX(ZERO,CRECI)

C    !---------------!         ...............! NEXT POCKET
      END DO                                  !
C    !---------------!                        !..............................!
C                                             !
           Betas=Betas/NPOCKET
           L2s=L2s/NPOCKET
           L4s=L4s/NPOCKET

c        !....................................!.......
          IF(Beta.eq.zero) THEN
              IF (DEVICE.EQ.1) WRITE(1,1400)
          ELSE
              IF (DEVICE.EQ.1) WRITE (1, 1550) BETAs, HREC,L2s,L4s
          END IF
c        !....................................!



C-------------------------------------------------------------------------

  88      IF ((ISYM.EQ.1).and.(ZO.NE.ZERO).and.
     +        (BEARING.NE.2) ) CALL TEMPASYM    ! => ISYM=0 temporal

C...............................................!
C Calculate coefficients for U&V&T first order eqns.
C...............................................!

           CALL UCOEFS               ! => compe1t.f
           CALL VCOEFS               !
           IF (ISOTH.NE.1) THEN      ! FOR Thermal SOLUTION
             CALL TCOEFS             ! => calcthert.f
           END IF                    !


C-------------------------------------------------------------------------
C         FORCE & MOMENT COEFFICIENTS FOR JOURNAL CENTER DISPLACEMENTS
C-------------------------------------------------------------------------

          IFM=0
          CALL PERTXY(YO,IFM,DEVICE,INERFREQ) ! Perturbs in X&Y direction

C-------------------------------------------------------------------------
C FORM Matrix of FLOWS: Crec Delta|j,k + Q1|k,j
C      and find INVERSE MATRIX
C-------------------------------------------------------------------------

 4411   IF (NPOCKET.EQ.0) GOTO 4412

          DO I=1, NPOCKET
              QRECREC(I,I)=QRECREC(I,I)+CREC(I)
              FXcom(I+1)=FXX(I+1)                  ! Save Component Force
              FYcom(I+1)=FYY(I+1)                  ! and Moments due to
              MXcom(I+1)=MXX(I+1)                  ! PR1=1+i1
              MYcom(I+1)=MYY(I+1)
          END DO
          CALL CINVERSE(NPOCKET,QRECREC,QINV)      ! ==> supportp.f
c        !...........................................!
          IF (INERFREQ.EQ.1) THEN
            DO I=1, NPOCKET
                QRECRECO(I,I)=QRECRECO(I,I)+CRECR(I)
            END DO
            CALL CINVERSE(NPOCKET,QRECRECO,QINVO)  !==> supportp.f
          END IF

          DO I=1, NPOCKET                     !
              QZERO(I)=QXZERO(I)
          END DO                              !
C-------------------------------------------------------------------------
C Calculate Coefficients: Kxx, Kyx, Cxx, Cyx, Mxx, Myx > PERTURB displac.
C                       Kaxx,Kayx, Caxx, Cayx,Maxx,Mayx> along X direction
C-------------------------------------------------------------------------
 4412      FXX(1)=XFXX                        !
           FYY(1)=XFYY                        !
           MXX(1)=XMXX                        !
           MYY(1)=XMYY                        !

C ------------------------------------------- ! --------------------------
          CALL KIJCIJ(FUNRX, L11, QINV, NPOCKET)
C-------------------------------------------- ! --------------------------
          !write (6,*) FREQU
          KXXD=K11*CSTIF                      !  FORCE coefficients
          !WRITE(6,*) KXXD
          KYXD=K12*CSTIF                      !  Stiffness in [N/m]
          CXXD=C11*CDAMP                      !  Damping in [Ns/m]
          CYXD=C12*CDAMP                      !  due to X displacement                    !
          MXXD=ZERO                           !
          MYXD=ZERO
c        !....................................!
          KmXXD=KM11*CSTIF*Radius             !  MOMENT coefficients
          KmYXD=KM12*CSTIF*Radius             !  Stiffness in [N]
          CmXXD=CM11*CDAMP*Radius             !  Damping in [N-s]
          CmYXD=CM12*CDAMP*Radius             !  due to X displacement
          MmXXD=zero                          !
          MmYXD=zero                          !
c        !....................................!

        !IF (FREQ.EQ.0) THEN
        !    INERFREQ=1
        !END IF

        !WRITE(6,*) INERFREQ
        !INERFREQ=1
        !WRITE(6,*) INERFREQ
                                 !
C        !........................            !
          IF (INERFREQ.EQ.1) THEN             ! at frequency=0.0
C        !........................            !
             IF (NPOCKET.GT.0) THEN           !
              DO I=2, NPOCKET+1               !
                  FXX( I )= FXXO( I )         !
                  FYY( I )= FYYO( I )         !
                  MXX( I )= MXXO( I )         !
                  MYY( I )= MYYO( I )         !
                  QZERO(I-1)=QXZEROO(I-1)     !
              END DO                          !
             END IF                           !
              FXX(1)=XFXXO                    !
              FYY(1)=XFYYO                    !
              MXX(1)=XMXXO                    !
              MYY(1)=XMYYO                    !

C            !.......................................!
              CALL KIJCIJ(FUNRX,zero, QINVO, NPOCKET)

              !WRITE(6,*)FREQ2
              !FREQ2 = (RPM/60)*(RPM/60)
              !WRITE(6,*)FREQ2
C            !.......................................!

              KdumD=K11*CSTIF                 ! FORCE coefficients
              KdumC=K12*CSTIF                 ! for frequency w=0.
c             CXXDO=C11*CDAMP                 ! due to X displacement
c             CYXDO=C12*CDAMP                 !
              MXXD=(KdumD-KXXD)/FREQ2         ! Inertia coefficients in N-s2/m
              MYXD=(KdumC-KYXD)/FREQ2         !
              !WRITE(6,*) MXXD
              KXXD=KdumD                      !
              KYXD=KdumC                      !
c            !................................!
              KdumD =KM11*CSTIF*Radius        ! MOMENT coefficients
              KdumC =KM12*CSTIF*Radius        ! for frequency w=0.
c             CmXXDO=CM11*CDAMP*Radius        ! due to X displacement
c             CmYXDO=CM12*CDAMP*Radius        !
              MmXXD=(KdumD -KmXXD)/FREQ2      ! Inertia coefficients in N-s2
              MmYXD=(KdumC -KmYXD)/FREQ2      !
              KmXXD=KdumD                     !
              KmYXD=KdumC                     !

         IF (NPOCKET.GT.0) THEN                !.....
          DO I=2, NPOCKET+1                    !
              FXX(I)=FXcom(I)                  ! RE-STORE Component Force
              FYY(I)=FYcom(I)                  ! and Moments due to
              MXX(I)=MXcom(I)                  ! PR1=1+i1
              MYY(I)=MYcom(I)                  !
          END DO                               !
         END IF                                !
C        !.......!                            !
          END IF ! -------------------------- ! Inerfreq=1
C        !.......!                            !
C
C                                             !
          IF (ICASE.EQ.1) THEN                !
c        !.....................!              !
              KYYD=KXXD                       !
              KXYD=-KYXD                      ! for concentric position
              CYYD=CXXD                       ! Ex=Ey=0.
              CXYD=-CYXD                      !
              MYYD=MXXD                       ! coefficients are skew-symmetric
              MXYD=-MYXD                      !
              KmYYD=KmXXD                     !
              KmXYD=-KmYXD                    ! for concentric position
              CmYYD=CmXXD                     ! Ex=Ey=0.
              CmXYD=-CmYXD                    ! AND YO=0  !#### CHECK
              MmYYD=MmXXD                     !
              MmXYD=-MmYXD                    !
              GOTO 402                        !=> Angle Perturbations
          END IF                              !
C        !......................!             !
C
C-------------------------------------------------------------------------
C Calculate Coefficients: Kxy, Kyy, Cxy, Cyy, Mxy, Myy > PERTURB displac.
C                       Kaxy,Kayy, Caxy, Cayy,Maxy,Mayy> along Y direction
C-------------------------------------------------------------------------
          IF (NPOCKET.GT.0) THEN
           DO I=1, NPOCKET
              QZERO(I)=QYZERO(I)
           END DO
          END IF

          FXX(1)=YFXX                         !
          FYY(1)=YFYY                         !
          MXX(1)=YMXX                         !
          MYY(1)=YMYY                         !
C ------------------------------------------- ! --------------------------
          CALL KIJCIJ(FUNRY, L11, QINV, NPOCKET)
C-------------------------------------------- ! --------------------------

          KXYD=K11*CSTIF                      !  FORCE coefficients
          KYYD=K12*CSTIF                      !  Stiffness in [N/m]
          CXYD=C11*CDAMP                      !  Damping in [Ns/m]
          CYYD=C12*CDAMP                      !  due to displacement in Y
          MYYD=zero                           !
          MXYD=zero                           !
c         ....................................!
          KmXYD=KM11*CSTIF*Radius             ! Dimensional MOMENT coefficients
          KmYYD=KM12*CSTIF*Radius             !  Stiffness in [N]
          CmXYD=CM11*CDAMP*Radius             !  Damping in [N-s]
          CmYYD=CM12*CDAMP*Radius             !  due to displacement in Y
          MmYYD=zero                          !
          MmXYD=zero                          !
C                                             !
C        !........................!           !
          IF (INERFREQ.EQ.1) THEN ! --------- ! Inerl=1 & w><0
C        !........................!           !
             IF (NPOCKET.GT.0) THEN           !
              DO I=2, NPOCKET+1               !
                  FXX( I )= FXXO( I )         ! X force components
                  FYY( I )= FYYO( I )         ! Y force components
                  MXX( I )= MXXO( I )         ! X moment components
                  MYY( I )= MYYO( I )         ! Y moment components
                  QZERO(I-1)=QYZEROO(I-1)     ! Flow:  Qzero-Qo1
              END DO                          !
             END IF                           !
              FXX(1)=YFXXO                    !
              FYY(1)=YFYYO                    !
              MXX(1)=YMXXO                    !
              MYY(1)=YMYYO                    !

C            !........................................!
              CALL KIJCIJ(FUNRY,zero, QINVO, NPOCKET)
C            !........................................!

              KdumC=K11*CSTIF                 ! FORCE coefficients
              KdumD=K12*CSTIF                 ! for frequency w=0.
c             CXYDO=C11*CDAMP                 ! due to Y displacement
c             CYYDO=C12*CDAMP                 !
              MXYD=(KdumC-KXYD)/FREQ2         ! Inertia coeffs. in N-s2/m=Kg
              MYYD=(KdumD-KYYD)/FREQ2         !
              KXYD=KdumC                      !
              KYYD=KdumD                      !
c         ....................................!
              KdumC=KM11*CSTIF*Radius         ! MOMENT coefficients
              KdumD=KM12*CSTIF*Radius         ! for frequency w=0.
c             CmXYDO=CM11*CDAMP*Radius        ! due to X displacement
c             CmYYDO=CM12*CDAMP*Radius        !
              MmXYD=(KdumC-KmXYD)/FREQ2       ! Inertia coeffs. in N-s2
              MmYYD=(KdumD-KmYYD)/FREQ2       !
              KmXYD=KdumC                     !
              KmYYD=KdumD                     !

         IF (NPOCKET.GT.0) THEN                !
          DO I=2, NPOCKET+1                    ! ....
              FXX(I)=FXcom(I)                  ! RE-STORE Component Force
              FYY(I)=FYcom(I)                  ! and Moments due to
              MXX(I)=MXcom(I)                  ! PR1=1+i1
              MYY(I)=MYcom(I)
          END DO                               !
         END IF
C        !.......!                            !
          END IF ! -------------------------- ! Inerfreq=1
C        !.......!                            !
C ........................................... ! .............................
C                                             !

 402      CONTINUE

c........................................................................
c FOR TILT-PADS we calculate here the force coefficients for pad rotations
c........................................................................
c    Rsinpk = R x sin(TETApivotk),  Rcospk = R x cos(TETApivotk)

c    !.................................!.......
      IF (TILTPAD.EQ.1) THEN
c    !.................................!.......
      KdXk=Rsinpk*KXXD - Rcospk*KYXD   ! moment stiffness due to X disp [N]
      KdYk=Rsinpk*KXYD - Rcospk*KYYD   ! moment stiffness due to Y disp [N]
      KXdk=RsinpK*KXXD - Rcospk*KXYD   ! force X stiffness due to rotation [N]
      KYdk=Rsinpk*KYXD - Rcospk*KYYD   ! force Y stiffness due to rotation [N]
      Kddk=Rsinpk*KXdk - Rcospk*KYdk   ! moment stiffness due to rotation [N-m]

      CdXk=Rsinpk*CXXD - Rcospk*CYXD   ! moment damping due to X disp [Ns]
      CdYk=Rsinpk*CXYD - Rcospk*CYYD   ! moment damping due to Y disp [Ns]
      CXdk=RsinpK*CXXD - Rcospk*CXYD   ! force X damping due to rotation [Ns]
      CYdk=Rsinpk*CYXD - Rcospk*CYYD   ! force Y damping due to rotation [Ns]
      Cddk=Rsinpk*CXdk - Rcospk*CYdk   ! moment damping due to rotation [N-m-s]
      GOTO 407                         !
c    !.................................!.......
      END IF
c    !.................................!.......


C-------------------------------------------------------------------------
C        FORCE AND MOMENT COEFFICIENTS FOR JOURNAL ANGULATIONS
C-------------------------------------------------------------------------
c       !.....................................!    no moment coefficients
          IF (IMOMFLAG.EQ.0) GOTO 407         !
          IF (MODEL.EQ.1) GOTO 407            ! => for double row HJB
c       !.....................................!    no moment coefficients

          CALL MODGS(YO)                      ! modify GUh & GVh coeffs.

          IFM=1
          CALL PERTXY(YO,IFM,DEVICE,inerfreq) ! Perturbs in X&Y direction
          YREC=zero                           ! Recess at center !######
C                                             !
C.............................................! .........................

C-------------------------------------------------------------------------
C Calculate Coefficients: Kxdx, Kydx, Cxdx, Cydx, Mxdx, Mydx > PERTURB angle
C                       Kaxdx,Kaydx, Caxdx, Caydx,Maxdx,Maydx> about -X dir.
C-------------------------------------------------------------------------
          IF (NPOCKET.GT.0) THEN              !
          DO I=1, NPOCKET                     !
              QZERO(I)=QXZERO(I)
          END DO                              !
          END IF
           FXX(1)=XFXX                        !
           FYY(1)=XFYY                        !
           MXX(1)=XMXX                        !
           MYY(1)=XMYY                        !

C ------------------------------------------- ! --------------------------
          CALL KIJCIJ(FUNRY,(YREC-YO)*L11, QINV, NPOCKET)
C-------------------------------------------- ! --------------------------
          CSTIF=CSTIF*Radius
          CDAMP=CDAMP*Radius

          KXXA=-K11*CSTIF                     !  FORCE coefficients
          KYXA=-K12*CSTIF                     !  Stiffness in [N]
          CXXA=-C11*CDAMP                     !  Damping in [N-s]
          CYXA=-C12*CDAMP                     !  due to -X rotation
          MXXA=zero                           !
          MYXA=zero                           !
c        !....................................!
          KmXXA=-KM11*CSTIF*Radius            !  MOMENT coefficients
          KmYXA=-KM12*CSTIF*Radius            !  Stiffness in [N-m]
          CmXXA=-CM11*CDAMP*Radius            !  Damping in [N-s-m]
          CmYXA=-CM12*CDAMP*Radius            !  due to -X rotation
          MmXXA=zero                          !
          MmYXA=zero                          !
c        !....................................!

C                                             !
C        !........................            !
          IF (INERFREQ.EQ.1) THEN             ! at frequency=0.0
C        !........................            !
          IF (NPOCKET.GT.0) THEN              !
              DO I=2, NPOCKET+1               !
                  FXX( I )= FXXO( I )         !
                  FYY( I )= FYYO( I )         !
                  MXX( I )= MXXO( I )         !
                  MYY( I )= MYYO( I )         !
                  QZERO(I-1)=QXZEROO(I-1)     !
              END DO                          !
          END IF                              !
              FXX(1)=XFXXO                    !
              FYY(1)=XFYYO                    !
              MXX(1)=XMXXO                    !
              MYY(1)=XMYYO                    !

C            !.......................................!
              CALL KIJCIJ(FUNRY,zero, QINVO, NPOCKET)
C            !.......................................!

              KdumD=-K11*CSTIF                ! FORCE coefficients
              KdumC=-K12*CSTIF                ! for frequency w=0.
c             CXXDO=-C11*CDAMP                ! due to -X rotation
c             CYXDO=-C12*CDAMP                !
              MXXA=(KdumD-KXXA)/FREQ2         ! Inertia coefficients in N-s2
              MYXA=(KdumC-KYXA)/FREQ2         !
              KXXA=KdumD                      !
              KYXA=KdumC                      !
c            !................................!
              KdumD=-KM11*CSTIF*Radius        ! MOMENT coefficients
              KdumC=-KM12*CSTIF*Radius        ! for frequency w=0.
c             CmXXDO=-CM11*CDAMP*Radius       ! due to -X rotation
c             CmYXDO=-CM12*CDAMP*Radius       !
              MmXXA=(KdumD -KmXXA)/FREQ2      ! Inertia coefficients in N-s2-m
              MmYXA=(KdumC -KmYXA)/FREQ2      !
              KmXXA=KdumD                     !
              KmYXA=KdumC                     !

          IF (NPOCKET.GT.0) THEN              !
          DO I=2, NPOCKET+1
              FXX(I)=FXcom(I)                  ! RE-STORE Component Force
              FYY(I)=FYcom(I)                  ! and Moments due to
              MXX(I)=MXcom(I)                  ! PR1=1+i1
              MYY(I)=MYcom(I)
          END DO
          END IF
C        !.......!                            !
          END IF ! -------------------------- ! Inerfreq=1
C        !.......!                            !
C
C                                             !
          IF (ICASE.EQ.1) THEN                !
c        !.....................!              !
              KYYA=KXXA                       !
              KXYA=-KYXA                      ! for concentric position
              CYYA=CXXA                       ! Ex=Ey=0.
              CXYA=-CYXA                      !
              MYYA=MXXA                       ! coefficients are skew-symmetric
              MXYA=-MYXA                      !
              KmYYA=KmXXA                     !
              KmXYA=-KmYXA                    ! for concentric position
              CmYYA=CmXXA                     ! Ex=Ey=0.
              CmXYA=-CmYXA                    ! AND YO=0  !#### CHECK
              MmYYA=MmXXA                     !
              MmXYA=-MmYXA                    !
              GOTO 407                        ! => PRINT Coefficients
          END IF                              !
C        !......................!             !
C

C-------------------------------------------------------------------------
C Calculate Coefficients: Kxdy, Kydy, Cxdy, Cydy, Mxdy, Mydy > PERTURB angle
C                       Kaxdy,Kaydy, Caxdy, Caydy,Maxdy,Maydy> about Y dir.
C-------------------------------------------------------------------------

          IF (NPOCKET.GT.0) THEN
          DO I=1, NPOCKET
              QZERO(I)=QYZERO(I)
          END DO
          END IF

          FXX(1)=YFXX                         !
          FYY(1)=YFYY                         !
          MXX(1)=YMXX                         !
          MYY(1)=YMYY                         !

C
C ------------------------------------------- ! --------------------------
          CALL KIJCIJ(FUNRX, (YREC-YO)*L11, QINV, NPOCKET)
C-------------------------------------------- ! --------------------------

          KXYA=K11*CSTIF                      !  FORCE coefficients
          KYYA=K12*CSTIF                      !  Stiffness in [N]
          CXYA=C11*CDAMP                      !  Damping in [N-s]
          CYYA=C12*CDAMP                      !  due to rotation in Y
          MYYA=zero                           !
          MXYA=zero                           !
c         ....................................!
          KmXYA=KM11*CSTIF*Radius             ! Dimensional MOMENT coefficients
          KmYYA=KM12*CSTIF*Radius             !  Stiffness in [N-m]
          CmXYA=CM11*CDAMP*Radius             !  Damping in [N-s-m]
          CmYYA=CM12*CDAMP*Radius             !  due to rotation in Y
          MmYYA=zero                          !
          MmXYA=zero                          !
C                                             !
C        !........................!           !
          IF (INERFREQ.EQ.1) THEN ! --------- ! Inerl=1 & w><0
C        !........................!           !
          IF (NPOCKET.GT.0) THEN              !
              DO I=2, NPOCKET+1               !
                  FXX( I )= FXXO( I )         ! X force components
                  FYY( I )= FYYO( I )         ! Y force components
                  MXX( I )= MXXO( I )         ! X moment components
                  MYY( I )= MYYO( I )         ! Y moment components
                  QZERO(I-1)=QYZEROO(I-1)     ! Flow:  Qzero-Qo1
              END DO                          !
          END IF                              !....
              FXX(1)=YFXXO                    !
              FYY(1)=YFYYO                    !
              MXX(1)=YMXXO                    !
              MYY(1)=YMYYO                    !

C            !........................................!
              CALL KIJCIJ(FUNRX,zero, QINVO, NPOCKET)
C            !........................................!

              KdumC=K11*CSTIF                 ! FORCE coefficients
              KdumD=K12*CSTIF                 ! for frequency w=0.
c             CXYDO=C11*CDAMP                 ! due to Y rotation
c             CYYDO=C12*CDAMP                 !
              MXYA=(KdumC-KXYA)/FREQ2         ! Inertia coeffs. in N-s2
              MYYA=(KdumD-KYYA)/FREQ2         !
              KXYA=KdumC                      !
              KYYA=KdumD                      !
c         ....................................!
              KdumC=KM11*CSTIF*Radius         ! MOMENT coefficients
              KdumD=KM12*CSTIF*Radius         ! for frequency w=0.
c             CmXYDO=CM11*CDAMP*Radius        ! due to X displacement
c             CmYYDO=CM12*CDAMP*Radius        !
              MmXYA=(KdumC-KmXYA)/FREQ2       ! Inertia coeffs. in N-s2-m
              MmYYA=(KdumD-KmYYA)/FREQ2       !
              KmXYA=KdumC                     !
              KmYYA=KdumD                     !
          IF (NPOCKET.GT.0) THEN              !
          DO I=2, NPOCKET+1
              FXX(I)=FXcom(I)                  ! RE-STORE Component Force
              FYY(I)=FYcom(I)                  ! and Moments due to
              MXX(I)=MXcom(I)                  ! PR1=1+i1
              MYY(I)=MYcom(I)
          END DO
          END IF
C        !.......!                            !
          END IF ! -------------------------- ! Inerfreq=1
C        !.......!                            !
C ........................................... ! .............................
C                                             !
 407  CONTINUE
      ISYM=ISYMold
      ICASE=2
C###  CALL PRINTCOEF(DEVICE)                  ! PRINTS Kij, Cij, Mij

C!............................................! OUTPUT FORMATS
C
 1190 FORMAT (' ',3X, 'APPROXIMATE only, Reduced damping factor',
     + ' is lower than 1, WARNING for TOTAL LOSS of DAMPING',/,
     + 4X,50('.'))

 1200 FORMAT (' ', 3X, 'APPROXIMATE: Break frequency=', E11.4E2,
     + 'Hz & Reduced DAMP factor=', E11.4E2,
     + /,10X,'-->  APPROXIMATE ONLY Whirl ratio=', E11.4E2,/,
     + 1X, 39('. '))

 1300 FORMAT (' ', 3X, 'Excitation frequency (w) =', E12.5E2, ' Hz',/,
     + 4X, 'Freq. parameter: w u (R/C)**2/(Ps-Pa) =', E12.5E2, /,
     + 4X, 'Squeeze Reynolds # Res : (p/u) w C**2 =', E12.5E2, /,
     + 1X, 79('.'))

 1400 FORMAT (' ',3X, 'Beta=0:, NO COMPRESSIBILITY EFFECTS.',/,
     +        ' ',79('.'))

 1550 FORMAT (' ',3X, 'Beta=', E12.5E2,'m2/N; Hrec:', E12.5E2, 'm',
     + 1X, ',L2=', E12.5E2,1X, ',L4=',E12.5E2,/,' ',79('.'))

 1600 FORMAT (' ', 3X, 'Hrec=', E12.5E2, 'm  Sigma=', E12.5E2,
     +        '  L2=', E12.5E2, '  L3=', E12.5E2)

 2200 FORMAT (' ', 79(':'))

      END

C
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Pertxy                                                      **
C **                                                                         **
C **  PERTXY:  Finds force, moment & flow components for first order soln.   **
C **           for displacement and angular rotations                        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PERTXY(YO,IFM,DEVICE,INERFREQ)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PERPARM/ FXX(MAXNPOCKP1), FYY(MAXNPOCKP1),
     +                 MXX(MAXNPOCKP1), MYY(MAXNPOCKP1),
     +                 QZERO(MAXNPOCK), QRECREC(MAXNPOCK, MAXNPOCK)
      COMMON /PERPARM0/ FXXO(MAXNPOCKP1), FYYO(MAXNPOCKP1),
     +                  MXXO(MAXNPOCKP1), MYYO(MAXNPOCKP1),
     +                  QZEROO(MAXNPOCK), QRECRECO(MAXNPOCK, MAXNPOCK)
      COMMON /XPERTURB/ XFXX, XFYY, XMXX, XMYY, QXZERO(MAXNPOCK)
      COMMON /YPERTURB/ YFXX, YFYY, YMXX, YMYY, QYZERO(MAXNPOCK)
      COMMON /XPERTURB0/ XFXXO, XFYYO, XMXXO, XMYYO, QXZEROO(MAXNPOCK)
      COMMON /YPERTURB0/ YFXXO, YFYYO, YMXXO, YMYYO, QYZEROO(MAXNPOCK)
      COMMON /PERTREC/ TZERO(MAXNPOCK), TRECREC(MAXNPOCK, MAXNPOCK),
     +                 TZEROO(MAXNPOCK),TRECRECO(MAXNPOCK,MAXNPOCK)
      COMMON /QRECT1/ CRECT(MAXNPOCK),CRECTR(MAXNPOCK)
      COMMON /Sxy/ Sro(Maxnpock),Srec(Maxnpock,Maxnpock)

      COMMON /FLOW01/ QO1(MAXNPOCK), QINO1, QOUTO1
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL

      DOUBLE COMPLEX FXX, FYY, MXX, MYY, QZERO, QRECREC,
     +               FXXO,FYYO,MXXO,MYYO,QZEROO,QRECRECO,
     +               XFXX, XFYY, XMXX, XMYY, QXZERO,
     +               XFXXO,XFYYO,XMXXO,XMYYO,QXZEROO,
     +               YFXX, YFYY, YMXX, YMYY, QYZERO,
     +               YFXXO,YFYYO,YMXXO,YMYYO,QYZEROO,Sro,Srec,
     +               tzero,trecrec,tzeroo,trecreco,CRECT,CRECTR
      DOUBLE PRECISION FREQU, SIGMA, L1, RES, QO1, QINO1, QOUTO1,
     +                 COSXP, SINXP, COSXU, SINXU

      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        IFULL, ICASE, NCASE


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE COMPLEX Isigma
      DOUBLE PRECISION YO
      INTEGER I, J, DEVICE, INERFREQ, IFM, IX, IY

C ----------------------------------------------------------------------------
C --  PERTXY code                                                           --
C ----------------------------------------------------------------------------
C Perturbation along X/Y direction : Displacement (IFM=0) or Angle (IFM=1)
C..........................................................................
      Isigma = DCMPLX(0.0D0, Sigma)

c    !-----------------------!
      IF (IFM.EQ.0) THEN
c    !..............................................! X Displacement, IFM=0
      IX=1
      CALL FXYQXY(COSXP,COSXU,RES,SIGMA,YO,IFM,DEVICE,
     +            inerfreq,icase,IX)
c    !..............................................!
c       ! save flow matrix components due to PRec=1+1i
c       !.............................................!
         IF (NPOCKET.GT.0) THEN
          DO J=1, NPOCKET
            DO I=1, NPOCKET
              QRECREC(I,J)=QRECREC(I,J)+Isigma*Srec(I,J)
     +                    +CRECT(I)*TRECREC(i,j)      ! w/ thermal effects
            END DO
          END DO
         END IF
c       !.............................................!


c    !..............................................! Fp=COS(Xp), Fu=COS(Xu)
      ELSE                                          ! .......................!
c    !..............................................! -X axis Rotation, IFM=1
      IX=2
      CALL FXYQXY(SINXP,SINXU,RES,SIGMA,YO,IFM,DEVICE,
     +            inerfreq,icase,IX)
c    !..............................................! Fp=SIN(Xp), Fu=SIN(Xu)
      END IF
c    !-----------------------!

      IF (NPOCKET.GT.0) THEN
      DO J=1, NPOCKET                               !
         QXZERO(J)=QZERO(J)+DCMPLX(Qo1(J),0.0D0)+Isigma*Sro(j)
     +     +CRECT(J)*TZERO(J)                       ! w/thermal effects
      END DO                                        !
      END IF

      XFXX=FXX(1)
      XFYY=FYY(1)
      XMXX=MXX(1)
      XMYY=MYY(1)

      IF (INERFREQ.EQ.1) THEN                   ! -----------------------
c    !........................!                 !
       IF (NPOCKET.GT.0) THEN
        DO J=1, NPOCKET                         !
         DO I=1, NPOCKET                        ! w/thermal effects
            QRECRECO(I,J)=QRECRECO(I,J)+CRECTR(I)*TRECRECO(I,J)
         END DO
         QXZEROO(J)=QZEROO(J)+DCMPLX(QO1(J),0.0D0)
     +           +CRECTR(J)*TZEROO(J)           ! w/thermal effects
        END DO                                  !
       END IF

          XFXXO=FXXO(1)
          XFYYO=FYYO(1)
          XMXXO=MXXO(1)
          XMYYO=MYYO(1)
      END IF                                    ! INERL=1 & w=0 ------------
c    !........................!                 !


      IF (ICASE.EQ.1) RETURN                    ! Perturb only along X dir


C ................................................. ! ........................
C Perturbation along Y direction: Displacement or Angle
C ............................................................................

c    !-----------------------!
      IF (IFM.eq.0) THEN
c    !..............................................! Y Displacement, IFM=0
      IY=2
      CALL FXYQXY(SINXP,SINXU,RES,SIGMA,YO,IFM,DEVICE,
     +            inerfreq,icase,IY)
c    !..............................................! Fp=SIN(Xp), Fu=SIN(Xu)
      ELSE                                          ! .......................!
c    !..............................................! Y axis Rotation, IFM=1
      IY=2
      CALL FXYQXY(COSXP,COSXU,RES,SIGMA,YO,IFM,DEVICE,
     +            inerfreq,icase,IY)
c    !..............................................! Fp=COS(Xp), Fu=COS(Xu)
      END IF
c    !-----------------------!


      IF (NPOCKET.GT.0) THEN
       DO J=1, NPOCKET
          QYZERO(J)=QZERO(J)+DCMPLX(Qo1(J),0.0D0)+Isigma*Sro(J)
     +       +CRECT(J)*TZERO(J)                     ! w/thermal effects
       END DO
      END IF

      YFXX=FXX(1)
      YFYY=FYY(1)
      YMXX=MXX(1)
      YMYY=MYY(1)

          IF (INERFREQ.EQ.1) THEN                   ! ------------------------
c        !........................!                 !
          IF (NPOCKET.GT.0) THEN
           DO J=1, NPOCKET
              QYZEROO(J)=QZEROO(J)+DCMPLX(QO1(J),0.0D0)
     +           +CRECTR(J)*TZEROO(J)               ! w/thermal effects
           END DO
          END IF

           YFXXO=FXXO(1)
           YFYYO=FYYO(1)
           YMXXO=MXXO(1)
           YMYYO=MYYO(1)

          END IF !.... ............................ ! INER=1 & w=0--------------

C ........................................................................
  100 FORMAT (' ', 3X, '$ -> ',
     +        'FIND: FORCE coefficients for PERTURB in X dir.')

  200 FORMAT (' ', 3X, '$ -> ',
     +        'FIND: FORCE coefficients for PERTURB in Y dir.')

      END

C
C
C *****************************************************************************
C **                                                                         **
C **  Subroutine Fxyqxy                                                      **
C **                                                                         **
C **  FXYQXY:  Main routine for solution algorithm of first order equations. **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FXYQXY(FP,FU,RES,SIGMA,YO,IFM,DEVICE,
     +                  inerfreq,icase,IXY)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /FLOW01/ QO1(MAXNPOCK), QINO1, QOUTO1
      COMMON /So/ Soo(MAXNPOCK)

      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /RECES1/ PREC(MAXNPOCK), TREC(MAXNPOCK),
     +                QREC(MAXNPOCK), QIN, QOUT
      COMMON /S1/ S11(MAXNPOCK)
      COMMON /PERPARM/ FXX(MAXNPOCKP1), FYY(MAXNPOCKP1),
     +                 MXX(MAXNPOCKP1), MYY(MAXNPOCKP1),
     +                 QZERO(MAXNPOCK), QRECREC(MAXNPOCK, MAXNPOCK)
      COMMON /PERPARM0/ FXXO(MAXNPOCKP1), FYYO(MAXNPOCKP1),
     +                  MXXO(MAXNPOCKP1), MYYO(MAXNPOCKP1),
     +                  QZEROO(MAXNPOCK), QRECRECO(MAXNPOCK, MAXNPOCK)
      COMMON /PERTREC/ TZERO(MAXNPOCK), TRECREC(MAXNPOCK, MAXNPOCK),
     +                 TZEROO(MAXNPOCK),TRECRECO(MAXNPOCK,MAXNPOCK)
      COMMON /Sxy/ Sro(MAXNPOCK),Srec(MAXNPOCK,MAXNPOCK)

c.............................................................................
      INTEGER ITER, ITMAX, ITPMAX, ISYM, ICSTEP,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL,
     +        ISOLN, SVERB,DVERB, BEEP, BEARING

      DOUBLE PRECISION Soo, QO1, QINO1, QOUTO1

      DOUBLE COMPLEX U1,V1,P1,T1, S11, PREC,TREC, QREC, QIN, QOUT,
     +               FXX, FYY, MXX, MYY, QZERO, QRECREC,
     +               FXXO,FYYO,MXXO,MYYO,QZEROO,QRECRECO, Sro, Srec,
     +               tzero,trecrec,tzeroo,trecreco

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER I,J,JP1,ITSUM,DEVICE,K,L,INERFREQ,II,KK,Icase,jmax,
     +        Kstart, IFM, IXY
      DOUBLE PRECISION ANGK, COSA, SINA, Res0, Sig0, YO,
     +                 PI, RES, FU(*), FP(*), SIGMA, ANGLE
      DOUBLE COMPLEX QDUMY(MAXNPOCK),TDUMY(MAXNPOCK),
     +               Sdumy(MAXNPOCK), ONE, ZERO

C ----------------------------------------------------------------------------
C --  FXYQXY code                                                           --
C ----------------------------------------------------------------------------
C Basis for the method ON HYDROSTATIC BEARINGS:
C First order flow field is decomposed into the following components, say
C Px = PXo + SUM (PXr pxr), etc,  r=1,2....,NPOCKET
C where Pxo satisfies the non-homogeneous form of the first order eqns.
C       with recess pressures = 0+i0 and hx,hy,hdx,hdy <> 0  (ISOLN=1)
C       Pxr satisfies the homogeneous form of the first order equations
c       with recess pressure pxr=1+i1 and hx=hy=hdx=hdy = 0  (ISOLN=0)
c Thus the rec. pressure components need only be calculated once  !!!!
C ----------------------------------------------------------------------------
c
      PI=DACOS(-1.0D0)
      ZERO=(0.0D0, 0.0D0)
c     .................................!
      IF (Sigma.EQ.(0.0D0)) THEN
         ONE=(1.0D0,0.0D0)
      ELSE
         ONE=(1.0D0,1.0D0)
      END IF
c     .................................!
      Res0=0.0D0
      Sig0=0.0D0
c     .................................!
      Kstart=ISYM+(ISYM-1)*NYI
      IF (BEARING.EQ.2) Kstart=1
c     .................................!
      ITSUM=0                                       !
C ................................................. ! ........................

      CALL FLOWO1(FU, FP, YO, IFM)                  ! Calc Pert flow from 0 Soln
C                                                   ! INT(RHOo H1 Vo dTau)
C ................................................. ! ........................

      IF (NPOCKET.GT.0) THEN
      DO I=1, NPOCKET                               ! IF Hydrostatic bearing:
          PREC(I)=ZERO                              ! SET ZERO for
          TREC(I)=ZERO                              ! first order recess pressures
      END DO                                        ! and temperatures
      END IF                                        !............................

      DO I=1, NXT                                   !
          DO J=Kstart, NYI                          !
              V1(I, J)=ZERO                         !
              U1(I, J)=ZERO                         !
              P1(I, J)=ZERO                         !
              T1(I, J)=ZERO                         !
          END DO                                    !
      END DO                                        !
c    ...............................................!

c   ------------                                    !.........................
      ISOLN=1                                       ! First component, freq<>0
c   ------------                                    !.........................
c   * Calculate flow on film lands and impedances with Precj = (0.0, 0.0)
      CALL SLAND1(FU,FP,YO,RES,SIGMA,DEVICE,IFM)
      CALL FORCE1(FXX(1),FYY(1),MXX(1),MYY(1),YO,IFM)
C     ..............................................!
      ITSUM=ITSUM+ITER                              !

c    !.......................!
      IF (NPOCKET.GT.0) THEN
c    !.......................!
        IF ((ISYM.eq.1).and.(IFM.eq.1)) THEN        ! X or Y axis rotation
          DO I=1, NPOCKET                           ! in a symmetric bearing
           QREC(I)=zero
           TREC(I)=zero
           S11(I)=zero
          END DO
        END IF

        DO I=1, NPOCKET                             ! and save flows due only
          QZERO(I)=QREC(I)                          !   to zero order component
          TZERO(I)=TREC(I)                          !
          Sro(I)=S11(i)+DCMPLX(Soo(I),0.0D0)        !
        END DO                                      !

c    !.......................!
      END IF
c    !.......................!

C     ..............................................!
      IF (SVERB.ge.1) THEN                          !
          WRITE (6, 300) 0, FXX(1), FYY(1)          !
      END IF                                        !
      IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN      !
       WRITE (1, 300) 0, FXX(1), FYY(1)             !
      END IF                                        !
C     ..............................................!



c    !........................!                     !
      IF (INERFREQ.EQ.1) THEN ! ------------------- ! Inerfreq=1
c    !........................!                     !

          DO II=1, NXT                              !
              DO KK=Kstart, NYI                     ! extract REAL part of
                  V1(II, KK)=DREAL(V1(ii,kk))       ! perturbed fields and
                  U1(II, KK)=DREAL(U1(ii,kk))       ! Solve for w=0. (Sigma=Res=0)
                  P1(II, KK)=DREAL(P1(ii,kk))       !
                  T1(II, KK)=DREAL(T1(ii, kk))      !
              END DO                                ! STATIC DISP & FORCES
          END DO                                    !

          IF (NPOCKET.GT.0) THEN                    !
            DO I=1, NPOCKET                         !
               TREC(I)=DREAL(TREC(I))               ! First order recess temps.
            END DO                                  ! Precj=0 ==> 0th components
c                                                   !
            CALL SETPREC1                           !
          END IF                                    !
c        ...........................................!
          CALL SLAND1(FU, FP,YO,Res0,Sig0,DEVICE,IFM)
          CALL FORCE1(FXXO(1),FYYO(1),MXXO(1),MYYO(1),YO,IFM)
C        ...........................................!
          ITSUM=ITSUM+ITER                          !

          IF (NPOCKET.GT.0) THEN
           IF ((ISYM.eq.1).and.(IFM.eq.1)) THEN
              DO I=1, Npocket
                QREC(I)=zero
                TREC(I)=zero
                S11(I)=zero
              END DO
           END IF
           DO I=1, NPOCKET
              QZEROO(I)=QREC(I)
              TZEROO(I)=TREC(I)
           END DO
          END IF

C     ..............................................!

C         ..........................................!
          IF (SVERB.ge.1) THEN                      !
              WRITE (6, 300) 0, FXXO(1), FYYO(1)    !
          END IF                                    !
          IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN  !
              WRITE (1, 300) 0, FXXO(1), FYYO(1)    !
          END IF                                    !
      END IF ! ------------------------------------ ! Inerfreq=1



C ................................................. ! ........................
c                                                   ! for seals & hydrodynamic
      IF ((IXY.eq.2).OR.(NPOCKET.EQ.0)) GOTO 705    ! bearings => RETURN
c                                                   !
c ..................................................! ........................
c       FOR HYDROSTATIC BEARINGS:
c 	Component solutions due to Recess Pressure need to be calculated
c 	ONLY ONCE for the first-order solution.
c 	These calculations are done IF: IXY=1 and IFM=0, i.e.
c 	for the perturbation in X displacement
c ..................................................! ........................

c    ---------                                      !
      ISOLN=0                                       ! non-homog components
c    ---------                                      ! will set PREC=1+i1

      IF (ICASE.EQ.1) THEN                          !....
        jmax=1                                      ! centered HJB symmetric
      ELSE                                          !.....
        jmax=NPOCKET                                ! off centered solution
      END IF                                        !....


c ::::::::::::::::::::::::::::::::::::::::::::::::::! Calcule flow & forces
C                                                   ! Find Uj,Pj,Vj
C   !::::::::::::::!                                !--------------------
      DO j=1, jmax                                  !......................
c   !:::::::::::::::!                               !
                                                    !
          DO I=1, NPOCKET                           !
              PREC(I)=ZERO                          !    On each recess:
              TREC(I)=ZERO                          !
          END DO                                    !    find solution to linear

          DO I=1, NXT                               !    system of equations
              DO L=Kstart, NYI                      !    with Precj = 1 + i1
                  V1(I, L)=ZERO                     !    for Isoln=0
                  U1(I, L)=ZERO                     !
                  P1(I, L)=ZERO                     !
                  T1(I, L)=ZERO                     !
              END DO                                !
          END DO                                    !
c        !!!!!!!!!!!!!!!!                           !
          PREC(J)=ONE                               !-> nonhomog recess press
c        !!!!!!!!!!!!!!!!                           !
c     *   set boundary conditions and solve on film lands
c         FU=FP=0 on first order equations

          CALL SETPREC1                             !
          CALL SLAND1(FU,FP,YO,RES,Sigma,DEVICE,IFM)
C                                                   !
C        ...........................................!

          DO I=1, NPOCKET                           ! Determines coeffs of matrix
              QRECREC(I, J)=QREC(I)/ONE             ! of flow rates due to Prec
              TRECREC(I, J)=TREC(I)/ONE             !
              Srec(i,j)=S11(I)/ONE                  !
          END DO                                    !

          JP1=J+1                                   !

      CALL FORCE1(FXX(JP1),FYY(JP1),MXX(JP1),MYY(JP1),YO,IFM)
C     ..............................................! Calculate Force Compons.
          FXX(JP1)=FXX(JP1)/ONE                     !
          FYY(JP1)=FYY(JP1)/ONE                     !
          MXX(JP1)=MXX(JP1)/ONE                     !
          MYY(JP1)=MYY(JP1)/ONE                     !
C         ..........................................!
          IF (SVERB.ge.1) THEN                      !
              WRITE (6, 300) J, FXX(JP1), FYY(JP1)  !
          END IF                                    !
          IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN  !
              WRITE (1, 300) J, FXX(JP1), FYY(JP1)  !
          END IF                                    !
C         ..........................................!

          ITSUM=ITSUM+ITER                          !

C        !........................!                 ! Inerfreq=1
          IF (INERFREQ.EQ.1) THEN ! --------------- !
C        !........................!                 ! Inerfreq=1
              PREC(J)=(1.0D0, 0.0D0)                ! w=0, Sigma=Res=0.0
              DO II=1, NXT                          !
                  DO KK=Kstart, NYI                 ! get REAL part of fields
                     V1(II, KK)=DREAL(V1(ii, kk)/ONE)
                     U1(II, KK)=DREAL(U1(ii, kk)/ONE)
                     P1(II, KK)=DREAL(P1(ii, kk)/ONE)
                     T1(II, KK)=DREAL(T1(ii, kk)/ONE)
                  END DO                            !
              END DO                                !

              DO I=1, NPOCKET                       !
                 TREC(I)=DREAL(TREC(I)/ONE)         !
              END DO                                !

              CALL SETPREC1                         !
              CALL SLAND1(FU,FP,YO,Res0,Sig0,DEVICE,IFM)
C             ......................................!

              DO I=1, NPOCKET                       !
                  QRECRECO(I, J)=QREC(I)            !
                  TRECRECO(I, J)=TREC(I)            !
              END DO                                !

              JP1=J+1                               !

      CALL FORCE1(FXXO(JP1),FYYO(JP1),MXXO(JP1),MYYO(JP1),YO,IFM)
C     ..............................................!

C             ......................................!
              IF (SVERB.ge.1) THEN                  !
                  WRITE (6, 300) J, FXXO(JP1), FYYO(JP1)
              END IF                                !
              IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN
                  WRITE (1, 300) J, FXXO(JP1), FYYO(JP1)
              END IF                                !
C             ......................................!

              ITSUM=ITSUM+ITER                      !

          END IF ! -------------------------------- ! Inerfreq=1
c       !........!
C
C !::::::::::::::!                                  !--------------------
      END DO                                        ! j=1,2,.. jmax
C !::::::::::::::!                                  !--------------------


      IF (ICASE.EQ.2) GOTO 705    !=> offcentered HJB

c IN a rotationally symmetric HJB with its journal at centered
c position, we need perturb only ONE recess and
c calculate the first order fields. Fields for other perturbed recesses
c are identical except that are shifted an angle = to angle between
c recesseses. Very nice idea !
c  !................................................!ROTATE flows & forces
c                                                   !for concentric case

          ANGLE=2.0D0*PI/NPOCKET                    ! Angle between recesses
          ANGK=ANGLE                                !

          DO J=2, NPOCKET                           ! Find forces in other reces
              COSA=DCOS(ANGLE)                      ! since Pk = Pj(0-Angle)
              SINA=DSIN(ANGLE)                      !   by rotation
              FXX(J+1)=FXX(2)*COSA-FYY(2)*SINA      !
              FYY(J+1)=FXX(2)*SINA+FYY(2)*COSA      !
              MXX(J+1)=MXX(2)*COSA-MYY(2)*SINA      !
              MYY(J+1)=MXX(2)*SINA+MYY(2)*COSA      !

              IF (INERFREQ.EQ.1) THEN ! ----------- !
                  FXXO(J+1)=FXXO(2)*COSA-FYYO(2)*SINA
                  FYYO(J+1)=FXXO(2)*SINA+FYYO(2)*COSA
                  MXXO(J+1)=MXXO(2)*COSA-MYYO(2)*SINA
                  MYYO(J+1)=MXXO(2)*SINA+MYYO(2)*COSA
              END IF ! ---------------------------- ! Inerfreq=1

              ANGLE=ANGLE+ANGK                      !

          END DO                                    ! ........................
                                                    !
          DO J=2, NPOCKET                           ! Switches flow rates due to
              DO I=1, NPOCKET                       !   Precj
                  QDUMY(I)=QRECREC(I, J-1)          ! ........................
                  TDUMY(I)=TRECREC(I, J-1)          !
                  Sdumy(I)=Srec(i,j-1)              !
              END DO                                !
              K=NPOCKET                             !
              DO I=1, NPOCKET                       !
                  QRECREC(I, J)=QDUMY(K)            !
                  TRECREC(I, J)=TDUMY(K)            !
                  Srec(i,j)=Sdumy(k)                !
                  K=I                               !
              END DO                                !
          END DO                                    !
                                                    !
          IF (INERFREQ.EQ.1) THEN ! --------------- !
              DO J=2, NPOCKET                       !
                  DO I=1, NPOCKET                   !
                      QDUMY(I)=QRECRECO(I, J-1)     !
                      TDUMY(I)=TRECRECO(I, J-1)     !
                  END DO                            !
                  K=NPOCKET                         !
                  DO I=1, NPOCKET                   !
                      QRECRECO(I, J)=QDUMY(K)       !
                      TRECRECO(I, J)=TDUMY(K)       !
                      K=I                           !
                  END DO                            !
              END DO                                !
          END IF ! -------------------------------- ! End of concentric sol
C                                                   ! ::::::::::::::::::::::::


 705   CONTINUE
C##    WRITE (6, 200) ITSUM
       IF (DEVICE.EQ.1) WRITE (1, 200) ITSUM
       ITER=ITSUM                                   ! = Total # of land steps
C
C-----------------------------------------------------------------------

  200 FORMAT (' ',  '$ -> ',
     +        'FIND: FORCE coefficients for X/Y dir.',
     +        ' Perturb. in ', I3, ' Steps')

  300 FORMAT (' ', 3X, 'j:', I2, 'FX_:', 2(E12.5E2, 1X), 2X,
     +        'FY_:', 2(E12.5E2, 1X))

  400 FORMAT (' ', 3X, 'FOR CONCENTRIC SOLUTION: ',
     +        'Rotate forces & Switch flow rates.')

      END

C

C *****************************************************************************
C **                                                                         **
C **  Subroutine Kijcij                                                      **
C **                                                                         **
C **  KIJCIJ:  extracts dimensionless stiffness and damping coefficients.    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE KIJCIJ(FR, L1, QINV, NR)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PERPARM/ FXX(MAXNPOCKP1), FYY(MAXNPOCKP1),
     +                 MXX(MAXNPOCKP1), MYY(MAXNPOCKP1),
     +                 QZERO(MAXNPOCK), QRECREC(MAXNPOCK, MAXNPOCK)
      COMMON /COEFS/ K11, K12, C11, C12, KM11,KM12,CM11,CM12

      DOUBLE PRECISION K11, K12, C11, C12, KM11,KM12,CM11,CM12
      DOUBLE COMPLEX FXX, FYY, MXX, MYY, QZERO, QRECREC
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE COMPLEX QINV(MAXNPOCK, MAXNPOCK),RHS(MAXNPOCK),
     +               FX, FY, MX, MY, PREC1
      DOUBLE PRECISION FR(*), L1, DZERO
      INTEGER I, J, NR
C ----------------------------------------------------------------------------
C --  KIJCIJ code                                                           --
C ----------------------------------------------------------------------------
c NOTE:  QINV: inverse matrix of flows calculated on sub COEFFIC
c =====
      FX=FXX(1)
      FY=FYY(1)
      MX=MXX(1)
      MY=MYY(1)

      IF (NR.EQ.0) GOTO 77                         !==> no pockets on bearing/pad
      DZERO=0.0D0                                  !
      DO J=1, NR                                   ! Nr=Npocket
          RHS(J)=-QZERO(J)-DCMPLX(DZERO, L1*FR(J)) ! Qzero=Qzero+Qo1
      END DO                                       ! ...........................

c    !..................!                          !
      DO I=1, NR                                   ! ...........................
c    !..................!                          !

          PREC1=(0.0D0,0.0D0)                      ! Calculate first order
          DO J=1, NR                               ! recess pressures
              PREC1=PREC1+QINV(I, J)*RHS(J)        !
          END DO                                   ! ...........................

          FX=FX+FXX(I+1)*PREC1                     ! Calculate dimensionless
          FY=FY+FYY(I+1)*PREC1                     ! Stiffness and damping
          MX=MX+MXX(I+1)*PREC1
          MY=MY+MYY(I+1)*PREC1

C          write (6,122) I,Prec1,FXX(I+1),FYY(I+1),
C     +                          MXX(I+1),MYY(I+1)

c    !..................!                          !
      END DO                                       ! I=1,NPOCKET
c    !..................!                          !

      !write (6,121) FX,FY,MX,MY

 77   K11=-DREAL(FX)                               ! = Kxx or Kxy FORCES
      K12=-DREAL(FY)                               ! = Kyx or Kyy
      C11=-DIMAG(FX)                               ! = Cxx or Cxy
      C12=-DIMAG(FY)                               ! = Cyx or Cyy
      KM11=-DREAL(MX)                              ! = Kxx or Kxy MOMENTS
      KM12=-DREAL(MY)                              ! = Kyx or Kyy
      CM11=-DIMAG(MX)                              ! = Cxx or Cxy
      CM12=-DIMAG(MY)                              ! = Cyx or Cyy
      RETURN

 120  format ('I, Prec, Fx,Fy,Mx,My:')
 121  format (' ',27X,4(E11.4E2,',',E11.4E2,1X))
 122  format (' ',I2,5(E11.4E2,',',E11.4E2,1X))
      END

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c analysis and program by Luis San Andres
c last updated on 9/13/95
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
