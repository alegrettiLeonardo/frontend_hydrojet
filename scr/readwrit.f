c 7/7/95 updated for angle injection in HJBS :
c PATCH used for compatibility with OLD HYDRO versions:
c PRjet(MAXNPOCK,MAXNPOCK+2) superseeds Prise(MAXNPOCK,-MAXNYI:MAXNYI)
c where PRjet(I, j=1,NPC) is circumf. pressure profile within recess.
c NOTE: PRECup(I)=PRjet(I,1), PRECdo=PRjet(I,NPC); I=1,NPOCKET
c
c ------------------------------------------------
c 6/20/95 updated for HCELL depth  READ/WRITE
c ------------------------------------------------
c 6/02/94 updated for radial heat flow parameters
c ------------------------------------------------
c
c #####   ######    ##    #####   #    #  #####      #     #####  ######   #####
c #    #  #        #  #   #    #  #    #  #    #     #       #    #          #
c #    #  #####   #    #  #    #  #    #  #    #     #       #    #####      #
c #####   #       ######  #    #  # ## #  #####      #       #    #          #
c #   #   #       #    #  #    #  ##  ##  #   #      #       #    #          #
c #    #  ######  #    #  #####   #    #  #    #     #       #    ######     #
C
C #    #   #   #  #####   #####    ####        #  ######   #####
C #    #    # #   #    #  #    #  #    #       #  #          #
C ######     #    #    #  #    #  #    #       #  #####      #
C #    #     #    #    #  #####   #    #       #  #          #
C #    #     #    #    #  #   #   #    #  #    #  #          #
C #    #     #    #####   #    #   ####    ####   ######     #

C  hydrojet.f Copyright Dr. Luis SanAndres TexasA&MUniversity / 1995
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR III
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C *****************************************************************************
C **                                                                         **
C **  Subroutine Readata                                                     **
C **                                                                         **
C **  READATA:  Reads from data file :  input data and U, V, P arrays.       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE READATA(FILE, DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /PEDGE/ Prise(MAXNPOCK,-MAXNYI:MAXNYI)
      COMMON /SPLDATA/Z(Nsl),Cl(Nsl),BCL(Nsl),CCL(Nsl),DCL(Nsl),Nj
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD,AR,  HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR

      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PADORIF(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)

      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ Ewx, Ewy, Ewear, Betaw, Iwear
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr
      COMMON /LOSPAD/ LOSleadP, KLOSpad
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FORCE0/ FFACTOR, FXO, FYO, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /TGROOVE/ LEMDA,DELTA

      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /THERMAL/ ALFT,UC,TC,EC
      COMMON /OILCOEF/ Talpha
      COMMON /THERMID/ ISOTH
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /LIQUID/ Tempk,Vsound, IF, IL
      COMMON /CONT/ IFF
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /GUESPI/ IGUESP
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)

c     .............................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION TORR,U, V, P, T,ALFT,UC,TC,EC,CPS,THS,BETAKS,
     +                 PREC, QREC,TREC, PRISE, QIN, QOUT, QFACTOR,
     +                 Z,Cl,Bcl,Ccl,Dcl
      DOUBLE PRECISION X1, X1r, LREC, LPAD,LENGTHR,LENGTHL,
     +                 PADORIF, ROTPAD, INERPAD,KSTPAD,CDAPAD,
     +                 LEMDA, DELTA
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,
     +                 EXO,EYO,AXO,AYO,ZO,Ewx, Ewy, Ewear, Betaw,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr
      DOUBLE PRECISION cinlet, cexit,ClearO,ClearR,ClearL,YR,YL,
     +                 Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 PRCOEFC,PRCOEFS,PLCOEFC, PLCOEFS
      DOUBLE PRECISION FREQU, SIGMA, L1, RES,HRECU, VSUP, BETA,
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD,KROTPAD,CROTPAD
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,
     +                 ALFP,FFACTOR, FXO, FYO, TO, TOR,MFACTOR,MX,MY,
     +                 PRATIO, CORIF,SMASS,MPEPS, PREPS, MMP, SFLOW,
     +                 AMOD, BMOD, RUGR, RUGS, EXPO, TEMPK, VSOUND,
     +                 LOSLEADP, KLOSPAD,AC, ETA, RELAXH,
     +                 PRELOAD, PRELOADB, OFFSET, ROTDEL,Xpivotk,
     +                 QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL ,Talpha,
     +                 TSHAFT, TSTATOR,
     +                 TBOUT, THERMALK, ROUTER, HKB,HCELL,HCDIM,
     +                 ANGLEJ, LOCJET, CJET, DPJET,
     +                 PRECdo, PRECup, PRjet

      INTEGER ICASE, NCASE, Iwear,IFF, IF, IL, IFULL, NPAD, NREC,
     +        Nj, ITER, ITMAX, ITPMAX, ISYM, ICSTEP, ISOTH, IGUESP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,LIFT,
     +        NPOCKET, NLC, NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,TILT,
     +        IPRuni, IPLuni, NPRcs, NPLcs,BEARING, TILTPAD,NLOBES
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*40 FILE
      CHARACTER*2 VALK
      CHARACTER*7 SEED
      CHARACTER*9 TEMPFILE
      INTEGER DEVICE, I, J, IJ, IOS, K, JSTART

C ----------------------------------------------------------------------------
C --  READATA code                                                          --
C ----------------------------------------------------------------------------
      SEED='TEMPPAD'

c    !---------------------------------------------! OPEN DATA FILE
      OPEN (UNIT=2, FILE=FILE, STATUS='OLD', IOSTAT=IOS, ERR=1000)
c    !---------------------------------------------! OPEN DATA FILE

      READ (2, 20, IOSTAT=IOS, ERR=1100) TITLE
   20 FORMAT (60A)
      READ (2, 30, IOSTAT=IOS, ERR=1100) DDATE
   30 FORMAT (10A)
      READ (2,*,IOSTAT=IOS,ERR=1100) NPAD,NPOCKET,NLC,NPC,NLA,NPA
      READ (2,*,IOSTAT=IOS,ERR=1100) NJ,INERL,INERP,ITURB,ICAV,
     +                               ISOTH,MODEL
      READ (2,*,IOSTAT=IOS,ERR=1100) IF,IL,IWEAR,ISYM,ICSTEP,NLOBES
      READ (2,*,IOSTAT=IOS,ERR=1100) ITMAX, ITPMAX, IFULL, BEARING,
     +                               TILTPAD
      TILT=TILTPAD

C--------------------------------------------------------
C CHECK:  to make sure...
c--------------------------------------------------------
      NXI=NLC+NPC-2
      NXT=NXI*NPOCKET+NLC
      IF (IFULL.EQ.1) NXT=NXI*NPOCKET+1
      IF (NPOCKET.EQ.0) NXT=NLC
      NYI=NLA+NPA
      IF (BEARING.GE.2) NYI=NLA
C
      IF (NJ.GT.NSL) THEN
         write (6,*) 'This input file violates the upper limit for NSL,
     + MAX number of axial coordinates. TO raise this limt, edit file
     + params.f and recompile the entire program'
       CALL BEEPER
       CLOSE (UNIT=2, IOSTAT=IOS, ERR=1300)
       CALL PAUSE
       RETURN

      ELSE IF (NPAD.GT.MAXNPAD) THEN
         write (6,*)'This input file violates the upper limit for NPAD,
     + MAX number of PADS. TO raise this limt, edit file
     + params.f and recompile the entire program'
       CALL BEEPER
       CLOSE (UNIT=2, IOSTAT=IOS, ERR=1300)
       CALL PAUSE
       RETURN

      ELSE IF (NPOCKET.GT.MAXNPOCK) THEN
          WRITE (6, *) 'This input file violates the upper limit for NPO
     +CKET:'
          WRITE (6, 75) NPOCKET
   75     FORMAT (' ', 'NPOCKET=', I3)
          WRITE (6, 76) MAXNPOCK
   76     FORMAT (' ', 'The maximum NPOCKET allowed in this version of h
     +ydroflex is ', I3, '.')
          WRITE (6, *) 'To raise this limit, edit file params.f   and re
     +compile the entire program.'
          WRITE (6, *) 'Aborting read.  The previous data has been corru
     +pted by this error.'
          CALL BEEPER
          CLOSE (UNIT=2, IOSTAT=IOS, ERR=1300)
          CALL PAUSE
          RETURN

      ELSE IF (NXT.GT.MAXNXT) THEN
          WRITE (6, *) 'This input file violates the upper limit for NXT
     +:'
          WRITE (6, 71) NXT
   71     FORMAT (' ', 'NXT=(NPC+NLC-2)*NPOCKET+NLC=', I3)
          WRITE (6, 72) MAXNXT
   72     FORMAT (' ', 'The maximum NXT allowed in this version of hydro
     +seal is ', I3, '.')
          WRITE (6, *) 'To raise this limit, edit file params.f   and re
     +compile the entire program.'
          WRITE (6, *) 'Aborting read.  The previous data has been corru
     +pted by this error.'
          CALL BEEPER
          CLOSE (UNIT=2, IOSTAT=IOS, ERR=1300)
          CALL PAUSE
          RETURN

      ELSE IF (NYI.GT.MAXNYI) THEN
          WRITE (6, *) 'This input file violates the upper limit for NYI
     +:'
          WRITE (6, 73) NYI
   73     FORMAT (' ', 'NYI=NLA+NPA=', I3)
          WRITE (6, 74) MAXNYI
   74     FORMAT (' ', 'The maximum NYI allowed in this version on hydro
     +seal is ', I3, '.')
          WRITE (6, *) 'To raise this limit, edit file params.f   and re
     +compile the entire program.'
          WRITE (6, *) 'Aborting read.  The previous data has been corru
     +pted by this error.'
          CALL BEEPER
          CLOSE (UNIT=2, IOSTAT=IOS, ERR=1300)
          CALL PAUSE
          RETURN
      END IF
C
C--------------------------------------------------------
C END OF CHECK:  to make sure...
c--------------------------------------------------------

      READ (2,*,IOSTAT=IOS,ERR=1100) (NREC(K), K=1, NPAD)
      READ (2,*,IOSTAT=IOS,ERR=1100) (X1(K),K=1,NPAD),   ! PAD: start angle
     +                             (LPAD(K),K=1,NPAD)    !     angular length
      READ (2,*,IOSTAT=IOS,ERR=1100) (ROTPAD(K),K=1,NPAD)!     rotations

      IF (TILT.EQ.1) THEN
      READ (2,*,IOSTAT=IOS,ERR=1100) (INERPAD(K),K=1,NPAD)! pad Inertia
      READ (2,*,IOSTAT=IOS,ERR=1100) (KSTPAD(K),K=1,NPAD) ! pad ROT Stiffness
      READ (2,*,IOSTAT=IOS,ERR=1100) (CDAPAD(K),K=1,NPAD) ! pad ROT Damping
      END IF

      DO K=1, NPAD
        NPOCKET=NREC(K)
        IF (NPOCKET.GE.1) THEN                                   ! PAD:
          READ (2,*,IOSTAT=IOS,ERR=1100) (X1r(J,K),J=1,NPOCKET), ! Start Angle
     +                                  (Lrec(J,K),J=1,NPOCKET), ! Angular Length
     +                               (PaDorif(J,K),J=1,NPOCKET)  ! Orifice Diams.
        END IF
      END DO

      READ (2, *, IOSTAT=IOS, ERR=1100) ALFU,ALFP,ALFT,SFLOW,SMASS,MMP
      READ (2, *, IOSTAT=IOS, ERR=1100) LOSXSIxu, LOSXSIxd, LOSXSIyl,
     +                                  LOSXSIyr, ALPHA, PRATIO
      READ (2, *, IOSTAT=IOS, ERR=1100) CLEAR, CINLET, CEXIT
      READ (2, *, IOSTAT=IOS, ERR=1100) (Z(J), J=1, NJ)
      READ (2, *, IOSTAT=IOS, ERR=1100) (CL(J), J=1, NJ)
      IF (ICSTEP.eq.1) THEN
        READ(2,*, IOSTAT=IOS, ERR=1100) ClearO,ClearL,ClearR,YR,YL
      END IF

      READ (2, *, IOSTAT=IOS, ERR=1100) DIAM, LENGTH,AR, HREC, VSUP,
     +                                  LENGTHR
      LENGTHL=LENGTH-LENGTHR
      READ (2, *, IOSTAT=IOS, ERR=1100) EMUS,RHOS,CPS,THS,BETAKS,
     +                                  EMUle,RHOle,EMUri,RHOri
      READ (2, *, IOSTAT=IOS, ERR=1100) P1PROP, EMU1,RHO1,
     +                                  P2PROP, EMU2,RHO2
      READ (2, *, IOSTAT=IOS, ERR=1100) RPM,FREQU,PS,PLeft,Pright,PC
      READ (2, *, IOSTAT=IOS, ERR=1100) Cleft,Cright,CD,DORIF,BETA
      READ (2, *, IOSTAT=IOS, ERR=1100) RUGR, RUGS, AMOD, BMOD, EXPO
      READ (2, *, IOSTAT=IOS, ERR=1100) EXO, EYO, AXO, AYO, ZO
      READ (2, *, IOSTAT=IOS, ERR=1100) TEMPK, VSOUND, LOSLEADP
      READ (2, *, IOSTAT=IOS, ERR=1100) PRELOAD, OFFSET
      READ (2, *, IOSTAT=IOS, ERR=1100) EWX, EWY, EWEAR, BETAW

cc    OIL Temperature-viscosity coeff.
      READ (2, *, IOSTAT=IOS, ERR=1100) Talpha

      PRELOADB=PRELOAD
      NXI=NLC+NPC-2

      CALL ZERO0
      CALL ZERO1

      IFF=IF
c     ----------------------
      IF (IF.LE.4) THEN    !
         CALL FDATA(IF)    ! => initialize mipropst.f
      END IF               !
c     ----------------------


C     ................!..............................................!
      DO K=1, NPAD
C     ................! re-store data for all bearing pads
      TEMPFILE=SEED//VALK(K)
c
      NPOCKET=NREC(K)
      NXT=NXI*NPOCKET+NLC
      IF (IFULL.EQ.1) NXT=NXI*NPOCKET+1
      IF (NPOCKET.EQ.0) NXT=NLC

c     ...........................................................
      OPEN (UNIT=42, FILE=TEMPFILE, STATUS='UNKNOWN', IOSTAT=IOS,
     +      ERR=1099)
      WRITE (6, 77) K, TEMPFILE, FILE
c     ...........................................................

      WRITE (42, *, IOSTAT=IOS, ERR=1199) NPOCKET, NXT, TILT
c ... READ PAD Forces, Torque, Moments and flows
      READ (2, *, IOSTAT=IOS, ERR=1100) TOR,FXO,FYO,MX,MY,QIN,QOUT
      WRITE(42,*, IOSTAT=IOS, ERR=1199) TOR,FXO,FYO,MX,MY,QIN,QOUT
c ... READ RECESS pressures, temperatures, flowrates, Pressure rise, torque
      IF (NPOCKET.GT.0) THEN
        READ (2, *, IOSTAT=IOS, ERR=1100) (PREC(I), I=1, NPOCKET)
        WRITE(42,*, IOSTAT=IOS, ERR=1199) (PREC(I), I=1, NPOCKET)
        READ (2, *, IOSTAT=IOS, ERR=1100) (QREC(I), I=1, NPOCKET)
        WRITE(42,*, IOSTAT=IOS, ERR=1199) (QREC(I), I=1, NPOCKET)
        READ (2, *, IOSTAT=IOS, ERR=1100) (TREC(I), I=1, NPOCKET)
        WRITE(42,*, IOSTAT=IOS, ERR=1199) (TREC(I), I=1, NPOCKET)
        READ (2, *, IOSTAT=IOS, ERR=1100) (TORR(I), I=1, NPOCKET)
        WRITE(42,*, IOSTAT=IOS, ERR=1199) (TORR(I), I=1, NPOCKET)

        DO I=1, NPOCKET
         READ (2, *, IOSTAT=IOS, ERR=1100)(Prise(I,J),J=-NPA,NPA)    !##PATCH
         IJ=-NPA
         DO J=1,NPC
          PRjet(I,J)=Prise(I,IJ)
          IJ=IJ+1
         END DO
         WRITE(42, *,IOSTAT=IOS, ERR=1199)(PRjet(I,J),J=1,NPC)
         PRECup(I)=PRjet(I,1)             ! upstream recess edge P
         PRECdo(I)=PRjet(I,NPC)           ! downstream recess edge P
        END DO

        WRITE(42,*, IOSTAT=IOS, ERR=1199) (PRECdo(I),I=1,NPOCKET)
        WRITE(42,*, IOSTAT=IOS, ERR=1199) (PRECup(I),I=1,NPOCKET)

      END IF              !=====> for HJBS only  42=TEMP*, 2=DATAfile

c NOTE: PRECup(I)=PRjet(I,1), PRECdo=PRjet(I,NPC); I=1,NPOCKET


      JSTART=-NYI
      IF (BEARING.EQ.2) JSTART=0

c ....READ FLOW Field pressure, temperature, velocities (U,V)
      DO I=1, NXT
       DO J= JSTART, NYI, 1
        READ (2, *,IOSTAT=IOS,ERR=1100) P(I,J),T(I,J),U(I,J),V(I,J)
        WRITE(42,*,IOSTAT=IOS,ERR=1199) P(I,J),T(I,J),U(I,J),V(I,J)
       END DO
      END DO
c ... FOR PAD ARC read leading & trailing edge flow and temperatures
      IF (IFULL.EQ.0) THEN
       READ (2,*, IOSTAT=IOS,ERR=1100) QLEAD(K),QTRAIL(K),TQTRAIL(K)
       WRITE(42,*,IOSTAT=IOS,ERR=1199) QLEAD(K),QTRAIL(K),TQTRAIL(K)
       READ (2,*, IOSTAT=IOS,ERR=1100) (TLEAD(K,J),J=-NYI,NYI)
       WRITE(42,*,IOSTAT=IOS,ERR=1199) (TLEAD(K,J),J=-NYI,NYI)
       READ (2,*, IOSTAT=IOS,ERR=1100) (TRAIL(K,J),J=-NYI,NYI)
       WRITE(42,*,IOSTAT=IOS,ERR=1199) (TRAIL(K,J),J=-NYI,NYI)
      END IF

c ... READ force/moment coefficients due to displacements/rotations
      READ (2, *, IOSTAT=IOS, ERR=1100)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      WRITE (42, *, IOSTAT=IOS, ERR=1199)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

c    !......................!
      IF (TILT.EQ.1) THEN
c    !......................! Reads Pad coefficients due to pad rotations
      IPAD=INERPAD(K)
      Xpivotk=(OFFSET*LPAD(K)+X1(K))*DACOS(-1.0D0)/180.0D0
      RSINPK=(DIAM/2.0D0)*DSIN(Xpivotk)
      RCOSPK=(DIAM/2.0D0)*DCOS(Xpivotk)

      READ  (2,  *, IOSTAT=IOS, ERR=1100)
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      WRITE (42, *, IOSTAT=IOS, ERR=1199)
C#   +      RSINPK, RCOSPK, IPAD,ROTPAD(K),KSTPAD(K),CDAPAD(K),
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
c    !......................!
      END IF
c    !......................!
      CALL CPARAM(RPM,PS,PA)     ! check initial parameters & miprops values
      CALL MINMAXP    ! find MAX & MIN Pressures
      CALL ADD0       ! ADD FORCES AND COEFFICIENTS
      CALL ADD1       !........
c     ...........................................................
      CLOSE (UNIT=42, IOSTAT=IOS, ERR=1299)
c     ...........................................................

C     ................!..............................................!
      END DO          ! K=1, 2,... , NPAD
C     ................!..............................................!

c     .......................! Non Uniform Exit Side pressures.
      READ (2, *, IOSTAT=IOS, ERR=1100) IPRuni, IPLuni, NPRcs, NPLcs
      IF (IPRuni.NE.1) THEN
       READ (2,*) (PRCOEFC(I), I=0, NPRCS)
       READ (2,*) (PRCOEFS(I), I=1, NPRCS)
      END IF
      IF (IPLuni.NE.1) THEN
       READ (2,*) (PLCOEFC(I), I=0, NPLCS)
       READ (2,*) (PLCOEFS(I), I=1, NPLCS)
      END IF
c     .......................! Non Uniform Exit Side pressures.


      RES=RHO*2.0*FREQU*DACOS(-1.0D0)*CLEAR*CLEAR/EMU     ! SQUEEZE FILM RE#

      IFF=IF

c ......................................!..................................
      IF (BEARING.EQ.2) THEN            ! TYP discharge pressure
        PA=Pright                       ! for seals
      ELSE                              !........... &
        PA= DMIN1(Pleft, Pright)        ! for HJBs and plain bearings
      END IF                            !
c ......................................!..................................


c     ----------------------
      IF (IF.GT.4) goto 102     ! IF=5,6,7, and 12
c     ----------------------

c     CHECK Properties at Supply and Sump Conditions
c           for cryogenic liquids.
c     ..............................................
c##   CALL FDATA(IF)    ! => initialize mipropst.f
      CALL CHECKPROPS   ! => initopst.f
c     ..............................................

 102  CONTINUE

c     ...........................................................
      READ (2, *, IOSTAT=IOS, ERR=1100) TSHAFT,TSTATOR
      READ (2, *, IOSTAT=IOS, ERR=1100) IGUESP
c     ...........................................................
      IF (IFULL.EQ.0) THEN
        READ (2,*,IOSTAT=IOS, ERR=1100) LEMDA   ! pad temperature mixing coefficient
      ELSE
         LEMDA=0.0D0
      END IF
c     ...........................................................

c     ...........................................................
       READ (2,*,IOSTAT=IOS, ERR=1500) AC, ETA, RELAXH, LIFT ! COMPLIANCE coefficients
c     ...........................................................
       IF (ISOTH.GE.4) THEN       ! read parameters for radial heat flow
         READ (2, *) TBOUT, THERMALK, ROUTER
       END IF
c     ...........................................................
       HCELL=0.0D0                                 ! STATOR SURFACE
       READ (2, *,IOSTAT=IOS, ERR=1600) HCELL      ! CELL depth [m]
c     ...........................................................

c     ...........................................................
      ANGLEJ=0.0D0                                  !  0deg, radial injection
      LOCJET=0.50                                   !  center of recess
      IF (BEARING.EQ.1) THEN                        ! orifice angle and location
       READ (2, *,IOSTAT=IOS, ERR=1610) ANGLEJ, LOCJET
      END IF
c     ...........................................................

c    !---------------------------------------------! CLOSE DATA FILE
 1205 CLOSE (UNIT=2, IOSTAT=IOS, ERR=1200)
c    !---------------------------------------------! CLOSE DATA FILE

 1210 ITER=ITMAX
      RETURN


C Error handling:
C --------------------------------------------------------------------
 77   FORMAT ('$',' For Pad#',I2,' Transfer TO TEMP File:',
     +        A9,' FROM INPUT DATA file:', A40)


 1000 WRITE (6, *) 'Error opening file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting read.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=2)
      RETURN
 1099 WRITE (6, *) 'Error opening TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)
      GOTO 1205

 1100 WRITE (6, *) 'Error during read:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting read.  The previous input data has been cor
     +rupted by this error.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=2)
      RETURN
 1199 WRITE (6, *) 'Error during WRITE TO TEMPORAL DATA FILE:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting WRITE.  The previous input data has been cor
     +rupted by this error.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=2)
      GOTO 1205

 1200 WRITE (6, *) 'Error closing file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Read completed.  Errors may occur if you attempt to
     +read or write again.'
      CALL BEEPER
      CALL PAUSE
      GOTO 1210
 1299 WRITE (6, *) 'Error closing file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'WRITE completed. Errors may occur if you attempt to
     +read or write again.'
      CALL BEEPER
      CALL PAUSE
      GOTO 1205

 1300 WRITE (6, *) 'ERROR closing file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Errors may occur if you attempt to read or write aga
     +in.'
      CALL BEEPER
      CALL PAUSE
      RETURN

 1500 WRITE (6, *) 'INFORMATION on BEARING COMPLIANCE not on FILE'
      WRITE (6, *) 'SET DEFAULT values: ZERO'
      AC=0.0D0
      ETA=0.0D0
      RELAXH=0.0D0
      LIFT=0
      CLOSE(UNIT=2)
      ITER=ITMAX
      RETURN

 1600 WRITE (6, *) 'INFORMATION on CELL DEPTH not on FILE'
      WRITE (6, *) 'SET DEFAULT value= ZERO [m]'
      HCELL=0.0D0
      HCDIM=0.0D0   ! HCELL/CLEAR
      CLOSE(UNIT=2)
      ITER=ITMAX
      RETURN

 1610 WRITE (6, *) 'INFORMATION on ORIFICE ANGLE and location'
      WRITE (6, *) 'NOT FOUND. SET DEFAULT value: radial at center'
      ANGLEJ=0.0D0
      LOCJET=0.50
      CLOSE(UNIT=2)
      ITER=ITMAX
      RETURN


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Writedata                                                   **
C **                                                                         **
C **  WRITEDATA: Writes input data and U, V, and P fields into file.         **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE WRITEDATA(FILE, DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /PEDGE/ PRISE(MAXNPOCK, -MAXNYI:MAXNYI)
      COMMON /SPLDATA/Z(Nsl),CL(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl),Nj
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD,AR,HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR

      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)

      COMMON /QIOPAD/ QLEAD(MAXNPAD),QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAW, IWEAR
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr
      COMMON /LOSPAD/ LOSleadP, KLOSpad
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FORCE0/ FFACTOR, FXO, FYO, TO, TOR
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /TGROOVE/ LEMDA,DELTA

      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /THERMAL/ ALFT,UC,TC,EC
      COMMON /OILCOEF/ Talpha
      COMMON /THERMID/ ISOTH
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /GUESPI/ IGUESP
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
c     .............................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION TORR,U, V, P, T,ALFT,UC,TC,EC,CPS,THS,BETAKS,
     +                 Z, Cl, Bcl, Ccl, Dcl,
     +                 PREC, QREC,TREC, PRISE, QIN, QOUT, QFACTOR
      DOUBLE PRECISION X1, X1r, LREC, LPAD,LENGTHL, LENGTHR,
     +                 PaDorif, ROTPAD, INERPAD,KSTPAD,CDAPAD,
     +                 LEMDA, DELTA
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,
     +                 EXO,EYO,AXO,AYO,ZO,EWX,EWY,EWEAR,BETAW,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA
      DOUBLE PRECISION LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr,
     +                 Cinlet, Cexit,ClearO,ClearR,ClearL,YR,YL,
     +                 Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 PRCOEFC,PRCOEFS,PLCOEFC,PLCOEFS
      DOUBLE PRECISION FREQU, SIGMA, L1, RES,HRECU, VSUP, BETA,
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION RSINPK, RCOSPK, IPAD,KROTPAD,CROTPAD,
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,
     +                 ALFP, FFACTOR, FXO, FYO, TO, TOR, MFACTOR, MX,MY,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 Tempk, Vsound, LOSLEADP, KLOSPAD,
     +                 PRELOAD, PRELOADB, OFFSET, ROTDEL, AC,ETA,RELAXH,
     +                 QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL ,Talpha,
     +                 TSHAFT, TSTATOR, TBOUT, THERMALK, ROUTER, HKB,
     +                 HCELL, HCDIM,ANGLEJ, LOCJET, CJET, DPJET,
     +                 PRECdo, PRECup, PRjet

      INTEGER If, IL, IWEAR, NJ, ICASE, NCASE, NREC, NPAD,
     +        ITER, ITMAX, ITPMAX, ISYM, ICSTEP, IFULL,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        IPRuni, IPLuni, NPRcs, NPLcs, BEARING,
     +        NLOBES, TILTPAD ,TILT, LIFT, ISOTH ,IGUESP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*40 FILE
      CHARACTER*2 VALK
      CHARACTER*7 SEED
      CHARACTER*9 TEMPFILE
      INTEGER DEVICE, I, J, IJ, K, IOS, JSTART, TILTK
      DOUBLE PRECISION PDUMY(MAXNPOCK, -MAXNYI:MAXNYI)
C ----------------------------------------------------------------------------
C --  WRITEDATA code                                                        --
C ----------------------------------------------------------------------------
      SEED='TEMPPAD'
      NPOCKET=0
      DO K=1, NPAD
        NPOCKET=MAX0(NPOCKET,NREC(K))
      END DO

c    !----------------------------------------------------------------!
      OPEN (UNIT=2, FILE=FILE, STATUS='UNKNOWN', IOSTAT=IOS, ERR=1000)
c    !----------------------------------------------------------------!

      WRITE (2, 20, IOSTAT=IOS, ERR=1100) TITLE
   20 FORMAT (60A)
      WRITE (2, 30, IOSTAT=IOS, ERR=1100) DDATE
   30 FORMAT (10A)

      WRITE (2,*,IOSTAT=IOS,ERR=1100) NPAD,NPOCKET,NLC,NPC,NLA,NPA
      WRITE (2,*,IOSTAT=IOS,ERR=1100) NJ,INERL,INERP,ITURB,ICAV,
     +                                ISOTH,MODEL
      WRITE (2,*,IOSTAT=IOS,ERR=1100) IF,IL,IWEAR,ISYM,ICSTEP, NLOBES
      WRITE (2,*,IOSTAT=IOS,ERR=1100) ITMAX, ITPMAX, IFULL,BEARING,
     +                                TILTPAD
      WRITE (2,*,IOSTAT=IOS,ERR=1100) (NREC(K), K=1, NPAD)

      WRITE (2,*,IOSTAT=IOS,ERR=1100)(X1(K),K=1,NPAD),     ! PAD: start angle
     +                             (LPAD(K),K=1,NPAD)      !     angular length
      WRITE (2,*,IOSTAT=IOS,ERR=1100) (ROTPAD(K),K=1,NPAD) !     rotations

      IF (TILTPAD.EQ.1) THEN
      WRITE (2,*,IOSTAT=IOS,ERR=1100) (INERPAD(K),K=1,NPAD)! pad Inertia
      WRITE (2,*,IOSTAT=IOS,ERR=1100) (KSTPAD(K),K=1,NPAD) ! pad ROT Stiffness
      WRITE (2,*,IOSTAT=IOS,ERR=1100) (CDAPAD(K),K=1,NPAD) ! pad ROT Damping
      END IF

      DO K=1, NPAD
        NPOCKET=NREC(K)
        IF (NPOCKET.GE.1) THEN                                   ! PAD: RECESS
         WRITE(2,*,IOSTAT=IOS,ERR=1100) (X1r(J,K),J=1,NPOCKET),  ! Start Angle
     +                                  (Lrec(J,K),J=1,NPOCKET), ! Angular Length
     +                               (PaDorif(J,K),J=1,NPOCKET)  ! Orifice Diams.
        END IF
      END DO

      WRITE (2, *, IOSTAT=IOS, ERR=1100) ALFU,ALFP,ALFT,SFLOW,SMASS,MMP
      WRITE (2, *, IOSTAT=IOS, ERR=1100) LOSXSIxu, LOSXSIxd, LOSXSIyl,
     +                                   LOSXSIyr, ALPHA, PRATIO
      WRITE (2, *, IOSTAT=IOS, ERR=1100) CLEAR, CINLET, CEXIT
      WRITE (2, *, IOSTAT=IOS, ERR=1100) (Z(J), J=1, NJ)
      WRITE (2, *, IOSTAT=IOS, ERR=1100) (CL(J), J=1, NJ)
      IF (ICSTEP.eq.1) THEN
        WRITE(2,*, IOSTAT=IOS, ERR=1100) ClearO,ClearL,ClearR,YR,YL
      END IF

      WRITE (2, *, IOSTAT=IOS, ERR=1100) DIAM, LENGTH,AR, HREC, VSUP,
     +                                   LENGTHR
      WRITE (2, *, IOSTAT=IOS, ERR=1100) EMUS,RHOS, CPS,THS,BETAKS,
     +                                   EMUle,RHOle,EMUri,RHOri
      WRITE (2, *, IOSTAT=IOS, ERR=1100) P1PROP, EMU1,RHO1,
     +                                   P2PROP, EMU2,RHO2
      WRITE (2, *, IOSTAT=IOS, ERR=1100) RPM,FREQU,PS,PLeft,Pright,PC
      WRITE (2, *, IOSTAT=IOS, ERR=1100) Cleft,Cright,CD,DORIF,BETA
      WRITE (2, *, IOSTAT=IOS, ERR=1100) RUGR, RUGS, AMOD, BMOD, EXPO
      WRITE (2, *, IOSTAT=IOS, ERR=1100) EXO, EYO, AXO, AYO, ZO
      WRITE (2, *, IOSTAT=IOS, ERR=1100) TEMPK, VSOUND, LOSLEADP
      WRITE (2, *, IOSTAT=IOS, ERR=1100) PRELOAD, OFFSET

      IF (IWEAR.EQ.0) THEN
         EWEAR=0.0D0
         BETAW=0.0D0
      END IF
      WRITE (2, *, IOSTAT=IOS, ERR=1100) EWX, EWY, EWEAR, BETAW

      WRITE (2, *, IOSTAT=IOS, ERR=1100) Talpha

C     ................!..............................................!
      DO I=1,MAXNPOCK              ! PATCH: DUMY array needed
       DO J=-NYI, NYI, 1           ! for compatibility with
        PDUMY(I,J)=0.0D0           ! past releases of
       END DO                      ! hydro program
      END DO                       ! 7/95 by LSA

C     ................!..............................................!
      DO K=1, NPAD
C     ................! re-store data for all bearing pads from TEMPfile
      TEMPFILE=SEED//VALK(K)

c    !..........................................................!
      OPEN (UNIT=42, FILE=TEMPFILE, STATUS='UNKNOWN', IOSTAT=IOS,
     +      ERR=1099)
c    !..........................................................!
      WRITE (6, 77) K, TEMPFILE, FILE

      READ (42, *, IOSTAT=IOS, ERR=1100) NPOCKET, NXT, TILTK
      READ (42, *, IOSTAT=IOS, ERR=1199) TOR,FXO,FYO,MX,MY,QIN,QOUT

      WRITE (2, *, IOSTAT=IOS, ERR=1100) TOR,FXO,FYO,MX,MY,QIN,QOUT
      IF (NPOCKET.GE.1) THEN
        READ (42, *, IOSTAT=IOS, ERR=1199) (PREC(I), I=1, NPOCKET)
        WRITE (2, *, IOSTAT=IOS, ERR=1100) (PREC(I), I=1, NPOCKET)
        READ (42, *, IOSTAT=IOS, ERR=1199) (QREC(I), I=1, NPOCKET)
        WRITE (2, *, IOSTAT=IOS, ERR=1100) (QREC(I), I=1, NPOCKET)
        READ (42, *, IOSTAT=IOS, ERR=1199) (TREC(I), I=1, NPOCKET)
        WRITE (2, *, IOSTAT=IOS, ERR=1100) (TREC(I), I=1, NPOCKET)
        READ (42, *, IOSTAT=IOS, ERR=1199) (TORR(I), I=1, NPOCKET)
        WRITE (2, *, IOSTAT=IOS, ERR=1100) (TORR(I), I=1, NPOCKET)
        DO I=1, NPOCKET
         READ (42, *,IOSTAT=IOS, ERR=1199)(PRjet(I,J),J=1,NPC)
         IJ=-NPA
         DO J=1,NPC
            PDUMY(I,IJ)=PRJET(I,J)
            IJ=IJ+1
         END DO
         WRITE (2, *, IOSTAT=IOS,ERR=1100)(PDUMY(I,J),J=-NPA,NPA)
        END DO
      END IF   ! NPOCKET.GE.1 (HJB), UNIT 42=TEMP*, 2=DATA file

        READ (42, *, IOSTAT=IOS, ERR=1199)(PRECdo(I),I=1,NPOCKET)
C##     WRITE( 2,*, IOSTAT=IOS, ERR=1100) (PRECdo(I),I=1,NPOCKET)
        READ (42, *, IOSTAT=IOS, ERR=1199)(PRECup(I),I=1,NPOCKET)
C##     WRITE( 2,*, IOSTAT=IOS, ERR=1100) (PRECup(I),I=1,NPOCKET)


      JSTART=-NYI
      IF (BEARING.EQ.2) JSTART=0       !for seal 360 deg

      DO I=1, NXT
       DO J= JSTART, NYI, 1
        READ (42,*,IOSTAT=IOS,ERR=1199) P(I,J),T(I,J),U(I,J),V(I,J)
        WRITE (2,*,IOSTAT=IOS,ERR=1100) P(I,J),T(I,J),U(I,J),V(I,J)
       END DO
      END DO

      IF (IFULL.EQ.0) THEN
       READ (42,*,IOSTAT=IOS,ERR=1199) QLEAD(K),QTRAIL(K),TQTRAIL(K)
       WRITE(2,*, IOSTAT=IOS,ERR=1100) QLEAD(K),QTRAIL(K),TQTRAIL(K)
       READ (42,*,IOSTAT=IOS,ERR=1199) (TLEAD(K,J),J=-NYI,NYI)
       WRITE(2,*, IOSTAT=IOS,ERR=1100) (TLEAD(K,J),J=-NYI,NYI)
       READ (42,*,IOSTAT=IOS,ERR=1199) (TRAIL(K,J),J=-NYI,NYI)
       WRITE(2,*, IOSTAT=IOS,ERR=1100) (TRAIL(K,J),J=-NYI,NYI)
      END IF

      READ (42, *, IOSTAT=IOS, ERR=1199)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      WRITE (2, *, IOSTAT=IOS, ERR=1100)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

c    !......................!
      IF (TILTK.EQ.1) THEN
c    !......................!
      READ (42, *, IOSTAT=IOS, ERR=1199)
C#   +          RSINPK, RCOSPK, IPAD,ROTDEL,KROTPAD,CROTPAD,
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      WRITE (2,  *, IOSTAT=IOS, ERR=1100)
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

c    !......................!
      END IF
c    !......................!

c    !..........................................................!
      CLOSE (UNIT=42, IOSTAT=IOS, ERR=1299)
c    !......................................! CLOSE TEMPK file .!

C     ................!..............................................!
      END DO          ! K=1, 2,... , NPAD
C     ................!..............................................!

      WRITE (2, *, IOSTAT=IOS, ERR=1100) IPRuni, IPLuni, NPRcs, NPLcs
      IF (IPRuni.NE.1) THEN
       WRITE (2,*) (PRCOEFC(I), I=0, NPRCS)
       WRITE (2,*) (PRCOEFS(I), I=1, NPRCS)
      END IF
      IF (IPLuni.NE.1) THEN
       WRITE (2,*) (PLCOEFC(I), I=0, NPLCS)
       WRITE (2,*) (PLCOEFS(I), I=1, NPLCS)
      END IF

      WRITE (2, *, IOSTAT=IOS, ERR=1100) TSHAFT,TSTATOR
      WRITE (2, *, IOSTAT=IOS, ERR=1100) IGUESP
c     ...........................................................
      IF (IFULL.EQ.0) THEN
        WRITE (2,*,IOSTAT=IOS, ERR=1100) LEMDA   ! pad temperature mixing coefficient
      END IF
c     ...........................................................
      WRITE (2,*) AC, ETA, RELAXH, LIFT          ! compliance coeffs.
c     ...........................................................
      IF (ISOTH.GE.4) THEN                       ! radial heat flow
       WRITE (2, *) TBOUT, THERMALK, ROUTER      ! parameters
      END IF                                     !
c     ...........................................................
      WRITE (2, *) HCELL                         ! CELL depth [m]
c     ...........................................................
      IF (BEARING.EQ.1) THEN
       WRITE (2, *) ANGLEJ, LOCJET               ! orifice angle and location
      END IFc     ...........................................................
      WRITE (2, *)
c     ...........................................................

c    !---------------------------------------------! CLOSE OUTPUT FILE
 1205 CLOSE (UNIT=2, IOSTAT=IOS, ERR=1200)
c    !---------------------------------------------! CLOSE OUTPUT FILE


 1210 RETURN
C.................................................
C
C Error handling:
C.................................................
 77   FORMAT ('$',' For Pad#',I2,'  Transfer from File:', A9,
     +        ' to DATA file:', A40)

 1000 WRITE (6, *) 'Error opening file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=2)
      RETURN

 1099 WRITE (6, *) 'Error opening TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)
      GOTO 1205

 1100 WRITE (6, *) 'Error during write:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.  The file is being deleted.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=2, STATUS='DELETE')
      RETURN

 1199 WRITE (6, *) 'Error during READ FROM TEMPORAL FILE:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.  The file is being deleted.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)   !, STATUS='DELETE')
      GOTO 1205

 1200 WRITE (6, *) 'Error closing file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Write completed.  Errors may occur if you attempt to
     + read or write again.'
      CALL BEEPER
      CALL PAUSE
      GOTO 1210

 1299 WRITE (6, *) 'Error closing TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'READ completed.  Errors may occur if you attempt to
     + read or write again.'
      CALL BEEPER
      CALL PAUSE
      GOTO 1205

      END
C
C *****************************************************************************
C **                                                                         **
C **  Function ValK                                                          **
C **                                                                         **
C **  VALK: Changes an integer number into a character string.               **
C **                                                                         **
C *****************************************************************************

      CHARACTER*2 FUNCTION VALK(K)

      IMPLICIT NONE


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER  K
      CHARACTER*2 DIGIT(20)
C ----------------------------------------------------------------------------
C --  VALK code                                                             --
C ----------------------------------------------------------------------------
      DATA DIGIT /'1.','2.','3.','4.','5.','6.','7.','8.','9.','10',
     +            '11','12','13','14','15','16','17','18','19','20'/

      VALK=DIGIT(K)

      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C *****************************************************************************
C **                                                                         **
C **  Subroutine Writetemp                                                   **
C **                                                                         **
C **  WRITEtemp: Writes U,V,P,T fields and results into temporal data file   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE WRITETEMP(K)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)
      COMMON /PEDGE/ PRISE(MAXNPOCK, -MAXNYI:MAXNYI)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)

      COMMON /FORCE0/ FFACTOR, FXO, FYO, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /BTYPE/ BEARING
      COMMON /GUESPI/ IGUESP
c     .............................................................
      DOUBLE PRECISION U,V,P,PREC,QREC,TREC,PRISE,QIN,QOUT,QFACTOR,
     +                 FXO,FYO,TO,TOR,FFACTOR,MX,MY,MFACTOR ,T,TORR,
     +                 PRECdo, PRECup, PRJET
      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION RSINPK, RCOSPK, IPAD,PRELOAD, OFFSET,ROTDEL,
     +                 KROTPAD, CROTPAD, KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      DOUBLE PRECISION QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL
      INTEGER IFULL,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        BEARING ,TILT, IGUESP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*2 VALK
      CHARACTER*7 SEED
      CHARACTER*9 TEMPFILE
      INTEGER  I,J,K, IOS, JSTART

C ----------------------------------------------------------------------------
C --  WRITETEMP code                                                        --
C ----------------------------------------------------------------------------
      SEED='TEMPPAD'
      TEMPFILE=SEED//VALK(K)
      JSTART=-NYI
      IF (BEARING.EQ.2) JSTART=0

      OPEN (UNIT=42, FILE=TEMPFILE, STATUS='UNKNOWN', IOSTAT=IOS,
     +      ERR=1099)

      WRITE (42, *, IOSTAT=IOS, ERR=1100) NPOCKET, NXT, TILT
      WRITE (42, *, IOSTAT=IOS, ERR=1100) TOR,FXO,FYO,MX,MY,QIN,QOUT

      IF (NPOCKET.GE.1) THEN
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (PREC(I), I=1, NPOCKET)
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (QREC(I), I=1, NPOCKET)
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (TREC(I), I=1, NPOCKET)
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (TORR(I), I=1, NPOCKET)
        DO I=1, NPOCKET
         WRITE (42, *, IOSTAT=IOS, ERR=1100)(PRjet(I,J),J=1,NPC)
        END DO
      END IF
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (PRECdo(I),I=1,NPOCKET)
        WRITE (42, *, IOSTAT=IOS, ERR=1100) (PRECup(I),I=1,NPOCKET)

      DO I=1, NXT
       DO J= JSTART, NYI, 1
        WRITE (42,*,IOSTAT=IOS,ERR=1100) P(I,J),T(I,J),U(I,J),V(I,J)
       END DO
      END DO

      IF (IFULL.EQ.0) THEN
        WRITE (42,*,IOSTAT=IOS,ERR=1100) QLEAD(K),QTRAIL(K),TQTRAIL(K)
        WRITE (42,*,IOSTAT=IOS,ERR=1100) (TLEAD(K,J),J=-NYI,NYI)
        WRITE (42,*,IOSTAT=IOS,ERR=1100) (TRAIL(K,J),J=-NYI,NYI)
      END IF

      WRITE (42, *, IOSTAT=IOS, ERR=1100)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

      IF (TILT.EQ.1) THEN
      WRITE (42, *, IOSTAT=IOS, ERR=1100)
c#   +          RSINPK, RCOSPK, IPAD,ROTDEL,KROTPAD,CROTPAD,
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      END IF
c     ...........................................................
      WRITE (42, *, IOSTAT=IOS, ERR=1100) IGUESP

c     ...........................................................
      CLOSE (UNIT=42, IOSTAT=IOS, ERR=1299)
c     ...........................................................

 1205 CONTINUE


 1210 RETURN
C.................................................
C
C Error handling:
C.................................................
 77   FORMAT ('$',' For Pad#',I2,' WRITE DATA TO TEMP File:',
     +            A9)


 1099 WRITE (6, *) 'Error opening TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)
      RETURN

 1100 WRITE (6, *) 'Error during WRITE TO  TEMPORAL FILE:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting write.  The file is being deleted.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)   !, STATUS='DELETE')
      RETURN


 1299 WRITE (6, *) 'Error closing TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'WRITE completed. Errors may occur if you attempt to
     + read or write again.'
      CALL BEEPER
      CALL PAUSE
      RETURN

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine  Readtemp                                                   **
C **                                                                         **
C **   READtemp: Reads U,V,P fields and results from temporal data file      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE READTEMP(K)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)
      COMMON /PEDGE/ PRISE(MAXNPOCK, -MAXNYI:MAXNYI)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /FORCE0/ FFACTOR, FXO, FYO, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /BTYPE/ BEARING
      COMMON /GUESPI/ IGUESP
c     .............................................................
      DOUBLE PRECISION U,V,P,PREC,QREC,TREC,PRISE,QIN,QOUT,QFACTOR,
     +                 FXO,FYO,TO,TOR,FFACTOR,MX,MY,MFACTOR ,T,TORR,
     +                 PRECdo, PRECup, PRjet
      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION RSINPK, RCOSPK, IPAD,PRELOAD, OFFSET,ROTDEL,
     +                 KROTPAD, CROTPAD, KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      DOUBLE PRECISION QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL
      INTEGER IFULL,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        BEARING ,TILT, IGUESP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*2 VALK
      CHARACTER*7 SEED
      CHARACTER*9 TEMPFILE
      INTEGER  I,J,K, IOS, NPOCKETD, NXTD, JSTART,  TILTDUMY

C ----------------------------------------------------------------------------
C --  READTEMP code                                                        --
C ----------------------------------------------------------------------------
C NOTE: Here we need not to upset the values NPOCKET & NXT and make sure
C###### that NYI, NPA are the same as with main call

      JSTART=-NYI
      IF (BEARING.EQ.2) JSTART=0

      SEED='TEMPPAD'
      TEMPFILE=SEED//VALK(K)

      OPEN (UNIT=42, FILE=TEMPFILE, STATUS='UNKNOWN', IOSTAT=IOS,
     +      ERR=1099)

      READ (42, *, IOSTAT=IOS, ERR=1100) NPOCKETD, NXTD, TILTDUMY
      READ (42, *, IOSTAT=IOS, ERR=1100) TOR,FXO,FYO,MX,MY,QIN,QOUT
      IF (NPOCKETD.GE.1) THEN
        READ (42, *, IOSTAT=IOS, ERR=1100) (PREC(I), I=1, NPOCKETD)
        READ (42, *, IOSTAT=IOS, ERR=1100) (QREC(I), I=1, NPOCKETD)
        READ (42, *, IOSTAT=IOS, ERR=1100) (TREC(I), I=1, NPOCKETD)
        READ (42, *, IOSTAT=IOS, ERR=1100) (TORR(I), I=1, NPOCKETD)
        DO I=1, NPOCKETD
         READ (42, *, IOSTAT=IOS, ERR=1100)(PRjet(I,J),J=1,NPC)
        END DO
      END IF
        READ (42, *, IOSTAT=IOS, ERR=1100) (PRECdo(I),I=1,NPOCKETD)
        READ (42, *, IOSTAT=IOS, ERR=1100) (PRECup(I),I=1,NPOCKETD)

      DO I=1, NXT
        DO J= JSTART, NYI, 1
          READ (42,*,IOSTAT=IOS,ERR=1100) P(I,J),T(I,J),U(I,J),V(I,J)
        END DO
      END DO

      IF (IFULL.EQ.0) THEN
        READ (42,*,IOSTAT=IOS,ERR=1100) QLEAD(K),QTRAIL(K),TQTRAIL(K)
        READ (42,*,IOSTAT=IOS,ERR=1100) (TLEAD(K,J),J=-NYI,NYI)
        READ (42,*,IOSTAT=IOS,ERR=1100) (TRAIL(K,J),J=-NYI,NYI)
      END IF

      READ (42, *, IOSTAT=IOS, ERR=1100)
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

      IF (TILT.EQ.1) THEN
      READ (42, *, IOSTAT=IOS, ERR=1100)
c#   +           RSINPK, RCOSPK, IPAD, ROTDEL,KROTPAD,CROTPAD,
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      END IF


c     ...........................................................
      READ (42, *, IOSTAT=IOS, ERR=1100) IGUESP

c     ...........................................................
      CLOSE (UNIT=42, IOSTAT=IOS, ERR=1299)
c     ...........................................................

 1205 CONTINUE

 1210 RETURN
C.................................................
C
C Error handling:
C.................................................
 77   FORMAT ('$',' For Pad#',I2,' READ DATA from TEMP File:',
     +            A9)


 1099 WRITE (6, *) 'Error opening TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting READ.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)
      RETURN

 1100 WRITE (6, *) 'Error during READ FROM  TEMPORAL FILE:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'Aborting READ.  The file is being deleted.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)   !, STATUS='DELETE')
      RETURN


 1299 WRITE (6, *) 'Error closing TEMPORAL DATA file:'
      CALL DECODIOS(IOS)
      WRITE (6, *) 'READ completed. Errors may occur if you attempt to
     + read or write again.'
      CALL BEEPER
      CALL PAUSE
      RETURN

      END

C *****************************************************************************
C LAST revised by 6/28/94 by Dr. Luis San Andres
C READ/WRITE HCELL on 6/20/95
C READ/WRITE angled jet injection parameters 7/7/95
C *****************************************************************************