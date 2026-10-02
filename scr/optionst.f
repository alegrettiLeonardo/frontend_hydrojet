c 7/3/95 updated for angle injection in HJBS
c ------------------------------------------------
c 6/20/95 UPDATE for CELL DEPTH effect
c -----------------------------------------------------------------------------
c SUBS MinmaxP & MinmaxT are located at END of this file
c
c SEARCH for >> COMPLI on bearing compliance effects
c
c  ####   #####    #####     #     ####   #    #   ####    #####          ######
c #    #  #    #     #       #    #    #  ##   #  #          #            #
c #    #  #    #     #       #    #    #  # #  #   ####      #            #####
c #    #  #####      #       #    #    #  #  # #       #     #     ###    #
c #    #  #          #       #    #    #  #   ##  #    #     #     ###    #
c  ####   #          #       #     ####   #    #   ####      #     ###    #
c
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
C
C Length=LengthL+LengthR
C
C Pright & Cright at Y=LengthR
C
C     |----------------------------------------------| ^
C ^   |             Br                               | |
C |y  |            <--->                             | | LengthR
C | P |   -----    ----- ^  -----    -----    -----  | |
C | G |   |   |    |   | |  |   |    |   |    |   |  | v
C |.R | ..|   |    |   | |Ar|   |    |   |    |   |  | ----> x
C   O |   |   |    |   | |  |   |    |   |    |   |  | ^
C   O |   -----    ----- v  -----    -----    -----  | | LenghtL
C   V |                                              | |
C   E |----------------------------------------------| v
C
C    Xlead              Hydrostatic PAD           Xtrail
C
C Pleft & Cleft at Y=-LengthL

C *****************************************************************************
C **                                                                         **
C **  Subroutine Kinput                                                      **
C **                                                                         **
C **  KINPUT:  Enters dimensional input data from keyboard.                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE KINPUT(RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NRECESS
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /THERMID/ ISOTH

      COMMON /CONT/ IFF
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /BTYPE/ BEARING
c     .............................................................

      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION  EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                  FREQU, SIGMA, L1, RES, TEMPK, VSOUND,
     +                  LENGTHL, LENGTHR,TSHAFT, TSTATOR,
     +                  CLEAR, DIAM, LENGTH, LD, AR, HREC

      INTEGER ICASE, NRECESS, ISOTH,
     +        INERL,INERP,ITURB,INTER,ICAV,MODEL,
     +        IFF, IF, IL, ITER, ITMAX, ITPMAX,BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*40 DUMYC
      CHARACTER*60 NEWTITLE
      CHARACTER*8 BCASE
      CHARACTER*1 YN
      CHARACTER*50 THERM
      DOUBLE PRECISION DUMY
      INTEGER DUMYI, IOS, imenu, imodel, BEARINGo
C ----------------------------------------------------------------------------
C --  KINPUT code                                                           --
C ----------------------------------------------------------------------------
C Unit 42 is used as a temporal place for translations between formats.
C.........................................................................

      WRITE (6, *) ' '
      OPEN (UNIT=42,  STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)

C..............................................................

      write (6, *) '$ ENTER INPUT DATA for HYDROFLEX '
      write (6, *) '--------------------------------- '
      write (6, *) ' '

 777  WRITE (6, *) '-------------------------------------------'
      WRITE (6, *) 'ENTER Numerical parameters:'
      WRITE (6, *) 'To accept the default value, press <RETURN>'
      WRITE (6, *) '-------------------------------------------'
      WRITE (6, *) ' '

C..............................................................

C DEFAULTS:
      itmax = MAX0(itmax,199)     ! Maximum number of its on lands
      itpmax=10                   ! Maximum number of its on rec/pressures
      BEARINGo=BEARING

C..............................................................
C  >>> THERMAL CASES are:
      IF(ISOTH.EQ. 1) THERM='Isothermal fluid film(T=Tin=Constant)'
      IF(ISOTH.EQ. 2) THERM='Adiabatic journal & Iso-bearing (Tb=Ts)'
      IF(ISOTH.EQ. 3) THERM='Iso-journal (Tj=Ts) & Adiabatic bearing'
      IF(ISOTH.EQ. 4) THERM='Adiabatic journal & Bearing radial heat'
      IF(ISOTH.EQ. 5) THERM='Iso-journal(Tj=Ts)& Bearing radial heat'
      IF(ISOTH.EQ. 0) THERM='Iso-journal & bearing: ( Tj=Ts, Tb=Ts)'
      IF(ISOTH.EQ.-1) THERM='Adiabatic bounding Surfaces (Qb=Qj=0)'
C..............................................................

 601  FORMAT (' ', 75('-'),/,3X,
     +       'SELECT BEARING TYPE [Def:',I2,']',/,3X,
     +       '(1) HYDROSTATIC PAD BEARING, (2) ANNULAR SEAL',/,3X,
     +       '(3) CYLINDRICAL JOURNAL/ TILT PAD BEARING',/,1X,
     +        75('.'))

      CALL BEEPER
      WRITE (6,601) BEARING
      CALL ENTERINT(BEARING)

      WRITE (6,2144) ISOTH, THERM
 2144 FORMAT(' ',3X,'ISOTH=',I2,' THERMAL Case:',A50,/,
     +       ' ',75('.'))

C  !.....................................!
      IF (BEARING.NE.BEARINGo) THEN
C  !.....................................!
      WRITE (6,602)
 602  FORMAT(2X, 38('. '),/,
     +       ' $ WARNING: You have selected a different BEARING type',
     +       ' than value on memory.',/,
     +       ' $ PLEASE remember to review ALL INPUT DATA on this',
     +       '   section before you',/,
     +       ' $ attempt to',
     +       ' perform any calculations with OPTS(31-34)',/,
     +       2X, 38('. '))
      CALL BEEPER
C  !.....................................!
      END IF
C  !.....................................!

C  Resets entrance inertia effects for BEARINGS
C  !.....................................!
       IF (BEARING.GE.3) THEN
C  !.....................................! BEARINGS
          BCASE='BEARING'
          BEARING=3
          INERP=0                        ! no inertia at Z=0
          MODEL=2
          IF (LENGTHR.EQ.LENGTH) THEN
              LENGTHR=LENGTH/2.0D0
              LENGTHL=LENGTHR
          END IF

C  !.....................................!
      ELSE IF (BEARING.EQ.2) THEN
C  !.....................................! ANNULAR SEAL
          BCASE='SEAL'
          INERP=1                        ! YES inertia at Z=0
          MODEL=2
          LENGTHR=LENGTH
          LENGTHL=0.0D0

C  !.....................................!
      ELSE IF (BEARING.LE.1) THEN
C  !.....................................! HYDROSTATIC BEARING
          BEARING=1
          INERP=1                        ! YES inertia at Z=REC edge

         IF (MODEL.EQ.2) IMODEL=1
         IF (MODEL.EQ.1) IMODEL=2
      WRITE (6, *) '-----------------------------------------'
      WRITE (6, *) 'SELECT MODEL: 1=SINGLE ROW, 2=DOUBLE ROW '
      WRITE (6, *) '-----------------------------------------'
      CALL BEEPER
      WRITE (6, 3456) IMODEL
 3456 FORMAT ('$', 'ENTER MODEL [Default =',I2,']: ')

      CALL ENTERINT(IMODEL)

          IF (IMODEL.LE.1) THEN
             MODEL=2
             BCASE='1row HJB'
          ELSE IF (IMODEL.GE.2) THEN
             MODEL=1
             BCASE='2row HJB'
          END IF

C    !............................!
      END IF
C    !............................!

c ............................ MENU .....................................

 999  write (6, 998) BCASE,BCASE
 998  format (' ',50('='),/,' INPUT/MODIFY DATA MENU: Select Choice '
     + 'for:',A8,/,' ',50('.'),/,
     + '(1) ANALYSIS TITLE/DATE',/,
     + '(2) ',A8,':GEOMETRY: L,D, PADS & RECESS, CLEARANCE,',/,
     + '(3) ROUGHNESS Surface/ CELL depth values,',/,
     + '(4) Journal Eccentricity  and Misalignment Angle,',/,
c    + '(5) Fluid Type, Temperature & B factor,',/,
     + '(6) OPERATING CONDITIONS: RPM, Psupply, Psump, ',/,
     + '    Temperature, FLUID TYPE & Properties, PRATIO',/,
     + '(7) Cd, Loss & Seal Coefficients, Swirl velocity,',/,
     + '    BEARING compliance and loss factor,',/,
     + '(8) Numerical convergence factors: Iters., relax params,',/,
     + '(9) MOODYs Friction Factors formulae coeffs: A, B, EXPO,',/,
     + '(10) Options for Thermal Analysis in Fluid Film,',/,
     + '(11) ACCEPT current values & EXIT,',/, ' ',40('='))

c ............................ MENU .....................................

      Imenu=1
      READ(5, 1000) DUMYC
      IF (DUMYC.NE.' ') THEN
          WRITE (42, 1010, IOSTAT=IOS, ERR=1234) DUMYC
          BACKSPACE (UNIT=42)
          READ (42, *, IOSTAT=IOS, ERR=1234) DUMYI
          Imenu=IABS(DUMYI)

      ELSE
          WRITE (6, *) ' PLEASE Select an option'
          GOTO 999
      END IF

      IF (Imenu.eq.0)  goto 912      ! SET Inerp & Inerl FLAGS
      IF (Imenu.gt.11) goto 999      !==> MENU

      go to (901,902,903,904,905,906,907,908,909,910,911), Imenu


C *****************************************************************************
C  SELECTION (1) TITLE and DATE
C *****************************************************************************

  901 WRITE (6, *) 'TITLE (CHAR*60): Title of data.'
      WRITE (6, 10) TITLE
C........................................................
   10 FORMAT (' ', '[hydrojet: ', 60A)
      write (6, 11)
   11 FORMAT ('$','Input ANALYSIS TITLE:')

      READ (5, 15) NEWTITLE
   15 FORMAT (60A)
      IF (NEWTITLE.NE.' ') THEN
          TITLE=NEWTITLE
      END IF
      WRITE (6, *) ' '
C............................................................
      WRITE (6, *) 'DATE (CHAR*10): Date of data.'
      WRITE (6, 20) DDATE
   20 FORMAT ('$', 'Input date: MM/DD/YY[DEF:',A10,']: ')
      READ (5, 25) DDATE
   25 FORMAT (10A)
C#      IF (DDATE.EQ.' ') THEN
C#          CALL DATE(DDATE)
C#      END IF
      WRITE (6, *) ' '

      GOTO 999

C *****************************************************************************
C SELECTION (2): BEARING GEOMETRY
C *****************************************************************************

 902  write (6, *) '---------------------------------------------'
      write (6, *) '(2) INPUT: Bearing PAD DIMENSIONs:'
      write (6, *) '---------------------------------------------'

      CALL INPUTPADS(RPM,PS,PA)    !=> on program calcmeshp.f

      goto 999

C *****************************************************************************
C SELECTION (3): RELATIVE SURFACE ROUGHNESS & CELL DEPTH COEFFICIENTS
C *****************************************************************************

 903  CALL INPUTRGRS    !=> turbcoefs.f
      CALL INPUTCELL    !=> inputp.f     ! 6/20/95
      GOTO 999
C
C *****************************************************************************
C
C SELECTION (4): journal center coordinates and misalignment
C                angles:
C *****************************************************************************

 904  CALL INPUTEXEY    !=> program inputp.f

      GOTO 999


C *****************************************************************************
C
C SELECTION (5): Physical properties and operating conditions
C                at Inlet (Supply) and Discharge (Exit)
C *****************************************************************************

 905  continue

      WRITE(6,1020)
      WRITE(6,1030)
 1020 FORMAT('  SELECT A FLUID BY ENTERING THE CORRESPONDING NUMBER',/)
 1030 FORMAT('  1=PARA HYDROGEN,  2=NITROGEN,  3=OXYGEN,  4=METHANE'
     C ,/,   '  5=AIR,  6=WATER,  7=OIL,      12=OTHER',/,
     C       '                          unknown to MIPROPS')
      CALL BEEPER
      write (6, *) '--------------------------------------------'
      write (6, *) '(5):SELECT/ENTER  FLUID TYPE: Default is ', IF
      write (6, *) '--------------------------------------------'

      CALL ENTERINT(IF)

      IF (IF.LT.1) IF=12
      IF (IF.GT.8) IF=12

      IFF=IF

  148 write (6, *) 'ENTER Supply Fluid Temperature (degree Kelvin)'
      write (6, *) '----------------------------------------------'
      write (6, *) ' '

      WRITE (6, *) 'TEMPK (REAL*8):Supply Fluid Temp (T) [degK]'
      WRITE (6, 149) TEMPK
  149 FORMAT ('$', 'ENTER Tempk [Default (T)= ',E12.5E2,']: ')

      CALL ENTERVAL(TEMPK)

      TEMPK=DABS(TEMPK)

c ...............................................................


      IF (IF.EQ.12.OR.IF.EQ.7) THEN  !=> UNKNOWN FLUID or OIL
c    !....................!

        CALL INPUTPROPS   ! => on program inputp.f

      ELSE IF (IF.GT.4) THEN
c    !....................! => Air, Water, and Oil

        GOTO 999

      ELSE
c    !....................!  CRYOGENIC
        CALL FDATA(IF)    ! => initialize miprops.f
        CALL CHECKPROPS   ! => on program initopsp.f

c    !....................!
      END IF

c ...............................................................
c
c  FOR: IF=1 LH2, IF=2 N2, IF=3 LO2, IF=4 CH4
c       check if Ps,Pa,Pc & T are within range & calculate properties
c

      GOTO 999


C *****************************************************************************
C
C SELECTION (6): BEARING OPERATING CONDITIONS
C
C *****************************************************************************
 906  CALL INPUTOPS    !=> ON PROGRAM inputp.f

C...........................................! SET:
      frequ=rpm/60.                         ! excitation frequency
      if (rpm.le.0.0) then                  ! w=rotational speed,
         frequ=100.                         ! but if rpm=0.
      end if                                ! take w=100Hz
C...........................................!.....................

      goto 905        !=> INPUT FLUID TYPE & PROPERTIES


C *****************************************************************************
C
C SELECTION (7): ENTRANCE LOSS PARAMETERS
C
C *****************************************************************************

 907  CALL INPUTLOSS    ! => program inputp.f

      CALL INPUTCOMPL   ! => program foilsubs.f COMPLIANCE COEFFS

      GOTO 999


C *****************************************************************************
C
C SELECTION (8): PARAMETERS FOR NUMERICAL CONVERGENCE
C
C *****************************************************************************

 908  CALL INPUTCONV     ! => program inputp.f

      GOTO 999

C *****************************************************************************
C
C SELECTION (9): COEEFICIENTS for MOODY's friction factor formula
C
C *****************************************************************************

  909 CALL INPUTMOODY    ! => program turbcoefs.f

      GOTO 999

C *****************************************************************************
C
C SELECTION (10): Options for THERMOHYDRODYNAMIC Analysis
C                 and SHAFT/STATOR Temperature.
C
C *****************************************************************************

  910 CALL INPUTTHD     ! => program inputp.f

      GOTO 999


C *****************************************************************************
C
C SELECTION (0): SET INERTIA on LANDS/RECESS EDGES
C
C *****************************************************************************
C
C SET for HYDROJET PROGRAM    ->  ALWAYS
C
C     INERL = 1;    YES fluid inertia on lands
C
C     INERP = 1;    YES fluid inertia at recess edges
C

C...............................................................
C Logical Parameters (1: Yes, 0: No)
C...............................................................

 912   WRITE (6, *) 'Logical parameters:  Y=yes, N=no.'
       CALL BEEPER

       IF (INERL.EQ.1) THEN
          WRITE (6, 264) 'Y'
       ELSE
          WRITE (6, 264) 'N'
       END IF
  264  FORMAT ('$', 'Enter INERL [Default ',A1,']: ')
       READ (5, 265) YN
  265  FORMAT (A1)
       IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
           INERL=1
       ELSE IF ((YN.EQ.'N').OR.(YN.EQ.'n')) THEN
          INERL=0
       END IF
C
       IF (BEARING.EQ.3) goto 999

       IF (INERP.EQ.1) THEN
           WRITE (6, 270) 'Y'
       ELSE
           WRITE (6, 270) 'N'
       END IF
  270  FORMAT ('$', 'Enter INERP [Default ',A1,']: ')
       READ (5, 265) YN
       IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
           INERP=1
       ELSE IF ((YN.EQ.'N').OR.(YN.EQ.'n')) THEN
           INERP=0
       END IF
C
C      INERL=1
C      INERP=1

       goto 999

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 911   CONTINUE	                ! Calculate initial parameters
       CALL CPARAM(RPM,PS,PA)		!==> initopst.f

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 1000 FORMAT (40A)
 1010 FORMAT (' ', 40A)

      CLOSE (UNIT=42)
      RETURN
C---------------------------------------------------------------
C Error in I/O to unit 42.
C
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) '$ I/O ERROR: RETURNS to main loop.'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=42)
      RETURN

C
      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Newdef                                                      **
C **                                                                         **
C **  NEWDEF:  Writes new DEFAULT.DAT file.                                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE NEWDEF

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAW, IWEAR
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr
      COMMON/ LOSPAD/ LOSleadP, KLOSPAD
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /spldata/Z(nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl),Nj
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES

      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /OILCOEF/ Talpha
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /THERMID/ ISOTH
      COMMON /TGROOVE/ LEMDA,DELTA
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET

      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c     .........................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION X1,LPAD,X1R,LREC,PADORIF,ROTPAD,
     +                 INERPAD, KSTPAD, CDAPAD,
     +                 AC, RELAXH, ETA
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD,AR,HREC,KLOSPAD,
     +                 EXO, EYO,EWX,EWY,EWEAR,BETAW,AXO,AYO,ZO,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr,
     +                 CINLET, CEXIT, ClearO,ClearR,ClearL,YR,YL,
     +                 z,cl,bcl,ccl,dcl,Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU, BETU,ALFP,
     +                 FREQU, SIGMA,L1,RES,AMOD,BMOD,RUGR,RUGS,EXPO,
     +                 HRECU,VSUP,BETA,TEMPK,VSOUND,LENGTHR,LENGTHL,
     +                 PRATIO,CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 PRCOEFC,PRCOEFS,PLCOEFC,PLCOEFS,
     +                 LOSLeadP, PRELOAD, OFFSET,ROTDEL,PRELOADB,
     +                 ALFT, UC, TC, Ec ,CPS,THS,BETAKS ,Talpha,
     +                 LEMDA, DELTA, TSHAFT, TSTATOR,
     +                 TBOUT, THERMALK, ROUTER, HKB,HCELL,HCDIM,
     +                 ANGLEJ, LOCJET, CJET, DPJET

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL, LIFT, ISOTH,
     +        NPAD, NREC, IFULL, BEARING, NLOBES, TILTPAD,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ITER, ITMAX, ITPMAX, IWEAR,IPRuni,IPLuni,NPRcs,NPLcs,
     +        Nj, ICASE, Ncase, IF, IL, ISYM, ICSTEP, J, K, NPOCK

C ----------------------------------------------------------------------------
C --  NEWDEF code                                                           --
C ----------------------------------------------------------------------------
      OPEN (UNIT=2, FILE='DEFAULT.DAT', STATUS='UNKNOWN', ERR=10)
      WRITE (2, 100) TITLE
  100 FORMAT (A60)
      WRITE (2, 110) DDATE
  110 FORMAT (A10)
      WRITE (2, *) NPAD,NPOCKET,NLC,NPC,NLA,NPA,IFULL,
     +             NJ,MODEL,NLOBES, TILTPAD
      WRITE (2, *) (NREC(K),K=1,NPAD)                   ! # or recess/pad
      WRITE (2, *) CINLET,CEXIT,DIAM,LENGTH,AR,HREC,VSUP,LENGTHR
      WRITE (2, *) (X1(K),K=1,NPAD),(LPAD(K),K=1,NPAD)  ! loc. and length of pads
      WRITE (2, *) (ROTPAD(K),K=1,NPAD)                 ! rotations of tilt-pads
      IF (TILTPAD.EQ.1) THEN
      WRITE (2, *) (INERPAD(K),K=1,NPAD)                ! pad INERTIA moments
      WRITE (2, *) (KSTPAD(K),K=1,NPAD)                 ! pivot stiffness
      WRITE (2, *) (CDAPAD(K),K=1,NPAD)                 ! pivot damping
      END IF

      DO K=1, NPAD
        NPOCK=NREC(K)
        IF (NPOCK.GT.0) THEN
           WRITE (2,*)  (X1r(J,K),J=1,NPOCK),           ! loc. relative /pad
     +                 (Lrec(J,K),J=1,NPOCK),           ! length of recess
     +              (PaDorif(J,K),J=1,NPOCK)            ! recess orifices.
        END IF
      END DO
      WRITE (2, *) (Z(J),J=1, NJ)                       ! C(Z) Coordinates &
      WRITE (2, *) (CL(J),J=1, NJ)                      !      values
      WRITE (2, *) EXO,EYO,AXO,AYO,ZO,EWX,EWY
      WRITE (2, *) EMUS,RHOS,EMUA,RHOA,BETA
      WRITE (2, *) CPS,THS,BETAKS ,Talpha
      WRITE (2, *) RPM,PS,PA,PC
      WRITE (2, *) Pleft, RHOle, EMUle, Cleft
      WRITE (2, *) Pright, RHOri, EMUri, Cright
      WRITE (2, *) P1PROP,P2PROP,EMU1,EMU2,RHO1,RHO2
      WRITE (2, *) CD,DORIF,LOSXSI,ALPHA,PRATIO
      WRITE (2, *) LOSXSIxu,LOSXSIxd,LOSXSIyl,LOSXSIyr,LOSLEADP
      WRITE (2, *) TEMPK, VSOUND, RUGR,RUGS, AMOD, BMOD, EXPO
      WRITE (2, *) INERL,INERP,ITMAX,ITPMAX,IF,IL,IWEAR,ICSTEP,
     +             BEARING ,ISOTH
      WRITE (2, *) ALFU,ALFP,ALFT,MPEPS,PREPS,SFLOW
      WRITE (2, *) PRELOAD, OFFSET
      IF (ICSTEP.eq.1) THEN
        write (2,*) ClearO,ClearR,ClearL,YR,YL
      END IF

      write (2,*) IPRuni, IPLuni, NPRcs, NPLcs

      IF (IPRuni.NE.1) THEN
        write (2, *) PRCOEFC(0)
        DO j=1, NPRcs
          write(2,*) PRCOEFC(j),PRCOEFS(j)
        END DO
      END IF

      IF (IPLuni.NE.1) THEN
        write (2, *) PLCOEFC(0)
        DO j=1, NPLcs
          write(2,*) PLCOEFC(j),PLCOEFS(j)
        END DO
      END IF

      IF (IFULL.EQ.0) THEN     ! save mixing coefficient
        write (2, *) LEMDA     ! for groove temperature
      END IF                   ! updates

      WRITE (2, *) AC, ETA, RELAXH, LIFT  !write COMPLIANCE params.

      WRITE (2, *) TSHAFT, TSTATOR        !shaft/bearing temperatures

      IF (ISOTH.GE.4) THEN                ! parameters for radial heat
       WRITE (2, *) TBOUT, THERMALK, ROUTER
      END IF

      WRITE (2, *) HCELL                  ! cell depth 6/20/95

      IF (BEARING.EQ.1) THEN
       WRITE (2, *) ANGLEJ, LOCJET        ! orifice angle & location
      END IF
      WRITE (2, *)

      CLOSE (UNIT=2)

      WRITE (6, *) 'Default data written to file DEFAULT.DAT.'
      RETURN

   10 WRITE (6, *) 'Error accessing file DEFAULT.DAT.'
      WRITE (6, *) 'Write cancelled.'
      CALL BEEPER
      CALL PAUSE

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Input                                                       **
C **                                                                         **
C **  INPUT:  Read DEFAULT.DAT or internal data set.                         **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUT

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
C ----------------------------------------------------------------------------
      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD,SPEEDS,FECC
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAW, IWEAR
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,
     +           RPM(MAXNXT),PS(MAXNXT),PA(MAXNXT),WX(MAXNXT),WY(MAXNXT)
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr
      COMMON /LOSPAD/ LOSLeadP,KLOSPad
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /spldata/Z(nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl),Nj
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES

      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /OILCOEF/ Talpha
      COMMON /THERMID/ ISOTH
      COMMON /TGROOVE/ LEMDA,DELTA
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /CONT/ IFF
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c     .........................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION X1, X1r, LPAD, Lrec,PaDorif, ROTPAD,
     +                 INERPAD, KSTPAD, CDAPAD, AC, ETA, RELAXH
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD,AR,HREC,KLOSPad,
     +                 EXO, EYO,EWX,EWY,EWEAR,BETAW,AXO,AYO,ZO,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr,
     +                 CINLET, CEXIT, ClearO,ClearR,ClearL,YR,YL,
     +                 z,cl,bcl,ccl,dcl,Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU, BETU,ALFP,
     +                 FREQU, SIGMA,L1,RES,AMOD,BMOD,RUGR,RUGS,EXPO,
     +                 HRECU,VSUP,BETA,TEMPK,VSOUND,LENGTHR,LENGTHL,
     +                 PRATIO,CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 LOSLeadP, PRELOAD, PRELOADB, OFFSET, ROTDEL,
     +                 ALFT, UC, TC, Ec, CPS,THS,BETAKS ,Talpha,
     +                 TBOUT, THERMALK, ROUTER, HKB, HCELL, HCDIM,WX,WY

      DOUBLE PRECISION BR,PRCOEFC,PRCOEFS,PLCOEFC,PLCOEFS,CONVER,
     + 			LEMDA, DELTA, TSHAFT, TSTATOR,
     +                  ANGLEJ, LOCJET, CJET, DPJET,ZERO

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL, ISOTH,
     +        NPAD, NREC, IFULL, BEARING, LIFT,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ITER, ITMAX, ITPMAX, IWEAR, IFF,
     +        Nj, j, Icase, Ncase, IF, IL, ISYM, ICSTEP, NPOCK,K,
     +        IPRuni, IPLuni, NPRcs, NPLcs, NLOBES, TILTPAD, SPEEDS,FECC

C ----------------------------------------------------------------------------
C --  INPUT code                                                            --
C ----------------------------------------------------------------------------
      CONVER=180.0D0/DACOS(-1.0D0)   ! RAD ==> DEG

      OPEN (UNIT=2, FILE='DEFAULT.txt', STATUS='OLD', ERR=10)

      READ (2, 100) TITLE
  100 FORMAT (A60)
      READ (2, 110) DDATE
  110 FORMAT (A10)

      READ (2, *) NPAD,NPOCKET,NLC,NPC,NLA,NPA,IFULL,
     +            NJ,MODEL,NLOBES,TILTPAD
    !'\ NPAD: # of pads,
    !'\ NPOCKET: # pockets/recesses in HJB, =0 for BEARING=2 or 3
    !'\ Circumference NLC: # grids on land, NPC: on pockets (<>0 for BEARING=1)
    !'\ Axial         NLA: # grids """""  , NPA: on pockets (<>0 """")
    !'\ IFULL=1 (360 deg bearing/seal), =0 (pad bearing)
    !'\ NJ =     2 (def) #### # of splines for clearance fit = Nj_spline
    !'\ MODEL = 1 (one row of pockets) BEARING=1 only
    !'\ NLOBES = 1,2....
    !'\ TILTPAD = 1 if pads tilting, 0: fixed pads.

      READ(2,*) INERL,INERP,ITMAX,ITPMAX,IF,IL,
     +          IWEAR,ICSTEP,BEARING ,ISOTH, LIFT

    !'\ ------------------------------------------------------------
    !'\ INERL: (1) YES, inertia on lands, (0) NO, inertia on lands,
    !'\ INERP: (1) YES, INERTIA at entrance, (0) NO
    !'\ ITMAX: Max # iterations on film lands (equations)
    !'\ ITPMAX: Max # iterations for orifice pressure equations (BEARING=1 HJB only)
    !'\ IF: Fluid type
    !'\ IL: ?????  = 1 ???? (will be reset by hydrojet02 if cryogenic)
    !'\ IWEAR=0   (not worn bearing/seal)
    !     iwear = 0
    !'\ ICSTEP=0  (NOT stepped AXIAL clearance)
    !     icstep = 0
    !'\ BEARING=2 (SEAL), (1)= HJB, (3)= hydrodynamic bearing
    !'\ ISOTH=1   (ISOTHERMAL model)
    !'\ LIFT 0(N) or 1(Y) coefficient for foil bearing dettachment (set =0)

      READ(2,*) IPRuni,IPLuni,NPRcs,NPLcs

    !'\ ------------------------------------------------------------
    !'\  conditions for uniform/nonuniform LEFT and RIGHT exit pressures
    !'\ ------------------------------------------------------------
    !'\ IPRuni , IPLuni, =1 (uniform pressures) on right and left sides
    !'\ NPRcs, NPLcs  (number of Fourier coefficients) >=0

      READ (2, *) CINLET,CEXIT,DIAM,LENGTH,AR,HREC,VSUP,LENGTHR
      LENGTHL=LENGTH-LENGTHR

    !'\ ------------------------------------------------------------
    !'\ Inlet and Outlet Clearances, Diameter
    !'\ Length, Axial length of recess, recess depth
    !'\  Vsup: Orifice supply volume,
    !'\  LENGTH= LEFT + RIGHT, IF RIGHT = LENGTH THEN SEAL
    !'\                           RIGHT = LENGTH/2 THEN HYDRODYN BEARING
    !'\ LENGTHR determines origin of axial coordinate.

      IF ((BEARING.EQ.1).OR.(BEARING.EQ.3)) THEN
        DO K = 1, NPAD
            READ (2, *) NREC(K), X1(K), LPAD(K), ROTPAD(K)
        END DO
      END IF

    !'\ ------------------------------------------------------------
    !'\ PADS # pockets, location, extent and rotation
    !'\ ------------------------------------------------------------
    !'\F      READ (2,*) NREC(K), X1(K), LPAD(K), ROTPAD(K)
    !'\F      K=1, NPAD
    !'\F      NREC(K)  pad # or recess/pad
    !'\F      (X1(K),(LPAD(K),:leading edge and length of bearing PADS(degrees)
    !\F       (ROTPAD(K): pad rotations (RADIANS)

    !\## STORE PADS of same length

      IF (TILTPAD.EQ.1) THEN !SET PADS properties
        DO K=1, NPAD
            READ(2,*) INERPAD(K), KSTPAD(K) , CDAPAD(K)
        END DO
      END IF
    !'\ ------------------------------------------------------------
    !'\  IF TILTPAD=1, SET PADS properties
    !'\ ------------------------------------------------------------
    !'\                INERPAD(K),K=1,NPAD  [kg m2]
    !'\                KSTPAD(K),K=1,NPAD  ROT STIFFNESS [Nm/rad]
    !'\                CDAPAD(K),K=1,NPAD  ROT DAMPING   [Nms/rad]

      DO K=1, NPAD
        NPOCKET=NREC(K)
        IF ((NPOCKET > 0).AND.(BEARING.EQ.1)) THEN
            DO J=1,NPOCKET
                READ (2,*) X1r(J,K),Lrec(J,K),PaDorif(J,K)
            END DO
        END IF
      END DO
    !'\ ------------------------------------------------------------
    !'\ IF PAD BEARING WITH POCKETS ONLY
    !'\ ------------------------------------------------------------
    !'\##?? RECESS loc. relative to pad, length of recess, recess orifice diam.

	!READ(finp,*) ZERO,CINLET
	!READ(finp,*) LENGTH,CEXIT

      DO J=1,NJ
        READ (2,*) Z(J),CL(J)
      END DO

      !READ(2,*) ZERO,CINLET
      !READ(2,*) LENGTH,CEXIT

      !READ(2,*) (Z(J),J=1, NJ)
      !READ(2,*) (CL(J),J=1, NJ)

	!'\ ------------------------------------------------------------
	!'\  (Z(J),J=1, NJ) : Axial coordinates for clearances (inlet & exit)
	!'\  (CL(J),J=1, NJ): values of Radial clearance (inlet & exit)
	!'\ ------------------------------------------------------------

      READ(2,*) EXO , EYO, AXO, AYO, ZO, EWX, EWY

	!'\ ------------------------------------------------------------
	!'\  EXO, EYO: rotor eccentricities along X & Y directions [m]
	!'\  AXo, AYO: rotor axis rotations about X & Y axis [rads]
	!'\  ZO: location of center of rotations respect to Z=0 axis [m]
	!'\  EWX, EWY: coordinates for worn seal eccentricity ( NOT WORKING)

      READ (2, *) EMUS,RHOS,EMUA,RHOA,BETA

	!'\ ------------------------------------------------------------
	!'\ EMUS, RHOS: viscosity and density at SUPPLY pressure  [Pa.sec, kg/m3]
	!'\ EMUA, RHOA: """              ""   at DISCHARGE pressure
	!'\ BETA: compressibility coefficient [1/Pa]
	!'/### LSA MUST UPDATE FOR GAS

      READ (2, *) CPS,THS,BETAKS ,Talpha

	!'\ Fluid material properties: CPS,THS,BETAKS ,Talpha
	!'\ ------------------------------------------------------------
	!'/ Talpha (1/degK) = temp_visc_coef
	!'/ for formula: VISC = VIS_atTR EXP(-Talpha(Temp-TR))

      READ (2, *) PC !VERIFICAR RPM,PS,PA

	!'\ ------------------------------------------------------------
	!'\ operating conds: ROTOR RPM, PS: SUPPLY and PA: DISCHARGE Pressures [Pa]
	!'\                             PC: cavitation pressure [Pa]

	!IF (BEARING.EQ.2) THEN !FOR SEAL
	!   READ(finp,*) Pleft, RHOS, EMUS, Cleft
	!ELSE
	!    READ(finp,*) Pleft, RHOA, EMUA, Cleft
	!END IF

      READ (2, *) Pleft, RHOle, EMUle, Cleft

!'\ Left end:  Fluid props and end seal =Pleft, RHOle, EMUle, Cleft

      READ (2, *) Pright, RHOri, EMUri, Cright
	!'\ Right end:  Fluid props and end seal = Pright, RHOright, EMUright, Cright

      READ (2, *) P1PROP,P2PROP,EMU1,EMU2,RHO1,RHO2

	!'\ for FLUID=12 BAROTROPIC: fluid properties at two pressures

      READ (2, *) CD,DORIF,LOSXSI,ALPHA,PRATIO

	!'\ CD: orifice discharge coefficient [empirical],
	!'\ DORIF: orifice diameter [m], reference only
	!'\ LOSXSI: inlet (pressurized) loss coefficient [empirical]
	!'\ ALPHA: Inlet swirl ratio (fraction of rotor speed)
	!'\ PRATIO: guess for inlet pressure ratio IN SEAL and HJB

      READ (2, *) LOSXSIxu,LOSXSIxd,LOSXSIyl,LOSXSIyr,LOSLEAdP

	!'\ LOSXSIxu,LOSXSIxd: pockets loss coefficients circumferential upstream/downstream
	!'\ LOSXSIyl,LOSXSIyr: pockets loss coefficients axial (left and right)
	!'\                    for seal: LOSXSIyl needed
	!'\ LOSLEAdP: loss coefficient for pad bearings only (IFULL=0)
	!'\ Oper conditions, roughness conditions and coefficients

      READ (2, *) TEMPK, VSOUND, RUGR,RUGS, AMOD, BMOD, EXPO

	!'\ ------------------------------------------------------------
	!'\  TEMPK: Suply temperature,
	!'\  VSOUND: Sound Speed [m/s] [NOT important - code evaluates value]
	!'\  RUGR,RUGS: roughness ratios % rotor and stator
	!'\  Moodys friction factor coefficients


!'\  Numerical solution: under-relaxation coefficients & error criteria:
      READ (2, *) ALFU,ALFP,ALFT,MPEPS,PREPS,SFLOW

	!'\ ------------------------------------------------------------
	!'\  ALFU , ALFP, ALFT
	!'\  land mass flow, pressure and orifice flow: MPEPS,PREPS,SFLOW


      ZERO = 1E-20

      If (BEARING.EQ.2) Then
        If (NLOBES.EQ.1) Then
            READ(2,*) ZERO, OFFSET
        Else
            READ(2,*) PRELOAD, OFFSET
        End If

      Else
        READ(2,*) PRELOAD, OFFSET
      End If

	!'\ ------------------------------------------------------------
	!'\  Preload and Offset for seal/bearing
	!'\ ------------------------------------------------------------
	!'\ preload = 0 if seal is cylindrical (360 degrees)
	!'\##?? IFULL=1 ??? BETTER

c ..............................................! non uniform exit pressures


      IF (IFULL.EQ.0) THEN     ! read mixing coefficient
        read (2, *) LEMDA      ! for groove temperature
      END IF                   ! updates

	!'\ ------------------------------------------------------------
	!'\ For pad bearings, print value of thermal mixing coefficient
	!'\                   at pad leading edge (Lambda)
	!'\ ------------------------------------------------------------

c     .................................... bearing compliance coefficients
      read (2,*) AC, ETA, RELAXH    !READ COMPLIANCE params.

c     ....................................! shaft/bearing surface temperatures

      read (2,*) TSHAFT, TSTATOR          !
c     ....................................! params. for radial heat flow
	!'\ ------------------------------------------------------------
	!'\   TEMPERATURES SHAFT and STATOR
	!'\ ------------------------------------------------------------
	!'\ Needed for thermal flow models, in particular with radial heat flow
	! shaft/bearing surface temperatures



      IF (ISOTH.GE.4) THEN
       read (2,*) TBOUT, THERMALK, ROUTER
      END IF
c     ....................................!.................................
	!'\ ------------------------------------------------------------
	!'\  IF ISOTH=4, SET parameters for RADIAL HEAT FLOW:
	!'\  TBOUT (bearing temperature -OUTER surface),
	!'\  THERMALK (bearing thermal conductivity)
	!'\  ROUTER (outer radius)
	!'\ ------------------------------------------------------------


      HCELL=0.0D0                         ! zero [m]
      read (2, *) HCELL                   ! cell depth parameter [m]
c     ......................................................................
	!'\ ------------------------------------------------------------
	!'\ HJB only: SET feed jet angle and location (anglej , locjet)
	!'\ ------------------------------------------------------------

      ANGLEJ=0.0D0                        ! 0 deg. radial injection
      LOCJET=0.5                          ! orifiIce at middle of recess
      IF (BEARING.EQ.1) THEN
       read(2,*) ANGLEJ, LOCJET           ! orifice angle and location
      END IF
c     ....................................!.................................

      IF (BEARING.GE.2) THEN    ! => FOR SEALS AND BEARINGS
         NYI=NLA
         NXT=NLC
      END IF


*      IF (FLUIDTYPE.EQ.5) Then
*        READ(finp,*) gas_constant, compressibility_factor_supply,
*        + compressibility_factor_exit   !verificar parametros
*      End If


      READ(2,*) SPEEDS,FECC

	!'==> NS: NUMBER OF SPEEDS
	!'==> BCASE: IF BCASE = 1 - FIXED eccentricity - IF BCASE = 2 - FIND eccentricity



      DO K=1, SPEEDS

        READ(2,*) RPM(K),PS(K), PA(K), WX(K), WY(K) !WX, WY

      END DO

c   !------------------------------------------------!

      CLOSE (UNIT=2)
      WRITE (6, *) '      Default data read from file DEFAULT.TXT.'

c   !------------------------------------------------!

      IFF=IF
      HRECU=HREC
      FREQU=RPM(1)/60.0
      IF (RPM(1).EQ.(0.0D0))  FREQU=100.0
c ...................................................................

      IF (IF.EQ.12) THEN                ! calculate props for linear fluid
       RHOs=RHO2+(RHO1-RHO2)*(PS(1)-P2prop)/(P1prop-P2prop)
       RHOle=RHO2+(RHO1-RHO2)*(Pleft-P2prop)/(P1prop-P2prop)
       RHOri=RHO2+(RHO1-RHO2)*(Pright-P2prop)/(P1prop-P2prop)
       RHOa=RHO2+(RHO1-RHO2)*(PA(1)-P2prop)/(P1prop-P2prop)

       EMUs=EMU2+(EMU1-EMU2)*(PS(1)-P2prop)/(P1prop-P2prop)
       EMUle=EMU2+(EMU1-EMU2)*(Pleft-P2prop)/(P1prop-P2prop)
       EMUri=EMU2+(EMU1-EMU2)*(Pright-P2prop)/(P1prop-P2prop)
       EMUa=EMU2+(EMU1-EMU2)*(PA(1)-P2prop)/(P1prop-P2prop)
       RETURN
      END IF

      IF (IF.GT.4) RETURN! => FOR AIR, WATER, AND OIL

C .......................!...........................................!
      CALL FDATA(IF)     ! => initialize miprops.f
      CALL CHECKPROPS    ! => FOR CRYOGENIC LIQUIDS IF=1,2,3,4
C .......................!...........................................!

      RETURN

c ======================================================================!

   10 WRITE (6, *) '$ Error accessing file DEFAULT.DAT.'
      WRITE (6, *) '$ Loading internal default data set.'
	TITLE='DEFAULT DATA water HJB'! <60A
      DDATE='01/20/96'
      MODEL=2               ! 2: SINGLE ROW, 1: DOUBLE ROW
      BEARING=1             ! hydrostatic bearing
      NPAD=1                ! # pads on bearing
      NREC(NPAD)=5          ! # recesses/pad
      NPOCKET=5             ! # of pockets/row  on mesh
      NLC=5                 ! # grid points, circumferential on land
      NPC=5           ! ODD ! "   "    "   , circumferential on pocket
      NLA=5                 ! "   "    "   , axial on land (w/o pocket)
      NPA=3                 ! "   "    "   , axial on pocket side
	  NLOBES=1              ! cylindrical HJB
      PRELOADB=0.0D0        ! [m] bearing preload
      OFFSET=0.5D0          ! pad pivot offset
      PRELOAD=PRELOADB      !
      TILTPAD=0             ! 0: fized pads, 1: movingpads.
C                           !......................................
      CINLET=125D-6         ! INLET Radial clearance (m)
      CEXIT=CINLET          ! EXIT radial clearance (m)
      ICSTEP=0              ! no STEP Bearing
      ClearO=CINLET         !
      ClearR=CEXIT          !
      ClearL=CEXIT          !
                            !
      NJ=1                  ! STRAIGHT Clearance
      Z(1)=0.0D0            !
      CL(1)=CINLET          !
      DIAM=0.7645D-1        ! Journal diameter (m)
      LENGTH=0.762D-1       ! Journal length (m)
      LENGTHR=LENGTH/2.0D0  !
      LENGTHL=LENGTHR       !
      AR=0.27026D-1         ! Axial length of recess (m)
      BR=0.27026D-1         ! Circumferential lenght of recess (m)
      X1(1)=-0.010426136    ! Pad Location Circ. [m]
      X1r(1,1)=0.020852     ! 1 recess Location rel to pad [m]
      Lrec(1,1)=BR          ! 1 recess circ. length [m]
      ROTDEL=0.0D0          ! fixed pad
      ROUTER=DIAM/2.0       ! bearing outer radius

      X1(1)=X1(1)*(2.0/DIAM)*Conver         ! in DEG
      LPAD(1)=360.0D0                       ! in DEG
      X1r(1,1)=X1r(1,1)*(2.0/DIAM)*Conver   ! in DEG
      Lrec(1,1)=Lrec(1,1)*(2.0/DIAM)*Conver ! in DEG

      DO J=2, NPOCKET
        X1r(J,1)=X1r(J-1,1)+LPAD(1)/NPOCKET
        Lrec(J,1)=Lrec(1,1)
      END DO                ! J=2, NPOCKET REC Loc & Length

      HCELL=0.0D0           ! zero cell depth (m)
      HREC=254D-6           ! Recess depth (m) only
      HRECU=HREC            !
      VSUP=0.12899E-06      ! Volume of orifice supply line
      YL=-LENGTHL           ! Step Location at left side of HJB
      YR=LENGTHR                !  ''   ''      at right side
C  !........................!......................................
                            ! Equilibrium or initial condition
      EXO=0.0D0             ! S-S or Initial X eccentricity
      EYO=0.0D0             !  "  "     "    Y      "
      AXO=0.0D0             ! Angular misalignment about X axis
      AYO=0.0D0             ! Angular misalignment about Y axis
      ZO=0.0D0              ! Journal Axis CROSSES Y axis [m]
c  !........................!.......................................
      Iwear=0               ! WEAR on bearing 1:YES, 0:NO
      EWX=0.0D0             ! > MIN (Clearance)
      EWY=0.0D0             !
C  !........................! Surface roughness parameters
      RUGR=0.27000E-02      ! Relative to TYP Clearance
      RUGS=RUGR             !

C	                      ! FLUID
C  !........................! Physical properties and conditions
                            ! at PSupply
      EMUS=0.49294E-03      ! Absolute viscosity (N-s/m2)
      RHOS=986.26           ! Density (kg/m3)
      CPS=4182.D0           ! Specific heat in J/kg K (water)
      BETAKS=2.07167D-4     ! Therm. expansion coef.(1/K)
      THS=0.597D0           ! Therm. conductivity (W/m K)
      Talpha=0.0303D0       ! Temperature-viscosity coeff.(1/K-deg)
                            ! at Pdischarge=Pa
      EMUA=0.49294E-03      ! Absolute viscosity (N-s/m2)
      RHOA=982.98           ! Density (kg/m3)
C  !........................!..
      THERMALK=10.0*THS     ! bearing thermal conductivity (W/m K)
                            !
      TEMPK=328.3           ! degK, fluid temperature
      TSHAFT=TEMPK          ! shaft temperature   [K]
      TSTATOR=TEMPK         ! bearing temperature [K]
      TBOUT=TEMPK           ! bearing outer temperature [K]
      VSOUND=1400           ! sonic speed (m/s)
      IL=1                  !
      IF=12                 ! incompressible/linear model
      LEMDA=0.0             ! no thermal mixing at pad grooves

C  !........................! Operating conditions
                            !.....................
      RPM(1)=10200.0D0         ! Rotating speed (rpm)
C
      FREQU=RPM(1)/60.0
      IF (RPM(1).EQ.(0.0D0)) FREQU=100.0
C
      PS(1)=6.877D+6           ! Pressure supply (Pa)
      PA(1)=0.0*6894.757       ! External pressure (Pa)
      PC=PA(1)                 ! Cavitation pressure (Pa)
      BETA=0.48500E-09      ! Compressibility parameter (m2/N)
      CD=0.67               ! Orifice discharge coefficient (Tentative)
      DORIF=.00249          ! Orifice diameter (m)          (Tentative)
      ANGLEJ=45.0D0         ! 0 deg. radial injection
      LOCJET=0.5            ! orifice at middle of recess
      DO J=1, NPOCKET       !
        PaDorif(J,1)=DORIF  !
      END DO                !..................

      LOSXSI=0.0D0          ! Inertia entrance loss coefficient k=(1+Losxsi)/2
      LOSXSIxu=LOSXSI       !
      LOSXSIxd=LOSXSI       !
      LOSXSIyr=LOSXSI       !
      LOSXSIyl=LOSXSI       !
      LOSLeadP=0.64D00      ! PAD: leading edge recovery factor
      ALPHA=0.50D0          ! Swirl fraction at inlet
      PRATIO=0.267D0        ! Pressure ratio for E=0.0
                            !
                            !
c !.........................!.......................................
      ISYM=1                ! Symmetric HJB
      PLeft=PA(1)              ! LEFT Exit Pressure
      PRight=PA(1)             ! RIGHT Exit Pressure
      RHOle=RHOA            ! Props at Left
      EMUle=EMUA            !
      RHOri=RHOA            ! Props at right
      EMUri=EMUA            !.....
      Cleft=0.0D0           ! Exit seal coeffs.
      Cright=0.0D0          !
c !.........................!.......................................
      P1prop=PS(1)             ! props at two different pressures/
      P2prop=PA(1)             !
      EMU1=EMUS             !
      EMU2=EMUA             !
      RHO1=RHOS             !
      RHO2=RHOA             !
c !.........................!.......................................
      IPRuni=1              ! uniform exit pressures
      IPLuni=1              !
C !.........................! Logical Parameters (1: Yes, 0: No)
      INERL=1               ! Fluid inertia effects on lands
      INERP=1               ! Fluid inertia offects on recess edges
      ISOTH=1               ! Isothermal fluid as default process
                            !
C !.........................! Parameters
      ITMAX=99              ! Maximum # of iterations for convergence on lands
      ITPMAX=10             ! Maximum # of iterations for convergence on recess
      ALFU=.90              ! Under_relaxation for moments equations
      ALFP=.60              ! Under_relaxation for pressure correction
      ALFT=.90              ! Under_relaxation for energy equation
      MPEPS=.01             ! Convergence criteria for SUM(Mass source) on lands
      SFLOW=.006            ! Convergence criteria for flow on recess
c !.........................!
      AMOD=0.001375D+0      ! Moodys formulae coefs.
      BMOD=5.0D+5           !
      EXPO=1.0/3.00D0       !
c !.........................!......................................
      DO K=1, MAXNPAD       ! zero data for fixed pads
        ROTPAD(K)=0.0D0     ! pad rotation (rads)
        INERPAD(K)=0.0D0    ! pad inertia  (kg.m2)
        KSTPAD(K)=0.0D0     ! pad rotational stiffness (Nm/rad)
        CDAPAD(K)=0.0D0     ! pad rotational damping (Nm.s/rad)
      END DO                !
c !.........................!......................................
      AC=0.0D0              ! compliance characteristics
      ETA=0.0D0             ! of bearing surface
      RELAXH=1.0D0          !
      LIFT=0                !
c !.........................!......................................
      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Load                                                        **
C **                                                                         **
C **  LOAD:  Given external load calculate equilibrium eccentricity.         **
C **         FOR fixed pad bearings only.                                    **
C **                                                                         **
C *****************************************************************************

C *****************************************************************************
C See SUB Loadtilt on tiltsubs.f for tilt-pad bearing load routine
C *****************************************************************************

      SUBROUTINE LOAD(DEVICE,RPM,PS,PA,WX,WY)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /PARAM2/ EXO, EYO
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAw, IWEAR
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /FREQ/ FREQU,SIGMA,L1,RES,JCASE,NCASE
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT

      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /STIAT/ KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT

      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /IDTHDPAD/ ITPAD
      COMMON /GUESPI/ IGUESP
c     .........................................................
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC, P,
     +                 EXO, EYO, EWX,EWY, EWEAR,BETAw,
     +                 EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 FXT,FYT,TOT,MXT,MYT,QINT,QOUTT,
     +                 PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 FREQU,SIGMA,L1,RES, AC, ETA, PBACK,
     +                 KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,
     +                 KmXYDT,KmYXDT,KXXAT,KYYAT,KXYAT,KYXAT,
     +                 KmXXAT,KmYYAT,KmXYAT,KmYXAT,
     +                 PRELOAD, OFFSET, ROTDEL

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        NPAD, NREC, IFULL, LIFT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        ITER, ITMAX, ITPMAX, SVERB,DVERB, BEEP,
     +        JCASE,NCASE, IWEAR, ISYM, ICSTEP ,ITPAD,IGUESP


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      REAL TIME
      DOUBLE PRECISION EMAX, DEMAX, EPSLOAD, WX, WY, W, FREQUold,
     +                 ECO, EXOLD, EYOLD, ERRX, ERRY, DET,
     +                 DEX, DEY, DECC, ECC, ZERO, DAMP, DUMYR, MAXIN,
     +                 RPM,PS,PA
      INTEGER I, J, ITMAXU, ITPMAU, IFR, ITLOAD, ITLDMAX,
     +        CHO, DEVICE, DUMYI, IFILE, KPAD, IOS
      CHARACTER FILE*80
      CHARACTER YN*1
C ----------------------------------------------------------------------------
C --  LOAD code                                                             --
C ----------------------------------------------------------------------------
      OPEN (UNIT=42, STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)

      NCASE=1                          ! FOR one recess depth
                                       !.........
      ZERO=0.0D0                       !
      EMAX=0.99999999                  ! Maximum allowed eccentricity
      EMAX=EMAX*(1.0D0-PRELOAD/CLEAR)  ! based on preload value

      DEMAX=0.001                      ! Maximum change on Ex & Ey for conv
      EPSLOAD=.010                     ! Convergence load criteria of 1.0 %
      ITLDMAX=30                       ! Maximum number of iterations for load
                                       !
      !WX=0.0D0                         ! Dumy load values
      !WY=0.0D0                         !
      DEX=0.0D0                        !
      DEY=0.0D0                        !
      ITLOAD=0                         !
      DUMYI=0                          ! Dumy OUTPUT Device
      FREQU=RPM/60.0
      FREQUold=FREQU                   !...................!

C...................................................
C Input:
C...................................................

*  100 WRITE (6, 110)
*  110 FORMAT ('$', 'ENTER Load(Wx) in [N] along X dir.: ')
*      CALL ENTERVAL(WX)
*
*      WRITE (6, 120)
*  120 FORMAT ('$', 'ENTER Load(Wy) in [N] along Y dir.: ')
*      CALL ENTERVAL(WY)
*
*      WRITE (6, 130)
*  130 FORMAT ('$', 'INPUT Max # Iters for load calculation (<=10): ')
*      CALL ENTERINT(ITLDMAX)
*
*      IF (ITLDMAX.LT.1) ITLDMAX=10
*      IF (ITLDMAX.GT.10) ITLDMAX=10


      WRITE (6, 67) WX,WY
      IF (DEVICE.EQ.1) THEN
          WRITE (1, 67) WX,WY
      END IF

C
C                                      !.......................
      W=DSQRT(WX*WX+WY*WY)             ! Total load
      IF (W.LT.EPSLOAD) THEN           !.......................
          WRITE (6, *) 'LOAD IS TOO SMALL, TRY AGAIN'
          !GOTO 100                     !
      END IF                           !

      !WRITE (6, 90) WX, WY
      !WRITE (6, *) 'SELECT  one of:'
      !WRITE (6, *) ' (1) Use a guess at Ex=Ey=0 to start soln'
      !WRITE (6, *) ' (2) Use current solution read from DATA file'
      !WRITE (6, *) '     as start guess'
      !write (6, *) '     DATA for guess must be already on memory.'
      !write (6, *) '   '
      !write (6, *) '     ENTER Any number >2, to EXIT'


      CHO=1
      !WRITE (6, 51)
   !51 FORMAT ('$', 'SELECT: ')
      !CALL ENTERINT(CHO)

      !CLOSE(UNIT=42)

      IF ((CHO.LT.1).OR.(CHO.GT.2)) THEN
          RETURN
      ELSE IF (CHO.EQ.1) THEN
          GOTO 400
      ELSE
          GOTO 500
      END IF
C......................................................................
CINPUT                                 !
C                                      ! Read input data from sub. INPUT
  400 EXO=ZERO                         ! for concentric position
      EYO=ZERO                         !......................
      CALL CPARAM(RPM,PS,PA)                      ! Calculate parameters for calc.
      CALL XYDATA(0,RPM,PS,PA)                   ! set and save mesh on PLOTMESH
      CALL ZERO0                       ! zeroes total forces and coeffics.
      CALL ZERO1                       ! ....................

      ITMAXU=ITMAX                     !
      ITPMAU=ITPMAX                    !
      FREQU=ZERO                       ! at zero whirl frequency
c     .................................!..........................
      IGUESP=0                         !
      DUMYR=MAXIN(P)                   ! Determines if P field
      IF (DUMYR.LT.(1.0D-3)) THEN      ! is different from zero.
          IGUESP=1                     ! =1: Call Subroutine GUESP
      END IF            !..............!..........................!

C## NOTE:                              ! Read the following notes!
C## ITPAD=0 means groove edge temperature=Tsupply or Tambient and
C## there is no dependence b/w any two pads. Otherwise, energy balance
C## needs to be considered at the groove and edge temperatures must be
C## updated.
C## Subroutine CALCPADT is designed to update Tedges at Ecc=0.
C##                                    ! Delete following line if needed.
c##   ITPAD=0


c     .................................!..........................
      IF (ITPAD.EQ.1) THEN             ! THD solution of grooved
        CALL CALCPADT(RPM,PS,PA)                  ! pad bearings to get the
      END IF
c     .................................!..........................
      !KPAD=1
c    !................!                !..................!
      DO KPAD=1, NPAD
c    !................!                ! SWEEP on all pads
        CALL XYDATA(KPAD,RPM,PS,PA)
        CALL FILMH(EXO, EYO, NXT, NYI) ! no wear on bearing !###
        CALL SETBC(RPM,PS,PA)                     ! Sets boundary condition


        IF (ITPAD.EQ.1) THEN           ! Get the converged
           CALL READTEMP(KPAD)         ! THD solutions
           CALL FILMC                  !
        ELSE                           !
           IF ((KPAD.EQ.1).AND.(IGUESP.EQ.1)) THEN
             CALL GUESP(DUMYI,RPM,PS,PA)         ! Calculates initial guess for pratio
           END IF
           CALL SOLVE(DUMYI,RPM,PS,PA)           ! SOLVE & calculate:
        END IF                         !

C       CALL PRINTDUVP(DUMYI, 0)       ! forces & torque
        CALL FORCE(DUMYI)              !
        CALL TORQUE(DUMYI)


                 !
        JCASE=1                        !
        IFILE=3                        !
        IGUESP=0                       !###...................
        CALL INPUT1(IFILE,DUMYI,RPM,PS,PA)       ! finds stiffness coeffs.
        CALL WRITETEMP(KPAD)           ! & store in temporal file
        CALL ADD0                      !
        CALL ADD1
                     ! & sums forces & coeffs.
c    !................!                ! K=1,2,..., NPAD
      END DO                           !
c    !................!                ! K=1,2,..., NPAD
c     TIME=SECNDS(0.0)                 ! t=0.0
      GOTO 600                         ! GOTO Start
C                                      !------------------------------------
C......................................!
C Read file with guessed solution:     ! CHO=2

  500 ITMAX=MAX0(ITMAX,99)             !
      ITPMAX=MAX0(ITPMAX,10)           !
      CALL CPARAM(RPM,PS,PA)                      ! Calculate parameters for calc.
      CALL XYDATA(0,RPM,PS,PA)                   ! and generates mesh


      JCASE=1                          !
      FREQU=ZERO                       !
      ITMAXU=ITMAX                     !
      ITPMAU=ITPMAX                    !
      IGUESP=0                         !

      IF ((EXO.GT.ZERO).OR.(EYO.GT.ZERO)) JCASE = 2

      IF (FREQUold.GT.ZERO) THEN       ! => calculate Kij's at w=0
        WRITE(6,*) '$ -----------------------------------------------'
        WRITE(6,*) '$ Calculate start force coeffs. at zero frequency'
        WRITE(6,*) '$ -----------------------------------------------'
        GOTO 520
      END IF

      IF ((KXXDT.GT.ZERO).AND.(KYYDT.GT.ZERO)) GOTO 600

 520  IFILE=1
      CALL ZERO0
      CALL ZERO1
      DO KPAD=1, NPAD
          CALL XYDATA(KPAD,RPM,PS,PA)
          CALL READTEMP(KPAD)
          CALL TEMPROPS
          CALL INPUT1(IFILE, DEVICE,RPM,PS,PA)
          CALL ADD0
          CALL ADD1
          CALL WRITETEMP(KPAD)

      END DO



 550  GOTO 600                         ! GOTO Start


C
C                                      !------------------------------
C--------------------------------------! START of LOAD Calculations  !
C                                      !------------------------------
C                                      !
  600 ITLOAD=ITLOAD+1                 !................................

      !WRITE (6,*) EXO,EYO
      !WRITE (6,*) EXO,EYO
      IF (ITLOAD.GT.ITLDMAX) THEN      ! Maxim. # of iterations exceeded
          GOTO 900                     ! -> Nonconverged
      END IF                           !................................

      ECO=DSQRT(EXO*EXO+EYO*EYO)


      EXOLD=EXO
      EYOLD=EYO

      IFR=0 !............................................................

 615  ERRX=(WX+FXT)                     ! Error in load calculation
      ERRY=(WY+FYT)                     !

          WRITE (6, 80) EXO, WX, FXT, ERRX, EYO, WY, FYT, ERRY
          WRITE (6, 85) KXXDT,KYYDT,KYXDT,KXYDT
      IF (DEVICE.EQ.1) THEN
          WRITE (1, 80) EXO, WX, FXT, ERRX, EYO, WY, FYT, ERRY
          WRITE (1, 85) KXXDT,KYYDT,KYXDT,KXYDT
      END IF



      IF (DABS(ERRX).LT.EPSLOAD*W) THEN ! Check if error in load is
          IFR=IFR+1                     ! less than convergence criteria
      END IF                            !
      IF (DABS(ERRY).LT.EPSLOAD*W) THEN !
          IFR=IFR+1                     !
      END IF                            !....................................

      IF ((CHO.eq.2).and.(IFR.eq.2)) IFR=0

      IF (IFR.GE.2) THEN
          ITLOAD=ITLOAD-1              !..........................
          GOTO 800                     ! Converged TO given load
      END IF                           !.........................

C
C                                         !........................
      DET=KXXDT*KYYDT-KYXDT*KXYDT         ! Determinant of stiffness matrix
      IF (DABS(DET).LT.(0.01D0)) THEN     !........
          GOTO 1200                       ! LOAD TO LARGE
      END IF                              !.........................

      DEX=(KYYDT*ERRX-KXYDT*ERRY)/DET/CLEAR ! Variation in eccentricity
      DEY=(KXXDT*ERRY-KYXDT*ERRX)/DET/CLEAR !
                                       !..........................

  700 EXO=EXOLD+DEX                    ! Improved eccentricy
      EYO=EYOLD+DEY                    ! components Ex & Ey
                                       !..........................

      DECC=DSQRT(DEX*DEX+DEY*DEY)      ! Total change in ECC
      ECC=DSQRT(EXO*EXO+EYO*EYO)       ! & new value of Ecc

      WRITE (6, 70) ITLOAD, EXO, EYO, ECC
      IF (DEVICE.eq.1)  WRITE (1, 70) ITLOAD, EXO, EYO, ECC

C
C
      IF (CHO.EQ.2) GOTO 721           !

      IF (DABS(DECC).LE.DEMAX) THEN    !
          IFR=2                        ! Check for maximum difference
          GOTO 615                     ! in eccentricity between 2 its.
      END IF                           !...............................

 721  CHO=1

      IF ((ECC.GT.EMAX).AND.(AC.EQ.0.0D0)) THEN
         DAMP=(1.0D0+(ECO-1.0D0)*EXP(ECO-ECC)-ECO)/(ECC-ECO)
         DEX=DEX*DAMP                  !
         DEY=DEY*DAMP                  ! DAMP large calculated Ex &Ey
         GOTO 700                      ! for rigid bearing surface
      END IF                           !...............................

      ITMAX=ITMAXU
      ITPMAX=ITPMAU


C## NOTE:                              ! Read the following notes!
C## ITPAD=0 means groove edge temperature=Tsupply or Tambient and
C## there is no dependence b/w any two pads. Otherwise, energy balance
C## needs to be considered at the groove and edge temperatures be updated.
C## Subroutine CALCPADT is designed to update Tedges at Ecc>0.
C##                                    ! Delete following line if needed.

      CALL ZERO0
      CALL ZERO1

c     .................................!..........................
      IF (ITPAD.EQ.1) THEN             ! THD solution of grooved
        CALL CALCPADT(RPM,PS,PA)                  ! pad bearings to get the
      END IF                           ! thermal solution
c     .................................!..........................
                                       !
C   !--------------------!
      DO KPAD=1, NPAD    ! SWEEP OVER ALL PADS:
C   !--------------------!
          CALL XYDATA(KPAD,RPM,PS,PA)
          CALL READTEMP(KPAD)
          CALL FILMH(EXO, EYO, NXT, NYI)

          IF (ITPAD.EQ.0) THEN         ! =1: Solved in Sub. CALCPADT
             CALL SOLVE(DUMYI,RPM,PS,PA)         ! SOLVE & calculate
          ELSE                         !
             CALL FILMC                !
          END IF                       !

          CALL FORCE(DUMYI)
                     !
          CALL TORQUE(DUMYI)
          !WRITE(6,*) 'BEFORE ADD0:',FXT            !
          CALL ADD0                    ! Adds forces & moments
          !WRITE(6,*) 'AFTER ADD0:',FXT
          JCASE=2                      !
          IFILE=3                      !
          CALL INPUT1(IFILE,DUMYI,RPM,PS,PA)     ! finds stiffness coeffs.

          CALL ADD1                    ! Add stiffness coefficients
          CALL WRITETEMP(KPAD)         ! -> temporal storage

C   !--------------------!
        END DO           ! KPAD=1, NPAD
C   !--------------------!


C------------------------!.............!...................................!
         GOTO 600                      ! G0TO Start
C------------------------!.............!...................................!



  800 CONTINUE                               !-----------------------------!
      IF (DEVICE.eq.1) THEN                  ! Convergence achieved
          WRITE (1, 47) ITLOAD               ! ....................
          WRITE (1, 60) EXO, EYO, WX, WY     ! Wx + Fx = 0
      END IF                                 ! Wy + Fy = 0
      WRITE (6, 47) ITLOAD                   !
      WRITE (6, 60) EXO, EYO, WX, WY         !.........................

      !WRITE (61, *) EXO, EYO

      DUMYI=1                                ! Print Table of results
      CALL PRINT0PAD(DEVICE)		     ! .......
      CALL PRINT1PAD(DEVICE,RPM,PS,PA)                 !
      CALL BEEPER                            !=> Calculate force coefs.
      GOTO 1700                              !...at other frequency.
c                                            !.........................


c 900 TIME=SECNDS(TIME)                      !........................
  900 IF (DEVICE.eq.1) THEN                  ! NOT Converged
          WRITE (1, 57) ITLOAD               !........................
          WRITE (1, 60) EXO, EYO, WX, WY     !
      END IF                                 !
      WRITE (6, 57) ITLOAD                   !
      WRITE (6, 60) EXO, EYO, WX, WY         !
 1000 RETURN                                 !........................

 1100 WRITE (6, 10)                          !
      WRITE (6, 20)                          !........................
      WRITE (6, 10)                          ! Eccentricity > EMAX
      IF (DEVICE.EQ.1) THEN                  !
          WRITE (1, 10)                      ! program aborts EXEC
          WRITE (1, 20)                      !
          WRITE (1, 10)                      !
      END IF                                 !
      CALL BEEPER                            !
      CALL PAUSE                             !
      RETURN                                 !.......................

 1200 WRITE (6, 30)                          !......................
      IF (DEVICE.EQ.1) THEN                  ! Determinant of stiffness
          WRITE (1, 30)                      ! coefficients is too small
      END IF                                 ! -> STATICALLY UNSTABLE
      CALL BEEPER                            !    PAD or BEARING
      CALL PAUSE                             !
      RETURN                                 !.........................

c------------------------!...................!.........................!
c CALCULATE FORCE COEFFICIENTS at given frequency
c------------------------!...................!.........................!

 1700 CONTINUE

c    !...........................................................!
c     INPUT FREQUENCY FOR CALCULATION of COEFFICIENTS
c    !................................................!..........!
      OPEN (UNIT=42,  STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)
          YN='Y'
          !WRITE (6, 1520)
 !1520     FORMAT ('$', 'Do you wish to calculate dynamic ',
!     +            'coefficients [Def: Y] (Y/N): ')
     	  !READ (5, 1521) YN
 !1521     FORMAT (1A)
          !YN = 'Y'

          IF ((YN.eq.'y').OR.(YN.eq.'Y')) goto 1522

          GOTO 1235             !==> close unit and exit
c
c     !.........................................!
c     ! CALCULATE COEFFICIENTS at FREQUENCY W
c     !.........................................!
 1522     FREQU=RPM/60.0D0
          IF (RPM.EQ.0.0D0 )FREQU=100.0D0
          !WRITE (6, 1523) FREQU
 !1523 FORMAT ('$', 'ENTER Excitation Frequency(w) in Hz [Default (w)=',
 !    + E12.5E2,']: ')
      !CALL ENTERVAL(FREQU)
      !CLOSE(42)


c    !.................................!
      IF (FREQU.EQ.0.0D0) RETURN       ! load calcs. were made at zero freq.
c    !.................................!

      CALL ZERO0         !
      CALL ZERO1         ! zeroes global force coefficients
      JCASE=2            ! Perturb in X & Y directions
      IFILE=1            ! ............................
C   !--------------------!
      DO KPAD=1, NPAD    ! SWEEP OVER ALL PADS:
C   !--------------------!
       CALL XYDATA(KPAD,RPM,PS,PA)                       ! generates Mesh
       CALL READTEMP(KPAD)                     ! reads input data
       CALL ADD0                               ! Adds forces & moments
       JCASE=2                                 ! ...
       CALL INPUT1(IFILE,DEVICE,RPM,PS,PA)               ! find stiffness coeffs.
       CALL ADD1                               ! Adds stiffness coefficients
       CALL WRITETEMP(KPAD)                    ! -> temporal storage
C   !--------------------!
       END DO            ! KPAD=1, NPAD
C   !--------------------!

C PRINTS TOTAL FORCES AND COEFFICIENTS FOR NPADS:
C -----
C       IF (NPAD.GT.1) THEN
        CALL PRINT0PAD(DEVICE)
        CALL PRINT1PAD(DEVICE,RPM,PS,PA)



C      END IF

C------------------------!.............!...................................!
       RETURN                          ! To MAIN on hydroflex.f
C------------------------!.............!...................................!

C!...........................................!..........................
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) '$ I/O ERROR on SUB: LOAD, RETURN to MAIN.'
      CALL BEEPER
 1235 CLOSE (UNIT=42)
      RETURN

C!...........................................!..........................
C
   10 FORMAT (' ', 79('.'))
   67 FORMAT (' ', 79('.'),/,3X,'SPECIFIED LOAD: WX=',E12.5E2,
     + ' N',2X, 'WY=', E12.5E2,' N',/,1X, 79('.'))
   20 FORMAT (' ', 'SUB LOAD: NEW ECC>EMAX, PROGRAM ABORTS EXECUTION')
   30 FORMAT (' ', 3X, 70('.'), /, ' ',
     +    'DETERMINANT of Kij=0, STATICALLY UNSTABLE BEARING OR PAD',
     +     /, ' ', 70('.'))
   40 FORMAT (' ', 3X, 'LOAD Convergence in Itload:', I3, 3X,
     +        'Time:', F9.2, ' sec.')
   50 FORMAT (' ', 3X, 'NO CONVERGENCE in Itload:', I3, 3X,
     +        'Time:', F9.2, ' sec.')
   47 FORMAT (' ', 79('='),/,3X,
     +'LOAD Convergence in ', I3, ' iterations')
   57 FORMAT (' ', 3X, 'NO CONVERGENCE in Itload:', I3)
   60 FORMAT (' ', 2X, 'Journal EXo:', F7.5,
     +        2X, 'EYo:', F7.5, 2X, 'for Load WX:',
     +        E11.4E2, 2X, 'WY:', E11.4E2, /, ' ', 79('='))
   70 FORMAT (' ', /,/, ' LOAD: Iter=', I3, ' NEW EXo:', F7.5,
     +2X,'EYo:',F7.5,2X,'Ecc:',F7.5,/)
   80 FORMAT (' ',79('.'),/,
     +        ' ',3X,'EX=', F6.4,3X,'WX=',E10.3E2,'+',
     +               'FXT=',E10.3E2,3X,
     +            'errorX=WX+FX',E10.3E2,/,
     +        ' ',3X,'EY=', F6.4,3X,'WY=',E10.3E2,'+',
     +               'FYT=',E10.3E2,3X,
     +            'errorY=WY+FY',E10.3E2)
   85 FORMAT (' ',3X,'KXX=',E11.4E2,2X,'KYY=',E11.4E2,2X,
     +             'KYX=',E11.4E2,2X,'KXY=',E11.4E2,2X,'[N/m]',/,
     +              ' ',79('.'))
   90 FORMAT (' ',39('= '), /, 4X, 'DATA: LOAD WX=', E12.5E2,
     +        '[N],   WY=', E12.5E2, ' [N]', /, 39('='))

 111  FORMAT (' +',77('-'),'+',
     +       /,' |',3X,'RESULTS FOR BEARING PAD #',
     +        I2,46X,' |',/,' |',77('.'),'|')

 112  FORMAT (' +',77('-'),'+',
     +       /,' |',3X,'COEFFICIENTS FOR BEARING PAD #',
     +        I2,' AT FREQ:',E11.4E2,'Hz',
     +        19X,' |',/,' |',77('.'),'|')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Design                                                      **
C **                                                                         **
C **  DESIGN:  Calculate orifice DIAMETERS  given Pratio.                    **
c **           Only for NO WEAR condition        Iwear=0                     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE DESIGN(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK),TREC(MAXNPOCK),QREC(MAXNPOCK),
     +               QIN,QOUT,QFACTOR
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /Pedge/ Pedrise(MAXNPOCK,-MAXNYI:MAXNYI)

      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA!,RPM,PS,PA
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING, /SWITCH/ IPROP

      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /QSIDE0/ QSIDE(MAXNPOCK)
      COMMON /TORREC/ TORR(MAXNPOCK)
c............................................................................

      DOUBLE PRECISION U,V,P,RHOP,EMUP,PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 EXO, EYO, AXO, AYO, ZO, Pedrise,
     +                 DIAORIF, CORIF,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW,     +                 FREQU, SIGMA, L1, RES,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HREC,
     +                 T,CK,BETAK,THC ,HCB,HCJ ,QSIDE,TORR
      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPAD, NREC, IFULL, BEARING,
     +        ITER, ITMAX, ITPMAX, ICASE, NCASE, ISYM,ICSTEP,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IPROP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION QR,AREACD,PI,PR,DELTA,CUORIF,Qright,Denr,Emur
      DOUBLE PRECISION FACTOR ,Cpr,Betakr,Thr
      INTEGER DEVICE,NPock,NXTPock, I, J, K, JR, JL, JS, DUMYI, KSTART

C ----------------------------------------------------------------------------
C --  DESIGN code                                                           --
C ----------------------------------------------------------------------------
C only for Full 2PIXDIAM HJB: 1 Pad with NPOCKET => IFULL =1
C
      PI=DACOS(-1.0D0)
      EXO=0.0D0                 ! For concentric/ not misaligned operation
      EYO=0.0D0                 !
      AXO=0.0D0                 !
      AYO=0.0D0                 !
      DUMYI=0                   ! Dumy Device

      ITMAX=MAX0(ITMAX,99)     !
      ITPMAX=MAX0(ITPMAX,10)   !
c     TIMES=SECNDS(0.0)        !
      NPock=NPOCKET            !SAVE bearing params.
      NXTPock=NXT              !
C..............................!..... On single recess/1pad of 360DEG length
C                                     only need to calculate ONE recess.
      IF ((IFULL.EQ.1).AND.(NPOCKET.GT.0).AND.(ICASE.EQ.1)) THEN
          NPOCKET=1
          NXT=NXI+1
      END IF

C..............................!...........................................
      IPROP=1                  ! for properties calculations
      CALL TEMPROPS            !
      CALL SLAND(DEVICE)       ! Iterative solution of U,V,P on flow region
      ITMAX=ITER               ! ..........................................
      NPOCKET=NPock            !

c..............................!.
      IF ((IFULL.EQ.1). AND.(NPOCKET.GT.1).AND.(ICASE.EQ.1) ) THEN
c..............................!.
      KSTART=ISYM+(ISYM-1)*NYI
      NXT=NXTPock

      DO I=2, NPOCKET          !..............................................
          QREC(I)=QREC(1)      ! ORDER:
          TREC(I)=TREC(1)      !
          TORR(I)=TORR(1)      !
          QSIDE(I)=QSIDE(1)    !
          JL=(I-1)*NXI+2       ! Copy U, V, & P fields to other recesses
          JR=I*NXI+1           !
          JS=1                 !..............................................
          DO J=JL, JR
              JS=JS+1
              DO K= KSTART, NYI
                  P(J, K)=P(JS, K)
                  U(J, K)=U(JS, K)
                  V(J, K)=V(JS, K)
                  T(J, K)=T(JS, K)
                  HCB(J,  K)=HCB(JS,  K)
                  HCJ(J,  K)=HCJ(JS,  K)
                  Rhop(j, k)=Rhop(js, k)
                  Emup(j, k)=Emup(js, k)
                  CK(J,   K)=CK(JS,   K)
                  BETAK(J,K)=BETAK(JS,K)
                  THC(J,  K)=THC(JS,  K)
              END DO
          END DO
          DO K=KSTART, NYI
            Pedrise(I,K)=Pedrise(1,K)
          END DO
      END DO

      QR=0.0D0
      PR=0.0D0
      DO I=1, NPOCKET
          QR=QR+QREC(I)
          PR=PR+PREC(I)
      END DO

      Qright=NPOCKET*QIN       !
      QIN  =QR                 !Inlet & Outlet Bearing flows
      QOUT =NPOCKET*QOUT       !.....
      QR=QR/NPOCKET            ! Average flow rate
      PR=PR/NPOCKET            ! Average pressure at recess

c>>>  UPDATE WITH COMPLIANCE of BEARING rest of land/pocket
      CALL FILMH(EXO,EYO,NXT,NYI) !## update film thickness with
      CALL FILMC                  !## compliance effect
c..............................!
      ELSE
c   !..........................!
        QRIGHT=QIN

c   !..........................!
      END IF                   ! 1 PAD OF 360DEG LENGTH
c..............................!.............................................

 300   IF (NPOCKET.EQ.0) GOTO 400

       FACTOR=CLEAR*CLEAR*CLEAR
       FACTOR=FACTOR*DSQRT(RHO*DABS(PS-PA)/2.0D0)
       FACTOR=EMU*CD*PI/(4.0D0*FACTOR)

c    !.........................!.............................!
      DO I=1, NPOCKET
c    !.........................! Calculate Orifice Diameters
       Pr=PREC(I)              ! Calculate density at Precess
       CALL LOCPROPS(Denr,Emur,Cpr,Betakr,Thr,Pr,TREC(I))


       CORIFS=QREC(I)/DSQRT(DABS(1.0D0-PR)*Denr)  ! Orifice Coefficient
c                                                 !...........................
       DORIF=DSQRT(CORIFS/FACTOR)                 ! Selected Orifice diameter [m]

       DIAORIF(I)=DORIF
       CORIF(I)=CORIFS

       IF ((IFULL.EQ.1).AND.(I.GT.1)) GOTO 350

       WRITE (6, 100) DORIF, CD, CORIFS
       IF (DEVICE.eq.1) WRITE (1, 100) DORIF, CD, CORIFS
c    !.........................!.............................!
 350  END DO
c    !.........................!.............................!

 400  IF (ISYM.EQ.0) QIN=QRIGHT

C........................................!...........................



C........................................!...........................
      PRATIO=PR
C                                !.........................
C#    CALL CPARAM                ! Update Orifice Coefficients.
      CALL PRINTDUVP(DEVICE,0  ) ! Print table of results
      CALL FORCE(DEVICE  )       ! Calculate forces Fx & Fy
      CALL TORQUE(DEVICE  )      ! Calculate torque
      CALL PRINTF(DEVICE,1)      !
C                                !.................. ...! one recess depth
C     ICASE=1                    ! X perturbation only
      DUMYI=3                    ! Ifile=3

      IF ((BEARING.GE.2).AND.(IFULL.EQ.1)) THEN
         ICASE=1
      END IF

C..............................! CALCULATE Dynamic coeffs.
      CALL INPUT1(DUMYI, DEVICE,RPM,PS,PA)
      CALL PRINTCOEF(DEVICE)
C..............................!

C.............................................! PRINT FORMATS

  100 FORMAT (' ','+' 77('-'),'+', /, 1X, 'CALCULATED:',
     +        'Dorif=', E11.4E2, 'm, for Cd=',
     +        F7.5, ' -> Orifice Param. =', E11.4E2, /,
     +        1X, '+' 77('-'),'+')

  200 FORMAT (' ', 3X, 'Time of execution:', E12.5E2, ' secs')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcpr                                                      **
C **                                                                         **
C **  CALCPR:  Calculate Pratio given orifice parameters.                    **
c **           only for NO WEAR condition  Iwear=0                           **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCPR(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK),TREC(MAXNPOCK),QREC(MAXNPOCK),
     +               QIN,QOUT,QFACTOR
      COMMON /Pedge/ Pedrise(MAXNPOCK,-MAXNYI:MAXNYI)
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /BTYPE/ BEARING
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /IDTHDPAD/ ITPAD

      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /QSIDE0/ QSIDE(MAXNPOCK)
      COMMON /TORREC/ TORR(MAXNPOCK)
c............................................................................

      DOUBLE PRECISION U,V,P,RHOP,EMUP,PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 EXO, EYO, AXO, AYO, ZO, Pedrise,
     +                 DIAORIF, CORIF,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 FREQU, SIGMA, L1, RES,
     +                 T,CK,BETAK,THC ,HCB,HCJ ,QSIDE,TORR
      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPAD, NREC, IFULL, BEARING, ITPAD,
     +        ITER, ITMAX, ITPMAX, ICASE, NCASE, ISYM,ICSTEP,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION QR, AREACD, PI, Denr, Emur, Qright, PRAVE,
     +                 Cpr,Betakr,Thr
      INTEGER DEVICE, NPock, NXTPock, I, J, K, JR, JL, JS, DUMYI, KSTART
C ----------------------------------------------------------------------------
C --  CALCPR code                                                           --
C ----------------------------------------------------------------------------
      PI=DACOS(-1.0D0)

      EXO=0.0D0                  ! For concentric/not misaligned operation
      EYO=0.0D0                  !
      AXO=0.0D0                  !
      AYO=0.0D0                  !
      DUMYI=0                    ! Dumy Device=Screen

c     TIMES=SECNDS(0.0)
      ITMAX=MAX0(ITMAX,99)
      ITPMAX=MAX0(ITPMAX,10)

      NPock=NPOCKET           !SAVE input values
      NXTPock=NXT             !
C.............................!...............................................

      IF ((IFULL.EQ.1).AND.(ICASE.EQ.1).AND.(NPOCKET.GT.0)) THEN
          NPOCKET=1                               ! solve just on one rec/land
          NXT=NXI+1                               ! for rotationally symmetric
      END IF                                      ! 360 deg HJB

C                            !................................................
      IF (ITPAD.EQ.0) THEN   ! =1: Solved in Sub. CALCPADT
         CALL SOLVE(DUMYI,RPM,PS,PA)   ! Iterative solution of U, V, & P on flow region.
      END IF                 !
C                            !................................................
      NPOCKET=NPock
      KSTART=ISYM+(ISYM-1)*NYI

c..............................!.
      IF ((IFULL.EQ.1).AND.(ICASE.EQ.1).AND.(NPOCKET.GT.1) ) THEN
c..............................!.
      NXT=NXTPock            ! for 2PI long PAD with NPOCKETs
      PRATIO=PREC(1)         !

      DO I=2, NPOCKET        !........................................
          QREC(I)=QREC(1)    !
          PREC(I)=PRATIO     !
          TREC(I)=TREC(1)      !
          TORR(I)=TORR(1)      !
          QSIDE(I)=QSIDE(1)    !

          JL=(I-1)*NXI+2     ! Copy U, V, & P fields to other recesses
          JR=I*NXI+1         !........................................
          JS=1
          DO J=JL, JR
              JS=JS+1
              DO K= KSTART, NYI
                  P(J, K)=P(JS, K)
                  U(J, K)=U(JS, K)
                  V(J, K)=V(JS, K)
                  T(J, K)=T(JS, K)
                  HCB(J,  K)=HCB(JS,  K)
                  HCJ(J,  K)=HCJ(JS,  K)
                  Rhop(j, k)=Rhop(js, k)
                  Emup(j, k)=Emup(js, k)
                  CK(J,   K)=CK(JS,   K)
                  BETAK(J,K)=BETAK(JS,K)
                  THC(J,  K)=THC(JS,  K)
              END DO
          END DO
          DO K=KSTART, NYI
            Pedrise(I,K)=Pedrise(1,K)
          END DO
      END DO

C   !..........................! Account
      QR=0.0D0                 ! for the total flow for the
      DO I=1, NPOCKET          ! concentric HJB cases
          QR=QR+QREC(I)        !
      END DO                   !
      Qright=NPOCKET*QIN       !
      QIN =QR                  ! Inlet & Outlet Bearing flows
      QOUT=NPOCKET*QOUT        !
      QR=QR/NPOCKET            ! Average flow rate
C   !..........................!

      I=1
      WRITE (6, 100) I, DORIF, CD, PRATIO
      IF (DEVICE.eq.1) WRITE (1, 100) I, DORIF, CD, PRATIO

c>>>  UPDATE WITH COMPLIANCE of BEARING rest of land/pocket
      CALL FILMH(EXO,EYO,NXT,NYI)  !## update film thickness
      CALL FILMC                   !## with compliance effect

C    !.........................!
      ELSE
C    !.........................! PAD Bearing

      IF (NPOCKET.GE.1) THEN

        PRATIO=0.0D0
        DO I=1, NPOCKET
          PRATIO=PRATIO+PREC(I)
           WRITE (6, 100) I, DIAORIF(I), CD, PREC(I)
          IF (DEVICE.eq.1) THEN
           WRITE (1, 100) I, DIAORIF(I), CD, PREC(I)
          END IF
        END DO
        PRATIO=PRATIO/NPOCKET

      END IF

      QRIGHT=QIN
      ICASE=2

C    !.........................!
      END IF                   ! 1PAD OF 360DEG LENGTH
C    !.........................!

C............................!...........................................
c300 TIMES=SECNDS(TIMES)

 300  QR=0.0D0

C     AREACD=(CD*PI*DORIF*DORIF/4.0D0)
      IF (ISYM.EQ.0) QIN=QRIGHT

C............................................................

      WRITE (6, 150)
C.................................!..................................
C                                 !
      CALL PRINTDUVP(DEVICE,0  )  ! Print Results
      CALL FORCE(DEVICE  )        ! Calculates forces Fx, Fy & Torque
      CALL TORQUE(DEVICE )        !...................................
      CALL PRINTF(DEVICE,1)       !

C     ICASE=1                     ! X perturbation only
      NCASE=1                     ! one recess depth
      DUMYI=3                     ! Ifile=3 & calculate DYN Coeffics.
      IF ((BEARING.GE.2).AND.(IFULL.EQ.1)) THEN
         ICASE=1
      END IF

      CALL INPUT1(DUMYI,DEVICE,RPM,PS,PA)  ! Calculates coefficients
      CALL PRINTCOEF(DEVICE)      !

C.................................!..................................


C..................................................... ! PRINT FORMATS

  100 FORMAT (' ', '+' 77('-'),'+', /, 2X, 'REC:',I2,2X,
     +        'Dorif=', E11.4E2, 'm,  Cd=',
     +        F6.4, '==> PRATIO=(Pr-Pa)/(Ps-Pa)=',E12.5E2 )
  150 FORMAT (' ', '+' 77('-'),'+')

  200 FORMAT (' ', 3X, 'Time of execution:', E12.5E2, ' secs')

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Input1                                                      **
C **                                                                         **
C **  INPUT1:  SET parameters for first order solution.                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUT1(IFILE,DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /FREQ/ FREQU, SIGMA, L11, RES, ICASE, NCASE
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR

      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
C    ..................................................................
      DOUBLE PRECISION EXO, EYO, AXO, AYO, ZO,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 FREQU, SIGMA, L11, RES,EWX,EWY,EWEAR,BETAW

      INTEGER ITER, ITMAX, ITPMAX, IWEAR, IFULL,ICASE, NCASE,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        IPRuni, IPLuni, NPRcs, NPLcs, BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      CHARACTER*40  DUMYC
      DOUBLE PRECISION RELAXP,RELAXU,ALFUU,ALFPU,BETUU,
     +                 MPEPSU, EO, DUMY, ZERO
      INTEGER IFILE, DEVICE, I, ITMAXU, IOS, ITMAXP
      REAL TIM, TIME

C ----------------------------------------------------------------------------
C --  INPUT1 code                                                           --
C ----------------------------------------------------------------------------
      ZERO=0.0D0

c...............................................
 7    GOTO (100, 200, 300, 400), IFILE
c...............................................

C !.........................................! IFILE=1
  100 IF (IWEAR.EQ.1) THEN
c    !----------------------------!
        CALL HWEAR(EXO,EYO,NXT,NYI)
c    !.......!
      ELSE
c    !.......!
        CALL FILMH(EXO,EYO,NXT,NYI)
c    !.......!
      END IF
c    !----------------------------!

c >> FOR COMPLIANT BEARINGS update film thickness

      CALL FILMC

c >>  FOR COMPLIANT BEARINGS update film thickness

      GOTO 250

c...........................................!
  200 FREQU=RPM/60.0D0                      ! IFILE=2
      IF (RPM.EQ.zero )   FREQU=100.0D0     ! Set Excitation frequency

  250 IF (ICASE.EQ.2) GOTO 300
      IF ((EXO.NE.zero).OR.(EYO.NE.zero)) THEN
         ICASE=2
      ELSE
         IF ((AXO.NE.zero).OR. (AYO.NE.zero)) ICASE=2
      END IF

c..............................................! IFILE=3
  300 NCASE=1                                  ! JUST ONE FREQUENCY & Ifile=3
      ITMAXU=ITER
      ITMAXP=ITMAX
      MPEPSU=MPEPS
                                               !
c..............................................! IFILE=4
 400  EO=DSQRT(EXO*EXO+EYO*EYO)                ! start first-order solution
      IF (BEARING.EQ.1) THEN  !................!............
        ITMAX=MAX0(ITMAXU,40)                  !RESET CONVERGENCE
        DUMY=0.02D0                            !PARAMS
        MPEPS=DMAX1(DUMY,Eo/20.0D0)            !
      ELSE !...................................!............
        ITMAX=MAX0(ITMAXU,99)
      END IF   !...............................!............
C     HRECU=HREC                               !

      IF (IWEAR.EQ.1) ICASE=2                        ! wear on HJB surface
      IF ((IPRuni.NE.1).OR.(IPLuni.NE.1)) ICASE=2    ! non-uniform exit Press.

c ::::::::::::::::::::::::::::::::             !.....................
      CALL COEFFIC(DEVICE,RPM,PS,PA)                   !calculate dyn coeffs

c ::::::::::::::::::::::::::::::::             !....................

c     TIME=SECNDS(TIME)                        !
      ITMAX=ITMAXP                             !
      ITER=ITMAXU                              ! RESTORE PARAMETERS
      MPEPS=MPEPSU                             !......................

C.......................................... ...! PRINT FORMATS

 1000 FORMAT (' ', 3X, 45('-'), ' Iter:', I5, '  Time:', E12.5E2, ' s')

      END


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


C *****************************************************************************
C **                                                                         **
C **  Subroutine Inpwear                                                     **
C **                                                                         **
C ** INPWEAR Coordinates for WEAR on BEARING according to Scharrer form.     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPWEAR

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAW, IWEAR

      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD,AR,  HREC,
     +                 EWX,EWY,EWEAR,BETAW, CINLET, CEXIT,
     +                 DUMY, CMIN
      INTEGER IWEAR, IOS
      CHARACTER*1 YN
      CHARACTER*40 DUMYC

C ----------------------------------------------------------------------------
C -- INPWEAR code                                                           --
C ----------------------------------------------------------------------------
C Unit 42 is used as temporal storage  for translations between formats.
C.........................................................................

      WRITE (6, *) ' '
      OPEN (UNIT=45,  STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)
C.........................................................................

       IF (IWEAR.EQ.1) THEN
          WRITE (6, 264) 'Y'
       ELSE
          WRITE (6, 264) 'N'
       END IF

  264 FORMAT ('$', 'DO YOU WANT WEAR on BEARING SURFACE,'
     + , '  Enter Y or N [Default ',A1,']: ')

       READ (5, 265) YN
  265  FORMAT (A1)

       IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
           IWEAR=1
       ELSE IF ((YN.EQ.'N').OR.(YN.EQ.'n')) THEN
           IWEAR=0
           Ewx=0.0D0
           Ewy=0.0D0
           CLOSE(UNIT=45)
           RETURN
       END IF

       Iwear=1
       CMIN=CLEAR     !## DMIN1(CINLET,CEXIT)

       WRITE (6, 267) CLEAR, Cinlet, Cexit

  267 FORMAT (' ', 72('.'),/,3X,
     + 'WEAR VECTOR EWear MUST BE LARGER THAN MINIMUM ',
     + ' CLEARANCE in BEARING',/, 3X,
     + ' MINIMUM CLEARANCE =', E12.5E2, '[m]',/,3X,
     + 'AT Inlet, Clearance=', E12.5E2, '[m]',/,3X,
     + 'AT Exit,  Clearance=', E12.5E2, '[m]',/,1X,72('.'))
      CALL BEEPER


  277 continue

      WRITE (6, *) 'EWX (REAL*8):X dir. Wear Coordinate, [m]'
      WRITE (6, 90) EWX

   90 FORMAT ('$', 'ENTER EWX [Default= ',E12.5E2,' m]: ')
      CALL ENTERVAL(EWX)
      WRITE (6, *) ' '

      WRITE (6, *) 'EWY (REAL*8):Y dir. Wear Coordinate, [m]'
      WRITE (6, 91) EWY

   91 FORMAT ('$', 'ENTER EWY [Default= ',E12.5E2,' m]: ')
      CALL ENTERVAL(EWY)
      WRITE (6, *) ' '


      EWEAR=DSQRT(EWX*EWX+EWY*EWY)


      IF (EWEAR.LT.CMIN) THEN

         write (6,92) EWX,EWY,EWEAR
         CALL Beeper
         EWX=CMIN
         EWY=0.0D0
         GOTO 277

      END IF

      CLOSE (UNIT=45)


c..........................................................

   92 FORMAT (' ',72('.'),/,
     + ' ERROR: Magnitude of WEAR Vector is SMALLER THAN',
     + ' Minimum Clearance, ',/,1X,
     + 'Ewx=', E12.5E2,' Ewy=',E12.5E2, '-> Ewear=', E12.5E2,
     + ' [m] < Clearance',/,1X,
     + 'REVISE YOUR INPUT DATA/SEE MANUAL',/,1X, 72('.'))


 1000 FORMAT (40A)
 1010 FORMAT (' ',40A)

      RETURN

C---------------------------------------------------------------
C Error in I/O to unit 45.
C
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) 'ERROR Detected: RETURN TO MENU'
      Write (6, *) '::::::::::::::::::::::::::::::'
      CALL BEEPER
      CALL PAUSE
      CLOSE (UNIT=45)
      RETURN

      END
c
c
C *****************************************************************************
C **                                                                         **
C **  SUBROUTINE MinMaxP                                                     **
C **                                                                         **
C **  MINMAXP: Find minimum and maximum pressures in field                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE MINMAXP

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DEN2,VIS2,PSA,PA,
     +                 DEN12P12, VIS12P12, P2
      COMMON /COMPRE/ L4, PCAV
      COMMON /PMINMAX/ PMIN, PMAX

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION P, RHOTYP,EMUTYP,DEN2,VIS2,PSA,PA,
     +                 DEN12P12, VIS12P12, P2, L4, PCAV, PMIN, PMAX
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER I, J, IStart, NN, K, BEARING


C ----------------------------------------------------------------------------
C --  MINMAXP code                                                          --
C ----------------------------------------------------------------------------
      Istart=ISYM+(ISYM-1)*NYI
      NN=NXT
      IF (BEARING.EQ.2) Istart=1    ! for annular seals

 17   PMIN=1.0D+30
      PMAX=-PMIN
      DO I=1, NN
          DO J= Istart, NYI, 1
              IF (P(I,J).GT.PMAX) THEN
                  PMAX=P(I,J)
              ELSE IF (P(I,J).LT.PMIN) THEN
                  PMIN=P(I,J)
              END IF
          END DO
      END DO

      IF (PMIN.LT.PCAV) PMIN=PCAV

      PMAX=PMAX*PSA+PA      !PA  ! max. pressure in MPA
      PMIN=PMIN*PSA+PA      !PA  ! min. pressure in MPA


      END

C *****************************************************************************
C **                                                                         **
C **  SUBROUTINE MinMaxT                                                     **
C **                                                                         **
C **  MINMAXT: Find minimum and maximum temperatures in field                **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE MINMAXT(TMIN,TMAX,TEA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION T
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION TMIN,TMAX,TEA
      INTEGER I, J, IDUM, IStart, NN, K, BEARING

C ----------------------------------------------------------------------------
C --  MINMAXT code                                                          --
C ----------------------------------------------------------------------------
      Istart=ISYM+(ISYM-1)*NYI
      NN=NXT
      IDUM=ISYM
      IF (BEARING.EQ.2) THEN  !=> for annular seals
          Istart=1
          IDUM=1
      END IF

 17   TMIN=1.0D+30
      TMAX=-TMIN
      DO I=1, NN
          DO J= Istart, NYI, 1
              IF (T(I,J).GT.TMAX) THEN
                  TMAX=T(I,J)
              ELSE IF (T(I,J).LT.TMIN) THEN
                  TMIN=T(I,J)
              END IF
          END DO
      END DO

c    Average temperature at exit plane Y=Lr and -Ll

        TEA=0.0D0
        DO I=1, NXT
           TEA=TEA+T(I,NYI)+(1-IDUM)*T(I,-NYI)
        END DO
        TEA=TEA/(NXT)/(2-IDUM)

      END

C ----------------------------------------------------------------------------
C LAST REVISED 2/18/94 by Dr. Luis San Andres texas A&M University
C 6/20/95  Read/Write HCELL to DEFAULT.DAT file
C 7//3/95  READ/WRITE LOCJET.ANGLEJET to DEFAULT.DAT file
C ----------------------------------------------------------------------------