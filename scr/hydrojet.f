C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LAST revised 01/02/96 - 9/12/95  by Dr. Luis San Andres
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    #    #   #   #  #####   #####    ####        #  ######   #####
C    #    #    # #   #    #  #    #  #    #       #  #          #
C    ######     #    #    #  #    #  #    #       #  #####      #
C    #    #     #    #    #  #####   #    #       #  #          #
C    #    #     #    #    #  #   #   #    #  #    #  #          #
C    #    #     #    #####   #    #   ####    ####   ######     #
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    MAIN BODY OF PROGRAM HYDROJET
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C
C  hydrojet.f Copyright Dr. Luis San Andres TexasA&MUniversity / 1995
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR III
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

      PROGRAM HYDROJET

C
c               C......................................C
c               C  TEXAS A&M UNIVERSITY                C
c               C MECHANICAL ENGINEERING DEPARTMENT    C
c               C......................................C
c               C ANALYSIS & PROGRAM BY                C
c               C DR. LUIS SANANDRES, Associate Prof.  C
c               C 12/31/95                             C
c               C......................................C
c
c
c
c--------------------------------------------------------------------
c hydrojet has the following fortran programs:
c .........
c       hydrojet.f
c	addcoefs.f     compe2t.f      guespp.f       params.f       supportp.f
c	calcmesht.f    fileplot.f                    pjet.f         tiltsubst.f
c	calcpadt.f     filmwt.f       initopst.f     prntrest.f     turbcoefN.f
c	calcsolnt.f    findprop.f     inputp.f       readwritet.f   turbcoefs.f
c	calcthert.f    firstt.f       mipropst.f     setpexit.f
c	compe1t.f      foilsubs.f     optionst.f     spline.f
c
c and DATA files: DEFAULT.DAT, HJBHELP.DAT
c                 O2.COF, PH2.COF, N2.COF, METH.COF
c--------------------------------------------------------------------c
C..............................................................
C  >>> THERMAL CASES are:
C     ISOTH.EQ. 1) 'Isothermal fluid film (T=Ts=Constant)'
C     ISOTH.EQ. 2) 'Adiabatic journal & Isothermal stator (Tb)'
C     ISOTH.EQ. 3) 'Isothermal journal (Tj) & Adiabatic stator'
C     ISOTH.EQ. 4) 'Adiabatic journal & Bearing radial heat flow'
C     ISOTH.EQ. 5) 'Isothermal journal (Tj) & Bearing radial heat flow'
C     ISOTH.EQ. 0) 'Isothermal journal (Tj) & stator: (Tb)'
C     ISOTH.EQ.-1) 'Adiabatic bounding Surfaces (Qb=Qj=0)'
C..............................................................


      IMPLICIT NONE

      INCLUDE 'params.f'
C
C params.f contains parameters for dimensioning of all arrays
C          and it is used by most *.f programs.

C
C...............................................................
C Character variables
C...............................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
C................................................................
C Double precision variables
C................................................................
C
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,
     +                 EXO,EYO,AXO,AYO,ZO,EWX,EWY,EWEAR,BETAW,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr,
     +                 CINLET, CEXIT,Pleft, Pright, Cleft, Cright,
     +                 Ple, Pri, Csel, Cser,LENGTHR,LENGTHL,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 CPS,THS,BETAKS, TRECL,TRECR,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2,AC, ETA, RELAXH,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 LOSleadP, KLOSpad, DYVK, DYPK, SVNK, SVSK,
     +                 PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 AMOD, BMOD, RUGR, RUGS, EXPO,
     +                 PRECL, PRECR,ClearO,ClearR,ClearL,YR,YL,
     +                 FFACTOR, MFACTOR, FX, FY, MX, MY, TO, TOR,
     +                 FXT,FYT,TOT,MXT,MYT,QINT,QOUTT,
     +                 YB, YT, YTIC, XL, XR, XTIC, LEMDA, DELTA,
     +                 QIN, QOUT, QFACTOR,  PMIN,PMAX,
     +                 PRELOAD, OFFSET, ROTDEL, PRELOADB,
     +                 RSINPK, RCOSPK, IPAD,KROTPAD, CROTPAD,
     +                 ANGLEJ, LOCJET, CJET, DPJET


C
C................................................................
C Double Precision arrays
C................................................................
C
      DOUBLE PRECISION XP, XU, YP, YV, DXP, DXU, SUW, SUE,
     +                 HP, HU, HV,  HPO, HUO, HVO,
     +                 Z,CL,BCL,CCL,DCL, DYP, DYV, SVN, SVS,
     +                 COSXP, SINXP, COSXU, SINXU,
     +                 U, V, P, Uor, Vor, DIAORIF,CORIF,
     +                 PROLD, QOLD, POLD, UOLD, VOLD,
     +                 DU, DV, DVV, A, B, C, D,
     +                 PREC, QREC,TREC, ASPE, L4R,
     +                 PRCOEFC, PRCOEFS,PLCOEFC,PLCOEFS,
     +                 X1,LPad, X1r,Lrec,PADorif, ROTPAD,
     +                 INERPAD, KSTPAD, CDAPAD,
     +                 PLEAD,PTRAIL,QLEAD,QTRAIL,TQTRAIL,TLEAD,
     +                 TRAIL,QSIDE,QSIDEH,TORR,TROLD,TOLD,
     +                 PRECdo,PRECup,PRjet,WX,WY
C
C................................................................
C Integer variables
C................................................................
C
      INTEGER IWEAR,NJ, ISYM, ICSTEP, NPAD,  NREC, TILTPAD,
     +        JUMIN, JUMAX, JUSTART, JUSTOP, LEFTBC, TILT,
     +        JVMIN, JVMAX, JVSTART, JVSTOP, RIGHTBC,
     +        JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, LIFT,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        K, KM1, KP1, ITER, ITMAX, ITPMAX, IPROP,
     +        SVERB, DVERB, BEEP, ICHEC1,ICHEC2, KPAD,NLOBES,
     +        IPRuni, IPLuni, NPRcs, NPLcs,IMOMflag, BEARING,
     +        SPEEDS,FECC,WWW

C
C................................................................
C   COMMON BLOCKS
C................................................................
C Character commons.
C................................................................
C
      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
C
C................................................................
C Double precision array commons.
C................................................................
C
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD,SPEEDS,FECC
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)
      COMMON /PIOPAD/ PLEAD(-MAXNYI:MAXNYI), PTRAIL(-MAXNYI:MAXNYI)
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HOFILM/ HPO(MAXNXT,-MAXNYI:MAXNYI),
     +                HUO(MAXNXT,-MAXNYI:MAXNYI),
     +                HVO(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /SPLDATA/Z(NSL),CL(NSL),BCL(NSL),CCL(NSL),DCL(NSL), NJ
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)

      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /PDUMY/ PROLD(MAXNPOCK), TROLD(MAXNPOCK), QOLD(MAXNPOCK),
     +               POLD(MAXNXT,-MAXNYI:MAXNYI),
     +               UOLD(MAXNXT,-MAXNYI:MAXNYI),
     +               VOLD(MAXNXT,-MAXNYI:MAXNYI),
     +               TOLD(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DPUV/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /QSIDE0/ QSIDE(MAXNPOCK)
      COMMON /QSIDE0H/ QSIDEH(MAXNPOCK)
      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /RECASP/ ASPE(MAXNPOCK)
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)
      COMMON /RECCOM/ L4R(MAXNPOCK)
C
C................................................................
C Integer commons.
C................................................................
C
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /PADK/ KPAD
      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /FLMOM / IMOMFLAG
      COMMON /BTYPE/ BEARING
      COMMON /GUESPI/ IGUESP
      COMMON /SWITCH/ IPROP
C
C................................................................
C Double precision commons.
C................................................................
C
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,
     +         RPM(MAXNXT),PS(MAXNXT),PA(MAXNXT),WX(MAXNXT),WY(MAXNXT)
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr
      COMMON /LOSPAD/ LOSleadP, KLOSpad
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PRELR/ Ple, Pri, Csel, Cser
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR
      COMMON /FORCE0/ FFACTOR, FX, FY, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /YMAXMIN/ YB, YT, YTIC
      COMMON /XMAXMIN/ XL, XR, XTIC
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /TGROOVE/ LEMDA,DELTA
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET


c...............................................................
c commons for liquid properties > INTERFACE with MIPROPS
c...............................................................
      DOUBLE PRECISION G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,
     +                 TUL,TLL,PUL,DCC, TEMPK, VSOUND,
     +                 EM,EOK,RM,TC,DC,X,PCo,SIG, BX,PX,
     +                 CP, RC , AIRMU, AGAMA, PR ,KS, Talpha
      INTEGER IF,IFF,IL

      DIMENSION G(32),VP(9)	!==> arrays for miprops

      COMMON/LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /AIR/ CP, RC , AIRMU, AGAMA, PR ,KS
      COMMON /OILCOEF/ Talpha
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/CONT/IFF
      COMMON/CRIT/EM,EOK,RM,TC,DC,X,PCo,SIG
      COMMON/DIEL/BX(6),PX(6)


C................................................................
C Commons for FIRST  order solution
C................................................................
      INTEGER ICASE, NCASE, ISOLN
C................................................................
      DOUBLE PRECISION APU, APUI, AWU, AEU, ASU, ANU, GUO, GUV,
     +                 APV, APVI, AWV, AEV, ASV, ANV, GVO, GVU,
     +                 GUPR, GUPI, GUT, GVPR, GVPI, GVT,
     +                 RHOP,EMUP,Drhop, Demup,Drhot, Demut,
     +                 QO1, QINO1, QOUTO1, Soo,
     +                 FREQU, SIGMA, L1, RES,
     +                 HREC1, VSUP, BETA,L4, Pcavi
      DOUBLE PRECISION K11, K12, C11, C12, KM11,KM12,CM11,CM12,
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION
     +         KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT,
     +         CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT,
     +         MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT,
     +         KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT,
     +         CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT,
     +         MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

C...............................................................
      DOUBLE COMPLEX U1, V1, P1, DU1, DV1, DVV1,
     +               A1, B1, C1, D1, PREC1,TREC1,QREC1,QIN1,QOUT1,
     +               PRECL1, PRECR1,S11, Sro, Srec,
     +               FXX, FYY, MXX, MYY, QZERO, QRECREC,
     +               FXXO,FYYO,MXXO,MYYO,QZEROO,QRECRECO,
     +               XFXX, XFYY, XMXX, XMYY, QXZERO,
     +               XFXXO,XFYYO,XMXXO,XMYYO,QXZEROO,
     +               YFXX, YFYY, YMXX, YMYY, QYZERO,
     +               YFXXO,YFYYO,YMXXO,YMYYO,QYZEROO
C
C................................................................
C COMMON BLOCKS for FIRST Order solution
C................................................................
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /SOLN/ ISOLN
      COMMON /COMPRE/ L4, Pcavi
C................................................................
      COMMON /UCOEF/
     + APU(MAXNXT,-MAXNYI:MAXNYI), AWU(MAXNXT,-MAXNYI:MAXNYI),
     + AEU(MAXNXT,-MAXNYI:MAXNYI), ASU(MAXNXT,-MAXNYI:MAXNYI),
     + ANU(MAXNXT,-MAXNYI:MAXNYI), GUO(MAXNXT,-MAXNYI:MAXNYI),
     + GUV(MAXNXT,-MAXNYI:MAXNYI), APUI(MAXNXT,-MAXNYI:MAXNYI),
     + GUPR(MAXNXT,-MAXNYI:MAXNYI),GUPI(MAXNXT,-MAXNYI:MAXNYI),
     + GUT(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /VCOEF/
     + APV(MAXNXT,-MAXNYI:MAXNYI),AWV(MAXNXT,-MAXNYI:MAXNYI),
     + AEV(MAXNXT,-MAXNYI:MAXNYI),ASV(MAXNXT,-MAXNYI:MAXNYI),
     + ANV(MAXNXT,-MAXNYI:MAXNYI),GVO(MAXNXT,-MAXNYI:MAXNYI),
     + GVU(MAXNXT,-MAXNYI:MAXNYI),APVI(MAXNXT,-MAXNYI:MAXNYI),
     + GVPR(MAXNXT,-MAXNYI:MAXNYI),GVPI(MAXNXT,-MAXNYI:MAXNYI),
     + GVT(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /TCOEF/
     + GTU(MAXNXT, -MAXNYI:MAXNYI), GTV(MAXNXT, -MAXNYI:MAXNYI),
     + GTP(MAXNXT, -MAXNYI:MAXNYI), GTH(MAXNXT, -MAXNYI:MAXNYI),
     + GTPI(MAXNXT,-MAXNYI:MAXNYI), APTI(MAXNXT,-MAXNYI:MAXNYI),
     + APT(MAXNXT, -MAXNYI:MAXNYI), AWT(MAXNXT, -MAXNYI:MAXNYI),
     + AET(MAXNXT, -MAXNYI:MAXNYI), AST(MAXNXT, -MAXNYI:MAXNYI),
     + ANT(MAXNXT, -MAXNYI:MAXNYI), BT1(MAXNXT, -MAXNYI:MAXNYI),
     + BT2(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FLOW01/ QO1(MAXNPOCK), QINO1, QOUTO1
      COMMON /So/ Soo(MAXNPOCK)
      COMMON /RECPAR/ HREC1, VSUP, BETA
      COMMON /COEFS/ K11, K12, C11, C12, KM11,KM12,CM11,CM12
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /DAMPT/ CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT
      COMMON /INERT/ MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT
      COMMON /STIAT/ KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT
      COMMON /DAMAT/ CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT
      COMMON /INEAT/ MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
C....................................................................
C Common blocks for complex data arrays
C....................................................................
      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /DPUV1/ DU1(MAXNXTP2), DV1(MAXNXTP2), DVV1(MAXNXTP2)
      COMMON /TDMA1/ A1(MAXNXTP2),B1(MAXNXTP2),C1(MAXNXTP2),
     +               D1(MAXNXTP2)
      COMMON /RECES1/ PREC1(MAXNPOCK),TREC1(MAXNPOCK),
     +                QREC1(MAXNPOCK),QIN1, QOUT1
      COMMON /PRECES1/ PRECL1, PRECR1
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
      COMMON /S1/ S11(MAXNPOCK)
      COMMON /Sxy/ Sro(MAXNPOCK), Srec(MAXNPOCK,MAXNPOCK)

C................................................................
C      Variables for thermal analysis
C................................................................
C
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd
      COMMON /IDTHDPAD/ ITPAD
      COMMON /DBTA/ DBKP(MAXNXT,-MAXNYI:MAXNYI),
     +              DBKT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DCPS/ DCPP(MAXNXT,-MAXNYI:MAXNYI),
     +              DCPT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DTHC/  DKP(MAXNXT,-MAXNYI:MAXNYI),
     +               DKT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ TK(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /THERMAL/ ALFT, UC, TCT, Ec
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /TISOBJD/ TBJ,TBS
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /TMAXMIN/ TKMAX,TKMIN

      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /THERMID/ ISOTH

      DOUBLE PRECISION DBKP, DBKT ,DCPP, DCPT ,DKP,  DKT,
     +                 CK,BETAK,HCB,HCJ,THC,TK, TKMAX,TKMIN,
     +                 ALFT, UC, TCT, Ec,
     +                 GTU, GTV,  GTP, GTH, GTPI,
     +                 APT, APTI, AWT, AET, AST, ANT, BT1, BT2,
     +                 TSUMP,TSUMPd,RHOSU,RHOSd,
     +                 TSHAFT, TSTATOR, TBJ, TBS,
     +                 TBOUT, THERMALK, ROUTER, HKB
      DOUBLE COMPLEX   T1
      INTEGER ISOTH , ITPAD ,IGUESP
C
C................................................................
C
C
C
C
C................................................................
C LOCAL Variables for main program
C................................................................
C
      CHARACTER*80 FILE, DUMPF, SCRATCH
      CHARACTER*1 YN
      CHARACTER*40 DUMYC
      DOUBLE PRECISION DUMYR,MAXIN
      INTEGER DEVICE,CHOICE, I,J,DUMYI,DUMP,IOS,ISYMold,LOSTAT
C
C......................................................................
C Initial definitions.
C......................................................................
C Moody's coefficients Amod, Bmod, Expo > TURBULENCE MODEL DEFAULT
      AMOD=0.001375D+00
      BMOD=5.00D+05
      EXPO=1.0D0/3.00D0
C.....................................................................
C Data constants for air
      CP=1004.D0            ! Specific heat (J.kg/degK)
      RC=287.D0             ! Gas constant [m3/s2/K]
      AIRMU=1.79D-5         ! viscosity standard at T=288.2K [N/m2.s]
      AGAMA=1.4D0           ! Ratio of specific heats
      PR=0.71D0             ! Prandtl Number
      KS=0.026D0            ! Heat conductivity (W/m.K)
C.....................................................................
      LEMDA=0.0D0           ! mixing temperature pad coefficient
C.....................................................................
      IPROP=1                            ! Yes, update properties
      SVERB=0                            ! 0: Low verbosity for screen
      DVERB=0                            ! Low verbosity for dump file
      BEEP=0                             ! Beep turned on
      NCASE=1                            ! 1 recess depth only
C......................................................................
C     INITIALIZE MAIN ARRAYS
C    !======================!
      DO I=1, MAXNXT
          DO J= -MAXNYI,MAXNYI,1         ! TEMPORARY ONLY >>>>>>>>>
              U(I, J)=0.0D0              ! circumferential velocity
              V(I, J)=0.0D0              ! axial velocity
              P(I, J)=0.0D0              ! fluid pressure field
              TK(I,J)=1.0D0              ! fluid film temperature
              HCB(I,J)=0.D0              ! Heat transfer coeff. to bearing
              HCJ(I,J)=0.D0              ! Heat transfer coeff. to journal
              Drhop(I,J)=0.D0            ! dRHO/dP
              T1(I, J)=(0.0D0,0.0D0)     ! First-order fluid film temp.
          END DO                         ! TEMPORARY ONLY >>>>>>>>>
      END DO

      DO I=1, MAXNPOCK                   ! On MAX # of pockets
          PREC(I)=0.0D0                  ! Recess pressures
          QREC(I)=0.0D0                  ! Recess flow rates
          TREC(I)=1.0D0                  ! Recess fluid temperature
          TORR(I)=0.0D0                  ! Torque over recess area
          TREC1(I)=(0.0D0,0.0D0)         ! 1st-order recess fluid temp.
          DO J= 1,NPC            	 !
            PRjet(I,J)=0.0D0             ! recess rise circ. pressures
          END DO			 !
      END DO                             !

      DO I=1, MAXNPAD                    ! On MAX # of bearing pads
         DO J=-MAXNYI,MAXNYI,1           !
             TLEAD(I,J)=1.0D0            ! Pad leading edge temperature
             TRAIL(I,J)=1.0D0            ! Pad trailing `'  `'
         END DO
         QTRAIL(I) =0.0D0                ! Flow rate at trailing edge
         TQTRAIL(I)=0.0D0                ! Flow rate x Temp at `'  `'
      END DO

      DO I=1, MAXNYI 			 ! Coefficients of exit (discharge)
       PRCOEFC(I)=0.0D0			 ! non-uniform pressures
       PRCOEFS(I)=0.0D0
       PLCOEFC(I)=0.0D0
       PLCOEFS(I)=0.0D0
      END DO
       PRCOEFC(0)=0.0D0
       PLCOEFC(0)=0.0D0
c !.........................!......................................
      DO K=1, MAXNPAD       ! zero data for fixed pads
        ROTPAD(K)=0.0D0     ! no pad rotation [rads]
        INERPAD(K)=0.0D0    ! null pad inertia [kg m2]
        KSTPAD(K)=0.0D0     ! null pad rotational stiffness [Nm/rad]
        CDAPAD(K)=0.0D0     ! null pad rotational damping [Nms/rad]
      END DO
c !.........................!......................................
C
C................................................................
C Print opening, read default data and ask for dump file.
C................................................................

      write (6, 175)
  175 format(/,/,/,/)
      write (6, 171)
  171 format (4X,60('-'))
      write (6, 172)
  172 format (4X,'|', 24X,'HYDROJET  ',24X,'|',/,4X,'|',58X,'|')
      write (6, 171)


      CALL DISCLOSURE

      !CALL BEEPER


      write (6, 175)
      write (6, 175)
      write (6, 171)
      write (6, 172)
      write (6, 171)
      write (6,174)
  174 format (5X,
     +'Numerical Soln. of Pressure, Flow, Film Forces and Force',/,5X,
     +'Coefficients for ANGLED INJECTION HYDROSTATIC BEARINGS'  ,/,5X,
     +'DAMPER SEALS, (flex) TILT PAD and SIMPLE FOIL BEARINGS',/,/,5X,
     +'ADIABATIC SU00000000000RFACES or CONSTANT BEARING',1X,
     +'& JOURNAL TEMPERATURES,',/,5X,
     +'Uses MIPROPS for LIQUID Oxygen, Hydrogen & Nitrogen')
      write (6, 171)
      write (6, 175)



C........................................................................
C READ INPUT DATA from DEFAULT.DAT file or internal set
C........................................................................

 177  CALL INPUT     ! Input default data from DEFAULT.DAT or internal set.
      CALL CPARAM(RPM(1),PS(1),PA(1))    ! Calculate parameters for default data.
      CALL XYDATA(DEVICE,RPM(1),PS(1),PA(1)) ! 0: creates mesh for DEFAULT.DAT
      CALL PRINPUT(DEVICE,RPM(1),PS(1),PA(1)) !Print Parameters


        OPEN(UNIT=61, FILE='RESULTS.txt',
     +   STATUS='UNKNOWN', ACTION='WRITE')

      DO WWW=1,SPEEDS
        !WRITE (61,*) W,SPEEDS
        WRITE (6,250)
        WRITE (6,260) RPM(WWW), RPM(WWW)/60,PS(WWW),PA(WWW),
     +            WX(WWW), WY(WWW)

        IF (TILTPAD.EQ.1) THEN
            CALL LOADTILT(DEVICE,RPM(WWW),PS(WWW), PA(WWW), WX(WWW),
     +                    WY(WWW))
        ELSE

            CALL LOAD(DEVICE,RPM(WWW),PS(WWW), PA(WWW), WX(WWW),
     +             WY(WWW))
        END IF



        IF (ISOTH.EQ.1) THEN
            TKMAX = TEMPK
            TKMIN = TEMPK
        END IF



        WRITE (61, *) RPM(WWW), EXO, EYO,FXT,FYT,
     +                TOT,QINT,PMAX/1E5,PMIN/1E5,
     +                KXXD,KYYD,KXYD,KYXD,
     +                CXXD,CYYD,CXYD,CYXD,
     +                MXXD,MYYD,MXYD,MYXD,
     +                TKMAX,TKMIN


      !IF (ISOTH.NQ.1) THEN
      !      WRITE (6, *) HJ()
      !END IF

        !WRITE (6, *) RPM(W), KXXD, KYXD, KYYD, KXYD
        !WRITE (6, *) CXXD, CYXD, CYYD, CXYD, EXO, EYO
        !WRITE (6, *) MXXD, MYXD, MYYD, MXYD!, KEQ, WHIRL
        !WRITE (6, *) FXT,FYT
        !WRITE (6, *) FX,FY
        !WRITE (6, *) PMIN,PMAX
        !WRITE (6, *) TLEAD,TRAIL
        !WRITE (6, *) PMIN,PMAX
        !WRITE (6, *) TW, TE


        !WRITE (6,250)
        !WRITE (61,*) W,SPEEDS
        !STOP

      END DO
      CLOSE(UNIT=61)

      WRITE(6,255)
      WRITE(6,*) 'EXIT hydrojet'
      WRITE(6,255)
        !PAUSE
      STOP
      !CLOSE
c INPUT ==> options.f ; CPARAM ==> initops.f ; XYDATA ==> calcmesh.f
 250  FORMAT (' ', 3X, 120('-'))
 255  FORMAT ('', 3X, 120('-'))
 260  FORMAT(' ',3X,'BCASE=1',' RPM:',F5.0, ' ',
     +   'FREQ:',F5.2, '[Hz]',' PS,PA:', E12.3E2,' ', E12.3E2, '[N/m2]',
     +    ' WX,WY:' ,E12.3E2, E12.3E2, '[N]')

C........................................................................
C CREATES OUTPUT DUMP FILE
C........................................................................

      CALL BEEPER

 1010 WRITE (6, *) ' '
      WRITE (6, 10)
   10 FORMAT ('$', 'Do you wish to have output dumped to file? (Y/N) ')
      READ (5, 20, IOSTAT=IOS, ERR=1000) YN
   20 FORMAT (A)
      IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
 1110     WRITE (6, 30)
   30     FORMAT ('$', 'Dump file name ? SPECIFY TXT extension:')
          READ (5, 40, IOSTAT=IOS, ERR=1100) DUMPF
   40     FORMAT (80A)
          OPEN (UNIT=1, FILE=DUMPF, STATUS='UNKNOWN', IOSTAT=IOS,
     +          ERR=1200)
C##          CALL DATE(DDATE)
          SCRATCH='Log file from program hydrojet'
C##       SCRATCH='Log file from program hydrojet, '//DDATE//'.'
          WRITE (1, *, IOSTAT=IOS, ERR=1200) SCRATCH
          DEVICE=1
          DUMP=1
      ELSE
          DEVICE=0
          DUMP=0
      END IF


C ................................................................
C REQUESTs USER to SELECT CALCULATION OF MOMENT COEFFICIENTS
C........................................................................

      CALL BEEPER
      YN='N'
      WRITE (6, 11)
   11 FORMAT ('$', 'OPTION: Moment Coefficients [Def=N] (Y/N) ?')
      READ (5, 20, IOSTAT=IOS, ERR=1000) YN
      IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
         IMOMFLAG=1
      ELSE
         IMOMFLAG=0
      END IF


C--------------------------------------------------------------------
C Main Menu
C--------------------------------------------------------------------
C
      CALL BEEPER
  100 WRITE (6, *) ' '
      WRITE (6, *) 'hydrojet: SELECT one of:'
c                  ---------------------------
      WRITE (6, *) '1. Obtain INPUT data:'
      WRITE (6, *) '  (11) Enter/Change INPUT DATA from keyboard.'
      WRITE (6, *) '  (12) Read INPUT DATA & Results from DATA file.'
      WRITE (6, *) '  (13) Read default INPUT DATA from DEFAULT.DAT'

      WRITE (6, *) '2. STORE data:'
      WRITE (6, *) '  (21) Write INPUT & Results to DATA file.'
      WRITE (6, *) '  (22) Make new INPUT default file DEFAULT.DAT'

      WRITE (6, *) '3. CALCULATE SOLN.: Flow,Forces & Dynamic Coeffs:'

c   !.....................................................................!
      IF (BEARING.EQ.1) THEN     ! hydrostatic bearing
c   !.....................................................................!
      WRITE (6, *) '  (31) Ecc=0, GIVEN Orifice Diameter, FIND' ,
     +  ' Recess Pressures'
      WRITE (6, *) '  (32) Ecc=0, GIVEN Recess Pressure, FIND ',
     +  ' Orifice diameters'
      WRITE (6, *) '   ==> (31,32) only for rotationally symmetric HJB'

c   !.....................................................................!
      ELSE                       ! seal or pad journal bearing
c   !.....................................................................!
      WRITE (6, *) '  (31) Ecc=0, Centered Position & NO Misalignment'
c   !.....................................................................!
      END IF
c   !.....................................................................!

      WRITE (6, *) '  (33) Off-centered position ( Ex & Ey <> 0), FIND S
     +olution'
      WRITE (6, *) '  (34) Calculate s-s eccentricity for static load.'
      WRITE (6, *) '  (35) Calculate dynamic coefficients given zeroth
     +order'
      WRITE (6, *) '       solution saved on a DATA File'

      WRITE (6, *) '4. PRINT DATA or RESULTS'
      WRITE (6, *) '  (42) View data/results/save for plot'

      WRITE (6, *) '5. Miscellaneous'
      WRITE (6, *) '  (52) HELP,  (54) Change verbosity and beep.'

      IF (DUMP.eq.1) THEN
          IF (DEVICE.eq.1) THEN
      WRITE (6, *) '  (55) Deactivate DUMP file. (56) Examine DUMP file'
          ELSE
      WRITE (6, *) '  (55) Activate DUMP file, (56) Examine DUMP file'
          END IF
      END IF

      WRITE (6, *) '6. END Execution of hydrojet'
      WRITE (6, *) '  (61) EXIT.--------------------------------------'



C.............................................................
C ENTER CHOICE
C.....................
      WRITE (6, *) ' '
      WRITE (6, 50)
   50 FORMAT ('$', 'SELECT: ')
      READ (5, *, IOSTAT=IOS, ERR=1300) CHOICE
C...............................................................
C
C
C OPTIONS FOR PROGRAM CALCULATIONS
C
C...............................................................
C Option (11): Enter/Change Input params from KBD
C................................................

c     -------------------------------!
      IF (CHOICE.EQ.11) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
          WRITE (1, *) 'Option (11): Enter INPUT data from keyboard.'
          WRITE (1, *) '--------------------------------------------'
          END IF

          CALL KINPUT(RPM(1),PS(1),PA(1))        ! => optionst.f
          CALL CPARAM(RPM(1),PS(1),PA(1))        ! => initopst.f
          CALL XYDATA(0,RPM(1),PS(1),PA(1))     ! => calcmeshp.f

C.................................................................
C Option (12): Read INPUT Parameters from a DATA File
C......................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.12) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
           WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (12): Read INPUT data from a file.'
           WRITE (1, *) '--------------------------------------------'
          END IF
          WRITE (6, 60)
   60     FORMAT ('$', 'Name of data file: ')
          READ (5, 40, IOSTAT=IOS, ERR=1400) FILE

          CALL READATA(FILE, DEVICE,RPM(1),PS(1),PA(1))    ! => readwritet.f
          CALL CPARAM(RPM(1),PS(1),PA(1))        		! => initopst.f
          CALL XYDATA(0,RPM(1),PS(1),PA(1))     		! => calcmeshp.f
          ICHEC2=1

C....................................................................
C Option (13): Reads INPUT from DEFAULT.DAT
C...........................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.13) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
           WRITE (1, *) '---------------------------------------------'
           WRITE (1, *) 'Option (13): Read INPUT data from DEFAULT.DAT'
           WRITE (1, *) '--------------------------------------------'
          END IF

          CALL INPUT	     ! => optionst.f
          CALL CPARAM        ! => initopst.f
          CALL XYDATA(0)     ! => calcmeshp.f




C....................................................................
C Option (21): Stores Data + Results into DATA File
C..................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.21) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
             WRITE (1, *) 'Option (21): Write INPUT data/Results to a D
     +ATA File'
          WRITE (1, *) '--------------------------------------------'

          END IF
          WRITE (6, 60)
          READ (5, 40, IOSTAT=IOS, ERR=1400) FILE

          CALL WRITEDATA(FILE, DEVICE)        !=> readwritet.f

C....................................................................
C Option (22): Creates NEW DEFAULT.DAT
C..........................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.22) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (22):Write INPUT Data -> DEFAULT.DAT'
          WRITE (1, *) '--------------------------------------------'
          END IF

          CALL NEWDEF 			!=> optionst.f


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  CALCULATIONS
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  Option (31):
C: GIVEN AN ORIFICE DIAMETER, FIND RECESS PRESSURE RATIO at
C: JOURNAL CONCENTRIC POSITION
C: THIS options serves a) for design of HJBS b) to provide start guess
c  in bearings and seals
C.......................................................................!
C WARNING:
C For tilt-pad bearings this routine DOES NOT BALANCE THE PAD MOMENTS !!!
C.......................................................................!

c    ------------------------------
      ELSE IF (CHOICE.EQ.31) THEN
c    ------------------------------

          IF ((DEVICE.EQ.1).AND.(BEARING.EQ.1)) THEN
          WRITE (1, *) '--------------------------------------------'
             WRITE (1, *) 'Option (31): Ecc=0, Given Orifice Diameter,
     + Calculate Recess Pressure ratio'
          WRITE (1, *) '--------------------------------------------'
          ELSE IF ((DEVICE.EQ.1).AND.(BEARING.GT.1)) THEN
          WRITE (1, *) '--------------------------------------------'
             WRITE (1, *) 'Option (31): Ecc=0, Concentric Solution'
          WRITE (1, *) '--------------------------------------------'
          END IF

c      ............................................................
          IF (Iwear.eq.1) THEN
 3197     FORMAT ('BEARING WEAR OPTION NOT AVAILABLE')
             WRITE(6,3197)
             IF (Device.eq.1) WRITE(1,3197)
             CALL BEEPER
             CALL PAUSE
             GOTO 100
          END IF
c      ............................................................
c      Works only for discharge pressures equal to a constant
c      not suitable for non-uniform (periodic) exit pressures
c      since bearing rotational symmetry is lost.

          IF ((IPRuni.ne.1).OR.(IPLuni.ne.1)) THEN
             IPRuni=1
             IPLuni=1
             WRITE (6,3198)
             CALL BEEPER
             IF (Device.eq.1) WRITE (1,3198)
          END IF
c      ............................................................

 3198     FORMAT (' ',30('- '),/, 4X,'W A R N I N G',/,4X,
     +        'OPTION Calculates FLOW for UNIFORM EXIT PRESSURES ONLY',
     +  /,4X, ' SET Default EQUAL to Mean Exit Pressure and',
     +        ' CONTINUE',
     +  /,1X,' ',30('- ') )

c      ............................................................

          EXO=0.0D0                         ! For concentric journal position
          EYO=0.0D0                         ! ......................
          AXO=0.0D0                         ! and no journal misalignment
          AYO=0.0D0  			    !

          CALL ZERO0                        ! zeroes forces & coeffs.
          CALL ZERO1                        ! ...............=> addcoefs.f
          IPROP=1                           ! yes to evaluate props

c     ............................................................
c      IGUESP=0  Use current P,U,V&T fields as initial guess for soln.
c      IGUESP=1  Guess P,U,V&T fields from guespp.f to initiate soln.
c     ............................................................


       IF (ITPAD.EQ.0) THEN   ! for Isothermal fluid film(T=Tin=Constant)
        IGUESP=1              !--->>>
       ELSE            !......! other THD models
        WRITE (6, *) 'Choose one of:'
        WRITE (6, *) ' (1) Start with guessed P,U,V & T fields'
        WRITE (6, *) ' (0) Use current P,U,V&T fields as start'
        WRITE (6, *) '     TEMP files need to exist'
	WRITE (*,*) 'OPTION = (?) '             ! =0: No call of GUESP
	READ (*,*)  IGUESP                      ! =1: Call Subroutine GUESP
        IF (IGUESP.LT.0) IGUESP=0
        IF (IGUESP.GT.1) IGUESP=1
       END IF                 !........................................!

C## .................................................................
C## IMPORTANT NOTES:
C## ITPAD is first declared on initopsp.f
C##
C## ITPAD=0 means groove temperatures = Tsupply = Tsump and
C## there is no coupling between any two adjacent pads.
c## Otherwise, heat carry-over energy balance in groove is needed
C## for leading pad edge temperatures to be updated.
C## Subroutine CALCPADT is designed to update Tedges at Ecc=0.
C##
C## ITPAD=0 IS ALWAYS ZERO FOR ISOTHERMAL SOLUTION (ISOTH=1)
C## .................................................................

C       !:::::::::::::::::::::!             !
          IF (ITPAD.EQ.1) THEN              !== FIND SOLUTION INCLUDING
             CALL CALCPADT(RPM(1),PS(1),PA(1))		    !   THERMAL EFFECTS ON
          END IF                            !== GROOVES
C       !:::::::::::::::::::::!             !   ==> calcpadt.f


C       --------------------------------------------------------------
C        FOR ITPAD=0:  LH2 or if ISOTH=1 (isothermal) conditions; THEN:
C       --------------------------------------------------------------

C       !................!..................!
        DO KPAD=1, NPAD			    ! FOR EACH PAD on BEARING
C       !................!..................!
          IF (IFULL.EQ.0) THEN              ! => Pad Bearing
           WRITE (6,7979) KPAD
           IF (DEVICE.EQ.1) WRITE (1,7979) KPAD
          END IF                            !...................!

          CALL XYDATA(KPAD)                 ! generates mesh & film H
          CALL FILMH(EXO, EYO, NXT, NYI)    !
          CALL FILMC                        !
          CALL SETBC(RPM,PS,PA)                        ! Set known Bound. Conds => calcsont.f
c        !--------------------!             !......................!
          IF (ITPAD.EQ.1) THEN              ! Get the converged
c        !--------------------!             !
             CALL READTEMP(KPAD)            ! THD solutions
c        !--------------------!             !......................!
          ELSE                              !
c        !--------------------!             !......................!
             IF (KPAD.EQ.1.AND.IGUESP.GT.0) THEN
                CALL GUESP(DEVICE,RPM,PS,PA)         ! Guesp a Start Pfield
                CALL FILMC                  ! and modify film thickness
             END IF                         !
c        !--------------------!             !......................!
          END IF                            !
c        !--------------------!             !......................!

          CALL CALCPR(DEVICE)               ! SOLVE FLOW PROBLEM

          CALL WRITETEMP(KPAD)              ! & store in temporal file
          CALL ADD0                         !
          CALL ADD1                         ! & adds forces & coeffs.

C       !................!..................!
        END DO                              ! KPAD=1,... NPAD
C       !................!..................!

       ICHEC2=1
c     !---------------------------------! FOR all PADS
       IF ((NPAD.GT.1).OR.(TILTPAD.EQ.1)) THEN
          CALL PRINT0PAD(DEVICE)        !==> addcoefs.f
          CALL PRINT1PAD(DEVICE,RPM,PS,PA)        !==> ''''
       END IF				!
c     !---------------------------------!
       CALL PAUSE			!==> supportp.f

C 	XYDATA               ==> calcmeshp.f
C       FILMH                ==> filmf.f
C       FILMC                ==> foilsubs.f
C       CALCPR               ==> optionst.f
C 	READTEMP, WRITETEMMP ==> readwritep.f
C 	ADD0, ADD1           ==> addcoefs.f
C
C..............................................................
C  OPTION (32):
C: GIVEN A RECESS PRESSURE RATIO FIND ORIFICE DIAMETER
C...................................................................
C Option 32 is used for design purposes only in hydrostatic bearings
c Given the recess pressure, hydrojet calculates the orifice
c diameters for each bearing recess.
C RECESS PRESSURE MUST BE THE SAME FOR ALL RECESSES
c...................................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.32) THEN
c     --------------------------------!
          IF (BEARING.GT.1) GOTO 100  !=>seals and pad bearings

          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (32): Ecc=0, Given Recess Pressure, C
     +alculate Orifice Diameter'
          WRITE (1, *) '--------------------------------------------'
          END IF

c      ............................................................
          IF (Iwear.eq.1) THEN       ! NO WEAR OPTION Available
             WRITE(6,3197)
             IF (Device.eq.1) THEN
               WRITE(1,3197)
             END IF
             CALL BEEPER
             CALL PAUSE
             GOTO 100
          END IF
c      ............................................................
c      Works only for discharge pressures equal to a constant
c      not suitable for non-uniform (periodic) exit pressures
c      since rotational symmetry is lost.

          IF ((IPRuni.ne.1).OR.(IPLuni.ne.1)) THEN
             IPRuni=1
             IPLuni=1
             WRITE (6,3198)
             CALL BEEPER
             IF (Device.eq.1) THEN
             WRITE (1,3198)
             END IF
          END IF
c      ............................................................

          EXO=0.0                         !for Concentric position
          EYO=0.0                         !.......................
          AXO=0.0D0			  ! w/o misalignment
          AYO=0.0D0  			  !
          CALL ZERO0		          ! zeroes forces and
          CALL ZERO1			  ! force coefficients => addcoefs.f

          IPROP=1                           ! yes to evaluate props
          IGUESP=1                          ! start with guess fields

C       !................!..................!
          DO KPAD=1, NPAD		    ! On all bearing pads
C       !................!..................!
          IF (IFULL.EQ.0) THEN              ! => Pad Bearing
           WRITE (6,7979) KPAD
           IF (DEVICE.EQ.1) WRITE (1,7979) KPAD
          END IF                            !.......................!

            CALL XYDATA(KPAD)		    ! generate pad mesh => calcmeshp.f
            CALL FILMH(EXO, EYO, NXT,NYI)   ! & film thickness
            CALL FILMC                      ! w/ bearing compliance
            CALL SETBC(RPM,PS,PA)                      ! Set known bound. conds.

          IF (KPAD.EQ.1) THEN               !.......................!
            CALL GUESP(DEVICE,RPM,PS,PA)              ! Guesp a Start Pfield
            CALL FILMC                      !
          END IF                            !.......................!

            CALL DESIGN(DEVICE)             ! SOLVE flow field => optionst.f

          IF (NPOCKET.GT.0) THEN            ! HJBs
            DO J=1, NPOCKET                 !....................
               PaDorif(J,KPAD)=Diaorif(J)   ! save orifice diams.
            END DO		            !....................
          END IF                            !

            CALL WRITETEMP(KPAD)            ! & store in temporal file => readwritep.f
            CALL ADD0                       !
            CALL ADD1                       ! & sums forces & coeffs.

C       !................!..................!
          END DO                            ! K=1,2,..., NPAD
C       !................!..................!

C 	XYDATA               ==> calcmesht.f
C       FILMH                ==> filmwt.f
C       FILMC                ==> foilsubs.f
C       DESIGN               ==> optionst.f
C 	READTEMP, WRITETEMMP ==> readwritet.f
C 	ADD0, ADD1           ==> addcoefs.f

       ICHEC2=1
c     !---------------------------------! FOR all PADS
       IF ((NPAD.GT.1).OR.(TILTPAD.EQ.1)) THEN
          CALL PRINT0PAD(DEVICE)        !==> addcoefs.f
          CALL PRINT1PAD(DEVICE)        !==> ''''
       END IF				!
c     !---------------------------------!
       CALL PAUSE			!==> supportp.f

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  Option (33):
C: FOR OFFCENTER POSITION, CALCULATE FIELDS AND FORCES
C.....................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.33) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '---------------------------------------'
              WRITE (1, *)
     +   'Option (33): Ecc>0., Solution for fixed Eccentricity'
          WRITE (1, *) '---------------------------------------'
          END IF
C                          ! IMPORTANT:    ! Input Data + guess have
C                                          ! already been read by (12)


          IPROP=1                          ! yes to evaluate props
          IGUESP=0                         ! TEMP files already existent

c        !:::::::::::::::::::::::::!
          IF (TILTPAD.EQ.1) THEN
c        !:::::::::::::::::::::::::!
          LOSTAT=0
          CALL TILTECC(DEVICE,LOSTAT)   ! => on tiltsubs.f
          GOTO  1530

c        !:::::::::::::::::::::::::!
          END IF
c        !:::::::::::::::::::::::::!

c        !.................................!
          ISYMold=ISYM                     ! Save current value of SYMMETRY
c                                          !................................!
c        !.................................! CHANGE Journal Displacements
          CALL ASKITER                     ! Exo, Eyo, Axo, Ayo, Zo
c        !.................................! ==> initopst.f

c        !.................................! REQUEST User IFF Force Coeffics.
          YN='Y'			   ! will be calculated
          WRITE (6, 1520)
 1520     FORMAT ('$', 'Do you wish to calculate dynamic ',
     +            'coefficients [Def: Y] (Y/N): ')
     	  READ (5, 1521) YN
 1521     FORMAT (1A)
c        !.................................! REQUEST Ends

 1525     CALL ZERO0                        !=> addcoefs.f
          CALL ZERO1			    ! zeroes forces and coefficients

C## .................................................................
C## IMPORTANT NOTES:
C## ITPAD is first declared on initopsp.f
C##
C## ITPAD=0 means groove temperatures = Tsupply = Tsump and
C## there is no coupling between any two adjacent pads.
c## Otherwise, heat carry-over energy balance in groove is needed
C## for leading pad edge temperatures to be updated.
C## Subroutine CALCPADT is designed to update Tedges at Ecc=<>0
C##
C## ITPAD=0 IS ALWAYS ZERO FOR ISOTHERMAL SOLUTION (ISOTH=1)
C## .................................................................


C       !:::::::::::::::::::::!             !
          IF (ITPAD.EQ.1) THEN              !== FIND SOLUTION INCLUDING
             CALL CALCPADT(RPM(1),PS(1),PA(1))		    !   THERMAL EFFECTS ON
          END IF                            !== GROOVES
C       !:::::::::::::::::::::!             !   ==> calcpadt.f


C      --------------------------------------------------------------
C      FOR ITPAD=0:  LH2 or if ISOTH=1 (isothermal) conditions; THEN:
C      --------------------------------------------------------------

C   !--------------------!
 1510  DO KPAD=1, NPAD   ! SWEEP OVER ALL PADS:
C   !--------------------!
          IF (IFULL.EQ.0) THEN                  ! => Pad Bearing
           WRITE (6,7979) KPAD                  !
           IF (DEVICE.EQ.1) WRITE (1,7979) KPAD !
          END IF                                !.................!

          CALL XYDATA(KPAD)                     !==> calcmeshp.f
          CALL READTEMP(KPAD)                   !==> readwritep.f
c                                               !    READ initial fields

          IF ((ISYMold.eq.1).AND.(ISYM.eq.0)) THEN   ! Set UVP mirror image to
	     CALL TEMPASYM                           ! start solution for
          END IF                                     ! asymmetric HJB

          IF (IWEAR.EQ.1) THEN                       ! Calculate film
             CALL HWEAR(EXO, EYO, NXT, NYI)          ! thickness due to
          ELSE                                       ! journal displacements
             CALL FILMH(EXO, EYO, NXT, NYI)          !==> filmf.f
          END IF                                     !......................!

          DUMYI=0


C IF ITPAD=1: Solution was already found in Sub. CALCPADT
C             THIS is for the thermal solution case.

c       !------------------------!             !.....................!
          IF (ITPAD.EQ.0) THEN                 ! ITPAD=0
c       !------------------------!             !.....................!

          DUMYR=MAXIN(P)                       ! Determines if P field
          IF (DUMYR.LT.(1.0D-3)) THEN          ! is zero.
              CALL SETBC(RPM,PS,PA)                       ! start with a Pguess,
              CALL GUESP(DEVICE,RPM,PS,PA)               ! otherwise Use
          END IF                               ! current soln.
c        !...................!                 !.....................!
          CALL SOLVE(DEVICE,RPM,PS,PA)                   ! SOLVE for flow field
c        !...................!                 !.....................!

c       !------------------------!             !.....................!
          ELSE
c       !------------------------!             !.....................!
          CALL FILMC                           ! upadte film thickness
c       !------------------------!             !.....................!
          END IF 			       ! ITPAD=1
c       !------------------------!             !.....................!

          CALL PRINTDUVP(DEVICE, 0)            ! get forces & torque
          CALL FORCE(DEVICE)                   !
          CALL TORQUE(DEVICE)                  !........................
          CALL PRINTF(DEVICE,1)                ! prints results for pad
          CALL ADD0                            ! Adds forces & moments

C 	SETBC, SOLVE, FORCE, TORQUE ==> calcsolnt.f
C 	PRINTF, PRINTDUVP    ==> prntrest.f
C       GUESP                ==> guespp.f
C       FILMC                ==> foilsubs.f

C       !......................................................!
c        CHECK whether Coefficients are to be calculated or not

c       !......................................!
     	  IF ((YN.EQ.' ').OR.(YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
c       !......................................!
              DUMYI=3
              ICASE=2
              IF (RPM(1).EQ.(0.0D0)) THEN
                    FREQU=100.00D0  ! HZ
              ELSE
                    FREQU=RPM(1)/60.0D0
              END IF

c          !...................................!
              CALL INPUT1(DUMYI, DEVICE,RPM(1),PS(1),PA(1))       ! solves for coefficients
              CALL PRINTCOEF(DEVICE)           ! prints them OUT
              CALL WRITETEMP(KPAD)             ! stores in TEMP file
              CALL ADD1                        ! and adds them to others
              ICHEC2=1                         !
c       !......................................!
      	  ELSE
c       !......................................! No calculation
              CALL WRITETEMP(KPAD)	       ! of force coefficients
              ICHEC2=0                         !
c       !......................................!
     	  END IF
c       !......................................! End of Check

C	 INPUT1               ==> optionst.f => firstt.f (coeffic.f)
C	 PRINTCOEF            ==> prntrest.f
C	 READTEMP, WRITETEMMP ==> readwritet.f

C   !--------------------!
        END DO           ! KPAD=1, NPAD
C   !--------------------!

C PRINTS TOTAL FORCES AND COEFFICIENTS FOR NPADS:
C -----
c     !---------------------------------! FOR all PADS
       IF ((NPAD.GT.1).OR.(TILTPAD.EQ.1)) THEN
          CALL PRINT1PAD(DEVICE,RPM,PS,PA)       !==> addcoefs.f
       IF (ICHEC2.EQ.1) CALL PRINT1PAD(DEVICE)
       END IF				!
c     !---------------------------------!
 1530  CALL PAUSE			!==> supportp.f
C
C................................................................
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Option (34)
C Calculate Exo,Eyo for given s-s load
C.......................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.34) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (34): Calculate s-s eccentricity for
     + a given static LOAD'
          WRITE (1, *) '--------------------------------------------'
          END IF

          IPROP=1                  ! yes to evaluate props
c        !:::::::::::::::::::::::::!
          IF (TILTPAD.EQ.1) THEN   !
c        !:::::::::::::::::::::::::!
          CALL LOADTILT(DEVICE,RPM(1),PS(1),PA(1))    ! => on tiltsubst.f
c        !:::::::::::::::::::::::::!
          ELSE
c        !:::::::::::::::::::::::::!
          CALL LOAD(DEVICE)        ! => on optionst.f
c        !:::::::::::::::::::::::::!
          END IF
c        !:::::::::::::::::::::::::!

          ICHEC2=1
          CALL PAUSE		   !=> on supportp.f


C.....................................................................
C Option (35)
C Calculate first order solution given a zeroth order soln.
C............................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.35) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (35): Calculate dynamic coefficients
     +given the s-s soln.'
          WRITE (1, *) '--------------------------------------------'
          END IF

          WRITE (6, 60)                     ! Read File Name
          READ (5, 40, IOSTAT=IOS, ERR=1400) FILE
      OPEN (UNIT=2, FILE=FILE, STATUS='OLD', IOSTAT=IOS, ERR=1777)
          CLOSE (UNIT=2)

          IPROP=1                            ! yes to evaluate props
          CALL READATA(FILE, DEVICE,RPM(1),PS(1),PA(1))         ! => readwritet.f
          CALL CPARAM(RPM(1),PS(1),PA(1))			     ! => initopst.f
          CALL XYDATA(0,RPM(1),PS(1),PA(1))		     ! => calcmeshp.f
          CALL ZERO0        		     ! => addcoefs.f
          CALL ZERO1                         !    '''''

          FREQU=RPM(1)/60.0D0
          IF (RPM(1).EQ.0.0D0 )FREQU=100.0D0

      WRITE (6, 6099) FREQU
 6099 FORMAT ('$', 'ENTER Excitation Frequency(w) in Hz [Default (w)=',
     + E12.5E2,']: ')
      READ (5, 6199) DUMYC
 6199 FORMAT (40A)
      IF (DUMYC.NE.' ') THEN
          OPEN (UNIT=45, STATUS='SCRATCH')
          WRITE (45, 6299) DUMYC
 6299     FORMAT (' ', 40A)
          BACKSPACE (UNIT=45)
          READ (45, *) FREQU
          CLOSE (UNIT=45)
      END IF

      DUMYI=1                           ! Ifile=1

C     !.................................!
       DO KPAD=1, NPAD
C     !.................................!
          IF (IFULL.EQ.0) THEN          ! => Pad Bearing
           WRITE (6,7979) KPAD
           IF (DEVICE.EQ.1) WRITE (1,7979) KPAD
          END IF                        !.....................!

          CALL XYDATA(KPAD)             ! forms pad mesh
          CALL READTEMP(KPAD)           ! reads zeroth order solution
          ICASE=2                       ! X&Y perturbation.
          CALL INPUT1(DUMYI, DEVICE,RPM(1),PS(1),PA(1))    ! first order soln.
          CALL PRINTCOEF(DEVICE)        ! calc. pad coefficients
          CALL ADD0                     ! add coefficients
          CALL ADD1                     !
          CALL WRITETEMP(KPAD)          ! transfer to TEMP files
C     !.................................!
       END DO
C     !.................................! KPAD=1, NPAD

C	 XYDATA               ==> calcmeshp.f
C	 READTEMP, WRITETEMMP ==> readwritet.f
C	 INPUT1               ==> optionst.f => firstt.f (coeffic.f)
C	 PRINTCOEF            ==> prntrest.f
C	 ADD0, ADD1           ==> addcoefst.f

       ICHEC2=1
c     !---------------------------------! FOR all PADS
       IF ((NPAD.GT.1).OR.(TILTPAD.EQ.1)) THEN
          CALL PRINT0PAD(DEVICE)        !==> addcoefs.f
          CALL PRINT1PAD(DEVICE,RPM,PS,PA)        !==> ''''
       END IF				!
c     !---------------------------------!
       CALL PAUSE			!==> supportp.f

C-------------------------------------------------------------------C

C................................................................
C Option (42): print results on monitor
C................................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.42) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Choice (42): View data and/or results'
          WRITE (1, *) '--------------------------------------------'
          END IF

          CALL DISPLAY(DEVICE,ICHEC2,RPM(1),PS(1),PA(1))		!=> prntrest.f

C................................................................
C Option (52): HELP
C................................................................
C
c     -------------------------------!
      ELSE IF (CHOICE.EQ.52) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Choice (52): HELP'
          WRITE (1, *) '--------------------------------------------'
          END IF
          CALL HELP			        !=> supportp.f


C................................................................
C Option (54): Verbosity
C................................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.54) THEN
c     -------------------------------!
          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Option (54): Modify verbosity and beep.'
          WRITE (1, *) '--------------------------------------------'
          END IF
          WRITE (6, *) ' '
          WRITE (6, *) 'Choices are:'
          WRITE (6, *) ' (0) Low: Gives only essential information.'
          WRITE (6, *) ' (1) Medium: Gives short progress reports during
     + execution.'
          WRITE (6, *) ' (2) High: Gives complete information during exe
     +cution.'
          WRITE (6, *) ' '
 1610     IF (SVERB.gt.1) THEN
              IF (SVERB.EQ.1) THEN
                  WRITE (6, 301)
  301             FORMAT ('$Screen verbosity (Default: Medium): ')
              ELSE
                  WRITE (6, 302)
  302             FORMAT ('$Screen verbosity (Default: High): ')
              END IF
          ELSE
              WRITE (6, 303)
  303         FORMAT ('$Screen verbosity (Default: Low): ')
          END IF
          READ (5, *, IOSTAT=IOS, ERR=1600) SVERB
          IF ((SVERB.LT.0).OR.(SVERB.GT.2)) THEN
              SVERB=1
          END IF
          IF (SVERB.EQ.2) THEN
              SVERB=3
          END IF
 1910     IF (DVERB.gt.1) THEN
              IF (DVERB.EQ.1) THEN
                  WRITE (6, 311)
  311             FORMAT ('$Dump verbosity (Default: Medium): ')
              ELSE
                  WRITE (6, 312)
  312             FORMAT ('$Dump verbosity (Default: High): ')
              END IF
          ELSE
              WRITE (6, 313)
  313         FORMAT ('$Dump verbosity (Default: Low): ')
          END IF
          READ (5, *, IOSTAT=IOS, ERR=1900) DVERB
          IF ((DVERB.LT.0).OR.(DVERB.GT.2)) THEN
              DVERB=0
          END IF
          IF (DVERB.EQ.2) THEN
              DVERB=3
          END IF
 1710     IF (BEEP.eq.1) THEN
              WRITE (6, 321)
  321         FORMAT ('$Beep (Default: Yes): ')
          ELSE
              WRITE (6, 322)
  322         FORMAT ('$Beep (Default: No): ')
          END IF
          READ (5, 20, IOSTAT=IOS, ERR=1700) YN
          IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
              BEEP=1
          ELSE
              BEEP=0
          END IF

C................................................................
C Option (55): Work with DUMP File
C...................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.55) THEN
c     -------------------------------!
          IF (DUMP.eq.1) THEN
              IF (DEVICE.eq.1) THEN
                  DEVICE=0
                  WRITE (1, *) 'Option (55): Deactivate DUMP file.'
                  WRITE (6, *) 'DUMP file deactivated.'
              ELSE
                  DEVICE=1
                  WRITE (1, *) 'Option (55): Activate DUMP file.'
                  WRITE (6, *) 'DUMP file activated.'
              END IF
          ELSE
              WRITE (6, *) 'No DUMP file was requested.'
          END IF



       ELSE IF (CHOICE.EQ.56) THEN

          IF (DUMP.eq.1) THEN
              I=0
              IF (DEVICE.eq.1) THEN
                  WRITE (1, *) 'Option (56): Examine DUMP file.'
              END IF
              WRITE (1, *) 'End of dump file.'
              REWIND (UNIT=1)
  210         READ (1, 220, IOSTAT=IOS, ERR=1800) SCRATCH
  220         FORMAT (80A)
              I=I+1
              IF (I.EQ.25) THEN
                  I=1
                  CALL PAUSE
              END IF
              WRITE (6, 222) SCRATCH		 ! output to screen
  222			FORMAT(80A)
              IF (SCRATCH.EQ.' End of dump file.') THEN
                  GOTO 200
              END IF
              GOTO 210
  200         BACKSPACE (UNIT=1)
              CALL PAUSE
          ELSE
              WRITE (6, *) 'No dump file was requested.'
          END IF

C................................................................
C Option (61): EXIT, end of execution
C................................................................

c     -------------------------------!
      ELSE IF (CHOICE.EQ.61) THEN
c     -------------------------------!

          IF (DEVICE.eq.1) THEN
          WRITE (1, *) '--------------------------------------------'
              WRITE (1, *) 'Choice (61): EXIT hydrojet'
          WRITE (1, *) '--------------------------------------------'
          END IF
          WRITE (6, *)
     +    'End of <hydrojet> for flex-pad and compliant bearings'
          CALL BEEPER
          IF (DEVICE.eq.1) THEN
              CLOSE (UNIT=1)
          END IF
          STOP

c.................................................................

c     -------------------------------!
      ELSE
c     -------------------------------!

        WRITE (6, *) 'Enter the number of your choice, or 52 for help'


c     -------------------------------!
      END IF
c     -------------------------------! END OF OPTIONS

c.................................................................
c.................................................................
c...............................................................

      GOTO 100
c.................................................................
c.................................................................
c.................................................................
C
 7979  FORMAT (' +',77('-'),'+',
     +       /,' |',3X,'RESULTS FOR BEARING PAD #',
     1        I2,46X,' |',/,' |',77('.'),'|')
C
C
C................................................................
C
C Error handling
C
C................................................................


C1000 - ERROR IN Y/N DUMP FILE OPTION
 1000 IF (IOS.EQ.64) THEN
          WRITE (6, *) 'Please enter a Y or N.'
          GOTO 1010
      ELSE
          CALL DECODIOS(IOS)
          WRITE (6, *) 'Exiting HYDROJET'
          STOP
      END IF

C1100 - ERROR IN RETRIVING DUMP FILE NAME
 1100 CALL DECODIOS(IOS)
      WRITE (6, *) 'Please reenter name.'
      GOTO 1110

C1200 - ERROR OPENING OR MAINTAINING DUMP FILE
 1200 CALL DECODIOS(IOS)
      WRITE (6, *) 'Returning to Y/N choice.'
      GOTO 1010

C1300 - ERROR GETTING CHOICE
 1300 IF (IOS.EQ.64) THEN
          WRITE (6, *) 'Enter the number of your choice, or 41 for help.
     +'
      ELSE
          CALL DECODIOS(IOS)
          CALL PAUSE
      END IF
      GOTO 100

C1400 - ERROR GETTING FILE NAME
 1400 WRITE (6, *) 'Error in file name input:'
      CALL DECODIOS(IOS)
      CALL PAUSE
      GOTO 100

C1777 - ERROR OPENING FILE NAME
 1777  WRITE (6, *) 'Error in OPENING file :'
      CALL DECODIOS(IOS)
      CALL PAUSE
      GOTO 100


C1600 - ERROR IN SCREEN VERBOSITY INPUT READ
 1600 IF (IOS.NE.64) THEN
          CALL DECODIOS(IOS)
      END IF
      WRITE (6, *) 'Please enter a number, 0, 1, or 2.'
      GOTO 1610

C1700 - ERROR IN BEEP INPUT READ
 1700 CALL DECODIOS(IOS)
      WRITE (6, *) 'Please enter a Y or N.'
      GOTO 1710

C1800 - ERROR IN DUMP FILE READ
 1800 CALL DECODIOS(IOS)
      WRITE (6, *) 'Returning to main menu.'
      CALL PAUSE
      GOTO 100

C1900 - ERROR IN DUMP VERBOSITY INPUT
 1900 IF (IOS.NE.64) THEN
          CALL DECODIOS(IOS)
      END IF
      WRITE (6, *) 'Please enter a number, 0, 1, or 2.'
      GOTO 1910

      END


c==============================================================================
c --------------------------------------------------------------------------- C

      SUBROUTINE DISCLOSURE

c --------------------------------------------------------------------------- C
c==============================================================================

      write (6,173)

  173 format (5X,
     +'    The Texas A&M University System (TAMU) warrants',       /,5X,
     +'that all software conforms to applicable TAMU published',   /,5X,
     +'specifications prevailing at the time of shipment. Except', /,5X,
     +'for this warranty,TAMU disclaims all warranties in regard', /,5X,
     +'to product/services including all implied warranties of',   /,5X,
     +'merchantibility and fitness for a particular purpose; and', /,5X,
     +'the stated express warranty is in lieu of all obligations', /,5X,
     +'or liabilities on the part of TAMU for damages including,',/,5X,
     +'but not limited to special, indirect, or consequential ',/,5X,
     +'damages arising out of/or in connection with the use or',/,5X,
     +'performance of the products or services provided.',      /,5X,
     +'(c)Copyright 1995 Texas A&M University System/',            /,5X,
     +'   Dr. Luis San Andres,    ALL RIGHTS RESERVED.',/,/,/)

      !CALL PAUSE

      END

c
c==============================================================================
c
c main program for hydrojet : by Drs. Luis San Andres 1995
c based on hydroflext.f & hydrosealt.f programs
c
C==============================================================================
C