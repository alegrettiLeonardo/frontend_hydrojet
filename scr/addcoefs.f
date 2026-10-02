c   ##    #####   #####    ####    ####   ######  ######   ####           ######
c  #  #   #    #  #    #  #    #  #    #  #       #       #               #
c #    #  #    #  #    #  #       #    #  #####   #####    ####           #####
c ######  #    #  #    #  #       #    #  #       #            #   ###    #
c #    #  #    #  #    #  #    #  #    #  #       #       #    #   ###    #
c #    #  #####   #####    ####    ####   ######  #        ####    ###    #

c #    #   #   #  #####   #####    ####   ######  #       ######  #    #
c #    #    # #   #    #  #    #  #    #  #       #       #        #  #
c ######     #    #    #  #    #  #    #  #####   #       #####     ##
c #    #     #    #    #  #####   #    #  #       #       #         ##
c #    #     #    #    #  #   #   #    #  #       #       #        #  #
c #    #     #    #####   #    #   ####   #       ######  ######  #    #

C addcoefs.f > hydroflex code Drs. Luis San Andres, TexasA&MUniv. 1994
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR I
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C *****************************************************************************
C **                                                                         **
C **  Subroutine ZERO0: zeroes TOTAL PADS forces and flow rates              **
C **                                                                         **
C *****************************************************************************
      SUBROUTINE ZERO0

      IMPLICIT NONE
C .....................................................................
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      DOUBLE PRECISION  PMIN, PMAX, FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
C .....................................................................
C      DATA FXT,FYT,TOT,MXT,MYT,QINT,QOUTT/7*0.0D0/
C F77 > DATA CAN NOT be used for more than one initialization !!!!
C ------------------------------------------------------------------------

      FXT=0.0D0
      FYT=0.0D0
      TOT=0.0D0
      MXT=0.0D0
      MYT=0.0D0
      QINT=0.0D0
      QOUTT=0.D0
      PMIN=0.0D0
      PMAX=0.0D0

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine ZERO1: zeroes TOTAL PADS dynamic force coefficients         **
C **                                                                         **
C *****************************************************************************
      SUBROUTINE ZERO1

      IMPLICIT NONE
C .....................................................................
      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /DAMPT/ CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT
      COMMON /INERT/ MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT
      COMMON /STIAT/ KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT
      COMMON /DAMAT/ CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT
      COMMON /INEAT/ MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
      DOUBLE PRECISION
     +         KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT,
     +         CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT,
     +         MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT,
     +         KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT,
     +         CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT,
     +         MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
C .....................................................................
c      DATA     KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT,
c     +         CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT,
c     +         MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT,
c     +         KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT,
c     +         CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT,
c     +         MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
c     +         /48*0.0D0/

      KXXDT=0.0D0
      KYYDT=0.0D0
      KXYDT=0.0D0
      KYXDT=0.0D0
      CXXDT=0.0D0
      CYYDT=0.0D0
      CXYDT=0.0D0
      CYXDT=0.0D0
      MXXDT=0.0D0
      MYYDT=0.0D0
      MXYDT=0.0D0
      MYXDT=0.0D0
      KXXAT=0.0D0
      KYYAT=0.0D0
      KXYAT=0.0D0
      KYXAT=0.0D0
      CXXAT=0.0D0
      CYYAT=0.0D0
      CXYAT=0.0D0
      CYXAT=0.0D0
      MXXAT=0.0D0
      MYYAT=0.0D0
      MXYAT=0.0D0
      MYXAT=0.0D0
      KmXXDT=0.0D0
      KmYYDT=0.0D0
      KmXYDT=0.0D0
      KmYXDT=0.0D0
      CmXXDT=0.0D0
      CmYYDT=0.0D0
      CmXYDT=0.0D0
      CmYXDT=0.0D0
      MmXXDT=0.0D0
      MmYYDT=0.0D0
      MmXYDT=0.0D0
      MmYXDT=0.0D0
      KmXXAT=0.0D0
      KmYYAT=0.0D0
      KmXYAT=0.0D0
      KmYXAT=0.0D0
      CmXXAT=0.0D0
      CmYYAT=0.0D0
      CmXYAT=0.0D0
      CmYXAT=0.0D0
      MmXXAT=0.0D0
      MmYYAT=0.0D0
      MmXYAT=0.0D0
      MmYXAT=0.0D0
C     .................!

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine ADD0: SUMS PAD FORCES & MOMENTS & FLOWS                     **
C **                                                                         **
C *****************************************************************************
      SUBROUTINE ADD0

      IMPLICIT NONE

      INCLUDE 'params.f'
C .....................................................................
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK),
     +               QREC(MAXNPOCK),QIN,QOUT,QFACTOR
      COMMON/FORCE0/ FFACTOR,FX,FY,TO,TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL

      DOUBLE PRECISION PMIN,PMAX,FFACTOR,FX,FY,TO,TOR,MFACTOR,MX,MY,
     +                 PREC,QREC,TREC,QIN,QOUT,QFACTOR,
     +                 FXT,FYT,TOT,MXT,MYT,QINT,QOUTT

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
C ------------------------------------------------------------------------
      QINT=QINT+QIN*QFACTOR
      QOUTT=QOUTT+QOUT*QFACTOR
      FXT=FXT+FX
      FYT=FYT+FY
      MXT=MXT+MX
      MYT=MYT+MY
      TOT=TOT+TO*TOR

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine ADD1: SUMS PAD FORCES & MOMENTS COEFFICIENTS                **
C **                                                                         **
C *****************************************************************************
      SUBROUTINE ADD1

      IMPLICIT NONE

      INCLUDE 'params.f'

C .....................................................................
      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /DAMPT/ CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT
      COMMON /INERT/ MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT
      COMMON /STIAT/ KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT
      COMMON /DAMAT/ CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT
      COMMON /INEAT/ MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT

      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT

      DOUBLE PRECISION
     +         KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT,
     +         CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT,
     +         MXXDT,MYYDT,MXYDT,MYXDT,MmXXDT,MmYYDT,MmXYDT,MmYXDT,
     +         KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT,
     +         CXXAT,CYYAT,CXYAT,CYXAT,CmXXAT,CmYYAT,CmXYAT,CmYXAT,
     +         MXXAT,MYYAT,MXYAT,MYXAT,MmXXAT,MmYYAT,MmXYAT,MmYXAT
      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION RSINPK, RCOSPK, IPAD
      INTEGER TILT
C ---------------------------------------------------------------------
C FOR TILT PAD bearings reduce the force coefficients with frequency w
c
      IF (TILT.EQ.1) THEN
          CALL REDUCECOEFS
          RETURN
      END IF

C .....................................................................
C ------------------------------------------------------------------------

c    !-----------------------!
C     for fixed-pad bearings:
c    !-----------------------!
C     .................! FORCE COEFS. DUE TO DISPACEMENTS
      KXXDT=KXXDT+KXXD
      KYYDT=KYYDT+KYYD
      KXYDT=KXYDT+KXYD
      KYXDT=KYXDT+KYXD
      CXXDT=CXXDT+CXXD
      CYYDT=CYYDT+CYYD
      CXYDT=CXYDT+CXYD
      CYXDT=CYXDT+CYXD
      MXXDT=MXXDT+MXXD
      MYYDT=MYYDT+MYYD
      MXYDT=MXYDT+MXYD
      MYXDT=MYXDT+MYXD
C     .................! FORCE COEFS. DUE TO ROTATIONS
      KXXAT=KXXAT+KXXA
      KYYAT=KYYAT+KYYA
      KXYAT=KXYAT+KXYA
      KYXAT=KYXAT+KYXA
      CXXAT=CXXAT+CXXA
      CYYAT=CYYAT+CYYA
      CXYAT=CXYAT+CXYA
      CYXAT=CYXAT+CYXA
      MXXAT=MXXAT+MXXA
      MYYAT=MYYAT+MYYA
      MXYAT=MXYAT+MXYA
      MYXAT=MYXAT+MYXA
C     .................! MOMENT COEFS. DUE TO DISPACEMENTS
      KmXXDT=KmXXDT+KmXXD
      KmYYDT=KmYYDT+KmYYD
      KmXYDT=KmXYDT+KmXYD
      KmYXDT=KmYXDT+KmYXD
      CmXXDT=CmXXDT+CmXXD
      CmYYDT=CmYYDT+CmYYD
      CmXYDT=CmXYDT+CmXYD
      CmYXDT=CmYXDT+CmYXD
      MmXXDT=MmXXDT+MmXXD
      MmYYDT=MmYYDT+MmYYD
      MmXYDT=MmXYDT+MmXYD
      MmYXDT=MmYXDT+MmYXD
C     .................! MOMENT COEFS. DUE TO ROTATIONS
      KmXXAT=KmXXAT+KmXXA
      KmYYAT=KmYYAT+KmYYA
      KmXYAT=KmXYAT+KmXYA
      KmYXAT=KmYXAT+KmYXA
      CmXXAT=CmXXAT+CmXXA
      CmYYAT=CmYYAT+CmYYA
      CmXYAT=CmXYAT+CmXYA
      CmYXAT=CmYXAT+CmYXA
      MmXXAT=MmXXAT+MmXXA
      MmYYAT=MmYYAT+MmYYA
      MmXYAT=MmXYAT+MmXYA
      MmYXAT=MmYXAT+MmYXA
C     .................!
      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine PRINT0PAD: PRINTS TOTAL FORCES , FLOW & MOMENTS             **
C **                                                                         **
C *****************************************************************************
C
      SUBROUTINE PRINT0PAD(DEVICE)
C
C----------------------------------------------------------------------
C
C PRINT0PAD: Prints total flow, forces, moments & torque
C
C .....................................................................
      IMPLICIT NONE

      INCLUDE 'params.f'

      COMMON /RESULTS0/ FX,FY,TO,MX,MY,QIN,QOUT
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION LOAD,ANGLE, FX,FY,TO,MX,MY,QIN,QOUT,
     +                 Qkgsec, Qkgmin, ZERO

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        DEVICE, ISYM, ICSTEP, Iwrite, IU, BEARING

C ------------------------------------------------------------------------
      ZERO=0.0D0

      LOAD=DSQRT(FX*FX+FY*FY)

C .............................................!
C FIND ANGLE BETWEEN FLUID FORCE F AND X-AXIS
C .............................................!

      IF (FX.NE.ZERO) THEN

          ANGLE=DATAN(DABS(FY/FX))*180.0D0/DACOS(-1.0D0)

          IF (FX.LT.ZERO) THEN
              IF (FY.LT.ZERO) THEN
                  ANGLE=ANGLE
              ELSE
                  ANGLE=-ANGLE
              END IF
          ELSE
              IF (FY.LT.ZERO) THEN
                  ANGLE=-ANGLE+180.0D0
              ELSE
                  ANGLE=ANGLE+180.0D0
              END IF
          END IF


      ELSE

          IF (FY.GE.ZERO) THEN
              ANGLE=-90.0D0
          ELSE
              ANGLE=90.0D0
          END IF

      END IF

C.......................................!.................
      Qkgsec=QIN
      IF ((ISYM.eq.0).AND.(BEARING.EQ.1))  Qkgsec=QOUT   ! Qout=Qleft+Qright
      IF (MODEL.eq.1) Qkgsec=2.0D0*Qkgsec                ! 2 row HJB
      Qkgmin=Qkgsec*60.0D0

      Iwrite=1
      IU=6
      CALL BEEPER
C.......................................!.!Qright=Qin,Qleft=Qo-Qright.
  1   WRITE (IU, 590)
      WRITE (IU, 600)
      WRITE (IU, 666)
      WRITE (IU, 600)
c   !.....................................!
      IF (BEARING.EQ.3) THEN
        WRITE(IU,501) Qkgsec, Qkgmin
        IF (IFULL.NE.1) WRITE (IU,502) QOUT
      ELSE
        WRITE (IU, 503) Qkgsec, Qkgmin
      END IF
      IF ((BEARING.EQ.1).AND.(ISYM.EQ.0).AND.(MODEL.EQ.2)) THEN
        WRITE (IU, 550) (QOUT-QIN), QIN
      END IF
C.......................................!.................
      WRITE (IU, 600)
      WRITE (IU, 100) FX, FY, LOAD, ANGLE
      IF (MODEL.EQ.2) WRITE (IU, 105) MX, MY
      WRITE (IU, 600)
      WRITE (IU, 110) TO



c .......................................................................
      IF (Iwrite.eq.2) RETURN
      IF (DEVICE.EQ.1) THEN
          Iwrite=2
          IU=1
          GOTO 1
      END IF

C.......................................!.............................

  666 FORMAT (' | hydroflex: TOTAL RESULTS FOR BEARING:',39X,'|')

  110 FORMAT (' |    Torque on Film Lands=', E12.5E2, ' N-m', 36X, '|',
     +        /, ' +', 77('-'), '+')

  200 FORMAT (' ', '+', 77('-'), '+')

  100 FORMAT (' |    Fx=', E12.5E2, ' N     Fy=', E12.5E2, ' N',
     +        34X, '|', /, ' |    Load=', E12.5E2,
     + ' N   Load Angle = ',
     +        F8.3, ' deg', 26X, '|')
  105 FORMAT (' |    Mx=', E12.5E2, ' Nm    My=', E12.5E2, ' Nm',
     +        33X, '|')

  400 FORMAT (' ', '|     MAX Rise in Pressure within recess: ',
     +        E12.5E2, 24X, '|')

  501 FORMAT (' |    Side Mass Flows=', E12.5E2,' Kg/s = ',E12.5E2,
     +        ' Kg/min', 18X, '|')
  502 FORMAT (' |    Lead edge Flows=', E12.5E2,' Kg/s   ',
     +                   37X, '|')
  503 FORMAT (' |    Mass Flows=', E12.5E2, ' Kg/s = ', E12.5E2,
     +        ' Kg/min', 23X, '|')

  550 FORMAT (' |    (Left Flow)=', E12.5E2, ' Kg/s  (Right Flow)=',
     +        E12.5E2, ' Kg/s',12X,'|')
  555 FORMAT (' |    Mass Flow/Rhos=',E12.5E2, ' lt/min', 38X, ' |')

  590 FORMAT (' ')
  600 FORMAT (' +', 77('-'), '+')

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine PRINT1PAD: PRINTS TOTAL FORCE & MOMENT COEFFICIENTS         **
C **                                                                         **
C **                                                                         **
C *****************************************************************************
          SUBROUTINE PRINT1PAD(DEVICE,RPM,PS,PA)
C
      IMPLICIT NONE
c.........................................................................
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /RECPAR/ HREC, VSUP, BETA
      COMMON /STIFT/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPT/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERT/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIAT/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMAT/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INEAT/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /FLMOM/ IMOMFLAG

      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION frequ,sigma,l1,res,hrec,vsup,beta,AXO,AYO,ZO,
     +                emu,rho,rpm,ps,pa,pc,cd,dorif,losxsi,alpha,
     +                RSINPK, RCOSPK, IPAD

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,ISYM,ICSTEP,
     +        DEVICE, ICASE, NCASE, INERFREQ, IMOMFLAG, TILT

C ------------------------------------------------------------------------
C    local variables
C ------------------------------------------------------------------------
      DOUBLE PRECISION omega,freq,Keq,I1,I2,I4,phi02,dumy,whirl,PI,
     +                 whirl1, whirl2,L,Lold,EPS
      INTEGER UN, Iset, Iwfr,I, INERLDUM
c.........................................................................
c calculation of whirl frequency ratio for freq=omega (synchrpnous excita)
c.........................................................................
      PI=DACOS(-1.0D00)
      Iwfr=0
      EPS=1.0D-5
      INERLDUM = 0   ! 9/21/93

c hydroflex calculates impedance coefficients and therefore does not provide
c           values of inertia force coefficients.

      IF (RPM.EQ.(0.0D0)) GOTO 100

      omega=rpm/60.0D0
      freq=frequ*2.0D0*PI

c    !............................................!
      IF (DABS(frequ-omega).GT.(1.0D0)) GOTO 100
c    !............................................!
        Iwfr=1
        Keq=(Kxxd*Cyyd+Cxxd*Kyyd-Cyxd*Kxyd-Cxyd*Kyxd)/(Cxxd+Cyyd)

        dumy=((Keq-Kxxd)*(Keq-Kyyd)-Kxyd*Kyxd)
        dumy=dumy/(Cxxd*Cyyd-Cxyd*Cyxd)
        phi02=dumy/freq/freq


c      !............................................!
        IF ((INERLDUM.EQ.1).AND.(TILT.EQ.0)) THEN
c      !............................................!
          I1=(Cyxd*Mxyd+Cxyd*Myxd)/(Cxxd+Cyyd)
          I2=(Kxyd*Myxd+Kyxd*Mxyd-I1*(Kxxd+Kyyd)+2*Keq*I1)
          I2=I2/(Cxxd*Cyyd-Cxyd*Cyxd)
          I4=(I1*I1-Mxyd*Myxd)/(Cxxd*Cyyd-Cxyd*Cyxd)
          I4=I4*freq*freq
          dumy=1.0-I2

          i=0
          L=Phi02
  17      Lold=L
          L=(Phi02+L*L*I4)/dumy
          i=i+1
          IF (i.gt.10) goto 18
          IF (DABS(L-Lold).LT.EPS) goto 18
          goto 17

  18      IF (L.LT.0.0D0) THEN
            whirl=0.0D0
          ELSE
            whirl=+DSQRT(+L)
          END IF


c      !.....................!
        ELSE                 ! INERL=0
c      !.....................!

          IF (phi02.lt.0.0D0) THEN
             whirl=0.0D0
          ELSE
             whirl=+DSQRT(+Phi02)
          END IF


c      !.....................!
        END IF
c      !.....................!



c ...............................................................!


 100  UN=6
      INERFREQ=0
      IF ((INERLDUM.EQ.1).AND.(FREQU.NE.(0.0))) INERFREQ=1
      IF (TILT.EQ.1) INERFREQ=0
      !INERFREQ=0        !####

c     We have set INERFREQ=0 always, Then force coefficients
c     Kij & Cij really represent the real part and imaginary part of
c     the complex dynamic force impedance.


      Iset=1
      INERLDUM = 1

 101  WRITE (Un,1800) frequ,Res
c     ............................. FORCE COEFFICs. due to DISPLACEMENTS
      WRITE (Un, 1900) KXXD, KYXD, KYYD, KXYD
      IF (FREQU.GT.(0.0)) WRITE (Un, 2000) CXXD, CYXD, CYYD, CXYD
      IF (INERFREQ.EQ.1) WRITE (Un, 2100) MXXD, MYXD, MYYD, MXYD
      IF (Iwfr.EQ.1) WRITE (Un, 1999) KEQ, WHIRL
      WRITE (Un, 2200)




      !IF (FREQU.GT.0) WRITE (61, *) RPM, KXXD, KYXD, KYYD, KXYD,CXXD, CYXD, CYYD, CXYD

c     ...........................................................!
      IF (TILT.EQ.1) GOTO 44
c     ...........................................................!
      IF (MODEL.eq.1) GOTO 44
c     ...........................................................!

c     ............................. MOMENT COEFFICs. due to DISPLACEMENTS
 11   IF ( (ISYM.EQ.1).AND.(ZO.EQ.0.0D0) ) GOTO 33
      WRITE (Un, 1901) KmXXD, KmYXD, KmYYD, KmXYD
      IF (FREQU.GT.(0.0)) WRITE (Un, 2001) CmXXD, CmYXD, CmYYD, CmXYD
      IF (INERFREQ.EQ.1) WRITE (Un, 2101) MmXXD, MmYXD, MmYYD, MmXYD
      WRITE (Un, 2200)

      IF (IMOMFLAG.EQ.0) GOTO 44
c     ............................. FORCE COEFFICs. due to ANGLE DISP
 22   WRITE (Un, 1902) KXXA, KYXA, KYYA, KXYA
      IF (FREQU.GT.(0.0)) WRITE (Un, 2002) CXXA, CYXA, CYYA, CXYA
      IF (INERFREQ.EQ.1) WRITE (Un, 2102) MXXA, MYXA, MYYA, MXYA
      WRITE (Un, 2200)
c     ............................. MOMENT COEFFICs. due to ANGLE DISP
 33   IF (IMOMFLAG.EQ.0) GOTO 44
      WRITE (Un, 1903) KmXXA, KmYXA, KmYYA, KmXYA
      IF (FREQU.GT.(0.0)) WRITE (Un, 2003) CmXXA, CmYXA, CmYYA, CmXYA
      IF (INERFREQ.EQ.1) WRITE (Un, 2103) MmXXA, MmYXA, MmYYA, MmXYA
      WRITE (Un, 2200)
c     ...................................................................
 44   IF (Iset.eq.2) RETURN
c     ......................!
      IF (DEVICE.eq.1) THEN
          Iset=2
          UN=1
          GOTO 101
      END IF

c.........................................................................

 1800 format (' ', 3X,'frequency(w):',E12.5E2,'Hz',
     + 2X,'Squeeze Reynolds#(Res):(p/u)wC**2 =', E11.4E2,
     + /, 1X, 79('.'))

 1900 FORMAT (' ',3X,'FORCE Coefficients due to Displacements:',/,
     +        4X,'Kxx=', E11.4E2, 2X, 'Kyx=', E11.4E2,
     +        2X, 'Kyy=', E11.4E2, 2X,'Kxy=', E11.4E2, ' (N/m)')
 2000 FORMAT (' ',3X,  'Cxx=', E11.4E2, 2X, 'Cyx=', E11.4E2,
     +        2X, 'Cyy=', E11.4E2, 2X,'Cxy=', E11.4E2, ' (Ns/m)')
 2100 FORMAT (' ',3X,  'Mxx=', E11.4E2, 2X, 'Myx=', E11.4E2,
     +        2X, 'Myy=', E11.4E2, 2X,'Mxy=', E11.4E2, ' (Kg)')

 1901 FORMAT (' ',3X,'MOMENT Coefficients due to Displacements:',/,
     +        4X,'Kxx=', E11.4E2, 2X, 'Kyx=', E11.4E2,
     +        2X, 'Kyy=', E11.4E2, 2X,'Kxy=', E11.4E2, ' (N)')
 2001 FORMAT (' ',3X,  'Cxx=', E11.4E2, 2X, 'Cyx=', E11.4E2,
     +        2X, 'Cyy=', E11.4E2, 2X,'Cxy=', E11.4E2, ' (N-s)')
 2101 FORMAT (' ',3X,  'Mxx=', E11.4E2, 2X, 'Myx=', E11.4E2,
     +        2X, 'Myy=', E11.4E2, 2X,'Mxy=', E11.4E2, ' (N-s2)')

 1902 FORMAT (' ',3X,'FORCE Coefficients due to ANGLE Rots. :',/,
     +        4X,'Kxx=', E11.4E2, 2X, 'Kyx=', E11.4E2,
     +        2X, 'Kyy=', E11.4E2, 2X,'Kxy=', E11.4E2, ' (N)')
 2002 FORMAT (' ',3X,  'Cxx=', E11.4E2, 2X, 'Cyx=', E11.4E2,
     +        2X, 'Cyy=', E11.4E2, 2X,'Cxy=', E11.4E2, ' (N-s)')
 2102 FORMAT (' ',3X,  'Mxx=', E11.4E2, 2X, 'Myx=', E11.4E2,
     +        2X, 'Myy=', E11.4E2, 2X,'Mxy=', E11.4E2, ' (N-s2)')

 1903 FORMAT (' ',3X,'MOMENT Coefficients due to ANGLE Rots. :',/,
     +        4X,'Kxx=', E11.4E2, 2X, 'Kyx=', E11.4E2,
     +        2X, 'Kyy=', E11.4E2, 2X,'Kxy=', E11.4E2, ' (N-m)')
 2003 FORMAT (' ',3X,  'Cxx=', E11.4E2, 2X, 'Cyx=', E11.4E2,
     +        2X, 'Cyy=', E11.4E2, 2X,'Cxy=', E11.4E2, ' (N-m-s)')
 2103 FORMAT (' ',3X,  'Mxx=', E11.4E2, 2X, 'Myx=', E11.4E2,
     +        2X, 'Myy=', E11.4E2, 2X,'Mxy=', E11.4E2, ' (N-m-s2)')

 1999 FORMAT (' ',3X,'Keq=',E11.4E2,'N/m;  WFR=', E12.5E2 )
 2200 FORMAT (' ',  79('_'))

      END

c.........................................................................
C modified 2/15/94  BY DR. LUIS SANANDRES / TEXAS A&M University
c.........................................................................