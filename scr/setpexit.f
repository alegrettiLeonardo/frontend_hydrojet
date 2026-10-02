c  ####   ######   #####  #####   ######  #    #     #     #####          ######
c #       #          #    #    #  #        #  #      #       #            #
c  ####   #####      #    #    #  #####     ##       #       #            #####
c      #  #          #    #####   #         ##       #       #     ###    #
c #    #  #          #    #       #        #  #      #       #     ###    #
c  ####   ######     #    #       ######  #    #     #       #     ###    #

c #    #   #   #  #####   #####    ####   ######  #       ######  #    #
c #    #    # #   #    #  #    #  #    #  #       #       #        #  #
c ######     #    #    #  #    #  #    #  #####   #       #####     ##
c #    #     #    #    #  #####   #    #  #       #       #         ##
c #    #     #    #    #  #   #   #    #  #       #       #        #  #
c #    #     #    #####   #    #   ####   #       ######  ######  #    #

C setpexit.f > hydroflex.f Dr. Luis San Andres, TexasA&MUniv. 1994
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR I
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C *****************************************************************************
C **                                                                         **
C **  Subroutine SetPright                                                   **
C **                                                                         **
C **  SETPright:Imposses non-uniform pressure at bearing right end  y=LenghtR**
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETPright

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)

      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c    !.................................................................!
      DOUBLE PRECISION P, XP, XU, YP, YV,
     +                 PRCOEFC, PRCOEFS, Pexit,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2 

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        IPRuni, IPLuni, NPRcs, NPLcs, I, J, BEARING

C ----------------------------------------------------------------------------
C --  SETPright code                                                        --
C     Fourier series expansion of discharge pressure at right plane y=L/D
C ----------------------------------------------------------------------------
           DO J=1, NXT                 ! set nonuniform pressure at exit
              Pexit=PRCOEFC(0)
              DO i=1, NPRCs
                Pexit=Pexit+PRCOEFC(i)*DCOS(i*XP(J))+
     +                      PRCOEFS(I)*DSIN(i*XP(J))
              END DO

              P(J,NYI)=(Pexit-Pa)/PSA

           END DO

           END


C *****************************************************************************
C **                                                                         **
C **  Subroutine SetPleft                                                    **
C **                                                                         **
C **  SETPleft :Imposses non-uniform pressure at bearing left end y=-LengthL **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETPleft

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)

      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c    !.................................................................!
      DOUBLE PRECISION P, XP, XU, YP, YV,
     +                 PLCOEFC, PLCOEFS, Pexit,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        IPRuni, IPLuni, NPRcs, NPLcs, I, J, BEARING 

C ----------------------------------------------------------------------------
C --  SETPleft code                                                        --
C     Fourier series expansion of discharge pressure at left plane 
C ----------------------------------------------------------------------------
           DO J=1, NXT                 ! set nonuniform pressure at exit
              Pexit=PLCOEFC(0)
              DO i=1, NPLcs
                Pexit=Pexit+PLCOEFC(i)*DCOS(i*XP(J))+
     +                      PLCOEFS(I)*DSIN(i*XP(J))
              END DO

              P(J,-NYI)=(Pexit-Pa)/PSA

           END DO

           END

C ----------------------------------------------------------------------------

C *****************************************************************************
C **                                                                         **
C **  Subroutine ReadPright                                                  **
C **                                                                         **
C **  ReadPright:reads Fourier pressure coefficients at right end y=LenghtR  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ReadPright(PRIGHT)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)

      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c    !.................................................................!
      DOUBLE PRECISION PRCOEFC, PRCOEFS, PRIGHT 

      INTEGER I, J, IPRuni, IPLuni, NPRcs, NPLcs, BEARING 

C ----------------------------------------------------------------------------
C --  ReadPright code                                                      --
C     Fourier series expansion of discharge pressure at right plane 
C ----------------------------------------------------------------------------
      CALL BEEPER
      WRITE (6, *) '--------------------------------------------------'
      WRITE (6, *) 'DISHARGE Pressure at Right END, Y=LengthR         '
      WRITE (6, *) '--------------------------------------------------'
      WRITE (6, *) 'SELECT IPRuni: 1= Uniform Pressure 0=Nonuniform   '
      WRITE (6, *) '--------------------------------------------------'
      CALL BEEPER

      WRITE (6, 17) IPRuni
 17   FORMAT ('$', 'ENTER IPRuni [Default =',I2,']: ')
      CALL ENTERINT(IPRUNI)

c   !------------------------------------------!
      IF (IPRUNI.EQ.1 ) THEN                   ! Uniform Discharge Pressure
c   !------------------------------------------!
      WRITE (6, *) 'PR (REAL*8):RIGHT Discharge pressure(PR), [N/m2]'
      WRITE (6, 191) PRIGHT

  191 FORMAT ('$', 'ENTER PR [Default(PR)= ',E14.7E2,' N/m2]: ')

      CALL ENTERVAL(PRIGHT)

c   !------------------------------------------!
      ELSE
c   !------------------------------------------! Non Uniform Discharge Pressure

      WRITE (6,*) '------------ NON UNIFORM PRESSURE ----------'
      WRITE (6,*) 'Fourier Coefficients for Pressure Discharge '
      WRITE (6,*) 'function equal to:                          '
      WRITE (6,*) 'PR = PRo  + SUM {PRci cos(X) + PRsi sin(X)},'
      WRITE (6,*) 'i=1,2,....   NPRcs <= NYI'

      WRITE (6, 27) NPRcs
 27   FORMAT ('$', 'ENTER NPRcs [Default =',I2,']: ')
      CALL ENTERVAL(NPRCS)

      IF (NPRcs.GT.MAXNYI) NPRcs=MAXNYI
      CALL BEEPER

      WRITE (6, 28) PRCOEFC(0)
 28   FORMAT ('$', 'ENTER PRo [Default(PRo)= ',E14.7E2,' N/m2]: ')
      CALL ENTERVAL(PRIGHT)
      PRCOEFC(0)=PRIGHT

      IF (NPRcs.eq.0) THEN
         IPRuni=1
         GOTO 77
      END IF

      DO I=1, NPRcs

        WRITE (6,29) I, I
 29     FORMAT ('$','ENTER cos and sin pressure coefs: ',
     +              'PRc(', I2 ,'), PRs(', I2 ,')' )

        READ (5, *) PRCOEFC(I),PRCOEFS(I)

      END DO

      WRITE (6,*) 
      WRITE (6,*) '------------------------------------------'
      WRITE (6,*) 'RIGHT Fourier Pressure Coefficients are:'
      WRITE (6,41) PRIGHT, NPRcs 
      DO I=1, NPRcs
       WRITE (6,42) I, PRCOEFC(I),PRCOEFS(I)
      END DO
      WRITE (6,*) '-------------------------------------------'

 41   FORMAT (' ',3X,   'PRo=',E12.5E2,'N/m2  NPRcs=',I2)
 42   FORMAT (' ',I2,1X,'PRc=',E12.5E2,2X,'PRs=',E12.5E2,'N/m2')     
c   !------------------------------------------!
 77   END IF
c   !------------------------------------------! Non Uniform Discharge Pressure


      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine ReadPleft                                                  **
C **                                                                         **
C **  ReadPright:reads Fourier pressure coefficients at left end Y=-LengthL  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ReadPLeft(PLeft)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)

      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING
c    !.................................................................!
      DOUBLE PRECISION PLCOEFC, PLCOEFS, PLeft 

      INTEGER I, J, IPRuni, IPLuni, NPRcs, NPLcs, BEARING

C ----------------------------------------------------------------------------
C --  ReadPleft code                                                  --
C     Fourier series expansion of discharge pressure at left plane 
C ----------------------------------------------------------------------------
      CALL BEEPER
      WRITE (6, *) '--------------------------------------------------'
      WRITE (6, *) 'DISHARGE Pressure at Left END, Y=-LengthL         '
      WRITE (6, *) '--------------------------------------------------'
      WRITE (6, *) 'SELECT IPLuni: 1= Uniform Pressure 0=Nonuniform   '
      WRITE (6, *) '--------------------------------------------------'
      CALL BEEPER

      WRITE (6, 17) IPLuni
 17   FORMAT ('$', 'ENTER IPRuni [Default =',I2,']: ')
      CALL ENTERINT(IPLUNI)

c   !------------------------------------------!
      IF (IPLUNI.EQ.1 ) THEN                   ! Uniform Discharge Pressure
c   !------------------------------------------!
      WRITE (6, *) 'PL (REAL*8):LEFT Discharge pressure(PL), [N/m2]'
      WRITE (6, 191) PLEFT

  191 FORMAT ('$', 'ENTER PL [Default(PL)= ',E14.7E2,' N/m2]: ')

      CALL ENTERVAL(PLEFT)

c   !------------------------------------------!
      ELSE
c   !------------------------------------------! Non Uniform Discharge Pressure

      WRITE (6,*) '------------ NON UNIFORM PRESSURE ----------'
      WRITE (6,*) 'Fourier Coefficients for Pressure Discharge '
      WRITE (6,*) 'function equal to:                          '
      WRITE (6,*) 'PL = PLo  + SUM {PLci cos(X) + PLsi sin(X)},'
      WRITE (6,*) 'i=1,2,....   NPLcs <= NYI'

      WRITE (6, 27) NPLcs
 27   FORMAT ('$', 'ENTER NPLcs [Default =',I2,']: ')
      CALL ENTERINT(NPLCS)

      IF (NPLcs.GT.MAXNYI) NPRcs=MAXNYI
      CALL BEEPER

      WRITE (6, 28) PLCOEFC(0)
 28   FORMAT ('$', 'ENTER PLo [Default(PLo)= ',E14.7E2,' N/m2]: ')
      CALL ENTERVAL(PLEFT)

      PLCOEFC(0)=PLEFT

      IF (NPLcs.eq.0) THEN
         IPLuni=1
         GOTO 77
      END IF

      DO I=1, NPLcs

        WRITE (6,29) I, I
 29     FORMAT ('$','ENTER cos and sin pressure coefs: ',
     +              'PLc(', I2 ,'), PLs(', I2 ,')' )

        READ (5, *) PLCOEFC(I),PLCOEFS(I)

      END DO

      WRITE (6,*) 
      WRITE (6,*) '------------------------------------------'
      WRITE (6,*) 'LEFT Fourier Pressure Coefficients are:'
      WRITE (6,41) PLEFT, NPLcs 
      DO I=1, NPLcs
       WRITE (6,42) I, PLCOEFC(I),PLCOEFS(I)
      END DO
      WRITE (6,*) '-------------------------------------------'
      CALL BEEPER

 41   FORMAT (' ',3X,   'PLo=',E12.5E2,'N/m2  NPLcs=',I2)
 42   FORMAT (' ',I2,1X,'PLc=',E12.5E2,2X,'PLs=',E12.5E2,'N/m2')     
c   !------------------------------------------!
 77   END IF
c   !------------------------------------------! Non Uniform Discharge Pressure


      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine WritePright                                                 **
C **                                                                         **
C **  WritePright:prints Fourier pressure coefficients at right end y=LENGTHR**
C **                                                                         **
C *****************************************************************************

      SUBROUTINE WritePright(UNIT)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PRexit/ PRCOEFC(0:MAXNYI), PRCOEFS(1:MAXNYI)

      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
c    !.................................................................!
      DOUBLE PRECISION PRCOEFC, PRCOEFS

      INTEGER I, UNIT, IPRuni, IPLuni, NPRcs, NPLcs 

C ----------------------------------------------------------------------------
C     Fourier series expansion of discharge pressure at right plane 
C ----------------------------------------------------------------------------

      IF (IPRuni.eq.1) RETURN

      WRITE (UNIT,15)
 15   FORMAT(' ',3X,
     +       '--NON UNIFORM RIGHT END PRESSURE(Y=LengthR)------',/,
     +       4X,'Fourier Coefficients for Pressure Discharge ',
     +       'function equal to:',/,3X,                          
     +       'PR = PRo  + SUM {PRci cos(X) + PRsi sin(X)} ',2X,
     +       'i=1,2,....   NPRcs ',/,4X,20('. '))

      WRITE (UNIT,41) PRCOEFC(0), NPRcs 
      DO i=1, NPRcs
       WRITE (UNIT,42) I, PRCOEFC(I),PRCOEFS(I)
      END DO
      WRITE (UNIT,44)

c   !------------------------------------------! Non Uniform Discharge Pressure

 41   FORMAT (' ',3X,   'PRo=',E12.5E2,'N/m2  NPRcs=',I2)
 42   FORMAT (' ',I2,1X,'PRc=',E12.5E2,2X,'PRs=',E12.5E2,'N/m2')     
 44   FORMAT (' ', 3X, 77('-'))

      RETURN


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine WritePleft                                                  **
C **                                                                         **
C **  WritePleft:prints Fourier pressure coefficients at left end y=-LengthL **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE WritePleft(UNIT)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PLexit/ PLCOEFC(0:MAXNYI), PLCOEFS(1:MAXNYI)

      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
c    !.................................................................!
      DOUBLE PRECISION PLCOEFC, PLCOEFS

      INTEGER I, UNIT, IPRuni, IPLuni, NPRcs, NPLcs 

C ----------------------------------------------------------------------------
C     Fourier series expansion of discharge pressure at right plane y=L/D
C ----------------------------------------------------------------------------

      IF (IPLuni.eq.1) RETURN

      WRITE (UNIT,15)
 15   FORMAT(' ',3X,
     +       '--NON UNIFORM LEFT END PRESSURE(Y=-LengthL)----',/,
     +       4X,'Fourier Coefficients for Pressure Discharge ',
     +       'function equal to:',/,3X,                          
     +       'PL = PLo  + SUM {PLci cos(X) + PLsi sin(X)} ',2X,
     +       'i=1,2,....   NPLcs ',/,4X,20('. '))

      WRITE (UNIT,41) PLCOEFC(0), NPLcs 
      DO i=1, NPLcs
       WRITE (UNIT,42) I, PLCOEFC(I),PLCOEFS(I)
      END DO
      WRITE (UNIT,44)

c   !------------------------------------------! Non Uniform Discharge Pressure

 41   FORMAT (' ',3X,   'PLo=',E12.5E2,'N/m2  NPLcs=',I2)
 42   FORMAT (' ',I2,1X,'PLc=',E12.5E2,2X,'PLs=',E12.5E2,'N/m2')     
 44   FORMAT (' ', 3X, 77('-'))

      RETURN

      END

C ----------------------------------------------------------------------------
C setpexit.f: Last revised 12/31/94  by Dr. Luis San Andres
C ----------------------------------------------------------------------------