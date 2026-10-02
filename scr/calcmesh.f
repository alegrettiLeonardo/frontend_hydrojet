c
c  ####     ##    #        ####   #    #  ######   ####   #    #   #####
c #    #   #  #   #       #    #  ##  ##  #       #       #    #     #
c #       #    #  #       #       # ## #  #####    ####   ######     #
c #       ######  #       #       #    #  #            #  #    #     #     ###
c #    #  #    #  #       #    #  #    #  #       #    #  #    #     #     ###
c  ####   #    #  ######   ####   #    #  ######   ####   #    #     #     ###
C
c #    #   #   #  #####   #####    ####   ######  #       ######  #    #
c #    #    # #   #    #  #    #  #    #  #       #       #        #  #
c ######     #    #    #  #    #  #    #  #####   #       #####     ##
c #    #     #    #    #  #####   #    #  #       #       #         ##
c #    #     #    #    #  #   #   #    #  #       #       #        #  #
c #    #     #    #####   #    #   ####   #       ######  ######  #    #

C calcmesht.f > hydroflex code Drs. Luis San Andres, TexasA&MUniv. 1994
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR I
c Technical monitor: Mr. James Walker, NASA Lewis Research Center
c
C
C FOR PAD:
C
C Length=LengthL+LengthR

C .....................................................
C Pright & Cright  at Y=LengthR
C
C ->X1<................ Lpad ......................>
C   ------------------------------------------------
C   ^        -> Lrec                              ^
C   |y      X1r <--->                              | LengthR
C   |  -----    ----- ^  -----    -----    -----   |
C   |  |   |    |   | |  |   |    |   |    |   |   v
C   |..|   |    |   | |Ar|   |    |   |    |   |   ----> x
C   |  |   |    |   | |  |   |    |   |    |   |   ^
C   |  -----    ----- v  -----    -----    -----   | LenghtL
C   |                                              v
C   |------------------------------------------------
C
C Pleft & Cleft  at Y=LengthL

C NOTE: for 2PI (360) circular bearing, ALWAYS start at a recess EDGE

C INPUTPADS: User input data for bearing geometry :

C DIAMeter, axial LENGTH, Lr, Ll=L-Lr, recess axial length AR,
C Number of pads NPAD, X1(Npad) = position of leading edge of pad [m]
C                      LPAD(Npad) = circ. length of pad [m]
C
C FOR each Pad, k=1, Npad, INPUT
C  Nrec(k)=NPOCKET = Number of recesses/pockets on pad,
C  X1rec(i,k) = position of leading edge of recess RELATIVE to PAD Leading edge
C  Lrec(i,k) = cir. length of recess on pad [m]; i=1,2...,NPOCKET
C  PaDorif(i,k) = orifice diameter [m], i=1,2,... NPOCKET
C
C .............................................................................
C

C *****************************************************************************
C **                                                                         **
C **  Subroutine InputPads                                                   **
C **                                                                         **
C **  INPUTPADS: KBD INPUT of PAD BEARING GEOMETRY & DIMENSIONS              **
C **                                                                         **
C *****************************************************************************
C
      SUBROUTINE INPUTPADS(RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /TGROOVE/ LEMDA,DELTA

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c     ................................................................
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD, AR, HREC,
     +                 X1, LPAD, X1r,Lrec, HRECU, VSUP, BETA,
     +                 LENGTHL, LENGTHR, PaDorif,ROTPAD,
     +                 PRELOAD, PRELOADB, OFFSET, ROTDEL,INERPAD,
     +                 KSTPAD, CDAPAD, LEMDA, DELTA,RPM,PS,PA

      INTEGER NPAD, NREC, NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI,
     +        NYI, NXT, ISYM, ICSTEP, IFULL, INERL, INERP, ITURB,
     +       INTER, ICAV, MODEL,BEARING, NLOBES, TILTPAD

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION ZERO, R, LMAX,LPmax, Conver, Lrmax, CROTPAD,
     +                 X1pmi,X1rmi,X2,X2r,SLPAD,SLREC,IPAD, KROTPAD
      INTEGER I, J, II, JJ, N, K, POINTER, KKK, Ndumy
      CHARACTER*1 YN
      CHARACTER*30 BCASE
C ----------------------------------------------------------------------------
C --  INPUTPADS code                                                      --
C ----------------------------------------------------------------------------
C Unit 42 is used as temporal storage for translations between formats.
C.........................................................................
      ZERO=0.0D0
      LMAX=360.00D0   ! degrees MAX PAD angular LENGTH

C    !.........................!
      IF (BEARING.EQ.1) THEN   !==> FOR HYDROSTATIC BEARINGS
C    !.........................!
        IF (MODEL.EQ.2) BCASE='1 ROW HYDROSTATIC BEARING'
        IF (MODEL.EQ.1) BCASE='2 ROW HYDROSTATIC BEARING'
C    !.........................!
      ELSE IF (BEARING.EQ.2) THEN   !==> FOR ANNULAR SEAL
C    !.........................!
        BCASE='ANNULAR SEAL'
        NPAD=1
        X1(1)=ZERO
        LPAD(1)=360.0D0
        NREC(1)=0
        SLPAD=360.0D0
        NLA=NYI
        NPA=1
        NLC=NXT
        NPC=0
        MODEL=2
C    !.........................!
      ELSE IF (BEARING.EQ.3) THEN   !==> FOR PLAIN BEARING
C    !.........................!
        BCASE='CYLINDRICAL BEARING'
        DO J=1, NPAD
          NREC(J)=0
        END DO
        NLA=NYI
        NPA=1
        NLC=NXT
        NPC=0
        MODEL=2
C    !.........................!
      END IF
C    !.........................!



C .......................................! MODIFY INPUT DATA
    1 CONTINUE
C .......................................! MODIFY INPUT DATA
      SLPAD=ZERO
      DO K=1, NPAD
        SLPAD=SLPAD+LPAD(K)
      END DO
C ................................... PRINT DEFAULT VALUES
      CALL BEEPER
      IF (TILTPAD.EQ.1) THEN
        WRITE (6,7011) BCASE
      ELSE
        WRITE (6,7012) BCASE
      END IF

      WRITE (6,7021) DIAM, LENGTH, LENGTHR, LENGTHL

      IF (BEARING.EQ.1) WRITE (6,7022) AR


c    !........................!
      IF (TILTPAD.EQ.1) THEN
c    !........................!
      WRITE (6,7023) NPAD,SLPAD
        DO K=1, NPAD
        WRITE (6,7031) K, X1(K),LPAD(K),ROTPAD(K),
     +                 INERPAD(K),KSTPAD(K),CDAPAD(K)
        END DO
        WRITE (6,7033) NPAD, PRELOAD, OFFSET
c    !........................!
      ELSE
c    !........................!
        WRITE (6,7024) NPAD,SLPAD
        DO K=1, NPAD
        WRITE (6,7032) K, X1(K),LPAD(K)
        END DO
        WRITE (6,7034) NLOBES, PRELOAD, OFFSET
c    !........................!
      END IF
c    !........................!



      DO K=1, NPAD
        NPOCKET=NREC(K)
        IF (NPOCKET.GT.0) THEN
          WRITE (6,704) K
          DO I=1, NPOCKET
          WRITE (6,705) I,X1r(I,K),LRec(I,K),PaDorif(I,K)
          END DO
        ELSE
          WRITE (6,7040) K
        END IF
      END DO
      WRITE (6,706)

      IF (BEARING.EQ.1) THEN
        WRITE (6,7081) NLA, NPA,  NLC, NPC
      ELSE
        WRITE (6,7082) NYI, NXT
        AR=ZERO ! POCKET AXIAL LENGTH
      END IF

      WRITE (6,709)
      POINTER = 1
C.............................................................
  17  WRITE (6, *) 'ENTER SELECTED CHOICE:'
      CALL ENTERINT(POINTER)
      IF (POINTER.LE.0) GOTO 1
      IF (POINTER.GT.9) POINTER = 9
      GO TO (90,100,105,115,1205,80,111,140,999), POINTER

      POINTER=9
      RETURN

C.............................................................
C CLEARANCE function and recess depth/supply line:
c Bearing clearances for non-uniform bearing
C C:C(Y): STRAIGHT, TAPERED, NON-UNIFORM via SPLINE fit
C................................................

   80 POINTER=6
      CALL INSPL     ! => on kbdinput.f

      IF (BEARING.GT.1) THEN
         HREC=0.0D0
         VSUP=0.0D0
         GOTO 1
      END IF
C..........................................................
C
C RECESS DEPTH
C.................................

      WRITE (6, *) 'HREC (REAL*8): Recess depth (Hrec), [m]'
      WRITE (6, 81) HREC

  81  FORMAT ('$', 'ENTER Hrec [Default (Hrec)= ',E14.7E2,']: ')
      CALL ENTERVAL(HREC)
      HREC=DABS(HREC)
      HRECU=HREC
      WRITE (6, *) ' '

C...................................................................
C
C VSUP: Volume of supply line to recess
C.......................................

      WRITE (6, *) 'VSUP (REAL*8): Volume of orifice supply line,[m3]'
      WRITE (6, 82) VSUP
  82  FORMAT ('$', 'ENTER Vsup [Default (Vsup)= ',E14.7E2,']: ')
      CALL ENTERVAL(VSUP)
      VSUP=DABS(VSUP)
      GOTO 1


C.............................................................
C DIAMETER:
C...................................
  90  POINTER=1
      CALL BEEPER
      WRITE (6, *) 'DIAM (REAL*8): Journal diameter(D), [m]'
      WRITE (6, 91) DIAM
   91 FORMAT ('$', 'ENTER Diam [Default (D)= ',E12.5E2,']: ')
      CALL ENTERVAL(DIAM)
      DIAM=DABS(DIAM)

      IF (DIAM.EQ.ZERO) THEN
        WRITE (6,*) 'ERROR: DIAMETER = ZERO [m] REVISE DATA'
        GOTO 90
      END IF

      LMAX=360.00D0   !MAX PAD LENGTH in degrees
      R=DIAM/2.0D0
      WRITE (6, *) ' '
      GOTO 1
C.................................................................
C BEARING AXIAL LENGTH
C...................................
  100 POINTER=2
      CALL BEEPER
      WRITE (6, *) 'LENGTH (REAL*8): BEARING Axial length (L), [m]'

      IF (MODEL.eq.1) THEN
        WRITE(6, *) 'FOR DOUBLE ROW HJB, ENTER HALF THE BEARING LENGTH'
      END IF

      WRITE (6, 101) LENGTH
  101 FORMAT ('$', 'ENTER Length [Default (L)=',E12.5E2,']: ')
      CALL ENTERVAL(LENGTH)
      LENGTH=DABS(LENGTH)

      IF (LENGTH.EQ.ZERO) THEN
        WRITE (6,*) 'ERROR: LENGTH = ZERO [m] REVISE DATA'
        GOTO 100
      END IF
      WRITE (6, *) ' '
      LD=LENGTH/DIAM
C .............................................................
      IF (BEARING.EQ.2) THEN    ! FOR ANNULAR SEAL
         LENGTHR=LENGTH
         LENGTHL=0.0D0
         ISYM=0
         GOTO 1
      END IF

      WRITE (6, *) '$ IS BEARING SYMMETRIC (?) ENTER Y or N'
      YN='Y'
      READ (5, 265) YN
 265  FORMAT (A1)

      IF ((YN.EQ.'Y').OR.(YN.EQ.'y').OR.(YN.EQ.' ')) THEN
           LENGTHL=LENGTH/2.0D0
           LENGTHR=LENGTHL
           WRITE (6, *) ' '
           GOTO 1
      END IF
C .............................................................

 102  WRITE (6,103) LENGTHR
 103  FORMAT ('$', 'ENTER RIGHT Length FROM ORIFICE CENTER[(Lr)='
     +             ,E12.5E2,']:')
      CALL ENTERVAL(LENGTHR)
      LENGTHR=DABS(LENGTHR)

      LENGTHL=LENGTH-LENGTHR

          IF (LENGTHR.GE.LENGTH) THEN
              WRITE(6,*) 'ERROR: LR>LENGTH OF BEARING; TRY AGAIN'
              LENGTHR=LENGTH/2.0D0
              GOTO 100
          END IF

      GOTO 1

C................................................................
C RECESS AXIAL LENGTH
C......................
  105 POINTER=3
      CALL BEEPER
      IF (BEARING.GE.2) THEN   ! FOR SEAL AND HYDRODYNAMIC BEARING
         AR=ZERO
         GOTO 1
      END IF

      WRITE (6, *) 'AR (REAL*8): Axial WIDTH of recess, [m]'
      WRITE (6, 106) AR
  106 FORMAT ('$', 'ENTER Ar [Default (Ar)= ',E14.7E2,']: ')
      CALL ENTERVAL(AR)

      IF (AR.EQ.ZERO) THEN
        WRITE (6,*) 'ERROR: RECESS LENGTH=ZERO[m] REVISE DATA'
        CALL BEEPER
        GOTO 105
      END IF
C##      IF (MODEL.EQ.2) THEN
       IF((LENGTHR.LE.(AR/2.0D0)).OR.(LENGTHL.LE.(AR/2.0D0))) THEN
         WRITE(6,*) 'ERROR: AR/2 > LENGTHL OR LENGTHR ! TRY AGAIN'
         CALL BEEPER
         GOTO 105
       END IF
C##      END IF

      GOTO 1

C..............................................................
C
C NPAD : PAD CIRCUMFERENTIAL ARC LENGTH AND LOCATION
C..............................................................
  115 POINTER=4

C    !.........................!
      IF (BEARING.EQ.2) THEN   !==> FOR ANNULAR SEAL
C    !.........................!
        WRITE (6,*)' FOR ANNULAR SEAL: AUTO DEFAULT IS 1-PAD OF',
     +             ' 360DEG AND START AT 0DEG'
        CALL BEEPER
        NPAD=1
        X1(1)=0.0D0
        LPAD(1)=360.0D0
        SLPAD=360.0D0
        GOTO 223
C    !.........................!
      END IF
C    !.........................!

      WRITE (6, *) 'NPAD (INT*4): Number of PADS on Bearing'
      WRITE (6, 116) NPAD
      CALL BEEPER
  116 FORMAT ('$', 'ENTER NPAD [Default NPAD= ',I2,']: ')
      CALL ENTERINT(NPAD)

      NPAD=MAX0(NPAD,1)


      IF (NPAD.GT.MAXNPAD) THEN
          WRITE(6,*) 'ERROR Number of PADS TOO LARGE ',
     +               'for this version of HYDROPAD'
          WRITE(6,1160) MAXNPAD
 1160     FORMAT (' ','MAX(NPAD)=', I3)
          NPAD=MAXNPAD
          CALL BEEPER
          GOTO 115
      END IF

      LMAX=360.0D0    !max bearing circ.length in degrees

      WRITE (6, *) ' '



C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  117 SLPAD=ZERO

c ==================== !...............................
      DO K=1, NPAD
c ==================== !
C     PAD LEADING EDGE Location and  CIRCUMFERENTIAL LENGTH
C     ..............................
  118 WRITE (6,1118) K
 1118 FORMAT ('$','ENTER/MODIFY DATA for PAD#:',I2,
     + ' LEADING EDGE AND ARC LENGTH [DEGREES]',/)

      CALL BEEPER

  119 WRITE (6, 1119) K,X1(K)
 1119 FORMAT ('$', 'ENTER FOR PAD#', I2,
     +          ' LEADING EDGE ANGLE [Def=',E12.5E2,'(DEG)]:')

      CALL ENTERVAL(X1(K))

      IF ((K.GT.1).AND.(X1(K).LE.X2)) THEN
	 WRITE(6,*)'ERROR: PAD Leading EDGE OVERLAPS Previous PAD'
         X1(K)=X2
         CALL BEEPER
         GOTO 119
      END IF

  120 WRITE (6, 121) K, LPAD(K)
  121 FORMAT ('$', 'ENTER FOR PAD#', I2,
     +          ' ARC LENGTH [Def=',E12.5E2,'(DEG)]:')

      CALL ENTERVAL(LPAD(K))
      LPAD(K)=DABS(LPAD(K))


      IF (LPAD(K).EQ.ZERO) THEN
	 WRITE(6,*)'ERROR: PAD Length is ZERO DEGREES'
         CALL BEEPER
         GOTO 120
      END IF
      WRITE (6, *) ' '

      X2=X1(K)+LPAD(K)                   ! Loc. of Trailing Edge of PAD
      SLPAD=SLPAD+LPAD(K)                ! SUM Pad Lengths

      IF (SLPAD.GT.LMAX*1.01) THEN
         WRITE(6,*) 'ERROR: SUM PAD LENGTHS > 360 degs'
         SLPAD=SLPAD-LPAD(K)
         CALL BEEPER
         GOTO 120
      END IF
c ==================== !
      END DO           ! K=1, NPAD       !.......................!
c ==================== !

      IF (SLPAD.GT.LMAX) THEN
         WRITE(6,*) 'ERROR: SUM PAD LENGTHS > 360.0 DEG'
         CALL BEEPER
         GOTO 117
      END IF


c.................................................................!
c LOBES and PRELOAD coefficients for pad bearings
c.................................................................!

 223  CONTINUE


c    !............................!
      IF (NPAD.GT.1) THEN
c    !............................!
        NLOBES=NPAD
c    !............................!
      ELSE IF (NPAD.EQ.1) THEN
c    !............................!
        IF (LPAD(1).LT.LMAX) THEN

           NLOBES=1

        ELSE

        WRITE (6, 224) NLOBES
 224    FORMAT ('$', 'ENTER Number of LOBES, [Def=',I2,']')
        CALL ENTERINT(NLOBES)
        NLOBES=MAX0(NLOBES,1)

        END IF

        WRITE (6,2241)
 2241   FORMAT ('$',' Nlobes option is in effect only for',/,
     +          '     360deg bearings',/,30('!'))
        CALL BEEPER

c    !............................!
      END IF
c    !............................!

c.................................................................
c    MULTIPLE PADS: Fixed Pads or Lobes, Tilt Pads
c.................................................................

         WRITE (6, 225) PRELOAD
 225     FORMAT ('$', 'ENTER PAD PRELOAD: ',E12.5E2,' (m)]:')
         CALL ENTERVAL(PRELOAD)
         PRELOADB=PRELOAD

         WRITE (6, 227) OFFSET
 227     FORMAT ('$', 'ENTER PAD OFFSET factor: ',E12.5E2,']:')
         CALL ENTERVAL(OFFSET)

c PRELOAD: center of curvature of all pads lies on a preload circle
c          of radius rp=PRELOAD [m] centered on bearing center
c OFFSET: (Angle_pivot - leading_pad_angle)/pad size
c       = 0.50 for pivot at mid-point of pad
c
c.............................................................
c for fixed 1-pad of 360 deg extent then
c.............................................................

         IF ((NPAD.EQ.1).AND.(SLPAD.EQ.LMAX)) THEN
            TILTPAD=0
            ROTPAD(1)=0.0D0 ! zero pad rotation for 360 deg bearing
            GOTO 229
         END IF

c.............................................................
c      FIXED or TILT / FLEX  PADS
c.............................................................

         IF (TILTPAD.EQ.1) THEN
            YN='Y'
         ELSE
            YN='N'
         END IF

         WRITE (6, 228) YN
 228     FORMAT ('$',' ARE PADS ABLE TO TILT ? Y or N ?',
     +               ' Def:[',A1,']')
         READ (5, 265) YN

         IF ((YN.EQ.' ').AND.(TILTPAD.EQ.1)) YN='Y'
         IF ((YN.EQ.' ').AND.(TILTPAD.EQ.0)) YN='N'
         IF ((YN.EQ.'Y').OR.(YN.EQ.'y')) THEN
              TILTPAD=1
         ELSE
              TILTPAD=0
         END IF


c       !...........................! .................!
 229     IF (TILTPAD.LE.0) THEN
c       !...........................! FIXED PAD BEARING
            TILTPAD=0
            ROTDEL=0.0D0
            DO K=1, NPAD
               ROTPAD(K)=ROTDEL     ! zero pad rotations
            END DO
c       !...........................!
          ELSE IF (TILTPAD.GT.0) THEN
c       !...........................! READ Pad rotations
             TILTPAD=1

             Ndumy=1                !# SET Ndumy=NPAD for values of
c                                   !# pad Inertia and rotational stiffness/
c                                   !# damping different for each pad

c            .............................!
             DO K=1, NPAD                 ! NPAD
c            .............................!
             WRITE (6,230) K, ROTPAD(K)
 230         FORMAT ('$','ENTER PAD ', I2,
     +                  ' rotation angle in RADIANS, Def[',
     +                   E12.5E2,']')
             ROTDEL=ROTPAD(K)
             CALL ENTERVAL(ROTDEL)
             ROTPAD(K)=ROTDEL
c            .............................!
             END DO
c            .............................!
             CALL BEEPER

c            .............................!
             DO K=1, Ndumy                ! PAD Mass Moment of INERTIA
c            .............................!
             WRITE (6,231) K, INERPAD(K)
 231         FORMAT ('$','ENTER PAD ', I2,
     +                  ' MASS Moment of INERTIA, Def[',
     +                   E12.5E2,'(kg.m2)]')
             IPAD=INERPAD(K)
             CALL ENTERVAL(IPAD)
             INERPAD(K)=DABS(IPAD)
c            .............................!
             END DO
c            .............................!

c            .............................!
             DO K=1, Ndumy                ! PAD Rotational Stiffness
c            .............................!
             WRITE (6,232) K, KSTPAD(K)
 232         FORMAT ('$','ENTER PAD ', I2,
     +                  ' Rotational Stiffness, Def[',
     +                   E12.5E2,'(N.m/rad)]')
             KROTPAD=KSTPAD(K)
             CALL ENTERVAL(KROTPAD)
             KSTPAD(K)=DABS(KROTPAD)
c            .............................!
             END DO
c            .............................!

c            .............................!
             DO K=1, Ndumy                ! PAD Rotational Damping
c            .............................!
             WRITE (6,233) K, CDAPAD(K)
 233         FORMAT ('$','ENTER PAD ', I2,
     +                  ' Rotational Damping, Def[',
     +                   E12.5E2,'(N.m.s/rad)]')
             CROTPAD=CDAPAD(K)
             CALL ENTERVAL(CROTPAD)
             CDAPAD(K)=DABS(CROTPAD)
c            .............................!
             END DO
c            .............................!

c            COPY values to other pads ......
             IF (Ndumy.EQ.1) THEN
              DO K=2, NPAD
               INERPAD(K)=IPAD
               KSTPAD(K)=KROTPAD
               CDAPAD(K)=CROTPAD
              END DO
             END IF
c           !................................!


c       !...........................! TILTPAD=1
          END IF
c       !...........................! TILT-PAD BEARING


c PRELOAD: center of curvature of all pads lies on a preload circle
c          of radius rp=PRELOAD [m] centered on bearing center
c OFFSET: (Angle_pivot - leading_pad_angle)/pad size
c       = 0.50 for pivot at mid-point of pad
c ROTDEL: angle of rotation of pad about pivot in RADIANS

      GOTO 1

c .................................................................
c ................... PAD :  RECESS DATA  .........................
c .................................................................

 1205 POINTER=5

      IF (BEARING.GE.2) THEN    ! Annular seal and plain bearing
        DO K=1, NPAD            ! NO recesses on pads.
          NREC(K)=0
        END DO
        GOTO 1
      END IF

c ==================== !
      DO K=1, NPAD
c ==================== !

      WRITE (6,1210) K
 1210 FORMAT ('$','ENTER/MODIFY RECESS DATA for PAD:',I2)
      CALL BEEPER
C     ..............................................................
C      NPOCKET = Number of pockets on K PAD
C     .................
  122 WRITE (6, 123) K, NREC(K)
      NPOCKET=NREC(K)
  123 FORMAT ('$', 'ENTER Number of RECESSES for PAD ',I2,
     +            ' [DEF Nrec=',I2,']:')

      CALL ENTERINT(NPOCKET)
      NREC(K)=NPOCKET

      IF (NPOCKET.GT.MAXNPOCK) THEN
          WRITE(6,*) 'ERROR Number of recesses TOO LARGE for this',
     +               'version of HYDROFLEX'
          NREC(K)=MAXNPOCK
          CALL BEEPER
          GOTO 122
      END IF

      IF (NPOCKET.EQ.0) GOTO 177


C....................................................................
C     POCKET LEADING EDGE Location relative to PAD leading edge
C                    and  CIRCUMFERENTIAL LENGTH
C..............................
      SLREC=ZERO
      X2r=LPAD(K)  ! length of K-PAD in degrees

      WRITE (6, 1235) K,LPAD(K)
 1235 FORMAT('$ FOR PAD # ',I2,' of LENGTH=', E12.5E2,'[DEG]',/,
     +     1X,' ENTER/MODIFY: RECESS Edge POSITION [DEG]',
     +       ' RELATIVE to PAD and RECESS Arc Length [DEG]',/)

c    !---------------!
      DO I=1, NPOCKET
c    !---------------!
  124 WRITE (6, 125) I,K,X1r(I,K)
  125 FORMAT ('$', 'ENTER FOR Recess', I2,' OF PAD#', I2,
     +          ' EDGE POSITION [Def= ',E12.5E2,'(DEG)]:')

      CALL ENTERVAL(X1R(I,K))
      X1r(I,K)=DABS(X1R(I,K))

      IF (I.eq.1) THEN
          X1rmi=ZERO   ! X1(K)
      ELSE
          X1rmi=X1r(I-1,K)+LREC(I-1,K)
      END IF

      IF ((X1r(I,K).LE.X1rmi).OR.(X1r(I,K).GE.X2r))  THEN
        write (6,*) 'WARNING: POCKET OUT OF PAD/OVERLAPS PREVIOUS',
     +              ' RECESS'
        X1r(I,K)=X1rmi
        CALL BEEPER
        GOTO 124
      END IF

  126 WRITE (6, 127) K,I,LREC(I,K)
  127 FORMAT ('$', 'ENTER FOR PAD#', I2,' RECESS:', I2,
     +            '  ARC Length [Def= ',E12.5E2,'(DEG)]:')
      CALL ENTERVAL(LREC(I,K))
      LREC(I,K)=DABS(LREC(I,K))

      WRITE (6, *) ' '

      SLREC=SLREC+LREC(I,K)

      IF (X1r(I,K)+LREC(I,K).GE.X2r)  THEN
        write (6,*) 'WARNING: RECESS EDGE at END/OVERLAPS PAD'
        SLREC=SLREC-LREC(I,K)
        CALL BEEPER
        GOTO 126
      END IF

      IF (LREC(I,K).GT.LPAD(K)) THEN
        WRITE (6,*) 'ERROR: RECESS ARC LENGTH GREATER THAN PAD ',
     +              'LENGTH=' , LPAD(K), '(DEG)'
        SLREC=SLREC-LREC(I,K)
        CALL BEEPER
        GOTO 126
      END IF

      IF (SLREC.GT.LPAD(K)) THEN
         WRITE(6,*) 'ERROR: SUM(REC LENGTHS) GREATER THAN PAD ',
     +              'LENGTH=' , LPAD(K), '(DEG)'
        SLREC=SLREC-LREC(I,K)
        CALL BEEPER
        GOTO 126
      END IF

  128 WRITE (6, 129) K,I, PADORIF(I,K)
  129 FORMAT ('$', 'ENTER FOR PAD#', I2,' RECESS:', I2,
     +            ' ORIFICE Diameter [Def= ',E12.5E2,'(m)]:')

      CALL ENTERVAL(PADORIF(I,K))
      PADORIF(I,K)=DABS(PADORIF(I,K))

      WRITE (6, *) ' '

c    !---------------!
      END DO         ! I=1,2, NPOCKET on Pad K
c    !---------------!

c ==================== !

 177  END DO           ! K=1, 2,  NPAD
c ==================== !
      GOTO 1


C......................................................................
C GRID POINTS IN AXIAL DIRECTION:
C................................
  111 POINTER=7
      CALL BEEPER

C   !..............................!...................!
      IF (BEARING.EQ.1) THEN
C   !..............................! RECESSED BEARING
      WRITE (6, *) 'NLA (INT*4): Number of grid points, axial on land (w
     +/o recess).'
      WRITE (6, 112) NLA
  112 FORMAT ('$', 'ENTER Nla [Default Nla=',I2,']: ')
      CALL ENTERINT(NLA)
      WRITE (6, *) 'NPA (INT*4): Number of grid points, axial on 1/2 rec
     +ess side.'
      WRITE (6, 113) NPA
  113 FORMAT ('$', 'ENTER Npa [Default Npa=',I2,']: ')
      CALL ENTERINT(NPA)

      IF ((NLA+NPA).GT.MAXNYI) THEN
          WRITE(6,*) 'ERROR Number of AXIAL GRID POINTS TOO LARGE ',
     +               'for this version of HYDROFLEX'
          WRITE(6,114) MAXNYI
  114     FORMAT (' ','MAX(NLA+NPA)=', I3)
          NLA=MAXNYI/2
          NPA=MAXNYI/2
          GOTO 111
      END IF

C   !..............................!...................!
      ELSE
C   !..............................! SEAL or PLAIN BEARING
      WRITE (6, *) 'NYI (INT*4): Number of AXIAL grid points'
      CALL BEEPER
      WRITE (6, 1121) NYI
1121  FORMAT ('$', 'ENTER Nyi [Default Nyi=',I2,']: ')
      CALL ENTERINT(NYI)

      IF (NYI.GT.MAXNYI) THEN
          WRITE(6,*) 'ERROR Number of AXIAL GRID POINTS TOO LARGE ',
     +               'for this version of HYDROFLEX'
          WRITE(6,1141) MAXNYI
1141      FORMAT (' ','MAX(NYI)=', I3)
          NYI=MAXNYI
          GOTO 111
      END IF

      NPA=1
      NLA=NYI
C   !..............................!...................!
      END IF
C   !..............................! SEAL or PLAIN BEARING

      WRITE (6, *) ' '
      GOTO 1


C.....................................................................
C NUMBER OF GRID POINTS CIRCUMFERENTIAL
C.....................................................................

  140 CALL BEEPER
      POINTER=8
C    !..........................!..................................!
      IF (BEARING.EQ.1) THEN    ! FOR RECESSED HYDROSTATIC BEARING
C    !..........................!..................................!
      WRITE (6, *) 'NLC (INT*4): Number of grid points, circumferential
     +on land.'
      WRITE (6, 141) NLC
  141 FORMAT ('$', 'ENTER Nlc [Default Nlc= ',I2,']: ')
      CALL ENTERINT(NLC)

      IF (NPC.LT.3) NPC=3

  142 WRITE (6, *) 'NPC (INT*4): Number of grid points, circumferential
     +on recess side; must be ODD'
      WRITE (6, 143) NPC
  143 FORMAT ('$', 'ENTER Npc [Default Npc=',I2,']: ')
      CALL ENTERINT(NPC)

      IF (MOD(NPC, 2).EQ.1) goto 144  !##?? JMOD function is not in CRAY
          WRITE (6, *) 'Npc must be ODD.  Please enter a new value.'
          NPC=NPC-1
          GOTO 142

C NOTE:
C  JMOD (A,2) = GIVES REMAINDER OF DIVISION BETWEEN INTEGER*4 A AND 2
C               IF = 1 THEN A IS ODD
C               IF = 0 THEN A IS EVEN
C  IN OTHER SYSTEMS USE IMOD OR MOD FUNCTION

  144 CONTINUE
      IF ((NPC+NLC-2)*NPOCKET+NLC.GT.MAXNXT) THEN
          write(6,*) 'ERROR Number of CIRC. GRID POINTS TOO LARGE ',
     +               'for this version of HYDROFLEX'
          write (6,145) MAXNXT
  145 FORMAT (' ','MAX({NPC+NLC-2}*NPOCKET+NLC)=', I3)
          GOTO 140
      END IF


C    !..........................!..................................!
      ELSE     ! FOR ANNULAR SEAL OR PLAIN BEARING
C    !..........................!..................................!
      WRITE (6, *) 'NXT (INT*4): Number of circumferential grid points'
      WRITE (6, 1411) NXT
 1411 FORMAT ('$', 'ENTER Nxt [Default Nxt= ',I2,']: ')
      CALL ENTERINT(NXT)

      IF (NXT.GT.MAXNXT) THEN
          write(6,*) 'ERROR Number of CIRC. GRID POINTS TOO LARGE ',
     +               'for this version of HYDROFLEX'
          write (6,1451) MAXNXT
          GOTO 140
      END IF

      NLC=NXT
      NPC=0

 1451 FORMAT (' ','MAX(NXT)=', I3)

C    !..........................!..................................!
      END IF
C    !..........................!..................................!


      WRITE (6, *) ' '
      GOTO 1

C
C ================================================ BUILDS MESH
C                                 AND SAVES IT ON PLOTMESH FILE
 999  POINTER=9
      LD=LENGTH/DIAM
      R=DIAM/2.0D0
      LMAX=360.0D0
      KKK=0
C...................................

      CALL XYDATA(KKK,RPM,PS,PA)
C................................... BUILDS MESH
C..........................................................
C
C READ MIXING COEFFICIENT FOR PAD BEARING
C..........................................................
      IF (BEARING.EQ.2) RETURN

      IF (IFULL.EQ.0) THEN
C    !.....................! FOR PAD BEARING
      CALL BEEPER
      WRITE (6, *) '==========================================='
      WRITE (6, *) '$ ENTER MIXING COEFFICIENT AT PAD GROOVES  '
      WRITE (6, *) '$ LEMDA=0, NO MIXING OF FLUID AT GROOVE'
      WRITE (6, *) '$ LEMDA=1=MAX, ALL HOT FLUID ENTERS PAD'
      WRITE (6, *) '$ AT LEADING EDGE'
      WRITE (6, 1981) LEMDA
 1981 FORMAT ('$', 'ENTER MIXING COEFFICIENT= [',F7.4,']: ')
      CALL ENTERVAL(LEMDA)
      LEMDA=DABS(LEMDA)
      LEMDA=DMIN1( 1.0D0, LEMDA)
      WRITE (6, *) ' '
C    !...................!
      END IF
C    !...................!


 1000 FORMAT (40A)
 1010 FORMAT (' ', 40A)

      RETURN
C
C..............................................................
 601  FORMAT (' ', 75('-'),/,3X,
     +       'SELECT BEARING TYPE [Def:',I2,']',/,3X,
     +       '(1) HYDROSTATIC PAD BEARING, (2) ANNULAR SEAL',/,3X,
     +       '(3) CIRCULAR JOURNAL PAD BEARING',  /,3X,
     +        75('.'))

7011  FORMAT (' ', 75('-'),/,3X,'GEOMETRICAL CHARACTERISTICS OF:'
     +       ,' TILT PAD ', 30A,/,1X,75('-'),/,3X,
     +       'SELECT PARAMETERS TO CHANGE FROM LIST BELOW',/)
7012  FORMAT (' ', 75('-'),/,3X,'GEOMETRICAL CHARACTERISTICS OF:'
     +       ,' FIXED PAD ', 30A,/,1X,75('-'),/,3X,
     +       'SELECT PARAMETERS TO CHANGE FROM LIST BELOW',/)

7021  FORMAT (1X,'(1) DIAMETER =   ',E12.5E2,'[m]',/,
     +        1X,'(2) AXIAL LENGTH=',E12.5E2,'[m]=',
     +        E12.5E2,'(Lr)+',E12.5E2,'(Ll) [m]')
7022  FORMAT (1X,'(3) RECESS AXIAL L =',E12.5E2,'[m]')
7023  FORMAT (1X,'(4) PADS Number =',I3,1X,'Total X-length=',
     +        E12.5E2,'[DEG]',/,5X,
     + '#  Position   Arc Length',1X,'Tilt angle', 2X,
     + 'Inertia',4X,'Krot',6X,'Crot',/,5X,
     + '-    [DEG]      [DEG]   ',1X,'[radians] ', 3X,
     + '(kg.m2)',3X,'N-m/rad',3X,'N-m.s/rad')

7024  FORMAT (1X,'(4) PADS Number =',I3,1X,'Total X-length=',
     +        E12.5E2,'[DEG]',/,5X,
     +        '#  Position      Arc Length',5X,'FIXED-RIGID PAD',
     + /,5X,  '-    [DEG]        [DEG]    ',5X,'---------------')


7031  FORMAT (4X,I2,6(1X,E10.3E2) )
7032  FORMAT (4X,I2,2(1X,E12.5E2) )

7033  FORMAT (4X,'NPADs =',I3,2X,'PRELOAD:',E12.5E2,'(m)',
     +        2X,'OFFSET:',F5.3)
7034  FORMAT (4X,'NLOBES=',I3,2X,'PRELOAD:',E12.5E2,'(m)',
     +        2X,'OFFSET:',F5.3)
7035  FORMAT (4X,'PAD Inertia Moments (kg.m2):',
     +        5(E10.4E2,1X))

 704  FORMAT (1X,'(5) PAD:',I2,
     +        4X,'RECESS  Position       Length         Orifice',/,
     +       15X,'        Angle[DEG]     Angle[DEG]     do (m)' )

7040  FORMAT (1X,'(5) PAD:',I2,' NO RECESS ')

 705  FORMAT (16X, I2,1X, 3(2X,E12.5E2))
 706  FORMAT (1X,'(6) CLEARANCE FUNCTION and RECESS DEPTH')
7081  FORMAT (1X,'(7) GRID Points AXIAL NLA=',I3,' NPA=',I3,/
     +        1X,'(8) GRID Points CIRC. NLC=',I3,' NPC=',I3)
7082  FORMAT (1X,'(7) GRID Points AXIAL NYI=',I3,/
     +        1X,'(8) GRID Points CIRC. NXT=',I3)
 709  FORMAT (1X,'(9) RETURN',/,' ',75('.'))


      END

C
C *****************************************************************************
C **                                                                         **
C **  Subroutine Xydata                                                      **
C **                                                                         **
C **  XYDATA: Builds vector arrays for P, U, and V control volumes.          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE XYDATA(KKK,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /RECASP/ ASPEC(MAXNPOCK)
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
c     ................................................................
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

c     ................................................................
      DOUBLE PRECISION XP, XU, YP, YV, DXP, DXU, SUW, SUE,
     +                 DYP, DYV, SVN, SVS,
     +                 COSXP, SINXP, COSXU, SINXU, ASPEC
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD, AR, HREC,ROTPAD,
     +                 X1, LPAD, X1r,Lrec, LENGTHR, LENGTHL,
     +                 PaDorif, DIAORIF, CORIF, FREQU,SIGMA,
     +                 L1, RES, PRELOAD, OFFSET,ROTDEL,INERPAD,
     +                 KSTPAD, CDAPAD, IPAD, RSINPK, RCOSPK,
     +                 KROTPAD, CROTPAD

      INTEGER NPAD, NREC, NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI,
     +        NYI, NXT, ISYM, ICSTEP, IFULL, ICASE, NCASE, BEARING,
     +        TILTPAD, TILT

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION ZERO, R,LMAX,LPmax, Lrmax,Xrmax,HGROV, HPOCK,
     +       DYPU, DYLU, DXL(MAXNPOCKP1),DXR(MAXNPOCKP1),RPM,PS,PA,
     +                 X2r, XPI, XUI,DYLUR,DYLUL, CONVER, Xpivot
      INTEGER I,J,II,JJ,N,K,COUNTER,KP,KK,IOS,idumy,jdumy,Istart,KKK,
     +        KSTART,KSTEP,ISET
C ----------------------------------------------------------------------------
C --  XYDATA code                                                           --
C ----------------------------------------------------------------------------
C Note: dimensionless coordinates divided by RADIUS = DIAM/2.0D0

      ZERO=0.0D0
      LMAX=360.0D0
      CONVER=DACOS(-1.0D0)/180.0D0       ! DEG to RADIANS.
      R=DIAM/2.0D0
      LD=LENGTH/DIAM
C ----------------------------------------------------------------------------
C Generate Mesh Vectors in Axial Direction
C ----------------------------------------------------------------------------
      TILT=TILTPAD

      IF (BEARING.GE.2) GOTO 120       ! => SEALS AND PLAIN BEARINGS
c     ..................................................!


C -------------------------
C FOR HYDROSTATIC BEARINGS
C -------------------------

  100 NYI=NLA+NPA                           ! Number of P-CV axially
      NPAP1=NPA+1

      DYPU=(AR/R)/(2*NPA-1)                 ! Y size on pocket side
      DYLUR=(LengthR-AR/2.0D0-R*DYPU/2.0D0) ! Y size on right(top) land
      DYLUR=2.0D0*DYLUR/(2.0D0*NLA-1)/R

      DYLUL=(LengthL-AR/2.0D0-R*DYPU/2.0D0) ! Y size on left(bottom) land
      DYLUL=2.0D0*DYLUL/(2.0D0*NLA-1)/R
      DYLU=DMIN1(DYLUL,DYLUR)

      IF ((DYPU.LE.ZERO).OR.(DYLU.LE.ZERO))THEN
        WRITE (6,615) DYPU,DYLU
        CALL BEEPER
        RETURN
      END IF

      KSTART=1
      KSTEP=1
      DYLU=DYLUR
      ISET=1
C ............................,,,,,......... ! Generate DY vectors
  110 DYV(KSTART)=DYPU/2.0D0                 !
      DYP(KSTART)=DYPU                       !
      DO I=KSTART+KSTEP,KSTEP*NPA,KSTEP      ! Y size P and U CVS
          DYP(I)=DYPU                        ! DYP(I) = DYU(I)
          DYV(I)=DYPU                        ! Y size V CV
      END DO                                 !
      DYV(KSTEP*NPAP1)=(DYPU+DYLU)/2.0D0     ! at pocket boundary
      DYP(KSTEP*NPAP1)=DYLU                  ! irregular CV
      DO I=KSTEP*(NPA+2),KSTEP*NYI,KSTEP     !
          DYP(I)=DYLU                   ! on extended lands
          DYV(I)=DYLU                   ! to bearing
      END DO                            ! outlet
      DYP(KSTEP*NYI)=DYLU/2.0D0         !

      YV(KSTEP)=0.0D0                             ! from bottom (y=0)
      YP(KSTEP)=KSTEP*DYPU/2.0D0                  ! to outlet (y=LR/R)
      SVN(KSTEP)=0.0D0                            !
      SVS(KSTEP)=0.0D0                            !
      DO I=2*KSTEP, KSTEP*NYI,KSTEP               ! where
          II=I-KSTEP                              !  YP(I) = YU(I)
          YV(I)=YV(II)+Kstep*DYP(II)              !
          YP(I)=YP(II)+Kstep*DYV(I)               !
          SVN(I)=DYP(I)/DYV(I)/2.0D0              ! North weight coef.
          SVS(I)=DYP(II)/DYV(I)/2.0D0             ! South weight coef.
      END DO                                !
      SVS(KSTEP)=0.5D00                     ! c###
      SVN(KSTEP)=0.5D00
      SVN(KSTEP*NYI)=2.0D0*SVN(KSTEP*NYI)

      GOTO 160

C-------------------------------------
C FOR ANNULAR SEALS AND PLAIN BEARINGS
C-------------------------------------

 120  DYLUR=LENGTHR/R/(NYI-1)
      DYLUL=LENGTHL/R/(NYI-1)
      KSTEP=1
      DYLU=DYLUR
      ISET=1

 130  DYP(KSTEP)=DYLU/2.0D00

      IF (DYLU.LE.ZERO)THEN
        WRITE (6,616) DYLU
        CALL BEEPER
        RETURN
      END IF

      DYV(KSTEP)=0.0D0
      YP(KSTEP) =0.0D0
      YV(KSTEP) =0.0D0
      SVN(KSTEP)=0.5D00
      SVS(KSTEP)=0.5D00

      DO I=KSTEP+KSTEP, KSTEP*NYI, KSTEP
         DYP(I)=DYLU
         DYV(I)=DYLU
         II=I-KSTEP
         YV(I)=YV(II)+KSTEP*DYP(II)
         YP(I)=YP(II)+KSTEP*DYV(I)
         SVN(I)=DYP(I)/DYV(I)/2.0D0              ! North weight coef.
         SVS(I)=DYP(II)/DYV(I)/2.0D0             ! South weight coef.
      END DO

      SVS(2*KSTEP)=2.0D0*SVS(2*KSTEP)
      DYP(KSTEP*NYI)=DYLU/2.0D0
C     SVN(KSTEP*NYI)=SVN(KSTEP*NYI)

      IF (BEARING.EQ.2) GOTO 200   ! => FOR ANNULAR SEAL


C .........................................................................C

 160  IF (ISET.EQ.2) GOTO 200
c   !.........................!
      IF (ISYM.EQ.0) THEN
c   !.........................!
         KSTART=-1
         KSTEP=-1
         DYLU=DYLUL
         ISET=2
         IF (BEARING.EQ.1) GOTO 110
         IF (BEARING.EQ.3) GOTO 130
c   !.........................!
      ELSE
c   !.........................!
        DO I= 1, NYI
         YP(-I)=-YP(I)
         YV(-I)=-YV(I)
         DYP(-I)=DYP(I)
         DYV(-I)=DYV(I)
         SVN(-I)=SVN(I)
         SVS(-I)=SVS(I)
        END DO
c   !.........................!
      END IF
c   !.........................!


 200     YP(0)=0.0D0
         YV(0)=0.0D0
         DYP(0)=0.0D0
         DYV(0)=0.0D0
         SVN(0)=0.0D0
         SVS(0)=0.0D0


C ===============================================================
C Construct MESH Vectors in X-Circumferential Direction
C
C  !======================!
      IF (KKK.EQ.0) THEN
C  !======================!
C     ..................................! CHECK FOR PADS & RECESS ARC LENGTHS
C     NOTE: here comparisons are performed in DEG
C
      LMAX=360.0D0      !degrees
      LPmax=0.0D0
      Xrmax=ZERO
c    !......................!
      DO K=1, NPAD
c    !......................! SWEEP on all PADS.
        LPmax=LPmax+LPAD(K)                           ! Total length of bearing/pads

        NPOCKET=NREC(K)
        IF (NPOCKET.eq.0) goto 12

        Lrmax=0.0D0
        DO I=1, NPOCKET
           Lrmax=Lrmax+Lrec(I,K)                      ! Total length of recesses on pad
        END DO

        Xrmax=Lrec(NPOCKET,K)+X1r(NPOCKET,K)+X1(K)    ! Xrmax:  ABSOLUTE Location of
c                                                     ! RECESS trailing edge
        IF (Lrmax.GT.LPAD(K)*1.01) THEN
           WRITE(6,611) K
           CALL BEEPER
           RETURN
        END IF

  12  CONTINUE
c    !......................!
      END DO
c    !......................!

      IF (LPmax.GT.361.0D0) THEN        !MAX ARC PAD LENGTHS = 360.0 deg
          WRITE (6,613) LPmax
          CALL BEEPER
          RETURN
      END IF
      Xrmax=Xrmax-X1(1)

      IFULL=0    ! IFULL=1: circular 360 DEG HJB,ONE PAD with NPOCKET

c     ..................................................!
      IF (BEARING.EQ.1) THEN       ! => HJB
c     ..................................................!
      IF ( (NPAD.EQ.1).AND.
     +     (DABS(LPmax/LMAX-1.0D0).LE.0.001D0).AND.
     +     (DABS(Xrmax/LMAX-1.0D0).LE.0.001D0) ) IFULL=1

         Jdumy=2*NYI+1+2
c     ..................................................!
      ELSE IF (BEARING.EQ.2) THEN  ! => SEAL
c     ..................................................!
         Jdumy=NYI+1
         IFULL=1
c     ..................................................!
      ELSE IF (BEARING.EQ.3) THEN  ! => BEARING PAD
c     ..................................................!
      IF ( (NPAD.EQ.1).AND.
     +     (DABS(LPmax/LMAX-1.0D0).LE.0.001D0)) IFULL=1
         Jdumy=2*NYI
c     ..................................................!
      END IF
c     ..................................................!

      Idumy=0
      KP=0
c    !...................!
 400     KP=KP+1         ! SWEEP ON PADS -> Find # of X-Grid points
c    !...................!

      NPOCKET=NREC(KP)
      IF (NPOCKET.EQ.0) GOTO 410

      Idumy=Idumy+(NLC+NPC)*NPOCKET

      IF (IFULL.EQ.0) Idumy=Idumy+NLC

      GOTO 412

  410 Idumy=Idumy+NLC   ! NO POCKET ON PAD:

  412 IF (KP.EQ.NPAD) GOTO 477

      IF ((KP.GE.1).AND.(KP.LT.NPAD)) Idumy=Idumy+2
c    !...............!
       GOTO 400      ! -------> NEXT PAD KP=1,2,3,... NPAD
c    !...............!

 477  WRITE (6,*)
  !23  FORMAT (' ','  *transfer bearing/mesh data to PLOTMESH File')

      OPEN (UNIT=2, FILE='PLOTMESH', STATUS='UNKNOWN', IOSTAT=IOS,
     +       ERR=1000)
      WRITE ( 2, *, IOSTAT=IOS, ERR=1200) Idumy, Jdumy

      K=1

C  !======================!
      ELSE
C  !======================!
      K=KKK               ! Pad Number
      NPOCKET=NREC(KKK)
      IF (NPOCKET.GT.0) THEN
       DO I=1, NPOCKET
        Diaorif(I)=PaDorif(I,KKK)     ! transfer Orifice diameter data
       END DO
      END IF
      ROTDEL=ROTPAD(KKK)              ! pad rotation angle in [rads]

      CALL CPARAM(RPM,PS,PA)		      !==> initopst.f (INITIAL CALCS.)

c   !.....................! for 360 deg HJB check if all rec+land are equal
      IF ((IFULL.EQ.1).AND.(NPOCKET.GT.0)) THEN
        ICASE=1           ! for perturbations in X&Y dir.
        Lrmax=Lrec(1,K)+X1r(1,K)
        DO I=2, NPOCKET
            Xrmax=Lrec(I,K)+X1r(I,K)-(I-1)*Lrmax !length of recesses + land
            IF ( DABS(Xrmax/Lrmax-1.0D0).GE.0.001D0) ICASE=2
        END DO
c   !.....................!
      ELSE
c   !.....................!
        ICASE=2
        IF ((IFULL.EQ.1).AND.(BEARING.GE.2)) ICASE=1
c   !.....................!
      END IF
c   !.....................!

C  !======================!
      END IF
C  !======================!
C

C ====================  MESH GENERATION =========================== C
C NOTE:  DX and X vectors and arrays are in radians.
C
C ..................................... ! Generate DX vextors
  300 CONTINUE

      NXI=NLC+NPC-2                     ! # of PCV on recess+land

      NPOCKET=NREC(K)

C    !..................................!
      IF (NPOCKET.EQ.0) THEN
C    !..................................!
           NXT=NLC
           DXL(1)=LPAD(K)*Conver/(NXT-1)     !#d# size
           GOTO 310
C    !..................................!
      ELSE
C    !..................................!
           NXT=NXI*NPOCKET+NLC
           IF (IFULL.EQ.1) NXT=NXT-NLC+1
C    !..................................!
      END IF
C    !..................................!



      DXL(1)=(X1r(1,K))*Conver/(NLC-1)       !#d# DX size on left of FIRST recess
      X2r=(X1r(1,K)+Lrec(1,K))               ! ABSOLUTE Location of downstream
c                                            ! recess edge

      IF (DXL(1).LE.ZERO) THEN
          print *, 'DXL(1)=', DXL(1), 'PAD#:', K
          GOTO 997
      END IF

      IF (NPOCKET.EQ.1) GOTO 302

c    !..................................!
       DO i=2, NPOCKET
c    !..................................!
         DXL(i)=(X1r(i,K)-X2r)*Conver/(NLC-1) !#d# DX size between recesses
         X2r=(X1r(i,K)+Lrec(i,K))             !#d# LOC of REC downstream edge

         IF (DXL(i).LE.ZERO) THEN
          print *, 'DXL(i)=', DXL(i), 'i:', i, 'PAD#:', K
          GOTO 997
         END IF
c    !..................................!
       END DO
c    !..................................!


c    !..................................!
 302  i=NPOCKET
c    !..................................!


c    !......................! #d# DX size between last recess edge and boundary
      IF (IFULL.EQ.1) THEN
c    !......................! 360 HJB
         DXL(i+1)=DXL(1)
c    !......................!
      ELSE
C    !......................! PAD HJB
         DXL(i+1)=(LPAD(K)-X2r)*Conver/(NLC-1)  !
c    !......................!
      END IF
c    !......................!


      IF (DXL(i+1).LE.ZERO) THEN
         print *, 'DXL(i)=', DXL(i+1), 'i+1', i+1, 'PAD#:', K
         GOTO 997
      END IF


c    !..................................!#d# DX size on recess X-boundary
 304  DO i=1, NPOCKET
c    !..................................!
c#d#   DXR(i)=(2.0*Lrec(i,K)/R-DXL(i)-DXL(i+1))/2.0/(Npc-2)
       DXR(i)=(2.0*Lrec(i,K)*Conver-DXL(i)-DXL(i+1))/2.0/(Npc-2)

       IF (DXR(i).LE.ZERO) THEN
          print *, 'DXR(i)=', DXR(i), 'i', i, 'PAD#:', K
          GOTO 997
       END IF
c    !..................................!
      END DO
c    !..................................!


c    !..................................!
  306 DO J=1, NPOCKET                   !
c    !..................................!
          JJ=(J-1)*NXI                  ! Control volumes
          DO I=1, NLC                   ! on inter-recess lands
              II=JJ+I                   ! between pockets
              DXP(II  )=DXL(J)          ! X size P & V control vols.
              DXU(II  )=DXL(J)          ! DXP(I) = DXV(I)
          END DO                        ! I=1, NLC

          DXU(JJ+NLC)=(DXL(J)+DXR(J))/2.0D0 ! on pocket boundary

          DO I=NLC+1, NXI               !
              II=JJ+I                   ! Control volumes for
              DXP(II  )=DXR(J)          ! lands on top of
              DXU(II  )=DXR(J)          ! pockets
          END DO                        ! I=NLC+1, NXI
          DXU(JJ+NXI)=(DXR(J)+DXL(J+1))/2.0D0
        ASPEC(J)=Lrec(J,K)*Conver/2.0D0 ! 1/2 recess length aspect ratio
c    !..................................!
      END DO                            ! J=1, NPOCKET
c    !..................................!



c    !..................................!
      IF (IFULL.EQ.1) THEN              ! 360 deg HJB 1PAD
c    !..................................!
          DXP(NXT)=DXP(1)
          DXU(NXT)=DXU(1)
c    !..................................!
      ELSE
c    !..................................! PAD BEARING
          JJ=NPOCKET*NXI
          J=NPOCKET+1
          DO I=1, NLC
            II=JJ+I
            DXP(II)=DXL(J)
            DXU(II)=DXL(J)
          END DO
c#?       NXT=II
c    !..................................!
      END IF
c    !..................................!

      GOTO 312

c    !..................................! NO NPOCKETS on PAD.
  310 DO II=1, NXT
        DXP(II)=DXL(1)
        DXU(II)=DXL(1)
      END DO
c    !..................................!

C     ................................. ! Generate X vectors
  312 XP(1)=X1(K)*Conver
      XU(1)=XP(1)+DXP(1)/2.0D0
      DO I=2, NXT
         II=I-1
         XP(I)=XP(II)+DXU(II)
         XU(I)=XU(II)+DXP(I)
      END DO

C    !............................!TRIG Functions for Film Thickness in Circular Journal                                                           --
        DO I=1, NXT-1
          XPI=XP(I)
          XUI=XU(I)
          COSXP(I)=DCOS(XPI)
          COSXU(I)=DCOS(XUI)
          SINXP(I)=DSIN(XPI)
          SINXU(I)=DSIN(XUI)
          SUW(I)=DXP(I)/DXU(I)/2.0D0           ! West weight coef.
          SUE(I)=DXP(I+1)/DXU(I)/2.0D0         ! East weight coef.
        END DO                                 !

c    !............................!............!..................!
        IF (IFULL.EQ.1) THEN                   !1PAD of 360 Length
c    !............................!
         COSXP(NXT)=COSXP(1)                   ! for periodicity
         COSXU(NXT)=COSXU(1)                   ! condition
         SINXP(NXT)=SINXP(1)                   !....................
         SINXU(NXT)=SINXU(1)
         SUW(NXT)=SUW(1)
         SUE(NXT)=SUE(1)
c    !............................!
        ELSE
c    !............................!
          XPI=XP(NXT)
          XUI=XU(NXT)
          COSXP(NXT)=DCOS(XPI)
          COSXU(NXT)=DCOS(XUI)
          SINXP(NXT)=DSIN(XPI)
          SINXU(NXT)=DSIN(XUI)
          SUW(NXT)=DXP(NXT)/DXU(NXT)/2.0D0
          SUE(NXT)=1.0D0-SUW(NXT)
          SUW(1)=1.0D0-SUE(1)
c    !............................!
        END IF
c    !............................!


C !------------------------------! KKK=KPAD (current pad)
      IF (KKK.GT.0) THEN
C !------------------------------!
c      Allocate values of Tilt-Pad parameters
       IF (TILTPAD.EQ.1) THEN                  ! Calculate location
         Xpivot=OFFSET*(XP(NXT)-XP(1))+XP(1)   ! of K-pad pivot in rads.
         RSINPK=R*DSIN(Xpivot)                 ! moment arms
         RCOSPK=R*DCOS(Xpivot)                 ! ......
         IPAD=INERPAD(KKK)                     ! pad mass moment of inertia
         KROTPAD=KSTPAD(KKK)                   ! pad rotational stiffness
         CROTPAD=CDAPAD(KKK)                   ! pad rotational damping
       END IF

       RETURN
C !------------------------------!
       END IF
C !------------------------------! KKK=KPAD (current pad)



C !------------------------------!   KKK=0  write mesh to datafile

  500 CONTINUE                   !-> STORE MESH DATA in PLOTMESH File
      KP=K
      HGROV=-2.0D0               !## assumed groove depth

      IF (BEARING.EQ.2) GOTO 666 ! ==> FOR ANNULAR SEAL/BEARING

c     ---------------
      IF (KP.GT.1) THEN
          I=1
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HGROV
          end do
        IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HGROV
          end do
        END IF

          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HGROV
          end do
      END IF      ! at leading edge of pad K+1: Groove depth
c     ---------------

      IF (NPOCKET.EQ.0) GOTO 610

c     ---------------
  606 DO K=1, NPOCKET
c     ---------------
        DO II=1, NLC
              I=(K-1)*Nxi+II
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
        END DO                    ! II=1, NLC

        KK=(K-1)*Nxi+Nlc-1

        DO II=1, NPC
          i=KK+II
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), -1.0D0  ! -HRR
          end do
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
        END DO       ! II=1, NPC
c    !---------------!..............
      END DO         ! K=1, NPOCKET
c    !---------------!..............


c    .....................!
      IF (IFULL.EQ.0) THEN
c    .....................!
        JJ=NPOCKET*NXI

        DO II=1, NLC
          i=JJ+II
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
        END DO                    ! II=1, NLC
c    ......................!
      END IF
c    ......................!

      GOTO 612

c    ----------------------!
  610 DO I=1, NXT          ! NO POCKET ON PAD:
c    ----------------------!
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
       IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
       END IF
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
       END DO              ! II=1, NLC
c    ----------------------!
       I=NXT

C     ................................. !
  612 CONTINUE

       IF (KP.EQ.NPAD) GOTO 777

c    ----------------------!
       IF (KP.LT.NPAD) THEN
c    ----------------------!
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), HGROV
          end do
        IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), HGROV
          end do
        END IF
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), HGROV
          end do
c    ----------------------!
       END IF              ! at trailing edge of pad K
c    ----------------------!


       K=KP
       K=K+1
       GOTO 300



c    ----------------------!
 666  DO I=1, NXT          ! ANNULAR SEAL
c    ----------------------!
          do j=0, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(i),YP(j), 0.0D0
          end do
       END DO              ! II=1, NXT
c    ----------------------!

c    !...............!
 777   CONTINUE
       CLOSE (2)
       WRITE (6,25)
c    !...............! KKK=0

 25   FORMAT(' ',3X,'XYDATA >> PLOTMESH File completed')


      RETURN

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  997 WRITE (6, 998)
  998 FORMAT ('$',' ERROR: Size of DX is zero or negative, RESULTS',
     +        ' will be in error: REVISE DATA')
      CALL BEEPER
      RETURN

 1000 WRITE (6, *) 'Error opening file:'
      WRITE (6, *) 'Aborting write.'
      CLOSE (UNIT=2)
      STOP

 1100 WRITE (6, *) 'Error on file name:'
      WRITE (6, *) 'Aborting program'
      CLOSE (UNIT=2)
      STOP


 1200 WRITE (6, *) 'Error during write:'
      WRITE (6, *) 'Aborting write.  The file is being deleted.'
      CLOSE (UNIT=2, STATUS='DELETE')
      STOP

c ......................................!.........................
 611  FORMAT (' ','ERROR: Total Recess Length on Pad', I3,
     +        ' is larger than Pad Length',/,1X,
     +        'PROGRAM STOPS: Review Input Data')
 613  FORMAT (' ','ERROR: Total Length of Pads=', E12.5E2,
     +        '[m] is larger than allowed 2PIxDIAM',/,1X,
     +        'PROGRAM STOPS: Review Input Data')
 615  FORMAT (' ','ERROR: AXIAL length of U or V CVs is zero',
     +        ' or negative, DYPU,DYLU=', 2(E12.5E2,','), /,1X,
     +        '   PROGRAM will give erroneous results')
 616  FORMAT (' ','ERROR: AXIAL length of U or V CVs is zero',
     +        ' or negative, DYLU=',E12.5E2,',', /,1X,
     +        '   PROGRAM will give erroneous results')

  221 FORMAT ( ' ', 70('='),/,1X,'PAD No=', I3,/,1X,
     +        'j  Xp=Xv    Xu   Dxp=Dxv   Dxu   Suw    Sue ',/,
     +        '---------------------------------------------')
  222 FORMAT ( I3, 6(1X,F6.3) )


  223 FORMAT (/,'  k  Yp=Yu    Yv   Dyp=Dyu  Dyv    Svs     Svn',
     +        /,'----------------------------------------------')
  123 FORMAT ( I3, 6(1X, F6.3))

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Minmax                                                      **
C **                                                                         **
C **  MINMAX:                                                                **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE MINMAX(X, XMIN, XMAX, XTIC, N)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION X(*), XMIN, XMAX, XTIC
      INTEGER N, I

C ----------------------------------------------------------------------------
C --  MINMAX code                                                           --
C ----------------------------------------------------------------------------

      XMAX=-1.00D+30                          !
      XMIN=1.00D+30                           !
      DO I=1, N                               !
          XMAX=DMAX1(XMAX, X(I))              !
          XMIN=DMIN1(XMIN, X(I))              !
      END DO                                  !
      XTIC=DABS(XMAX-XMIN)/2.0D0              !
      IF (XTIC.EQ.(0.0D0)) THEN               !
          WRITE (6, *) 'Tic marc = 0, STOP.'  !
          CALL BEEPER                         !
          STOP                                !
      END IF                                  !

      END
c
c-----------------------------------------------------------------------------
C LAST REVIEWED 2/17/94 BY DR LUIS SAN ANDRES AT TAMU
c-----------------------------------------------------------------------------