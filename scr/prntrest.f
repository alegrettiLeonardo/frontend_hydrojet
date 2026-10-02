c 7/3/95 updated for angle injection in HJBS
c -------------------------------------------------------
c 6/20/95 Updated for cell depth effect COMMON/HONEY
c -------------------------------------------------------
c 6/28/94 for radial heat flow parameters
c -------------------------------------------------------
c
C prntrest.f: program for display/print of calculated results
C
C
C #####   #####   #    #   #####  #####   ######   ####    #####          ######
C #    #  #    #  ##   #     #    #    #  #       #          #            #
C #    #  #    #  # #  #     #    #    #  #####    ####      #            #####
C #####   #####   #  # #     #    #####   #            #     #     ###    #
C #       #   #   #   ##     #    #   #   #       #    #     #     ###    #
C #       #    #  #    #     #    #    #  ######   ####      #     ###    #
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

C *****************************************************************************
C **                                                                         **
C **  Subroutine Display                                                     **
C **                                                                         **
C **  DISPLAY:  Calculates Xydata and displays input data.                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE DISPLAY(DEVICE,Ichec2,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARPAD/ X1(MAXNPAD),LPAD(MAXNPAD),
     +                X1r(MAXNPOCK,MAXNPAD),Lrec(MAXNPOCK,MAXNPAD),
     +                PaDorif(MAXNPOCK,MAXNPAD)
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /INERTPAD/ INERPAD(MAXNPAD)
      COMMON /PARAPAD/ KSTPAD(MAXNPAD), CDAPAD(MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /PARAM2/ EXO, EYO

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL

      DOUBLE PRECISION X1,LPAD,X1r,Lrec,PaDorif,
     +                 EXO,EYO,EWX,EWY,EWEAR,BETAW,
     +                 RSINPK,RCOSPK, IPAD, CROTPAD,KROTPAD,
     +                 INERPAD, KSTPAD, CDAPAD, ROTPAD,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HREC,
     +                 PRELOAD, OFFSET,ROTDEL

      INTEGER NPAD,NREC,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,
     +        NXI,NYI,NXT,IFULL, TILT, TILTPAD
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER CHOICE, DEVICE, ichec2,ichec1, ICASE, IWEAR, KPAD
      DOUBLE PRECISION R, CONV, Xpivot,RPM,PS,PA

C ----------------------------------------------------------------------------
C --  DISPLAY code                                                          --
C ----------------------------------------------------------------------------

         R=DIAM/2.0D0                  ! Bearing Radius
         CONV=DACOS(-1.0D0)/180.0D0    ! PI/180. Conversion deg to radians

      ICASE = 2
c...................................................!
c MENU:
c...................................................!


  200 WRITE (6, *) 'SELECT  one of:'
      WRITE (6, *) ' (1) Print parameters.'
      WRITE (6, *) ' (2) Print X, Y and H vectors for mesh.'
      WRITE (6, *) ' (3) Summary of Results'
      WRITE (6, *) ' (4) Print MAX/MIN Speeds, Reynolds & Mach#s'
      WRITE (6, *) ' (5) List of dimensionless P,U,V, & T fields'
      WRITE (6, *) ' (6) Create files for P, T, and H-Film Plots'
      WRITE (6, *) ' (7) Continue.'
c     WRITE (6, *) ' (8) Print CURRENT Fluid Properties'

      WRITE (6, 50)
   50 FORMAT ('$', 'Choice is: ')
  100 READ (5, *) CHOICE
      IF ((CHOICE.LT.1).OR.(CHOICE.GT.8)) THEN
          GOTO 200
      END IF
C...................................................!
C OPTION 1: Prints Input Data & Parameters
C...................................................!

      IF (CHOICE.EQ.1) THEN

          CALL PRINPUT(DEVICE,RPM,PS,PA)

C...................................................!
C OPTION 2: Prints X&Y vectors for mesh             !
C...................................................!

      ELSE IF (CHOICE.EQ.2) THEN

          CALL DISPVECTOR(DEVICE,RPM,PS,PA)
C...................................................!
C OPTION 3: PRINT SUMMARY OF RESULTS                !
C           for each PAD
C...................................................!

      ELSE IF (CHOICE.EQ.3) THEN

         CALL ZERO0
         CALL ZERO1

         DO KPAD=1, NPAD
c       !.........................!
           CALL READTEMP(KPAD)    ! Reads DATA from FILE

           IF (TILTPAD.EQ.1) THEN                       ! Calculate location
            Xpivot=(OFFSET*LPAD(KPAD)+X1(KPAD))*CONV    ! of K-pad pivot in rads.
            RSINPK=R*DSIN(Xpivot)                       ! moment arms
            RCOSPK=R*DCOS(Xpivot)                       ! ......
            ROTDEL=ROTPAD(KPAD)                         ! angular pad rotation
            IPAD=INERPAD(KPAD)                          ! pad mass moment of inertia
            KROTPAD=KSTPAD(KPAD)                        ! pad rotational stiffness
            CROTPAD=CDAPAD(KPAD)                        ! pad rotational damping
           END IF                                       !<<< TILT-PAD

           WRITE (6,111) KPAD
           IF (DEVICE.EQ.1) WRITE (1,111) KPAD

           CALL PRINTDUVP(DEVICE, 0)    ! Prints Flow Rates & Pressures
           ICHEC1=1

           CALL PRINTF(DEVICE,ICHEC1)   ! Prints Forces

           IF (ICHEC2.EQ.1) CALL PRINTCOEF(Device)

           CALL ADD0               ! Add forces/moments
           CALL ADD1               ! Add force coefficients

         END DO
c       !..........................!

c       !..........................! Prints TOTAL forces and
         IF ((NPAD.GT.1).OR.(TILT.EQ.1)) THEN       ! coefficients.
           CALL PRINT0PAD(DEVICE)
           IF (ICHEC2.EQ.1) CALL PRINT1PAD(DEVICE,RPM,PS,PA)
         END IF

C...................................................!
C OPTION 4: PRINT MAX fluid SPeeds & Reynolds Numbers
C...................................................!

      ELSE IF (CHOICE.EQ.4) THEN
         CALL CPARAM(RPM,PS,PA)
         DO KPAD=1, NPAD
           WRITE (6,111) KPAD
           CALL XYDATA(KPAD,RPM,PS,PA)
           CALL READTEMP(KPAD)
           IF (IWEAR.eq.1) THEN
              CALL HWEAR(EXO, EYO, NXT, NYI)
           ELSE
              CALL FILMH(EXO, EYO, NXT, NYI)
           END IF
           CALL FILMC                        ! ### update for compliance
           CALL PRINTUVP(DEVICE,0)
           CALL FINPROP(3,DEVICE,RPM,PS,PA)            ! ICASE=2 print
         END DO

C...................................................!
C OPTION 5: LIST of dimensionless values of U,V,P
C...................................................!

      ELSE IF (CHOICE.EQ.5) THEN

         DO KPAD=1, NPAD
c       !.........................!
           CALL XYDATA(KPAD,RPM,PS,PA)      !
           CALL READTEMP(KPAD)    ! Reads Results/Data
           WRITE (6,111) KPAD
           IF (IWEAR.eq.1) THEN
              CALL HWEAR(EXO, EYO, NXT, NYI)
           ELSE
              CALL FILMH(EXO, EYO, NXT, NYI)
           END IF
           CALL FILMC                        ! ### update for compliance
           IF (DEVICE.EQ.1) WRITE (1,111) KPAD
           CALL PRINTFUVP(DEVICE,0)
           CALL PFIELDS(DEVICE)
           CALL PAUSE
         END DO
c       !.........................!


C...................................................!
C OPTION 6: CREATE Pres, Temp & Film Data files for plot
C...................................................!

      ELSE IF (CHOICE.EQ.6) THEN

          CALL FILEPLOT(RPM,PS,PA)

C...................................................!
C OPTION 7: EXITS Menu
C...................................................!

      ELSE IF (CHOICE.EQ.7) THEN

          RETURN

C...................................................!
C OPTION 8: LIST of dimensionless Fluid Properties
C...................................................!
C Only lists props for last(current) pad of array, otherwise
c it takes too much time for a cryogen.
c FOR DEBUGGING PURPOSES ONLY

      ELSE IF (CHOICE.EQ.8) THEN

           CALL PROPERTIES(DEVICE)
C...................................................!

      END IF
c     END of MENU Selection

C...................................................!

      GOTO 200

 111  FORMAT (' ',/,' +',77('-'),'+',/,' |',3X,
     1         'RESULTS FOR BEARING PAD #',I2,46X,' |')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Dispvector                                                  **
C **                                                                         **
C **  DISPVECTOR:  Prints X, Y, and H vectors for mesh.                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE DISPVECTOR(DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT),DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)

      COMMON /SPLDATA/Z(Nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl),Nj
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c.........................................................................
      DOUBLE PRECISION XP, XU, YP, YV,DYP, DYV, SVN, SVS,
     +                 DXP, DXU, SUW, SUE,HP, HU, HV,
     +                 EXO,EYO,EWX,EWY,EWEAR,BETAW,Z,Cl,Bcl,Ccl,Dcl,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HREC,
     +                 ClearO,ClearR,ClearL,YR,YL,AXO,AYO,ZO
      INTEGER Nj,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        IWEAR, ISYM, ICSTEP, NPAD, NREC, IFULL,BEARING


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION Rj, Xsip,Xsiv, Cp, Cv,MAXHP,MINHP,
     +                 MAXIN, MININ, Px,RPM,PS,PA
      INTEGER Device,J, K, Kstart, U, Iwrite, KPAD

C ----------------------------------------------------------------------------
C --  DISPVECTOR code                                                       --
C ----------------------------------------------------------------------------
      Rj=DIAM/2.0D0
      Iwrite=1
      U=6

 1    CONTINUE

      WRITE (U, 1100)
      WRITE (U, 1110) CLEAR, EXO, EYO
      IF (IWEAR.EQ.1) THEN
         WRITE (U, 1120) EWX,EWY,EWEAR,BETAW
      END IF
      IF ( (AXO.NE.0.0D0).OR. (AYO.NE.0.0D0) ) THEN
         WRITE (U,1100)
         WRITE (U,1130) AXO,AYO, ZO
      END IF

      WRITE (U,1100)

      DO KPAD=1, NPAD
c    !-----------------! Sweep over all Pads
      NPOCKET=NREC(KPAD)
      CALL XYDATA(KPAD,RPM,PS,PA)
      CALL READTEMP(KPAD)    ! Reads Results/Data

      IF (IWEAR.eq.1) THEN
         CALL HWEAR(EXO, EYO, NXT, NYI)
      ELSE
         CALL FILMH(EXO, EYO, NXT, NYI)
      END IF

      CALL FILMC                        ! ### update for compliance

      MAXHP=MAXIN(HP)
      MINHP=MININ(HP)
C ......................................!
      IF (BEARING.EQ.1) THEN
        WRITE (U,1005) KPAD, NPOCKET, NLC, NPC-2, NXI, NXT
      ELSE
        WRITE (U,1006) KPAD, NXT
      END IF

      WRITE (U, 300)
c ..................................... !Write X, DX vectors at P&U nodes
c                                       ! and Film at inlet(z=0) and exits
c
C    !..................................!
      IF ((ISYM.eq.1).OR.(BEARING.EQ.2)) THEN
C    !..................................!
         Kstart=1
         WRITE (U,400)
         DO J=1, NXT
           WRITE (U, 500) J, XP(J), XU(J), DXP(J), DXU(J),
     +              HV(J,1), HP(J,NYI),SUW(J), SUE(J)
         END DO

C    !..................................!
      ELSE
C    !..................................!

         Kstart=0
         WRITE (U,401)
         DO J=1, NXT
           WRITE (U, 501) J, XP(J), XU(J), DXP(J), DXU(J),
     +              HV(J,1), HP(J,-NYI),HP(J,NYI),SUW(J), SUE(J)
         END DO
C    !..................................!
      END IF
C    !..................................!

      WRITE (U,1100)
      WRITE (U, 650) MAXHP, MINHP


c    !-----------------!      DONE on Sweep over all Pads
       END DO
c    !-----------------!      KPAD=1, 2, ...., NPAD

      WRITE (U,1100)
      IF (BEARING.EQ.1) THEN
         WRITE (U, 2005) NLA, NPA, NYI
      ELSE
         WRITE (U, 2006) NYI
      END IF

      WRITE (U, 700)
      WRITE (U, 800)



c    !----------------------------------!@
      DO K= Kstart, NYI
c    !----------------------------------!@
           Xsip=YP(k)*Rj
           Xsiv=YV(k)*Rj

C       ................................!
        IF (ICSTEP.eq.1) THEN
c       ................................! STEP CLEARANCE BEARING
           IF (Xsip.gt.YR) THEN
              Cp=ClearR/CLEAR
           ELSE
              Cp=Clearo/CLEAR
           END IF
           IF (Xsiv.gt.YR) THEN
              Cv=ClearR/CLEAR
           ELSE
              Cv=Clearo/CLEAR
           END IF

C        ...............................!
         ELSE
c        ...............................! CONTINUOUS CLEARANCE C(Y)

           call seval(Nj,Xsip,Z,CL,BCL,CCL,DCL,Cp)
           call seval(Nj,Xsiv,Z,CL,BCL,CCL,DCL,Cv)
           Cp=Cp/CLEAR
           Cv=Cv/CLEAR

C        ...............................!
         END IF
c        ...............................!

          WRITE (U, 900) K, YP(K), YV(K), DYP(K), DYV(K),
     +                   SVS(K), SVN(K),Cp,Cv
c    !----------------------------------!@
      END DO
c    !----------------------------------! K=Kstart, Nyi

      IF (BEARING.EQ.2) GOTO 666        ! FOR ANNULAR SEAL



c    !======================!
      IF (ISYM.eq.0) THEN
c    !======================!
      Px=1.0D0
      IF (Nj.eq.2) Px=-1.0D0      !#### TAPERED SYMMETRIC

c      !----------------------------------!
        DO K= -1, -NYI, -1
c      !----------------------------------!

c        ...............................!
          IF (ICSTEP.eq.1) THEN
c       ................................! STEP CLEARANCE BEARING
             Xsip=YP(k)*Rj
             Xsiv=Yv(k)*Rj
             IF (Xsip.LT.YL) THEN
                Cp=ClearL/CLEAR
             ELSE
                Cp=ClearO/CLEAR
             END IF
             IF (Xsiv.LT.YL) THEN
                Cv=ClearL/CLEAR
             ELSE
                Cv=ClearO/CLEAR
             END IF

c        ...............................!
           ELSE
c        ...............................! CONTINUOUS CLEARANCE C(Y)
             Xsip=YP(k)*Rj
             Xsiv=Yv(k)*Rj
             call seval(Nj,Px*Xsip,Z,CL,BCL,CCL,DCL,Cp)
             call seval(Nj,Px*Xsiv,Z,CL,BCL,CCL,DCL,Cv)
             Cp=Cp/CLEAR
             Cv=Cv/CLEAR

c        ...............................!
           END IF
c         ...............................!

          WRITE (U, 900) K, YP(K), YV(K), DYP(K), DYV(K),
     +                   SVS(K), SVN(K),Cp,Cv

c      !----------------------------------!
        END DO
c      !----------------------------------! K=-1,Nyi,-1

c    !======================!
      END IF                              ! ASYMMETRIC BEARING
c    !======================!


 666  WRITE (U,1111)

      IF (Iwrite.eq.2) RETURN

c    !--------------------------!
      IF (DEVICE.eq.1) THEN
c    !--------------------------! --> DUMP.FILE
          U=1
          Iwrite=2
          GOTO 1
      END IF
c    !--------------------------! --> DUMP.FILE


C........................................... FORMAT statements
 1111 FORMAT (' ', 3X, 76('='),
     +        /,4X,'NOTE: '
     +         'AXIAL Clearance values DO NOT reflect',
     +        /,10X, 'variations with PAD preload',
     +        /,4X, 76('='))

 1005 FORMAT (' ', 3X, 'PAD #:',  I2, ' with Nrecess:', I2, /, 4X,
     +        'Number grid points(CVs) circumferential:', /, 4X,
     +        'Inter-recess, Nlc:', I2, 3X, ', on recess Npc-2:', I2,
     +        ' Total, Nxi:', I3, 2X, 'Nxt:', I3)
 2005 FORMAT (' ', 3X, 'Number grid points(CVs) axially:', /, 4X,
     +        'Inter-recess, Nla:', I2, 3X, ', Above Npa:',
     +        I2, 3X, 'Total, Nyi:', I3)

 1006 FORMAT (' ', 3X, 'PAD #:', I2, 1X,
     +        'Number grid points(CVs) circumferential: Nxt:', I3)
 2006 FORMAT (' ', 3X, 'Number grid points(CVs) axial : Nyi:', I3)


  300 FORMAT (' ', /, 10X,
     +        'nonDIM X vectors, DeltaX sizes and film thickness',
     +        ' at Inlet & Exit')
  400 FORMAT (' ', 3X, 'j', 2X, 'Xp=Xv', 6X, 'Xu', 3X,
     +        'Dxp=Dxv', 4X, 'Dxu', 5X, 'Ho', 5X, 'Hr', 7X,
     +        'Suw', 4X, 'Sue', /, 3X, 75('-'))
  500 FORMAT (' ', 2X, I2, 2X, 8(F6.3, 2X))
  401 FORMAT (' ', 3X, 'j', 2X, 'Xp=Xv', 6X, 'Xu', 3X,
     +        'Dxp=Dxv', 4X, 'Dxu', 5X, 'Ho', 5X, 'Hl', 5X,
     +        'Hr',5X, 'Suw', 4X, 'Sue', /, 3X, 79('-'))
  501 FORMAT (' ', 2X, I2, 2X, 9(F6.3, 2X))
  600 FORMAT (' ', 2X, I2, 2X, F6.3, 10X, F6.3, 10X, F6.3)
  700 FORMAT (' ', /, 10X, 'nonDIM Axial Y vectors, DeltaY sizes',
     +                     ' and Clearance/C*' )
  800 FORMAT (' ', 3X, 'k', 2X, 'Yp=Yu', 6X, 'Yv', 3X, 'Dyp=Dyu',
     +        4X, 'Dyv', 4X, 'Svs', 5X, 'Svn', 4X,
     +            'Cp/C*',4X,'Cv/C*',/, 3X, 70('-'))
  650 FORMAT (' ', 3X, 'Film Thickness at Pnodes, MAX:', E12.5E2,
     +         3X, 'MIN:', E12.5E2,/)
  900 FORMAT (' ', 2X, I2, 2X, 8(F6.3, 2X))
 1000 FORMAT (' ', 2X, I2, 2X, F6.3, 10X, F6.3, 10X, F6.3)
 1100 FORMAT (' ', 2X, 76('-'))
 1110 FORMAT (' ',2X, 'TYP Clearance C*=', E12.5E2,'[m]',
     +       3X, 'Eccentricity, Ex/C*=',F6.4,' Ey/C*=',F6.4)
 1120 FORMAT (' ', 2X, 76('.'),/,3X,'WEAR on BEARING, Ewx=',
     +        E12.5E2,'  Ewy=', E12.5E2, ', Ewear=',
     +        E12.5E2,' [m]',/,3X,'Wear Angle relative to -Xaxis:',
     +        F9.4,' deg')
 1130 FORMAT (' ', 2X, 'Misalign Journal Angles AX=',E11.4E2,
     +        ', AY=', E11.4E2,' RAD, ZO=', E11.4E2,'[m]')

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Pfields                                                     **
C **                                                                         **
C **  PFIELDS:  Print P, U, V T & H fields.                                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PFIELDS(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION EXO, EYO, AXO, AYO, ZO, P, T, U, V,
     +                 PRECdo,PRECup,PRjet,
     +                 XP, XU, YP, YV, HP, HU, HV
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL,BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER DEVICE, I, J, K, KMAX, WID, NXX, IU, Iwrite ,
     +        IMI,IMA, II, IJ, IDUMY
C ----------------------------------------------------------------------------
C --  PFIELDS code                                                          --
C ----------------------------------------------------------------------------
      WID=10-NYI+6
      IF (WID.LT.7) WID=7

      Iwrite=1
      IU=6

 1    IF (ISYM.eq.1) THEN
         WRITE(IU,21)
      ELSE
         WRITE(IU,20)
      END IF

      WRITE (IU,410)

C    !....................................................!
C     Calculates Index Numbers for PRINT OUT
C    !....................................................!

      IF (IFULL.EQ.1) THEN
c    !..........................! circular 2PI bearing
      IF ((EXO.eq.(0.0)).AND.(EYO.eq.(0.0))) THEN
         IF ( (AXO.eq.(0.0)).AND.(AYO.eq.(0.0)) ) THEN
           NXX=NXI
           KMAX=1
         ELSE
           NXX=NXT
           KMAX=NPOCKET
         END IF
      ELSE
         NXX=NXT
         KMAX=NPOCKET
      END IF
c    !..........................! Pad bearing
      ELSE
          NXX=NXT
          KMAX=NPOCKET+1
      END IF
c    !..........................!



      GOTO (701,702, 703) , BEARING

C    !------------------------!
C     IF (BEARING.EQ.1) THEN  ! FOR HYDROSTATIC BEARING
C    !------------------------!
C    !....................................................!
C       PRINTS PRESSURE FIELD
C    !....................................................!

 701  WRITE (IU, 40)
      IF (NPOCKET.EQ.0) GOTO 702

      DO K=1, KMAX
         IMI=(K-1)*NXI+1
         IMA=IMI-1+NLC
         DO I=IMI,IMA
             WRITE (IU, 101) XP(I), (P(I,J),J=1,NYI)
             IF (ISYM.eq.0) THEN
               WRITE(IU,102) XP(I), (P(I,J),J=-1,-NYI,-1)
             END IF
         END DO
         IF (K.GT.NPOCKET) GOTO 77

         IMI=(K-1)*NXI+NLC
         IMA=K*NXI+1
         II=0
         DO I=IMI, IMA
          II=II+1
             WRITE (IU, 101) XP(I), ( PRjet(K,II),J=1,NPA-1),
     +                              ( P( I,J), J=NPA,NYI   )
          IF (ISYM.eq.0) THEN
             WRITE (IU, 102) XP(I), ( PRjet(K,II),J=-1,-NPA+1,-1),
     +                              ( P( I,J), J=-NPA,-NYI,-1)
          END IF
         END DO
 77   CONTINUE
      END DO

C    !....................................................!
C       PRINTS HP - FILM FIELD
C    !....................................................!
      WRITE (IU, 451)
      DO I=1, NXX
          WRITE (IU, 101)  XP(I), (HP(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I), (HP(I,J),J=-1,-NYI,-1)
          END IF
      END DO

C    !....................................................!
C       PRINTS X- VELOCITY FIELD
C    !....................................................!
      WRITE (IU, 50)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (U(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I), (U(I,J),J=-1,-NYI,-1)
          END IF

      END DO

C    !....................................................!
C       PRINTS Y-VELOCITY FIELD
C    !....................................................!
      WRITE (IU, 60)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (V(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I),(V(I,J),J=-1,-NYI,-1)
          END IF
      END DO

C    !....................................................!
C       PRINTS TEMPERATURE FIELD
C    !....................................................!
      WRITE (IU, 70)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (T(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I),(T(I,J),J=-1,-NYI,-1)
          END IF
      END DO

C    !....................................................!
      GOTO 707

C    !------------------------------!
C      ELSE IF (BEARING.EQ.2) THEN  ! FOR SEALS
C    !------------------------------!

 702  WRITE (IU, 40)
      NXX=NXT
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (P(I, J), J=0, NYI)
      END DO

      WRITE (IU, 451)
      DO I=1, NXX
          WRITE (IU, 101)  XP(I), (HP(I, J), J=1, NYI)
      END DO

      WRITE (IU, 50)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (U(I, J), J=1, NYI)
      END DO

      WRITE (IU, 60)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (V(I, J), J=1, NYI)
      END DO

      WRITE (IU, 70)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (T(I, J), J=1, NYI)
      END DO
      GOTO 707
C    !....................................................!

C    !------------------------------!
C      ELSE IF (BEARING.EQ.3) THEN  ! FOR PLAIN BEARINGS
C    !------------------------!

 703  WRITE (IU, 40)
      NXX=NXT
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (P(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I), (P(I,J),J=-1,-NYI,-1)
          END IF
      END DO

      WRITE (IU, 451)
      DO I=1, NXX
          WRITE (IU, 101)  XP(I), (HP(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I), (HP(I,J),J=-1,-NYI,-1)
          END IF
      END DO

      WRITE (IU, 50)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (U(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I),(U(I,J),J=-1,-NYI,-1)
          END IF

      END DO

      WRITE (IU, 60)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (V(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I),(V(I,J),J=-1,-NYI,-1)
          END IF
      END DO

      WRITE (IU, 70)
      DO I=1, NXX
          WRITE (IU, 101) XP(I), (T(I, J), J=1, NYI)
          IF (ISYM.eq.0) THEN
             WRITE(IU,102) XP(I),(T(I,J),J=-1,-NYI,-1)
          END IF
      END DO

      GOTO 707
C    !....................................................!


 707  IF (Iwrite.eq.2) RETURN

c    !----------------------------! --> DUMP.FILE
      IF (DEVICE.eq.1) THEN
c    !----------------------------! --> DUMP.FILE
         Iwrite=2
         IU=1
         GOTO 1
      END IF
c    !----------------------------! --> DUMP.FILE

c   10 FORMAT (' ', I2, '  ', 9(F<WID>.<WID-3>, ' '))
c  101 FORMAT ('R: ', F7.3, '  ', 9(F<WID>.<WID-3>, ' '))
c  102 FORMAT ('L: ', F7.3, '  ', 9(F<WID>.<WID-3>, ' '))
   10  FORMAT (' ', I2, '  ', 9(F7.5, ' '))
  101 FORMAT ('R: ', F7.3, '  ', 9(F7.5, ' '))
  102 FORMAT ('L: ', F7.3, '  ', 9(F7.5, ' '))
   20 FORMAT (' ', 'ASYMMETRIC BEARING:',/,1X,20('='))
   21 FORMAT (' ', ' SYMMETRIC BEARING:',/,1X,20('='))
   40 FORMAT (' ', 'Dimensionless Pressures : P Array:')
  451 FORMAT (' ', 'Dimensionless Film thickness: HP array:')
  452 FORMAT (' ', 'Dimensionless Film thickness: HU array:')
  453 FORMAT (' ', 'Dimensionless Film thickness: HV array:')
   50 FORMAT (' ', 'Dimensionless Circ. Vel. : U Array:')
   60 FORMAT (' ', 'Dimensionless Axial Vel. : V Array:')
   70 FORMAT (' ', 'Dimensionless Temperature : T Array:')
  410 FORMAT (' ','Zero-th order Pressure and Velocities in',
     +' dimensionless form',/,' x=X/R :::::::::: FLOW FIELDs',
     +' in axial direction',/,'        R: right side, L: '
     +' left side of bearing',/,' ',78('.'))


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Prinput                                                     **
C **                                                                         **
C **  PRINPUT:  Prints input data parameters                                 **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PRINPUT(DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TITLE/ TITLE
      COMMON /DATES/ DDATE
      COMMON /SPLDATA/Z(Nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl),Nj
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
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /LOSPAD/ LOSLeadP, KLOSPad
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl,LOSXSIyr
      COMMON /PARAM4/ CINLET, CEXIT
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PAA,
     +                 DEN12P12, VIS12P12, P2
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /LIQUID/ TEMPK,VSOUND, IF, IL
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /TGROOVE/ LEMDA,DELTA
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /OILCOEF/ Talpha
      COMMON /HONEY/ HCELL, HCDIM
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      COMMON /THERMAL/ ALFT, UC, TC,Ec
      COMMON /THERMID/ ISOTH
c     ......................................................................
      CHARACTER*60 TITLE
      CHARACTER*10 DDATE
      DOUBLE PRECISION Z,Cl,Bcl,Ccl,Dcl
      DOUBLE PRECISION X1, LPad, X1r, Lrec, PaDorif, ROTPAD, INERPAD,
     +                 KSTPAD, CDAPAD, AC, ETA, RELAXH
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,LENGTHR,LENGTHL,
     +                 EXO, EYO, AXO, AYO, ZO, EWX,EWY,EWEAR,BETAW,
     +                 LOSLeadP, KLOSPad,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr,
     +                 CINLET, CEXIT,ClearO,ClearR,ClearL,YR,YL,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PAA,
     +                 DEN12P12, VIS12P12, P2

      DOUBLE PRECISION Pleft, Pright, Cleft, Cright,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 HRECU, VSUP, BETA,AMOD, BMOD, RUGR, RUGS, EXPO,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 Tempk, Vsound, ALFT, UC, TC,Ec,CPS,THS,BETAKS,
     +                 PRELOAD, PRELOADB, ROTDEL, OFFSET ,Talpha,
     +                 LEMDA, DELTA, TSHAFT, TSTATOR,
     +                 TBOUT, THERMALK, ROUTER, HKB, HCELL, HCDIM,
     +                 ANGLEJ, LOCJET, CJET, DPJET

      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, ISOTH,
     +        Nj,ITER,ITMAX,ITPMAX,IF, IL, IWEAR, ISYM, ICSTEP,
     +        NPAD, NREC, IFULL, BEARING, NLOBES, TILTPAD, LIFT

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER J, DEVICE, Iwrite, U, KPAD
      DOUBLE PRECISION AOCD, ECC, LDIFF
      CHARACTER*8 FLUID
      CHARACTER*40 BCASE
      CHARACTER*50 THERM
      CHARACTER*50 OILALFA
C ----------------------------------------------------------------------------
C --  PRINPUT code                                                          --
C ----------------------------------------------------------------------------
      LDIFF=DABS(LENGTHR/LENGTH-0.5D0)
      AOCD=DACOS(-1.0D0)*DORIF*DORIF*CD/4.0D0
      ECC = DSQRT(EXO*EXO+EYO*EYO)

      IF (IF.eq.1) fluid='HYDROGEN'
      IF (IF.eq.2) fluid='NITROGEN'
      IF (IF.eq.3) fluid='OXYGEN'
      IF (IF.eq.4) fluid='METHANE'
      IF (IF.eq.5) fluid='AIR'
      IF (IF.eq.6) fluid='WATER'
      IF (IF.eq.7) THEN
         fluid='OIL'
         OILALFA='Temperature-viscosity coefficient of oil'
      END IF

      IF(ISOTH.EQ. 1) THERM='Isothermal fluid film(T=Tin=Constant)'
      IF(ISOTH.EQ. 2) THERM='Adiabatic journal & Iso-bearing (Tb=Ts)'
      IF(ISOTH.EQ. 3) THERM='Iso-journal (Tj=Ts) & Adiabatic bearing'
      IF(ISOTH.EQ. 4) THERM='Adiabatic journal & Bearing radial heat'
      IF(ISOTH.EQ. 5) THERM='Iso-journal(Tj=Ts)& Bearing radial heat'
      IF(ISOTH.EQ. 0) THERM='Iso-journal & bearing: ( Tj=Ts, Tb=Ts)'
      IF(ISOTH.EQ.-1) THERM='Adiabatic bounding Surfaces (Qb=Qj=0)'

      IF (BEARING.EQ.1) BCASE='1 RECESS ROW HJB'
      IF (BEARING.EQ.2) BCASE='ANNULAR SEAL'
      IF (BEARING.EQ.3) BCASE='CYLINDRICAL BEARING'

      Iwrite=1
      U=6
c    !.................................!..................................!   1  CONTINUE
c    !.................................!..................................!
      WRITE (U, 240)                   ! PRINT title/date/header
      WRITE (U, 210) TITLE             !
      WRITE (U, 211) DDATE             !
      WRITE (U, 240)                   !
c    !.................................!..................................!
      IF (MODEL.eq.1) THEN             ! two recess row symmetric  bearing
        BCASE='2 RECESS ROW HJB'       !
        WRITE (U, 2111) BCASE          !
      ELSE                             !
        IF (ISYM.eq.1) THEN            ! SYM/ASYM bearing 1 recess row
         WRITE (U, 2113) BCASE         !
        ELSE                           !
         IF (BEARING.eq.2) THEN        !
           WRITE(U, 2115) BCASE        !
         ELSE                          !
           WRITE (U, 2114) BCASE       !
         END IF                        !
        END IF                         !
      END IF                           !
c    !.................................!.................................!

c---- Print FLUID TYPE

      IF (IF.EQ.12) THEN                                ! => unknown fluid
        WRITE(U,213)					!
      ELSE IF (IF.EQ.5.OR.IF.EQ.6.OR.IF.EQ.7) THEN      !=> air, water, oil
        WRITE(U,2141) IF, FLUID				!
      ELSE 						!
        WRITE(U,214) IF, FLUID				!=> cryogens
      END IF						!....................!

      IF ((BEARING.EQ.1).AND.(IF.EQ.5)) THEN		! WARNING
         WRITE( U,2133)					! for AIR hydrostatic bearing
      END IF						!....................!

      WRITE (U,2144) ISOTH, THERM

      WRITE (U,240)

c---- TYPical Clearance and HJB dimensions

      WRITE (U, 1001) CLEAR, DIAM, LENGTH
      IF (BEARING.EQ.1) WRITE (U,1002) AR
      WRITE (U, 240)

      IF ((BEARING.NE.2).AND.(LDIFF.GT.0.01D0)) THEN
        WRITE (U, 2011) LENGTHR, LENGTHL
        WRITE (U,240)
      END IF

c---- PAD & RECESS Location and Angular LENGTHs

c    !.....................!...................................!
      DO KPAD=1, NPAD
c    !.....................! sweep over PADS
        NPOCKET=NREC(KPAD)
        WRITE (U, 2001) KPAD, X1(KPAD),LPAD(KPAD)

        IF (NPOCKET.GT.0) THEN

          WRITE (U, 2002)     ! recess: position, length, orifice dia.
          DO J=1, NPOCKET
             WRITE (U, 2003) J, X1r(J,KPAD), Lrec(J,KPAD),
     +                          PaDorif(J,KPAD)
          END DO

          WRITE (U, 1201) CD,ANGLEJ, LOCJET ! orifice parameters

        ELSE

          WRITE (U, 2004)
        END IF

        WRITE (U, 2005)
c    !.....................!
      END DO
c    !.....................! sweep over PADS ..................!

      IF (IFULL.EQ.1) WRITE (U, 2006)    !=> cylindrical 2PI bearing

c---- preload and offset factors for pad bearings
c---- tilt-pad bearing parameters

      IF (TILTPAD.EQ.1) THEN
        WRITE (6,7032) NLOBES, PRELOAD, OFFSET
        WRITE (U,7033) (ROTPAD(J), J=1,NPAD)    ! pad angle rotation (rads)
        WRITE (U,7035) (INERPAD(J),J=1,NPAD)    ! pad mass moment of inertia
        WRITE (U,7036) (KSTPAD(J), J=1,NPAD)    ! pad rotational stiffness
        WRITE (U,7037) (CDAPAD(J),J=1,NPAD)     ! pad rotational damping
      ELSE
        WRITE (6,7034) NLOBES, PRELOAD, OFFSET
      END IF

      WRITE (U,240)

c---- Clearance function for taper/step/spline bearings

      IF (ICSTEP.eq.1) GOTO 11

c---- TAPERED BEARING

      IF (Nj.eq.2) then
         write (U, 240)
         write (U, 201) cinlet, cexit
         write (U, 240)
      END IF
c
c---- NON Uniform Bearing clearance

      IF (Nj.gt.2) then
          write (U, 240)
          write (U,202) Nj
          do j=1,Nj
            write (U,203) Z(j),Cl(j)
          end do
          write (U,240)
      end if
      GOTO 12

c---- STEP BEARING
 11   IF ((ISYM.eq.1).OR.(BEARING.EQ.2)) THEN
         write (U,240)
         write (U,205) ClearO, ClearR, YR
      ELSE
         write (U,240)
         write (U,206) ClearO, ClearR, YR, ClearL, YL
      END IF
         write (U,240)

 12    CONTINUE

c------ Recess Depth & Supply Orifice line volume HJB

       IF (BEARING.EQ.1) THEN
          WRITE (U, 190) HREC, VSUP
          WRITE (U,240)
       END IF

c ----- Roughness parameters

      IF ((RUGR.GT.(0.0)).OR.(RUGS.GT.(0.0))) THEN
          WRITE (U, 401) RUGS, RUGR
          WRITE (U, 240)
      END IF

c ----  Cell depth parameter

      IF (HCELL.GT.(0.0D0)) THEN
          WRITE (U, 402) HCELL, HCDIM
          WRITE (U, 240)
      END IF

C.............................
      WRITE (U,180) EXO, EYO, ECC
      WRITE (U, 240)
C......................................
      IF ( (AXO.NE.0.0D0).OR. (AYO.NE.0.0D0) ) THEN
         WRITE (U,1130) AXO,AYO,ZO
      ELSE
         WRITE (U,1131) ZO
      END IF
      WRITE (U, 240)
C.............................
      IF (IWEAR.EQ.1) THEN
        WRITE (U, 1120) EWX,EWY,EWEAR,BETAW
        WRITE (U, 240)
      END IF
c.............................................
      WRITE (U, 160) RPM, RPM*DACOS(-1.0D0)/30.0D0
      WRITE (U, 240)
c............................................. KNOWN PRESSURES
      IF (BEARING.EQ.1) THEN
         WRITE (U, 1501) PS, Pleft, Pright, PC
      ELSE IF (BEARING.EQ.2) THEN
         WRITE (U, 1502) PS, Pright, PC
      ELSE IF (BEARING.EQ.3) THEN
         WRITE (U, 1503) Pleft, Pright, PC
      END IF
c...............................................

      WRITE (U, 240)
      CALL WritePright(U)
      CALL WritePleft(U)
      IF (BEARING.EQ.1) THEN
         WRITE (U, 121) PRATIO
         WRITE (U, 240)
      END IF
c................................................. SUPPLY TEMPERATURE
      WRITE (U, 2150) TEMPK, Tempk*9.0D0/5.0D0
c................................................. SHAFT/STATOR TEMPERATURE
      IF(ISOTH.EQ. 2) THEN
        WRITE (U, 2152) TSTATOR,TSTATOR*1.8D0
        WRITE(U,240)
      ELSE IF (ISOTH.EQ.3) THEN
        WRITE (U, 2151) TSHAFT, TSHAFT*1.8D0
        WRITE(U,240)
      ELSE IF (ISOTH.EQ.0) THEN
        WRITE (U, 2152) TSTATOR,TSTATOR*1.8D0
        WRITE (U, 2151) TSHAFT, TSHAFT*1.8D0
        WRITE(U,240)
      ELSE IF (ISOTH.EQ.4) THEN
        WRITE (U, 2154) TBOUT,TBOUT*1.8D0, THERMALK, ROUTER
      ELSE IF (ISOTH.EQ.5) THEN
        WRITE (U, 2151) TSHAFT, TSHAFT*9.0D0/5.0D0
        WRITE (U, 2154) TBOUT,TBOUT*1.8D0, THERMALK, ROUTER
      END IF
c.............................................

      IF (IF.eq.7) THEN                        ! FLUID=OIL
         WRITE (U,240)
         WRITE (U,2145) Talpha, OILALFA
         WRITE (U,240)
      END IF

      IF (BEARING.LE.2) THEN
          WRITE (U, 170) EMUs, RHOs
      END IF

      IF ((ISYM.eq.1).OR.(BEARING.EQ.2)) THEN
        WRITE (U,171) EMUri, RHOri
      ELSE
        WRITE (U,172) EMUle, RHOle
        WRITE (U,171) EMUri, RHOri
      END IF
        WRITE (U, 240)
c.............................................

      IF (BEARING.EQ.1) THEN
        WRITE (U, 250) BETA
        WRITE (U, 1300) LOSXSIxu,LOSXSIxd,LOSXSIyl,LOSXSIyr,ALPHA
      ELSE IF (BEARING.EQ.2) THEN
        WRITE (U, 1301) LOSXSIyr, ALPHA
      END IF

c.............................................
      IF (IFULL.EQ.0) THEN
        WRITE (U, 1302) LOSleadP, LEMDA
        WRITE (U, 240)
      END IF
c.............................................

      IF (BEARING.EQ.2) THEN
         WRITE (U, 1311) Cright
      ELSE IF (BEARING.EQ.1) THEN
         WRITE (U, 1312) Cleft, Cright
      END IF
c      WRITE (U, 240)
c.............................................

      IF (BEARING.EQ.1) THEN 				 !=> HJB
        WRITE (U, 1101)  ITMAX, ITPMAX, ALFU, ALFP ,ALFT
        WRITE (U, 1005)  MPEPS, SFLOW
      ELSE    				                 !=> SEAL/BEARING
        WRITE (U, 1102)  ITMAX, ALFU, ALFP ,ALFT
        WRITE (U, 1006)  MPEPS, SFLOW
      END IF
      WRITE (U, 240)

c.............................................

      WRITE (U, 80) AMOD, BMOD, EXPO
      WRITE (U, 240)

c.............................................
      IF (BEARING.EQ.1) THEN
        WRITE (U, 2301) NLC, NPC-2, NXI, NXT
        WRITE (U, 2201) NLA, NPA, NYI
        WRITE (U, 240)
        WRITE (U, 3001) LENGTH/DIAM, AR/LENGTH,
     +              CLEAR*2.0D0/DIAM, HREC/CLEAR
c    +         DBLE(NPOCKET)*BR/DIAM/DACOS(-1.0D0)
c##     WRITE (U, 601) CORIF
c.............................................
      ELSE
        WRITE (U, 2302) NXT, NYI
        WRITE (U, 240)
        WRITE (U, 3002) LENGTH/DIAM,CLEAR*2.0D0/DIAM

      END IF
c.............................................
      WRITE (U, 602) SPEED, RENC
      WRITE (U, 50) REP, REY
      WRITE (U, 603) PSA/1.0D6,PAA/1.0D6
      WRITE (U, 240)
C...........................................................
      IF (AC.GT.0.0D0) THEN		! bearing compliance
         WRITE (U,615) AC, ETA, LIFT    ! parameters
         WRITE (U,240)                  !
      END IF			        !
C...........................................................


      IF (Iwrite.eq.2) RETURN
C............................................................
      IF (DEVICE.eq.1) THEN
C !-----------------------! -> DUMP FILE
         Iwrite=2
         U=1
         GOTO 1
      END IF
C !-----------------------!

C
C-------------------------------------------------------! PRINT FORMATS
C
  615 FORMAT(' ',3X,'COMPLIANT BEARING, ac=',F8.5,2X,'loss ',
     +       'factor n=',F8.5,2X,'lift:',I2)
 2150 FORMAT(' ',3X,'FLUID Supply Temperature   TS=',F9.4, ' K =',
     +           2X,F9.4,' R',/,4X,76('.'))
 2151 FORMAT(' ',3X,'SHAFT(Journal) Temperature TJ=',F9.4, ' K =',
     +           2X,F9.4,' R')
 2152 FORMAT(' ',3X,'STATOR(Bearing)Temperature TB=',F9.4, ' K =',
     +           2X,F9.4,' R')
 2154 FORMAT(' ',3X,'BEARING Outer Temperature TBout=',F9.4,' K =',
     +           2X,F9.4,' R',
     +        /, 10X,'Thermal Conductivity:', E11.4E2,'[watt/mK]'
     +         ,  2X,'Outer Radius:', F9.4, '[m]' )
  213 FORMAT(' ',3X,'IF=12, FLUID Type: UNKNOWN',2X,
     +       ':Uses a linear function for properties')
 2133 FORMAT(' ',3X,'WARNING: Orifice eqns. DO NOT account for'
     +       ' AIR expansion factor')
  214 FORMAT(' ',3X,'IF=',I2,' FLUID Type:',A8,2X,
     +      ':Uses MIPROPS for properties as f(P,Ts)')
 2141 FORMAT(' ',3X,'IF=',I2,' FLUID Type:',A8)

 2144 FORMAT(' ',3X,'ISOTH=',I2,' THERMAL Case:',A50)
 2145 FORMAT(' ',3X,'Talpha=',F8.6,' (1/K-deg): ',A50)

  250 FORMAT (' ', 19X, 'Compressibility B=', E12.5E2,
     +        ' m2/N at Recess Pressure',/,' ', 3X,76('.'))
   10 FORMAT (' ', 3X, 'CONV FACTORS, Cflow:', E12.5E2, 'm3/s', X,
     +        'Cfor:', E12.5E2, 'N', X, 'CTor:', E12.5E2, 'Nm')
 3001 FORMAT (' ', 3X, 'L/D=', F8.3, 3X, 'Y=Ar/L=', F6.4,
     +        3X, 'C/R:', F8.6, 3X, ',  Hrec/C:', F7.4)
c     +        , 3X, 'X=Nrec x Br/(PIxD)=', F6.4)
 3002 FORMAT (' ', 3X, 'L/D=', F8.3, 3X,  'C/R:', F8.6)

  401 FORMAT (' ', 3X, 'Non smooth bearing, effective relative',
     + ' roughness:', /, 4X, 'Rugs(bearing)=', E12.5E2, 3X,
     +                       'Rugr(journal)=', E12.5E2)
  402 FORMAT (' ', 3X, 'CELL Depth, Hcell=', E12.5E2, '[m]',
     +             3X, 'Hcell/TypC=', F8.4)
   50 FORMAT (' ', 3X,
     + 'Poiseuille Reynolds # (Rep): (rho/emu)Psa c**3/uR = ',
     + E12.5E2, /, 4X,
     + 'Flow Reynolds # (Rep*) : Rep (C/R) =', E12.5E2)
  601 FORMAT (' ', 3X, 'AVE. Orifice Param: ',
     + 'Ao Cd u/C**3/SQR[p(Ps-Pa)/2] = ', E12.5E2)
  602 FORMAT (' ', 3X,
     + 'Speed Parameter: Omega EMU (R/C)**2/(Psa)=', E12.5E2, /,
     + 4X, 'Rec= (rho/emu) Omega R C =', E12.5E2)
  603 FORMAT (' ',3X,
     + 'dim. pressure p=(P-PA)/PSA where: PSA=', E12.5E2,' and',
     + ' PA=', E12.5E2, ' MPa')
   80 FORMAT (' ', 3X, 'FLOW TURBULENCE MODEL based on Moody', 1H',
     +        's Friction Factor Formulae with coeffs:', /, 4X,
     +       'A=',E12.5E2,', B=', E12.5E2,', EXPO=', E12.5E2)

 1005 FORMAT (' ', 3X, 'Conv. Crit. for Mass Sources (MPEPS):',
     + E12.5E2, /,
     +    4X, 'Conv. Crit. for recess flow  (SFLOW):', E12.5E2)
 1006 FORMAT (' ', 3X, 'Conv. Crit. for Mass Sources (MPEPS):',
     + E12.5E2, /,
     +    4X, 'Conv. Crit. for dim. pressures (PEPS):',E12.5E2)


 1101 FORMAT (' ', 3X, 'Itmax:', I4, 3X, 'Itpmax:', I4, 3X,
     + ', Under_Relax Param. Au=',F5.3,', Ap=',F5.3,', At=',F5.3)
 1102 FORMAT (' ', 3X, 'Itmax:', I4, 3X,
     + ', Under_Relax Param. Au=',F5.3,', Ap=',F5.3,', At=',F5.3)


 1201 FORMAT (/,' ', 3X, 'Cd = ORIFICE Coeffic:', F7.5,
     +        2X,'ANGLE:',F7.3,'deg, Location:',F7.3,' within REC')
 121  FORMAT (' ', 3X, 'Recess Pressure ratio (Pr-Pa)/(Ps-Pa):',
     + F7.5, '  at concentric position, Ecc=0')
 1300 FORMAT (' ', 3X, 'Recess Edge Loss Factors (XSI):',
     +       'Circ.-X:(u)', F7.3,' (d)', F7.3,/,
     +   35X,'Axial-Y:(L)', F7.3, ' (R)',F7.3,/, 4X,
     + 'Alpha:', F7.3 , ', entrance swirl velocity factor')
 1301 FORMAT (' ', 3X, 'Edge Loss Factor XSIy:',F7.3,/, 4X,
     + 'Alpha:', F7.3 , ', entrance swirl velocity factor')
 1302 FORMAT (' ', 3X, 'PAD Leading Edge Recovery Factor, LOSleadP:',
     +             F7.3,/,4X,
     +        'Pad-groove Heat Carryover coefficient Lemda:', F7.3)

 1311 FORMAT (' ', 3X,76('.'),/,4X, 'EXIT SEAL Coef: 'E12.5E2)
 1312 FORMAT (' ', 3X,76('.'),/,4X, 'EXIT SEAL Coefs, LEFT:',E12.5E2,
     +        '  RIGHT:',E12.5E2)

 1501 FORMAT (' ', 3X, 'PRESSURES, Psupply:', E12.5E2, 'N/m2,',
     +      /,' ', 3X,' Discharge PLEFT: ', E12.5E2,
     +                         ' PRIGHT:', E12.5E2, 'N/m2',
     +                           ' Pcav:', E11.5E2,'N/m2')
 1502 FORMAT (' ', 3X,          'Psupply: ', E12.5E2,
     +                         ' PRIGHT:', E12.5E2, 'N/m2',
     +                           ' Pcav:', E11.5E2,'N/m2')
 1503 FORMAT (' ', 3X,'Discharge PLEFT: ', E12.5E2,
     +                         ' PRIGHT:', E12.5E2, 'N/m2',
     +                           ' Pcav:', E11.5E2,'N/m2')


  160 FORMAT (' ', 3X, 'ROTATING Speed:', E12.5E2, ' rpm -> Omega =',
     +        E12.5E2, ' rad/s')
  170 FORMAT (' ', 3X, 'at Psupply(Ps):',1X,
     +        'Viscosity:', E12.5E2,
     +        ' Ns/m2 Density:', E12.5E2, ' Kg/m3')
  171 FORMAT (' ', 3X, 'at PRIGHT (PR):',1X,
     +        'Viscosity:', E12.5E2,
     +        ' Ns/m2 Density:', E12.5E2, ' Kg/m3')
  172 FORMAT (' ', 3X, 'at PLEFT  (PL):',1X,
     +        'Viscosity:', E12.5E2,
     +        ' Ns/m2 Density:', E12.5E2, ' Kg/m3')
  180 FORMAT (' ', 3X, 'nonDIM eccentricities:  Exo:', F7.5, 3X,
     +        'Eyo:', F7.5, 3X, 'Ecc:', F7.5)
  190 FORMAT (' ', 3X, 'Recess depth (Hrec):',E12.5E2,' m', 2X,
     +       'Vsup:', E12.5E2,' m3')
 1001 FORMAT (' ',3X,'Clearance(C):',E11.5E2, 'm, Diameter(D):',
     +        E11.5E2, 'm, Length(L):', E11.5E2, 'm')
 1002 FORMAT (' ', 32X,
     +        'Recess, Axial length (Ar):', E12.5E2, ' m  ')

 2011 FORMAT (' ',3X,'Length(R)=',E11.5E2,'m :',
     +               'Length(L)=',E11.5E2,'m')

 2000 FORMAT (' ',3X,'TILT PAD # =',I2, ' Position', 17X,
     +       'Angular',17X,'Rotation',
     +        /,  4X,'=========== ',3X, '  Angle =', E11.4E2,
     +        '[o]  Length=', E11.4E2, '[deg]',
     +        ' in RADS:', E11.4E2)
 2001 FORMAT (' ',3X,'PAD Number =',I3, ' Position', 17X,
     +       'Angular',
     +        /,  4X,'=========== ',3X, '  Angle =', E12.5E2,
     +        '[o]  Length=', E11.4E2, '[deg]')
 2002 FORMAT (15X,'RECESS  Position[o]    Length[o]     Orifice',/,
     +        15X,'          Angle          Angle       do (m)' )
 2003 FORMAT (18X, I2,1X, 3(2X,E12.5E2))
 7032 FORMAT (4X,'TILT  PADS; NPADS=',I3,2X,'PRELOAD:'
     +        ,E12.5E2,'(m)', 2X,'OFFSET:',F5.3)
 7037 FORMAT (4X,'Rot DAMPING (N.m.s):',5(E10.4E2,1X))
 7036 FORMAT (4X,'Rot STIFFN  (N.m/r):',5(E10.4E2,1X))
 7035 FORMAT (4X,'Pad INERTIA (kg.m2):',5(E10.4E2,1X))
 7033 FORMAT (4X,'Pad ROTATION(Rads): ',5(E10.4E2,1X))
 7034 FORMAT (4X,'FIXED PADS; NLOBES=',I3,2X,'PRELOAD:'
     +        ,E12.5E2,'(m)', 2X,'OFFSET:',F5.3)
 2004 FORMAT (' ',3X, 'NO RECESSes ON PAD ')
 2005 FORMAT (' ',3X,37('. '))
 2006 FORMAT (' ',3X,'NOTE: FULL CYLINDRICAL BEARING 360deg')
  201 FORMAT (' ', 3X, 'HJBearing with tapered  axial clearance',
     + /, 4X, 'Cinlet:', E12.5E2, 'm,   Cexit:', E12.5E2,'m')
  202 FORMAT (' ', 3X, 'HJBearing with non-uniform axial clearance',
     + /, 4X, 'Nj: Number of axial points for spline fit=', I3,/,
     +/, 4X, 'Axial Z' , 8X, 'Clearance', /,
     + 4X, '  [m]  ' , 8X, '  [m]', /,
     + 4X, '-------' , 8X, '---------')
  203 format (3X,E12.5E2,2X,E12.5E2)
  205 FORMAT (' ', 3X, 'STEP CLEARANCE: C=', E12.5E2,'m at y=0;',1X,
     +        'CR=',E12.5E2,'m at YR=',E12.5E2,'m')
  206 FORMAT (' ', 3X, 'STEP CLEARANCE: C=', E12.5E2,'m at y=0;',1X,
     +        'CR=',E12.5E2,'m at YR=',E12.5E2,'m',/,44X,
     +        'CL=',E12.5E2,'m at YL=',E12.5E2,'m')
  211 FORMAT (' ',3X, 'Date:',10A)
 2111 FORMAT (' ',3X, 40A, /,4X, 37(' .'))
 2113 FORMAT (' ',3X, 'SYMMETRIC ', 40A,/, 4X, 37(' .'))
 2114 FORMAT (' ',3X, 'ASYMMETRIC ', 40A,/, 4X, 37(' .'))
 2115 FORMAT (' ',3X, ' ', 40A,/, 4X, 37(' .'))

 2201 FORMAT (' ', 3X, 'Number grid points axially:', /, 4X,
     +        'Inter-recess, Nla:',I2,3X, ',  Above Npa:', I2, 3X,
     +        '-> Total, Nyi:', I3)
 2301 FORMAT (' ', 3X,
     +        'Number grid points circumferential:', /, 4X,
     +        'Inter-recess, Nlc:', I2, 3X, 'on recess, Npc-2:', I2,
     +        '   -> Total, Nxi:', I3, 2X, 'Nxt:', I3)
 2302 FORMAT (' ', 3X,
     +        'Number grid points circumferential: Nxt=', I3,
     +        3X, 'Axial Nyi=', I3)

 1120 FORMAT (' ', 3X, 76('.'),/,4X,'WEAR on BEARING, Ewx=',
     +        E12.5E2,', Ewy=', E12.5E2, ', Ewear=',
     +        E12.5E2,' [m]',/,4X,'Wear Angle relative to -Xaxis:',
     +        F9.4,' deg')
 1130 FORMAT (' ', 3X, 'Misalign Journal Angles AX=',E11.4E2,
     +        ', AY=', E11.4E2,' RAD, ZO=', E11.4E2,'[m]')
 1131 FORMAT (' ', 3X, 'Journal Angles AX=AY=0 rad, Pivot ZO='
     +        , E11.4E2,'[m]')
  240 FORMAT (' ', 3X, 76('-'))

  210 FORMAT (' ', 3X, 'hydrojet: SEALs/ HJBs/ Flex Pad JBs',
     +        ' : BEARINGs',3X,'(CR) Texas A&M Univ/94',/, 4X, 60A)
      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Printuvp                                                    **
C **                                                                         **
C **  PRINTUVP:  Prints summary of U, V, and P fields.                       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PRINTUVP(DEVICE, DEST)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/ RHO(MAXNXT,-MAXNYI:MAXNYI),
     +                EMU(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /THERMID/ ISOTH
c     .......................................................................
      DOUBLE PRECISION U, V,P,T,PREC, QREC, TREC, QIN, QOUT, QFACTOR,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PA,
     +                 DEN12P12, VIS12P12, P2,
     +                 ALFT, UC, TC, Ec,ck,betak,hcb,hcj,thc,rho,emu,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 CPS,THS,BETAKS,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HREC
      INTEGER NPOCKET,NLC,NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ISYM, ICSTEP, IFULL, BEARING ,isoth

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER K, J, DEVICE, DEST, IU, Iwrite
      DOUBLE PRECISION MAXIN, MININ ,HCC

C ----------------------------------------------------------------------------
C --  PRINTUVP code                                                         --
C ----------------------------------------------------------------------------
C DEST is the destination of the data: 0 for screen and dump, 1 for
C screen only, and 2 for dump only.  It allows the verbosity to work.

      Iwrite=1
      IU=6
 1    IF ((DEST.EQ.0).OR.(DEST.EQ.1)) THEN
          GOTO 10
      ELSE
          GOTO 20
      END IF


 10     IF (NPOCKET.EQ.0) GOTO 12
          DO K=1, NPOCKET
              WRITE (IU, 800) K, PREC(K), QREC(K)
          END DO
 12       CONTINUE
          HCC=RHOS*CPS*UC*CLEAR/DIAM*2.0D0
          WRITE (IU, 700) QIN, QOUT
          WRITE (IU, 900) MAXIN(P)*PSA+PA,MININ(P)*PSA+PA,PCAV*PSA+PA
          WRITE (IU, 1000) MAXIN(U)*UC, MININ(U)*UC
          WRITE (IU, 1100) MAXIN(V)*UC, MININ(V)*UC
          WRITE (IU, 1200) MAXIN(T)*TC, MININ(T)*TC
          WRITE (IU, 1300) MAXIN(RHO)*RHOS, MININ(RHO)*RHOS
          WRITE (IU, 1400) MAXIN(EMU)*EMUS, MININ(EMU)*EMUS
          WRITE (IU, 1500) MAXIN(CK)*CPS, MININ(CK)*CPS
          WRITE (IU, 1600) MAXIN(BETAK)/TC, MININ(BETAK)/TC
          WRITE (IU, 1700) MAXIN(THC)*THS, MININ(THC)*THS
          WRITE (IU, 1800) MAXIN(HCB)*HCC, MININ(HCB)*HCC
          WRITE (IU, 1900) MAXIN(HCJ)*HCC, MININ(HCJ)*HCC
          WRITE (IU, 2200)
          WRITE (IU, 2300) PSA, PA
          WRITE (IU, 2200)

      IF (Iwrite.eq.2) RETURN

 20   IF ((DEVICE.eq.1).AND. ((DEST.EQ.0).OR.(DEST.EQ.2)) ) THEN
          Iwrite=2
          IU=1
          goto 10
      END IF

      RETURN

  300 FORMAT (' ', 'Index j:', I3, 65('.'))
  700 FORMAT (' ', 3X, 'Inlet Flow: ', E12.5E2, 3X, 'Outlet Flow: ',
     +        E12.5E2)
  800 FORMAT (' ', 3X, 'Recess:', I2, 3X, '-> Pres:', E12.5E2,
     +        3X, 'Flow:', E12.5E2)
  900 FORMAT (' ', 3X, 'LANDS MaxP: ', E12.5E2, 3X, '       MinP: ',
     +        E12.5E2, 2X, 'Pcav:', E12.5E2,' (Pa)')
 1000 FORMAT (' ', 3X, '      MaxU: ', E12.5E2, 3X, '       MinU: ',
     +        E12.5E2,' (m/s)')
 1100 FORMAT (' ', 3X, '      MaxV: ', E12.5E2, 3X, '       MinV: ',
     +        E12.5E2,' (m/s)')
 1200 FORMAT (' ', 3X, '      MaxT: ', E12.5E2, 3X, '       MinT: ',
     +        E12.5E2' (K-deg)')
 1300 FORMAT (' ', 3X, '      MaxRo:', E12.5E2, 3X, '       MinRo:',
     +        E12.5E2,' (kg/m^3)')
 1400 FORMAT (' ', 3X, '      MaxMu:', E12.5E2, 3X, '       MinMu:',
     +        E12.5E2,' (N-s/m^2)')
 1500 FORMAT (' ', 3X, '      MaxCp:', E12.5E2, 3X, '       MinCp:',
     +        E12.5E2,' (J/kg-K)')
 1600 FORMAT (' ', 3X, '      MaxBt:', E12.5E2, 3X, '       MinBt:',
     +        E12.5E2,' (1/K)')
 1700 FORMAT (' ', 3X, '      MaxK: ', E12.5E2, 3X, '       MinK: ',
     +        E12.5E2,' (W/m-K)')
 1800 FORMAT (' ', 3X, '      MaxHb:', E12.5E2, 3X, '       MinHb:',
     +        E12.5E2,' (W/m^2-K)')
 1900 FORMAT (' ', 3X, '      MaxHj:', E12.5E2, 3X, '       MinHj:',
     +        E12.5E2,' (W/m^2-K)')
 2300 FORMAT (' | dim pressure=(P-PA)/PSA where, PSA=', E12.5E2,
     +         1X,'PA=',E12.5E2,'MPa')
 2200 FORMAT (' +', 77('-'), '+')

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Printfuvp                                                   **
C **                                                                         **
C **  PRINTFUVP:  Print summary of U, V, and P fields in final format.       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PRINTFUVP(DEVICE, DEST)

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
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
c     .................................................
      DOUBLE PRECISION U,V, P, T,PREC,TREC, QREC, QIN, QOUT, QFACTOR,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      INTEGER NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ISYM, ICSTEP, IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER K, J, DEVICE, DEST, IU, Iwrite
      DOUBLE PRECISION MAXIN, MININ

C ----------------------------------------------------------------------------
C --  PRINTFUVP code                                                        --
C ----------------------------------------------------------------------------
      Iwrite=1
      IU=6

      IF ((DEST.EQ.0).OR.(DEST.EQ.1)) THEN
          GOTO 10
      ELSE
          GOTO 20
      END IF


 10       WRITE (IU, 1200)
          IF (NPOCKET.EQ.0) GOTO 12
          DO K=1, NPOCKET
              WRITE (IU, 800) K, PREC(K), QREC(K)
          END DO
 12       WRITE (IU, 700) QIN, QOUT
          WRITE (IU, 900) MAXIN(P), MININ(P), PCAV
          WRITE (IU, 1100) MAXIN(V), MININ(V)
          WRITE (IU, 1000) MAXIN(U), MININ(U)
          WRITE (IU, 1300) MAXIN(T), MININ(T)
          WRITE (IU, 1200)

 20   IF (Iwrite.eq.2) RETURN

      IF ((DEVICE.eq.1). AND. ((DEST.EQ.0).OR.(DEST.EQ.2)) ) THEN
          IU=1
          Iwrite=2
          GOTO 10
      END IF

  700 FORMAT (' |    Inlet Flow:', E12.5E2, 3X, 'Outlet Flow:',
     +        E12.5E2, 23X, '|')
  800 FORMAT (' |    Recess:', I2, 3X, '-> Pres:', E12.5E2,
     +        3X, 'Flow:', E12.5E2, 21X, '|')
  900 FORMAT (' |    LANDS MaxP:', E12.5E2, 3X, 'MinP:',
     +        E12.5E2, 3X,'Pcav:', E12.5E2,10X, '|')
 1000 FORMAT (' |          MaxU:', E12.5E2, 3X, 'MinU:',
     +        E12.5E2, 30X, '|')
 1100 FORMAT (' |          MaxV:', E12.5E2, 3X, 'MinV:',
     +        E12.5E2, 30X, '|')
 1300 FORMAT (' |          MaxT:', E12.5E2, 3X, 'MinT:',
     +        E12.5E2, 30X, '|')
 1200 FORMAT (' +', 77('-'), '+')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Printduvp                                                   **
C **                                                                         **
C **  PRINTDUVP:  Prints summary of U, V, and P fields on dimensional form   **
C **              and final format.                                          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PRINTDUVP(DEVICE, DEST)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /SOURCEA/ PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /Pdisch/ Pleft, Pright, Cleft, Cright
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c     ..............................................................
      DOUBLE PRECISION U, V, P,PREC, QREC,TREC, QIN, QOUT, QFACTOR,
     +                 EXO, EYO, ECC, AXO,AYO,ZO
      DOUBLE PRECISION PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW ,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 Pleft, Pright, Cleft, Cright,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ITER, ITMAX, ITPMAX, ISYM, ICSTEP, IFULL, BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER K, J, DEVICE, DEST, Iwrite, IU

C ----------------------------------------------------------------------------
C --  PRINTDUVP code                                                        --
C ----------------------------------------------------------------------------
      IU=6
      Iwrite=1
      Ecc=DSQRT(EXO*EXO+EYO*EYO)

      IF ((DEST.EQ.0).OR.(DEST.EQ.1)) THEN
          GOTO 101
      ELSE
          GOTO 201
      END IF

 101      WRITE (IU,60)
          WRITE (IU,10) ITER, ITPMAX, ISYM
          WRITE (IU,20)EXO,EYO,ECC,AXO,AYO,ZO,RPM,PS,PLeft,Pright
          IF (BEARING.EQ.1) THEN
            WRITE (IU,30) DORIF, CD, CLEFT, CRIGHT
          END IF

          IF (NPOCKET.EQ.0) GOTO 102
          DO K=1, NPOCKET
              WRITE (IU, 40) K, PREC(K), PREC(K)*PSA+PA,
     +                      QREC(K)*QFACTOR
          END DO
 102      WRITE (IU, 60)

 201  IF (Iwrite.eq.2) RETURN

      IF ((DEVICE.eq.1) .AND. ((DEST.EQ.0).OR.(DEST.EQ.2)) ) THEN
         IU=1
         Iwrite=2
         GOTO 101
      END IF

C.................................................! PRINT FORMATS

   10 FORMAT (' |    SOLN FOUND IN Iter=', I4, ', Iterp=',
     +         I2, ' ISYM=', I2, 32X, '|')
   20 FORMAT (' ', '|', 4X, 'Ex=', F5.3, ',Ey=', F5.3,
     +        '; Ec=', F5.3, '; Ax=',E11.4E2,',Ay=',E11.4E2,
     +        ' ZO=',E11.4E2,'|', /,
     +        ' |    rpm=', E12.5E2, '  Ps=', E12.5E2,2X,
     +        'PL=',E12.5E2,' PR=',E12.5E2,' N/m2',2X,'|')
   30 FORMAT (' |    Dorif=', E12.5E2, ' m, Cd=', F7.4,'; Seals (L)=',
     +         E12.5E2,',(R)=', E12.5E2,'|',
     +         /, ' |',77X,'|')
   40 FORMAT (' |    Rec:', I2, 2X, '-> Prec: (', E12.5E2, ')',
     +        E12.5E2, ' N/m2', 2X, 'Flow:', E12.5E2, ' Kg/s |')
   60 FORMAT (' +', 77('-'), '+')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Printfuvp1                                                  **
C **                                                                         **
C **  PRINTFUVP1:  Print first order U, V, and P field in final format.      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PRINTFUVP1(DEVICE, DEST)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /RECES1/ PREC(MAXNPOCK), TREC(MAXNPOCK),
     +                QREC(MAXNPOCK), QIN, QOUT

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX


      DOUBLE PRECISION PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      DOUBLE COMPLEX PREC, QREC,TREC, QIN, QOUT
      INTEGER NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ITER, ITMAX, ITPMAX, IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER K, DEVICE, DEST, IU

C ----------------------------------------------------------------------------
C --  PRINTFUVP1 code                                                       --
C ----------------------------------------------------------------------------
      IF ((DEST.EQ.0).OR.(DEST.EQ.1)) THEN
          IU=6
          GOTO 10
      ELSE
      IF ((DEVICE.eq.1) .AND. ((DEST.EQ.0).OR.(DEST.EQ.2))) THEN
          IU=1
          GOTO 10
      END IF
      END IF
      RETURN

 10       WRITE (IU, 600)
          WRITE (IU, 100) ITER, SMASS, MMP
C          DO K=1, NPOCKET
C              WRITE (IU, 1000) K, PREC(K), QREC(K)
C          END DO
          WRITE (IU, 900) QIN, QOUT
          WRITE (IU, 600)

C
  100 FORMAT (' |    Iter=', I4, ', SUM(Sources):',
     +        E12.5E2, 3X, 'MAX(Mp):', E12.5E2, '   |')
  900 FORMAT (' |    Qin1= (', E12.5E2, 1X, E12.5E2, ')   Qout1= (',
     +        E12.5E2, 1X, E12.5E2, ')  |')
 1000 FORMAT (' |    REC:', I2, '  -> P:', E12.5E2, 1X, E12.5E2,
     +        '   Flow:', E12.5E2, 1X, E12.5E2, ' |')
  600 FORMAT (' +', 77('-'), '+')



      END

C
C *****************************************************************************
C **                                                                         **
C **  Subroutine PRINTF                                                      **
C **                                                                         **
C **  PRINTF  : Prints flow, forces, moments & torque FOR 1-PAD              **
C **                                                                         **
C *****************************************************************************
C
      SUBROUTINE PRINTF(DEVICE,ICHECK)

      IMPLICIT NONE

      INCLUDE 'params.f'

C .....................................................................
      COMMON /PARAM2/EXO,EYO
C     COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /FORCE0/FFACTOR,FX,FY,TO,TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/  KROTPAD, CROTPAD
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL

      COMMON /TMAXMIN/ TKMAX,TKMIN
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /THERMID/ ISOTH

      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION EMU, RHO, RPM, PS, PA, PC, CD,
     +                 DORIF, LOSXSI, ALPHA, MFACTOR,MX,MY,
     +                 FFACTOR,FX,FY,TO,TOR,LOAD,ANGLE,TORQ,
     +                 PREC,QREC,TREC,QIN,QOUT,QFACTOR,EXO,EYO,
     +                 PMIN, PMAX, TKMAX, TKMIN,ALFT,UC,TC,Ec,
     +                 RSINPK, RCOSPK, IPAD,
     +                 PRELOAD, OFFSET,ROTDEL, KROTPAD, CROTPAD


      DOUBLE PRECISION FACTOR, Ecc, Qkgsec, Qkgmin, Qave,
     +                 TEA, TMAX,TMIN, PADMOMENT
      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,ISOTH,I,
     +        DEVICE,ICHECK, ISYM, ICSTEP, Iwrite, IU,BEARING
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        IDUM, TILT

C .....................................................................
C find MIN and MAX Pressures

      CALL MINMAXP    !=> on optionst.f
c.....................................................!
C find MIN and MAX Temperatures

c    !..........................! Average Exit Temperature
      IF (ISOTH.EQ.1) THEN
c    !..........................! barotropic model.
        TKMAX=1.0D0
        TKMIN=1.0D0
c    !..........................!
      ELSE
c    !..........................! from THERMAL analysis
        CALL MINMAXT(TMIN,TMAX,TEA)      !=> on optionst.f
        TKMAX=DMAX1(TKMAX,TMAX)
        TKMIN=DMIN1(TKMIN,TMIN)
        TMAX=TMAX*TC            ! MAX, MIN on lands
        TMIN=TMIN*TC            ! Average Temp. at exit
        TEA=TEA*TC              ! planes
c    !..........................!
      END IF
c    !..........................!

C .....................................................................

      Ecc=DSQRT(EXO*EXO+EYO*EYO)

      IF (IFULL.EQ.0) GOTO 77

C   !....................................! FIND ANGLE FOR FLUID FORCES
      LOAD=DSQRT(FX*FX+FY*FY)
      IF (FX.NE.(0.0D0)) THEN

          ANGLE=DATAN(DABS(FY/FX))*180.0D0/DACOS(-1.0D0)
          IF (FX.LT.(0.0D0)) THEN
              IF (FY.LT.(0.0D0)) THEN
                  ANGLE=ANGLE
              ELSE
                  ANGLE=-ANGLE
              END IF
          ELSE
              IF (FY.LT.(0.0D0)) THEN
                  ANGLE=-ANGLE+180.0D0
              ELSE
                  ANGLE=ANGLE+180.0D0
              END IF
          END IF


      ELSE

          IF (FY.GE.(0.0D0)) THEN
              ANGLE=-90.0D0
          ELSE
              ANGLE=90.0D0
          END IF

      END IF

C.......................................!.! Qout=Qleft+Qright
 77   Qkgsec=QIN*QFACTOR ! seal, SYMM HJB or pad bearing

      IF ((BEARING.EQ.1).AND.(ISYM.eq.0)) Qkgsec=QOUT*QFACTOR

      IF (MODEL.eq.1) Qkgsec=2.0D0*Qkgsec  ! 2row, HJB

      Qkgmin=Qkgsec*60.0D0

      Iwrite=1
      IU=6
C.......................................!.!Qright=Qin,Qleft=Qo-Qright.
  1   CONTINUE
c   !...................................!..........................
      IF (BEARING.EQ.3) THEN            ! for journal pad bearings
c   !...................................!..........................
      WRITE (IU, 501) Qkgsec, Qkgmin    ! side flow
      IF (IFULL.NE.1) WRITE (IU, 502) QOUT*QFACTOR ! leading edge flow
c   !...................................!..........................
      ELSE                              ! for seals/HJBs
c   !...................................!..........................
      WRITE (IU, 503) Qkgsec, Qkgmin    ! Flow at Pressure side
c   !...................................!..........................
      END IF                            !
c   !...................................!..........................

c     print side flows (left and right) for asymmetric HJBs/pad bearings
c     Qtotal=> Qout; QR => Qin,   QL=(Qtotal-QR)
      IF ((ISYM.eq.0).AND.(MODEL.EQ.2).AND.(BEARING.EQ.1)) THEN
        WRITE (IU, 550) (QOUT-QIN)*QFACTOR, QIN*QFACTOR
      END IF
c   !...................................!..........................
      WRITE (IU, 600)
      WRITE (IU, 100) FX, FY
      IF (IFULL.EQ.1) WRITE (IU, 101) LOAD, ANGLE
      IF (MODEL.eq.2) WRITE (IU, 105) MX, MY

c.......................................! Tilt Pad moment & rotation
      IF (TILT.EQ.1) THEN
         PADMOMENT=RSINPK*FX-RCOSPK*FY-KROTPAD*ROTDEL
         WRITE(IU,600)
         WRITE(IU,125) PADMOMENT, ROTDEL
      END IF
c..............................................................!
      WRITE (IU, 600)
      WRITE (IU,115) PMAX/1.0D6, PMIN/1.0D6
      IF (ISOTH.NE.1) THEN
       WRITE (IU, 120) TMAX, TMIN, TEA
      END IF
      IF (ICHECK.EQ.1) WRITE (IU, 110) TO*TOR
      WRITE (IU, 600)
c .......................................................................
      IF (Iwrite.eq.2) RETURN
      IF (DEVICE.eq.1) THEN
          Iwrite=2
          IU=1
          GOTO 1
      END IF
C.......................................!.............................
  120 FORMAT(' |    Temperatures MAX:', E12.5E2,' MIN:',E12.5E2,'K',
     +       2X,'Exit AVE:',E12.5E2,'K',2X,'|')

  110 FORMAT (' |    Torque on Film Lands=', E12.5E2, ' N-m', 36X, '|')
  115 FORMAT (' |    Pressure at film lands, MAX', E12.5E2, ' MPa',
     +        2X,'MIN=', E12.5E2, ' MPa', 8X, '|')
  125 FORMAT (' |    Tilt Pad Moment:',E12.5E2,'(Nm),  Pad angle ',
     +       ('rotation:',E11.4E2,'(rad)',3X,'|'))
  200 FORMAT (' ', '+', 77('-'), '+')

  100 FORMAT (' |    FX=', E12.5E2, ' N     FY=', E12.5E2, ' N',
     +        34X, '|')
  101 FORMAT (' |    Load=', E12.5E2,' N   Load Angle = ',
     +        F8.3, ' deg', 26X, '|')
  105 FORMAT (' |    MX=', E12.5E2, ' Nm    MY=', E12.5E2, ' Nm',
     +        33X, '|')

  400 FORMAT (' ', '|     MAX Rise in Pressure within recess: ',
     +        E12.5E2, 24X, '|')

  501 FORMAT (' |    Side Mass  Flow=', E12.5E2,' Kg/s = ',E12.5E2,
     +        ' Kg/min', 18X, '|')
  502 FORMAT (' |    Lead edge  Flow=', E12.5E2,' Kg/s   ',
     +                   37X, '|')
  503 FORMAT (' |    Mass  Flow=', E12.5E2, ' Kg/s = ', E12.5E2,
     +        ' Kg/min', 23X, '|')

  550 FORMAT (' |    (Left Flow)=', E12.5E2, ' Kg/s  (Right Flow)=',
     +        E12.5E2, ' Kg/s',12X,'|')
  555 FORMAT (' |    Mass Flow/Rhos=',E12.5E2, ' lt/min', 38X, ' |')

  600 FORMAT (' +', 77('-'), '+')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine PRINTCOEF                                                   **
C **                                                                         **
C **  PRINTCOEF  : Prints Dynamic force coefficientsin final format.         **
C **                                                                         **
C *****************************************************************************
          SUBROUTINE PRINTCOEF(Device)

      IMPLICIT NONE

      INCLUDE 'params.f'
c.........................................................................
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /PMINMAX/ PMIN, PMAX
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE
      COMMON /RECPAR/ HREC, VSUP, BETA
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD
      COMMON /INERC/ MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD
      COMMON /STIFA/ KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA
      COMMON /DAMPA/ CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA
      COMMON /INERA/ MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA

      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /FLMOM/ IMOMFLAG

      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD,
     +                 MXXD,MYYD,MXYD,MYXD,MmXXD,MmYYD,MmXYD,MmYXD,
     +                 KXXA,KYYA,KXYA,KYXA,KmXXA,KmYYA,KmXYA,KmYXA,
     +                 CXXA,CYYA,CXYA,CYXA,CmXXA,CmYYA,CmXYA,CmYXA,
     +                 MXXA,MYYA,MXYA,MYXA,MmXXA,MmYYA,MmXYA,MmYXA
      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD,
     +                 PRELOAD, OFFSET,ROTDEL
      DOUBLE PRECISION FREQU, SIGMA, L1,RES, HREC, VSUP, BETA,
     +                EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                AXO,AYO,ZO, PMIN,PMAX

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL, IMOMFLAG,
     +        DEVICE, ICASE, NCASE, INERFREQ, ISYM, ICSTEP, TILT,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        NPAD, NREC
c.........................................................................

      DOUBLE PRECISION omega,freq,Keq,I1,I2,I4,phi02,dumy,whirl,PI,
     +                 whirl1,whirl2,L,Lold,EPS
      INTEGER UN, Iset, Iwfr,I, INERLDUM
c.........................................................................
c calculation of whirl frequency ratio for freq=omega
c.........................................................................
      Un=6                             ! OUTPUT Unit
      PI=DACOS(-1.0D00)
      EPS=1.0D-05
      INERLDUM=0    !INERL
      Iwfr=0
      Iset=1

c    !.................................!.........
      IF (PMAX.EQ.PC) THEN             ! PAD IS UNLOADED
        WRITE(6,121)                   ! All coefficients = 0
        IF (DEVICE.EQ.1) WRITE(1,121)
        RETURN
      END IF
c    !.................................!.......
      IF (TILT.EQ.1) THEN
c    !.................................!.......
      INERFREQ=0                       !
      GOTO 101                         !
c    !.................................!.......
      END IF
c    !.................................!.......

c............................................................................
c     FOR 360 deg bearings here we calculate the whirl frequency ratio
c............................................................................

c    !....................................!
      IF (NPAD.GT.1)           GOTO 100   !=> NO WFR calculations
c    !....................................!

c    !....................................!
      IF (RPM.EQ.(0.0D0)) GOTO 100        !=> NO WFR caclulations
c    !....................................!

      omega=rpm/60.0D0
      freq=frequ*2.0D0*PI

      IF (DABS(frequ-omega).le.(1.0D0)) THEN
c    !....................................!
        Iwfr=1
        Keq=(Kxxd*Cyyd+Cxxd*Kyyd-Cyxd*Kxyd-Cxyd*Kyxd)/(Cxxd+Cyyd)

        dumy=((Keq-Kxxd)*(Keq-Kyyd)-Kxyd*Kyxd)
        dumy=dumy/(Cxxd*Cyyd-Cxyd*Cyxd)
        phi02=dumy/freq/freq


c      !.....................!
        IF (INERLDUM.EQ.1) THEN
c      !.....................!
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
            whirl=0.0D0 ! -DSQRT(-L)
          ELSE
            whirl=DSQRT(L)
          END IF


c      !.....................!
        ELSE                 ! INERL=0
c      !.....................!

          IF (phi02.lt.0.0D0) THEN
             whirl=-0.0D0  ! DSQRT(-Phi02)
          ELSE
             whirl=+DSQRT(+Phi02)
          END IF

c      !.....................!
        END IF
c      !.....................!


      END IF      ! cpm/60 = freq => rotordynamic coefficients
c ...............................................................!


 100  INERFREQ=0
      IF ((INERL.EQ.1).AND.(FREQU.NE.(0.0))) INERFREQ=1
c.............................................!..............................!
      INERFREQ=0                              ! 9/21/93 Calculate coeffs at w always
c.............................................!..............................!
c INERFREQ=0 always, Then force coefficients at frequency w
c Kij & Cij represent the real part and imaginary part of
c           the complex dynamic force impedance.
c     .............................!.............................!

 101  WRITE (Un,1800) frequ,Res
c     ............................. FORCE COEFFICs. due to DISPLACEMENTS
      WRITE (Un, 1900) KXXD, KYXD, KYYD, KXYD
      IF (FREQU.GT.(0.0)) WRITE (Un, 2000) CXXD, CYXD, CYYD, CXYD
      IF (INERFREQ.EQ.1) WRITE (Un, 2100) MXXD, MYXD, MYYD, MXYD
      IF (Iwfr.eq.1) WRITE (Un, 1999) KEQ, WHIRL
      WRITE (Un, 2200)
c     ................................................................

      IF (TILT.EQ.1) THEN             ! write pad-rotation coeffs.
        WRITE (Un,3900) ROTDEL, KdXk,KdYk,KXdk,KYdk,Kddk
        IF (FREQU.GT.0.0D0) WRITE (Un,3910) CdXk,CdYk,CXdk,CYdk,Cddk
        WRITE (Un, 2200)
        GOTO 34
      END IF

c     ................................................................
      IF (MODEL.eq.1) GOTO 34
c     ................................................................

c     ............................. MOMENT COEFFICs. due to DISPLACEMENTS
 11   IF ( (ISYM.EQ.1).AND.(ZO.EQ.0.0D0) ) GOTO 33
      WRITE (Un, 1901) KmXXD, KmYXD, KmYYD, KmXYD
      IF (FREQU.GT.(0.0)) WRITE (Un, 2001) CmXXD, CmYXD, CmYYD, CmXYD
      IF (INERFREQ.EQ.1) WRITE (Un, 2101) MmXXD, MmYXD, MmYYD, MmXYD
      WRITE (Un, 2200)

c     ................................................................
      IF (IMOMFLAG.EQ.0) GOTO 34
c     ............................. FORCE COEFFICs. due to ANGLE DISP
 22   WRITE (Un, 1902) KXXA, KYXA, KYYA, KXYA
      IF (FREQU.GT.(0.0)) WRITE (Un, 2002) CXXA, CYXA, CYYA, CXYA
      IF (INERFREQ.EQ.1) WRITE (Un, 2102) MXXA, MYXA, MYYA, MXYA
      WRITE (Un, 2200)

c     ................................................................
 33   IF (IMOMFLAG.EQ.0) GOTO 34
c     ............................. MOMENT COEFFICs. due to ANGLE DISP
      WRITE (Un, 1903) KmXXA, KmYXA, KmYYA, KmXYA
      IF (FREQU.GT.(0.0)) WRITE (Un, 2003) CmXXA, CmYXA, CmYYA, CmXYA
      IF (INERFREQ.EQ.1) WRITE (Un, 2103) MmXXA, MmYXA, MmYYA, MmXYA
      WRITE (Un, 2200)

C ------------------------------------------------------------------------
 34   IF (Iset.eq.2) RETURN
      IF (DEVICE.eq.1) THEN            !==> DUMP FILE
          Iset=2
          UN=1
          GOTO 101		       !==> DUMP FILE
      END IF



C ------------------------------------------------------------------------
c.........................................................................

 1800 format (' ', 3X,'frequency(w):',E12.5E2,'Hz',
     + 2X,'Squeeze Reynolds#(Res):(p/u)wC**2 =', E11.4E2,
     + /, 1X, 79('.'))
 3900 FORMAT (' ',3X,'PAD Rotation Coefficients:', 3X,
     +        'ANGLE (rotation):',E10.3E2,'(Rads)',/,
     +        4X,'Kdx=', E10.3E2, 1X, 'KdY=', E10.3E2,
     +        1X,'Kxd=', E10.3E2, 1X, 'Kyd=', E10.3E2, '(N)',/,
     +        4X,'Kdd=', E10.3E2,'(N.m)')
 3910 FORMAT (4X,'Cdx=', E10.3E2, 1X, 'CdY=', E10.3E2,
     +        1X,'Cxd=', E10.3E2, 1X, 'Cyd=', E10.3E2, '(N.s)',/,
     +        4X,'Cdd=', E10.3E2,'(N.m.s)')

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

 121  FORMAT (' ',3X,'PAD IS UNLOADED, PMAX=PC, FORCE COEFFICIENT',
     +               'S ARE NULL',/,1X,79('_'))

 2200 FORMAT (' ',  79('_'))

      END

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Last revised 7/3/95 by Dr. Luis San Andres
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::