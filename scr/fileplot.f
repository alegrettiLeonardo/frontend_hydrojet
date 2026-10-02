C......................................................
C UPDATED for recess/pressure with jet: 7/6/95
C......................................................
C Updated for CALL to FILMC where compliance is updated.
C ......................................................
C Need to OPEN a file to save deformed bearing surface
C ......................................................
C
C
C ######     #    #       ######  #####   #        ####    #####          ######
C #          #    #       #       #    #  #       #    #     #            #
C #####      #    #       #####   #    #  #       #    #     #            #####
C #          #    #       #       #####   #       #    #     #     ###    #
C #          #    #       #       #       #       #    #     #     ###    #
C #          #    ######  ######  #       ######   ####      #     ###    #

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
C **  Subroutine Fileplot                                                    **
C **                                                                         **
C **  SAVES pressure,temperature and film thickness on data files for plot   **
C **                                                                         **
C *****************************************************************************
C Assumed PGROOVE: Pressure between pads is equal to zero.


      SUBROUTINE Fileplot(RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /XYVEC/   XP(MAXNXT), XU(MAXNXT),
     +                 YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /HFILM/   HP(MAXNXT,-MAXNYI:MAXNYI),
     +                 HU(MAXNXT,-MAXNYI:MAXNYI),
     +                 HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU, BETU, ALFP
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,
     +                 DEN12P12, VIS12P12, P2
c     ................................................................
      COMMON /PARAM2/ EXO, EYO
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /TMAXMIN/ TKMAX,TKMIN
      COMMON /THERMID/ ISOTH
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
c     ................................................................
      DOUBLE PRECISION XP,XU,YP,YV,HP,HU,HV,P,T,PREC,TREC,QREC,
     +                 PRECdo,PRECup,PRjet
      DOUBLE PRECISION EXO, EYO,EWX,EWY,EWEAR,BETAW, QIN,QOUT,
     +                 QFACTOR,REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,
     +                 BETU, ALFP, TEMPK, VSOUND,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,
     +                 DEN12P12, VIS12P12, P2 ,TKMAX,TKMIN,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HREC

      INTEGER NPAD, NREC, NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI,
     +        NYI,NXT,ISYM,ICSTEP,IFULL,IWEAR,BEARING,ISOTH ,IF,IL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER I,J,II,JJ,N,K,COUNTER,KP,KK,IOS,idumy,jdumy,Istart,KKK,
     +        Jstart
      DOUBLE PRECISION ZERO,HGROV, HPOCK, PGROOVE, Pres , FACTOR,
     +       DTMAX, MAXIN, DT(MAXNXT,-MAXNYI:MAXNYI),RPM,PS,PA

C ----------------------------------------------------------------------------
C --  Fileplot code                                                        --
C ----------------------------------------------------------------------------
C
C Structure of PLOTPRES file is as follows
C ----------------------------------------
C line # 1: N1, N2, PSA, PA
C following lines:
C   { [(X, Y, P)i,j],  j=1, N2 } i=1, N1
C
C where N1: # points in circumferential direction (X)
C       N2: # points in axial direction (Y)
C       Xij  array of points in circ. direction
C       Yij  array of points in axial direction
C       Pij  dimensionless pressure defined as P=(Pres-PA)/PSA where
C            PA and PSA are in MPa (Mega-pascals)
C
C Groove pressures are assumed as zero
C
C ...........................................................................
C Structure for PLOTFILM and PLOTMESH files is identical
C ...........................................................................
C Structure of PLOTTEMP file is as follows
C ----------------------------------------
C line # 1: N1, N2, Tsupply
C following lines:
C   { [(X, Y, T)i,j],  j=1, N2 } i=1, N1
C
C where N1: # points in circumferential direction (X)
C       N2: # points in axial direction (Y)
C       Xij  array of points in circ. direction
C       Yij  array of points in axial direction
C       Tij  dimensionless temperature defined as
C            T= FACTOR*(T-Tsupply)/Tsupply
C            where Tsupply [degK] & FACTOR = 10 (arbitrary)  Tempk=Tsupply
C
C ...........................................................................

      FACTOR=10.0D0      ! dumy for temperature diffference amplification
C ...........................................................................

      ZERO=0.0D0
      CALL CPARAM(RPM,PS,PA)

C # of points in axial direction is equal to :

      IF (BEARING.EQ.1) THEN           !=> HJB
         Jdumy=2*Nyi+1+2               !
         CALL PJET  ! generate Precess !........ ##7/7/95
      ELSE IF (BEARING.EQ.2) THEN      !=> SEAL
         Jdumy=Nyi+1                   !
      ELSE IF (BEARING.EQ.3) THEN      !=>BEARING
         Jdumy=2*Nyi                   !
      END IF                           !........................!

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

  412 IF (KP.EQ.NPAD) THEN
          GOTO 477      ! => start transfer
      ELSE
         Idumy=Idumy+2  ! 2 points for groove pressures L&R
      END IF
c    !...............!
       GOTO 400      ! -------> NEXT PAD KP=1,2,3,... NPAD
c    !...............!

c ----------------------------------------------------------------c
 477  CONTINUE  ! Transfer P & H data
c ----------------------------------------------------------------c
      WRITE (6,23)
  23  FORMAT (' ','  *transfer  pressure data to PLOTPRES File',/,
     +        ' ','  *  ..   temperature data to PLOTTEMP File',/,
     +        ' ','  *and film thickness data to PLOTFILM File')

      OPEN (UNIT=2, FILE='PLOTPRES', STATUS='UNKNOWN', IOSTAT=IOS,
     +       ERR=1000)
      WRITE ( 2, *, IOSTAT=IOS, ERR=1200) Idumy, Jdumy , PSA, PA

      OPEN (UNIT=3, FILE='PLOTFILM', STATUS='UNKNOWN', IOSTAT=IOS,
     +       ERR=1000)
      WRITE ( 3, *, IOSTAT=IOS, ERR=1200) Idumy, Jdumy, CLEAR

      OPEN (UNIT=4, FILE='PLOTTEMP', STATUS='UNKNOWN', IOSTAT=IOS,
     +       ERR=1000)
      WRITE (4,*,IOSTAT=IOS,ERR=1200)Idumy, Jdumy, Tempk, Factor

c ----------------------------------------------------------------c

      PGROOVE=ZERO     ! assumed pressure between pads.
      KP=1

c    !...........................!

  300 CALL XYDATA(KP,RPM,PS,PA)            !-> builds X&Y mesh arrays
      CALL READTEMP(KP)          !-> reads UVP values
      IF (IWEAR.eq.1) THEN       !-> builds film thickness vectors
         CALL HWEAR(EXO, EYO, NXT, NYI)
      ELSE
         CALL FILMH(EXO, EYO, NXT, NYI)
      END IF

c###----------------------------! update for compliance coef.
      CALL FILMC

c####---------------------------! update for compliance coef.

      Jstart=ISYM+(ISYM-1)*NYI   !
         DO J=Jstart,NYI,1       ! Relative temperature rise
            DO I=1,NXT           !
               DT(I,J)=(T(I,J)-1.0D0)*FACTOR
            END DO               !
         END DO                  !

c    !...........................!
      IF (ISYM.EQ.1) THEN        !-> image data for pad left side
c    !...........................!   on symmetric bearing
         DO I=1, NXT
            DO J=1, NYI
              P(I,-J)=P(I,J)
              DT(I,-J)=DT(I,J)
              HP(I,-J)=HP(I,J)
            END DO
              P(I,0)=P(I,1)
              DT(I,0)=DT(I,1)
              HP(I,0)=HP(I,1)
         END DO


c    !...........................! ISYM = 1
      END IF
c    !...........................!


c    !...........................!
  500 CONTINUE                   !-> STORE MESH DATA in PLOTMESH File
c    !...........................!

      IF (BEARING.EQ.2) GOTO 666 ! ==> FOR ANNULAR SEAL


c    !...........................! left side of PAD
      IF (KP.GT.1) THEN
c    !...........................! left side of PAD
          I=1
          do j=-Nyi, -Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), PGROOVE
           write(3,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HP(I,J)
           write(4,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), DT(I,J)
          end do
      IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), PGROOVE
           write(3,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HP(I,J)
           write(4,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), DT(I,J)
          end do
      END IF
          do j=Npa, Nyi, 1
           write(2,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), PGROOVE
           write(3,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), HP(I,J)
           write(4,*, IOSTAT=IOS, ERR=1200) XP(I),YP(j), DT(I,J)
          end do
      END IF      ! at leading edge of pad K+1: Groove P
c    !...........................! left side of PAD

      IF (NPOCKET.EQ.0) GOTO 610

c    !---------------!
  606 DO K=1, NPOCKET
c    !---------------!
c !->   on land between rec/grooves ....................................!
        DO II=1, NLC
          I=(K-1)*Nxi+II
          do j=-Nyi, -Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

          do j=-Npa, Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

          do j=Npa, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

        END DO                    ! II=1, NLC

c  --  !-> on recesses & land above ....................................!
        KK=(K-1)*Nxi+Nlc-1
        DO II=1, NPC
          I=KK+II
          do j=-Nyi, -Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

          do j=-Npa, Npa, 1
            Pres=PRjet(K,II)     ! recess pressure field
            write(2,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),DMAX1(Pres,Pcav)
            write(3,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),HP(I,J)
            write(4,*,IOSTAT=IOS,ERR=1200) XP(I),YP(j),DT(I,J)
          end do

          do j=Npa, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

        END DO       ! II=1, NPC

c    !---------------!..............
      END DO         ! K=1, NPOCKET
c    !---------------!..............


c    .....................!
      IF (IFULL.EQ.0) THEN
c    .....................! for PAD BEARING: downstream of recess
        JJ=NPOCKET*NXI

        DO II=1, NLC
          i=JJ+II

          do j=-Nyi, -Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

          do j=-Npa, Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

          do j=Npa, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

        END DO                    ! II=1, NLC
c    ......................!
      END IF
c    ......................!

      GOTO 612

c    ............................!
  610 DO I=1, NXT                ! NO POCKET ON PAD:
c    ............................!
          do j=-Nyi, -Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

      IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do
      END IF

          do j=Npa, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),DMAX1(P(i,j),Pcav)
           write(3,*,IOSTAT=IOS,ERR=1200)XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200)XP(I),YP(j),DT(I,J)
          end do

       END DO                    ! II=1, NLC

       I=NXT
C     ................................. !
  612 CONTINUE

       IF (KP.EQ.NPAD) GOTO 777

c    .........................!
       IF (KP.LT.NPAD) THEN
c    .........................! at trailing edge of pad K
c      PGROOVE=P(I,Nyi)       !####

          do j=-Nyi, -Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),PGROOVE
           write(3,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200) XP(I),YP(j),DT(I,J)
          end do

      IF (BEARING.EQ.1) THEN
          do j=-Npa, Npa, 1
           write(2,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),PGROOVE
           write(3,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200) XP(I),YP(j),DT(I,J)
          end do
      END IF

          do j=Npa, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),PGROOVE
           write(3,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),HP(I,J)
           write(4,*,IOSTAT=IOS,ERR=1200) XP(I),YP(j),DT(I,J)
          end do
c    .........................!
       END IF                 ! at trailing edge of pad K
c    .........................!


       KP=KP+1
       GOTO 300

c    ----------------------!
  666 DO I=1, NXT          ! ANNULAR SEAL
c    ----------------------!
          do j=0, Nyi, 1
           write(2,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),P(i,j)
           write(3,*,IOSTAT=IOS,ERR=1200) XP(i),YP(j),HP(i,j)
           write(4,*,IOSTAT=IOS,ERR=1200) XP(I),YP(j),DT(I,J)
          end do
       END DO              ! II=1, NXT
c    ----------------------!

c    !...............!
 777   CONTINUE
       CLOSE (2)
       CLOSE (3)
       CLOSE (4)
       WRITE (6,25)
c    !...............! KKK=0

 25   FORMAT(' ',3X,'PRESSURE(X,Y) >> PLOTPRES File completed',/,
     +       ' ',3X,'Tfluid(X,Y)   >> PLOTTEMP File completed',/,
     +       ' ',3X,'FILM H(X,Y)   >> PLOTFILM File completed')


      RETURN

c    .........................!
 1000 WRITE (6, *) 'Error opening file:'
      WRITE (6, *) 'Aborting write.'
      CLOSE (UNIT=2)
      CLOSE (UNIT=3)
      CLOSE (UNIT=4)
      STOP

c    .........................!
 1100 WRITE (6, *) 'Error on file name:'
      WRITE (6, *) 'Aborting program'
      CLOSE (UNIT=2)
      CLOSE (UNIT=3)
      CLOSE (UNIT=4)
      STOP


c    .........................!
 1200 WRITE (6, *) 'Error during write:'
      WRITE (6, *) 'Aborting write.  The file is being deleted.'
      CLOSE (UNIT=2, STATUS='DELETE')
      CLOSE (UNIT=3, STATUS='DELETE')
      CLOSE (UNIT=4, STATUS='DELETE')
      STOP

c ......................................!.........................

      END

c======================================================================
c hydrojet.f program / Dr. Luis San Andres/ Texas A&M University 1995
c last updated 7/6/95
c======================================================================