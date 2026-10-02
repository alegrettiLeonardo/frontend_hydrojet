C
C SUB  calcpadt.f  4/6/94
C I have merged original SUBS calcpadt.f and eccpadt.f of hydrosealt
c into one single routine which is more efficient
c
c  ####     ##    #        ####   #####     ##    #####    #####          ######
c #    #   #  #   #       #    #  #    #   #  #   #    #     #            #
c #       ######  #       #       #####   ######  #    #     #     ###    ###
c #    #  #    #  #       #    #  #       #    #  #    #     #     ###    #
c  ####   #    #  ######   ####   #       #    #  #####      #     ###    #
c
c hydroflex.f    Copyright Luis SanAndres/TexasA&MUniversity/1994

c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR II
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C
C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcpadt                                                    **
C **                                                                         **
C **  CALCPADT: Thermohydrodynamic calculation of pad bearings               **
C **            Called in main routine HYDROFLEX.F at options 31 and 33      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCPADT(RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /PARAM2/ EXO, EYO
      COMMON /SOURCEA/ PRATIO,CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /PADK/ KPAD
      COMMON /GUESPI/ IGUESP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL, P,
     +                 EXO,EYO,EWX,EWY,EWEAR,BETAW,
     +                 PRATIO,CORIFS,SMASS,MPEPS,PREPS,MMP,SFLOW
      DOUBLE PRECISION TSUMP,TSUMPd,RHOSU,RHOSd
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        NPAD,NREC, ISYM, ICSTEP, KPAD ,IWEAR, IGUESP, BEARING

C................................................................
C     Local Variables
C................................................................
C
      DOUBLE PRECISION DUMYR,MAXIN, Tdiff, DTconv,RPM,PS,PA
      INTEGER JPAD, J,Kstart, IGCONV,KSKIP(MAXNPAD) ,ISYMold, ITER,
     +        Icenter, IGUESPold

C ----------------------------------------------------------------------------
C --  CALCPADT code                                                          --
C ----------------------------------------------------------------------------
c     ........................................................................
c      IGUESP=0  Use current P,U,V&T fields as initial guess for soln.
c                These fields are read from TEMP files (already existing)
c                IGUESP=0 from OPTION (33)
c                IGUESP=0 or 1 from choices in OPTION (31)
c      IGUESP=1  Guess P,U,V&T fields from guespp.f to initiate soln.
c                since there may not be TEMP data files with guess fields
c                IGUESP=0 from OPTION (33)
c                IGUESP=0 or 1 from choices in OPTION (31)
c     ........................................................................
c                                     !
      DTconv= 0.5D0/TSUMP             ! 0.5deg K for convergence between
c                                     ! two consecutive PAD edge temperatures
c     ................................! ........................................
      IGCONV=0                        ! Groove temperature convergence parameter
      IGUESPold=IGUESP                ! ........................................
      ISYMold=ISYM                    ! save current status on bearing symmetry

c     ................................! RESET indicators on pad calculation
      DO J=1,NPAD                     !
         KSKIP(J)=0                   ! =1: Skip solution of pad (SOLVE)
      END DO                          !
C    !................................!
      Icenter=1                       ! IF EXO=EYO=0 , THEN Icenter=0 and most
      IF ((EXO.EQ.0.0D0).AND.(EYO.EQ.0.0D0)) Icenter=0
c    !................................! likely there are no TEMP files available



      ITER=0                          !=============================!
c  !==================================! START THERMAL LOOP
  100 KPAD=0			      !=============================!
      JPAD=NPAD                       ! Upstream pad when Kpad=1
c                                     !.............................!
C  !..................................!
  200 KPAD=KPAD+1                     !
C  !..................................! FOR KPAD=1,2,...., NPAD
c                                     !.............................!

cNOTE!. . . . . . .. . . . . . . . . . . .! Current pad has leading
      IF (KSKIP(KPAD).EQ.KPAD) THEN       ! edge temp. determined
         IGCONV=IGCONV+1                  ! from within pad flow
         GOTO 600                         ! Thus, there is no need
      END IF                              ! to recalculate flow field.
c    !. . . . . . .. . . . . . . . . . . .!............................!

      CALL XYDATA(KPAD,RPM,PS,PA)               ! generates mesh and sets current ISYM
c                                     ! ==> calcmesht.f
c                                     ! ...............................!
      IF (IGUESP.EQ.0) THEN           ! READ current flow fields for
        CALL READTEMP(KPAD)           !  KPAD saved on TEMP file
      END IF                          !   (returns with IGUESP=0)
c                                     ! ==> readwritet.f
c                                     ! ...............................!
      IF ((ISYMold.eq.1).AND.         ! IF BEARING Modified on ASKITER then
     +    (ISYM.eq.0)) THEN           ! Set UVPT fields mirror image to
 	      CALL TEMPASYM               ! start solution for
      END IF                          ! asymmetric bearings
c     ISYM=0 (ASYMMETRIC), ISYM=1 (SYMMETRIC) geometry
c     ==> TEMPASYM resides on calcsolnt.f
c                                     !
      Kstart=ISYM+(ISYM-1)*NYI        !..........................!
c                                     !
      CALL FILMH(EXO, EYO, NXT, NYI)  ! generate film thickness
c                                     ! ==> filmwt.f
      IF (IGUESP.EQ.1) THEN           !..........................!
          CALL SETBC(RPM,PS,PA)                  ! Set known Boundary Conds.
          CALL GUESP(0,RPM,PS,PA)               ! and guess a Pfield
      END IF                          !..........................!
C     ==>                         SETBC: calcsolnt.f, GUESP: guespp.f

C.....................................!..................................!
C A concentric JB w/o preload can not be calculated
C with present model, i.e. since there is no axial flow to carry away the
C heat then temperature should keep increasing w/o bounds. Theoretically
c this is correct. However, a thermal wedge is created which forces the
c hot fluid out of the bearing. Fresh fluid then replenishes the bearing.
c This phenomena can not be modeled here. Fortunately, a concentric JB
c is a rarity since it does not have load and shows very poor force coeffs.
c besides being unstable at all speeds. (LSA 11/8/93)
C.....................................!..................................!

  300 CONTINUE

C    !................................!
      IF ((BEARING.EQ.3).AND.(Icenter.EQ.0)) THEN
C    !................................! For Concentric JB, no axial flow
         DO J=Kstart,NYI              ! causes infinite temperature rise !
            TLEAD(KPAD,J)=TSUMPd      ! Use approximate formula==>
         END DO                       ! Note: Results are meaningless
C    !................................!
      END IF                          !
C    !................................!

C     KSKIP=0 , inlet pad temperatures depend on Tupstream & Tgroove
C     KSKIP>0 , inlet pad temperatures depend on Tdownstream (pad)
C    !.....................................!
C##   IF (KSKIP(KPAD).EQ.0) THEN           !
C    !.....................................!
C                                          ! SOLVE: U,V,P,T on KPAD
        CALL SOLVE(0,RPM,PS,PA)                      !
C                                          ! ==> calcsolnt.f
C    !.....................................!
C##   END IF				   !
C    !.....................................!

      IF ((BEARING.EQ.3).AND.(Icenter.EQ.0)) GOTO 400
C     FOR centered JB do not update pad inlet temps.

C    !.....................................!........................!
C     IF Flow at leading edge of pad is POSITIVE, then fluid enters
c     pad. Upstream of edge (at groove separating pads),
c     thermal fluid mixing occurs between hot fluid
c     leaving upstream pad (at trailing egde) and cold fluid feed
c     at groove. Here, the entrance (inlet) temperature to pad
c     is calculated (updated) using SUB GROOVETEMP

c    !................................!.............................!
      IF (QLEAD(KPAD).GT.0.0D0) THEN  ! Update groove edge temps.
c    !................................!.............................!

        CALL GROOVETEMP(JPAD,KPAD,Kstart,Tdiff)   !==> calcthert.f
c
c      CHECK for temperature convergence on pad
c      Tdiff=|TLEADnew - TLEADold|
c      here we choose convergence to be 0.5 deg K ~ 1 deg F
c      since most thermometers are that accurate

        IF (Tdiff.LE.DTconv) IGCONV=IGCONV+1
c                                     !
        KSKIP(KPAD)=0                 !
        GOTO 500                      !

c    !................................!.............................!
      END IF                          !
c    !................................!.............................!


c     IF QLEAD<=0 then flow leaves pad and leading edge temperature
c     is obtained from the solution within the film lands by
c     upwinding to the land internal points

 400    KSKIP(KPAD)=KPAD              ! Thus, Skip call of SOLVE
        IGCONV=IGCONV+1               ! and No update of groove temps.
c                                     ! i.e. convergence achieved on KPAD

c    !................................!.................................!
 500  IGUESP=0                        ! Redefine IGUESP for next iteration
C    !................................!
      IF (NPAD.GT.1) THEN             ! for multiple pad bearing
        CALL WRITETEMP(KPAD)          ! store fields in temporal file
      END IF                          !  including IGUESP
C    !................................!

C    !================================!
 600  IF (KPAD.EQ.NPAD) THEN
C    !================================!
          ISYMold=ISYM                ! RESET value of SYM parameter
          ITER=ITER+1                 ! ---> get ready for new iteration
          WRITE (6,177) ITER          ! write thermal loop iteration #

C       !.............................!.................................!
          IF (IGCONV.LT.NPAD) THEN    ! Groove temps. DONOT converge
C       !.............................!.................................!
           IGCONV=0                   ! Ready for next KPAD do-loop
c                                     !
           IF (NPAD.EQ.1) THEN        ! Single pad groove bearing
             GOTO 300                 !
           ELSE                       ! Multipad bearing
             GOTO 100                 ! New iteration loop with
           END IF                     ! updated groove temperatures
C       !.............................!.................................!
        ELSE                          !
C       !.............................!.................................!
c                                     !
           IF (NPAD.EQ.1) THEN        ! for Single pad bearing
              CALL WRITETEMP(KPAD)    ! store fields in TEMP1 file
           END IF                     !.......................

           RETURN                     ! Groove temps CONVERGE
C       !.............................! ==> goto MAIN program
        END IF                        !
C       !.............................!.................................!

C    !================================!
      ELSE			      ! KPAD=1,2,... NPAD-1
C    !================================!

        JPAD=KPAD                     ! Ready for next KPAD

        IF (ITER.EQ.0.AND.IGUESPold.EQ.1) THEN   ! for centered bearing
               IGUESP=1                          ! from opt 31 since there
        END IF                                   ! are no TEMP data files
c                                                ! on initial iteration

        GOTO 200                      ! => calculate solution
C                                     !    for next pad
C    !================================!
      END IF                          !
C    !================================!

  177 FORMAT ('$ THERMAL LOOP, ITERATION #', I4, ' COMPLETED')

      END

C ----------------------------------------------------------------------------
C --  created 4/6/94 by Dr. Luis San Andres                                --
C I have merged SUBS CALCPADT and ECCPADT into a single one : CALCPADT
C ----------------------------------------------------------------------------