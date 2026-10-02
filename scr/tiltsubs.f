
c  #####     #    #        #####   ####   #    #  #####    ####           ######
c    #       #    #          #    #       #    #  #    #  #               #
c    #       #    #          #     ####   #    #  #####    ####           #####
c    #       #    #          #         #  #    #  #    #       #   ###    #
c    #       #    #          #    #    #  #    #  #    #  #    #   ###    #
c    #       #    ######     #     ####    ####   #####    ####    ###    #

c hydrojet.f: copyright Dr. L. San Andres/ Texas A&M University 1993.
c


c
c 4/27/94 introduced outher thermal loop for thermal effects
c
c
c 2/15/94 correct convergence to loaded pad is attained IFF
c         the pad moment stiffness Kddk+Krotdel > 0.
c
c
c 1/6/94: need to change EMAX on LOAD sub since there may be cases
c         where e/c can be larger than (1-preload/C) !!!

c IMPORTANT: FOR tilt-pad bearings
c            stiffnesses vary rapidly, so small increments
c           in ecc. produce large changes in load. Thus, load
c           routine more than often overshoots to very large
c           eccentricities. Need to improve this.
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine TILTECC                                                     **
C **                                                                         **
C ** TILTECC: Given journal eccentricity find load components for            **
C **          tilt-pad bearing and such that Pad moments = zero              **
C **                                                                         **
C *****************************************************************************
C *****************************************************************************

      SUBROUTINE TILTECC(DEVICE,LOSTAT,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PARAM2/ EXO, EYO
      COMMON /FORCE0/ FFACTOR,FX,FY,TO,TOR
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd
      COMMON /FREQ/ FREQU,SIGMA,L1,RES,JCASE,NCASE
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /PMINMAX/ PMIN, PMAX

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD
      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/  KROTPAD, CROTPAD
      COMMON /RESFXY / FXTmom, FYTmom, DEX, DEY

      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /IDTHDPAD/ ITPAD
      COMMON /GUESPI/ IGUESP
      COMMON /SWITCH/ IPROP
      COMMON /PADK/ KPAD
      COMMON /BTYPE/ BEARING
c     .........................................................
      DOUBLE PRECISION QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL,
     +                 TSUMP,TSUMPd,RHOSU,RHOSd

      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,
     +                 EXO, EYO, FFACTOR, FX, FY, TO, TOR,
     +                 FXT,FYT,TOT,MXT,MYT,QINT,QOUTT,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 FREQU,SIGMA,L1,RES
      DOUBLE PRECISION ROTPAD, PRELOAD, OFFSET, ROTDEL,
     +                 KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD, KROTPAD, CROTPAD,
     +                 PMIN,PMAX,FXTmom, FYTmom, DEX, DEY
      DOUBLE PRECISION XP,XU,YP,YV

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        NPAD, NREC, IFULL, TILTPAD, TILT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        ITER, ITMAX, ITPMAX, SVERB,DVERB, BEEP,
     +        JCASE,NCASE, IWEAR, ISYM, ICSTEP,
     +        ITPAD, IGUESP, IPROP, KPAD, BEARING


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION MomK, MomKold, Momdiff, Rotold,ZERO, Delta,
     +                 Rdelmax, Deldiff, Rotoldd, FREQUB,
     +                 ROTtrail, ROTlead, Rdelta, Rdelinc,
     +                 Merror, Fnorm, deltaold,Sign, Mdif,
     +                 M, Xtrail, Xlead, Xpivot, KROTdel
      DOUBLE PRECISION Tdiff,DUMYR,MAXIN,DTconv

      INTEGER Count, CountMAX,CHO,DEVICE,DUMYI,IFILE, IOS,
     +        ILOAD, ILOADo, CountMAXX, LOSTAT, IL10, Imax, Icl,
     +        JPAD, ITERTH, IGUESPold, ISYMold, Kstart, IGconv,
     +        J, ITERTHmax, KSKIP(MAXNPAD), ITHpad(MAXNPAD),
     +        ITbal(MAXNPAD)

      CHARACTER FILE*80
      CHARACTER*1 YN

C ----------------------------------------------------------------------------
C --  TILTECC  code                                                         --
C ----------------------------------------------------------------------------
C     TILTECC: Only for SYMMETRIC/ALIGNED TILT PAD BEARINGS since we have not
c              included forces due to misalignment or bearing asymmetry.
c
c     KROTPAD and CROTPAD are the rotational stiffness (Nm/rad)
c                     and damping coefficient (Nm.sec/rad) for flexure pad
c
C ----------------------------------------------------------------------------
c
c    !......................................................!
c    ! LOSTAT=1 , TILTECC called from LOADTILT subroutine OPT(44)
c    ! LOSTAT=0 , TILTECC called from Main Program, OPT (33)
c    !......................................................!

      ZERO=0.0D0                       !
      ISYMold=ISYM                     ! save current status on bearing symmetry


c   !..................................!................................!
      IF (LOSTAT.EQ.0) THEN            ! Given eccentricity find solution
c   !..................................!.................. OPT 33    ...!
          YN='Y'
          WRITE (6, 11)
 11       FORMAT ('$', 'Do you wish to calculate dynamic ',
     +            'coefficients [Def: Y] (Y/N): ')
     	  READ (5, 12) YN
 12       FORMAT (1A)

          FREQU=RPM/60.0D0                    	! HZ, for synchronous whirl
          IF (RPM.EQ.0.0D0 ) FREQU=100.0D0      ! HZ, for HJBS, no rotation

          IF ((YN.eq.'N').OR.(YN.eq.'n')) THEN
             FREQU=0.0D0
          ELSE
           OPEN (UNIT=42, STATUS='SCRATCH',IOSTAT=IOS)
C          Unit 42 is used as a dumping place for translations between formats.
           WRITE (6, 6099) FREQU
           CALL ENTERVAL(FREQU)
           CLOSE (UNIT=42)
          END IF

          FREQUB=FREQU
          CALL ASKITER                 ! Set Ex, Ey
c   !..................................!................................!
      ELSE                             ! Given load find eccentricity
c   !..................................!................................!
      FXTmom=ZERO                      ! for load calculations
      FYTmom=ZERO                      ! FiTmom due to unbalanced pads
c   !..................................!................................!
      END IF                           !
c   !..................................!................................!

c     ........................................................................
c                                     !
      DTconv= 0.5D0/TSUMP             ! 0.5deg K for convergence between
c                                     ! two consecutive PAD edge temperatures
c                                     ! TSUMP=TSUPPLY
c     ................................! ........................................
      IGCONV=0                        ! Groove temperature convergence parameter
      IGUESPold=IGUESP                ! ........................................
c     ................................! RESET thermal convergence
      DO J=1,NPAD                     !       indicator
         KSKIP(J)=0                   ! =1: Skip solution of pad (SOLVE)
         ITHpad(J)=0                  !
      END DO                          !
C    !................................!

c   !... SET some Maximum values for # its. and error criteria .........!
c NOTE:        These % values should be declared on params.f for
c              improved I/O performance of program
c
      Rdelmax=0.15D0*CLEAR/(DIAM/2.0D0)! max pad rotation 15% clearance/step
      Momdiff=0.001D0                  ! MAX moment difference (N.m)
      Deldiff=0.02D0*Rdelmax           ! MAX 0.2% change in max pad rotation
      Merror=Momdiff                   !
      Imax=3                           ! MAX # of searches IF pad is unloaded
                                       !.........
      CountMAXX=15                     ! MAX#  of iterations for Pad
c                                      !       moment balance

c   !..................................!.Control statements ......!
      CountMAX =CountMAXX              !
      DUMYI=0                          ! Dumy Device
      NCASE=1                          ! FOR one recess depth only


c    !..................................!........................!
      IF ((EXO.EQ.ZERO).AND.(EYO.EQ.ZERO)) THEN
c    !..................................!........................!
        ROTDEL=Rdelmax                  ! set pad rotation to 15%
       IF (KROTPAD.EQ.ZERO) THEN        ! of pad clearance
        DO KPAD=1, NPAD                 ! as initial guess only
         ROTPAD(KPAD)=ROTDEL            !
        END DO                          !
       END IF                           !
c    !..................................!........................!
      ELSE                              !
c    !..................................!........................!
        Count=0                         ! check here if all pads
        DO KPAD=1, NPAD                 ! have some rotation.
         IF (ROTPAD(KPAD).EQ.ZERO) Count=Count+1
        END DO                          !........................!
        IF ((KROTPAD.EQ.ZERO).AND.      !
     +            (Count.EQ.NPAD)) THEN ! IF all pads have zero
          DO KPAD=1, NPAD               ! rotation then, give an
             ROTPAD(KPAD)=RDELMAX       ! initial (guessed) rotation
          END DO                        !
        END IF                          ! only for tilting pad JB
c    !..................................!........................!
      END IF
c    !..................................!........................!

C......................................................................!
C Input: Calculations performed at FREQU=zero whirl frequency to determine
C        static (equilibrium) eccentricity and pad rotations.
C......................................................................!


 33   FREQU=ZERO                       ! set whirl frequency = 0 HZ
C                                      !------------------------------
C--------------------------------------! START of ECCXY Calculations  !
C                                      !------------------------------
C                        !.............!...............!
      CALL ZERO0         ! ZERO static load and
      CALL ZERO1         !      global force coefficients for
      JCASE=2            ! Perturbations in X & Y directions
      IFILE=3            ! ............................!

c .........................................................................!
c The algorithm is based on increasing rotations until pad is loaded ILOAD=1
c ROTDEL is a pad angle value which loads the pad,
c Then, algorithm works until the pad moment is set to null or close to zero.
c .........................................................................!
c Based on FIXED eccentricities EX,EY, a N-R algorithm shows that the
c improved pad rotation angle is given by:
c
c RDel(new)=RDel(old) + [Momk(old)-Kp RDel(old)]/[Kp+Kdd(old)]
c
c where Kp=KROTPAD is the rotational stiffness of a flex tilt-pad and
c       Kdd is the rotational angular stiffness of the fluid film.
c
c .........................................................................!
c AT least two iterations are always required on the pads to check trends
c    on loading conditions
c .........................................................................!

      ITERTH=0           !
      ITERTHmax=10       ! Maximum # of iterations for thermal loop

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 77   CONTINUE           ! OUTER LOOP ON THERMAL EFFECTS
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

      JPAD=NPAD

      ITERTH=ITERTH+1

C   !--------------------!------------------------------------------!
      DO KPAD=1, NPAD    ! SWEEP OVER ALL PADS:
C   !--------------------!------------------------------------------!
          Count=0        ! reset all counters
          ILOAD=0
          IL10=0
          Icl=0

c      !........................................! message to USER
          IF (SVERB.GT.1) WRITE(6,110) KPAD
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) WRITE(1,110) KPAD
c      !........................................! message to USER

c      !........................................!...................!
          CALL XYDATA(KPAD,RPM,PS,PA)                     ! Generates KPAD Mesh
c      !........................................!...................!
          CALL READTEMP(KPAD)                   ! Reads input data
c      !........................................!...................!

c NOTE: ISYM=0 (ASYMMETRIC), ISYM=1 (SYMMETRIC) geometry
c                                     ! ...............................!
      IF ((ISYMold.eq.1).AND.         ! IF BEARING SYM Modified on ASKITER
     +    (ISYM.eq.0)) THEN           ! Then set UVPT fields mirror image to
 	  CALL TEMPASYM               ! start solution for
      END IF                          ! asymmetric bearings
c                                 ==>   TEMPASYM resides on calcsolnt.f
c                                     !
      Kstart=ISYM+(ISYM-1)*NYI        !..........................!
c
c     !........................................!...................!


c      !........................................!...................!
          IF (LOSTAT.EQ.1) THEN                 ! find new pad rotation
c        !......................................! OPTION (34)
            CALL FINDROT(DEX,DEY,ROTDEL)        ! for static load calculation
            ROTPAD(KPAD)=ROTDEL                 !
c   >>      FINDROT is in tiltsubs.f program and determines:
c...>>      ROTDEL(new)=ROTDEL(old)+(Mbalan-KdXk*DEX-KdYk*DEY)/(Kddk+Kpad)
c        !......................................!....................!
          ELSE                                  !
c        !......................................!....................!
            ROTDEL=ROTPAD(KPAD)                 ! initial rotation
c        !......................................!....................!
          END IF                                !
c      !........................................!....................!

c        !..........................................................!
c        ! find maximum allowed angles of rotation for current pad  !
c        ! ONLY valid for uniform clearance tilt-pads w/o misalignment
c        !..........................................................!
c
          M=PRELOAD/CLEAR                       ! pad preload dimensionless
c SET:                                          ! ..
          Xtrail=XP(NXT)                        ! leading and trailing edge
          Xlead=XP(1)                           ! pad angles [rad]
          Xpivot=OFFSET*(Xtrail-Xlead)+Xlead    !..........................!

c AT trailing end leading pad edges:

          ROTtrail=(1.0D0+EXO*DCOS(Xtrail)+EYO*DSIN(Xtrail) -
     +              M*DCOS(Xtrail-Xpivot))*CLEAR*2.0D0/DIAM/
     +              DSIN(Xtrail-Xpivot)
          ROTlead =(1.0D0+EXO*DCOS(Xlead )+EYO*DSIN(Xlead ) -
     +              M*DCOS(Xlead-Xpivot ))*CLEAR*2.0D0/DIAM/
     +              DSIN(Xlead-Xpivot)
c        !..........................................................!
c         ROTtrail(>0) and ROTlead(<0) are the rotations needed for the
c         KPAD to CONTACT the journal for the current values of eccentricity
c         at the trailing and leading edges of the pads, respectively.
c        !..........................................................!
c         PAD rotation ROTdel should be between ROTtrail and ROTlead
c
          Rdelta=(ROTtrail-ROTlead)/COUNTmax        ! rotation angle step
          Deldiff=0.0005*DABS(ROTtrail-ROTlead)     ! max allowed difference
          ILOADo=-1

c        !.....................................!....................!
c         ITERATIVE PROCESS TO STATICALLY BALANCE PAD starts here
c        !.....................................!....................!
 35       Count=Count+1
c        !.....................................!....................!
          IF (ROTDEL.LT.ROTlead) THEN          ! IF New Rotation causes
             ROTDEL=ROTlead +Rdelta            ! PAD to contact at leading
             ILOADo=ILOAD                      ! edge then reduce angle
          END IF                               !....................!
c                                              !
          IF (ROTDEL.GT.ROTtrail) THEN         ! IF New Rotation causes
             ROTDEL=ROTtrail-Rdelta            ! PAD to contact at trailing
             ILOADo=ILOAD                      ! edge then reduce angle.
          END IF                               !.....................!

c        !.....................................!....................!
c NOW     Build film thickness array, solve zeroth-order flow field and
c         calculate pad fluid film forces:
c
c     !........................................!...................!
          CALL FILMH(EXO, EYO, NXT, NYI)       ! ==> film.f
c     !........................................!...................!
          CALL SOLVE(DUMYI,RPM,PS,PA)                    ! ==>  calcsolnp.f
c     !........................................!...................!
          CALL FORCE(DUMYI)                    !
c     !........................................!...................!
c NOTE the pad moment is equal to Momk =
c      Momk = Moment from fluid film forces - elastic moment =
c      Momk = R x Fperp. - Krotpad x Rotdel:

          Momk=RSINPK*FX-RCOSPK*FY-KROTPAD*ROTDEL

c and the fluid film force normal to pad (along pivot) is:

          Fnorm=-RCOSPK*FX-RSINPK*FY           !=Fxsi = Pad normal force x R
          Fnorm=DABS(Fnorm)                    !

c      Fnorm & Fperp are the fluid film forces normal and perp. to pad.
c     !........................................!...................!
c
c    !==================================!////////////////////////////////////
C          IF ISOTH=1 (isothermal case) or IF fluid is LH2 then
C    ITPAD=0   Tlead=Tsump (no update of pad leading edge temperature)
C
c         FOR Thermal cases and IF IFULL=0  (pad bearing < 360 deg)
c    ITPAD=1   Update pad leading edge temperature.
C
C    ITPAD is defined on CPARAM of initopst.f
c    !==================================!////////////////////////////////////

c    !................................!.............................!
      IF ((ITPAD.EQ.1).AND.(QLEAD(KPAD).GT.0.0D0)) THEN
c    !................................!.............................!
c     IF Flow enters leading edge of KPAD, then calculate Leading edge
c     temperature from thermal mixing with upstream JPAD

        CALL GROOVETEMP(JPAD,KPAD,Kstart,TDIFF)   !==> calcthert.f
c
c      CHECK for temperature convergence on pad   Tdiff=|TLEADnew - TLEADold|
c      here we choose convergence to be 0.5 deg K ~ 1 deg F
c      since most thermometers are that accurate

        IF (TDIFF.LE.DTconv) THEN     !Tdiff less than 0.5 deg K
           IGCONV=IGCONV+1            !
           ITHpad(KPAD)=KPAD          !
        ELSE                          !.............................!
           ITHpad(KPAD)=0             ! Tdiff larger than 0.5 deg K
        END IF                        !
           KSKIP(KPAD)=0              !.............................!

c    !................................!.............................!
      ELSE                            ! all other cases
c    !................................!.............................!
c     IF QLEAD<=0 then flow leaves pad and leading edge temperature
c     is obtained from the solution within the film lands by
c     upwinding to the land internal points

        ITHpad(KPAD)=KPAD             !
        KSKIP(KPAD)=KPAD              ! Thus, Skip call of SOLVE
        IGCONV=IGCONV+1               ! and No update of groove temps.
c                                     ! i.e. convergence achieved
c    !................................!.............................!
      END IF                          !
c    !................................!.............................!


c     !........................................!...................!
c      IF PMAX=PMIN=PC means that pad is unloaded and fluid film
c                      forces, fluid induced pad moment, and all force
c         coefficients are equal to zero. In this case we simply
c         force the pad to rotate until it becomes loaded.
c
c      THIS argument is valid for a conventional tilt-pad. However,
c      for a flex-pad with KROTPAD >0, this means that there will be
c      an elastic moment.
c      In this case, the only valid solution should be that
c      as long as the pad is unloaded, then the pad rotation = 0 !!
c     !........................................!...................!
c
c UNLOADED PAD   >>>> ILOAD=0
c ::::::::::::
c        !............................!        !......................!
           IF ((PMAX.EQ.PC).AND.(DABS(Momk).LT.Merror)) THEN
c        !............................!        !......................!
             Momk=-KROTPAD*ROTdel              ! increase pad rotation
             ILOAD=0                           ! with a fixed angle
             ROTold=ROTdel                     ! until pad is loaded.
             ROTdel=ROTdel+Rdelta              ! ....................!

             KROTdel=KROTPAD                   ! since Kddk=0

             IF (KROTPAD.GT.ZERO) ROTdel=ZERO  !==> flexpad

             IF (ROTdel.GE.ROTtrail) THEN !........
                ROTdel=ROTdel-Rdelta      ! Do not allow to touch trail egde
                GOTO 776                  ! No convergence
             END IF                       !........
             Delta=ZERO                        !

c        !............................!        !......................!
           ELSE                                ! Loaded PAD
c        !............................!        !......................!
c LOADED PAD   >>> ILOAD=1                     ! at zero frequency
c ::::::::::				       ! find static stiffness
             ILOAD=1                           ! (first order soln.)
             FREQU=ZERO                        !
             CALL INPUT1(IFILE,DUMYI,RPM,PS,PA)          !==>optionsp.f==> firstp.f
             ROTold=ROTdel                     !

             KROTdel=KROTPAD+Kddk              !Pad rotational stiffness
c                                              !
c          Pad angle increment is: 	       !>>
           IF (KROTdel.GT.ZERO) THEN           !
             Delta=Momk/KROTdel                !
           ELSE                                !IF KROTdel<0 then we
             Delta=Rdelta                      !may converge to unloaded
           END IF                              !pad

c         ..........................>>>>       !
          IF (DABS(Momk).LT.Merror) GOTO 777   ! => convergence
c         ..........................>>>>       !
c					       !
c         SET new pad rotation                 !
             ROTdel=ROTold+Delta               ! Recall that here ex & ey
             ROToldd=ROTold                    ! are fixed.
c        !............................!        !......................!
           END IF                              !
c        !............................!        !......................!
c         Delta above denotes the increase in pad rotation angle
c         using N-R scheme to balance pad.
c        !.....................................!......................!
c
c In order for pad to be STABLE and converge to actual pad rotation
c the total pad moment stiffness MUST BE
c             KROTdel=KROTPAD+Kddk>=0.0
c otherwise angle Delta (increment) is in the wrong direction !.
c
c
c NOTE that for flex-pad bearings (KROTPAD>0) IF PMIN.EQ.PMAX then algorithm
c      for LOADED pad returns:
c      Delta==Momk/(KROTPAD+Kddk)=-KROTPAD*ROTDEL/KROTPAD=ROTDEL
c      Hence, new ROTDEL=ROTDEL-Delta=0.0
c
c        !.....................................!......................!
c          DISPLAY convergence in Screen or DUMP file
c        !.....................................!......................!
c
          IF (SVERB.GT.1) THEN
             WRITE (6,115) Count,ILOAD,ILOADo,ROTold,Delta,Momk,Fnorm
          END IF
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
             WRITE (1,115) Count,ILOAD,ILOADo,ROTold,Delta,Momk,Fnorm
          END IF
c        !.....................................!......................!


c         !=======================!            !----------------------!
            IF (COUNT.EQ.1) THEN               !
c         !=======================!            ! direct flow of
             Momkold=Momk                      ! calculations
             Deltaold=Delta                    !
             ILOADo=ILOAD                      ! At least 2 iterations
             GOTO 35                           ! needed
c         !=======================!            !----------------------!
            ELSE IF (COUNT.LT.COUNTMAX) THEN   !...! Count>1
c         !=======================!            !----------------------!

c           !...........................................!...=> PAD BALANCED
c###         IF (DABS(Momk-Momkold).LT.Momdiff) GOTO 777
c           !...........................................!...=> converged ...!
c###         IF (DABS(Delta).LT.Deldiff) GOTO 777
c           !...........................................!...=> converged ...!

             IF (DABS(ROTdel-ROTold).LT.Deldiff) GOTO 777
c           !...........................................!...=> converged ...!

c              Deldiff = 0.0005 (ROTtrail-ROTlead)  => 0.05 % difference

c HERE we propose to check the convergence to a solution by looking
c      at how the pad behaves between 2 iterations.

c           !................................................!
             IF ((ILOADo.EQ.1).AND.(ILOAD.EQ.1)) THEN
c            pad has been loaded on last 2 steps: continue search
c           !................................................!
                  ILOADo=ILOAD                 ! Reset values for
                  Momkold=Momk                 ! pad moment and
                  Deltaold=Delta               ! last angle
                  Icl=0                        !
                  GOTO 35                      !==> START over
c           !................................................!
             ELSE IF ((ILOADo.EQ.1).AND.(ILOAD.EQ.0)) THEN
c            pad was loaded and became unloaded: refine search
c           !................................................!
c            Pad rotation angle which caused pad to be loaded is ROToldd
c            current ROTdel=ROToldd+Delta caused pad to become unloaded
c      THEN, use oldd value and refine search with new pad angle
c            given as ROTdel=ROToldd+Delta/2.0D0

C##               ILOADo=ILOAD     ! Keep previous step as LOADED condition
C##               Momkold=Momkold  !

                  Delta=Deltaold/2.0D0       ! TO refine search USE a
c                                            ! SECANT like scheme.
                  Sign=+1.0D0                ! i.e. divide interval

                  IF (Icl.GE.Imax) GOTO 775   ! ==> give up refinement

                  ROTdel=ROToldd+Sign*Delta

                  Deltaold=Delta

                  Icl=Icl+1

                  GOTO 35                      !==> START over

c                 IF refinement does not provide a measureable difference
c                 then leave after a few its. (Imax=3)
c           !................................................!
             ELSE IF ((ILOADo.EQ.0).AND.(ILOAD.EQ.0)) THEN
c            pad is unloaded on 2 consecutive steps
c            continue increasing the pad rotation angle
c           !................................................!

                  ILOADo=ILOAD
                  Momkold=Momk
                  Deltaold=Delta
                  Icl=0
                  GOTO 35                      !==> START over
c           !................................................!
             ELSE IF ((ILOADo.EQ.0).AND.(ILOAD.EQ.1)) THEN
c            pad was unloaded and now is loaded, continue the search
c           !................................................!
c            use at least two iterations to make sure that
c            pad has been loaded sufficiently
c            otherwise may converge to unloaded eccentricity

             IL10=1    !### 2/15 /94 skip over rotation

               IF (IL10.EQ.0) THEN
                 IL10=1
                 ILOADo=0
                 ROTold=ROTdel          !==> increase ROTdel
                 ROTdel=ROTdel+Rdelta   !    by a fixed amount
                 IF (ROTdel.GE.ROTtrail) THEN
                        ROTdel=ROTdel-Rdelta
                        GOTO 776
                 END IF
                 Delta=ZERO
               ELSE
                 IL10=0
                 ILOADo=ILOAD
                 Momkold=Momk
                 Deltaold=Delta
               END IF

               Icl=0
               GOTO 35                      !==> START over



             END IF

c         !=======================!            !----------------------!
            ELSE IF (COUNT.EQ.COUNTMAX) THEN   ! not converged
c         !=======================!            !----------------------!
                  GOTO 776
c         !=======================!            !----------------------!
            END IF
c         !=======================!            !----------------------!



c
c .........................................................................!
 775      WRITE (6, 665) KPAD, COUNT            ! Failed refinement
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
                 WRITE (1, 665) KPAD, COUNT
          END IF
          ITbal(KPAD)=0
          GOTO 779
c .........................................................................!
 776      ROTDEL=ROTold                         ! NO CONVERGENCE
          WRITE (6, 666) KPAD, COUNT            ! Failed
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
                 WRITE (1, 666) KPAD, COUNT
          END IF
          ITbal(KPAD)=0
          GOTO 779

c .........................................................................!
c  CONVERGENCE:   Moment difference is null or almost zero for current pad
c .........................................................................!

 777      CONTINUE
          ITbal(KPAD)=1
          IF (SVERB.GT.1) THEN
             WRITE (6,115) Count,ILOAD,ILOADo,ROTdel,Delta,Momk,Fnorm
          END IF
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
             WRITE (1,115) Count,ILOAD,ILOADo,ROTdel,Delta,Momk,Fnorm
          END IF
c .........................................................................!

 779   ROTPAD(KPAD)=ROTDEL                     !==> save pad rotation

c     !........................................!............................!
c      CHECK HERE FOR CONVERGENCE ON PAD LEADING EDGE TEMPERATURES.
c     !........................................!............................!

       IF (ITHpad(KPAD).EQ.0   ) THEN          ! Groove T = Leading edge T
         Count=0                               ! have not converged
         IF (ITERTH.GE.ITERTHmax) GOTO 844
         ITERTH=ITERTH+1                       ! start OVER AGAIN
         GOTO 35                               !
       END IF                                  !


         JPAD=KPAD                             !==> get ready for next pad
         ITERTH=1                              !
         GOTO 855
c     !........................................!............................!
 844   WRITE (6,845) KPAD, ITERTH
 845   FORMAT (3X,'PAD #:',I3,' Thermal Analysis did not CONVERGE after'
     +        ,1X,I3,' STEPS',/,3X,70('.'))
       GOTO 866
c     !........................................!............................!
 855   IF (ITPAD.EQ.1) THEN
         WRITE (6,856) KPAD, ITERTH
 856     FORMAT (3X,'PAD #:',I3,' Thermal Analysis DID CONVERGE after'
     +        ,1X,I3,' STEPS',/,3X,70('.'))
       END IF
c     !........................................!.......................!

 866     JPAD=KPAD                             ! Current pad -> previous pad
         ITERTH=0                              !
         CALL TORQUE(DEVICE)                   ! Calculate drag torque
         CALL ADD0                             ! add forces and torque
c     !........................................!.......................!
       IF (LOSTAT.EQ.1) THEN                   ! for LOAD calculations
c     !......................!                 !.......................!
         IF ((PMAX.EQ.PC).AND.(KROTPAD.EQ.0.0D0)) GOTO 888  !=> unloaded pad
         Delta=Momk/(KROTPAD+Kddk)             ! Momk=Mk - Kpad x Del
         FXTmom=FXTmom-KXdk*Delta              ! Residual forces induced
         FYTmom=FYTmom-KYdk*Delta              ! by unbalanced pads.
         GOTO 888                              !

c     !......................!                 !.......................!
       END IF                                  ! LOSTAT=1
c     !......................!                 !.......................!

c     !......................!                 !.......................!
       IF (ITPAD.EQ.0) THEN                    !ISOTHERMAL ANALYSIS
c     !......................!                 !.......................!
       WRITE(DEVICE,111) KPAD                  ! Print results
       CALL PRINTDUVP(DEVICE, 0)               ! and Calculate force
       IF (FREQUB.EQ.ZERO) GOTO 877            ! coefficients
       FREQU=FREQUB                            ! at frequency specified
       CALL INPUT1(IFILE,DEVICE,RPM,PS,PA)               ! if different from
 877   CALL PRINTCOEF(DEVICE)                  ! zero
c     !......................!                 !.......................!
       END IF                                  ! ITPAD=0
c     !......................!                 !.......................!

 888   CALL ADD1                               ! Add force coefficients
c                                              !
       CALL WRITETEMP(KPAD)                    ! -> temporal storage

C   !--------------------!------------------------------------------!
       END DO            ! KPAD=1, NPAD
C   !--------------------!------------------------------------------!

c   ----------------------------------------------------------------
c   ISOTHERMAL ANALYSIS RETURN TO MAIN CALL
       IF (ITPAD.EQ.0) GOTO 955
c   ----------------------------------------------------------------



c   ============================================================
c   THERMOHYDRODYNAMIC ANALYSIS proceeds with Iterative Solution
c   ============================================================
c

       IGconv=0          ! Check for convergence on all pad leading edge
c     !..................! temperatures.
       DO KPAD=1, NPAD
         IF (ITHpad(KPAD).EQ.KPAD) IGconv=IGconv+1
       END DO
c     !..................!

       IF (IGconv.EQ.NPAD) GOTO 933

       WRITE (6,901)
 901   FORMAT (/,3X,'NEW THERMAL LOOP:', 50('*'),/)
       GOTO 77



c     !........................................!.......................!
C       ONCE THERMAL SOLUTION IS CONVERGED and PADS ARE BALANCED
c       THEN CALCULATE FORCE COEFFICIENTS AT REQUESTED FREQUENCY
c     !........................................!.......................!

 933  CONTINUE

c     !........................................!.......................!
       IF (LOSTAT.EQ.1) RETURN                 ! for LOAD OPTION
c     !......................!                 !.......................!

      FREQU=FREQUB

      CALL ZERO0         !
      CALL ZERO1         ! zeroes global force coefficients
      JCASE=2            ! Perturb in X & Y directions
      IFILE=1            ! ............................


C   !--------------------!
      DO KPAD=1, NPAD    ! SWEEP OVER ALL PADS:
C   !--------------------!
       WRITE(DEVICE,111) KPAD                  !
       CALL XYDATA(KPAD,RPM,PS,PA)                       ! generates Mesh
       CALL READTEMP(KPAD)                     ! reads input data
       CALL ADD0                               ! Adds forces & moments
       CALL PRINTDUVP(DEVICE, 0)               ! and prints
       CALL PRINTF(DEVICE,1)                   !
       IF (FREQU.NE.0.0D0) THEN                !
         CALL INPUT1(IFILE,DEVICE,RPM,PS,PA)             ! find force coeffs.
       END IF
       CALL PRINTCOEF(DEVICE)                  ! & prints them OUT
       CALL ADD1                               ! Adds stiffness coefficients
       CALL WRITETEMP(KPAD)                    ! -> temporal storage
C   !--------------------!
       END DO            ! KPAD=1, NPAD
C   !--------------------!


C PRINTS TOTAL FORCES AND COEFFICIENTS FOR NPADS:
C -----
 955    IF (LOSTAT.EQ.1) RETURN                 ! for LOAD OPTION
        CALL PRINT0PAD(DEVICE)
        CALL PRINT1PAD(DEVICE,RPM,PS,PA)

C------------------------!.............!...................................!
       RETURN                          ! To MAIN on hydroflex.f
C------------------------!.............!...................................!


C!...........................................!..........................
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) '$ I/O ERROR on SUB: TILTECC, RETURN to MAIN.'
      CALL BEEPER
      CLOSE (UNIT=42)
      RETURN

C!...........................................!..........................
C
 110  FORMAT (///,'+',77('-'),'+',
     +       /,4X,'START ITERATIVE PROCEDURE TO BALANCE PAD',
     +            ' MOMENT for PAD #:', I3,/,
     +            '+',77('-'),'+')

 111  FORMAT (' +',77('-'),'+',
     +       /,' |',3X,'RESULTS FOR BEARING PAD #',
     1        I2,46X,' |',/,' |',77('.'),'|')
 112  FORMAT (//,'PAD:',I3,'  ITERATION #', I3, //)

 115  FORMAT (4X,'#:',I2,' IL:',I2,1X, I2,' rot:',E11.4E2,
     +        ' DD:',E11.4E2, ' Mk:',E11.4E2,' Fxi:',E10.3E2)

 116  FORMAT (/,3X,'PAD:',I2,' It:',I2,/, 3X,
     +                     'Moment(I):',E11.4E2, 2X,
     +                     'Moment(I-1):',E11.4E2,/,3X,
     +                        'Rotpad(I):',E11.4E2, 2X,
     +                     'Rotpad(I-1):',E11.4E2, 2X,
     +                     'Max diff:', E11.4E2)
 665   FORMAT(' ',3X,30('* '),/,4X,
     + 'Failed to improve search on PAD',I3,' after ', I3, 'its.',
     + /,4X,30('* '))
 666   FORMAT(' ',3X,30('* '),/,4X,
     + 'NO convergence on PAD',I3,' after ', I3, 'its.',
     + /,4X,30('* '))

 6099 FORMAT ('$', 'ENTER Excitation Frequency(w) in Hz [Default (w)=',
     + E12.5E2,']: ')

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Loadtilt                                                    **
C **                                                                         **
C **  LOADTILT:  Given external load calculate equilibrium eccentricity.     **
C **         for tilt-pad bearing                                            **
C *****************************************************************************
C *****************************************************************************

      SUBROUTINE LOADTILT(DEVICE,RPM,PS,PA,WX,WY)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PARAM2/ EXO, EYO
      COMMON /WEAR/ EWX, EWY, EWEAR, BETAw, IWEAR
      COMMON /FORCE0/ FFACTOR,FX,FY,TO,TOR
      COMMON /RESULTS0/ FXT,FYT,TOT,MXT,MYT,QINT,QOUTT
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /FREQ/ FREQU,SIGMA,L1,RES,JCASE,NCASE
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT
      COMMON /PMINMAX/ PMIN, PMAX

      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /STIAT/ KXXAT,KYYAT,KXYAT,KYXAT,KmXXAT,KmYYAT,KmXYAT,KmYXAT
      COMMON /DAMPT/ CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT

      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /ROTAPAD/ ROTPAD(MAXNPAD), TILTPAD


      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/  KROTPAD, CROTPAD
      COMMON /RESFXY / FXTmom, FYTmom, DEXC, DEYC


      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
c     .........................................................
      DOUBLE PRECISION ROTPAD, PRELOAD, OFFSET, ROTDEL
      DOUBLE PRECISION CLEAR,DIAM,LENGTH,LD,AR,HREC,
     +                 EXO, EYO, EWX,EWY, EWEAR,BETAw,
     +                 FFACTOR, FX, FY, TO, TOR,
     +                 FXT,FYT,TOT,MXT,MYT,QINT,QOUTT,
     +                 EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 PRATIO,CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 FREQU,SIGMA,L1,RES,AC, ETA, PBACK,
     +                 KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,
     +                 CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,
     +                 CmXYDT,CmYXDT,
     +                 KmXYDT,KmYXDT,KXXAT,KYYAT,KXYAT,KYXAT,
     +                 KmXXAT,KmYYAT,KmXYAT,KmYXAT

      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD, KROTpad, CROTpad,
     +                 PMIN,PMAX

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        NPAD, NREC, IFULL, TILTPAD, TILT, LIFT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        ITER, ITMAX, ITPMAX, SVERB,DVERB, BEEP,
     +        JCASE,NCASE, IWEAR, ISYM, ICSTEP


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      REAL TIME
      DOUBLE PRECISION EMAX, DEMAX, EPSLOAD, WX, WY, W, FREQUold,
     +                 ECO, EXOLD, EYOLD, ERRX, ERRY, DET,
     +                 DEX, DEY, DECC, ECC, ZERO, DAMP,
     +                 FXTmom, FYTmom, Momk, DEXC, DEYC,
     +                  RPM,PS,PA

      INTEGER I, J, ITMAXU, ITPMAU, IFR, ITLOAD, ITLDMAX,
     +        CHO, DEVICE, DUMYI, IFILE, KPAD, IOS, LOSTAT
      CHARACTER FILE*80, YN*1

C ----------------------------------------------------------------------------
C --  LOADTILT code for calculation of eccentricity on TILT-PAD bearings
C ----------------------------------------------------------------------------
      OPEN (UNIT=42, STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)

      NCASE=1                          ! FOR one recess depth
                                       !.........
      ZERO=0.0D0                       !
      EMAX=0.99999999                  ! Maximum allowed eccentricity
      EMAX=1.20000000
      EMAX=EMAX*(1.0D0-PRELOAD/CLEAR)  ! based on preload value

      DEMAX=0.0001                     ! Maximum change on Ex & Ey for conv
      EPSLOAD=.010                     ! Convergence load criteria of 0.5%
      ITLDMAX=10                       ! Maximum number of iterations for load
      FREQUold=FREQU                   !....................................
c                                      !
      WX=0.0D0                         ! Dumy load values
      WY=0.0D0                         !
      DEX=0.0D0                        !
      DEY=0.0D0                        !
      ITLOAD=0                         !
      DUMYI=0                          ! Dumy Device
C--------------------------------------!--------------------------!
C NOTE: In tilt pad bearings, stiffness coefficients at zero
c       frequency may be quite different than those at synchronous
c       frequency especially for low eccentricities (large Sommerfeld
c       numbers). Thus, user needs to make sure a good initial guess
c       is provided for starting LOAD routine with choice (2)
c       3/16/93 LSA
C--------------------------------------!--------------------------!

C......................................!.......................
c                                      ! calculations are performed
      FREQU=ZERO                       ! at zero whirl frequency
C......................................!.......................      LOSTAT=1                         ! status for TILTECC routine
C
C......................................!.......................
  100 WRITE (6, 110)
  110 FORMAT ('$', 'ENTER Load(Wx) in [N] along X dir.: ')
      CALL ENTERVAL(WX)

      WRITE (6, 120)
  120 FORMAT ('$', 'ENTER Load(Wy) in [N] along Y dir.: ')
      CALL ENTERVAL(WY)

      WRITE (6, 130)
  130 FORMAT ('$', 'INPUT Max # Iters for load calculation (<=10): ')
      CALL ENTERINT(ITLDMAX)

      IF (ITLDMAX.LT.1) ITLDMAX=10
      IF (ITLDMAX.GT.10) ITLDMAX=10


      WRITE (6, 67) WX,WY
      IF (DEVICE.EQ.1) THEN
          WRITE (1, 67) WX,WY
      END IF

C                                      !.......................
      W=DSQRT(WX*WX+WY*WY)             ! Total load
      IF (W.LT.EPSLOAD) THEN           !.......................
          WRITE (6, *) 'LOAD IS TOO SMALL, TRY AGAIN'
          GOTO 100                     !
      END IF                           !

c##   WRITE (6, 90) WX, WY
      WRITE (6, *) 'SELECT  one of:'
      WRITE (6, *) ' (1) Use a guess at Ex=Ey=0 to start soln'
      WRITE (6, *) ' (2) Use current solution read from DATA file'
      WRITE (6, *) '     as start guess'
      write (6, *) '     DATA for guess must be already on memory.'
      write (6, *) '   '
      write (6, *) '     ENTER Any number >2, to EXIT'

      WRITE (6, 51)
   51 FORMAT ('$', 'SELECT: ')
      CALL ENTERINT(CHO)
      CLOSE(UNIT=42)

      JCASE=2    ! for perturbations in X & Y directions

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
      CALL CPARAM(RPM,PS,PA)          ! Calculate parameters for calc.
      CALL XYDATA(0,RPM,PS,PA)                   ! set and save mesh on PLOTMESH
      CALL ZERO0                       ! zeroes total forces and coeffics.
      CALL ZERO1                       ! ....................

      ITMAXU=ITMAX                     !
      ITPMAU=ITPMAX                    !

      FXTmom=ZERO
      FYTmom=ZERO

c    !..................................!........................!
      IF (TILTPAD.EQ.1) THEN            ! for tilt pad bearings
        ROTDEL=0.25D0*CLEAR/(DIAM/2.0D0)! set pad rotaion to 25%
        DO KPAD=1, NPAD                 ! of pad clarance
         ROTPAD(KPAD)=ROTDEL            ! initial guess only
        END DO                          !
      END IF                            !
c    !..................................!........................!

c    !................!                !..................!
      DO KPAD=1, NPAD                  !
c    !................!                ! SWEEP on all pads
        CALL XYDATA(KPAD,RPM,PS,PA)              !..................!
        CALL FILMH(EXO, EYO, NXT, NYI) ! no wear on bearing
        CALL SETBC(RPM,PS,PA)                     ! Sets boundary condition
      IF (KPAD.EQ.1) CALL GUESP(DUMYI) ! Calculates initial guess for pratio
        CALL SOLVE(DUMYI,RPM,PS,PA)              ! SOLVE & calculate:
C       CALL PRINTDUVP(DUMYI, 0)       ! forces & torque
C##     WRITE(6,111) KPAD              !
        CALL FORCE(DUMYI)              !Dumyi=0 device
        CALL TORQUE(DUMYI)             !
        IFILE=3                        !
        ROTDEL=ROTPAD(KPAD)            !
        CALL INPUT1(IFILE,DUMYI,RPM,PS,PA)       ! finds stiffness coeffs.
        CALL ADD0                      !
        CALL ADD1                      ! & sums forces & coeffs.
c     !........................................!...................!
       IF (TILTPAD.EQ.1) THEN                  ! calculate loads due
          IF ((PMAX.EQ.PC).AND.(KROTPAD.EQ.0.0D0)) GOTO 450
c                                              ! to non-zero Moments
          Momk=RSINPK*FX-RCOSPK*FY -           ! in tilt-pad bearings
     +         KROTPAD*ROTDEL                  !
          FXTmom=FXTmom-KXdk*Momk/(Kddk+KROTPAD)
          FYTmom=FYTmom-KYdk*Momk/(Kddk+KROTPAD)
       END IF                                  ! PMIN=PC=PMAX pad unloaded
c     !........................................!...................!
 450    CALL WRITETEMP(KPAD)                   ! & store in temporal file
c    !.........................................!...................!
      END DO                           !
c    !................!                ! K=1,2,..., NPAD

c     TIME=SECNDS(0.0)                 ! t=0.0
      GOTO 600                         ! GOTO Start
C                                      !------------------------------------
C......................................!
C Read file with guessed solution:     ! CHO=2
C......................................!......
  500 ITMAX=MAX0(ITMAX,99)             !
      ITPMAX=MAX0(ITPMAX,10)           !
      ITMAXU=ITMAX                     !
      ITPMAU=ITPMAX                    !

      CALL CPARAM(RPM,PS,PA)                      ! Calculate parameters for calc.
      CALL XYDATA(0,RPM,PS,PA)                   ! and generates mesh


      IF (FREQUold.GT.ZERO) THEN       ! => calculate Kij's at w=0
        WRITE(6,*) '$ -----------------------------------------------'
        WRITE(6,*) '$ Calculate start force coeffs. at zero frequency'
        WRITE(6,*) '$ -----------------------------------------------'
        GOTO 525
      END IF

      IF ((CXXDT.NE.ZERO).AND.(CYYDT.NE.ZERO)) GOTO 525

      IF ((KXXDT.NE.ZERO).AND.(KYYDT.NE.ZERO)) GOTO 600

c...................................................................c
c for tilt-pads if damping coefs. are different from zero
c we need to simply add pad coeffs. at zero freq. from data
c stored from file. (Of course this would be valid iff land
c inertia is not accounted, i.e. INERL=0)
c 3/16/93
c...................................................................c

 525  IFILE=1

      CALL ZERO0
      CALL ZERO1
      FXTmom=ZERO
      FYTmom=ZERO

c     !........................................!...................!
      DO KPAD=1, NPAD
c     !........................................!...................!
          CALL XYDATA(KPAD,RPM,PS,PA)
          CALL READTEMP(KPAD)
          CALL TEMPROPS
          CALL MINMAXP
c#        WRITE (6,111) KPAD
          CALL INPUT1(IFILE, DEVICE,RPM,PS,PA)
          CALL ADD0
          CALL ADD1
c     !........................................!...................!
       IF (TILTPAD.EQ.1) THEN                  ! calculate loads due
          IF ((PMAX.EQ.PC).AND.(KROTPAD.EQ.0.0D0)) GOTO 530
c                                              ! to non-zero Moments
          Momk=RSINPK*FX-RCOSPK*FY  -          ! for tilt-pads
     +         KROTPAD*ROTDEL                  !
          FXTmom=FXTmom-KXdk*Momk/(Kddk+KROTPAD)
          FYTmom=FYTmom-KYdk*Momk/(Kddk+KROTPAD)
       END IF                                  !PMIN=PC=PMAX pad unloaded
c     !........................................!...................!
 530     CALL WRITETEMP(KPAD)
c     !........................................!...................!
      END DO
c    !................!                ! K=1,2,..., NPAD

 550  GOTO 600                         ! GOTO Start


C
C                                      !------------------------------
C--------------------------------------! START of LOAD Calculations  !
C                                      !------------------------------
C                                      !
  600 ITLOAD=ITLOAD+1                  !................................

      ECO=DSQRT(EXO*EXO+EYO*EYO)
      EXOLD=EXO
      EYOLD=EYO

      IFR=0 !............................................................

 615  ERRX=(WX+FXT+FXTmom)              ! Error in load calculation
      ERRY=(WY+FYT+FYTmom)              !


      WRITE (6,82) EXO,WX,FXT,FXTmom,ERRX,
     +             EYO,WY,FYT,FYTmom,ERRY
      WRITE (6,85) KXXDT,KYYDT,KYXDT,KXYDT

      IF (DEVICE.EQ.1) THEN
          WRITE (1,82) EXO,WX,FXT,FXTmom,ERRX,
     +                 EYO,WY,FYT,FYTmom,ERRY
          WRITE (1, 85) KXXDT,KYYDT,KYXDT,KXYDT
      END IF


      IF (ITLOAD.GT.ITLDMAX) THEN       ! Maxim. # of iterations exceeded
          GOTO 900                      ! -> Nonconverged
      END IF                            !................................

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
          GOTO 1200                       ! STATICALLY UNSTABLE bearing
      END IF                              !.........................

      DEX=(KYYDT*ERRX-KXYDT*ERRY)/DET/CLEAR ! Variation in eccentricity
      DEY=(KXXDT*ERRY-KYXDT*ERRX)/DET/CLEAR !
                                       !..........................

  700 EXO=EXOLD+DEX                    ! Improved eccentricy
      EYO=EYOLD+DEY                    ! components Ex & Ey
                                       !..........................

      DECC=DSQRT(DEX*DEX+DEY*DEY)      ! Total change in ECC
      ECC=DSQRT(EXO*EXO+EYO*EYO)       ! & new value of Ecc

          WRITE (6, 10)
          WRITE (6, 70) ITLOAD, EXO, EYO, ECC
      IF (DEVICE.EQ.1) THEN
          WRITE (1, 10)
          WRITE (1, 70) ITLOAD, EXO, EYO, ECC
      END IF

C
C
      IF (CHO.EQ.2) GOTO 721           !

      IF (DABS(DECC).LE.DEMAX) THEN    !
          IFR=2                        ! Check for maximum difference
          GOTO 615                     ! in eccentricity between 2 its.
      END IF                           !...............................

 721  CHO=1

      IF ((ECC.GT.EMAX).AND.(AC.EQ.ZERO)) THEN
         DAMP=(EMAX+(ECO-EMAX)*EXP(ECO-ECC)-ECO)/(ECC-ECO)
         DEX=DEX*DAMP                  !
         DEY=DEY*DAMP                  ! DAMP large calculated Ex &Ey
         GOTO 700                      !
      END IF                           !...............................

      ITMAX=ITMAXU
      ITPMAX=ITPMAU

      FXTmom=ZERO
      FYTmom=ZERO

c    !........................................!...................!
       DEXC=DEX*CLEAR                         ! Delta changes in
       DEYC=DEY*CLEAR                         ! journal eccentricity
c     !.......................................!...................!
      CALL TILTECC(DUMYI,LOSTAT,RPM,PS,PA)              ! Balance Pads for given Ecc
c     !.......................................!...................!


C------------------------!.............!...................................!
         GOTO 600                      ! G0TO Start
C------------------------!.............!...................................!





c 800 TIME=SECNDS(TIME)                      !-------------------------
  800 CONTINUE                               !
      IF (DEVICE.eq.1) THEN                  ! Convergence achieved
          WRITE (1, 47) ITLOAD               ! ....................
          WRITE (1, 60) EXO, EYO, WX, WY     ! Wx + Fx = 0
      END IF                                 ! Wy + Fy = 0
      WRITE (6, 47) ITLOAD                   !
      WRITE (6, 60) EXO, EYO, WX, WY         !.........................

      DUMYI=1                                ! Print Table of results
      CALL PRINT0PAD(DEVICE)                 !
      CALL PRINT1PAD(DEVICE,RPM,PS,PA)                 !
c                                            !.........................
      WRITE (6, 83)


c------------------------!...................!.........................!
c CALCULATE FORCE COEFFICIENTS at given frequency
c------------------------!...................!.........................!

c    !...........................................................!
c     INPUT FREQUENCY FOR CALCULATION of COEFFICIENTS
c    !................................................!..........!
      OPEN (UNIT=42,  STATUS='SCRATCH', ERR=1234,
     +      IOSTAT=IOS)
          YN='Y'
          WRITE (6, 1520)
 1520     FORMAT ('$', 'Do you wish to calculate dynamic ',
     +            'coefficients [Def: Y] (Y/N): ')
     	  READ (5, 1521) YN
 1521     FORMAT (1A)

          IF ((YN.eq.'y').OR.(YN.eq.'Y')) goto 1522
          CLOSE(42)
          RETURN

 1522     FREQU=RPM/60.0D0
          IF (RPM.EQ.0.0D0 )FREQU=100.0D0
          WRITE (6, 1523) FREQU
 1523 FORMAT ('$', 'ENTER Excitation Frequency(w) in Hz [Default (w)=',
     + E12.5E2,']: ')
      CALL ENTERVAL(FREQU)
      CLOSE(42)
c    !.................................!
      IF (FREQU.EQ.0.0D0) RETURN       ! load calcs. were made at 0freq.
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
       JCASE=2                                 !.....
       CALL INPUT1(IFILE,DEVICE,RPM,PS,PA)               ! find stiffness coeffs.
       CALL ADD1                               ! Adds stiffness coefficients
       CALL WRITETEMP(KPAD)                    ! -> temporal storage
C   !--------------------!
       END DO            ! KPAD=1, NPAD
C   !--------------------!

C PRINTS TOTAL FORCES AND COEFFICIENTS FOR NPADS:
C -----
        CALL PRINT0PAD(DEVICE)
        CALL PRINT1PAD(DEVICE,RPM,PS,PA)

C------------------------!.............!...................................!
       RETURN                          ! To MAIN on hydroflex.f
C------------------------!.............!...................................!


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
      WRITE (6, 10)                          ! Eccentricity > 1.0
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

C!...........................................!..........................
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) '$ I/O ERROR on SUB: LOAD, RETURN to MAIN.'
      CALL BEEPER
      CLOSE (UNIT=42)
      RETURN

C!...........................................!..........................
C
   10 FORMAT (' ', 79('.'))
   67 FORMAT (' ', 79('.'),/,3X,'SPECIFIED LOAD: WX=',E12.5E2,
     + ' N',2X, 'WY=', E12.5E2,' N',/,1X, 79('.'))
   20 FORMAT (' ', 'SUB LOAD: NEW ECC>1, PROGRAM ABORTS EXECUTION')
   30 FORMAT (' ', 3X, 70('.'), /, ' ',
     +    'DETERMINANT of Kij=0, STATICALLY UNSTABLE BEARING ',
     +    'OR PAD', /, ' ', 70('.'))
   40 FORMAT (' ', 4X, 'LOAD Convergence in Itload:', I3, 3X,
     +        'Time:', F9.2, ' sec.')
   50 FORMAT (' ', 3X, 'NO CONVERGENCE in Itload:', I3, 3X,
     +        'Time:', F9.2, ' sec.')
   47 FORMAT (' ', 79('='),/,3X,
     +'LOAD Convergence in ', I3, ' iterations')
   57 FORMAT (' ', 3X, 'NO CONVERGENCE in Itload:', I3)
   60 FORMAT (' ', 39('= '), /, 4X, 'Eccentricity EXo:', F7.5,
     +        3X, 'EYo:', F7.5, /, 4X, 'for Load WX:',
     +        E12.5E2, 3X, 'WY:', E12.5E2, /, ' ', 79('='))
   70 FORMAT (' ', 39('. '),/,
     +        ' LOAD: Iter=', I3, ' NEW Exo:', F7.5,
     +        2X,'Eyo:',F7.5,2X,'Ecc:',F7.5,/,' ',39('. '))
   82 FORMAT (' ',79('.'),/,
     +        ' ',3X,'EX=', F6.4,3X,'WX=',E10.3E2,'=',
     +               'FXT=',E10.3E2,1X,'+ FXTm=',E10.3E2,2X,
     +            'errorX=',E10.3E2,/,
     +        ' ',3X,'EY=', F6.4,3X,'WY=',E10.3E2,'=',
     +               'FYT=',E10.3E2,1X,'+ FYTm=',E10.3E2,2X,
     +            'errorY=',E10.3E2)

   83 FORMAT (' ',79('/')
     +    ,/,4X,'NOTE: Save Results in Data File and use OPT (35)'
     +    ,/,10X,'to find force coefficients at other frequencies'
     +    ,/,' ', 79('/'))

   85 FORMAT (' ',3X,'KXX=',E12.5E2,2X,'KYY=',E12.5E2,2X,
     +             'KYX=',E12.5E2,2X,'KYY=',E12.5E2,2X,'[N/m]',/,
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

C *****************************************************************************
C **                                                                         **
C **  Subroutine FINDROT     for TILT PAD BEARINGS                           **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FINDROT(DEX,DEY,ROTDEL)

      IMPLICIT NONE

      INCLUDE 'params.f'
C ------------------------------------------------------------------------
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD,AR,  HREC
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /FORCE0/ FFACTOR,FX,FY,TO,TOR
      COMMON /PMINMAX/ PMIN, PMAX

      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD

      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/  KROTPAD, CROTPAD

      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD,AR,  HREC,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD

      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      DOUBLE PRECISION FFACTOR,FX,FY,TO,TOR,PMIN,PMAX,
     +                 RSINPK, RCOSPK, IPAD, KROTPAD, CROTPAD
      INTEGER TILT

C -----------------------------------------------------------------
C    LOCAL VARIABLES
C -----------------------------------------------------------------
      DOUBLE PRECISION MomentK, DEX,DEY, ROTDEL, DELTA, DELnew,
     +                 DD, KDD, ZERO, Mbalan
      INTEGER UNLOAD
C -----------------------------------------------------------------
C Find new value of rotation of pad about pivot at frequency w=0
C
C ------------------------------------------------------------------------
      ZERO=0.0D0
c    !..............................!
C     Fluid induced Moment at Pad K is equal to:
C     Mk = R { sin(Tetapivotk) x FXk - cos(Tetapivotk) x FYk }
c    !..............................!
      Momentk=RSINPK*FX - RCOSPK*FY !
c                                   !.............................!
      UNLOAD=0                      !
      CALL MINMAXP                  ! find out if pad is unloaded
      IF (PMAX.eq.PC) UNLOAD=1      !.............................!
c                                   ! IF PAD is unloaded then
c                                   ! all pad force coefficients = 0
c                                   ! in particular, Kddk = 0
c                                   !.............................!
      KDD = Kddk+ KROTPAD           ! total rotational stiffness
      Mbalan=MomentK-KROTPAD*ROTDEL ! Mk- Kpad x Del = 0 for equilibrium

c ....change in pad rotation angle (rads) is given as


c    !..............................!
      IF ((UNLOAD.EQ.1).AND.(KROTPAD.EQ.ZERO)) THEN
c    !..............................! for unloaded pad.
         Delta=0.0D0
c    !..............................!
      ELSE
c    !..............................! for loaded pad
         Delta=(Mbalan-KdXk*DEX-KdYk*DEY)/KDD
c    !..............................!
      END IF
c    !..............................!

c .... improved pad rotation angle (rads) is equal to

      ROTDEL=ROTDEL+Delta           ! in radians

c     Check if new pad rotation is within bounds is performed on
c     TILTECC subprogram.

c    !..................................................!

      END


C *****************************************************************************
C *****************************************************************************
C **                                                                         **
C **  Subroutine REDUCEs COEFFICIENTS for TILT PAD BEARINGS                  **
C **                                                                         **
C *****************************************************************************
      SUBROUTINE REDUCECOEFS

      IMPLICIT NONE

      INCLUDE 'params.f'
C ------------------------------------------------------------------------
      COMMON /FORCE0/ FFACTOR,FX,FY,TO,TOR
      COMMON /PMINMAX/ PMIN, PMAX

      COMMON /STIFT/ KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT
      COMMON /DAMPT/ CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT

      COMMON /STIFF/ KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD
      COMMON /DAMPI/ CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD

      COMMON /TILTCOE/ KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk

      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /ROTPAD/ KROTPAD, CROTPAD
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASE, NCASE

      DOUBLE PRECISION FFACTOR,FX,FY,TO,TOR,PMIN,PMAX,
     +         KXXDT,KYYDT,KXYDT,KYXDT,KmXXDT,KmYYDT,KmXYDT,KmYXDT,
     +         CXXDT,CYYDT,CXYDT,CYXDT,CmXXDT,CmYYDT,CmXYDT,CmYXDT

      DOUBLE PRECISION KXXD,KYYD,KXYD,KYXD,KmXXD,KmYYD,KmXYD,KmYXD,
     +                 CXXD,CYYD,CXYD,CYXD,CmXXD,CmYYD,CmXYD,CmYXD

      DOUBLE PRECISION KdXk,KdYk,KXdk,KYdk,Kddk,
     +                 CdXk,CdYk,CXdk,CYdk,Cddk,
     +                 RSINPK, RCOSPK, IPAD, KROTPAD, CROTPAD,
     +                 FREQU, SIGMA, L1, RES
      INTEGER ICASE, NCASE, TILT

C -----------------------------------------------------------------
C    LOCAL VARIABLES
C -----------------------------------------------------------------
      DOUBLE PRECISION FREQ, ZERO, RZDD,IZDD
      DOUBLE COMPLEX ZXX,ZXY,ZYX,ZYY,ZXD,ZDX,ZYD,ZDY,ZDD,
     +               ZXXR,ZXYR,ZYXR,ZYYR

C -----------------------------------------------------------------
c    Calculate pad reduced force coefficients for harmonic
c    motion exp(iwt) with frequency w
c    Ipad : mass moment of inertia of k-pad in [kg.m2]
c    Krotpad: pad rotational stifness of k-pad in [N.m/rad]
c    Crotpad: pad rotational damping of k-pad in [N.m.s/rad]
C ------------------------------------------------------------------------

      FREQ=2.0D0*FREQU*DACOS(-1.0D0)   ! excit. frequency in rad/s
      ZERO=0.0D0

c    !...............................!
c     IF PAD forces=0 then NO impedance coefficients
c        reduced coefficients = 0.0D0
c    !...............................!

c##     IF ((FX.EQ.0.0D0).and.(FY.EQ.0.0D0)) RETURN

        IF (DABS(PMIN-PMAX).LT.1.0D0) THEN       ! Tilt pad is unloaded
          RETURN                     ! ................
        END IF                       ! does not provide any coefficients
c                                    ! i.e. all Kij = 0

c    !...............................!
C     set fluidic impedances: Z = K + i w C
c    !...............................!
      ZXX=DCMPLX(KXXD,FREQ*CXXD)
      ZYY=DCMPLX(KYYD,FREQ*CYYD)
      ZXY=DCMPLX(KXYD,FREQ*CXYD)
      ZYX=DCMPLX(KYXD,FREQ*CYXD)

      ZXD=DCMPLX(KXDK,FREQ*CXDK)
      ZYD=DCMPLX(KYDK,FREQ*CYDK)
      ZDX=DCMPLX(KDXK,FREQ*CDXK)
      ZDY=DCMPLX(KDYK,FREQ*CDYK)

      RZDD=KDDK+KROTPAD-FREQ*FREQ*IPAD
      IZDD=FREQ*(CDDK+CROTPAD)
      ZDD=DCMPLX(RZDD,IZDD)


c    !................................!
      IF ((RZDD.EQ.ZERO).AND.(IZDD.EQ.ZERO)) THEN
       WRITE (6,123)
       STOP
      END IF

 123  FORMAT ('$', 2X,50('='),/,3X,
     +        'ERROR: Pad Moment Impedance Zddk=0+i0',/,3X,
     +        '       PROGRAM aborts execution',/,3X,
     +        50('='))


c    !................................! Zdd includes Pad Inertia

C     reduced impedances for frequency w are equal to:
c    !................................................!
      ZXXR=ZXX-ZXD*ZDX/ZDD
      ZXYR=ZXY-ZXD*ZDY/ZDD
      ZYXR=ZYX-ZYD*ZDX/ZDD
      ZYYR=ZYY-ZYD*ZDY/ZDD

C .....................................................................!
C     reduced TILT-PAD BEARING force coefficients at freq w are equal to:
c    !.................................................................!

        KXXDT=KXXDT+DREAL(ZXXR)
        KYYDT=KYYDT+DREAL(ZYYR)
        KXYDT=KXYDT+DREAL(ZXYR)
        KYXDT=KYXDT+DREAL(ZYXR)

      IF (FREQ.NE.ZERO) THEN
        CXXDT=CXXDT+IMAG(ZXXR)/FREQ
        CYYDT=CYYDT+IMAG(ZYYR)/FREQ
        CXYDT=CXYDT+IMAG(ZXYR)/FREQ
        CYXDT=CYXDT+IMAG(ZYXR)/FREQ
      END IF

c    !..................................................!

C      WRITE (6,99) ZXXR, ZYYR, ZXYR, ZYXR,
C     +             ZXXR*ZYYR-ZXYR*ZYXR

 99   FORMAT ( 3X,
     +        'Pad reduced Impedances:',/,3X,
     +        'ZXX:',(E11.5E2,'+i',E11.5E2),/,3X,
     +        'ZYY:',(E11.5E2,'+i',E11.5E2),/,3X,
     +        'ZXY:',(E11.5E2,'+i',E11.5E2),/,3X,
     +        'ZYX:',(E11.5E2,'+i',E11.5E2),3X,
     +        'DET:',(E11.5E2,'+i',E11.5E2),/,3X)


      END

C *****************************************************************************
C ANALYSIS AND PROGRAM BY DR. LUIS SAN ANDRES
C NASA PROJECT 1993, TEES 42490
C *****************************************************************************
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012