c --------------------------------------------------------------------------
c 10/4-12/95 modified torque,edgev,calcv to include recess edge
c         circumferential velocity (UREC). This is determined by continuity
c         on sub pjet.f
c---------------------------------------------------------------------
c 7/5/95 started modifications for angled injections BEARING=1 : HJB
c   SLAND, CALCU, EDGEU, FORCE, CALCV
c   INCLUDE PJET routine
c---------------------------------------------------------------------

C
C  ####     ##    #        ####    ####    ####   #       #    #  #######
C #    #   #  #   #       #    #  #       #    #  #       ##   #     #
C #       #    #  #       #        ####   #    #  #       # #  #     #
C #       ######  #       #            #  #    #  #       #  # #     #     ###
C #    #  #    #  #       #    #  #    #  #    #  #       #   ##     #     ###
C  ####   #    #  ######   ####    ####    ####   ######  #    #     #     ###

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
C
c---------------------------------------------------------------------
c 6/28/94 This VERSION is the ORIGINAL for the FOIL PAPER
c       i.e. does not use average axial pressure on foil deformation
c---------------------------------------------------------------------
c 6/10/94 added FX&YBACK for air bearing due to PBACK on FORCE
c---------------------------------------------------------------------
c 6/1/94  removed call to SURFACET since TB & TJ are constants
c         removed calls to arrays TB&TJ in SOLVE, SLAND, TEMPASYM
c         since they are not needed.
c---------------------------------------------------------------------
c 3/4/94 update for compliance calculations
c        added CALL PMINMAX on sub FORCE
c---------------------------------------------------------------------
C *****************************************************************************
C **                                                                         **
C **  Subroutine Solve                                                       **
C **                                                                         **
C **  SOLVE:  Main routine for solution algorithm.                           **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SOLVE(DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PDUMY/ PROLD(MAXNPOCK), TROLD(MAXNPOCK), QOLD(MAXNPOCK),
     +               POLD(MAXNXT,-MAXNYI:MAXNYI),
     +               UOLD(MAXNXT,-MAXNYI:MAXNYI),
     +               VOLD(MAXNXT,-MAXNYI:MAXNYI),
     +               TOLD(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /DORIFS/ DIAORIF(MAXNPOCK), CORIF(MAXNPOCK)
      COMMON /RECCOM/ L4(MAXNPOCK)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      COMMON /SOURCEA/ PRATIO, CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /SWITCH/ IPROP
      COMMON /PADS/ NPAD, NREC(MAXNPAD)

c     ............................................................
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC, Ec

      COMMON /THERMID/ ISOTH

c     ............................................................
      DOUBLE PRECISION U, V, P, PROLD, QOLD, POLD, UOLD, VOLD,
     +                 PREC, TREC, QREC, QIN, QOUT, QFACTOR,
     +                 DIAORIF, CORIF, L4,RPM,PS,PA,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIFS,SMASS, MPEPS, PREPS, MMP, SFLOW
      DOUBLE PRECISION CK,BETAK,THC,T,ALFT, UC, TC, Ec ,
     +                 TROLD,TOLD ,TLEAD,TRAIL

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,IFULL,BEARING,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT, IPROP,
     +        ITER, ITMAX, ITPMAX, ISYM, ICSTEP, SVERB, DVERB, BEEP
      INTEGER ISOTH ,NPAD, NREC

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION RHS(MAXNPOCK), RHOR(MAXNPOCK),
     +                 SS(MAXNPOCK, MAXNPOCK), SINV(MAXNPOCK, MAXNPOCK),
     +                 COR2, QORIF, DIFLOW, DPRI, PRECI, SIG,
     +                 SUM, PNEW, RHO, EMU, Deltap,zero,
     +                 CKR,BTR,THR, BKR(MAXNPOCK),
     +                 SSPI(MAXNPOCK),SSTI(MAXNPOCK)
      INTEGER DEVICE, II, JSTART, I, J, K, IPR, ITERP, ITSUM
C     REAL TIME

C ----------------------------------------------------------------------------
C --  SOLVE code                                                            --
C ----------------------------------------------------------------------------
      zero=0.0D0
      JSTART=ISYM+(ISYM-1)*NYI

      IPROP=1                      ! =1: Call MIPROPS--update fluid properties
      Deltap=0.001D0               ! Small change in dimensionless pressure
      ITERP=0                      !
      ITSUM=0                      !
      IPR=0                        !
C     TIME=SECNDS(0.0)             ! Starting time of solution
                                   !...................................
      CALL SETBC(RPM,PS,PA)                   ! Set known boundary conditions
c                                  !...................................

c .................................!...................................
      IF (ISOTH.EQ.1) THEN         ! For Isothermal solution
c    !.....................!       !.....................
         DO J=Jstart,NYI,1         !
            DO I=1,NXT             !
               T(I,J)=1.0D0        ! Uniform fluid temperature
            END DO                 !
         END DO                    !
         IF (IFULL.EQ.0) THEN      ! On Pad groove bearings
            DO I=1,NPAD            ! (IFULL=0)
              DO J=Jstart,NYI,1    ! set temperatures at pad adges =
                 TLEAD(I,J)=1.0D0  ! supply temperature.
                 TRAIL(I,J)=1.0D0  !
              END DO               !
            END DO                 !
         END IF                    !
c    !.....................!       !.....................
      END IF                       ! ISOTH=1
c    !.....................!       !.....................
C..................................! ...................................
                                   !
      CALL TEMPROPS                ! Calculates props at P nodes
C..................................! ...................................

c :::::::::::::::::::::::::::::::  ! .................................
C FOR SEALS AND PAD BEARINGS withot RECESSES go directly to solve
c		 equations on film lands.

      IF ((NPOCKET.EQ.0).OR.(BEARING.GE.2)) THEN
        CALL SLAND(DEVICE)         !
        RETURN                     !
      END IF                       !
                                   !
c :::::::::::::::::::::::::::::::  ! BEARING=1     ....................
                                   ! FOR RECESS HYDROSTATIC BEARING
 77   CONTINUE                     !...................................
c                                  ! .......................................
      CALL SETPREC                 ! Sets edge & rec pressures
                                   !
C ................................ ! ........................................
      CALL SLAND(DEVICE)           ! Solve Land Flow with given Prec
C ................................ ! ........................................
      ITSUM=ITSUM+ITER             ! Total # Iterations on land flow soln.
                                   !
                                   !
C   !..................!           !
      DO I=1, NPOCKET              !
C   !..................!           !
          SIG=1.0D0                !
          IF (PREC(I).GE.(1.0D0)) SIG=-1.0D0

          QOLD(I)=QREC(I)          ! Save original recess flow,
          PROLD(I)=PREC(I)         ! recess pressure, and
          TROLD(I)=TREC(I)         ! Recess temperature before perturbation

          CALL LOCPROPS(RHO,EMU,CKR,
     +    BTR,THR,PREC(I),TREC(I)) ! fluid properties at Precess
          RHOR(I)=RHO              !
          BKR(I)=BTR               ! Volumetric expansion coefficient (non-d)

          QORIF=CORIF(I)*DSQRT(RHO*DABS(1.0D0-PREC(I))) ! Inlet orifice flow
          DIFLOW=QREC(I)-QORIF*SIG                      ! Flow difference

          CALL LOCPROPS(RHO,EMU,CKR,                    ! Deltap=0.001 of (Ps-Pa)
     &    BTR,THR,PREC(I)+Deltap,TREC(I))               ! B = (1/rho) drho/dP
          L4(i)=(rho-rhor(I))/rhor(i)/Deltap            ! Compressibility factor

         !.........................!
          IF(QORIF.EQ.zero) QORIF=0.01*DABS(QIN)
         !..........................!

c        => small flow difference-> converges
          IF (DABS(DIFLOW).LT.(SFLOW*QORIF)) IPR=IPR+1

         !.........................! TRACE Calculations for DEBUGGING
          IF (SVERB.EQ.1) THEN     !
              WRITE (6, 500) ITERP, I, PREC(I), TREC(I),
     +              QREC(I),DIFLOW*100.0/QORIF, IPR
          END IF                   !
          IF ((DEVICE.EQ.1).AND.(DVERB.EQ.3)) THEN
              WRITE (1, 500) ITERP, I, PREC(I), TREC(I),
     +              QREC(I),DIFLOW*100.0/QORIF, IPR
          END IF                   !
         !.........................!....................................!

C   !.................!            !
      END DO                       ! i=1, npocket .............................
C   !.................!            !
                                   !
      IF ((IPR.GE.NPOCKET).AND.    ! CONVERGENCE of PREC pressures
     +    (ITERP.GE.1)) THEN       !
          GOTO 700 		   ! => RETURN
      ENDIF                        !...............................

      IF (ITERP.GE.ITPMAX) GOTO 900! =>MAX # its. is exceeded
                                   !...............................

C --------------------             ! **************************************
                                   !
  600 ITERP=ITERP+1                ! Start iterations on rec press
C                                  ! **************************************
      DO I=1, NPOCKET              !
          DO J=1, NPOCKET          !
              SS(I, J)=zero        ! Zero Jacobian matrix of flow derivatives
          END DO                   !
      END DO                       !
                                   !
      DO I=1, NXT                  !
          DO J= JSTART, NYI        !
              POLD(I, J)=P(I, J)   ! Save original pressure, velocity, and
              UOLD(I, J)=U(I, J)   ! temperature fields before perturbation
              VOLD(I, J)=V(I, J)   !
              TOLD(I, J)=T(I, J)   !
          END DO                   !
      END DO                       !
c                                  ! Props. change little during perturbation.
      IPROP=0                      ! 0==> No update of properties during
c                                  ! recess perturbation process--save time !
C  !................!              !
      DO I=1, NPOCKET              ! ......... PERTURB Recess press
C  !................!              !
          DPRI=.02D0*PREC(I)       ! Delta_P for ith-recess

          IF (DPRI.EQ.zero) DPRI=0.01D0

          PRECI=PROLD(I)+DPRI      ! Perturb ith Precess

          IF (PRECI.LT.PCAV) PRECI=PCAV

          PREC(I)=PRECI            ! ::::::::::::::::::::::::::::::::::::::::
                                   !
          CALL SETPREC             !     Reset recess pressure
                                   !     and
          CALL SLAND(DEVICE)       !     Solve FLOW ON LANDS
                                   !
          ITSUM=ITSUM+ITER         ! ::::::::::::::::::::::::::::::::::::::::
                                   !
          DO J=1, NPOCKET          ! Calculate Flow Derivatives
             SS(J, I)=((QREC(J))**2-(QOLD(J))**2)/DPRI
          END DO                   !
                                   !
          PREC(I)=PROLD(I)         ! get ready for next pocket

          SIG=1.0D0                !
          IF (QOLD(I).LT.zero) SIG=-1.0D0

          COR2=CORIF(I)*CORIF(I)   ! Orifice coefficient ^ 2

          RHS(I)=COR2*Rhor(I)*(1.0D0-L4(i)*Prold(i)*(1.0D0-Prold(i))
     +    - BKR(I)*(TREC(I)-TROLD(I)) )  ! Density (mass) change due to
     +    - Sig*Qold(i)*Qold(i)          ! temperature variation in the recess

          SSTI(I)=-COR2*Rhor(I)
     +    *BKR(I)*(TREC(I)-TROLD(I))    ! BK=volumetric expansion factor

          SSPI(I)= COR2*Rhor(i)*( 1.0D0-L4(i)*( 1.0D0-Prold(i) ) )

          DO J=1, NPOCKET          !Restore old Trecess field
             TREC(J)=TROLD(J)      !
          END DO                   !......

          DO K=1, NXT              !......
              DO J=Jstart,NYI      !
                P(K, J)=POLD(K, J) ! Restore old P, U, V, and T
                U(K, J)=UOLD(K, J) ! fields.
                V(K, J)=VOLD(K, J) !
                T(K, J)=TOLD(K, J) !
              END DO               !
          END DO                   !
C  !.................!             !----
      END DO                       ! END of recess pressure perturbations
C  !.................!             !----

c                                  !........................................!
      IPROP=1                      ! Restore value=1. Note: only during
c                                  ! the recess perturbation, no props update

C  !................ !             !
      DO I=1, NPOCKET              ! ........................................
C  !.................!             ! Prepare system of equations
          SUM=zero                 !-------!
          DO J=1, NPOCKET                  !
              SUM=SUM+SS(I, J)*PROLD(J)    ! [ S ] {Pold}
          END DO                           !
          RHS(I)=RHS(I)+SUM                !
c                                          !
          SS(i,i)=SS(i,i)+SSPI(I)+SSTI(I)  !
c                                          !
C  !.................!                     !
      END DO                       ! .......................................
C  !.................!             !

c     .....................................................

      CALL INVERSE(NPOCKET, SS, SINV) ! Get Inverse jacobian matrix
c     ..................................................... ==> supportp.f

      IPR=0                        !
C    !...............!             !
      DO I=1, NPOCKET              ! .......................................
C    !...............!             !
          PNEW=zero                !
          DO J=1, NPOCKET          !        Calculate improved
              PNEW=PNEW+SINV(I,J)*RHS(J) ! new recess pressures
          END DO                   !
          PNEW=DABS(PNEW)          !
          PREC(I)=DMAX1(PNEW,Pcav) !

          DPRI=DABS(Prec(I)-Prold(I))
          if (DPRI.LT.PREPS) IPR=IPR+1
          PRold(I)=Prec(I)         !
C    !...............!             !
      END DO
C    !...............!             !

      IF (IPR.LT.NPOCKET) IPR=0     ! Reset convergence counter
C
C    !------------------------------!
C    !                              !
          GOTO 77                   !...... START New iteration
C    !                              !
C    !------------------------------!


C                                  !
C ................................ ! Convergence on lands & recesses
  700 CONTINUE
        !.............................!==> TO DUMP file
         IF (DEVICE.eq.1) THEN        !
           WRITE (1, 801) ITERP, itsum
         END IF                       !
        !.............................!

      ITER=ITSUM                   !
      RETURN                       !
C ---------------------------------!-----------------------------

  900 WRITE (6, 1000)              ! Maximum number of its exceeded
      IF (DEVICE.eq.1) WRITE (1, 1000)
      ITER=ITSUM                   !
      RETURN                       !
C ---------------------------------!------------------------------

  500 FORMAT (' ', 'Ip:', I2, ' I:', I2, ', Prec:', E12.5E2,
     +        ', Trec:', E12.5E2, ', Qrec:', E12.5E2,
     +        ' %DifQ:', E12.5E2, '  Ipr:', I2)

  800 FORMAT (' ', 'Convergence in Iterp:', I3, '   EXEC Time:'
     +        E12.5E2, ' secs.')
  801 FORMAT (' ',79('-'),/,4X,
     + 'ZEROTHth O SOLN CONVERGENCE:', I3,
     + 'its on REC press, &', I4, 'TOTAL its on lands')


 1000 FORMAT (' ',72('-'),/,3X,'WARNING: Number of Its.(Iterp)',
     + 1X,'exceeded' ,/,3X,'RESULTS MAY BE ERRONEOUS',/,3X,
     +    'Start with a different initial guess.',/,
     + 1X,72('-'))


      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Setprec                                                     **
C **                                                                         **
C **  SETPREC:  Resets recess pressures on pockets after iteration.          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETPREC

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARRAY/ P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /SWITCH/ IPROP

      DOUBLE PRECISION RHOP,EMUP,P,T, PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 CK,BETAK,THC
      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL,
     +        INERL,INERP,ITURB,INTER,ICAV,MODEL,ISYM,ICSTEP,IPROP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION Rhor, EMur ,CKR,BTR,THR
      INTEGER IFF, IP, JSTA, JSTO, J, K, KK, Kstart
C ----------------------------------------------------------------------------
C --  SETPREC code                                                          --
C ----------------------------------------------------------------------------

      IFF= INERP
      KK = NPA -IFF
      Kstart=ISYM + (ISYM-1)*KK

C    !.............................!....................! IP loop
       DO IP=1, NPOCKET            ! Update Recess pressure, etc.
c     !.................!          ! ............................!
          JSTA=(IP-1)*NXI+NLC+IFF  ! pressure=PREC at interior
          JSTO=IP*NXI+1-IFF        ! and edge boundaries

          CALL LOCPROPS(RHOR,EMUR,CKR,   ! evaluate fluid props at
     +       BTR,THR,PREC(IP),TREC(IP))  ! Prec, Trec
          DO J=JSTA, JSTO          ! sweep axial and circ.
              DO K= Kstart, KK, 1  ! Left & Right side of HJB recess
                  P(J, K)=PREC(IP) !
                  Rhop(j,k)=Rhor   !
                  EMup(j,k)=EMur   !
                  CK(j,k)  =CKR    ! Dimensionless specific heat,
                  BETAK(j,k)=BTR   ! Volumetric expansion factor,
                  THC(j,k) =THR    ! and Thermal conductivity
              END DO               !
          END DO                   !
c       !.....................!    !
c     !.................!          !
       END DO                      !
C    !.............................!................! END OF Pocket sweep

      IF (INERP.EQ.1) RETURN       ! No edge inertia
      IF (IFULL.EQ.0) RETURN       ! pad bearing < 360 degress

      Kstart=Kstart+(-1+ISYM)

      DO K=Kstart, NPA             ! periodicity condition
          P(1, K)=P(NXT, K)        ! for cylindrical bearing
          Rhop(1,K)=Rhop(NXT,K)    ! of length 360 deg
          Emup(1,K)=Emup(NXT,K)    !
          CK(1,k)=CK(Nxt,k)        !
          BETAK(1,k)=BETAK(Nxt,k)  !
          THC(1,k)=THC(Nxt,k)      !
          T(1, k)=T(Nxt, k)        !
      END DO                       !.........................!

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Setbc                                                       **
C **                                                                         **
C **  SETBC:  Imposses known boundary conditions on closure of domain.       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETBC(RPM, PS, PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PIOPAD/ PLEAD(-MAXNYI:MAXNYI), PTRAIL(-MAXNYI:MAXNYI)

      COMMON /PreLR/ Ple, Pri, Csel, Cser
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      COMMON /LOSPAD/ LOSleadP, KLOSPad
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2
      COMMON /PARAM3/ EMU, RHO, PC, CD, DORIF, LOSXSI,ALPHA

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c    !.................................................................!
      DOUBLE PRECISION XP,XU,YP,YV,PLEAD,PTRAIL, U, V, P,
     +                 Ple, Pri, Csel, Cser, zero,LOSleadP, KLOSPad,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL,
     +        IPRuni, IPLuni, NPRcs, NPLcs, BEARING
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION a1u,b1u,c1u,a2u,b2u, Yri,Yle, DY,
     +                 a2d,b2d, Pgrooveu, Pgrooved
      INTEGER I, J, K, Kstart, Kstop,JSTART,JSTOP,Kmin,Kmax

C ----------------------------------------------------------------------------
C --  SETBC code                                                            --
C ----------------------------------------------------------------------------
C Approximate BC refer to those not actually needed on the analysis and
C calculated in the solution procedure.
C These are given here so as not to start the solution with unknown values
C ----------------------------------------------------------------------------

      Kstop=NPA                        ! recess edge pointer
      ZERO=0.0D0

c    !---------------------------!     !.......................
      IF (ISYM.EQ.1) THEN              !
c    !---------------------------!     ! for symmetric HJB

      DO J=1, NXT                      !
          V(J, 1)=ZERO                 ! Zero axial velocity at y=0
          U(J, NYI)=SWIRL              ! => Approximate only at y=LR
      END DO                           !
      Kstart = 1

c    !.................................!.................................!
        IF (Cser.eq.zero) THEN         ! No discharge seal at right plane
c    !.................................! set pressure at Y=LengthR
         IF (IPRuni.eq.1) THEN

          DO J=1, NXT
            P(J,NYI)=Pri               ! set uniform pressure at exit
          END DO

         ELSE
           CALL SETPright              ! set nonuniform pressure at exit
         END IF

c    !.................................!.................................!
        END IF			       ! No end seal at Y=LengthR
c    !.................................!.................................!

      Kmin=1
      Kmax=NYI

c    !---------------------------!     !.................................!
      ELSE                             !
c    !---------------------------!     ! for ASYMMETRIC BEARING, ISYM=0

      DO J=1, NXT                      !
          U(J,-NYI)=SWIRL              ! Approximate only
          U(J, NYI)=SWIRL              !  ''''' at both ends  Y=-LL, LR
      END DO                           !

c      !....................................! No discharge right seal
        IF ((Cser.eq.zero).AND.(MODEL.eq.2)) THEN
c      !....................................! asymmetric single recess row
         IF (IPRuni.eq.1) THEN
          DO J=1, NXT                  ! set pressure at y=LengthR
            P(J,NYI)=Pri               ! = Pexit for NOseal & SINGLE ROW
          END DO
         ELSE
            CALL SETPright             ! non uniform pressure
         END IF
c    !.................................!.................................!
        END IF			       ! NO END SEAL at Y=LENGTHR
c    !.................................!.................................!

c      !....................................! No discharge left seal
        IF ((BEARING.NE.2).AND.(Csel.eq.zero)) THEN  !
c      !....................................! asymmetric single recess row
         IF (IPLuni.eq.1) THEN
          DO J=1, NXT                  ! set uniform pressure at y=-LengthL
            P(J,-NYI)=Ple
          END DO
         ELSE
            CALL SETPleft              ! non uniform pressure
         END IF
c    !.................................!.................................!
        END IF			       ! NO END SEAL at Y=-LENGTHL
c    !.................................!.................................!


      IF (BEARING.EQ.1) THEN       !==> FOR HJB
          Kstart=-NPA
          Kmin=-NYI
          Kmax=NYI

      ELSE IF (BEARING.EQ.2) THEN  !==> FOR SEAL
          Kstart=1
          Kstop=NYI
          Kmin=Kstart

      ELSE IF (BEARING.EQ.3) THEN  !==> plain BEARING
          Kmin=-NYI
          Kmax=NYI

      END IF

c    !---------------------------!
      END IF
c    !---------------------------!     ! ISYM ..........................

      JSTART=0

      IF (NPOCKET.EQ.0) GOTO 555       !==> seal or plain bearing

c    !--------------------!            ! AT RECESS BOUNDARIES: Approximate only
      DO I=1, NPOCKET                  ! ....................................
c    !--------------------!            ! Set V=0 on X-edges, U=SWIRL on Y-edges
          JSTART=JSTART+1              ! On circ-X recess edges:
          JSTOP=JSTART+NLC-1           !
          DO K=Kstart, Kstop           !
              V(JSTART, K)=zero        ! Zero axial velocity at recess
              V(JSTOP, K)=zero         ! edges in axial direction
          END DO                       !
          JSTART=JSTOP                 !
          JSTOP=I*NXI                  !
          DO J=JSTART, JSTOP           ! Circumferential velocity U equal
              DO K=Kstart,Kstop        ! to
                  U(J, K)=SWIRL        ! Swirl velocity at Y-recess edge
              END DO                   !
              V(J,Kstop)=zero          !C
              IF (ISYM.EQ.0) V(J,Kstart)=zero  ! C??? really NEEDed ?
          END DO                       ! J between recesses
          JSTART=JSTOP                 !
c    !----------------------!          ! ---------------!
      END DO                           ! I on pockets
c    !----------------------!          ! ---------------!

c    !.......................!.........!.....................................
 555  IF (IFULL.EQ.1) THEN             ! SET PERIODICITY on single pad 360deg
c    !.......................!         !....................................

c      !......................!        !..................!
       IF (BEARING.EQ.1) THEN          ! FOR HJB:
c      !......................!        !..................!
        DO K=Kstart, Kstop             ! Last P, U & V CVs
           V(NXT, K)=zero              ! V=0 on pocket side
        END DO                         !

c      !......................!        !..................!
       ELSE IF (BEARING.EQ.2) THEN     ! FOR SEAL
c      !......................!        !..................!

        DO J=1, NXT		       ! Set entrance pressure
          P(J,0)=(PS-PATYP)/PSA        ! and inlet swirl
C##       P(J,1)=P(J,1)
          U(J,1)=SWIRL
        END DO

c      !......................!        !..................!
       ELSE IF (BEARING.EQ.3) THEN
c      !......................!        !..................!
C         NO B.C. at Y=0 PLANE
c      !......................!        !..................!
       END IF
c      !......................!        !..................!

c    !.......................!         !--------------------------!
      ELSE                   !.IFULL=0 ! SET V=0 ON EDGES OF PAD
c    !.......................!           Pad with axial grooves

        DO K=Kmin, Kmax
           V(1,K)=ZERO
           V(NXT,K)=ZERO
        END DO                         !......................................!
c                                      ! Groove RAM Pressures based on
c                                      ! NON uniform exit left&right pressures
      Yri=YP(NYI)                      !......................................!
      Yle=YP(-NYI)
      DY=Yri-Yle                       ! dim. axial length
      a1u=-4.0D0/DY/DY                 !
      b1u=-a1u*(Yri+Yle)               ! coefs. of quadratic for leading
      c1u=a1u*Yri*Yle                  ! edge pressure ram effect.

      IF (ISYM.EQ.1) THEN              ! for symmetric HJB
        a2u=0.0D0                      !......................
        a2d=0.0D0                      ! Ple=Pri
      ELSE                !............!......................
        a2u=(P(1,NYI)-P(1,-NYI))/DY    ! non-uniform side pressures
        a2d=(P(NXT,NYI)-P(NXT,-NYI))/DY! (Pri-Ple)/DY
      END IF              !............!......................
        b2u=P(1,NYI)-a2u*Yri           !
        b2d=P(NXT,NYI)-a2d*Yri         !

      DO K=Kmin,Kmax                    ! On -Ll < Y < Lr set pressures
        PGrooved=a2d*YP(K)+b2d          ! linear variation between left&right
        P(NXT,K)=PGrooved               !
        PTRAIL(K)=PGrooved              ! pressure at trailing edge (downstream)
        PGrooveu=a2u*YP(K)+b2u          !
        PGrooveu=PGrooveu+KLOSPAD*(YP(K)*(a1u*YP(K)+b1u)+c1u)
        P(1,K) =PGrooveu                ! pressure at leading edge (upstream)
        PLEAD(K)=PGrooveu               ! with RAM effect
      END DO                            !....

c    !.......................!          !.....................................!
      END IF			        ! IFULL=0 PAD BEARING
c    !.......................!          !.....................................!

c .......................................!....................................
c NOTE: leading edge RAM Pressure = KLOSPAD
c where KLOSpad= k/2 x RHO x (OmegaxR/2)^2, TYP k=0.60
c actual is Pram = KLOSpad - Rey x u^2 where u is mean velocity on groove
c according to Burton and Carper, 1967.
c .......................................!....................................

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Temprops                                                    **
C **                                                                         **
C ** Temprops: calculate props & store in arrays, temporal storage           **
C **                                                                         **
C *****************************************************************************


      SUBROUTINE TEMPROPS

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ TK(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /SWITCH/ IPROP

      DOUBLE PRECISION P, RHOP, EMUP, CK,BETAK,THC,TK
      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1,NXI,NYI,NXT,IFULL,
     +        J,K, ISYM, ICSTEP, Kstart, Kend, BEARING,IPROP

C ----------------------------------------------------------------------------
c Calculate fluid props at pressure nodes & store in arrays
c.............................................................................
c SUB LOCPROPS is on initopst.f file and directs to miprops.f program

      IF (IPROP.EQ.0) RETURN       !==> no update of fluid props.

      Kend=NYI
      Kstart=ISYM+(ISYM-1)*NYI
      IF (BEARING.EQ.2) Kstart=0   !=> SEAL

      DO K=Kstart, Kend
         DO  J=1, NXT
             CALL LOCPROPS(RHOP(J,K),EMUP(J,K),CK(J,K),
     +       BETAK(J,K),THC(J,K),P(J,K),TK(J,K))
         END DO
      END DO

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Sland                                                       **
C **                                                                         **
C **  SLAND:  Iterative solution of flow field on land-sills regions.        **
C **                                                                         **
C *****************************************************************************


      SUBROUTINE SLAND(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /DPUV/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /PIOPAD/ PLEAD(-MAXNYI:MAXNYI), PTRAIL(-MAXNYI:MAXNYI)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)
      COMMON /RECASP/ ASPEC(MAXNPOCK), /RECCOM/ L4R(MAXNPOCK)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR
      COMMON /PRELR/ Ple,Pri, Csel, Cser
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPE, HRECD
      COMMON /LOSPAD/ LOSleadP, KLOSPad
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP
      COMMON /COMPLIA/ AC, NETA, PBACK, LIFT
      COMMON /TDIFMAX/ DELTAT
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd

      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /JET/ ANGLEJ, XJET, CJET, DPJET
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)

      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /BTYPE/ BEARING
      COMMON /SWITCH/ IPROP
      COMMON /PADK/ KPAD
      COMMON /PADS/ NPAD, NREC(MAXNPAD)
      COMMON /THERMID/ ISOTH

c.........................................................................

      DOUBLE PRECISION HP, HU, HV,
     +                 DYP, DYV, SVN, SVS, U, V, P, RHOP, EMUP,
     +                 DU, DV, DVV,A,B,C,D,PLEAD,PTRAIL,
     +                 PREC,TREC, QREC, QIN, QOUT, QFACTOR, ASPEC,
     +                 L4R, CK,BETAK,THC,T ,ALFT,UC,TC,EC ,TB,TJ
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPE, HRECD,
     +                 PRECL, PRECR,TRECL, TRECR, Ple,Pri, Csel, Cser,
     +                 DXP, DXU, SUW, SUE, DYVK, DYPK, SVNK, SVSK,
     +                 HCB,HCJ,TLEAD,TRAIL ,TSUMP,TSUMPd,RHOSU,RHOSd
      DOUBLE PRECISION PGROOVE, PRAM, LOSleadP, KLOSPad,
     +                 PMAX, PMAXP, PEPS, AC, NETA, PBACK, DELTAT,
     +                 ANGLEJ, XJET, CJET, DPJET,
     +                 PRECdo, PRECup, PRjet

      INTEGER ITER, ITMAX, ITPMAX,IFULL,ISOTH, KPAD,NPAD,NREC,
     +        JUMIN, JUMAX, JUSTART, JUSTOP,
     +        JVMIN, JVMAX, JVSTART, JVSTOP,
     +        JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, LIFT,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        K, KM1, KP1, SVERB, DVERB, BEEP, ISYM, ICSTEP,
     +        LEFTBC,RIGHTBC, BEARING, COUNTP, MAXCOUNTP, IPROP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION UW, SOLD, relaxp, relaxu, zero, C1U, C2U,C1V,C2V,
     +                 RHOW,EMUW,HEDGE,ETA,RENR,ALFAR, KLOSY,relaxt,
     +                 Inflow, Factor, Psump, Cseal, Sgn ,L4,MACH2,
     +                 RHOr, EMUr, RHOem, EMUem, Pem, Pemold,
     +                 Rec ,PS(MAXNXT),CKR,BTR,THR,PEPSold, SMASSc,
     +                 DELTATmax

      INTEGER Kstart, Kend, Kstep, KVm1, Icase, Iprint, KM1U,
     +        Iconv, Inerldum, I, J, JP1, DEVICE, DUMYI, NYIP1, Kr,KV,
     +        Countpold, ITstop, ITdelta, IPROPold, CHPROP,
     +        JT,KEM1,KEM2,KEM3
C ----------------------------------------------------------------------------
C --  SLAND code                                                            --
C ----------------------------------------------------------------------------
C LEFTBC & RIGHTBC: 0 (IFULL=1) FLOW SPECIFIED
C                   1 GROOVE AT BOUNDARY, 2 RECESS AT BOUNDARY
C ............................................................................
      ITdelta=30          ! maximum allowed # of iterations to reduce error
      ITstop=ITdelta      ! ...............................................!
      DELTATmax=1.0D0/TC  ! maximum difference in temperature between iterations
c   !......................................................................!
      ZERO=0.0D0
c   !........................!                !............................!
c     CHPROPMAX                               ! evaluate props every CHPROPmax its.
c                                             ! given in params.f
c                                             !............................!
      IPROPold=IPROP                          ! IPROP=1 evaluate props
c                                             !      =0 skip calculation
c   !........................!                !............................!
      C1U=0.50                                ! Us = C1 US + C2 UP
      C2U=0.50                                ! COEFS FOR INTERPOLATION OF US
      C1V=0.50                                ! Vs = C1 VS + C2 VP
      C2V=0.50                                ! COEFS FOR INTERPOLATION OF VS
c   !........................!                !............................!

c NOTES:  ....................................................................
c Convergence for HJBS based on orifice flows,
c   """""     for SEALS and BEARINGS based on Pressure diferences
c.............................................................................


      Kstart=ISYM+(ISYM-1)*NYI

c .................................!...................................
      IF (ISOTH.EQ.1) THEN         ! Isothermal MODEL
         DO J=Kstart,NYI,1         ! ..................................
            DO I=1,NXT             !
               T(I,J)=1.0D0        ! Set Temperature = Tsupply
            END DO                 !
         END DO                    !
         IF (IFULL.EQ.0) THEN      ! Pad groove bearings
            DO I=1,NPAD            !
              DO J=Kstart,NYI,1    !
                 TLEAD(I,J)=1.0D0  !
                 TRAIL(I,J)=1.0D0  !
              END DO               !
            END DO                 !
         END IF                    !
      END IF                       !
C..................................! ...................................

c    !.........................!
      IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN
c    !.........................!  for hydrostatic bearing

       IF (ISOTH.EQ.1) THEN        !....................
         DO I=1,NPOCKET            ! Isothermal Model
            TREC(I)=1.0D0          ! Temp at recess = Tsupply
         END DO                    !
       END IF                      !....................
c     !............................!
       IF (MODEL.EQ.1) THEN        !=> FOR DOUBLE ROW HJBS
         NYIP1=NYI+1
         DO J=1, NXT
           HV(J,NYIP1)=zero
           RHOP(J,NYIP1)=zero
         END DO
       END IF
c     !........................!

c     ..................................! ANNULAR SEAL
      ELSE IF (BEARING.EQ.2) THEN
C     ..................................! IFULL=1 always (360 deg)
       Maxcountp=(NXT-2+IFULL)*(NYI-1)  ! # of pressure unknowns.
       IF (Cser.gt.zero) Maxcountp=Maxcountp+(NXT-2+IFULL)

C     ..................................!
      ELSE IF (BEARING.EQ.3) THEN
C     ..................................! HYDRODYNAMIC BEARING
       Maxcountp=(NXT-2+IFULL)*(NYI-1)
       IF (ISYM.EQ.0) Maxcountp=2*Maxcountp
c    !.........................!
      END IF
c    !.........................! BEARING TYPE

c     ...................................!..........................
      REC=REY/EC                         ! Reynolds No./Eckert No.
c     ...................................!..........................

c     ================== #################### ==========================
C     UPDATE FILM THICKNESS DUE TO PRESSURE_COMPLIANCE BEARING

      CALL FILMC     ! => on foilsubs.f

c     ================== #################### ==========================

c     ...................................!..........................
      PEPS=SFLOW
      PEPSold=PEPS                       !
      INERLdum=INERL			 !
      relaxp=alfp                        ! save nominal values of
      relaxu=alfu		         ! convergence parameters
      relaxt=alft                        !
      Iconv=0                            !
c     ...................................!..........................
      iprint=0                           ! set = 1 for TRACE in calculations
c     ...................................!..........................
      CHPROP=0                           !
      PMAX=1.0D0                         !
      ITER=0                             ! Start of iterative process

C    !------------------!                !-------------------------------
  100 ITER=ITER+1                        ! Zeroes calculation vectors
C    !------------------!                !-------------------------------
      SMASS=zero                         !
      MMP=zero                           !
      DELTAT=zero                        !
      PMAX=DMAX1(1.0D0,PMAX)             !
      PEPS=SFLOW*PMAX                    ! convergence as % of MAX Presure
      PMAX=zero                          !
      PMAXP=zero                         !
      COUNTP=0                           !

        IF (CHPROP.EQ.CHPROPMAX) THEN    !.... Time to evaluate fluid props.
           CHPROP=0			 !
           IPROP=IPROPold		 ! =1 YES, 0=NO
        ELSE                     !.......!
           CHPROP=CHPROP+1		 !
           IPROP=0                       ! 0=NO upgrade of fluid props.
        END IF                           !................................!

c    !..................!                !..................
      Icase=1                            ! Right Side of HJB
      Kstep=1                            !..................
      Psump=Pri                          ! Discharge pressure
      Cseal=Cser                         !  & end seal param.
      KLOSY=KLOSYr                       !  & entrance edge coeff.
      Sgn=1.0D0                          ! direction of integration
c    !..................!                !..................

C    !...................................!
C      IF (SVERB.EQ.3) WRITE (6, 10) ITER             !
C   10 FORMAT (' ', 'Solve on film lands iteration: ', I3)
C    !...................................!

 200  CONTINUE

      DO I=1, NXT+2                      ! .........................
          DVV(I)=zero                    ! Zero calculation vectors
          A(I)=zero                      ! for TDMA algorithm
          B(I)=zero                      !
          C(I)=zero                      !
          D(I)=zero                      !
          DU(I)=zero                     !
          DV(I)=zero                     !
      END DO                             !----!



C    !-----------------------!...........! SYMMETRIC BEARING
      IF (ISYM.EQ.1) THEN
C    !-----------------------!...........! SYMMETRIC BEARING
        Kstart=Kstep         !  =+/- 1

        IF ((NPOCKET.EQ.0).OR.(BEARING.EQ.3)) THEN ! => SYM JOURNAL BEARING
           DO J=1, NXT
C             DVV(J)=ZERO
              V(J,0)=ZERO
              V(J,1)=ZERO
           END DO
          KM1=Kstart      ! => start at Y=0 CV
          K=KM1
          DYPK=DYP(K)
          KV=K
          GOTO 650        !==> calculate solution on extended lands
        END IF

C    !-----------------------!...........! ASYMMETRIC BEARING
      ELSE
C    !-----------------------!...........! ASYMMETRIC BEARING
c     for 1st V-CV on centerline, y=0, and asymmetric BEARING
c                                        ! set initial counters
       IF (BEARING.EQ.1) THEN            !=> HJB ......
         K=Kstep                         !  1,-1
         KM1=-Kstep                      ! -1, 1
         KVM1=-2*Kstep                   ! -2,
         KP1=2*Kstep                     !  2, -2
         Dyvk=2.0D0*DYV(Kstep)           !
       ELSE IF (BEARING.EQ.2) THEN       ! => SEAL .....
         KM1=Kstep                       ! 1
         K=KM1+Kstep                     ! 2
         KP1=K+Kstep                     ! 3
         KVM1=KM1                        ! 1
         Dyvk=DYV(K)                     !
         C1V=1.0D0                       !
         C2V=0.0D0                       !
       ELSE IF (BEARING.EQ.3) THEN       !=> JOURNAL BEARING
         KM1=-Kstep                      ! -1, 1
         K=KM1-Kstep                     ! -2, 2
         KP1=K-Kstep                     ! -3, 3
         KVM1=-Kstep                     ! -1, 1
         Dyvk=DYV(K)                     !
         C1V=1.0D0                       !
         C2V=0.0D0                       !
         Sgn=-Sgn                        !
       END IF                            !.................
         SvnK=SVN(K)
         Svsk=SVS(K)

C       FIND AXIAL VELOCITY at Y=0 plane
c       LEFTBC,RIGHTBC = 0 periodicity, 1: groove, 2: recess edge
c      !..................................!
        IF (NPOCKET.GT.0) THEN
c      !..................................! Hydrostatic Bearing
        LEFTBC=1+IFULL
        RIGHTBC=2
        DO I=1, NPOCKET+1-IFULL
           IF (I.GT.NPOCKET) RIGHTBC=1+IFULL
           CALL SETLIM(I,1)               ! INTER=1
           CALL CALCV(Sgn,Klosy,Cseal,Psump,C1V,C2V,Kstep,KVm1)
           LEFTBC=RIGHTBC
        END DO
        Kstart=-Kstep
c      !..................................!
        ELSE
c      !..................................! SEALS or JOURNAL BEARINGS   ###
        LEFTBC=1-IFULL
        RIGHTBC=LEFTBC
        INTER=0
        CALL SETLIM(I,INTER)
        CALL CALCV(Sgn,Klosy,Cseal,Psump,C1V,C2V,Kstep,KVM1)
        IF (IFULL.EQ.1) V(NXT,K)=V(1,K)    ! inforce periodicity
        Kstart=-Kstep
c      !..................................!
        END IF
c      !..................................!


c      !.............................!
        IF (BEARING.EQ.2) THEN       ! => SEAL
c      !.............................!
           C1V=0.50D0
           C2V=C1V
           Kstart=1
           DO J=1, NXT
             DVV(J)=DV(J)
           END DO
           KM1=Kstart
           K=KM1+Kstep
           DYPK=DYP(K)
           KV=K
           GOTO 655 !==> INTER=0 goto extended lands
c      !.............................!
        ELSE IF (BEARING.EQ.1) THEN  ! => HJB
c      !.............................!
C##??? if npocket=0 then KM1,K,KV, DYPK = ???????
           DO J=1, NXT
              DVV(J)=DV(J)           ! Assure a unique V(y=0)
              V(J,0)=0.5D0*(V(J,1)+V(J,-1))
              V(J, 1) = V(J,0)
              V(J,-1) = V(J,0)
           END DO

c      !.............................!
        ELSE IF (BEARING.EQ.3) THEN  !=> PLAIN BEARING
c      !.............................!
            KV=2
            DYPK=DYP(1)+DYP(-1)
            SVSK=DYP(-1)/DYPK
            SVNK=DYP( 1)/DYPK        !: coeffs. for linear interpolation
            DO J=1, NXT
              DVV(J)=DV(J)
              V(J, Kstep)=V(J,KV)*SVSK+V(J,-KV)*SVNK ! at middle of first CV
            END DO

              KM1=-2*Kstep           ! -2, 2
              K=Kstep                !  1,-1
              Sgn=-Sgn
              KV=KM1
              GOTO 655  !==> INTER=0 GOTO extended lands

c      !.............................!
        END IF
c      !.............................!


C    !-----------------------!...........! ON FIRST control volume
      END IF                             ! of
C    !-----------------------!...........! ASYMMETRIC BEARING

C    !...................................!....................................
      IF (NPOCKET.EQ.0) GOTO 650         ! => treat as seals or JBs
C    !...................................!....................................

c     ========................................................................
      INTER=1                            ! ON LANDs BETWEEN RECESSES OR GROOVES
c     ========...........................!------------------------------------

C    !...................................! SET boundary conditions
      Kend=Kstep*NPA                     ! FOR RECESS (POCKETS) BEARINGS
C    !...................................! SET boundary conditions

c     !:::::::::::! determine fluid pressures and circumferential
       CALL PJET  ! velocities whitin recesses. 10/95
c     !:::::::::::! ==> pjet.f
c                                        !
      IF (IFULL.EQ.1) THEN               ! LEFT BDY: A recess
          LEFTBC=2                       ! NPOCKET RECESS EDGE
          PRECL=PRECdo(NPOCKET)          ! (downstream)
          TRECL=TREC(NPOCKET)            ! Recess temp. from Sub. FLOWT
          J=NPOCKET                      !
      ELSE                        !......!....................
          LEFTBC=1                       ! LEFT BDY: A groove (LEADING edge)
          PRECL=PLEAD(KSTEP)             ! Set Pad RAM pressure Rise
          TRECL=TLEAD(KPAD,KSTEP)        ! and temperature
          DO K=Kstart,Kend               !
            P(1,K)=PLEAD(K)              !
            P(NXT,K)=PTRAIL(K)           !
            T(1,K)=TLEAD(KPAD,K)         ! Tlead & Trail from Sub GROOVETEMP
            T(NXT,K)=TRAIL(KPAD,K)       ! on calctherm.f
          END DO                         !
      END IF                      !......!....................

c     !...............!                  !
  300  DO I=1, NPOCKET+1-IFULL           !
C     !...............!                  !
        CALL SETLIM(I, INTER)            !

        IF (ISOTH.NE.1) THEN             !
          DO JT=JPMIN, JPMAX             ! for Sub. CALCT b/w the recess land
            PS(JT)=P(JT,Kstart)          ! PS=Psouth when K=Kstart
          END DO                         !
        ENDIF                            ! NOTE:
c                                        ! Recess pressures got in sub SOLVE
        IF (I.LE.NPOCKET) THEN           !
          RIGHTBC=2                      ! RIGHT BDY: An I-th recess or pocket
          PRECR=PRECup(I)                ! (upstream edge)
          TRECR=TREC(I)                  !
        ELSE                             !......
          RIGHTBC=1                      ! RIGHT BDY: A PAD groove
          PRECR=PTRAIL(KSTEP)            ! P and T at GROOVE
          TRECR=TRAIL(KPAD,KSTEP)        !
        END IF                           !...................................

          K=0                            !
          KM1=Kstart                     ! 1 if ISYM=1, -1/+1 if ISYM=0
c       :::::::::::::::::::::::::::::::::! ...................................

  400     K=K+Kstep                      ! SWEEP AXIALLY, K varies
          KP1=K+Kstep                    !
          DYPK=DYP(K)                    ! Y size of P & U CVs
C                                        ! ...................................
          IF (LEFTBC.EQ.1) THEN          ! BDY: PAD groove at LEFT side
              PRECL=PLEAD(K)             !
              TRECL=TLEAD(KPAD,K)        !
          END IF                         !
          IF (RIGHTBC.EQ.1) THEN         ! BDY: PAD Groove at Right side
             PRECR=PTRAIL(K)             !
             TRECR=TRAIL(KPAD,K)         !
          END IF                         !
c     !..................................! ...................................
          CALL CALCU(Sgn,C1U,C2U,K)      ! Solve U-eqn on k row -> U(j, k)
c     !..................................! ...................................
          IF (K.eq.Kend)  goto 500       ! --> NEXT recess
                                         !
          KM1=K                          !
          K=KM1+Kstep                    !
          KP1=K+Kstep                    !
          DYVK=DYV(K)                    ! Y size of Vcv
          SVNK=SVN(K)                    ! weight coeffs for V-eqn
          SVSK=SVS(K)                    !
                                         !
c     !..................................! ...................................
          CALL CALCV(Sgn,KLOSY,Cseal,Psump,C1V,C2V,Kstep,KM1)
c     !..................................! ...................................
c                                        ! Solve V-eqn on k row -> V(j, k)
          KP1=K                          !
          K=KM1                          !
          KM1=K-Kstep                    !
          IF (K.eq.Kstep) KM1=Kstart
                                         !
c     !..................................! ...................................
          CALL CALCP(Sgn,K)              ! Correct P(j, k), U(j, k), V(j, k)
          IF (ISOTH.NE.1) THEN           !         and V(j, k+1) on k-row
             CALL CALCTM(PS,REC,Sgn,K)   ! Solve T-equation land b/w recesses
          ENDIF                          !
c     !..................................! ...................................

          KM1=K                          !
          GOTO 400                       !==> AXIAL sweep
c                                        !....................................
  500     LEFTBC=2                       ! SET BDY: recess edge on left side
          IF (I.LE.NPOCKET) THEN         !
             PRECL=PRECdo(I)             !
             TRECL=TREC(I)               !
          END IF                         !
          J=I                            !...................................
C    !.............!                     ! NEXT LAND/Recess, I=1,2,....
      END DO                             ! ...................................
C    !.............!                     !
C                                        ! Set velocity at y=0
C.!......................................! ASYM HJB: Set V=AVE(Vleft+Vright)
 550  IF (ISYM.eq.0) THEN                !
        DO J=1 , NXT
           V(J,0)=0.5D0*(V(J,1)+V(J,-1))
           V(J,1)=V(J,0)
           V(J,-1)=V(J,0)
        END DO
      END IF !...........................!....................................!


C     -----------------------------------!.......................................
  600 INTER=0                            ! SWEEP on extended LAND above REcesses
C     -----------------------------------!.......................................
      CALL SETLIM(I, INTER)              ! RESets limits for sweeps
      IF (IFULL.EQ.1) THEN               ! C###
          JVMAX=NXT                      !
          JVSTOP=2                       !
      END IF                             !...........................

      K=Kstep*NPAP1                      !
      KM1=Kstep*NPA                      !
      KP1=K+Kstep                        !
      DYVK=DYV(K)                        ! Y size of Vcv above recess edge
      DYPK=DYP(KM1)                      ! X size of Pcv on recess edge
      SVNK=SVN(K)                        ! weight coeffecients
      SVSK=SVS(K)                        !

      LEFTBC=1-IFULL                     ! SET BDY Type: 0= periodicity
      RIGHTBC=LEFTBC                     ! or 1=groove
c    !...................................!.......................................!
      CALL CALCV(Sgn,KLOSY,Cseal, Psump, C1V,C2V, Kstep,KM1)
c    !...................................! Solve V-eqn above rec edge V(j,npap1)
      IF (IFULL.EQ.1) THEN               !......
        V(1,KM1)=V(NXT,KM1)              ! KM1=+/- NPA ------ and
        P(1,KM1)=P(NXT,KM1)              ! Inforce periodicity conditions
        T(1,KM1)=T(NXT,KM1)              !
        V(1,K)=V(NXT,K)                  ! K=+/- NPA+1
      END IF                             !.............

C     -----------------------------------!.......................................
      INTER=1                            ! On lands between recess edges
C     -----------------------------------! K=+/- NPA ...........................
      K=Kstep*NPA                        !
      KP1=K+Kstep                        !
      KM1=K-Kstep                        !
      DYPK=DYP(K)                        ! Y size of Pcv

      IF (IFULL.EQ.1) THEN               ! SET Left boundary: RECESS
        LEFTBC=2                         !
        PRECL=PRECdo(NPOCKET)            !
        TRECL=TREC(NPOCKET)              !
      ELSE                               !.....
        LEFTBC=1                         ! SET Left boundary: GROOVE
        PRECL=PLEAD(K)                   !
        TRECL=TLEAD(KPAD,K)              !
      END IF                             !....

        RIGHTBC=2                        ! Right boundary is a recess
c    !...................!               !
      DO I=1, NPOCKET                    ! SWEEP between recesses
c    !...................!               ! ......
          PRECR=PRECup(I)                ! Pressure at right boundary : recess
          TRECR=TREC(I)                  ! Temperat at right boundary : recess
          CALL SETLIM(I, INTER)          ! Correct pressure and flow
          CALL CALCP(Sgn,K)              ! at edge of inter-recess lands correct
c                                        ! U(j,npa),P(j,npa),V(j,npa),V(j,npap1)
          IF (ISOTH.NE.1) THEN           !
             CALL CALCTM(PS,REC,Sgn,K)   ! find land Temperature  K=+/- NPA
          ENDIF                          !
c                                        !
          PRECL=PRECdo(I)                !
          TRECL=TRECR                    !
          LEFTBC=2                       !
c    !...................!               ! ......
      END DO                             ! I=1, NPOCKET
c    !...................!               ! ......

      IF (IFULL.EQ.1) THEN               ! FOR 360 deg hydrostatic bearing
        V(NXT,K)=V(1,K)                  ! Inforce periodicity
        U(NXT,K)=U(1,K)                  ! K=+/- NPA
        P(NXT,K)=P(1,K)                  !
        T(NXT,K)=T(1,K)                  !
        V(NXT,KP1)=V(1,KP1)              !
      ELSE                               !....................................!
        PRECL=PRECdo(NPOCKET)            ! Sweep on last land from recess to
        TRECL=TREC(NPOCKET)              ! PAD groove:
        RIGHTBC=1                        !
        PRECR=PTRAIL(K)                  !
        TRECR= TRAIL(KPAD,K)             !
        CALL SETLIM(NPOCKET+1,INTER)     !
        CALL CALCP(Sgn,K)                ! Correct Pressure & Veloc. Fields
        IF (ISOTH.NE.1) THEN             !
           CALL CALCTM(PS,REC,Sgn,K)     ! find land temperature K=+/- NPA
        ENDIF                            !
      END IF                             !....................................!

      DO J=1, NXT                        ! Save pressure coeffs for extended
          DVV(J)=DV(J)                   !  land correction
      END DO                             !....................................!

      C1U=0.25D0                         ! ......... coefficients for
      C2U=0.75D0                         !           Us interpoltion AT K=NPA+1

      Kstart=Kstep*NPA                   ! resets axial counters
      KM1=Kstep*NPA                      !
      K=KM1+Kstep                        !
      DYPK=DYP(K)                        !
      KV=K                               !

C   !------------------------------------!-----------------------------------!
  650 INTER=0                            ! On extended lands above recesses
C   !------------------------------------!-----------------------------------!
      LEFTBC=1-IFULL                     ! BDY COND: =1 a groove
      RIGHTBC=LEFTBC                     ! =0 periodicity
      CALL SETLIM(I,INTER)               ! RESET circumferential counters

C----!
 655  Kend=Kstep*NYI                     ! == axial exit plane
C----!
C        !...............................! for thermal analysis
          IF (ISOTH.NE.1) THEN           !
          DO JT=JPMIN, JPMAX             ! for CALCT on 360-deg film land
            PS(JT)=P(JT,K  )*SVS(K)      ! PS=Psouth
     +            +P(JT,KM1)*SVN(K)      !
          END DO                         !
          ENDIF                          !
C        !...............................!

 700  CONTINUE

c    !.........................!         !
      IF ( K.eq.Kend ) THEN              ! AT Exit plane of Bearing
c    !.........................!         !

c       !-----------------------------------------------!
C        FOR TWO recess ROW BEARING (right side)
c       !-----------------------------------------------!
         IF ((Icase.eq.1).AND.(MODEL.eq.1) ) THEN
            KP1=NYIP1
            DO J=1, NXT
             V(J,KP1)=-V(J,K)            ! Symmetry BC.
             U(J,KP1)=zero               ! dumy value to set AN UN = zero
            END DO

c         !--------------------------!
            CALL CALCU(Sgn,C1U,C2U,K)
c         !--------------------------!

            DO J=1, NXT
               DV(J)=zero                ! To set ANP=0
               V(J,KP1)=zero             !        Fnp=zero
            END DO
            DV(NXT+1)=zero
            DV(NXT+2)=zero

            IF (IFULL.EQ.1) U(JUSTOP,K)=U(JUMIN,K)! Periodicity condition

c          !--------------------------!
            CALL CALCP(Sgn,K)
c          !--------------------------!

         ELSE
c       !---------------------------------------!
C        FOR SINGLE recess ROW BEARING OR SEAL
c       !---------------------------------------!
            KP1=Kend                            ! MODEL=2 = single row HJB
c          !--------------------------!
            CALL CALCU(Sgn,C1U,C2U,K)
c          !--------------------------!
            IF (IFULL.EQ.1) U(JUSTOP,K)=U(JUMIN,K) ! Periodicity condition

         END IF

         GOTO 800              !==> update inlet values

c    !.........................!
      ELSE                     ! |K| < |KEND|
c    !.........................!  AT internal control volumes on extended lands

         KP1=K+Kstep
c       !--------------------------!     !.................................!
         CALL CALCU(Sgn,C1U,C2U,KV)      ! Solve U-eqn on K row -> U(j, k)
c       !--------------------------!     !.................................!
      IF (KPAD.EQ.0) KPAD=NPAD
      IF (IFULL.EQ.1) THEN               ! 360 deg bearing:
          U(JUSTOP,K)=U(JUMIN,K)         ! Periodicity condition
      ELSE                               ! --
          P(1,K)=PLEAD(K)                ! pad bearing: take upstream (leading edge)
          T(1,K)=TLEAD(KPAD,K)           ! groove values
      END IF                             !
c    !.........................!         !
      END IF                   !         ! K=Kstart, ... Kend-1
c    !.........................!         !

      C1U=0.5D0
      C2U=C1U
      KM1U = KM1

      KM1=K                              !
      K=KM1+Kstep                        !
      DYVK=DYV(K)                        ! Y size at V cv &
      SVNK=SVN(K)                        ! Weight coeffs.
      SVSK=SVS(K)                        !
c    !.........................!         !
      IF (K.eq.Kend) THEN      !         ! K=NPa+1, ... Kend-1
c    !.........................!         !

         IF ((Icase.eq.1).AND.(MODEL.eq.1) ) THEN  ! Right Side of 2 row HJB
            KP1=NYIP1
            DO J=1, NXT
             V(J,KP1)=-V(J,K)            ! Symmetry BC.
            END DO
         ELSE
            KP1=Kend                     ! Model=2: single row HJB
         END IF

c    !.........................!
      ELSE
c    !.........................!
          KP1=K+Kstep                    !
c    !.........................!         !K=NPa+1, ... Kend-1
      END IF                   !         !
c    !.........................!         !

c                                        ! Solve V-eqn on k row -> V(j, k)
c       !--------------------------!     !.................................!
      CALL CALCV(Sgn,KLOSY,Cseal,PSump,C1V,C2V,Kstep,KM1)
c       !--------------------------!     !.................................!

      IF (IFULL.EQ.1) V(NXT,K)=V(1,K)    ! and inforce periodicity

      KP1=K                              !
      K=KM1                              !
      KM1=KM1U ! K-Kstep                 !

c                                        ! Correct P(j,k), U(j,k), V(j,kv), V(j,kp1)
c       !--------------------------!     !.................................!
      CALL CALCP(Sgn,KV)
c       !--------------------------!     !.................................!
      IF (ISOTH.NE.1) THEN               !
         CALL  CALCT(PS,REC,Sgn,KV)      ! Solve T-equation: T(J,K)
      ENDIF                              !
c       !--------------------------!     !.................................!


      KM1=K                !===> sweep next K-row
      K=K+Kstep
      KV=K
      DYPK=DYP(K)

      C1V=0.5D0
      C2V=C1V

c  !::::::::::::::::::::::!              !
      GOTO 700                           !=> NEXT K axial row
C                                        !
C :::::::::::::::::::::::::::::::::::::: !::::::::::::::::::::::::::::::::::
C AFTER solving land flow equations:

  800 CONTINUE
c     FOR THERMAL MODEL: SET EXIT Temperature and PROPERTIES
      IF (ISOTH.NE.1) THEN               ! Exit temperatures (K=NYI)
        KEM1=Kend-Kstep                  ! Calculated by the zero-
        KEM2=KEM1-Kstep                  ! temperature gradient
        KEM3=KEM2-Kstep                  ! at discharge ends
        DO J=1,NXT                       !
           T(J,Kend)=(3.D0*T(J,KEM1)-    !
     +                 T(J,KEM2))/2.D0   ! From Roach's Book,1976,P189-191
        END DO                           !

       IF (IPROP.EQ.1) THEN              ! Set Properties at Exit plane
        DO J=1,NXT                       !
         CALL LOCPROPS(RHOP(J,Kend), EMUP(J,Kend),CK(J,Kend),
     +   BETAK(J,Kend),THC(J,Kend),P(J,Kend),T(J,Kend))
        END DO                           !
       ENDIF                             !
      ENDIF                              ! ISOTH <> 1
C ...................................... ! ----------------------------

  810 CONTINUE

c    !...................................!
      IF ((INERP.EQ.0).OR.(NPOCKET.EQ.0)) goto 850
c    !...................................!
C     FOR HJBS: Pressure CORRECTIONS AT RECESS EDGES: CIRCUMFERENTIAL
C                                        ! If INERP=1 for HJB then:
      IF (IFULL.EQ.1) THEN               ! specify pressure at boundaries
        LEFTBC=2                         !
        PRECL=PRECdo(NPOCKET)            ! LEFTBC=2 a pocket at left edge
        TRECL=TREC(NPOCKET)              !
        J=NPOCKET                        !
        ASPE=ASPEC(J)                    !
      ELSE                               !...............................
        LEFTBC=1                         ! LEFTBC=1 a groove at left edge
        J=1                              !
      END IF                             !...............................

c    !...................................!
      DO I=1, NPOCKET+1-IFULL            ! Modify edge pressures & inlet Us
C    !........................!          !.................................!
          CALL SETLIM(I, 1)              ! SET counters
c                                        ! ....
          IF (I.LE.NPOCKET) THEN         ! SET BDY conditions
             RIGHTBC=2                   ! =2 a recess edge at right edge
             PRECR=PRECup(I)             !
             TRECR=TREC(I)               !
          ELSE                           ! ..
             RIGHTBC=1                   ! =1 a groove at right edge
          END IF                         !.................................!

c        !...............................!
          DO K=Kstep,Kstep*(NPA-1),Kstep ! AXIAL SWEEP
c        !...............................!
              IF (LEFTBC.EQ.1) THEN      ! LEFT BDY: a groove
                 PRECL=PLEAD(K)          ! pad leading edge
                 TRECL=TLEAD(KPAD,K)     !
              END IF                     !.....
              IF (RIGHTBC.EQ.1) THEN     ! RIGHT BDY: a groove
                 PRECR=PTRAIL(K)         ! pad trailing edge
                 TRECR=TRAIL(KPAD,K)     !
              END IF                     !
c            !-------------------!       !.............................
              CALL EDGEU(Sgn,J,K)        ! Find DeltaP at recess edge

c        !...--------------------........!.............................
          END DO                         !sweep K=...-> +/- (NPA-1)
c        !...............................!

          LEFTBC=2                       ! LEFT BDY: A recess edge
          PRECL=PRECdo(I)                !
          TRECL=TRECR                    !
          IF (I.LE.NPOCKET)ASPE=ASPEC(I) !
          J=I                            !............................!

C    !........................!          !.................................!
      END DO                             ! END SWEEP on POCKET edges
C    !........................!          !.1, NPOCKET+1-IFULL..............!


  850 CONTINUE

c    !...................................!
      IF ((IFULL.EQ.1).AND.(NPOCKET.GT.0)) THEN
c    !...................................!
      DO K=Kstep,Kstep*NPA,Kstep         ! Inforce periodicity condition
          V(NXT, K)=V(1, K)              ! along left & right cut boundary
          U(NXT, K)=U(1, K)              ! FOR 360 deg HJB
          P(NXT, K)=P(1, K)              !
          T(NXT, K)=T(1, K)              !
          Rhop(NXT,K)=Rhop(1,K)          !
          Emup(NXT,K)=Emup(1,K)          !
          CK(NXT, K)=CK(1, K)            !
          THC(NXT, K)=THC(1, K)          !
          BETAK(NXT, K)=BETAK(1, K)      !
      END DO                             !
c    !...................................! POR PAD bearing
      ELSE IF (IFULL.EQ.0) THEN          ! advect downstream last U vel
c    !...................................!
       DO K=Kstep, Kstep*NYI, Kstep      !
         U(NXT,K)=U(JUMAX,K)             !
         T(NXT,K)=TRAIL(KPAD,K)          ! Trailing edge temperatures
       END DO                            !
       P(1,Kstep*NPA)=PLEAD(Kstep*NPA)   !  ###????
       T(1,Kstep*NPA)=TLEAD(KPAD,Kstep*NPA)
c    !...................................!
      END IF
c   !....................................!................................!


 860  CONTINUE   ! ==> GET READY to modify entrance AXIAL (edge) pressures

c    !.....................................!
      IF (BEARING.EQ.1) THEN               ! => HJB RECESS EDGES
c    !.....................................!
        SVNK=SVN(Kstep*NPAP1)
        SVSK=SVS(Kstep*NPAP1)
        DYPK=DYP(KSTEP*NPA)
c    !.....................................!
      ELSE IF (BEARING.EQ.2) THEN          ! => SEAL INLET
c    !.....................................!
        SVNK=SVN(2)
        SVSK=SVS(2)
        DYPK=DYP(1)
c    !.....................................!
      ELSE IF (BEARING.EQ.3) THEN          ! => BEARING
c    !.....................................!
        GOTO 880
c    !.....................................!
      END IF
c   !......................................! Modify Pedge & Vinlet

 870  CONTINUE

c   !....................................!................................!
      CALL EDGEV(Sgn,DYPK,SVNK,SVSK,KLOSY,INERP,Kstep)
c   !....................................!................................!

c                                        ! Modify P(NYI) for end seal cond.
 880  IF (Cseal.GT.zero) THEN
        CALL ENDSEAL(Sgn,DYP(Kstep*NYI),Psump,Cseal,alfp,Kstep)
      END IF
c   !....................................!................................!

c    !============================!.............................!
      IF (AC.GT.ZERO) THEN        ! COMPLAINT SURFACE BEARING
c    !............................!.............................!
c     AT exit side of bearing Y=L, we need to provide values of film
c     thickness due to pressure deformation,
c     Assume pad leading edge is fixed (i.e. it does not deform)
c            pad trailing edge deforms with upstream value.
c    !............................! at bearing exit plane,  Y=Lright
c##   KU=NYI-1                    ! foil deformation should occur
      KEND=Kstep*NYI              ! but since p=0 here, we take
      HP(1,KEND)=HV(1,KEND)       ! the first value of film thickness
      DO J=1, NXT-1               ! with deformation
        JP1=J+1                   !
        HP(JP1,KEND)=HV(JP1,KEND) !
        HU(J,KEND)=SUE(J)*HP(J,KEND)+SUW(J)*HP(JP1,KEND)
      END DO                      !
        HU(NXT,KEND)=HP(NXT,KEND) !
c    !............................!
      END IF                      ! for bearing compliance effects
c    !============================!.............................!

C   !.......................!............!......................!
c     Direct for symmetry/asymmetry on bearing
C   !.......................!............!......................!

C   !.......................!
      IF (ISYM.eq.0) THEN   !            !=> FOR asymmetric HJB
c   !.......................!            !......................!
         Factor=1.0D0
         IF (BEARING.EQ.2) GOTO 895  !=> SEAL

         IF (BEARING.EQ.3) THEN      !=> BEARING at Y=0
            KV=2
            DYPK=DYP(1)+DYP(-1)
            SVSK=DYP(-1)/DYPK
            SVNK=DYP( 1)/DYPK  ! coeffs. for linear interpolation
            DO j=1, NXT
               P(j,-Kstep)=P(j,Kstep)   ! SET unique fields
               T(j,-Kstep)=T(j,Kstep)   ! SET unique fields
               U(j,-Kstep)=U(j,Kstep)   ! at Y=0
               V(j, Kstep)=V(j,KV)*SVSK+V(J,-KV)*SVNK ! at middle of first CV
               V(j,-Kstep)=V(j, Kstep)
            END DO
         END IF

         IF (Icase.eq.2) GOTO 890

         Sgn=-1.0D0
         Kstep=-1
         Psump=Ple                       ! LEFT side of HJB
         Cseal=Csel                      ! Discharge pressures & end seal
         KLOSY=KLOSYl                    ! & entrance edge coeff.
         Icase=2
         GOTO 200
C   !.......................!
       ELSE                              ! symmetric bearing ISYM=1
C   !.......................!
         Factor=2.0D0
         GOTO 895
C   !.......................!
       END IF
c   !.......................!


c    !..................................! ASYMMETRIC BEARING=3
 890  IF (BEARING.EQ.3) THEN
c    !..................................! unify Flow fields
        DO j=1, NXT
          P(j, 0)=P(j,1)
          T(j, 0)=T(j,1)
          U(j, 0)=U(j,1)
          V(j, 0)=V(j,1)
          CALL LOCPROPS(Rhop(j,0),Emup(j,0),CK(J,0),
     +    BETAK(J,0),THC(J,0),P(j,0),T(J,0))
          Rhop(j, 1) = Rhop(j,0)
          Rhop(j,-1) = Rhop(j,0)
          Emup(j, 1) = Emup(j,0)
          Emup(j,-1) = Emup(j,0)
          CK(j, 1)   =   CK(j,0)
          BETAK(j,1) =BETAK(j,0)
          THC(j, 1)  =  THC(j,0)
          CK(j,-1)   =   CK(j,0)
          BETAK(j,-1)=BETAK(j,0)
          THC(j,-1)  =  THC(j,0)
        END DO
        GOTO 895
c    !..................................!
      END IF
c    !..................................!

c   !....................................! HYDROSTATIC BEARING AT center line:
 892    DO j=1, NXT
          P(j,0)=0.5D0*(P(j,1)+P(J,-1))
          T(j,0)=0.5D0*(T(j,1)+T(J,-1))
          U(j,0)=0.5D0*(U(j,1)+U(J,-1))
          V(j,0)=0.5D0*(V(j,1)+V(J,-1))
          IF (IPROP.EQ.1) THEN
            CALL LOCPROPS(RHOP(J,0),EMUP(J,0),CK(J,0),
     +      BETAK(J,0),THC(J,0),P(J,0),T(J,0))
          END IF
        END DO

c   !-----------------------!            !..............................

 895  CONTINUE                           ! Calculates :
c   !-----------------------!            !..............................
      IF (ISOTH.NE.1) THEN               !
         CALL TORQUER                    ! TorquexU over recesses for DTr
      END IF                             ! &
c   !.......................!  ..........! flows on recesses
      CALL FLOWT(REC,Inflow)             ! and recess temperatures as well
c   !-----------------------!            ! as side flows ............


C........................................!.........................

C      if (iprint) then                   ! REMOVE comment lines to have
C      call Pfields(DEVICE)               ! TRACE on calculations
C      print *, 'enter iprint'
C      read *, iprint
C      if (iprint.eq.99) STOP
C      end if                             !^^^^^^^^^^^^^^^^^^ TRACE

C    !...................................! ==> PRINT CONVERGENCE
      IF (SVERB.EQ.1) THEN               !     OF SOLUTION
          WRITE (6, 1200) ITER,SMASS,QIN,QOUT,
     +                    PMAX,PMAXP,Countp,MAXCOUNTP

      ELSE IF (SVERB.EQ.3) THEN          !
          DUMYI=0                        !
          CALL PRINTUVP(DUMYI, 1)        !
      END IF                             !

      IF ((DEVICE.eq.1).AND.(DVERB.GE.1)) THEN
              WRITE (1, 1200) ITER, SMASS, QIN, QOUT,
     +                    PMAX,PMAXP,Countp,MAXCOUNTP
         IF (DVERB.EQ.3) CALL PRINTUVP(DEVICE, 2)
      END IF                             !
C    !...................................!............................!
                                         !
c    !-----------------------------------!---------------------------!

C##   IF ((ISOTH.EQ.1).OR.(ISOTH.EQ.-1)  ) GOTO 900

C     ISOTH=0, Isothermal stator and journal
C     ISOTH=2, Adiabatic journal, Isothermal stator
C     ISOTH=3, Adiabatic stator, Isothermal journal
C     ISOTH=1, no thermal solution, ISOTH=-1, adiabatic solution
c    !-----------------------------------!---------------------------!
C##   CALL SURFACET                      ! Bush and Journal Surface Temperatures
C     This call is not needed for constant surface temperatures.
c    !-----------------------------------!---------------------------!

 900  CONTINUE

C    !...................................!......................
C     CHECKS CONVERGENCE FOR BEARINGS
C    !...................................!......................
      IF (ITER.EQ.1) THEN                !
            SOLD=SMASS                   !
            GOTO 100                     ! 2 its minimum
C    !...................................!......................
      ELSE IF (ITER.GE.ITMAX) THEN       !
            GOTO 1000                    ! => NO CONVERGENCE
C    !...................................!......................
      ELSE IF (ITER.EQ.10) THEN          !
            SMASSc=SMASS                 ! save mass error
C    !...................................!.......................
      ELSE IF (ITER.EQ.ITstop) THEN      ! C##
        IF (DABS(SMASS/SMASSc).GT.10) GOTO 1600
            SMASSc=SMASS                 !
            ITstop=ITstop+ITdelta        !=> divergence occurs.
      END IF
C    !...................................!......................

C     MASS error has increased 100 times !!!!
c     ----------------------------------- ABORT calculations
      IF (SOLD.EQ.0) SOLD=1
      IF ((SMASS/SOLD).GT.1.0D2) GOTO 1600      !==> divergence   !###????

c    !-----------------------------------!..........................!
      IF (ISOTH.NE.1) THEN               ! Convergence on Thermal field
c    !....................!              ! requires aat least 2 its. on
c                                        ! properties
        IF (ITER.LT.2*CHPROPMAX+1) GOTO 955
        IF (DELTAT.GT.DELTATmax) GOTO 955
      END IF                             ! 1 deg K max. difference
c    !-----------------------------------!..........................!

c    !.......................!
      IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN
c    !.......................!           ! => HJB convergence in lands
       IF (SMASS.LE.MPEPS*DABS(Inflow)/Factor) GOTO 1100
c
c     Here it is based on % of total mass flow = MPEPS
c    !.......................!
      ELSE
c    !.......................!           !=> SEAL or PLAIN  BEARING

       IF (COUNTP.GE.Maxcountp) GOTO 1100 ! Convergence on Pressure

c    !.......................!
      END IF
c    !.......................!

 955  IF (Iconv.GE.3) THEN               !==> error grows
        IF((SMASS-SOLD).GT.ZERO) GOTO 1570
      ELSE                               !    reduce relax params.
        Iconv=Iconv+1                    !..........
      END IF                             !
c                                        !
      SOLD=SMASS                         ! save last values
      COUNTPold=COUNTP                   !
C    !...................................!......................
                                         !

C    !...................................!......................
C   !==========!
 977   GOTO 100                           !->  NEW iteration
C   !==========!



c !......................................!......................
 1570 Iconv=0                            !
      IF (Alfp.LE.0.25D0) THEN           ! Reduce relaxation params.
        Alfp=0.25D0                      ! Min values are
      ELSE                               ! alfp=alfu=0.25
        Alfp=Alfp**1.25                  !..........................!
        Iconv=1
      END IF

      IF (Alfu.LE.0.25D0) THEN
        Alfu=0.25D0
      ELSE
        Alfu=Alfu**1.25
        Iconv=2
      END IF
      Betu=1.0D0-Alfu

      IF (Iconv.GE.1) THEN
           ITER=ITER-1
           Iconv=0
      END IF

C#      write(6,*)  'Changed relax Params: Ap & Au:', Alfp, Alfu,
C#     +             ' on iter:', iter
      GOTO 100                           ! -> START OVER

c !......................................!......................
c NO CONVERGENCE>                        !
c !......................................!......................
 1600 CONTINUE
      write (6,*) '$ Source mass error has increased 10 fold on last'
      write (6,*) '$ ', ITdelta, ' its, ==> PROGRAM aborts'
      write (6,*) 'Smass:',Smass, '  Sold:', Sold
      GOTO 1100
c........................................!..........................

 1000 WRITE (6, *) 'WARNING:  MAXIMUM NUMBER OF ITERATIONS ON LANDS EXCE
     +EDED., ITMAX:', ITER
c........................................!.............................!
c
c !......................................!......................
c CONVERGENCE>                           !
c !......................................!......................
 1100 CONTINUE                           ! Convergence or ITER>ITMAX
c........................................!...........................      IF ((Alfp.LE.0.25D0).OR.(Alfu.LE.0.25D0)) THEN
        write(6,*) 'Relax parameters AP & AU too low = 0.25,',
     +             '==> RESULTS may be in error'
      END IF
c........................................!...........................
      alfp=relaxp                        ! reset under-relax values
      alfu=relaxu                        !.....
      alft=relaxt			 !
      betu=1.0D0-Alfu                    !
      IPROP=IPROPold		         !

C##      IF (ISOTH.NE.1) THEN
C##        print *, 'SLAND================================='
C##        print *, 'MAX Tdifference:', DELTAT*TC, ' degK'
C##        print *, 'ITER:', ITER
C##      END IF

      RETURN
c........................................!...........................

 1200 FORMAT (' ','It:',I3,' SUM(Ms):',E10.3E2,
     +        ' Qi:',E9.2E2,' Qo:',E9.2E2, ' PMAX:',
     +        E9.2E2,' Ppmax:',E9.2E2,' #Ps:',I3,
     +        '=>',I3)

 1300 FORMAT (' ', 'Convergence achieved in ', I3, ' iterations.')


C.........................................!.........................

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Setlim                                                      **
C **                                                                         **
C **  SETLIM:  Sets lower and upper limits for control volumes.              **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETLIM(I, INTER)

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL

      INTEGER JUMIN, JUMAX, JUSTART, JUSTOP
      INTEGER JVMIN, JVMAX, JVSTART, JVSTOP
      INTEGER JPMIN, JPMAX, JPSTART, JPSTOP
      INTEGER NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      INTEGER I, INTER
C ----------------------------------------------------------------------------
C --  SETLIM code                                                           --
C ----------------------------------------------------------------------------

      IF (INTER.EQ.1) THEN        ! Inter=1, on INTER_RECESS lands
C                                 !................................
          JVMIN=(I-1)*NXI+2       ! min: leftmost point
          JVMAX=JVMIN+NLC-3       ! max: rightmost point
          JVSTART=JVMIN-1         ! start: loc. left BC
          JVSTOP=JVMAX+1          ! stop: loc. right BC
          JUMIN=JVMIN-1           !
          JUMAX=JVMAX             !
          JUSTOP=JVSTOP           !
          IF (I.EQ.1) THEN        ! first land/recess
             IF (IFULL.EQ.1) THEN
                 JUSTART=NXT-1
             ELSE
                 JUSTART=JUMIN
             END IF
          ELSE
             JUSTART=JUMIN-1
          END IF
C###      IF (I.GT.NPOCKET) JUSTOP=JUMAX  !### check here for i=npocket+1

          JPMIN=JVMIN             ! = JUMIN+1
          JPMAX=JVMAX             ! = JUMAX
          JPSTART=JVSTART         ! = JUSTART+1
          JPSTOP=JVSTOP           ! = JUSTOP
C                                 !.................................
      ELSE                        ! Inter=0, on extended lands
C                                 !.................................
       IF (IFULL.EQ.1) THEN       !=> 360 deg pad
          JUMIN=1                 ! min: leftmost point
          JUMAX=NXT-1             ! max: rightmost point
          JUSTART=NXT-1           ! start: loc. left BC
          JUSTOP=NXT              ! stop: loc. right BC
          JPMIN=JUMIN             !
          JPMAX=JUMAX             !
          JPSTART=JUSTART         !
          JPSTOP=JUSTOP           !
          JVMIN=1                 !
          JVSTART=NXT-1           !
          JVMAX=JUMAX             !
          JVSTOP=JUSTOP           !
      ELSE                        !=> pad < 360deg
          JUSTART=1
          JUMIN=1
          JUMAX=NXT-1
          JUSTOP=NXT-1
          JPMIN=2
          JPMAX=NXT-1
          JPSTART=1
          JPSTOP=NXT
          JVMIN=JPMIN
          JVMAX=NXT-1
          JVSTART=1
          JVSTOP=NXT
      END IF
C                                 !
      END IF                      !.................................

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcu                                                       **
C **                                                                         **
C **  CALCU:  Finds U velocity from marching equation:                       **
C **    Apu Up = Aep Ue + Awu Uw + Asu Us + Dpu (Pp-Pe) + Spu + (1-a)Apu Up  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCU(Dir, C1, C2, KV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DPUV/ BU(MAXNXTP2), BV(MAXNXTP2), BVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR

      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /LRBOUND/ LEFTBC,RIGHTBC
      COMMON /BTYPE/ BEARING, /SWITCH/ IPROP
c........................................................................

      DOUBLE PRECISION DXP, DXU,SUW, SUE, HP, HU, HV,
     +                 U, V, P, RHP, EMP,
     +                 BU, BV, BVV, A, B, C, D
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD,
     +                 DYVK, DYPK, SVNK, SVSK, PRECL,PRECR,TRECL,TRECR
      INTEGER JUMIN,JUMAX,JUSTART,JUSTOP,JPMIN,JPMAX,JPSTART,JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,IFULL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        K, KM1, KP1, ISYM, ICSTEP, LEFTBC, RIGHTBC, BEARING,IPROP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYPRE, SANB, UW, UE, FWU, AWUO, AEU,
     +                 VS, VN, VP, KXP, SPU, SPS, DPU, FEU,
     +                 FSU, FNU, ANU, AWU, ASU, AEUC, CORRE,
     +                 KXJ,KXB, AWUC, CORRW, ZERO, FWUC, FEUC,
     +                 FNKX, C1, C2, US, Dir,
     +                 RHOP,EMUP,RHON,EMUN,RHOS,EMUS,RHOW,EMUW,
     +                 RHOE,EMUE,RHSE,EMSE,RHNE,EMNE,SUEJ,SUWJ,
     +                 PP,HEDGE,ETA,AX,RENH, HS, HN,
     +                 RHOR,EMUR,CKR,BTR,THR
      INTEGER J, JJ, JJJ, JJMAX, JP1, INERLDUM, ICASE,KK,KV
C ----------------------------------------------------------------------------
C--  CALCU code                                                            --
C ----------------------------------------------------------------------------
c Dir: flow direction, IF Dir=1.0   RIGHT SIDE of HJB
c                            =-1.0  LEFT SIDE of HJB
C ............................................................................
C KV=K except when K=+/- 1 for asymmetric bearing=3

      KK=IABS(K)
      JJJ=2-JUMIN                                !
      JJMAX=JUMAX+JJJ
      DYPRE=DYPK*REY                             ! Dypk: Size of Ucv
      INERLDUM=INERL                             ! :::::::::::::::::::::::::::

      J=JUMIN
      ZERO=0.0D0
      SANB=ZERO
C                          ! ------------------- ! Boundary conditions
  100 IF (INTER.EQ.1) THEN                       ! On inter-recess lands:
c    !.....................!                     !
          RHOE=Sue(J)*Rhp(J,K)+Suw(J)*Rhp(J+1,K)
          FWU=DYPRE*RHOE*HU(J,K)*U(J,K)          ! =Feu at left edge
          RHOW=Rhp(J,K)
          UW=FWU/DYPRE/RHOW/HP(J,K)              !west U vel. at left edge
          FWUC=FWU                               !
C                                                !
          RHOW=Sue(JUMAX)*Rhp(JUMAX,K)+Suw(JUMAX)*Rhp(JPSTOP,K)
          FEUC=DYPRE*RHOW*HU(JUMAX,K)*U(JUMAX,K) ! = Fwu at right edge
          RHOE=Rhp(JPSTOP,K)
          UE=FEUC/DYPRE/RHOE/HP(JPSTOP,K)        !
          IF ((KK.EQ.NPA).AND.(LEFTBC.EQ.2)) INERL=0  !NO EDGE inertia
C     !....................!                     !-----------
       ELSE                                      ! BCs On extended lands:
c     !....................!                     !....
          UW=U(JUSTART, K)                       ! Left U velocity=BC
          UE=U(JUSTOP, K)                        ! Right U velocity=BC
          RHOW=Rhp(JUMIN,K)
          FWU=DYPRE*RHOW*HP(JUMIN,K)*(UW+U(JUMIN,K))/2.0D0
c                                                ! west flow on Ucv
C     !.....................!                    !
       END IF                                    ! --------------------------
C     !.....................!
      FWU=FWU*INERL                              ! West flow on Ucv
      AWUO=DMAX1(FWU,ZERO)                       !
      AEU=ZERO                                   !
C :::::::::::::::::::::::::::::::::::::::::::::: ! --------------------------

      ICASE=1
      J=JUMIN-1

C    !....................!                      !
 701  J=J+1

      IF (J.EQ.JUMAX) GOTO 702

      GOTO 800

C    !....................!                      !

C                                                ! on last Ucv of edge k=Npa
 702  IF ((KK.EQ.NPA).AND.(RIGHTBC.EQ.2)) THEN
C    !....................!                      ! neglect EDGE fluid inertia
          INERL=0                                ! and correct at volume
          D(JJ)=D(JJ)-Betu*B(JJ)*U(J-1,K)        !
          B(JJ)=B(JJ)+A(JJ)/ALFU                 !    interface j=jumax-1
          BU(J-1)=B(JJ)-A(JJ)-SANB               !
          A(JJ)=ZERO                             !
          D(JJ)=D(JJ)+Betu*B(JJ)*U(J-1,K)        !
          J=JUMAX                                !
      END IF                                     !
C    !....................!                      !
      ICASE=2                                    !
C     .....................!                     !

 800  JP1=J+1                                    ! Sweep from left -> right
      JJ=J+JJJ                                   ! J+2-JUMIN
      SUEJ=SUE(J)                                !
      SUWJ=SUW(J)                                !

      VS=V(J,KV)*SUEJ+V(JP1,KV)*SUWJ             ! South vel=Vp for upwind
      VN=V(J,KP1)*SUEJ+V(JP1,KP1)*SUWJ           ! North vel=Vp for upwind
      VP=(VS+VN)/2.0D0                           !

      Rhop=Suej*Rhp(J,K)+Suwj*Rhp(JP1,K)
      Emup=Suej*Emp(J,K)+Suwj*Emp(JP1,K)

      KXP=FNKX(U(J,K),VP,HU(J,K),SPEED,REP,RHOP,
     +         EMUP,KXJ,KXB)                     ! Shear parameter
      SPU=KXP*DXU(J)*DYPK/HU(J,K)                ! Source term, friction
      SPS=(KXJ/KXP)*SPU                          !
      DPU=HU(J,K)*DYPK                           !
      D(JJ)=DPU*(P(J,K)-P(JP1,K))+SPS*MSPEED     ! = Dpv (Ps-Pa)+U/2 Sps

      B(JJ)=SPU/ALFU                             ! Spu/Alfu
      C(JJ)=ZERO                                 !
      A(JJ)=ZERO                                 !
      BU(J)=B(JJ)

C     !....................!                     ! Account for Inertia
      IF (INERL.EQ.1) THEN                       ! Effects on lands:
C     !....................!

          AWU=DMAX1(+FWU,ZERO)                   !At west Ucv face

c      !..........................................................!
          RHOE=Rhp(JP1,K)                        ! EAST Ucv face
          FEU=DYPRE*RHOE*HP(JP1,K)*(U(J,K)+U(JP1,K))/2.0D0
          AEU=DMAX1(-FEU,ZERO)

          HS=HV(J,KV)*SUEJ+HV(JP1,KV)*SUWJ       !At south Ucv face
          RHOS=(SUEJ*(Rhp(J,KM1)+Rhp(J,K))+
     +          SUWJ*(Rhp(JP1,KM1)+Rhp(JP1,K)))/2.0D0
          FSU=Dir*REY*RHOS*DXU(J)*HS*VS
          ASU=DMAX1(+FSU, ZERO)

          HN=HV(J,KP1)*SUEJ+HV(JP1,KP1)*SUWJ     !At north Ucv face
          RHON=(SUEJ*(Rhp(J,K)+Rhp(J,KP1))+
     +          SUWJ*(Rhp(JP1,K)+Rhp(JP1,KP1)))/2.0D0
          FNU=Dir*REY*RHON*DXU(J)*HN*VN          !
          ANU=DMAX1(-FNU, ZERO)                  !

          SANB=AEU+AWU+ASU+ANU                   ! SUM(Anb)
          C(JJ)=-AWU                             ! Coef. of Uw
          A(JJ)=-AEU                             ! Coef. of Ue
          B(JJ)=B(JJ)+SANB/ALFU                  ! Coef. of Up = Apu/Alfu

          US=C1*U(J, KM1)+C2*U(J, K)             !
          D(JJ)=D(JJ)+ASU*US+ANU*U(J,KP1)        ! RHS coef.
          FWU=FEU                                !
          BU(J)=B(JJ)-SANB
      ELSE
          INERL=INERLDUM                         ! RESET value
      END IF                                     ! Inerl=1
C    !......................!                    !                           !

      D(JJ)=D(JJ)+BETU*B(JJ)*U(J, K)             ! RHS+Underrelax, Betu=1-Alfu


      GOTO (701,200), ICASE                      !

C  ..............................................!.......................
                                                 !
  200 IF ( (INTER.EQ.0). OR .
     +     ((INERP+INERL).EQ.0). OR.
     +     ( IABS(K).EQ.NPA ) ) GOTO 500         ! => SOLVE U EQUATION

C  !.............................................! CORRECT rec edge flows
C                                                !
C 210  IF (RIGHTBC.EQ.2) THEN
      JJ=JJMAX                                   ! ---------------- Last Ucv
      AEUC=DMAX1(-FEUC,ZERO)                     !
      CORRE=(AEUC-AEU)/ALFU                      !
      B(JJ)=B(JJ)+CORRE                          ! coef. of Up
      D(JJ)=D(JJ)+BETU*CORRE*U(JUMAX, K)         !
      BU(JUMAX)=BU(JUMAX)+BETU*CORRE             ! for SIMPLEC procedure
C      END IF

C 220  IF (LEFTBC.EQ.2) THEN
      JJ=2 ! jumin+jjj                           ! ---------------- First Ucv
      AWUC=DMAX1(FWUC,ZERO)                      ! Upwind coef.
      CORRW=(AWUC-AWUO)/ALFU                     ! correction of Uw
      B(JJ)=B(JJ)+CORRW                          ! coef. of Up
      D(JJ)=D(JJ)+BETU*CORRW*U(JUMIN, K)         !
      BU(JUMIN)=BU(JUMIN)+BETU*CORRW             ! for SIMPLEC procedure
C      END IF

C                                                ! ::::::::::::::::::::::::::
C                                                !
C ! -------------------------------------------- ! On inter-recess lands set
C                                                ! Bernoulli effects on recess
  300 IF (INERP.EQ.0) GOTO 400                   ! ..........................

      IF ((UW.GT.ZERO).AND.(LEFTBC.EQ.2)) THEN   ! Flow out of pocket at
          JJ=JUMIN+JJJ                           ! left boundary with recess
          HEDGE=HRECD+HP(JUMIN,K)                !  JPSTART=JUMIN
          ETA=HP(JUMIN,K)/HEDGE                  !
          RHOW=Rhp(JUMIN,K)                      !
          CALL LOCPROPS(RHOR,EMUR,CKR,           ! Rhoe-, Pe-
     +                 BTR,THR,PRECL,TRECL)      ! Te- = Trecl
          AX=KLOSXd*RHOW*(1.0D0-(ETA*RHOW/RHOR)**2)

          DPU=DYPK*HU(JUMIN,K)                   ! .....
          AWUC=AWUC-AX*DPU*UW                    ! Edge inertia effect on Uw
          D(JJ)=D(JJ)+DPU*(PRECL-P(JUMIN,K))     ! Effect of Uinf velocity
          BU(JUMIN)=BU(JUMIN)+AX*DPU*UW          ! for SIMPLEC procedure

      END IF                                     ! ..........................


      IF ((UE.LT.ZERO).AND.(RIGHTBC.EQ.2)) THEN  ! Flow out of pocket at
          JJ=JJMAX                               ! right boundary with recess
          HEDGE=HRECD+HP(JPSTOP,K)               ! JPSTOP=JUSTOP
          ETA=HP(JPSTOP,K)/HEDGE                 !
          RHOE=Rhp(JPSTOP,K)
          CALL LOCPROPS(RHOR,EMUR,CKR,           ! Rhoe-, Pe-
     +                  BTR,THR,PRECR,TRECR)     ! Te- = Trecr
          AX=KLOSXu*RHOE*(1.0D0-(ETA*RHOE/RHOr)**2)

          DPU=DYPK*HU(JUMAX,K)                   !
          AEUC=AEUC+AX*DPU*UE                    ! Edge inertia effect on Ue
          D(JJ)=D(JJ)+DPU*(P(JPSTOP,K)-PRECR)    ! Effect of Uinf velocity
          BU(JUMAX)=BU(JUMAX)-AX*DPU*UE          ! for SIMPLEC procedure

      END IF                                     ! ..........................
                                                 ! k=Npa:Corner> Use prev value
C                                                ! ::::::::::::::::::::::::::
C  !.............................................!
C                                                !
  400 A(JJMAX)=-AEUC                             ! Corrected east coeff.
      C(2)=-AWUC                                 ! Corrected west coeff.

C                                                ! ::::::::::::::::::::::::::
C .............................................. ! Solution of U-equations by
  500 CALL TDMA(UW, UE, JJMAX)                   ! TDMA algorithm
C .............................................. ! ..........................
C                                                !
  600 DO J=JUMIN, JUMAX                          ! B(jj)=Usol(j)
          JJ=J+JJJ                               !
          U(J, K)=B(JJ)                          ! Order solution & set
          BU(J)=HU(J,K)*DYPK/BU(J)               ! coeffics. for Peqn.
      END DO                                     !
C................................................!...........................
      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcv                                                       **
C **                                                                         **
C **  CALCV:  Finds V velocity from marching equation:                       **
C **       Apu Vp = Aep Ve + Awu Vw + Asu Vs + Dpu (Ps-Pp) + (1-a) Apu Vp*   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCV(Dir,KLOSY, Cseal, Pexit, C1,C2,Kstep, KVM1 )

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT),DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /DPUV/ BU(MAXNXTP2), BV(MAXNXTP2), BVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /URECJET/UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK

c     ................................................................
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

c     ................................................................

      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 U, V, P, RHP, EMP,
     +                 BU, BV, BVV, A, B, C, D,
     +                 PRECdo, PRECup, PRjet, UREC
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC,ASPEC,HRECD,
     +                 DYVK, DYPK, SVNK, SVSK

      INTEGER JVMIN, JVMAX, JVSTART, JVSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,IFULL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1,NXI,NYI,NXT,
     +        K, KM1, KP1, ISYM, ICSTEP, BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYVRE, SANB, VW, VE, UW, FWV, UE, UP,
     +                 KYP, SPV, DPV, FEV, FSV, AEV, AWV, ASV,
     +                 ANV, DYP, DIV, FEW, VS, DDXH, FSVO, ASVO,
     +                 CORRS, FORCEP, ZERO, FNV, FACTOR,
     +                 FNKY, Cseal, Pexit,HW,HE, ETA,AY,VN, KLOSY, Dir,
     +                 RHOP,EMUP,RHOS,EMUS,RHON,EMUN,RHOW,EMUW,
     +                 RHSW,EMSW,RHSE,EMSE,RHOE,EMUE,C1,C2,
     +                 CKR,BTR,THR,Uedge

      INTEGER J,JJ,JJJ,JJMAX,JP1,I,IJ,JL,JR,JMAX,KK,KKP1,ICASE,Kstep,
     +        KVm1,II

C ----------------------------------------------------------------------------
C --  CALCV code                                                            --
C ----------------------------------------------------------------------------
c Dir: flow director: = 1.0 Right side ob HJB, 0<y -> <L/D
c                     =-1.0 Left side of HJB, 0: y -> -L/D
c............................................................................

      JJJ=2-JVMIN                                 !
      JJMAX=JVMAX+JJJ                             !
      DYVRE=DYVK*REY                              ! Dyvk = Size of Vcv
      SANB=0.0D0                                  !
      ZERO=0.0D0                                  !
C                                                 ! .........................
 100  VW=V(JVSTART, K)                            ! Left V velocity = BC
      VE=V(JVSTOP, K)                             ! right V velocity = BC
C                                                 !
      UW=U(JVSTART, K)*SVSK+U(JVSTART, KM1)*SVNK  ! West velocity on Vcv
      HW=HU(JVSTART,K)*SVSK+HU(JVSTART,KM1)*SVNK  !

      J=JVSTART
      RHOW=SVSK*(Rhp(J,K)*Sue(J)  +Rhp(JVMIN,K)*Suw(J))+
     +     SVNK*(Rhp(J,KM1)*Sue(J)+Rhp(JVMIN,KM1)*Suw(J))
      FWV=DYVRE*RHOW*HW*UW*INERL                  ! West flow on Vcv

C                                                 ! .........................
C ......................                          ! On standard V equation
      DO J=JVMIN, JVMAX                           ! .........................
C ......................                          !
          JJ=J+JJJ                                ! SWEEP from left -> right
          C(JJ)=ZERO                              ! on V control volume
          A(JJ)=ZERO                              !
          UE=U(J, K)*SVSK+U(J, KM1)*SVNK          ! East velocity
          UP=(UE+UW)/2.0D0                        ! U velocity at center Vcv
          RHOP=Rhp(J,K)*SVSK+Rhp(J,KM1)*SVNK
          EMUP=Emp(J,K)*SVSK+Emp(J,KM1)*SVNK

          KYP=FNKY(UP,V(J,K),HV(J,K),SPEED,REP,RHOP,EMUP)

          SPV=KYP*DXP(J)*DYVK/HV(J,K)             ! Source term, friction
          DPV=HV(J,K)*DXP(J)                      ! Area of pressure action

          B(JJ)=SPV/ALFV                          ! Alfv=Underrelaxation par
          D(JJ)=Dir*DPV*(P(J,KM1)-P(J,K))         ! = Dpv (Ps-Pp)
          BV(J)=B(JJ)
C                                                 !
C .............................                   ! Account for inertia
          IF (INERL.EQ.1) THEN                    ! Effects on lands:
C .............................                   !.................
                                                  !ON Vcv south face
              VS=C2*V(J,K)+C1*V(J,KVM1)           !
              RHOS=Rhp(J,KM1)
              FSV=Dir*RHOS*REY*DXP(J)*HP(J,KM1)*VS
              ASV=DMAX1(+FSV,ZERO)                !


              VN=(V(J,K)+V(J,KP1))/2.0D0          !ON Vcv north face
              RHON=Rhp(J,K)
              FNV=Dir*RHON*REY*DXP(J)*HP(J,K)*VN
              ANV=DMAX1(-FNV,ZERO)                !

                                                  !......
              HE=HU(J,K)*SVSK+HU(J,KM1)*SVNK      ! ON Vcv east face

              IF (J.EQ.JVMAX) THEN
                   JP1=JVSTOP
              ELSE
                   JP1=J+1
              END IF

              RHOE=SVSK*(Rhp(J,K)*Sue(J)+Rhp(JP1,K)*Suw(J))+
     +             SVNK*(Rhp(J,KM1)*Sue(J)+Rhp(JP1,KM1)*Suw(J))
              FEV=DYVRE*UE*HE*RHOE                !
              AEV=DMAX1(-FEV,ZERO)                !
                                                  !
              AWV=DMAX1(+FWV,ZERO)                !

              SANB=AEV+AWV+ASV+ANV                ! SUM(Anb)

              C(JJ)=-AWV                          ! coef. of Vw
              A(JJ)=-AEV                          ! coef. of Ve
              B(JJ)=B(JJ)+SANB/ALFV               ! coef. of Up = Apv/Alfv
              D(JJ)=D(JJ)+ASV*V(J,KVM1)+ANV*V(J,KP1) ! RHS coef=Dpv(Ps-Pe)+Asv Vs

              FWV=FEV                             !
              BV(J)=B(JJ)-SANB                    !=> for SIMPLEc
c.............................                    !.........
          END IF                                  ! Inerl=1
C .............................                   !

          D(JJ)=D(JJ)+BETV*B(JJ)*V(J, K)          ! RHS+Underrelax, Betv=1-Alfv

          UW=UE                                   !
C .............................                   !
      END DO                                      ! :::::::::::::::::::::::::
C .............................                   !



C    !................................................!
  200 IF (((INERP+INERL).EQ.0).OR.(BEARING.GE.3)) THEN
C    !................................................!
C     No edge inertia effects or PLAIN BEARING
c
                 GOTO 777 ! -> SOLVE equations

C FOR HYDROSTATIC BEARINGS:

C    !................................................!
      ELSE                                            ! FOR HYDROSTATIC BEARING
     +IF ((BEARING.EQ.1).AND.(IABS(K).EQ.NPAP1).AND.  !--------------------------
     +    (NPOCKET.GT.0)) THEN                        ! MODIFY Eqns when REC edge
C    !................................................!is at bottom of Vcv

c          KK=Kstep*NPA = KM1
c          KKP1=KK+Kstep = K

          DYP=DYPK/2.0D0                          ! Y size of Pcv at edge
          DIV=2.0D0                               !
C        !................                        !SWEEP ACROSS RECESSES
          DO I=1, NPOCKET                         !:::::::::::::::::::::
C        !................                        !
              II=0                                !
              JMAX=I*NXI+1                        !
 611          J=(I-1)*NXI+NLC                     ! On left corner of recess
              II=II+1                             !
              Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KM1))
              FEW=DYPK*Uedge*(HU(J,KM1)-HP(J,KM1))
              ICASE=1
              GOTO 614


 612          J=J+1                               ! at interior V control vs.
              IF(J.EQ.JMAX) GOTO 613
              II=II+1
              Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KM1))
              FEW=DYP*Uedge*(HU(J,KM1)-HU(J-1,KM1))
              DIV=1.0D0
              ICASE=2
              GOTO 614

 613          DIV=2.0D0                            ! At right corner of recess
              II=II+1
              Uedge=UREC(I,II)*(1.0D0+HRECD/HP(J,KM1))
              FEW=DYPK*Uedge*(HP(J,KM1)-HU(J-1,KM1))
              ICASE=3

 614          JJ=J+JJJ                            ! Correct Flow & Coefs.
              RHOS=Rhp(J,KM1)                     !
              FEW=FEW*RHOS                        !

              VS=(V(J,KM1)+V(J,K))/2.0D0          ! WRONG Vsouth * Flow
              FSVO=Dir*RHOS*DXP(J)*HP(J,KM1)*VS   !
              ASVO=DMAX1(FSVO*INERL*REY,ZERO)     !

              RHON=Rhp(j,k)*Svsk+Rhp(j,km1)*Svnk
              FNV=RHON*DXP(J)*HV(J,K)*V(J,K)      ! North Vflow on 1/2Pcv

              FSV=FNV+Dir*FEW                     ! CORRECT South flow
              VS=FSV/RHOS/DXP(J)/HP(J,KM1)        !
              ASV=DMAX1(REY*FSV*Dir/DIV,ZERO)     !

              CORRS=(ASV-ASVO)/ALFV               ! Correction on south coefs.
              B(JJ)=B(JJ)+CORRS                   ! node
              BV(J)=BV(J)+CORRS*BETV              ! Simplec procedure
              FORCEP=DXP(J)*HV(J,K)*Dir*(P(J,KM1)-P(J,K))  ! = Dpv (Ps-Pp)
              D(JJ)=BETV*B(JJ)*V(J,K)             ! RHS side = Under-relax V

c         !.........................................!.....................
           IF ((INERP.EQ.1).AND.((Dir*VS).GT.ZERO)) THEN
c         !.........................................!
                  RHOP=RHP(JMAX-1,KM1-Kstep)        !
                  ETA=HP(J,KM1)/(HP(J,KM1)+HRECD)   ! IF inertia at edge.
                  AY=KLOSY*RHOS*(1.0D0-(ETA*RHOS/RHOP)**2)
                  FACTOR=AY*DXP(J)*HV(J,K)*VS*Dir   !
                  ASV=ASV-FACTOR                    !
                  BV(J)=BV(J)+FACTOR                !
                  FORCEP=DXP(J)*HV(J,K)*Dir*(PRjet(I,II)-P(J,K)) !Dpv (Prec-Pp)
           END IF                                   ! Inerp=1
c         !.........................................!.....................
              VN=(V(J,K)+V(J,KP1))/2.0D0            !

              FNV=Dir*REY*Rhp(J,K)*DXP(J)*HP(J,K)*VN*INERL
              ANV=DMAX1(-FNV,ZERO)                !

              D(JJ)=D(JJ)+ASV*VS/DIV+FORCEP+ANV*V(J,KP1)

          GOTO (612, 612, 615) , ICASE            !

 615      CONTINUE                                !

C      !.......................!                  !------
          END DO                                  ! i=1, npocket -----------
C      !.......................!                  !------

        IF (IFULL.EQ.1) THEN
          B(2)=B(JJMAX)                           ! for 1st Vcv copy from
          D(2)=D(JJMAX)                           ! last CV; j=2=Jvmin+jjj
          BV(JVMIN)=BV(JVMAX)                     ! jvmin=1; jvmax=nxt
        END IF
        GOTO 777

C      END IF                                     ! k=Npa+1
C  !-------------------------!                    ! FOR HYDROSTATIC BEARIN



C FOR SEALS:

C   !.............................................!
      ELSE IF ((BEARING.EQ.2).AND.(K.EQ.2)) THEN  !Modify eqns. at SEAL
C   !.............................................!inlet

      DO J=JVMIN, JVMAX                           !
c   !..........................!                  !..........................

              JJ=J+JJJ                            ! Correct Flow & Coefs.
              Vs=V(j,1)                           !

              RHOS=RHP(j,1)                       !
              FSV=Rey*RHOS*Dxp(J)*HP(J,1)*Vs      ! south flow
              FSVO=FSV*INERL                      !
              ASV=DMAX1(FSV,ZERO)                 !
              ASVO=DMAX1(Fsvo,ZERO)               !
              CORRS=(ASV-ASVO)/ALFV               ! Correction on south coefs.
              B(JJ)=B(JJ)+CORRS                   ! node
              BV(J)=BV(J)+CORRS*BETV              ! Simplec procedure
              FORCEP=DXP(J)*HV(J,K)*(P(J,1)-P(J,K))  ! = Dpv (Ps-Pp)
              D(JJ)=BETV*B(JJ)*V(J,1)             ! RHS side = Under-relax V

         IF ((INERP.eq.1).AND.(VS.GT.ZERO)) THEN  !............. Correct:
                  AY=KLOSY*RHOS                   !
                  FACTOR=AY*DXP(J)*HV(J,K)*VS     !
                  ASV=ASV-FACTOR                  !
                  BV(J)=BV(J)+FACTOR              !
                  FORCEP=DXP(J)*HV(J,K)*(P(J,0)-P(J,K)) ! Dpv (Prec-Pp)
         END IF                                   ! Inerp=1

              VN=(V(J,K)+V(J,KP1))/2.0D0          !
              FNV=REY*RHP(j,k)*DXP(J)*HP(J,K)*VN*INERL
              ANV=DMAX1(-FNV,ZERO)                !
              D(JJ)=D(JJ)+ASV*VS+FORCEP+ANV*V(J,KP1)

          END DO                                  !
C      !.......................!                  !------


C   !.............................................!AT SEAL
       END IF
C   !.............................................!inlet




c ................................................!............................!
c                                                 ! Modify eqn. for end seal
 777  CONTINUE

      IF ((Cseal.GT.ZERO).AND.((K*Kstep).EQ.NYI)) THEN
c    !................................................!

      DO J=JVMIN, JVMAX
         JJ=J+JJJ
         DPV=HV(J,K)*DXP(J)
         VN=(V(J,K)+V(J,KP1))/2.00D0
         IF ((Dir*VN).gt.zero) THEN
           Factor=DPV*Cseal*Rhp(J,K)*VN*Dir
           D(JJ)=D(JJ)+DPV*Dir*(P(J,K)-Pexit)-Factor*VN
C##        BV(J)=BV(J)+Factor
         END IF
      END DO

      END IF
c    !..............!   Cseal>0 & K=Nyi

C                                                 !
C ............................................... ! Solution of V-equations by
  350 CALL TDMA(VW, VE, JJMAX)                    ! TDMA algorithm
C ............................................... ! .........................
  400 DO J=JVMIN, JVMAX                           ! Order V solution
          JJ=J+JJJ                                !
          V(J, K)=B(JJ)                           !
          BV(J)=HV(J,K)*DXP(J)/BV(J)              !
      END DO                                      !..........................

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcp                                                       **
C **                                                                         **
C **  CALCP:  Solves Pres Correc Eqn at inter-recess lands, Inter=1:         **
C **          Ap Pp = Ae Pe' + Aw Pw' + As Ps' + Mp                          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCP(Dir,KV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HOFILM/ HPO(MAXNXT,-MAXNYI:MAXNYI),
     +                HUO(MAXNXT,-MAXNYI:MAXNYI),
     +                HVO(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +                BT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DPUV/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT

      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /BTYPE/ BEARING, /SWITCH/ IPROP

      DOUBLE PRECISION DXP, DXU, HP, HU, HV, SUW, SUE, HPO, HUO, HVO,
     +                 DYP,DYV,SVN,SVS, U, V, P, RHP,EMP,
     +                 DU, DV, DVV,A, B, C, D
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 DYVK, DYPK, SVNK, SVSK,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 PMAX, PMAXP,PEPS, AC, ETA, PBACK,
     +                 CK,BT,THC,T

      INTEGER JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,IFULL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        K, KM1, KP1, LEFTBC, RIGHTBC, BEARING, COUNTP,
     +        MAXCOUNTP, IPROP, LIFT

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION FWP,AWP,DPU,DPV,FEP,FSP, PP, PU, PV,
     +                 FNP, MP,AEP, ASP, ANP, PW, PE, PPJK, FSNP,
     +                 RHOW,EMUW,RHOE,EMUE,RHOS,EMUS,RHON,EMUN,
     +                 SVNKM1,SVSKM1, DMP, Dir, Dirk,PPMAG,PABS
      INTEGER J, JJ, JJJ, JM1, JP1, JJMAX, KV

C ----------------------------------------------------------------------------
C --  CALCP code                                                            --
C ----------------------------------------------------------------------------
c Dir: flow director: =1.0 Right side ob HJB, 0<y -> <L/D
c                     =-1.0 Left side of HJB, 0: y -> -L/D
c............................................................................
C KV=K except when K=+/- 1 for asymmetric bearing=3

      JJJ=2-JPMIN                              !
      JJMAX=JPMAX+JJJ                          !
      JJ=JPMIN+JJJ                             !
      J=JPSTART                                !

      RHOW=Rhp(J,K)*Sue(J)+Rhp(JPMIN,K)*Suw(J)
      DPU=HU(J,K)*DYPK*RHOW                    !
      FWP=DPU*U(J, K)                          ! Flow on west side of Pcv

      IF (LEFTBC.EQ.0) THEN       !.........................!
         AWP=0.0D0                ! SPEC: FLOW
      ELSE                        !.........................!
         AWP=DPU*DU(J)            ! SPEC: PRESSURE
      END IF                      !.........................!

      Dirk=Dir
      IF ((BEARING.EQ.1).AND.(IABS(K).EQ.1)) Dirk=-Dir

      SVNKM1=SVN(K)                            !
      SVSKM1=SVS(K)                            !

C ......................                       !
      DO J=JPMIN, JPMAX                        !----------------------------
C ......................                       !
          JJ=J+JJJ                             ! Sweep from left-right on
                                               ! P cvs

          RHOE=Rhp(J,K)*Sue(J)+Rhp(J+1,K)*Suw(J)
          DPU=HU(J,K)*DYPK*RHOE                ! east U flow on Pcv
          FEP=DPU*U(J,K)                       !
          AEP=DPU*DU(J)                        !

          RHON=Rhp(J,KP1)*SVSK+Rhp(J,K)*SVNK   !
          DPV=HV(J,KP1)*DXP(J)*RHON            !
          FNP=Dir*DPV*V(J,KP1)                 ! North V flow on Pcv
          ANP=DPV*DV(J)                        !

          RHOS=Rhp(J,K)*SVSKM1+Rhp(J,KM1)*SVNKM1
          DPV=HV(J,KV)*DXP(J)*RHOS             !
          FSP=Dir*DPV*V(J,KV)                  ! South V flow on Pcv
          ASP=DPV*DVV(J)                       !

          MP=FWP-FEP+FSP-FNP                   ! Mass source on Pcv

          C(JJ)=-AWP                           ! coeffs. for TDMA soln.
          A(JJ)=-AEP                           !
          B(JJ)=AEP+AWP+ASP+ANP                !
          D(JJ)=MP                             ! P'(j, k-1)=0

          FWP=FEP                              !
          AWP=AEP                              !
          DMP=DABS(MP)                         !
          SMASS=SMASS+DMP                      ! Global mass source
          IF (DMP.GT.MMP) MMP=DMP              ! MAX(Local Mass source)
C
C......................                        !

      END DO                                   !
C ......................                       ! ::::::::::::::::::::::::::::
C                                              !
      IF (RIGHTBC.EQ.0) THEN                   ! specified flow on right Boundary
          B(JJMAX)=B(JJMAX)+A(JJMAX)           ! AEP=0
          A(JJMAX)=0.0D0                       !
      END IF                                   !.............................
                                               !
C                                              ! ............................
  100 PW=0.0D0 ! =Pcor(Jpstart)                ! Spec. pressure at left B
      PE=0.0D0 ! =Pcor(Jpstop)                 ! Spec. pressure at right B
C                                              ! ----------------------------
  200 CALL TDMA(PW, PE, JJMAX)                 ! Solve by TDMA algorithm
C                                              ! ----------------------------
      IF ((INTER.EQ.0).AND.(IFULL.EQ.1)) THEN  !
          B(JJMAX+1)=B(2)                      ! P'(Nxt)=P'(1)
      END IF                                   ! where Pcor(j)=B(jj)
C ............................................ !
C                                              ! ::::::::::::::::::::::::::::
  300 DO J=JPMIN, JPMAX                        ! CORRECT:
C    !....................!                    ! .......
          JJ=J+JJJ                             !
          PPJK=B(JJ)                           ! P'(j, k)
          P(J, K)=P(J, K)+ALFP*PPJK            ! Corrected P(j, k),
          U(J, K)=U(J, K)+(PPJK-B(JJ+1))*DU(J) ! Corrected U(j, k) velocity

          PPMAG=DABS(PPJK)                     !......
          IF (PPMAG.LE.PEPS) COUNTP=COUNTP+1   ! check for convergence
          IF (PPMAG.GT.PMAXP) PMAXP=PPMAG      ! on pressures
          PABS=DABS(P(J,K))                    !
          IF (PABS.GT.PMAX) PMAX=PABS          !......

          V(J,KV)=V(J,KV)-Dirk*PPJK*DVV(J)     !     "     south V velocity
          V(J,KP1)=V(J,KP1)+Dir*PPJK*DV(J)     !     "     north V velocity
          DVV(J)=DV(J)                         ! for use on next k

C    !....................!                    !
      END DO                                   ! sweep X-direction
C    !....................!                    !.
                                              ! ::::::::::::::::::::::::::::

C
C ............................................ ! ----------------------------

      IF (IPROP.EQ.1) THEN                     !
        DO J=JPMIN, JPMAX                      !----------------------------
         CALL LOCPROPS(RHP(J,K),EMP(J,K),CK(J,K),
     +   BT(J,K),THC(J,K),P(J,K),T(J,K))       ! Update of properties
        END DO                                 ! at the INTERIOR P-nodes.
      END IF                                   ! Properties at the edges
C ......................                       ! updated in EDGEU & EDGEV

      IF (LEFTBC.GT.0) THEN                    ! IF P SPECIFIED
          J=JPSTART                            ! Correct Ue velocity
          U(J, K)=U(J, K)-B(2)*DU(J)           ! Jpstart=Jumin
      END IF                                   ! ............................

      IF ((IFULL.EQ.1).AND.(INTER.EQ.0)) THEN  !................ 360 deg bearing
          P(JPSTOP, K  ) = P(JPMIN, K)         ! On extended lands
          V(JPSTOP, KV ) = V(JPMIN, KV)        ! jpstop=nxt, jpmin=1
          U(JPSTOP, K  ) = U(JPMIN, K)         ! Inforce periodicity
          V(JPSTOP, KP1) = V(JPMIN, KP1)       !
          Rhp(JPSTOP, K) = Rhp(JPMIN, K)       !
          Emp(JPSTOP, K) = Emp(JPMIN, K)       !
          ck(jpstop ,k)=ck(jpmin ,k)           !
          bt(jpstop ,k)=bt(jpmin ,k)           !
          thc(jpstop,k)=thc(jpmin,k)           !
      END IF                                   !........................

C .............................................!..............................
C     UPDATE FILM THICKNESS FOR COMPLIANT SURFACE BEARING
C .............................................!..............................
C  NO deformation occurs if P<Pback pressure of elastic matrix
C .............................................!......................c
c    !......................!                  !....................!
      IF (AC.GT.0.0D0) THEN !                  ! Calculate deformation
c    !......................!                  !....................!
C
      JM1=JPSTART
c    !....................!
      DO J=JPMIN, JPMAX
c    !....................! Sweep along X axis
        PP=P(J,K)         !##     DMAX1(P(J,K),PCAV)
       IF (PP.GT.PBACK) THEN
        HP(J,K)=HPO(J,K)+AC*(P(J,K)-PBACK)     ! at P - CV
       ELSE
        HP(J,K)=HPO(J,K)
       END IF


        PU=P(J,K)*Sue(J)+P(J+1,K)*Suw(j)       ! at U - CV
       IF (PU.GT.PBACK) THEN
        HU(J,K)=HUO(J,K)+AC*(PU-PBACK)
       ELSE
        HU(J,K)=HUO(J,K)
       END IF

        PV=P(J,K)*SVSKM1+P(J,KM1)*SVNKM1       ! at V - CV
       IF (PV.GT.PBACK) THEN
        HV(J,KV)=HVO(J,KV)+AC*(PV-PBACK)       ! DMAX1(PV,PCAV)   ! ..
       ELSE
        HV(J,KV)=HVO(J,KV)
       END IF

        PV=P(J,K)*SVSK+P(J,KP1)*SVNK           ! at V - CV
       IF (PV.GT.PBACK) THEN
        HV(J,KP1)=HVO(J,KP1)+AC*(PV-PBACK)     !    DMAX1(PV,PCAV)   ! ..
       ELSE
        HV(J,KP1)=HVO(J,KP1)
       END IF

        JM1=J
c    !....................!
      END DO
c    !....................! J=JPMIN, JPMAX

         IF ((IFULL.EQ.1).AND.(INTER.EQ.0)) THEN
           HP(JPSTOP,K)=HP(JPMIN,K)            ! periodicity
           HU(JPSTOP,K)=HU(JPMIN,K)            ! if foil does not lift
           HV(JPSTOP,KV)=HV(JPMIN,KV)          ! i.e. a continuous elastic shell
         ELSE IF (IFULL.EQ.0) THEN             ! .........
           HP(JPSTOP,K)=HU(JPMAX,K)            ! upwind deformation
           HU(JPSTOP,K)=HP(JPSTOP,K)           ! at trailing edge
           HV(JPSTOP,KV)=HP(J,K)*SVSKM1+HP(J,KM1)*SVNKM1
           HV(JPSTOP,KP1)=HP(J,K)*SVSK+HP(J,KP1)*SVNK      ! 9/30/93 !###
         END IF                                !
c    !......................!                  !....................!
      END IF                                   ! compliant surface bearing
c    !......................!                  !....................!

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Edgeu                                                       **
C **                                                                         **
C **  EDGEU:  Updates inlet U velocity and recess pressure due to inertia    **
C **          effects.                                                       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE EDGEU(Dir,IREC,K)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD

      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /LRBOUND/ LEFTBC, RIGHTBC, /SWITCH/ IPROP
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
c     ................................................................
      DOUBLE PRECISION DXP, DXU,SUW, SUE, HP, HU,HV,
     +                 U, V, P, RHP, EMP
      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC, ASPEC,HRECD,
     +                 PRECL, PRECR, CK,BETAK,THC,T,TRECL, TRECR

      INTEGER JUMIN, JUMAX, JUSTART, JUSTOP,
     +        JPMIN, JPMAX, JPSTART, JPSTOP, LEFTBC, RIGHTBC, IPROP,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION Zero, UW, UE,FWU,FEU,RHOE,EMUE,RHOW,EMUW, Dir,
     +                 RHOr, EMUr, RHOem, EMUem, Pep, Pem, Pemold,
     +                 HEDGE,ETA,AX,RENR,RENH,ALFAR,MACH2, L4,
     +                 CKR,BTR,THR
      INTEGER IREC, K, JP1, Kr, IPROPo

C ----------------------------------------------------------------------------
C --  EDGEU code: Boundary condition=2 ==>>> recess edge
C ----------------------------------------------------------------------------
C On interrecess lands, Inter=1 & Inerp=1         ! Update Ue & Uw velocities
C ................................................!...........................
      ZERO=0.0D0
      P(JPSTART,K)=PRECL                          ! Jumin=Jpstart
      P(JPSTOP, K)=PRECR                          ! Jumax=Jpstop

C    !..............................!             !
      IF ((LEFTBC.EQ.2).AND.(U(JUMIN,K).GT.ZERO))  THEN
C    !..............................!             !ON west side FLOW enters lands
      Pem=PRECL                     ! Pe-         !...........................
      RHOE=Rhp(JUMIN,K)*Sue(JUMIN)+Rhp(JUMIN+1,K)*Suw(JUMIN)
      FWU=RHOE*HU(JUMIN,K)*U(JUMIN,K)             ! *DYPK = FEU

      RHOW=RHP(JUMIN,K)
      EMUW=EMP(JUMIN,K)

      UW=FWU/RHOW/HP(JUMIN,K)

      HEDGE=HRECD+HP(JUMIN,K)
      ETA=HP(JUMIN,K)/HEDGE

      CALL LOCPROPS(RHOr,EMUr,CKR,BTR,THR,PRECL,TRECL)

      AX=KLOSXd*RHOw*(1.0D0-(ETA*RHOw/RHOr)**2)      ! RHOe+=RHOw

C     Calc. Pedge at entrance                    ! Pe+:  Pedge at entrance
      Pep =DMAX1(PRECL-AX*UW*UW,PCAV)            !       of film lands

 45   P(JPSTART,K)=Pep                           ! Pe+
                                                 ! Update properties at edges
      CALL LOCPROPS(Rhp(JPSTART,K),Emp(JPSTART,K),CK(JPSTART,K),
     +        BETAK(JPSTART,K),THC(JPSTART,K),Pep,T(JPSTART,K)) ! Te+ = Tedge from CALCTM

C    !..............................!
      END IF
C    !..............................! LEFTBC=rec edge & U>0



C    !..............................!             !...........................
      IF ((RIGHTBC.EQ.2).AND.(U(JUMAX,K).LT.ZERO)) THEN
C    !..............................!             ! ON EAST SIDE FLOW enters lands
      Pem=PRECR                     ! Pe-         !...........................
      CALL LOCPROPS(RHOem,EMUem,    ! RHOe-               ! Recess pressure
     +         CKR,BTR,THR,Pem,TRECR)                      ! and properties
      RHOW=Rhp(JUMAX,K)*Sue(JUMAX)+Rhp(JPSTOP,K)*Suw(JUMAX)
      FEU=RHOW*HU(JUMAX,K)*U(JUMAX,K)             ! *DYPK = FWU
      RHOE=Rhp(JPSTOP,K)                          !
      UE=FEU/RHOE/HP(JPSTOP,K)                    !

      HEDGE=HRECD+HP(JPSTOP,K)                    !
      ETA=HP(JPSTOP,K)/HEDGE                      !
      AX=KLOSXu*RHOE*(1.0D0-(ETA*RHOE/RHOem)**2)  !
C                                                 ! Calc. Pedge at entrance
      Pep=DMAX1(Pem-AX*UE*UE,PCAV)                !       of film lands
      P(JPSTOP,K)=Pep                             ! Update properties at edges

      CALL LOCPROPS(Rhp(JPSTOP,K),Emp(JPSTOP,K),CK(JPSTOP,K),
     +BETAK(JPSTOP,K),THC(JPSTOP,K),Pep,T(JPSTOP,K))

C    !..............................!
      END IF
C    !..............................! RIGHTBC=rec edge & U<0


      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Edgev                                                       **
C **                                                                         **
C **  EDGEV:  Updates V velocity and edge pressures due to inertia effects.  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE EDGEV(Dir, DYPK, SVNK, SVSK,KLOSY, INERP, Kstep)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HOFILM/ HPO(MAXNXT,-MAXNYI:MAXNYI),
     +                HUO(MAXNXT,-MAXNYI:MAXNYI),
     +                HVO(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /URECJET/UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP
      COMMON /COMPLIA/ AC, ETAC, PBACK, LIFT

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /BTYPE/ BEARING
c     .................................................................
      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 U, V, P, RHP, EMP, HPO, HUO, HVO,
     +                 PRECdo, PRECup, PRjet, UREC,
     +                 PMAX,PMAXP,PEPS, AC, ETAC, PBACK,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD,
     +                 CK,BETAK,THC,T

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,BEARING,
     +        COUNTP, MAXCOUNTP, LIFT
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION zero,DYPK, SVSK, SVNK, KLOSY, Dir,FNV,FSV,
     +                 DYP,FEW,VS,PEDGE,RHOR,EMUR,CKR,BTR,THR,
     +                 AY,ETA,RHOS,EMUS,RHON,EMUN,Pcorr,Uedge
      INTEGER INERP, I, J, IJ, JL, JR, Icase,jmax, KK, KKP1, Kstep,II
C ----------------------------------------------------------------------------
C --  EDGEV code                                                            --
C ----------------------------------------------------------------------------
c 	Dir: flow director: = 1.0 Right side ob HJB, 0< y -> <Lright/R
c                    	    =-1.0 Left side of HJB,  0< y -> -Lleft/R
c............................................................................
      ZERO=0.0D0

C    !-------------------------!               !----------------!
      IF (BEARING.EQ.1) THEN                   ! FOR HJBS
C    !-------------------------!               !----------------!

      KK=Kstep*NPA
      KKP1=KK+Kstep
      IF (INERP.EQ.0)  goto 400 ! EXIT: NO correction

      DYP=DYPK/2.0D0                           ! Size of 1/2 Pcv
C    !.........................!               !..........................
      DO I=1, NPOCKET                          ! Sweep from
C    !.........................!               ! left Vcv to right vcv
          IJ=(I-1)*NXI+NLC                     ! on i-th recess
          JMAX=I*NXI+1                         !
c                                              !....................
          ICASE=1                              ! AT interior recess
          J=IJ   !left node                    ! P nodes
          II=1
C        !.....................!
 101      J=J+1                                !evaluate Vs & Pedge
          IF (J.EQ.JMAX) GOTO 102              ! . .  . . . . . . .
          II=II+1                              ! at internal nodes
          Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KK))
          FEW=DYP*Uedge*(HU(J,KK)-HU(J-1,KK))  ! internal nodes
          GOTO 105                             ! ICASE=1

 201          IF(PEDGE.LE.PCAV) PEDGE=P(J,KK)
              P(J,KK)=PEDGE
              GOTO 101
C        !....................!

C                                              ! On right corner
 102     II=NPC                                !....................
         Uedge=UREC(I,II)*(1.0D0+HRECD/HP(J,KK))
         FEW=DYPK*Uedge*(HP(J,KK)-HU(J-1,KK))
         ICASE=2
         GOTO 105

 202     IF(PEDGE.eq.PCAV) PEDGE=P(J-1,KK)
         P(J,KK)=PEDGE

 103      J=IJ                                 ! On left corner:
          II=1                                 !...................
          Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KK))
          FEW=DYPK*Uedge*(HU(J,KK)-HP(J,KK))
          ICASE=3
          GOTO 105

 203      IF(PEDGE.EQ.PCAV) PEDGE=P(J+1,KK)
          P(J,KK)=PEDGE
          GOTO 204

C       !.......................................!...................
c       ! Correct Pedge
c       !.......................................!..................

 105          RHOS=Rhp(J,KK)
              RHON=RHOS*SVNK+Rhp(j,KKP1)*SVSK

              VS=(V(J,KKP1)*HV(J,KKP1)*(RHON/RHOS)+
     +            Dir*FEW/DXP(J))/HP(J,KK)      ! by continuity

              PEDGE=PRjet(I,II)                 ! = Pe-
              RHOR=RHP(JMAX-1,KK-Kstep)         ! Rhoe-
              IF (Dir*VS.GT.ZERO) THEN
                 ETA=HP(J,KK)/(HRECD+HP(J,KK))
                 AY=KLOSY*RHOS*(1.0D0-(ETA*RHOS/RHOR)**2)
                 PEDGE=PEDGE-AY*VS*VS
                 PEDGE=DMAX1(PEDGE,PCAV)
                 HP(J,KK)=HPO(J,KK)+AC*(PEDGE-PBACK)   !=> IFF COMPLIANCE
              END IF

              GOTO (201, 202, 203) , Icase
C       !........................................!...................

 204     CONTINUE

C    !........................!
      END DO                                   ! ----------------------------
C    !........................!                !

C                                              !
  300 CONTINUE                                 !

C    !.........................!               !..........................
      DO I=1, NPOCKET                          ! Find props at edges from
C    !.........................!               ! left Vcv to right vcv
          IJ=(I-1)*NXI+NLC                     ! on i-th recess
          JMAX=I*NXI+1                         !
c        !.....................!               !
          DO J=IJ,JMAX                         ! Update properties at edges
            CALL LOCPROPS(Rhp(j,KK),Emp(j,KK),CK(j,KK),
     &      BETAK(j,KK),THC(j,KK),P(j,KK),T(j,KK))
          END DO
c        !.....................!
      END DO
c   !..........................!
c


C    !-------------------------!               !----------------!
      ELSE IF (BEARING.EQ.2) THEN              ! FOR SEALS
C    !-------------------------!               !----------------!

C    !.........................!               !..........................
      DO J=2, NXT                              ! SWEEP from
C    !.........................!               ! left Vcv to right vcv
              RHOS=Rhp(j,1)
              RHON=Rhp(j,2)*Svsk+Rhos*Svnk
              FNV=RHON*Dxp(J)*Hv(j,2)*V(j,2)
              FEW=RHOS*Dypk*SWIRL*(Hu(j,1)-Hu(j-1,1))

              FSV=FNV+FEW
              VS=FSV/RHOS/Dxp(j)/Hp(j,1)
              PEDGE=P(j,0)                       != Pexternal
              Pcorr=0.0D0

              IF ((VS.GT.(0.0D0)).AND.(INERP.EQ.1)) THEN
c            !........................................!
                 PEDGE=PEDGE-KLOSY*RHOS*VS*VS
                 PEDGE=DMAX1(PEDGE,PCAV)
                 IF (PEDGE.GT.PMAX) PMAX=PEDGE
                 Pcorr=DABS(PEDGE-P(j,1))
                 IF (Pcorr.GT.PMAXP) PMAXP=Pcorr
              END IF
c           !..........................................!

         IF(Pcorr.LE.PEPS) COUNTP=COUNTP+1

         V(J,1)=VS
         P(J,1)=PEDGE

         CALL LOCPROPS(rhp(j,1),emp(j,1),CK(J,1),
     +   BETAK(J,1),THC(J,1),PEDGE,T(J,1))     ! Update properties at edges

         KK=1
         KKP1=2

C    !........................!
      END DO                                   ! ----------------------------
C    !........................!                !


C    !-------------------------!               !----------------!
      ELSE IF (BEARING.EQ.3) THEN              ! FOR JBS
C    !-------------------------!               !----------------!
      RETURN
C    !-------------------------!               !----------------!
      END IF
C    !-------------------------!               !----------------!

      IF (IFULL.EQ.1) THEN                      ! ............................
      Rhp(1,KK)=Rhp(NXT,KK)                     ! 1 PAD OF 360 degrees
      Emp(1,KK)=Emp(NXT,KK)                     !
      CK(1,KK)=CK(NXT,KK)                       !
      BETAK(1,KK)=BETAK(NXT,KK)                 !
      THC(1,KK)=THC(NXT,KK)                     !
      T(1,KK)=T(NXT,KK)                         !
      V(1,KK )=V(NXT,KK)                        ! 1st Vcv = Nxt(last) Vcv
      P(1,KK) =P(NXT,KK)                        ! set periodicity conditions
      V(1,KKP1)=V(NXT,KKP1)                     ! ----------------------------
      GOTO 500
      END IF

  400 IF (IFULL.EQ.1) THEN
        V(1,KK )=V(NXT,KK)                        ! 1st Vcv = Nxt(last) Vcv
        P(1,KK) =P(NXT,KK)                        ! set periodicity conditions
        T(1,KK) =T(NXT,KK)                        !
        V(1,KKP1)=V(NXT,KKP1)                     ! ----------------------------
      END IF

c    !........................! COMPLIANT effects at bearing inlet
  500 IF (AC.GT.0.0D0) THEN
c    !........................!
      DO J=1, NXT-1
        HP(j,KK)=HPO(j,KK)+AC*(P(j,KK)-PBACK)
        HV(j,KK)=HP(j,KK)
        HU(j,KK)=HUO(j,KK)+AC*( (P(j,KK)*Sue(j)+P(j+1,KK)*Suw(j))
     +                          -PBACK                           )
      END DO
       IF (IFULL.EQ.1) THEN
        HP(NXT,KK)=HP(1,KK)
        HV(NXT,KK)=HV(1,KK)
        HU(NXT,KK)=HU(1,KK)
       END IF
c    !........................!
      END IF
c    !........................! AC>0, COMPLIANCE bearing surface.

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Force                                                       **
C **                                                                         **
C **  FORCE:  Integrates pressure to find forces & moments                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FORCE(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/ P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)

      COMMON /PARAM1/ CLEAR,DIAM,LENGTH,LD,AR,HREC
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /FORCE0/ FFACTOR, FX, FY, TO, TOR
      COMMON /MOMENT0/ MFACTOR, MX,MY
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
C     ................................................................
      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HV, HU,
     +                 XP, XU, YP, YV, DYP, DYV, SVN, SVS,
     +                 COSXP, SINXP, COSXU, SINXU, P,
     +                 PREC,TREC, QREC, QIN, QOUT, QFACTOR,
     +                 PRECdo,PRECup,PRjet
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD, AR, HREC,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 FFACTOR, FX, FY, TO, TOR, MFACTOR, MX, MY,
     +                 EXO, EYO, AXO,AYO, ZO, TEMPK, VSOUND,
     +                 AC, ETA, PBACK
      INTEGER NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,ICHECK,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP, BEARING, LIFT, IF, IL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYVK, PP, PCL, PSL, AREA, PCR, PSR, PRECI,
     +                 LOAD, ANGLE, UR, YO,YY,ZERO
      DOUBLE PRECISION FXBACK,FYBACK

      INTEGER DEVICE, Kstart, Kend, Kstep, K0, Icase, ISYMold,
     +        I, J, JJ,  K, KM1, JMI, JMA, JP1, IJ
C ----------------------------------------------------------------------------
C --  FORCE code                                                            --
C ----------------------------------------------------------------------------
      !WRITE (6,*) 'BEGIN FORCE:', FX
      ZERO=0.0D0
      ISYMold=ISYM
      Icase=1
      YO=ZO*2.0D0/DIAM                  ! ZO/Radius
      FX=ZERO
      FY=ZERO
      MX=ZERO
      MY=ZERO

      IF ((ISYM.EQ.1).AND.(ZO.NE.ZERO)) THEN
          CALL TEMPASYM                        ! => ISYM=0 temporal
      END IF

      Kstep=1
      K0=ISYM
C    !...................................!
      IF (NPOCKET.GT.0) THEN             ! HJB
C    !...................................!
       CALL PJET  !SET Pressures on recesses !###
       GOTO 778
C    !...................................!
      ELSE                               ! SEAL or PAD JB
C    !...................................!
       Kstart=2
       GOTO 222
C    !...................................!
      END IF
C    !...................................!

c
c    !..........................!
 777  IF (NPOCKET.EQ.0) GOTO 222
c    !..........................!
 778  Kstart=Kstep
      Kend=Kstep*NPA

c    !...............                   !................
      DO I=1, NPOCKET-IFULL+1           ! Sweep on INTER-RECESS lands
C    !...............                   !...............
          KM1=K0                        ! NXI=NLC+NPC-2
          JMI=(I-1)*NXI+2               ! limits for pressures on lands
          JMA=JMI+NLC-2                 ! between recesses

c        !..........................!   !
          DO K=Kstart,Kend,Kstep        ! from bottom to top of recess
c        !..........................!   !
              DYVK=DYV(K)               !  Y size of VCV
              J=JMI-1                   !  first pressure index
              PP=DMAX1(PCAV, P(J, K))+DMAX1(PCAV, P(J, KM1))
              PCL=PP*COSXP(J)           ! left P*COS(Angle)
              PSL=PP*SINXP(J)           !  ''  P*SIN(Angle)
              YY=YV(K)-YO               !
c                                       !
              DO J=JMI, JMA             !
c            !...............!          !
                  AREA=DYVK*DXU(J-1)/4.0D0 ! Area of integration/4
                  PP=DMAX1(PCAV, P(J, K))+DMAX1(PCAV, P(J, KM1))
                  PCR=PP*COSXP(J)       ! right P*COS(Angle)
                  PSR=PP*SINXP(J)       !  ''   P*SIN(Angle)
                  FX=FX+(PCR+PCL)*AREA  ! X-force
                  FY=FY+(PSR+PSL)*AREA  ! Y-force
                  MY=MY+(PCR+PCL)*AREA*YY  ! Y-Moment
                  MX=MX-(PSR+PSL)*AREA*YY  ! X-Moment
                  PCL=PCR               !
                  PSL=PSR               !
              END DO                    !
c            !...............!          ! J=JMI,JMA
c                                       !
              KM1=K                     !
c                                       !
            IF (I.GT.NPOCKET) GOTO 111  !

              IJ=0
c                                       !JMA+NPC-2=I*NXI
c           !.....................!     ! .......................
              DO J=JMA, JMA+NPC-2       ! At interior of recess I
c           !.....................!     ! .......................
                JP1=J+1                 !
                AREA=DXU(J)*DYVK/2.0D0  ! Area of integration/2
                IJ=IJ+1
                FX=FX+AREA*(PRjet(I,IJ  )*COSXP(J)+
     +                      PRjet(I,IJ+1)*COSXP(JP1))
                FY=FY+AREA*(PRjet(I,IJ  )*SINXP(J)+
     +                      PRjet(I,IJ+1)*SINXP(JP1))
                MY=MY+AREA*(PRjet(I,IJ  )*COSXP(J)+
     +                      PRjet(I,IJ+1)*COSXP(JP1))*YY
                MX=MX-AREA*(PRjet(I,IJ  )*SINXP(J)+
     +                      PRjet(I,IJ+1)*SINXP(JP1))*YY
c           !.....................!     ! .......................
              END DO                    ! end recess contribution
c           !.....................!     ! .......................

c        !..........................!   !
 111      END DO                        ! k=1, Npa
c        !..........................!   !

c    !............                      !..............
      END DO                            ! i=1, Npocket
C    !.......                           !..............

      Kstart=Kstep*NPAP1
C----!

 222  Kend=Kstep*NYI
C----!
C    !......................!           !------------------

      DO K=Kstart,Kend,Kstep            ! On EXTENDED LANDS
C    !......................!           !------------------
          KM1=K-Kstep                   !
          DYVK=DYV(K)                   ! Y size VCV
          YY=YV(K)-YO                   ! moment arm
          J=1                           !
          PP=DMAX1(PCAV, P(J, K))+DMAX1(PCAV, P(J, KM1))
          PCL=PP*COSXP(J)               ! left P*COS(Angle)
          PSL=PP*SINXP(J)               !  ''  P*SIN(Angle)

c        !.....................!        !
          DO J=2, NXT                   !
c        !.....................!        !
              AREA=DYVK*DXU(J-1)/4.0D0  ! Area of integration/4
              PP=DMAX1(PCAV, P(J, K))+DMAX1(PCAV, P(J, KM1))
              PCR=PP*COSXP(J)           ! right P*COS(Angle)
              PSR=PP*SINXP(J)           !  ''   P*SIN(Angle)
              FX=FX+(PCR+PCL)*AREA      ! X-force
              FY=FY+(PSR+PSL)*AREA      ! Y-force
              MY=MY+(PCR+PCL)*AREA*YY   ! Y-Moment
              MX=MX-(PSR+PSL)*AREA*YY   ! X-Moment
              PCL=PCR                   !
              PSL=PSR                   !
          END DO                        ! j=2, Nxt
c        !.....................!        !
                                      !
      END DO                            ! k= Npap1,  Nyi
C    !....................              ! k=-Npap1, -Nyi
       !write(6,*)FX
       !write(6,*)FX

C --------------------------------------!----------------------------
C --------------------------------------!----------------------------

 888  IF (ISYM.EQ.0) THEN               ! for asymmetric HJB

c    !..................................!
      IF (BEARING.EQ.2) GOTO 999

      IF (Icase.eq.2) THEN

        IF (MODEL.eq.1) THEN              !
           FX=FX*2.0D0                     ! DOUBLE ROW HJB
           FY=FY*2.0D0                     !
           MX=0.0D0
           MY=0.0D0
        END IF

        GOTO 999

      END IF      ! ICASE=2

         Icase=2
         Kstep=-1
         Kstart=-2
         goto 777
c    !..................................!
      ELSE                              ! for symmetric HJB
c    !..................................!
         !WRITE (6,*) 'MIDDLE FORCE:', FX
         FX=FX*2.0D0
         FY=FY*2.0D0
         MX=0.0D0
         MY=0.0D0

c    !..................................!
      END IF
c    !..................................!

c.......................................!...................

 999  CONTINUE

c    !...................................!
c     Forces due to back pressure in foil or air bearing
c    !...................................!
      Kstart=-NYI
      IF (BEARING.EQ.2) Kstart=1
      FXback=-PBACK*(SINXP(NXT)-SINXP(1))*(YP(NYI)-YP(Kstart))
      FYback=+PBACK*(COSXP(NXT)-COSXP(1))*(YP(NYI)-YP(Kstart))
c    !...................................!
      IF ((AC.GT.ZERO).OR.(IF.EQ.5)) THEN
c    !...................................!
      FX=FX+FXBACK
      FY=FY+FYBACK
c    !...................................!
      END IF
c    !...................................!

      FX=FFACTOR*FX                     ! forces Fx&Fy in N
      FY=FFACTOR*FY                     !
      MX=MFACTOR*MX                     ! Moments Mx&My in N-m
      MY=MFACTOR*MY                     !
      ICHECK=0                          !
      ISYM=ISYMold                      ! RE-store original value

C .....................................................................
C     find MIN and MAX Pressures        ! on optionst.f
      CALL MINMAXP
C .....................................................................
      !WRITE (6,*) 'END FORCE:', FX

      END

C
C *****************************************************************************
C **                                                                         **
C **  Subroutine Torque                                                      **
C **                                                                         **
C **  TORQUE:  Calculates friction torque in film lands.                     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TORQUE(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /URECJET/ UREC(MAXNPOCK,MAXNPOCK+2)
      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC,HRECD
      COMMON /FORCE0/ FFACTOR, FX, FY, TO, TOR

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c     ................................................................
      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 DYP, DYV, SVN, SVS, U, V, P, RHP, EMP,
     +                 PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 PRECdo, PRECup,PRJET,UREC

      DOUBLE PRECISION REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 FFACTOR, FX, FY, TO, TOR
      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO, RUGRs,RUGSs

      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,IFULL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP, BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION VS, H, VN, VP, GAMA, PE, PP, DX,
     +                 CKR,BTR,THR,RHOR,EMUR,FNGAMA,
     +                 TORQ, RHOP,EMUP, Uedge
      INTEGER DEVICE, Icase, Kstart, Kend, Kendp1, Kstep, K0,
     +        I, J, K, KP1, JSTA, JSTO, JP1, IJ
C ----------------------------------------------------------------------------
C --  TORQUE code                                                           --
C ----------------------------------------------------------------------------
C JOURNAL SHEAR STRESS is equal to:
C TXY(h)=(H/2) dP/dx + GAMA
C where GAMA=EMUxSPEED/H for laminar flows,
c       GAMA=(EMU/H)x(U ks - (U-SPEED) kr)/4.0 for fully turbulent flows.
c
C function GAMA returns moving boundary shear stress as def. by HIRS
c modified recess torque by LSA on 10/12/95
C.............................................................................

      Icase=1
      Kstep=1
      TOR=0.0D0

c   !....................!
 777  Kstart=Kstep
      IF (NPOCKET.EQ.0) GOTO 222        !==> SEAL or JBearing

      Kend=Kstep*NPA
      Kendp1=Kend+kstep
c   !....................!              !:
      DO I=1, NPOCKET-IFULL+1           ! --> on lands between recesses
C   !....................!              !
          JSTA=(I-1)*NXI+1              !
          JSTO=JSTA+NLC-2               !=(I-1)*NXI+NLC-1

c        !.................!            !from
          DO J=JSTA, JSTO               !left to right on land
c        !.................!            !between recesses
              JP1=J+1                   !
              VS=V(J,Kstart)*SUE(J)+V(JP1,Kstart)*SUW(J)
              DX=DXU(J)
c                                       !
              DO K=Kstart,Kend,Kstep    ! from bottom to edge of recess
                  KP1=K+Kstep           !
                  VN=V(J, KP1)*SUE(J)+V(JP1, KP1)*SUW(J)
                  VP=(VN+VS)/2.0D0      ! V vel. at center of UCV

                  Rhop=Rhp(J,K)*Sue(J)+Rhp(JP1,K)*Suw(J)
                  EMup=Emp(J,K)*Sue(J)+Emp(JP1,K)*Suw(J)
                  H=HU(J,K)
                  GAMA=FNGAMA(U(J,K),VP,H,SPEED,REP,RHOP,EMUP)
                  PE=DMAX1(PCAV, P(JP1, K))
                  PP=DMAX1(PCAV, P(J, K)) !
                  TOR=TOR+(H*(PE-PP)/2.0D0/DX+GAMA)*DX*DYP(K)
                  VS=VN                 !
              END DO                    ! K=-1,-Npa,-1 LEFT side
C                                       ! K= 1, Npa, 1 RIGHT Side

C         !..............!              !
           END DO                       ! J=left to right on land
C         !..............!              !

         IF (I.GT.NPOCKET) GOTO 111     !-> goto extended lands

c

            CALL LOCPROPS(RHOR,EMUR,CKR,! fluid props. at recesses
     +      BTR,THR,PREC(I),TREC(I))    ! pressure & temperature

         JSTA=JSTO+1                    ! (I-1)*NXI+NLC
         JSTO=I*NXI                     !
         IJ=0

C     !..................!              ! K=Kend=+/- NPA
          DO J=JSTA, JSTO               ! On circumferential edge of
c     !..................!              ! recess
              JP1=J+1                   !  KP1=NPAP1
              H=HU(J,Kend)              !
              IJ=IJ+1                   !
              VP=V(J,Kendp1)*SUE(J)+V(JP1,Kendp1)*SUW(J)
              Rhop=Rhp(J,Kend)*Sue(J)+Rhp(JP1,Kend)*Suw(J)
              EMup=Emp(J,Kend)*Sue(J)+Emp(JP1,Kend)*Suw(J)
              Uedge=UREC(I,IJ)*(1.0D0+HRECD/H)
              GAMA=FNGAMA(Uedge,VP,H,SPEED,REP,RHOP,EMUP)
              TOR=TOR+(DXU(J)*DYP(Kend)/2.0D0)*( GAMA +
     +              H*(P(JP1,Kend)-P(J,Kend))/2.0D0/DXU(J) )

c           within recess assume no axial velocity / WITH smooth surface
c           and YES circumferential pressure gradient
            H=HU(J,Kend)+HRECD          ! for Half V-cv in the recess
            GAMA=FNGAMA(UREC(I,IJ),0.0D0,H,SPEED,REP,RHOR,EMUR)

            TOR=TOR+(DXU(J)*DYP(Kend)/2.0D0)*(GAMA +
     +                 H*(P(JP1,Kend)-P(J,Kend))/2.0D0/DXU(J))
C          !............................! Torque on recess areas
            DO K=Kstart,Kend-Kstep,Kstep! LET V=0 dP/dx=0
C          !............................!
               H=HU(J,K)+HRECD
               GAMA=FNGAMA(UREC(I,IJ),0.0D0,H,SPEED,REP,RHOR,EMUR)
               TOR=TOR+DXU(J)*DYP(K)*(GAMA +
     +                 H*(PRjet(I,IJ+1)-PRjet(I,IJ))/2.0D0/DXU(J)  )
C          !............................!
            END DO                      ! Torque over recess
C          !............................! includes pressure gradient

c     !..................!              !
          END DO                        !J=left to right on recess.
C     !...................!             !

c###vv remove when done
c       IF (I.EQ.1) THEN
c        WRITE (6,170)
c        WRITE (6,171)(PRJET(I,IJ),IJ=1,NPC)
c        WRITE (6,172)(UREC(I,IJ), IJ=1,NPC)
c       END IF
 170  FORMAT (1X,'## calcsolnt/SUB Torque:',20("."))
 171  FORMAT (1X,'PR:',7(G10.5,1X))
 172  FORMAT (1X,'UR:',7(G10.4,1X))
c###^^^^

 111   CONTINUE
C     !--------------------!            ! I=1,..NPOCKET+1-IFULL
       END DO                           ! Next recess+land region
C     !--------------------!            !

        Kstart=Kendp1        !  (NPA+1) or -(NPA+1)
c----!
 222    Kend=Kstep*NYI       ! 1:Right Side, -1:Left side
c----!
c    !.....................!
      DO J=1, NXT-1                     ! ---> Torque on EXTENDED lands
C    !....................!             !
          JP1=J+1                       ! -----------------------------
          VS=V(J,Kstart)*SUE(J)+V(JP1,Kstart)*SUW(J)
          DX=DXU(J)                     !
                                        !
c        !............................! !
          DO K=Kstart,Kend-Kstep,Kstep  ! at interior UCV of extented land
c        !............................! !
              KP1=K+Kstep               !
              VN=V(J, KP1)*SUE(J)+V(JP1, KP1)*SUW(J)
              VP=(VN+VS)/2.0D0
              H=HU(J,K)
              Rhop=Rhp(J,K)*Sue(J)+Rhp(JP1,K)*Suw(J)
              EMup=Emp(J,K)*Sue(J)+Emp(JP1,K)*Suw(J)
              GAMA=FNGAMA(U(J,K),VP,H,SPEED,REP,RHOP,EMUP)
              PE=DMAX1(PCAV,P(JP1,K))   !
              PP=DMAX1(PCAV,P(J,K))     !
              TOR=TOR+(H*(PE-PP)/2.0D0/DX+GAMA)*DX*DYP(K)
              VS=VN                     !
c        !............................! !
          END DO                        !
c        !............................! !K: axial sweep

c        At last U-1/2 cv:
          H=HU(J,Kend)
          Rhop=Rhp(J,Kend)*Sue(J)+Rhp(JP1,Kend)*Suw(J)
          Emup=Emp(J,Kend)*Sue(J)+Emp(JP1,Kend)*Suw(J)
c                                       ! At bearing exit
          GAMA=FNGAMA(U(J,Kend),VS,H,SPEED,REP,RHOP,EMUP)

          PE=DMAX1(PCAV,P(JP1,Kend))
          PP=DMAX1(PCAV,P(J,Kend))
          TOR=TOR+(H*(PE-PP)/2.0D0/DX+GAMA)*DX*DYP(Kend)

C     !....................!            !
      END DO                            !J: circumferential sweep
C    !..................................!......................


 888  IF (ISYM.EQ.0) THEN
c    !......................! ISYM=0, asymmetric bearing
        IF (BEARING.EQ.2) GOTO 999
        IF (Icase.EQ.2) GOTO 998
        Kstep=-1
        Icase=2
        GOTO 777

 998    IF (MODEL.EQ.1) THEN              !..............
          TOR=2.0D0*TOR                   ! 2-row HJB
        END IF                            !..............

c    !......................! ISYM=1, symmetric HJB
      ELSE
c    !......................!
        TOR=2.0D0*TOR

c    !......................!
      END IF
c    !......................!


 999    TORQ=TO*TOR                     ! TORQ: Local variable

C.......................................!..............................

  100 FORMAT (' |    Torque on Film Lands=', E12.5E2, ' N-m', 36X, '|',
     +        /, ' +', 77('-'), '+')
  200 FORMAT (' ', '+', 77('-'), '+')

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Tdma                                                        **
C **                                                                         **
C **  TDMA:  Solution of algebraic tridiagonal system of equations.          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TDMA(ZLEFT, ZRIGHT, JMAX)


      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)

      DOUBLE PRECISION A, B, C, D

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION ZLEFT, ZRIGHT, DEN
      INTEGER JMAX, J, JM1

C ----------------------------------------------------------------------------
C --  TDMA code                                                             --
C ----------------------------------------------------------------------------

      D(1)=ZLEFT                        ! = F(1)
      A(1)=0.0D0                        ! = E(1)
      DO J=2, JMAX                      ! 1st sweep:
          JM1=J-1                       ! Find recursive
          DEN=1.0D0/(C(J)*A(JM1)+B(J))  !
          A(J)=-A(J)*DEN                ! = E(J)
          D(J)=(D(J)-C(J)*D(JM1))*DEN   ! = F(J)
      END DO                            !
      B(JMAX+1)=ZRIGHT                  ! Set BC at right
      DO J=JMAX, 1, -1                  ! 2nd sweep:
          B(J)=A(J)*B(J+1)+D(J)         ! calculate soln.
      END DO                            ! stored in vector B

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Endseal                                                     **
C **                                                                         **
C **  ENDSEAL: Updates V velocity and discharge pressure due to end seal  .  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ENDSEAL(Dir,DYPK, Psump, Sealcoef,alfp, Kstep)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
c     .................................................................
      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 U, V, P, RHP, EMP,PMAX,PMAXP,PEPS
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        COUNTP, MAXCOUNTP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION zero,DYPK, Dir, Psump, Sealcoef,Pcorr,
     +                 Alfp, Rhoe,Rhow, Fe, Fw, Vend, Pend,
     +                 CK,BETAK,THC,T
      INTEGER J,KK, Kstep

C ----------------------------------------------------------------------------
C --  ENDSEAL CODE                                                         --
C ----------------------------------------------------------------------------
c Dir: flow director: =1.0 Right side ob HJB, 0<y -> <L/D
c                     =-1.0 Left side of HJB, 0: y -> -L/D
c............................................................................
      IF (IFULL.EQ.0) THEN
       PRINT *, 'calcsolnt> ENDSEAL ERROR !'
       PRINT *, 'End seals make sense only for 360deg bearings'
       PRINT *, 'PROGRAM STOPS - REVISE INPUT DATA'
       CALL BEEPER
       STOP
      END IF

      KK=Kstep*NYI
      ZERO=0.0D0

      J=NXT-1
      Rhow=Sue(J)*Rhp(J,KK)+Suw(J)*Rhp(1,KK)
      Fw=DYPK*Rhow*U(J,KK)*Hu(J,KK)

C    !.........................!               !..........................
      DO J=1, NXT-1                            ! Sweep from
C    !.........................!               ! left Vcv to right vcv
        Rhoe=Sue(j)*Rhp(j,kk)+Suw(j)*Rhp(j+1,kk)
        Fe=DYPK*Rhoe*U(j,KK)*Hu(j,KK)

c        Vend=(V(j,KK)*HV(j,KK)-Dir*(Fe-Fw)/Dxp(j)/Rhp(j,KK))/HP(j,KK)
         Vend=V(j,KK)

              IF (Dir*Vend.GT.zero) THEN
                 Pend=Psump+Sealcoef*Vend*Vend*Rhp(j,KK)
                 Pend=P(j,KK)*(1.0D0-Alfp)+Pend*Alfp
              ELSE
                 Pend=Psump
              END IF

              Pcorr=DABS(Pend-P(j,KK))
              P(j,KK)=Pend
              IF (Pcorr.le.PEPS) Countp=Countp+1

         Fw=Fe

         CALL LOCPROPS(Rhp(j,KK),Emp(j,KK),CK(J,KK),
     +   BETAK(J,KK),THC(J,KK),P(j,KK),T(J,KK))
C    !........................!
      END DO                                   ! ----------------------------
C    !........................!                !

      P(NXT,KK) =P(1,KK)                        ! set periodicity conditions
      Rhp(NXT,KK)=Rhp(1,KK)
      Emp(NXT,KK)=Emp(1,KK)
      CK(NXT,KK)=CK(1,KK)
      BETAK(NXT,KK)=BETAK(1,KK)
      THC(NXT,KK)=THC(1,KK)
      T(NXT,KK)=T(1,KK)

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine TEMPASYM                                                    **
C **                                                                         **
C **  Temporal Asymmetric HJB for case of pivot point not at Yo=0            **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TEMPASYM

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX,EWY,EWEAR,BETAW,IWEAR

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ISYM, ICSTEP
c     .................................................................
      DOUBLE PRECISION U,V,P, RHOP, EMUP, DRHOP, DEMUP,
     +       EXO,EYO,EWX,EWY,EWEAR,BETAW,AXO,AYO,ZO,
     +       CK,BETAK,THC,T,HCB,HCJ,DRHOT,DEMUT
      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ISYM, ICSTEP, IWEAR, IFULL


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER J,K,Kstart,Kend,Imod

C ----------------------------------------------------------------------------
c     .......................................................................
c     Extend U,V,P fields if TEMPORAL ASYMMETRY: ISYM=0, Yo<>0
c     .......................................................................

  71  ISYM=0
c     -------------!!!
c      write (6,72)
  72  format (' ',2X,'Pivot for moment coeffs. Yo<>0.0D0',/,
     +            3X,'SET Bearing Assymetry and Extend Mirror',/,
     +            3X,'Image of UVP fields',/,3X,50('.'))

C##      IF (IWEAR.eq.1) THEN
C##         CALL HWEAR(EXO, EYO, NXT, NYI)
C##      ELSE
         CALL FILMH(EXO, EYO, NXT, NYI)
C##      END IF


c    !................................. ! SWEEP in circumferential dir.
      DO J=1, NXT
c    !................................. ! ............................!
        DO K=1, NYI			! DIMENSIONLESS
           U(J,-K)=U(J,K)               ! circumferential velocity
           V(J,-K)=-V(J,K)		! axial velocity
           P(J,-K)=P(J,K)		! pressure
           T(J,-K)=T(J,K)		! temperature
           RHOP(J,-K)=RHOP(J,K)		! fluid density
           EMUP(J,-K)=EMUP(J,K)		! fluid viscosity
           Drhop(J,-K)=Drhop(J,K)	! dRHO/dP
           Demup(J,-K)=Demup(J,K) 	! dEMU/dP
           Drhot(J,-K)=Drhot(J,K)	! dRHOT/dT
           Demut(J,-K)=Demut(J,K)	! dEMU/dT
           CK(J,-K)=CK(J,K)		! fluid specific heat
           BETAK(J,-K)=BETAK(J,K)	! fluid thermal expansion coeff
           THC(J,-K)=THC(J,K)		! fluid thermal conductivity
           HCB(J,-K)=HCB(J,K)		! heat convection coeffs
           HCJ(J,-K)=HCJ(J,K)		! at bearing and journal surfaces
        END DO                          !..............................!
           U(J,0)=.5D0*(U(J,-1)+U(J,1)) ! AT Y=0, index 0 is average
           V(J,0)=0.0D0                 !         value from indexes -1 & 1
           P(J,0)=.5D0*(P(J,-1)+P(J,1))
           RHOP(J,0)=.5D0*(RHOP(J,-1)+RHOP(J,1))
           EMUP(J,0)=.5D0*(EMUP(J,-1)+EMUP(J,1))
           Drhop(J,0)=.5D0*(Drhop(J,-1)+Drhop(J,1))
           Demup(J,0)=.5D0*(Demup(J,-1)+Demup(J,1))
           T(J,0)=.5D0*(T(J,-1)+T(J,1))
           CK(J,0)=.5D0*(CK(J,-1)+CK(J,1))
           BETAK(J,0)=.5D0*(BETAK(J,-1)+BETAK(J,1))
           THC(J,0)=.5D0*(THC(J,-1)+THC(J,1))
           HCB(J,0)=.5D0*(HCB(J,-1)+HCB(J,1))
           HCJ(J,0)=.5D0*(HCJ(J,-1)+HCJ(J,1))
           Drhot(J,0)=.5D0*(Drhot(J,-1)+Drhot(J,1))
           Demut(J,0)=.5D0*(Demut(J,-1)+Demut(J,1))
c    !................................. ! ............................!
      END DO
c    !................................. ! ............................!


      END

C *****************************************************************************
c
c calcsolnt.f > hydrojet program by Luis San Andres  10/12/95
c
c *****************************************************************************


