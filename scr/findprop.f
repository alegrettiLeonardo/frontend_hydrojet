C ######     #    #    #  #####   #####   #####    ####   #####           ######
C #          #    ##   #  #    #  #    #  #    #  #    #  #    #          #
C #####      #    # #  #  #    #  #    #  #    #  #    #  #    #          #####
C #          #    #  # #  #    #  #####   #####   #    #  #####    ###    #
C #          #    #   ##  #    #  #       #   #   #    #  #        ###    #
C #          #    #    #  #####   #       #    #   ####   #        ###    #
C
C #    #   #   #  #####   #####    ####   ######  #       ######  #    #
C #    #    # #   #    #  #    #  #    #  #       #       #        #  #
C ######     #    #    #  #    #  #    #  #####   #       #####     ##
C #    #     #    #    #  #####   #    #  #       #       #         ##
C #    #     #    #    #  #   #   #    #  #       #       #        #  #
C #    #     #    #####   #    #   ####   #       ######  ######  #    #

C findprop.f > hydroflex code Drs. Luis San Andres, TexasA&MUniv. 1994
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR I
c Technical monitor: Mr. James Walker, NASA Lewis Research Center


C *****************************************************************************
C **                                                                         **
C **  Subroutine Finprop                                                     **
C **                                                                         **
C ** Finprop:  Calculates Final values of Rho&Emu and derivatives to P       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FINPROP(ICASE, DEVICE,RPM,PS,PA)

      IMPLICIT NONE

      INCLUDE 'params.f'

c----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT,QFACTOR
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DBTA/ DBKP(MAXNXT,-MAXNYI:MAXNYI),
     +              DBKT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DCPS/ DCPP(MAXNXT,-MAXNYI:MAXNYI),
     +              DCPT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DTHC/  DKP(MAXNXT,-MAXNYI:MAXNYI),
     +               DKT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ TK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TMAXMIN/ TKMAX,TKMIN
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /THERMID/ ISOTH

      COMMON /FACTORS/ REP,REYP,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /PARAM3/ EMU, RHO, PC, CD, DORIF, LOSXSI,ALPHA
      COMMON /PROPTYP/ RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2

      COMMON /RECPAR/Hrecu,Vsup,Beta
      COMMON /LIQUID/Tempk,Vsound,IF,IL

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /SWITCH/IPROP
c     ....................................................................
      DOUBLE PRECISION Prec,Qrec,Qin,Qout, Qfactor, Hp, Hu, Hv, U, V,
     +                 P,Rhop,EMup,Drhot,Demut,Drhop,Demup ,Trec,
     +                 DBKP, DBKT ,DCPP, DCPT ,DKP,  DKT,
     +                 CK,BETAK,THC,TK, ALFT, UC, TC, Ec, TKMAX,TKMIN


      DOUBLE PRECISION REP,REYP,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 Clear,Diam,Length,LD,Ar,Hrec,
     +                 RHOTYP,EMUTYP,DENA,VISA,PSA,PATYP,
     +                 DEN12P12, VIS12P12, P2,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 Hrecu,Vsup,Beta, Tempk, Vsound

      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, BEARING,ISOTH,
     +        IF, IL, DEVICE, SVERB, DVERB, BEEP, ISYM, ICSTEP, IPROP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DeltaP,DeltaT,Rhopp,Emupp, Vorif,Pro,
     +                 Rhomin,Rhomax,Emumin,EMumax, Mflow,
     +                 Xmin,Xmax,Vel, Up, Vp, zero,PI,
     +                 Vfactor,Velmin,Velmax,Vsmin,Vsmax,
     +                 Rey,Reymax,Reymin ,CP1,BK1,TH1,
     +                 Rhopt,Emupt,CPPT,BTAT,THT,CPPP,BTAP,THP,
     +                 TEA,TMAX,TMIN
      INTEGER Icsym, I,J,K,Kp1, Kstart,Kstep, jm1,ICASE,NYIP1 ,JREC,
     +        IDUM
C ----------------------------------------------------------------------------
C --  FINPROP code                                                          --
C ----------------------------------------------------------------------------
      IPROP=1                                       ! evaluate props
      PI=DACOS(-1.0D0)                              !.................
      Vfactor=Clear*Clear*2.0D0*(PSA)/EMU/DIAM      ! veloc. conversion factor
      zero=0.0D0

      Icsym=1
      Kstart=1
      Kstep=1

      DeltaP=0.0010D0
      DeltaT=0.0010D0
      call Locprops(Rhomin,Emumin,CP1,BK1,TH1,P(1,1),tk(1,1))
      Rhomax=Rhomin
      Emumax=Emumin

      Vsmin=Vsound
      Vsmax=Vsound

      jm1=1+IFULL*(NXT-2)

      Up=(U(1,1)+U(jm1,1))/2.0D0       ! Speed at Pnodes
      Vp=(V(1,1)+V(1,2))/2.0D0         !
      Vel=DSQRT(Up*Up+Vp*Vp)           ! ABS Fluid velocity
      Velmin=Vfactor*Vel
      Velmax=Velmin

      Reymin=Rep*Hp(1,1)*(Rhomin/Emumin)*Vel
      Reymax=Reymin

 66   jm1=1+IFULL*(NXT-2)

c   !......................!
      DO J=1, NXT
c   !......................!                   ! Sweep over all bearing surface
c                                              ! & calculate properties at
c                                              ! all pressure nodes
         DO k=Kstart, Kstep*NYI, Kstep
c       !.............................!

         CALL Locprops(Rhopp,Emupp,CP1,BK1,TH1,P(j,k),tk(j,k))
         Rhop(j,k)=Rhopp
         Emup(j,k)=Emupp
         CK(J,K)=CP1
         BETAK(J,K)=BK1
         THC(J,K)=TH1

         IF(Rhopp.gt.Rhomax) THEN
            Rhomax=Rhopp
         ELSE IF(Rhopp.lt.Rhomin) THEN
            Rhomin=Rhopp
         END IF

         IF(Emupp.gt.Emumax) THEN
            Emumax=Emupp
         ELSE IF(Emupp.lt.Emumin) THEN
            Emumin=Emupp
         END IF

c       !.......................................!
c        Determine MAx & MIn Fluid Speeds
c                  & Flow Reynolds numbers
c       !.......................................!
c        kp1=MIN0(k+1,Nyi)
         Kp1=K+Kstep
         IF ((Kstep*Kp1).gt.Nyi) Kp1=Kstep*Nyi

         Up=(U(j,k)+U(jm1,k))/2.0D0       ! Speeds at Pnodes
         Vp=(V(j,k)+V(j,kp1))/2.0D0       !
         Vel=DSQRT(Up*Up+Vp*Vp)           ! ABS Fluid velocity

         Rey=Rep*(Rhopp/Emupp)*Hp(j,k)*Vel

         Vel=Vfactor*Vel
         IF (Vel.gt.Velmax) THEN
             Velmax=Vel
         ELSE
             IF(Vel.lt.Velmin) Velmin=Vel
         END IF

         IF (Rey.gt.Reymax) THEN
               Reymax=Rey
         ELSE
             IF(Rey.lt.Reymin) Reymin=Rey
         END IF

c       !.......................................!
c        Determine max&min sonic speed
c       !.......................................!
         IF (Vsound.eq.zero) goto 78

c         IF (IF.GT.4)  GOTO 77

           IF (Vsound.gt.Vsmax) THEN
             Vsmax=Vsound
           ELSE
             IF(Vsound.lt.Vsmin) Vsmin=Vsound
           END IF

 77      IF ((Vel.ge.Vsound).AND.(ICASE.GT.2)) THEN
            WRITE(6,91) J,K,Vel,Vsound
             IF (DEVICE.eq.1) THEN
                WRITE (1,91) J,K, Vel, Vsound
             END IF
         END IF

c       !...............................................!
c       ! Calculate derivatives of density, viscosity
c         and other material properties w/respect to
c         pressure and temperature (dimensionless)
c       !...............................................!
 78      CONTINUE
         CALL LOCPROPS(Rhopt,Emupt,CPPT,BTAT,THT,
     +              P(j,k),TK(j,k)+DeltaT)              ! Rho(T+DT), Emu(T+DT)
         CALL LOCPROPS(Rhopp,Emupp,CPPP,BTAP,THP,
     +              P(j,k)+DeltaP,tk(j,k))              ! Rho(P+DP), Emu(P+DP)

         DBKP(J,K) =(BTAP-BETAK(J,K))/DELTAP
         DBKT(J,K) =(BTAT-BETAK(J,K))/DELTAT
         DCPP(J,K) =(CPPP-CK(J,K))/DELTAP
         DCPT(J,K) =(CPPT-CK(J,K))/DELTAT
         DKT(J,K)  =(THT -THC(J,K))/DELTAT
         DKP(J,K)  =(THP -THC(J,K))/DELTAP
         Drhop(j,k)=(Rhopp-Rhop(j,k))/Deltap
         Demup(j,k)=(Emupp-Emup(j,k))/Deltap
	     DRHOT(j,k)=(Rhopt-Rhop(j,k))/DeltaT
 	     DEMUT(j,k)=(Emupt-Emup(j,k))/DeltaT
c       !....................!
         END DO
c       !....................!     k=1,2.....,NYI*Kstep

       jm1=j

c   !........................!
      END DO
c   !........................!  SWEEP J=1, NXT



c...............................................................
c for asymmetric HJB, calculate on Left side of Bearing
c and update center values.
c...............................................................

c    !-------------------------------!
 79   IF (ISYM.eq.0) THEN
c    !-------------------------------!
         IF (BEARING.EQ.2) GOTO 81   ! => FOR SEAL

         IF (Icsym.eq.2) GOTO 80

            IF (MODEL.eq.1) THEN     ! ==> two row HJB
               NYIP1=NYI+1           !     at symmetry plane
               DO J=1, NXT           !     Z=Lengthr
                  Drhop(J,NYIP1)=zero
                  Demup(J,NYIP1)=zero
                  Drhot(J,NYIP1)=zero           ! ??C#
                  Demut(J,NYIP1)=zero           ! ??C#
                  Rhop(J,NYIP1)=zero
                  Emup(J,NYIP1)=zero
                  CK(J,NYIP1)=zero              ! ??C#
                  BETAK(J,NYIP1)=zero
                  THC(J,NYIP1)=zero
                  Hv(J,NYIP1)=zero
               END DO
            END IF                   ! ==> two row HJB
         Kstart=-1
         Kstep=-1
         Icsym=2
         GOTO 66
      ELSE
         GOTO 81
      END IF

 80   CONTINUE

c    !.................................................!
c     Find props and derivatives at Z=0 for asymmetric
c     bearing
c    !.................................................!

      DO J=1, NXT
         CALL LOCPROPS(Rhop(j,0),Emup(j,0),CK(J,0),
     +        BETAK(J,0),THC(J,0),P(j,0),TK(J,0))
         CALL LOCPROPS(Rhopt,Emupt,CPPT,BTAT,THT,
     +              P(j,0),TK(j,0)+DeltaT)              ! Rho(T+DT), Emu(T+DT)
         CALL LOCPROPS(Rhopp,Emupp,CPPP,BTAP,THP,
     +              P(j,0)+DeltaP,TK(j,0))              ! Rho(P+DP), Emu(P+DP)
         DBKP(J,0) =(BTAP-BETAK(J,0))/DELTAP
         DBKT(J,0) =(BTAT-BETAK(J,0))/DELTAT
         DCPP(J,0) =(CPPP-CK(J,0))/DELTAP
         DCPT(J,0) =(CPPT-CK(J,0))/DELTAT
         DKT(J,0)  =(THT -THC(J,0))/DELTAT
         DKP(J,0)  =(THP -THC(J,0))/DELTAP
         Drhop(j,0)=(Rhopp-Rhop(j,0))/Deltap
         Demup(j,0)=(Emupp-Emup(j,0))/Deltap
         Drhot(j,0)=(Rhopt-Rhop(j,0))/Deltat
         Demut(j,0)=(Emupt-Emup(j,0))/Deltat
      END DO

c.........................................
c AT Recess orifices, FIND Fluid Speed
c.........................................

 81   K=NPOCKET

      IF (K.EQ.0) GOTO 82     !==>  no orifice on bearing pad

      IF (ICASE.EQ.1) K=1     !==> rotationally symmetric HJB

      IF (Vsound.EQ.zero) GOTO 82
c    !---------------------------!
          IF (ICASE.GT.2) THEN
          WRITE(6,97)
          IF (DEVICE.EQ.1) WRITE(1,97)
          END IF

c    !....................!
      DO I=1,K
c    !....................!
        CALL Locprops(Rhopp,Emupp,CP1,BK1,TH1,Prec(I),Trec(I))
        Pro=Prec(I)*(Psa)+Pa
        Rhopp=Rhopp*RHO
        Vorif=DSQRT(2.0D0*DABS(Ps-Pro)/Rhopp)

        IF ((Vorif.GE.Vsound).AND.(ICASE.GT.2)) THEN
           WRITE(6,89) I, Vorif, Vsound
           IF (Device.eq.1) WRITE(1,89) I, Vorif, Vsound
        ELSE IF (ICASE.GT.2) THEN
           WRITE(6,90) I,Vorif, Vsound
           IF (Device.eq.1) WRITE (1,90) I, Vorif, Vsound
        END IF

      END DO
c    !....................! I=1, Npocket

c................................................................
c PRINT Min & Max values of land speed, reynolds numbers
c                 and liquid properties
c................................................................
 82   IF (ICASE.LE.2) GOTO 87

      WRITE(6,911) Reymax, Reymin
      WRITE(6,92) Velmax,Velmin
      IF (DEVICE.eq.1) THEN
        WRITE(1,911) Reymax,Reymin
        WRITE(1,92) Velmax,Velmin
      END IF

c    !................................!
      IF(Vsound.GT.zero) THEN
        WRITE(6,93) Vsmax,Vsmin
        IF (DEVICE.eq.1) WRITE(6,93) Vsmax,Vsmin
      END IF
c    !................................!,

      WRITE(6,94) Rhomax*RHO,Rhomin*RHO,
     +            Emumax*EMU,Emumin*EMU

      IF(DEVICE.eq.1) THEN
        WRITE(1,94) Rhomax*RHO,Rhomin*RHO,
     +              Emumax*EMU,Emumin*EMU
      END IF

c    !..........................! Average Exit Temperature
      IF (ISOTH.EQ.1) THEN
c    !..........................! barotropic model.
        TKMAX=1.0D0
        TKMIN=1.0D0
c    !..........................!
      ELSE
c    !..........................! from THERMAL analysis
        CALL MINMAXT(TMIN,TMAX,TEA)    !=> optionst.f
        TKMAX=DMAX1(TKMAX,TMAX)
        TKMIN=DMIN1(TKMIN,TMIN)
        TMAX=TMAX*TC            ! MAX, MIN on lands
        TMIN=TMIN*TC            ! Average Temp. at exit
        TEA=TEA*TC              ! planes

        WRITE (6, 120) TMAX, TMIN
        WRITE(6,99) TEA
        IF(DEVICE.eq.1) THEN
          WRITE (1, 120) TMAX, TMIN
          WRITE(1,99) TEA
        END IF
c    !..........................!
      END IF
c    !..........................!

 87   CONTINUE

c............................................. FORMAT STATEMENTS

 89   FORMAT(' ',77('.'),/,3X,'WARNING: At Recess(',I3,'):',/,
     +    3X,'Fluid Orifice Velocity is:', E12.5E2,
     +    ' m/s and LARGER THAN',/,
     +    3X,'Sonic Speed:', E12.5E2,/,' ',79('.'))

 90   FORMAT(' ',3X,'Rec:',I3,1X,'Vel. at orifice=', E12.5E2, 2X,
     +       'Sonic Speed=',E12.5E2,' m/s')

 91   FORMAT(' ',77('.'),/,3X,'WARNING: At Pnode (',I3,',',I3,'):',/,
     +    3X,'Fluid Velocity is:', E12.5E2,' m/s and LARGER THAN',/,
     +    3X,'Sonic Speed:', E12.5E2,/,' ',79('.'))
 911  FORMAT(' ',77('.'),/,3X,'FLOW REYNOLDS Number, MAX=', E12.5E2,
     +      '     MIN=',E12.5E2)
 92   FORMAT(' ',77('.'),/,3X,'Fluid Speed on Lands, MAX=', E12.5E2,
     +      'm/s, MIN=',E12.5E2,'m/s')
 93   FORMAT(' ',77('.'),/,3X,'Fluid Land SonicSpeed MAX=', E12.5E2,
     +      'm/s, MIN=',E12.5E2,'m/s')
 94   FORMAT(' ',77('.'),/,3X,'FLUID Density on Lands, MAX=',E12.5E2,
     +       3X,'MIN=',E12.5E2,' Kg/m3',/,3X,
     + '      Viscosity """     MAX=',E12.5E2,3X,'MIN=',E12.5E2,
     + ' Pa-s')
 95   FORMAT(' ',3X,'d(Density)/dPress:   Max=',E12.5E2,
     +                           3X,'Min=',E12.5E2)
 96   FORMAT(' ',3X,'d(Viscosity)/dPress: Max=',E12.5E2,
     +                           3X,'Min=',E12.5E2,/,
     +   1X,77('.'))
 97   FORMAT(' ',34('. '),'VALUES OF:')
 99   FORMAT(' ',79('.'),/,3X,'Average Exit Temperature=', E12.5E2,
     + ' K',/,' ', 79('.'))
 120  FORMAT(' ',79('.'),/,3X,'Max & Min Temperatures, MAX=',E12.5E2,
     +       3X,'MIN=',E12.5E2,' K-deg')

      END
c------------------------------------------------------------------c
c                                                                  c
c hydroflex.f :   revised 2/17/94  by Luis San Andres              c
c                                                                  c
c------------------------------------------------------------------c
