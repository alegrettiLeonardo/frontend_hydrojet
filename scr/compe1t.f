c --------------------------------------------------------------------------
c 10/4/95 modified VCALC for recess edge velocities
c         see UCALC for SWIRL ???
c --------------------------------------------------------------------------
c 9/13/95 updated for angled jet in hydrostatic bearings
c MUST CALL PJET1 on SLAND1
c --------------------------------------------------------------------------
c
c  ####    ####   #    #  #####   ######   ##      #####          ######
c #    #  #    #  ##  ##  #    #  #       # #        #            #
c #       #    #  # ## #  #    #  #####     #        #            #####
c #       #    #  #    #  #####   #         #        #     ###    #
c #    #  #    #  #    #  #       #         #        #     ###    #
c  ####    ####   #    #  #       ######  #####      #     ###    #
C
c
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

c --------------------------------------------------------------------------
c 6/20/95 coefficient for compliance of cell depth on rough seals
c is given as  ix(Res.Uo.HCdim)xRHOj and written also as
c              ix(Res.Uo.HCdim)x(drho/dp.Pj+drho/dT.Tj)
c              here I have neglected the effect of temperature
c              otherwise I will require a GVTimaginary coeffic.
c    RECALL: HCdim=Cell_depth/TYPclearance
c    coefficient GUPI and GVPI have been modified on UCALC & VCALC
c --------------------------------------------------------------------------
c 6/1/94 removed call to SURFACET1 < not needed  if TB&TJ are uniform
c --------------------------------------------------------------------------
c 3/8/94
c   Updated terms with compliance coefficient. UCALC, VCALC

c   Coefficient of dynamic pressure for real fluids with bearing compliance
c   effects is equal to:
c    Gxp = Gxp(real,imag) = Gxh ac/(1+in) + gxrho drho/dp + gxemu demu/dp
c
c    gur RHO1 + gumu EMU1 have been clustered as equal to
c
c    gup P1 + gut T1 = (gur drho/dp + gumu demu/dp) P1 +
c                      (gur drho/dT + gumu demu/dT) T1
c
c			in UCOEFS and VCOEFS subroutines.
c --------------------------------------------------------------------------
c  11/30/93       
c modified SLAND1 for 2 parallel recess row HJB so as not to calculate T
c field at symmetry line, otherwise gives erroneous results
c search for CTTT for changes made from original file
c
c --------------------------------------------------------------------------

C *****************************************************************************
C **                                                                         **
C **  Subroutine Setprec1                                                    **
C **                                                                         **
C **  SETPREC1:  Allocates perturbed recess pressures on P array.            **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETPREC1

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK), 
     +                QREC1(MAXNPOCK), QIN1, QOUT1

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE COMPLEX U1, V1, P1, PREC1,TREC1,QREC1,QIN1,QOUT1
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,
     +        IFULL,ICSTEP, BEARING
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER IP, J, JSTA, JSTO, K, Kstart
C ----------------------------------------------------------------------------
C --  SETPREC1 code                                                         --
C ----------------------------------------------------------------------------
      Kstart=ISYM+(ISYM-1)*NPA
      IF ((NPOCKET.EQ.0).OR.(BEARING.GE.2)) GOTO 12

          
C    !.................!             !--------
      DO IP=1, NPOCKET               ! UPDATE:
C    !.................!             !--------
                                     ! 
          JSTA=(IP-1)*NXI+NLC        ! Pressure=Prec at interior
          JSTO=IP*NXI+1              ! and edge point
                                     ! of pocket regions
          DO J=JSTA, JSTO            !
              DO K=Kstart, NPA       !
                  P1(J, K)=PREC1(IP) !
              END DO                 !
          END DO                     !
C    !.................!             !
      END DO                         !
C    !.................!             !

C    !-------------------------------! .............................
 11     IF (IFULL.EQ.1) THEN         ! PIxD FULL PAD 
          DO K=Kstart, NPA           ! periodicity condition
              P1(1, K)=P1(NXT, K)    ! .............................
          END DO                     !
        ELSE                         ! .............................
          Kstart=ISYM+(ISYM-1)*NYI   ! PAD BEARING
          DO K=Kstart, NYI           !
            P1(1 , K)=(0.0D0,0.0D0)  !
            P1(NXT,K)=(0.0D0,0.0D0)  !
          END DO                     !
        END IF                       !  
C    !-------------------------------!
c#        RETURN

C    !-------------------------------! SET P-BCS for SEAL/BEARING
 12   RETURN

      END 


C *****************************************************************************
C **                                                                         **
C **  Subroutine Ucoefs                                                      **
C **                                                                         **
C **  UCOEFS:  Calculates Au coefficients  from zero-th order solution.      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE UCOEFS

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /UCOEF/
     + APU(MAXNXT,-MAXNYI:MAXNYI), AWU(MAXNXT,-MAXNYI:MAXNYI),
     + AEU(MAXNXT,-MAXNYI:MAXNYI), ASU(MAXNXT,-MAXNYI:MAXNYI),
     + ANU(MAXNXT,-MAXNYI:MAXNYI), GUO(MAXNXT,-MAXNYI:MAXNYI),
     + GUV(MAXNXT,-MAXNYI:MAXNYI), APUI(MAXNXT,-MAXNYI:MAXNYI),
     + GUPR(MAXNXT,-MAXNYI:MAXNYI),GUPI(MAXNXT,-MAXNYI:MAXNYI),
     + GUT(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ISYM, ICSTEP
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION APU, AWU, AEU, ASU, ANU, GUO, GUV, APUI,
     +                 U,V, DYP, DYV, SVN, SVS, GUPR, GUPI, GUT

      INTEGER INERL,INERP,ITURB,INTER,ICAV,MODEL,ISYM,ICSTEP,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        IFULL, LEFTBC, RIGHTBC, BEARING 

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYPK, C1, C2, Dir, zero, xmin,xmax
      INTEGER J, K, KP1, IP, INERPL, DEVICE, Icsym, NYIP1,
     +        I, KM1, INERPDUM, INERLDUM, Kstart, Kend, Kstep, KV
C ----------------------------------------------------------------------------
C --  UCOEFS code                                                           --
C ----------------------------------------------------------------------------
      zero=0.0D0
      Kstep= ISYM+(ISYM-1)*NYI
      Kstart=ISYM+(ISYM-1)
      IF (BEARING.EQ.2) Kstep=1           ! => SEAL

C    !...........................................! ISYM=0/1

      INERLDUM=INERL                             ! ..........................
      INERPDUM=INERP                             ! SAVE inertia parameters
      IF (INERL.EQ.0) INERP=0                    !
cc                                               ! Set coefficients w/o
      INERPL=INERL+INERP                         ! inertia edge pdrop

C    !...........................................!
      DO I=1, NXT                                ! ..........................
          DO K=Kstep, NYI                        !
              APU(I, K)=0.0D0                    !
              AWU(I, K)=0.0D0                    ! zero array coefficients
              AEU(I, K)=0.0D0                    !
              ANU(I, K)=0.0D0                    !
              ASU(I, K)=0.0D0                    !
              APUI(I,K)=0.0D0                    !
          END DO                                 !
      END DO                                     !

c   !...............................!.start with right side  of bearing
      Kstep=1
      Icsym=1
      Dir=1.0D0
      C1=0.5D0
      C2=0.5D0

      IF (ISYM.EQ.1) GOTO 80


c   !-----------------------------! ISYM=0, ASYMMETRICAL BEARING

      IF (BEARING.EQ.2) THEN      !=> SEAL
         KM1=1                    ! 1
         K=KM1+Kstep              ! 2
         KV=K                     ! 2
         DYPK=DYP(K)              ! size U-CV at K
         GOTO 650
      END IF

 77   IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN  !=> HJB
         GOTO 90

      ELSE IF (BEARING.EQ.3) THEN !=> PLAIN BEARING
         KM1=-2*Kstep             ! -2, 2
         K=Kstep                  !  1, -1
         DYPK=DYP(K)+DYP(-K)      ! size of U-CV at Y=0
         KV=KM1                   ! -2, 2
         GOTO 650
      END IF     

c   !-----------------------------! ISYM=1, SYMMETRICAL BEARING

 80   IF ((NPOCKET.EQ.0).OR.(BEARING.EQ.3)) THEN
         KM1=Kstep      ! 1,-1   => start at Y=0
         K=KM1          ! 1,-1
         DYPK=DYP(K)    !
         KV=K           ! 1,-1
         GOTO 650
      END IF

C   !------------------------------------!-----------------------------------!
C   ! ON LANDs BETWEEN RECESSES           INTER=1
C   !------------------------------------!-----------------------------------!

      
 90   K=Kstep    
      INTER=1
      Kend=Kstep*NPA
      LEFTBC=1+IFULL                             !1: GROOVE, 2: RECESS
      RIGHTBC=2                                  !..........................
C    !...................!                       ! 
      DO IP=1, NPOCKET+1-IFULL                   ! Sweep on lands between
C    !...................!                       ! pockets: Inter=1
          CALL SETLIM(IP, 1)                     ! 
          IF (IP.GT.NPOCKET) RIGHTBC=1
          K=0                                    !
          KM1=Kstart                             !

  100     K=K+Kstep                              ! Dypk= Y size Ucv
          KP1=K+Kstep                            !
          DYPK=DYP(K)                            !

          CALL UCALC(DYPK, INERPL, INERLDUM,     ! Calcule Au coefficients
     +               K,KP1,KM1,K,Kstep,C1,C2,Dir)!    on internal lands

          IF (K.EQ.Kend) goto 200                !

          KM1=K                                  !

          GOTO 100                               !
C    !...................!     
  200 LEFTBC=2
      END DO                                     ! ...........................
C    !...................!                       !    
C                                                !
C     K=Kstep*NPA                                !           at k=Npa+1
      C1=0.25D0                                  ! ......... coefficients for
      C2=0.75D0                                  !           Us interpoltion
      KM1=Kstep*NPA                              !
      K=KM1+Kstep                                !
      DYPK=DYP(K)
      KV=K

C   !------------------------------------!-----------------------------------!
C   ! ON EXTENDED LAND ABOVE RECESSES    ! INTER=0
C   !------------------------------------!-----------------------------------!

  650 INTER=0                                    ! Sweep on extended lands
                                                 !
      CALL SETLIM(IP,INTER)                      ! Inter=0
      Kend=Kstep*NYI                             !
      LEFTBC=1-IFULL                             ! BC's at sides 
      RIGHTBC=LEFTBC
      

  700 CONTINUE

c     .....................!
      IF (K.eq.Kend) THEN
c     .....................!

         IF ((Icsym.eq.1).AND.(MODEL.eq.1) ) THEN  
c       !...................................! Right Side of 2row HJB
            KP1=NYI+1
            DO J=1, NXT
             V(J,KP1)=-V(J,K)            ! Symmetry BC.
             U(J,KP1)=zero               ! dumy value to set AN UN = zero
            END DO
c       !...................................!
         ELSE  
c       !...................................!
            KP1=Kend                     ! MODEL=2 ! single row HJB
c       !...................................!
         END IF
c       !...................................!

         CALL UCALC( DYPK, INERPL, INERLDUM,               
     +               K, KP1, KM1, K, Kstep, C1, C2, Dir )       
         GOTO 800
c     ......................!
      ELSE
c     ......................!
         KP1=K+Kstep

         CALL UCALC( DYPK, INERPL, INERLDUM,               
     +               K, KP1, KM1, KV, Kstep, C1, C2, Dir)       
         C1=0.5D0                                   ! ...........................
         C2=0.5D0                                   ! Us = C1 Us + C2 Up
 
         KM1 = K
         K = KM1 + Kstep
         KV = K
         DYPK = DYP(K)
         
         GOTO 700                                   !
c     ......................!
      END IF
c     .......................! K=Kend & MODELs=1 or 2
 



c !..............................................! Direct for HJB asymmetry

c    !.....................! ISYM=0       
  800 IF (ISYM.eq.0) THEN
c    !.....................!
        IF ((BEARING.EQ.2).OR.(Icsym.eq.2)) GOTO 900

             IF (MODEL.EQ.1) THEN     ! TWO recess ROW HJB
               DO J=1, NXT
                  ANU(J,Kend)=zero
                  GUV(J,Kend)=zero
               END DO
             END IF

        Kstep=-1
        Kstart=-Kstep
        Dir=-1.0D0
        Icsym=2
        GOTO 77            !-> left side of bearing
c    !.....................! ISYM=0       
      END IF
c    !.....................! ISYM=0       


  900 INERP=INERPDUM


      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Ucalc                                                       **
C **                                                                         **
C **  UCALC:  Subroutine Calcule of the BASIC subroutine UCOEFS has been     **
C **          made into a full-blown FORTRAN subroutine.                     **
C **                                                                         **
C *****************************************************************************
                           
      SUBROUTINE UCALC( DYPK, INERPL, INERLDUM,
     +                 K, KP1, KM1, KV, Kstep, C1, C2, Dir)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /UCOEF/
     + APU(MAXNXT,-MAXNYI:MAXNYI), AWU(MAXNXT,-MAXNYI:MAXNYI),
     + AEU(MAXNXT,-MAXNYI:MAXNYI), ASU(MAXNXT,-MAXNYI:MAXNYI),
     + ANU(MAXNXT,-MAXNYI:MAXNYI), GUO(MAXNXT,-MAXNYI:MAXNYI),
     + GUV(MAXNXT,-MAXNYI:MAXNYI), APUI(MAXNXT,-MAXNYI:MAXNYI),
     + GUPR(MAXNXT,-MAXNYI:MAXNYI), GUPI(MAXNXT,-MAXNYI:MAXNYI),
     + GUT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT), SUE(MAXNXT) 
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

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU, ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPE, HRECD
      COMMON /RECASP/ ASPEC(MAXNPOCK)
      COMMON /COMPLIA/ AC, ETAC, RELAXH, LIFT
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASED, NCASE
      COMMON /HONEY/ HCELL,HCDIM

      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /LRBOUND/ LEFTBC,RIGHTBC
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION APU, AWU, AEU, ASU, ANU, GUO, GUV, APUI,
     +                 DXP, DXU, HP, HU, HV, SUW, SUE,
     +                 U, V, P, GUPR, GUPI, GUT, RHOP,EMUP,ASPEC,
     +                 Drhop, Drhot, Demup, Demut,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPE, HRECD,
     +                 AC, ETAC, RELAXH,FREQU,SIGMA,L1,RES,
     +                 HCELL,HCDIM

      INTEGER JUMIN,JUMAX,JUSTART,JUSTOP,JPMIN,JPMAX,JPSTART,JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, 
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,
     +        IFULL, LEFTBC, RIGHTBC, BEARING, LIFT, ICASED, NCASE

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION zero,DYPK, C1, C2, Hs, Hn, Suej, Suwj,
     +                 rhopu,emupu, rhon,rhow,emuw,
     +                 rhos, rhoe, rhor, emue, 
     +                 drhodp, demudp, drhodt, demudt,
     +                 UW, UE, FWU, FEU, VS, VN, VP, DAREA,
     +                 GUOU, GUUU, GUVU,  FSU, FNU, SANB,
     +                 AEUC, CORRE, FWUC, AWUC,
     +                 DYPRE, UN, US, CORRW, FEUC, NUMC,
     +                 Eta, Ax, Renh, Guru, Gumu, dpx, Dir
      INTEGER J, JP1, INERLDUM, Icase,K,KP1,KM1,Kstep,KK,INERPL,
     +        KV 
C ----------------------------------------------------------------------------
C --  UCALC code                                                            --
C ----------------------------------------------------------------------------
c Dir: flow direction, IF Dir=1.0   RIGHT SIDE of HJB
c                            =-1.0  LEFT SIDE of HJB
c ............................................................................
C KV=K except when K=+/- 1 for asymmetric bearing=3
c ............................................................................
C COMPLIANCE EFFECT : Ac/(1+ i ETAC) * P1 
c ............................................................................
C
      zero=0.0D0
      KK=K*Kstep
      NUMC=AC/(1.0D0+ETAC*ETAC)               ! COMPLIANCE coefficient
      DYPRE=DYPK*REY                          ! :::::::::::::::::::::::::::::
                                              !
C    !......................!                 !......
      IF (INTER.EQ.1) THEN                    ! BCs on inter_recess lands:
C    !......................!                 !.....
          j=jumin                             ! LEFT SIDE:
          rhoe=Sue(j)*rhop(j,k)+Suw(j)*rhop(j+1,k)
          rhow=rhop(j,k)
          Fwu=Dypre*rhoe*Hu(j,k)*U(j,k)       !=Feu at left edge
          Uw=Fwu/dypre/rhow/hp(j,k)           ! West Uvel at left edge
          Fwuc=Fwu                            !...
                                              !
          j=jumax                             ! RIGHT SIDE
          rhow=Sue(j)*rhop(j,k)+Suw(j)*rhop(jpstop,k)
          rhoe=rhop(jpstop,k)
          Feuc=Dypre*rhow*Hu(j,k)*U(j,k)      !=Fwu at right edge
          Ue=Feuc/dypre/rhoe/hp(jpstop,k)     ! East Uvel at right edge
          Feu=Feuc*Inerl                      !
          IF ((KK.EQ.NPA).AND.(LEFTBC.EQ.2)) INERL=0       
          Fwu=Fwu*INERL                       !

C     !......................!                !...............
      ELSE                                    ! BCs on extended lands:
C     !......................!                !....
          UW=U(justart, k)                    ! left U velocity = BC
          rhow=rhop(jumin,k)
          Fwu=Dypre*Rhow*Hp(jumin,k)*(Uw+U(jumin,k))/2.0D0
C     !......................!                !
      END IF                                  !
C     !......................!                !......................

      j=jumin-1                               !
      Icase=1                                 !

  301 j=j+1                                   ! Sweep on Ucvs
      if (j.eq.jumax) goto 302                !
      goto 303                                !

  302 IF ((KK.EQ.NPA).AND.(RIGHTBC.EQ.2)) THEN  !.............................!
        INERL=0                               !Neglect inertia on last Ucv
        j=jumax-1                             !
        Apu(j,k)=Apu(j,k)-(Aeu(j,k)+Feu)/Alfu !
        Aeu(j,k)=zero                         !
        Apui(j,k)=zero                        !
        j=jumax                               !
      END IF           !......................!.............................! 
      Icase=2                                 !

C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  303   jp1=j+1                               ! SWEEP ON U CONTROL VOLUMES
        Suej=Sue(j)                           ! AND CALCULATE AUNB COEFFS.
        Suwj=Suw(j)                           !

        Vs=V(j,kv)*Suej+V(jp1,kv)*Suwj        ! Axial velocity at center of
        Vn=V(j,kp1)*Suej+V(jp1,kp1)*Suwj      ! U control volume
        Vp=0.5D0*(Vs+Vn)                      !

        rhopu=rhop(j,k)*Suej+rhop(jp1,k)*Suwj ! Density and viscosity at
        emupu=emup(j,k)*Suej+emup(jp1,k)*Suwj ! center of U control volume   

        drhodp=Drhop(j,k)*Suej+Drhop(jp1,k)*Suwj ! dRho/dP and dEmu/dP at center
        demudp=Demup(j,k)*Suej+Demup(jp1,k)*Suwj ! of U control volume

        drhodt=Drhot(j,k)*Suej+Drhot(jp1,k)*Suwj ! dRho/dT and dEmu/dT at center
        demudt=Demut(j,k)*Suej+Demut(jp1,k)*Suwj ! of U control volume

        dpx=(P(jp1,k)-P(j,k))/DXu(J)          ! dp/dx

        CALL CALCKX1( U(j,k), Vp, Hu(j,k), Speed, Rep,
     +      dpx,rhopu,emupu, guou,guuu,guvu, guru, gumu  )

        Darea=DXu(j)*Dypk                     !
        APU(j,k)=guuu*darea/alfu              ! Gxu: UP1/UP2 coefs.
        GUO(j,k)=-guou*Darea                  ! Gxh: H1/H2 coefs.
        GUV(j,k)=-guvu*Darea                  ! Gxv: V1/V2 coefs.
        
        GUPR(j,k)= GUO(j,k)*NUMC              ! Gxp: COMPLIANCE effect
        GUPI(j,k)=-GUPR(j,k)*ETAC -           ! 
     +            (RES*U(J,K)*HCDIM)*drhodp*DAREA  ! cell depth effect

        GUPR(j,k)=-(guru*drhodp+gumu*demudp)*Darea+GUPR(j,k)
        GUT (j,k)=-(guru*drhodt+gumu*demudt)*Darea  ! GxT: T1/T2 coeff.

c NOTES:
c    Coefficient of dynamic pressure for  bearings with compliance
c    effects is equal to:
c    Gxp = Gxp(real,imag) = Gxh ac/(1+in) + gxrho drho/dp + gxemu demu/dp
c 
c    ac/(1+in)= ac (1-in)/(1+n^2) = numc (1-in)
c
c NOTES: 6/20/95 coefficient for compliance of cell depth on rough seals
c is given as  ix(Res.Uo.HCdim)xRHOj and written also as
c              ix(Res.Uo.HCdim)x(drho/dp.Pj+drho/dT.Tj)
c              here I have neglected the effect of temperature
c              otherwise I will require a GVTimaginary coeffic.


C      !.......................               ! Account for inertia:   
          IF (INERL.EQ.1) THEN                ! Effects on lands:
C      !.......................               ! 
C           Imaginary part of APu coefficient
            APUI(j,k)=Rhopu*Darea*RES*Hu(j,k)/Alfu 

            AWU(j,k)=DMAX1(Fwu,zero)          !

            rhoe=rhop(jp1,k)
            Feu=Dypre*rhoe*Hp(jp1,k)*(U(j,k)+U(jp1,k))/2.0D0
            AEU(j,k)=DMAX1(-Feu,zero)

            Hs=Hv(j,kv)*suej+hv(jp1,kv)*Suwj  ! At south Ucv face
            rhos=((rhop(j,km1)+rhop(j,k))*Suej +
     +            (rhop(jp1,km1)+rhop(jp1,k))*Suwj)/2.0D0
            Fsu=Rey*rhos*DXu(j)*Hs*Vs*Dir     !
            ASU(j,k)=DMAX1(Fsu,zero)          !

            Hn=Hv(j,kp1)*suej+hv(jp1,kp1)*suwj! At north Ucv face
            rhon=((rhop(j,k)+rhop(j,kp1))*suej +
     +            (rhop(jp1,k)+rhop(jp1,kp1))*suwj)/2.0D0
            Fnu=Rey*rhon*DXu(j)*Hn*Vn*Dir    !
            ANU(j,k)=DMAX1(-Fnu,zero)         !

            Sanb=AEU(j,k)+AWU(j,k)+ASU(j,k)+ANU(j,k)+(Feu-Fwu)

            APU(j,k)=APU(j,k)+Sanb/Alfu       !

            Fwu=Feu                           !

            Un=0.50D0*(U(j,kp1)+U(j,k))       !
            Us=C1*U(j,km1)+C2*U(j,k)          !

            GUV(j,k)=GUV(j,k)+(Us-Un)*Dir*Dxu(j)*Hu(j,k)*Rey*Rhopu

C      !................!                     ! 
          ELSE
C      !................!                     ! 
          INERL=INERLDUM                      ! reset inertia parameter
C      !................!                     ! 
          END IF                              ! Inerl=1
C      !................!                     !

        GOTO (301, 500) , Icase

C --------------------------------------------------------------


  500 IF ( (INTER. EQ.0). OR.
     +     (INERPL.EQ.0). OR.
     +     (KK.EQ.NPA)       ) RETURN

C  ...........................................! Modify W&E flows on inter-rec
      Aeuc=DMAX1(-Feuc,zero)                  ! regions
      Corre=(Aeuc-AEU(jumax,k)+Feuc-Feu)/Alfu !
      APU(jumax,k)=APU(jumax,k)+Corre         ! Corrected APU & AEU on
      AEU(jumax,k)=Aeuc                       ! east boundary

      Awuc=DMAX1(+Fwuc,zero)                  !
      Fwu=Fwuc*INERL                          !
      Corrw=(Awuc-AWU(jumin,k)+Fwu-Fwuc)/Alfu !
      APU(jumin,k)=APU(jumin,k)+Corrw         ! Corrected APU & AWU on
      AWU(jumin,k)=Awuc                       ! west boundary

C ........................................... ! .............................
  600 IF (INERP.EQ.0) RETURN                   

      IF ((UW.GT.zero).AND.(LEFTBC.EQ.2))THEN ! Correct for Edge PRES drop
c    !.....................!                  ! on west (LEFT) boundary
      Eta=Hp(jumin,k)/(hrecd+hp(jumin,k))     !
      Rhow=rhop(jumin,k)    
      Emuw=emup(jumin,k)
      Rhor=Rhop(justart,K)                     
      Ax=klosxd*(1.0D0-(Eta*Rhow/Rhor)**2)*Rhow      

C## removed on 10/5/95      
C##      IF (SWIRL.GT.zero ) THEN   !...............! 
C##        Renh=Renc*Rhow*Hp(jumin,k)/Emuw
C##        Ax=Ax*(1.0D0+1.95/(Renh**.43))
C##      END IF                     !...............!

      AWU(jumin,k)=Awuc-2.*Ax*dypk*Hu(jumin,k)*Uw

      END IF                                  ! ....................... &
c   !.......................!                 !
 
      IF ((UE.LT.zero).AND.(RIGHTBC.EQ.2))THEN! Correct for Edge PRES drop
c   !........................!                ! on east (RIGHT) boundary
      Eta=Hp(jpstop,k)/(hrecd+hp(jpstop,k))   !
      Rhoe=rhop(jpstop,k)                     !
      Rhor=Rhop(JPSTOP+1,K)                   ! 
      Ax=klosxu*(1.0D0-(Eta*Rhoe/Rhor)**2)*Rhoe

      AEU(jumax,k)=Aeuc+2.*Ax*dypk*hu(jumax,k)*Ue
  
      END IF                                  !
c   !.......................!                 ! 
C                                             ! :::::::::::::::::::::::::::::
C      INERL=INERLDUM                         ! RESET VALUE

      END 


C *****************************************************************************
C **                                                                         **
C **  Subroutine Vcoefs                                                      **
C **                                                                         **
C **  VCOEFS:  Calculates Av coefficients from zero-th order solution.       **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE VCOEFS

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /VCOEF/
     + APV(MAXNXT,-MAXNYI:MAXNYI),AWV(MAXNXT,-MAXNYI:MAXNYI),
     + AEV(MAXNXT,-MAXNYI:MAXNYI),ASV(MAXNXT,-MAXNYI:MAXNYI),
     + ANV(MAXNXT,-MAXNYI:MAXNYI),GVO(MAXNXT,-MAXNYI:MAXNYI),
     + GVU(MAXNXT,-MAXNYI:MAXNYI),APVI(MAXNXT,-MAXNYI:MAXNYI),
     + GVPR(MAXNXT,-MAXNYI:MAXNYI),GVPI(MAXNXT,-MAXNYI:MAXNYI),
     + GVT(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /PreLR/ Ple,Pri, Csel, Cser

      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c     ...............................................................
      DOUBLE PRECISION APV, AWV, AEV, ASV, ANV, GVO, GVU, APVI,
     +                 GVPR, GVPI, GVT, DYP, DYV, SVN, SVS,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD,
     +                 Ple,Pri, Csel, Cser

      INTEGER JVMIN, JVMAX, JVSTART, JVSTOP,      
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, ISYM, ICSTEP,
     +        NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER Kstep,Kstart,Kend,KVM1,Icsym,J,I,K,KM1,KP1,
     +        INERPL,IP,INERPDUM,INERLDUM
      DOUBLE PRECISION zero, Cseal,KLOSY, SVNK, SVSK, DYVK, Dir,
     +                 C1V,C2V, xmin,xmax

C ----------------------------------------------------------------------------
C --  VCOEFS code                                                           --
C ----------------------------------------------------------------------------
      zero=0.0D0
      C1V=0.50                                ! Vs = C1 VS + C2 VP
      C2V=0.50                                ! COEFS FOR INTERPOLATION OF VS

      Kstart=ISYM+(ISYM-1)*NYI
      IF (BEARING.EQ.2) Kstart=2              !=> SEAL


      INERLDUM=INERL                          ! Save inerl
      INERPDUM=INERP                          ! Save inerp
      IF (INERL.EQ.0) INERP=0                 ! ..............................
c                                             ! Set Vcoefs without
      INERPL=INERP+INERL                      ! inertia effect at edge

      DO I=1, NXT                             ! ..............................
          DO K=Kstart,NYI                     !
              APV(I, K)=ZERO                  !
              AWV(I, K)=ZERO                  !
              AEV(I, K)=ZERO                  ! Zeroes matrix coeficients
              ANV(I, K)=ZERO                  !
              ASV(I, K)=ZERO                  !
              APVI(I,K)=ZERO                  !
          END DO                              !
      END DO                                  !


C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Kstep=1                            !
      Icsym=1                            ! AT: 
      Dir=1.0D0                          ! Right side of BEARING
      KLOSY=KLOSYr                       ! entrance coefficient &
      Cseal=Cser                         !  & end seal param.
c    !...................................!...........................

      IF (ISYM.EQ.1) GOTO 88

c    !...................................! ASYMMETRIC BEARING


       IF (BEARING.EQ.2) THEN            ! => SEAL
         KM1=Kstep                       ! 1
         K=KM1+Kstep                     ! 2
         KP1=K+Kstep                     ! 3
         KVM1=KM1                        ! 1
         Dyvk=DYV(K)
         C1V=1.0D0
         C2V=0.0D0
         IF (INERPDUM.EQ.1) THEN         ! 
            INERP=INERPDUM               !
            INERL=1                      ! => inertia of 1st CV
            INERPL=INERP+INERL           !    at seal inlet
         END IF
         GOTO 79 
       END IF

 77    CONTINUE

       IF (BEARING.EQ.1) THEN            !=> HJB
         K=Kstep                         !  1,-1
         KM1=-Kstep                      ! -1, 1
         KVM1=-2*Kstep                   ! -2, 2
         KP1=2*Kstep                     !  2, -2
         Dyvk=2.0D0*DYV(Kstep)           !
         INTER=1                         ! 
       ELSE IF (BEARING.EQ.3) THEN       !=> JOURNAL BEARING
         KM1=-Kstep                      ! -1, 1
         K=KM1-Kstep                     ! -2, 2
         KP1=K-Kstep                     ! -3, 3
         KVM1=-Kstep                     ! -1, 1
         Dyvk=DYV(K)
         C1V=1.0D0
         C2V=0.0D0
         Dir=-Dir  
       END IF                

 79      SvnK=SVN(K)                 
         Svsk=SVS(K)

c      !..................................!
         IF (NPOCKET.GT.0) THEN
c      !..................................!
        DO IP=1, NPOCKET+1-IFULL
           CALL SETLIM(IP, INTER)
           CALL VCALC(DYVK,INERPL,K,KP1,KM1,KVM1,SVNK,SVSK,
     +                Dir,KLOSY,Cseal, C1V, C2V)
        END DO
c      !..................................!
          ELSE
c      !..................................! SEAL/PLAIN BEARING
           INTER=0
           CALL SETLIM(IP, INTER)
           CALL VCALC(DYVK,INERPL,K,KP1,KM1,KVM1,SVNK,SVSK,
     +                Dir,KLOSY,Cseal, C1V, C2V)
c      !..................................!
         END IF
c      !..................................!

c      !.............................!
        IF (BEARING.EQ.2) THEN       ! => SEAL
c      !.............................!
           C1V=0.50D0
           C2V=C1V
           K=KM1+Kstep               !=2, start with K=K+1=3
           INERL=INERLDUM            ! Reset params at seal
           INERP=INERPDUM            ! inlet
           IF (INERL.EQ.0) INERP=0   !
           INERPL=INERL+INERP        ! 
           GOTO 655

c      !.............................!
        ELSE IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN  ! => HJB 
c      !.............................!
          INTER=1
          GOTO 99
C##??? if npocket=0 then K=???

c      !.............................!
        ELSE IF (BEARING.EQ.3) THEN  !=> PLAIN BEARING
c      !.............................!
              K=Kstep                !  start with K=K+Kstep=2
              Dir=-Dir
              GOTO 655   
c      !.............................!
        END IF
c      !.............................!


C    !-----------------------!...........! SYMMETRIC BEARING

 88    IF ((NPOCKET.EQ.0).OR.(BEARING.EQ.3)) THEN ! => SYM JOURNAL BEARING
          K=Kstep
          C1V=1.0D0    !###
          C2V=0.0D0    !###
          GOTO 650
       ELSE IF (BEARING.EQ.1) THEN
          INTER =1
          GOTO 99
       END IF
 

 99   K=Kstep  

      Kend=Kstep*NPA

c   !------------------------!.................!-------------------------! 
      DO IP=1, NPOCKET+1-IFULL                 ! Sweep on lands between
C   !------------------------!                 ! pockets: Inter=1
           CALL SETLIM(IP, INTER)              ! 
           KM1=Kstep                           !

 100       K=KM1+Kstep                         ! Note: K starts at k=2*Kstep
           KP1=K+Kstep                         !
           SVNK=SVN(K)                         !
           SVSK=SVS(K)                         ! Calculate Av coefficients
           DYVK=DYV(K)                         ! on internal lands

           CALL VCALC(DYVK,INERPL,K,KP1,KM1,KM1,SVNK,SVSK,
     +                Dir, KLOSY, Cseal, C1V, C2V)

           IF (K.EQ.Kend) goto 200             !
           KM1=K                               !
           GOTO 100                            !
 200  END DO                                   ! END of sweep ON
C    !----------------------!                  ! land between recesses.

      INTER=0                                  !
C                                              ! -------------------      
 250  CALL SETLIM(IP,INTER)                    
      IF (IFULL.EQ.1) THEN  
          JVMAX=NXT                          
          JVSTOP=2
      END IF

      KM1=Kstep*NPA                      ! ON V-CV just above recesses
      K=KM1+Kstep                        !
      KP1=K+Kstep                        ! 
      DYVK=DYV(K)                        ! 
      SVNK=SVN(K)                        ! 
      SVSK=SVS(K)                        !
      CALL VCALC(DYVK,INERPL,K,KP1,KM1,KM1,SVNK,SVSK,
     +              Dir, KLOSY, Cseal, C1V, C2V)

C   !------------------------------------!-----------------------------------!
C   ! ON EXTENDED LAND ABOVE RECESSES
C   !------------------------------------!-----------------------------------!

 650  INTER=0                            

      CALL SETLIM(IP,INTER)                     

C----!      
 655  Kend=Kstep*NYI
C----!      

 700  KM1=K                                   !
      K=KM1+Kstep                             !
      SVNK=SVN(K)                             !
      SVSK=SVS(K)                             ! Calcule Av coefficients
      DYVK=DYV(K)                             ! on extended lands

c     ......................!
      IF (K.eq.Kend) THEN
c     ......................! at exit/end sides of bearing

         IF ((Icsym.eq.1).AND.(MODEL.eq.1) ) THEN 
c       !...................................! Right Side of 2row HJB
            KP1=NYI+1
c       !...................................!
         ELSE  
c       !...................................!
            KP1=Kend                        ! MODEL=2 ! single row HJB
c       !...................................!
         END IF
c       !...................................!

         CALL VCALC(DYVK,INERPL,K,KP1,KM1,KM1,SVNK,SVSK,
     +              Dir, KLOSY, Cseal, C1V, C2V)

         GOTO 800
c     .......................!
      ELSE
c     .......................! at internal extended CVs

         KP1=K+Kstep                             
         CALL VCALC(DYVK,INERPL,K,KP1,KM1,KM1,SVNK,SVSK,
     +              Dir, KLOSY, Cseal, C1V, C2V)

         C1V=0.5D0
         C2V=C1V

         GOTO 700                                
c     .......................!
      END IF
c     .......................!


c.............................................!... GOTO Left side of HJB

C     !.............................! ASYMMETRIC BEARING
  800  IF (ISYM.eq.0) THEN
C     !.............................!

         IF ((BEARING.EQ.2).OR.(Icsym.eq.2)) GOTO 900

         Kstep=-1
         Dir=-1.0D0
         Cseal=Csel                      ! Left side end seal
         KLOSY=KLOSYl                    ! & entrance edge coeff.
         Icsym=2
         GOTO 77                         ! => LEFT SIDE of BEARING
C     !.............................!
       END IF
C     !.............................! ASYMMETRIC BEARING

c.............................................!............................!

  900  INERP=INERPDUM
      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Vcalc                                                       **
C **                                                                         **
C **  VCALC:  Subroutine Calcule of BASIC subroutine Vcoefs has been made    **
C **          into a full-blown FORTRAN subroutine.                          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE VCALC(DYVK,INERPL,K,KP1,KM1,KVM1,SVNK,SVSK,
     +                 Dir,KLOSY,Cseal, C1, C2)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT), SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /URECJET/UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /VCOEF/
     + APV(MAXNXT,-MAXNYI:MAXNYI),AWV(MAXNXT,-MAXNYI:MAXNYI),
     + AEV(MAXNXT,-MAXNYI:MAXNYI),ASV(MAXNXT,-MAXNYI:MAXNYI),
     + ANV(MAXNXT,-MAXNYI:MAXNYI),GVO(MAXNXT,-MAXNYI:MAXNYI),
     + GVU(MAXNXT,-MAXNYI:MAXNYI),APVI(MAXNXT,-MAXNYI:MAXNYI),
     + GVPR(MAXNXT,-MAXNYI:MAXNYI),GVPI(MAXNXT,-MAXNYI:MAXNYI),
     + GVT(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV, ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /COMPLIA/ AC, ETAC, RELAXH, LIFT
      COMMON /FREQ/ FREQU, SIGMA, L1, RES, ICASED, NCASE
      COMMON /HONEY/ HCELL,HCDIM

      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /BTYPE/ BEARING

c     ................................................................
      DOUBLE PRECISION APV, AWV, AEV, ASV, ANV, GVO, GVU, APVI,
     +                 GVPR, GVPI, GVT, DYP, DYV, SVN, SVS,
     +                 DXP, DXU, HP, HU, HV , UREC, SUW, SUE,
     +                 U, V, P, RHOP, EMUP,Drhop, Demup, Drhot, Demut,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD,
     +                 AC, ETAC, RELAXH, FREQU,SIGMA,L1,RES,
     +                 HCELL,HCDIM

      INTEGER JVMIN, JVMAX, JVSTART, JVSTOP,      
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        BEARING, LIFT, ICASED,NCASE
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYVRE, VEE(MAXNXTP2), VWW, UW, FWV, UE, DAREA,
     +                 FEV, FSV, FNV, SANB, DYPK, DYPP,
     +                 drhodp,demudp, drhodt, demudt, Uedge,
     +                 DIV, FEW, DPV, FSVO, VS, VN, ASVV, CORRS,
     +                 GVOU, GVUU, GVVU, GVRU, GVMU, dpy,
     +                 SVSK, SVNK, DUMY, DYVK, KYO, suej, suwj,
     +                 rhopv,emupv,rhos,emus,rhow, rhor, Ay, C1,C2,
     +                 KLOSY,Dir,rhon, rhoe, He, Hw,Eta,Cseal, NUMC

      INTEGER J,JP1, JL, JR, K, KM1,KVM1, KP1, I, INERPL, Icase, II
C ----------------------------------------------------------------------------
C --  VCALC code                                                            --
C ----------------------------------------------------------------------------
      DYVRE=DYVK*REY                            ! Calcule Av coefficients
C                                               ! :::::::::::::::::::::::::::
      DO J=JVMIN, JVMAX-1                       ! V velocities at Vcv fases
          VEE(J)=0.5D0*(V(J,K)+V(J+1, K))       ! ...........................
      END DO                                    !
      VWW=0.5D0*(V(JVSTART,K)+V(JVMIN,K))       ! <--- Vwest
      VEE(JVMAX)=0.5D0*(V(JVMAX,K)+V(JVSTOP,K)) ! ---> Veast
C                                               ! 
      j=jvstart                                 !west Uvel on first Vcv
      Uw=U(j,k)*Svsk+U(j,km1)*Svnk              !
      Hw=Hu(j,k)*Svsk+Hu(j,km1)*Svnk            !
      rhow=(rhop(j,k)*Sue(j)+rhop(jvmin,k)*Suw(j))*Svsk +
     +     (rhop(j,km1)*Sue(j)+rhop(jvmin,km1)*Suw(j))*Svnk 
      Fwv=dyvre*Hw*Uw*rhow*INERL                ! West Flow
C                                               ! 
      NUMC=AC/(1.0+ETAC*ETAC)                   ! compliance coefficient
c                                               ! with loss coefficient

C                                               ! ...........................
C !---------------------------!                 ! Sweep on V equation
      DO J=JVMIN, JVMAX                         ! ...........................
C !---------------------------!                 ! 

          UE=U(J, K)*SVSK+U(J, KM1)*SVNK        ! East U velocity on Vcv

          rhopv=rhop(j,k)*Svsk+rhop(j,km1)*Svnk ! density and viscosity at
          emupv=emup(j,k)*Svsk+emup(j,km1)*Svnk ! center of V control volume

          drhodp=drhop(j,k)*Svsk+drhop(j,km1)*Svnk  ! dRho/dP and DEmu/dP
          demudp=demup(j,k)*Svsk+demup(j,km1)*Svnk  ! at V-cv

          drhodt=drhot(j,k)*Svsk+drhot(j,km1)*Svnk  ! dRho/dT and DEmu/dT
          demudt=demut(j,k)*Svsk+demut(j,km1)*Svnk  ! at V-cv

          Dpy=Dir*(P(j,k)-P(j,km1))/DYVK        ! dp/dy
                                                ! ............. 1st shear coef
          CALL CALCKY1(0.5D0*(UE+UW), V(J, K), Hv(J,K), SPEED, REP, 
     +      dpy, rhopv, emupv, GVOU, GVUU, GVVU, GVRU, GVMU) 

          DAREA=DXP(J)*DYVK                     !
          APV(J, K)=GVVU*DAREA/ALFV             ! Gyv: coef. of V1 or V2
          GVO(J, K)=-GVOU*DAREA                 ! Gyh: coef. of H1 or H2
          GVU(J, K)=-GVUU*DAREA                 ! Gyv: coef. of U1 or U2

          GVPR(J,K)=GVO(J,K)*NUMC               ! Gyp: compliance effect
          GVPI(J,K)=-ETAC*GVPR(J,K) -            
     +               (RES*V(J,K)*HCDIM)*drhodp*DAREA  ! cell depth effect
          GVPR(J,K)=-(GVRU*drhodp+GVMU*demudp)*DAREA+GVPR(J,K)
          GVT (J,K)=-(GVRU*drhodt+GVMU*demudt)*DAREA   ! GyT coef. of T1 or T2

c NOTES:
c    Coefficient of dynamic pressure for bearing with compliance
c    effects is equal to:
c    Gyp = Gyp(real,imag) = Gyh ac/(1+in) + gyrho drho/dp + gyemu demu/dp
c
c    ac/(1+in)= ac (1-in)/(1+n^2) = numc (1-in)
c -
c 6/20/95 coefficient for compliance of cell depth on rough seals
c is given as  ix(Res.Vo.HCdim)xRHOj and written also as
c              ix(Res.Vo.HCdim)x(drho/dp.Pj+drho/dT.Tj)
c              here I have neglected the effect of temperature
c              otherwise I will require a GVTimaginary coeffic,


C                                               ! ---------------------------
C      ........................                 ! Account for inertia
          IF (INERL.EQ.1) THEN                  ! Effects on lands:
C      ........................                 !
C         IMAGINARY part ot APv coefficient:
          APVI(j,k)=Rhopv*RES*DAREA*Hv(j,k)/ALFV   

          Vs=C2*V(j,k)+C1*V(j,kvm1)             ! On Vcv south face
          Fsv=Dir*Rey*dxp(j)*Hp(j,km1)*Vs*rhop(j,km1)
          ASV(j,k)=DMAX1(+Fsv,0.0D0)            !....

          Vn=(V(j,k)+V(j,kp1))/2.0D0            ! On Vcv north face
          Fnv=Dir*Rey*dxp(j)*Hp(j,k)*Vn*rhop(j,k)   
          ANV(j,k)=DMAX1(-Fnv,0.0D0)            !....

          He=Hu(j,k)*Svsk+Hu(j,km1)*Svnk        ! On Vcv east face

          IF (j.eq.jvmax) THEN 
            jp1=jvstop                          !
          ELSE
            jp1=j+1                             !
          END IF

          Suej=sue(j)                           !
          Suwj=suw(j)                           !
          rhoe=(rhop(j,k)*Suej+rhop(jp1,k)*Suwj)*Svsk +
     +       (rhop(j,km1)*Suej+rhop(jp1,km1)*Suwj)*Svnk 
          Fev=dyvre*Ue*He*rhoe                  !
          AEV(j,k)=DMAX1(-Fev,0.0D0)            !

          AWV(j,k)=DMAX1(+Fwv,0.0D0)            ! 

          Sanb=AEV(j,k)+AWV(j,k)+ASV(j,k)+ANV(j,k)+ (Fnv-Fsv)    

          APV(J, K)=APV(J, K)+SANB/ALFV         ! coef. of Vp = Apv/alfv

           FWV=FEV                              !

           GVU(J,K)=GVU(J,K)+DYVRE*Hv(j,k)*(VWW-VEE(J))*Rhopv !  coef. of U1 or U2

           VWW=VEE(J)                           !
C     ............................              !
          END IF                                ! Inerl=1
C     ............................              !

          UW=UE                                 !
C   !.....................!                     !
      END DO                                    ! j = jvmin, jvmax

C --------------------------------------------- !


C    !................................................! 	
  200 IF (((INERPL).EQ.0).OR.(BEARING.GE.3)) THEN
C    !................................................! 	
C     No edge inertia effects or PLAIN BEARING
c
                 GOTO 777 ! => EXIT


C    !................................................! 	
      ELSE
     +IF ((BEARING.EQ.1).AND.(IABS(K).EQ.NPAP1).AND.
     +    (NPOCKET.GT.0)) THEN                        ! MODIFY Eqns when REC edge
C    !................................................!is at bottom of Vcv
C     FOR HYDROSTATIC BEARINGS:

          DYPK=DYP(KM1)                         ! ...........................
          DYPP=DYPK/2.0D0                       ! Y size of 1/2 Pcv at edge
          DIV=2.0D0                             !
C        !....................!                 !
          DO I=1, NPOCKET                       ! Correct at recess I:
C        !....................!                 !
              JR=I*NXI+1                        !
              J=(I-1)*NXI+NLC                   !
              II=0                              !
c                                               ! REC Left corner
              II=II+1                             !
              Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KM1))
              FEW=DYPK*Uedge*(HU(J,KM1)-HP(J,KM1))     
              Icase=1 
              goto 704                          !

 701          Icase=2   
              DIV=1.0D0                         !

 702          J=J+1                             ! at interior V-cvs
              IF(J.EQ.JR) goto 703              !
              II=II+1
              Uedge=UREC(I,II)*(1.0D0+HRECD/HU(J,KM1))
              FEW=DYPP*Uedge*(HU(J,KM1)-HU(J-1,KM1))
              goto 704

 703          Icase=3                           !
              DIV=2.0D0                         ! Rec right corner
              II=II+1
              Uedge=UREC(I,II)*(1.0D0+HRECD/HP(J,KM1))
              FEW=DYPK*Uedge*(HP(J,KM1)-HU(J-1,KM1))   

C                                               ! -----------------------------
 704          Vs=(V(j,k)+V(j,km1))/2.0D0        ! Corrects south coefficients
              rhos=rhop(j,km1)
              emus=emup(j,km1) 
              rhon=rhos*Svnk+rhop(j,k)*Svsk

              FSVO=Dir*REY*DXP(J)*HP(J,KM1)*Vs*rhos*INERL   ! wrong south flow

              Fnv=DXP(J)*Hv(J,K)*V(J,K)*Rhon
              Fsv=Fnv+Dir*Rhos*Few                          ! Actual South flow
              Vs=Fsv/DXP(J)/HP(J,KM1)/Rhos                  ! ''''' Vsouth

              Fsv=Dir*Fsv*Rey/Div                           ! corrected Sflow
              Asvv=DMAX1(Fsv,0.0D0)                         ! ...............
              Corrs=(Asvv-Asv(J,K)+Fsvo-Fsv)/Alfv

              APV(J, K)=APV(J, K)+CORRS          
              ASV(J, K)=ASVV                    

c            !.........................................!
              IF ((INERP.EQ.1).AND.(Dir*VS.GT.(0.0D0))) THEN ! Flow out of pocket
                  Eta=Hp(J,KM1)/(Hrecd+Hp(J,KM1))            ! set Bernoulli effect
                  Rhor=RHOP(JR-1,1)                          ! 
                  AY=KLOSY*RHOS*(1.0D0-(ETA*Rhos/Rhor)**2)   !         
                  ASV(J,K)=ASVV-2.0D0*AY*DXP(J)*HV(J,K)*Dir*VS 
              END IF                                         !Inerp=1 & Vs>0.
c            !.........................................!


          GOTO (701,702,705), Icase

C       !.......................................!.................. 
  705     CONTINUE                              ! with next recess

C       !.....................!                 !
          END DO                                ! -----------------------------
c       !.....................!                 !
    
         IF (IFULL.EQ.1) THEN
          APV(1, K)=APV(NXT, K)                 ! Periodicity condn.
          ASV(1, K)=ASV(NXT, K)                 !
         END IF


C   !.............................................!
      ELSE IF ((BEARING.EQ.2).AND.(K.EQ.2).AND.
     +          (INERP.EQ.1) ) THEN               !Modify eqns. at SEAL
C   !.............................................!inlet
C    FOR SEALS:

C                                               ! -----------------------------
      DO J=Jvmin, Jvmax                         !
c    !..................!                       !
        VS=V(j,1)                               !
        Rhos=Rhop(j,1)                          !
c            !.........................................!
              IF (VS.GT.(0.0D0)) THEN                  ! Inflow to Seal
                  AY=Klosy*Rhos
                  ASV(J,K)=ASV(j,k)-2.0*AY*DXP(J)*HV(J,2)*VS
              END IF                                   !Inerp=1 & Vs>0.
c            !.........................................!

c     !..............!
       END DO 
c     !..............! X-sweep on first CV

C   !.............................................!
      END IF                                      ! END corrections
C   !.............................................!


c ................................................!............................!
c                                                 ! Modify eqn. for end seal
 777  CONTINUE

c##      IF ((Cseal.gt.(0.0D0)).AND.(IABS(K).eq.NYI)) THEN
c##     !.................................................!
c##      DO J=JVMIN, JVMAX
c##       VN=V(J,K)
c##       IF ((Dir*VN).gt.(0.0D0)) THEN
c##        ANV(J,K)=ANV(J,K)-2.0D0*Hv(J,K)*DXP(J)*Cseal*Rhop(j,K)*Dir*VN
c##       END IF
c##      END DO
c##       ANV(JVSTOP,K)=ANV(J,JVMIN)
c##      END IF
c    !..............!   Cseal>0 & K=Nyi

      END 


C *****************************************************************************
C **                                                                         **
C **  Subroutine Sland1                                                      **
C **                                                                         **
C **  SLAND1:  Iterative solution of first order flow field on land-sills    **
C **           regions.                                                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SLAND1(FU,FP,YO,RES,SIGMA,DEVICE,IFM)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT), 
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK),QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /DPUV1/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA1/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK), 
     +                QREC1(MAXNPOCK), QIN1, QOUT1
      COMMON /RECASP/ ASPE(MAXNPOCK)
      COMMON /PRECES1/ PRECL1, PRECR1
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /THERMID/ ISOTH   
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /RECJET1/ PREC1do(MAXNPOCK),PREC1up(MAXNPOCK),
     +                 PR1jet(MAXNPOCK,MAXNPOCK+2)

      COMMON /SOURCEA/ PRATIO, CORIF, SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV, ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /LOSPAD/ LOSleadP, KLOSPad
      COMMON /PreLR/ Ple,Pri, Csel, Cser
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP

      COMMON /KVALB/ K, KM1, KP1
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /SOLN/ ISOLN
      COMMON /BTYPE/ BEARING
c     .................................................................
      DOUBLE PRECISION XP,XU,YP,YV,DYP, DYV, SVN, SVS, U, V,ASPE,
     +                 PREC,TREC,QREC,QIN,QOUT,QFACTOR,PRECL,PRECR,
     +                 DYVK, DYPK, SVNK, SVSK, TRECL,TRECR,
     +                 PRATIO, CORIF, SMASS, MPEPS,PREPS, MMP, SFLOW,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC, ASPEC, HRECD,
     +                 Ple,Pri, Csel, Cser,LOSleadP, KLOSPad,
     +                 PMAX, PMAXP, PEPS, ALFT,UC,TC,EC,
     +                 PRECdo, PRECup, PRjet


      DOUBLE COMPLEX U1, V1, P1, DU, DV, DVV,PREC1,QREC1,QIN1, QOUT1,
     +               A, B, C, D, PRECL1, PRECR1, zero, U1e, U1w,
     +               T1,TREC1, PREC1do, PREC1up, PR1jet

      INTEGER K, KM1, KP1, ITER, ITMAX, ITPMAX, ISOTH,
     +        JVMIN, JVMAX, JVSTART, JVSTOP, ISOLN,
     +        JUMIN, JUMAX, JUSTART, JUSTOP,
     +        JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, BEARING,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,IFULL,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP, LEFTBC, RIGHTBC,
     +        COUNTP, MAXCOUNTP ,KEM1,KEM2,KEM3 
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION YO, YPO, YVN, YVS, RES, SIGMA, Cseal,C1U,C2U,
     +                 Factor, KLOSY,QMAX, sold, relaxp, relaxv, Dir,
     +                 FU(MAXNXT),FP(MAXNXT),FUY(MAXNXT),FPY(MAXNXT),
     +                 FPYN(MAXNXT),FPYS(MAXNXT), relaxt 
      INTEGER Icsym, NYIP1, Kstep, Kstart, KVM1, Kend , KV, KM1U,
     +        Iconv, I, J, INERPDUM,INERLDUM, INERPL, DEVICE, IFM,
     +        ITdelta, ITstop 
      DOUBLE COMPLEX P1S(MAXNXT)
C ----------------------------------------------------------------------------
C --  SLAND1 code                                                           --
C ----------------------------------------------------------------------------
C LEFTBC & RIGHTBC: 0 (IFULL=1) FLOW SPECIFIED
C                   1 GROOVE AT BOUNDARY, 2 RECESS AT BOUNDARY
C ............................................................................
      zero=(0.0D0,0.0D0)    ! 0+i0
      relaxp=alfp           ! save nominal values of 
      relaxv=alfv           ! relaxation parameters
      relaxt=alft           !
      Iconv=0               !.....
      C1U=0.5D0             ! Us = C1 US + C2 UP
      C2U=0.5D0             ! coefficients for Us vel. interpolation
c                           !.....
      DO J=1, NXT           ! SAVE film thickness perturbation
        FUY(J)=FU(J)        ! vectors
        FPY(J)=FP(J)        !
        FPYN(J)=FP(J)       !
        FPYS(J)=FP(J)       !
      END DO                !...............................

c     ..................................! SET for TWO RECESS-ROW HJB
      IF (MODEL.EQ.1) THEN
C     ..................................! HJB 2ROW
         NYIP1=NYI+1                    ! symmetry line
         DO J=1, NXT                    ! axial flow=0
            V (J,NYIP1)=0.0D0           ! 
            V1(J,NYIP1)=zero            !
            P1(J,NYIP1)=zero            !
            T1(J,NYIP1)=zero            ! 
         END DO                         !
c     ..................................! for SEALS/JBEARINGS
      ELSE IF (BEARING.EQ.2) THEN       ! COUNT Number of unknowns
C     ..................................! (pressure)
       Maxcountp=(NXT-2+IFULL)*(NYI-1)  !      
       IF (Cser.GT.0.0D0) Maxcountp=Maxcountp+NXT-1
C     ..................................!
      ELSE IF (BEARING.EQ.3) THEN
C     ..................................!
       IF (ISYM.EQ.1) THEN
            Maxcountp=(NXT-2+IFULL)*(NYI-1)
       ELSE IF (ISYM.EQ.0) THEN
            Maxcountp=2*(NXT-2+IFULL)*(NYI-1)
       END IF
C     ..................................!
      END IF
c     ..................................! MODEL=2

c NOTES:  Convergence for HJBS based on orifice flows, 
c         """""     for SEALS and BEARINGS based on Pressure diffs.
c.............................................................................

      PEPS=SFLOW    !=> convergence param for BEARING=2 & 3
      PMAX=1.0D0
      ITER=0                             
      INERLDUM=INERL                     ! Save inertia parameters 
      INERPDUM=INERP                     ! ...................................
      INERPL=INERP+INERL                 ! Correct pedge1 at Edgeu1 & Edgev1
c                                        ! set & save parameters
C ............................................................................
      IF (INERL.EQ.0) INERP=0            !
C ............................................................................

  100 ITER=ITER+1                        ! Start of iterative process
      MMP=0.0D0                          ! --------------------------
      SMASS=0.0D0                        !
      PMAXP=0.0D0                        !
      PMAX=DMAX1(1.0D0,PMAX)             ! modifed
      PEPS=SFLOW*PMAX*BEARING**2         ! dynamic convergence param.

c NOTE: This means that for JBearings, PEPS=9xPEPS
c       i.e. an order of magnitude larger than specified by user
c       This may cause erroneous results if user sets 
c       PEPS=0.001 or larger.
c       
C##   PEPS=SFLOW*DMAX1(1.0D0,PMAX)       ! 
      PMAX=0.0D0                         !
      COUNTP=0                           !
C    !...................................!..................
      Icsym=1                            ! Right side OF BEARIN
      Kstep=1                            !..................
      Cseal=Cser                         !  exit end seal param.
      KLOSY=KLOSYr                       !  & entrance edge coeff.
      Dir=1.0D0                          !

  150 CONTINUE

      DO I=1, NXT+2                      ! .........................
          DVV(I)=zero                    ! Zero calculation vectors
          DV(I)= zero                    ! for TDMA algorithm
          DU(I)= zero                    !
          A(I)=  zero                    !
          B(I)=  zero                    !
          C(I)=  zero                    !
          D(I)=  zero                    !
      END DO                             !
                                         ! --------------------------------
      INTER=1                            ! =1 -> Inter-pockets region
c   

C    !-----------------------!...........! SYMMETRIC BEARING
      IF (ISYM.EQ.1) THEN
C    !-----------------------!...........! SYMMETRIC BEARING
        Kstart=Kstep         ! +/- 1

        IF ((NPOCKET.EQ.0).OR.(BEARING.EQ.3)) THEN ! => SYM JOURNAL BEARING
           DO J=1, NXT                
C             DVV(J)=0.0D0
              V1(J,0)=zero
              V1(J,1)=zero
           END DO
          KM1=Kstart      ! => start at Y=0 CV
          K=KM1
          DYPK=DYP(K)
          KV=K
          GOTO 650         !==> calculate solution on extended lands
        END IF

C    !-----------------------!...........! ASYMMETRIC BEARING
      ELSE
C    !-----------------------!...........! ASYMMETRIC BEARING
c     for 1st V-CV on centerline, y=0, and asymmetric BEARING
c                                        ! set initial counters
       IF (BEARING.EQ.1) THEN            !=> HJB
         K=Kstep                         !  1,-1
         KM1=-Kstep                      ! -1, 1
         KVM1=-2*Kstep                   ! -2, 
         KP1=2*Kstep                     !  2, -2
         Dyvk=2.0D0*DYV(Kstep)           !
       ELSE IF (BEARING.EQ.2) THEN       ! => SEAL
         KM1=Kstep                       ! 1
         K=KM1+Kstep                     ! 2
         KP1=K+Kstep                     ! 3
         KVM1=KM1                        ! 1
         Dyvk=DYV(K)
         IF (INERPDUM.EQ.1) THEN   !.....! Inertia at edge of seal
            INERP=INERPDUM
            INERL=1
         END IF                    !.....!.................
       ELSE IF (BEARING.EQ.3) THEN       !=> JOURNAL BEARING
         KM1=-Kstep                      ! -1, 1
         K=KM1-Kstep                     ! -2, 2
         KP1=K-Kstep                     ! -3, 3
         KVM1=-Kstep                     ! -1, 1
         Dyvk=DYV(K)                     !
         Dir=-Dir                        !  
       END IF                            !.......................
         SvnK=SVN(K)                     !
         Svsk=SVS(K)                     !
        
C       FIND AXIAL VELOCITY at Y=0 plane 
c       LEFTBC,RIGHTBC = 0 periodicity, 1: groove, 2: recess edge
c      !..................................!
        IF (NPOCKET.GT.0) THEN
c      !..................................! HJB BEARING=1
        LEFTBC=1+IFULL                       
        RIGHTBC=2
        DO I=1, NPOCKET+1-IFULL
           IF (I.GT.NPOCKET) RIGHTBC=1+IFULL
           CALL SETLIM(I,1)               ! INTER=1
           CALL CALCV1(FP, Res, Dir, Cseal, KVM1)
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
        CALL CALCV1(FP, Res, Dir, Cseal, KVM1)
        IF (IFULL.EQ.1) V1(NXT,K)=V1(1,K)    
        Kstart=-Kstep  				!## LSA 9/1/3/95
c      !..................................!
        END IF
c      !..................................!


c      !.............................!
        IF (BEARING.EQ.2) THEN       ! => SEAL
c      !.............................!
           Kstart=1
           DO J=1, NXT
             DVV(J)=DV(J)
           END DO
           KM1=Kstart
           K=KM1+Kstep
           DYPK=DYP(K)
           KV=K
           INERL=INERLDUM            ! Reset inertia at lands 
           GOTO 655                  !==> INTER=0 goto extended lands
c      !.............................!
        ELSE IF (BEARING.EQ.1) THEN  ! => HJB 
c      !.............................!
C##??? if npocket=0 then KM1,K,KV, DYPK = ???????
           DO J=1, NXT                      
              DVV(J)=DV(J)           ! Assure a unique V(y=0)
              V1(J,0)=0.5D0*(V1(J,1)+V1(J,-1))
              V1(J, 1) = V1(J,0)
              V1(J,-1) = V1(J,0)
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
              V1(J, Kstep)=V1(J,KV)*SVSK+V1(J,-KV)*SVNK ! at middle of first CV
            END DO

              KM1=-2*Kstep           ! -2, 2
              K=Kstep                !  1,-1
              Dir=-Dir
              KV=KM1
              GOTO 655    !==> INTER=0 GOTO extended lands

c      !.............................!
        END IF
c      !.............................!

C    !-----------------------!...........! FIRST volume 
      END IF                             ! of
C    !-----------------------!...........! ASYMMETRIC BEARING
      
c
C    !...................................!....................................
      IF (NPOCKET.EQ.0) GOTO 650         ! treat as PLAIN Bearing OR SEAL
C    !...................................!....................................
      
c    !........................! SET first-order recess-pressures
      CALL PJET               ! including angled-jet model
      CALL PJET1              ! 10/95
c    !........................!............................................!

c     ========................................................................
      INTER=1                            ! ON LANDs BETWEEN RECESSES OR GROOVES
c     ========...........................!------------------------------------

C    !...................................! SET boundary conditions
      Kend=Kstep*NPA                     ! FOR RECESS (POCKETS) BEARINGS
C    !...................................! SET boundary conditions at LEFT SIDE
      IF (IFULL.EQ.1) THEN               !
          LEFTBC=2                       ! RECESS EDGE NPOCKET
          PRECL1=PREC1do(NPOCKET)        !
      ELSE                          !....!........
          LEFTBC=1                       ! PAD LEADING EDGE
          PRECL1=ZERO                    ! P1GROOVE=0
      END IF                             ! i.e. P0groove = constant                                   
C    !...................................!....................................

C    !---------------!                   ! -----------------------------
  300 DO I=1, NPOCKET+1-IFULL            ! Sweep over inter-recess lands
C    !---------------!                   ! -----------------------------

          CALL SETLIM(I, INTER)          ! Sets limits for sweeps
C     ...................................! 
      IF (ISOTH.NE.1) THEN               !
          DO J=JPMIN,JPMAX               ! 
            P1S(J)=P1(J,Kstart)          ! for thermal soln.
          END DO                         !
      ENDIF                              !
C     ...................................! 
        IF (I.LE.NPOCKET) THEN           ! & Bond Cond. at RIGHT BDY
          RIGHTBC=2                      ! An I-th recess
          PRECR1=PREC1up(I)              ! (upstream edge)
        ELSE                             !......
          RIGHTBC=1                      ! RIGHT BDY: A PAD groove
          PRECR1=ZERO                    ! i.e. P0groove = constant 
        END IF                           !..........................

          K=0                            !
          KM1=Kstart                     ! 1 if ISYM=1, -1/+1 if ISYM=0
c       :::::::::::::::::::::::::::::::::! ...................................

  400     K=K+Kstep                      !
          KP1=K+Kstep                    !
          DYPK=DYP(K)                    ! ................................

          IF (K.EQ.Kend) THEN            ! BDYC's for U1 - equation
             U1w=(0.0D0,0.0D0)           ! 
             U1e=U1w                     ! at axial edge K=KEND
          ELSE                           !...
            U1w=U1(justart,k)            ! at edge K=1,..KEND-1
            U1e=U1(justop,k)             !
          END IF                         !.............

c                                        ! Solve U1 equation
c     !..................................! ...................................
          CALL CALCU1(FU,RES,Dir,C1U,C2U,U1e, U1w, K)        
c     !..................................! ...................................
          IF (K.EQ.Kend)  goto 500       !                 --> NEXT RECESS

          KM1=K                          !
          K=KM1+Kstep                    !
          KP1=K+Kstep                    !
          DYVK=DYV(K)                    ! Y size of Vcv
          SVNK=SVN(K)                    ! weight coeffs for V-eqn
          SVSK=SVS(K)                    ! 

c     !..................................! Solve V1 equation.................
          CALL CALCV1(FP, RES,Dir,Cseal,KM1)
c     !..................................! ...................................
          KP1=K                          !
          K=KM1                          ! ................................
          KM1=K-Kstep                    !
          IF(K.eq.Kstep) KM1=Kstart

c            !-- FOR JOURNAL rotations multiply by axial moment ARM
              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for journal rotations
                  YPO=YP(K)-YO
                  YVN=YV(KP1)-YO
                  YVS=YV(K)-YO
                  DO J=JVSTART,JVMAX     
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                    FPYN(J)=FP(J)*YVN
                    FPYS(J)=FP(J)*YVS
                  END DO
c            !--
              END IF                     
c            !--

c     !..................................! ...................................
c         ON K ROW: Correct P1,U1,V1(j,k),V1(j,k+1)
          CALL CALCP1(FUY,FPY,FPYN,FPYS,SIGMA,Dir,K) 
c     !..................................! ...................................
          IF (ISOTH.NE.1) THEN           ! Solution of 1st-order energy eq.
            CALL CALCT1(FP,RES,SIGMA,P1S,Dir,K)
          ENDIF                          !...................................

          KM1=K                          
          GOTO 400                       ! -> NEXT ROW

  500     LEFTBC=2
          PRECL1=PREC1do(I)                  
C    !---------------!                   ! Next land/recess ..............  
      END DO                             ! i=1,2,,, NPOCKET+1-IFULL
C    !---------------!                   !................................
C                                        ! Set velocity at y=0
C.!......................................! ASYM HJB: Set V=AVE(Vleft+Vright)
      IF (ISYM.eq.0) THEN     
         DO J=1, NXT
           V1(J, 0)=0.5D0*(V1(J,1)+V1(J,-1))
           V1(J, 1)=V1(J,0)
           V1(J,-1)=V1(J,0)
         END DO
      END IF!............................!.......................

C    ---------                           !...................................
  600 INTER=0                            ! At edge k=Npa+1 sweep on ext lands
C    ---------                           !
      CALL SETLIM(I, INTER)              ! RESets limits for sweeps
      IF (IFULL.EQ.1) THEN               ! C###
          JVMAX=NXT                      !    
          JVSTOP=2                       !
      END IF                             !...........................

      K=NPAP1*Kstep                      !
      KM1=NPA*Kstep                      !
      KP1=K+Kstep                        !
      DYVK=DYV(K)                        ! Y size of Vcv above recess edge
      SVNK=SVN(K)                        ! weight coefficients
      SVSK=SVS(K)                        ! .................................

      LEFTBC=1-IFULL                     ! SET BDY Type: 0= periodicity
      RIGHTBC=LEFTBC                     ! or 1=groove 
c    !...................................!.......................................!
      CALL CALCV1(FP,RES,Dir,Cseal,KM1)  ! Solve V1 equation at recess edge
c    !...................................!.......................................!
      IF (IFULL.EQ.1) THEN               !..................................
      V1(1,KM1)=V1(NXT,KM1)              ! inforce periodicity conditions
      P1(1,KM1)=P1(NXT,KM1)              ! for 360 deg. Bearing
      T1(1,KM1)=T1(NXT,KM1)              ! for 360 deg. Bearing
      V1(1,K)  =V1(NXT,K)                !
      END IF                             !.....

C     -----------------------------------!....................................... 
      INTER=1                            ! On lands between recess edges
C     -----------------------------------! K=+/- NPA ........................... 
      K=NPA*Kstep                        ! On inter-recess lands at recess edge
      KP1=K+Kstep                        ! .................................
      KM1=K-Kstep                        !
      DYPK=DYP(K)                        ! Y size of Pcv

      YPO=YP(K)-YO                       !.......
      YVN=YV(KP1)-YO                     ! axial moment arm for rotations
      YVS=YV(K)-YO                       !.......

      IF (IFULL.EQ.1) THEN               ! Pressure at left boundary
        LEFTBC=2                         ! a recess:
        PRECL1=PREC1do(NPOCKET)          !
        PRECL =PRECdo(NPOCKET)           !
        TRECL =TREC(NPOCKET)             !
      ELSE                               !...
        LEFTBC=1                         ! a pad groove
        PRECL1=ZERO                      !
      END IF                             !

c     ...................................!
      DO I=1, NPOCKET-IFULL+1            ! at edge of recess K=NPA
c     ...................................!
          CALL SETLIM(I, INTER)           

        IF (I.LE.NPOCKET) THEN
          RIGHTBC=2
          PRECR1=PREC1(I)
          PRECR =PREC(I)
          TRECR =TREC(I)
        ELSE
          RIGHTBC=1
          PRECR1=ZERO
        END IF

              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for rotations
                  DO J=JVSTART,JVMAX     !............
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                    FPYN(J)=FP(J)*YVN
                    FPYS(J)=FP(J)*YVS
                  END DO
              END IF                     !............
c                                        ! CORRECT at edge of inter-rec lands:
          CALL CALCP1(FUY,FPY,FPYN,FPYS,SIGMA,Dir,K) 
          IF (ISOTH.NE.1) THEN           ! Solution of 1st-order energy eq.
            CALL CALCT1(FP,RES,SIGMA,P1S,Dir,K)
          ENDIF                          !
C                                        ! U1,P1,T1,V1(j,npa)
        LEFTBC=2
        PRECL1=PREC1do(I)
        PRECL =PRECdo(I)                  
        TRECL =TRECR                  
c     ...................................!
      END DO                             ! I=1, NPOCKET-1+IFULL
c     ...................................!................................

      IF (IFULL.EQ.1) THEN               ! for 360 deg HJBearing
      V1(NXT,K)=V1(1,K)                  ! inforce periodicity
      U1(NXT,K)=U1(1,K)                  ! k=NPAxKstep
      P1(NXT,K)=P1(1,K)                  !
      T1(NXT,K)=T1(1,K)                  !
      V1(NXT,KP1)=V1(1,KP1)              !
      END IF        !....................!.....

      DO J=1, NXT                        ! Save pressure coeff. for extended
          DVV(J)=DV(J)                   !  land correction
      END DO                             !..........................

      C1U=0.25D0                         ! ......... coefficients for
      C2U=0.75D0                         !           Us interpoltion AT K=NPA+1

      KM1=Kstep*NPA                      !
      K=KM1+Kstep                        ! resets axial counters
      DYPK=DYP(K)                        !
      KV=K                               !

C   !------------------------------------!-----------------------------------!
  650 INTER=0                            ! On extended lands above recesses
C   !------------------------------------!-----------------------------------!
      LEFTBC=1-IFULL                     ! BDY COND: 0: continous 360 land
      RIGHTBC=LEFTBC                     ! 1: specified P at groove
      CALL SETLIM(I,INTER)               ! RESET circumferential counters

C----!      
  655 Kend=Kstep*NYI
C----!

C        !...............................! for thermal analysis
      IF (ISOTH.NE.1) THEN               !
      DO J=JPMIN,JPMAX                   ! 
        P1S(J)=P1(J,K  )*SVS(K)+P1(J,KM1)*SVN(K)           
      END DO                             !
      ENDIF                              !
C        !...............................!

  700 CONTINUE
      U1e=U1(justop ,K)    ! U1 - BC's
      U1w=U1(justart,K)
c    !.........................!         ! 
      IF ( K.eq.Kend ) THEN              ! AT Exit plane of Bearing
c    !.........................!         !

C        FOR TWO recess ROW BEARING (right side)
c       !-----------------------------------------------!
         IF ((Icsym.eq.1).AND.(MODEL.eq.1) ) THEN  
            KP1=NYIP1
         ELSE  
            KP1=Kend                    
         END IF

c         !--------------------------!     
         CALL CALCU1(FU, RES, Dir, C1U, C2U, U1e,U1w, K )
c         !--------------------------!     

         GOTO 750

c    !.........................!         
      ELSE                     ! |K| < |KEND|
c    !.........................!  AT internal control volumes on extended lands

         KP1=K+Kstep
c          !--------------------------!     
         CALL CALCU1(FU, RES, Dir, C1U, C2U, U1e, U1w, KV )        
c          !--------------------------!     

c     .....................!K=NPa+1, ... Kend-1
      END IF
c     .....................!

      C1U=0.5D0
      C2U=C1U
      KM1U = KM1
 
      KM1=K                              !
      K=KM1+Kstep                        !
      DYVK=DYV(K)                        ! Y size of Vcv
      SVNK=SVN(K)                        ! weight coeffs.
      SVSK=SVS(K)                        !

c     .....................!.......................!
      IF (K.eq.Kend) THEN
c     .....................!
         IF ((Icsym.eq.1).AND.(MODEL.eq.1) ) THEN  
            KP1=NYIP1
         ELSE  
            KP1=Kend                    
         END IF
c     .....................!
      ELSE
c     .....................!
         KP1=K+Kstep
c     .....................!
      END IF
c     .....................!

c       !--------------------------! Solve V1 equation
      CALL CALCV1(FP,RES,Dir,Cseal, KM1)         
c       !--------------------------!

      IF (IFULL.EQ.1) V1(JVSTOP, K)=V1(JVMIN, K)  ! inforce periodicity

c     ...................................!. Correct: UVP fields..!
      KP1=K                              !
      K=KM1                              ! 
      KM1=KM1U                           ! K-Kstep

c            !..                               for journal rotations
              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for rotations
                  YPO=YP(K)-YO
                  YVN=YV(KP1)-YO                  YVS=YV(K)-YO
                  DO J=JVMIN,JVSTOP      !............
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                    FPYN(J)=FP(J)*YVN
                    FPYS(J)=FP(J)*YVS
                  END DO
c            !..                               for rotations
              END IF                     
c            !..                               for rotations

c                                  ! Correct P1(j,k),U1(j,k),V1(j,kv),V1(j,kp1)
c       !--------------------------!     
      CALL CALCP1(FUY,FPY,FPYN,FPYS,SIGMA,Dir,KV)     
c       !--------------------------!     !.................................!
      IF (ISOTH.NE.1) THEN               !  Solution of 1st-order energy eq.
        CALL CALCT1(FP,RES,SIGMA,P1S,Dir,KV)    
      ENDIF                              !
c       !--------------------------!     !.................................!

      KM1=K                !===> ready for sweep next K-row
      K=K+Kstep 
      KV=K
      DYPK=DYP(K)
c##   C1V=0.5D0
c##   C2V=C1V
c  !::::::::::::::::::::::!              !

      GOTO 700                           ! => NEXT K axial CV 
C ...................................... ! ...................................

C :::::::::::::::::::::::::::::::::::::: !::::::::::::::::::::::::::::::::::
C AFTER solving land flow equations:
c     ...................................!....................!
  750 CONTINUE                           ! for 2row HJB: Solve for P on SYM Line

c     ...................................!....................!
      IF ( (Icsym.eq.1).AND.(MODEL.eq.1) ) THEN  
c    ....................................! MODEL= (1) = 2row HJB 
            DO J=1, NXT+2
             DV(J)=zero                  ! To set ANP=0
            END DO

              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for journal rotations
                  YPO=YP(K)-YO
                  YVN=YV(KP1)-YO
                  YVS=YV(K)-YO
                  DO J=JVMIN,JVSTOP      !............
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                    FPYN(J)=FP(J)*YVN
                    FPYS(J)=FP(J)*YVS
                  END DO
              END IF                     !............

c                                        ! Correct P,U,V at SYMMETRY side=RIGHT
        CALL CALCP1(FUY,FPY,FPYN,FPYS,SIGMA,Dir,K)   

CTTT        IF (ISOTH.NE.1) THEN             ! Solution of 1st-order energy eq.
CTTT          CALL CALCT1(FP,RES,SIGMA,P1S,Dir,K)  
CTTT        ENDIF                            !
c     ...................................!....................!
            DO J=1, NXT
              V1(J,NYIP1)=zero
              P1(J,NYIP1)=zero
              T1(J,NYIP1)=zero           ! 
            END DO
CTTT            GOTO 800
c    ....................................! 
      END IF
c    ....................................! MODEL= 1 = 2row HJB 

C ...................................... ! MODEL= 2 = 1 row HJB
      IF (ISOTH.NE.1) THEN               !
        KEM1=Kend-Kstep                  !
        KEM2=KEM1-Kstep                  !
        KEM3=KEM2-Kstep                  !
        DO J=1,NXT                       ! Exit temperatures
           T1(J,Kend)=(3.D0*T1(J,KEM1)-  !
     +                 T1(J,KEM2))/2.D0  ! From Roach's Book P189-191 
        END DO                           !
      ENDIF                              !
C ...................................... ! ----------------------------

  800 INERP=INERPDUM
      Kend=Kstep*NPA

c    !...................................!
      IF ((INERP.EQ.0).OR.(NPOCKET.EQ.0)) goto 850 
c    !...................................!
                         

      LEFTBC=1+IFULL
      IF (IFULL.EQ.1) THEN               ! specify pressure at boundaries
        PRECL1=PREC1do(NPOCKET)          ! LEFTBC=2 a pocket at left edge
        ASPEC=ASPE(J)                    !
      ELSE                               !...............................
         PRECL1=ZERO                     ! LEFTBC=1 a groove at left edge
      END IF                             !...............................

c   !....................................!................................! 
      DO I=1, NPOCKET+1-IFULL            ! Correct Uinlet & Pedge lateral
c   !....................................!................................! 
         CALL SETLIM(I, 1)               ! Precl1 | ............. | Precr1

          IF (I.LE.NPOCKET) THEN         ! SET BDY conditions
             RIGHTBC=2                   ! =2 a recess edge at right edge
             PRECR1=PREC1up(I)           !
          ELSE                           !..
             RIGHTBC=1                   ! =1 a groove at right edge
             PRECR1=ZERO                 !
          END IF                         !.................................!

c        !...............................!
          DO K=Kstep,Kend,Kstep          ! AXIAL sweep
c        !...............................!

c            !..                               for journal rotations
              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for rotations
                  YPO=YP(K)-YO
                  DO J=JVSTART,JVMAX     
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                  END DO
              END IF                     !............

c            !-------------------!       !.............................
              CALL EDGEU1(FUY,FPY,SIGMA,Dir,K,INERP,Kend)
c            !-------------------!       !.............................

c        !...............................!
          END DO                         !sweep K=...-> +/- (NPA-1)
c        !...............................!

          LEFTBC=2                       ! LEFT BDY: A recess edge
          PRECL1=PREC1do(I)              !
          IF (I.LE.NPOCKET) ASPEC=ASPE(I)
c    !...................................! 
      END DO                             ! i=1,...NPOCKET-1+IFULL
c    !...................................! 


 850  CONTINUE

c    !...................................!
      IF ((IFULL.EQ.1).AND.(NPOCKET.GT.0)) THEN ! => HJB 
c    !...................................!
      DO K=Kstep,Kend,Kstep              ! .....................................
          V1(NXT, K)=V1(1, K)            ! Inforce periodicity at right bound.
          U1(NXT, K)=U1(1, K)            !
          P1(NXT, K)=P1(1, K)            !
          T1(NXT, K)=T1(1, K)            !
      END DO                             ! ............
c    !...................................!
      ELSE IF (IFULL.EQ.0) THEN          ! Advect downstream last U vel
c    !...................................!
c         DO K=Kstep, Kend, Kstep         !##      
c            U1(NXT,K)=U1(NXT-1,K)        !NXT-1=JUMAX
c          END DO                         !##
c          P1(1,Kend)=ZERO                ! ????   
c    !...................................!
      END IF
c   !....................................!................................!



 860  CONTINUE      ! ==> GET READY to modify entrance AXIAL (edge) pressures

c    !.....................................!
      IF (BEARING.EQ.1) THEN               ! => HJB RECESS EDGES
c    !.....................................!
        SVNK=SVN(Kstep*NPAP1)              
        SVSK=SVS(Kstep*NPAP1)              
        DYPK=DYP(Kstep*NPA)                
c    !.....................................!
      ELSE IF (BEARING.EQ.2) THEN          ! => SEAL
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


c   !....................................!................................!
 870  CONTINUE
 
              IF ((IFM.EQ.1).AND.(ISOLN.EQ.1)) THEN        
c            !..                               for journal rotations
                  YPO=YP(Kend)-YO
                  YVN=YV(NPAP1*Kstep)-YO
                  DO J=1, NXT            
                    FUY(J)=FU(J)*YPO
                    FPY(J)=FP(J)*YPO
                    FPYN(J)=FP(J)*YVN
                  END DO
              END IF                    

c   !....................................!  Correct Vinlet & Pedge axial 
      CALL EDGEV1(FUY,FPY,FPYN,SIGMA,Dir,DYPK,
     +            KLOSY,SVNK,SVSK,Kstep) 
c   !....................................!................................!

c                                        ! Modify P1(NYI) for end seal cond.
 880  IF (Cseal.gt.(0.0D0)) THEN
        CALL ENDSEAL1(Dir,DYP(Kstep*NYI),Cseal,alfp,Kstep)
      END IF
c   !....................................!................................!


C   !.......................!............!......................!
c     Direct for symmetry/asymmetry on bearing
C   !.......................!............!......................!

c     ...................................!....................!
      IF (ISYM.eq.0) THEN
c    !..........................!     ! asymmetric BEARING
         Factor=1.0D0

         IF (BEARING.EQ.2) GOTO 900   !=> SEAL
         IF (BEARING.EQ.3) THEN       !=> BEARING at Y=0

            KV=2
            DYPK=DYP(1)+DYP(-1)
            SVSK=DYP(-1)/DYPK
            SVNK=DYP( 1)/DYPK  
            DO j=1, NXT
               P1(j,-Kstep)=P1(j,Kstep)      ! SET Unique fields
               T1(j,-Kstep)=T1(j,Kstep)      ! SET Unique fields
               U1(j,-Kstep)=U1(j,Kstep)      ! at Y=0
               V1(j, Kstep)=V1(j,KV)*SVSK+V1(J,-KV)*SVNK 
               V1(j,-Kstep)=V1(j, Kstep)
            END DO
         END IF

         IF (Icsym.eq.2) GOTO 890

         Dir=-1.0D0
         Kstep=-1
         Cseal=Csel                      ! Discharge pressures & end seal
         KLOSY=KLOSYl                    ! & entrance edge coeff.
         Icsym=2                         !
         GOTO 150                        ! --> SOlve Left side of bearing

c    !...........................!
      ELSE
c    !...........................! FOR SYMMETRIC BEARING
         Factor=2.0D0
         GOTO 900
c    !...........................!
      END IF
c    !...........................!


c    !..................................! ASYMMETRIC BEARING=3
 890  IF (BEARING.EQ.3) THEN
c    !..................................! 
        DO j=1, NXT
          P1(j, 0)=P1(j,1)
          T1(j, 0)=T1(j,1)
          U1(j, 0)=U1(j,1)
          V1(j, 0)=V1(j,1)
        END DO
        GOTO 900
c    !..................................! 
      END IF
c    !..................................! 

c   !....................................! HYDROSTATIC BEARING AT center line:
 892    DO j=1, NXT
        P1(J,0)=0.5D0*(P1(J,1)+P1(J,-1))
        T1(J,0)=0.5D0*(T1(J,1)+T1(J,-1))
        U1(J,0)=0.5D0*(U1(J,1)+U1(J,-1))
        V1(J,0)=0.5D0*(V1(J,1)+V1(J,-1))
      END DO

c    !...................................!................................

  900 CONTINUE

C........................................!
      IF (ISOTH.EQ.1) THEN               ! ISOTH=1: Isothermal cases
         CALL FLOWS1(QMAX)               ! Calculates flows on recesses
      ELSE
         CALL FLOWT1(QMAX,SIGMA,FP,FU)   ! Flows and energy balance on recesses
      ENDIF
C    !...................................!


 1120 CONTINUE


C    !...................................!................................
C                                        ! Print convergence steps
      IF (SVERB.EQ.1) THEN               !
          WRITE (6, 1400) ITER, SMASS, QMAX,
     +                    PMAX,PMAXP,Countp,MAXCOUNTP
      END IF                             !

      IF ((DEVICE.EQ.1).AND.(DVERB.EQ.1)) THEN 
          WRITE (1, 1400) ITER, SMASS, QMAX,
     +                    PMAX,PMAXP,Countp,MAXCOUNTP
      END IF                             !

C    !...................................!............................

C    !...................................!......................
C     CHECKS CONVERGENCE FOR BEARINGS
C    !...................................!......................
      IF(ITER.EQ.1) THEN                 !
            SOLD=SMASS                   !
            GOTO 100                     ! 2 its minimum
C    !...................................!......................
      ELSE IF (ITER.GE.ITMAX) THEN       !
            GOTO 1200                    ! => NO CONVERGENCE
C    !...................................!......................
C##   ELSE IF (ITER.EQ.10) THEN          !
C##         SMASSc=SMASS                 ! save mass error
C    !...................................!.......................
C??## ELSE IF (ITER.EQ.ITstop) THEN      !
C##        IF (DABS(SMASS/SMASSc).GT.10) GOTO 1600
C###            SMASSc=SMASS                 !
C###            ITstop=ITstop+ITdelta        !=> divergence occurs.           
C    !...................................!......................
      END IF
C    !...................................!......................


c    !...................................!................................
      IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN
c    !...................................!................................
        SMASS=Factor*SMASS/QMAX          !=> for HJB
        IF (SMASS.LE.MPEPS) GOTO 1300    !-> Converged
c    !...................................!................................
      ELSE IF (BEARING.GE.2) THEN
c    !...................................!................................
      IF (COUNTP.GE.Maxcountp) GOTO 1300 ! Convergence on Pressure
c    !...................................!................................
      END IF
C    !...................................!................................

C    !...................................!................................
      IF (Iconv.GE.3) THEN               ! ==> error grows
        IF((SMASS-SOLD).GT.0.0D0) GOTO 1570   
      ELSE                               !
        Iconv=Iconv+1                    !..........
      END IF                             !
C    !...................................!................................

      SOLD=SMASS			 ! save last values
C ...........................................................................

C   !==========!
 977  GOTO 100                           ! -> NEW ITERATION
C   !==========!


C ...........................................................................
 1200 CALL BEEPER                        !
C      WRITE (6, *) 'MAXIMUM NUMBER OF ITERs EXCEEDED, SLAND1 ABORTS'
C      write (6, *) 'DYNAMIC force coefficients may be in error'

      IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN        
      WRITE (1, *) 'MAXIMUM NUMBER OF ITERs EXCEEDED, SLAND1 ABORTS'
      write (1, *) 'DYNAMIC force coefficients may be in error'
      END IF                             !
      goto 1700                          !

 1577 CALL BEEPER                        !............................
      WRITE (6, 1578) iter               ! Error grows
      IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN  
      WRITE (1, 1578) iter               ! 
      END IF                             !
      GOTO 1700                          !..........................
                                         !
 1300 IF ((DEVICE.eq.1).AND.(DVERB.ge.1)) THEN
          WRITE (1, 1310) ITER           ! CONVERGENCE
      END IF                             !
      GOTO 1700 
C    !...................................!......................

 1570 Iconv=0                            !..........................!
      IF (Alfp.LE.0.25D0) THEN           ! Reduce relaxation params.
        Alfp=0.25D0                      ! Min values are
      ELSE                               ! alfp=alfu=0.25
        Alfp=Alfp**1.10                  !..........................!
        Alft=Alft**1.10
        ITER=ITER-1
      END IF
      GOTO 100                           !==> START OVER
c........................................!..........................

 1600 CONTINUE
      write (6,*) '$ Source mass error has increased 10 fold on last'
      write (6,*) '$ ', ITdelta, ' its, ==> PROGRAM aborts'
 
 
c........................................!..........................
 1700 alfv=relaxv                        ! re-store relax pars.
      IF (alfp.LE.0.25D0) THEN
        write(6,*) 'Relax parameter AP too low = 0.25,',
     +             '==> RESULTS may be in error'
      END IF
      alfp=relaxp                        !
      alft=relaxt
      betv=1.0D0-alfv                    ! 
      PEPS=SFLOW                         !
c........................................!...........................

c........................................!...........................
 1578 FORMAT (' ', 'Error started to grow at it:', I3,
     + 3X, 'SLAND1 aborts', /, 3X, 
     + 'DYNAMIC force coefficients ARE in ERROR')

 1310 FORMAT (' ', 'Convergence achieved in', I3, ' iterations.')

 1400 FORMAT (' ', 'It:', I3, 1X, 'SUM(Ms):', E10.3E2, 1X,
     +        ' QMAX:',E10.3E2, ' PMAX:',
     +        E10.3E2,' Ppmax:',E10.3E2,' #Ps:',I3,'=>',I3)

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine MODGS                                                       **
C **                                                                         **
C ** MODGS: Modifies GUO and GVO arrays for moment/rotation coefficients     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE MODGS(YO)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT), 
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
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

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
c     .................................................................
      DOUBLE PRECISION XP, XU, YP, YV,
     +                 APU, AWU, AEU, ASU, ANU, GUO, GUV, APUI,
     +                 APV, AWV, AEV, ASV, ANV, GVO, GVU, APVI,
     +                 GUPR, GUPI, GUT, GVPR, GVPI, GVT
      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        ISYM, ICSTEP, IFULL, BEARING


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION YO,Ycu,Ycv
      INTEGER J,K,Kstart,Kend

C ----------------------------------------------------------------------------
      Kend=NYI

      IF ((ISYM.eq.1).AND.(YO.eq.0.0D0) ) THEN
         Kstart=1
      ELSE
         Kstart=-NYI
      END IF

      IF (BEARING.EQ.2) Kstart=1

      DO K=Kstart,Kend
         Ycu=YP(K)-YO
         Ycv=YV(K)-YO
         DO J=1, NXT
           GUO(J,K)=Ycu*GUO(J,K)
           GVO(J,K)=Ycv*GVO(J,K)
         END DO
      END DO

      RETURN
      
      END

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Last revised 3/8/94 by Luis San Andres
C Modifications for TWO-CELL model on 6/20/95 by Luis San Andres
c revised for angled jet injection on 9/13/95 by LSA
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::