c -------------------------------------------------------
c 10/4/95 modified URECedge on EDGEV1 replaces SWIRL
c         on EDGEU1  removed SWIRL !!
c -------------------------------------------------------
c 9/15/95 modified for angled injection on 9/15/95
c check for SWIRL on EDGEV1, CALCV1
c -------------------------------------------------------
c 6/20/95 includes modifications for CELL depth compliance
c         effects
c -------------------------------------------------------
c INCLUDES modifications for BEARING COMPLIANCE EFFECTS
c -------------------------------------------------------
c                                        #####
c  ####    ####   #    #  #####   ###### #     #   #####          ######
c #    #  #    #  ##  ##  #    #  #            #     #            #
c #       #    #  # ## #  #    #  #####   #####      #            #####
c #       #    #  #    #  #####   #      #           #     ###    #
c #    #  #    #  #    #  #       #      #           #     ###    #
c  ####    ####   #    #  #       ###### #######     #     ###    #
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
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcu1                                                      **
C **                                                                         **
C **  CALCU1:  Finds U1 velocity from X perturbed momentum equation.         **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCU1(FPU, RES, Dir,C11,C22,Ue,Uw,KV)

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
     + GUPR(MAXNXT,-MAXNYI:MAXNYI),GUPI(MAXNXT,-MAXNYI:MAXNYI),
     + GUT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /DPUV1/ BUC(MAXNXTP2), BVC(MAXNXTP2), BVVC(MAXNXTP2)
      COMMON /T1ARRAY/ T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /TDMA1/ A1(MAXNXTP2), B1(MAXNXTP2),
     +               C1(MAXNXTP2), D1(MAXNXTP2)
      COMMON /PRECES1/ PRECL1, PRECR1

      COMMON /KVALB/ K, KM1, KP1
      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /LRBOUND/ LEFTBC,RIGHTBC
      COMMON /SOLN/ ISOLN

      DOUBLE PRECISION APU, AWU, AEU, ASU, ANU, GUO, GUV, APUI,
     +                 GUPR,GUPI, GUT, DXP, DXU, HP, HU,HV, SUW, SUE,
     +                 Drhop,Demup,Drhot,Demut, DYVK, DYPK, SVNK, SVSK,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      DOUBLE COMPLEX U1,V1,P1,T1,BUC,BVC,BVVC,A1,B1,C1,D1,PRECL1,PRECR1

      INTEGER K, KM1, KP1, ISOLN,
     +        JUMIN,JUMAX,JUSTART,JUSTOP,JPMIN,JPMAX,JPSTART,JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT
      INTEGER  LEFTBC, RIGHTBC, IFULL
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION RES, factor, FPU(*) ,
     +                 Suej, SUwj, C11, C22, DPU, SANB, Dir
      DOUBLE COMPLEX PP1, TP1, SO1U, VP1, US, UN, Uw, Ue
      INTEGER J, JJ, JJJ, JJMAX, JM1, JP1, KV

C ----------------------------------------------------------------------------
C --  CALCU1 code                                                           --
C ----------------------------------------------------------------------------
      JJJ=2-JUMIN                              !
      JJMAX=JUMAX+JJJ                          !

C    !-----------------!                       ! on STANDARD U equation
      DO J=JUMIN, JUMAX                        ! ::::::::::::::::::::::::::::
C    !-----------------!                       !
          JP1=J+1                              ! Sweep from left -> right
          JJ=J+JJJ                             ! J+2-Jumin
          Suej=SUE(J)                          !
          Suwj=SUW(J)                          !

          VP1=0.5D0*(V1(J,KV)+V1(J,KP1))*Suej  !
          VP1=VP1+0.5D0*(V1(JP1,KV)+V1(JP1,KP1))*Suwj ! V1 at center of Ucv

          PP1=Suej*P1(j,k)+Suwj*P1(jp1,k)      ! P1 at center CV
          TP1=Suej*T1(j,k)+Suwj*T1(jp1,k)      ! T1 at center CV

          SO1U=GUV(J,K)*VP1 +                  ! Source terms
     +         DCMPLX(GUPR(J,K),GUPI(J,K))*PP1+
     +         GUT(J,K)*TP1                    !

          IF (ISOLN.EQ.1) THEN                 ! -> Only for zero component
              SO1U=SO1U+GUO(J,K)*FPU(J)        !    solution
          END IF                               !

          DPU=HU(j,k)*DYPK                     ! ...........................
          D1(JJ)=SO1U+DPU*(P1(J,K)-P1(JP1,K))  ! Source + Pressure force
          C1(JJ)=-AWU(J, K)                    ! coef. of U1w
          A1(JJ)=-AEU(J, K)                    ! coef. of U1e

          B1(JJ)=DCMPLX(APU(J,K),APUI(J,K))    ! coef. of U1p

          US=C11*U1(J,KM1)+C22*U1(J,K)         ! South U veloc.

          UN=U1(J,KP1)

          D1(JJ)=D1(JJ)+ASU(J,K)*US+ANU(J,K)*UN

          D1(JJ)=D1(JJ)+BETU*B1(JJ)*U1(J, K)   !

          BUC(J)=DPU/B1(JJ)                    ! for SIMPLE procedure
      END DO                                   !
C -----------------------                      !......................!

      IF (INTER.EQ.0) GOTO 100    !=> for lands w/o recess

C    !............................................! correct edge pressures
      IF ((INERP.EQ.1).AND.(IABS(K).LT.NPA)) THEN
C    !............................................! correct edge pressures
        IF (RIGHTBC.EQ.2) THEN
          D1(JJMAX)=D1(JJMAX)+HU(JUMAX,k)*DYPK*(P1(JUSTOP,K)-PRECR1)
        END IF
        IF (LEFTBC.EQ.2) THEN
          D1(2)=D1(2)+HU(JUMIN,k)*DYPK*(PRECL1-P1(JUMIN, K))
        END IF
C    !............................................! correct edge pressures
      END IF
C    !............................................!

C ............................................ ! ............ Solve by .......
  100 CALL TDMACPLX(Uw,Ue, JJMAX)              ! TDMA Algorithm
C ............................................ ! .............................
  200 DO J=JUMIN, JUMAX                        ! B(jj) = Usol(j)
          JJ=J+JJJ                             !
          U1(J, K)=B1(JJ)                      ! Order solution
      END DO                                   ! .............................

      IF ((INTER.EQ.0).AND.(IFULL.EQ.1)) THEN  ! SET periodicity condition
          U1(JUSTOP, K)=U1(JUMIN, K)           ! on extended lands
      END IF                                   ! Justop=Nxt; Jumin=1
                                               ! -----------------------------
      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcv1                                                      **
C **                                                                         **
C **  CALCV1:  Finds V1 velocity for first order solution.                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCV1(FPV, RES, Dir, Cseal, KVM1 )

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /VCOEF/
     + APV(MAXNXT,-MAXNYI:MAXNYI),AWV(MAXNXT,-MAXNYI:MAXNYI),
     + AEV(MAXNXT,-MAXNYI:MAXNYI),ASV(MAXNXT,-MAXNYI:MAXNYI),
     + ANV(MAXNXT,-MAXNYI:MAXNYI),GVO(MAXNXT,-MAXNYI:MAXNYI),
     + GVU(MAXNXT,-MAXNYI:MAXNYI),APVI(MAXNXT,-MAXNYI:MAXNYI),
     + GVPR(MAXNXT,-MAXNYI:MAXNYI),GVPI(MAXNXT,-MAXNYI:MAXNYI),
     + GVT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFV,BETV,ALFP
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/ T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /DPUV1/ BUC(MAXNXTP2), BVC(MAXNXTP2), BVVC(MAXNXTP2)
      COMMON /TDMA1/ A1(MAXNXTP2), B1(MAXNXTP2),
     +               C1(MAXNXTP2), D1(MAXNXTP2)
      COMMON /RECJET1/ PREC1do(MAXNPOCK),PREC1up(MAXNPOCK),
     +                 PR1jet(MAXNPOCK,MAXNPOCK+2)

      COMMON /KVALB/ K, KM1, KP1
      COMMON /VVEC/ JVMIN, JVMAX, JVSTART, JVSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN
      COMMON /LRBOUND/ LEFTBC,RIGHTBC
      COMMON /BTYPE/ BEARING

c     ................................................................
      DOUBLE PRECISION APV, AWV, AEV, ASV, ANV, GVO, GVU, APVI,
     +                 GVPR, GVPI, GVT, DXP, DXU, HP, HU,HV, SUW, SUE,
     +                 Drhop,Demup,Drhot,Demut,
     +                 REP, REY, MSPEED, SPEED, SWIRL, PCAV,
     +                 ALFV, BETV, ALFP, DYVK, DYPK, SVNK, SVSK

      DOUBLE COMPLEX U1, V1, P1,T1,BUC, BVC, BVVC,
     +               A1, B1, C1,D1,PREC1do, PREC1up, PR1jet

      INTEGER K, KM1, KP1, ISOLN,
     +        JVMIN, JVMAX, JVSTART, JVSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,BEARING,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT
      INTEGER LEFTBC, RIGHTBC, IFULL,II
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION RES,factor, Dpv, Sanb, Dir, Cseal, FPV(*)
      DOUBLE COMPLEX PP1, TP1,UE1, UW1, UP1, SO1V, V1E, V1W
      INTEGER J, JJ, JJJ, JJMAX, JP1, I, IJ, JL, JR, KVM1

C ----------------------------------------------------------------------------
C --  CALCV1 code                                                           --
C ----------------------------------------------------------------------------
      IF (LEFTBC.EQ.0) THEN                    !SETS BC FOR V1 VELOCITY
          V1W=V1(JVSTART,K)                    ! BC=0, 360 deg pad
      ELSE                                     ! BC=1,2, spec P, V=0 normal
          V1W=(0.0D0,0.0D0)
      END IF

      IF (RIGHTBC.EQ.0) THEN
          V1E=V1(JVSTOP,K)
      ELSE
          V1E=(0.0D0,0.0D0)
      END IF

      JJJ=-JVMIN+2                             !
      JJMAX=JVMAX+JJJ                          !
C                                              ! .........................
      UW1=U1(JVSTART,K)*SVSK+U1(JVSTART,KM1)*SVNK ! <--- West U1 velocity


C                                              ! :::::::::::::::::::::::::
C   !------------------!                       ! On
      DO J=JVMIN, JVMAX                        ! V-CVs
C   !------------------!                       ! Sweep from left -> right
          JJ=J+JJJ                             !
          UE1=U1(J, K)*SVSK+U1(J, KM1)*SVNK    ! East U1 velocity
          UP1=(UE1+UW1)/2.0D0                  ! U1 velocity at center Vcv

          PP1=Svsk*P1(j,k)+Svnk*P1(j,km1)      ! P1 at center CV
          TP1=Svsk*T1(j,k)+Svnk*T1(j,km1)      ! T1 at center CV

          SO1V=GVU(J,K)*UP1  +                 ! Source term:
     +         DCMPLX(GVPR(J,K),GVPI(J,K))*PP1+
     +         GVT(J,K)*TP1                    !

          IF (ISOLN.EQ.1) THEN                 ! -> only for zero component
              SO1V=SO1V+GVO(J, K)*FPV(J)       !
          END IF                               !
C
          DPV=HV(J,K)*DXP(J)                   ! Area of pressure action

          D1(JJ)=SO1V+Dir*DPV*(P1(J,KM1)-P1(J,K))  ! Source + Pressure force

          C1(JJ)=-AWV(J, K)                    ! coef. of V1w
          A1(JJ)=-AEV(J, K)                    ! coef. of V1e

          B1(JJ)=DCMPLX(APV(J,K),APVI(J,K))    ! coef. of V1p

          D1(JJ)=D1(JJ)+ASV(J,K)*V1(J,KVM1)+ANV(J,K)*V1(J,KP1)

          D1(JJ)=D1(JJ)+BETV*B1(JJ)*V1(J,K)    !

          BVC(J)=DPV/B1(JJ)                    ! for SIMPLE procedure
          UW1=UE1                              !

      END DO                                   ! .........................
C   !----------------------------------------- ! END OF SWEEP
c
C    !................................................!
 200  IF ((INERP.EQ.0).OR.(BEARING.EQ.3)) THEN
C    !................................................!
C    PLAIN BEARINGS:
          GOTO 300

C    !................................................!
      ELSE
     +IF ((IABS(K).EQ.NPAP1).AND.(BEARING.EQ.1).AND.
     +    (NPOCKET.GT.0) ) THEN
c   !..........................................! Correct for edge pressures
C   HYDROSTATIC BEARING:                       !............................

c       !..................!                   !............................!
          DO I=1, NPOCKET                      ! SWEEP across recesses
c       !..................!                   !............................!
              II=1                             !
              JL=(I-1)*NXI+NLC                 ! from left corner of recess
              JR=I*NXI+1                       ! to right corner of recess
c           !. . . . . . . .!                  !
              DO J=JL, JR                      !
c           !. . . . . . . .!                  !
               JJ=J+JJJ                        !
               II=II+1                         !
               D1(JJ)=D1(JJ)+Dir*DXP(J)*HV(J,K)*(PR1jet(I,II)-P1(J,KM1))
c           !. . . . . . . .!                  !
              END DO                           !
c           !. . . . . . . .!                  !
c       !..................!                   !............................!
          END DO                               !
c       !..................!                   !............................!
          D1(2)=D1(JJMAX)                      !

C   !.............................................!
      ELSE IF ((BEARING.EQ.2).AND.(K.EQ.2)) THEN  !Modify eqns. at SEAL
C   !.............................................!inlet
C    ANNULAR SEALS:

         DO J=JVMIN, JVMAX
              JJ=J+JJJ
              D1(JJ)=D1(JJ)+DXP(J)*HV(J,K)*(-P1(J,1))
         END DO                                !

C    !................................................!
      END IF
C    !................................................!

c ................................................!............................!
 300  CONTINUE

c##      IF ((Cseal.gt.(0.0D0)).AND.(IABS(K).eq.NYI)) THEN
c##  !................................................!
c##      DO J=JVMIN, JVMAX
c##         JJ=J+JJJ
c##         DPV=HV(J,K)*DXP(J)
c##         D1(JJ)=D1(JJ)+HV(J,K)*DXP(J)*Dir*P1(J,K)    ! P1(J,exit)=(0,0)
c##      END DO
c##      END IF
c    !..............!   Cseal>0 & K=Nyi

C -------------------------------------------- ! ........ Solve by ........
 400  CALL TDMACPLX(V1W, V1E, JJMAX)           !         TDMA algorithm
c ---------------------------------------------!

      DO J=JVMIN, JVMAX                        ! ..........................
          JJ=J+JJJ                             !
          V1(J, K)=B1(JJ)                      ! Order solution
      END DO                                   ! --------------------------
C ---------------------------------------------!

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calcp1                                                      **
C **                                                                         **
C **  CALCP1:  Solves pressure correc equation for first order solution.     **
C **                                                                         **
C *****************************************************************************
c Modified on 9/17/93 to include variation of density: p=p*+p'

      SUBROUTINE CALCP1(FU, FP, FN, FS, SIGMA, Dir, KV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /SOURCEA/ PRATIO, CORIF, SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /KVALA/ DYVK,DYPK,SVNK,SVSK
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP
      COMMON /HONEY/ HCELL,HCDIM

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /DPUV1/ DU1(MAXNXTP2), DV1(MAXNXTP2), DVV1(MAXNXTP2)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /TDMA1/ A1(MAXNXTP2), B1(MAXNXTP2),
     +               C1(MAXNXTP2), D1(MAXNXTP2)

      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /KVALB/ K, KM1, KP1
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN
      COMMON /LRBOUND/ LEFTBC, RIGHTBC
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION DXP, DXU, HP, HU,HV,SUW,SUE,DYP,DYV,SVN,SVS,
     +                 U, V, P, Rhop, Emup, Drhop, Drhot,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 PRATIO, CORIF, SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 DYVK, DYPK, SVNK, SVSK,PMAX, PMAXP, PEPS,
     +                 AC, ETA, RELAXH, HCELL,HCDIM

      DOUBLE COMPLEX U1,V1,P1,T1,DU1, DV1, DVV1,A1, B1, C1, D1

      INTEGER ITER, ITMAX, ITPMAX, ISOLN, K, KM1, KP1,
     +        JPMIN, JPMAX, JPSTART, JPSTOP, LIFT,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT
      INTEGER IFULL, LEFTBC, RIGHTBC, BEARING,Countp,MAXCOUNTP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DPU, FWPO, SIGMA, DPV,Dir,PPmag,Pmag,
     +                 FEPO, FSNO, AMASS, DZERO, suej,suwj,
     +                 svskm1,svnkm1,rhow,rhoe,rhos,rhon,
     +                 Bw,Be,Bn,Bs,FE0,FW0,FS0,FN0,AW0,AE0,AS0,AN0,
     +                 FU(*), FP(*), FN(*), FS(*)
      DOUBLE COMPLEX AWP,AEP,FEP,FWP,PW1,PE1,PPJK,MASS,FSP,FNP,
     +                 Imag, ASP, ANP,zero, ACeta, PN1, PS1,
     +                 rhon1,rhos1,rhoe1,rhow1,rhop1,Isigma
      INTEGER J, JJ, JJJ, JJMAX, Jp1,KV, Dirk
C ----------------------------------------------------------------------------
C --  CALCP1 code                                                           --
C ----------------------------------------------------------------------------
C KV=K except when K=+/- 1 for asymmetric bearing=3

      zero=(0.0D0,0.0D0)
      Imag=(0.0D0,1.0D0)                       ! imaginary unit
      DZERO=0.0D0                              !

      ACeta=AC/DCMPLX(1.0D0,ETA)               ! COMPLIANCE FACTOR
      Isigma=DCMPLX(0.0D0,Sigma)               ! i x FREQUENCY PARAM

      svskm1=svs(k)                            !
      svnkm1=svn(k)                            !

      Dirk=Dir
      IF ((BEARING.EQ.1).AND.(IABS(K).EQ.1)) Dirk=-Dir

      JJJ=2-JPMIN                              !
      JJMAX=JPMAX+JJJ                          !

c    !.........................................!.............................
      j=jpstart                                !-> On west face Pcv
      rhow=rhop(j,k)*sue(j)+rhop(jpmin,k)*suw(j)
      PW1= P1(j,k)*sue(j) + P1(jpmin,k)*suw(j) ! P1

      rhow1=( Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k) )*sue(j) +
     +( Drhop(jpmin,k)*P1(jpmin,k)+Drhot(jpmin,k)*T1(jpmin,k))*suw(j)

      Dpu=Hu(j,k)*Dypk                         !
      Fw0=Dpu*rhow*U(j,k)                      ! zeroth order west flow
      Fwp=Dpu*(U1(j,k)*rhow+U(j,k)*rhow1) +    ! coeff of W(est) flow
     +    Dypk*rhow*ACeta*PW1*U(j,k)           ! with compliance
      Fwpo=Dypk*U(j,k)*Fu(j)*Rhow              !

      IF (LEFTBC.EQ.0) THEN       !.........................!
         AWP=ZERO                 ! SPEC: FLOW
         Bw=0.0D0                 !
      ELSE                        !.........................!
         AWP=Dpu*DU1(J)*rhow      ! SPEC: PRESSURE
         Bw=Drhop(jpmin,k)/rhop(jpmin,k)
      END IF                      !.........................!

C   !.........................!                ! ---------------------------
      DO J=JPMIN, JPMAX                        ! Sweep from left-right on
C   !.........................!                !
          JJ=J+JJJ                             ! Pcv's determines flow on
          jp1=J+1                              !............................
          suej=sue(j)
          suwj=suw(J)

      rhoe=rhop(j,k)*suej+rhop(jp1,k)*suwj     !..........................!
      PE1= P1(j,k)*suej+P1(jp1,k)*suwj         ! on E(ast) face of P-CV
      rhoe1=( Drhop(j,k)*P1(j,k)   +Drhot(j,k)*T1(j,k)     )*suej +
     +      (Drhop(jp1,k)*P1(jp1,k)+Drhot(jp1,k)*T1(jp1,k) )*suwj
      Be=Drhop(jp1,k)/rhop(jp1,k)

      Dpu=Hu(j,k)*Dypk                         !
      Fe0=Dpu*rhoe*U(j,k)                      ! zeroth order east flow
      Fep=Dpu*(U1(j,k)*rhoe+U(j,k)*rhoe1)+     ! Fe=DY He (U1*Rho+Uo*Rho1)e +
     +        Dypk*rhoe*ACeta*U(j,k)*PE1       !    DY Ac (Uo*Rho*P1)e/(1+in)

          Aep=Dpu*Du1(j)*rhoe
c                                              !...........................!
c                                              ! -> ON SOuth face Pcv
      rhos=rhop(j,km1)*Svnkm1+rhop(j,k)*Svskm1 !
      PS1 =P1(j,km1)*Svnkm1+P1(j,k)*Svskm1     !
      rhos1=(Drhop(j,km1)*P1(j,km1)+Drhot(j,km1)*T1(j,km1) )*Svnkm1
     +     +(Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k) )*Svskm1

c##   Bs=Drhop(j,km1)/rhop(j,km1)              !C## not really needed

      Dpv=Hv(j,KV)*Dxp(j)                      !
      Fs0=Dpv*rhos*V(j,KV)*Dir                 ! zeroth order south flow
      Fsp=Dpv*(V1(j,KV)*rhos+V(j,KV)*rhos1)+   ! Fs=Dx Hs (V1*Rho+Vo*Rho1)s
     +      Dxp(j)*rhos*ACeta*V(j,KV)*PS1      ! Dx Ac (Vo*Rho*P1)s/(1+in)

          Asp=Dpv*Dvv1(J)*rhos
c                                              !...........................!
c                                              ! -> ON north face Pcv
      rhon=rhop(j,k)*Svnk+rhop(j,kp1)*Svsk     !
      PN1=P1(j,k)*Svnk+P1(j,kp1)*Svsk          !
      rhon1=(Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k))*Svnk +
     +      (Drhop(j,kp1)*P1(j,kp1)+Drhot(j,kp1)*T1(j,kp1))*Svsk

c##   Bn=Drhop(j,kp1)/rhop(j,kp1)              !C## not really needed

      Dpv=Hv(j,kp1)*Dxp(j)                     !
      Fn0=Dpv*rhon*V(j,kp1)*Dir                ! zeroth order north flow
      Fnp=Dpv*(V1(j,kp1)*rhon+V(j,kp1)*rhon1)+ !Fn=Dx Hn (V1*Rho+Vo*Rho1)n
     +       Dxp(j)*rhon*ACeta*V(j,kp1)*PN1    !   Dx Ac (Vo*Rho*P1)n/(1+in)

          Anp=Dpv*Dv1(J)*rhon
c                                              !...........................!
c                                              ! Mass Balance on Pcv
          rhop1=Drhop(j,k)*P1(j,k)             !...........................!
     +         +Drhot(j,k)*T1(j,k)

          Mass=Fwp-Fep+(Fsp-Fnp)*Dir

c        !......................!
          IF (ISOLN.EQ.1) THEN                 ! ...........................
C        !......................!              ! Only for zeroth first comp.
            Fepo=Dypk*U(j,k)*Fu(j)*rhoe        ! DX H1e (URho)e
            Fsno=Dir*Dxp(j)*(FS(j)*V(j,k)*rhos-FN(j)*V(j,kp1)*rhon)

          Mass=Mass+Fwpo-Fepo+Fsno -Isigma*Dxp(j)*Dypk*
     + ( rhop(j,k)*(Fp(j)+ACeta*P1(j,k)) +
     +       rhop1*(Hp(j,k)+HCDIM      )               )

            Fwpo=Fepo                          !

          ELSE
c        !.......................!             ! Other components

          Mass=Mass-Isigma*Dxp(j)*Dypk*(
     +         rhop1*(Hp(j,k)+HCDIM) +
     +         rhop(j,k)*ACeta*P1(j,k)  )

          END IF                               ! ...........................
C        !......................!              ! Coeffs. for P1
          AW0=DMAX1(+FW0,DZERO)                ! AW upwind coeff.
          AE0=DMAX1(-FE0,DZERO)                ! AE upwind coeff.
          AS0=DMAX1(+FS0,DZERO)                ! AS upwind coeff.
          AN0=DMAX1(-FN0,DZERO)                ! AN upwind coeff.

          C1(JJ)=-AWP-AW0*Bw                   ! coeffs. for TDMA soln.
          A1(JJ)=-AEP-AE0*Be                   !
          B1(JJ)=AEP+AWP+ANP+ASP               !
     +   +(AE0+AW0+AN0+AS0)*Drhop(j,k)/rhop(j,k)+
     +    Isigma*Dxp(j)*Dypk*Drhop(j,k)*Hp(j,k)

          D1(JJ)=MASS                          ! P'(j, k-1)=P(j, k+1)=0

          FW0=FE0                              !for next j:
          Fwp=Fep                              !downstream -> upstream
          AWP=AEP                              !values
          Bw=Be                                !at west face
          PW1=PE1                              !

          AMASS=CDABS(MASS)                    !
          SMASS=SMASS+AMASS                    ! Global mass source
          IF (AMASS.GT.MMP) MMP=AMASS          ! MAX(Local mass source)

C   !............................!             !
      END DO                                   ! ...........................
C -------------------------------------------- ! j=jpmin,..,jpmax


      IF (RIGHTBC.EQ.0) THEN                   ! SPEC FLOW ON EAST BOUNDARY
        B1(JJMAX)=B1(JJMAX)+A1(JJMAX)          ! On extended lands, Ue spec
        A1(JJMAX)=ZERO                         ! & Aep = 0.0
      END IF
C                                              ! ...........................
  100 PW1=zero                                 ! Spec. pressure at left B
      PE1=zero                                 ! Spec. pressure at right B
C                                              ! ...........................
C ...................................          !
  200 CALL TDMACPLX(PW1, PE1, JJMAX)           ! Solve by TDMA algorithm
C ...................................          !

      IF ((INTER.EQ.0).AND.(IFULL.EQ.1)) THEN  ! P1'(Nxt)=P1'(1)
          B1(JJMAX+1)=B1(2)                    ! P1'(j) = B1(jj)
      END IF                                   !
C ............................................ ! ...........................


C ----------------------                       !
  300 DO J=JPMIN, JPMAX                        ! Correct U, V, & P fields
C ----------------------                       ! ::::::::::::::::::::::::

          JJ=J+JJJ                                 !
          PPJK=B1(JJ)                              ! P'(j, k)
          P1(J, K)=P1(J, K)+PPJK*ALFP              ! Corrected P1(j, k),
          U1(J, K)=U1(J, K)+(PPJK-B1(JJ+1))*DU1(J) ! Corr. U1(j, k) velocity
          V1(J,KV)=V1(J,KV)-PPJK*DVV1(J)*Dirk      ! Corrected south V1 velocity
          V1(J, KP1)=V1(J, KP1)+PPJK*DV1(J)*Dir    ! Corrected north V1 velocity
          DVV1(J)=DV1(J)                           ! For use on next K

          Ppmag=CDABS(PPJK)                    ! check convergence
          IF (Ppmag.gt.PMAXp) PMAXp=Ppmag      ! in pressures
          IF (Ppmag.LE.PEPS) Countp=Countp+1   !

          Pmag=CDABS(P1(J,K))
          IF (Pmag.GT.PMAX) PMAX=Pmag

C ----------------------                       !
      END DO                                   ! ...........................
C ----------------------                       !

      IF (LEFTBC.GT.0) THEN                    ! On inter-lands:
          J=JPSTART                            ! Correct Ue velocity
          U1(J, K)=U1(J, K)-B1(2)*DU1(J)       ! Jpstart=Jumin
      END IF                                   ! ...........................

      IF ((IFULL.EQ.1).AND.(INTER.EQ.0)) THEN  !
          P1(JPSTOP, K)=P1(JPMIN, K)           ! On extended lands
          V1(JPSTOP,KV)=V1(JPMIN, KV)          ! jpstop=Nxt, jpmin=1
          U1(JPSTOP, K)=U1(JPMIN, K)           ! Inforce periodicity
          V1(JPSTOP, KP1)=V1(JPMIN, KP1)       !
      END IF                                   ! ...........................
C..............................................!

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Tdmacplx                                                    **
C **                                                                         **
C **  TDMACPLX:  Solution of complex algebraic tridiagonal system of         **
C **             equations.                                                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TDMACPLX(ZLEFT, ZRIGHT, JMAX)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TDMA1/ A(MAXNXTP2), B(MAXNXTP2),
     +               C(MAXNXTP2), D(MAXNXTP2)
      DOUBLE COMPLEX A, B, C, D
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE COMPLEX DEN, ZLEFT, ZRIGHT
      INTEGER J, JM1, JMAX
C ----------------------------------------------------------------------------
C --  TDMACPLX code                                                         --
C ----------------------------------------------------------------------------

      D(1)=ZLEFT                          ! = F(1)
      A(1)=(0.0D0, 0.0D0)                 ! = E(1)
      DO J=2, JMAX                        ! 1st sweep:
          JM1=J-1                         ! Find recursive
          DEN=1.0/(C(J)*A(JM1)+B(J))      !
          A(J)=-A(J)*DEN                  ! = E(j)
          D(J)=(D(J)-C(J)*D(JM1))*DEN     ! = F(J)
      END DO                              !
      B(JMAX+1)=ZRIGHT                    ! Set BC at right
      DO J=JMAX, 1, -1                    ! 2nd sweep:
          B(J)=A(J)*B(J+1)+D(J)           ! Calculate soln.
      END DO                              ! stored in vector B

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Edgeu1                                                      **
C **                                                                         **
C **  EDGEU1:  Updates U1 and edge P1 due to inertia effects.                **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE EDGEU1(FU, FP, SIGMA,Dir, K, INERP, NPA)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /PRECES/ PRECL, PRECR
      COMMON /TRECES/ TRECL, TRECR

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /PRECES1/ PRECL1, PRECR1

      COMMON /UVEC/ JUMIN, JUMAX, JUSTART, JUSTOP
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /SOLN/ ISOLN
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPAA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /LRBOUND/ LEFTBC, RIGHTBC

      DOUBLE PRECISION DXP,DXU, SUE, SUW, HP, HU, HV, TRECL, TRECR,
     +                 U, V, RHOP,EMUP, Drhop,Demup , PRECL, PRECR,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC,HRECD,
     +                 Drhot,Demut
      DOUBLE COMPLEX U1, V1, P1,T1, PRECL1, PRECR1
      INTEGER JUMIN,JUMAX,JUSTART,JUSTOP,JPMIN,JPMAX,JPSTART,JPSTOP,
     +        NPOCKET,NLC,NPC,NLA,NPAA,NPAP1,NXI,NYI,NXT,IFULL,
     +        ISOLN, LEFTBC, RIGHTBC
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION FU(*), FP(*), Sigma, Renr, Alfar, Hedge,Dir,
     +                 rhmu, F1WR, F1WI, F1ER, F1EI, hrat,UW, UE,
     +                 zero,emue,emuw,rhoe,rhow,eta,ax,renh,rhorat,
     +                 rhoem, emuem ,CR,BR,TR

      DOUBLE COMPLEX  emu1,rhoe1,rhow1,imag, Pedge
      INTEGER J, JP1, K, NPA, INERP, KEDGE
C ----------------------------------------------------------------------------
C --  EDGEU1 code                                                           --
C ----------------------------------------------------------------------------
C NOTE: IF K=NPA, CORRECT P1 ON EDGEV1 FOR BOUNDARY WITH RECESS
C### WARNING Justart = JUMIN ????  for pad bearing if LBC=1  ????

       KEDGE=0
       IF (IABS(K).EQ.NPA) KEDGE=1
       zero=0.0D0
       Imag=(0.0D0,1.0D0)

c     -------------------------------
c ... AT Left boundary with a recess :
c     --------------------------------
 17   IF ((KEDGE.EQ.1).AND.(LEFTBC.EQ.2)) THEN
          U1(JUSTART,K)=(0.0D0,0.0D0)
          GOTO 18
      END IF

      j=jumin
      rhow=rhop(j,k)
      rhoe=(rhow+rhop(j+1,k))/2.0D0
      emuw=emup(j,k)
      rhow1=Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k)
      rhoe1=(rhow1+Drhop(j+1,k)*P1(j+1,k)
     +            +Drhot(j+1,k)*T1(j+1,k))/2.0D0
      emu1=Demup(j,k)*P1(j,k)+Demut(j,k)*T1(j,k)
      rhorat=rhoe/rhow
      hrat=hu(j,k)/hp(j,k)
      Uw=rhorat*U(j,k)*hrat                          !Uw veloc.

      U1(justart,k)=Rhorat*U1(j,k)*hrat +            ! West inlet U1 vel
     +    (U(j,k)*hrat*rhoe1-Uw*rhow1)/rhow +        ! .................
     +    Imag*Sigma*Dxp(j)*rhow1/rhow/2.0D0         !

      IF (ISOLN.EQ.1) THEN            !.......................!
         f1wr=(rhorat*U(j,k)*Fu(j)-Uw*Fp(j))/Hp(j,k)
         f1wi=sigma*fp(j)*dxp(j)/2.0D0/Hp(j,k)
         U1(justart,k)=U1(justart,k)+DCMPLX(f1wr,f1wi)
      END IF                          !.......................!

      IF ((LEFTBC.EQ.2).AND.(INERP.EQ.1).AND.(UW.GT.ZERO)) THEN
c    !................................................!
         rhmu=rhow/emuw                               !
         Hedge=hrecd+Hp(j,k)                          ! Pdrop at edge due
         Eta=Hp(j,k)/Hedge                            ! to Bernoulli effect
c        RHOem=Rhop(Justart,K)                        !
c        EMUem=Emup(Justart,K)                        !
         CALL LOCPROPS(rhoem,emuem,cr,br,TR,PRECL,TRECL) ! Rhoe-, Pe-
         Ax=klosxd*rhow*(1.0D0-(Eta*rhow/rhoem)**2)   !

C##         IF (SWIRL.GT.ZERO) THEN     !................! EDGEU1
C##            Renh=Renc*rhmu*hp(j,k)                    ! removed on 10/95
C##            Ax=Ax*(1.0D0+1.95D0/(Renh**0.43))         !
C##            Renr=Renc*rhoem*Hedge/emuem
C##            ALFAR=ASPEC*((Renr**.681)/7.752963)/Hedge/Hedge
C##            PRECL1=PRECL1+emu1*ALFAR*(1.0D0-Eta*(Rhow/Rhoem))*SWIRL
C##         END IF                      !................!
C##c NOTE:  here emu1 on PRECL1 should be emu1 evaluated at Prec  ! 10/23/91 !##

         P1(j,k)=PRECL1-Ax*Uw*(2.0D0*U1(justart,k)+Uw*rhow1/rhow)

c    !.....................
      ELSE
c    !....................

         P1(j,k)=PRECL1

      END IF !........................................! Inerp=1 & Uw>0
c    !.......!

c     -------------------------------
c ... AT Right boundary with a recess
c     -------------------------------
 18   IF ((KEDGE.EQ.1).AND.(RIGHTBC.EQ.2)) THEN
          U1(JUSTOP,K)=(0.0D0,0.0D0)
          RETURN
      END IF
      j=jumax
      rhoe=rhop(jpstop,k)
      rhow=(rhoe+rhop(j,k))/2.0D0
      emue=emup(jpstop,k)
      rhoe1=Drhop(justop,k)*P1(justop,k)+Drhot(justop,k)*T1(justop,k)
      emu1 =Demup(justop,k)*P1(justop,k)+Demut(justop,k)*T1(justop,k)
      rhow1=(rhoe1+Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k))/2.0D0
      rhorat=rhow/rhoe
      hrat=hu(j,k)/hp(jpstop,k)
      Ue=rhorat*U(j,k)*hrat                           !Fue=Fuw sinve Vo=0

      U1(justop,k)=Rhorat*U1(j,k)*hrat+               !east inlet U1 vel
     +         (U(j,k)*hrat*rhow1-Ue*rhoe1)/rhoe -   !...................
     +          Imag*Sigma*rhoe1*dxp(jpstop)/rhoe/2.0D0

      IF (ISOLN.EQ.1) THEN            !.......................!
         f1er=(rhorat*U(j,k)*Fu(j)-Ue*Fp(jpstop))/Hp(jpstop,k)
         f1ei=-sigma*fp(jpstop)*dxp(jpstop)/2.0D0/Hp(jpstop,k)
         U1(justop,k)=U1(justop,k)+DCMPLX(f1er,f1ei)
      END IF                          !.......................!

      IF ((RIGHTBC.EQ.2).AND.(INERP.EQ.1).AND.(UE.LT.ZERO)) THEN
c    !................................................!
c        RHOem=RHOP(JPSTOP+1,K)
c        EMUem=EMUP(JPSTOP+1,K)
         CALL LOCPROPS(rhoem,emuem,cr,br,
     +   TR,PRECR,TRECR)                              ! Rhoe-, Pe-
         rhmu=rhoe/emue                               ! Pdrop at edge due
         Eta=Hp(jpstop,k)/(hrecd+Hp(jpstop,k))        ! to Bernoulli effect
         Ax=klosxu*rhoe*(1.0D0-(Eta*rhoe/rhoem)**2)   !

         P1(jpstop,k)=PRECR1-Ax*Ue*(2.0D0*U1(justop,k) +
     +                   Ue*rhoe1/rhoe )
c    !...................!
      ELSE
c    !...................!
         P1(jpstop,k)=PRECR1
      END IF  !........................................! Inerp=1 & ue<0
c    !........!

 19   RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Edgev1                                                      **
C **                                                                         **
C **  EDGEV1:  Updates inlet V1 velocity and edge pressure P1 due to         **
C **           inertia effects.                                              **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE EDGEV1(FU,FP,FN,SIGMA,Dir,DYPK,KLOSY,Svnk,Svsk,Kstep)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT),DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /URECJET/ UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD
      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECJET1/ PREC1do(MAXNPOCK),PREC1up(MAXNPOCK),
     +                 PR1jet(MAXNPOCK,MAXNPOCK+2)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /SOLN/ ISOLN
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION DXP, DXU, HP, HU, HV,SUW, SUE,
     +                 U, V, RHOP, EMUP, Drhop,Demup,Drhot,Demut,
     +                 PREC,TREC,QREC,QIN,QOUT,UREC,
     +                 QFACTOR,PMAX,PMAXP,PEPS,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC, ASPEC, HRECD

      DOUBLE COMPLEX U1, V1, P1,T1,PREC1do,PREC1up,PR1jet

      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +     ISOLN, INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL,
     +     BEARING, Countp, MAXCOUNTP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYPK, Svnk, SVsk, SIGMA, FU(*), FP(*),FN(*),
     +                 F1RS, F1IS, DYP, FEW, FEW1, FACT1, VS,KLOSY,
     +                 rhoem,rhos,emus,rhon, eta, ay, hrat,
     +                 rhorat,Dir,Ppmag, APiner ,cr,br,tr, Uedge
      DOUBLE COMPLEX Imag, rhos1,rhon1,emu1, Piner
      INTEGER I, J, JM1, K, KP1, JR, Icase, Kstep,II
C ----------------------------------------------------------------------------
C --  EDGEV1 code                                                           --
C ----------------------------------------------------------------------------

C    !-------------------------!               !----------------!
      IF ((BEARING.EQ.1).AND.(NPOCKET.GT.0)) THEN !=> FOR HJBS
C    !-------------------------!               !----------------!

      DYP=DYPK/2.0D0                                 !
      K=Kstep*Npa                                    !
      KP1=K+Kstep                                    !
      imag=(0.0D0,1.0D0)                             !

c    !-----------------!                             !...............................
      DO I=1, NPOCKET                                !
C    !-----------------!                             ! Sweep from
C                                                    ! left to right Vcvs on rec edge
C                                                    ! ........................
          jr=I*NXI+1                                 !
          j=(i-1)*NXI+NLC                            !On left recess edge
          Uedge=UREC(I,1)*(1.0D0+HRECD/hu(j,k))
          fact1=dypk*Uedge/dxp(j)                    !..          few=fact1*(hu(j,k)-hp(j,k))                !
          few1=fact1*(fu(j)-fp(j))                   !
          Icase=1                                    !
          II=0                                       ! II=1,...NPC within pocket
          goto 704

 701      V1(j,K)=0.5D0*V1(j,K)                      !
          Icase=2                                    !

 702      j=j+1                                      !
          if(j.eq.jr) goto 703                       !

          jm1=j-1                                    ! along interior of reces
          Uedge=UREC(I,II+1)*(1.0D0+HRECD/hu(j,k))
          fact1=dyp*Uedge/dxp(j)                     !....
          few=fact1*(hu(j,k)-hu(jm1,k))              !
          few1=fact1*(fu(j)-fu(jm1))                 !
          goto 704                                   !

 703      jm1=j-1                                    !On right recess edge
          Uedge=UREC(I,NPC)*(1.0D0+HRECD/hp(j,k))    !
          fact1=dypk*Uedge/dxp(j)                    !...
          few=fact1*(hp(j,k)-hu(jm1,k))              !
          few1=fact1*(fp(j)-fu(jm1))                 !
          Icase=3

 704      rhos=rhop(j,K)                             !INTERPOLATE
          rhon=rhos*Svnk+rhop(j,KP1)*Svsk            !props. on CV faces

          rhos1=Drhop(j,K)*P1(j,K)+Drhot(j,K)*T1(j,K)
          rhon1=rhos1*Svnk+ (Drhop(j,KP1)*P1(j,KP1)+
     +                       Drhot(j,KP1)*T1(j,KP1) )*Svsk

          rhorat=rhon/rhos                           !
          hrat=Hv(j,KP1)/Hp(j,k)                     !
c
          few=few/Hp(j,K)                            ! MOD: 040491: CORRECTIONS
          Vs=V(j,KP1)*hrat*rhorat+Dir*few

          V1(j,k)=rhorat*hrat*V1(j,KP1) +            ! Pert. South Vel V1
     +            (rhon1*V(j,KP1)*hrat +
     +             rhos1*(-Vs+Dir*few+Imag*Sigma*dyp))/rhos

c        !.......................!
          IF (ISOLN.EQ.1) THEN
             f1rs=(Dir*few1+rhorat*V(j,KP1)*fn(j)-Vs*fp(j))/Hp(j,k)
             f1is=Sigma*Dyp*fp(j)/hp(j,k)
             V1(j,k)=V1(j,k)+DCMPLX(f1rs,f1is)
          END IF
c        !.......................!

          II=II+1                ! II=1,...NPC within pocket edge
c        !.................................!
          IF ((INERP.EQ.1).AND.(VS*DIR.GT.(0.0D0))) THEN
c        !.................................! CORRECT EDGE PRESSURE
c            rhoem=RHOP(jr-1,1)            ! find RHOe-, Pe-
             CALL LOCPROPS(rhoem,emus,CR,BR,TR,PREC(I),TREC(I))      !##should be PRjet

             Eta=hp(j,k)/(hp(j,k)+hrecd)
             Ay=Klosy*rhos*(1.0D0-(Eta*rhos/rhoem)**2)
             P1(j,K)=PR1jet(I,II)-Ay*Vs*(2.0D0*V1(j,k)+Vs*rhos1/rhos)
c        !..............!
          ELSE
c        !..............!
             P1(j,K)=PR1jet(I,II)
c        !..................................!
          END IF
         !..................................! PR1jet: perturbed recess pressure

          GOTO (701,702,705) , Icase

 705      V1(j,k)=0.5D0*V1(j,k)

c   !----------------------------!
      END DO
c   !----------------------------! END FOR HJB, K=1, NPOCKET




C    !-------------------------!               !----------------!
      ELSE IF (BEARING.EQ.2) THEN              ! FOR SEALS
C    !-------------------------!               !----------------!
       imag=(0.0D0,1.0D0)
       K=1
C     !----------------!
       DO J=2, Nxt                                     ! Find Inlet V1
c     !----------------!                               ! & correct P1

        jm1=J-1

        Few1=SWIRL*(Dypk/Dxp(j))*(Fu(j)-Fu(jm1))     !
        Few =SWIRL*(Dypk/Dxp(j))*(Hu(j,1)-Hu(jm1,1))/Hp(j,1)  !C### corr 04/91

        rhos=rhop(j,1)
        rhon=rhos*svnk+rhop(j,2)*svsk
        rhos1=Drhop(j,1)*P1(j,1)+Drhot(j,1)*T1(j,1)
        rhon1=rhos1*Svnk+(Drhop(j,2)*P1(j,2) +
     +                    Drhot(j,2)*T1(j,2)  )*Svsk
        rhorat=rhon/rhos
        hrat=Hv(j,2)/Hp(j,1)

        Vs=V(j,1)
        f1rs=(few1+rhorat*V(j,2)*fn(j)-Vs*fp(j))/Hp(j,1)
        f1is=Sigma*Dypk*fp(j)/hp(j,1)

        V1(j,k)=rhorat*hrat*v1(j,2) +
     +            (rhon1*V(j,2)*hrat  +
     +      rhos1*(-Vs+few+Imag*Sigma*dypk))/rhos +
     +      DCMPLX(f1rs,f1is)

c        !.................................!
          if ((inerp.eq.1).and.(Vs.gt.(0.0))) then    !.....................
c        !.................................!          !Correct Edge Pressure
             AY=Klosy*rhos                            !by Bernoulli effect

             Piner =   AY*Vs*(2.0*V1(j,k) +
     +                      Vs*rhos1/rhos)

             APiner=CDABS(Piner)

             Ppmag=CDABS(P1(j,1)+Piner)

             IF (Ppmag.LE.PEPS) Countp=Countp+1
             IF (PMAXp.LT.APiner) PMAXp=APiner

             IF (APiner.GT.PMAX) PMAX=APiner

             P1(J,1)= -Piner

c        !..............!
          else
c        !..............!

             P1(j,1)=(0.0D0,0.0D0)
             Countp=Countp+1

c        !..................................!
          end if
         !..................................!

c     !----------------!
       END DO
c     !----------------! => END FOR SEAL
       K=1
       KP1=2
C    !-------------------------!
      END IF
C    !-------------------------! END of edge corrections!



  300 IF (IFULL.EQ.1) THEN
      V1(1,K )=V1(NXT,K )                            ! 1st vcv = Nxt(last) Vcv
      P1(1,K )=P1(NXT,K )                            !  periodicity conditions
      V1(1,KP1)=V1(NXT,KP1)                          ! ........................
      END IF
      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Endseal1                                                    **
C **                                                                         **
C **  ENDSEAL1: Updates V1 velocity and discharge pres P1 due to end seal .  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ENDSEAL1(Dir,DYPK, Sealcoef,alfp, Kstep)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /PCOUNT/ PMAX, PMAXP, PEPS, COUNTP, MAXCOUNTP

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/ T1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN

      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 U, V, P, RHP, EMP, Drhop,Drhot,
     +                 PMAX, PMAXP, PEPS

      DOUBLE COMPLEX U1,V1,P1, T1

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,ISOLN,
     +        COUNTP, MAXCOUNTP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION zero,DYPK, Dir, Sealcoef,Ppcor, P1mag,
     +                 Alfp, Rhoe,Rhow, Fe, Fw, Vend
      DOUBLE COMPLEX Vend1, Pend1, Rho1
      INTEGER J,KK, Kstep
C ----------------------------------------------------------------------------
C --  ENDSEAL1 CODE                                                          --
C ----------------------------------------------------------------------------
c Dir: flow director: =1.0 Right side ob HJB, 0<y -> <L/D
c                     =-1.0 Left side of HJB, 0: y -> -L/D
c............................................................................
      IF (IFULL.NE.1) THEN
       PRINT *, 'CALCSOLNZ> ENDSEAL ERROR !'
       PRINT *, 'THIS PROGRAM ONLY UPDATED FOR 1 PAD FULL 2PI HJB'
       CALL BEEPER
       STOP
      END IF

      KK=Kstep*NYI
      zero=0.0D0

      j=NXT-1
      Rhow=Sue(j)*Rhp(j,kk)+Suw(j)*Rhp(1,kk)
      Fw=DYPK*Rhow*U(j,KK)*Hu(j,KK)

C    !.........................!               !..........................
      DO J=1, NXT-1                            ! Sweep from
C    !.........................!               ! left Vcv to right vcv
        Rhoe=Sue(j)*Rhp(j,kk)+Suw(j)*Rhp(j+1,kk)
        Fe=DYPK*Rhoe*U(j,KK)*Hu(j,KK)

c##        Vend=(V(j,KK)*HV(j,KK)-Dir*(Fe-Fw)/Dxp(j)/Rhp(j,KK))/HP(j,KK)

         Vend=V(j,KK)                          ! Advect exit velocities
         Vend1=V1(j,KK)                        ! ........................!
         Rho1 = Drhop(j,KK)*P1(j,KK)+Drhot(j,KK)*T1(j,KK)

          IF (Dir*Vend.GT.zero) THEN
               Pend1=Sealcoef*(2.0D0*Rhp(j,KK)*Vend1+
     +                         Rho1*Vend)*Vend
               Pend1=P1(j,KK)*(1.0D0-Alfp)+Pend1*Alfp
          ELSE
               Pend1=(0.0D0,0.0D0)
          END IF

               Ppcor=CDABS(Pend1-P1(j,KK))
               P1mag=CDABS(Pend1)
               IF (Ppcor.LT.PEPS) Countp=Countp+1
               IF (P1mag.GT.PMAX) PMAX=P1mag

               P1(j,KK)=Pend1

               Fw=Fe

C    !........................!
      END DO                                   ! ----------------------------
C    !........................!                !
      IF (IFULL.EQ.1) P1(NXT,KK) =P1(1,KK)     ! set periodicity conditions

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Flowo1                                                      **
C **                                                                         **
C **  FLOWO1:  Calculates components of 1st order flow independent of P1     **
C **           pressures.                                                    **
C **          = INT ( Rho0 V0 H1 dT)  &&  INT(Rho H1 dx dy )on recess border**
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FLOWO1(FU, FP, YO, IFM)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT), SUE(MAXNXT)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /QSIDE0H/ QSIDEH(MAXNPOCK)
      COMMON /FLOW01/ QO1(MAXNPOCK), QIN, QOUT
      COMMON /So/ Soo(MAXNPOCK)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION DXP, DXU, SUW, SUE, XP, XU, YP, YV,
     +                 DYP, DYV, SVN, SVS, U, V, RHOP, EMUP,
     +                 Soo,QO1, QIN, QOUT ,QSIDEH
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,
     +        ICSTEP,IFULL, BEARING

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION QR, svsk, svnk, factor, rho, Dir,YO,YPO,YVO
      DOUBLE PRECISION Ss, FU(*), FP(*)
      INTEGER I,J,K,JMIN,JMAX,jp1,nyim1,Kstep,Kend,Icsym,IFM
C ----------------------------------------------------------------------------
C --  FLOWO1 code                                                           --
c     Qo1= INT { RHOo H1 Vo.n dT } on (U,V) grid recess boundary
c     Soo= INT { RHOo H1 dArea} on area between recess and grid
C ----------------------------------------------------------------------------
c Dir: flow director: = 1.0 Right side ob HJB, 0<y -> <L/D
c                     =-1.0 Left side of HJB,  0: y -> -L/D
c .................................
      QOUT=0.0D0
      QIN=0.0D0

      Kstep=1
      Dir=1.0D0
      Icsym=1

      YPO=1.0D0
      YVO=YPO

c    !........................!
      IF (NPOCKET.GT.0) THEN
c    !........................! HJB
      DO I=1, NPOCKET
        QO1(I)=QIN
        Soo(I)=QIN
        QSIDEH(I)=QIN
      END DO
c    !........................!
      ELSE
c    !........................! SEAL/BEARING
        GOTO 88
c    !........................!
      END IF
c    !........................!

 77   Kend=Kstep*NPA

      Svsk=SVS(Kstep*NPAP1)
      Svnk=SVN(Kstep*NPAP1)

C    !................!                      ! Sweep over recesses
      DO I=1, NPOCKET                        ! -------------------
C    !................!
          Ss=0.0D0
          QR=0.0D0                           ! INLET flow
          JMIN=(I-1)*NXI+NLC                 ! for i-th recess
          JMAX=I*NXI+1                       ! is equal to
          K=Kend+Kstep                       !

          IF (IFM.eq.1) THEN
            YVO=YV(K)-YO
            YPO=YP(K)-YO
          END IF
C        !.................!                 !
          DO J=JMIN, JMAX                    ! Flow at top of recess = SUM(V Fp)
C        !.................!                 !
          rho=Rhop(j,Kend)*Svnk+Rhop(j,K)*Svsk
          QR=QR+DXP(J)*YVO*V(J,K)*FP(J)*rho*Dir

          Ss=Ss+Rhop(j,Kend)*YPO*Fp(j)*Dxp(j)*Dyp(Kend)/2.0D0
C        !.................!                 !
          END DO                             !
C        !.................!                 !

          QSIDEH(I)=QSIDEH(I)+QR             !
          factor=1.0D0
          Jmin=Jmin-1                        ! C

C        !.................!                 !
          DO K=Kstep,Kend,Kstep              ! (Right-Left) Circ. flow (U Fu)
C        !.................!                 !
             IF (K.eq.Kend) Factor=2.0D0
             IF (IFM.EQ.1) THEN
               YPO=YP(K)-YO
             END IF

             j=jmin
             rho=rhop(j,k)*sue(j)+rhop(j+1,k)*suw(j)
             Qr=Qr-dyp(K)*rho*YPO*Fu(j)*U(jmin,k)

             j=j+1
             Ss=Ss+Rhop(j,k)*YPO*Fp(j)*Dxp(j)*Dyp(k)/2.0D0/Factor

             j=jmax
          IF ((j.eq.NXT).AND.(IFULL.EQ.1)) THEN
             jp1=2
          ELSE
             jp1=j+1
          END IF

             rho=rhop(j,k)*sue(j)+rhop(jp1,k)*suw(j)
             Qr=Qr+Dyp(k)*rho*YPO*Fu(j)*U(jmax,k)

             Ss=Ss+Rhop(j,k)*YPO*Fp(j)*Dxp(j)*Dyp(k)/2.0D0/Factor

C        !.................!                 !
          END DO                             !
C        !.................!                 !
          QO1(I)=QO1(I)+QR                   ! = Flow at i-th recess
          Soo(I)=Soo(I)+Ss
          QIN=QIN+QR
C    !.....................!                 !
      END DO                                 ! .....................
C    !---------------------!                 ! i=1,2,.... Npocket
      GOTO 89

 88   IF (IFM.eq.1) YVO=YV(1)-YO

C    !.....................!                 !
      DO J=1, NXT-1                          ! INLET FLOW for seal
C    !.....................!                 !
          QIN=QIN+DXP(J)*YVO*FP(J)*Dir*V(J,1)*Rhop(J,1)
C    !.....................!                 !
      END DO                                 !
C    !.....................!                 !

 89   NYIm1=Kstep*(NYI-1)
      Kend=Kstep*NYI
      Svsk=SVS(Kend)
      Svnk=SVN(Kend)

      IF (IFM.eq.1) YVO=YV(Kend)-YO

C    !.....................!                 !
      DO J=1, NXT-1                          ! Outlet flow at y=L/D
C    !.....................!                 !
          rho=rhop(j,NYIm1)*Svnk+rhop(j,Kend)*Svsk
          QOUT=QOUT+DXP(J)*YVO*FP(J)*Dir*V(J,Kend)*Rho
C    !.....................!                 !
      END DO                                 !
C    !.....................!                 !

      IF (IFULL.EQ.1) GOTO  666              !=> 360 deg BEARING
c    !-------------------------!             !...................!
          j=NXT-1
          jp1=NXT
C        !.................!
          DO K=Kstep, Kend, Kstep
C        !.................!             !..............SIDE PAD FLOW
          IF (IFM.eq.1) YPO=YP(K)-YO
          rho=rhop(1,K)*Sue(1)+rhop(2,K)*Suw(1)
          QOUT=QOUT-DYP(K)*rho*YPO*FU(1)*U(1,K)
          rho=rhop(j,K)*Sue(j)+rhop(jp1,K)*Suw(j)
          QOUT=QOUT+DYP(K)*rho*YPO*FU(J)*U(J,K)
C        !.................!
          END DO
C        !.................!             !........................

C    !---------------------!
 666  IF (ISYM.EQ.0) THEN
c    !.....................!
       IF (BEARING.EQ.2) goto 99

       IF(Icsym.eq.2) goto 99
          Kstep=-1
          Dir=-1.0D0
          Icsym=2
          IF (NPOCKET.EQ.0) GOTO 89
          GOTO 77
c    !......................!
      ELSE
c    !......................!
         IF (IFM.eq.1) THEN
            Factor=0.0D0
         ELSE
            Factor=2.0D0
         END IF

         IF (NPOCKET.GT.0) THEN
             DO I=1, NPOCKET
              QO1(I)=Factor*QO1(I)
              Soo(I)=Factor*Soo(I)
              QSIDEH(I)=Factor*QSIDEH(I)
             END DO
         END IF

         QIN=Factor*QIN
         QOUT=Factor*QOUT
c    !.......................!
      END IF
c    !.......................!


 99   END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Flows1                                                      **
C **                                                                         **
C **  FLOWS1:  Calculates inlet and outlet flow for first order fields.      **
C **          = INT ( Rho0 V1 H0 + RHO1 V0 H0 + AC/(1+in) V0 P1 ) dT         **
C **          &  {Rho1 H0 + Rho0 Ac/(1+in) P1} dx dy on recess border        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FLOWS1( QABS)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /FLOW01/ QO1(MAXNPOCK), QINO1, QOUTO1
      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK),
     +                QREC1(MAXNPOCK), QIN1, QOUT1
      COMMON /S1/ S11(MAXNPOCK)
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

c     ................................................................
      DOUBLE PRECISION DXP, DXU, HP, HU, HV, SUW, SUE,
     +                 DYP, DYV, SVN, SVS, U, V, Drhot,Demut,
     +                 Rhop,Emup,Drhop,Demup, QO1, QINO1, QOUTO1,
     +                 AC, ETA, RELAXH
      DOUBLE COMPLEX U1,V1,P1,T1, PREC1,TREC1,QREC1,S11,QIN1,QOUT1
      INTEGER ISOLN,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,
     +        ICSTEP, IFULL, BEARING, LIFT
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION QABS, rho, svsk, svnk, factor, Dir
      DOUBLE COMPLEX   Ss, Qr,rho1, ACeta, PN1, PU1
      INTEGER I, J, JP1, K, JMIN, JMAX, NYIM1, Kstep, Kend, Icsym

C ----------------------------------------------------------------------------
C --  FLOWS1 code: calculates                                                           --
c     Qrec1= INT [{Ho(RHO1 Vo+RHOo V1)+RHOo Vo ACeta P1}.n dTr]
c     on (U,V) grid recess boundary
c     S11= INT { (RHO1 HO+RHOo ACeta P1) dArea} on area between recess and grid
C ----------------------------------------------------------------------------
      Nyim1=Nyi-1
      QIN1=(0.0D0, 0.0D0)
      QOUT1=QIN1
      QABS=0.0D0
      ACeta=AC/DCMPLX(1.0D0,ETA)         ! COMPLIANCE Ac/(1+in)

      Kstep=1
      Icsym=1
      Dir=1.0D0

c    !.......................!
      IF (NPOCKET.GT.0) THEN
c    !.......................!
          DO I=1, NPOCKET
          QREC1(I)=QIN1
          S11(I)=QIN1
          END DO
c    !.......................!
      ELSE
c    !.......................!
        GOTO 88
c    !.......................!
      END IF
c    !.......................!

 77   Kend=Kstep*NPA

      Svsk=SVS(Kstep*NPAP1)
      Svnk=SVN(Kstep*NPAP1)

C   !-------------------                      !
      DO I=1, NPOCKET                         ! ------------------------------
C   !-------------------                      !
          QR=(0.0D0, 0.0D0)                   ! Inlet flow
          Ss=(0.0D0, 0.0D0)

          JMIN=(I-1)*NXI+NLC                  ! for i-th Pocket
          JMAX=I*NXI+1                        !

          K=Kend+Kstep
C        !...............!                    !
          DO J=JMIN, JMAX                     ! Flow at top of rec = SUM(H V1)
C        !...............!                    !

          rho= rhop(j,Kend)*Svnk+rhop(j,k)*Svsk
          rho1=(Drhop(j,Kend)*P1(j,Kend)+Drhot(j,Kend)*T1(j,Kend))*Svnk
     +        +(Drhop(j, k)*P1(j,k)+Drhot(j, k)*T1(j,k) )*Svsk
          PN1=P1(j,Kend)*Svnk+P1(j,k)*Svsk

          QR=QR+DXP(J)*Dir*(Hv(J,k)*(V1(j,k)*rho+V(j,k)*rho1) +
     +                       rho*V(j,k)*ACeta*PN1)

          rho1=Drhop(j,Kend)*P1(j,Kend)+Drhot(j,Kend)*T1(j,Kend)

          Ss=Ss+0.5D0*Dxp(j)*Dyp(Kend)*
     +         ( rho1*Hp(j,Kend) + rho*ACeta*PN1)

          END DO                              ! DO J
C        !...............!                    !

          JMIN=JMIN-1                         ! Plus

          Factor=1.0D0
C        !...............!                    !
          DO K=Kstep,Kend,Kstep               ! (Right-Left) Circ. flow (H VU1)
C        !...............!                    !
          IF (K.eq.Kend) Factor=2.0D0

          j=jmin
          jp1=j+1
          rho=rhop(j,k)*sue(j)+rhop(jp1,k)*suw(j)
          rho1=(Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k) )*sue(j) +
     +         (Drhop(jp1,k)*P1(jp1,k)+Drhot(jp1,k)*T1(jp1,k))*suw(j)
          PU1=P1(j,k)*sue(j)+P1(jp1,k)*suw(j)
          Qr=Qr-Dyp(k)*(Hu(j,k)*(rho*U1(j,k)+rho1*U(j,k) ) +
     +                  rho*U(j,k)*ACeta*PU1 )

          j=jp1
          rho1=Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k)
          Ss=Ss+Dxp(j)*Dyp(k)/2.0D0/Factor*
     +         (rho1*Hp(j,k)+rhop(j,k)*ACeta*P1(j,k))

          j=jmax
          IF ((j.eq.NXT).AND.(IFULL.EQ.1)) THEN
             jp1=2
          ELSE
             jp1=j+1
          END IF

          rho=rhop(j,k)*sue(j)+rhop(jp1,k)*suw(j)
          rho1=(Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k))*sue(j) +
     +         (Drhop(jp1,k)*P1(jp1,k)+Drhot(jp1,k)*T1(jp1,k))*suw(j)
          PU1=P1(j,k)*sue(j)+P1(jp1,k)*suw(j)

          Qr=Qr+Dyp(k)*( Hu(j,k)*(rho*U1(j,k)+rho1*U(j,k)) +
     +                   rho*U(j,k)*ACeta*PU1 )

          rho1=Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k)
          Ss=Ss+Dxp(j)*Dyp(k)/2.0D0/Factor*
     +         (rho1*Hp(j,k)+rhop(j,k)*ACeta*P1(j,k))
C        !................!                   !
          END DO
C        !................!                   !

          QREC1(I)=QREC1(I)+QR                ! Perturbed flow at i-th recess
          S11(I)=S11(I)+ Ss
          QIN1=QIN1+QR
C    !------------------!                     !
      END DO                                  ! ------------------------------
C    !------------------!                     !
      GOTO 89



 88   CONTINUE
c    !..................!
      DO J=1, NXT-1                           ! INLET flow at Y=0
c    !.................!                      ! FOR SEAL
      rho = rhop(j,1)
      rho1= Drhop(j,1)*P1(j,1)+Drhot(j,1)*T1(j,1)
      QIN1=QIN1+Dxp(j)*Dir*( Hp(j,1)*(V1(j,1)*rho+V(j,1)*rho1) +
     +                       rho*V(j,1)*ACeta*P1(j,1) )

C    !.................!                      !
      END DO
C    !.................!                      !


 89   NYIM1=Kstep*(NYI-1)
      Kend=Kstep*NYI
      Svsk=SVS(Kend)
      Svnk=SVN(Kend)
c    !..................!
      DO J=1, NXT-1                           ! Outlet flow at L/R side
c    !.................!                      !
      rho = rhop(j,Kend)*Svsk+rhop(j,NYIm1)*Svnk
      rho1=(Drhop(j,Kend)*P1(j,Kend)+Drhot(j,Kend)*T1(j,Kend))*Svsk+
     + (Drhop(j,NYIm1)*P1(j,NYIm1)+Drhot(j,NYIm1)*T1(j,NYIm1))*Svnk
      PN1 = P1(j,Kend)*Svsk+P1(j,NYIm1)*Svnk
      QOUT1=QOUT1+Dxp(j)*Dir*(
     +            Hv(j,Kend)*(V1(j,Kend)*rho + V(j,Kend)*rho1 ) +
     +            rho*V(j,Kend)*ACeta*PN1 )

C    !.................!                      !
      END DO
C    !.................!                      !


c    !-------------------------!              !...................!
      IF (IFULL.EQ.1) GOTO  666               !=> 360deg bearing
c    !-------------------------!              !...................!

          j=NXT-1
          jp1=NXT
C        !.................!
          DO K=Kstep, Kend, Kstep
C        !.................!             !..............SIDE PAD FLOW
          rho=rhop(1,k)*Sue(1)+rhop(2,k)*Suw(1)
          rho1=(Drhop(1,k)*P1(1,k)+Drhot(1,k)*T1(1,k))*Sue(1) +
     +         (Drhop(2,k)*P1(2,k)+Drhot(2,k)*T1(2,k))*Suw(1)
          QOUT1=QOUT1-Dyp(k)*Hu(1,k)*(rho*U1(1,k)+rho1*U(1,k))

          rho=rhop(j,k)*Sue(j)+rhop(jp1,k)*Suw(j)
          rho1=(Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k) )*Sue(j) +
     +         (Drhop(jp1,k)*P1(jp1,k)+Drhot(jp1,k)*T1(jp1,k))*Suw(j)
          QOUT1=QOUT1+Dyp(k)*Hu(j,k)*(rho*U1(j,k)+rho1*U(j,k))
C        !.................!
          END DO
C        !.................!             !........................
C At pad inlet/outlet no compliance effect since P=Pback here, i.e
C        there is no pressure deformation P1=0


C    !---------------------!
 666  IF (ISYM.EQ.0) THEN
c    !.....................! ASYMMETRIC BEARING
       IF (BEARING.eq.2) GOTO 99

       IF(Icsym.eq.2) goto 99

          Kstep=-1
          Dir=-1.0D0
          Icsym=2
          IF (NPOCKET.EQ.0) GOTO 89
          GOTO 77
c    !......................!
      ELSE
c    !......................! SYMMETRIC BEARING
         Factor=2.0D0

         IF (NPOCKET.GT.0) THEN
           DO I=1, NPOCKET
            QREC1(I)=Factor*QREC1(I)
            S11(I)=Factor*S11(I)
           END DO
         END IF

         QIN1=Factor*QIN1
         QOUT1=Factor*QOUT1

c    !.......................!
      END IF
c    !.......................!

 99   QIN1=QIN1+QINO1*ISOLN                   !
      QOUT1=QOUT1+QOUTO1*ISOLN                ! ------------------------------

c    !..........................!
      IF (NPOCKET.GT.0) THEN
c    !..........................!
      DO I=1, NPOCKET
        QABS=QABS+CDABS(QREC1(I))+ISOLN*DABS(QO1(I)) ! Total absolute flow
      END DO
c    !..........................!
      ELSE
c    !..........................!
        QABS=CDABS(QIN1)        !ZABS on VAX
c    !..........................!
      END IF
c    !..........................!

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Force1                                                      **
C **                                                                         **
C **  FORCE1:  Integrates first-order pressure field to find forces.         **
C **           & moments                                                     **
C *****************************************************************************

      SUBROUTINE FORCE1(FX, FY, MX, MY, YO, IFM)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT), SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT),
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /PARRAY/ P(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECJET1/ PREC1do(MAXNPOCK),PREC1up(MAXNPOCK),
     +                 PR1jet(MAXNPOCK,MAXNPOCK+2)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBSYM/ ISYM,ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION DXP, DXU, SUW, SUE, DYP, DYV, SVN, SVS, P,
     +                 XP, XU, YP, YV, COSXP, SINXP, COSXU, SINXU,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      DOUBLE COMPLEX U1,V1,P1,PREC1do,PREC1up,PR1jet
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, IFULL, BEARING
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYVK, AREA, YY, YO
      DOUBLE COMPLEX FX, FY, MX, MY,PP,PCL,PSL,PCR,PSR,ZERO
      INTEGER I, J, K, K0, KM1, JMI, JMA, JP1,IJ,
     +        Kstart, Kend, Kstep, Icsym, IFM
C ----------------------------------------------------------------------------
C --  FORCE1 code  : Calculates Forces Fx&Fy and Moments Mx&My
C ----------------------------------------------------------------------------
      ZERO=(0.0D0, 0.0D0)
      Icsym=1
      Fx=zero
      Fy=zero
      Mx=zero
      My=zero

      !write(6,*) PP

      Kstep=1
      K0=ISYM

 777  IF (NPOCKET.EQ.0) THEN  ! => SEAL or BEARING
           Kstart=2*Kstep
           GOTO 222
      ELSE
           Kstart=Kstep
           Kend=Kstep*NPA
      END IF

C --------------------                             !
      DO I=1, NPOCKET-IFULL+1                      ! Sweep on inter-recess lands
C --------------------                             ! ---------------------------
          KM1=K0                                   !
          JMI=(I-1)*NXI+2                          ! limits for land pressures
          JMA=JMI+NLC-2                            ! betweeen recesses.
C                                                  !
c        !..........................!              !
          DO K=Kstart, Kend, Kstep                 ! from bottom to top of rec
c        !..........................!              !
              DYVK=DYV(K)                          !  Y size of Vcv
              J=JMI-1                              !  first pressure index
              PP=ZERO                              !
              IF (P(J, K).GT.PCAV) PP=P1(J,K)      !
              IF (P(J, KM1).GT.PCAV) PP=PP+P1(J,KM1)
              PCL=PP*COSXP(J)                      ! left P*COS(Angle)
              PSL=PP*SINXP(J)                      ! left P*SIN(Angle)
              YY=YV(K)-YO                          !

              DO J=JMI, JMA                        !
c            !...............!                     !
                  AREA=DYVK*DXU(J-1)/4.0D0         ! Area of integration/4.0
                  PP=ZERO                          !
                  IF (P(J, K).GT.PCAV) PP=P1(J,K)  !
                  IF (P(J, KM1).GT.PCAV) PP=PP+P1(J,KM1)
                  PCR=PP*COSXP(J)                  ! right P*COS(Angle)
                  PSR=PP*SINXP(J)                  ! right P*SIN(Angle)
                  FX=FX+(PCR+PCL)*AREA             ! X-force
                  FY=FY+(PSR+PSL)*AREA             ! Y-force
                  MY=MY+(PCR+PCL)*AREA*YY          ! Y-Moment
                  MX=MX-(PSR+PSL)*AREA*YY          ! X-Moment
                  PCL=PCR                          !
                  PSL=PSR                          !
              END DO                               !
c            !...............!                     ! J=JMI,JMA

              KM1=K

              IF (I.GT.NPOCKET) GOTO 111

              IJ=0
c            !...............!                     ! .......................
              DO J=JMA, I*NXI                      ! At interior of recess I
c            !...............!                     ! .......................
                JP1=J+1                            !
                AREA=DXU(J)*DYVK/2.0D0             ! Area of integration/2.0
                IJ=IJ+1
                FX=FX+AREA*(PR1jet(I,IJ  )*COSXP(J)+
     +                      PR1jet(I,IJ+1)*COSXP(JP1))
                FY=FY+AREA*(PR1jet(I,IJ  )*SINXP(J)+
     +                      PR1jet(I,IJ+1)*SINXP(JP1))
                MY=MY+AREA*(PR1jet(I,IJ  )*COSXP(J)+
     +                      PR1jet(I,IJ+1)*COSXP(JP1))*YY
                MX=MX-AREA*(PR1jet(I,IJ  )*SINXP(J)+
     +                      PR1jet(I,IJ+1)*SINXP(JP1))*YY
c            !...............!                     ! .......................
              END DO                               !
c            !...............!                     ! .......................

c       !......................!                   !
 111      END DO                                   ! k=1,2,.. Npa*Kstep
c       !......................!                    !

        END DO                                     ! i=1,2, Npocket
c   !--------------------------!                   !.......................

        Kstart=NPAP1*Kstep

C---!..........................!...................!....................!
 222    Kend=NYI*Kstep
C---!

C    !......................!           !------------------
      DO K=Kstart,Kend,Kstep            ! On EXTENDED LANDS
C    !......................!           !------------------
          KM1=K-Kstep                              !
          DYVK=DYV(K)                              ! Y size of Vcv
          YY=YV(K)-YO                              ! moment arm
          J=1                                      ! Start of integration
          PP=ZERO
          IF (P(J, K).GT.PCAV) PP=P1(J, K)         !
          IF (P(J, KM1).GT.PCAV) PP=PP+P1(J, KM1)  !
          PCL=PP*COSXP(J)                          ! left P*COS(Angle)
          PSL=PP*SINXP(J)                          ! left P*SIN(Angle)

c        !.....................!                   !
          DO J=2, NXT                              !
c        !.....................!                   !
              AREA=DYVK*DXU(J-1)/4.0D0             ! Area of integration/4.0
              PP=ZERO                              !
              IF (P(J, K).GT.PCAV) PP=P1(J, K)     !
              IF (P(J, KM1).GT.PCAV) PP=PP+P1(J, KM1)
              PCR=PP*COSXP(J)                      ! right P*COS(Angle)
              PSR=PP*SINXP(J)                      ! right P*SIN(Angle)
              FX=FX+(PCR+PCL)*AREA                 ! X-force
              FY=FY+(PSR+PSL)*AREA                 ! Y-force
              MY=MY+(PCR+PCL)*AREA*YY              ! Y-Moment
              MX=MX-(PSR+PSL)*AREA*YY              ! X-Moment
              PCL=PCR                              !
              PSL=PSR                              !
          END DO                                   ! j=2,3,.., Nxt
c        !.....................!                   !

      END DO                                       ! k=Kstart,... Nyi
c    !-------------------------!                   !

c..................................................!.........................!
      IF (ISYM.EQ.0) THEN
c    !.....................!         ! ASYMMETRIC HJB
       IF (BEARING.EQ.2) GOTO 999    !=> annular Seal

       IF (Icsym.eq.2) THEN
          IF (MODEL.EQ.1) THEN
             Fx=Fx*2.0D0             ! 2row HJB
             Fy=Fy*2.0D0
             Mx=zero
             My=zero
          END IF
          GOTO 999
       END IF

       Icsym=2
       Kstep=-1
       GOTO 777
c    !.....................!
      ELSE
c    !.....................!                       ! SYMMETRIC HJB       IF (IFM.eq.0) THEN  ! IFM=0 translations for Yo=0
          FX=FX*2.0D00
          FY=FY*2.0D00
          MX=zero
          MY=zero
        ELSE               ! IFM=1 rot. displc. for Yo=0
          FX=zero
          FY=zero
          MX=MX*2.0D0
          MY=MY*2.0D0
        END IF
c    !.....................!
      END IF
c    !.....................!


 999  END

C
C ================================================================!
C support subroutines:

      SUBROUTINE uminmax(x, xmin,xmax,jmin,jmax,kmin,kmax)

      IMPLICIT NONE

      INCLUDE 'params.f'

      integer jmin,jmax,kmin,kmax, j ,k
      double precision X(maxnxt,-maxnyi:maxnyi),xmin,xmax

      xmin=x(jmin,kmin)
      xmax=x(jmin,kmin)

      do j=jmin, jmax
        do k=kmin,kmax
          if (X(j,k).gt.Xmax) Xmax=X(j,k)
          if (X(j,k).lt.Xmin) Xmin=X(j,k)
        end do
      end do

      end


      SUBROUTINE vminmax(x, xmin,xmax,jmin,jmax,kmin,kmax)

      IMPLICIT NONE

      INCLUDE 'params.f'

      integer jmin,jmax,kmin,kmax, j ,k
      double precision X(maxnxt,-maxnyi:maxnyi),xmin,xmax

      xmin=x(jmin,kmin)
      xmax=x(jmin,kmin)
      do j=jmin, jmax
        do k=kmin,kmax
          if (X(j,k).gt.Xmax) Xmax=X(j,k)
          if (X(j,k).lt.Xmin) Xmin=X(j,k)
        end do
      end do

      end
c
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c compe2t.f last revised 3/11/94  by Dr. Luis SanAndres
c           updated for cell depth effect on 6/20/95 by LSA
c           modified for angled injection on 9/15/95
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::