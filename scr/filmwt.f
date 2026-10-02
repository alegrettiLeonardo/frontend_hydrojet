C NOTE:
c  For compliant bearings, two film thickness arrays are required:
C  COM HOFILM holds the conventional (undeformed) film thickness
c  COM  HFILM saves the deformed film thickness such that Hdef = H + a P
C
C ######     #    #       #    #  #    #   #####          ######
C #          #    #       ##  ##  #    #     #            #
C #####      #    #       # ## #  #    #     #            #####
C #          #    #       #    #  # ## #     #     ###    #
C #          #    #       #    #  ##  ##     #     ###    #
C #          #    ######  #    #  #    #     #     ###    #

c  hydroflex.f  Copyright LuisS anAndres / TexasA&MUniversity / 1994
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR II
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C *****************************************************************************
C **                                                                         **
C **  Subroutine Filmh                                                       **
C **                                                                         **
C **  FILMH: Determines film thickness at Xp and Xu coordinates.             **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FILMH(EXO, EYO, NXT, NYI)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT), 
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /HOFILM/ HP(MAXNXT,-MAXNYI:MAXNYI), 
     +                HU(MAXNXT,-MAXNYI:MAXNYI),
     +                HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /HFILM/  HPN(MAXNXT,-MAXNYI:MAXNYI), 
     +                HUN(MAXNXT,-MAXNYI:MAXNYI),
     +                HVN(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /SPLDATA/Z(Nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl), Nj

      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON /PADPOS/ PRELOAD, OFFSET,ROTDEL
      COMMON /LOBES/ PRELOADB, NLOBES
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /NODES/  NPOC,NLCd,NPCd,NLAd,NPAd,NPAP1d,NXId,NYId,
     +                NXTd,IFULL
c     ...............................................................
      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,HPN,HUN,HVN,
     +                 XP, XU, YP, YV, Z,Cl,Bcl,Ccl,Dcl,
     +                 COSXP, SINXP, COSXU, SINXU
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD, AR,  HREC,
     +                 AXo,AYo,Zo,ClearO,ClearR,ClearL,YR,YL, XJ,
     +                 PRELOAD, OFFSET, ROTDEL, PRELOADB, 
     +                 AC, ETA, PBACK, angleP, PI2
      INTEGER ISYM, ICSTEP, BEARING, NLOBES, IFULL, LIFT
      INTEGER NPOC,NLCd,NPCd,NLAd,NPAd,NPAP1d,NXId,NYId,NXTd  ! <- DUMY VALUES

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER NXT, NYI, I, J, NJ, K, Jstart,Jend,Jstep, Index, Iset,
     +        Ilarge
      DOUBLE PRECISION EXO,EYO,XSIP,XSIV,CP,CV,ZERO,R,PI,
     +                 DXOp,DYOp, DXOv,DYOv,Px,M,Size,Xpivot,Crot,
     +                 Coepadc,Coepads 
      DOUBLE PRECISION Xlobel(1:MAXNYI), Xlobet(1:MAXNYI),
     +          CosXlobe(1:MAXNYI),SinXlobe(1:MAXNYI)
C ----------------------------------------------------------------------------
C --  FILMH code                                                            --
C ----------------------------------------------------------------------------
c    dimensionless film thickness is equal to:
c
c    H= C(Z)/C* + { Ex+Ay(Z-Zo)/C*} cos(X) + { Ey-Ax(Z-Zo)/C*} sin(X)
c
c            - (rp/C*) cos(X-Xpivot) - (R/C*) Delta sin(X-Xpivot)
c
c    Y=Z, and 
c    where C(Z) pad machined clearance, C*= TYP clearance
c          rp   = C - Cm : pad preload  [m], Cm: bearing clearance
c          Delta tilt-pad rotation about pivot located at Xpivot.
c          Delta =0 for fixed pad bearings.         
c
c    and OFFSET=(Xpivot-Xleading)/Pad_size, TYP = 0.5
c    
c    ex,ey are journal center displacements [m]
c    dx,dy are journal axis rotations [rad] about Zo.
c 
c    Nlobes option only available for 360deg bearings.
c.............................................................................
c NOTE:
c HP, HU, HV contain film  thickness w/o compliance or cavitation effects
c HPN,HUN,HVN contain film thickness which will be modified for compliance/cavitation
c
c.............................................................................
      ZERO=0.0D0
      R=DIAM/2.0D0
      PI=DACOS(-1.0D0)
      PI2=2.0D0*PI
      Ilarge=0                                ! pointer for H < 0
      Jstart=1
      IF (ISYM.eq.0) Jstart=0

      IF (IFULL.EQ.1) GOTO 777                !==> 360 deg bearing

c.............................................!.....................!
c  BEARING PAD data needed
c.............................................! FOR BEARING TILT PADS
      M=PRELOAD/CLEAR                         ! dim. preload factor
      Size=XP(NXT)-XP(1)                      ! pad size in rads.
      Xpivot=OFFSET*Size+XP(1)                ! loc. of pivot in rads.
      Crot=R*ROTDEL/CLEAR                     ! dim. pad rotation
C##   Crot=0.0D0                              ! ==> for fixed pads.
      CoePadc=-M*DCOS(Xpivot)+Crot*DSIN(Xpivot)
      CoePads=-M*DSIN(Xpivot)-Crot*DCOS(Xpivot)  
c.............................................!...............................!


      IF (ICSTEP.eq.1) goto 20                !-=> STEP Clearance HJBearing


c     ........................................!..............................
c     RIGHT HALF of Bearing, 0=<y<=Lr/R
c     ........................................!..............................
        Px=1.0D0
        Jstart=1
        Jend=NYI
        Jstep=1
        Iset=1
c     !------------------------!              !..............................
 10    DO   j=Jstart,Jend,Jstep               ! Sweep in axial direction 
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]
 
        call seval(Nj,Px*Xsip,Z,CL,BCL,CCL,DCL,Cp) !clearance in [m] at the
        call seval(Nj,Px*Xsiv,Z,CL,BCL,CCL,DCL,Cv) !nodal points
        Cp=Cp/Clear                             
        Cv=Cv/Clear                              

        DXOp=EXO+AYO*(Xsip-ZO)/Clear+Coepadc
        DYOp=EYO-AXO*(Xsip-ZO)/Clear+Coepads
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear+Coepadc
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear+Coepads

c                                               !.............................
        DO   I=1, NXT                           
          HP(I,J)=CP+DXOp*COSXP(I)+DYOp*SINXP(I) 
          HU(I,J)=CP+DXOp*COSXU(I)+DYOp*SINXU(I)
          HV(I,J)=CV+DXOv*COSXP(I)+DYOv*SINXP(I)

          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)
          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
        END DO                                  !..sweep in circumferential direction

c     !------------------------!              !..............................
      END DO                                  ! sweep in axial direction
c     !------------------------!              !..............................

c     ........................................!..............................
      IF ((ISYM.EQ.1).OR.(BEARING.EQ.2)) THEN ! IF bearing is symetrical
         GOTO 99                              !    or an annular seal then EXIT
      ELSE IF (Iset.EQ.2) THEN                !
         GOTO 99                              !
      END IF                                  !
c     ........................................!..............................
c     LEFT HALF of Bearing, -Ll/R<y<0
c     ........................................!..............................
      Px=1.0D0
      IF (Nj.eq.2) Px=-1.0D0                  !=> tapered clearance both sides
      Jstart=-1
      Jend=-NYI
      Jstep=-1
      Iset=2
      GOTO 10

C:::::::::::::::::  STEPPED CLEARANCE BEARING Icstep=1 :::::::::::::::::::::::::
c     RIGHT HALF of Bearing, 0=<y<=Lr/R
c     ........................................!..............................

 20     Jend=NYI
        Jstep=1
        Px=1.0D0

c     !------------------------!              !..............................
        DO   j=Jstart,Jend,Jstep               ! Sweep in axial direction 
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]

        IF (Xsip.gt.YR) THEN                  !
           Cp=ClearR/Clear                    !   ----------. ClearO  
        ELSE                                  !             .--------- ClearR
           Cp=ClearO/Clear                    !   --------------------
        END IF                                !  y=0        YR        Lr

        IF (Xsiv.gt.YR) THEN
           Cv=ClearR/Clear
        ELSE
           Cv=ClearO/Clear
        END IF

        DXOp=EXO+AYO*(Xsip-ZO)/Clear+Coepadc
        DYOp=EYO-AXO*(Xsip-ZO)/Clear+Coepads
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear+Coepadc
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear+Coepads

        DO   I=1, NXT                           
          HP(I,J)=CP+DXOp*COSXP(I)+DYOp*SINXP(I) 
          HU(I,J)=CP+DXOp*COSXU(I)+DYOp*SINXU(I)
          HV(I,J)=CV+DXOv*COSXP(I)+DYOv*SINXP(I)
          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)
          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
        END DO             !......................! sweep circ. direction

c     !------------------------!              !..............................
      END DO				      ! sweep in axial direction
c     !------------------------!              !..............................

c     ........................................!..............................
      IF ((ISYM.EQ.1).OR.(BEARING.EQ.2)) THEN ! IF bearing is symetrical
         GOTO 99                              !    or an annular seal then EXIT
      END IF
c     ........................................!..............................

c     ........................................!..............................
c     LEFT HALF of Bearing, -Ll<y<0           ! Asymmetric Bearing
c     ........................................!..............................
      Jstart=-1
      Jstep=-1
      Jend=-NYI
      
c     !------------------------!              !..............................
       DO  j= Jstart,Jend,Jstep               ! Sweep in axial direction
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]

        IF (Xsip.gt.YL) THEN                  !
           Cp=ClearO/Clear                    !               .-------- ClearO  
        ELSE                                  ! ClearL ........
           Cp=ClearL/Clear                    !   ---------------------
        END IF                                !  -Ll          YL      y=0

        IF (Xsiv.gt.YL) THEN
           Cv=ClearO/Clear
        ELSE
           Cv=ClearL/Clear
        END IF
        DXOp=EXO+AYO*(Xsip-ZO)/Clear+Coepadc
        DYOp=EYO-AXO*(Xsip-ZO)/Clear+Coepads
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear+Coepadc
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear+Coepads

        DO   I=1, NXT                           
          HP(I,J)=CP+DXOp*COSXP(I)+DYOp*SINXP(I) 
          HU(I,J)=CP+DXOp*COSXU(I)+DYOp*SINXU(I)
          HV(I,J)=CV+DXOv*COSXP(I)+DYOv*SINXP(I)
          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)
          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
        END DO				      !.. sweep in circ. direction

c     !------------------------!              !..............................
      END DO
c     !------------------------!              !sweep in axial direction

      GOTO 99






C....................................................................
C....................................................................
C   FOR IFULL=1 ,   360 deg. bearing with NLOBES
C....................................................................
C....................................................................

 777  IF (NLOBES.EQ.0) NLOBES=1

      M=PRELOADB/CLEAR                        ! dimensionless preload factor
      DO K=1, NLOBES
         Xlobel(K)=2*(K-1)*PI/Nlobes          ! leading edge of lobe
         Xlobet(K)=2*K*PI/Nlobes              ! trailing edge of lobe
         anglep=(2*K-1)*PI/Nlobes             ! loc. of pivot: half-way
         COSXLobe(K)=-M*DCOS(anglep)          ! 
         SINXlobe(K)=-M*DSIN(anglep)          !
      END DO                                  !
c.............................................!...............................!

      IF (ICSTEP.eq.1) goto 40                !-=> STEP Clearance HJBearing

c     ........................................!..............................
c     RIGHT HALF of Bearing, 0=<y<=Lr/R
c     ........................................!..............................
      Jend=NYI
      Jstep=1
      Px=1.0D0
      Iset=1

c     !------------------------!              !..............................
 30    DO   j=Jstart,Jend,Jstep               ! sweep in axial direction
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]
 
        call seval(Nj,Px*Xsip,Z,CL,BCL,CCL,DCL,Cp) !clearance in [m] at the
        call seval(Nj,Px*Xsiv,Z,CL,BCL,CCL,DCL,Cv) !nodal points
        Cp=Cp/Clear                             
        Cv=Cv/Clear                              

        DXOp=EXO+AYO*(Xsip-ZO)/Clear
        DYOp=EYO-AXO*(Xsip-ZO)/Clear
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear

c      !....................!                   !.............................
        DO   I=1, NXT                           !
c      !....................!                   ! sweep in circ. direction.
          XJ=XP(I)
          IF (XJ.LT.0.0D0) XJ=XJ+PI2
          K=1
  311     IF ((XJ.LT.XlobeT(K)).AND.(XJ.GE.XlobeL(K))) THEN
           goto 312                          
          ELSE
           IF (K.GE.NLOBES) GOTO 312
           K=K+1
           GOTO 311
          END IF

  312     HP(I,J)=CP+(DXOp+CosXlobe(K))*COSXP(I)+
     +               (DYOp+SinXlobe(K))*SINXP(I) 
          HU(I,J)=CP+(DXOp+CosXlobe(K))*COSXU(I)+
     +               (DYOp+SinXlobe(K))*SINXU(I)
          HV(I,J)=CV+(DXOv+CosXlobe(K))*COSXP(I)+
     +               (DYOv+SinXlobe(K))*SINXP(I)
          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)

          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
c      !....................!                   
        END DO
c      !....................! I: sweep in circ. direction

c     !------------------------!              !..............................
      END DO
c     !------------------------!              ! J: Axial sweep ........
 
c     ........................................!..............................
      IF ((ISYM.EQ.1).OR.(BEARING.EQ.2)) THEN ! IF bearing is symmetrical
        GOTO 99                               !    or an annular seal then EXIT
      ELSE IF (Iset.EQ.2) THEN                !
        GOTO 99                               !
      END IF                                  !
c     ........................................!..............................
c     LEFT HALF of Bearing, -Ll/R<y<0
c     ........................................!..............................
      Px=1.0D0
      IF (Nj.eq.2) Px=-1.0D0                  !=> tapered clearance both sides
      Jstart=-1
      Jend=-NYI
      Jstep=-1
      Iset=2
      GOTO 30

C:::::::::::::::::  STEPPED CLEARANCE BEARING Icstep=1 :::::::::::::::::::::::::
c     ........................................!..............................
c     RIGHT HALF of Bearing, 0=<y<=Lr/R
c     ........................................!..............................

c     !------------------------!              !..............................
 40     DO   j=Jstart, NYI                    ! sweep in axial direction
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]

        IF (Xsip.gt.YR) THEN                  !
           Cp=ClearR/Clear                    !   ----------. ClearO  
        ELSE                                  !             .--------- ClearR
           Cp=ClearO/Clear                    !   --------------------
        END IF                                !  y=0        YR        Lr

        IF (Xsiv.gt.YR) THEN
           Cv=ClearR/Clear
        ELSE
           Cv=ClearO/Clear
        END IF

        DXOp=EXO+AYO*(Xsip-ZO)/Clear
        DYOp=EYO-AXO*(Xsip-ZO)/Clear
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear

c      !....................!                   !.............................
        DO   I=1, NXT                           !
c      !....................!                   ! sweep on circ. direction.
          XJ=XP(I)
          IF (XJ.LT.0.0D0) XJ=XJ+PI2
          K=1
  321     IF ((XJ.LT.XlobeT(K)).AND.(XJ.GE.XlobeL(K))) THEN
           goto 322                          
          ELSE
           IF (K.GE.NLOBES) GOTO 322
           K=K+1
           GOTO 321
          END IF

  322     HP(I,J)=CP+(DXOp+CosXlobe(K))*COSXP(I)+
     +               (DYOp+SinXlobe(K))*SINXP(I) 
          HU(I,J)=CP+(DXOp+CosXlobe(K))*COSXU(I)+
     +               (DYOp+SinXlobe(K))*SINXU(I)
          HV(I,J)=CV+(DXOv+CosXlobe(K))*COSXP(I)+
     +               (DYOv+SinXlobe(K))*SINXP(I)
          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)

          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
c      !....................!                   
        END DO
c      !....................! sweep in X-direction

c     !------------------------!              !..............................
      END DO				      ! sweep in axial direction
c     !------------------------!              !..............................

c     ........................................!..............................
      IF ((ISYM.EQ.1).OR.(BEARING.EQ.2)) THEN ! IF bearing is symmetrical
        GOTO 99                               !    or an annular seal then EXIT
      END IF                                  !
c     ........................................!..............................

c     ........................................!..............................
c     LEFT HALF of Bearing, -Ll/R<y<0         ! Asymmetric Bearing
c     ........................................!..............................

c     !------------------------!              !..............................
      DO   j= -1, -NYI, -1                    ! sweep in axial direction
c     !------------------------!              !..............................
        Xsip=YP(j)*R                          ! axial coordinate of Pnode [m]
        Xsiv=YV(j)*R                          ! """"             of Vnode [m]

        IF (Xsip.gt.YL) THEN                  !
           Cp=ClearO/Clear                    !               .-------- ClearO  
        ELSE                                  ! ClearL ........
           Cp=ClearL/Clear                    !   ---------------------
        END IF                                !  -Ll          YL      y=0

        IF (Xsiv.gt.YL) THEN
           Cv=ClearO/Clear
        ELSE
           Cv=ClearL/Clear
        END IF
        DXOp=EXO+AYO*(Xsip-ZO)/Clear
        DYOp=EYO-AXO*(Xsip-ZO)/Clear
        DXOv=EXO+AYO*(Xsiv-ZO)/Clear
        DYOv=EYO-AXO*(Xsiv-ZO)/Clear

c      !....................!                   !.............................
        DO   I=1, NXT                           !
c      !....................!                   ! sweep on circ. direction.
          XJ=XP(I)
          IF (XJ.LT.0.0D0) XJ=XJ+PI2
          K=1
  331     IF ((XJ.LT.XlobeT(K)).AND.(XJ.GE.XlobeL(K))) THEN
           goto 332                          
          ELSE
           IF (K.GE.NLOBES) GOTO 332
           K=K+1
           GOTO 331
          END IF

  332     HP(I,J)=CP+(DXOp+CosXlobe(K))*COSXP(I)+
     +               (DYOp+SinXlobe(K))*SINXP(I) 
          HU(I,J)=CP+(DXOp+CosXlobe(K))*COSXU(I)+
     +               (DYOp+SinXlobe(K))*SINXU(I)
          HV(I,J)=CV+(DXOv+CosXlobe(K))*COSXP(I)+
     +               (DYOv+SinXlobe(K))*SINXP(I)
          HPN(I,J)=HP(I,J)
          HUN(I,J)=HU(I,J)
          HVN(I,J)=HV(I,J)

          IF ((HP(I,J).LE.ZERO).OR.(HU(I,J).LE.ZERO).OR.
     +        (HV(I,J).LE.ZERO)) Ilarge=1
c      !....................!                   
        END DO
c      !....................!    sweep in X-direction                  

c     !------------------------!              !..............................
      END DO
c     !------------------------!              !.........sweep in axial direction



      GOTO 99

C....................................................................      
 444    write (6,*) 'ERROR on calculation of film thickness for LOBES'
        write (6,*) 'J:',I,'Xj:',Xj,'K:',K
        write (6,*) 'PRG ABORTS execution : STOP'
        RETURN 

C....................................................................      
  99    CONTINUE
C       CHECK if film thickness <=0 for rigid surface bearing
c      !.......................................................!
        IF ((AC.EQ.ZERO).AND.(ILARGE.EQ.1)) THEN
        WRITE (6,100)
 100    FORMAT('$ Eccentricity: Exo or Eyo TOO LARGE',
     +         '  Film thickness: H is negative or zero',/,
     +         '$ PROGRAM ABORTS execution on SUB FILMH')
        STOP
        END IF
c....................................................................

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Hwear                                                       **
C **                                                                         **
C **  FILMH: Determines clearance & wear at Xp and Xu coordinates.           **
C **                                                                         **
C *****************************************************************************
c REFERENCE: Scharrer J. et al., 1990
c            The effects of Wear on the Rotordynamic Coefficients of a Hydros-
c            tatic Journal Bearing
c            1990 ASME-STLE Tribology Conference
c..............................................................................

      SUBROUTINE HWEAR(EXO, EYO, NXT, NYI)

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /TRIGS/ COSXP(MAXNXT), SINXP(MAXNXT), 
     +               COSXU(MAXNXT), SINXU(MAXNXT)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /SPLDATA/ Z(Nsl),Cl(Nsl),Bcl(Nsl),Ccl(Nsl),Dcl(Nsl), Nj

      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /WEAR/ EWX,EWY,EWEAR,BETA, IWEAR 
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL

      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION DXP, DXU, SUW, SUE, HP, HU, HV,
     +                 XP, XU, YP, YV, COSXP, SINXP, COSXU, SINXU
      DOUBLE PRECISION CLEAR, DIAM, LENGTH, LD, AR, HREC
      DOUBLE PRECISION Z,Cl,Bcl,Ccl,Dcl,
     +                 EWX,EWY,EWEAR,BETA, AXO, AYO, ZO,   
     +                 ClearO,ClearR,ClearL,YR,YL

      INTEGER ISYM, ICSTEP,BEARING
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      INTEGER NXT, NYI, I, J, NJ, Iwear,  JSTART
      DOUBLE PRECISION EXO, EYO, XSIP, XSIV, CP, CV, ZERO,
     +                 Rj, Rbp,Rbv,Rw,Cpp,Cuu,Cvv, DXOp,DXOv,
     +                 Nump,Numv,Den,Anglep,Anglev, DYOp, DYOv,
     +                 Theta,L1p,L1v,L2p,L2v, PI, Factor,Px  
C ----------------------------------------------------------------------------
C --  HWEAR code                                                            --
C ----------------------------------------------------------------------------
      Iwear=0
      EWX=0.0D0
      EWY=0.0D0
      BETA=0.0D0
      write (6,77)
 77   format (3X,74('='),/,3X,
     + 'WEAR option not available in this version of hydroflex',/,3X,
     + 'Contact Dr. Luis San Andres at Texas A&M for an upgrade',/,3X,
     +  74('='))

      CALL FILMH(EXO, EYO, NXT, NYI)

      END 
c
c--------------------------------------------------------------------
c Last revision 2/24/94 by Dr. Luis SanAndres
c for film thickness in preloaded pad bearings.
c--------------------------------------------------------------------