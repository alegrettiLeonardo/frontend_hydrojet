c #####        #  ######   #####          ######
c #    #       #  #          #            #
c #    #       #  #####      #            #####
c #####        #  #          #     ###    #
c #       #    #  #          #     ###    #
c #        ####   ######     #     ###    #

c  hydrojet.f  Copyright Luis San Andres / TexasA&MUniversity / 1995
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR III
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

C These routines calculate pressure and circ. flow within recesses
c of a hydrostatic bearing with angled injection.
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Pjet                                                        **
C **                                                                         **
C **  Pjet:  Find 1-D pressure field within recess of HJB                    **
C **                                                                         **
C *****************************************************************************


      SUBROUTINE PJET

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
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT), 
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT), SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /Pedge/ PEDRISE(MAXNPOCK,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECASP/ ASPEC(MAXNPOCK), /RECCOM/ L4R(MAXNPOCK)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPE, HRECD
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW 

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP

      COMMON /JET/ ANGLEJ, XJET, CJETnom, DPJET
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)   
      COMMON /URECJET/UREC(MAXNPOCK,MAXNPOCK+2)

      DOUBLE PRECISION U,V,P,RHOP,EMUP, XP,XU,YP, YV, PEDRISE,
     +                 HP, HU, HV, 
     +                 DXP, DXU, SUW, SUE, DYP, DYV, SVN, SVS,
     +                 PREC,TREC,QREC,QIN,QOUT,QFACTOR, ASPEC, L4R,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPE, HRECD,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 ANGLEJ, XJET, CJETnom, DPJET,
     +                 PRECdo, PRECup,PRJET,UREC

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION  HEDGE,ETA,RENR,ALFAR,kR, MACH2,
     +                  RHO, RHOR, EMUR, CKR, BTR, THR,
     +                  bd, bu, Udown, Uup, DPdown, DPup,
     +                  X,XL,XR,XRML,Xo,DXL,DXR, PRdo,PRup, PP,
     +                  Flow, Ysize, FAC, QRL, QRR, QSIDE,CJET

      INTEGER DEVICE, K, Kstart, Kend, Kstep, I, IJ,
     +        J, JL, JR, JP1, JJ, Jorif

C ----------------------------------------------------------------------------
C --  PJET code                                                            --
C ----------------------------------------------------------------------------
c ..................................................................!
c determine circumferential velocities within pocket using GLOBAL
c CONTINUITY RELATIONSHIP: U = [Fe - (Fn-Fs)]/(RHO [HREC+H] Yr ]
c where Yr=Ysize is the effective axial length of pocket
c ...................................................................!
c Angled injection coefficient as defined on initopst.f
c CJETnom=PI*CD*DORIF**2*DSIN(ANGLEJ*PI/180.0D0)/(2.0D0*AR*CLEAR(1+HRECD))      


c                                     ! ON POCKET:
      Kend=NPA                        ! Number of grid points axial
      IF (ISYM.EQ.0) THEN             ! ..........................
         Kstart=-NPA                  ! Account for symmetric 
         FAC=1.0D0                    ! and asymmetric
         Ysize=YV(Kend)-YV(-Kend-1)   ! bearing geometry
      ELSE                            !
         Kstart=1                     !
         FAC=2.0D0                    !
         Ysize=Fac*YV(Kend+1)         !
      END IF                          !........................

         

c    !.............................!
      DO I=1, NPOCKET                   
c    !.............................!
         JL=(I-1)*NXI+NLC          ! NXI=NLC+NPC-2
         JR=I*NXI+1                !
         XR=XP(JR)                 ! coordinates of recess
         XL=XP(JL)                 ! ..
         XRML=XR-XL                ! x-size of recess
         DXL=XJET*XRML             ! upstream size -> orifice
         DXR=(1.0D0-XJET)*XRML     ! orifice -> downstream size
         Xo=XL+DXL  ! = XR-DXR     ! location of orifice
         ASPE=ASPEC(I)*2.0D0       ! size of recess ASPE=XRML ????
         bu=ASPE*XJET              ! x-length of upstream recess
         bd=ASPE*(1.0-XJET)        ! x-length of downstream recess


        QRL=0.0D0                  ! circ. flow on left & right
        QRR=0.0D0                  ! of pockets
        JL=JL-1                    ! ==========

C        !.................!             !........................         
          DO K=1,Kend, 1                 ! circumferential flows
C        !.................!             !........................
          J=JL                           ! At upstream edge (left)
          JP1=J+1
          RHO=Rhop(J,K)*SUE(J)+Rhop(JP1,K)*SUW(J)
          QRL=QRL+FAC*DYP(K)*RHO*HU(J,K)*U(J,K) 
        IF (ISYM.EQ.0) THEN
          RHO=Rhop(J,-K)*SUE(J)+Rhop(JP1,-K)*SUW(J)
          QRL=QRL+DYP(-K)*RHO*HU(J,-K)*U(J,-K) 
        END IF

          J=JR  			 ! At downtream edge (right)
          JP1=J+1                        !
          IF ((J.EQ.NXT).AND.(IFULL.EQ.1)) JP1=2
          RHO=Rhop(J,K)*SUE(J)+Rhop(JP1,K)*SUW(J)
          QRR=QRR+FAC*DYP(K)*RHO*HU(J,K)*U(J,K)
        IF (ISYM.EQ.0) THEN
          RHO=Rhop(J,-K)*SUE(J)+Rhop(JP1,-K)*SUW(J)
          QRR=QRR+DYP(-K)*RHO*HU(J,-K)*U(J,-K) 
        END IF
C        !..................!
          END DO
C        !..................!            ! K=1,... NPA

         IJ=0
         QSIDE=0.0D0 
         JL=JL+1                   !=(I-1)*NXI+NLC
         K=Kend+1                  ! 
c       !...............!          !
         DO J= JL, JR              ! X-sweep along recess
c       !...............!          ! NPC-1=JR-JL
c         find axial flow contribution:
          RHO=Rhop(J,Kend)*SVN(K)+Rhop(J,K)*SVS(K)
          QSIDE=QSIDE+HV(J,K)*DXP(J)*V(J,K)*RHO*Fac
        IF (ISYM.EQ.0) THEN
          RHO=Rhop(J,-Kend)*SVN(-K)+Rhop(J,-K)*SVS(-K)
          QSIDE=QSIDE-HV(J,-K)*DXP(J)*V(J,-K)*RHO*Fac
        END IF

           IJ=IJ+1                 ! AT X-coordinate
           X=XP(J)                 !
c         !........................!........................
           IF (X.LE.Xo) THEN       ! XL < X <= Xo (upstream)
c         !........................! upstream of orifice
            Flow=(QRL-QSIDE)
            Jorif=J
c         !........................! ........................
           ELSE IF (X.GT.Xo) THEN  ! Xo < X < XR (downstream)
c         !........................! downstream of orifice
            Flow=(QRL-QSIDE+QREC(I))
c         !........................!........................
           END IF
c         !........................!.........................
c         estimate circumferential velocity within pocket
c###      CALL LOCPROPS(RhoR,EMUR,CKR,BTR,THR,PRjet(I,IJ),TREC(I))  !<== correct one
          CALL LOCPROPS(RhoR,EMUR,CKR,BTR,THR,PREC(I),TREC(I))
          Urec(I,IJ)=Flow/RhoR/Ysize/(HRECd+HU(J,1))
c       !...............!          !
         END DO
c       !...............!          !

         Jorif=Jorif+1             ! approx. location of recess orifice

c      !...........................!........ DOWNSTREAM EDGE of RECESS
c  ::: DPdown: estimate pressure rise due to step/hydrodynamics
          IF (Urec(I,IJ).GT.0.0D0) THEN
             CALL LOCPROPS(RHOR,EMUR,CKR,BTR,THR,PREC(I),TREC(I))
             HEDGE=HRECD+HU(Jorif,1)              
             ETA=1.0D0/HEDGE                 !
             RENR=RENC*HEDGE*RHOR/EMUR       ! circumferential Re# at recess.
             kR=(RENR**0.681)/7.752963D+00   ! shear factor at recess.
             ALFAR=kR*EMUR*bd/(HEDGE*HEDGE)
             DPdown=ALFAR*(MSPEED-UREC(I,IJ))    !IJ=NPC
          ELSE
             DPdown=0.0D0
          END IF
c      !...........................!...............................!

c      !...........................!............UPSTREAM EDGE of RECESS
c  ::: DPup: estimate pressure rise/drop due to step/hydrodynamics
       IF (Urec(I,1).LT.0.0D0) THEN
             CALL LOCPROPS(RHOR,EMUR,CKR,BTR,THR,PREC(I),TREC(I))
             HEDGE=HRECD+HU(Jorif,1)              
             ETA=1.0D0/HEDGE                 !
             RENR=RENC*HEDGE*RHOR/EMUR       ! circumferential Re# at recess.
             kR=(RENR**0.681)/7.752963D+00   ! shear factor at recess.
             ALFAR=kR*EMUR*bu/(HEDGE*HEDGE)
             DPup=ALFAR*(MSPEED-UREC(I,1))
       ELSE
             DPup=0.0D0
       END IF

             DPup=0.0D0            !## based on exp. tests for turb. flows

c      !...........................!...............................!
         CJET=CJETnom*(1.0D0+HRECD)/(HU(Jorif,1)+HRECD)
         PRECdo(I)=PREC(I)-CJET*(1.0D0-PREC(I))*(1.0D0-XJET)+DPdown
         IF (PRECdo(I).LT.PCAV) PRECdo(I)=PCAV
         PRECup(I)=PREC(I)+CJET*(1.0D0-PREC(I))*XJET        -DPup
         IF (PRECup(I).LT.PCAV) PRECup(I)=PCAV
c      !...........................................................!



           IJ=0                    ! ESTIMATE pressures within
c       !...............!          ! recess: LINEAR Prise/drop
         DO J= JL, JR              ! X-sweep along recess
c       !...............!          ! NPC-1=JR-JL
           IJ=IJ+1
           X=XP(J)                 !
c         !........................!........................
           IF (X.LE.Xo) THEN       ! XL < X <= Xo (upstream)
c         !........................! upstream of orifice
             PP = (PRECup(I)*(Xo-X)+PREC(I)*(X-XL))/DXL
c         !........................! ........................
           ELSE IF (X.GT.Xo) THEN  ! Xo < X < XR (downstream)
c         !........................! downstream of orifice
             PP = ( PREC(I)*(XR-X)+PRECdo(I)*(X-Xo))/DXR
c         !........................!........................
           END IF
c         !........................!.........................
           PRjet(I,IJ)=PP
c       !...............!          !
         END DO
c       !...............!          !

         DO K=Kstart,NPA
           PEDRISE(I,K)=DPdown
         END DO

C    !......................!            !
      END DO                             ! NEXT Recess, I=1,2,...NPOCKET
C    !......................!            !


 171  FORMAT (1X,'PR:',7(G10.5,1X))
 172  FORMAT (1X,'UR:',7(G10.4,1X))
 173  FORMAT (1X,'Qr:',5(G10.4,1X))
	END

c============================================================================



C *****************************************************************************
C **                                                                         **
C **  Subroutine Pjet1                                                       **
C **                                                                         **
C **  Pjet1:  Find 1-D perturbed pressure field within recess of HJB         **
C **                                                                         **
C *****************************************************************************


      SUBROUTINE PJET1

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT), 
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECASP/ ASPEC(MAXNPOCK), /RECCOM/ L4R(MAXNPOCK)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPE, HRECD
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW 

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /VERB/ SVERB, DVERB, BEEP
      COMMON /HJBSYM/ ISYM, ICSTEP

      COMMON /JET/ ANGLEJ, XJET, CJETnom, DPJET
      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK), 
     +                QREC1(MAXNPOCK), QIN1, QOUT1
      COMMON /RECJET1/ PREC1do(MAXNPOCK),PREC1up(MAXNPOCK),
     +                 PR1jet(MAXNPOCK,MAXNPOCK+2)   

      DOUBLE PRECISION U,V,P,RHOP,EMUP, XP,XU,YP, YV, PEDRISE,
     +                 HP, HU, HV, 
     +                 PREC,TREC,QREC,QIN,QOUT,QFACTOR, ASPEC, L4R,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPE, HRECD,
     +                 PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 ANGLEJ, XJET, CJETnom, DPJET
      DOUBLE COMPLEX PREC1,TREC1,QREC1,QIN1,QOUT1,
     +               PREC1do, PREC1up,PR1JET

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        SVERB, DVERB, BEEP, ISYM, ICSTEP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION  HEDGE,ETA,RENR,ALFAR,kR, MACH2,
     +                  RHOR, EMUR, CKR, BTR, THR, CJET,
     +                  bd, bu, Udown, Uup, DPdown, DPup,
     +                  X,XL,XR,XRML,Xo,DXL,DXR, PRdo,PRup
      DOUBLE COMPLEX PP

      INTEGER DEVICE, I, IJ, J, JL, JR

C ----------------------------------------------------------------------------
C --  PJET1 code                                                            --
C ----------------------------------------------------------------------------

c    !.............................!
      DO I=1, NPOCKET                   
c    !.............................!
         JL=(I-1)*NXI+NLC          ! NXI=NLC+NPC-2
         JR=I*NXI+1                !
         XR=XP(JR)                 ! coordinates of recess
         XL=XP(JL)                 ! ..
         XRML=XR-XL                ! x-size of recess
         DXL=XJET*XRML             ! upstream size -> orifice
         DXR=(1.0D0-XJET)*XRML     ! orifice -> downstream size
         Xo=XL+DXL  ! = XR-DXR     ! location of orifice
         ASPE=ASPEC(I)*2.0D0       ! size of recess ASPE=XRML ????
         bu=ASPE*XJET              ! x-length of upstream recess
         bd=ASPE*(1.0-XJET)        ! x-length of downstream recess
         CJET=CJETnom  !*(1.0+HRECD)/(HU+HRECD)
         PREC1do(I)=PREC1(I)*(1.0D0+CJET*(1.0D0-XJET))
         PREC1up(I)=PREC1(I)*(1.0D0-CJET*XJET)
         IJ=0
c       !...............!          !
         DO J= JL, JR              ! X-sweep along recess
c       !...............!          ! NPC-1=JR-JL
           IJ=IJ+1
           X=XP(J)                 !
c         !........................!........................
           IF (X.LE.Xo) THEN       ! XL < X <= Xo (upstream)
c         !........................! upstream of orifice
             PP = (PREC1up(I)*(Xo-X)+PREC1(I)*(X-XL))/DXL
c         !........................! ........................
           ELSE IF (X.GT.Xo) THEN  ! Xo < X < XR (downstream)
c         !........................! downstream of orifice
             PP = ( PREC1(I)*(XR-X)+PREC1do(I)*(X-Xo))/DXR
c         !........................!........................
           END IF
c         !........................!.........................
           PR1jet(I,IJ)=PP
c       !...............!          !
         END DO
c       !...............!          !

c    !.............................!
      END DO 
c    !.............................!


	END

c=============================================================================
c hydrojet.f program: Luis San Andres 10/12/95
c=============================================================================

