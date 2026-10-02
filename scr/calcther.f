c===============================================================================
c 10/12/95 modified recess torque equation for better accuracy on
c          angled jet injections
c probably more work is needed for first order solution (FLOWT1 *TORRECJ)
c===============================================================================
c 6/28/94 modified for bearing radial heat flow
c
c         NOTE: RECESS ENERGY EQUATION DOES NOT ACCOUNT FOR
c               HEAT FLOW to bearing or journal surfaces, i.e.
c               it is insulated.
c .............................................................................
c 5/31/94 modified to include isothermal stator and journal
c         different from Tsupply
c         removed arrays of bearing and shaft temperatures since these
c         are constant.
c
c 5/11/94 modified to save maximum temperature difference
c         DELTAT=Tnew-Told for checking on SLAND
c                          as convergence parameter
c
C 4/26/94 modified GROOVE subroutine for improved evaluation of
c         leading edge temperature RMST renamed as TDIFF
c         other was misleading
c         original file saved as calcthert.f426
c 4/27/94 improved calculations in CALCTM for JJ
c
c ==============
c THD MODELS are:
c ==============
C      IF(ISOTH.EQ. 1) THERM='Isothermal fluid film(T=Tin=Constant)'         93
C      IF(ISOTH.EQ. 2) THERM='Adiabatic journal & Isothermal bearing (Tb)'  5/94
C      IF(ISOTH.EQ. 3) THERM='Isothermal journal (Tj) & Adiabatic bearing'  5/94
C      IF(ISOTH.EQ. 4) THERM='Adiabatic journal & Bearing radial heat'      6/94
C      IF(ISOTH.EQ. 5) THERM='Isothermal journal(Tj) & Bearing radial heat' 6/94
C      IF(ISOTH.EQ. 0) THERM='Isothermal journal(Tj) & bearing: (Tb)'       5/94
C      IF(ISOTH.EQ.-1) THERM='Adiabatic bounding Surfaces (Qb=Qj=0)'         93
c ==============
c
C  ####     ##    #        ####    #####  #    #  ######  #####    #####
C #    #   #  #   #       #    #     #    #    #  #       #    #     #
C #       #    #  #       #          #    ######  #####   #    #     #
C #       ######  #       #          #    #    #  #       #####      #     ###
C #    #  #    #  #       #    #     #    #    #  #       #   #      #     ###
C  ####   #    #  ######   ####      #    #    #  ######  #    #     #     ###
C
c hydroflex.f . Copyright Luis SanAndres / TexasA&MUniversity / 1994
c
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR II
c Technical monitor: Mr. James Walker, NASA Lewis Research Center

c========================================================================
c  4/4/94: modified CALCT1 for compliance of bearing surface in calculation
c          of first order solution for temperature.
c          AC/(1+ixETAC) *GTH * P1
c          modified FLOWT1 for boundary flows with compliance
c
c  On 11/05/93, LEMDA is input by user on inputpad routine
c  of calcmesht.f
c
c
c 11/05/93: I have modified the cavitation algorithm so as  to effectively
c           set all energy transfer and dissipation on cavitation zone, i.e.
c           it is an insulated zone
c
c========================================================================

C *****************************************************************************
C **                                                                         **
C **  Subroutine Calct                                                       **
C **                                                                         **
C **  CALCT:  Solves Energy equation                                         **
C **          on film lands with TW=T(1,K), TE=T(NXT,K)                      **
C **          Ap Tp = Ae Te + Aw Tw + As Ts+ An Tn + St                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCT(PS,REC,DIR,KV)

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
      COMMON /IOPROP/ RHS,EMS,RHA,EMA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HB(MAXNXT,-MAXNYI:MAXNYI),
     +               HJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC ,Ec
      COMMON /TISOBJD/ TBJ,TBS
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /DPUV/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /TDIFMAX/ DELTAT

      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /THERMID/ ISOTH

      DOUBLE PRECISION DXP, DXU, HP, HU, HV, SUW, SUE,
     +                 RHS,EMS,RHA,EMA,RHOle,EMUle,RHOri,EMUri,
     +                 CPS,THS,BETAKS, DYP,DYV,SVN,SVS,
     +                 U, V, P, RHP,EMP, DU, DV, DVV, A, B, C, D,
     +                 REP, REY, MSPEED, SPEED,SWIRL,PCAV,ALFU,BETU,
     +                 ALFP, DYVK, DYPK, SVNK, SVSK, PRATIO, CORIF,
     +                 SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 CK, BETAK, THC, T ,HB,HJ,TBJ,TBS,
     +                 ALFT, UC, TC, Ec ,DIR, DELTAT,
     +                 TBOUT, THERMALK, ROUTER, HKB

      INTEGER JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL ,ISOTH,
     +        NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        K, KM1, KP1 ,KV,IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER J, JJ, JJJ, JJMAX ,JP1
      DOUBLE PRECISION DPU, DPV, FWP, FEP, FSP, FNP,
     +                 RHOW,RHOE,RHOS,RHON,
     +                 SVNKM1,SVSKM1, BETAT, DELTA,
     +                 APT, AWT, AET, AST, ANT, SP1, SP2, STC,
     +                 PW,PE,TW,TE, SANB, UP,VP,UW,
     +                 ST1,ST2,ST3,KXP,KXJ,KXB,FNKX ,REC,REYCK ,ZERO,
     +                 PN(MAXNXT),PS(MAXNXT), CH,GAMA, HBEQUIV

C ----------------------------------------------------------------------------
C --  CALCT code                                                            --
C ----------------------------------------------------------------------------
C  Pw,Pe,Pn and Ps are interpolated pressures at faces of P control volume
C  T and P share the same contro volume
C  Ps is extracted on Sland of calcsolnt.f
C
C  Please note crudness of cavitation model, i.e. essentially we neglect
C  all mechanical energy being disposed as heat, i.e. cavitation region
C  should be a zone of constant temperature.
C ----------------------------------------------------------------------------
C     TBJ=TSHAFT/TC                            ! dimensionless bearing and
C     TBS=TSTATOR/TC                           ! journal temperatures.
c					       !.............def in initopst.f
      ZERO=0.0D0			       !
      BETAT=1.0D0-ALFT                         ! relaxation parameters

c					       !...........................!
      JJJ=2-JPMIN                              ! indexes for circ. sweeps
      JJMAX=JPMAX+JJJ                          !
      J=JPSTART                                !
c					       !...........................!

      RHOW=Rhp(j,k)*Sue(j)+Rhp(jpmin,k)*Suw(j) ! AT WEST most T node
      DPU=HU(J,K)*DYPK*RHOW                    !
      UW=U(J,K)				       !
      FWP=DPU*UW                               ! Flow on west side of Pcv

      Pw=P(J,K)*SUE(J)+P(JPMIN,K)*SUW(J)       ! Pw at J=JPMIN
      TW=T(JPSTART,K)	                       ! For TDMA Boundary Values
      TE=T(JPSTOP,K)	                       ! at Two Ends--from FLOWT of
                                               ! energy balance at recesses

      SVNKM1=SVN(K)                            !
      SVSKM1=SVS(K)                            !

      CH=2.D0*REY*(EMS*CPS/THS)**0.66667D0     ! Parameter for Heat coefs.
C ......................                       !
      DO J=JPMIN, JPMAX                        !----------------------------
C ......................                       ! Sweep from left-right on P cvs
          JJ=J+JJJ                             !
          JP1=J+1                              !
C                                              ! REC = REY/EC = Rep* /EC
          REYCK=REC*CK(J,K)                    ! Rep* x specific heat / EC
c                                              !
          AWT=DMAX1(+FWP,ZERO)*REYCK           ! Awt coeff.
                                               !
          Rhoe=Rhp(j,k)*Sue(j)+Rhp(jp1,k)*Suw(j)
          DPU=HU(J,K)*DYPK*RHOE                ! East U flow on Pcv
          FEP=DPU*U(J,K)                       !
          AET=DMAX1(-FEP,ZERO)*REYCK           !

          Rhon=Rhp(j,kp1)*Svsk+Rhp(j,k)*Svnk   !
          DPV=HV(J,KP1)*DXP(J)*RHON            !
          FNP=Dir*DPV*V(J,KP1)                 ! North V flow on Pcv
          ANT=DMAX1(-FNP,ZERO)*REYCK           !

          Rhos=Rhp(j,k)*Svskm1+Rhp(j,km1)*Svnkm1
          DPV=HV(J,KV)*DXP(J)*RHOS             !
          FSP=Dir*DPV*V(J,KV)                  ! South V flow on Pcv
          AST=DMAX1(+FSP,ZERO)*REYCK           !

          SANB=AET+AWT+AST+ANT                 !--------------------

          UP=(UW+U(J,K))/2.0D0                 ! Velocities at center
          VP=(V(J,K)+V(J,KP1))/2.0D0           ! of T-CV

          Pe=P(J,K)*SUE(J)+P(JP1,K)*SUW(J)     ! P on east/north faces
          Pn(J)=P(J,KP1)*SVS(KP1)+P(J,K)*SVN(KP1)
c      !.......................!               !............................!
          HB(J,K)=0.D0                         ! (1) Adiabatic or Isothermal
          HJ(J,K)=0.D0                         ! Heat transfer coefs.=0
c      !.......................!               !............................!
c                                              !
C.............................                 ! In full fluid film region
      IF (P(J,K).GE.PCAV) THEN                 ! Note: even P=Pcav, like the
C.............................                 ! concentric case, shear exits
c                                              !
          KXP=FNKX(UP,VP,HP(J,K),              ! Shear stress on journal
     +    SPEED,REP,RHP(J,K),EMP(J,K),KXJ,KXB) ! Dr. San Andres' formula
                                               !
          GAMA=(UP*KXB-(UP-SPEED)*KXJ)/HP(J,K)/4.0D0
c                                              !
          IF (KXP.EQ.(12.0D0*EMP(J,K))) THEN   !
             GAMA=SPEED*EMP(J,K)/HP(J,K)       ! Shear for laminar flow
          END IF                               !
C..............................................!
      HBequiv=0.0D0                            ! Heat transfer coefficients
c    !....................!                    ! ----------------------------
      IF (ISOTH.EQ.0) THEN                     ! (0) Tb & Tj constants: isothermal
c    !....................!                    ! bearing & journal surfaces
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HJ(J,K)=HB(J,K)*KXJ/KXB              !
          HBequiv=HB(J,K)                      !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.2) THEN                ! (2) Qj=0, Tb=Tc: isothermal
c    !....................!                    ! bearing and adiabatic
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   ! journal surface
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HBequiv=HB(J,K)                      !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.3) THEN                ! (3) Qb=0, Tj=Tc: isothermal
c    !....................!                    ! journal surface and adiabatic
          HJ(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   ! bearing
     +    *THC(J,K))**0.333333*KXJ/EMP(J,K)/HP(J ,K)/CH
          HBequiv=0.0D0                        !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.4) THEN                ! (4) Qj=0, adiabatic journal
c    !....................!                    ! & bearing radial heat flow
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HBequiv=HB(J,K)/(1.0D0+HB(J,K)/HKB)  !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.5) THEN                ! (5) Isothermal journal
c    !....................!                    ! & bearing radial heat flow
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HJ(J,K)=HB(J,K)*KXJ/KXB              !
          HBequiv=HB(J,K)/(1.0D0+HB(J,K)/HKB)  !
c    !....................!                    !
      ENDIF                                    !
c    !....................!....................!THERMAL CASES



c     FOR CASES 4 and 5: HKB represents the dimensionless heat transfer
c     coefficient due to conductivity on solid bearing, and equal to
c     HKB = [K/radius]/[REY x (CP EMU/CLEAR]* x ln(Router/Radius) ]
c     and where K is the bearing material conductivity coefficient  LAS 6/28/94
c         TBS = TBout/Ts, i.e. outer bearing temperature.
C..............................................!..............................c
C NOTE: These convection heat transfer coefficients are made dimensionless
C       according to         H [Wats/m2 degK]
c                      H =  -------------------------
c                           REY x [ Cp* EMU*/ CLEAR*]
c
c where REY = REp* = [ RHO/EMU x U x Clear**2. / R]*
c
c in general for turbulent flows:
c  H = 0.5 RHO x CP x  VELOCITY x FRICTION FAC / (Prandtl #)**(2/3)
c
c  where Prandtl # = PR = CP x EMU / K ; K: fluid thermal conductivity.
c
C..............................................!..............................c
                                               !
      SP1=REc*(HBequiv+HJ(J,K))*DXP(J)*DYPK    ! Surface Heat transfer term

      ST1=(KXP*(UP*UP+VP*VP)-UP*MSPEED*KXJ)    ! Friction disipation
     +   /HP(J,K)+GAMA*SPEED                   !

      ST3=HP(J,K)*MSPEED*(PE-PW)*DYPK          ! h Speed/2 dp/dx

      ST2=REc*(HBequiv*TBS+HJ(J,K)*TBJ)        ! Heat transfer to walls
c                                              ! Constant surface temperatures

      SP2=-BETAK(J,K)*HP(J,K)*(UP*(PE-PW)      ! Pressure extrution
     +   *DYPK+VP*(PN(J)-PS(J))*Dir*DXP(J))    ! mechanical work

      STC=(ST1+ST2)*DXP(J)*DYPK+ST3 +          ! Constant term on RHS
     +          DMAX1(-SP2,ZERO)*T(J,K)        ! Negative slope of source term

      APT=(SANB+SP1+DMAX1(SP2,ZERO))/ALFT      ! Main coefficient of Tp

C.............................                 !.....................
      ELSE                                     ! In cavitation region
C.............................                 ! Viscous shear and
C     set RHS of energy equation as zero       ! extrusion terms = 0
c                                              ! i.e. model as
      STC=0.0D0                                ! insulated region
      APT=SANB/ALFT                            ! (no heat flow)
C.............................                 !------------------------
      END IF                                   ! End of CAVITATION model
C.............................P>Pcav           !        APPROXIMATE only
C                                              !------------------------

      C(JJ)=-AWT                               ! coeffs. for TDMA soln.
      A(JJ)=-AET                               !
      B(JJ)= APT
      D(JJ)=STC+AST*T(J,KM1)+ANT*T(J,KP1)      ! Constant term on RHS
     +         +BETAT*B(JJ)*T(J,K)             !
C..............................................!------------------- T-EQ.
          FWP=FEP                              !
          UW=U(J,K)                            ! NEXT j STEP--Uw|j=Ue|j-1
	      PW=PE                                !
	      PS(J)=PN(J)                          ! Reserved for next K-row
C
C......................                        !
      END DO                                   ! end of j loop
C......................                        ! ::::::::::::::::::::::::::::
C                                              ! ----------------------------
      CALL TDMA(TW,TE,JJMAX)                   ! Solve by TDMA algorithm-T-EQ.
C ............................................ ! ----------------------------
C                                              ! ::::::::::::::::::::::::::::
  100 DO J=JPMIN, JPMAX                        ! CORRECT:
          JJ=J+JJJ                             ! Temperature field
          DELTA=ALFT*(B(JJ)-T(J,K))            ! using under-relaxation
          T(J,K)=T(J,K)+DELTA                  !
          DELTA=DABS(DELTA)                    !
          IF (DELTA.GT.DELTAT) DELTAT=DELTA    ! Save Max. Tdifference
      END DO                                   !
C ............................................ ! ----------------------------
C                                              !
      IF (IFULL.EQ.1.AND.INTER.EQ.0) THEN      ! 360 deg bearing
          T(JPSTOP, K)=T(JPMIN, K)             ! On extended lands
      END IF                                   !
c .............................................!.............................

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Calctm                                                      **
C **                                                                         **
C **  CALCTM: Solves Energy equation on film lands between recesses          **
C **          with TW=TRECL,TE=TRECR. Note Tleft & Tright at recess          **
C **          edges are unknowns and obtained from the solution:             **
C **          Ap Tp = Ae Te + Aw Tw + As Ts +An Tn + St                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCTM(PS,REC,DIR,KV)

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
      COMMON /IOPROP/ RHS,EMS,RHA,EMA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BETAKS
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HB(MAXNXT,-MAXNYI:MAXNYI),
     +               HJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC ,Ec
      COMMON /TISOBJD/ TBJ,TBS
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB
      COMMON /DPUV/ DU(MAXNXTP2), DV(MAXNXTP2), DVV(MAXNXTP2)
      COMMON /TDMA0/ A(MAXNXTP2), B(MAXNXTP2), C(MAXNXTP2), D(MAXNXTP2)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /KVALA/ DYVK, DYPK, SVNK, SVSK
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS, MPEPS, PREPS, MMP, SFLOW
      COMMON /TRECES/ TRECL, TRECR
      COMMON /TDIFMAX/ DELTAT

      COMMON /THERMID/ ISOTH
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /KVALB/ K, KM1, KP1


      DOUBLE PRECISION DXP, DXU, HP, HU, HV, SUW, SUE,
     +                 RHS,EMS,RHA,EMA,RHOle,EMUle,RHOri,EMUri,
     +                 CPS,THS,BETAKS, DYP,DYV,SVN,SVS,
     +                 U, V, P, RHP,EMP,DU, DV, DVV,A, B, C, D,
     +                 REP, REY, MSPEED, SPEED,SWIRL,PCAV,ALFU,BETU,
     +                 ALFP, DYVK, DYPK, SVNK, SVSK, PRATIO, CORIF,
     +                 SMASS, MPEPS, PREPS, MMP, SFLOW,
     +                 CK, BETAK, THC, T ,HB,HJ,TBJ,TBS,
     +                 ALFT, UC, TC, Ec ,TRECL,TRECR, DIR, DELTAT,
     +                 TBOUT, THERMALK, ROUTER, HKB

      INTEGER JPMIN, JPMAX, JPSTART, JPSTOP,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL ,ISOTH,
     +        NPOCKET, NLC, NPC,
     +        NLA, NPA, NPAP1, NXI, NYI, NXT,
     +        K, KM1, KP1 ,KV,IFULL

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER J, JJ, JJJ, JJMAX ,JP1
      DOUBLE PRECISION DPU, DPV, FWP, FEP, FSP, FNP, FWP0,
     +                 RHOW,RHOE,RHOS,RHON,
     +                 SVNKM1,SVSKM1, BETAT, DELTA,
     +                 APT, AWT, AET, AST, ANT, SP1, SP2, STC,
     +                 PW,PE,TW,TE, SANB, UP,VP,UW,
     +                 ST1,ST2,ST3,KXP,KXJ,KXB,FNKX ,REC,REYCK ,ZERO,
     +                 PN(MAXNXT),PS(MAXNXT), CH,GAMA, HBequiv

C ----------------------------------------------------------------------------
C --  CALCTM code                                                           --
C ----------------------------------------------------------------------------
      ZERO=0.0D0
      BETAT=1.0D0-ALFT                         ! relaxation parameters

      JJJ=3-JPMIN                              ! size of sweep counters
      JJMAX=JPMAX+JJJ                          ! and max. value
      J=JPSTART                                !

      RHOW=Rhp(j,k)*Sue(j)+Rhp(jpmin,k)*Suw(j)
      DPU=HU(J,K)*DYPK*RHOW                    !
      UW=U(J,K)
      FWP=DPU*UW                               ! Flow on west side of Pcv
      FWP0=FWP                                 ! Reserved for j=jpstart

      TW=TRECL                                 ! Recess temperatures are
      TE=TRECR                                 ! calculated in Sub. FLOWT
      Pw=P(J,K)*SUE(J)+P(JPMIN,K)*SUW(J)       ! Pw at J=JPMIN
                                               ! energy balance at recesses

      SVNKM1=SVN(K)                            !
      SVSKM1=SVS(K)                            !

      CH=2.D0*REY*(EMS*CPS/THS)**0.66667       !  Parameter for Heat coefs.
C ......................                       !
      DO J=JPMIN, JPMAX                        !----------------------------
C ......................                       ! Sweep from left-right on T cvs
          JJ=J+JJJ                             !----------------------------
          JP1=J+1                              !
          REYCK=REC*CK(J,K)                    ! Reynolds No. * specific heat
          AWT=DMAX1(+FWP,ZERO)*REYCK           ! Awt coeff.
                                               !
          Rhoe=Rhp(j,k)*Sue(j)+Rhp(jp1,k)*Suw(j)
          DPU=HU(J,K)*DYPK*RHOE                ! East U flow on Pcv
          FEP=DPU*U(J,K)                       !
          AET=DMAX1(-FEP,ZERO)*REYCK           !

          Rhon=Rhp(j,kp1)*Svsk+Rhp(j,k)*Svnk   !
          DPV=HV(J,KP1)*DXP(J)*RHON            !
          FNP=Dir*DPV*V(J,KP1)                 ! North V flow on Pcv
          ANT=DMAX1(-FNP,ZERO)*REYCK           !

          Rhos=Rhp(j,k)*Svskm1+Rhp(j,km1)*Svnkm1
          DPV=HV(J,KV)*DXP(J)*RHOS             !
          FSP=Dir*DPV*V(J,KV)                  ! South V flow on Pcv
          AST=DMAX1(+FSP,ZERO)*REYCK           !

          SANB=AET+AWT+AST+ANT                 !-------------------- T-EQ.

          UP=(UW+U(J,K))/2.0D0                 ! U,V velocities at center of
          VP=(V(J,K)+V(J,KP1))/2.0D0           ! T-CV
c                                              !
          Pe=P(J,K)*SUE(J)+P(JP1,K)*SUW(J)     ! Pe,Pn on faces of T-CV
          Pn(J)=P(J,KP1)*SVS(KP1)+P(J,K)*SVN(KP1)

c      !.......................!               !............................!
          HB(J,K)=0.D0                         ! (1) Adiabatic or Isothermal
          HJ(J,K)=0.D0                         ! Heat transfer coefs.=0
c      !.......................!               !............................!

c                                              !..............................!
C.............................                 ! In full fluid film region
      IF (P(J,K).GE.PCAV) THEN                 ! Note: even P=Pcav, like the
C.............................                 ! concentric case, shear exits
c                                              !
          KXP=FNKX(UP,VP,HP(J,K),              ! Shear stress on journal
     +    SPEED,REP,RHP(J,K),EMP(J,K),KXJ,KXB) ! Dr. San Andres' formula
c                                              !
          GAMA=(UP*KXB-(UP-SPEED)*KXJ)/HP(J,K)/4.0D0
c                                              !
          IF (KXP.EQ.(12.0D0*EMP(J,K))) THEN   !
             GAMA=SPEED*EMP(J,K)/HP(J,K)       ! Shear for laminar flow
          END IF                               !
C..............................................! -----------------------------
      HBequiv=0.0D0                            ! Heat transfer coefficients
c    !....................!                    ! -----------------------------
      IF (ISOTH.EQ.0) THEN                     ! (2) Tb;Tj=Tc: isothermal
c    !....................!                    ! bearing & journal surfaces
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J ,K)/CH
          HJ(J,K)=HB(J,K)*KXJ/KXB              !
          HBequiv=HB(J,K)                      !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.2) THEN                ! (3) Qj=0, Tb=Tc: isothermal
c    !....................!                    ! bearing & adiabatic
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   ! journal surface
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J ,K)/CH
          HBequiv=HB(J,K)                      !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.3) THEN                ! (4) Qb=0, Tj=Tc: isothermal
c    !....................!                    ! journal surface
          HJ(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   ! and adiabatic bearing surface
     +    *THC(J,K))**0.333333*KXJ/EMP(J,K)/HP(J ,K)/CH
          HBequiv=0.0D0                        !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.4) THEN                ! (4) Qj=0, adiabatic journal
c    !....................!                    ! & bearing radial heat flow
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HBequiv=HB(J,K)/(1.0D0+HB(J,K)/HKB)  !
c    !....................!                    !
      ELSE IF (ISOTH.EQ.5) THEN                ! (5) Isothermal journal
c    !....................!                    ! & bearing radial heat flow
          HB(J,K)=(EMP(J,K)*CK(J,K)*THC(J,K)   !
     +    *THC(J,K))**0.333333*KXB/EMP(J,K)/HP(J,K)/CH
          HJ(J,K)=HB(J,K)*KXJ/KXB              !
          HBequiv=HB(J,K)/(1.0D0+HB(J,K)/HKB)  !
c    !....................!                    !
      END IF                                   !
c    !....................!....................!THERMAL CASES
c     FOR CASES 4 and 5: HKB represents the dimensionless heat transfer
c     coefficient due to conductivity on solid bearing, and equal to
c     HKB = [K/radius]/[REY x (CP EMU/CLEAR]* x ln(Router/Radius) ]
c     and where K is the bearing material conductivity coefficient  LAS 6/28/94
c         TBS = TBout/Ts, i.e. outer bearing temperature.
C..............................................!..............................c
                                               !
      SP1=REc*(HBequiv+HJ(J,K))*DXP(J)*DYPK    ! Surface heat transfer

      ST1=(KXP*(UP*UP+VP*VP)-UP*MSPEED*KXJ)    ! Friction disipation
     +   /HP(J,K)+GAMA*SPEED                   !

      ST2=REc*(HBequiv*TBS+HJ(J,K)*TBJ)        ! for constant surface T's

      SP2=-BETAK(J,K)*HP(J,K)*(UP*(PE-PW)      ! Pressure extrution
     +   *DYPK+VP*(PN(J)-PS(J))*Dir*DXP(J))    ! term

      ST3=HP(J,K)*MSPEED*(PE-PW)*DYPK          ! h Speed/2  dP/dx

      STC=(ST1+ST2)*DXP(J)*DYPK+ST3            ! Constant term on RHS +
     +   +DMAX1(-SP2,ZERO)*T(J,K)              ! Negative slope of source term

      APT=(SANB+SP1+DMAX1(SP2,ZERO))/ALFT      ! Main coefficient

C.............................                 !.......................
      ELSE                                     ! In cavitation region
C.............................                 ! model as insulated region
      STC=0.0D0                                ! i.e. no heat flow
      APT=SANB/ALFT			       !
C.............................                 !--------------------------
      END IF                                   ! End of CAVITATION model
C.............................                 ! APPROXIMATE ONLY
c                                              !--------------------------

      C(JJ)=-AWT                               ! coeffs. for TDMA soln.
      A(JJ)=-AET                               !
      B(JJ)= APT
      D(JJ)=STC+AST*T(J,KM1)+ANT*T(J,KP1)      ! Constant term on RHS
     +         +BETAT*B(JJ)*T(J,K)             !
C..............................................!------------------- T-EQ.
          FWP=FEP                              !
          UW=U(J,K)                            ! NEXT j STEP--Uw|j=Ue|j-1
	  PW=PE                                !
	  PS(J)=PN(J)                          ! Reserved for next K-row
C
C......................                        !
      END DO                                   ! end of j loop
C ......................                       ! ::::::::::::::::::::::::::::

c Recess at right side of land                 !
      JJ=JJMAX+1                               ! Upwind scheme to determine
      AET=DMAX1(-FEP,ZERO)*REYCK               ! right edge temp.
      AWT=DMAX1( FEP,ZERO)*REYCK               ! assume FWP=FEP at jpstop
      A(JJ)=-AET                               ! where only advection terms
      C(JJ)=-AWT                               ! are considered to determine
      B(JJ)=(AET+AWT)/ALFT                     ! the recess edge temperatures
      D(JJ)=BETAT*B(JJ)*T(JPSTOP,K)            !

c Recess at left side of land                  !............................
      JJ=2                                     !
      AET=DMAX1(-FWP0,ZERO)*REC*CK(JPMIN,K)    ! left edge temp.
      AWT=DMAX1( FWP0,ZERO)*REC*CK(JPMIN,K)    !
      A(JJ)=-AET                               ! Therefore, T(1,K) & T(NXT,K)
      C(JJ)=-AWT                               ! are not B.C.'s, but calculated
      B(JJ)=(AET+AWT)/ALFT                     ! from the TDMA with JJMAX+1
      D(JJ)=BETAT*B(JJ)*T(JPSTART,K)           !
C ............................................ ! ----------------------------
C                                              !
      CALL TDMA(TW,TE,JJMAX+1)                 ! Solve by TDMA algorithm-T-EQ.
C ......................                       ! ::::::::::::::::::::::::::::

  100 DO J=JPSTART, JPSTOP                     ! Note here J=JPSTART, JPSTOP
          JJ=J+JJJ                             ! Correct Tfield
          DELTA=ALFT*(B(JJ)-T(J,K))            ! using under-relaxation
          T(J,K)=T(J,K)+DELTA                  !
          DELTA=DABS(DELTA)                    !
          IF (DELTA.GE.DELTAT) DELTAT=DELTA    !save Maximum Tdifference
      END DO                                   !
C ............................................ ! ----------------------------

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Flowt                                                       **
C **                                                                         **
C **  FLOWS:  Calculates inlet and outlet flow from hydrostatic bearing.     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FLOWT(REC,Inflow)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT), SUE(MAXNXT)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC, EC
      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /QSIDE0/ QSIDER(MAXNPOCK)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /QIOPAD/ QLEAD(MAXNPAD), QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV1, MODEL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      COMMON /PADK/ KPAD
      COMMON /THERMID/ ISOTH
c     ................................................................
      DOUBLE PRECISION DXP, DXU,SUW, SUE, HP,HU, HV,
     +                 DYP, DYV, SVN, SVS,U, V, P, T, RHP, EMP,
     +                 PREC, QREC,TREC, QIN, QOUT, QFACTOR,
     +                 REC,TORR,CK, BETAK, QSIDER,
     +                 ALFT, UC, TC, EC,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 QLEAD,QTRAIL,TQTRAIL,TLEAD,TRAIL
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL,BEARING ,KPAD,INERL,INERP,ITURB,INTER,ICAV1,MODEL,
     +        ISOTH
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION SVSK, SVNK, RHO, SGN,FACTOR,Inflow,Qright,
     +                 QRL(MAXNPOCK),QRR(MAXNPOCK),QLT(MAXNPOCK),
     +                 QRT(MAXNPOCK),TROLD(MAXNPOCK),DTRC,CPR,DUMY,
     +                  QsideB,Qsidej,QleadK,QtrailK ,RHO1,RHO2
      INTEGER I, J, JM1, JP1, K, JMIN, JMAX,Kstart,Kend,
     +        Icase,Kendp1,Isgn, KNY,KNYM1 ,IDCAV,JDCAV,ICAV(MAXNXT),JC
C ----------------------------------------------------------------------------
C -- FLOWT code                                                             --
C ----------------------------------------------------------------------------
  100 QIN=0.0D0                          ! zero variables
      QOUT=0.0D0			 !.............................!
      QsideB=0.0D0                       ! axial side flow
      Qlead(KPAD)=0.0D0                  ! Kpad leading edge circ. flow
      Qtrail(KPAD)=0.0D0                 ! Kpad trailing edge circ. flow
      TQtrail(KPAD)=0.0D0                ! Sum(rho h U T dy ) at trailing edge

      Kstart=1                           !.............................!
      Kend=NPA				 ! Indices for flow calculation
      Sgn=1.0D0				 ! on recess boundary
      Isgn=1				 !
      Icase=1                            ! TOP Boundary of recess
      KNY=NYI			         !.............................!

      IF (NPOCKET.GT.0) THEN             ! Calculation of recess temperatures
      DO I=1, NPOCKET                    !
          TROLD(I)=TREC(I)               ! For under-relaxation
          QSIDER(I)=0.0D0                !
          QRL(I)=0.0D0                   !
          QRR(I)=0.0D0                   !
          QLT(I)=0.0D0                   ! Sum(rho h u T dy ) at left edge
          QRT(I)=0.0D0                   ! Sum(rho h u T dy ) at right edge
      END DO
      END IF

 555  CONTINUE
      KNYM1=KNY-Isgn
c    !.................................................!
      IF ((BEARING.GE.2).OR.(NPOCKET.EQ.0)) GOTO 577   !=> SEAL & BEARING
c    !.................................................!

      Kendp1=Kend+Isgn
      SVSK=SVS(Kendp1)
      SVNK=SVN(Kendp1)
C    !.................!                 !
      DO I=1, NPOCKET                    !.................
C    !.................!                 !

          JMIN=(I-1)*NXI+NLC             ! for i-th pocket
          JMAX=I*NXI+1                   !
          K=Kendp1                       !

          DO J=JMIN, JMAX                ! flow at axial edge of recess
C        !.................!             !......................

          RHO=Rhp(J,Kend)*SVNK+Rhp(J,K)*SVSK

          QSIDER(I)=QSIDER(I)+HV(J,K)*
     +             DXP(J)*V(J,K)*RHO*Sgn ! = SUM(hxVxRho)

          END DO                         !
C        !.................!             !

          JMIN=JMIN-1                    !
C        !.................!             !........................

          DO K=Kstart,Kend, Isgn         ! (right-left) circ. flow
C        !.................!             !........................
          J=JMIN
          JP1=J+1
          RHO=Rhp(J,K)*SUE(J)+Rhp(JP1,K)*SUW(J)
          DUMY=DYP(K)*RHO*HU(J,K)*U(J,K) !
          QRL(I)=QRL(I)+DUMY             ! East(left) U flow into recess
          QLT(I)=QLT(I)+DUMY*T(JP1,K)    !
c                                        ! Recess left-edge flow * Temp
          J=JMAX			 !
          JP1=J+1                        !
          IF ((J.EQ.NXT).AND.(IFULL.EQ.1)) JP1=2
          RHO=Rhp(J,K)*SUE(J)+Rhp(JP1,K)*SUW(J)
          DUMY=DYP(K)*RHO*HU(J,K)*U(J,K) !
          QRR(I)=QRR(I)+DUMY             ! West(right)Uflow out of recess
          QRT(I)=QRT(I)+DUMY*T(J,K)      ! Recess right-edge flow * Temp
C        !..................!
          END DO
C        !..................!            ! K

C    !......................!            !
      END DO                             ! NEXT Recess, I=1,2,...NPOCKET
C    !......................!            !

      GOTO 588

C ---------------------------------------! SIDE FLOWS:

 577  IF (Icase.EQ.2) GOTO 588

C    !......................!            !..........................
      DO J=1, NXT-1                      ! INLET flow FOR SEAL at Y=0
C    !......................!            !..........................
        QIN=QIN+DXP(J)*HP(J,1)*V(J,1)*RHP(J,1)
      END DO                             !
C    !......................!            !

 588  CONTINUE

C    !......................!            !..........................
      DO J=1, NXT-1                      ! outlet flow at sides of bearing
C    !......................!            !..........................
C     here assume flow does not enter from exit ends (Y=L) if P<Pcav
        RHO=Rhp(J,KNYM1)*SVN(KNY)+Rhp(J,KNY)*Svs(KNY)
        IF (P(J,KNYM1).GE.PCAV) THEN     !
            QSIDEJ=DXP(J)*HV(J,KNY)*V(J,KNY)*RHO*Sgn
            QOUT=QOUT+QSIDEJ
            QSIDEb=QSIDEb+QSIDEJ
        END IF
      END DO                             !
C    !......................!            !..............

      IF (IFULL.EQ.1) GOTO 666           !=>  360-deg bearings or seals
C                                        !    no need for circ. pad flows

C ..............................................!.............................!
C ONLY FOR PAD BEARINGS:                        !
C ..............................................! FLOW AT LEFT&RIGHT
C                                               ! PAD BOUNDARIES
          J=NXT-1                               !
          JM1=J-1                               ! Note: Cavitation is not
          JP1=NXT                               ! considered for BEARING=1
          IDCAV=0                               ! =0: No cavitation
C    !......................!                   !
       IF (BEARING.EQ.3) THEN                   ! Grooved pad journal bearing
C    !......................!                   ! Determine cavitation boundaries
             ICAV(NXT)=0                        ! at bearing midplane Y=0
          DO JC=2,J                             !
             IF (P(JC,1).LT.PCAV) THEN          ! if P< Pcav => cavitated
                ICAV(JC)=JC                     !
                IDCAV=1                         !
             ELSE 			        ! J=NXT-1
                ICAV(JC)=0                      !
             END IF                             !
          END DO                                !
c        !.......................!              !.....................
          IF (IDCAV.EQ.1) THEN                  ! If cavitated, decide
c        !.......................!              ! on cavitation type
             DO JC=2,J                          !.....................
                IF (ICAV(JC).GT.0) THEN         ! (1) Diverge-converge type cavi.
                   IF (ICAV(JC+1).EQ.0) THEN    ! ICAV(JC)=0,2,3,...,JC,0,0,0,...
                      JDCAV=1                   !
                      GO TO 222                 !
                   END IF                       !
                ELSE IF (ICAV(JC+1).GT.0) THEN  ! (2) Converge-diverge type cavi.
                   JDCAV=JC+1                   ! ICAV(JC)=0,0,0,...,JC+1,JC+2,...
                   GO TO 222                    !
                END IF                          !
             END DO                             !
  222        CONTINUE                           !
c        !.......................!              ! cavitation type
          END IF                                ! end search
c        !.......................!              !
C    !......................!                   !
       END IF                                   ! End of grooved pad JBs
C    !......................!                   !
C        !.................!
          DO K=Kstart,KNY, Isgn                 !
C        !.................!                    !..................!
          RHO1=Rhp(1,K)*SUE(1)+Rhp(2,K)*SUW(1)  ! Density at leading edge

c        !......................!               !.......................!
          IF (IDCAV.EQ.0) THEN                  ! No cavitation
c        !......................!               ! on PAD  ..............!
             QleadK= DYP(K)*RHO1*HU(1,K)*U(1,K) !
             RHO2=Rhp(J,K)*SUE(J)+Rhp(JP1,K)*SUW(J)
             QtrailK=DYP(K)*RHO2*HU(J,K)*U(J,K) ! J=NXT-1
c        !......................!               !.......................!
          ELSE IF (JDCAV.EQ.1) THEN             ! First type cavitation
c        !......................!               !.......................!
             QleadK= DYP(K)*RHO1*HU(1,K)*MSPEED ! =shear flow approximately
             QtrailK=QleadK                     ! =leading edge flow `'
c        !......................!               !.......................!
          ELSE IF (JDCAV.GT.1) THEN             !
c        !......................!               !.......................!
             QleadK= DYP(K)*RHO1*HU(1,K)*U(1,K) ! Second type cavitation
             RHO2=Rhp(JDCAV,K)*SUE(JDCAV)       !
     +      +Rhp(JDCAV+1,K)*SUW(JDCAV)          !
             QtrailK=DYP(K)*RHO2*HU(JDCAV,K)    ! The flow at trailing edge
     +      *U(JDCAV,K)                         ! =flow at inception of cav.
c        !......................!               !.......................!
          END IF                                !
c        !......................!               !.......................!
          Qlead(KPAD)=Qlead(KPAD)+Qleadk        !
          Qtrail(KPAD)=Qtrail(KPAD)+QtrailK     !
          QOUT=QOUT-Qleadk+QtrailK              !
c 					        ! PAD
          TRAIL(KPAD,K)=T(J,K)                  ! Trailing temperature
          TQtrail(KPAD)=TQtrail(KPAD)           ! Sum(Q*T) at trailing edge
     +                 +QtrailK*TRAIL(KPAD,K)   !

          IF (QLEAD(KPAD).LE.0.0D0) THEN        ! Flow out of the pad
             TLEAD(KPAD,K)=T(2,K)               ! Then Leading edge temperature
          END IF                                ! = land temperature
C        !..................!
          END DO            ! K = Kstart,  KNY
C        !..................!

C The Product T * Q is the energy being convected at pad edges
c here we have considered the specific heat as constant across the edge
c

C...............................................!...........................!
 666  IF (ISYM.EQ.0) THEN
C    !........................!          ! ASYMMETRIC BEARING

         IF (BEARING.EQ.2) THEN          !=> SEAL
            Qright=QIN
            GOTO 778
         END IF

         IF (Icase.EQ.2) GOTO 777

         Kstart=-1
         Kend=-Npa
         Isgn=-1
         Sgn=-1.0D0
         KNY=-Nyi
         Icase=2
         Qright=QOUT
         goto 555

C    !........................! ISYM=0, ASYMMETRIC BEARING
      ELSE
C    !........................!
        IF (NPOCKET.GE.1) THEN
         DO I=1, NPOCKET
           QSIDER(I)=2.0D0*QSIDER(I)! For recess temperature: Trec
           QRL(I)=2.0D0*QRL(I)      ! Full recess energy balance
           QRR(I)=2.0D0*QRR(I)      ! for both symmetric (ISYM=1)
           QLT(I)=2.0D0*QLT(I)      ! and asymmetric (ISYM=0) cases
           QRT(I)=2.0D0*QRT(I)      !
         END DO
        END IF
         Qright=QOUT
         QIN=2.0D0*QIN
         QOUT=2.0D0*QOUT
         QsideB=2.0D0*QsideB
         Qlead(KPAD)=2.0D0*Qlead(KPAD)
         Qtrail(KPAD)=2.0D0*Qtrail(KPAD)
         TQtrail(KPAD)=2.0D0*TQtrail(KPAD)
C    !........................! ISYM=1, SYMMETRIC BEARING
      END IF
C    !........................!

 777  continue

C    !...................................!
      IF (NPOCKET.GT.0) THEN
C      !......................!          !
        DO I=1, NPOCKET                  !...................
          QREC(I)=QSIDER(I)-QRL(I)+QRR(I)! Flow at i-th recess
          QIN=QIN+QREC(I)                ! Inlet flow
        END DO

C       !......................!         !
        IF (ISOTH.EQ.1) GOTO 778         ! Uniform Recess temperature
C       !......................!         !


C Calculation of Recess temperature
C =================================      ! REC=REY/Eckert #
c      !----------------!                !...........................!
        DO I=1, NPOCKET                  ! DOES not include radial heat flow
c      !----------------!                !...........................!
         JP1=(I-1)*NXI+NLC+1             ! SIMPLE Recess energy balance eqn

         CPR=CK(JP1,1)                   ! fluid specific heat

         DTRC=SPEED/REC/CPR              ! rotor-speed !####

         TREC(I)=1.0D0                   != Supply temperature

         IF((QRL(I).GT.0.D0).and.(QRR(I).GT.0.D0)) THEN
            TREC(I)=(QREC(I)+QLT(I)+TORR(I)*DTRC)/(QSIDER(I)+QRR(I))
            TREC(I)=TROLD(I)+ALFT*(TREC(I)-TROLD(I))

         ELSE IF((QRL(I).LT.0.D0).and.(QRR(I).LT.0.D0)) THEN
            TREC(I)=(QREC(I)-QRT(I)+TORR(I)*DTRC)/(QSIDER(I)-QRL(I))
            TREC(I)=TROLD(I)+ALFT*(TREC(I)-TROLD(I))
         END IF
c      !----------------!                !...........................!
        END DO                           ! NEXT Recess, I=1,2,...NPOCKET
c      !----------------!                !...........................!

      CALL SETTREC                       ! Set recess temp & props

C    !...................................!
      END IF                             ! End of Recess temperature
C    !...................................! if NPOCKET > 0
 778  Inflow=QIN                         !
      IF (MODEL.EQ.1) THEN               ! ==> 2 row HJB
        Inflow=QIN                       !
        QIN=Qright                       !
      ELSE IF (ISYM.eq.0) THEN           !=> asymmetric bearing
        QIN=Qright                       !
        Inflow=QIN                       !
      END IF                             !
C    !...................................!

      IF (BEARING.GE.2) THEN
      IF (DABS(Inflow).LE.1.0D0) Inflow=1.0D0
        IF (BEARING.EQ.3) THEN   ! for journal pad bearing
           QIN=QsideB            ! set side flow to QIN
           QOUT=Qlead(KPAD)      ! & leading flow to QOUT ???
        END IF
      END IF

      END
c

C *****************************************************************************
C **                                                                         **
C **  Subroutine Settrec                                                     **
C **                                                                         **
C **  SETTREC:  Resets recess temperatures and properties on pockets.        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE SETTREC

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /PARRAY/ P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RHOEMU/ RHOP(MAXNXT, -MAXNYI:MAXNYI),
     +                EMUP(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +                BT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /LIQUID/ TEMPK, VSOUND,IF,IL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /SWITCH/ IPROP

      DOUBLE PRECISION RHOP,EMUP,P,PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 TEMPK, VSOUND ,CK,BT,THC ,T
      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT
      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL ,IF,IL,
     +        IFULL,ISYM, ICSTEP ,IPROP
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER IP, JSTA, JSTO, J, K, KK ,Kstart
      DOUBLE PRECISION Rhor, EMur ,CKR,BTR,THR

C ----------------------------------------------------------------------------
C --  SETTREC code                                                          --
C ----------------------------------------------------------------------------
      KK = NPA -1
      Kstart=ISYM + (ISYM-1)*KK

C    !.............................!....................!
      IF (IPROP.EQ.1) THEN         ! Update properties
C    !.............................!....................! IP loop
       DO IP=1, NPOCKET            ! Update
          JSTA=(IP-1)*NXI+NLC+1    ! temperature=TREC at interior
          JSTO=IP*NXI              ! of recesses

          CALL LOCPROPS(RHOR,EMUR,CKR,
     +    BTR,THR,PREC(IP),TREC(IP))

C        !.........................!
          DO J=JSTA, JSTO          ! Only interior points (w/o edges)
              DO K=Kstart, KK ,1   !
                  T(J, K)=TREC(IP) ! are at uniform temperature,
                  Rhop(j,k)=Rhor   ! pressures and properties.
                  EMup(j,k)=EMur   !
                  CK(j,k)  =CKR    ! specific heat
                  BT(j,k)  =BTR    ! expansion coefficient
                  THC(j,k) =THR    ! thermal conductivity
              END DO               !
              K=NPA                ! At top edge, pressure is not uniform
              T(J,K)=TREC(IP)      ! but temperature is.
              CALL LOCPROPS(RHOP(J,K),EMUP(J,K),CK(J,K),
     +        BT(J,K),THC(J,K),P(J,K),TREC(IP))
          END DO                   !
C        !.........................!
          IF (ISYM.EQ.0) THEN      !
          DO J=JSTA, JSTO          ! At bottom edge,
              K=-NPA               ! pressure is not uniform
              T(J,K)=TREC(IP)      ! but temperature is.
              CALL LOCPROPS(RHOP(J,K),EMUP(J,K),CK(J,K),
     +        BT(J,K),THC(J,K),P(J,K),TREC(IP))
          END DO                   !
          ENDIF                    !
C        !.........................!
                                   !
       END DO                      !
C                                  !
c                                  !
C    !.............................!....................!
      ELSE                         ! NO update of properties
C    !.............................!....................!
       DO IP=1, NPOCKET            ! Update
          JSTA=(IP-1)*NXI+NLC+1    ! temperature=TREC at interior
          JSTO=IP*NXI              ! recess
          DO J=JSTA, JSTO          ! Only interior points (w/o edges)
              DO K=Kstart, KK ,1   !
                  T(J, K)=TREC(IP) ! are at uniform temperature,
              END DO               !
              K=NPA                ! At top edge, pressure is not uniform
              T(J,K)=TREC(IP)      ! but temperature is.
          END DO                   !
C        !.........................!
          IF (ISYM.EQ.0) THEN      !
          DO J=JSTA, JSTO          ! At bottom edge,
              K=-NPA               ! pressure is not uniform
              T(J,K)=TREC(IP)      ! but temperature is.
          END DO                   !
          ENDIF                    !
C        !.........................!
                                   !
       END DO                      !
C    !.............................!....................! END OF IP
      END IF
C    !.............................!....................!

      END
c

C *****************************************************************************
C **                                                                         **
C **  Subroutine Torquer                                                     **
C **                                                                         **
C **  TORQUE:  Calculates drag friction torque over recess areas             **
C **           Used for calculation of recess temperatures on FLOWT          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TORQUER

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC, HRECD
      COMMON /TORREC/ TORR(MAXNPOCK) ! Torques over recess areas
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /RECJET/ PRECdo(MAXNPOCK),PRECup(MAXNPOCK),
     +                PRjet(MAXNPOCK,MAXNPOCK+2)
      COMMON /URECJET/ UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL

      DOUBLE PRECISION P, DXP, DXU, SUW, SUE, HP, HU, HV, TORR,
     +                 DYP, DYV, SVN, SVS,
     +                 REP, REY, MSPEED, SPEED,
     +                 SWIRL, PCAV, ALFU, BETU, ALFP, RHOP,EMUP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 PREC,TREC,QREC,QIN,QOUT,QFACTOR,
     +                 PRECdo, PRECup,PRJET,UREC

      INTEGER NPOCKET, NLC, NPC, NLA, NPA, NPAP1, NXI, NYI, NXT,
     +                 ISYM, ICSTEP ,IFULL
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER I,J,K,JSTA,JSTO,JP1,Kstart, Kstep,Kend,Kendm1,ICASE,
     +        IJ
      DOUBLE PRECISION H,DX,GAMA,FNGAMA,RHOR,EMUR,Uedge,
     +                 CKR,BTR,THR
C ----------------------------------------------------------------------------
C --  TORQUER code                                                          --
C ----------------------------------------------------------------------------
      DO I=1, NPOCKET                   ! TORQUER: calculates DRAG
         TORR(I)=0.0D0                  ! torque x fluid speed at rotor surface
      END DO                            ! within recess area
C    !..................................! i.e. Shear power on recess.

      Icase=1
      Kstep=1
c   !....................!

 777  Kstart=Kstep
      Kend=Kstep*NPA
      Kendm1=Kend-kstep
c   !....................!              !sweep on pockets
      DO I=1, NPOCKET                   !.................
C   !....................!              !

         JSTA=(I-1)*NXI+NLC             ! Left edge of i recess
         JSTO=I*NXI                     ! +1 =right `' `' `'
         IJ=0                           !
         CALL LOCPROPS(RHOR,EMUR,CKR,   ! fluid props. at recesses
     +      BTR,THR,PREC(I),TREC(I))    ! pressure & temperature

C     !..................!              ! K=Kend=+/- NPA
          DO J=JSTA, JSTO               !
c     !..................!              !
c           within recess assume no axial velocity / smooth surface
c           and YES circumferential pressure gradient
            IJ=IJ+1                     !
            JP1=J+1                     !
            DX=DXU(J)                   !
            H=HU(J,Kend)+HRECD          !On circumferential edge of
            Uedge=UREC(I,IJ)            !recess, 1/2 U-CV
c                                       !
            GAMA=FNGAMA(Uedge,0.0D0,H,SPEED,REP,RHOR,EMUR)
            TORR(I)=TORR(I)+(DX*DYP(Kend)/2.0D0)*(GAMA+
     +                 H*(P(JP1,Kend)-P(J,Kend))/2.0D0/DX)

C          !............................! within the recess areas
            DO K=Kstart,Kendm1,Kstep    ! .....................
               H=HU(J,K)+HRECD          ! Recess film thickness
               GAMA=FNGAMA(Uedge,0.0D0,H,SPEED,REP,RHOR,EMUR)
               TORR(I)=TORR(I)+DX*DYP(Kend)*(GAMA+
     +                 H*(PRjet(I,IJ+1)-PRjet(I,IJ))/2.0D0/DX  )
            END DO                      ! Torque over recess
C          !............................!

C    !....................!             !
         END DO                         !sweep along X on recess
C    !....................!             !

      END DO                            ! Next Pocket
C    !--------------------!             !----------------------!

      IF (ISYM.EQ.0) THEN
c    !......................! ISYM=0, asymmetric bearing
        IF (Icase.EQ.2) RETURN
        Kstep=-1
        Icase=2
        GOTO 777

      ELSE
c    !......................! ISYM=1, symmetric bearing

        DO I=1, NPOCKET
           TORR(I)=TORR(I)*2.0D0
        END DO

      END IF
c    !......................!

      END




C *****************************************************************************
C **                                                                         **
C **  Subroutine Groovetemp                                                  **
C **                                                                         **
C **  GROOVETEMP:Energy Balance Equation at Groove to Find Groove-Edge Or    **
C **             Pad Leading and Trailing Temperatures--                     **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE GROOVETEMP(JPAD,KPAD,Kstart,TDIFF)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /QIOPAD/ QLEAD(MAXNPAD),QTRAIL(MAXNPAD),TQTRAIL(MAXNPAD)
      COMMON /TIOPAD/ TLEAD(MAXNPAD,-MAXNYI:MAXNYI),
     +                TRAIL(MAXNPAD,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /SUMPCOND/ TSUMP,TSUMPd,RHOSU,RHOSd
      COMMON /TGROOVE/ LEMDA,DELTA

      DOUBLE PRECISION Qlead,Qtrail,TQtrail,Tlead,Trail, T,
     +                 TSUMP,TSUMPd,RHOSU,RHOSd ,LEMDA,DELTA
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        ISYM,ICSTEP

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION TDIFF,TLEADnew, TLEADold
      INTEGER Kstart,K, JPAD,KPAD
C ----------------------------------------------------------------------------
C --  GROOVETEMP code                                                       --
C ----------------------------------------------------------------------------
C  This routine is only called from CALCPADT
c  when Qlead(KPAD) > 0, i.e. if flow enters the current pad
c  ...........................................................................
c  The notation is the following at a groove between pads
C
C
C         Leading edge =========== JPAD=========== Trailing edge  --> X
C
C
C         Leading edge =========== KPAD=========== Trailing edge  --> X
C
C
C                                         | Tsump
C                                         v  Qgroove
C   JPAD trailing                                           KPAD leading
C        edge ======  ---> Qtrail (Jpad)    ---> Qlead(KPAD)  edge =======
C                          Trail                 Tlead
C
C    Qtrail(Jpad) >0 : flow leaves JPAD,    <0: flow enters JPAD
C    Qlead (Kpad) >0 : flow enters KPAD,    <0: flow leaves KPAD
C
C  IF QLEAD >0 then some flow mixing and thermal mixing with groove feed flow
c              and upstream flow and temperature must occur at the
c              groove. The KPAD leading edge temperature (Tlead) is a function of
c              the flow rates, the sump temperature and the trailing edge
c              temperature (Trail)
C
C
C  IF QLEAD < 0, then KPAD leading edge temperature is obtained from solution
C                of thermal field within pad (SUB FLOWT)
C
C
C Case 1:
C  IF Qtrail<0 and Qlead>0 then fresh fluid enters both pad edges with
C                          a Temperature equal to Tsump
C   i.e.  Trail=Tlead=Tsump  (1)
C
C Case 2:
C  IF Qtrail>0 and Qlead>0, gives by conservation of flow
C
C     (2)  Qlead = Lembda*Qtrail + Qgroove, where Qgroove is feeding flow at Tsump
C
C  and     Lembda coefficient denoting % of mixing occuring
C
C  then a simple balance of thermal energies leaving JPAD, [Lembda x Qtrail x Trail]
C  plus energy IN from feeding groove Qgroove*Tsump = (Qlead-Lembda*Qtrail)*Tsump
c  should be equal to the energy INTO KPAD: Qlead*Tlead.
c
C  For equation (2) to make sense, Qlead > Lembda*Qtrail, i.e. Qgroove >0
C
C  Otherwise if Qlead< Lembda*Qtrail, then
c
C  Case 3:
c  i.e. Qgroove <=0 which means that there is no mixing and most of hot fluid
c                   from the upstream pad just enters the leading edge of KPAD
c                   at a temperature equal to the exit temperature of the
c                   trailing edge JPAD.
c
c=============================================================================
c     LEMDA= 0 to 1.0                          ! -> set on calcmesht.f by USER
c============================================================================

C
      TLEADold=TLEAD(KPAD,Kstart)              ! OLD loading edge temperature
c                                              !
      TDIFF=0.0D0                              ! (Tnew-Told) at leading edge


c  CASE  #1; Qtrail<0 and Qlead>0
c     ............................................................
      IF (QTRAIL(JPAD).LE.0.0D0) THEN          !
c     ............................................................
        DO K=Kstart,NYI                        ! Flow into both
          TLEAD(KPAD,K)=Tsumpd                 ! pads from feeding source
          TRAIL(JPAD,K)=Tsumpd                 ! (groove) at Tsump
        END DO                                 !
        TLEADnew=Tsumpd                        !
c     ............................................................
      ELSE                                     ! Flow leaving up-pad (JPAD)
c     ............................................................


c   CASE #2: Qlead > Lembda x Qtrail, i.e. Qgroove=QIN >0
c                                          some fresh fluid at Tsump enters KPAD
C       !......................!               !
        IF (QLEAD(KPAD).GE.LEMDA*QTRAIL(JPAD)) THEN
C       !......................!               ! Mixing occurs
c                                              !
         TLEADnew=(LEMDA*TQTRAIL(JPAD)+        ! Formula given by J.Mitsui
     +   (QLEAD(KPAD)-LEMDA*QTRAIL(JPAD))      ! 1983,JLT; Lemda=0.4 to 1.0,
     +   *Tsumpd)/QLEAD(KPAD)                  ! typical value=0.6 to 0.8
c                                              !
C       !......................!               !
        ELSE                                   !
C       !......................!               !........................!

c  CASE #3: Qlead < Lembda x Qtrail, i.e. Qgroove=QIN <0
c                                         NO thermal mixing occurs
c                                              !
          TLEADnew=0.0D0                       ! find average trailing edge
          DO K=Kstart,NYI                      ! temperature and set it
            TLEADnew=TLEADnew+TRAIL(JPAD,K)    ! equal to new
          END DO                               ! leading edge temperature
            TLEADnew=TLEADNEW/(NYI-1+Kstart)   !
C       !......................!               !
        END IF
C       !......................!               !

         DO K=Kstart,NYI                       ! SET Average uniform temperature
           TLEAD(KPAD,K)=TLEADnew              ! at the leading pad edge
         END DO                                ! equal to a constant

c     ............................................................
      END IF
c     ............................................................
c                                              !
         TDIFF=DABS(TLEADnew-TLEADold)         ! Temperature DIFFerence
c                                              !  (New-Old) at KPAD leading edge
c                                              !.............................!
      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Tcoefs--for solution of the first-order energy equation     **
C **                                                                         **
C **  TCOEFS:  Calculates At coefficients from Uo, Vo, Po and To solutions   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TCOEFS

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /BTYPE/ BEARING
      COMMON /HJBSYM/ISYM, ICSTEP
      COMMON /LRBOUND/ LEFTBC,RIGHTBC

      COMMON /TCOEF/
     + GTU(MAXNXT, -MAXNYI:MAXNYI), GTV(MAXNXT, -MAXNYI:MAXNYI),
     + GTP(MAXNXT, -MAXNYI:MAXNYI), GTH(MAXNXT, -MAXNYI:MAXNYI),
     + GTPI(MAXNXT,-MAXNYI:MAXNYI), APTI(MAXNXT,-MAXNYI:MAXNYI),
     + APT(MAXNXT, -MAXNYI:MAXNYI), AWT(MAXNXT, -MAXNYI:MAXNYI),
     + AET(MAXNXT, -MAXNYI:MAXNYI), AST(MAXNXT, -MAXNYI:MAXNYI),
     + ANT(MAXNXT, -MAXNYI:MAXNYI), BT1(MAXNXT, -MAXNYI:MAXNYI),
     + BT2(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)

      DOUBLE PRECISION APT,AWT,AET,AST,ANT, APTI,BT1,BT2,
     +                 GTU,GTV,GTP,GTH,GTPI,
     +                 DYP,DYV,SVN,SVS , U, V, P, T
      INTEGER INERL, INERP, ITURB, INTER, ICAV, MODEL,BEARING,
     +        JPMIN, JPMAX, JPSTART, JPSTOP,ISYM, ICSTEP,IFULL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI ,NXT, NYI,
     +        LEFTBC, RIGHTBC

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DYPK,Dir,zero, xmin,xmax,PS(MAXNXT),TS(MAXNXT)
      INTEGER I,J,K,KM1,KP1,IP,Icsym, Kstart, Kend, Kstep, KV
C ----------------------------------------------------------------------------
C --  TCOEFS code                                                           --
C ----------------------------------------------------------------------------
      zero=0.0D0
      Kstep= ISYM+(ISYM-1)*NYI
      Kstart=ISYM+(ISYM-1)
      IF (BEARING.EQ.2) Kstep=1                ! => SEAL

C   !..........................................! ISYM=0/1
      DO J=Kstep, NYI                          !
          DO I=1, NXT                          !
              APT(I, J)=zero                   !
              APTI(I,J)=zero                   !
              AWT(I, J)=zero                   !
              AET(I, J)=zero                   !
              ANT(I, J)=zero                   !
              AST(I, J)=zero                   !
              BT1(I, J)=zero                   !
              BT2(I, J)=zero                   !
          END DO                               !
      END DO                                   !

c   !...............................!.start with right side  of bearing
      Kstep=1
      Icsym=1
      Dir=1.0D0

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
      LEFTBC=1+IFULL                             ! =1: GROOVE, =2: RECESS
      RIGHTBC=2                                  !..........................

C    !...................!                       !
      DO IP=1, NPOCKET+1-IFULL                   ! Sweep on lands between
C    !...................!                       ! pockets: Inter=1
          CALL SETLIM(IP, 1)                     !

          DO J=JPMIN,JPMAX
            PS(J)=P(J,Kstart)                    ! Psouth at Kstart
            TS(J)=T(J,Kstart)                    ! Tsouth at Kstart
          END DO

          IF (IP.GT.NPOCKET) RIGHTBC=1
          K=0                                    !
          KM1=Kstart                             !

  100     K=K+Kstep                              ! Dypk= Y size Ucv
          KP1=K+Kstep                            !
          DYPK=DYP(K)                            !

          CALL TCALC(DYPK,K,KP1,KM1,K,PS,TS,Dir) ! Calculate At coefficients

          IF (K.EQ.Kend) goto 200                !

          KM1=K                                  !

          GOTO 100                               !
C    !...................!
  200 LEFTBC=2
      END DO                                     ! ...........................
C    !...................!                       !
C                                                !
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

      DO J=JPMIN, JPMAX                          ! for TCALC on extended land
        PS(J)=P(J,K)*SVS(K)+P(J,KM1)*SVN(K)      ! Ps=Psouth
        TS(J)=T(J,K)*SVS(K)+T(J,KM1)*SVN(K)      ! Ts=Tsouth
      ENDDO

      Kend=Kstep*NYI                             !
      LEFTBC=1-IFULL                             ! BC's at sides
      RIGHTBC=LEFTBC


  700 CONTINUE

c     .....................!
      IF (K.eq.Kend) THEN
c    !.....................!

c    !................!
       IF ((Icsym.EQ.1).AND.(MODEL.EQ.1)) THEN  ! 2row HJB
c    !................!
           KP1=NYI+1
           DO J=1, NXT
             V(J,KP1)=-V(J,K)                   ! symmetry line
           END DO
           CALL TCALC(DYPK,K,KP1,KM1,K,PS,TS,Dir) ! Calcule At coefficients
c    !................!
       END IF
c    !................! 1 row HJB
           GOTO 800

c    !.....................!
       ELSE
c    !.....................! K < KEND  (within film lands)

         KP1=K+Kstep

         CALL TCALC(DYPK,K,KP1,KM1,KV,PS,TS,Dir) ! Calcule At coefficients

         KM1 = K
         K = KM1 + Kstep
         KV = K
         DYPK = DYP(K)

         GOTO 700
c    !.....................!
       END IF
c    !.....................! K Counter on film lands


c !..............................................! Direct for HJB asymmetry

c    !.....................! ISYM=0
  800 IF (ISYM.eq.0) THEN
c    !.....................!
        IF ((BEARING.EQ.2).OR.(Icsym.eq.2)) RETURN

        IF (MODEL.EQ.1) THEN     ! TWO recess ROW HJB
          DO J=1, NXT
             ANT(J,Kend)=zero
             GTV(J,Kend)=zero
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


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Tcalc--Called by TCOEFS                                     **
C **                                                                         **
C **  TCALC:  Calculates At coefficients from Uo, Vo, Po and To solutions    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TCALC(DYPK,K,KP1,KM1,KV,PS,TS,Dir)

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
     +               HV(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +             BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TCOEF/
     + GTU(MAXNXT, -MAXNYI:MAXNYI), GTV(MAXNXT, -MAXNYI:MAXNYI),
     + GTP(MAXNXT, -MAXNYI:MAXNYI), GTH(MAXNXT, -MAXNYI:MAXNYI),
     + GTPI(MAXNXT,-MAXNYI:MAXNYI), APTI(MAXNXT,-MAXNYI:MAXNYI),
     + APT(MAXNXT, -MAXNYI:MAXNYI), AWT(MAXNXT, -MAXNYI:MAXNYI),
     + AET(MAXNXT, -MAXNYI:MAXNYI), AST(MAXNXT, -MAXNYI:MAXNYI),
     + ANT(MAXNXT, -MAXNYI:MAXNYI), BT1(MAXNXT, -MAXNYI:MAXNYI),
     + BT2(MAXNXT, -MAXNYI:MAXNYI)

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

      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /TISOBJD/ TBJ,TBS
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB

      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /THERMID/ ISOTH

      DOUBLE PRECISION APT,AWT,AET,AST,ANT, APTI,BT1,BT2,
     +                 GTU,GTV,GTP,GTH,GTPI,
     +                 CK, BETAK, HCB, HCJ ,THC ,DKP,DKT,
     +                 DRHOP,DRHOT,DEMUP,DEMUT,DBKP,DBKT,DCPP,DCPT,
     +                 DXP, DXU, HP, HU, HV, SUW, SUE,
     +                 DYP,DYV,SVN,SVS, TBS, TBJ,
     +                 AMOD, BMOD, RUGR, RUGS, EXPO,
     +                 U, V, P, T, RHP, EMP,
     +                 ALFT, UC, TC, Ec ,
     +                 REP, REY, MSPEED, SPEED,
     +                 SWIRL, PCAV, ALFU, BETU, ALFP ,Dir,
     +                 TBOUT, THERMALK, ROUTER, HKB

      INTEGER JPMIN, JPMAX, JPSTART, JPSTOP,ISOTH,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL, KV

C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION DX,DYPK,DAREA, U0,V0, BT,UC0,UC1,REC,REYCK,
     +                 FWP, FEP, FSP, FNP,
     +                 RHOW, RHOE, RHOS, RHON, ZERO,
     +                 PW,PE,TW,TE, SANB,
     +                 PN(MAXNXT),PS(MAXNXT),TN(MAXNXT),TS(MAXNXT),
     +                 GTUU,GTVV,GTM,GTR,GTPP,GTT,GTHH,GCP,GTB,GTK,
     +                 H, REPT, DPX,DPY,DTX,DTY, KXO, KXR, KXS,
     +                 UR, UR2, RS, RR, CR, CS, BR, BS, FRO, FSO,
     +                 DUMY, GR, GS, FR1, FS1, CCR, CCS,
     +                 rho, mu0, CP, expm1 ,UW ,QB0,QJ0,QS0, Den
      INTEGER J, JJ, JJJ, JJMAX, JP1, K, KM1, KP1

C ----------------------------------------------------------------------------
C --  TCALC code                                                            --
C ----------------------------------------------------------------------------
C No cavitation model is accounted on the AT coefficients
C turbulence coefficients based on Moody's equation

      Den = 1.0D0
      ZERO=0.0D0
      REC=REY/EC                               ! Reynolds No./Eckert No.

      JJJ=2-JPMIN                              !
      JJMAX=JPMAX+JJJ                          !
      J=JPSTART                                !

 100  RHOW=Rhp(j,k)*Sue(j)+Rhp(jpmin,k)*Suw(j) !
      FWP=HU(J,K)*DYPK*RHOW*U(J,K)             ! Flow on west side of Pcv
      Pw=P(J,K)*SUE(J)+P(JPMIN,K)*SUW(J)       ! Pw at J=JPMIN
      Tw=T(J,K)*SUE(J)+T(JPMIN,K)*SUW(J)       ! Tw at J=JPMIN
      Uw=U(J,K)
C ......................                       !
      DO J=JPMIN, JPMAX                        !----------------------------
C ......................                       !

          JJ=J+JJJ                             ! Sweep from left-right on
          JP1=J+1                              ! P cvs
          DX=DXP(J)
          DAREA=DX*DYPK                        !
          CP=CK(J,K)                           ! Dimensionless specific heat
          REYCK=REC*CP                         ! Rey*Cp/Ec

          AWT(J,K)=DMAX1(+FWP,ZERO)*REYCK      ! Awt coeff.

          Rhoe=Rhp(j,k)*Sue(j)+Rhp(jp1,k)*Suw(j)
          FEP=HU(J,K)*DYPK*RHOE*U(J,K)         !
          AET(J,K)=DMAX1(-FEP,ZERO)*REYCK      ! Aet coeff.

          Rhon=Rhp(j,kp1)*Svs(kp1)+Rhp(j,k)*Svn( kp1)
          FNP=HV(J,KP1)*DX*RHON*V(J,KP1)*Dir   ! North V flow on Pcv
          ANT(J,K)=DMAX1(-FNP,ZERO)*REYCK      ! Ant coeff.

          Rhos=Rhp(j,k)*Svs(k)+Rhp(j,km1)*Svn(k)
          FSP=HV(J,KV)*DX*RHOS*V(J,KV)*Dir     ! South V flow on Pcv
          AST(J,K)=DMAX1(+FSP,ZERO)*REYCK      ! Ast coeff.

      SANB=AET(J,K)+AWT(J,K)+AST(J,K)+ANT(J,K) !
      APT(J,K)=SANB                            !

C ----------------------------------------------------------------
C     Find the gama coefficients for the 1st order energy equation
C ----------------------------------------------------------------

      PE=P(J,K)*SUE(J)+P(JP1,K)*SUW(J)         ! FIND Pe
      TE=T(J,K)*SUE(J)+T(JP1,K)*SUW(J)         ! FIND Te
      PN(J)=P(J,KP1)*SVS(KP1)+P(J,K)*SVN(KP1)  ! Pn on P-C.V.
      TN(J)=T(J,KP1)*SVS(KP1)+T(J,K)*SVN(KP1)  ! Tn on P-C.V.
      DPX=(PE-PW)/DX                           ! Pressure gradient in x
      DTX=(TE-TW)/DX                           ! Temperature gradient in x
      DPY=(PN(J)-PS(J))*Dir/DYPK               ! Pressure gradient in y
      DTY=(TN(J)-TS(J))*Dir/DYPK               ! Temperature gradient in y
      U0=(Uw+U(J,K))/2.0D0                     ! U-velocity at P-node
      V0=(V(J,K)+V(J,KP1))/2.0D0               ! V-velocity at P-node
      MU0=EMP(J,K)
      RHO=RHP(J,K)
      H=HP(J,K)
      BT=BETAK(J,K)*T(J,K)	               ! betat*tk, not one variable!
C .............................................!
C         Find KXO,KXR
C .............................................!
      UC0=U0*U0+V0*V0+U0*MSPEED                ! For bulk flow model find
      UC1=MSPEED*MSPEED-U0*SPEED               ! shear coefficients
      UR=U0-SPEED                              ! based on Moody's friction factor
      UR2=UR*UR                                !
      REPT=REP*RHO/MU0                         !
      RS=H*DSQRT(U0*U0+V0*V0)*REPT             ! Stator Reynolds number/rep
      RR=H*DSQRT(UR2+V0*V0)*REPT               ! Rotor Reynolds number/rep
      CR=RUGR*10000.0D0/H                      ! Rotor roughness coeff.
      CS=RUGS*10000.0D0/H                      ! Stator roughness coeff.
      BR=BMOD/RR                               !
      BS=BMOD/RS                               !
      FRO=AMOD*(1.0D0+(CR+BR)**EXPO)           ! Zeroth rotor friction factor
      FSO=AMOD*(1.0D0+(CS+BS)**EXPO)           !        stator friction factor
      KXR=RR*FRO                               !
      KXS=RS*FSO                               !
      KXO=0.5D0*(KXS+KXR)                      ! Zeroth turbulent shear coeffs
      DUMY=12.0D0                              !
      KXO=DMAX1(DUMY, KXO)                     ! Select laminar or turbulent
      KXR=DMAX1(DUMY, KXR)                     !
      KXS=DMAX1(DUMY, KXS)                     !

c    !=========================================!..............................!
      IF (KXO.GT.DUMY) THEN                    ! Turbulent flow condition
c    !=========================================!..............................!
C        Find CCR,CCS,FR1,FS1,GR,GS
         expm1=1.D0/EXPO-1.D0                            !
         GR=-AMOD*EXPO/(DABS(FRO/AMOD-1.0D0))**expm1     !............
         GS=-AMOD*EXPO/(DABS(FSO/AMOD-1.0D0))**expm1     !
         FR1=0.5D0*REPT*REPT*H*(FRO+GR*BR)/RR            !
         FS1=0.5D0*REPT*REPT*H*(FSO+GS*BS)/RS            !
         CCR=0.5D0*(CR*RR+BMOD)*GR                       !
         CCS=0.5D0*(CS*RS+BMOD)*GS                       !
C    !...................................................!
C        Gama coefficients for 1st order energy equation
C    !...................................................!
         GTUU=REYCK*RHO*H*DTX-BT*H*DPX                   ! Gama coefficient
     +       -MU0/H*(KXO*(2.D0*U0+MSPEED)-SPEED*KXR)     ! for Uj
     +       -MU0*UC0*(UR*FR1+U0*FS1)-2.D0*MU0*UC1*UR*FR1!
         GTVV=REYCK*RHO*H*DTY-BT*H*DPY                   ! for Vj
     +       -MU0*V0*(2.D0*KXO/H+UC0*(FR1+FS1)           !
     +       +2.D0*UC1*FR1)                              !
         GTM =BMOD/H*(UC0*(GR+GS)/2.D0+UC1*GR)           ! for MUj
         GTR =REYCK*H*(U0*DTX+V0*DTY)                    ! for RHOj
     +       -MU0/RHO*GTM-MU0/RHO/H*(UC0*KXO+UC1*KXR)    !
         GTB =-H*T(J,K)*(U0*DPX+V0*DPY)                  ! for BETAKj
         GTK =0.0D0                                      ! for Thermal conductivity
         GCP =REC*RHO*H*(U0*DTX+V0*DTY)                  ! for Cpj
         GTHH=REYCK*RHO*(U0*DTX+V0*DTY)-BT*(U0*DPX+V0*DPY)
     +       -MSPEED*DPX-MU0/H/H*(UC0*(CCR+CCS)          ! for Hj
     +       +UC1*2.D0*CCR)                              !
c                                                        !
C       !.....................!..........................! NON ADIABATIC CASE
         IF (ISOTH.NE.-1) THEN                           ! Perturbation of
C       !.....................!                          ! heat transfer coef.
C        Heat Transfer To Stator & Journal: Qb0 & Qj0    !
C       !................................................!
          IF (ISOTH.GE.4) THEN     !............! FOR BEARING RADIAL HEAT FLOW
            Den=1.0D0+HCB(J,K)/HKB
          ELSE
            Den=1.0D0
          END IF                   !..........................................!

C##       QB0=HCB(J,K)*(T(J,K)-TB(J,K))/Den/Den          ! removed for simplicity
C##       QJ0=HCJ(J,K)*(T(J,K)-TJ(J,K))                  ! with uniform temps.

          QB0=(HCB(J,K)*(T(J,K)-TBS)/Den)/Den            ! Heat transfer to bearing
          QJ0=HCJ(J,K)*(T(J,K)-TBJ)                      ! and to journal surface

          QS0=QB0+QJ0                                    !QS0 equivalent

c         actual QB0 is  =  HCB (T-TB)/Den
C       !................................................!
C         Modify the gama coefficients due to Hbj and Hjj ! NON ADIABATIC CASE
C       !................................................!  & Turbulent flow

          GTUU=GTUU+REC*2.D0*H*(FS1*U0*QB0/KXS           ! These terms come from
     +        +FR1*UR*QJ0/KXR)                           ! the perturbation of
          GTVV=GTVV+REC*2.D0*H*V0*(FS1*QB0/KXS           ! the heat transfer
     +        +FR1*QJ0/KXR)                              ! coefficients.
          GTM =GTM -REC/MU0*(2.D0*QS0/3.D0               !
     +        +BMOD*(GS*QB0/KXS+GR*QJ0/KXR))             !
          GTR =GTR +REC*(QS0+BMOD*(GS*QB0/KXS            !
     +        +GR*QJ0/KXR))/RHO                          !
          GTK =2.D0*REC*QS0/3.D0/THC(J,K)                ! THC: THERMAL CONDUCTIVITY
          GCP =GCP +REC*QS0/CP/3.D0                      !
          GTHH=GTHH+REC*2.D0*(CCS*QB0/KXS+CCR*QJ0/KXR)/H !

C       !.....................!                          ! NON-ADIABATIC CASE
         END IF                                          !
C       !.....................!                          !

c    !=========================================!..............................!
      ELSE                                               ! LAMINAR FLOW KX=12
c    !=========================================!..............................!
c                                                        ! Laminar flow case
      GTUU=REYCK*RHO*H*DTX-BT*H*DPX                      !
     +    -MU0/H*(2.D0*U0+MSPEED-SPEED)*12.0D0           ! KXR=KXS=12
      GTVV=REYCK*RHO*H*DTY-BT*H*DPY                      !
     +    -MU0*V0*24.0D0/H                               !
      GTM =-(UC0+UC1-SPEED*SPEED/6.0D0)*12.D0/H          ! Based on rotor shear
      GTR =REYCK*H*(U0*DTX+V0*DTY)                       !
      GTB =-H*T(J,K)*(U0*DPX+V0*DPY)                     !
      GTK =0.0D0                                         !
      GCP =REC*RHO*H*(U0*DTX+V0*DTY)                     !
      GTHH=REYCK*RHO*(U0*DTX+V0*DTY)-BT*(U0*DPX+V0*DPY)  !
     +    -MSPEED*DPX+MU0/H/H*(UC0+UC1)*12.0D0           !

C    !.....................!.............................! NONADIABATIC CASES
        IF (ISOTH.NE.-1) THEN                            ! Perturbation of
C    !.....................!                             ! heat transfer coef.
C       Find Heat Transfer To Stator & Journal: Qb0&Qj0  !
C    !...................................................!

          IF (ISOTH.GE.4) THEN     !............! FOR BEARING RADIAL HEAT FLOW
            Den=1.0D0+HCB(J,K)/HKB
          ELSE
            Den=1.0D0
          END IF                   !..........................................!

C##       QB0=HCB(J,K)*(T(J,K)-TB(J,K))/Den/Den          ! removed for simplicity
C##       QJ0=HCJ(J,K)*(T(J,K)-TJ(J,K))                  ! with uniform temps.


          QB0=(HCB(J,K)*(T(J,K)-TBS)/Den)/Den            ! Heat transfer to bearing
          QJ0=HCJ(J,K)*(T(J,K)-TBJ)                      ! and to journal surface
          QS0=QB0+QJ0                                    !

C    !...................................................!
C       Modify the gama coefficients due to Hbj and Hjj
C    !...................................................! Laminar flow case
          GTM =GTM +REC*QS0/MU0/3.D0                     ! These terms come from
          GTK =2.D0*REC*QS0/3.D0/THC(J,K)                ! the perturbation of
          GCP =GCP +REC*QS0/CP/3.D0                      ! the heat transfer
          GTHH=GTHH-REC*QS0/H                            ! coefficients.
                                                         !
C    !.....................!                             !
        END IF                                           !
C    !.....................!.............................!

c    !=========================================!..............................!
      END IF                      ! End of flow condition TURBULENT or LAMINAR
c    !=========================================!..............................!


C    !...................................................!...
      GTPP=GTR*DRHOP(J,K)+GTM*DEMUP(J,K)                 ! Gama coeffs. for Pj
     +    +GTB*DBKP(J,K)+GCP*DCPP(J,K)+GTK*DKP(J,K)      !
c                                                        !...
      GTT =GTR*DRHOT(J,K)+GTM*DEMUT(J,K)                 ! Gamma coeffs for Tj
     +    +GTB*DBKT(J,K)+GCP*DCPT(J,K)+GTK*DKT(J,K)      !
     +    +REC*(HCB(J,K)/Den+HCJ(J,K))                   !
     +    -BETAK(J,K)*H*(U0*DPX+V0*DPY)                  !C    !...................................................!

	  GTU(J, K) =-GTUU*DAREA                         ! These coeffs. are used
          GTV(J, K) =-GTVV*DAREA                         ! in CALCT1 where the 1st
          GTP(J, K) =-GTPP*DAREA                         ! order energy eq. is
	  GTPI(J,K) = BT*H*DAREA                         ! solved.
          GTH(J, K) =-GTHH*DAREA                         !

 	  BT1(J, K) =(BT*U0+MSPEED)*H*DYPK               !
	  BT2(J, K) = BT*V0*H*DX                         !
	  APT(J ,K) = APT(J,K)+GTT*DAREA                 !
	  APTI(J,K) = RHO*H*CP/EC*DAREA                  !
C    !...................................................!

          FWP=FEP                              ! CALCULATE THE PARAMETERS for
	  PW=PE                                ! NEXT STEP--Pw|i=Pe|i-1
          TW=TE
          Uw=U(J,K)
	  PS(J)=PN(J)
	  TS(J)=TN(J)
C..............................!               !
      END DO                                   ! J=JPMIN, JPMAX
C------------------------------!               !

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Calct1                                                      **
C **                                                                         **
C **  CALCT:  Solve the First Order Energy Eqn in Terms of Temperature       **
C **          Ap Tp = Ae Te + Aw Tw + As Ts + An Tn + S                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CALCT1(FP,RES,SIGMA,P1S,Dir,KV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                         --
C ----------------------------------------------------------------------------

      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /COMPLIA/ AC, ETAC, RELAXH, LIFT
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB

      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PARRAY/ P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /TDMA1/A1(MAXNXTP2),B1(MAXNXTP2),C1(MAXNXTP2),D1(MAXNXTP2)
      COMMON /UVP1/ U1(MAXNXT,-MAXNYI:MAXNYI),
     +              V1(MAXNXT,-MAXNYI:MAXNYI),
     +              P1(MAXNXT,-MAXNYI:MAXNYI)

      COMMON /TCOEF/
     + GTU(MAXNXT, -MAXNYI:MAXNYI), GTV(MAXNXT, -MAXNYI:MAXNYI),
     + GTP(MAXNXT, -MAXNYI:MAXNYI), GTH(MAXNXT, -MAXNYI:MAXNYI),
     + GTPI(MAXNXT,-MAXNYI:MAXNYI), APTI(MAXNXT,-MAXNYI:MAXNYI),
     + APT(MAXNXT, -MAXNYI:MAXNYI), AWT(MAXNXT, -MAXNYI:MAXNYI),
     + AET(MAXNXT, -MAXNYI:MAXNYI), AST(MAXNXT, -MAXNYI:MAXNYI),
     + ANT(MAXNXT, -MAXNYI:MAXNYI), BT1(MAXNXT, -MAXNYI:MAXNYI),
     + BT2(MAXNXT, -MAXNYI:MAXNYI)

      COMMON /SOLN/ ISOLN
      COMMON /THERMID/ ISOTH
      COMMON /FLAGS/ INERL, INERP, ITURB, INTER, ICAV, MODEL
      COMMON /KVALB/ K, KM1, KP1
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /PRES/ JPMIN, JPMAX, JPSTART, JPSTOP

      DOUBLE PRECISION REP, REY, MSPEED, SPEED, P,
     +                 SWIRL, PCAV, ALFU, BETU, ALFP,
     +                 AC, ETAC, RELAXH,TBOUT, THERMALK, ROUTER,
     +                 APT,AWT,AET,AST,ANT, APTI,BT1,BT2, HKB,
     +                 GTU,GTV,GTP,GTH,GTPI ,HCB,HCJ,
     +                 ALFT, UC, TC,Ec, SIGMA, RES,
     +                 DXP,DXU,SUE,SUW,DYP,DYV,SVN,SVS ,T,Dir

      INTEGER K, KM1, KP1 ,KV,
     +        JPMIN, JPMAX, JPSTART, JPSTOP, IFULL,
     +        INERL, INERP, ITURB, INTER, ICAV, MODEL,
     +        NPOCKET, NLC, NPC,NLA, NPA, NPAP1, NXI,
     +        NXT, NYI, ISOTH ,ISOLN, LIFT

      DOUBLE COMPLEX A1,B1,C1,D1, T1,U1,V1,P1
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------

      INTEGER J, JJ, JJJ, JJMAX, JP1
      DOUBLE PRECISION FP(*) ,REC, BETAT
      DOUBLE COMPLEX T1W,T1E, U1W,P1W,P1E ,UP1,VP1, SO1T
      DOUBLE COMPLEX P1s(MAXNXT),P1n(MAXNXT) ,IMAG, ACeta

C ----------------------------------------------------------------------------
C --  CALCT1 code                                                           --
C ----------------------------------------------------------------------------
C Added COMPLIANCE EFFECT : Ac/(1+ i ETAC) * GTH * P1
c ............................................................................

      REc=REY/Ec                               ! Reynolds No./Eckert No.
      IMAG=(0.0D0,1.0D0)                       ! imaginary unit
      BETAT=1.0D0-ALFT                         ! relaxT factor
      ACeta=AC/(1.0D0+IMAG*ETAC)               !........................!

      J=JPSTART
      JJJ=2-JPMIN                              !
      JJMAX=JPMAX+JJJ                          !

      T1w=T1(Jpstart,k)
      T1e=T1(Jpstop,k)

      U1W=U1(JPSTART,K)
      P1w=P1(J,K)*SUE(J)+P1(JPMIN,K)*SUW(J)    ! P1w at J=JPMIN
C ......................                       !
      DO J=JPMIN, JPMAX                        !----------------------------
C ......................                       !
      JP1=J+1
      JJ=J+JJJ                                 ! Sweep from left-right on
      P1E=P1(J,K)*SUE(J)+P1(JP1,K)*SUW(J)      ! FIND Pe, Pn on P-C.V.
      P1N(J)=P1(J,KP1)*SVS(KP1)+P1(J,K)*SVN(KP1)

      UP1=(U1W+U1(J,K))/2.D0
      VP1=(V1(J,KV)+V1(J,KP1))/2.D0

      A1(JJ)=-AET(J,K)
      C1(JJ)=-AWT(J,K)
      B1(JJ)=DCMPLX(APT(J,K),APTI(J,K)*Res)/ALFT

      SO1T=GTU(J,K)*UP1+GTV(J,K)*VP1+          ! Source terms
     +(GTP(J,K)+IMAG*SIGMA*GTPI(J,K)+          !
     + GTH(J,K)*ACeta)*P1(J,K)+                !
     +ISOLN*GTH(J,K)*FP(J)                     ! Isoln=0 ==> Homogeneous Eqs.

CNOTE:
C##  ++REc*(HCB(J,K)/(1.0D0+HCB(J,K)/HKB)*TB1(J,K)+
C##  +      HCJ(J,K)*TJ1(J,K) )*DXP(J)*DYP(K)
CNOTE TB1=TJ1=0 perturbed bearing and shaft temperatures  6/1/94

      D1(JJ)=SO1T+BT1(J,K)*(P1E-P1W)+BT2(J,K)
     +      *(P1N(J)-P1S(J))*Dir
      D1(JJ)=D1(JJ)+AST(J,K)*T1(J,KM1)+ANT(J,K)
     +    *T1(J,KP1)+BETAT*B1(JJ)*T1(J,K)

      U1W=U1(J,K)
      P1W=P1E                                  ! NEXT STEP--Vw|i=Ve|i-1
      P1S(J)=P1N(J)

C......................                        !
      END DO                                   !
C ......................                       ! ::::::::::::::::::::::::::::
C                                              ! ----------------------------
  200 CALL TDMACPLX(T1W, T1E, JJMAX)           ! Solve by TDMA algorithm
C                                              ! ----------------------------
C ............................................ !
C                                              ! ::::::::::::::::::::::::::::
  300 DO J=JPMIN, JPMAX                        ! CORRECT:
C    !...................!                     ! no under-relaxation on first
        JJ=J+JJJ                               ! order temperatures.
        T1(J,K)=B1(JJ)                         !
      END DO                                   !
C    !...................!                     !
C                                              !
C ............................................ ! ----------------------------
C                                              !
      IF (IFULL.EQ.1.AND.INTER.EQ.0) THEN      !
        T1(JPSTOP, K)=T1(JPMIN, K)             ! On extended lands
      ENDIF                                    !
c .............................................!.............................

      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Flowt1                                                      **
C **                                                                         **
C **  FLOWT1:  Calculates inlet and outlet flow for first order fields.      **
C **          = INT ( Rho0 V1 H0 + RHO1 V0 H0 + AC/(1+in) V0 P1 ) dT         **
C **          &  {Rho1 H0 + Rho0 Ac/(1+in) P1} dx dy on recess border        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FLOWT1(QABS,SIGMA,FP,FU)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /XYVEC/ XP(MAXNXT), XU(MAXNXT),
     +               YP(-MAXNYI:MAXNYI),YV(-MAXNYI:MAXNYI)
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT),SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /HFILM/ HP(MAXNXT,-MAXNYI:MAXNYI),
     +               HU(MAXNXT,-MAXNYI:MAXNYI),
     +               HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /UVARRAY/ U(MAXNXT,-MAXNYI:MAXNYI),
     +                 V(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/  T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/  RHOP(MAXNXT,-MAXNYI:MAXNYI),
     +                 EMUP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP1/ CK(MAXNXT,-MAXNYI:MAXNYI),
     +            BETAK(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DRho/ Drhop(MAXNXT,-MAXNYI:MAXNYI),
     +              Drhot(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DEmu/ Demup(MAXNXT,-MAXNYI:MAXNYI),
     +              Demut(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /DCPS/ DCPP(MAXNXT,-MAXNYI:MAXNYI),
     +              DCPT(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RECES/ PREC(MAXNPOCK), TREC(MAXNPOCK), QREC(MAXNPOCK),
     +               QIN, QOUT, QFACTOR
      COMMON /URECJET/ UREC(MAXNPOCK,MAXNPOCK+2)

      COMMON /FLOW01/ QO1(MAXNPOCK), QINO1, QOUTO1
      COMMON /QSIDE0H/ QSIDEH(MAXNPOCK)
      COMMON /TORREC/ TORR(MAXNPOCK)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /FACTOR2/ KLOSXu,KLOSXd,KLOSYl,KLOSYr, RENC,ASPEC, HRECD
      COMMON /RECPAR/ HREC, VSUP, BETA
      COMMON /RECASP/ ASPE(MAXNPOCK)
      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HRECC
      COMMON /PARAM2/ EXO, EYO
      common /PARAM4/ CINLET, CEXIT
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      COMMON /COMPLIA/ AC, ETA, RELAXH, LIFT

      COMMON /UVP1/ U1(MAXNXT, -MAXNYI:MAXNYI),
     +              V1(MAXNXT, -MAXNYI:MAXNYI),
     +              P1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /T1ARRAY/  T1(MAXNXT, -MAXNYI:MAXNYI)
      COMMON /RECES1/ PREC1(MAXNPOCK), TREC1(MAXNPOCK),
     +                QREC1(MAXNPOCK), QIN1, QOUT1
      COMMON /S1/ S11(MAXNPOCK)

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /SOLN/ ISOLN
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

c     ................................................................
      DOUBLE PRECISION DXP, DXU, HP, HU, HV, SUW, SUE, DCPP,DCPT,
     +                 DYP, DYV, SVN, SVS, U, V, T, Drhot,Demut,
     +                 XP,XU,YP,YV,
     +                 Rhop,Emup,CK,BETAK, Drhop,Demup,
     +                 PREC, QREC,TREC, QIN, QOUT, QFACTOR,UREC,
     +                 QO1, QINO1, QOUTO1 , QSIDEH, TORR,
     +                 REP, REY, MSPEED, SPEED,
     +                 SWIRL, PCAV, ALFU, BETU, ALFP,
     +                 KLOSXu,KLOSXd,KLOSYl,KLOSYr,RENC,ASPEC,HRECD,
     +                 HREC, VSUP, BETA, ASPE,
     +                 CLEAR, DIAM, LENGTH, LD, AR, HRECC,
     +                 EXO, EYO, ALFT,UC,TC,EC,
     +                 CINLET, CEXIT, AC, ETA, RELAXH

      DOUBLE COMPLEX U1,V1,P1,T1, PREC1,TREC1,QREC1,S11,QIN1,QOUT1
      INTEGER ISOLN,NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,
     +        ICSTEP, IFULL, BEARING, LIFT
C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION QABS,rho,svsk,svnk,factor,Dir,FPR,FP(*),FU(*),
     +                 TORC,DTRR,DTRM,DTRH, XR,HR,VR,SIGMA, BR,
     +                 CPR,RHOR,EMUR,DUMYR,URECave
      DOUBLE COMPLEX   Ss, Qr,rho1, Imag,ZERO, RHOR1,EMUR1,
     +                 CPR1,DUMYI, ACeta, PN1, PU1,
     +                 QL1T0(MAXNPOCK),QR1T0(MAXNPOCK),QSIDE1(MAXNPOCK)
      INTEGER I,J,JP1,K, JMIN,JMAX, NYIM1, Kstep, Kstart,Kend,Icsym,IJ

C ----------------------------------------------------------------------------
C --  FLOWT1 code: calculates                                               --
c     Qrec1= INT [{Ho(RHO1 Vo+RHOo V1)+RHOo Vo ACeta P1}.n dTr]
c     on (U,V) grid recess boundary
c     S11= INT { (RHO1 HO+RHOo ACeta P1) dArea} on area between recess and grid
C ----------------------------------------------------------------------------
      Nyim1=Nyi-1
      Imag=(0.0D0, 1.0D0)
      ZERO=(0.0D0, 0.0D0)
      QIN1=ZERO
      QOUT1=QIN1
      QABS=0.0D0
      ACeta=AC/DCMPLX(1.0D0,ETA)              ! COMPLIANCE Ac/(1+in)
      DTRR=0.0D0                              ! Perturbation coeffs. of
      DTRM=0.0D0                              ! temperature due to recess
      DTRH=0.0D0                              ! torque

      Kstep=1
      Icsym=1
      Dir=1.0D0

c    !.......................!
      IF (NPOCKET.GT.0) THEN
c    !.......................!
          DO I=1, NPOCKET
             QREC1(I)=QIN1
             S11(I)=QIN1
             QL1T0(I)=ZERO                    ! 1st-order Inlet left flow * To
             QR1T0(I)=ZERO                    !  `'     Outlet right flow * To
             QSIDE1(I)=ZERO                   ! 1st order top flow (with hj)
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
          QR=ZERO                             ! Inlet flow
          Ss=ZERO

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

          QSIDE1(I)=QSIDE1(I)+QR              ! 1st order top flow (with hj)

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

          DUMYI=Dyp(k)*(Hu(j,k)*(rho*U1(j,k)+rho1*U(j,k) ) +
     +                  rho*U(j,k)*ACeta*PU1 )

          DUMYR=DYP(K)*RHO*U(J,K)
          QR=QR-DUMYI

          QL1T0(I)=QL1T0(I)+(DUMYI             ! 1st-order Inlet left flow * To
     +            +DUMYR*ISOLN*FU(J))*T(JP1,K) !

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

          DUMYI=Dyp(k)*( Hu(j,k)*(rho*U1(j,k)+rho1*U(j,k)) +
     +                   rho*U(j,k)*ACeta*PU1 )
          DUMYR=DYP(K)*RHO*U(J,K)
          QR=QR+DUMYI
          QR1T0(I)=QR1T0(I)+(DUMYI            ! 1st-order outlet flow * To
     +            +DUMYR*ISOLN*FU(J))*T(J,K)  !

          rho1=Drhop(j,k)*P1(j,k)+Drhot(j,k)*T1(j,k)
          Ss=Ss+Dxp(j)*Dyp(k)/2.0D0/Factor*
     +         (rho1*Hp(j,k)+rhop(j,k)*ACeta*P1(j,k))
C        !................!                   !
          END DO                              !.....................!
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

C At pad inlet/outlet: no compliance effect since P=Pback here, i.e
C        there is no pressure deformation P1=0


 666  CONTINUE

C    !---------------------!
      IF (ISYM.EQ.0) THEN
c    !.....................! ASYMMETRIC BEARING
       IF (BEARING.eq.2) GOTO 99

       IF(Icsym.eq.2) goto 98

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
            QREC1(I)= Factor*QREC1(I)
            QSIDE1(I)=Factor*QSIDE1(I)
            QL1T0(I)= Factor*QL1T0(I)
            QR1T0(I)= Factor*QR1T0(I)
            S11(I)=Factor*S11(I)
           END DO
         END IF

         QIN1=Factor*QIN1
         QOUT1=Factor*QOUT1

c    !.......................!
      END IF
c    !.......................!

 98   CONTINUE


C    !........................................! 1st order recess energy balance
      IF (NPOCKET.GT.0) THEN                  ! FOR HYDROSTATIC POCKET BEARINGS
C    !........................................!
C##   NOTE: No bearing compliance effects accounted : COMPLIANCE

         Kstart=ISYM+(ISYM-1)*(NPA-1)         !
c    !.......................!
         DO I=1, NPOCKET                      !...................
c    !.......................!                ! NXI=NLC+NPC-2
          JMIN=(I-1)*NXI+NLC                  ! JMAX-JMIN=NPC-1
          JMAX=I*NXI+1                        !

          BR=ASPE(I)*DIAM                     ! Cir. length of recess
          XR=(XP(JMAX)+XP(JMIN))/2.0D0        ! Find recess volume
          HR=CINLET*(1.0D0+EXO*DCOS(XR)+EYO*DSIN(XR))  !### no misalignment
          VR=(AR*BR*(HREC+HR)+VSUP)           ! Half of the dimensionless
     +       *4.0D0/DIAM/DIAM/CLEAR           ! recess volume

          JP1=JMIN+1
          CPR=CK(JP1,Kstep)                   ! fluid Properties in recesses
          RHOR=RHOP(JP1,Kstep)                ! Note: both pressure
          EMUR=EMUP(JP1,Kstep)                ! and temperature in
          RHOR1=DRHOP(JP1,Kstep)*PREC1(I)     ! the recess are considered
     +         +DRHOT(JP1,Kstep)*TREC1(I)     ! uniform
          CPR1 =DCPP(JP1 ,Kstep)*PREC1(I)     !
     +         +DCPT(JP1 ,Kstep)*TREC1(I)     !
          EMUR1=DEMUP(JP1,Kstep)*PREC1(I)     !
     +         +DEMUT(JP1,Kstep)*TREC1(I)     !
c -  base analysis on average recess speed:   ! 10/4/95 by LSA
          URECave=0.0D0                       !
          DO IJ=1,NPC-1                       !
            URECave=URECave+UREC(I,IJ)        !
          END DO                              !
          URECave=URECave/(NPC-1)             !
          TORC=URECave*EC/REY                 !...... 10/4/95
c                                             !
          HR=(HR+HREC)/CLEAR                  !
c                                             !
          FPR=(FP(JP1)+FP(JMAX-1))*ISOLN/2.D0 !
C        !....................................!
          CALL TORRECJ(URECave,HR,SPEED,REP,RHOR,
     +                 EMUR,DTRR,DTRM,DTRH)   ! Include the recess temp-rise
C        !....................................!
          DUMYR=VR*RHOR*SIGMA                 ! Recess 1st order temperature
          TREC1(I)=(QREC1(I)+ISOLN*QO1(I)     ! is uniform
     +    -(QSIDE1(I)+ISOLN*QSIDEH(I))        !
     +     *TREC(I)+QL1T0(I)-QR1T0(I)         !
     +    +(DTRR*RHOR1/RHOR+DTRM*EMUR1/EMUR   ! Note: All the recess
     +     +DTRH*FPR/HR-CPR1/CPR)*TORR(I)     ! variables are ontained from   !######
     +     *TORC/CPR-imag*VR*SIGMA*TREC(I)    ! the whole recess, not
     +     *RHOR1)/DCMPLX(QREC(I),DUMYR)      ! half of it, even symm.

          DO J=JMIN, JMAX                     ! 1st order Temps at top edge
             T1(J,NPA)=TREC1(I)               ! of recess
          END DO                              !...
          IF (ISYM.EQ.0) THEN                 !
           DO J=JMIN, JMAX                    ! 1st order Temps at bottom edge
             T1(J,-NPA)=TREC1(I)              ! of recess
           END DO                             !
          END IF                              !....
          DO K=Kstart, NPA-1                  ! 1st order Temps at side edges
             T1(JMIN,K)=TREC1(I)              ! of recess
             T1(JMAX,K)=TREC1(I)              !
          END DO                              ! END OF
C    !........................................! 1st order recess energy balance
         END DO                               ! NEXT Recess, I=1,2,...NPOCKET
C    !........................................!

c    !.......................!                !
      END IF
c    !.......................!                ! ONLY FOR HYDROSTATIC BEARING

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
        QABS=CDABS(QIN1)        !
c    !..........................!
      END IF
c    !..........................!

      END

C    !........................................!


C *****************************************************************************
C **                                                                         **
C **  Subroutine Torrecj                                                     **
C **                                                                         **
C **  TORRECJ: Calculates shear factors Ksj and coeffs. for recess temp-rise **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE TORRECJ(U,H,SP,REPT,RHO,EMU,DTRR,DTRM,DTRH)

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /MOODY/ AMOD, BMOD, RUGR, RUGS, EXPO

      DOUBLE PRECISION AMOD, BMOD, RUGR, RUGS, EXPO

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION U,H,SP,REPT,rho,emu,DTRR,DTRM,DTRH,REP,
     +                 RS,CS,BS,FSO,KSO,DUMY,GS,CCS,expm1
C ----------------------------------------------------------------------------
C --  TORRECJ code                                                          --
C ----------------------------------------------------------------------------
C no roughness effects on recess surface, i.e. smooth surface
c no axial flow velocity

      REP=REPT*RHO/EMU                           !
      RS=DABS(H*U*REP)                           ! Stator Reynolds number
      CS=0.0d0 ! RUGS*10000.0D0/H                ! Stator roughness coeff.
      BS=BMOD/RS                                 !
      FSO=AMOD*(1.0D0+(CS+BS)**EXPO)             ! Stator friction factor
      KSO=RS*FSO                                 ! Zeroth turbulent shear coefs
      DUMY=12.0D0                                !
      KSO=DMAX1(DUMY, KSO)                       ! Select Laminar or turbulent

      IF (KSO.GT.DUMY) THEN                      ! Turbulent flow condition
          expm1=1.0D0/expo-1.0D0                 !
          GS=-AMOD*EXPO/(DABS(FSO/AMOD-1.0D0))**expm1
          CCS=GS*(RS*CS+BMOD)/2.0D0

          DTRM=-GS*BMOD/KSO                      ! Coeffs. of emuj
          DTRR=1.0D0-DTRM                        ! Coeffs. of rhoj
          DTRH=2.0D0*CCS/KSO			 ! Coeffs. of Hj
          RETURN                                 !
      END IF                                     !.....................

      DTRM= 1.0D0                                ! Laminar flow case
      DTRR= 0.0D0
      DTRH=-1.0D0

      RETURN

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine PROPERTIES                                                  **
C **                                                                         **
C **  PFIELDS:  Print entire Bt, Cp, Mu, Rho and K fields.  This routine     **
C **            is added for debugging purposes.                             **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PROPERTIES(DEVICE)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /PARAM2/ EXO, EYO
      COMMON /PROP1/  CK(MAXNXT,-MAXNYI:MAXNYI),
     +              BETA(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP2/ HCB(MAXNXT,-MAXNYI:MAXNYI),
     +               HCJ(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /PROP3/ THC(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /RHOEMU/RHP(MAXNXT,-MAXNYI:MAXNYI),
     +               EMP(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /TARRAY/ T(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /THERMAL/ ALFT, UC, TC,Ec

      INTEGER NPOCKET, NLC, NPC,
     +                NLA, NPA, NPAP1, NXI, NYI, NXT ,IFULL,NXX
      DOUBLE PRECISION CK,BETA,HCB,HCJ,THC,RHP,EMP,
     +                 T, UC, TC,Ec, ALFT ,EXO,EYO

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER DEVICE, I, J, Ifile, U,  WID

C ----------------------------------------------------------------------------
C --  PROPERTIES code                                                       --
C ----------------------------------------------------------------------------
      Ifile=1

      IF (NYI.LE.6) THEN                         !
          WID=10                                 ! One of the few bits of
      ELSE IF (NYI.LE.7) THEN                    ! cleverness FORTRAN allows...
          WID=9                                  !
      ELSE IF (NYI.LE.8) THEN                    !
          WID=8                                  !
      ELSE                                       !
          WID=6                                  !
      END IF

      IF ((EXO.eq.(0.0)).AND.(EYO.eq.(0.0))) THEN
         NXX=NXI
      ELSE
         NXX=NXT
      END IF
      U=6                                        ! Unit=6 (screen)

  17  CONTINUE
      WRITE (U,410)
 410  FORMAT (' ','FLUID Properties in dimensionless form',/,
     +        ' ',70('.'))
      WRITE (U, 40)                              !
      DO I=1, NXX                                !
          WRITE (U, 10) I, (BETA(I, J), J=1, NYI)! compressibility parameter
      END DO                                     !
      WRITE (U, 50)                              !
      DO I=1, NXX                                !
          WRITE (U, 10) I, (CK(I, J), J=1, NYI)  ! specific heat
      END DO                                     !
      WRITE (U, 60)                              !
      DO I=1, NXX                                !
          WRITE (U, 10) I, (EMP(I, J), J=1, NYI) ! viscosity
      END DO                                     !
      WRITE (U, 70)                              !
      DO I=1, NXX                                !
          WRITE (U, 10) I, (RHP(I, J), J=1, NYI) !density
      END DO                                     !
      WRITE (U, 80)                              !
      DO I=1, NXX                                !
          WRITE (U, 10) I, (THC(I, J), J=1, NYI) !thermal conductivity
      END DO                                     !
      WRITE (U, 90)                              !
      DO I=1, NXX                                !convection
          WRITE (U, 10) I, (HCB(I, J), J=1, NYI) !heat transfer coeff.
      END DO                                     !to bearing
      WRITE (U, 100)                             !
      DO I=1, NXX                                !convection
          WRITE (U, 10) I, (HCJ(I, J), J=1, NYI) !heat transfer coeff.
      END DO                                     !to journal
c                                                !
c  !.............................................!.................
      IF (Ifile.eq.2) RETURN                     !
      IF (DEVICE.EQ.1) THEN                      !
         U=1                                     ! Unit=1 (DUMP file)
         Ifile=2                                 !
         GOTO 17                   !=> print     !
      END IF                                     !
c  !.............................................!.................
      RETURN
c  !.............................................!.................

   10 FORMAT (' ', I2, ' : ', 9(F7.5, ' '))
c   10 FORMAT (' ', I2, ' : ', 9(F<WID>.<WID-2>, ' '))
C   10 FORMAT (' ', I2, ':',9(E12.5E2,' '))

   40 FORMAT (' ', 'Beta Array: Dimensionless Fluid Volumetric',
     +      /,' ', '            Expansion Coefficient')

   50 FORMAT (' ', 'Cp Array: Dimensionless Specific Heat')

   60 FORMAT (' ', 'Mu Array: Dimensionless Viscosity')

   70 FORMAT (' ', 'Rho Array: Dimensionless Density')

   80 FORMAT (' ', 'K Array: Dimensionless Fluid heat Conductivity')

   90 FORMAT (' ', 'Hcb Array: Dimensionless Turbulent Heat Transfer',
     +      /,' ', '           Coefficient--to Bearing Surface')

  100 FORMAT (' ', 'Hcj Array: Dimensionless Turbulent Heat Transfer',
     +      /,' ', '           Coefficient--to Journal Surface')

      END

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Last revised 10/12/95 by Dr. Luis San Andres TexasA&M University
c
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::