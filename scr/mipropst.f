c #    #     #    #####   #####    ####   #####    ####    #####          ######
c ##  ##     #    #    #  #    #  #    #  #    #  #          #            #
c # ## #     #    #    #  #    #  #    #  #    #   ####      #            #####
c #    #     #    #####   #####   #    #  #####        #     #     ###    #
c #    #     #    #       #   #   #    #  #       #    #     #     ###    #
c #    #     #    #       #    #   ####   #        ####      #     ###    #
c
c for hydroflex.f program  1994
c
c  SEARCH for CVAX to find changes neeeded for VAX SYSTEMS
c  FUNCTION CPI uses triple precision
c
C------------------------------------------------------------------------|
c      PROGRAM MIPROPS                                                   |
c -----------------------------------------------------------------------+
c MIPROPS program has been modified by Dr. Luis SanAndres on 03/90       |
c         for implicit/active use in compressible liquid, hybrid         |
c         bearing calculations.                                          |
c         This version calculates Density, Viscosity, specific heat      |
c         and sound speed GIVEN the absolute liquid pressure & Temp      |
c         ------- INPUT and OUTPUT are in SI (metric units)              |
c------------------------------------------------------------------------+ 
c 4/17/90 Added FUNCTION PTVISC for evaluation of viscosity of LH2       |
c........................................................................+
c
c
c            **************************************************'
c            *       NBS STANDARD REFERENCE DATABASE 12       *'   
c            *      THERMOPHYSICAL PROPERTIES OF FLUIDS       *'
c            *                  MIPROPS 1986                  *'
c            *                                                *'
c            *                   WRITTEN BY                   *'
c            *               ROBERT D. McCARTY                *'
c            **************************************************')
c            *             THERMOPHYSICS DIVISION             *'
c            *        CENTER FOR CHEMICAL ENGINEERING         *'
c            *          NATIONAL BUREAU OF STANDARDS          *'
c            *             BOULDER,COLORADO 80303             *'
c            *                                                *'
c            *                 DISTRIBUTED BY                 *'
c            *       OFFICE OF STANDARD REFERENCE DATA        *'
c            *          NATIONAL BUREAU OF STANDARDS          *'
c            *          GAITHERSBURG, MARYLAND 20899          *'
c            *                                                *'
c            *               COPYRIGHT 1986 BY                *'
c            *           U.S. SECRETARY OF COMMERCE           *'
c            *   ON BEHALF OF THE UNITED STATES GOVERNMENT    *'
c            **************************************************')

c   SELECTED  FLUIDs for HYDRO program are
c  
c            PARA HYDROGEN,  Choice IF=1, DATA File: PH2.COF
c            NITROGEN,       CHoice IF=2,   "" ""    N2.COF
c            OXYGEN,         Choice IF=3,   "" ""    O2.COF
c            METHANE,        Choice IF=4,   "" ""    METH.COF
c
c            OTHER,          IF=12
c                            Uses linear interpolation supply&sump




C....................................................................
C
      SUBROUTINE REPRO(P,D,T,Dk,Vk,W,TH,CPP)
c....................................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/CONT/IF
      COMMON/CRIT/EM,EOK,RM,TC,DC,X,PC,SIG
      COMMON/CPID/GI(11),GH(11),GL(11)
      COMMON/DIEL/ BX(6),PX(6)
      COMMON/HAN/CR,TCI 
c...................................................................

      Dk=D*EM                        ! Density in Kg/m3
      W=SOUND(D,T)                   ! Calculate sonic speed
      CPP=CP(D,T)*1000.D0/EM         ! Specific Heat Cp in J/Kg-K
c
c     Calculate  fluid Viscosity
c     ----------------------------

      GO TO (130,150,150,150) , IF

      return

c...................................................
c    FOR: IF= 1 LH2
c..................................................

  130 Vk=PTVISC(P,T)                 ! Viscosity in PA-SE+6 of H2  
      TH=0.09174D0   ! APPROXIMATELY ! Conductivity in W/M-K of H2
      return

c....................................................
c     FOR: IF=2 N2, IF=3 O2, IF=4 Methane
c.....................................................
     
  150 VK=VISC(D,T)/1.0D+6            ! Conductivity in W/M-K
      TH=THERM(D,T)       
      return

      END 


c...............................................................

      SUBROUTINE FDATA(IF)
c.....................................................................
c READs DATA file containing DATA for Fluid TYPE: IF

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION G(32),VP(9),GV(9),GT(9),FV(4),FT(4),EV(8),ET(8),A(20)
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/SEN/BETA,XO,DELTA,E1,E2,AGAM
      COMMON/CRIT/EM,EOK,RM,TC,DC,X,PC,SIG
      COMMON/SATC/A,DTPV,EG
      COMMON/CPID/GI(11),GH(11),GL(11)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      COMMON/ISP/N,NW,NWW
      COMMON/FIXPT/T0,S0,H0
      COMMON/DIEL/BX(6),PX(6)
      COMMON/B/GMW,SEOK,SSIG
      COMMON/A/SE1,G1,B1,DE,BK,D1,XZ,ZZ,X1,X2,X3,X4

c....................................................................

      WRITE(6, 1201)
 1201 FORMAT(' ',/,
     1 '     **************************************************'
     1/'     *       NBS STANDARD REFERENCE DATABASE 12       *'
     4/'     *      THERMOPHYSICAL PROPERTIES OF FLUIDS       *'
     5/'     *                  MIPROPS 1986                  *'
     6/'     *                                                *'
     7/'     *                   WRITTEN BY                   *'
     9/'     *               ROBERT D. McCARTY                *'
     +/'     **************************************************')
      

      GO TO (3,5,7,4,6,1,8,2,9,10,11,12),IF
c....................................................................
c IF=6, ARGON

    1 WRITE(6,200)
  200 FORMAT(' ',40('-'),/,'IF=6 -> FLUID: ARGON'/,1X, 40('-'),
     1/' THE TEMPERATURE RANGE FOR ARGON IS 83.8 TO 400 K'
     1/' (-308.8 TO 260 F) WITH PRESSURES TO 100 MPA (14504 PSIA)')
      OPEN(44,FILE='ARGON.COF',STATUS='UNKNOWN')
      N=0
      NW=0
      NWW=0
      GO TO 50
c...................................................................
c IF=8, ETHYLENE

    2 WRITE(6,201)
  201 FORMAT(' ',40('-'),/,'IF=8 -> FLUID: ETHYLENE',/,1X,40('-'),
     1/ ' THE TEMPERATURE RANGE FOR ETHYLENE IS 104 TO 400 K'
     1/' (-272.4 TO 260 F) WITH PRESSURES TO 40 MPA (5801 PSIA)')
      OPEN(44,FILE='C2H4.COF',STATUS='UNKNOWN')
      N=0
      NW=0
      NWW=0
      GO TO 50
c..................................................................
c IF=1, paraHYDROGEN

    3 WRITE(6,202)
  202 FORMAT(' ',40('-'),/,'IF=1 -> FLUID:HYROGEN',/,1X,40('-'),
     1/' THE TEMPERATURE RANGE FOR HYDROGEN IS 13.8 TO 400 K'
     1/' (-434.8 TO 260 F) WITH PRESSURES TO 120 MPA (17404 PSIA)')
      OPEN(44,FILE='PH2.COF',STATUS='UNKNOWN')
      N=1
      NW=0
      NWW=0
      GO TO 50
c....................................................................
c IF=4, METHANE

    4 WRITE(6,203)
  203 FORMAT(' ',40('-'),/,'IF=4 -> FLUID: METHANE',/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR METHANE IS 90.68 TO 600 K'
     1/' (-296.45 TO 620 F) WITH PRESSURES TO 200 MPA (29007 PSIA)')
      OPEN(44,FILE='METH.COF',STATUS='UNKNOWN')
      N=0
      NW=1
      NWW=0
      GO TO 50
c...................................................................
c IF=2 NITROGEN

    5 WRITE(6,204)
  204 FORMAT(' ',40('-'),/,'IF=2 -> FLUID:NITROGEN',/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR NITROGEN IS 63.15 TO 1900 K'
     1/' (-346 TO 2960 F) WITH PRESSURES TO 1000 MPA (145034 PSIA)')
      OPEN(44,FILE='N2.COF',STATUS='UNKNOWN')
      N=0
      NW=0
      NWW=0
      GO TO 50
c.....................................................................
c IF=5 NITROGEN TRIFLUORIDE

    6 WRITE(6,205)
  205 FORMAT(' ',40('-'),/,'IF=5 -> FLUID:NITROGEN FLUORIDE',
     1/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR NITROGEN TRIFLUORIDE IS'
     1,/' 66.36 TO 500 K (-340.2 T0 440 F) AND  50 MPA (7252 PSIA)')
      OPEN(44,FILE='NF3.COF',STATUS='UNKNOWN')
      N=0
      NW=0
      NWW=0
      GO TO 50
C.....................................................................
C IF=3 OXYGEN

    7 WRITE(6,206)
  206 FORMAT(' ',40('-'),/,'IF=3 -> FLUID:OXYGEN',/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR OXYGEN IS 54.359 TO 400 K'
     1/' (-361.8 TO 260 F) WITH PRESSURES TO 120 MPA (17404 PSIA)')
      OPEN(44,FILE='O2.COF',STATUS='UNKNOWN')
      N=0
      NW=0
      NWW=0
      GO TO 50
C....................................................................
C IF=7 ETHANE

    8 WRITE(6,207)
  207 FORMAT(' ',40('-'),/,'IF=7 -> FLUID:ETHANE',/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR ETHANE IS 90.35 TO 600 K'
     1/' (-297 TO 620 F) WITH PRESSURES TO 70 MPA (10153 PSIA)')
      OPEN(44,FILE='C2H6.COF',STATUS='UNKNOWN')
      N=0
      NW=1
      NWW=0
      GO TO 50
C....................................................................
C IF=9 PROPANE

    9 WRITE(6,208)
  208 FORMAT(' ',40('-'),/,'IF=9 -> FLUID:PROPANE',/,1X,40('-'),
     +/' THE RANGE OF TEMPERATURE FOR PROPANE IS 85.47 TO 600 K'
     1/' (-305.8 TO 620 F) WITH PRESSURES TO 100 MPA (14504 PSIA)')
      OPEN(44,FILE='C3H8.COF',STATUS='UNKNOWN')
      N=0
      NW=1
      NWW=0
      GO TO 50
C....................................................................
C IF=10 ISOBUTANE

   10 WRITE(6,209)
  209 FORMAT(' ',40('-'),/,'IF=10 -> FLUID:ISOBUTANE',/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR ISO BUTANE IS'
     1/' 113.55 TO 600 K (-255.2 TO 620 F)'
     1/' WITH PRESSURES TO 35 MPA (5076 PSIA)')
      OPEN(44,FILE='ISOB.COF',STATUS='UNKNOWN')
      N=0
      NW=1
      NWW=0
      GO TO 50
C....................................................................
C IF=11 NORMAL BUTANE

   11 WRITE(6,210)
  210 FORMAT(' ',40('-'),/,'IF=11 -> FLUID: NORMAL BUTANE',
     1/,1X,40('-'),
     1/' THE RANGE OF TEMPERATURE FOR NORMAL BUTANE IS'
     1/' 134.68 TO 500 K (-216.9 TO 440 F)'
     1/' WITH PRESSURES TO 70 MPA (10153 PSIA)')
      OPEN(44,FILE='NORB.COF',STATUS='UNKNOWN')
      N=0
      NW=1
      NWW=0
      GO TO 50
C....................................................................
C IF=12 STOP EXECUTION

   12 STOP
c....................................................................


   50 READ(44,100)SIG,XO,BETA,DELTA,E1,E2,AGAM
      READ(44,100)EM,EOK,RM,TC,DC,X,PC
  100 FORMAT(7E12.6)
      DO 60 I=1,32
   60 READ(44,101)G(I)
  101 FORMAT(3D20.13)
      DO 70 I=1,20
   70 READ(44,101)A(I)
      DO 75 I=1,11
   75 READ(44,101)GI(I),GH(I),GL(I)
      DO 80 I=1,9
   80 READ(44,101)VP(I),GV(I),GT(I)
      TCC=VP(8)
      PTP=VP(9)
      TTP=VP(7)
      DO 90 I=1,8
   90 READ(44,101)EV(I),ET(I)
      DO 110 I=1,4
  110 READ(44,101)FV(I),FT(I)
      READ(44,101)DTP,DTPV
      READ(44,101)T0,S0,H0
      READ(44,102)R,GAMMA,TUL,TLL,PUL,DCC,PCC
      DO 95  I=1,6
   95 READ(44,101)BX(I),PX(I)
      IF(NW.EQ.0)GO TO 96
      READ(44,100)SE1,G1,B1,DE,BK,D1
      READ(44,100)XZ,ZZ,X1,X2,X3,X4
      READ(44,100)GMW,SEOK,SSIG,EG
      EM=GMW
   96 CONTINUE
c
      CLOSE(44,STATUS='KEEP')


      TTPF=(TTP-273.15D0)*1.8D0+32.D0
      TCCF=(TCC-273.15D0)*1.8D0+32.D0
      PCCP=PCC*14.6959D0/.101325D0
      WRITE(6,1065) TCC,TCCF
 1065 FORMAT (/' CRITICAL TEMPERATURE IS ',F7.2,' K  (',F8.2,' F).')
      WRITE(6,1075) PCC,PCCP
 1075 FORMAT (' CRITICAL PRESSURE IS   ',F6.3,' MPA (',F6.1,' PSIA).'/)
  102 FORMAT(F10.8,E14.8,3F8.2,2F8.4)

      RETURN
      END

c----------------------------------------------------------------------

      SUBROUTINE LIMITS(P,T,IL)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      IF(P .GT. PUL)GO TO 10
      IF(T .GT. TUL.OR.T .LT. TLL)GO TO 12
      PM=PMELT(T) 
      IF(P .GT. PM) GO TO 20
      IL=1
      RETURN

   10 PULF=(PUL/.101325D0)*14.6959D0
      WRITE(*,11)PUL,PULF
   11 FORMAT(' THE INPUT PRESSURE IS OUT OF THE RANGE OF THIS EQUATION '
     1/' THE PRESSURE MUST BE BELOW ',F6.0,' MPA OR ',F7.0,' PSIA')
      write (6,*) 'Liquid Pressure is:', P, ' MPA'
      IL=0
      RETURN

   12 TLLF= (TLL-273.15D0)*1.8D0+32.D0
      TULF= (TUL-273.15D0)*1.8D0+32.D0
      WRITE(*,13) TLL,TUL,TLLF,TULF 
   13 FORMAT(' THE INPUT TEMPERATURE IS OUT OF RANGE' 
     A /' THE RANGE FOR THIS EQUATION IS  ',F6.2,' K   TO ',F6.0,' K',/,
     B 27X,' OR ',F8.2,' F   TO ',F6.0,' F')
      write (6,*) 'Liquid Temperature is:', T, ' K'
      IL=0
      RETURN

   20 TM=TMELT(P) 
      TF=(TM-273.15D0)*1.8D0+32.D0
      WRITE(*,22) TM,TF 
   22 FORMAT(' SOLID PHASE DETECTED.',/,' FOR THIS PRESSURE,  TEMP' 
     A ' SHOULD EXCEED ',F8.3,' K, OR',F9.3,' F') 
      IL=0
      END 

c-----------------------------------------------------------
      SUBROUTINE PROPS(PP,DD,TT,K)

c..........................................................
C THE 32 TERM EQUATION OF STATE, INPUT IS DENSITY(MOLES/L),
C TEMPERATURE(K), OUTPUT (PP) IS PRESSURE(MPA),OR DP/DD IN
C LITER-MPA/MOLE OR DP/DT MPA/K OR S,H,OR CV AT ONE LIMIT OF
C INTEGRATION

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION X(33) 
      DIMENSION B(33),G(32),VP(9) 
      EQUIVALENCE (B,X) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      DATA M/32/
      D=DD
      P=PP
      T=TT
      GM=GAMMA
      D2=D*D
      D3=D2*D 
      D4=D3*D 
      D5=D4*D 
      D6=D5*D 
      D7=D6*D 
      D8=D7*D 
      D9=D8*D 
      D10=D9*D
      D11=D10*D 
      D12=D11*D 
      D13=D12*D 
      TS=DSQRT (T) 
      T2=T*T
      T3=T2*T 
      T4=T3*T 
      T5=T4*T 
      F=DEXP (GM*D2) 
      GO TO (100,200,300,400,500,600),K
C     ENTRY PRESS
  100 B( 1)=D2*T
      B( 2)=D2*TS 
      B( 3)=D2
      B( 4)=D2/T
      B( 5)=D2/T2 
      B( 6)=D3*T
      B( 7)=D3
      B( 8)=D3/T
      B( 9)=D3/T2 
      B(10)=D4*T
      B(11)=D4
      B(12)=D4/T
      B(13)=D5
      B(14)=D6/T
      B(15)=D6/T2 
      B(16)=D7/T
      B(17)=D8/T
      B(18)=D8/T2 
      B(19)=D9/T2 
      B(20)=D3*F/T2 
      B(21)=D3*F/T3 
      B(22)=D5*F/T2 
      B(23)=D5*F/T4 
      B(24)=D7*F/T2 
      B(25)=D7*F/T3 
      B(26)=D9*F/T2 
      B(27)=D9*F/T4 
      B(28)=D11*F/T2
      B(29)=D11*F/T3
      B(30)=D13*F/T2
      B(31)=D13*F/T3
      B(32)=D13*F/T4
      P=0.D0 
      DO 101 I=1,M
      P=P+B(I)*G(I) 
  101 CONTINUE
      P=P+R*D*T 
      PP=P
      RETURN

C     ENTRY DPDD
  200 F1=2.D0*F*GM*D
      F21=3.D0*F*D2 +F1*D3 
      F22=5.D0*F*D4 +F1*D5 
      F23=7.D0*F*D6 +F1*D7 
      F24=9.D0*F*D8 +F1*D9 
      F25=11.D0*F*D10+F1*D11
      F26=13.D0*F*D12+F1*D13
      B( 1)=2.D0*D*T
      B( 2)=2.D0*D*TS 
      B( 3)=2.D0*D
      B( 4)=2.D0*D/T
      B( 5)=2.D0*D/T2 
      B( 6)=3.D0*D2*T 
      B( 7)=3.D0*D2 
      B( 8)=3.D0*D2/T 
      B( 9)=3.D0*D2/T2
      B(10)=4.D0*D3*T 
      B(11)=4.D0*D3 
      B(12)=4.D0*D3/T 
      B(13)=5.D0*D4 
      B(14)=6.D0*D5/T 
      B(15)=6.D0*D5/T2
      B(16)=7.D0*D6/T 
      B(17)=8.D0*D7/T 
      B(18)=8.D0*D7/T2
      B(19)=9.D0*D8/T2
      B(20)=F21/T2
      B(21)=F21/T3
      B(22)=F22/T2
      B(23)=F22/T4
      B(24)=F23/T2
      B(25)=F23/T3
      B(26)=F24/T2
      B(27)=F24/T4
      B(28)=F25/T2
      B(29)=F25/T3
      B(30)=F26/T2
      B(31)=F26/T3
      B(32)=F26/T4
      P=0.D0 
      DO 201 I=1,M
  201 P=P+B(I)*G(I) 
      P=P+R*T 
      PP=P
      RETURN

C     ENTRY DPDT
  300 X( 1)=D2
      X( 2)=D2/(2.D0*TS)
      X( 3)=0.D0
      X( 4)=-D2/T2
      X( 5)=-2.D0*D2/T3 
      X( 6)=D3
      X( 7)=0.D0
      X( 8)=-D3/T2
      X( 9)=-2.D0*D3/T3 
      X(10)=D4
      X(11)=0.D0
      X(12)=-D4/T2
      X(13)=0.D0
      X(14)=-D6/T2
      X(15)=-2.D0*D6/T3 
      X(16)=-D7/T2
      X(17)=-D8/T2
      X(18)=-2.D0*D8/T3 
      X(19)=-2.D0*D9/T3 
      X(20)=-2.D0*D3*F/T3 
      X(21)=-3.D0*D3*F/T4 
      X(22)=-2.D0*D5*F/T3 
      X(23)=-4.D0*D5*F/T5 
      X(24)=-2.D0*D7*F/T3 
      X(25)=-3.D0*D7*F/T4 
      X(26)=-2.D0*D9*F/T3 
      X(27)=-4.D0*D9*F/T5 
      X(28)=-2.D0*D11*F/T3
      X(29)=-3.D0*D11*F/T4
      X(30)=-2.D0*D13*F/T3
      X(31)=-3.D0*D13*F/T4
      X(32)=-4.D0*D13*F/T5
      P=0.D0 
      DO 301 I=1,M 
  301 P=P+G(I)*X(I) 
      PP=P+R*D
      RETURN

C     ENTRY DSDN
C     PARTIAL OF ENTROPY WITH 
C     RESPECT TO THE G COEFFICIENTS 
C     S=S0-R*LOGF(D*R*T/P0)+(DSDN(D)-DSDN(0))*1000.DO +CPOS(T)
  400 G1=F/(2.D0*GM)
      G2=(F*D2-2.D0*G1)/(2.D0*GM) 
      G3=(F*D4-4.D0*G2)/(2.D0*GM) 
      G4=(F*D6-6.D0*G3)/(2.D0*GM) 
      G5=(F*D8-8.D0*G4)/(2.D0*GM) 
      G6=(F*D10-10.D0*G5)/(2.D0*GM) 
      X( 1)=-D
      X( 2)=-D/(2.D0*TS)
      X( 3)=0.D0
      X( 4)=+D/T2 
      X( 5)=2.D0*D/T3 
      X( 6)=-D2/2.D0
      X( 7)=0.D0
      X( 8)=D2/(2.D0*T2)
      X( 9)=D2/T3 
      X(10)=-D3/3.D0
      X(11)=0.D0
      X(12)=D3/(3.D0*T2)
      X(13)=0.D0
      X(14)=D5/(5.D0*T2)
      X(15)= 2.D0*D5/(5.D0*T3)
      X(16)=D6/(6.D0*T2)
      X(17)=D7/(7.D0*T2)
      X(18)=2.D0*D7/(7.D0*T3) 
      X(19)=D8/(4.D0*T3)
      X(20)=2.D0*G1/T3
      X(21)=3.D0*G1/T4
      X(22)=2.D0*G2/T3
      X(23)=4.D0*G2/T5
      X(24)=2.D0*G3/T3
      X(25)=3.D0*G3/T4
      X(26)=2.D0*G4/T3
      X(27)=4.D0*G4/T5
      X(28)=2.D0*G5/T3
      X(29)=3.D0*G5/T4
      X(30)=2.D0*G6/T3
      X(31)=3.D0*G6/T4
      X(32)=4.D0*G6/T5
      P=0.D0
      DO 401 I=1,M 
  401 P=P+G(I)*X(I) 
      PP=P
      RETURN

C     ENTRY DUDN
C     TERMS NEEDED FOR ENTHALPY CALCULATION 
C     H=H0+(T*DSDN(D)-DSDN(0))*1000.D0+(DUDN(D-DUDN(0))*1000.D0+CPOH(T) 
C     +(P/D-R*T)*1000.D0
  500 G1=F/(2.D0*GM)
      G2=(F*D2-2.D0*G1)/(2.D0*GM) 
      G3=(F*D4-4.D0*G2)/(2.D0*GM) 
      G4=(F*D6-6.D0*G3)/(2.D0*GM) 
      G5=(F*D8-8.D0*G4)/(2.D0*GM) 
      G6=(F*D10-10.D0*G5)/(2.D0*GM) 
      X( 1)=D*T 
      X( 2)=D*TS
      X( 3)=D 
      X( 4)=D/T 
      X( 5)=D/T2
      X( 6)=D2*T/2.D0 
      X( 7)=D2/2.D0 
      X( 8)=D2/(2.D0*T) 
      X( 9)=D2/(2.D0*T2)
      X(10)=D3*T/3.D0 
      X(11)=D3/3.D0 
      X(12)=D3/(3.D0*T) 
      X(13)=D4/4.D0 
      X(14)=D5/(5.D0*T) 
      X(15)=D5/(5.D0*T2)
      X(16)=D6/(6.D0*T) 
      X(17)=D7/(7.D0*T) 
      X(18)=D7/(7.D0*T2)
      X(19)=D8/(8.D0*T2)
      X(20)=G1/T2 
      X(21)=G1/T3 
      X(22)=G2/T2 
      X(23)=G2/T4 
      X(24)=G3/T2 
      X(25)=G3/T3 
      X(26)=G4/T2 
      X(27)=G4/T4 
      X(28)=G5/T2 
      X(29)=G5/T3 
      X(30)=G6/T2 
      X(31)=G6/T3 
      X(32)=G6/T4
      P=0.D0
      DO 501 I=1,M
  501 P=P+G(I)*X(I) 
      PP=P
      RETURN

C     ENTRY TDSDT 
C     TEMP. TIMES THE PARTIAL OF
C     ENTROPY WITH RESPECT TO TEMP. 
C     CV=CV0+(TDSDN(/)-TDSDN(D))*1000.D0
  600 G1=F/(2.D0*GM)
      G2=(F*D2-2.D0*G1)/(2.D0*GM) 
      G3=(F*D4-4.D0*G2)/(2.D0*GM) 
      G4=(F*D6-6.D0*G3)/(2.D0*GM) 
      G5=(F*D8-8.D0*G4)/(2.D0*GM) 
      G6=(F*D10-10.D0*G5)/(2.D0*GM) 
      X(1)=0.D0
      X( 2)=-D/(4.D0*TS)
      X(3)=0.D0
      X( 4)=2.D0*D/T2 
      X( 5)=6.D0*D/T3 
      X(6)=0.D0
      X(7)=0.D0
      X( 8)=D2/T2 
      X( 9)=3.D0*D2/T3
      X(10)=0.D0
      X(11)=0.D0
      X(12)=(2.D0*D3)/(3.D0*T2) 
      X(13)=0.D0
      X(14)=(2.D0*D5)/(5.D0*T2) 
      X(15)=(6.D0*D5)/(5.D0*T3) 
      X(16)=D6/(3.D0*T2)
      X(17)=(2.D0*D7)/(7.D0*T2) 
      X(18)=(6.D0*D7)/(7.D0*T3) 
      X(19)=(3.D0*D8)/(4.D0*T3) 
      X(20)=6.D0*G1/T3 
      X(21)=12.D0*G1/T4 
      X(22)=6.D0*G2/T3 
      X(23)=20.D0*G2/T5 
      X(24)=6.D0*G3/T3 
      X(25)=12.D0*G3/T4 
      X(26)=6.D0*G4/T3 
      X(27)=20.D0*G4/T5 
      X(28)=6.D0*G5/T3 
      X(29)=12.D0*G5/T4 
      X(30)=6.D0*G6/T3 
      X(31)=12.D0*G6/T4 
      X(32)=20.D0*G6/T5
      P=0.D0
      DO 601 I=1,M 
  601 P=P+G(I)*X(I) 
      PP=P
      END
c------------------------------------------------------------------
 
      DOUBLE PRECISION FUNCTION FINDTV(POBS)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  ITERATES VAPOR PRESS EQN TO FIND TEMP(K), FOR INPUT OF PRESS(MPA). 
C     GIVEN AN INPUT PRESSURE(MPA)
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      T=VP(8)
      DO 7 I=1,10 
      P=VPN(T)
      IF(DABS (P-POBS)-.000001D0*POBS)8,8,6
    6 CONTINUE
      CORR=(POBS-P)/DPDTVP(T,P)
    7 T=T+CORR
    8 CONTINUE
      FINDTV=T
      END 
      DOUBLE PRECISION FUNCTION CV(D,T)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  CALCULATES CV(J/(MOL*K)).  INPUT DENS(MOL/L) AND TEMP(K).
      DATA R/8.31434D0/ 
      DD=D
      TT=T
      CALL PROPS(CD,DD,TT,6)
      DD=0.0D0
      CALL PROPS(C0,DD,TT,6)
      CV=CPI(TT,1)+(C0-CD)*1000.D0
      CV=CV-R 
      END
c------------------------------------------------------------------
 
      FUNCTION FINDD(P,T)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      TT=T
      IF(TT.GT.VP(8)*.99999D0)GO TO 100 
      IF( P.GT.VPN(TT))GO TO 101
      DD=SATV(TT) 
      GO TO 102 
  100 PC=PCC
      X=(1.1D0/(9.D0*PC))*P+.7D0/9.D0 
      DD=P/(R*T*X)
      IF(P/PC.GT.2D0.AND.T/VP(8).LT.2.5D0)DD=DTP
      GO TO 102 
  101 DD=SATL(TT) 
  102 CONTINUE
      DO 10 I=1,50
      IF(DD.LE.0.0D0.OR.DD.GT.50.D0)GO TO 11
      CALL PROPS(PP,DD,TT,1)
      IF(PP.LE.0.0D0)GO TO 11 
      P2=PP 
      IF(DABS (P-P2)-1.D-7*P)20,20,1 
    1 CALL PROPS(PP,DD,TT,2) 
      DP=PP 
      CORR=(P2-P)/DP
      IF(DABS (CORR)-1.D-7*DD)20,20,10 
   10 DD=DD-CORR
   11 CALL REGULA(P,DD,T) 
   20 FINDD=DD
      END 
c------------------------------------------------------------------

      SUBROUTINE REGULA(PI,DD,TT)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  ITERATES EQN OF STATE FOR DENSITY WHEN SUBPROG FINDD FAILS.
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      T=TT
      P=PI
      D2=0.0D0
      IF(T.LT.TCC)GO TO 10
      D0=DCC*TCC/T
      GO TO 20
   10 PP=VPN(T) 
      IF(P.GT.PP)GO TO 15 
      D0=SATV(T)
      DO 11 I=1,150 
      CALL PROPS(P0,D0,T,1) 
      IF(P0.GE.P)GO TO 12 
   11 D0=D0+.0001D0*D0
      GO TO 42
   12 D1=D0 
   13 CALL PROPS(P1,D1,T,1) 
      IF(P1.LT.P)GO TO 14 
      IF(D1.LE..1D0*PTP)GO TO 42
      D0=D1 
      Z=(P1-P)/P
      IF(Z.LT..1D0)Z=.1D0 
      IF(Z.GT..9D0)Z=.9D0 
      D1=D1-Z*D1
      GO TO 13
   14 CALL PROPS(P0,D0,T,1) 
      DO 140 I=1,50 
      D=D1
      P3=P1 
      IF(DABS(P-P1).LT..00001D0*P)GO TO 40 
      P2=P-P1 
      D1=D1+(D1-D0)*P2/(P1-P0)
      IF(DABS(D-D1).LE..00001D0*D)GO TO 40 
      IF(DABS(P-P1).LT..005D0*P)D2=FINDM(P,T,D1) 
      IF(D2.GT.0.0D0.AND.D2.LT.50.D0)D1=D2
      D2=0.0D0
      CALL PROPS(P1,D1,T,1) 
      IF(P0.GT.P.AND.P1.GT.P)GO TO 120
      IF(P0.LT.P.AND.P1.LT.P)GO TO 120
      GO TO 140 
  120 P0=P3 
      D0=D
  140 CONTINUE
      GO TO 41
   15 D0=SATL(T)
      DO 16 I=1,10
      CALL PROPS(P0,D0,T,1) 
      IF(P0.LE.P)GO TO 17 
   16 D0=D0-.0001D0*D0
      GO TO 42
   17 D1=D0 
   18 CALL PROPS(P1,D1,T,1) 
      IF(D1.GE.50.D0)GO TO 42 
      IF(P1.GT.P)GO TO 14 
      D0=D1 
      Z=(P-P1)/P
      Z=Z*10.D0
      IF(T/TCC.LT..6D0)Z=1.D0 
      IF(Z.LT.1.D0)Z=1.D0 
      IF(Z.GT.9.D0)Z=9.D0 
      D1=D1+.01D0*D1*Z
      GO TO 18
   20 CALL PROPS(P0,D0,T,1) 
      IF(P.LE.P0)GO TO 30 
      D1=D0 
   21 CALL PROPS(P1,D1,T,1) 
      IF(P1.GE.P)GO TO 14 
      IF(D1.GE.50.D0)GO TO 42 
      D0=D1 
      Z=(P-P1)/P
      Z=Z*10.D0
      IF(Z.LT.1.D0)Z=1.D0 
      IF(Z.GT.9.D0)Z=9.D0 
      D1=D1+.1D0*D1*Z 
      GO TO 21
   30 D1=D0 
   31 CALL PROPS(P1,D1,T,1) 
      IF(P1.LE.P)GO TO 14 
      IF(D1.LE..1D0*PTP)GO TO 42
      D0=D1 
      Z=(P1-P)/P
      Z=Z*10.D0
      IF(Z.LT.1.D0)Z=1.D0 
      IF(Z.GT.9.D0)Z=9.D0 
      D1=D1-.1D0*D1*Z 
      GO TO 31
   40 DD=D1 
      RETURN


   41 WRITE(*,101)P,T,D 
  102 FORMAT(' REGULA FAILED AT P=',E11.4,' AND T=',F7.2) 
  101 FORMAT(' DENSITY ITERATION FAILED AT P=',F7.2,' AND T=',F7.2, 
     1/' DENSITY RETURNED IS',E17.8)
      RETURN
   42 WRITE(*,102)P,T 
      END 

c------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION CP(D,T)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  CALCULATES CP(J/(MOL*K)).  INPUT DENS(MOL/L), TEMP(K). 
      CVEE=CV(D,T)
      CALL PROPS(DPT,D,T,3)
      CALL PROPS(DPD,D,T,2)
      CP=CVEE+(T/(D**2)*(DPT**2)/DPD)*1000.D0 
      END 

c-----------------------------------------------------------------

      DOUBLE PRECISION FUNCTION DPDTVP(TT,P)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  CALCULATES THE DERIVATIVE OF PRESSURE WITH RESPECT TO TEMPERATURE
C  AT SATURATION.  INPUT IS TEMP(K),  OUTPUT IS DPDT(MPA/K).
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      T=TT
      IF(TT.GT.VP(8))GO TO 1
      X=(1.D0-VP(7)/T)/(1.D0-VP(7)/VP(8)) 
      DXDT=(VP(7)/T**2)/(1.D0-VP(7)/VP(8))
      DPDT=VP(1)*DXDT+2.D0*VP(2)*X*DXDT+VP(3)*3.D0*X**2*DXDT+VP(5)* 
     1((1.D0-X)**VP(6))*DXDT+VP(5)*X*((1.D0-X)**(VP(6)-1.D0))*VP(6)
     2*(-DXDT)
      DPDTVP=DPDT*P
      RETURN
    1 DPDTVP=0.0D0
      END 
c------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION FINDM(P,T,DD)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  SOLVES FOR DENSITY(MOL/L) BY ITERATION.  INPUT IS PRESSURE(MPA), 
C  TEMPERATURE(K), AND A STARTING VALUE OF DENSITY.  THIS FCN IS AN 
C  ALTERNATIVE FOR  FUNCTION FINDD. 
      TT=T
      D=DD
      DO 10 I=1,50
      CALL PROPS(PP,D,TT,1)
      P2=PP 
      IF(DABS(P-P2)-1.D-7*P)20,20,1 
    1 CALL PROPS(PP,D,TT,2) 
      DP=PP 
      CORR=(P2-P)/DP
      IF(DABS(CORR)-1.D-7*D)20,20,10
   10 D=D-CORR
      FINDM=FINDD(P,T)
      RETURN
   20 FINDM=D
      END 
c------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ENTHAL(P,D,T)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      R= .00831434D0
      DD=D
      TT=T
      CALL PROPS(SD,DD,TT,4) 
      CALL PROPS(UD,DD,TT,5) 
      DD=0.D0
      CALL PROPS(S0,DD,TT,4) 
      CALL PROPS(U0,DD,TT,5) 
      ENTHAL=T*(SD-S0)*1000.D0+(UD-U0)*1000.D0+CPI(T,3)+(P/D-R*T)*1.D+3
      END 
c--------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ENTROP(D,T)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  CALCULATES ENTROPY(J/(MOL-K), FROM INPUT OF DENSITY(MOL/L) AND TEMP(K).
      R=  .00831434D0 
      P0=  .101325D0
      DD=D
      TT=T
      CALL PROPS(SD,DD,TT,4) 
      DD=0
      CALL PROPS(S0,DD,TT,4) 
      ENTROP=(SD-S0)*1000.D0-R*DLOG(D*R*T/P0)*1000.D0+CPI(T,2) 
      END 
c------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION SATL(TT)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  CALCULATES DENSITY(MOL/L) OF SATURATED LIQUID.  INPUT IS TEMP(K).
      DIMENSION A(20) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/SATC/A,DTPV,EG
      COMMON/ISP/N,NW,NWW
      IF(NW.EQ.1)GO TO 30
      T=TT
      K=14
      KK=7
   10 IF(T.GE.TCC*.99999D0)GO TO 20 
      ITT=TCC 
      IF(ITT+1-T.LT.1.D0)T=ITT
      X=(T-TCC)/(TTP-TCC) 
      D=A(K)*DLOG(X)
      DO 11 I=2,KK
      K=K+1 
      MM=I
      IF(MM.GE.5)MM=MM+1
   11 D=D+A(K)*(1.D0-X**((MM-5)/3.D0))
      IF(K.LT.14)GO TO 12 
      D=DCC+DEXP(D)*(DTP-DCC)
      GO TO 13
   12 D=DCC+DEXP(D)*(DTPV-DCC) 
   13 SATL=D
      IF(ITT+1-TT.LT.1)SATL=D-(D-DCC)*(TT-T) 
      RETURN
   20 SATL=DCC
      RETURN
   30 CALL SSATL(DL,TT)
      SATL=DL
      END
c-----------------------------------------------------------------

      DOUBLE PRECISION FUNCTION SATV(TT)
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  CALCULATES DENSITY(MOL/L) OF SATURATED VAPOUR.  INPUT IS TEMP(K).
      DIMENSION A(20) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
c###  COMMON/SATC/A,DTPV1,EG          ! Original using DTPV1
      COMMON/SATC/A,DTPV,EG           ! Changed by Zhou Yang using DTPV      COMMON/ISP/N,NW,NWW
      IF(NW.EQ.1)GO TO 30
      K=1 
      KK=13 
      T=TT
   10 IF(T.GE.TCC*.99999D0)GO TO 20 
      ITT=TCC 
      IF(ITT+1-T.LT.1.D0)T=ITT
      X=(T-TCC)/(TTP-TCC) 
      D=A(K)*DLOG(X)
      DO 11 I=2,KK
      K=K+1 
      MM=I
      IF(MM.GE.5)MM=MM+1
   11 D=D+A(K)*(1.D0-X**((MM-5)/3.D0))
      IF(K.LT.14)GO TO 12 
      D=DCC+DEXP(D)*(DTP-DCC)
      GO TO 13
   12 D=DCC+DEXP(D)*(DTPV-DCC) 
   13 SATV=D
      IF(ITT+1-TT.LT.1) SATV=D-(D-DCC)*(TT-T) 
      RETURN
   20 SATV=DCC
      RETURN
   30 CALL SSATV(DV,TT)
      SATV=DV
      END 
c.......................................................................
c
c
      DOUBLE PRECISION FUNCTION PTVISC(PRESPA,TEMPK)
c
c......................................................................
c FOR LH2 only, (IF=1), delivered by Rocketdyne on 4/4/90
c......................................................................
   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   

      DIMENSION PS(20),TS(20),JP(21),MX(21),LOC(21),BP(21),DP(21),BT(21)
     &,DT(21),V(536)

      DIMENSION AA( 66),AB( 67),AC( 67),AD( 66),AE( 67),AF( 66),AG( 66) 
     & ,AH( 66),AI(  5) 

      EQUIVALENCE( V,AA),( V(  67),AB),( V( 134),AC),( V( 201),AD)  
     &    ,( V( 267),AE),( V( 334),AF),( V( 400),AG),( V( 466),AH)  
     &    ,( V( 532),AI)

      DATA PS/1.022,2.,4.,8.,14.,25.,43.,69.,99.,128.,151.,165.,176.,   
     &182.,185.,186.5,187.25,187.46875,187.506,187.6385/

      DATA TS/24.845,27.07,29.81,33.07,36.18,39.96,44.12,48.33,51.97,54.
     &79,56.72,57.80,58.57,58.99,59.18,59.29,59.34,59.353,59.356,59.4/  

      DATA LOC/1,21,37,46,55,90,112,134,155,176,204,216,276,296,320,356,
     &376,436,480,510,528/  

      DATA JP/4,4,3,3,5,3,3,3,3,4,3,6,4,6,4,5,5,11,6,3,3/   

      DATA MX/2,2,1,1,3,1,1,1,1,2,1,4,2,4,2,3,3,9,4,1,1/

      DATA BP/0.   ,0.,-1., 0.  ,1.,-10., 0. ,1469.6,1469.6,1469.6,1469.
     &6,0.,0.,0.,587.84,0.,293.92,190.,190.,190.,190./  

      DATA DP/1000., 1000.,3.,1500.,1.,20.,1500.,1763.52,1763.52,1175.68
     &,1763.52,293.92,73.48,293.92,293.92,146.96,73.48,10.,20.,50.,50./ 

      DATA BT/180.,500.,2000.,2000.,3000.,3000.,3000.,30.6,41.4,63.,126.
     &,27.,36.,126.,54.,99.,59.4,59.4,64.8,72.,81./ 

      DATA DT/80.,500.,500.,500.,500.,500.,500.,1.8,3.6,9.,18.,3.6,18., 
     &18.,9.,9.,3.6,1.8,1.8,1.8,9./ 

      DATAAA/2.354D-11,2.797D-11,3.278D-11,3.762D-11,3.1D-11,3.452D-11,3
     &.812D-11,4.172D-11,3.754D-11,3.99D-11,4.22D-11,4.451D-11,4.351D-11
     &,4.498D-11,4.638D-11,4.779D-11,4.906D-11,4.986D-11,5.06D-11,5.138D
     &-11,4.906D-11,4.986D-11,5.06D-11,5.138D-11,7.852D-11,7.688D-11,7.5
     &43D-11,7.418D-11,1.033D-10,1.002D-10,9.74D-11,9.495D-11,1.256D-10,
     &1.213D-10,1.174D-10,1.14D-10,1.256D-10,1.256D-10,1.256D-10,1.463D-
     &10,1.463D-10,1.462D-10,1.657D-10,1.657D-10,1.657D-10,1.256D-10,1.1
     &93D-10,1.14D-10,1.463D-10,1.384D-10,1.317D-10,1.657D-10,1.564D-10,
     &1.485D-10,1.657D-10,1.657D-10,1.657D-10,1.657D-10,1.657D-10,1.843D
     &-10,1.843D-10,1.842D-10,1.842D-10,1.842D-10,2.027D-10,2.024D-10/  

      DATAAB/2.023D-10,2.023D-10,2.022D-10,2.223D-10,2.215D-10,2.21D-10,
     &2.208D-10,2.206D-10,2.441D-10,2.427D-10,2.418D-10,2.412D-10,2.407D
     &-10,2.617D-10,2.635D-10,2.634D-10,2.631D-10,2.626D-10,2.647D-10,2.
     &742D-10,2.782D-10,2.801D-10,2.811D-10,1.658D-10,1.657D-10,1.655D-1
     &0,1.844D-10,1.842D-10,1.84D-10,2.026D-10,2.021D-10,2.018D-10,2.221
     &D-10,2.201D-10,2.194D-10,2.446D-10,2.394D-10,2.377D-10,2.681D-10,2
     &.608D-10,2.576D-10,2.782D-10,2.821D-10,2.794D-10,0.,1.657D-10,1.56
     &4D-10,1.485D-10,1.842D-10,1.736D-10,1.645D-10,2.02D-10,1.901D-10,1
     &.799D-10,2.197D-10,2.06D-10,1.948D-10,2.38D-10,2.217D-10,2.094D-10
     &,2.58D-10,2.373D-10,2.238D-10,2.799D-10,2.534D-10,2.385D-10,0./   

      DATAAC/1.838D-10,2.745D-10,0.,1.645D-10,2.572D-10,0.,1.496D-10,2.3
     &41D-10,3.211D-10,1.381D-10,2.145D-10,2.94D-10,1.266D-10,1.958D-10,
     &2.659D-10,1.182D-10,1.804D-10,2.432D-10,1.1D-10,1.674D-10,2.237D-1
     &0,1.1D-10,1.674D-10,2.237D-10,9.723D-11,1.467D-10,1.938D-10,8.717D
     &-11,1.308D-10,1.717D-10,7.898D-11,1.184D-10,1.547D-10,7.217D-11,1.
     &084D-10,1.413D-10,6.587D-11,9.78D-11,1.297D-10,6.091D-11,9.124D-11
     &,1.209D-10,6.091D-11,8.175D-11,1.02D-10,1.209D-10,5.159D-11,7.068D
     &-11,8.807D-11,1.043D-10,4.432D-11,6.179D-11,7.734D-11,9.191D-11,3.
     &911D-11,5.524D-11,6.936D-11,8.264D-11,3.553D-11,5.021D-11,6.317D-1
     &1,7.54D-11,3.28D-11,4.627D-11,5.823D-11,6.956D-11,3.134D-11/  

      DATAAD/4.364D-11,5.458D-11,6.503D-11,3.05D-11,4.65D-11,6.177D-11,2
     &.946D-11,4.211D-11,5.525D-11,2.962D-11,4.035D-11,5.165D-11,3.034D-
     &11,3.954D-11,4.874D-11,1.29D-10,1.48D-10,1.677D-10,1.866D-10,2.058
     &D-10,2.25D-10,1.033D-10,1.183D-10,1.335D-10,1.493D-10,1.656D-10,1.
     &825D-10,8.563D-11,9.818D-11,1.107D-10,1.233D-10,1.363D-10,1.496D-1
     &0,7.246D-11,8.358D-11,9.433D-11,1.05D-10,1.158D-10,1.267D-10,6.196
     &D-11,7.234D-11,8.199D-11,9.135D-11,1.006D-10,1.1D-10,5.305D-11,6.3
     &24D-11,7.223D-11,8.073D-11,8.902D-11,9.723D-11,4.487D-11,5.552D-11
     &,6.423D-11,7.216D-11,7.974D-11,8.717D-11,3.647D-11,4.869D-11,5.747
     &D-11,6.505D-11,7.214D-11,7.898D-11,6.15D-12,4.232D-11,5.16D-11/   

      DATAAE/5.906D-11,6.576D-11,7.217D-11,6.56D-12,3.651D-11,4.645D-11,
     &5.353D-11,6.013D-11,6.587D-11,6.3D-12,7.49D-12,5.8D-12,0.,9.3D-12,
     &1.004D-11,1.209D-11,2.554D-11,1.198D-11,1.249D-11,1.318D-11,1.434D
     &-11,1.442D-11,1.483D-11,1.53D-11,1.587D-11,1.667D-11,1.701D-11,1.7
     &38D-11,1.78D-11,1.877D-11,2.006D-11,2.192D-11,2.457D-11,2.739D-11,
     &3.019D-11,2.075D-11,2.185D-11,2.325D-11,2.513D-11,2.729D-11,2.946D
     &-11,2.263D-11,2.36D-11,2.473D-11,2.617D-11,2.785D-11,2.962D-11,2.4
     &28D-11,2.532D-11,2.625D-11,2.744D-11,2.88D-11,3.029D-11,5.444D-11,
     &6.191D-11,6.881D-11,7.543D-11,4.154D-11,4.921D-11,5.559D-11,6.091D
     &-11,3.164D-11,3.988D-11,4.613D-11,5.159D-11,2.517D-11,3.308D-11/  

      DATAAF/3.919D-11,4.432D-11,2.234D-11,2.867D-11,3.427D-11,3.911D-11
     &,2.131D-11,2.623D-11,3.099D-11,3.553D-11,2.113D-11,2.509D-11,2.9D-
     &11,3.28D-11,2.141D-11,2.463D-11,2.792D-11,3.134D-11,2.192D-11,2.45
     &7D-11,2.739D-11,3.019D-11,1.557D-11,1.635D-11,1.739D-11,1.901D-11,
     &2.131D-11,1.667D-11,1.738D-11,1.826D-11,1.948D-11,2.113D-11,1.774D
     &-11,1.839D-11,1.916D-11,2.014D-11,2.141D-11,1.877D-11,1.941D-11,2.
     &006D-11,2.099D-11,2.192D-11,3.651D-11,3.9D-11,4.148D-11,4.397D-11,
     &4.645D-11,2.868D-11,3.35D-11,3.675D-11,3.938D-11,4.154D-11,2.076D-
     &11,2.766D-11,3.178D-11,3.483D-11,3.716D-11,1.748D-11,2.262D-11,2.7
     &22D-11,3.068D-11,3.329D-11,1.649D-11,2.D-11,2.377D-11,2.712D-11/  

      DATAAG/3.D-11,1.614D-11,1.857D-11,2.156D-11,2.448D-11,2.72D-11,1.6
     &12D-11,1.79D-11,2.025D-11,2.27D-11,2.517D-11,1.625D-11,1.766D-11,1
     &.95D-11,2.154D-11,2.37D-11,1.649D-11,1.766D-11,1.913D-11,2.086D-11
     &,2.268D-11,1.676D-11,1.774D-11,1.896D-11,2.044D-11,2.202D-11,1.707
     &D-11,1.79D-11,1.892D-11,2.019D-11,2.157D-11,1.739D-11,1.82D-11,1.9
     &01D-11,2.016D-11,2.131D-11,2.484D-11,2.818D-11,2.97D-11,3.086D-11,
     &3.181D-11,3.261D-11,3.335D-11,3.401D-11,3.461D-11,3.52D-11,3.573D-
     &11,1.431D-11,1.56D-11,1.812D-11,2.297D-11,2.59D-11,2.757D-11,2.885
     &D-11,2.987D-11,3.076D-11,3.154D-11,3.225D-11,1.368D-11,1.429D-11,1
     &.511D-11,1.626D-11,1.791D-11,2.015D-11,2.259D-11,2.458D-11/   

      DATAAH/2.611D-11,2.728D-11,2.831D-11,1.35D-11,1.394D-11,1.438D-11,
     &1.508D-11,1.578D-11,1.689D-11,1.8D-11,1.949D-11,2.099D-11,2.246D-1
     &1,2.394D-11,1.35D-11,1.438D-11,1.578D-11,1.8D-11,2.099D-11,2.394D-
     &11,1.347D-11,1.412D-11,1.503D-11,1.634D-11,1.815D-11,2.032D-11,1.3
     &53D-11,1.404D-11,1.471D-11,1.562D-11,1.683D-11,1.833D-11,1.363D-11
     &,1.405D-11,1.459D-11,1.526D-11,1.615D-11,1.725D-11,1.377D-11,1.414
     &D-11,1.458D-11,1.512D-11,1.582D-11,1.664D-11,1.377D-11,1.484D-11,1
     &.664D-11,1.393D-11,1.486D-11,1.635D-11,1.409D-11,1.489D-11,1.61D-1
     &1,1.427D-11,1.499D-11,1.604D-11,1.445D-11,1.509D-11,1.6D-11,1.464D
     &-11,1.523D-11,1.604D-11,1.464D-11,1.523D-11,1.604D-11,1.563D-11/  

      DATAAI/1.606D-11,1.658D-11,1.662D-11,1.697D-11,1.736D-11/ 

c.......................................................................
c UNITS: PRESPA in MPA, TEMPK in deg K
c
c convert to british units: PRES in PSI, TEMP in deg R
c.........................

      PRES=PRESPA*1000.0D0/6.8947572D0
      TEMP=TEMPk*9.0D0/5.0D0

      P=PRES
      IF(P.LT.1.0D0) P=1.0D0
      T=TEMP
      IF(T.LT.180.0D0) GO TO 7
      IF(T.GE.500.0D0) GO TO 1
      N=1   
      GO TO 30  
    1 IF(T.GE.3000.0D0) GO TO 4   
      IF(T.GE.2000.0D0) GO TO 2   
      N=2   
      GO TO 30  
    2 IF(P.GT.5.0D0) GO TO 3  
      N=3   
      GO TO 30  
    3 N=4   
      GO TO 30  
    4 IF(T.GE.6000.0D0) T=5999.99999D0
      IF(P.GE.5.0D0) GO TO 5  
      N=5   
      GO TO 30  
    5 IF(P.GE.30.0D0) GO TO 6 
      N=6   
      GO TO 30  
    6 N=7   
      GO TO 30  
    7 IF(P.LT.1469.6D0) GO TO 12  
      IF(T.GE.  63.0D0) GO TO 10  
      IF(T.GE.41.4D0)   GO TO 9   
      N=8   
    8 TZ=24.84D0+0.00317D0*P
      IF(T.LT.TZ) T=TZ  
      GO TO 30  
    9 N=9   
      GO TO 30  
   10 IF(T.GE.117.0D0) GO TO 11   
      N=10  
      GO TO 30  
   11 N=11  
      GO TO 30  
   12 IF(T.GE.59.4D0)  GO TO 17   
      N=12  
      IF(P.GE.187.6385D0) GO TO 8 
      DO 13 I=2,20  
      IF(P-PS(I))15,14,13   
   13 CONTINUE  
   14 TL=TS(I)  
      GO TO 16  
   15 TL=TS(I-1)+(TS(I)-TS(I-1))*(P-PS(I-1))/(PS(I)-PS(I-1))
   16 IF(T.LE.TL) GO TO 8   
      N=13  
      GO TO 30  
   17 IF(T.LT.126.0D0) GO TO 18   
      N=14  
      GO TO 30  
   18 IF(P.LT.587.84D0) GO TO 19  
      N=15  
      GO TO 30  
   19 IF(T.LT.99.0D0) GO TO 20
      N=16  
      GO TO 30  
   20 IF(P.GE.190.0D0)GO TO 21
      N=13  
      GO TO 30  
   21 IF(P.LT.293.92D0) GO TO 22  
      N=17  
      GO TO 30  
   22 IF(T.GE.72.0D0) GO TO 24
      IF(T.GE.64.8D0) GO TO 23
      N=18  
      GO TO 30  
   23 N=19  
      GO TO 30  
   24 IF(T.GE.81.0D0) GO TO 25
      N=20  
      GO TO 30  
   25 N=21  
   30 FP=(P-BP(N))/DP(N)
      IP=FP 
      IF(IP.GT.MX(N)) IP=MX(N)  
      F=FP-IP   
      FP=1.0-F  
      FT=(T-BT(N))/DT(N)
      IT=FT 
      FF=FT-IT  
      FT=1.0-FF 
      I=IT*JP(N)+IP+LOC(N)  
      J=I+JP(N) 
      PTVISC=FP*FT*V(I)+F*FT*V(I+1)+FP*FF*V(J)+F*FF*V(J+1)  
c.........................................................
c Viscosity units are in lbf-hr/ft^2
c
c convert to Pa-s
c.........................................................

      PTVISC=PTVISC*172368.93600    

      RETURN

      END 

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Modified to SI units by LSA on 4/17/90
c----------------------------------------------------------
 
 
C --  Here is the begin of the old MIPROP2.FOR
c---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION SOUND(D,T)
c....................................................................
c calculates speed of sound (m/s): Input is density(mol/l) and T(K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      COMMON/CRIT/ W, EOK, RM, TC, DC, X , PC, SIG

      CALL PROPS(DP,D,T,2) 
      SOUND=((CP(D,T)/CV(D,T))*DP*1000000.D0/W)**.5D0 
      END 


c------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VISC(DD,T)
c..................................................................
c calculates viscosity (micro Pa s): INPUT is density(mol/l) & T(K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG
      COMMON/ISP/N,NW,NWW

      IF(NW.EQ.1)GO TO 10
      D=DD*GMW/1000.D0
      VISC=DILV(T)+FDCV(D,T)+EXCESV(D,T)
      RETURN

   10 VISC=VISCE(DD,T)
      END 

c-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION THERM(DD,T)
c................................................................
c Returns thermal conductivity (W/(M K)). INPUT: Density(mol/L), T(K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/ISP/N,NW,NWW 
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG

      IF(NW.EQ.1)GO TO 10
      D=DD*GMW/1000.D0
      CR=CRITC(D,T)
      THER=DILT(T)+FDCT(D,T)+EXCEST(D,T)+CR 
      TCI=THER-CR 
      THERM=THER
      RETURN

   10 THERM=THERME(DD,T)
      END 

c----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION EXCESV(D,T)
c...............................................................
c     calculates excess viscosity

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8) 
      R2=D**(.5D0)*((D-EV(8))/EV(8))
      R=D**(.1D0) 
      X=EV(1)+EV(2)*R2+EV(3)*R+EV(4)*R2/(T*T)+EV(5)*R/T**(1.5D0)+EV(6)/T
     1+EV(7)*R2/T 
      X1=EV(1)+EV(6)/T
      EXCESV= DEXP(X)- DEXP(X1) 
      END
c----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION EXCEST(D,T)
c...............................................................
c     calculates excess thermal conductivity

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      R2=D**(.5D0)*((D-ET(8))/ET(8))
      R=D**(.1D0) 
      X=ET(1)+ET(2)*R2+ET(3)*R+ET(4)*R2/(T*T)+ET(5)*R/T**(1.5D0)+ET(6)/T
     1+ET(7)*R2/T 
      X1=ET(1)+ET(6)/T
      EXCEST= DEXP(X)- DEXP(X1) 
      END 
c----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FDCV(D,T)
c...............................................................
c  First density correction for viscosity

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8) 
      FDCV=(FV(1)+FV(2)*(FV(3)-DLOG(T/FV(4)))**2)*D 
      END

c--------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FDCT(D,T)
c.............................................................
c First density correction for thermal cond.

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG

      FDCT=(FT(1)+FT(2)*(FT(3)-DLOG(T/FT(4)))**2)*D 

      END 

c------------------------------------------------------------
      FUNCTION CRITC(D,T)
c............................................................
C  CALCULATES CRITICAL ENHANCEMENT FOR THERM. COND.
C  INPUT UNITS ARE G/CC, K,  OUTPUT IS W/(M*K).

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG
      COMMON/CHECK/ DELD,DELT,DSTAR,TSTAR
      COMMON/ISP/N,NW,NWW 

      AV=6.0225D+23
      BK=1.38054D-16
      DELD=DABS (D-DC)/DC
      DELT=DABS (T-TC)/TC 
C  CALCULATE DISTANCE PARAMETER 
      R=(RM**2.5D0)*(D**0.5D0)*(AV/GMW)**0.5D0
      R=R*(EOK**0.5D0)*X/(T**0.5D0) 
      RRR=R 
C   GENERAL  EQUATION 
      DX=D*1000.D0/GMW 
C  DX IN MOL/L,  D IN G/CM3.
      CALL PROPS(DPT,DX,T,3) 
C  DPDT IN MPA/K. 
      DPT=DPT*1.0D+7
C  DPDT NOW IN DYNES/(CM2*K)
      CALL PROPS(DPD,DX,T,2) 
C  DPDD IN L*MPA/MOL. 
      DPD=DPD*1.0D+7*1000.D0/GMW
C  DPDD NOW IN DYNE*CM/G. 
      IF( DPD.LT.0.0D0) DPD=1.D0 
   94 VIS=VISC(DX,T)*(1.0D-05)
C  VISCOSITY NOW IS G/(CM*S). 
      IF(DELD.GT.0.25D0)GO TO 10
    8 IF(DELT.GT.0.025D0)GO TO 10
    9 COMPRES=SENG(D,T) 
      GO TO 12
   10 COMPRES=1.D0/(D*DPD)**0.5D0
   12 EX=BK*T**2*(DPT**2)*COMPRES 
      EXB=R*((BK*T)**0.5D0)*(D**0.5D0)*((AV/GMW)**0.5D0)
      CRIT=EX/(EXB*6.D0*3.14159D0*VIS) 
C   THERMAL COND, CRIT, IS IN ERG/(CM*SEC*K)
C   PUT IN DAMPING FACTOR 
      BDD=((D-DC)/DC)**4
      BTT=((T-TC)/TC)**2
      BXX= -18.66D0*BTT - 4.25D0*BDD
      IF(BXX.LT.-1.D+2) BXX= -1.D+2 
      FACT= DEXP( BXX )
C     FACT=DEXP (-18.66D0*BTT - 4.25D0*BDD)
      DELC=CRIT*FACT
      CRITC=DELC/100000.D0
C   THERMAL COND, CRITC, IS NOW IN W/(M*K)
      AKT=COMPRES*COMPRES 
      EPSI=R*R*BK*T*(AV*D/GMW)*AKT
      EPSI=EPSI**0.5D0
C   CALC  CP-CV 
      CPCV=T*(DPT**2)*AKT/D 
      END 

c------------------------------------------------------
      DOUBLE PRECISION FUNCTION SENG(D,T)
c......................................................
c SCALED EQUATION OF STATE for CRITICAL REGION

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/CRIT/ GMW, EOK, RM, TC, DC, X , PC, SIG
      COMMON/SEN/BETA,XO,DELTA,E1, E2, AGAM 
      COMMON/CHECK/ DELD,DELT,DSTAR,TSTAR

      DSTAR= D/DC
      TSTAR=T/TC
      BETO=1.D0/BETA
      XX=DELT/DELD**BETO
      AG=AGAM-1.D0 
      BET2= 2.0D0*BETA 
      AGB=AG/BET2 
      DEL1=DELTA-1.D0
      AGBB=(AG-BET2)/BET2
      XXO=(XX+ XO)/XO 
      XXB=XXO**BET2 
      BRAK=1.D0 + E2*XXB 
      BRAK1=BRAK**AGB
      H=E1*XXO*BRAK1
      HPRIM=(E1/XO)*BRAK1 + (AG/XO)*E1*E2*(XXB)*(BRAK**AGBB)
      RCOM=(DELD**DEL1)*(DELTA*H - (XX/BETA)*HPRIM  ) 
      RCOMP=1.D0/(RCOM*DSTAR**2) 
      RCM=RCOMP/(PC*1.0D+7) 
C  RCM IN CM2/DYNE, PC IN MPA 
      RCM=RCM**0.5D0
       SENG=RCM 
      END 

c---------------------------------------------------------
      DOUBLE PRECISION FUNCTION DILV(T)
c.........................................................
c DILUTE GAS VISCOSITY

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8) 

      SUM=0.0D0 
      TF=T**(1.D0/3.D0) 
      TFF=T**(-4.D0/3.D0) 
      DO 10 I=1,9 
      TFF=TFF*TF
   10 SUM=SUM+GV(I)*TFF 
      DILV=SUM
      END

c---------------------------------------------------
      DOUBLE PRECISION FUNCTION DILT(T)
c...................................................
c DILUTE GAS THERMAL CONDUCTIVITY

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)

      TF=T**(1.D0/3.D0) 
      TFF=T**(-4.D0/3.D0) 
      SUM=0 
      DO 20 I=1,9 
      TFF=TFF*TF
   20 SUM=SUM+GT(I)*TFF 
      DILT=SUM
      END 

c---------------------------------------------------
      DOUBLE PRECISION FUNCTION FINDP(D,T)
c..................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      DD=D
      TT=T
      IF(TT.LT.TCC)GO TO 10 
    1 CALL PROPS(PP,DD,TT,1)
      FINDP=PP
      RETURN

   10 P=VPN(TT) 
      DV=FINDD(P-.0001D0,TT)
      DL=FINDD(P+.0001D0,TT)
      IF(DD.LE.DV.OR.DD.GE.DL)GO TO 1 
      WRITE(*,100)DV,DL,DD
      CALL PROPS(PP,DV,TT,1)
      FINDP=PP
      D=DV
  100 FORMAT(' THE STATE POINT YOU HAVE SPECIFIED CORRESPONDS TO A '
     1/' DENSITY IN THE LIQUID VAPOR COEXISTENCE REGION'
     2/' THE DENSITY OF THE SATURATED VAPOR IS ',F6.4,' MOLES/LITER'
     3/' THE DENSITY OF THE SATURATED LIQUID IS ',F8.4,' MOLES/LITER' 
     4/' AND THE INPUT DENSITY IS ',F8.4,' MOLES/LITER' 
     5/' SATURATED VAPOR IS ASSUMED') 
      END 

c------------------------------------------------------
      DOUBLE PRECISION FUNCTION FINDT(P,D)
c......................................................
C  RETURNS TEMPERATURE(K), FROM THE 32-TERM MBWR EQN OF STATE.
C  INPUT IS PRESSURE(MPA) AND DENSITY(MOL/L).

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC

      PP=P
      DD=D
      IF(P.GE.PCC)GO TO 1 
      TSAT=FINDTV(PP) 
      DV=FINDD(PP-.00001D0,TSAT)
      DL=FINDD(PP+.0001D0,TSAT) 
      IF(DD.GT.DV.AND.DD.LT.DL)GO TO 30 
      TT=TSAT 
      GO TO 2 

    1 TT=TCC

    2 DO 10 I=1,10
      CALL PROPS(P2,DD,TT,1)
      IF(DABS(PP-P2)-1.D-7*PP)20,20,11 
   11 CALL PROPS(DP,DD,TT,3) 
      CORR=(P2-PP)/DP 
      IF(DABS(CORR)-1.D-5)20,20,10 
   10 TT=TT-CORR
   20 FINDT=TT
      RETURN
   30 FINDT=TSAT
      D=DV
      WRITE(*,100)DV,DL,DD
  100 FORMAT(' THE STATE POINT YOU HAVE SPECIFIED CORRESPONDS TO' 
     1/' A DENSITY IN THE LIQUID VAPOR COEXISTENCE REGION'
     2/' DENSITY OF THE SATURATED VAPOR IS',F8.4,' MOLES/LITER' 
     3/' DENSITY OF THE SATURATED LIQUID IS',F8.4,' MOLES/LITER'
     4/' INPUT DENSITY IS',F8.4,' MOLES/LITER'
     5/' SATURATED VAPOR CONDITIONS ARE ASSUMED') 
      END 

c---------------------------------------------------------
      DOUBLE PRECISION FUNCTION FDIEL(P,D,T)
c........................................................
c DIELECTRIC Constant. INPUT P(MPA), D(MOL/L), T(K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/DIEL/BX(6),PX(6)
      COMMON/ISP/N,NW,NWW

      IF(NW.EQ.1)GO TO 1 
      CM= BX(1)+ BX(2)*D+ BX(3)*D**2+ BX(4)*D**3+ BX(5)*P+ BX(6)*T
      FDIEL=(1.D0+2.D0*D*CM)/(1.D0-D*CM)
      RETURN

    1 FDIEL=SDIEL(P,D,T) 
      END 

c---------------------------------------------------------
      DOUBLE PRECISION FUNCTION CPI(T,K)
c.........................................................

C  CALCULATES SPECIFIC HEAT, ENTROPY, AND ENTHALPY FOR THE IDEAL GAS. 
C  OUTPUT IS IN J/(MOL*K), FOR CP AND S,  AND J/MOL FOR H.
C  HYDROGEN, (N=1), IS TREATED AS A SPECIAL CASE AS THE COEFF. FOR
C  CP ARE IN THREE TEMPERATURE RANGES.  T < 40 K,  40 < T < 140 K,
C  AND T > 140 K. 

CVAX  CHANGES FOR VAX SYSTEM ARE NOTED WITH CVAX comment

      IMPLICIT REAL*8(A-D)
      IMPLICIT REAL*8(F-H)
      IMPLICIT REAL*8(O-T)
      IMPLICIT REAL*8(V-Z)
      IMPLICIT INTEGER*4(I-N)

CVAX  REAL*16 U, EU
      REAL*8  U,EU    

      COMMON/CPID/G(11),GH(11),GL(11)
      COMMON/ISP/N,NW,NWW
      COMMON/H2/FHI,FSI

      IF(N.NE.1)GO TO 15
      TX1=140.D0
      TX2=40.D0
      DO 110 J=1,11
  110 G(J)=GH(J)
      IF(T.LT.140.D0)GO TO 130
      GO TO 180
  130 CALL SHI(TX1)
      G(10)=G(10)+FHI
      G(11)=G(11)+FSI
      DO 140 J=1,8
  140 G(J)=GL(J)
      CALL SHI(TX1)
      G(10)=G(10)-FHI
      G(11)=G(11)-FSI
      IF(T.LT.40.D0)GO TO 160
      GO TO 180
  160 CALL SHI(TX2)
      G(10)=G(10)+FHI
      G(11)=G(11)+FSI
      DO 170 J=1,8
  170 G(J)=0.0D0
      G(4)=2.5000315D0
      CALL SHI(TX2)
      G(10)=G(10)-FHI
      G(11)=G(11)-FSI
  180 CONTINUE
   15 U=G(9)/T
      EU=DEXP(U)            ! on VAX change to => 
CVAX  EU=QEXP(U)      
      TS=1.D0/T/T/T/T
      GO TO (20,40,55),K
   20 CPI=G(8)*U*U*EU/(EU-1.D0)/(EU-1.D0)
      DO 25 I=1,7 
      TS=TS*T 
   25 CPI=CPI+G(I)*TS 
      CPI=CPI*8.31434D0 
      RETURN

   40 CPI=G(8)*(U/(EU-1.D0)-DLOG(1.D0-1.D0/EU)) 
     1-G(1)*TS*T/3.D0-G(2)*TS*T*T/2.D0-G(3)/T+G(4)*DLOG(T)+G(5)*T
     2+G(6)*T*T/2.D0+G(7)*T**3/3.D0
CVAX  on VAX systems change  statement 40 to:
CVAX40 CPI=G(8)*(U/(EU-1.D0)-QLOG(1.D0-1.D0/EU)) 
CVAX 1-G(1)*TS*T/3.D0-G(2)*TS*T*T/2.D0-G(3)/T+G(4)*DLOG(T)+G(5)*T
CVAX 2+G(6)*T*T/2.D0+G(7)*T**3/3.D0
      CPI=CPI*8.31434D0+G(11) 

      RETURN
   55 CPI=G(8)*U*T/(EU-1.D0)-G(1)/(2.D0*T*T)-G(2)/T+G(3)*DLOG(T)+G(4)*T 
     1+G(5)*T*T/2.D0+G(6)*T**3/3.D0+G(7)*T**4/4.D0
      CPI=CPI*8.31434D0+G(10) 
      END 

c--------------------------------------------------------
      SUBROUTINE SHI(T) 
c.......................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/CPID/G(11),GH(11),GL(11)
      COMMON/H2/FHI,FSI

    1 U=G(9)/T
      EU=DEXP(U) 
      GHI=G(8)*U*T/(EU-1.D0)-G(1)/(2.D0*T*T)-G(2)/T+G(3)*DLOG(T)+G(4)*T 
     A +G(5)*T*T/2.D0+G(6)*T**3/3.D0+G(7)*T**4/4.D0
      FHI=GHI*8.31434D0
      U=G(9)/T
      EU=DEXP(U) 
      TS=1.D0/T**4
      GHS= G(8)*(U/(EU-1.D0)-DLOG(1.D0-1.D0/EU))- 
     A G(1)*TS*T/3.D0-G(2)*TS*T*T/2.D0-G(3)/T+G(4)*DLOG(T)+G(5)*T+
     B G(6)*T*T/2.D0+G(7)*T**3/3.D0 
      FSI=GHS*8.31434D0 

      END 

c------------------------------------------------------
      DOUBLE PRECISION FUNCTION PMELT(T)
c......................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  COMPUTES MELTING PRESSURE(MPA) FOR INPUT TEMPERATURE(K). 
      COMMON/DIEL/BX(6),PX(6) 
      COMMON/ISP/N,NW,NWW 

      IF(N.EQ.1)GO TO 20
   10 PMELT= PX(1)+ PX(2)*T**PX(3)
      RETURN

   20 IF(T.LT.22.D0)GO TO 10
   30 PMELT= PX(4)+ PX(5)*T**PX(6)
      END 

c------------------------------------------------------
      FUNCTION TMELT(P)
c......................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C  COMPUTES MELTING TEMPERATURE(K) FOR INPUT PRESSURE(MPA)
      COMMON/DIEL/BX(6),PX(6) 
      COMMON/ISP/N,NW,NWW 

      IF(N.EQ.1)GO TO 20
   10 TMELT=((P-PX(1))/PX(2))**(1.D0/PX(3)) 
      RETURN

   20 IF(P.LT.31.64D0)GO TO 10
   30 TMELT=((P-PX(4))/PX(5))**(1.D0/PX(6)) 
      END 

c-----------------------------------------------------
      DOUBLE PRECISION FUNCTION VPN(TT)
c.....................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C  CALCULATES VAPOR PRESSURE(MPA), INPUT IS TEMP(K).
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      T=TT
      X=(1.D0-VP(7)/T)/(1.D0-VP(7)/VP(8)) 
      VPN=VP(9)*DEXP (VP(1)*X+VP(2)*X*X+VP(3)*X**3+VP(4)*X**4+VP(5)*X* 
     1(1.D0-X)**VP(6))
      END
c-----------------------------------------------------
      SUBROUTINE SSATL(D1,T1)
c.....................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
C            SATURATED LIQUID AND VAPOR DENSITIES.
C            LIQUID DENSITIES ARE FROM R.D. MCCARTY.
C            VAPOR DENSITIES ARE FROM R.D. GOODWIN. 
      DIMENSION G(32),VP(9),A(20) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/SATC/ A,DTPV,EG
      IF(T1 .GE. TCC) GO TO 20
      X1=(T1-TCC)/(TTP-TCC) 
      Y1=A(7)*DLOG(X1)+ A(8)*(1.D0-1.D0/X1)+
     A A(9)*(1.D0-X1**(-2.D0/3.D0))+ A(10)*(1.D0-X1**(-1.D0/3.D0))+ 
     B A(11)*(1.D0-X1**( 1.D0/3.D0))+ A(12)*(1.D0-X1**( 2.D0/3.D0))+
     C A(13)*(1.D0-X1)
      D1=DCC+ (DTP-DCC)*DEXP(Y1)
      RETURN
   20 D1=DCC
      END

c--------------------------------------------------------
      SUBROUTINE SSATV(D1,T1)
c........................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9),A(20) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/SATC/ A,DTPV,EG
      TT=T1 
      IF(TT .GE. TCC) GO TO 20
      YN=DLOG(DCC/DTPV) 
      X1=(TCC-TT)/(TCC-TTP) 
      Z1=TCC*X1/TT
      Y1= A(1)*Z1+ A(2)*X1**EG+ A(3)*X1+ A(4)*X1**(4.D0/3.D0) 
     A   +A(5)*X1**(5.D0/3.D0)+ A(6)*X1*X1
      D1=DCC*DEXP(-YN*Y1)
      RETURN
   20 D1=DCC
      END 

c--------------------------------------------------------
      DOUBLE PRECISION FUNCTION SDIEL(P,D,T)
c........................................................
c     Dielectric constant. INPUT P(MPA), D(MOL/L) & T(K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION G(32),VP(9) 
      COMMON/DIEL/BX(6),PX(6) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC

      CM= BX(1)+ BX(2)*D+ BX(3)*D**2+ BX(4)*DLOG(1.D0+TCC/T)+ 
     A BX(5)*P
      SDIEL=(1.D0+2.D0*D*CM)/(1.D0-D*CM) 
      END 

c---------------------------------------------------------
      DOUBLE PRECISION FUNCTION THERME(DD,TT)
c.........................................................
c     THERMAL Conductivity (W/M*K)

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET
         TI  = 1.D0/TT 
         TRM0= ET(1)+ ET(2)*TI+ ET(3)*TI*TI 
         TRM1= ET(4)+ ET(5)*TI+ ET(6)*TI*TI 
         TRM2= ET(7)+ ET(8)*TI
      BACKG  = TCOND0(TT)+(TRM0+TRM1*DD)*DD/(1.D0+TRM2*DD)
      THERME = BACKG + TCRIT(DD,TT) 
      END 
c--------------------------------------------------------
      DOUBLE PRECISION FUNCTION TCOND0(TT)
c........................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
C            THERMAL CONDUCTIVITY (W/(M*K)), TT(K). 
C            LOW DENSITY LIMIT. 
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET 
      COMMON/B/ GMW,EOK,SIG
C            R IS GAS CONSTANT IN J/(MOL*K).
      R= 8.31434D0
      CON1= 15.D0*R/4.D0
      CON2= 2.D0*CON1/3.D0
      CP0= CPI(TT,1)
C            CP0 IS SPECIFIC HEAT IN J/(MOL*K), 
C            VSCTY0 IS VISCOSITY IN MICRO-PA*S. 
      ETA0   = VSCTY0(TT,0.0D0)/1000000.D0
      TS     = TT/EOK 
      YC     = (GT(1)+ GT(2)/TS)*(CP0- CON2)
      TCOND0 = 1000.D0*ETA0*(CON1+ YC)/GMW
C            FACTOR OF 1000 CONVERTS G/MOL TO KG/MOL. 
      END 

c------------------------------------------------------
      DOUBLE PRECISION FUNCTION TCRIT(DD,TT) 
c.........................................................
C            CRITICAL ENHANCEMENT W/(M*K), INPUT MOL/L AND K.

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9) 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC
      COMMON/A/ E1,G1,B1,DE,BK,D1,XZ,ZZ,X1,X2,X3,X4 
      COMMON/B/ GMW,EOK,SIG
      PC1=PCC*1.D+6 
      DELD=DABS(DD-DCC)/DCC
      DELT=DABS(TT-TCC)/TCC
      FACT= X1*DELT**4.D0 + X2*DELD**4.D0 
      IF(FACT  .GT.  100.D0) FACT=100.D0
      DFACT=DEXP(-FACT)
      RSTAR=DD/DCC
C            CONVERTING MICRO-PA*S TO PA*S. 
      VIS= 1.0D-06*VISCE(DD,TT) 
      CALL PROPS(DPD,DD,TT,2)
C            CONVERTING M-PA TO PA. 
      DPD=DPD*1.D+6 
      CALL PROPS(DPT,DD,TT,3)
      DPT=DPT*1.D+6 
      IF(DELD .NE. 0.0D0)GO TO 20
C            CRITICAL ISOCHORE. 
   10 BGAM=XZ**G1/D1*((1.D0+E1)/E1)**((G1-1.D0)/(2.D0*B1))
      CHISTAR=BGAM*(DELT)**(-G1)
      GO TO 50
   20 IF(DELD .LE. 0.25D0.AND.DELT .LT. 0.03D0) GO TO 30
      GO TO 40
C            CRITICAL REGION
   30 XX=DELT/DELD**(1.D0/B1)
      Y=(XX+XZ)/XZ
      TOP=DELD**(-G1/B1)*((1.D0+E1)/(1.D0+E1*Y**(2.D0*B1)))**((G1-1.D0
     1)/(2.D0*B1))
      DIV=D1*(DE+(Y-1.D0)*(DE-1.D0/B1+E1*Y**(2.D0*B1))/(1.D0+E1*Y**(2.
     1D0*B1)))
      CHISTAR=TOP/DIV 
      GO TO 50
C            NON CRITICAL REGION
   40 CHISTAR=PC1*DD/(DCC**2*DPD) 
   50 CHI=CHISTAR**X3 
      UPPER=X4*BK/PC1*(TT*DPT/RSTAR)**2*CHI*DFACT 
      SSENG=UPPER/(ZZ*6.D0*3.14159D0*VIS)
      TCRIT=SSENG
      END

c------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VISCE(DS,TS) 
c.................................................................
      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N) 
      DIMENSION G(32),VP(9)
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET 
      COMMON/DATA/G,R,GAMMA,VP,DTP,PCC,PTP,TCC,TTP,TUL,TLL,PUL,DCC

      ETA0= VSCTY0(TS,DS) 
      TRM1= EV(3) + EV(4)*TS**(-3.D0/2.D0)
      TRM2= EV(5) + EV(6)/TS + EV(7)*TS**(-2.D0)
      TRMX= DEXP(EV(1) + EV(2)/TS) 
      R1  = DS**0.1D0 
      R2  = ((DS-DCC)/DCC)*DS**0.5D0
      VISCE = TRMX*(DEXP(TRM1*R1 + TRM2*R2) - 1.D0) + ETA0 
      END 


c----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VSCTY0(TX,DX)
c................................................................
C            VISCOSITY THROUGH LINEAR TERM IN DENSITY IN MICRO-PA*S.
C            DENSITY IS IN MOL/L, AND TEMP IN K.

      IMPLICIT REAL*8(A-H)
      IMPLICIT REAL*8(O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION C0(9)
      DIMENSION GV(9),GT(9),FV(4),FT(4),EV(8),ET(8)
      COMMON/DATA1/GV,GT,FV,FT,EV,ET 
      COMMON/B/ GMW,EOK,SIG 
      DATA C0/-3.0328138281D+00, 1.6918880086D+01,-3.7189364917D+01,
     *         4.1288861858D+01,-2.4615921140D+01, 8.9488430959D+00,
     *        -1.8739245042D+00, 2.0966101390D-01,-9.6570437074D-03/
      TS = TX / EOK 
      TY = 1.D0 / TS 
      TZ= TS**(1.D0/3.D0) 
      E0=0.D0
         DO 200 J=1,9 
         E0=E0+C0(J)*TY 
  200    TY=TY*TZ 
      OM22 = 1.D0 / E0 
C            ETA0 IS VISCOSITY AT THE LOW DENSITY LIMIT.
      ETA0 = 2.6693D0 *DSQRT(GMW * TX) / (SIG*SIG * OM22) 
      ETA1= FV(1) + FV(2)*(FV(3)-DLOG(TX/FV(4)))**2 
      VSCTY0=(ETA0+ETA1*DX) 
      END 

c-----------------------------------------------------------------
c end of mipropst.f  LAST Revised 12/31/93
c-----------------------------------------------------------------