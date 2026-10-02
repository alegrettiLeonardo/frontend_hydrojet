C Updated for angled injection 7/3/95
C
C    #    #    #  #####   #    #   #####  #####           ######
C    #    ##   #  #    #  #    #     #    #    #          #
C    #    # #  #  #    #  #    #     #    #    #          #####
C    #    #  # #  #####   #    #     #    #####    ###    #
C    #    #   ##  #       #    #     #    #        ###    #
C    #    #    #  #        ####      #    #        ###    #
C
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

C SUBS for inputdata from keyboard
C
C *****************************************************************************
C **                                                                         **
C **  Subroutine INSPL                                                       **
C **                                                                         **
C **  INSPL: KBD INPUT of CLEARANCE GEOMETRY & AXIAL LOCATIONS               **
C **                                                                         **
C *****************************************************************************
c

      SUBROUTINE INSPL

      IMPLICIT NONE

c ----------------------------------------------------------------------------
c

      INCLUDE 'params.f'

c.................................................................
c USES spline.f => program CAN BE COPIED from COMMON directory
c.................................................................

      COMMON /PARAM1/ CLEAR, DIAM, LENGTH, LD, AR, HREC
      COMMON /HBLEN/ LENGTHL, LENGTHR
      COMMON/ PARAM4/ Cinlet,Cexit
      COMMON /HJBSTEP/ ClearO,ClearR,ClearL,YR,YL
      COMMON/spldata/Z(Nsl),CL(Nsl),BCL(Nsl),CCL(Nsl),DCL(Nsl),NJ
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION clear,diam,length,ld,ar,hrec,cinlet,cexit,
     +                 dumy, dumyz, dumycl, z,cl, bcl,ccl,dcl,
     +                 lhalf,zval, f,ClearO,ClearR,ClearL,YR,YL,
     +                 lengthr, lengthl, zero

      INTEGER ISYM, ICSTEP,Icase, nj,kcl,j,ICLEAR, BEARING
      CHARACTER*1 dumyc

c................................................................
      zero=0.0D0
      if (nj.lt.1) nj=1
c................................................................
c Nj: Number of points for spline evaluation
c ----------
      lhalf=length/2.0D0

c     ........................!
      IF (ICSTEP.eq.1) THEN
c     ........................! STEP CLEARANCE BEARING
      IF (BEARING.EQ.2) THEN
      write (6, *) ' '
      write (6, *) 'STEP CLEARANCE SEAL:'
      write (6, *) ' .-------. ClearO'
      write (6, *) '         .------- ClearR'
      write (6, *) '---------------------'
      write (6, *) '0       YR          L'
      ELSE
      write (6, *) ' '
      write (6, *) 'STEP CLEARANCE BEARING:'
      write (6, *) '             .-------. ClearO'
      write (6, *) ' ClearL -----.       .------- ClearR'
      write (6, *) '        ---------------------'
      write (6, *) '    -Ll     YL   0   YR      Lr'
      END IF
      write (6, 243) length, lengthl, lengthr
      write (6, 253) ClearO, YR, ClearR, YL, ClearL

c     .......................!
      ELSE
c     .......................! Continuos clearance bearing

      write (6, *) ' '
      write (6, *) 'NJ (INT*4): Number of axial coordinates for'
      write (6, *) '           continuous clearance calculation'
      write (6, *) '------------- DEFAULT VALUES --------------'
      write (6, 243) length, lengthl, lengthr
      write (6,21) Nj
      do j=1,nj                        ! Print default  values
         write (6,23) j, z(j), cl(j)   !
      end do                           !

      END IF
c     .................................!

      write (6, 240)
      call beeper
c............................................................
      dumyc='Y'
      write (6, 25)
      CALL BEEPER
      read (5,26) dumyc
      IF ((dumyc.eq.'N').or.(dumyc.eq.'n')) THEN
          GOTO 1
      ELSE
          GOTO 28
      END IF
c...........................................................

 1    CONTINUE
      ICLEAR=1
      write (6,*) 'SELECT (1) for CONTINUOUS CLEARANCE'
      write (6,*) '       (2) for STEP CLEARANCE'

      CALL ENTERINT(ICLEAR)

      IF (ICLEAR.eq.2) GOTO 31    ! => STEP CLEARANCE

c ::::::::::::::::::::: CONTINUOUS CLEARANCE BEARING ::::::::::::::::::::

 111  write (6, 2) NJ
      ICSTEP=0

  2   format ('$', /,' ','ENTER NJ [Default Nj= ',I2, ']: ',
     + /, ' ENTER : 1 for uniform clearance,' ,/,
     +    '         2 for tapered clearance (linear),',
     + /, '        >2 for nonuniform clearance (spline fit)',
     + /, '   NOTE: Tapered only for symmetric HJB'
     + /, ' .................................................')

      CALL ENTERINT(nj)

c     ............................... Check for Upper limits
           if (nj.gt.nsl) then
              write (6, 75) nj, nsl
              call beeper
              nj=nsl
              goto 1
           end if
      lhalf=lengthr*1.01
c.................................................................
c ENTER Axial coordinates & Clearance values
c.................................................................

      WRITE(6,10) NJ
   10 FORMAT(/,2X,'Number of axial coordinates selected ',I4)

      j=0
      Icase=1
 11   j=j+1
      if(j.gt.Nj) goto 20
      goto 12


 20   write (6,21) Nj
      do j=1,nj                        ! Print typed values
         write (6,23) j, z(j), cl(j)   !
      end do                           !
      write (6,240)
c................................................................
c CHECK IF Input DATA satisfies USER
c......................................
      dumyc='Y'
      write (6, 25)
  25  format (' ',/, 1X,  40('.'),
     + /,' IS THE DATA for AXIAL CLEARANCE OK,',
     + /,' IF satisfied then press (Y or ENTER)',
     +/,2X,'Otherwise press (N or n)',/, 1X,40('.') )
      CALL BEEPER
      read (5,26) dumyc
  26  format (A1)


      if ((dumyc.eq.'N').or.(dumyc.eq.'n')) THEN
c    !............................................!

           write(6,*) 'INPUT Number of data point to be changed'
           call enterint(j)
           IF(j.lt.1) j=1
           IF(j.gt.Nj) j=Nj
           Icase=2
           goto 12
      else
           ICSTEP=0
           goto 28
      end if
c...............................................................
 12   WRITE (6, 241) j

         read(5, *) dumyz,dumycl

         if (j.eq.1) then
            if (Nj.eq.2) THEN
              dumy=-1.0D-7      ! TAPERED SYMMETRIC  !###
            else
              dumy=-Lengthl*1.01
            end if
         else
            dumy=Z(j-1)
         end if

         if(dumyz.lt.dumy) then
             write(6,*) 'WARNING: -------------------------------------'
             write(6,*) 'Axial Y values need to be strictly increasing,'
         if (Nj.eq.2) THEN
             write(6,*) 'Tapered SYMMETRIC: range between [0, Lr]'
         else
             write(6,*) 'in range [-Ll,Lr]'
         end if
             write(6,*) 'TYPE again correct value'
             write(6,*) '----------------------------------------------'
             call beeper
             goto 12
         else
             Z(j)=dumyz
         end if

         if (Z(j).gt.lhalf) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Axial Y value larger than Lr=',Lengthr
             write(6, *) 'TYPE a value less than Lr    '
             write(6, *) '--------------------------------------'
             call beeper
             goto 12
         end if

         if(dumycl.le.(0.0)) then
             write(6,*) 'WARNING: ------------------------------'
             write(6,*) 'Clearance value can not be negative'
             write(6,*) 'TYPE again correct value'
             write(6,*) '---------------------------------------'
             call beeper
             goto 12
         else
           Cl(j)=DABS(Dumycl)
         end if


         goto (11,20) , Icase

C
C...........................................
c calculate spline coefficients for continuous clearance

 28   CONTINUE

      CALL SPLINE(NJ,Z,CL,BCL,CCL,DCL)

      RETURN


c ICLEAR=2 :::::::::::::::::: STEP CLEARANCE BEARING :::::::::::::


 31      write (6,250)         ! Read ClearO
         CALL ENTERVAL(ClearO)
         ClearO=DABS(ClearO)

         if (ClearO.eq.zero) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Negative CLEARANCE at Bearing Center  '
             write(6, *) '--------------------------------------'
             call beeper
             goto 31
         end if
         lhalf=LengthR

 32      write (6,251)         ! read YR, ClearR
         read (5,*) YR, ClearR
         ClearR=DABS(ClearR)

         if (YR.gt.lhalf) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Axial Y value larger than Lr=',Lengthr
             write(6, *) 'TYPE a value within range ]0,Lr]     '
             write(6, *) '--------------------------------------'
             call beeper
             goto 32
         end if

         if (ClearR.eq.zero) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Negative CLEARANCE at RIGHT SIDE '
             write(6, *) '--------------------------------------'
             call beeper
             goto 32
         end if


 33      IF (BEARING.EQ.2) THEN    ! => for annular seal
               ClearL=ClearO
               YL=Lengthl
               goto 34
         END IF

         write (6,252)         ! Read YL, ClearL
         read (5,*) YL, ClearL
         lhalf=Lengthl
         ClearL=DABS(ClearL)

         if (YL.lt.-lhalf) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Axial Y value smaller than -Ll=-',Lengthl
             write(6, *) 'TYPE a value within range ]0,-Ll]    '
             write(6, *) '--------------------------------------'
             call beeper
             goto 33
         end if
         if (ClearL.eq.zero) then
             write(6, *) 'WARNING:------------------------------'
             write(6, *) 'Negative CLEARANCE at LEFT SIDE '
             write(6, *) '--------------------------------------'
             call beeper
             goto 32
         end if

 34      write (6, 253) ClearO, YR, ClearR, YL, ClearL
         ICSTEP=1
         dumyc='Y'
         write (6, 25)
         CALL BEEPER
         read (5,26) dumyc

      if ((dumyc.eq.'N').or.(dumyc.eq.'n')) GOTO 31

c......................! ICSTEP=1 STEP CLEARANCE HJB


c-----------------------------------------------------------------
 250  format ( ' ', 'INPUT clearance [m] at Y=0 ')
 251  format ( ' ', 'INPUT Axial coord. [m] and clearance [m] for',
     + ' RIGHT STEP')
 252  format ( ' ', 'INPUT Axial coord. [m] and clearance [m] for',
     + ' LEFT STEP')
 253  format (' ', 15('.'),'ICSTEP=1 STEP CLEARANCE',15('.'),/,3X,
     + 'at Y=0.0(center) ,    ClearO=', E12.5E2, 'm',/,3X,
     + 'at YR=',E12.5E2,'m : ClearR=', E12.5E2, 'm', /, 3X,
     + 'at YL=',E12.5E2,'m : ClearL=', E12.5E2, 'm', /, 1X,60('.'))

 240  format ( ' ', 40('.'))
 241  format ( ' ', 'INPUT Axial coord. [m] and clearance [m] for',
     + 1X,I3,' data point')
 243  format ( ' ', 'Bearing AXIAL Length L=', E12.5E2, ' m',
     + /, 'Axial values to describe clearance are on range:',
     + /,'-',E12.5E2,'=-Ll <= y <= Lr:', E12.5E2, 'm', /, 40('.'))
 23   format (' ',I3, E12.5E2, 2X, E12.5E2)
 21   format (' ',' Number of axial coordinates Nj:', I3,/, 3X,
     + 'j  ', 4X, 'Axial Y' , 4X, 'Clearance',/, 3X,
     + '   ', 4X, '  [m]  ' , 4X, '  [m]    ',/, /, 3X,
     + 40('.') )
 75   format (' ','You have violated the upper limited for Nj=',
     + I3,/, 'The maximum Nj allowed in this version is Nj=Nsl=',
     + I3,/, 'To raise this limit, edit params.f , change Nls value,'
     +,/, 'and recompile entire program',
     +/, 'PLEASE RE-ENTER the value')

 1000 format(A40)
 1010 format(' ',A40)

c...............................................................

      END


C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputcell                                                   **
C **                                                                         **
C **  INPUTcell: keyboard input changes in eccentricity and misalignments    **
C **                                                                         **
C *****************************************************************************
C 6/20/95

      SUBROUTINE INPUTCELL

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /HONEY/ HCELL, HCDIM
      DOUBLE PRECISION HCELL, HCDIM
C ----------------------------------------------------------------------------
C --  INPUTCELL code                                                       --
C ----------------------------------------------------------------------------

 904  write (6, *) '---------------------------------------------'
      write (6, *) '(3)INPUT: CELL DEPTH for bearing surface [m]'
      write (6, *) '---------------------------------------------'
      CALL BEEPER

c##   WRITE (6, *) 'HCELL (REAL*8): CELL DEPTH in meters.'
      WRITE (6, 131) HCELL

  131 FORMAT ('$', 'ENTER CELL depth [m]= ',E12.5E2,']: ')

      CALL ENTERVAL(HCELL)
      HCELL=DABS(HCELL)

      END
C..................................................................



C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputexey                                                   **
C **                                                                         **
C **  INPUTEXEY: keyboard input changes in eccentricity and misalignments    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTEXEY

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /TILTPAD/ RSINPK, RCOSPK, IPAD, TILT
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO

      DOUBLE PRECISION EXO,EYO,AXO,AYO,ZO, DUMY, ECC
      DOUBLE PRECISION RSINPK, RCOSPK, IPAD
      INTEGER  TILT

C ----------------------------------------------------------------------------
C --  INPUTEXEY code                                                       --
C ----------------------------------------------------------------------------

 904  write (6, *) '---------------------------------------------'
      write (6, *) '(4)INPUT: Coordinates of journal eccentricity'
      write (6, *) '---------------------------------------------'
      CALL BEEPER

      WRITE (6, *) 'EXO (REAL*8): Dimensionless X eccentricity.'
      WRITE (6, 131) EXO

  131 FORMAT ('$', 'ENTER Exo [Default (Exo)= ',E12.5E2,']: ')

      CALL ENTERVAL(EXO)
C..................................................................

      WRITE (6, *) 'EYO (REAL*8): Dimensionless Y eccentricity.'
      WRITE (6, 132) EYO

  132 FORMAT ('$', 'ENTER Eyo [Default (Eyo)= ',E12.5E2,']: ')

      CALL ENTERVAL(EYO)

          ECC=DSQRT(EXO*EXO+EYO*EYO)        !
          IF (ECC.GE.(1.0D0)) THEN          ! ECC>1
            WRITE (6,179) ECC               !
            CALL BEEPER                     !
c###        GOTO 904                        !
          END IF                            !
                                            !
  179     FORMAT (' ','MAGNITUDE OF DIMENSIONLESS ECCENTRICITY'
     +' ECC=', F7.5, ' IS >=1.0, CAUTION NEEDED')

c ................................................................
c    No Journal misalignment effects for Tilt Pads
c    .............................................
      IF (TILT.EQ.1) THEN
         AXO=0.0D0
         AYO=0.0D0
         ZO=0.0D0
         RETURN
      END IF
c ................................................................
      write (6, *) '--------------------------------------------'
      write (6, *) '(4)INPUT:Angular Misalignment in RADIANS '
      write (6, *) '--------------------------------------------'
      CALL BEEPER

      WRITE (6, *) 'AXO (REAL*8): Misalignment Angle about X axis'
      WRITE (6, 133) AXO

  133 FORMAT ('$', 'ENTER AXO [Default (AXo)= ',E12.5E2,']: ')

      CALL ENTERVAL(AXO)
C.................
      WRITE (6, *) 'AYO (REAL*8): Misalignment Angle about Y axis'
      WRITE (6, 134) AYO
  134 FORMAT ('$', 'ENTER AYO [Default (AYo)= ',E12.5E2,']: ')

      CALL ENTERVAL(AYO)

c.................
      WRITE (6, *) 'ZO: Location Journal AXIS crosses Z axis [m]'
      WRITE (6, 135) ZO
  135 FORMAT ('$', 'ENTER ZO [Default (Zo)= ',E12.5E2,' m]: ')

      CALL ENTERVAL(ZO)


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputops                                                    **
C **                                                                         **
C **  INPUTOPS: keyboard input changes in rpm and pressures                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTOPS

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA!,RPM,PS,PA
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /PDISCH/ Pleft, Pright, Cleft, Cright
      COMMON /IPresLR/ IPRuni, IPLuni, NPRcs, NPLcs
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 PRATIO, CORIF, SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 Pleft, Pright, Cleft, Cright
      INTEGER BEARING, IPRuni, IPLuni, NPRcs, NPLcs
C ----------------------------------------------------------------------------
C --  INPUTOPS code                                                       --
C ----------------------------------------------------------------------------
 906  CONTINUE

      write (6, *) '-------------------------------------------'
      write (6, *) '(6)INPUT:operating parameters:RPM,Ps,PL,PR'
      write (6, *) '-------------------------------------------'
      CALL BEEPER

      WRITE (6, *) 'RPM (REAL*8): Journal rotating speed, [rpm]'
      WRITE (6, 170) RPM

  170 FORMAT ('$', 'ENTER Rpm [Default (rpm)= ',E12.5E2,']: ')

      CALL ENTERVAL(RPM)
      RPM=DABS(RPM)
C
C...............................
C SUPPLY PRESSURE
C...............................


C FOR HYDROSTATIC BEARING OR ANNULAR SEAL

      IF (BEARING.EQ.3) GOTO 192

      WRITE (6, *) 'PS (REAL*8): Pressure supply (Ps), [N/m2]'
      WRITE (6, 180) PS

  180 FORMAT ('$', 'ENTER Ps [Default (Ps)= ',E14.7E2,']: ')

      CALL ENTERVAL(PS)
      PS=DABS(PS)

      WRITE (6, *) ' '
C.................................
C
C EXTERNAL DISCHARGE PRESURES:
C.................................
       IF (BEARING.EQ.2) THEN !=> ANNULAR SEAL
          PLEFT=PS
          IPLUNI=1
          GOTO 194
       END IF

  192 CALL ReadPleft(PLeft)

  194 CALL ReadPright(PRight)

      IF (BEARING.EQ.2) THEN       ! seal
          PA=PRight
      ELSE IF (BEARING.EQ.1) THEN  ! HJB
          PA=DMIN1(Pleft,Pright)
      ELSE IF (BEARING.EQ.3) THEN  ! JB
          PA=DMIN1(Pleft,Pright)
          PS=PA
      END IF

C.............................................................
C
C CAVITATION PRESSURE
C................................

      WRITE (6, *) 'PC (REAL*8): Cavitation pressure(Pc), [N/m2]'
      WRITE (6, 200) PC

  200 FORMAT ('$', 'ENTER Pc [Default (Pc)= ',E14.7E2,']: ')

      CALL ENTERVAL(PC)

      CALL BEEPER
C
C RECESS PRESSURE RATIO
C..................................

      IF (BEARING.GE.2) GOTO 250

      WRITE (6, *) 'PRATIO (REAL*8): Desired or guessed pressure ratio.'
      WRITE (6, *) '                 PRATIO = (Prec-Pa)/(Psupply-Pa)'
      write (6, *) '                 at the concentric position.'
      WRITE (6, 245) PRATIO

  245 FORMAT ('$', 'ENTER Pratio [Default = ', E12.5E2,']: ')

      CALL ENTERVAL(PRATIO)
      PRATIO=DABS(PRATIO)

  250 CALL BEEPER


      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputconv                                                   **
C **                                                                         **
C **  INPUTCONV: keyboard input changes in iterations and convergence params **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTCONV

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /SOURCEA/ PRATIO, CORIF,SMASS,MPEPS,PREPS,MMP,SFLOW
      COMMON /SOURCEB/ ITER, ITMAX, ITPMAX
      COMMON /BTYPE/ BEARING
      COMMON /THERMAL/ ALFT, UC, TC, Ec
      DOUBLE PRECISION PRATIO, CORIF, SMASS,MPEPS,PREPS,MMP,SFLOW,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 ALFT, UC, TC, Ec

      INTEGER ITER, ITMAX,ITPMAX,BEARING
C ----------------------------------------------------------------------------
C --  INPUTRCONV code                                                       --
C ----------------------------------------------------------------------------

 908  write (6, *) '--------------------------------------------------'
      write (6, *) '(8)INPUT: Numerical Convergence Parameters'
      write (6, *) '--------------------------------------------------'
      CALL BEEPER

      WRITE (6, *) 'ITMAX (INT*4): MAX # of iterations for convergence
     +on lands.'
      WRITE (6, 280) ITMAX

  280 FORMAT ('$', 'ENTER Itmax [Default Itmax= ',I3,']: ')

      CALL ENTERINT(ITMAX)

C......................................................................
C
C ITPMAX: for HJB only
C.....................
      IF (BEARING.EQ.1) THEN
      WRITE (6, *) 'ITPMAX (INT*4): MAX # of iterations for convergence
     +on Precess'
      WRITE (6, 290) ITPMAX

  290 FORMAT ('$', 'ENTER Itpmax [Default Itpmax= ',I2,']: ')

      CALL ENTERINT(ITPMAX)
      END IF
C.............................................................
C
C ALFU: under-relax parameter for convergence of momentum eqs.
C       typically large for large journal speeds
C.........................

      WRITE (6, *) 'ALFU (REAL*8): Under_relaxation for moment equations
     +(Au); 0 < Au < 1'
      WRITE (6, 300) ALFU

  300 FORMAT ('$', 'ENTER Alfu [Default (Au)= ',E12.5E2,']: ')

      CALL ENTERVAL(ALFU)
      ALFU=DABS(ALFU)
      ALFU=DMIN1(ALFU,1.0D0)

C..............................................................
C
C ALFP: Under-relax parameter for pressure correction equation
C       typically low =0.5 for pure hydrostatic operation
C....................................

      WRITE (6, *) 'ALFP (REAL*8): Under_relaxation for pressure equatio
     +n (Ap); 0< Ap <1'
      WRITE (6, 310) ALFP

  310 FORMAT ('$', 'ENTER Alfp [Default (Ap)= ',E12.5E2,']: ')

      CALL ENTERVAL(ALFP)
      ALFP=DABS(ALFP)
      ALFP=DMIN1(ALFP,1.0D0)

C................................................................
C
C ALFT: Under-relax parameter for energy equation: typically =0.8
C................................................................

      WRITE (6, *) 'ALFT (REAL*8): Under_relaxation for energy equatio
     +n (At); 0< At <1'
      WRITE (6, 311) ALFT

  311 FORMAT ('$', 'ENTER Alft [Default (At)= ',E12.5E2,']: ')

      CALL ENTERVAL(ALFT)
      ALFT=DABS(ALFT)
      ALFT=DMIN1(ALFT,1.0D0)

C.......................................................................
C
C MPEPS: Max Error in CV mass/Bearing Inflow
C.................................

      WRITE (6, *) 'MPEPS (REAL*8): Convergence criteria for sum of mass
     + error flows on lands; a'
      WRITE (6, *) '                percent of total bearing inlet flow'
      WRITE (6, 320) MPEPS

  320 FORMAT ('$', 'ENTER Mpeps [Default (Mpeps)= ',E12.5E2,']: ')

      CALL ENTERVAL(MPEPS)
      MPEPS=DABS(MPEPS)

      PREPS=.0005D0   ! difference in DIM recess pressures.


C......................................................................
C
C SFLOW: MAX Difference Rec Flows/ Rec flow on 2 interations
C.....................................
      IF (BEARING.EQ.1) THEN  ! => FOR HJBS
C   !.........................!
      WRITE (6, *) 'SFLOW (REAL*8): Convergence criteria for flow on rec
     +ess'
      write (6, *) 'SFLOW: MAX percent recess flow difference on 2 iter'
      WRITE (6, 340) SFLOW


  340 FORMAT ('$', 'ENTER Sflow [Default (Sflow)= ',E12.5E2,']: ')

C   !.........................!
      ELSE
C   !.........................!

      WRITE (6, *) 'PEPS (REAL*8): Convergence criteria for pressure cal
     +culations'
      WRITE (6, 350) SFLOW
  350 FORMAT ('$', 'ENTER PEPS [Default (Peps)= ',E12.5E2,']: ')

C   !.........................!
      END IF
C   !.........................!

      CALL ENTERVAL(SFLOW)
      SFLOW=DABS(SFLOW)

c...............................................................
      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputprops                                                  **
C **                                                                         **
C **  INPUTPROPS: keyboard input for fluid props at 2 pressures              **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTPROPS

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PROPS12/ P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2
      COMMON /IOPROP/ RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                CPS,THS,BKS
      COMMON /PDISCH/ Pleft, Pright, Cleft, Cright
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /RECPAR/ HRECU, VSUP, BETA
      COMMON /LIQUID/ TEMPK, Vsound, IF, IL
      COMMON /OILCOEF/ Talpha
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION P1PROP,P2PROP,RHO1,RHO2,EMU1,EMU2,
     +                 RHOS,EMUS,RHOA,EMUA,RHOle,EMUle,RHOri,EMUri,
     +                 Pleft, Pright, Cleft, Cright,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 HRECU, VSUP, BETA ,CPS,THS,BKS, Talpha,
     +                 TEMPK, Vsound
      INTEGER BEARING  ,IF, IL
C ----------------------------------------------------------------------------
C --  INPUTPROPS code                                                       --
C ----------------------------------------------------------------------------
c     IF=12, unknown fluid type for MIPROPS,
c            USER Inputs fluid properties at two different pressures
c     ..........................................................

c    !======================!
      IF (IF.EQ.12) THEN
c    !======================! unknown barotropic fluid

      write (6, *) 'IF=12 ENTER Fluid Viscosity/Density at 2 pressures'
      write (6, *) '--------------------------------------------------'
  150 write (6, *) 'FLUID Props. as linear functions of Pressure'
      write (6, *) '--------------------------------------------------'

      CALL BEEPER
      write (6, *) ' '
      write (6, *) 'ENTER a value of pressure, say P1 [N/m2]'
      WRITE (6, 180) P1prop
  180 FORMAT ('$', 'ENTER Press#1 [Default (P1)= ',E14.7E2,' N/m2]:')
      CALL ENTERVAL(P1prop)
      P1prop=DABS(P1prop)

      write (6, *) ' '

  250 WRITE (6, *) 'EMU1 (REAL*8):Absolute viscosity (u) at P1, [Ns/m2]'
      WRITE (6, 350) EMU1

  350 FORMAT ('$', 'ENTER Emu1 [Default (u)= ',E12.5E2,']: ')
      CALL ENTERVAL(EMU1)
      EMU1=DABS(EMU1)

  260 WRITE (6, *) 'RHO1 (REAL*8) Density (p) at P1, [kg/m3]'
      WRITE (6, 360) RHO1

  360 FORMAT ('$', 'ENTER Rho1 [Default (p)= ',E14.7E2,']: ')
      CALL ENTERVAL(RHO1)
      RHO1=DABS(RHO1)

c...............................................................

      CALL BEEPER
      write (6, *) ' '
      write (6, *) 'ENTER a value of pressure, say P2 [N/m2]'
      write (6, *) 'IF P2=P1 then fluid is incompressible & isoviscous'
      write (6, *) '--------------------------------------------------'
      write (6, *) ' '

      WRITE (6, 400) P2prop
  400 FORMAT ('$', 'ENTER Press#2 [Default (P2)= ',E14.7E2,' N/m2]:')
      CALL ENTERVAL(P2prop)
      P2prop=DABS(P2prop)

c    !.............................!
      IF (P2prop.eq.P1prop) THEN
c    !.............................!
         write (6, 420)
  420    format (' ',10X,40('-'),/,11X,
     +           'P2=P1 for properties, THEN =>',/,11X,
     +           'fluid is incompressible and isoviscous',/11X,40('-'))
         write (6, *) ' '
         CALL BEEPER
         RHO2=RHO1
         RHOS=RHO1
         RHOle=RHOS
         RHOri=RHOS
         RHOA=RHOS
         EMU2=EMU1
         EMUS=EMU1
         EMUle=EMUS
         EMUri=EMUS
         EMUA=EMUS
         BETA=0.0D0
         RETURN
c    !.............................!
      END IF
c    !.............................!

  450 WRITE (6, *) 'EMU2 (REAL*8):Absolute viscosity (u) at P2, [Ns/m2]'
      WRITE (6, 550) EMU2

  550 FORMAT ('$', 'ENTER Emu2 [Default (u)= ',E12.5E2,']: ')

      CALL ENTERVAL(EMU2)
      EMU2=DABS(EMU2)

  460 WRITE (6, *) 'RHO2 (REAL*8):Density (p) at P2, [kg/m3]'
      WRITE (6, 560) RHO2

  560 FORMAT ('$', 'ENTER Rho2 [Default (p)= ',E14.7E2,']: ')

      CALL ENTERVAL(RHO2)
      RHO2=DABS(RHO2)

c     ...... calculate props. at PS, Pleft,Pright from linear relations

      RHOs=RHO2+(RHO1-RHO2)*(Ps-P2prop)/(P1prop-P2prop)
      RHOle=RHO2+(RHO1-RHO2)*(Pleft-P2prop)/(P1prop-P2prop)
      RHOri=RHO2+(RHO1-RHO2)*(Pright-P2prop)/(P1prop-P2prop)
      RHOa=RHO2+(RHO1-RHO2)*(Pa-P2prop)/(P1prop-P2prop)

      EMUs=EMU2+(EMU1-EMU2)*(Ps-P2prop)/(P1prop-P2prop)
      EMUle=EMU2+(EMU1-EMU2)*(Pleft-P2prop)/(P1prop-P2prop)
      EMUri=EMU2+(EMU1-EMU2)*(Pright-P2prop)/(P1prop-P2prop)
      EMUa=EMU2+(EMU1-EMU2)*(Pa-P2prop)/(P1prop-P2prop)

c    !======================!
      ELSE IF (IF.EQ.7) THEN
c    !======================!Oil properties

      write (6, *) '--------------------------------------------------'
      write (6, *) 'IF=7 ENTER OIL Properties at SUPPLY Condition with'
      write (6, *) 'P=Psupply and T=To=Tsupply(Tempk).  Note that only'
      write (6, *) 'The viscosity varies:   MU=MUo*EXP[-Talpha*(T-To)]'
      write (6, *) '--------------------------------------------------'
      write (6, *) 'Note:  Following values are provided for reference'
      write (6, *) '--------------------------------------------------'
      write (6, *) 'EMUS=0.0172D0   ! Viscosity at Supply P,T (Ns/m^2)'
      write (6, *) 'RHOS=865.12D0   ! Density at Supply P & T (kg/m^3)'
      write (6, *) 'CPS=2000.D0     ! Specific  heat  of  oil (J/kg.K)'
      write (6, *) 'BetaTs=2.0D-4   ! Volumetric expansion coeff.(1/K)'
      write (6, *) 'THS=0.15D0      ! Oil thermal conductivity (W/m.K)'
      write (6, *) 'Talpha=0.0303D0 ! Temperature-viscosity coef.(1/K)'
      write (6, *) 'BETA=4.148D-10  ! Compressibility  factor  (m^2/N)'
      write (6, *) '--------------------------------------------------'

      CALL BEEPER

c      write (6, *) ' '

 3250 WRITE (6, *) 'EMUS (REAL*8): Viscosity (u) at T=Tempk, [Ns/m2]'
      WRITE (6,3350) EMUS

 3350 FORMAT ('$', 'ENTER Emus [Default (u)= ',E12.5E2,']: ')
      CALL ENTERVAL(EMUS)
      EMUS=DABS(EMUS)

 1260 WRITE (6,*)'RHOS (REAL*8) Density (p) at Supply Cond., [kg/m3]'
      WRITE (6,1360) RHOS

 1360 FORMAT ('$', 'ENTER Rhos [Default (p)= ',E14.7E2,']: ')
      CALL ENTERVAL(RHOS)
      RHOS=DABS(RHOS)

 1450 WRITE (6, *) 'CPS (REAL*8):Specific Heat (Cp), [J/kg.K]'
      WRITE (6,1550) CPS

 1550 FORMAT ('$', 'ENTER CPS [Default (Cp)= ',E12.5E2,']: ')

      CALL ENTERVAL(CPS)
      CPS=DABS(CPS)

 1460 WRITE (6, *) 'BetaTs (REAL*8): Thermal Expansion Coefficient (Bt),
     +[1/K-deg]'
      WRITE (6,1560) BKS

 1560 FORMAT ('$', 'ENTER BetaTs [Default (Bt)= ',E14.7E2,']: ')

      CALL ENTERVAL(BKS)
      BKS=DABS(BKS)

 1461 WRITE (6, *) 'THS (REAL*8): Thermal Conductivity (K), [W/m.K]'
      WRITE (6,1561) THS

 1561 FORMAT ('$', 'ENTER THS [Default (K)= ',E14.7E2,']: ')

      CALL ENTERVAL(THS)
      THS=DABS(THS)

 1250 WRITE (6, *) 'Talpha (REAL*8): Temp.-viscosity coefficient,[1/K]'
      WRITE (6,1350) Talpha

 1350 FORMAT ('$', 'ENTER Talpha [Default = ',E12.5E2,']: ')
      CALL ENTERVAL(Talpha)
      Talpha=DABS(Talpha)

      END IF
C...................................................................
C
C Compressibility parameter Beta=1./Bulk modulus
C.................................
C##      IF (BEARING.GE.3) RETURN

      WRITE (6, *) 'BETA (REAL*8): Compressibility parameter(B),[m2/N]'
      WRITE (6, 201) BETA

  201 FORMAT ('$', 'ENTER Beta [Default (B)= ',E14.7E2,']: ')

      CALL ENTERVAL(BETA)
      BETA=DABS(BETA)

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputloss                                                   **
C **                                                                         **
C **  INPUTLOSS: keyboard input for edge loss coefficients                  **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTLOSS

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /LOSPAR/ LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr
      COMMON /LOSPAD/ LOSleadP, KLOSpad
      COMMON /PARAM3/ EMU,RHO,PC,CD,DORIF,LOSXSI,ALPHA,RPM,PS,PA
      COMMON /PDISCH/ Pleft, Pright, Cleft, Cright
      COMMON /JET/ ANGLEJ, LOCJET, CJET, DPJET
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION LOSXSIxu, LOSXSIxd, LOSXSIyl, LOSXSIyr,
     +                 LOSleadP, KLOSpad,
     +                 EMU,RHO,RPM,PS,PA,PC,CD,DORIF,LOSXSI,ALPHA,
     +                 Pleft, Pright, Cleft, Cright,
     +                 ANGLEJ, LOCJET, CJET, DPJET

      INTEGER BEARING
C ----------------------------------------------------------------------------
C --  INPUTLOSS code                                                       --
C ----------------------------------------------------------------------------

 907  write (6, *) '---------------------------------------------'
      write (6, *) '(7)INPUT: Cd, Loss and SEAL Parameters       '
      write (6, *) '---------------------------------------------'
      CALL BEEPER

      IF (BEARING.GE.2) THEN
         CD=0.0D0
         ANGLEJ=0.0D0
         LOCJET=0.5D0
         LOSXSIxu=0.0D0
         LOSXSIxd=0.0D0
         GOTO 330
      END IF
C..............................................................
C ORIFICE DISCHARGE COEFFICIENT, ORIFICE ANGLE & LOCATION
C..............................................................

      WRITE (6, *) 'CD (REAL*8): Orifice discharge coefficient.'
      WRITE (6, 210) CD

  210 FORMAT ('$', 'ENTER Cd [Default (Cd)=',F7.4,']: ')

      CALL ENTERVAL(CD)
      CD=DABS(CD)

      WRITE (6, *) 'ANGLEJ (REAL*8): Orifice Supply Angle:'
      WRITE (6, *) '0 deg: RADIAL, +90: TANGENT against shaft rotation'
      WRITE (6, *) '     -90 <= ANGLEJ <= +90 degrees'
      WRITE (6, 211) ANGLEJ

cvvvv added on 7/3/95
  211 FORMAT ('$', 'ENTER Angle Jet[Default=',F7.3,'(degrees):]')

      CALL ENTERVAL(ANGLEJ)
      IF (ANGLEJ.GT.90.0D0)  ANGLEJ=90.0D0
      IF (ANGLEJ.LT.-90.0D0) ANGLEJ=-90.0D0

      WRITE (6, *) 'LOCJET (REAL*8) circumferential location of orifice'
      WRITE (6, *) 'relative to RECESS, i.e 0 = left edge (upstream)'
      WRITE (6, *) '0.5 = middle, 1.0 = right edge (downstream)'
      WRITE (6, 212) LOCJET

  212 FORMAT ('$', 'ENTER Location [Default=',F7.3,' :]')

      CALL ENTERVAL(LOCJET)
      LOCJET=DABS(LOCJET)
      IF (LOCJET.GT.1.0D0) LOCJET=1.0D0

c^^^ added on 7/3/95

      CALL BEEPER

C...............................................................
C
C LOSXSI: Added Inertial entrance loss coefficient
C         Value used is (1+LOSXSI)/2.
C.................................................
      write (6, *) ' Used Recess edge entrance loss is k=(1+Xsi)/2'
      WRITE (6, *) ' '
      WRITE (6, *) 'XSIx,XSIy (REAL*8): Circumf. & Axial loss coeffs.'
      write (6, *) ' '

      WRITE (6, 229) LOSXSIxu
  229 FORMAT ('$', 'ENTER XSIxUPSTREAM[Default(Xsixu)= ',F7.3,']: ')

      CALL ENTERVAL(LOSXSIxu)

      WRITE (6, 230) LOSXSIxd
  230 FORMAT ('$', 'ENTER XSIxDOWNSTREAM[Default(Xsix)= ',F7.3,']: ')

      CALL ENTERVAL(LOSXSIxd)

C ..............................................................!
  330 CONTINUE

      IF (BEARING.EQ.3) THEN   ! FOR PLAIN BEARING, INERP=0
         LOSXSIyl=0.0D0
         LOSXSIyr=0.0D0
         CLEFT=0.0D0
         CRIGHT=0.0D0
         GOTO 433
      ELSE IF (BEARING.EQ.2) THEN  ! FOR SEAL
         LOSXSIyl=0.0D0
         WRITE (6,231) LOSXSIyr
  231 FORMAT ('$', 'ENTER LOSS COEFFICIENT Xsi[Default= ',F7.3,']: ')
         CALL ENTERVAL(LOSXSIyr)
         GOTO 333
      END IF
C ..............................................................!

      WRITE (6, 232) LOSXSIyl
  232 FORMAT ('$', 'ENTER LEFT XSIy AXIAL[Default= ',F7.3,']: ')

      CALL ENTERVAL(LOSXSIyl)

  332 WRITE (6, 233) LOSXSIyr
  233 FORMAT ('$', 'ENTER RIGHT XSIy AXIAL[Default= ',F7.3,']: ')

      CALL ENTERVAL(LOSXSIyr)
      CALL BEEPER
      CALL PAUSE

  333 WRITE (6, *) '.......................................... '
      write (6, *) 'ENTER: EXIT Discharge Seal Coefficients    '
      write (6, *) '       Cseal=0.0 no seal                   '
      write (6, *) '       ....................................'

c    !.....................!
      IF (BEARING.EQ.2) THEN   !=> annular seal
c    !.....................!
          CLEFT=0.0D0
c    !.....................!
      ELSE
c    !.....................!
      WRITE (6, 243) CLEFT
  243 FORMAT ('$', 'ENTER LEFT SEAL [Default= ',F9.3,']: ')
      CALL ENTERVAL(CLEFT)
c    !.....................!
      END IF
c    !.....................!

      WRITE (6, 244) CRIGHT
  244 FORMAT ('$', 'ENTER RIGHT SEAL [Default= ',F9.3,']: ')

      CALL ENTERVAL(CRIGHT)

      WRITE (6, *) '.......................................... '

C ..............................................................!

  433 LOSXSI=LOSXSIyR
C

C.................................................................
C
C Recovery factor for leading edge RAM Pressure=
C    (1/2)LOSleadPxRHO(OmegaxR/2)**2
c
C.................................................................
      IF (BEARING.EQ.2) THEN
         LOSLEADP=0.0D0
      ELSE
C .........#### HERE ONLY IF IFULL=0 FOR PAD BEARING

         CALL BEEPER
      WRITE (6, *) ' PAD LEADING EDGE: RECOVERY FACTOR (LOSleadP)'
      WRITE (6, *) ' for inlet pressure a fraction of (OMEGAxR/2)**2'
      WRITE (6, *) ' ==> ONLY FOR PAD BEARINGS < 360 deg'

      WRITE (6, 250) LOSleadP
  250 FORMAT ('$', 'ENTER LOSleadP [Default = ',F7.3,']: ')

          CALL ENTERVAL(LOSLEADP)

          CALL BEEPER

      END IF

C...............................................................
C
C ALPHA: SWIRL FRACTION
C...........................

      IF (BEARING.EQ.3) THEN  !no swirl effects for JBs.
         ALPHA=0.50
         CALL BEEPER
         RETURN
      END IF


      WRITE (6, *) 'ALPHA (REAL*8): Swirl fraction at entrance.'
      write (6, *) 'ALPHA = Umean circ./(Omega x R)'
      WRITE (6, 260) ALPHA

  260 FORMAT ('$', 'ENTER Alpha [Default = ',F6.3,']: ')

      CALL ENTERVAL(ALPHA)

      CALL BEEPER
      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputthd                                                    **
C **                                                                         **
C **  INPUTTHD: Options for Thermal Analysis in Fluid Film                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INPUTTHD

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /LIQUID/ TEMPK, VSOUND, IF, IL
      COMMON /TISOBJ/ TSHAFT, TSTATOR
      COMMON /THERMID/ ISOTH
      COMMON /CONT/ IFF
      COMMON /IDTHDPAD/ ITPAD
      COMMON /RADHEAT/ TBOUT, THERMALK, ROUTER, HKB

      DOUBLE PRECISION TEMPK, VSOUND, TSHAFT, TSTATOR,
     +                 TBOUT, THERMALK, ROUTER, HKB

      INTEGER ISOTH, IFF ,ITPAD, IF, IL

C ----------------------------------------------------------------------------
C --  INPUTTHD code                                                         --
C ----------------------------------------------------------------------------

 910  CONTINUE
      IF(IFF.EQ.12) THEN
         ISOTH=1
         ITPAD=0       ! Isothermal pad bearings
         write (6, *) '*****************************************'
         write (6, *) '*  For linear property fluids, ISOTH=1  *'
         write (6, *) '*  TSHAFT=TSTATOR=TSUPPLY               *'
         write (6, *) '*****************************************'
         TSHAFT=TEMPK
         TSTATOR=TEMPK
         CALL BEEPER
         CALL PAUSE
         RETURN
      END IF

      write (6, *) 'ISOTH = 1;    Isothermal (Barotropic) fluid film (T
     &=Ts)'
      write (6, *)
      write (6, *) 'ISOTH =-1;    Adiabatic Bounding Surfaces (Qb=Qj=0)'
      write (6, *) 'ISOTH = 0;    Isothermal journal & stator: Tj ; Tb'
      write (6, *) 'ISOTH = 2;    Adiabatic journal (Qj=0) & Isothermal
     &stator (Tb)'
      write (6, *) 'ISOTH = 3;    Isothermal journal(Tj) & Adiabatic
     & stator (Qb=0)'
      write (6, *) 'ISOTH = 4;    Adiabatic journal (Qj=0) & Radial heat
     & flow on bearing'
      write (6, *) 'ISOTH = 5;    Isothermal journal (Tj) &  Radial heat
     & flow on bearing'

      write (6, *) '---------------------------------------------------'
      write (6, *) '(1) For liquid hydrogen (LH2), it is recommended to'
      write (6, *) '    perform a  BAROTROPIC  fluid  calculation first'
      write (6, *) '    (ISOTH=1), and  then use the  results as a good'
      write (6, *) '    starting guess for the THERMAL analysis;'
      write (6, *) '(2) For HJBs with high  pressure gradients or large'
      write (6, *) '    mass flow rates, the ADIABATIC (ISOTH=-1) model'
      write (6, *) '    is an excellent model.'
      write (6, *) '---------------------------------------------------'
      WRITE (6, 400) ISOTH
  400 FORMAT ('$', 'Enter ISOTH [Default ISOTH= ',I2,']: ')
      CALL ENTERINT(ISOTH)
      CALL BEEPER
C##   CALL PAUSE

      IF ((ISOTH.LT.-1).or.(ISOTH.GT.5)) THEN
         GOTO 910
      END IF

c    !--------------------------- ISOTHERMAL ANALYSIS
      IF (ISOTH.EQ.1) THEN
c    !.....................!
         TSHAFT=TEMPK             ! NOT REALLY NEEDED
         TSTATOR=TEMPK            ! temporal only

c    !--------------------------- ADIABATIC JOURNAL / ISOTHERMAL STATOR
      ELSE IF (ISOTH.EQ.2) THEN
c    !.....................!
         WRITE (6,171) TEMPK
         WRITE (6,172) TSTATOR
         CALL ENTERVAL(TSTATOR)
         TSTATOR=DABS(TSTATOR)
         TSHAFT=TEMPK            ! temporal only

c    !--------------------------- ADIABATIC STATOR / ISOTHERMAL JOURNAL
      ELSE IF (ISOTH.EQ.3) THEN
c    !.....................!
         WRITE (6,171) TEMPK
         WRITE (6,173) TSHAFT
         CALL ENTERVAL(TSHAFT)
         TSHAFT=DABS(TSHAFT)
         TSTATOR=TEMPK            ! temporal only

c    !--------------------------- ADIABATIC JOURNAL / RADIAL HEAT FLOW STATOR
      ELSE IF (ISOTH.EQ.4) THEN
c    !.....................!
         WRITE (6,171) TEMPK
         WRITE (6,174) THERMALK   ! bearing solid thermal conductivity
         CALL ENTERVAL(THERMALK)
         THERMALK=DABS(THERMALK)
         WRITE (6,175) ROUTER     ! bearing or pad outer radius ROUTER > R
         CALL ENTERVAL(ROUTER)
         ROUTER=DABS(ROUTER)
         WRITE (6,176) TBOUT      ! bearinng outer temperature
         CALL ENTERVAL(TBOUT)
         TBOUT=DABS(TBOUT)
         TSHAFT=TEMPK             ! temporal only
         TSTATOR=TEMPK            ! temporal only

c    !--------------------------- ISOTHERMAL JOURNAL / RADIAL HEAT FLOW STATOR
      ELSE IF (ISOTH.EQ.5) THEN
c    !.....................!
         WRITE (6,171) TEMPK
         WRITE (6,173) TSHAFT     ! journal temperature
         CALL ENTERVAL(TSHAFT)
         TSHAFT=DABS(TSHAFT)
         WRITE (6,174) THERMALK   ! bearing solid thermal conductivity
         CALL ENTERVAL(THERMALK)
         THERMALK=DABS(THERMALK)
         WRITE (6,175) ROUTER     ! bearing or pad outer radius ROUTER > R
         CALL ENTERVAL(ROUTER)
         ROUTER=DABS(ROUTER)
         WRITE (6,176) TBOUT      ! bearinng outer temperature
         CALL ENTERVAL(TBOUT)
         TBOUT=DABS(TBOUT)
         TSTATOR=TEMPK            ! temporal only

c    !--------------------------- ISOTHERMAL JOURNAL & BEARING
      ELSE IF (ISOTH.EQ.0) THEN
c    !.....................!
         WRITE (6,171) TEMPK
         WRITE (6,172) TSTATOR
         CALL ENTERVAL(TSTATOR)
         TSTATOR=DABS(TSTATOR)
         WRITE (6,173) TSHAFT
         CALL ENTERVAL(TSHAFT)
         TSHAFT=DABS(TSHAFT)
c    !--------------------------- ADIABATIC JOURNAL & BEARING SURFACES
      ELSE IF (ISOTH.EQ.-1) THEN
c    !.....................!
         TSHAFT=TEMPK            ! TEMPORAL ONLY, CHANGED ON CALCULATIONS
         TSTATOR=TEMPK
c    !---------------------------
      END IF
c    !.....................!

      RETURN



 171  FORMAT ('$', 'SUPPLY TEMPERATURE= ',F9.3,' K]: ',/)
 172  FORMAT ('$', 'ENTER STATOR or BEARING TEMPERATURE [Def= ',F9.3,
     +             ' K]: ')
 173  FORMAT ('$', 'ENTER SHAFT or  JOURNAL TEMPERATURE [Def= ',F9.3,
     +             ' K]: ')
 174  FORMAT ('$', 'ENTER BEARING or PAD MATERIAL CONDUCTIVITY [Def= ',
     +             E12.5E4,' (Watts/(m . degK)]: ')
 175  FORMAT ('$', 'ENTER BEARING or PAD OUTER RADIUS [Def= ',
     +             F9.3,' (m)]: ')
 176  FORMAT ('$', 'ENTER BEARING or PAD OUTER TEMPERATURE [Def= ',
     +             F9.3,' K]: ')
      END

C *****************************************************************************

      SUBROUTINE ENTERVAL(VALUE)

C *****************************************************************************

      IMPLICIT NONE
      DOUBLE PRECISION VALUE, DUMY
      CHARACTER*40 DUMYC
      INTEGER IOS
c .............................................................................
c Reads from KBD a double precision number
c .............................................................................
C Unit 42 is used as a temporal place for translations between formats.
C.........................................................................

      READ (5, 1000) DUMYC
      IF (DUMYC.NE.' ') THEN
          WRITE (42, 1010, IOSTAT=IOS, ERR=1234) DUMYC
          BACKSPACE (UNIT=42)
          READ (42, *, IOSTAT=IOS, ERR=1234) DUMY
          VALUE=DUMY
      END IF
      WRITE (6, *) ' '

c...............................................................

      RETURN
c...............................................................
 1000 FORMAT (40A)
 1010 FORMAT (' ', 40A)
C                                           !
C                                           !
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) 'I/O ERROR ON INPUT DATA (unit 42)'
      RETURN

c...............................................................
      END

C *****************************************************************************

      SUBROUTINE ENTERINT(VALUE)

C *****************************************************************************

      IMPLICIT NONE
      CHARACTER*40 DUMYC
      INTEGER IOS, VALUE, DUMYI
c .............................................................................
c Reads from KBD an integer number
c .............................................................................
C Unit 42 is used as a temporal place for translations between formats.
C.........................................................................

      READ (5, 1000) DUMYC
      IF (DUMYC.NE.' ') THEN
          WRITE (42, 1010, IOSTAT=IOS, ERR=1234) DUMYC
          BACKSPACE (UNIT=42)
          READ (42, *, IOSTAT=IOS, ERR=1234) DUMYI
C###      VALUE=IABS(DUMYI)
          VALUE=DUMYI
      END IF
      WRITE (6, *) ' '

c...............................................................

      RETURN
c...............................................................
 1000 FORMAT (40A)
 1010 FORMAT (' ', 40A)
C                                           !
C                                           !
 1234 CALL DECODIOS(IOS)
      WRITE (6, *) 'I/O ERROR ON INPUT DATA (unit 42)'
      RETURN

c...............................................................
      END
C *****************************************************************************
C LAST REVISED 2/18/94 BY Dr. Luis San Andres
C              7/3/95  for angled orifice injection
C *****************************************************************************