c  ORIGINAL:
c  This version uses local pressure for deformation on foil element
c  -----------------------------------------------------------------
c  9/30/93: modified to deform trailing edge of foil pad 
c  9/7/93: Set h at exit from upstream value, i.e. h(z=L)=h(z=L-Dz)
c  9/2/93: gone back to basics, give deformation only for P > Pcav
c  started july, 1993
c
c ######   ####      #    #        ####   #    #  #####    ####           ######
c #       #    #     #    #       #       #    #  #    #  #               #
c #####   #    #     #    #        ####   #    #  #####    ####           #####
c #       #    #     #    #            #  #    #  #    #       #   ###    #
c #       #    #     #    #       #    #  #    #  #    #  #    #   ###    #
c #        ####      #    ######   ####    ####   #####    ####    ###    #
c
c
c  hydrofoil.f program. Copyright LuisSanAndres/Texas A&M University 1994.
c
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Filmc                                                       **
C **                                                                         **
C **  FILMC: Updates film thickness due to local compliance effect.          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE FILMC

      IMPLICIT NONE

      INCLUDE 'params.f' 

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /DXVEC/ DXP(MAXNXT), DXU(MAXNXT), SUW(MAXNXT),SUE(MAXNXT)
      COMMON /DYVEC/ DYP(-MAXNYI:MAXNYI), DYV(-MAXNYI:MAXNYI),
     +               SVN(-MAXNYI:MAXNYI), SVS(-MAXNYI:MAXNYI)
      COMMON /PARRAY/  P(MAXNXT,-MAXNYI:MAXNYI)
      COMMON /HOFILM/ HPO(MAXNXT,-MAXNYI:MAXNYI), 
     +                HUO(MAXNXT,-MAXNYI:MAXNYI),
     +                HVO(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /HFILM/  HP(MAXNXT,-MAXNYI:MAXNYI), 
     +                HU(MAXNXT,-MAXNYI:MAXNYI),
     +                HV(MAXNXT,-MAXNYI: MAXNYI)
      COMMON /FACTORS/ REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT
      COMMON /NODES/ NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBSYM/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

C ----------------------------------------------------------------------------
C --  Global variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION DXP, DXU, SUW, SUE,DYP, DYV, SVN, SVS, 
     +                 P, HPO, HUO, HVO, HP, HU, HV,
     +                 REP,REY,MSPEED,SPEED,SWIRL,PCAV,ALFU,BETU,ALFP,
     +                 AC, ETA, PBACK
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL,
     +        ISYM, ICSTEP, BEARING, LIFT
C ----------------------------------------------------------------------------
C --  Local variable declarations                                          --
C ----------------------------------------------------------------------------
      DOUBLE PRECISION PP, PU, PV, PCCV
      INTEGER J, JP1, JM1, K, KM1, Kstart, Kend, Icav, Ilarge
C ----------------------------------------------------------------------------
C --  FILMC code                                                            --
C ----------------------------------------------------------------------------
c    General film thickness for pad with elastic matrix is equal to:
c
c    H/C = Ho/C  + ac (P-Pback)
c
c    where ac is a dimensionless compliance coefficient.
c
c    the local bearing surface deformation is equal to:
c
c    d =(P-Pback) / Ke 
c
c    with Ke = stiffness of elastic foundation / unit area [N/m/m2]
c
c    here Pback is the pressure beneath the compliance surface = Charact. pressure
c
c    since p=(P-P*)/Psa  then, in dimensionless form, the local deflection is
c
c    d/C* = { PSA X/ Ke x C* } p = ac (p-pback)
c  
c    ON Initopsp.f : SUB CPARAM I have set PBACK=PA always !
c
c    AND:   
c
c    Ho = C(Z) + { ex+dy(Z-Zo) } cos(X) + { ey-dx(Z-Zo) } sin(X)
c
c            - rp cos(X-Xpivot) - R x Delta sin(X-Xpivot)
c
c    where C(Z) pad machined clearance.
c          rp   = C - Cm : pad preload  [m], Cm: bearing clearance
c          Delta tilt-pad rotation about pivot located at Xpivot.
c          
c    and OFFSET=(Xpivot-Xleading)/Pad_size, TYP = 0.5
c    
c    ex,ey are journal center displacements [m]
c    dx,dy are journal axis rotations [rad] about Zo.
c
c    USER Provides the value of compliance coefficient ac.
c
c    NOTE that if Ps & P* or PSA change, SO DOES ac in this formulation.
c    however ac is INPUT by the user
c
c    IMPORTANT Limitation
c    In this version, IF P <=Pback then as usual we 
c    do not allow deformation of the foil bearing surface. IN this case, 
c    the foil should lift and conform to the runner 
c    surface. However, we do not have a model for this yet.
c
c    This is a very simple model and it will need a lot of improvement
c    in the future.
c
c    A NOTE on cavitation:
c    In principle, cavitation could be formulated such that when
c    P<=Pcav then H=Hupstream in order to satisfy mass continuity
c    This approach does not work very well since it leads to
c    spurius oscillations on the Pfield and lack of convergence.
C
C ----------------------------------------------------------------------------
C --  FILMC code                                                            --
C ----------------------------------------------------------------------------
C NOTE: On first call to film.f HP=HPO, etc.

      
C   !............................................ ! Rigid surface bearing
      IF (AC.EQ.0.0D0) RETURN                     ! NO need to deform
C   !............................................ ! bearing surface
C                                                 ! ::::::::::::::::::::::::::::
      Ilarge=0              ! indicator of H < 0.0D0

      Kstart=1
      IF (ISYM.EQ.0) Kstart=-NYI   
      IF (BEARING.EQ.2) Kstart=1

c    !..............................! AT P control volumes
      DO K= Kstart, NYI  
        DO J=1, NXT
          PP=P(J,K)       
         IF (PP.GT.PBACK) THEN
          HP(J,K)=HPO(J,K)+AC*(PP-PBACK)
          IF (HP(J,K).LE.0.0D0) Ilarge=1
         ELSE
          HP(J,K)=HPO(J,K)                
         END IF
        END DO
      END DO

c    !.............................! AT U control volumes
      DO K = Kstart, NYI
        DO J=1, NXT-1
          JP1=J+1
          PU=Sue(J)*P(J,K)+Suw(JP1)*P(JP1,K)
         IF (PU.GT.PBACK) THEN
          HU(J,K)=HUO(J,K)+AC*(PU-PBACK)
          IF (HU(J,K).LE.0.0D0) Ilarge=1
         ELSE
          HU(J,K)=HUO(J,K)      
         END IF
        END DO
        IF (IFULL.EQ.1) HU(NXT,K)=HU(1,K)
      END DO

c    !............................! AT V control volumes
      DO K = 2, NYI
        KM1=K-1
        DO J=1, NXT
          PV=Svs(K)*P(J,K)+Svn(KM1)*P(J,KM1)
         IF (PV.GT.PBACK) THEN
          HV(J,K)=HVO(J,K)+AC*(PV-PBACK)
          IF (HV(J,K).LE.0.0D0) Ilarge=1
         ELSE
          HV(J,K)=HVO(J,K)         
         END IF
        END DO
      END DO

c    !............................!.............................!
c     AT exit side of bearing, we need to provide values of film
c     thickness due to pressure deformation,
c     Assume pad leading edge is not fixed (i.e. it does deform)
c            pad trailing edge deforms with upstream value.
c    !............................! at bearing exit plane,  Y=Lright
c##   KU=NYI-1                    ! foil deformation should occur
      KEND=NYI                    ! but since p=0 here, we take
      HP(1,KEND)=HV(1,KEND)       ! the first value of film thickness
      DO J=1, NXT-1               ! with deformation
        JP1=J+1                   !
        HP(JP1,KEND)=HV(JP1,KEND) ! 
        HU(J,KEND)=SUE(J)*HP(J,KEND)+SUW(J)*HP(JP1,KEND)
      END DO                      !
        HU(NXT,KEND)=HP(NXT,KEND) !
c    !............................!.............................!

c    !...........................! FOR COMPLIANT SEAL or PAD BEARING
      IF ((BEARING.EQ.2).OR.     ! at bearing centerline Y=0    
     +    ((BEARING.EQ.3).AND.(ISYM.EQ.1))) THEN
         DO J=1, NXT             ! ...........................
           HV(J,1)=HP(J,1)       ! Yv = Yp at centerline 
         END DO                  ! first V control volume
         GOTO 99
      END IF

c    !..........................! FOR ASYMMETRIC BEARING
      IF (ISYM.EQ.0) THEN 
c    !..........................! find film thickness at left side Y<0
      DO K=-2, -NYI, -1
        KM1=K+1
        DO J=1, NXT
          PV=Svs(K)*P(J,K)+Svn(K)*P(J,KM1)          
         IF (PV.GT.PBACK) THEN
          HV(J,K)=HVO(J,K)+AC*(PV-PBACK)
          IF (HV(J,K).LE.0.0D0) Ilarge=1
         ELSE
          HV(J,K)=HVO(J,K)            
         END IF
        END DO
      END DO

        DO J=1, NXT               ! at centerline Y=0
          HV(J,-1) = HP(J,-1)     ! ### only if no misaligment
          HV(J, 1) = HP(J, 1)
        END DO

      KEND=-NYI                   ! exit plane deformation
      HP(1,KEND)=HV(1,KEND)       ! 
      DO J=1, NXT-1               ! at Y=-LengthL
        JP1=J+1                   !
        HP(JP1,KEND)=HV(JP1,KEND) ! 
        HU(J,KEND)=SUE(J)*HP(J,KEND)+SUW(J)*HP(JP1,KEND)
      END DO                      !
        HU(NXT,KEND)=HP(NXT,KEND) !
c    !..........................!
      ELSE
c    !..........................! at symmetry line
         DO J=1, NXT            ! IF NOT Misaligned ONLY !###
           HV(J,1)=HP(J,1)
         END DO
c    !..........................!
      END IF
c    !..........................!

  99    CONTINUE
c    !............................!.............................!
c     AT PAD Trailing edge assume deformation equals 
c     to upstream value
c    !...........................! 9/30/93 only IF NO MISALIGNMENT
      IF (IFULL.EQ.0) THEN
c    !...........................!
      Kstart=1
      IF (ISYM.EQ.0) Kstart=-NYI   
      IF (BEARING.EQ.2) Kstart=1
      Kend=NYI
      JM1=NXT-1
      DO K=Kstart, Kend
        HP(NXT,K)=HU(JM1,K)
        HU(NXT,K)=HP(NXT,K)
      END DO
      KM1=1
      J=NXT
      DO K=1, Kend
        HV(NXT,K)=HP(J,K)*SVS(K)+HP(J,KM1)*SVN(K)
        KM1=K
      END DO
      IF (ISYM.EQ.0) THEN 
        KM1=-1
        DO K=-1, -KEND, -1
          HV(NXT,K)=HP(J,K)*SVS(K)+HP(J,KM1)*SVN(K)
          KM1=K
        END DO
        HV(NXT,0)=0.5D0*(HV(NXT,1)+HV(NXT,-1))
      END IF   
c    !...........................!
      END IF
c    !...........................! PAD compliant bearing 9/30/93
c    !............................!.............................!
C....................................................................      
C       CHECK if film thickness <=0 for compliant surface bearing
c      !.......................................................!
        IF (ILARGE.EQ.1) THEN
        WRITE (6,100)
 100    FORMAT('$ Journal travel too large and surface',
     +         ' deformation not enough,',/,
     +         '$ Film thickness: H is negative or zero',/,
     +         '$ PROGRAM ABORTS execution on SUB FILMC')
        CALL PAUSE
        END IF
c....................................................................

      END



C *****************************************************************************
C **                                                                         **
C **  Subroutine Inputcompl                                                  **
C **                                                                         **
C **  INPUTCOMPL: keyboard input of bearing compliance parameters            **
C **                                                                         **
C *****************************************************************************
 
      SUBROUTINE INPUTCOMPL

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /COMPLIA/ AC, ETA, PBACK, LIFT
      COMMON /BTYPE/ BEARING

      DOUBLE PRECISION AC, ETA, PBACK
      INTEGER BEARING, LIFT
C ----------------------------------------------------------------------------
C --  INPUTCOMPL code                                                       --
C ----------------------------------------------------------------------------
 906  CONTINUE

      write (6, *) '-----------------------------------------------'
      write (6, *) '(6)INPUT:COMPLIANCE COEFFICIENTS (dimensionless)'
      write (6, *) '-----------------------------------------------'
      write (6, *) 'Local deformation u [m] = (P-Pback) / Ke '
      write (6, *) 'A local compliance coeffient is needed as'
      write (6, *) 'ac = PSA/ (Ke x C* ) ,              where'
      write (6, *) '     PSA is the characteristic pressure,'
      write (6, *) '     Ke = local stiffness / unit area [N/m3],'
      write (6, *) 'C* = TYP Clearance, Pback: foil back pressure'
      write (6, *) '...............................................'
      CALL BEEPER

      WRITE (6, *) 'DIM Compliance Coefficient ac'
      WRITE (6, 170) AC

  170 FORMAT ('$', 'ENTER ac [Default (ac)= ',F8.4,']: ')

      CALL ENTERVAL(AC)
      AC=DABS(AC)
C
      WRITE (6, *) 'DIM LOSS FACTOR ETA'
      WRITE (6, 172) ETA

  172 FORMAT ('$', 'ENTER ETA [Default (eta)= ', F8.4,']: ')

      CALL ENTERVAL(ETA)
      ETA=DABS(ETA)

      WRITE (6,173)
  173 FORMAT ('$', 40('.'),/,
     +        '$ FOIL Lift option and under-relaxation not ready yet',
     +       /,'$', 40('.'))

      PBACK=0.0D0
      LIFT=0

      END
     
C ----------------------------------------------------------------------------
C --  foilsubs.f                                                          --
C ----------------------------------------------------------------------------