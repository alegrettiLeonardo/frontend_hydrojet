C
C  ####   #    #  #####   #####    ####   #####    #####  #####           ######
C #       #    #  #    #  #    #  #    #  #    #     #    #    #          #
C  ####   #    #  #    #  #    #  #    #  #    #     #    #    #          #####
C      #  #    #  #####   #####   #    #  #####      #    #####    ###    #
C #    #  #    #  #       #       #    #  #   #      #    #        ###    #
C  ####    ####   #       #        ####   #    #     #    #        ###    #
C
C
C SUPPORT programs for I/O (Input/Output) operation with hydroxxxx codes.
C FROM Library of programs of Tribology Group
C 
c hydrosealt.f  Copyright Luis SanAndres and Zhou Yang/TexasA&MUniversity/1993
c
C *****************************************************************************
C **                                                                         **
C **  Subroutine Beeper                                                      **
C **                                                                         **
C **  BEEPER:  Send a beep to the terminal if the beep flag is set.          **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE BEEPER

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------

      COMMON /VERB/ SVERB, DVERB, BEEP

      INTEGER SVERB, DVERB, BEEP

C ----------------------------------------------------------------------------
C --  BEEPER code                                                           --
C ----------------------------------------------------------------------------

      IF (BEEP.eq.1) THEN
          WRITE (6, 1) 7
      END IF

    1 FORMAT ('$', A)

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Zero                                                        **
C **                                                                         **
C **  ZERO:  Zero array z of length jmax.                                    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ZERO(Z, JMAX)

      IMPLICIT NONE


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER JMAX
      DOUBLE PRECISION Z(JMAX)
      INTEGER I

C ----------------------------------------------------------------------------
C --  ZERO code                                                             --
C ----------------------------------------------------------------------------

      DO I=1, JMAX
          Z(I)=0.0D0
      END DO

      END


C *****************************************************************************
C **                                                                         **
C **  Function Maxin                                                         **
C **                                                                         **
C **  MAXIN:  Find the maximum element in a two dimensional array.           **
C **                                                                         **
C *****************************************************************************

      DOUBLE PRECISION FUNCTION MAXIN(ARRAY)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING

      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL, BEARING 

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION EXO,EYO, ECC,AXO, AYO, ZO, 
     +                 ARRAY(MAXNXT,-MAXNYI:MAXNYI)                 
      INTEGER I, J, NN, K, ISTART

C ----------------------------------------------------------------------------
C --  MAXIN code                                                            --
C ----------------------------------------------------------------------------
      ECC=DSQRT(EXO*EXO+EYO*EYO)
      Istart=ISYM+(ISYM-1)*NYI
      NN=NXT

      IF (BEARING.EQ.2) THEN
          Istart=1
      ELSE IF (BEARING.EQ.1) THEN
        IF (IFULL.EQ.0) THEN
         NN=NXT
         GOTO 17
        END IF

        IF (ECC.GT.(0.01D0)) THEN
          NN=NXT
        ELSE
         NN=NXI
         IF ((AXO.NE.0.0D0).AND.(AYO.NE.0.0D0)) THEN
             NN=NXT
         END IF 
        END IF
      END IF


 17   MAXIN=-1.00D+30

      DO I=1, NN
          DO J=ISTART, NYI, 1
              MAXIN=DMAX1(MAXIN, ARRAY(I,J))
          END DO 
      END DO

      END

C *****************************************************************************
C **                                                                         **
C **  Function Minin                                                         **
C **                                                                         **
C **  MININ:  Find the minimum element in a two dimensional array.           **
C **                                                                         **
C *****************************************************************************

      DOUBLE PRECISION FUNCTION MININ(ARRAY)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  COMMON variable declarations                                          --
C ----------------------------------------------------------------------------
      COMMON /PARAM2/ EXO, EYO
      COMMON /ALIGNM/ AXO, AYO, ZO
      COMMON /NODES/  NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,IFULL
      COMMON /HJBsym/ ISYM, ICSTEP
      COMMON /BTYPE/ BEARING
      INTEGER NPOCKET,NLC,NPC,NLA,NPA,NPAP1,NXI,NYI,NXT,ISYM,ICSTEP,
     +        IFULL 

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      DOUBLE PRECISION EXO,EYO, ECC,AXO, AYO, ZO, 
     +                 ARRAY(MAXNXT,-MAXNYI:MAXNYI)               
      INTEGER I, J, IStart, NN, K, BEARING

C ----------------------------------------------------------------------------
C --  MININ code                                                            --
C ----------------------------------------------------------------------------
      ECC=DSQRT(EXO*EXO+EYO*EYO)
      Istart=ISYM+(ISYM-1)*NYI
      NN=NXT
      IF (BEARING.EQ.2) THEN
         Istart=1
      ELSE IF (BEARING.EQ.1) THEN
        IF (IFULL.EQ.0) THEN
         NN=NXT
         GOTO 17
        END IF

        IF (ECC.GT.(0.01D0)) THEN
         NN=NXT
        ELSE
         NN=NXI
         IF ((AXO.NE.0.0D0).AND.(AYO.NE.0.0D0)) THEN
             NN=NXT
        END IF 
       END IF
      END IF

 17   MININ=1.0D+30
      DO I=1, NN
          DO J= Istart, NYI, 1
              MININ=DMIN1(MININ, ARRAY(I,J))
          END DO
      END DO

      END

C *****************************************************************************
C **                                                                         **
C **  Subroutine Echelon                                                     **
C **                                                                         **
C **  ECHELON:  Reduce an augmented matrix to reduced-row echelon form.      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE ECHELON(A, B, C, AR1, AR2)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER A, B, C
      DOUBLE PRECISION AR1(MAXNPOCK, MAXNPOCK), AR2(MAXNPOCK, MAXNPOCK)
      DOUBLE PRECISION T
      INTEGER I, J, K, L, CJ

C ----------------------------------------------------------------------------
C --  ECHELON code                                                          --
C ----------------------------------------------------------------------------

      CJ=1
      DO I=1, B
          J=CJ
          DO WHILE ((J.LT.A).AND.(AR1(J, I).EQ.(0.0)))
              J=J+1
          END DO
          IF (J.LE.A) THEN
              IF (J.NE.CJ) THEN
                  DO K=1, B
                      T=AR1(J, K)
                      AR1(J, K)=AR1(CJ, K)
                      AR1(CJ, K)=T
                      T=AR2(J, K)
                      AR2(J, K)=AR2(CJ, K)
                      AR2(CJ, K)=T
                  END DO
                  J=CJ
              END IF
              CJ=CJ+1
              T=1/AR1(J, I)
              DO K=1, B
                  AR1(J, K)=AR1(J, K)*T
                  AR2(J, K)=AR2(J, K)*T
              END DO
              DO K=1, A
                  IF ((K.NE.J).AND.(AR1(K, I).NE.(0.0))) THEN
                      T=-AR1(K, I)
                      DO L=1, B
                          AR1(K, L)=AR1(K, L)+T*AR1(J, L)
                          AR2(K, L)=AR2(K, L)+T*AR2(J, L)
                      END DO
                  END IF
              END DO
          END IF
      END DO

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Inverse                                                     **
C **                                                                         **
C **  INVERSE:  Invert a matrix by the echelon trick.                        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE INVERSE(A, AR1, INV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER A
      DOUBLE PRECISION AR1(MAXNPOCK, MAXNPOCK), INV(MAXNPOCK, MAXNPOCK)
      DOUBLE PRECISION SCRATCH(MAXNPOCK, MAXNPOCK)
      INTEGER I, J

C ----------------------------------------------------------------------------
C --  INVERSE code                                                          --
C ----------------------------------------------------------------------------

      DO I=1, A
          DO J=1, A
              IF (I.EQ.J) THEN
                  INV(I, J)=1.0D+00
              ELSE
                  INV(I, J)=0.0D+00
              END IF
              SCRATCH(I, J)=AR1(I, J)
          END DO
      END DO
      CALL ECHELON(A, A, A, SCRATCH, INV)
      
      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Help                                                        **
C **                                                                         **
C **  HELP:  Provide access to file HJBHELP.DAT during the execution of     **
C **         HYDROSEAL.                                                      **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE HELP

      IMPLICIT NONE

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      CHARACTER*79 LINE, TOPIC, SCR1

C ----------------------------------------------------------------------------
C --  HELP code                                                             --
C ----------------------------------------------------------------------------

      OPEN (UNIT=10, FILE='HJBHELP.TXT', STATUS='OLD', ERR=10)
      
  100 FORMAT (A80)
  
  200 READ (10, 100) SCR1
      IF (SCR1.NE.'OPENING') THEN
          GOTO 200
      END IF

  300 READ (10, 100) LINE
      IF (LINE.NE.'CLOSE') THEN
          WRITE (6, 400) LINE
          GOTO 300
      END IF
  400     FORMAT (' ', A79)

      WRITE (6, *) 'Topics:'

  500 READ (10, 100) SCR1
      IF (SCR1(1:6).EQ.'TOPIC') THEN
          WRITE (6, 400) SCR1(7:79)
          GOTO 500
      ELSE IF (SCR1(1:7).NE.'ENDHELP') THEN
          GOTO 500
      END IF
      WRITE (6, *) ' '
      REWIND (UNIT=10)
  600 FORMAT ('$', 'Topic, or <ENTER> to exit: ')
  900 WRITE (6, 600)
      READ (5, 100) SCR1
      IF (SCR1.EQ.' ') THEN
          WRITE (6, *) 'Returning from help.'
          RETURN
      END IF
      TOPIC='TOPIC '//SCR1
  700 READ (10, 100) LINE
      IF (LINE.EQ.TOPIC) THEN
          WRITE (6, *) ' '
  800     READ (10, 100) LINE
          IF (LINE.EQ.'PAUSE') THEN
 1000         FORMAT ('$', 'Press <Return> to continue.')
              WRITE (6, 1000)
              READ (5, 100) SCR1
              GOTO 800
          ELSE IF (LINE.NE.'ENDTOPIC') THEN
              WRITE (6, *) LINE
              GOTO 800
          END IF
          WRITE (6, *) ' '
          WRITE (6, *) 'Press <ENTER> to return to MAIN HELP menu'
          READ (5, 100) SCR1
          REWIND (UNIT=10)
          GOTO 200
      ELSE IF (LINE.EQ.'ENDHELP') THEN
          WRITE (6, *) 'Topic not found.  Please enter the topic name ex
     +actly as it appears, '
          WRITE (6, *) 'starting  with a blank & using capital letters'
          REWIND (UNIT=10)
          GOTO 900
      ENDIF
      GOTO 700

   10 WRITE (6, *) 'File HJBHELP.DAT could not be found.'
      RETURN
                         
      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Pause                                                       **
C **                                                                         **
C **  PAUSE:  Prompt the user and wait for a <Return> key.                   **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE PAUSE

      IMPLICIT NONE


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      CHARACTER*1 X

C ----------------------------------------------------------------------------
C --  PAUSE code                                                            --
C ----------------------------------------------------------------------------

      WRITE (6, 10)
   10 FORMAT ('$', 'Press <Return> to continue.')
      READ (5, 20) X
   20 FORMAT (1A)

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Decodios                                                    **
C **                                                                         **
C **  DECODIOS:  Decode the VAX's IOSTAT regester and report the error       **
C **             to the user.                                                **
C **                                                                         **
C *****************************************************************************
                             
      SUBROUTINE DECODIOS(I)

      IMPLICIT NONE


C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER I

C ----------------------------------------------------------------------------
C --  DECODIOS code                                                         --
C ----------------------------------------------------------------------------

      WRITE (6, 10) I
   10 FORMAT (' ', 'Fortran I/O Error: IOSTAT=', I4)
      IF (I.EQ.-1) THEN
          WRITE (6, *) 'End of file.'
      ELSE IF (I.EQ.0) THEN
          WRITE (6, *) 'No error:  What are we doing in this routine?'
      ELSE IF (I.EQ.21) THEN
          WRITE (6, *) 'Vax error FOR-F-DUPFILSPE'
          WRITE (6, *) 'Duplicate file specification.'
      ELSE IF (I.EQ.22) THEN
          WRITE (6, *) 'Vax error FOR-F-INPRECTOO'
          WRITE (6, *) 'Input record too long.'
      ELSE IF (I.EQ.24) THEN
          WRITE (6, *) 'Vax error FOR-F-ENDDURREA'
          WRITE (6, *) 'End of file during read.'
      ELSE IF (I.EQ.29) THEN
          WRITE (6, *) 'Vax error FOR-F-FILNOTFOU'
          WRITE (6, *) 'File not found.'
      ELSE IF (I.EQ.30) THEN
          WRITE (6, *) 'Vax error FOR-F-OPEFAI'
          WRITE (6, *) 'Open failure.'
      ELSE IF (I.EQ.34) THEN
          WRITE (6, *) 'Vax error FOR-F-UNIALROPE'
          WRITE (6, *) 'Unit already open.'
      ELSE IF (I.EQ.38) THEN
          WRITE (6, *) 'Vax error FOR-F-ERRDURWRI'
          WRITE (6, *) 'Error during write.'
      ELSE IF (I.EQ.39) THEN
          WRITE (6, *) 'Vax error FOR-F-ERRDURREA'
          WRITE (6, *) 'Error during read.'
      ELSE IF (I.EQ.42) THEN
          WRITE (6, *) 'Vax error FOR-F-NO_SUCDEV'
          WRITE (6, *) 'No such device.'
      ELSE IF (I.EQ.43) THEN
          WRITE (6, *) 'Vax error FOR-F-FILNAMSPE'
          WRITE (6, *) 'File name specification error.'
      ELSE IF (I.EQ.63) THEN
          WRITE (6, *) 'Vax error FOR-F-OUTCONERR'
          WRITE (6, *) 'Output conversion error.'
      ELSE IF (I.EQ.64) THEN
          WRITE (6, *) 'Vax error FOR-F-INPCONERR'
          WRITE (6, *) 'Input conversion error.'
      ELSE
          WRITE (6, *) 'Error is not in my dictionary.'
      END IF

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Cechelon                                                    **
C **                                                                         **
C **  CECHELON:  Reduce a complex matrix to reduced-row echelon form.        **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CECHELON(A, B, C, AR1, AR2)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER A, B, C
      DOUBLE COMPLEX AR1(MAXNPOCK, MAXNPOCK), AR2(MAXNPOCK, MAXNPOCK)
      DOUBLE COMPLEX T
      INTEGER I, J, K, L, CJ

C ----------------------------------------------------------------------------
C --  CECHELON code                                                         --
C ----------------------------------------------------------------------------

      CJ=1
      DO I=1, B
          J=CJ
          DO WHILE ((J.LT.A).AND.(AR1(J, I).EQ.(0.0, 0.0)))
              J=J+1
          END DO
          IF (J.LE.A) THEN
              IF (J.NE.CJ) THEN
                  DO K=1, B
                      T=AR1(J, K)
                      AR1(J, K)=AR1(CJ, K)
                      AR1(CJ, K)=T
                      T=AR2(J, K)
                      AR2(J, K)=AR2(CJ, K)
                      AR2(CJ, K)=T
                  END DO
                  J=CJ
              END IF
              CJ=CJ+1
              T=(1.0D+00, 0.0D+00)/AR1(J, I)
              DO K=1, B
                  AR1(J, K)=AR1(J, K)*T
                  AR2(J, K)=AR2(J, K)*T
              END DO
              DO K=1, A
                  IF ((K.NE.J).AND.(AR1(K, I).NE.(0.0))) THEN
                      T=-AR1(K, I)
                      DO L=1, B
                          AR1(K, L)=AR1(K, L)+T*AR1(J, L)
                          AR2(K, L)=AR2(K, L)+T*AR2(J, L)
                      END DO
                  END IF
              END DO
          END IF
      END DO

      END 

C *****************************************************************************
C **                                                                         **
C **  Subroutine Cinverse                                                    **
C **                                                                         **
C **  CINVERSE:  Invert a complex matrix.                                    **
C **                                                                         **
C *****************************************************************************

      SUBROUTINE CINVERSE(A, AR1, INV)

      IMPLICIT NONE

      INCLUDE 'params.f'

C ----------------------------------------------------------------------------
C --  Local variable declarations                                           --
C ----------------------------------------------------------------------------

      INTEGER A
      DOUBLE COMPLEX AR1(MAXNPOCK, MAXNPOCK), INV(MAXNPOCK, MAXNPOCK)
      DOUBLE COMPLEX SCRATCH(MAXNPOCK, MAXNPOCK)
      INTEGER I, J

C ----------------------------------------------------------------------------
C --  CINVERSE code                                                         --
C ----------------------------------------------------------------------------
            
      DO I=1, A
          DO J=1, A
              IF (I.EQ.J) THEN
                  INV(I, J)=(1.0D0, 0.0D0)
              ELSE
                  INV(I, J)=(0.0D0, 0.0D0)
              END IF
              SCRATCH(I, J)=AR1(I, J)
          END DO
      END DO
      CALL CECHELON(A, A, A, SCRATCH, INV)
      
      END 
c
c----------------------------------------------------------------------------
c last revised 12/31/93 by Dr. Luis San Andres at Texas A&M University
c----------------------------------------------------------------------------
