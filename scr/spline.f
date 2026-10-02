c
c  ####   #####   #          #    #    #  ######          ######
c #       #    #  #          #    ##   #  #               #
c  ####   #    #  #          #    # #  #  #####           #####
c      #  #####   #          #    #  # #  #        ###    #
c #    #  #       #          #    #   ##  #        ###    #
c  ####   #       ######     #    #    #  ######   ###    #
c
c #    #   #   #  #####   #####    ####   ######  #       ######  #    #
c #    #    # #   #    #  #    #  #    #  #       #       #        #  #
c ######     #    #    #  #    #  #    #  #####   #       #####     ##
c #    #     #    #    #  #####   #    #  #       #       #         ##
c #    #     #    #    #  #   #   #    #  #       #       #        #  #
c #    #     #    #####   #    #   ####   #       ######  ######  #    #

c...................................................
C spline.f > hydroflex.f / Dr. Luis San Andres, TexasA&MUniv. 1994
C
c NASA Grant NAG3-1434 "Thermohydrodynamic Analysis of Cryogenic Liquid
c                       Turbulent Flow Fluid Film Bearings" YEAR II
c Technical monitor: Mr. James Walker, NASA Lewis Research Center
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


      SUBROUTINE SPLINE(N,X,Y,B,C,D)

c...............................................................
C
C SPLINE: Evaluates spline coefficients for Y(X) function
C
C source code from COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS
C               BY G. FORSYTHE, M. MALCOLM. C. MOLER,
C                  PRENTICE-HALL , 1977
C
C...............................................................
C The coefficients B(i),C(i),D(i) are computed for a cubic inter-
c polating spline:
c s(X)= Y(i)+ B(i)*(X-Xi) + C(i)*(X-Xi)**2 + D(i)*(X-Xi)**3
c on interval
c                   Xi.LE.X.LE.Xi+1
c INPUT: N number of data points, N.GE.2
c        X the absissas of the knots in increasing order
c        Y the ordinates of the knots.
c.................................................................

      INCLUDE 'params.f'

      DOUBLE PRECISION T, X(Nsl),Y(Nsl),B(Nsl),C(Nsl),D(Nsl)
      INTEGER N, kst, NM1,IB,I
C
c      IF(KST.NE.0) GO TO 60      !Kst=1 -> linear interpolation
c.....................................................
      DO Kst=1, Nsl
         B(kst)=0.0D0
         C(Kst)=0.0D0
         D(Kst)=0.0D0
      END DO


      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 50
C
c ...... Set up tridiagonal system
c        B=diagonal, D=offdiagonal, C=RHS
c
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 10 I=2,NM1
        D(I)=X(I+1)-X(I)
        B(I)=2.*(D(I-1)+D(I))
        C(I+1)=(Y(I+1)-Y(I))/D(I)
        C(I)=C(I+1)-C(I)
   10 CONTINUE
C
c End conditions, third derivatives at X(1) & X(N)
c obtained by divided differences.
c
      B(1)=-D(1)
      B(N)=-D(N-1)
      C(1)=0.0
      C(N)=0.0
      IF(N.EQ.3) GO TO 15
      C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1)=C(1)*D(1)**2/(X(4)-X(1))
      C(N)=-C(N)*D(N-1)**2/(X(N)-X(N-3))
C
c ... forward elimination
c
   15 DO 20 I=2,N
        T=D(I-1)/B(I-1)
        B(I)=B(I)-T*D(I-1)
        C(I)=C(I)-T*C(I-1)
   20 CONTINUE
C
c .... back substitution
c
      C(N)=C(N)/B(N)
      DO 30 IB=1,NM1
        I=N-IB
        C(I)=(C(I)-D(I)*C(I+1))/B(I)
   30 CONTINUE
C
c .... C(i) -> Sigma
c
      B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.*C(N))
      DO 40 I=1,NM1
        B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.*C(I))
        D(I)=(C(I+1)-C(I))/D(I)
        C(I)=3.*C(I)
   40 CONTINUE
      C(N)=3.*C(N)
      D(N)=D(N-1)
      RETURN
c
C....................................................
c N. LT. 3 case
c
   50 B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0
      D(1)=0.0
      B(2)=B(1)
      C(2)=0.0
      D(2)=0.0
      RETURN
c....................................................
c Linear interpolation case
c
   60 CONTINUE

      DO 70 I=1,N
      B(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      C(I)=0.0
   70 D(I)=0.0
      B(N)=B(NM1)
      RETURN
c....................................................
c

      END

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::


      SUBROUTINE SEVAL(N,U,X,Y,B,C,D,F)

c........................................................
c
c SEVAL: evaluates the cubic spline function
c
c Seval=Y(I)+B(I)*(U-X(I)) +C(I)*(U-X(I))**2+D(I)*(U-X(I))**3
c
c.........................................................

      include 'params.f'

      double precision DX, U,F, X(Nsl),Y(Nsl),B(Nsl),C(Nsl),D(Nsl)
      integer N, i,j, k
c........................................................
      DATA I/1/
      IF(I.GE.N) I=1
      IF(U.LT.X(I)) GO TO 10
      IF(U.LE.X(I+1)) GO TO 30
C................................ Start Binary search
c
   10 I=1
      J=N+1
   20 K=(I+J)/2.
      IF(U.LT.X(K)) J=K
      IF(U.GE.X(K)) I=K
      IF(J.GT.I+1) GO TO 20
C
C............................... evaluate spline
   30 DX=U-X(I)
      F=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
c      FD=B(I)+(2.*C(I)+3.*D(I)*DX)*DX     ! first derivative
c      FDD=2.*C(I)+6.*D(I)*DX              ! second derivative
      RETURN
      END
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  LAST revised 09/14/93 LuisSanAndres TexasA&M University
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::