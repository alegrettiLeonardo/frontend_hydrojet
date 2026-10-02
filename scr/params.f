C HYDROFLEXT program
c
C params.f: Maximum dimensions for arrays
C
C Maximum NUmber of pads on bearing surface

      INTEGER MAXNPAD
      PARAMETER (MAXNPAD=20)

C Maximum Number of pockets or recesses/pad

      INTEGER MAXNPOCK
      PARAMETER (MAXNPOCK=8)
C     
C Maximum number of grid poinst in circumferential direction/pad is
C MAXNXT: (Npc+Nlc-2)*Npockets+Nlc for 1pad of length L<2*PI*D
c       or(Npc+Nlc-2)*Npockets+1   for 1pad/2PID length    , Npc must be odd

      INTEGER MAXNXT
      PARAMETER (MAXNXT=85)
C                                     
C Maximum number of grid points in axial direction is
C MAXNYI: Npa+Nla+1 

      INTEGER MAXNYI
      PARAMETER (MAXNYI=31)
C
C
C maximum number of iterations/steps for evaluation of fluid 
c properties in subroutine SOLVE

      INTEGER CHPROPMAX
      PARAMETER (CHPROPMAX=2)  
c CHPROPMAX =1 or 2 for highly compressible fluids like LH2 or AIR
c           =5 or 10 for almost incompressible fluids like water, LOx, LN2, etc.
c
c program speed of execution depends on this number.
c
C
c maximum number of interpolation points for
c clearance evaluation

      INTEGER NSL
      PARAMETER (NSL=10)

c........
      INTEGER MAXNPOCKP1
      PARAMETER (MAXNPOCKP1=MAXNPOCK+1)
      INTEGER MAXNXTP2
      PARAMETER (MAXNXTP2=MAXNXT+2)
C
C End of params.f 12/19/94 by Dr. Luis San Andres.
c