C
      PROGRAM MAIN
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C .. INPUT VARIABLES
C
C    P      PRESSURE                             IN PA  
C    TL     TEMPERATURE OF THE LIQUID            IN K
C
C .. OUTPUT VARIABLES
C
C    ROL    DENSITY OF THE LIQUID                IN KG/(M**3)
C    CPL    SPECIFIC HEAT OF THE LIQUID          IN J/(KG K)
C    VISL   DYNAMIC VISCOSITY OF THE LIQUID      IN N S/(M**2)
C    CL     THERMAL CONDUCUTIVITY OF THE LIQUID  IN W/(M K)
C

      INCLUDE 'IOUNITS'
C
      OPEN(UNIT=IMOUT,FILE='../output/steam.dat')
C
      P  = 1.550D7
      TL = 559.150
C
      CALL PF1STEAM( P, TL, ROL, CPL, VISL, CL )
C
      WRITE(6,1000)     P,TL,ROL,CPL,VISL,CL
      WRITE(IMOUT,1000) P,TL,ROL,CPL,VISL,CL
C
      CLOSE( UNIT=IMOUT)
C
 1000 FORMAT(1H ,1X,'P   = ',1P,E12.5,2X,'PA      ',/,
     1           2X,'TL  = ',1P,E12.5,2X,'K       ',/,
     2           2X,'ROL = ',1P,E12.5,2X,'KG/M3   ',/,
     3           2X,'CPL = ',1P,E12.5,2X,'J/(KG K)',/,
     4           2X,'VISL= ',1P,E12.5,2X,'N S/M2  ',/,
     5           2X,'CL  = ',1P,E12.5,2X,'W/(M K) ')
C
      STOP
      END
