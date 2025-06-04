C*DK SIGMA                                                              
      FUNCTION SIGMA( P, TSAT )                                         
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C     SURFACE TENSION OF WATER.                                         
C     E. SCHMIDT.  PROPERTIES OF WATER AND STEAM IN SI UNITS -- 1979.   
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C                                                                       
C*CA FUNCTION                                                           
C*CC FUNCTION                                                           
      INCLUDE 'FUNCTION'                                                
C                                                                       
      TT    = DMAX1((647.15D0-TSAT)/647.15D0,0.1D0)                     
      SIGMA = 0.2358D0*STAR(TT,1.256D0)*(1.0D0-0.625D0*TT)              
      RETURN                                                            
      END                                                               
crc@PROCESS XOPT(NOAMOVE) OPT(2)                                        
