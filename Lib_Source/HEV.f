C*DK HEV                                                                
      FUNCTION HEV( TEMP )                                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     THIS FUNCTION CALCULATES THE HEAT OF EVAPORATION OF LIQUID        
C     CORRESPONDING TO A GIVEN TEMPERATURE FOR LOW PRESSURES            
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C                                                                       
      A   = 3180619.59D0                                                
      B   = 2470.2120D0                                                 
      T   = DMAX1(TEMP,CEOSLP(32))                                      
      HEV = A-B*T                                                       
      RETURN                                                            
      END                                                               
