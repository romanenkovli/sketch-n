C*DK SATDER                                                             
      FUNCTION SATDER( PRES, TEMP )                                     
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     THIS FUNCTION CALCULATES THE DERIVATIVE OF SATURATION TEMPERATURE 
C     OF VAPOR WITH RESPECT TO PRESSURE                                 
C                                                                       
C     EQUATION OF STATE CONSTANTS ARRAY                                 
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C                                                                       
      IF( TEMP.GE.CEOSLP(21) ) GO TO 10                                 
      P      = DMAX1(PRES,CEOSLP(36))                                   
      T      = DMAX1(TEMP,CEOSLP(32))                                   
      SATDER = CEOSLP(12)*T*T/(P*HEV(T))                                
      RETURN                                                            
   10 CONTINUE                                                          
      IF( TEMP.GE.CEOSLP(40) ) GO TO 20                                 
      SATDER = CEOS2*(TEMP-CEOS3)/PRES                                  
      RETURN                                                            
   20 CONTINUE                                                          
      IF( TEMP.GE.CEOSLP(38) ) GO TO 30                                 
      SATDER = -TEMP**2/(PRES*(-8529.6481905883D0+2333338.6556656D0/    
     1         TEMP))                                                   
      RETURN                                                            
   30 CONTINUE                                                          
      SATDER = 2.0304886238506D-04*TEMP**2/PRES                         
      RETURN                                                            
      END                                                               
