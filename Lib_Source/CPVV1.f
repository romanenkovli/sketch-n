C*DK CPVV1                                                              
      FUNCTION CPVV1( T, P, PA )                                        
      IMPLICIT REAL*8(A-H,O-Z)                                          
C*CA CONTRLLR                                                           
C*CC CONTRLLR                                                           
      INCLUDE 'CONTRLLR'                                                
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
      EQUIVALENCE (CPAX,CEOSLP(22))                                     
C        SPECIFIC HEAT OF WATER VAPOR AS A FUNTION OF TEMPERATURE       
C        AND PRESSURE.                                                  
C                                                                       
C          TEMPERATURE      T     IN      (K)                           
C          PRESSURE         P     IN      (N/M**2)                      
C          SPECIFIC HEAT    CPV   IN      (J/KG/K)                      
C                                                                       
C                                                                       
      DATA C1,C2,C3,C4,C5,C6/1688.35968D0,0.6029856D0,                  
     *     482.0979623D0,2.95317905D+7,1.8D0,460.D0/                    
C*CA FUNCTION                                                           
C*CC FUNCTION                                                           
      INCLUDE 'FUNCTION'                                                
      IF( IEOS.NE.0 ) GO TO 1                                           
      TB    = T*C5-C6                                                   
      P3    = P*P*P                                                     
      CPVV1 = C1+C2*T+C3*P/STAR(TB,2.4D0)+C4*P3/TB**9                   
      IF( PA.GE.0.1D0 ) GO TO 2                                         
      RETURN                                                            
    1 CPVV1 = CPAX                                                      
      RETURN                                                            
    2 CPVV1 = ((P-PA)*CPVV1+PA*CPAX)/P                                  
      RETURN                                                            
      END                                                               
