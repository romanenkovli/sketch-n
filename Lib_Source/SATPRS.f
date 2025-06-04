C*DK SATPRS                                                             
      FUNCTION SATPRS( TEMP )                                           
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     THIS FUNCTION CALCULATES SATURATION PRESSURE OF VAPOR             
C     CORRESPONDING TO A GIVEN TEMPERATURE                              
C                                                                       
C     EQUATION OF STATE CONSTANTS ARRAY                                 
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C                                                                       
C*CA FUNCTION                                                           
C*CC FUNCTION                                                           
      INCLUDE 'FUNCTION'                                                
C                                                                       
      IF( TEMP.GE.CEOSLP(21) ) GO TO 10                                 
      T      = DMAX1(TEMP,CEOSLP(32))                                   
      SATPRS = 24821.0D0*STAR(T/338.0D0,-5.3512D0)*EXP(20.387D0*(T-     
     1         338.0D0)/T)                                              
      SATPRS = DMAX1(SATPRS,CEOSLP(36))                                 
      RETURN                                                            
   10 CONTINUE                                                          
      IF( TEMP.GE.CEOSLP(40) ) GO TO 20                                 
      SATPRS = (1.0D0/AEOS14)*STAR((TEMP-CEOS3)/CEOS1,1.0D0/CEOS2)      
      RETURN                                                            
   20 CONTINUE                                                          
      IF( TEMP.GE.CEOSLP(38) ) GO TO 30                                 
      SATPRS = 7.2166948490268D+11*EXP((-8529.6481905883D0+             
     1         1166669.3278328D0/TEMP)/TEMP)                            
      RETURN                                                            
   30 CONTINUE                                                          
      SATPRS = CEOSLP(37)*EXP(7.6084086799277D0-4924.9229385171D0/TEMP) 
      RETURN                                                            
      END                                                               
