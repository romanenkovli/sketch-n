C*DK SATTMP                                                             
      FUNCTION SATTMP( PRES )                                           
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     THIS FUNCTION CALCULATES SATURATION TEMPERATURE OF VAPOR          
C     CORRESPONDING TO A GIVEN PRESSURE                                 
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
      IF( PRES.GE.CEOSLP(20) ) GO TO 10                                 
      P      = DMAX1(PRES,CEOSLP(36))                                   
      SATTMP = CEOSLP(1)/(CEOSLP(2)*DLOG(P/CEOSLP(11))+CEOSLP(3))       
      HFGREF = HEV(SATTMP)                                              
      PSREF  = SATPRS(SATTMP)                                           
      SATTMP = SATTMP/(1.0D0-CEOSLP(12)*SATTMP*DLOG(P/PSREF)/HFGREF)    
      SATTMP = DMAX1(SATTMP,CEOSLP(32))                                 
      RETURN                                                            
   10 CONTINUE                                                          
      IF( PRES.GE.CEOSLP(39) ) GO TO 20                                 
      SATTMP = CEOS1*STAR(AEOS14*PRES,CEOS2)+CEOS3                      
      RETURN                                                            
   20 CONTINUE                                                          
      IF( PRES.GE.CEOSLP(37) ) GO TO 30                                 
      PLOG   = DLOG(PRES)                                               
      SATTMP = (4264.8240952941D0+SQRT(-13666986.708428D0+              
     1         1166669.3278328D0*PLOG))/(27.304833093884D0-PLOG)        
      RETURN                                                            
   30 CONTINUE                                                          
      SATTMP = 4924.9229385171D0/(24.520401414546D0-DLOG(PRES))         
      RETURN                                                            
      END                                                               
