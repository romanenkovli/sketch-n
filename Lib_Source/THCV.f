C*DK THCV                                                               
      FUNCTION THCV( H, P, ROV, TV, PA )                                
      IMPLICIT REAL*8(A-H,O-Z)                                          
C*CA CONTRLLR                                                           
C*CC CONTRLLR                                                           
      INCLUDE 'CONTRLLR'                                                
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
      EQUIVALENCE (THAX,CEOSLP(13))                                     
C                                                                       
C   THERMAL CONDUCTIVITY OF STEAM AS FUNCTION OF PRESSURE AND ENTHALPY  
C        PRESSURE  P          IN (N/SQ.M)                               
C        ENTHALPY  H          IN (J/KG)                                 
C        DENSITY      ROV        IN (KG/CU.M)                           
C        TEMPERATURE  TV         IN (K)                                 
C        THERMAL CONDUCTIVITY IN (W/M K)                                
C              REQUIRES POLY                                            
C                                                                       
      DIMENSION A(4) ,B(3)                                              
      DATA A/1.76D-2, 5.87D-5, 1.04D-7, -4.51D-11/, B/1.0351D-4,        
     *     .4198D-6, -2.771D-11/, C/2.1482D5/                           
      POLY2(A1,A2,A3,X) = A1+X*(A2+X*A3)                                
      POLY3(A1,A2,A3,A4,X) = A1+X*(A2+X*(A3+X*A4))                      
C                                                                       
      IF( IEOS.NE.0 ) GO TO 5                                           
      TDIFF = DMAX1(0.1D0,TV-CEOSLP(34))                                
      THCV  = DMAX1(0.0001D0,POLY3(A(1),A(2),A(3),A(4),TDIFF)+ROV*      
     1        (POLY2(B(1),B(2),B(3),TDIFF)+ROV*C*EXP(-4.2D0*            
     2        DLOG(TDIFF))))                                            
      IF( PA.LT.1.D0 ) GO TO 10                                         
      THCV = ((P-PA)*THCV+PA*THAX)/P                                    
   10 CONTINUE                                                          
      RETURN                                                            
    5 THCV = THAX                                                       
      RETURN                                                            
      END                                                               
crc@PROCESS XOPT(NOAMOVE) OPT(2)                                        
