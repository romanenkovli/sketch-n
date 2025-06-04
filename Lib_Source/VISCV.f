C*DK VISCV                                                              
      FUNCTION VISCV( H, P, ROV, TV, PA )                               
      IMPLICIT REAL*8(A-H,O-Z)                                          
C*CA CONTRLLR                                                           
C*CC CONTRLLR                                                           
      INCLUDE 'CONTRLLR'                                                
C                                                                       
C                                                                       
C   DYNAMIC VISCOSITY OF STEAM AS FUNCTION OF PRESSURE AND ENTHALPY     
C        PRESSURE  P       IN (N/SQ.M)                                  
C        ENTHALPY  H       IN (J/KG)                                    
C         DENSITY   ROV      IN (KG/CU.M)                               
C         TEMPERATURE T      IN (K)                                     
C        DYNAMIC VISCOSITY IN (N S/SQ.M.)                               
C              REQUIRES POLY                                            
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
      DIMENSION A(3) ,F(4) ,G(4)                                        
      DATA A/3.53D-8, 6.765D-11, 1.021D-14/, B/.407D-7/, C/8.04D-6/,    
     *     D/1.858D-7/, E/5.9D-10/,                                     
     *     F/-.2885D-5, .2427D-7, -.67893333333D-10, .6317037037D-13/,  
     *     G/.176D3, -.16D1, .48D-2, -.47407407407D-5/                  
C                                                                       
      POLY2(A1,A2,A3,X) = A1+X*(A2+X*A3)                                
      POLY3(A1,A2,A3,A4,X) = A1+X*(A2+X*(A3+X*A4))                      
C  POLY4 IS VISCOSITY FUNCTION FOR HYDROGEN;2ND ORDER POLYNOMIAL FIT TO 
C  REFERENCE: HANDBOOK OF THERMODYNAMIC TABLES AND CHARTS,              
C             KUZMAN RAZNJEVIC, HEMISPHERE PUBLISHING CORP.             
      POLY4(TH) = 4.175D-06+1.588D-08*TH+7.6705D-13*TH*TH               
      IF( IEOS.NE.0.AND.IGAS.EQ.1 ) GO TO 50                            
      IF( IEOS.NE.0 ) GO TO 65                                          
      T  = TV-273.15D0                                                  
      RO = ROV                                                          
      V1 = B*T+C                                                        
      IF( T.LE.300.D0 ) GO TO 10                                        
      IF( T.GE.375.D0 ) GO TO 20                                        
      VISCV = V1+(POLY3(F(1),F(2),F(3),F(4),T)+POLY3(G(1),G(2),G(3),    
     1        G(4),T)*POLY2(A(1),A(2),A(3),RO))*RO                      
C                                                                       
      GO TO 70                                                          
10    VISCV = V1-RO*(D-E*T)                                             
      IF( VISCV.LT.1.D-7 ) VISCV = 1.D-7                                
C                                                                       
      GO TO 70                                                          
                                                                        
   20 VISCV = V1+POLY2(A(1),A(2),A(3),RO)*RO                            
C                                                                       
      GO TO 70                                                          
  50  T = TV-273.15D0                                                   
      IF( T.GT.229.D0 ) GO TO 60                                        
      VISCV = 17.08D-6+T*(5.927D-8-8.14D-11*T)                          
      RETURN                                                            
  60  VISCV = 17.35D-6+T*(4.193D-8-1.09D-11*T)                          
      RETURN                                                            
   65 VISCV = POLY4(TV)                                                 
      RETURN                                                            
  70  IF( PA.LT.0.1D0 ) RETURN                                          
      IF( IGAS.EQ.1 ) GO TO 80                                          
      VISCA = POLY4(TV)                                                 
      GO TO 100                                                         
  80  T = TV-273.15D0                                                   
      IF( T.GT.229.D0 ) GO TO 90                                        
      VISCA = 17.08D-6+T*(5.927D-8-8.14D-11*T)                          
      GO TO 100                                                         
  90  VISCA = 17.35D-6+T*(4.193D-8-1.09D-11*T)                          
 100  VISCV = ((P-PA)*VISCV+PA*VISCA)/P                                 
      RETURN                                                            
      END                                                               
