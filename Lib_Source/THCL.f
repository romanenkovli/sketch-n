C*DK THCL                                                               
      FUNCTION THCL( H )                                                
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C                                                                       
C   THERMAL CONDUCTIVITY OF WATER AS FUNCTION OF PRESSURE AND ENTHALPY  
C        ENTHALPY  H          IN (J/KG)                                 
C        THERMAL CONDUCTIVITY IN (W/M K)                                
C              REQUIRES POLY                                            
C                                                                       
      DIMENSION A(4)                                                    
      DATA  A/5.73738622D-01,2.536103551D-01,-1.45468269D-01,           
     *      1.387472485D-02/, HO/5.815D05/                              
      POLY3(A1,A2,A3,A4,X) = A1+X*(A2+X*(A3+X*A4))                      
C                                                                       
      XI   = H/HO                                                       
      THCL = POLY3(A(1),A(2),A(3),A(4),XI)                              
C                                                                       
      RETURN                                                            
      END                                                               
