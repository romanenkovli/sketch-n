C*DK VISCL                                                              
      FUNCTION VISCL( H, P )                                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C                                                                       
C   DYNAMIC VISCOSITY OF WATER AS FUNCTION OF PRESSURE AND ENTHALPY     
C        PRESSURE  P       IN (N/SQ.M)                                  
C        ENTHALPY  H       IN (J/KG)                                    
C        DYNAMIC VISCOSITY IN (N S/SQ.M)                                
C              REQUIRES POLY                                            
C                                                                       
      DIMENSION A(5) ,B(4) ,D(5) ,E(4) ,F(4)                            
      DATA A/1.299470229D-3,-9.264032108D-4, 3.81047061D-4,             
     *     -8.219444458D-5, 7.022437984D-6/, HO/8.581289699D-6/,        
     *     CON/4.265884D4/, ECON/5.53588D4/, EHO/6.484503981D-6/,       
     *     B/-6.5959D-12, 6.763D-12, -2.88825D-12, 4.4525D-13/,         
     *     D/3.026032306D-4, -1.836606896D-4, 7.567075775D-5,           
     *     -1.647878879D-5, 1.416457633D-6/, HOO/3.892077365D-6/,       
     *     E/1.4526052612D-3, -6.9880084985D-9, 1.5210230334D-14,       
     *     -1.2303194946D-20/, H1/.276D6/, H2/.394D6/, CN/4.014676D5/,  
     *     F/-3.8063507533D-11, 3.9285207677D-16, -1.2585799292D-21,    
     *     1.2860180788D-27/, PI/6.894575293D5/                         
C                                                                       
      POLY3(A1,A2,A3,A4,X) = A1+X*(A2+X*(A3+X*A4))                      
      POLY4(A1,A2,A3,A4,A5,X) = A1+X*(A2+X*(A3+X*(A4+X*A5)))            
      XI  = (H-CON)*HO                                                  
      ETA = (H-ECON)*EHO                                                
      IF( H.GT.H1 ) GO TO 40                                            
C                                                                       
      VISCL = POLY4(A(1),A(2),A(3),A(4),A(5),XI)-POLY3(B(1),B(2),B(3),  
     1        B(4),ETA)*(P-PI)                                          
      RETURN                                                            
40    IF( H.GT.H2 ) GO TO 50                                            
C                                                                       
      VISCL = POLY3(E(1),E(2),E(3),E(4),H)+POLY3(F(1),F(2),F(3),F(4),H)*
     1        (P-PI)                                                    
      RETURN                                                            
C                                                                       
   50 CONTINUE                                                          
      VISCL = POLY4(D(1),D(2),D(3),D(4),D(5),(H-CN)*HOO)                
C                                                                       
      RETURN                                                            
      END                                                               
