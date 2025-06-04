C*DK CPLL                                                               
      FUNCTION CPLL( H, P )                                             
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C                                                                       
C                                                                       
C   SPECIFIC HEAT OF LIQUID WATER AS A FUNCTION OF                      
C   ENTHALPY AND PRESSURE                                               
C            PRESSURE    P   IN  (N/SQ.M)                               
C            ENTHALPY    H   IN  (J/KG)                                 
C            SPEC. HEAT  CPL IN  (J/(KG K))                             
C                        REQUIRES NO SUBPROGRAMS  (S)                   
C                                                                       
      DATA B0,B1/2.394907D-04,-5.196250D-13/                            
      DATA C0,C1/1.193203D-11,2.412704D-18/                             
      DATA D0,D1/-3.944067D-17,-1.680771D-24/                           
C                                                                       
      Z    = H*(H*(D0+D1*P)+(C0+C1*P))+B0+B1*P                          
      CPLL = 1.0D0/Z                                                    
C                                                                       
      RETURN                                                            
      END                                                               
