C*DK FPROP                                                              
      SUBROUTINE FPROP                                                  
     1  ( P    ,EL   ,EV   ,ROL  ,ROV  ,TL   ,TV   ,TSAT ,HFG  ,CPL  ,  
     2    CPV  ,VISL ,VISV ,CL   ,CV   ,SIG  ,NC   ,PA   ,ROA  ,TSSN ,  
     3    EA   ,DER   )                                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C     CALCULATES VALUES FOR FLUID ENTHALPY, TRANSPORT PROPERTIES,       
C     AND SURFACE TENSION.                                              
C                                                                       
C*CA DIMNSION                                                           
C*CC DIMNSION                                                           
      INCLUDE 'DIMNSION'                                                
      DIMENSION P(NC)    ,EL(NC)   ,EV(NC)   ,ROL(NC)  ,ROV(NC)  ,      
     1          TL(NC)   ,TV(NC)   ,TSAT(NC) ,HFG(NC)  ,CPL(NC)  ,      
     2          CPV(NC)  ,VISL(NC) ,VISV(NC) ,CL(NC)   ,CV(NC)   ,      
     3          SIG(NC)  ,PA(NC)   ,ROA(NC)  ,TSSN(NC) ,EA(NC)   ,      
     4          DER(1)                                                  
C                                                                       
      IF( NC.LE.0 ) RETURN                                              
CAD +                                                                   
C
CITJ+ (11/11/98)
C
CCC   CALL STIMER( 11, 0, 0, 'FPROP' )                                  
C
CITJ- (11/11/98)
C
CAD -                                                                   
      DO 100 I = 1, NC                                                  
         IPOINT = NTHM*(I-1)                                            
C                                                                       
         HL = EL(I)+P(I)/ROL(I)                                         
         HV = EV(I)+P(I)/ROV(I)                                         
C                                                                       
         HFG(I) = DER(IPOINT+10)-DER(IPOINT+11)                         
C                                                                       
C<<ADD>>                                                                
         IF( TL(I).GT.TSAT(I) ) HL = DER(IPOINT+11)                     
         CPL(I) = CPLL(HL,P(I))                                         
         CPV(I) = CPVV1(TV(I),P(I),PA(I))                               
C                                                                       
         VISL(I) = VISCL(HL,P(I))                                       
         VISV(I) = VISCV(HV,P(I),ROV(I),TV(I),PA(I))                    
C                                                                       
         CL(I) = THCL(HL)                                               
         CV(I) = THCV(HV,P(I),ROV(I),TV(I),PA(I))                       
C                                                                       
         SIG(I) = SIGMA(P(I),TSAT(I))                                   
C                                                                       
  100 CONTINUE                                                          
CAD +                                                                   
C
CITJ+ (11/11/98)
C
CCC   CALL TTIMER( 11, 0, 0 )                                           
C
CITJ- (11/11/98)
C
CAD -                                                                   
c      write(*,*) 'CPL =', (CPL(i),i=1,NC)
c      write(*,*) 'HL =', HL
c      write(*,*) 'P(I) =', (p(i), i=1, NC)
c      pause
      
      
      RETURN                                                            
      END                                                               
