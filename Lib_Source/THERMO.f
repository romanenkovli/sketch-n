C*DK THERMO                                                             
      SUBROUTINE THERMO                                                 
     1  ( P      ,EL     ,EV     ,TL     ,TV     ,TSAT   ,ROL    ,      
     2    ROV    ,PA     ,ROVA   ,TSSN   ,EVA    ,DR     ,IOP    ,      
     3    JSTART ,NCELLS  )                                             
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION P(1)    ,EL(1)   ,EV(1)   ,TL(1)   ,TV(1)   ,TSAT(1) ,  
     1          ROL(1)  ,ROV(1)  ,PA(1)   ,ROVA(1) ,TSSN(1) ,EVA(1)  ,  
     2          DR(1)                                                   
C                                                                       
C                                                                       
C     THERMODYNAMIC PROPERTIES OF WATER                                 
C                                                                       
C     INPUT VARIABLES                                                   
C        1. P      PRESSURE                                             
C        2. TL     LIQUID TEMPERATURE                                   
C        3. TV     VAPOR TEMPERATURE                                    
C        4. PA     PARTIAL PRESSURE OF THE NONCONDENSIBLE               
C        5. IOP    OPTION SELECTOR - NOT IN PRESENT VERSION             
C                                                                       
C     OUTPUT VARIABLES                                                  
C        1. EL     LIQUID INTERNAL ENERGY                               
C        2. EV     VAPOR (STEAM AND NONCONDENSABLE MIXTURE) INTERNAL    
C                  ENERGY                                               
C        3. TSAT   SATURATION TEMPERATURE CORRESPONDING TO THE TOTAL    
C                  PRESSURE                                             
C        4. ROL    LIQUID DENSITY                                       
C        5. ROV    VAPOR (STEAM AND NONCONDENSABLE MIXTURE) DENSITY     
C        6. ROVA   DENSITY OF THE NONCONDENSABLE                        
C        7. TSSN   SATURATION TEMPERATURE CORRESPONDING TO THE STEAM    
C                  PARTIAL PRESSURE                                     
C        8. EVA    INTERNAL ENERGY OF THE NONCONDENSABLE                
C        9. DTSDP  DERIVATIVE OF TSAT WRT PRESSURE                      
C       10. DELDP  DERIVATIVE OF EL WRT PRESSURE                        
C       11. DEVDP  DERIVATIVE OF STEAM INTERNAL ENERGY WRT STEAM        
C                  PARTIAL PRESSURE                                     
C       12. DELDT  DERIVATIVE OF EL WRT TL                              
C       13. DEVDT  DERIVATIVE OF STEAM INTERNAL ENERGY WRT TV           
C       14. DROLP  DERIVATIVE OF ROL WRT PRESSURE                       
C       15. DROVP  DERIVATIVE OF STEAM DENSITY WRT STEAM PARTIAL        
C                  PRESSURE                                             
C       16. DROLT  DERIVATIVE OF ROL WRT TL                             
C       17. DROVT  DERIVATIVE OF STEAM DENSITY WRT TV                   
C       18. HVST   SATURATED STEAM ENTHALPY (PSTEAM,TSSN)               
C       19. HLST   SATURATED LIQUID ENTHALPY (P,TSSN)                   
C       20. DHVSP  DERIVATIVE OF HVST WRT STEAM PARTIAL PRESSURE        
C       21. DHLSP  DERIVATIVE OF HLST WRT PRESSURE                      
C       22. DTSSP  DERIVATIVE OF TSSN WRT STEAM PARTIAL PRESSURE        
C       23. DEVAT  DERIVATIVE OF EVA WRT TV                             
C       24. DEVAP  DERIVATIVE OF EVA WRT PA                             
C       25. DRVAP  DERIVATIVE OF ROVA WRT PA                            
C       26. DRVAT  DERIVATIVE OF ROVA WRT TV                            
C                                                                       
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C*CA CONTRLLR                                                           
C*CC CONTRLLR                                                           
      INCLUDE 'CONTRLLR'                                                
C*CA DIMNSION                                                           
C*CC DIMNSION                                                           
      INCLUDE 'DIMNSION'                                                
C*CA NMFAIL                                                             
C*CC NMFAIL                                                             
      INCLUDE 'NMFAIL'                                                  
C*CA FIXEDLT                                                            
C*CC FIXEDLT                                                            
C
CITJ+
C
CCC   INCLUDE 'FIXEDLT'                                                 
C
CITJ-
C
C*CA IOUNITS                                                            
C*CC IOUNITS                                                            
      INCLUDE 'IOUNITS'                                                 
C                                                                       
      DIMENSION ALE(8) ,BLE(8) ,CLE(8) ,DLE(8)                          
      DIMENSION AVE(11) ,BVE(11) ,CVE(11) ,DVE(11)                      
      DIMENSION AVG(11) ,BVG(11) ,CVG(11) ,DVG(11)                      
      DIMENSION ACP(10) ,BCP(10) ,CCP(10) ,DCP(10)                      
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C     --------------------------PRELIMINARIES---------------------------
C                                                                       
C                                                                       
C                                                                       
C     DATA FOR RAW CONSTANTS USED IN FITS                               
C                                                                       
      DATA ALE(1)/-1.1436668993222D+06/,BLE(1)/ 4.1868000000000D+03/,   
     *     CLE(1)/ 0.D0               /,DLE(1)/ 0.D0               /    
      DATA ALE(2)/ 8.0957542810383D+06/,BLE(2)/-5.7008855264640D+04/,   
     *     CLE(2)/ 1.3443632119671D+02/,DLE(2)/-9.7879669155946D-02/    
      DATA ALE(3)/-1.9373932457007D+06/,BLE(3)/ 9.7492797103351D+03/,   
     *     CLE(3)/-1.3299615999876D+01/,DLE(3)/ 1.0879999999922D-02/    
      DATA ALE(4)/-5.3245827703670D+06/,BLE(4)/ 2.9179372045334D+04/,   
     *     CLE(4)/-5.0452192000967D+01/,DLE(4)/ 3.4560000000583D-02/    
      DATA ALE(5)/-6.3583523639930D+07/,BLE(5)/ 3.2873715263424D+05/,   
     *     CLE(5)/-5.6371182000208D+02/,DLE(5)/ 3.2760000000116D-01/    
      DATA ALE(6)/-6.6239163195929D+09/,BLE(6)/ 3.1605562257270D+07/,   
     *     CLE(6)/-5.0263730855532D+04/,DLE(6)/ 2.6650075114186D+01/    
      DATA ALE(7)/-5.4759091078157D+09/,BLE(7)/ 2.4635618770681D+07/,   
     *     CLE(7)/-3.6931079506707D+04/,DLE(7)/ 1.8454719393083D+01/    
      DATA ALE(8)/-7.1536399439453D+07/,BLE(8)/ 3.0560801674842D+05/,   
     *     CLE(8)/-4.2424553999630D+02/,DLE(8)/ 1.9719999999823D-01/    
      DATA AVE(1)/ 2.4949771766385D+06/,BVE(1)/ 2.0855856331827D-01/,   
     *     CVE(1)/-1.3553894579716D-07/,DVE(1)/ 2.8522684989198D-14/    
      DATA AVE(2)/ 2.5600870370371D+06/,BVE(2)/ 3.1086111111026D-02/,   
     *     CVE(2)/-6.8988888888580D-09/,DVE(2)/ 4.3203703703379D-16/    
      DATA AVE(3)/ 2.5915500000006D+06/,BVE(3)/ 8.7749999997567D-03/,   
     *     CVE(3)/-1.7499999999663D-09/,DVE(3)/ 4.2999999998503D-17/    
      DATA AVE(4)/ 2.6606000000024D+06/,BVE(4)/-1.3545000000581D-02/,   
     *     CVE(4)/ 6.4250000004682D-10/,DVE(4)/-4.2100000001248D-17/    
      DATA AVE(5)/ 3.8201600000097D+06/,BVE(5)/-2.3019900000170D-01/,   
     *     CVE(5)/ 1.4068900000098D-08/,DVE(5)/-3.1786000000187D-16/    
      DATA AVE(6)/-1.2103411633350D+08/,BVE(6)/ 1.8018803375785D+01/,   
     *     CVE(6)/-8.7442426507726D-07/,DVE(6)/ 1.4091076856088D-14/    
      DATA AVE(7)/ 2.2000000000000D+06/,BVE(7)/ 0.D0               /,   
     *     CVE(7)/ 0.D0               /,DVE(7)/ 0.D0               /    
      DATA AVE(8)/ 2.2000000000000D+06/,BVE(8)/ 0.D0               /,   
     *     CVE(8)/ 0.D0               /,DVE(8)/ 0.D0               /    
      DATA AVE(9)/ 2.2000000000000D+06/,BVE(9)/ 0.D0               /,   
     *     CVE(9)/ 0.D0               /,DVE(9)/ 0.D0               /    
      DATA AVE(10)/ 2.2000000000000D+06/,BVE(10)/ 0.D0               /, 
     *     CVE(10)/ 0.D0               /,DVE(10)/ 0.D0               /  
      DATA AVE(11)/ 2.2000000000000D+06/,BVE(11)/ 0.D0               /, 
     *     CVE(11)/ 0.D0               /,DVE(11)/ 0.D0               /  
      DATA AVG(1)/ 1.0666845123419D+00/,BVG(1)/ 2.8310838172462D-08/,   
     *     CVG(1)/-2.1151097428905D-14/,DVG(1)/ 4.7404001285964D-21/    
      DATA AVG(2)/ 1.0735412407407D+00/,BVG(2)/ 2.6518055555551D-09/,   
     *     CVG(2)/-6.3461111111128D-16/,DVG(2)/ 3.9824074074117D-23/    
      DATA AVG(3)/ 1.0777730000000D+00/,BVG(3)/-2.4300000008021D-11/,   
     *     CVG(3)/-7.1979999998378D-17/,DVG(3)/ 4.8799999990422D-25/    
      DATA AVG(4)/ 1.0851130000007D+00/,BVG(4)/-1.9307000001824D-09/,   
     *     CVG(4)/ 8.9100000014826D-17/,DVG(4)/-3.8960000003946D-24/    
      DATA AVG(5)/ 1.1639800000015D+00/,BVG(5)/-1.6338350000254D-08/,   
     *     CVG(5)/ 9.5856000001448D-16/,DVG(5)/-2.1194000000274D-23/    
      DATA AVG(6)/ 3.8898867259868D+00/,BVG(6)/-3.8595945559811D-07/,   
     *     CVG(6)/ 1.7476370114910D-14/,DVG(6)/-2.6377008249858D-22/    
      DATA AVG(7)/ 2.7168710524682D+00/,BVG(7)/-2.2832718294604D-07/,   
     *     CVG(7)/ 1.0417331983836D-14/,DVG(7)/-1.5842822199773D-22/    
      DATA AVG(8)/ 3.9749829999964D+00/,BVG(8)/-3.0657099999960D-07/,   
     *     CVG(8)/ 1.0637899999985D-14/,DVG(8)/-1.2257999999981D-22/    
      DATA AVG(9)/ 1.2946929999997D+00/,BVG(9)/-2.4834999999979D-08/,   
     *     CVG(9)/ 7.8979999999944D-16/,DVG(9)/-8.0799999999948D-24/    
      DATA AVG(10)/ 1.0590519999963D+00/,BVG(10)/-2.4615999996941D-09/, 
     *     CVG(10)/ 8.8399999991573D-17/,DVG(10)/-8.0799999992269D-25/  
      DATA AVG(11)/ 1.1430199999838D+00/,BVG(11)/-7.7095999988588D-09/, 
     *     CVG(11)/ 1.9335999997331D-16/,DVG(11)/-1.4639999997924D-24/  
      DATA ACP(1)/-7.9678485852270D+02/,BCP(1)/ 2.8187658437259D+01/,   
     *     CCP(1)/-1.0180624999920D-01/,DCP(1)/ 1.2499999999912D-04/    
      DATA ACP(2)/-9.7082632232795D+02/,BCP(2)/ 2.8324981030402D+01/,   
     *     CCP(2)/-9.7656200001157D-02/,DCP(2)/ 1.1600000000110D-04/    
      DATA ACP(3)/-1.6649701690752D+03/,BCP(3)/ 3.3159363169596D+01/,   
     *     CCP(3)/-1.0861179999898D-01/,DCP(3)/ 1.2399999999915D-04/    
      DATA ACP(4)/-6.1420486441088D+03/,BCP(4)/ 6.3630987079837D+01/,   
     *     CCP(4)/-1.7762319999965D-01/,DCP(4)/ 1.7599999999975D-04/    
      DATA ACP(5)/-8.2289951961933D+04/,BCP(5)/ 5.3773958896061D+02/,   
     *     CCP(5)/-1.1612491999609D+00/,DCP(5)/ 8.5599999997375D-04/    
      DATA ACP(6)/-6.5842104212475D+05/,BCP(6)/ 3.7934294783212D+03/,   
     *     CCP(6)/-7.2924928000022D+00/,DCP(6)/ 4.7040000000014D-03/    
      DATA ACP(7)/ 3.4561620732510D+05/,BCP(7)/-2.2129380791446D+02/,   
     *     CCP(7)/-2.4524285999925D+00/,DCP(7)/ 3.1479999999958D-03/    
      DATA ACP(8)/ 1.9798369474597D+06/,BCP(8)/-1.4782551342826D+04/,   
     *     CCP(8)/ 3.1656481897637D+01/,DCP(8)/-2.0843356864237D-02/    
      DATA ACP(9)/-9.6249385211359D+07/,BCP(9)/ 4.3633668884423D+05/,   
     *     CCP(9)/-6.5887615106930D+02/,DCP(9)/ 3.3146147264269D-01/    
      DATA ACP(10)/-1.1074934463333D+07/,BCP(10)/ 4.8073794630970D+04/, 
     *     CCP(10)/-6.9212173247881D+01/,DCP(10)/ 3.3091693999800D-02/  
      DATA C26/0.30D0/                                                  
      DATA CK0/-8.329595D-04/,CK2/-2.245825D-17/,CK4/-1.450382D-06/     
      DATA A11/ .10000887519691D-02/,A12/ .76916250454393D+03/          
      DATA A13/ .13001153775598D-02/                                    
      DATA A15/ .15448787270199D-02/                                    
      DATA DELC0/ 3.737452204297500D+03/                                
      DATA DELC1/ 7.868071497565012D+00/                                
      DATA DELC2/-4.511456897485311D-02/                                
      DATA DELC3/ 1.204554788927213D-04/                                
      DATA DELD0/-.81454700000000D+05/,DELD1/ .89319600000000D+03/      
      DATA DELD2/-.31234800000000D+01/,DELD3/ .37040880000000D-02/      
      DATA DELE0/-.26221567700000D+08/,DELE1/ .22589733400000D+06/      
      DATA DELE2/-.64870195500000D+03/,DELE3/ .62113375200000D+00/      
C                                                                       
C     DEFINITIONS OF COMBINATIONS OF CONSTANTS                          
C     (AS INITIALIZED IN DATA STATEMENTS ABOVE)                         
C                                                                       
C         A11=2.0*C26/(CEOSLP(16)*CEOSLP(12))=2.0/CEOSLP(23)            
C         A12=1.0/A13=CEOSLP(4)/2.0                                     
C         A13=A11*(1.0+C26)=2.0/CEOSLP(4)                               
C         AEOS14=1.0/C28                                                
C         C26=CEOSLP(16)-1.0                                            
C                                                                       
C                                                                       
C     -------------------------EXECUTION LOGIC--------------------------
C                                                                       
CAD +                                                                   
C
CITJ+ (11/13/98)
C
CCC   CALL STIMER( 10, 0, 0, 'THERMO' )                                 
C
CITJ- (11/13/98)
C
CAD -                                                                   
C                                                                       
      DO 300 J = JSTART, NCELLS                                         
         LDTSDP = NTHM*(J-1)+1                                          
         LDELDP = LDTSDP+1                                              
         LDEVDP = LDELDP+1                                              
         LDELDT = LDEVDP+1                                              
         LDEVDT = LDELDT+1                                              
         LDROLP = LDEVDT+1                                              
         LDROVP = LDROLP+1                                              
         LDROLT = LDROVP+1                                              
         LDROVT = LDROLT+1                                              
         LHVST  = LDROVT+1                                              
         LHLST  = LHVST+1                                               
         LDHVSP = LHLST+1                                               
         LDHLSP = LDHVSP+1                                              
         LDTSSP = LDHLSP+1                                              
         LDEVAT = LDTSSP+1                                              
         LDEVAP = LDEVAT+1                                              
         LDRVAP = LDEVAP+1                                              
         LDRVAT = LDRVAP+1                                              
         IF( P(J).GE.CEOSLP(30).AND.P(J).LE.CEOSLP(31) ) GO TO 5        
         P(J) = DMIN1(CEOSLP(31),DMAX1(CEOSLP(30),P(J)))                
         IF( IFTP.EQ.1 ) GO TO 5                                        
C
CITJ+ (11/13/98)
C
CCC      CALL ERROR( 2, 32H*THERMO* PRESSURE LIMIT EXCEEDED, 4 )        
         WRITE(IMOUT,8001)
 8001    FORMAT(1H ,1X,'*THERMO* PRESSURE LIMIT EXCEEDED' )        
         NUM = 999
C
CITJ- (11/13/98)
C
         WRITE(IMOUT,18) NUM ,J                                         
    5    CONTINUE                                                       
C                                                                       
C                                                                       
C     CALCULATE SATURATION PROPERTIES                                   
C                                                                       
         PS         = P(J)-PA(J)                                        
         TSAT(J)    = SATTMP(P(J))                                      
         DR(LDTSDP) = SATDER(P(J),TSAT(J))                              
         IF( PA(J).LT.1.D-5 ) GO TO 9                                   
         TSSN(J)    = SATTMP(PS)                                        
         DR(LDTSSP) = SATDER(PS,TSSN(J))                                
         GO TO 10                                                       
    9    TSSN(J)    = TSAT(J)                                           
         DR(LDTSSP) = DR(LDTSDP)                                        
   10    IF( EV(J).EQ.-1.0D0.AND.IEOS.EQ.0 ) TV(J) = TSAT(J)            
         IF( EL(J).EQ.-1.0D0.AND.IEOS.EQ.0 ) TL(J) = TSAT(J)            
         IF( EV(J).EQ.-2.D0.AND.IEOS.EQ.0 ) TV(J) = TSSN(J)             
         IF( EL(J).EQ.-2.D0.AND.IEOS.EQ.0 ) TL(J) = TSSN(J)             
         IF( IEOS.NE.0 ) GO TO 15                                       
         IF( TV(J).GE.CEOSLP(34).AND.TV(J).LE.CEOSLP(35) ) GO TO 15     
         TV(J) = DMIN1(CEOSLP(35),DMAX1(TV(J),CEOSLP(34)))              
         IF( IFTP.EQ.1 ) GO TO 15                                       
C
CITJ+ (11/13/98)
C
CCC      CALL ERROR( 2, 34H*THERMO* VAPOR TEMP LIMIT EXCEEDED, 4 )      
         WRITE(IMOUT,8002)
 8002    FORMAT(1H ,1X,'*THERMO* VAPOR TEMP LIMIT EXCEEDED' )      
         NUM = 999
C
CITJ- (11/13/98)
C
         WRITE(IMOUT,18) NUM ,J                                         
   15    IF( TL(J).GE.CEOSLP(32).AND.TL(J).LE.CEOSLP(33) ) GO TO 20     
         TL(J) = DMIN1(CEOSLP(33),DMAX1(TL(J),CEOSLP(32)))              
         IF( IFTP.EQ.1 ) GO TO 20                                       
         IF( IEOS.NE.0 ) GO TO 20                                       
C
CITJ+ (11/13/98)
C
CCC      CALL ERROR( 2, 35H*THERMO* LIQUID TEMP LIMIT EXCEEDED, 4 )     
         WRITE(IMOUT,8002)
 8003    FORMAT(1H ,1X,'*THERMO* LIQUID TEMP LIMIT EXCEEDED' )     
         NUM = 999
C
CITJ- (11/13/98)
C
         WRITE(IMOUT,18) NUM ,J                                         
   18 FORMAT(15X,9HCOMPONENT,I4,2H, ,4HCELL,I4)                         
   20    CONTINUE                                                       
C                                                                       
C     CALCULATE LIQUID PROPERTIES                                       
C                                                                       
C                                                                       
C        1. INTERNAL ENERGY AND ITS DERIVATIVES                         
C                                                                       
         II = IDINT(DMIN1(DMAX1((TL(J)-373.15D0)/50.0D0,0.0D0)+1.0D0,   
     1        7.0D0))                                                   
         IF( TL(J).GT.645.15D0 ) II = II+1                              
         EL(J)      = ALE(II)+TL(J)*(BLE(II)+TL(J)*(CLE(II)+TL(J)*      
     1                DLE(II)))                                         
         DR(LDELDT) = BLE(II)+TL(J)*(2.0D0*CLE(II)+TL(J)*3.0D0*DLE(II)) 
         PSL        = SATPRS(TL(J))                                     
         DPS        = 1.0D0/SATDER(PSL,TL(J))                           
         EXPPS      = EXP(CK4*PSL)                                      
         DR(LDELDP) = CK0*(1.0D0-EXPPS)+CK2*PSL*PSL                     
         ELP        = (P(J)-PSL)*DR(LDELDP)                             
         ERT        = DPS*(CK0*(-1.0D0+EXPPS*(1.0D0-CK4*(P(J)-PSL)))+   
     1                CK2*(2.0D0*P(J)*PSL-3.0D0*PSL*PSL))               
         EL(J)      = EL(J)+ELP                                         
         DR(LDELDT) = DR(LDELDT)+ERT                                    
         II         = IDINT(DMIN1(DMAX1((TSSN(J)-373.15D0)/50.0D0,      
     1                0.0D0)+1.0D0,7.0D0))                              
         IF( TSSN(J).GT.645.15D0 ) II = II+1                            
         ELSAT  = ALE(II)+TSSN(J)*(BLE(II)+TSSN(J)*(CLE(II)+TSSN(J)*    
     1            DLE(II)))                                             
         DELSAT = BLE(II)+TSSN(J)*(2.0D0*CLE(II)+TSSN(J)*3.0D0*DLE(II)) 
         CALL RHOLIQ( P(J), TSSN(J), ROLST, DRLSDP, DRLSDT )            
         RROLST     = 1.0D0/ROLST                                       
         DR(LHLST)  = ELSAT+P(J)*RROLST                                 
         DR(LDHLSP) = DELSAT*DR(LDTSDP)+RROLST-P(J)*RROLST*RROLST*      
     1                (DRLSDT*DR(LDTSDP)+DRLSDP)                        
C                                                                       
C        2. DENSITY AND ITS DERIVATIVES                                 
C                                                                       
         CALL RHOLIQ( P(J), TL(J), ROL(J), DR(LDROLP), DR(LDROLT) )     
C                                                                       
C                                                                       
C                                                                       
C     CALCULATE STEAM PROPERTIES                                        
C                                                                       
C                                                                       
C                                                                       
C     PROPERTIES AT SATURATION                                          
C                                                                       
C     -----SPECIFIC HEAT AND ITS DERIVATIVE                             
         II = IDINT(DMIN1(DMAX1((TSSN(J)-273.15D0)/50.0D0,0.0D0)+1.0D0, 
     1        9.0D0))                                                   
         IF( TSSN(J).GT.647.3D0 ) II = II+1                             
         CPS   = ACP(II)+TSSN(J)*(BCP(II)+TSSN(J)*(CCP(II)+TSSN(J)*     
     1           DCP(II)))                                              
         DPCPS = (BCP(II)+TSSN(J)*(2.0D0*CCP(II)+TSSN(J)*3.0D0*         
     1           DCP(II)))*DR(LDTSSP)                                   
C     -----INTERNAL ENERGY, ENTHALPY, AND THEIR DERIVATIVES             
         IF( PS.GT.5.0D+05 ) GO TO 140                                  
         DR(LHVST)  = CEOSLP(8)+CEOSLP(7)*(TSSN(J)-CEOSLP(5))+          
     1                HEV(TSSN(J))                                      
         DR(LDHVSP) = (CEOSLP(7)-2470.2120D0)*DR(LDTSSP)                
         ES         = DR(LHVST)-CEOSLP(12)*TSSN(J)                      
         DPES       = DR(LDHVSP)-CEOSLP(12)*DR(LDTSSP)                  
         GAMS       = DR(LHVST)/ES                                      
         GAMSM      = GAMS-1.0D0                                        
         DPGAMS     = (DR(LDHVSP)-GAMS*DPES)/ES                         
         GO TO 200                                                      
  140    CONTINUE                                                       
         II = IDINT(DMIN1(PS/50.0D+05+1.0D0,9.0D0))                     
         IF( PS.GT.20.0D+05 ) II = II+1                                 
         IF( PS.GT.220.0D+05 ) II = II+1                                
         ES         = AVE(II)+PS*(BVE(II)+PS*(CVE(II)+PS*DVE(II)))      
         DPES       = (BVE(II)+PS*(2.0D0*CVE(II)+PS*3.0D0*DVE(II)))     
         GAMS       = AVG(II)+PS*(BVG(II)+PS*(CVG(II)+PS*DVG(II)))      
         DPGAMS     = BVG(II)+PS*(2.0D0*CVG(II)+PS*3.0D0*DVG(II))       
         GAMSM      = GAMS-1.0D0                                        
         DR(LHVST)  = GAMS*ES                                           
         DR(LDHVSP) = GAMS*DPES+ES*DPGAMS                               
  200    CONTINUE                                                       
C                                                                       
C     PROPERTIES AT ACTUAL STEAM TEMPERATURE                            
C                                                                       
         IF( IEOS.NE.0 ) GO TO 280                                      
         IF( TV(J).GT.TSSN(J) ) GO TO 210                               
         DR(LDEVDT) = CPS/CEOSLP(16)                                    
         DE         = DR(LDEVDT)*(TV(J)-TSSN(J))                        
         EV(J)      = ES+DE                                             
         DR(LDEVDP) = DPES+DE*DPCPS/CPS-DR(LDEVDT)*DR(LDTSSP)           
         GO TO 220                                                      
  210    CONTINUE                                                       
         T1         = 1.0D0/(A11*CPS-1.0D0)                             
         BETA       = TSSN(J)*TSSN(J)*(1.0D0-T1*T1)                     
         T2         = TSSN(J)*T1                                        
         DE         = A12*(TV(J)-TSSN(J)+SQRT(TV(J)*TV(J)-BETA)-T2)     
         EV(J)      = ES+DE                                             
         CAPK       = A13*DE+TSSN(J)+T2                                 
         DBETAP     = 2.0D0*(BETA*DR(LDTSSP)+A11*DPCPS*T2**3)/TSSN(J)   
         DCAPKP     = -A13*DPES+(1.0D0+T1)*DR(LDTSSP)-A11*T1*T2*DPCPS   
         T3         = 1.0D0-BETA/(CAPK*CAPK)                            
         DR(LDEVDT) = CEOSLP(4)/T3                                      
         DR(LDEVDP) = -0.5D0*(T3*DCAPKP+DBETAP/CAPK)*DR(LDEVDT)         
  220    CONTINUE                                                       
         T4         = 1.0D0/(GAMSM*ES+C26*DE)                           
         ROV(J)     = PS*T4                                             
         DRVDE      = -ROV(J)*C26*T4                                    
         DR(LDROVT) = DRVDE*DR(LDEVDT)                                  
         DR(LDROVP) = T4-ROV(J)*(ES*DPGAMS+(GAMSM-C26)*DPES)*T4+DRVDE*  
     1                DR(LDEVDP)                                        
         IF( PS.GT.0.0D0.AND.ROV(J).LE.0.0D0 ) GO TO 230                
         IF( ROV(J).LT.(0.999D0*ROL(J)) ) GO TO 280                     
         ROV(J)     = 0.999D0*ROL(J)                                    
         DR(LDROVT) = 0.999D0*DR(LDROLT)                                
         DR(LDROVP) = 0.999D0*DR(LDROLP)                                
         GO TO 280                                                      
  230    CONTINUE                                                       
         ROV(J)     = PS/(CEOSLP(12)*TV(J))                             
         DR(LDROVT) = -ROV(J)/TV(J)                                     
         DR(LDROVP) = ROV(J)/PS                                         
C                                                                       
C                                                                       
C                                                                       
C     CALCULATE AIR PROPERTIES                                          
C                                                                       
C                                                                       
C                                                                       
  280    EVA(J)     = CEOSLP(17)*TV(J)                                  
         DR(LDEVAT) = CEOSLP(17)                                        
         DR(LDEVAP) = 0.0D0                                             
         DR(LDRVAP) = 1.0D0/(CEOSLP(25)*TV(J))                          
         ROVA(J)    = DR(LDRVAP)*PA(J)                                  
         DR(LDRVAT) = -CEOSLP(25)*ROVA(J)*DR(LDRVAP)                    
C                                                                       
C                                                                       
C                                                                       
C     CALCULATE AIR-STEAM MIXTURE PROPERTIES                            
C                                                                       
C                                                                       
C                                                                       
         IF( IEOS.EQ.0 ) GO TO 281                                      
         EV(J)      = EVA(J)                                            
         ROV(J)     = DR(LDRVAP)*PS                                     
         DR(LDEVDT) = DR(LDEVAT)                                        
         DR(LDEVDP) = DR(LDEVAP)                                        
         DR(LDROVP) = DR(LDRVAP)                                        
         DR(LDROVT) = -CEOSLP(25)*ROV(J)*DR(LDRVAP)                     
  281    CONTINUE                                                       
         EV(J)  = (EV(J)*ROV(J)+EVA(J)*ROVA(J))                         
         ROV(J) = ROV(J)+ROVA(J)                                        
         EV(J)  = EV(J)/ROV(J)                                          
  300 CONTINUE                                                          
CAD +                                                                   
C
CITJ+ (11/13/98)
C
CCC   CALL TTIMER( 10, 0, 0 )                                           
C
CITJ- (11/13/98)
C
CAD -                                                                   
      RETURN                                                            
      END                                                               
