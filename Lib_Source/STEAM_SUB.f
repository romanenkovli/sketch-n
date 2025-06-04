      subroutine JAERI_STEAM( NIN, TZ, PZ, RO_OUT , CP_OUT,
     &           Viscos_OUT, Conduct_Out)
	                                                                          
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*30    COM(21)                                                   
      CHARACTER*10   UNIT(21)                                                   
      DIMENSION IQUICK(9), fdata(21), ZTAB(19)                                          
                                                                                
      DATA IQUICK(1) / 1 /                                                      
      DATA IQUICK(2) / 2 /                                                      
      DATA IQUICK(3) / 3 /                                                      
      DATA IQUICK(4) / 4 /                                                      
      DATA IQUICK(5) / 16/                                                      
      DATA IQUICK(6) / 17/                                                      
      DATA IQUICK(7) / 18/                                                      
      DATA IQUICK(8) / 20/                                                      
      DATA IQUICK(9) / 21/                                                      
                                                                                
      DATA COM( 1) /'Pressure                      '/                           
      DATA COM( 2) /'Temperature                   '/                           
      DATA COM( 3) /'Density                       '/                           
      DATA COM( 4) /'Enthalpy                      '/                           
      DATA COM( 5) /'Entolopy                      '/                           
      DATA COM( 6) /'Internal Energy               '/                           
      DATA COM( 7) /'Isobaric Specific Heat  CP    '/                           
      DATA COM( 8) /'         Specific Heat  CV    '/                           
      DATA COM( 9) /'Squared Sound Velocity C**2   '/                           
      DATA COM(10) /'[D(Enthalpy)/D(Density)| P    '/                           
      DATA COM(11) /'1/V[D(Volume)/D(Temp)| P      '/                           
      DATA COM(12) /'CP / CV                       '/                           
      DATA COM(13) /'-1/V [D(Volume)/D(Press.)| T  '/                           
      DATA COM(14) /'Specific GIBBS Function       '/                           
      DATA COM(15) /'Specific HELMHOLTZ Function   '/                           
      DATA COM(16) /'Dynamic Viscocity             '/                           
      DATA COM(17) /'Kinematic Viscocity           '/                           
      DATA COM(18) /'Thermal Conductivity          '/                           
      DATA COM(19) /'Thermal Diffusivity           '/                           
      DATA COM(20) /'Prendtle Number               '/                           
      DATA COM(21) /'Surface Tension               '/                           
                                                                                
      DATA UNIT( 1) /' MPa      '/                                              
      DATA UNIT( 2) /' K        '/                                              
      DATA UNIT( 3) /' kG/m3    '/                                              
      DATA UNIT( 4) /' J/kG     '/                                              
      DATA UNIT( 5) /' J/kG/K   '/                                              
      DATA UNIT( 6) /' J/kG     '/                                              
      DATA UNIT( 7) /' J/kG/K   '/                                              
      DATA UNIT( 8) /' J/kG/K   '/                                              
      DATA UNIT( 9) /' k2/s2    '/                                              
      DATA UNIT(10) /' J m3/kg2 '/                                              
      DATA UNIT(11) /' 1/K      '/                                              
      DATA UNIT(12) /' -        '/                                              
      DATA UNIT(13) /' m2/N     '/                                              
      DATA UNIT(14) /' J/kG     '/                                              
      DATA UNIT(15) /' J/kG     '/                                              
      DATA UNIT(16) /' kG/m2/s  '/                                              
      DATA UNIT(17) /' m2/s     '/                                              
      DATA UNIT(18) /' W/m/K    '/                                              
      DATA UNIT(19) /' m2/s     '/                                              
      DATA UNIT(20) /' -        '/                                              
      DATA UNIT(21) /' N/m      '/                                              
                                                                                
c  100 write(6,'(1h1)')                
c      write(6,*) ' Mode Select   K = C + 273.15,  Pa = 1.0197D-5 kgf/m2'                 
c      write(6,*) '   0 : P(MPa) & T(K) : subcooled liquid or superheated 
c     &                steam '                 
c      write(6,*) '   1 : T(K)          : saturated steam'          
c      write(6,*) '   2 : P(MPa)        :           steam'          
c      write(6,*) '   3 : T(K)          : saturated liquid'         
c      write(6,*) '   4 : P(MPa)        :           liquid'         
c      WRITE(6,*) '  '                                                      
c      READ(5,*) NIN                                                             
                                                                                
      PZ = PZ*1.0D6                                                             
      CALL STEAMS(PZ,TZ,NIN,NOUT,ZTAB)                                          


C                                       WRITE(6,*) ' CHECK AFTER CALL PZ,TZ,NIN',PZ,TZ,NIN                       
      P = PZ / 1.0D6                                                            
      T = TZ             
	                                                       
      ZTAB(1) = 1.0 / ZTAB(1)                                                   

                                                                            
      DO 150 I = 1,21                                                           
      IF(I.EQ.1) FDATA(1) = P                                                   
      IF(I.EQ.2) FDATA(2) = T                                                   
      IF(I.GE.3) THEN                                                           
         J = I - 2                                                              
         FDATA(I) = ZTAB(J)                                                     
      ENDIF                                                                     
  150 CONTINUE                                                                  

                                                                                
c      WRITE(6,*) ' Detailed Data ?  YES => 1,  Simple-Set > others'          
c      READ(5,*) I                                                               
c      IF(I.eq.1) GO TO 400                                                      
                                                                                
c      DO 200 I = 1,9                                                            
c      II = IQUICK(I)                                                            
c 200  WRITE(6,10) COM(II),FDATA(II),UNIT(II)
c      go to 450
 
c 400  continue                                                                                
c      DO 300 I = 1,21                                                           
c  300 WRITE(6,10) COM(I),FDATA(I),UNIT(I)                                       
                                                                                
c  450 WRITE(6,*) '  '            
c      WRITE(6,*) ' Continue ?   YES => 1 ,  NO => others'            
c      READ(5,*) I                                                               
c      IF(I.EQ.1) GO TO 100                                                      
                                                                                
   10 FORMAT(1h ,A30,1PE12.5,A10)                
   
      RO_OUT = FDATA(3)                               
	CP_OUT = FDATA(7)
	Viscos_Out = FDATA(16)
	Conduct_Out = FDATA(18)

      return
      END                                                                       
C
C=================================================================================================
C
      SUBROUTINE STEAMS(P,T,NIN,NOUT,ZTAB)                                      
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION ZTAB(19),ZKTAB(18),NTAB(4)                                      
C                                                                               
      PK = P * 1.0197D-5                                                        
      TK = T - 273.15                                                           
      NTAB(2) = 0                                                               
      NTAB(4) = 1                                                               
C                                                                               
      IF(NIN.EQ.1) THEN                                                         
            PK = PSAT(TK)                                                       
            P = PK * 9.8067D4                                                   
            NTAB(2) = 1                                                         
      ENDIF                                                                     
C                                                                               
      IF(NIN.EQ.2) THEN                                                         
            TK = TSAT(PK)                                                       
            T = TK + 273.15                                                     
            NTAB(2) = 1                                                         
      ENDIF                                                                     
C                                                                               
      IF(NIN.EQ.3) THEN                                                         
            PK = PSAT(TK)                                                       
            P = PK*9.8067D4                                                     
            NTAB(2) = 2                                                         
      ENDIF                                                                     
C                                                                               
      IF(NIN.EQ.4) THEN                                                         
            TK = TSAT(PK)                                                       
            T = TK + 273.15                                                     
            NTAB(2) = 2                                                         
      ENDIF                                                                     
C                                                                               
      CALL STEAMF(PK,TK,ZKTAB(1),ZKTAB(2),ZKTAB(3),ZKTAB(4),                    
     1            ZKTAB(5),ZKTAB(6),ZKTAB(7),ZKTAB(8),                          
     1            ZKTAB(9),ZKTAB(10),ZKTAB(11),ZKTAB(12),                       
     1            ZKTAB(13),ZKTAB(14),ZKTAB(15),ZKTAB(16),                      
     1            ZKTAB(17),ZKTAB(18),NTAB(1),NTAB(2),NTAB(3),NTAB(4) )         
C                                                                               
C                                                                               
      DO 100 I = 1,18                                                           
      ZTAB(I) = ZKTAB(I)                                                        
  100 CONTINUE                                                                  
C                                                                               
C                                                                               
C   SPECIFIC VOLUME     M3/KG => M3/KG                                          
C                                                                               
      ZTAB(1) = ZKTAB(1) * 1.00000                                              
C                                                                               
C   SPECIFIC ENTHALPY   KCAL/KG => J/KG                                         
C                                                                               
      ZTAB(2) = ZKTAB(2) * 4.1860D3                                             
C                                                                               
C   SPECIFIC ENTROPY    KCAL/KG C => J/KG K                                     
C                                                                               
      ZTAB(3) = ZKTAB(3) * 4.1860D3                                             
C                                                                               
C   SPECIFIC ENERGY     KCAL/KG => J/KG                                         
C                                                                               
      ZTAB(4) = ZKTAB(4) * 4.1860D3                                             
C                                                                               
C   SPECIFIC HEAT   CP    KCAL/KG C => J/KG K                                   
C                                                                               
      ZTAB(5) = ZKTAB(5) * 4.1860D3                                             
C                                                                               
C   SPECIFIC HEAT   CV    KCAL/KG C => J/KG K                                   
C                                                                               
      ZTAB(6) = ZKTAB(6) * 4.1860D3                                             
C                                                                               
C   SQUARED LIGHT VELOSITY  M2/SEC2 => M2/SEC2                                  
C                                                                               
      ZTAB(7) = ZKTAB(7) * 1.00000                                              
C                                                                               
C   DH/DRO     KCAL M3/KG2 => J M3/KG2                                          
C                                                                               
      ZTAB(8) = ZKTAB(8) * 4.1860D3                                             
C                                                                               
C                       1/K => 1/K                                              
C                                                                               
      ZTAB(9) = ZKTAB(9) * 1.00000                                              
C                                                                               
C   GAMMA                     =>                                                
C                                                                               
      ZTAB(10) = ZKTAB(10) * 1.00000                                            
C                                                                               
C   COMPRESION FACTOR    M2/N => M2/N                                           
C                                                                               
      ZTAB(11) = ZKTAB(11) * 1.00000                                            
C                                                                               
C   GIBBS FREE ENERGY     KCAL/KG => J/KG                                       
C                                                                               
      ZTAB(12) = ZKTAB(12) * 4.1860D3                                           
C                                                                               
C   HELMHOLTZ FREE ENERGY  M3/KG => M3/KG                                       
C                                                                               
      ZTAB(13) = ZKTAB(13) * 4.1860D3                                           
C                                                                               
C   VISCOSITY     MICRO POISE => KG/M2 S                                        
C                                                                               
      ZTAB(14) = ZKTAB(14) * 1.0D-7                                             
C                                                                               
C   KINEMATIC VISCOSITY    M2/S => M2/S                                         
C                                                                               
      ZTAB(15) = ZKTAB(15) * 1.00000                                            
C                                                                               
C   THERMAL CONDUCTIVITY  KCAL/M H C => W/M C                                   
C                                                                               
      ZTAB(16) = ZKTAB(16) * 1.16280                                            
C                                                                               
C   TERMAL DUFISITIVITY    M2/H => M2/S                                         
C                                                                               
      ZTAB(17) = ZKTAB(17) / 3600.00                                            
C                                                                               
C   PRANDTLE NUMBER                   =>                                        
C                                                                               
      ZTAB(18) = ZKTAB(18) * 1.00000                                            
C                                                                               
C   SURFACE TENSION                   =>                                        
C                                                                               
      CALL SYGCAL( T,ZTAB(19) )                                                 
C                                                                               
C                                                                               
      NOUT=0                                                                    
C                                                                               
      IF(ZKTAB(14).EQ.0.0.AND.ZKTAB(16).EQ.0.0)  THEN                           
             NOUT=3                                                             
             CALL VISCAL (T,1.0/ZTAB(1),VIS)                                    
               DVIS = VIS*ZTAB(1)                                               
             CALL CONCAL (T,1.0/ZTAB(1),CON)                                    
               TDIF = CON*ZTAB(1)/ZTAB(5)                                       
             PR = DVIS / TDIF                                                   
C                                                                               
             ZTAB(14) = VIS                                                     
             ZTAB(15) = DVIS                                                    
             ZTAB(16) = CON                                                     
             ZTAB(17) = TDIF                                                    
             ZTAB(18) = PR                                                      
             GO TO 1000                                                         
      ENDIF                                                                     
C                                                                               
C                                                                               
      IF(ZKTAB(16).EQ.0.0) THEN                                                 
             NOUT = 2                                                           
             CALL CONCAL(T,1.0/ZTAB(1),CON)                                     
             TDIF = CON*ZTAB(1)/ZTAB(5)                                         
             PR = ZTAB(15)/TDIF                                                 
C                                                                               
             ZTAB(16)=CON                                                       
             ZTAB(17)=TDIF                                                      
             ZTAB(18)=PR                                                        
             GO TO 1000                                                         
       ENDIF                                                                    
C                                                                               
C                                                                               
      IF(ZKTAB(14).EQ.0.0) THEN                                                 
             NOUT = 1                                                           
             CALL VISCAL(T,1.0/ZTAB(1),VIS)                                     
             DVIS = VIS*ZTAB(1)                                                 
             PR = DVIS/ZTAB(17)                                                 
C                                                                               
             ZTAB(14)=VIS                                                       
             ZTAB(15)=DVIS                                                      
             ZTAB(18)=PR                                                        
             GO TO 1000                                                         
       ENDIF                                                                    
C                                                                               
C                                                                               
 1000 CONTINUE                                                                  
C                                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE CONCAL ( T,DEN,CON )                                           
C                                                                               
C     CONDUCTIVITY CALCULATION                                                  
C                                                                               
C                              REF.1980 JSME STEAM TABLE                        
C                                                                               
C        SI UNIT     DEN : DENSITY     KG/M3                                    
C                      T ; TEMPERATURE     K                                    
C                    CON ; THERMAL CONDUCTIVITY  W/M K                          
C                                                                               
      IMPLICIT REAL*8  (A-H,O-Z)                                                
C                                                                               
      DATA  TSTA / 647.30 /                                                     
      DATA  ROSTA / 317.7 /                                                     
      DATA  BB1 / -1.71587D-1 /                                                 
      DATA  BB2 / 2.39219 /                                                     
      DATA  C1 / 6.42857D-1 /                                                   
      DATA  C2 / -4.11717 /                                                     
      DATA  C3 / -6.17937 /                                                     
      DATA  C4 / 3.08976D-3 /                                                   
      DATA  C5 / 8.22994D-2 /                                                   
      DATA  C6 / 1.00932D1 /                                                    
      DATA  A0 / 1.02811D-2 /                                                   
      DATA  A1 / 2.99621D-2 /                                                   
      DATA  A2 / 1.56146D-2 /                                                   
      DATA  A3 / -4.22464D-3 /                                                  
      DATA  B0 / -3.97070D-1 /                                                  
      DATA  B1 / 4.00302D-1 /                                                   
      DATA  B2 / 1.06 /                                                         
      DATA  D1 / 7.01309D-2 /                                                   
      DATA  D2 / 1.18520D-2 /                                                   
      DATA  D3 / 1.69937D-3 /                                                   
      DATA  D4 / -1.02 /                                                        
C                                                                               
C                                                                               
      DTSTA = ABS( T / TSTA - 1.0 ) + C4                                        
      TRAT = T / TSTA                                                           
      IF( TRAT.GE.1.0 ) THEN                                                    
           S = DTSTA**(-1.0)                                                    
           ELSE                                                                 
           S = C6*DTSTA*(-0.6)                                                  
      ENDIF                                                                     
C                                                                               
      Q = 2.0 + C5*DTSTA**(-0.6)                                                
      R = Q + 1.0                                                               
C                                                                               
      RAM0 = TRAT**0.5 *( A0 + A1*TRAT + A2*TRAT**2 + A3*TRAT**3 )              
      RAMBAR = B0 + B1*DEN/ROSTA + B2*DEXP( BB1*(DEN/ROSTA+BB2)**2 )            
      DELRAM = ( D1*(1/TRAT)**10+D2 ) * (DEN/ROSTA)**1.8 *                      
     1         DEXP( C1*( 1.0-(DEN/ROSTA)**2.8 ) ) +                            
     2         D3*S*(DEN/ROSTA)**Q*DEXP( Q/R*(1.0-(DEN/ROSTA)**R) )             
     3       + D4*DEXP( C2*TRAT**1.5 + C3*(ROSTA/DEN)**5 )                      
C                                                                               
C                                                                               
      CON = RAM0 + RAMBAR + DELRAM                                              
C                                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE VISCAL( T,DEN,VIS )                                            
C                                                                               
C        PROGRAM FOR VISCOSITY CALCULATION                                      
C                                                                               
C                                                                               
C           T IN K   DEN IN KG/M3  VIS IS PA.S                                  
C                                  REFFERENCE    1980 JSME STEAM TABLE          
C                                                                               
      IMPLICIT  REAL*8  (A-H,O-Z)                                               
C                                                                               
      DIMENSION A(4),B(6,5)                                                     
C                                                                               
      DATA   TSTA /647.27/                                                      
      DATA   DENSTA /317.763/                                                   
      DATA   A /  0.0181583,  0.0177624,  0.0105287, -0.0036744 /               
      DATA   B /                                                                
     * 0.501938,  0.162888, -0.130356,  0.907919, -0.551119,  0.146543,         
     * 0.235622,  0.789393,  0.673665,  1.207552,  0.0670665, -0.084337,        
     *-0.274637, -0.743539, -0.959456, -0.687343, -0.497089,  0.195286,         
     * 0.145831,  0.263129,  0.347247,  0.213486,  0.100754, -0.032932,         
     *-0.0270448,-0.0253093,-0.0267758,-0.0822904, 0.0602253,-0.0202595/        
C                                                                               
C                                                                               
      SUM = 0.0                                                                 
C                                                                               
      DO 100 I = 1,4                                                            
C                                                                               
      SUM = SUM + A(I) * ( TSTA / T )**(I-1)                                    
C                                                                               
  100 CONTINUE                                                                  
C                                                                               
      VIS0 = DSQRT( T/TSTA ) / SUM                                              
C                                                                               
      SUM = 0.0                                                                 
C                                                                               
      DO 200 I = 1,6                                                            
C                                                                               
      DO 200 J = 1,5                                                            
C                                                                               
      SUM = SUM + B(I,J) * ( TSTA / T - 1.0 )**(I-1) *                          
     *                     ( DEN / DENSTA - 1.0 )**(J-1)                        
C                                                                               
  200 CONTINUE                                                                  
C                                                                               
      VIS = VIS0 * DEXP( DEN / DENSTA * SUM )                                   
      VIS = 1.0D-6 * VIS                                                        
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE  SYGCAL(T,SYG)                                                 
C                                                                               
C       SURFACE TENSION CALCULATION                                             
C                                                                               
C                     T IN K    ;  SYG IN N/M2                                  
C                                                                               
      REAL*8 SYG,T,TC                                                           
C                                                                               
      DATA  TC / 647.15 /                                                       
C                                                                               
      SYG = 0.2358 * (1.0-T/TC)**1.256*( 1.0-0.625*(1.0-T/TC) )                 
C                                                                               
      RETURN                                                                    
      END                                                                       
                                                                                
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
      SUBROUTINE STEAMF (P,T,SPVOL,ENTAL,ENTRO,EI,CP,CV,CC,HRHO,BETHA,          
     1                   GAMMA,AKAP,G,F,VIS,DVIS,TCON,A,PR,NREGON,              
     2                   NSTWAT,IPR,IZV)                                        
C                                                                               
C     NREGON=1 -- SUBREGION 1 (SUBROUTINE SUB1)                                 
C     NREGON=2 -- SUBREGION 2 (SUBROUTINE SUB2)                                 
C     NREGON=3 -- SUBREGION 3 (SUBROUTINE SUB3)                                 
C     NREGON=4 -- SUBREGION 4 (SUBROUTINE SUB4)                                 
C                                                                               
C     NSTWAT=0 -- COMPRESSUED WATER OR SUPERHEATED STEAM                        
C     NSTWAT=1 -- SATURATED STEAM                                               
C     NSTWAT=2 -- SATURATED WATER                                               
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      PB=0.980665*P                                                             
      BET=PB/221.2                                                              
      TK=T+273.15                                                               
      SIT=TK/647.3                                                              
      BTK=BETK(SIT)                                                             
      CALL SUBBET(SIT,BTL,BETLD,BETLDD)                                         
      BT2=4.520795660                                                           
      SIT1=0.9626911787                                                         
      SIT2=1.333462073                                                          
      SIT3=1.657886606                                                          
      SITT=273.15/647.3                                                         
      IF ((SIT.GE.SITT.AND.SIT.LE.SIT1).AND.(NSTWAT.EQ.2)) GO TO 1              
      IF ((SIT.GE.SITT.AND.SIT.LE.SIT1).AND.(NSTWAT.EQ.1)) GO TO 2              
      IF ((SIT.GT.SIT1.AND.SIT.LE.1.0).AND.(NSTWAT.EQ.1)) GO TO 3               
      IF ((SIT.GT.SIT1.AND.SIT.LE.1.0).AND.(NSTWAT.EQ.2)) GO TO 4               
      IF ((SIT.GE.SITT.AND.SIT.LE.SIT1).AND.(BET.GT.BTK.AND.BET.LE.BT2)         
     1   ) GO TO 1                                                              
      IF ((SIT.GE.SITT.AND.SIT.LE.SIT1).AND.(BET.GE.0.0.AND.BET.LT.BTK))        
     1   GO TO 2                                                                
      IF ((SIT.GT.SIT1.AND.SIT.LT.SIT2).AND.(BET.GE.0.0.AND.BET.LE.BTL))        
     1   GO TO 2                                                                
      IF ((SIT.GE.SIT2.AND.SIT.LE.SIT3).AND.(BET.GE.0.0.AND.BET.LE.BT2))        
     1   GO TO 2                                                                
      IF ((SIT.GT.SIT1.AND.SIT.LT.1.0).AND.(BET.GT.BTL.AND.BET.LT.BTK))         
     1   GO TO 3                                                                
      IF ((SIT.GE.1.0.AND.SIT.LT.SIT2).AND.(BET.GT.BTL.AND.BET.LE.BT2))         
     1   GO TO 3                                                                
      IF ((SIT.GT.SIT1.AND.SIT.LT.1.0).AND.(BET.GT.BTK.AND.BET.LE.BT2))         
     1   GO TO 4                                                                
      GO TO 5                                                                   
    1 NREGON=1                                                                  
      CALL SUB1 (SIT,BET,ZET,CHI,SIG,EPS,FAI,FAV,RHRHO,RCC,RBETHA,RKAP,         
     1           RGAM,IZV)                                                      
      GO TO 6                                                                   
    2 NREGON=2                                                                  
      CALL SUB2 (SIT,BET,ZET,CHI,SIG,EPS,FAI,FAV,RHRHO,RCC,RBETHA,RKAP,         
     1           RGAM,IZV)                                                      
      GO TO 6                                                                   
    3 NREGON=3                                                                  
      CALL SUB3 (SIT,BET,PSI,CHI,SIG,EPS,FAI,FAV,RHRHO,RCC,RBETHA,RKAP,         
     1           RGAM,IZV)                                                      
      GO TO 6                                                                   
    4 NREGON=4                                                                  
      CALL SUB4 (SIT,BET,PSI,CHI,SIG,EPS,FAI,FAV,RHRHO,RCC,RBETHA,RKAP,         
     1           RGAM,IZV)                                                      
      GO TO 6                                                                   
    5 STOP 'STEAMF (P,T,SPVOL,ENTAL,ENTRO,CP,CV,G,F,VIS,DVIS,TCON,A,PR,N        
     1REGON,NSTWAT)'                                                            
    6 SPVOL=CHI*0.00317                                                         
      ENTAL=EPS*70120.4*0.238846*0.001                                          
      ENTRO=SIG*108.3275143*0.238846*0.001                                      
      IF(IZV.EQ.0)RETURN                                                        
      TC1=647.3                                                                 
      PC1=22120000.                                                             
      VC1=.00317                                                                
      CP=FAI*108.3275143*0.238846*0.001                                         
      IF(IZV.EQ.-1)RETURN                                                       
      CV=FAV*108.3275143*0.238846*0.001                                         
      CC=RCC*PC1*VC1                                                            
      HRHO=RHRHO*VC1*VC1*PC1*.238846*.001                                       
      BETHA=RBETHA/TC1                                                          
      GAMMA=RGAM                                                                
      AKAP=RKAP/PC1                                                             
      EI=ENTAL-SPVOL*P*0.234228*100.                                            
      F=EI-TK*ENTRO                                                             
      G=ENTAL-TK*ENTRO                                                          
      CALL VISCON(P,T,SPVOL,VISMP,TCONMW,NREGON,NSTWAT)                         
      VIS=VISMP                                                                 
      DVIS=VIS*SPVOL*1.000000E-07                                               
      TCON=0.859845*0.001*TCONMW                                                
C     WRITE(6,10) NREGON                                                        
   10 FORMAT(I10)                                                               
C     RETURN                                                                    
  101 CONTINUE                                                                  
      A=TCON*SPVOL/CP                                                           
      IF(DABS(A).GE.1.0D-20) GO TO 20                                           
      IPR=0                                                                     
      PR=0.                                                                     
      RETURN                                                                    
   20 IPR=1                                                                     
      PR=DVIS*3600.0/A                                                          
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION PSAT(T)                                         
C                                                                               
C     PRESSURE OF SATURATED STEAM AS A FUNCTION OF TEMPERATURE                  
C                                                                               
C     PSAT = SATURATION PRESSURE (KG/CM2)                                       
C     T    = TEMPERATURE (OC) (0-374.15 OC)                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      SIT=(T+273.15)/647.3                                                      
      BTK=BETK(SIT)                                                             
      PSAT=BTK*221.2/0.980665                                                   
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION TSAT(P)                                         
C                                                                               
C     TEMPERATURE OF SATURATED STEAM AS A FUNCTION OF PRESSURE                  
C                                                                               
C     TSAT = SATURATION TEMPERATURE (OC)                                        
C     P    = PRESSURE (KG/CM2) (0.006228-225.65 KG/CM2)                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      BET=P*0.980665/221.2                                                      
      TSAT=BETKP(BET)*647.3-273.15                                              
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION TSATA(P)                                        
C                                                                               
C     TEMPERATURE OF SATURATED STEAM AND WATER AS A FUNCTION OF PRESSURE        
C                                                                               
C     TSATA = SATURATED TEMPERATURE (OC)                                        
C                                                                               
C     P = PRESSURE (KG/CM2)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      IF(P.GE.0.1) GO TO 1                                                      
      C0=-1.4027184                                                             
      C1=2.9622845                                                              
      C2=-1.1582744                                                             
      C3=3.1894969E-1                                                           
      C4=-5.5043472E-2                                                          
      C5=5.6724437E-3                                                           
      C6=-3.1813374E-4                                                          
      C7=7.4526790E-6                                                           
      XS=100.                                                                   
      YS=0.1                                                                    
      X=P*XS                                                                    
      Y=(C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*C7)))))))/YS                  
      GO TO 5                                                                   
    1 IF(P.GE.1.5) GO TO 2                                                      
      C0=1.3909155                                                              
      C1=4.4918984                                                              
      C2=-1.7145823                                                             
      C3=3.9849678E-1                                                           
      C4=-5.0037269E-2                                                          
      C5=3.1144544E-3                                                           
      C6=-7.4511981E-5                                                          
      XS=10.                                                                    
      YS=0.1                                                                    
      X=P*XS                                                                    
      Y=(C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*C6))))))/YS                         
      GO TO 5                                                                   
    2 IF(P.GE.10.) GO TO 3                                                      
      C0=5.4115033E-1                                                           
      C1=7.1971995E-1                                                           
      C2=-4.0200320E-1                                                          
      C3=1.7558544E-1                                                           
      C4=-5.3055241E-2                                                          
      C5=1.0852398E-2                                                           
      C6=-1.4695072E-3                                                          
      C7=1.2579400E-4                                                           
      C8=-6.1513953E-6                                                          
      C9=1.3072177E-7                                                           
      XS=1.                                                                     
      YS=0.01                                                                   
      X=P*XS                                                                    
      Y=(C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*(C8+X*C9))))))))        
     1)/YS                                                                      
      GO TO 5                                                                   
    3 IF(P.GE.100.) GO TO 4                                                     
      C0=1.0164701                                                              
      C1=1.3449215                                                              
      C2=-9.0888377E-1                                                          
      C3=4.6834233E-1                                                           
      C4=-1.6377515E-1                                                          
      C5=3.8669958E-2                                                           
      C6=-6.1522423E-3                                                          
      C7=6.4864736E-4                                                           
      C8=-4.3375412E-5                                                          
      C9=1.6638736E-6                                                           
      C10=-2.7863220E-8                                                         
      XS=0.1                                                                    
      YS=0.01                                                                   
      X=P*XS                                                                    
      Y=(C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*(C8+X*(C9+X*C10)        
     1)))))))))/YS                                                              
      GO TO 5                                                                   
    4 C0=2.0679536                                                              
      C1=8.4703533E-1                                                           
      C2=1.6565001                                                              
      C3=-3.1971445                                                             
      C4=2.6767174                                                              
      C5=-1.2202267                                                             
      C6=2.9377612E-1                                                           
      C7=-2.9311442E-2                                                          
      XS=0.01                                                                   
      YS=0.01                                                                   
      X=P*XS                                                                    
      Y=(C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*C7)))))))/YS                  
    5 TSATA=Y                                                                   
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE SUBBET(SIT,BETL,BETLD,BETLDD)                                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DOUBLE PRECISION L0,L1,L2                                                 
C                                                                               
      L0=1.574373327D+01                                                        
      L1=-3.417061978D+01                                                       
      L2=1.931380707D+01                                                        
      BETL=L0+SIT*(L1+SIT*L2)                                                   
      BETLD=L1+2.0*L2*SIT                                                       
      BETLDD=2.0*L2                                                             
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C     THE OPTION 'FORT='AUTODBL(DBLPAD4)' ' SHOULD BE USED                      
C     IN COMPILING THIS PROGRAM BY THE 'FORTRAN-HE COMPILER.                    
C     OPTION DOUBLE                                                             
      SUBROUTINE SUB1 (SIT,BET,ZETA,CHI1,SIG1,EPS1,FAI1,FAV1,RHRHO1,            
     1                 RCC1,RBETHA,RKAP1,RGAM1,IZV)                             
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
C     SUBREGION 1                                                               
C                                                                               
      A0=6.824687741D+03                                                        
      A1=-5.422063673D+02                                                       
      A2=-2.096666205D+04                                                       
      A3=3.941286787D+04                                                        
      A4=-6.733277739D+04                                                       
      A5=9.902381028D+04                                                        
      A6=-1.093911774D+05                                                       
      A7=8.590841667D+04                                                        
      A8=-4.511168742D+04                                                       
      A9=1.418138926D+04                                                        
      A10=-2.017271113D+03                                                      
      A11=7.982692717D+00                                                       
      A12=-2.616571843D-02                                                      
      A13=1.522411790D-03                                                       
      A14=2.284279054D-02                                                       
      A15=2.421647003D+02                                                       
      A16=1.269716088D-10                                                       
      A17=2.074838328D-07                                                       
      A18=2.174020350D-08                                                       
      A19=1.105710498D-09                                                       
      A20=1.293441934D+01                                                       
      A21=1.308119072D-05                                                       
      A22=6.047626338D-14                                                       
      AS1=8.438375405D-01                                                       
      AS2=5.362162162D-04                                                       
      AS3=1.720000000D+00                                                       
      AS4=7.342278489D-02                                                       
      AS5=4.975858870D-02                                                       
      AS6=6.537154300D-01                                                       
      AS7=1.150000000D-06                                                       
      AS8=1.510800000D-05                                                       
      AS9=1.418800000D-01                                                       
      AS10=7.002753165D+00                                                      
      AS11=2.995284926D-04                                                      
      AS12=2.040000000D-01                                                      
      Y=1.0-AS1*SIT**2-AS2*SIT**(-6)                                            
      Z=Y+DSQRT(AS3*Y**2-2.0*AS4*SIT+2.0*AS5*BET)                               
      ZETA=A0*SIT*(1.0-DLOG(SIT))+A1+SIT*(A2+SIT*(A3+SIT*(A4+SIT*(A5            
     1     +SIT*(A6+SIT*(A7+SIT*(A8+SIT*(A9+SIT*A10))))))))+A11*(17.0/          
     2     29.0*Z-17.0/12.0*Y)*Z**(12./17.)+(A12+A13*SIT+A14*SIT**2+A15*        
     3     (AS6-SIT)**10+A16*(AS7+SIT**19)**(-1))*BET-(AS8+SIT**11)**           
     4     (-1)*BET*(A17+BET*(A18+BET*A19))-A20*SIT**18*(AS9+SIT**2)*           
     5     ((AS10+BET)**(-3)+AS11*BET)+A21*(AS12-SIT)*BET**3+A22*SIT**          
     6     (-20)*BET**4                                                         
      CHI1=A11*AS5*Z**(-5./17.)+(A12+A13*SIT+A14*SIT**2+A15*(AS6-SIT)**         
     1     10+A16*(AS7+SIT**19)**(-1))-(AS8+SIT**11)**(-1)*(A17+2.0*A18*        
     2     BET+3.0*A19*BET**2)-A20*SIT**18*(AS9+SIT**2)*(-3.0*(AS10+BET)        
     3     **(-4)+AS11)+3.0*A21*(AS12-SIT)*BET**2+4.0*A22*SIT**(-20)*           
     4     BET**3                                                               
      YS=-2.0*AS1*SIT+6.0*AS2*SIT**(-7)                                         
      SIG1=A0*DLOG(SIT)-A2-SIT*(2.0*A3+SIT*(3.0*A4+SIT*(4.0*A5+SIT*(5.0         
     1     *A6+SIT*(6.0*A7+SIT*(7.0*A8+SIT*(8.0*A9+SIT*9.0*A10)))))))           
     2     +A11*((5.0/12.0*Z-(AS3-1.0)*Y)*YS+AS4)*Z**(-5./17.)+(-A13-2.*        
     3     A14*SIT+10.0*A15*(AS6-SIT)**9+19.0*A16*(AS7+SIT**19)**(-2)*          
     4     SIT**18)*BET-11.0*(AS8+SIT**11)**(-2)*SIT**10*(A17*BET+A18*          
     5     BET**2+A19*BET**3)+A20*SIT**17*(18.0*AS9+20.0*SIT**2)*((AS10         
     6     +BET)**(-3)+AS11*BET)+A21*BET**3+20.0*A22*SIT**(-21)*BET**4          
      EPS1=ZETA+SIT*SIG1                                                        
      IF(IZV.EQ.0)RETURN                                                        
      YSS=-2.0*AS1-42.0*AS2*SIT**(-8)                                           
      ZS=YS+(AS3*Y*YS-AS4)*(AS3*Y**2-2.0*AS4*SIT+2.0*AS5*BET)**(-0.5)           
      ZSS=YSS+AS3*(AS3*Y**2-2.0*AS4*SIT+2.0*AS5*BET)**(-0.5)*(YS**2+Y*          
     1    YSS)-(AS3*Y*YS-AS4)**2*(AS3*Y**2-2.0*AS4*SIT+2.0*AS5*BET)**           
     2    (-1.5)                                                                
      ZESS=-A0*SIT**(-1)+1.0*2.0*A3+SIT*(2.0*3.0*A4+SIT*(3.0*4.0*A5+SIT*        
     1     (4.0*5.0*A6+SIT*(5.0*6.0*A7+SIT*(6.0*7.0*A8+SIT*(7.0*8.0*A9+         
     2     SIT*8.*9.*A10))))))+A11*((12./29.*Z-Y)*(Z**(-5./17.)*ZSS-5./1        
     3     7.*Z**(-22./17.)*ZS**2)+(24./29.*ZS-2.*YS)*Z**(-5./17.)*ZS+(         
     4     17./29.*ZSS-17./12.*YSS)*Z**(12./17.))+BET*(2.0*A14+90.0*A15*        
     5     (AS6-SIT)**8+722.0*A16*SIT**36*(AS7+SIT**19)**(-3)-342.0*A16*        
     6     SIT**17*(AS7+SIT**19)**(-2))-(242.0*SIT**20*(AS8+SIT**11)**          
     7     (-3)-110.0*SIT**9*(AS8+SIT**11)**(-2))*(A17*BET+A18*BET**2+          
     8     A19*BET**3)-A20*SIT**16*(306.0*AS9+380.0*SIT**2)*((AS10+BET)         
     9     **(-3)+AS11*BET)+420.0*A22*SIT**(-22)*BET**4                         
      FAI1=-SIT*ZESS                                                            
      IF(IZV.EQ.-1)RETURN                                                       
      ZB=AS5/(Z-Y)                                                              
      ZSB=-0.5*AS5/(Z-Y)**3*(2.0*AS3*Y*YS-2.0*AS4)                              
      ZBB=-AS5**2/(Z-Y)**3                                                      
      ZESB=-                                                                    
     A     5.0/17.0*A11*ZB*Z**(-22.0/17.0)*(Z*YS+(AS3-1.0)*Y*YS-AS4)+           
     1     (A13+2.0*A14*SIT-10.0*A15*(AS6-SIT)**9-19.0*A16*(AS7+SIT**19)        
     2     **(-2)*SIT**18)+11.0*(AS8+SIT**11)**(-2)*SIT**10*(A17+2.0*A18        
     3     *BET+3.0*A19*BET**2)-2.0*A20*SIT**17*(9.0*AS9+10.0*SIT**2)*(         
     4     -3.0*(AS10+BET)**(-4)+AS11)-3.0*A21*BET**2-80.0*A22*SIT**            
     5     (-21)*BET**3                                                         
      ZEBB=-5.0/17.0*A11*AS5*Z**(-22.0/17.0)*ZB-(AS8+SIT**11)**(-1)*(2.0        
     1     *A18+6.0*A19*BET)-12.0*A20*SIT**18*(AS9+SIT**2)*(AS10+BET)**         
     2     (-5)+6.0*A21*(AS12-SIT)*BET+12.0*A22*SIT**(-20)*BET**2               
      FAV1=FAI1+SIT*ZESB**2/ZEBB                                                
      CHI1SQ=CHI1*CHI1                                                          
      ZEBBSQ=ZEBB*ZEBB                                                          
      RHRHO1=CHI1SQ*SIT*ZESS/ZESB                                               
      ZESBSQ=ZESB*ZESB                                                          
      RCC1=-CHI1SQ/(ZEBB-ZESBSQ/ZESS)                                           
      RBETHA=ZESB/CHI1                                                          
      RKAP1=-ZEBB/CHI1                                                          
      RGAM1=FAI1/FAV1                                                           
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE SUB2 (SIT,BET,ZETB,CHI2,SIG2,EPS2,FAI2,FAV2,RHRHO2,            
     1                  RCC2,RBETHA,RKAP2,RGAM2,IZV)                            
C                                                                               
C     SUBREGION 2                                                               
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION B(8,3),BS(8,2)                                                  
C                                                                               
      B0= 1.683599274D+01                                                       
      B01=2.856067796D+01                                                       
      B02=-5.438923329D+01                                                      
      B03=4.330662834D-01                                                       
      B04=-6.547711697D-01                                                      
      B05=8.565182058D-02                                                       
      B90=1.936587558D+02                                                       
      B91=-1.388522425D+03                                                      
      B92=4.126607219D+03                                                       
      B93=-6.508211677D+03                                                      
      B94=5.745984054D+03                                                       
      B95=-2.693088365D+03                                                      
      B96=5.235718623D+02                                                       
      BB=7.633333333D-01                                                        
      B(1,1)=6.670375918D-02                                                    
      B(1,2)=1.388983801D+00                                                    
      B(1,3)=0.                                                                 
      B(2,1)=8.390104328D-02                                                    
      B(2,2)=2.614670893D-02                                                    
      B(2,3)=-3.373439453D-02                                                   
      B(3,1)=4.520918904D-01                                                    
      B(3,2)=1.069036614D-01                                                    
      B(3,3)=0.                                                                 
      B(4,1)=-5.975336707D-01                                                   
      B(4,2)=-8.847535804D-02                                                   
      B(4,3)=0.                                                                 
      B(5,1)=5.958051609D-01                                                    
      B(5,2)=-5.159303373D-01                                                   
      B(5,3)=2.075021122D-01                                                    
      B(6,1)=1.190610271D-01                                                    
      B(6,2)=-9.867174132D-02                                                   
      B(6,3)=0.                                                                 
      B(7,1)=1.683998803D-01                                                    
      B(7,2)=-5.809438001D-02                                                   
      B(7,3)=0.                                                                 
      B(8,1)=6.552390126D-03                                                    
      B(8,2)=5.710218649D-04                                                    
      B(8,3)=0.                                                                 
      BS(1,1)=0.                                                                
      BS(1,2)=0.                                                                
      BS(2,1)=0.                                                                
      BS(2,2)=0.                                                                
      BS(3,1)=0.                                                                
      BS(3,2)=0.                                                                
      BS(4,1)=0.                                                                
      BS(4,2)=0.                                                                
      BS(5,1)=0.                                                                
      BS(5,2)=0.                                                                
      BS(6,1)=4.006073948D-01                                                   
      BS(6,2)=0.                                                                
      BS(7,1)=8.636081627D-02                                                   
      BS(7,2)=0.                                                                
      BS(8,1)=-8.532322921D-01                                                  
      BS(8,2)=3.460208861D-01                                                   
      TC1=647.3                                                                 
      PC1=22120000.0                                                            
      VC1=0.00317                                                               
      R1=461.51                                                                 
      AI1=R1*TC1/(PC1*VC1)                                                      
      X=DEXP(BB*(1.0-SIT))                                                      
      CALL SUBBET(SIT,BETL,BETLD,BETLDD)                                        
      S01=0.0                                                                   
      S02=0.0                                                                   
      S03=0.0                                                                   
      S04=0.0                                                                   
      S05=0.0                                                                   
      S06=0.0                                                                   
      S07=0.0                                                                   
      S08=0.0                                                                   
      S09=0.0                                                                   
      S10=0.0                                                                   
      DO 4 I=1,8                                                                
      BX=0.0                                                                    
      ZBX=0.0                                                                   
      BZX=0.0                                                                   
      BXA=0.0                                                                   
      BXB=0.0                                                                   
      BXC=0.0                                                                   
      NNI=NN(I)                                                                 
      NLI=NL(I)                                                                 
      DO 2 J=1,NNI                                                              
      SBX=0.0                                                                   
      XBX=0.0                                                                   
      NZIJ=NZ(I,J)                                                              
      FNZIJ=FLOAT(NZIJ)                                                         
      DO 1 K=1,NLI                                                              
      NXIK=NX(I,K)                                                              
      FNXIK=FLOAT(NXIK)                                                         
      SBX=SBX+BS(I,K)*X**NXIK                                                   
      XBX=XBX+FNXIK*BS(I,K)*X**NXIK                                             
    1 CONTINUE                                                                  
      BX=BX+B(I,J)*X**NZIJ                                                      
      ZBX=ZBX+FNZIJ*B(I,J)*X**NZIJ                                              
      BZX=BZX+B(I,J)*(1.0+FNZIJ*BB*SIT)*X**NZIJ                                 
      BXA=BXA+B(I,J)*X**NZIJ*(FNZIJ-XBX/(BET**(2-I)+SBX))                       
      BXB=BXB+B(I,J)*X**NZIJ*(-FNZIJ*(BET**(2-I)+SBX)+2.0*XBX)                  
      BXC=BXC+B(I,J)*X**NZIJ*FLOAT(2-I)*(FLOAT(1-I)*BET**(-I)*(BET**            
     1    (2-I)+SBX)-2.0*FLOAT(2-I)*BET**(2*(1-I)))                             
    2 CONTINUE                                                                  
      IF (I.GE.6) GO TO 3                                                       
      S01=S01+BET**I*BX                                                         
      S03=S03+FLOAT(I)*BET**(I-1)*BX                                            
      S05=S05+BET**I*ZBX                                                        
      S07=S07+FLOAT(I)*BET**(I-1)*ZBX                                           
      S09=S09+FLOAT(I*(I-1))*BET**(I-2)*BX                                      
      GO TO 4                                                                   
    3 S02=S02+BX/(BET**(2-I)+SBX)                                               
      S04=S04+FLOAT(I-2)*BET**(1-I)*BX/(BET**(2-I)+SBX)**2                      
      S06=S06+BXA/(BET**(2-I)+SBX)                                              
      S08=S08+FLOAT(2-I)*BET**(1-I)*BXB/(BET**(2-I)+SBX)**2/(BET**(2-I)+        
     1    SBX)                                                                  
      S10=S10+BXC/(BET**(2-I)+SBX)**2/(BET**(2-I)+SBX)                          
    4 CONTINUE                                                                  
      ZETB=AI1*SIT*DLOG(BET)+B0*SIT*(1.0-DLOG(SIT))+B01+SIT*(B02+SIT*           
     1 (B03+SIT*(B04+SIT*B05)))-S01-S02+BET*(BET/BETL)**10*(B90+X*(B91+         
     2     X*(B92+X*(B93+X*(B94+X*(B95+X*B96))))))                              
      CHI2=AI1*SIT/BET-S03-S04+11.0*(BET/BETL)**10*(B90+X*(B91+X*(B92+X*        
     1     (B93+X*(B94+X*(B95+X*B96))))))                                       
      BBLD=10.0*BETLD/BETL                                                      
      SIG2=-AI1*DLOG(BET)+B0*DLOG(SIT)-B02-SIT*(2.0*B03+SIT*(3.0*B04+SIT        
     1     *4.0*B05))-BB*S05-BB*S06+BET*(BET/BETL)**10*(BBLD*B90+(BBLD+         
     2     BB)*B91*X+(BBLD+2.0*BB)*B92*X**2+(BBLD+3.0*BB)*B93*X**3+(BBLD        
     3     +4.0*BB)*B94*X**4+(BBLD+5.0*BB)*B95*X**5+(BBLD+6.0*BB)*B96*X         
     4     **6)                                                                 
      EPS2=ZETB+SIT*SIG2                                                        
           IF(IZV.EQ.0)RETURN                                                   
      BX61=B(6,1)*X**12*BET**4                                                  
      BX62=B(6,2)*X**11*BET**4                                                  
      BX71=B(7,1)*X**24*BET**5                                                  
      BX72=B(7,2)*X**18*BET**5                                                  
      BX81=B(8,1)*X**24*BET**6                                                  
      BX82=B(8,2)*X**14*BET**6                                                  
      BSX61=BS(6,1)*X**14*BET**4                                                
      BSX71=BS(7,1)*X**19*BET**5                                                
      BSX81=BS(8,1)*X**54*BET**6                                                
      BSX82=BS(8,2)*X**27*BET**6                                                
      F60=BX61+BX62                                                             
      F61=-BB*(12.0*BX61+11.0*BX62)                                             
      F62=BB**2*(144.0*BX61+121.0*BX62)                                         
      F70=BX71+BX72                                                             
      F71=-BB*(24.0*BX71+18.0*BX72)                                             
      F72=BB**2*(576.0*BX71+324.0*BX72)                                         
      F80=BX81+BX82                                                             
      F81=-BB*(24.0*BX81+14.0*BX82)                                             
      F82=BB**2*(576.0*BX81+196.0*BX82)                                         
      G60=1.0+BSX61                                                             
      G61=-14.0*BB*BSX61                                                        
      G62=196.0*BB**2*BSX61                                                     
      G70=1.0+BSX71                                                             
      G71=-19.0*BB*BSX71                                                        
      G72=361.0*BB**2*BSX71                                                     
      G80=1.0+BSX81+BSX82                                                       
      G81=-BB*(54.0*BSX81+27.0*BSX82)                                           
      G82=BB**2*(2916.0*BSX81+729.0*BSX82)                                      
      ZESS=-B0*SIT**(-1)+2.0*B03+6.0*B04*SIT+12.0*B05*SIT**2-(169.0*B(1,        
     1  1)*X**13+9.0*B(1,2)*X**3)*BB**2*BET-(324.0*B(2,1)*X**18+4.0*            
     2     B(2,2)*X**2+B(2,3)*X)*BB**2*BET**2-(324.0*B(3,1)*X**18+100.*         
     3     B(3,2)*X**10)*BB**2*BET**3-(625.0*B(4,1)*X**25+196.0*B(4,2)*         
     4     X**14)*BB**2*BET**4-(1024.0*B(5,1)*X**32+784.0*B(5,2)*X**28          
     5     +576.0*B(5,3)*X**24)*BB**2*BET**5-(F62/G60-(2.0*F61*G61+F60*         
     6     G62)/G60**2+2.0*F60*G61**2/G60**3)-(F72/G70-(2.0*F71*G71+F70         
     7     *G72)/G70**2+2.0*F70*G71**2/G70**3)-(F82/G80-(2.0*F81*G81+F80        
     8     *G82)/G80**2+2.0*F80*G81**2/G80**3)                                  
      ZESS=ZESS+(BET/BETL)**11*((110.0*BETLD**2/BETL-10.0*BETLDD)*(B90+         
     1     X*(B91+X*(B92+X*(B93+X*(B94+X*(B95+X*B96))))))+20.0*BB*BETLD*        
     2     X*(B91+X*(2.0*B92+X*(3.0*B93+X*(4.0*B94+X*(5.0*B95+X*6.0*B96)        
     3     ))))+BB**2*BETL*X*(B91+X*(4.0*B92+X*(9.0*B93+X*(16.0*B94+X*(         
     4     25.0*B95+X*36.0*B96))))))                                            
      FAI2=-SIT*ZESS                                                            
      IF(IZV.EQ.-1)RETURN                                                       
      ZESB=AI1/BET+BB*(S07+S08)-11.0*(BET/BETL)**10*(B90*BBLD+X*(B91*(          
     1     BBLD+BB)+X*(B92*(BBLD+2.0*BB)+X*(B93*(BBLD+3.0*BB)+X*(B94*(          
     2     BBLD+4.0*BB)+X*(B95*(BBLD+5.0*BB)+X*B96*(BBLD+6.0*BB)))))))          
      ZEBB=-AI1*SIT/BET**2-S09+S10+110.0/BETL*(BET/BETL)**9*(B90+X*(B91+        
     1     X*(B92+X*(B93+X*(B94+X*(B95+X*B96))))))                              
      FAV2=FAI2+SIT*ZESB**2/ZEBB                                                
      CHI2SQ=CHI2*CHI2                                                          
      ZESBSQ=ZESB*ZESB                                                          
      RHRHO2=CHI2SQ*SIT*ZESS/ZESB                                               
      RCC2=-CHI2SQ/(ZEBB-ZESBSQ/ZESS)                                           
      RBETHA=ZESB/CHI2                                                          
      RKAP2=-ZEBB/CHI2                                                          
      RGAM2=FAI2/FAV2                                                           
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE SUB3 (SIT,BET,PSIC,CHI,SIG3,EPS3,FAI3,FAV3,RHRHO3,RCC3,        
     1                 RBETHA,RKAP3,RGAM3,IZV)                                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
C     SUBREGION 3                                                               
C                                                                               
      C00=-6.839900000D+00                                                      
      C01=-1.722604200D-02                                                      
      C02=-7.771750390D+00                                                      
      C03=4.204607520D+00                                                       
      C04=-2.768070380D+00                                                      
      C05=2.104197070D+00                                                       
      C06=-1.146495880D+00                                                      
      C07=2.231380850D-01                                                       
      C08=1.162503630D-01                                                       
      C09=-8.209005440D-02                                                      
      C010=1.941292390D-02                                                      
      C011=-1.694705760D-03                                                     
      C012=-4.311577033D+00                                                     
      C11=7.086360850D-01                                                       
      C12=1.236794550D+01                                                       
      C13=-1.203890040D+01                                                      
      C14=5.404374220D+00                                                       
      C15=-9.938650430D-01                                                      
      C16=6.275231820D-02                                                       
      C17=-7.747430160D+00                                                      
      C21=-4.298850920D+00                                                      
      C22=4.314305380D+01                                                       
      C23=-1.416193130D+01                                                      
      C24=4.041724590D+00                                                       
      C25=1.555463260D+00                                                       
      C26=-1.665689350D+00                                                      
      C27=3.248811580D-01                                                       
      C28=2.936553250D+01                                                       
      C31=7.948418420D-06                                                       
      C32=8.088597470D+01                                                       
      C33=-8.361533800D+01                                                      
      C34=3.586365170D+01                                                       
      C35=7.518959540D+00                                                       
      C36=-1.261606400D+01                                                      
      C37=1.097174620D+00                                                       
      C38=2.121454920D+00                                                       
      C39=-5.465295660D-01                                                      
      C310=8.328754130D+00                                                      
      C40=2.759717760D-06                                                       
      C41=-5.090739850D-04                                                      
      C50=2.106363320D+02                                                       
      C60=5.528935335D-02                                                       
      C61=-2.336365955D-01                                                      
      C62=3.697071420D-01                                                       
      C63=-2.596415470D-01                                                      
      C64=6.828087013D-02                                                       
      C70=-2.571600553D+02                                                      
      C71=-1.518783715D+02                                                      
      C72=2.220723208D+01                                                       
      C73=-1.802039570D+02                                                      
      C74=2.357096220D+03                                                       
      C75=-1.462335698D+04                                                      
      C76=4.542916630D+04                                                       
      C77=-7.053556432D+04                                                      
      C78=4.381571428D+04                                                       
      SX=SIT-1.0                                                                
      SX1=1.0/SIT                                                               
      CH0=CHI3(SIT,BET)                                                         
      N=0                                                                       
    1 CX=1.0/CH0                                                                
      F3=-(C01-CX**2*(C02+CX*(2.0*C03+CX*(3.0*C04+CX*(4.0*C05+CX*(5.0*          
     1   C06+CX*(6.0*C07+CX*(7.0*C08+CX*(8.0*C09+CX*(9.0*C010+CX*10.0*          
     2   C011)))))))))+C012*CX)-(C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*            
     3   C14+CX*(4.0*C15+CX*5.0*C16))))+C17*CX)*SX-(C21-CX**2*(C22+CX*          
     4   (2.0*C23+CX*(3.0*C24+CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))           
     5   +C28*CX)*SX**2-(C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(            
     6   4.0*C35+CX*(5.0*C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))          
     7   +C310*CX)*SX**3+5.0*C41*CH0**(-6)*SIT**(-23)*SX-6.0*CH0**5*            
     8   SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))-BET                 
      F3D=-CX**3*(2.0*C02+CX*(6.0*C03+CX*(12.0*C04+CX*(20.0*C05+CX*(30.0        
     1    *C06+CX*(42.0*C07+CX*(56.0*C08+CX*(72.0*C09+CX*(90.0*C010+CX*         
     2    110.0*C011)))))))))+C012*CX**2-(CX**3*(2.0*C12+CX*(6.0*C13+CX*        
     3    (12.0*C14+CX*(20.0*C15+CX*30.0*C16))))-C17*CX**2)*SX-(CX**3*(         
     4    2.0*C22+CX*(6.0*C23+CX*(12.0*C24+CX*(20.0*C25+CX*(30.0*C26+CX*        
     5    42.0*C27)))))-C28*CX**2)*SX**2-(CX**3*(2.0*C32+CX*(6.0*C33+CX*        
     6    (12.0*C34+CX*(20.0*C35+CX*(30.0*C36+CX*(42.0*C37+CX*(56.0*C38+        
     7    CX*72.0*C39)))))))-C310*CX**2)*SX**3-30.0*C41*CX**7*SIT**(-23)        
     8    *SX-30.0*CH0**4*SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64        
     9    ))))                                                                  
      CHI=CH0-F3/F3D                                                            
      IF (DABS(CH0/CHI-1.0).LE.1.0D-07) GO TO 2                                 
      N=N+1                                                                     
      IF(N.EQ.50) WRITE(6,10) N,CH0,CHI                                         
   10 FORMAT(10X,'SUB3',I5,2E20.7)                                              
      IF (N.GE.50) GO TO 2                                                      
      CH0=CHI                                                                   
      GO TO 1                                                                   
    2 CONTINUE                                                                  
      CX=1.0/CHI                                                                
      PSIC=C00+C01*CHI+CX*(C02+CX*(C03+CX*(C04+CX*(C05+CX*(C06+CX*(C07+         
     1     CX*(C08+CX*(C09+CX*(C010+CX*C011)))))))))+C012*DLOG(CHI)+(C11        
     2     *CHI+CX*(C12+CX*(C13+CX*(C14+CX*(C15+CX*C16))))+C17*DLOG(CHI)        
     3     )*SX+(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(C25+CX*(C26+CX*C27)        
     4     ))))+C28*DLOG(CHI))*SX**2+(C31*CHI+CX*(C32+CX*(C33+CX*(C34+CX        
     5     *(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)))))))+C310*DLOG(CHI))*         
     6     SX**3+(C40+C41*CHI**(-5))*SIT**(-23)*SX+C50*SIT*DLOG(SIT)+CHI        
     7     **6*SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))+SX*           
     8     (C70+SX*(C71+SX*(C72+SX*(C73+SX*(C74+SX*(C75+SX*(C76+SX*(C77+        
     9     SX*C78))))))))                                                       
      SIG3=-(C11*CHI+CX*(C12+CX*(C13+CX*(C14+CX*(C15+CX*C16))))                 
     1     +C17*DLOG(CHI)+C50)-2.0*(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(        
     2     C25+CX*(C26+CX*C27)))))+C28*DLOG(CHI))*SX-3.0*(C31*CHI+CX*           
     3     (C32+CX*(C33+CX*(C34+CX*(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)         
     4     ))))))+C310*DLOG(CHI))*SX**2+(C40+C41*CHI**(-5))*(22.0*SIT**         
     5     (-23)-23.0*SIT**(-24))-C50*DLOG(SIT)+CHI**6*SX1**3*(2.0*C60+         
     6     SX1*(3.0*C61+SX1*(4.0*C62+SX1*(5.0*C63+SX1*6.0*C64))))-C70-          
     7     SX*(2.0*C71+SX*(3.0*C72+SX*(4.0*C73+SX*(5.0*C74+SX*(6.0*C75+         
     8     SX*(7.0*C76+SX*(8.0*C77+SX*9.0*C78)))))))                            
      EPS3=PSIC+SIT*SIG3+CHI*BET                                                
      IF(IZV.EQ.0)RETURN                                                        
      PSICD1=2.0*(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(C25+CX*(C26+CX*C27        
     1       )))))+C28*DLOG(CHI))+6.0*(C31*CHI+CX*(C32+CX*(C33+CX*(C34+         
     2       CX*(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)))))))+C310*DLOG(CHI        
     3       ))*SX+(C40+C41*CHI**(-5))*(506.0*SIT-552.0)*SIT**(-25)+C50*        
     4       SX1+CHI**6*SX1**4*(6.0*C60+SX1*(12.0*C61+SX1*(20.0*C62+SX1*        
     5       (30.0*C63+SX1*42.0*C64))))+2.0*C71+SX*(6.0*C72+SX*(12.0*C73        
     6       +SX*(20.0*C74+SX*(30.0*C75+SX*(42.0*C76+SX*(56.0*C77+SX*72.        
     7       0*C78))))))                                                        
      PSICD2=C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*C14+CX*(4.0*C15+CX*5.0*         
     1       C16))))+C17*CX+2.0*(C21-CX**2*(C22+CX*(2.0*C23+CX*(3.0*C24         
     2       +CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))+C28*CX)*SX+3.0*(          
     3       C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(4.0*C35+CX*(5.0*        
     4       C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))+C310*CX)*SX**        
     5       2-5.0*C41*CHI**(-6)*(23.0-22.0*SIT)*SIT**(-24)-6.0*CHI**5*         
     6       SX1**3*(2.0*C60+SX1*(3.0*C61+SX1*(4.0*C62+SX1*(5.0*C63+SX1*        
     7       6.0*C64))))                                                        
      PSICD3=-F3D                                                               
      CHISQ=CHI*CHI                                                             
      RHRHO3=CHISQ*SIT*(-PSICD1*PSICD3/PSICD2+PSICD2)                           
      RCC3=CHISQ*(PSICD3-PSICD2*PSICD2/PSICD1)                                  
      RKAP3=1./CHI/PSICD3                                                       
      RBETHA=-RKAP3*PSICD2                                                      
      FAV3=-SIT*PSICD1                                                          
      FAI3=FAV3+SIT*PSICD2**2/PSICD3                                            
      RGAM3=FAI3/FAV3                                                           
    3 RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE SUB4 (SIT,BET,PSID,CHI,SIG4,EPS4,FAI4,FAV4,RHRHO4,             
     1                 RCC4,RBETHA,RKAP4,RGAM4,IZV)                             
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
C     SUBREGION 4                                                               
C                                                                               
      C00=-6.839900000D+00                                                      
      C01=-1.722604200D-02                                                      
      C02=-7.771750390D+00                                                      
      C03=4.204607520D+00                                                       
      C04=-2.768070380D+00                                                      
      C05=2.104197070D+00                                                       
      C06=-1.146495880D+00                                                      
      C07=2.231380850D-01                                                       
      C08=1.162503630D-01                                                       
      C09=-8.209005440D-02                                                      
      C010=1.941292390D-02                                                      
      C011=-1.694705760D-03                                                     
      C11=7.086360850D-01                                                       
      C012=-4.311577033D+00                                                     
      C12=1.236794550D+01                                                       
      C13=-1.203890040D+01                                                      
      C14=5.404374220D+00                                                       
      C15=-9.938650430D-01                                                      
      C16=6.275231820D-02                                                       
      C17=-7.747430160D+00                                                      
      C21=-4.298850920D+00                                                      
      C22=4.314305380D+01                                                       
      C23=-1.416193130D+01                                                      
      C24=4.041724590D+00                                                       
      C25=1.555463260D+00                                                       
      C26=-1.665689350D+00                                                      
      C27=3.248811580D-01                                                       
      C28=2.936553250D+01                                                       
      C31=7.948418420D-06                                                       
      C32=8.088597470D+01                                                       
      C33=-8.361533800D+01                                                      
      C34=3.586365170D+01                                                       
      C35=7.518959540D+00                                                       
      C36=-1.261606400D+01                                                      
      C37=1.097174620D+00                                                       
      C38=2.121454920D+00                                                       
      C39=-5.465295660D-01                                                      
      C310=8.328754130D+00                                                      
      C40=2.759717760D-06                                                       
      C41=-5.090739850D-04                                                      
      C50=2.106363320D+02                                                       
      C60=5.528935335D-02                                                       
      C61=-2.336365955D-01                                                      
      C62=3.697071420D-01                                                       
      C63=-2.596415470D-01                                                      
      C64=6.828087013D-02                                                       
      C70=-2.571600553D+02                                                      
      C71=-1.518783715D+02                                                      
      C72=2.220723208D+01                                                       
      C73=-1.802039570D+02                                                      
      C74=2.357096220D+03                                                       
      C75=-1.462335698D+04                                                      
      C76=4.542916630D+04                                                       
      C77=-7.053556432D+04                                                      
      C78=4.381571428D+04                                                       
      D30=-1.717616747D+00                                                      
      D31=3.526389875D+00                                                       
      D32=-2.690899373D+00                                                      
      D33=9.070982605D-01                                                       
      D34=-1.138791156D-01                                                      
      D40=1.301023613D+00                                                       
      D41=-2.642777743D+00                                                      
      D42=1.996765362D+00                                                       
      D43=-6.661557013D-01                                                      
      D44=8.270860589D-02                                                       
      D50=3.426663535D-04                                                       
      D51=-1.236521258D-03                                                      
      D52=1.155018309D-03                                                       
      SIT1=9.626911787D-01                                                      
      SX=SIT-1.0                                                                
      SX1=1.0/SIT                                                               
      SX2=1.0-SIT1                                                              
      Y=(1.0-SIT)/(1.0-SIT1)                                                    
      CH0=CHI4(SIT,BET)                                                         
      N=0                                                                       
    1 CX=1.0/CH0                                                                
      F3=-(C01-CX**2*(C02+CX*(2.0*C03+CX*(3.0*C04+CX*(4.0*C05+CX*(5.0*          
     1   C06+CX*(6.0*C07+CX*(7.0*C08+CX*(8.0*C09+CX*(9.0*C010+CX*10.0*          
     2   C011)))))))))+C012*CX)-(C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*            
     3   C14+CX*(4.0*C15+CX*5.0*C16))))+C17*CX)*SX-(C21-CX**2*(C22+CX*          
     4   (2.0*C23+CX*(3.0*C24+CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))           
     5   +C28*CX)*SX**2-(C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(            
     6   4.0*C35+CX*(5.0*C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))          
     7   +C310*CX)*SX**3+5.0*C41*CH0**(-6)*SIT**(-23)*SX-6.0*CH0**5*            
     8   SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))                     
      F3D=-CX**3*(2.0*C02+CX*(6.0*C03+CX*(12.0*C04+CX*(20.0*C05+CX*(30.0        
     1    *C06+CX*(42.0*C07+CX*(56.0*C08+CX*(72.0*C09+CX*(90.0*C010+CX*         
     2    110.0*C011)))))))))+C012*CX**2-(CX**3*(2.0*C12+CX*(6.0*C13+CX*        
     3    (12.0*C14+CX*(20.0*C15+CX*30.0*C16))))-C17*CX**2)*SX-(CX**3*(         
     4    2.0*C22+CX*(6.0*C23+CX*(12.0*C24+CX*(20.0*C25+CX*(30.0*C26+CX*        
     5    42.0*C27)))))-C28*CX**2)*SX**2-(CX**3*(2.0*C32+CX*(6.0*C33+CX*        
     6    (12.0*C34+CX*(20.0*C35+CX*(30.0*C36+CX*(42.0*C37+CX*(56.0*C38+        
     7    CX*72.0*C39)))))))-C310*CX**2)*SX**3-30.0*C41*CX**7*SIT**(-23)        
     8    *SX-30.0*CH0**4*SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64        
     9    ))))                                                                  
      F4=F3-BET+Y**3*CX**2*(D31+CX*(2.0*D32+CX*(3.0*D33+CX*4.0*D34)))           
     1   +Y**4*CX**2*(D41+CX*(2.0*D42+CX*(3.0*D43+CX*4.0*D44)))-Y**32*(         
     2   D51+2.0*D52*CH0)                                                       
      F4D=F3D-Y**3*CX**3*(2.0*D31+CX*(6.0*D32+CX*(12.0*D33+CX*20.0*D34))        
     1    )-Y**4*CX**3*(2.0*D41+CX*(6.0*D42+CX*(12.0*D43+CX*20.0*D44)))         
     2    -Y**32*2.0*D52                                                        
      CHI=CH0-F4/F4D                                                            
      IF (DABS(CH0/CHI-1.0).LE.1.0D-07) GO TO 2                                 
      N=N+1                                                                     
      IF(N.EQ.50) WRITE(6,10) N,CH0,CHI                                         
   10 FORMAT(10X,'SUB4',I5,2E20.7)                                              
      IF (N.GE.50) GO TO 2                                                      
      CH0=CHI                                                                   
      GO TO 1                                                                   
    2 CONTINUE                                                                  
      CX=1.0/CHI                                                                
      PSIC=C00+C01*CHI+CX*(C02+CX*(C03+CX*(C04+CX*(C05+CX*(C06+CX*(C07+         
     1     CX*(C08+CX*(C09+CX*(C010+CX*C011)))))))))+C012*DLOG(CHI)+(C11        
     2     *CHI+CX*(C12+CX*(C13+CX*(C14+CX*(C15+CX*C16))))+C17*DLOG(CHI)        
     3     )*SX+(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(C25+CX*(C26+CX*C27)        
     4     ))))+C28*DLOG(CHI))*SX**2+(C31*CHI+CX*(C32+CX*(C33+CX*(C34+CX        
     5     *(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)))))))+C310*DLOG(CHI))*         
     6     SX**3+(C40+C41*CHI**(-5))*SIT**(-23)*SX+C50*SIT*DLOG(SIT)+CHI        
     7     **6*SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))+SX*           
     8     (C70+SX*(C71+SX*(C72+SX*(C73+SX*(C74+SX*(C75+SX*(C76+SX*(C77+        
     9     SX*C78))))))))                                                       
      SIG3=-(C11*CHI+CX*(C12+CX*(C13+CX*(C14+CX*(C15+CX*C16))))                 
     1     +C17*DLOG(CHI)+C50)-2.0*(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(        
     2     C25+CX*(C26+CX*C27)))))+C28*DLOG(CHI))*SX-3.0*(C31*CHI+CX*           
     3     (C32+CX*(C33+CX*(C34+CX*(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)         
     4     ))))))+C310*DLOG(CHI))*SX**2+(C40+C41*CHI**(-5))*(22.0*SIT**         
     5     (-23)-23.0*SIT**(-24))-C50*DLOG(SIT)+CHI**6*SX1**3*(2.0*C60+         
     6     SX1*(3.0*C61+SX1*(4.0*C62+SX1*(5.0*C63+SX1*6.0*C64))))-C70-          
     7     SX*(2.0*C71+SX*(3.0*C72+SX*(4.0*C73+SX*(5.0*C74+SX*(6.0*C75+         
     8     SX*(7.0*C76+SX*(8.0*C77+SX*9.0*C78)))))))                            
      BET3=F3                                                                   
      EPS3=PSIC+SIT*SIG3+CHI*BET3                                               
      PSICD1=2.0*(C21*CHI+CX*(C22+CX*(C23+CX*(C24+CX*(C25+CX*(C26+CX*C27        
     1       )))))+C28*DLOG(CHI))+6.0*(C31*CHI+CX*(C32+CX*(C33+CX*(C34+         
     2       CX*(C35+CX*(C36+CX*(C37+CX*(C38+CX*C39)))))))+C310*DLOG(CHI        
     3       ))*SX+(C40+C41*CHI**(-5))*(506.0*SIT-552.0)*SIT**(-25)+C50*        
     4       SX1+CHI**6*SX1**4*(6.0*C60+SX1*(12.0*C61+SX1*(20.0*C62+SX1*        
     5       (30.0*C63+SX1*42.0*C64))))+2.0*C71+SX*(6.0*C72+SX*(12.0*C73        
     6       +SX*(20.0*C74+SX*(30.0*C75+SX*(42.0*C76+SX*(56.0*C77+SX*72.        
     7       0*C78))))))                                                        
      PSICD2=C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*C14+CX*(4.0*C15+CX*5.0*         
     1       C16))))+C17*CX+2.0*(C21-CX**2*(C22+CX*(2.0*C23+CX*(3.0*C24         
     2       +CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))+C28*CX)*SX+3.0*(          
     3       C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(4.0*C35+CX*(5.0*        
     4       C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))+C310*CX)*SX**        
     5       2-5.0*C41*CHI**(-6)*(23.0-22.0*SIT)*SIT**(-24)-6.0*CHI**5*         
     6       SX1**3*(2.0*C60+SX1*(3.0*C61+SX1*(4.0*C62+SX1*(5.0*C63+SX1*        
     7       6.0*C64))))                                                        
      PSICD3=-F3D                                                               
      FAV3=-SIT*PSICD1                                                          
      FAI3=FAV3+SIT*PSICD2**2/PSICD3                                            
      PSID=Y**3*(D30+CX*(D31+CX*(D32+CX*(D33+CX*D34))))+Y**4*(D40+CX*(          
     1     D41+CX*(D42+CX*(D43+CX*D44))))+Y**32*(D50+CHI*(D51+CHI*D52))         
      SIG4=SIG3+(3.0*Y**2*(D30+CX*(D31+CX*(D32+CX*(D33+CX*D34))))+4.0*          
     1     Y**3*(D40+CX*(D41+CX*(D42+CX*(D43+CX*D44))))+32.0*Y**31*(D50+        
     2     CHI*(D51+CHI*D52)))/SX2                                              
      EPS4=EPS3+PSID+(SIG4-SIG3)*SIT+(BET-BET3)*CHI                             
      IF(IZV.EQ.0)RETURN                                                        
      FAV4=SIT*(-PSICD1-1.0/SX2**2*(6.0*Y*(D30+CX*(D31+CX*(D32+CX*(D33+         
     1     CX*D34))))+12.0*Y**2*(D40+CX*(D41+CX*(D42+CX*(D43+CX*D44))))+        
     2     992.0*Y**30*(D50+CHI*(D51+CHI*D52))))                                
      FAI4=FAV4+SIT*(PSICD2+1.0/SX2*(3.0*Y**2*CX**2*(D31+CX*(2.0*D32+CX*        
     1     (3.0*D33+CX*4.0*D34)))+4.0*Y**3*CX**2*(D41+CX*(2.0*D42+CX*(          
     2     3.0*D43+CX*4.0*D44)))                                                
     A     -32.*Y**31*(D51+2.0*D52*CHI)                                         
     B                          ))**2/(PSICD3+Y**3*CX**3*(2.0*D31+CX*(6.        
     3     *D32+CX*(12.0*D33+CX*20.0*D34)))+Y**4*CX**3*(2.0*D41+CX*(6.0         
     4     *D42+CX*(12.0*D43+CX*20.0*D44)))+Y**32*2.0*D52)                      
      CHISQ=CHI*CHI                                                             
      RHRHO4=CHISQ*SIT*(-PSICD1*PSICD3/PSICD2+PSICD2)                           
      RCC4=CHISQ*(PSICD3-PSICD2*PSICD2/PSICD1)                                  
      RKAP4=1./CHI/PSICD3                                                       
      RBETHA=-RKAP4*PSICD2                                                      
      RGAM4=FAI4/FAV4                                                           
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      SUBROUTINE VISCON(P,T,SPVOL,VIS,TCON,NREGON,NSTWAT)                       
C                                                                               
C     VISCOSITY AND THERMAL CONDUCTIVITY OF STEAM AND WATER                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,K8,K9                               
C                                                                               
      K1=-7.691234564D+00                                                       
      K2=-2.608023696D+01                                                       
      K3=-1.681706546D+02                                                       
      K4=6.423285504D+01                                                        
      K5=-1.189646225D+02                                                       
      K6=4.167117320D+00                                                        
      K7=2.097506760D+01                                                        
      K8=1.0D+09                                                                
      K9=6.0D+00                                                                
      S=1.0-(T+273.15)/647.3                                                    
      BEKB=DEXP(S*(K1+S*(K2+S*(K3+S*(K4+S*K5))))/((T+273.15)/647.3*(1.0+        
     1     K6*S+K7*S*S))-S/(K8*S*S+K9))*221.2                                   
      TT=273.16                                                                 
      TC=647.30                                                                 
      P2=1000.0                                                                 
      TK=T+273.15                                                               
      PB=0.980665*P                                                             
      DENS=0.001/SPVOL                                                          
      IF (((NREGON.EQ.1.AND.PB.EQ.BEKB).OR.(NREGON.EQ.4.AND.PB.EQ.BEKB))        
     1   .AND.(T.GE.0.0.AND.T.LE.300.0)) GO TO 1                                
      IF (((NREGON.EQ.1.AND.NSTWAT.EQ.2).OR.(NREGON.EQ.4.AND.NSTWAT.EQ.2        
     1   )).AND.(T.GE.0.0.AND.T.LE.300.0)) GO TO 1                              
      IF ((T.GE.0.0.AND.T.LE.300.0).AND.(PB.GE.BEKB.AND.PB.LE.800.0).AND        
     1   .(NREGON.EQ.1.OR.NREGON.EQ.4)) GO TO 2                                 
      IF ((T.GE.100.0.AND.T.LE.700.0).AND.(P .EQ.1.0)) GO TO 3                  
      IF (((NREGON.EQ.2.AND.PB.EQ.BEKB).OR.(NREGON.EQ.3.AND.PB.EQ.BEKB))        
     1   .AND.(T.GE.100.0.AND.T.LE.300.0)) GO TO 4                              
      IF (((NREGON.EQ.2.AND.NSTWAT.EQ.1).OR.(NREGON.EQ.3.AND.NSTWAT.EQ.1        
     1   )).AND.(T.GE.100.0.AND.T.LE.300.0)) GO TO 4                            
      IF ((T.GE.100.0.AND.T.LE.300.0).AND.(PB.GE.1.0.AND.PB.LE.BEKB))           
     1   GO TO 4                                                                
      IF ((T.GE.375.0.AND.T.LE.700.0).AND.(PB.GE.1.0.AND.PB.LE.800.0))          
     1   GO TO 5                                                                
      VIS=0.0                                                                   
      GO TO 6                                                                   
    1 VIS=241.4*10.0**(247.8*(1.0/(TK-140.0)))                                  
      GO TO 6                                                                   
    2 VIS=241.4*10.0**(247.8*(1.0/(TK-140.0)))*(1.0+(PB-BEKB)*1.046             
     1    *1.0E-6*(TK-305.0))                                                   
      GO TO 6                                                                   
    3 VIS=0.407*T+80.4                                                          
      GO TO 6                                                                   
    4 VIS=0.407*T+80.4-DENS*(1858.0-5.90*T)                                     
      GO TO 6                                                                   
    5 VIS=0.407*T+80.4+DENS*(353.0+DENS*(676.5+DENS*102.1))                     
    6 IF (((NREGON.EQ.1.AND.PB.EQ.BEKB).OR.(NREGON.EQ.4.AND.PB.EQ.BEKB))        
     1   .AND.(T.GE.0.0.AND.T.LE.350.0)) GO TO 7                                
      IF (((NREGON.EQ.1.AND.NSTWAT.EQ.2).OR.(NREGON.EQ.4.AND.NSTWAT.EQ.2        
     1   )).AND.(T.GE.0.0.AND.T.LE.350.0)) GO TO 7                              
      IF ((NREGON.EQ.1.OR.NREGON.EQ.4).AND.((T.GE.0.0.AND.T.LE.350.0)           
     1   .AND.(PB.GE.BEKB.AND.PB.LE.500.0))) GO TO 7                            
      IF ((T.GE.100.0.AND.T.LE.700.0).AND.(P .EQ.1.0)) GO TO 8                  
      IF (PB.LE.175.) TS=BETKP(PB/221.2)                                        
      IF ((T.GE.TS.AND.T.LE.700.0).AND.(PB.GT.1.0.AND.PB.LE.175.0))             
     1   GO TO 9                                                                
      IF ((T.GE.400.0.AND.T.LE.700.0).AND.(PB.GT.175.0.AND.PB.LE.225.0))        
     1   GO TO 9                                                                
      IF ((T.GE.425.0.AND.T.LE.700.0).AND.(PB.GT.225.0.AND.PB.LE.275.0))        
     1   GO TO 9                                                                
      IF ((T.GE.450.0.AND.T.LE.700.0).AND.(PB.GT.275.0.AND.PB.LE.350.0))        
     1   GO TO 9                                                                
      IF ((T.GE.500.0.AND.T.LE.700.0).AND.(PB.GT.350.0.AND.PB.LE.450.0))        
     1   GO TO 9                                                                
      IF ((T.GE.550.0.AND.T.LE.700.0).AND.(PB.GT.450.0.AND.PB.LE.500.0))        
     1   GO TO 9                                                                
      TCON=0.                                                                   
      GO TO 10                                                                  
    7 TX=TK/273.15                                                              
      TCON=-922.47+TX*(2839.5+TX*(-1800.7+TX*(525.77-TX*73.440)))+(PB-          
     1     BEKB)*(-0.94730+TX*(2.5186+TX*(-2.0012+TX*0.51536)))+(PB-BEKB        
     2     )**2*(1.6563E-3+TX*(-3.8929E-3+TX*(2.9323E-3-TX*7.1693E-4)))         
      GO TO 10                                                                  
    8 TCON=17.6+T*(5.87E-2+T*(1.04E-4-T*4.51E-8))                               
      GO TO 10                                                                  
    9 DENS=0.001/SPVOL                                                          
      TCON=17.6+T*(5.87E-2+T*(1.04E-4-T*4.51E-8))+DENS*(103.51+T*(0.4198        
     1     -T*2.771E-5)+DENS*2.1482E14/T**4.2)                                  
   10 RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION BETK(SIT)                                       
C                                                                               
C     REDUCED PRESSURE OF SATURATED STEAM AND WATER AS A FUNCTION OF RED        
C                                                                               
C     TEMPERATURE                                                               
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,K8,K9                               
C                                                                               
      K1=-7.691234564D+00                                                       
      K2=-2.608023696D+01                                                       
      K3=-1.681706546D+02                                                       
      K4=6.423285504D+01                                                        
      K5=-1.189646225D+02                                                       
      K6=4.167117320D+00                                                        
      K7=2.097506760D+01                                                        
      K8=1.0D+09                                                                
      K9=6.0D+00                                                                
      S=1.0-SIT                                                                 
      X=S*(K1+S*(K2+S*(K3+S*(K4+S*K5))))/(SIT*(1.0+K6*S+K7*S*S))-S/(K8*S        
     1  *S+K9)                                                                  
      BETK=DEXP(X)                                                              
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION BETKP(BET)                                      
C                                                                               
C     REDUCED TEMPERATURE OF SATURATED STEAM AND WATER AS A FUNCTION OF         
C                                                                               
C     PRESSURE                                                                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,K8,K9                               
C                                                                               
      K1=-7.691234564D+00                                                       
      K2=-2.608023696D+01                                                       
      K3=-1.681706546D+02                                                       
      K4=6.423285504D+01                                                        
      K5=-1.189646225D+02                                                       
      K6=4.167117320D+00                                                        
      K7=2.097506760D+01                                                        
      K8=1.0D+09                                                                
      K9=6.0D+00                                                                
      N=0                                                                       
      X0=(TSATA(BET*221.2/0.980665)+273.15)/647.3                               
    1 S=1.0-X0                                                                  
      F=DEXP(S*(K1+S*(K2+S*(K3+S*(K4+S*K5))))/(X0*(1.0+K6*S+K7*S*S))            
     1  -S/(K8*S*S+K9))-BET                                                     
      FD=(((-K1-S*(2.0*K2+S*(3.0*K3+S*(4.0*K4+S*5.0*K5))))*X0*(1.0+K6*S         
     1   +K7*S*S)-S*(K1+S*(K2+S*(K3+S*(K4+S*K5))))*(1.0+K6*(1.0-2.0*X0)         
     2   +K7*S*(1.0-3.0*X0)))/(X0*(1.0+K6*S+K7*S*S))**2-(K8*S*S-K9)/            
     3   (K8*S*S+K9)**2)*(F+BET)                                                
      BETKP=X0-F/FD                                                             
      IF (DABS(X0/BETKP-1.0).LE.1.0D-07) GO TO 2                                
      N=N+1                                                                     
      IF (N.GE.50) GO TO 2                                                      
      X0=BETKP                                                                  
      GO TO 1                                                                   
    2 RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION BET3(SIT,CH0)                                   
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      C00=-6.839900000D+00                                                      
      C01=-1.722604200D-02                                                      
      C02=-7.771750390D+00                                                      
      C03=4.204607520D+00                                                       
      C04=-2.768070380D+00                                                      
      C05=2.104197070D+00                                                       
      C06=-1.146495880D+00                                                      
      C07=2.231380850D-01                                                       
      C08=1.162503630D-01                                                       
      C09=-8.209005440D-02                                                      
      C010=1.941292390D-02                                                      
      C011=-1.694705760D-03                                                     
      C012=-4.311577033D+00                                                     
      C11=7.086360850D-01                                                       
      C12=1.236794550D+01                                                       
      C13=-1.203890040D+01                                                      
      C14=5.404374220D+00                                                       
      C15=-9.938650430D-01                                                      
      C16=6.275231820D-02                                                       
      C17=-7.747430160D+00                                                      
      C21=-4.298850920D+00                                                      
      C22=4.314305380D+01                                                       
      C23=-1.416193130D+01                                                      
      C24=4.041724590D+00                                                       
      C25=1.555463260D+00                                                       
      C26=-1.665689350D+00                                                      
      C27=3.248811580D-01                                                       
      C28=2.936553250D+01                                                       
      C31=7.948418420D-06                                                       
      C32=8.088597470D+01                                                       
      C33=-8.361533800D+01                                                      
      C34=3.586365170D+01                                                       
      C35=7.518959540D+00                                                       
      C36=-1.261606400D+01                                                      
      C37=1.097174620D+00                                                       
      C38=2.121454920D+00                                                       
      C39=-5.465295660D-01                                                      
      C310=8.328754130D+00                                                      
      C40=2.759717760D-06                                                       
      C41=-5.090739850D-04                                                      
      C50=2.106363320D+02                                                       
      C60=5.528935335D-02                                                       
      C61=-2.336365955D-01                                                      
      C62=3.697071420D-01                                                       
      C63=-2.596415470D-01                                                      
      C64=6.828087013D-02                                                       
      C70=-2.571600553D+02                                                      
      C71=-1.518783715D+02                                                      
      C72=2.220723208D+01                                                       
      C73=-1.802039570D+02                                                      
      C74=2.357096220D+03                                                       
      C75=-1.462335698D+04                                                      
      C76=4.542916630D+04                                                       
      C77=-7.053556432D+04                                                      
      C78=4.381571428D+04                                                       
      SX=SIT-1.0                                                                
      SX1=1.0/SIT                                                               
    1 CX=1.0/CH0                                                                
      BET3=                                                                     
     O   -(C01-CX**2*(C02+CX*(2.0*C03+CX*(3.0*C04+CX*(4.0*C05+CX*(5.0*          
     1   C06+CX*(6.0*C07+CX*(7.0*C08+CX*(8.0*C09+CX*(9.0*C010+CX*10.0*          
     2   C011)))))))))+C012*CX)-(C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*            
     3   C14+CX*(4.0*C15+CX*5.0*C16))))+C17*CX)*SX-(C21-CX**2*(C22+CX*          
     4   (2.0*C23+CX*(3.0*C24+CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))           
     5   +C28*CX)*SX**2-(C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(            
     6   4.0*C35+CX*(5.0*C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))          
     7   +C310*CX)*SX**3+5.0*C41*CH0**(-6)*SIT**(-23)*SX-6.0*CH0**5*            
     8   SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))                     
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION BET4(SIT,CH0)                                   
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      C00=-6.839900000D+00                                                      
      C01=-1.722604200D-02                                                      
      C02=-7.771750390D+00                                                      
      C03=4.204607520D+00                                                       
      C04=-2.768070380D+00                                                      
      C05=2.104197070D+00                                                       
      C06=-1.146495880D+00                                                      
      C07=2.231380850D-01                                                       
      C08=1.162503630D-01                                                       
      C09=-8.209005440D-02                                                      
      C010=1.941292390D-02                                                      
      C011=-1.694705760D-03                                                     
      C11=7.086360850D-01                                                       
      C012=-4.311577033D+00                                                     
      C12=1.236794550D+01                                                       
      C13=-1.203890040D+01                                                      
      C14=5.404374220D+00                                                       
      C15=-9.938650430D-01                                                      
      C16=6.275231820D-02                                                       
      C17=-7.747430160D+00                                                      
      C21=-4.298850920D+00                                                      
      C22=4.314305380D+01                                                       
      C23=-1.416193130D+01                                                      
      C24=4.041724590D+00                                                       
      C25=1.555463260D+00                                                       
      C26=-1.665689350D+00                                                      
      C27=3.248811580D-01                                                       
      C28=2.936553250D+01                                                       
      C31=7.948418420D-06                                                       
      C32=8.088597470D+01                                                       
      C33=-8.361533800D+01                                                      
      C34=3.586365170D+01                                                       
      C35=7.518959540D+00                                                       
      C36=-1.261606400D+01                                                      
      C37=1.097174620D+00                                                       
      C38=2.121454920D+00                                                       
      C39=-5.465295660D-01                                                      
      C310=8.328754130D+00                                                      
      C40=2.759717760D-06                                                       
      C41=-5.090739850D-04                                                      
      C50=2.106363320D+02                                                       
      C60=5.528935335D-02                                                       
      C61=-2.336365955D-01                                                      
      C62=3.697071420D-01                                                       
      C63=-2.596415470D-01                                                      
      C64=6.828087013D-02                                                       
      C70=-2.571600553D+02                                                      
      C71=-1.518783715D+02                                                      
      C72=2.220723208D+01                                                       
      C73=-1.802039570D+02                                                      
      C74=2.357096220D+03                                                       
      C75=-1.462335698D+04                                                      
      C76=4.542916630D+04                                                       
      C77=-7.053556432D+04                                                      
      C78=4.381571428D+04                                                       
      D30=-1.717616747D+00                                                      
      D31=3.526389875D+00                                                       
      D32=-2.690899373D+00                                                      
      D33=9.070982605D-01                                                       
      D34=-1.138791156D-01                                                      
      D40=1.301023613D+00                                                       
      D41=-2.642777743D+00                                                      
      D42=1.996765362D+00                                                       
      D43=-6.661557013D-01                                                      
      D44=8.270860589D-02                                                       
      D50=3.426663535D-04                                                       
      D51=-1.236521258D-03                                                      
      D52=1.155018309D-03                                                       
      SIT1=9.626911787D-01                                                      
      SX=SIT-1.0                                                                
      SX1=1.0/SIT                                                               
      SX2=1.0-SIT1                                                              
      Y=(1.0-SIT)/(1.0-SIT1)                                                    
    1 CX=1.0/CH0                                                                
      F3=-(C01-CX**2*(C02+CX*(2.0*C03+CX*(3.0*C04+CX*(4.0*C05+CX*(5.0*          
     1   C06+CX*(6.0*C07+CX*(7.0*C08+CX*(8.0*C09+CX*(9.0*C010+CX*10.0*          
     2   C011)))))))))+C012*CX)-(C11-CX**2*(C12+CX*(2.0*C13+CX*(3.0*            
     3   C14+CX*(4.0*C15+CX*5.0*C16))))+C17*CX)*SX-(C21-CX**2*(C22+CX*          
     4   (2.0*C23+CX*(3.0*C24+CX*(4.0*C25+CX*(5.0*C26+CX*6.0*C27)))))           
     5   +C28*CX)*SX**2-(C31-CX**2*(C32+CX*(2.0*C33+CX*(3.0*C34+CX*(            
     6   4.0*C35+CX*(5.0*C36+CX*(6.0*C37+CX*(7.0*C38+CX*8.0*C39)))))))          
     7   +C310*CX)*SX**3+5.0*C41*CH0**(-6)*SIT**(-23)*SX-6.0*CH0**5*            
     8   SX1**2*(C60+SX1*(C61+SX1*(C62+SX1*(C63+SX1*C64))))                     
      BET4=  F3+Y**3*CX**2*(D31+CX*(2.0*D32+CX*(3.0*D33+CX*4.0*D34)))           
     1   +Y**4*CX**2*(D41+CX*(2.0*D42+CX*(3.0*D43+CX*4.0*D44)))-Y**32*(         
     2   D51+2.0*D52*CH0)                                                       
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION CHI2(SIT,BET)                                   
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION B(8,3),BS(8,2)                                                  
C                                                                               
      B90=1.936587558D+02                                                       
      B91=-1.388522425D+03                                                      
      B92=4.126607219D+03                                                       
      B93=-6.508211677D+03                                                      
      B94=5.745984054D+03                                                       
      B95=-2.693088365D+03                                                      
      B96=5.235718623D+02                                                       
      BB=7.633333333D-01                                                        
      B(1,1)=6.670375918D-02                                                    
      B(1,2)=1.388983801D+00                                                    
      B(1,3)=0.                                                                 
      B(2,1)=8.390104328D-02                                                    
      B(2,2)=2.614670893D-02                                                    
      B(2,3)=-3.373439453D-02                                                   
      B(3,1)=4.520918904D-01                                                    
      B(3,2)=1.069036614D-01                                                    
      B(3,3)=0.                                                                 
      B(4,1)=-5.975336707D-01                                                   
      B(4,2)=-8.847535804D-02                                                   
      B(4,3)=0.                                                                 
      B(5,1)=5.958051609D-01                                                    
      B(5,2)=-5.159303373D-01                                                   
      B(5,3)=2.075021122D-01                                                    
      B(6,1)=1.190610271D-01                                                    
      B(6,2)=-9.867174132D-02                                                   
      B(6,3)=0.                                                                 
      B(7,1)=1.683998803D-01                                                    
      B(7,2)=-5.809438001D-02                                                   
      B(7,3)=0.                                                                 
      B(8,1)=6.552390126D-03                                                    
      B(8,2)=5.710218649D-04                                                    
      B(8,3)=0.                                                                 
      BS(1,1)=0.                                                                
      BS(1,2)=0.                                                                
      BS(2,1)=0.                                                                
      BS(2,2)=0.                                                                
      BS(3,1)=0.                                                                
      BS(3,2)=0.                                                                
      BS(4,1)=0.                                                                
      BS(4,2)=0.                                                                
      BS(5,1)=0.                                                                
      BS(5,2)=0.                                                                
      BS(6,1)=4.006073948D-01                                                   
      BS(6,2)=0.                                                                
      BS(7,1)=8.636081627D-02                                                   
      BS(7,2)=0.                                                                
      BS(8,1)=-8.532322921D-01                                                  
      BS(8,2)=3.460208861D-01                                                   
      TC1=647.3                                                                 
      PC1=22120000.0                                                            
      VC1=0.00317                                                               
      R1=461.51                                                                 
      AI1=R1*TC1/(PC1*VC1)                                                      
      X=DEXP(BB*(1.0-SIT))                                                      
      CALL SUBBET(SIT,BETL,BETLD,BETLDD)                                        
      S03=0.0                                                                   
      S04=0.0                                                                   
      DO 4 I=1,8                                                                
      BX=0.0                                                                    
      NNI=NN(I)                                                                 
      NLI=NL(I)                                                                 
      DO 2 J=1,NNI                                                              
      SBX=0.0                                                                   
      NZIJ=NZ(I,J)                                                              
      DO 1 K=1,NLI                                                              
      NXIK=NX(I,K)                                                              
      SBX=SBX+BS(I,K)*X**NXIK                                                   
    1 CONTINUE                                                                  
      BX=BX+B(I,J)*X**NZIJ                                                      
    2 CONTINUE                                                                  
      IF (I.GE.6) GO TO 3                                                       
      S03=S03+FLOAT(I)*BET**(I-1)*BX                                            
      GO TO 4                                                                   
    3 CONTINUE                                                                  
      S04=S04+FLOAT(I-2)*BET**(1-I)*BX/(BET**(2-I)+SBX)**2                      
    4 CONTINUE                                                                  
      CHI2=AI1*SIT/BET-S03-S04+11.0*(BET/BETL)**10*(B90+X*(B91+X*(B92+X*        
     1     (B93+X*(B94+X*(B95+X*B96))))))                                       
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION CHI3(SIT,BET)                                   
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DOUBLE PRECISION L0,L1,L2                                                 
C                                                                               
      L0=15.74373327                                                            
      L1=-34.17061978                                                           
      L2=19.31380707                                                            
      SIT1=0.9626911787                                                         
      BET1=0.7475191707                                                         
      BET2=4.52079566                                                           
      BTK=BETK(SIT)                                                             
      L0=L0-BET                                                                 
      IF ((BET.GT.BET1.AND.BET.LE.1.0).OR.(SIT.GT.SIT1.AND.SIT.LE.1.0))         
     1  GO TO 1                                                                 
      IF (BET.GT.1.0.AND.SIT.GT.SIT1) GO TO 2                                   
    1 CH0=CHI2(SIT,BET)                                                         
      GO TO 6                                                                   
    2 CH0=0.001368/0.00317                                                      
      SITL=(-L1+DSQRT(L1**2-4.0*L0*L2))/(2.0*L2)                                
    3 CH2=CHI2(SITL,BET)                                                        
      IF (CH2.GT.0.) GO TO 4                                                    
      SITL=SITL+0.1                                                             
      GO TO 3                                                                   
    4 BET0=BET3(SIT,CH0)-BET                                                    
    5 CH1=CH0+0.0001/0.00317                                                    
      BET1=BET3(SIT,CH1)-BET                                                    
      IF (BET0*BET1.LE.0.) GO TO 6                                              
      IF (CH1.GE.CH2) GO TO 6                                                   
C                                                                               
      CH0=CH1                                                                   
      BET0=BET1                                                                 
      GO TO 5                                                                   
    6 CHI3=CH0                                                                  
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      DOUBLE PRECISION FUNCTION CHI4(SIT,BET)                                   
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      SIT1=0.96269117                                                           
      CH0=0.001318/0.00317                                                      
      CH2=0.002448/0.00317                                                      
      BET0=BET4(SIT,CH0)-BET                                                    
    1 CH1=CH0+0.0001/0.00317                                                    
      BET1=BET4(SIT,CH1)-BET                                                    
      IF (BET0*BET1.LT.0.) GO TO 2                                              
      IF (CH1.GT.CH2) GO TO 2                                                   
C                                                                               
      CH0=CH1                                                                   
      BET0=BET1                                                                 
      GO TO 1                                                                   
    2 CHI4=CH0                                                                  
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      FUNCTION NL(I)                                                            
C                                                                               
      NL=0                                                                      
      IF (I.EQ.6.OR.I.EQ.7) NL=1                                                
      IF (I.EQ.8) NL=2                                                          
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      FUNCTION NN(I)                                                            
C                                                                               
      NN=2                                                                      
      IF (I.EQ.2.OR.I.EQ.5) NN=3                                                
      RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      FUNCTION NX(I,K)                                                          
C                                                                               
      NX=0                                                                      
      IF (K-1) 1,1,2                                                            
    1 IF (I.EQ.6) NX=14                                                         
      IF (I.EQ.7) NX=19                                                         
      IF (I.EQ.8) NX=54                                                         
      GO TO 3                                                                   
    2 IF (I.EQ.8) NX=27                                                         
    3 RETURN                                                                    
      END                                                                       
C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
C                                                                               
      FUNCTION NZ(I,J)                                                          
C                                                                               
      NZ=0                                                                      
      IF (J-2) 1,2,3                                                            
    1 IF (I.EQ.1) NZ=13                                                         
      IF (I.EQ.2.OR.I.EQ.3) NZ=18                                               
      IF (I.EQ.4) NZ=25                                                         
      IF (I.EQ.5) NZ=32                                                         
      IF (I.EQ.6) NZ=12                                                         
      IF (I.EQ.7.OR.I.EQ.8) NZ=24                                               
      GO TO 4                                                                   
    2 IF (I.EQ.1) NZ=3                                                          
      IF (I.EQ.2) NZ=2                                                          
      IF (I.EQ.3) NZ=10                                                         
      IF (I.EQ.4.OR.I.EQ.8) NZ=14                                               
      IF (I.EQ.5) NZ=28                                                         
      IF (I.EQ.6) NZ=11                                                         
      IF (I.EQ.7) NZ=18                                                         
      GO TO 4                                                                   
    3 IF (I.EQ.2) NZ=1                                                          
      IF (I.EQ.5) NZ=24                                                         
    4 RETURN                                                                    
      END                                                                       
