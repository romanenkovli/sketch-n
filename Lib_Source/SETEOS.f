C*DK SETEOS                                                             
      SUBROUTINE SETEOS                                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     THIS SOUBOUTINE SETS THE EQUATION OF STATE CONSTANTS              
C                                                                       
C-----------------------------COMMON BLOCKS-----------------------------
C                                                                       
C*CA TSATCN                                                             
C*CC TSATCN                                                             
      INCLUDE 'TSATCN'                                                  
C                                                                       
C------------------------------NOMENCLATURE-----------------------------
C                                                                       
C     AIRMOL  = MOLECULAR WEIGHT OF AIR                                 
C     CEOSLP( 1) = 1ST COEFF OF SAT VAP TEMP FUNCTION                   
C     CEOSLP( 2) = 2ND COEFF OF SAT VAP TEMP FUNCTION                   
C     CEOSLP( 3) = 3RD COEFF OF SAT VAP TEMP FUNCTION                   
C     CEOSLP( 4) = VAPOR SPECIFIC HEAT AT CONSTANT VOLUME               
C     CEOSLP( 5) = REFERENCE TEMPERATURE                                
C     CEOSLP( 6) = VAPOR SPECIFIC INTERNAL ENERGY AT REF. TEMP.         
C     CEOSLP( 7) = LIQUID SPECIFIC HEAT AT CONSTANT VOLUME              
C     CEOSLP( 8) = LIQUID SPECIFIC INTERNAL ENERGY AT REF. TEMP.        
C     CEOSLP( 9) = MICROSCOPIC DENSITY OF LIQUID                        
C     CEOSLP(10) = HEAT OF EVAPORATION AT REFERENCE TEMPERATURE         
C     CEOSLP(11) = REFERENCE PRESSURE                                   
C     CEOSLP(12) = GAS CONSTANT FOR VAPOR                               
C     CEOSLP(13) = AIR THERMAL CONDUCTIVITY                             
C     CEOSLP(14) = LIQUID THERMAL CONDUCTIVITY                          
C     CEOSLP(15) = SQUARE OF RECIPROCAL SOUND SPEED FOR LIQUID          
C     CEOSLP(16) = VAPOR GAMMA, RATIO OF SPECIFIC HEATS                 
C     CEOSLP(17) = AIR SPECIFIC HEAT AT CONSTANT VOLUME                 
C     CEOSLP(18) = AIR GAMMA, RATIO OF SPECIFIC HEATS                   
C     CEOSLP(19) = RATIO OF GAS COMPONENT 2 TO COMPONENT 1              
C     CEOSLP(20) = UPPER LIMIT ON PRESSURE FOR LOW-PRESSURE PROPERTIES  
C     CEOSLP(21) = UPPER LIMIT ON TEMP. FOR LOW-PRESSURE PROPERTIES     
C     CEOSLP(22) = AIR SPECIFIC HEAT AT CONSTANT PRESSURE               
C     CEOSLP(23) = VAPOR SPECIFIC HEAT AT CONSTANT PRESSURE             
C     CEOSLP(24) = LIQUID SPECIFIC HEAT AT CONSTANT PRESSURE            
C     CEOSLP(25) = GAS CONSTANT FOR AIR                                 
C     CEOSLP(26) = LIQUID SPECIFIC ENTHALPY AT REFERENCE TEMPERATURE    
C     CEOSLP(27) = VAPOR SPECIFIC ENTHALPY AT REFERENCE TEMPERATURE     
C     CEOSLP(28) = (VAP. GAS CONST. - AIR GAS CONST.) / VAP. GAS CONST. 
C     CEOSLP(29) = DATUM TEMPERATURE FOR THE ENTHALPY OF LIQUID         
C     CEOSLP(30) = MINIMUM ALLOWABLE PRESSURE                           
C     CEOSLP(31) = MAXIMUM ALLOWABLE PRESSURE                           
C     CEOSLP(32) = MINIMUM ALLOWABLE LIQUID TEMPERATURE                 
C     CEOSLP(33) = MAXIMUM ALLOWABLE LIQUID TEMPERATURE                 
C     CEOSLP(34) = MINIMUM ALLOWABLE VAPOR TEMPERATURE                  
C     CEOSLP(35) = MAXIMUM ALLOWABLE VAPOR TEMPERATURE                  
C     CEOSLP(36) = MINIMUM SATURATION PRESSURE                          
C                  (CORRESPONDING TO MINIMUM LIQUID TEMPERATURE)        
C     CEOSLP(37) = CRITICAL PRESSURE                                    
C     CEOSLP(38) = CRITICAL TEMPERATURE                                 
C     CEOSLP(39) = LOWER LIMIT ON SATURATION PRESSURE FOR THE           
C                  HIGH-PRESSURE EXPRESSION                             
C     CEOSLP(40) = LOWER LIMIT ON SATURATION TEMPERATURE FOR THE        
C                  HIGH-PRESSURE EXPRESSION                             
C     GASCON     = UNIVERSAL GAS CONSTANT                               
C     IGAS       = 1 FOR AIR, 2 FOR HYDROGEN                            
C     VAPMOL     = MOLECULAR WEIGHT OF VAPOR                            
C                                                                       
C----------------------------EXECUTION LOGIC----------------------------
C                                                                       
      IGSET     = 0                                                     
      AEOS14    = 1.0D-05                                               
      CEOS1     = 117.8D0                                               
      CEOS2     = 0.223D0                                               
      CEOS3     = 255.2D0                                               
      CEOSLP(1) = -2263.0D0                                             
      CEOSLP(2) = 0.434D0                                               
      CEOSLP(3) = -6.064D0                                              
      GO TO 10                                                          
C                                                                       
C     ENTRY POINT SETIG FOR RESETTING THE INERT GAS OPTION              
C                                                                       
      ENTRY SETIG                                                       
      IGSET = 1                                                         
   10 CONTINUE                                                          
      GASCON     = 6.022169D+26*1.380622D-23                            
      AIRMOL     = 28.967D0                                             
      CEOSLP(22) = 1004.832D0                                           
      CEOSLP(13) = 0.0228D0                                             
      IF( IGAS.EQ.1 ) GO TO 100                                         
      AIRMOL     = 2.016D0                                              
      CEOSLP(22) = 14533.2D0                                            
      CEOSLP(13) = 0.3094D0                                             
  100 CONTINUE                                                          
      CEOSLP(25) = GASCON/AIRMOL                                        
      CEOSLP(17) = CEOSLP(22)-CEOSLP(25)                                
      CEOSLP(18) = CEOSLP(22)/CEOSLP(17)                                
CC    CEOSLP(28)=(CEOSLP(12)-CEOSLP(25))/CEOSLP(12)                     
      IF( IGSET.EQ.0 ) GO TO 200                                        
      RETURN                                                            
C                                                                       
C     CONTINUE THE REMAINING INITIALIZATION IF ENTRY POINT  IS NOT SETIG
C                                                                       
  200 CONTINUE                                                          
      CEOSLP(5)  = 273.15D0                                             
      CEOSLP(29) = 273.16D0                                             
      CEOSLP(30) = 1.0D0                                                
      CEOSLP(31) = 450.0D+5                                             
      CEOSLP(32) = CEOSLP(5)                                            
      CEOSLP(33) = 713.94025779311D0                                    
      CEOSLP(34) = CEOSLP(5)                                            
      CEOSLP(35) = 3000.0D0                                             
      CEOSLP(36) = 610.8D0                                              
      CEOSLP(37) = 221.2D+5                                             
      CEOSLP(38) = 647.3D0                                              
      CEOSLP(39) = 139.69971285053D+5                                   
      CEOSLP(40) = 609.62462615967D0                                    
      VAPMOL     = 18.016D0                                             
      CEOSLP(12) = GASCON/VAPMOL                                        
      CEOSLP(16) = 1.3D0                                                
      CEOSLP(4)  = CEOSLP(12)/(CEOSLP(16)-1.0D0)                        
      CEOSLP(23) = CEOSLP(16)*CEOSLP(4)                                 
CC    CEOSLP(28)=(CEOSLP(12)-CEOSLP(25))/CEOSLP(12)                     
      CEOSLP(24) = 4186.800D0                                           
      CEOSLP(14) = 0.65141D0                                            
      CEOSLP(7)  = CEOSLP(24)                                           
      CEOSLP(9)  = 990.0D0                                              
      CEOSLP(15) = 0.0D0                                                
CC    CEOSLP(19)=NOT USED ANY MORE                                      
      CEOSLP(20) = 9.056466D+4                                          
      CEOSLP(21) = 370.4251D0                                           
      CEOSLP(11) = 100000.0D0                                           
      CEOSLP(10) = HEV(CEOSLP(5))                                       
      CEOSLP(8)  = -611.2D0*0.0010002D0+CEOSLP(7)*(CEOSLP(5)-CEOSLP(29))
      CEOSLP(26) = CEOSLP(24)*(CEOSLP(5)-CEOSLP(29))                    
      CEOSLP(27) = CEOSLP(26)+CEOSLP(10)                                
      CEOSLP(6)  = CEOSLP(27)-CEOSLP(12)*CEOSLP(5)                      
      RETURN                                                            
      END                                                               
