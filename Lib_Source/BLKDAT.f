C*DK BLKDAT                                                             
      BLOCK DATA
C
      IMPLICIT REAL*8(A-H,O-Z)
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
C*CA IOUNITS                                                            
C*CC IOUNITS                                                            
      INCLUDE 'IOUNITS'                                                 
C                                                                       
      DATA IEOS  /  0 /
      DATA NTHM  /  0 /
      DATA IFTP  /  0 /
      DATA IGAS  /  0 /
      DATA ILIQ  /  0 /
      DATA IMOUT / 60 /
C                                                                       
      END                                                               
