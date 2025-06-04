C
      SUBROUTINE PF1STEAM( P, TL, ROL, CPL, VISL, CL )
C
      PARAMETER ( IDRMAX=100 )
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'IOUNITS'
C
      DIMENSION DR(IDRMAX), DRS(IDRMAX)
C
      DO 10 I = 1, IDRMAX
        DR (I) = 0.0D0
        DRS(I) = 0.0D0
   10 CONTINUE
C
      ISAT = 0
      TV   = TL
C
      ILAST = 0
C
      PMIN   = 1.0D4
      PMAX   = 2.212D7
      TLMIN  = 273.15D0
      TLMAX  = 647.8D0
      TVMIN  = 273.15D0
      TVMAX  = 1000.0D0
C
      JSTART = 1
      NCELLS = 1
      IOP    = 0
      NC     = 1
C
      PA     = 0.0D0
C
      IMAX   = 1

      DO 100 I = 1, IMAX
C
        EL  = 0.0D0
        EV  = 0.0D0
C
        CALL SETEOS
C
        CALL THERMO( P      , EL     , EV     , TL     , TV     ,
     1               TSAT   , ROL    , ROV    , PA     , ROA    ,
     2               TSSN   , EA     , DR     , IOP    , JSTART ,
     3               NCELLS  )
C
        IF( ISAT.EQ.0 ) GO TO 50
C
        TL   = TSAT
        TV   = TSAT
C
        ELS  = 0.0D0
        EVS  = 0.0D0
C
        CALL THERMO( P      , ELS    , EVS    , TL     , TV     ,
     1               TSAT   , ROLS   , ROVS   , PA     , ROAS   ,
     2               TSSN   , EAS    , DRS    , IOP    , JSTART ,
     3               NCELLS  )
C
        EL  = ELS
        EV  = EVS
        EA  = EAS
C
        ROL = ROLS
        ROV = ROVS
        ROA = ROAS
C
        DO 30 KK = 1, IDRMAX
          DR(KK) = DRS(KK)
   30   CONTINUE
C
   50   CONTINUE
C
        CALL FPROP  (P      , EL     , EV     , ROL    , ROV    ,
     1               TL     , TV     , TSAT   , HFG    , CPL    ,
     2               CPV    , VISL   , VISV   , CL     , CV     ,
     2               SIG    , NC     , PA     , ROA    , TSSN   ,
     3               EA     , DR     )
C
  100 CONTINUE
C
      RETURN
      END
