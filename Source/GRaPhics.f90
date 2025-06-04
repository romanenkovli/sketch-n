      SUBROUTINE GRaPhics_Write_INI_Data
!=====================================================================!
! Output Parameters into .GRF file                                    !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

      OPEN(io_unit,file='Output/SKETCH.grf', status='unknown', &
       FORM= 'UNFORMATTED')
     
      CALL GRaPhics_Output_Parameters(io_unit)
      CALL GRaPhics_Output_GeoMTData(io_unit)
      IF(FILE_CD.NE."".AND.XS_model.EQ."TABLE") THEN
!        write(*,*) "writing steady-state distributions"
          CALL GraPhics_Output_ST_Dist(io_unit)
      END IF

      CLOSE(io_unit)

      RETURN
      END


      SUBROUTINE GRaPhics_Write_Distr_Data(time)
!=====================================================================!
! Output Parameters into .GRF file                                    !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      REAL time

      OPEN(io_unit,file='Output/SKETCH.grf', status='OLD', &
       FORM= 'UNFORMATTED', ACCESS ='APPEND')
        CALL GRaPhics_Output_TRN_Distr(io_unit, time)

      CLOSE(io_unit)

      RETURN
      END




      SUBROUTINE GRaPhics_Output_Parameters(unit)
!=====================================================================!
! Output Parameters into .GRF file                                    !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER unit
! Local:
      INTEGER nl, n_st_dist, n
      LOGICAL flag_dist_reactor 

      CALL OUTput_Write_Header_unformatted(Unit)

      WRITE(unit) N_LINE_PROBLEM_TITLE 
      DO nl = 1, N_LINE_PROBLEM_TITLE
         WRITE(unit) PROBLEM_TITLE(nl)
      END DO

! Geometry Description
      WRITE(unit) GMT_CRD_TYPE

! GeoMTry Module                                                      !
      WRITE(unit) N_POLY, NH, NZR, NZ, NXR, NYR, NX, NY, &
               NCHM, NDD,N_BUNDLE_TYPE
! TH_Model                                                             !
      WRITE(unit) NP_Reactor_Core, NZR_Core

! XS Module (Cross Section & Neutron Kinetics Constant)                !
      WRITE(unit) NNODE, NG, MD, N_FEEDBACK

! CRD control rod Module                                               !
      WRITE(unit) NN_CRod, NN_CRod_Comp, NN_CRod_El, NN_CRod_Type, &
       NN_Crod_Bundle

! Fuel RoD (FRD) Heat Conduction Model                                 !
      WRITE(unit) NN_FRD_FUEL, NN_FRD_CLAD, NN_FRD_TOTAL

! PVM Interface Module for TRAC                                        !
      WRITE(unit) NN_RT_HC_TRAC, NN_RT_FD_TRAC, NN_Z_HC_TRAC, &
       NN_Z_FD_TRAC 

! REACTOR TYPE
      WRITE(unit) REACTOR_TYPE

! Number of Steady-State Distributions from File_CD
      WRITE(unit) n_out_st_scal
! N_out_st_scal == 0
!      DO n = 1, n_out_st_scal
!         WRITE(unit) name_st_scal(n)
!      END DO

      WRITE(unit) n_out_st_vec
! N_out_st_vec == 0
!      DO n = 1, n_out_st_vec
!         WRITE(unit) dim_st_vec(n)
!         WRITE(unit) name_st_vec(n)
!      END DO

      IF(File_CD.NE."".AND.XS_model.EQ."TABLE") THEN
          n_st_dist = N_OUT_RG_DIST
      ELSE
          n_st_dist = 0
      END IF

      WRITE(unit) n_st_dist

      DO n = 1,  n_st_dist
         WRITE(unit) name_st_dist(n)
      END DO

! ALL steady-state sdistributions given for reactor core   
      flag_dist_reactor = .False. 
      DO n = 1, n_st_dist
         WRITE(unit) flag_dist_reactor
      END DO

! Number of Internal T/H Distributions
      WRITE(unit) n_out_tr_scal

      DO n = 1, n_out_tr_scal
         WRITE(unit) name_tr_scal(n)
      END DO

      WRITE(unit) n_out_tr_vec
      DO n = 1, n_out_tr_vec
         WRITE(unit) dim_tr_vec(n)
         WRITE(unit) name_tr_vec(n)
      END DO

! N_TR_DIST already computed
!      n_tr_dist = 1 + N_FEEDBACK 

!      IF(TH_Model.EQ."Internal".OR.TH_Model.EQ."SKAZKA") THEN
!          n_tr_dist = n_tr_dist +  N_OUT_TH_DIST 
!      END IF

!      IF(Problem_Type.EQ."Burnup") THEN
!         N_TR_DIST  = N_TR_DIST + 1 + N_ISOTOPE
! + burnup + isotopes
!      END IF 
!      n_tr_dist = n_tr_dist + NG
! neutron flux

      WRITE(unit) n_tr_dist


      DO n = 1,  n_tr_dist
         WRITE(unit) name_tr_dist(n)
      END DO

! the first 1 + N_FEEDBACK + N_OUT_TH_DIST distributions are for the core 
      flag_dist_reactor = .False. 
      DO  n = 1, n_tr_dist-NG
         WRITE(unit) flag_dist_reactor
      END DO
! neutron flux distribution is given for the reactor
      flag_dist_reactor = .True. 
      DO  n = n_tr_dist-NG+1, n_tr_dist
         WRITE(unit) flag_dist_reactor
      END DO

      RETURN
      END


      SUBROUTINE GRaPhics_Output_GeoMTData(unit)
!=====================================================================!
! Output Geometry Model Data into .GRF file                           !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER unit

! Local: 
      INTEGER nlx, nly, nlz
      INTEGER np, ns
      INTEGER ib, ie, ir
      INTEGER n1, nd
      INTEGER nc

! Output:  
! Core Numbering
      WRITE(unit) (Numb_Reactor_Core(nc), nc = 1, NP_Reactor_Core)

! Channels Coordinates:
      WRITE(unit) ((N_Coord(np, nd), np=1, N_POLY), nd=1,2)

! Reactor Geometry
      WRITE(unit) ((npoly(nly,nlx), nly=1, NYR), nlx=1, NXR),&
                 (Nxr_B_Reactor(nly), nly=1, NYR),&
                 (Nxr_E_Reactor(nly), nly=1, NYR),&
                 Nyr_B_Reactor, Nyr_E_Reactor, &
                 Nxr_Max_Reactor, Nxr_B_Min_Reactor


! Reactor Core Geometry
      WRITE(unit) ((Index_Core(nly,nlx), nly=1, NYR), nlx=1, NXR),&
                 (Nxr_B_Core(nly), nly=1, NYR), &
                 (Nxr_E_Core(nly), nly=1, NYR),&
                 Nyr_B_Core, Nyr_E_Core, &
                 Nxr_Max_Core, Nxr_B_Min_Core,&
                 NZR_Core_Beg, NZR_Core_End


! Spatial Mesh
      WRITE(unit) (hx(nlx), nlx = 1,NXR),&
                 (hy(nly), nly = 1, NYR),&
                 (hz(nlz), nlz = 1, NZR)
! Fine Axial Mesh
      WRITE(unit) (npz(ns), ns = 1, NZR),&
                 (hzt(n1), n1 =1, NZ)

! Reactor & Core Volumes & Nodal Volumes
      WRITE(unit) v_reactor, v_core, &
        ((vol_ass(np, ns), np=1, N_POLY), ns=1, NZR)

! Control RoD Description
      WRITE(unit) &
     ((Mat_Com_Rod(ir,ie), ir = 1,NN_CRod), ie =1,NN_CRod_El), &
     ((nrods(ir, ib), ir = 1,NN_CRod), ib = 1, NN_Crod_Bundle),&
     (h_rod_el(ie), ie = 1, NN_CRod_El)

      RETURN
      END

      SUBROUTINE GRaPhics_Output_TRN_Distr(unit, time)

!=====================================================================!
! Output Geometry Model Data into .GRF file                           !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER unit
      REAL time

! Local: 
      INTEGER ir, np, ns, i, n
! Local
      REAL sec_to_ms
      PARAMETER (sec_to_ms = 1.E+3)

! TIME 
      WRITE(unit) time

! SCALAR: 
!    Reactivity
!      write(*,*) 'bet, react =', bet, react
!      write(*,*) 'p_total, react/bet =', p_total, REAL(react/bet)
      WRITE(unit) p_total, REAL(react/bet), &
!                     dt_sketch*SEC_to_MS, dt_trac*SEC_to_MS,&
                     dt_sketch, dt_trac,&
                     data_scalar_skazka(1), data_scalar_skazka(2),&
                     k_ef, pow_ax_offset, conc_ax_offset_i, &
                     conc_ax_offset_xe, p_av_reactor

! ARRAYS: 
!    Control RoD Positions
      WRITE(unit) (zrods(ir), ir = 1, NN_CRod)
!      WRITE(*,*) '(zrods(ir), ir = 1, NN_CRod)'
!      WRITE(*,*) (zrods(ir), ir = 1, NN_CRod)

! TRANSIENT DISTRIBUTION:
!    Power 
      WRITE(unit) ((p_col(np,ns), np=1, N_POLY),ns=1,NZR)
      
!    Feedback 
      DO i = 1, N_FEEDBACK
        WRITE(unit) ((fdback_col(np,ns, i), np=1, N_POLY),ns=1,NZR)
      END DO
!    T/H Distributions
      IF(TH_Model.EQ."Internal".OR.TH_Model.EQ."SKAZKA") THEN
         DO i = 1, N_OUT_TH_DIST 
!            WRITE(*,*) 'WRITING T/H Distributions i =', i
          WRITE(unit) ((dist_th_col(np,ns,i),np=1,N_POLY),ns=1,NZR)
         END DO
      END IF

      IF(Problem_Type.EQ."Burnup") THEN

         WRITE(unit) ((brn_col(np,ns),np=1,N_POLY),ns=1,NZR)

         CALL Xe_Sm_Compute_Distribution 
         DO i = 1, N_ISOTOPE
            WRITE(unit)((conc_isotope_col(np,ns,i),np=1,N_POLY),ns=1,NZR)
         END DO
      END IF 

!    Neutron Flux 
      DO n = 1, NG
         WRITE(unit) ((dist_flux(np,ns,n), np=1, N_POLY),ns=1,NZR)
      END DO

      RETURN 
      END


      SUBROUTINE GRaPhics_Output_st_Dist(unit)

!=====================================================================!
! Output Geometry Model Data into .GRF file                           !
!                  Vyachreslav Zimin (c) April 1 1998                 !
!              vgzimin@mail.ru                  !
!=====================================================================!
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER unit

! Local: 
      INTEGER np, ns, i

! Ringhals Distributions
      DO i = 1, N_OUT_RG_DIST 
         WRITE(unit) ((dist_rg_col(np,ns, i), np=1, N_POLY),ns=1,NZR)
      END DO

      RETURN 
      END
