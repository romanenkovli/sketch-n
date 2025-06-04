MODULE IO_Files
   INTEGER, PARAMETER :: I_UNIT = 2, O_UNIT = 3, IO_UNIT = 4
   CHARACTER(LEN=80) :: File_GRF, File_Inp, File_Out, &
       File_TR_Out, FIle_TM_Out, File_ST_Out
END MODULE IO_Files

MODULE Formats

   IMPLICIT NONE
   CHARACTER(LEN=15) :: fmt_tr_r, fmt_tr_i, fmt_ti_r, fmt_ti_i, &
                        fmt_st_r, fmt_st_i
   INTEGER :: n_digit_tr_r, n_digit_tr_i, &
              n_digit_ti_r, n_digit_ti_i, &
              n_digit_st_r, n_digit_st_i

   CHARACTER(LEN=7), PARAMETER :: fmt_r_def = "(E12.5)", &
                                  fmt_i_def = "(I5)"
   INTEGER                     :: n_digit_r_def,  n_digit_i_def
   INTEGER, PARAMETER :: n_data_in_line=10
   CHARACTER(LEN=3)   :: TM_3D_Output
   LOGICAL Flag_Core_Absorber_Included

END MODULE Formats


MODULE Model_Parameters

! HEADER and PROBLEM TITLE
   INTEGER ::  N_LINE_HEADER, N_LINE_PROBLEM_TITLE
   CHARACTER(LEN=80), ALLOCATABLE :: HEADER(:), PROBLEM_TITLE(:)

! Geometry Description ! XYZ or HEX-Z
   CHARACTER(LEN=4) gmt_crd_type
          
! GeoMTry Module                                                      !
   INTEGER :: N_POLY, NH, NZR, NZ, NXR, NYR, NX, NY, &
                NCHM, NDD,N_BUNDLE_TYPE
! TH_Model                                                             !
   INTEGER :: NP_Reactor_Core, NZR_Core

! XS Module (Cross Section & Neutron Kinetics Constant)                !
   INTEGER :: NNODE, NG, MD, N_FEEDBACK

! CRD control rod Module                                               !
   INTEGER :: NN_CRod, NN_CRod_Comp, NN_CRod_El, &
     NN_CRod_Type, NN_Crod_Bundle

! Fuel RoD (FRD) Heat Conduction Model                                 !
   INTEGER :: NN_FRD_FUEL, NN_FRD_CLAD, NN_FRD_TOTAL

! PVM Interface Module for TRAC                                        !
   INTEGER :: NN_RT_HC_TRAC, NN_RT_FD_TRAC, NN_Z_HC_TRAC, &
      NN_Z_FD_TRAC 

! REACTOR TYPE
   CHARACTER*3 :: REACTOR_TYPE

! Number of Steady-State and Transient Distributions

   INTEGER :: N_ST_SCAL, N_ST_VEC, N_ST_DIST,       &
                    N_TR_DIST, N_TR_SCAL, N_TR_VEC
!                    N_TR_DIST_REACTOR, N_ST_DIST_REACTOR  

   INTEGER, ALLOCATABLE :: DIM_ST_VEC(:), DIM_TR_VEC(:)

   CHARACTER(LEN = 80), ALLOCATABLE ::                          &
                    NAME_ST_SCAL(:), NAME_ST_VEC(:), NAME_ST_DIST(:), &
                    NAME_TR_SCAL(:), NAME_TR_VEC(:), NAME_TR_DIST(:)             

END MODULE Model_Parameters

MODULE Model_Description

! REACTOR CORE DATA

  INTEGER, ALLOCATABLE :: Numb_Reactor_Core(:)

  INTEGER, ALLOCATABLE :: N_Coord(:,:)

  INTEGER, ALLOCATABLE ::  npoly(:,:),        &
             Nxr_B_Reactor(:),                      &
             Nxr_E_Reactor(:)                     

  INTEGER ::                                  &
             Nyr_B_Reactor, Nyr_E_Reactor,          &
             Nxr_Max_Reactor, Nxr_B_Min_Reactor

! Reactor Core Geometry
  INTEGER, ALLOCATABLE :: Index_Core(:, :),   &
              Nxr_B_Core(:),                        &
              Nxr_E_Core(:)
  INTEGER ::                                  &
              Nyr_B_Core, Nyr_E_Core,               &
              Nxr_Max_Core, Nxr_B_Min_Core,         &
              NZR_Core_Beg, NZR_Core_End

! Spatial Mesh
  REAL, ALLOCATABLE ::                        &
                               hx(:),               &
                               hy(:),               &
                               hz(:)

! Fine Axial Mesh
  INTEGER, ALLOCATABLE ::                     &                  
                                npz(:)
  REAL, ALLOCATABLE  ::                       &
                                hzt(:)

! Reactor & Core Volumes & Nodal Volumes
  REAL :: v_reactor, v_core
  REAL, ALLOCATABLE :: vol_ass(:, :)

! Control RoD Description
  INTEGER, ALLOCATABLE ::                     &                  
                                Mat_Com_Rod(:,:),   &
                                nrods(:,:)
  REAL, ALLOCATABLE :: h_rod_el(:)

  CONTAINS 

     INTEGER FUNCTION  nlz_core(nlz)
!=====================================================================!
! convert reactor coordinates into core coordinates (z)               !
!=====================================================================!
         implicit none
! Input:
         integer nlz
! Output:
!        integer nlz_core

         nlz_core = nlz - NZR_Core_Beg + 1

      RETURN
         END FUNCTION  nlz_core

      SUBROUTINE Convert_Coord_Core_xy(nlx, nly, nlx_core, nly_core)
!=====================================================================!
! convert reactor coordinates into core coordinates (x, y)            !
!=====================================================================!
         implicit none
! Input:
         integer nlx, nly
! Output:
         integer nlx_core, nly_core

         nlx_core = nlx - NXR_B_Core(nly) + 1
         nly_core = nly - NYR_B_Core      + 1

      RETURN
         END SUBROUTINE Convert_Coord_Core_xy

END MODULE Model_Description


MODULE Distributions


   USE Model_Parameters

   INTEGER, PARAMETER :: NN_TIME_STEP = 10000

   REAL :: time_sk(NN_TIME_STEP) 

   INTEGER :: N_TIME_STEP 

   REAL, ALLOCATABLE :: scal_tr(:)

   REAL, ALLOCATABLE :: vect_tr(:,:)

   REAL, ALLOCATABLE :: dist_st(:, :, :),   &
                                dist_st_mm(:,:)
   INTEGER, ALLOCATABLE ::                  &
                              k_dist_st_mm(:,:,:)
   LOGICAL, ALLOCATABLE :: flag_reactor_st_dist(:)


   REAL, ALLOCATABLE :: dist_tr(:, :, :),   &
                                dist_tr_mm(:,:)
   INTEGER, ALLOCATABLE ::                  &
                              k_dist_tr_mm(:,:,:)

   LOGICAL, ALLOCATABLE :: flag_reactor_tr_dist(:)

END MODULE Distributions

MODULE Dist_Time_Int

   USE Model_Parameters

   REAL ::  time_output

   INTEGER :: n_time_data  ! number of the data 2 or 1 

   INTEGER :: n_time_int

   REAL :: xsi_int

   REAL :: time_ti(2) 

   REAL, ALLOCATABLE :: scal_ti(:)

   REAL, ALLOCATABLE :: vect_ti(:,:)

      REAL, ALLOCATABLE :: dist_ti(:,:,:),   &
                                dist_ti_mm(:,:)
   INTEGER, ALLOCATABLE ::                  &
                              k_dist_ti_mm(:,:,:)
   

END MODULE Dist_Time_Int

MODULE TR_Output
   
   IMPLICIT NONE

   INTEGER, PARAMETER :: NN_OUT_TR = 100
   INTEGER      :: n_out_tr
   CHARACTER(LEN=1) :: type_out_tr(NN_OUT_TR)
   INTEGER          :: i_out_tr(NN_OUT_TR),     &
                             k_out_tr(3, NN_OUT_TR),  &
                             d_out_tr(NN_OUT_TR)
   INTEGER, PARAMETER :: NN_IDN_TYPE_OUT_TR = 5

   CHARACTER(LEN=2), DIMENSION(NN_IDN_TYPE_OUT_TR), PARAMETER :: &
                IDN_TYPE_OUT_TR = (/ "TS", "TV", "TD", "TL", "TM" /)

END MODULE TR_Output


MODULE TI_Output
   
   IMPLICIT NONE

   INTEGER, PARAMETER :: NN_OUT_TI = 100
   INTEGER      :: n_out_ti
   CHARACTER(LEN=1) :: type_out_ti(NN_OUT_TI)
   INTEGER          :: i_out_ti(NN_OUT_TI),     &
                             d_out_ti(NN_OUT_TI)
   REAL             :: scale_factor_ti(NN_OUT_TI)

   INTEGER, PARAMETER :: NN_IDN_TYPE_OUT_TI = 3
   CHARACTER(LEN=2), DIMENSION(NN_IDN_TYPE_OUT_TI), PARAMETER :: &
                IDN_TYPE_OUT_TI = (/ "MS", "MV", "MD" /)

END MODULE TI_Output

MODULE ST_Output
   
   IMPLICIT NONE

   INTEGER, PARAMETER :: NN_OUT_ST = 100
   INTEGER      :: n_out_st
   CHARACTER(LEN=1) :: type_out_st(NN_OUT_ST)
   INTEGER          :: i_out_st(NN_OUT_st),     &
                             d_out_st(NN_OUT_st)
   REAL             :: scale_factor_st(NN_OUT_st)

   INTEGER, PARAMETER :: NN_IDN_TYPE_OUT_st = 1
   CHARACTER(LEN=2), DIMENSION(NN_IDN_TYPE_OUT_st), PARAMETER :: &
                IDN_TYPE_OUT_st = (/ "SD" /)

END MODULE ST_Output
      
