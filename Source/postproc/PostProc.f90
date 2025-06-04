PROGRAM GRF_PostPocessing

!   WRITE(*,*) 'Input Namelist read'
   CALL INPut_Read_Namelist

!   WRITE(*,*) 'Input Geometry Data'
   CALL Input_Geometry_Data

!   WRITE(*,*) 'Read Write'
   CALL Read_Write_TR_Data


END PROGRAM GRF_PostPocessing

SUBROUTINE Input_Geometry_Data

   USE IO_Files
   USE Model_Parameters, ONLY : N_ST_DIST
   USE Distributions, ONLY : n_time_step, time_sk

   IMPLICIT NONE

   INTEGER :: ios
   write(*,*) '======================',File_GRF
   OPEN(i_unit, FILE = File_GRF, Status='OLD', &
      ACTION = 'READ', FORM = 'UNFORMATTED') 
   ! OPEN(io_unit,file=File_GRF, status='OLD', &
   !  FORM= 'UNFORMATTED')

   OPEN(o_unit, FILE = File_OUT, Status='UNKNOWN', &
      ACTION = 'WRITE') 

   CALL Read_Model_parameters(i_unit)
   CALL Allocate_Model_Description
   CALL Read_Model_Description(i_unit)

   CALL Write_Headers(o_unit)


   IF(File_Inp.NE."") THEN
      OPEN(io_unit, FILE = File_Inp, STATUS = 'OLD', ACTION = 'READ')
         CALL Read_TR_Input_Data(io_unit)
      CLOSE(io_unit)

      OPEN(io_unit, FILE = File_Inp, STATUS = 'OLD', ACTION = 'READ')
         CALL Read_TI_Input_Data(io_unit)
      CLOSE(io_unit)

      OPEN(io_unit, FILE = File_Inp, STATUS = 'OLD', ACTION = 'READ')
         CALL Read_ST_Input_Data(io_unit)
      CLOSE(io_unit)

      
      CALL WRITE_TR_Input_Data(o_unit)
      CALL WRITE_TI_Input_Data(o_unit)
      CALL WRITE_ST_Input_Data(o_unit)

   END IF

   CLOSE(o_unit) 

RETURN
END SUBROUTINE Input_Geometry_Data

SUBROUTINE Read_Model_parameters(i_unit)

   USE Model_Parameters

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i_unit

   INTEGER ::  n

! SKETCH Header
   READ(i_unit) N_LINE_HEADER

   IF (.NOT. ALLOCATED(HEADER)) ALLOCATE ( HEADER(N_LINE_HEADER) )

   DO n = 1, N_LINE_HEADER
       READ(i_unit) HEADER(n)
   END DO

! Problem Title
   READ(i_unit) N_LINE_PROBLEM_TITLE

   IF (.NOT. ALLOCATED(PROBLEM_TITLE)) ALLOCATE ( PROBLEM_TITLE( N_LINE_PROBLEM_TITLE ) )

   DO n = 1, N_LINE_PROBLEM_TITLE
       READ(i_unit) PROBLEM_TITLE(n)
   END DO
  READ(i_unit) GMT_CRD_TYPE

! GeoMTry Module                                                      !
  READ(i_unit) N_POLY, NH, NZR, NZ, NXR, NYR, NX, NY, &
              NCHM, NDD,N_BUNDLE_TYPE

! TH_Model                                                             !
   READ(i_unit) NP_Reactor_Core, NZR_Core

! XS Module (Cross Section & Neutron Kinetics Constant)                !
   READ(i_unit) NNODE, NG, MD, N_FEEDBACK

! CRD control rod Module                                               !
   READ(i_unit) NN_CRod, NN_CRod_Comp, NN_CRod_El, NN_CRod_Type, &
       NN_Crod_Bundle

! Fuel RoD (FRD) Heat Conduction Model                                 !
   READ(i_unit) NN_FRD_FUEL, NN_FRD_CLAD, NN_FRD_TOTAL

! PVM Interface Module for TRAC                                        !
   READ(i_unit) NN_RT_HC_TRAC, NN_RT_FD_TRAC, NN_Z_HC_TRAC, &
       NN_Z_FD_TRAC 

! REACTOR TYPE
   READ(i_unit) REACTOR_TYPE
! Steady-State Data

   READ(i_unit) N_ST_SCAL


   IF (.NOT. ALLOCATED(Name_ST_SCAL))  ALLOCATE (Name_ST_SCAL(N_ST_SCAL))


      DO n = 1,  N_ST_SCAL
         READ(i_unit) name_st_scal(n)
   END DO

   READ(i_unit) N_ST_VEC

   IF (.NOT. ALLOCATED(Name_St_VEC)) ALLOCATE (Name_St_VEC(N_ST_VEC))
   IF (.NOT. ALLOCATED(Dim_ST_VEC))  ALLOCATE (Dim_ST_VEC(N_ST_VEC))

      DO n = 1,  n_st_vec
         READ(i_unit) dim_st_vec(n)
         READ(i_unit) name_st_vec(n)
   END DO

   READ(i_unit) N_ST_DIST

   IF (.NOT. ALLOCATED(Name_St_Dist)) ALLOCATE (Name_St_Dist(N_ST_DIST))

      DO n = 1,  n_st_dist
         READ(i_unit) name_st_dist(n)
   END DO

! Transient Data

   READ(i_unit) N_TR_SCAL

   IF (.NOT. ALLOCATED(Name_TR_SCAL)) ALLOCATE (Name_TR_SCAL(N_TR_SCAL))

      DO n = 1,  N_TR_SCAL
         READ(i_unit) name_tr_scal(n)
   END DO

   READ(i_unit) N_TR_VEC

   IF (.NOT. ALLOCATED(Name_TR_VEC)) ALLOCATE (Name_TR_VEC(N_TR_VEC))
   IF (.NOT. ALLOCATED(Dim_TR_VEC))  ALLOCATE (Dim_TR_VEC(N_TR_VEC) )

      DO n = 1,  n_TR_vec
         READ(i_unit) dim_TR_vec(n)
         READ(i_unit) name_TR_vec(n)
   END DO
!
   READ(i_unit) N_TR_DIST

   IF (.NOT. ALLOCATED(Name_TR_Dist)) ALLOCATE (Name_TR_Dist(N_TR_DIST))

      DO n = 1,  n_tr_dist
         READ(i_unit) name_tr_dist(n)
   END DO

RETURN
END SUBROUTINE Read_Model_parameters

SUBROUTINE Allocate_Model_Description

  USE Model_Parameters, ONLY : NXR, NYR, NZR, N_POLY, NZ,          &
                               NN_CRod,NN_CRod_El, NN_Crod_Bundle, &
                               NP_Reactor_Core
  USE Model_Description

  IMPLICIT NONE

! Core Numbering
   ALLOCATE( Numb_Reactor_Core(NP_Reactor_Core) )

! Channel Coordinates (x,y)
   ALLOCATE( N_Coord(N_POLY, 2) )   

! Reactor Core Geometry
   ALLOCATE( npoly(NYR,NXR),                        &
             Nxr_B_Reactor(NYR),                    &
             Nxr_E_Reactor(NYR) )                     

! Reactor Core Geometry
   ALLOCATE (Index_Core(NYR, NXR),                  &
              Nxr_B_Core(NYR),                      &
              Nxr_E_Core(NYR) )
! Spatial Mesh
   ALLOCATE  ( hx(NXR), hy(NYR), hz(NZR) )

! Fine Axial Mesh
   ALLOCATE ( npz(NZR), hzt(NZ) )

! Reactor & Core Volumes & Nodal Volumes
  ALLOCATE ( vol_ass(N_POLY, NZR) )

! Control RoD Description
  ALLOCATE ( Mat_Com_Rod(NN_CRod,NN_CRod_El), &
             nrods(NN_CRod, NN_Crod_Bundle),  &
             h_rod_el(NN_CRod_El) )

RETURN
END SUBROUTINE Allocate_Model_Description

SUBROUTINE Read_Model_Description(unit)

  USE Model_Parameters, ONLY : NXR, NYR, NZR, N_POLY, NZ,          &
                               NN_CRod,NN_CRod_El, NN_Crod_Bundle, &
                               NP_Reactor_Core
  USE Model_Description

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: unit

  INTEGER nc
  INTEGER nd
  INTEGER nlx, nly, nlz
  INTEGER np, ns
  INTEGER ib, ie, ir
  INTEGER n1

! Core Numbering 
  READ(unit) (Numb_Reactor_Core(nc), nc = 1, NP_Reactor_Core)
! Channels Coordinates
  READ(unit) ((N_Coord(np, nd), np=1, N_POLY), nd=1,2)
! Reactor Geometry
  READ(unit) ((npoly(nly,nlx), nly=1, NYR), nlx=1, NXR),          &
                  (Nxr_B_Reactor(nly), nly=1, NYR),               &
                  (Nxr_E_Reactor(nly), nly=1, NYR),               &
                  Nyr_B_Reactor, Nyr_E_Reactor,                   &
                  Nxr_Max_Reactor, Nxr_B_Min_Reactor


! Reactor Core Geometry
  READ(unit) ((Index_Core(nly,nlx), nly=1, NYR), nlx=1, NXR),     &
                  (Nxr_B_Core(nly), nly=1, NYR),                  &
                  (Nxr_E_Core(nly), nly=1, NYR),                  &
                  Nyr_B_Core, Nyr_E_Core,                         &
                  Nxr_Max_Core, Nxr_B_Min_Core,                   &
                  NZR_Core_Beg, NZR_Core_End


! Spatial Mesh
  READ(unit) (hx(nlx), nlx = 1,NXR),                              &
                  (hy(nly), nly = 1, NYR),                        &
                  (hz(nlz), nlz = 1, NZR)
! Fine Axial Mesh
   READ(unit) (npz(ns), ns = 1, NZR),                             &
                  (hzt(n1), n1 =1, NZ)

! Reactor & Core Volumes & Nodal Volumes
   READ(unit) v_reactor, v_core,                                  &
         ((vol_ass(np, ns), np=1, N_POLY), ns=1, NZR)

! Control RoD Description
   READ(unit)                                                     &
      ((Mat_Com_Rod(ir,ie), ir = 1,NN_CRod), ie =1,NN_CRod_El),   &
      ((nrods(ir, ib), ir = 1,NN_CRod), ib = 1, NN_Crod_Bundle),  &
      (h_rod_el(ie), ie = 1, NN_CRod_El)

!   WRITE(*,*)                                                     &
!      (h_rod_el(ie), ie = 1, NN_CRod_El)
!   PAUSE



RETURN
END SUBROUTINE Read_Model_Description


SUBROUTINE INPut_Read_Namelist
!=============================================================================!
!        Input from the NAMELIST "PostProc.INI"                               !
!                   Vyachreslav Zimin (c) 25 July 2000                        !
!               vgzimin@mail.ru                         !
!=============================================================================!

   USE IO_Files
   USE Formats
   USE Dist_Time_Int, ONLY : time_output

   IMPLICIT NONE
! Local
   INTEGER ios


   NAMELIST /PostProc/  File_grf, File_Inp, File_Out,  File_TR_Out, &
                        File_TM_Out, File_ST_Out, fmt_tr_r, fmt_tr_i,   &
                        time_output, fmt_ti_r, fmt_ti_i
!=============================================================================!
! File_grf - *.grf file generated by SKETCH.exe         "Output/SKETCH.grf"   !
! File_Inp - ASCII Input File for PostProcessor         ""                    !
! File_Out - ASCII Output File                          "Output/PostProc.lst" !
! File_TR_Out - ASCII Output File                       "Output/PProc_TR.lst" !
! File_TM_Out - ASCII Output File                       "Output/PProc_TM.lst" !
! File_ST_Out - ASCII Output File                       "Output/PProc_ST.lst" !
! fmt_tr_r - Output Format for Transient Real Data                 "(E12.5)"  !
! fmt_tr_i - Output Format for Transient Integer Data              "(I5)"     !
! fmt_ti_r - Output Format for Time Moment Real Data               "(F6.3)"   !
! fmt_ti_i - Output Format for Time_Moment Integer Data            "(I5)"     !
! fmt_st_r - Output Format for Time Moment Real Data               "(F6.3)"   !
! fmt_st_i - Output Format for Time_Moment Integer Data            "(I5)"     !
! Time_Output - Output Time for Distributions                      "0.0"      !
! TM_3D_Output - Write 3D Distributions                            "YES"      !
!=============================================================================!

   File_grf = "Output/SKETCH.grf"
   File_Inp = ""
   File_Out =    "Output/PostProc.lst"
   File_TR_Out = "Output/PProc_TR.lst"
   File_TM_Out = "Output/PProc_TM.lst"
   File_ST_Out = "Output/PProc_ST.lst"
   fmt_tr_r = "(E12.5)"
   fmt_tr_i = "(I5)"
   fmt_ti_r = "(F6.3)"
   fmt_ti_i = "(I5)"
   fmt_st_r = "(F6.3)"
   fmt_st_i = "(I5)"
   time_output = 0.
   TM_3D_Output = "YES"

   OPEN(i_unit, FILE ='Input/PostProc.INI', STATUS='OLD')
      READ(i_unit, NML = PostProc)
   CLOSE(i_unit)

   CALL get_n_digit_in_format(fmt_tr_r, n_digit_tr_r)
   CALL get_n_digit_in_format(fmt_tr_i, n_digit_tr_i)

   CALL get_n_digit_in_format(fmt_ti_r, n_digit_ti_r)
   CALL get_n_digit_in_format(fmt_ti_i, n_digit_ti_i)

   CALL get_n_digit_in_format(fmt_st_r, n_digit_st_r)
   CALL get_n_digit_in_format(fmt_st_i, n_digit_st_i)

   CALL get_n_digit_in_format(fmt_r_def, n_digit_r_def)
   CALL get_n_digit_in_format(fmt_i_def, n_digit_i_def)

   
RETURN
END SUBROUTINE INPut_Read_Namelist      



SUBROUTINE OUTput_Write_Header(unit)
!=====================================================================!
! Writes Header for the SKETCH output files                           !
! 22.II.1999 (c) Slava                                                ! 
!=====================================================================!
IMPLICIT NONE
INTEGER unit

CALL OUTPut_Write_Separator_Bold(unit)

WRITE(unit, '(A)')                                       &
      "!             PostProcessing of the SKETCH-N"//   &
      "   Results                         !"    
WRITE(unit,'(A)')                                        &
      "!              Version 1.0 (c)"//                 &
      " Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !" 
CALL OUTPut_Write_Separator_Bold(Unit)

RETURN
END SUBROUTINE OUTput_Write_Header

SUBROUTINE OUTPut_Write_Separator_Bold(unit)
!=====================================================================!
! Writes a separator line  for the SKETCH output files                !
! 22.II.1999 (c) Slava                                                ! 
!=====================================================================!
   IMPLICIT NONE
   INTEGER :: unit

   CHARACTER*80 string

   WRITE(string, '( "!",78("="), "!" )' )

   WRITE(unit, '(A)') string

RETURN
END SUBROUTINE OUTPut_Write_Separator_Bold

SUBROUTINE OUTPut_Write_Separator(unit)
!=====================================================================!
! Writes a separator line  for the SKETCH output files                !
! 22.II.1999 (c) Slava                                                ! 
!=====================================================================!
   IMPLICIT NONE
   INTEGER :: unit

   CHARACTER*80 string

   WRITE(string, '( "!",78("-"), "!" )' )

   WRITE(unit, '(A)') string

RETURN
END SUBROUTINE OUTPut_Write_Separator


SUBROUTINE Write_Headers(unit)

   USE IO_Files
   USE Model_Parameters
!   USE Distributions, ONLY : time_sk, n_time_step

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: unit

   INTEGER :: n

! SKETCH PostProcessing Title

   CALL OUTput_Write_Header(unit)

   WRITE(unit, '(/, A, A)' )                               &
" Input  SKETCH   *.grf File : ", TRIM(File_grf)         
   WRITE(unit, '(A, A)' )                                  &
" Input  PostProc *.dat File : ", TRIM(File_inp)        
   WRITE(unit, '(A, A)' )                                  &
" Output PostProc *.dat File : ", TRIM(File_out)

! SKETCH Title

   DO n = 1, N_LINE_HEADER
       WRITE(unit, '(A)' ) HEADER(n)
   END DO

! Problem Title

   DO n = 1, N_LINE_PROBLEM_TITLE
       WRITE(unit, '(A)' ) PROBLEM_TITLE(n)
   END DO

! Distribution Titles
   WRITE(unit,'(/,A)')                &
"Data Available in File_GRF :"

   WRITE(unit,'(/,A)')                &
"  Steady-State Data :"

   WRITE(unit,'(/,A)')                &
"     Scalar Data :"
   DO n = 1, N_ST_SCAL
         WRITE(unit, '(5x, I3, A)') n, "- "//TRIM( NAME_ST_SCAL(n) )
   END DO
   WRITE(unit,'(/,A)')                &
"     Vector Data :"
   DO n = 1, N_ST_VEC
         WRITE(unit, '(5x, I3, A, A, I8)') & 
           n, "- "//TRIM( NAME_ST_VEC(n) ), "; Dimension =", DIM_ST_VEC(n)
   END DO

   WRITE(unit,'(/,A)')                &
"     Distributions :"
   DO n = 1, N_ST_DIST
         WRITE(unit, '(5x, I3, A)') & 
           n, "- "//TRIM( NAME_ST_DIST(n) )
   END DO
   
   WRITE(unit,'(/,A)')                &
"  Transient Data :"

   WRITE(unit,'(/,A)')                &
"     Scalar Data :"
   DO n = 1, N_TR_SCAL
         WRITE(unit, '(5x, I3, A)') n, "- "//TRIM( NAME_TR_SCAL(n) )
   END DO
   WRITE(unit,'(/,A)')                &
"     Vector Data :"
   DO n = 1, N_TR_VEC
         WRITE(unit, '(5x, I3, A, A, I8)') & 
           n, "- "//TRIM( NAME_TR_VEC(n) ), "; Dimension =", DIM_TR_VEC(n)
   END DO

   WRITE(unit,'(/,A)')                &
"     Distributions :"
   DO n = 1, N_TR_DIST
         WRITE(unit, '(5x, I3, A)') & 
           n, "- "//TRIM( NAME_TR_DIST(n) )
   END DO


RETURN
END SUBROUTINE Write_Headers   

SUBROUTINE Write_Time_Step_Data(unit)

   USE Distributions, ONLY : time_sk, n_time_step
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: unit
   INTEGER :: n

! Time Step Informations
   WRITE(unit,'(/,A,I5)')                &
         "Number of Time Steps in File_GRF :", n_time_step 
   WRITE(unit,'(/,A)')                &
         "Output in the Time Moments :"
   WRITE(unit,'(10x, 5E12.5)') (time_sk(n), n=1, n_time_step)
   
   RETURN
END SUBROUTINE Write_Time_Step_Data



SUBROUTINE Allocate_Distributions

   USE Distributions
   INTEGER n, MAX_DIM_TR_VEC

! Steady-State Distributions

   IF(N_ST_DIST /= 0 ) THEN
      ALLOCATE ( dist_st(0:N_POLY,0:NZR, N_ST_DIST),  &
                      dist_st_mm(-3:3, N_ST_DIST) )
      ALLOCATE ( k_dist_st_mm(3, -3:3, N_ST_DIST) )
   END IF

! Transient Scalar
   ALLOCATE( scal_tr(N_TR_SCAL) )

! Transient Vectors
   MAX_DIM_TR_VEC = N_POLY*NZR
   ALLOCATE( vect_tr(MAX_DIM_TR_VEC, N_TR_VEC)  )

! Transient  Distributions
   IF(N_TR_DIST /= 0 ) THEN
      ALLOCATE ( dist_tr(0:N_POLY,0:NZR, N_TR_DIST),  &
                      dist_tr_mm(-3:3, N_TR_DIST) )
      ALLOCATE ( k_dist_tr_mm(3, -3:3, N_TR_DIST) )
   END IF


RETURN
END SUBROUTINE Allocate_Distributions

SUBROUTINE Read_ST_Dist(unit)

   USE Model_Parameters
   USE Distributions

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: unit

   INTEGER :: i, np, ns

   DO i = 1, N_ST_DIST 
      READ(unit) ((dist_st(np, ns, i), np=1, N_POLY),ns=1,NZR)
   END DO

RETURN 
END SUBROUTINE Read_ST_Dist


SUBROUTINE Read_TR_Dist(unit, time, ios)

   USE IO_Files, ONLY : o_unit, File_TM_Out, File_TR_Out
   USE Model_Parameters
   USE Model_Description, ONLY : NZR_Core_Beg, NZR_Core_End, index_core
   USE Dist_Time_Int, ONLY : time_output
   USE TI_Output, ONLY : N_OUT_TI
   USE Dist_Time_Int, ONLY : n_time_data, time_ti, xsi_int, n_time_int  
   USE Distributions

   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: unit
   INTEGER, INTENT(OUT) :: ios
   REAL,    INTENT(OUT) :: time

! Local: 
   INTEGER :: ir
   INTEGER :: np, ns, i, n
   REAL, SAVE :: time_prev
   LOGICAL, SAVE :: Output_Done
   DATA Output_Done /.FALSE./

! TIME 
   READ(unit, iostat=ios ) time

! Output of the Distributions in the time interval
   IF(N_OUT_TI /= 0) THEN
      IF( time_output <= time .AND. (.NOT.Output_Done) ) THEN
         CALL Allocate_Dist_Time_Interval
         IF( n_time_step /= 1 ) THEN
            n_time_data = 2
            time_ti(1) = time_prev
            time_ti(2) = time
            xsi_int = ( time_output - time_ti(1) )/ ( time_ti(2)-time_ti(1) )
            n_time_int = n_time_step
            CALL Save_Dist_Time_Interval
         END IF ! n_time_step /= 1
      END IF ! time_output < = time
   END IF ! N_OUT_TI /= 0

   IF(ios /= 0 ) RETURN
   WRITE(*,*) 'Time =', Time
! SCALAR: 
!    Reactivity
!   WRITE(*,*) 'Transient Scalar'
   READ(unit) (scal_tr(n), n =1, N_TR_SCAL)
!   WRITE(*,*) '(scal_tr(n), n =1, N_TR_SCAL)=', (scal_tr(n), n =1, N_TR_SCAL)
!   pause

! ARRAYS: 
!    Control RoD Positions
!   WRITE(*,*) 'Transient Vectors'
   DO i = 1, N_TR_VEC
      READ(unit) (vect_tr(n,i), n=1, DIM_TR_VEC(i))
!      WRITE(*,*) 'i, DIM_TR_VEC(i) (vect_tr(n,i), n=1, DIM_TR_VEC(i))', i, DIM_TR_VEC(i)
!      WRITE(*,*) (vect_tr(n,i), n=1, DIM_TR_VEC(i))
!      PAUSE
   END DO

! DISTRIBUTION:
!   WRITE(*,*) 'Transient Distributions'
   DO i = 1, N_TR_DIST
!      WRITE(*,*) 'dist, i =', i
      READ(unit) ((dist_tr(np,ns,i), np=1, N_POLY),ns=1,NZR)
!      WRITE(*,*) ((dist_tr(np,ns,i), np=1, N_POLY),ns=1,NZR)
!      PAUSE
   END DO      

   time_prev = time

! Computing Average & Maximum Values
      DO i = 1, N_TR_DIST

         CALL Compute_3D_Average(dist_tr(0,0,i), NZR_Core_Beg, NZR_Core_End, &
            Index_Core, dist_tr_mm(3,i), dist_tr_mm(-3,i),                   &
            k_dist_tr_mm(1,3,i), k_dist_tr_mm(1,-3,i) )      
         CALL Compute_2D_Average(dist_tr(0,0,i), NZR_Core_Beg, NZR_Core_End, &
            Index_Core, dist_tr_mm(2,i), dist_tr_mm(-2,i),                   &
            k_dist_tr_mm(1,2,i), k_dist_tr_mm(1,-2,i) )      
         CALL Compute_1D_Average(dist_tr(0,0,i), NZR_Core_Beg, NZR_Core_End, &
            Index_Core, dist_tr_mm(1,i), dist_tr_mm(-1,i),                   &
            k_dist_tr_mm(1,1,i), k_dist_tr_mm(1,-1,i) )      
      END DO
! Output of the Distributions in the time moment
!      write(*,*) 'N_OUT_TI=', N_OUT_TI

      IF(N_OUT_TI /= 0) THEN
         IF( time_output <= time.AND. (.NOT.Output_Done)) THEN
            IF( n_time_step /= 1 ) THEN
               CALL Interpolate_TR_Data
            ELSE
               n_time_data = 1
               n_time_int = 1
               time_ti(1) = time 
               CALL Save_Dist_Time_Interval
            END IF
            OPEN(o_unit, FILE = File_TM_Out, STATUS = "UNKNOWN", ACTION = "WRITE")
               CALL Write_TI_Data(o_unit)
            CLOSE(o_unit)
            Output_Done = .TRUE.
         END IF ! time_output <= time.AND. (.NOT.Output_Done)
      END IF

      RETURN 
      END

SUBROUTINE Read_TR_Input_Data(unit)

   USE TR_Output

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   ios = 1
   n_out_tr = 0

   DO 
      READ(unit, '(A)', IOSTAT = ios) string
      IF(ios /= 0) EXIT

      input_string = .FALSE.
      DO n = 1, NN_IDN_TYPE_OUT_TR
         IF(string(1:2) == IDN_TYPE_OUT_TR(n) ) THEN
            n_out_tr = n_out_tr + 1
            READ(string(2:2), '(A1)') type_out_tr(n_out_tr)
            input_string = .TRUE.
            EXIT
         END IF
      END DO

      IF(input_string) THEN
         READ(string(3:5),'(I3)') i_out_tr(n_out_tr)

         SELECT CASE(type_out_tr(n_out_tr))
            CASE('M', 'L') 
               READ(string(6:80), *) k_out_tr(1,n_out_tr)
               d_out_tr(n_out_tr) = ABS( k_out_tr(1,n_out_tr) )
            CASE('D')
               d_out_tr(n_out_tr) = 3
            CASE('V')
               d_out_tr(n_out_tr) = 1
            CASE('S')
               d_out_tr(n_out_tr) = 0
         END SELECT
         
         IF((type_out_tr(n_out_tr) /= "M")                           &
             .AND.(type_out_tr(n_out_tr) /= "L")) THEN               
            READ(string(6:80), *)                                    &
               ( k_out_tr(nd, n_out_tr), nd=1, d_out_tr(n_out_tr) )
         END IF
      END IF

   END DO

RETURN
END


SUBROUTINE Read_ti_Input_Data(unit)

   USE TI_Output

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios, ios_sc_fc
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   ios = 1
   n_out_ti = 0

   DO 
      READ(unit, '(A)', IOSTAT = ios) string
      IF(ios /= 0) EXIT

      input_string = .FALSE.
      DO n = 1, NN_IDN_TYPE_OUT_ti
         IF(string(1:2) == IDN_TYPE_OUT_ti(n) ) THEN
            n_out_ti = n_out_ti + 1
            READ(string(2:2), '(A1)') type_out_ti(n_out_ti)
            input_string = .TRUE.
            EXIT
         END IF
      END DO

      IF(input_string) THEN
         READ(string(3:5),'(I3)') i_out_ti(n_out_ti)

         SELECT CASE(type_out_ti(n_out_ti))
            CASE('D')
               d_out_ti(n_out_ti) = 3
            CASE('V')
               d_out_ti(n_out_ti) = 1
            CASE('S')
               d_out_ti(n_out_ti) = 0
         END SELECT
        
        READ(string(6:80), FMT=*, IOSTAT=ios_sc_fc)  &
           scale_factor_ti(n_out_ti)
! Default values = 1.0
        IF( ios_sc_fc /= 0) THEN
           scale_factor_ti(n_out_ti) = 1.0
        END IF

!        WRITE(*,*) 'scale_factor_ti(n_out_ti) =', &
!        scale_factor_ti(n_out_ti) 


      END IF
   END DO

RETURN
END

SUBROUTINE Read_st_Input_Data(unit)

   USE ST_Output

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios, ios_sc_fc
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   ios = 1
   n_out_st = 0

   DO 
      READ(unit, '(A)', IOSTAT = ios) string
      IF(ios /= 0) EXIT

      input_string = .FALSE.
      DO n = 1, NN_IDN_TYPE_OUT_st
         IF(string(1:2) == IDN_TYPE_OUT_st(n) ) THEN
            n_out_st = n_out_st + 1
            READ(string(2:2), '(A1)') type_out_st(n_out_st)
            input_string = .TRUE.
            EXIT
         END IF
      END DO

      IF(input_string) THEN
         READ(string(3:5),'(I3)') i_out_st(n_out_st)

         SELECT CASE(type_out_st(n_out_st))
            CASE('D')
               d_out_st(n_out_st) = 3
         END SELECT
        
        READ(string(6:80), FMT=*, IOSTAT=ios_sc_fc)  &
           scale_factor_st(n_out_st)
! Default values = 1.0
        IF( ios_sc_fc /= 0) THEN
           scale_factor_st(n_out_st) = 1.0
        END IF

!        WRITE(*,*) 'scale_factor_st(n_out_st) =', &
!        scale_factor_st(n_out_st) 


      END IF
   END DO

RETURN
END
            
               
SUBROUTINE Write_TR_Input_Data(unit)

   USE TR_Output

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   WRITE(unit, '(A, I4)' ) "Number of the Output Values vs Time =",   &
      n_out_tr
   DO n = 1, N_OUT_TR
      IF( (type_out_tr(n) == "M").OR.(type_out_tr(n) == "L") ) THEN
         WRITE(unit, '("T",A1,I3,3I4)')                                 &
         type_out_tr(n),  i_out_tr(n),                                  &
         (k_out_tr(nd, n), nd=1, 1 )
      ELSE
         WRITE(unit, '("T",A1,I3,3I4)')                                 &
         type_out_tr(n),  i_out_tr(n),                                  &
         (k_out_tr(nd, n), nd=1, d_out_tr(n) )
      END IF ! 
   END DO

RETURN
END

SUBROUTINE Write_TI_Input_Data(unit)

   USE TI_Output
!   USE Dist_Time_Int

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   WRITE(unit, '(A, I4)' ) "Number of the Output Data in the Time Moment =", &
      n_out_ti

   DO n = 1, N_OUT_TI
      WRITE(unit, '("M",A1,I3, E12.5)')                                 &
      type_out_ti(n),  i_out_ti(n), scale_factor_ti(n)
   END DO

RETURN
END


SUBROUTINE Write_ST_Input_Data(unit)

   USE ST_Output
!   USE Dist_Time_Int

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN)  :: unit
! Local:
   INTEGER :: ios
   CHARACTER(LEN=80) :: string
   LOGICAL :: input_string
   INTEGER n, nd

   WRITE(unit, '(A, I4)' ) "Number of the Steady-State Output Data =", &
      n_out_st

   DO n = 1, N_OUT_st
      WRITE(unit, '("S",A1,I3, E12.5)')                                 &
      type_out_st(n),  i_out_st(n), scale_factor_st(n)
   END DO

RETURN
END


SUBROUTINE Compute_2D_Average(var_col, NZR_Beg, NZR_End, &
      Index_Map, var_2D_max, var_2D_min, k_2D_max, k_2D_min )
!=====================================================================!
! Computing 2D Average, Max-Min Values and their Locations            !
! 27.VII.2000 (c) Slava                                               ! 
!=====================================================================!
   USE Model_Parameters
   USE Model_Description

   IMPLICIT NONE
      
! Input:
   REAL, INTENT(INOUT) :: var_col(0:N_POLY, 0:NZR)
   INTEGER, INTENT(IN) :: NZR_Beg, NZR_End, Index_Map(NYR, NXR)
! Output:
!     REAL var_col(1:N_POLY, 0) - 2D distribution      
   REAL, INTENT(OUT) :: var_2D_max, var_2D_min 
   INTEGER, INTENT(OUT) :: k_2D_max(2), k_2D_min(2)
! Local Variables
   INTEGER :: n1, np, nlx, nly
      REAL vol, vol_chan

   var_2D_max = 0. ! Max
   var_2D_min = 1.E+30

   DO np = 1, N_POLY

            var_col(np,0) = 0.
         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
! Check if we are inside the Reactor Core 

         IF(Index_Map(nly, nlx) .NE. 0) THEN

            vol_chan = 0.
            DO n1 = NZR_Beg, NZR_End
               vol = vol_ass(np, n1)              
               var_col(np,0) = var_col(np,0) + var_col(np,n1)*vol
                  vol_chan = vol_chan + vol
            END DO ! NZ

            var_col(np,0) = var_col(np,0)/vol_chan
               IF(var_col(np,0) .GT. var_2D_max) THEN
                  var_2D_max = var_col(np,0)
               k_2D_max(1) = nlx
               k_2D_max(2) = nly
            END IF
               IF(var_col(np,0) .LT. var_2D_min) THEN
                  var_2D_min = var_col(np,0)
               k_2D_min(1) = nlx
               k_2D_min(2) = nly
            END IF
         END IF ! Index_Map
      END DO ! NH

RETURN
END SUBROUTINE Compute_2D_Average

SUBROUTINE Compute_1D_Average(var_col, NZR_Beg, NZR_End, &
     Index_Map, var_1D_max, var_1D_min, k_1D_max, k_1D_min)
!=====================================================================!
! Computing 2D Average, Max-Min Values and their Locations            !
! 27.VII.2000 (c) Slava                                               ! 
!=====================================================================!
   USE Model_Parameters
   USE Model_Description

   IMPLICIT NONE      
! Input:
   REAL, INTENT(INOUT) :: var_col(0:N_POLY, 0:NZR)
   INTEGER, INTENT(IN) :: NZR_Beg, NZR_End, Index_Map(NYR, NXR)
! Output:
   REAL, INTENT(OUT)  :: var_1D_max, var_1D_min
      INTEGER, INTENT(OUT) :: k_1D_max(1), k_1D_min(1)
! Local Variables
   INTEGER n1, k, np, kt, ns, nlx, nly
      REAL vol, vol_plane

     var_1D_Max = 0. 
        var_1D_Min = 1.E+30
      

      DO n1 = NZR_Beg, NZR_End
            var_col(0, n1) = 0.
         vol_plane = 0.
         DO np = 1, N_POLY
            nlx = N_Coord(np, 1)
            nly = N_Coord(np, 2)
! Check if we are inside the Reactor Core 
           IF(Index_Map(nly, nlx) .NE. 0) THEN
              vol = vol_ass(np, n1)              
              var_col(0, n1) = var_col(0, n1) + var_col(np,n1)*vol
                 vol_plane = vol_plane + vol
            END IF ! Index_Map
         END DO ! NH

         var_col(0, n1) = var_col(0, n1)/vol_plane

         IF(var_col(0, n1) .GT. var_1D_max) THEN
            var_1D_max = var_col(0, n1)
         k_1D_max = n1
      END IF

         IF(var_col(0, n1) .LT. var_1D_min) THEN
            var_1D_min = var_col(0, n1)
         k_1D_min = n1
      END IF

      END DO ! NZ

RETURN
END SUBROUTINE Compute_1D_Average
         

SUBROUTINE Compute_3D_Average(var_col, NZR_Beg, NZR_End, &
      Index_Map, var_3D_max, var_3D_min, k_3D_max, k_3D_min )
!=====================================================================!
! Computing 3D Average, Max-Min Values and their Locations            !
! 27.VII.2000 (c) Slava                                               ! 
!=====================================================================!
   USE Model_Parameters
   USE Model_Description

   IMPLICIT NONE

! Input:
   REAL,    INTENT(INOUT) :: var_col(0:N_POLY, 0:NZR)
   INTEGER, INTENT(IN)    :: NZR_Beg, NZR_End, Index_Map(NYR, NXR)
! Output:
   REAL,    INTENT(OUT)  :: var_3D_max, var_3D_min
      INTEGER, INTENT(OUT)  :: k_3D_max(3), k_3D_min(3)
! Local Variables
   INTEGER :: n1, np, nlx, nly
      REAL vol, vol_core


   var_3D_max = 0. ! Max
   var_3D_min = 1.E+30
   vol_core = 0.
   var_col(0,0) = 0.

      DO np = 1, N_POLY

         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
! Check if we are inside the Reactor Core 

         IF(Index_Map(nly, nlx) .NE. 0) THEN

            DO n1 = NZR_Beg, NZR_End
               vol = vol_ass(np, n1)              
               var_col(0,0) = var_col(0,0) + var_col(np,n1)*vol
               vol_core = vol_core + vol
                  IF(var_col(np,n1) .GT. var_3D_max) THEN
                     var_3D_max = var_col(np,n1)
                  k_3D_max(1) = nlx
                  k_3D_max(2) = nly
                  k_3D_max(3) = n1
               END IF
                  IF(var_col(np,n1) .LT. var_3D_min) THEN
                     var_3D_min = var_col(np,n1)
                  k_3D_min(1) = nlx
                  k_3D_min(2) = nly
                  k_3D_min(3) = n1
               END IF
            END DO ! NZ
         END IF ! Index_Map
     END DO ! N_POLY

     var_col(0,0) = var_col(0,0)/vol_core

RETURN
END SUBROUTINE Compute_3D_Average


SUBROUTINE Read_Write_TR_Data

   USE IO_Files
   USE Model_Parameters
   USE Model_Description
   USE Distributions
   USE Dist_Time_Int
   USE TI_Output
   USE ST_Output, ONLY : n_out_st

   IMPLICIT NONE

   INTEGER :: ios, i, n

   REAL time_tmp

   CALL Allocate_Distributions

!   WRITE(*,*) ' N Steady-State Distributions =', N_ST_DIST
   IF(N_ST_DIST /= 0) THEN
!      WRITE(*,*) 'Reading Steady-State Distributions'
      CALL Read_ST_Dist(i_unit)
      IF(N_OUT_ST /= 0) THEN
        OPEN(o_unit, FILE = File_ST_Out, STATUS = "UNKNOWN", ACTION = "WRITE")
            CALL Write_ST_Data(o_unit)   
        CLOSE(o_unit)
      END IF
   END IF


! Starting Here

! Deleting the old output
!   write(*,*) 'deleting the old output'

   OPEN(o_unit, FILE = File_TR_Out, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(o_unit, *)
   CLOSE(o_unit)
   
!   ios = 0
   n_time_step = 0
   DO ! EXIT from the loop when ios /= 0) using EXIT
      n_time_step = n_time_step + 1
!      WRITE(*,*) 'Read IN'
      CALL Read_TR_Dist(i_unit,  time_sk(n_time_step), ios)
!      WRITE(*,*) 'Read OUT, time =', time_sk(n_time_step)
!      WRITE(*,*) 'ios = ', ios
      IF( ios /= 0) EXIT
!      write(*,*) 'writing the new output'
      OPEN(o_unit, FILE = File_TR_Out, STATUS = "OLD", ACTION = "WRITE", &
             POSITION="APPEND" )
!      write(*,*) 'open done'
      CALL Write_TR_Dist(o_unit, time_sk(n_time_step))
      CLOSE(o_unit)
   END DO
   n_time_step = n_time_step - 1

   OPEN(o_unit, FILE = File_OUT, Status='UNKNOWN', &
      ACTION = 'WRITE', POSITION = 'APPEND') 
      CALL Write_Time_Step_Data(o_unit)
   CLOSE(o_unit)

   RETURN
END SUBROUTINE Read_Write_TR_Data

SUBROUTINE Write_TR_Dist(o_unit, time )

   USE Model_Description, ONLY : npoly
   USE Distributions
   USE TR_Output      
   USE Formats

   IMPLICIT NONE

! Input:
   INTEGER, INTENT(IN) :: o_unit
   REAL, INTENT(IN)    :: time
! Local
   CHARACTER(LEN = 1000) :: STRING
   INTEGER :: n_digit, ios
   INTEGER :: n_str_start, n_str_end
   CHARACTER(LEN = 15) fmt_out
   INTEGER n, np, nd
   REAL VALUE

   STRING = ""

! WRITING TIME
   n_str_start = 1
   n_digit = n_digit_tr_r
   n_str_end = n_str_start + n_digit - 1
   WRITE(string(n_str_start:n_str_end), FMT = fmt_tr_r) time
   n_str_start = n_str_end + 1

! Output Values
   DO n = 1, N_OUT_TR
      
      fmt_out = ""
      IF(type_out_tr(n) == "L") THEN
         n_digit = n_digit_tr_i*ABS(k_out_tr(1, n))

         WRITE(fmt_out(1:5), '(A1, I1, A1, I1, A1)') &
            "(", ABS(k_out_tr(1,n)), "I", n_digit_tr_i, ")"
      ELSE
         n_digit = n_digit_tr_r
         fmt_out = fmt_tr_r
      END IF

      SELECT CASE (type_out_tr(n))
         CASE("S") ! Salacr
            value = scal_tr( i_out_tr(n) )            
!            write(*,*) 'n =', "i_out_tr(n) =", i_out_tr(n), &
!               "scal_tr( i_out_tr(n) ) =", scal_tr( i_out_tr(n) )
!            pause
         CASE("V") 
            value = vect_tr(k_out_tr(1, n), i_out_tr(n))
         CASE("D")
             IF( (k_out_tr(2,n)/=0).AND.(k_out_tr(1,n)/=0) ) THEN
                np = npoly(k_out_tr(2,n), k_out_tr(1,n))
             ELSE
                np = 0
             END IF
             value = dist_tr(np, k_out_tr(3,n), i_out_tr(n) )
         CASE("M") 
              value = dist_tr_mm(k_out_tr(1, n), i_out_tr(n))
       END SELECT
       
       n_str_end = n_str_start + n_digit - 1
       IF( type_out_tr(n) == "L") THEN
          WRITE(string(n_str_start:n_str_end), FMT = fmt_out)                &
          (k_dist_tr_mm(nd,k_out_tr(1,n),i_out_tr(n)),nd=1,ABS(k_out_tr(1,n)))
       ELSE
          WRITE(string(n_str_start:n_str_end), FMT = fmt_out) value
       END IF
       n_str_start = n_str_end + 1
    END DO

    WRITE(o_unit,'(A)') TRIM(string)
RETURN
END SUBROUTINE Write_TR_Dist


REAL FUNCTION R_Interpolate_1D( xsi, y1, y2)
!=====================================================================!
!        Lagrange one-dimensional Linear Intrerpolation               !
!  (c) Slava 27.V.1999                                                !
!=====================================================================!
      IMPLICIT NONE
! Input:
      REAL y1, y2, xsi
! Output:
      R_Interpolate_1D = (1. - xsi)*y1 + xsi*y2

      RETURN
END FUNCTION R_Interpolate_1D


SUBROUTINE Allocate_Dist_Time_Interval

   USE Dist_Time_Int

   INTEGER n, MAX_DIM_TR_VEC

! Transient Scalar
   ALLOCATE( scal_ti(N_TR_SCAL) )

! Transient Vectors
   MAX_DIM_TR_VEC = N_POLY*NZR
   ALLOCATE( vect_ti(MAX_DIM_TR_VEC, N_TR_VEC)  )

! Transient  Distributions
   IF(N_TR_DIST /= 0 ) THEN
      ALLOCATE ( dist_ti(0:N_POLY,0:NZR, N_TR_DIST),  &
                      dist_ti_mm(-3:3, N_TR_DIST) )
      ALLOCATE ( k_dist_ti_mm(3, -3:3, N_TR_DIST) )
   END IF


RETURN
END SUBROUTINE Allocate_Dist_Time_Interval


SUBROUTINE Save_Dist_Time_Interval

   USE Dist_Time_Int
   USE Distributions

!Local
   INTEGER n, MAX_DIM_TR_VEC

! Transient Scalar
   scal_ti(:) = scal_tr(:)

! Transient Vectors
   vect_ti(:, :) = vect_tr(:, :)

! Transient  Distributions
   dist_ti(:, :, :) = dist_tr(:, :, :)
   dist_ti_mm(:, :) = dist_tr_mm(:, :)
   k_dist_ti_mm(:, :, :) = k_dist_tr_mm(:, :, :)

RETURN
END SUBROUTINE Save_Dist_Time_Interval


SUBROUTINE Interpolate_TR_Data

   USE Distributions
   USE Dist_Time_Int
   
!Local
   INTEGER i, n1, np, k, nd
! External Function
   REAL R_Interpolate_1D

! Transient Scalar
   DO n = 1, N_TR_SCAL
      scal_ti(n) = R_Interpolate_1D         &
         (xsi_int,scal_ti(n), scal_tr(n)) 
   END DO
! Transient Vectors
   DO n = 1, N_TR_VEC
      DO k = 1, DIM_TR_VEC(n)
         vect_ti(k,n) = R_Interpolate_1D     &
           (xsi_int,vect_ti(k,n),vect_tr(k,n)) 
      END DO
   END DO
   DO n = 1, N_TR_DIST
! Transient Distributions
      DO n1 = 0, NZR
         DO np = 0, N_POLY
         dist_ti(np,n1,n) = R_Interpolate_1D     &
           (xsi_int,dist_ti(np,n1,n),dist_tr(np,n1,n)) 
         END DO
      END DO
! Maximum & Minimum Distribution Values
      DO nd = 1, 3
         dist_ti_mm(nd,n) = R_Interpolate_1D     &
           (xsi_int,dist_ti_mm(nd,n),dist_ti_mm(nd,n))          
      END DO
   END DO

RETURN
END SUBROUTINE Interpolate_TR_Data


SUBROUTINE Write_TI_Data(unit)

   USE IO_Files
   USE Formats
   USE Model_Parameters
   USE Dist_Time_Int
   USE TI_Output

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: unit

   INTEGER :: n, NN

   INTEGER      :: n_digit_out
   CHARACTER*15 :: fmt_out
   REAL         :: scale_factor_dst  

! SKETCH PostProcessing Title

   CALL OUTput_Write_Header(unit)

   WRITE(unit, '(/, A, A)' )                               &
" Input  SKETCH   *.grf File : ", TRIM(File_grf)         
   WRITE(unit, '(A, A)' )                                  &
" Input  PostProc *.dat File : ", TRIM(File_inp)        

! SKETCH Title

   DO n = 1, N_LINE_HEADER
       WRITE(unit, '(A)' ) HEADER(n)
   END DO

! Problem Title

   DO n = 1, N_LINE_PROBLEM_TITLE
       WRITE(unit, '(A)' ) PROBLEM_TITLE(n)
   END DO

   WRITE(unit, '(A, E12.5)') "Output of the data in time moment =", &
         time_output
   IF(n_time_data == 2) THEN
      WRITE(unit, '(A, 2E12.5, A)' ) "Linear interpolation in interval [", &
          time_ti(1), time_ti(2), "]"
      WRITE(unit, '(A, F8.5)') "Interpolation factor =", xsi_int
   END IF

   DO n = 1, N_OUT_TI

      CALL OUTput_Write_Separator(unit)

      SELECT CASE(type_out_ti(n))
         CASE( "S" ) 
!            WRITE(*,*) 'FMT_TI_R =', fmt_ti_r

            CALL check_real_format(1,scal_ti(i_out_ti(n)),scale_factor_ti(n), &
               n_digit_ti_r, fmt_ti_r,  n_digit_out, fmt_out)
            CALL Write_Scalar_ti(unit, scale_factor_ti(n), fmt_out, &
                    i_out_ti(n) )

         CASE( "V" ) 

            CALL check_real_format(Dim_TR_Vec(n), vect_ti(1, i_out_ti(n)),    &
                 scale_factor_ti(n), n_digit_ti_r, fmt_ti_r,  n_digit_out,    &
                 fmt_out)
            CALL Write_Vector_ti(unit, scale_factor_ti(n), fmt_out, &
                    i_out_ti(n) )

         CASE( "D" ) 

            NN = (N_POLY+1)*(NZR+1)

            scale_factor_dst = scale_factor_ti(n)
            IF(scale_factor_ti(n).EQ.0) THEN
              scale_factor_dst = 1./dist_ti(0,0, i_out_ti(n))
            END IF
            CALL check_real_format(NN, dist_ti(0,0, i_out_ti(n)),             &
                 scale_factor_dst, n_digit_ti_r, fmt_ti_r,  n_digit_out,      &
                 fmt_out)
            CALL Write_Dist_TM(unit, Name_TR_Dist(i_out_ti(n)),              &
                 dist_ti(0,0, i_out_ti(n)), dist_ti_mm(-3, i_out_ti(n)),     &
                 k_dist_ti_mm(1, -3, i_out_ti(n)), scale_factor_dst,         &
                 n_digit_out, fmt_out)

      END SELECT
   END DO

!   CALL OUTput_Write_Separator(unit)

RETURN
END SUBROUTINE Write_TI_Data

SUBROUTINE Write_Scalar_ti(unit, scale_factor, fmt_out,  n )

   USE Formats, ONLY : fmt_r_def
   USE Dist_Time_Int
! Input:
   INTEGER, INTENT(IN) :: unit, n
   CHARACTER*(*), INTENT(IN)  :: fmt_out
   REAL, INTENT(IN) :: scale_factor
! Local:
   REAL :: value 
   LOGICAL :: error_write

   WRITE(unit, '(A)') TRIM(Name_TR_Scal(n))
   WRITE(unit, '(A, E12.5)') "Scaling Factor =", scale_factor

   WRITE(unit, fmt=fmt_out) scal_ti(n)*scale_factor
!   WRITE(*, fmt=fmt_checked) scal_ti(n)*scale_factor

RETURN
END SUBROUTINE Write_Scalar_ti
      
SUBROUTINE Write_Vector_ti(unit, scale_factor, fmt_out,  n )

   USE Dist_Time_Int
   USE Formats, ONLY : fmt_r_def, n_data_in_line

! Input:
   INTEGER, INTENT(IN) :: unit, n
   CHARACTER*(*), INTENT(IN)  :: fmt_out
   REAL, INTENT(IN) :: scale_factor
! Local


   INTEGER N_Out, N_Left, N_Cycle, nc, no
   CHARACTER*15 fmt_vec
      CHARACTER*2  char_n_out
   INTEGER ios_vec_out

   WRITE(unit, '(A)') TRIM(Name_TR_Vec(n))
   WRITE(unit, '(A, I8)') "Dimension of the Vector = ", Dim_TR_VEc(n)
   WRITE(unit, '(A, E12.5)') "Scaling Factor =", scale_factor
! CHecking OUtput   


   N_Cycle = Dim_TR_Vec(n)/n_data_in_line + 1
   N_Left = MOD(Dim_TR_Vec(n), n_data_in_line)
       IF(N_Left.EQ.0) THEN
      N_Cycle = N_Cycle - 1
      N_Left = n_data_in_line
   END IF
   N_Out = n_data_in_line
   WRITE(fmt_vec, FMT='(A,I2,A)') "(", N_Out, TRIM(fmt_out(2:) )
   DO nc = 1, N_Cycle
      IF(nc.EQ.N_Cycle) THEN 
         N_out = N_Left
         fmt_vec = ""
         WRITE(fmt_vec, FMT='(A,I2,A)') "(", N_Out, TRIM(fmt_out(2:) )
      END IF  
      WRITE(unit, FMT = fmt_vec)                            &
            (vect_ti(no+(nc-1)*n_data_in_line, n)*scale_factor,no=1, N_OUT)
!      WRITE(*, FMT = fmt_vec)                            &
!            (vect_ti(no+(nc-1)*n_data_in_line, n)*scale_factor,no=1, N_OUT)
   END DO

RETURN
END SUBROUTINE Write_Vector_ti

SUBROUTINE get_n_digit_in_format(fmt_out, n_digit)

IMPLICIT NONE

!Input:
   CHARACTER*(*), INTENT(IN) :: fmt_out
!Output:
   INTEGER, INTENT(out) :: n_digit
!Local:
   CHARACTER(LEN=2) :: char_digit
   CHARACTER(LEN=4) :: fmt_char_digit
   INTEGER :: n_fmt, n_char, ios, i_char

   n_fmt=1
   n_char=1
   DO WHILE ( (fmt_out(n_fmt:n_fmt) /= ".") .AND. &
              (fmt_out(n_fmt:n_fmt) /= ")")      )
      READ(fmt_out(n_fmt:n_fmt), FMT='(I1)', IOSTAT=ios) i_char
      IF(ios == 0) THEN
         WRITE(char_digit(n_char:n_char), FMT='(I1)') i_char
         n_char = n_char + 1
      END IF
      n_fmt = n_fmt + 1
   END DO
   n_char = n_char - 1
   WRITE(fmt_char_digit, '(A,I1,A)' ) "(I", n_char, ")"
   READ(char_digit, FMT=fmt_char_digit) n_digit
RETURN
END SUBROUTINE get_n_digit_in_format


SUBROUTINE check_real_format(nn, value, scale_factor, &
               n_digit_in, fmt_in, n_digit_out, fmt_out)

USE Formats, ONLY : fmt_r_def, n_digit_r_def

IMPLICIT NONE

INTEGER,                  INTENT(IN)  :: NN
     REAL,                INTENT(IN)  :: scale_factor, value(NN)
  INTEGER,                INTENT(IN)  :: n_digit_in
CHARACTER*(*),            INTENT(IN)  :: fmt_in
  INTEGER,                INTENT(OUT)  :: n_digit_out
CHARACTER*(*),            INTENT(OUT)  :: fmt_out
! Local
INTEGER :: ios, n
CHARACTER(LEN=n_digit_in) :: string

fmt_out = TRIM(fmt_in)
n_digit_out = n_digit_in
   DO n = 1, NN
      WRITE(string, FMT = fmt_in, IOSTAT = ios) value(n)*scale_factor
      IF(ios /= 0) THEN
         fmt_out = TRIM(fmt_r_def)
         n_digit_out = n_digit_r_def
         EXIT
      END IF
   END DO

RETURN
END 


SUBROUTINE OUTput_Distrb_Summary( unit, title, scale_factor,          &
   dist, dist_mm, k_dist_mm, fmt_out)
!=====================================================================!
! Writing Distribution Summary into "SKETCH.lst"                      !
! (c) Slava 21.VII.2000 JAERI                                         !
!=====================================================================!
      USE Model_Parameters, ONLY : N_POLY, NZR
      USE Model_Description
      USE Formats, ONLY : fmt_r_def, fmt_i_def
      IMPLICIT NONE

! Input:
      INTEGER,       INTENT(IN) :: unit 
      CHARACTER*(*), INTENT(IN) :: title
         REAL,          INTENT(IN) :: scale_factor, dist(0:N_POLY, 0:NZR), &
                                      dist_mm(-3:3) 
      INTEGER,       INTENT(IN) :: k_dist_mm(3, -3:3)
      CHARACTER*(*), INTENT(IN) :: fmt_out
!LOCAL Valus
      INTEGER n, nlx_core, nly_core
! Module function
!      INTEGER nlz_core

      WRITE(unit, '(/,A)') title

      WRITE(unit, FMT='(A)', ADVANCE = 'NO')        &
         "     Scaling Factor                        : "
      WRITE(unit, FMT=fmt_r_def) scale_factor


      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
         "     Average Value                         : " 
      WRITE(unit, FMT=fmt_out)                      &
         dist(0,0)*scale_factor

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
          "     3D Maximum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
           dist_mm(3)*scale_factor
      WRITE(unit,'(A, 3I3)')                              &
          "     Reactor Location of 3D Maximum (X,Y,Z):", &
             (k_dist_mm(n,3), n=1,3)
      CALL Convert_Coord_Core_xy                                 &
              (k_dist_mm(1,3), k_dist_mm(2,3), nlx_core, nly_core)
      WRITE(unit,'(A, 3I3)')                                     &
          "     Core    Location of 3D Maximum (X,Y,Z):",        &
           nlx_core, nly_core, nlz_core( k_dist_mm(3,3) )

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
      "     3D Minimum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
          dist_mm(-3)*scale_factor
      WRITE(unit,'(A, 3I3)')                        &
     "     Reactor Location of 3D minimum (X,Y,Z):",(k_dist_mm(n,-3), n=1,3)
      CALL Convert_Coord_Core_xy                    &
           (k_dist_mm(1,-3), k_dist_mm(2,-3),       &
            nlx_core, nly_core)
      WRITE(unit,'(A, 3I3)')                         &
      "     Core    Location of 3D minimum (X,Y,Z):",&
            nlx_core, nly_core, nlz_core( k_dist_mm(3,-3) )

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
          "     2D Maximum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
           dist_mm(2)*scale_factor
      WRITE(unit,'(A, 2I3)')                          &
      "     Reactor Location of 2D Maximum (X,Y)  :", &
            k_dist_mm(1,2), k_dist_mm(2,2)
      CALL Convert_Coord_Core_xy                      &
           (k_dist_mm(1,2), k_dist_mm(2,2),           &
           nlx_core, nly_core)
      WRITE(unit,'(A, 2I3)')                          &
      "     Core    Location of 2D Maximum (X,Y)  :", &
            nlx_core, nly_core

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     & 
      "     2D Minimum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
        dist_mm(-2)*scale_factor
      WRITE(unit,'(A, 2I3)')                           &
      "     Reactor Location of 2D minimum (X,Y)  :",  &
            k_dist_mm(1,-2), k_dist_mm(2,-2)
      CALL Convert_Coord_Core_xy                       &
           (k_dist_mm(1,-2), k_dist_mm(2,-2),          &
         nlx_core, nly_core)
      WRITE(unit,'(A, 2I3)')                           &
      "     Core    Location of 2D minimum (X,Y)  :",  &
            nlx_core, nly_core

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
          "     1D Maximum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
           dist_mm(1)*scale_factor
      WRITE(unit,'(A,6x,I3)')                          &
      "     Reactor Location of the 1D Maximum (Z):",  &
            k_dist_mm(1,1)
      WRITE(unit,'(A,6x,I3)')                          &
      "     Core    Location of the 1D Maximum (Z):",  &
            nlz_core (k_dist_mm(1,1))

      WRITE(unit, FMT='(/, A)', ADVANCE = 'NO')     &
      "     1D Minimum Value                      : "
      WRITE(unit, FMT=fmt_out)                      &
        dist_mm(-1)*scale_factor
      WRITE(unit,'(A,6x,I3)')                        &
      "     Reactor Location of the 1D minimum (Z):",&
            k_dist_mm(1,-1)
      WRITE(unit,'(A,6x,I3)')                        &
      "     Core    Location of the 1D minimum (Z):",&
            nlz_core (k_dist_mm(1,-1))

      CALL OUTput_Write_Separator(unit)

	RETURN
END SUBROUTINE OUTput_Distrb_Summary

SUBROUTINE  Write_Dist_TM(unit, name_dist,                     &
                 dist, dist_mm, k_dist_mm, scale_factor,       &
                 n_digit_out, fmt_out)

   USE Model_Parameters, ONLY : N_POLY, NZR, NXR, NYR
   USE Model_Description
   USE Formats, ONLY : TM_3D_Output
   IMPLICIT NONE
   
   INTEGER,       INTENT(IN) :: unit, n_digit_out
   CHARACTER*(*), INTENT(IN) :: name_dist, fmt_out
   REAL,          INTENT(IN) :: dist(0:N_POLY, 0:NZR), dist_mm(-3:3), scale_factor
   INTEGER,       INTENT(IN) :: k_dist_mm(3,-3:3)
  
! LOCAL
   INTEGER                    :: n, ns
   CHARACTER(LEN=80)          :: Header
      CHARACTER(LEN=n_digit_out) :: val_char(0:N_POLY)
   CHARACTER*5  fmt_char


   WRITE(unit,'(A)') name_dist
   WRITE(unit, '(A, E12.5)') "Scaling Factor =", scale_factor
   CALL OUTput_Write_Separator(unit)


   Header = " Distribution Summary:" 
   CALL OUTput_Distrb_Summary( unit, Header, scale_factor,          &
            dist, dist_mm, k_dist_mm, fmt_out)

   Header = " 2D Distribution     :" 
   DO n=1, N_POLY
      WRITE(val_char(n), FMT=fmt_out ) dist(n,0)*scale_factor
   END DO
   CALL OUT_Write_2D_Map(N_POLY, NXR, NYR, NXR_B_Min_Core,               &
         NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,   &
         index_core, Header, unit, val_char, n_digit_out)

   CALL OUTput_Write_Separator(unit)

   WRITE(unit,'(A, /)') "     1D Distribution"

   DO ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End
          WRITE(unit, FMT= "(1x, I3,: "//fmt_out(2:)  ) &
!          nlz_core( ns ) , dist(0,ns)*scale_factor
           ns , dist(0,ns)*scale_factor
   END DO

   CALL OUTput_Write_Separator(unit)

   IF(TM_3D_Output == "YES") THEN
     Header = " 3D Distribution     :" 
     CALL Write_3D_Dist(Header, unit, dist(0,0), scale_factor, &
        n_digit_out, fmt_out)
   END IF

   RETURN
END SUBROUTINE  Write_Dist_TM


SUBROUTINE OUT_Write_2D_Map(N_POLY, NXR, NYR, NXR_Beg_Min, NXR_Max, &
         NXR_Beg, NXR_End, NYR_Beg, NYR_End, npoly,              &
         Header_Map, io_out, Value, n_digits) 
!=====================================================================!
!             Write Character Array Map (Core or Reactor)             !
!                   Vyachreslav Zimin (c) 31 May 2000                 !
!               vgzimin@mail.ru                 !
!=====================================================================!
      IMPLICIT NONE
! Input:
      INTEGER N_POLY, NYR, NXR
      INTEGER NXR_Beg_Min, NXR_Max, NXR_Beg(NYR), NXR_End(NYR), &
           NYR_Beg, NYR_End, npoly(NYR, NXR)
      CHARACTER*(*) value(N_POLY) ! Output Value
      INTEGER io_out
      CHARACTER*(*) Header_Map
      INTEGER N_Digits


! Local Varia      
      CHARACTER*3 val_fmt
      INTEGER  N_Position, N_Blank, nly, nlx
      CHARACTER*5 C_Blank, C_Position
      CHARACTER*2 Char_Digits
      CHARACTER*100 fmt

!      CALL OUTput_Write_Separator(io_out)
      WRITE(io_out,'(/,A)') Header_Map
!      CALL OUTput_Write_Separator(io_out)

!      read(val_fmt, '(1X, I1)') N_Digits
      IF(N_Digits < 10) THEN
        WRITE(val_fmt, '(A,I1)') "A", N_digits
      ELSE
        WRITE(val_fmt, '(A,I2)') "A", N_digits
      END IF

      N_Position = NXR_Max - NXR_Beg_Min + 1

      WRITE(C_Position, '(I5)') N_Position


!         write(Char_Digits, '(I1)') N_Digits
      IF(N_Digits < 10) THEN
        WRITE(Char_Digits, '(I1)') N_digits
      ELSE
        WRITE(Char_Digits, '(I2)') N_digits
      END IF
      fmt = '(/,5x,'//C_Position//'I'//Char_Digits//')'
      WRITE(io_out,fmt)                          &
                 (nlx,  nlx = Nxr_Beg_Min, Nxr_Max)

      WRITE(Char_Digits, '(I2)') N_digits-1
  
      fmt = '(5x,'//C_Position//'(1x,'//Char_Digits//'("-")))'
      WRITE(io_out, fmt) 
 
         DO nly = Nyr_Beg, Nyr_End
            N_Blank = Nxr_Beg(nly) - Nxr_Beg_Min
            N_Position = Nxr_End(nly) - Nxr_Beg(nly) + 1
            WRITE(C_Position, '(I5)') N_Position
            WRITE(C_Blank, '(I5)') N_Blank*N_Digits+1

            fmt ='(I3,":",'//C_Blank//'x,'//C_Position//val_fmt//')'
            WRITE(io_out, fmt)                       &
                 nly,(value(npoly(nly,nlx)),         &
                     nlx=NXR_Beg(nly), NXR_End(nly))


         END DO ! nly

      WRITE(io_out,*)
!      CALL OUTput_Write_Separator(io_out)

      RETURN
END SUBROUTINE OUT_Write_2D_Map

SUBROUTINE Write_3D_Dist(Header, unit, dist, scale_factor, &
              n_digit_out, fmt_out)

   USE Model_Parameters, ONLY : N_POLY, NZR, NP_Reactor_Core
   USE Model_Description, ONLY : Numb_Reactor_Core, &
                                 NZR_Core_Beg, NZR_Core_End
   USE Formats, ONLY : n_data_in_line
   IMPLICIT NONE

! Input:
   INTEGER,       INTENT(IN) :: unit, n_digit_out
   CHARACTER*(*), INTENT(IN) :: Header, fmt_out
   REAL,          INTENT(IN) :: dist(0:N_Poly, 0:NZR), scale_factor

! Local:
   INTEGER :: n_out, n_left, n_cycle, nc, no, ns, np, n, i
   CHARACTER*100 fmt_wr
   CHARACTER*2 char_n_digit, char_n_out, char_n_digit1

   WRITE(unit,'(A)') Header

   WRITE(char_n_digit, '(I2)') n_digit_out
   WRITE(char_n_digit1, '(I2)') n_digit_out - 1

   N_Cycle = NP_Reactor_Core/n_data_in_line + 1
   N_Left = MOD(NP_Reactor_Core, n_data_in_line)
       IF(N_Left.EQ.0) THEN
      N_Cycle = N_Cycle - 1
      N_Left = n_data_in_line
   END IF
   N_Out = n_data_in_line

   WRITE(char_n_out, '(I2)') n_out

   DO nc = 1, N_Cycle
      IF(nc.EQ.N_Cycle) THEN 
         N_out = N_Left
         char_n_out = ""
         WRITE(char_n_out, '(I2)') n_out
      END IF  

      fmt_wr = "(/, A5,"//char_n_out//"I"//char_n_digit//")"
      WRITE(unit, FMT = fmt_wr)  " CHAN",                     &
                  (no+(nc-1)*n_data_in_Line, no=1, N_OUT)
!      WRITE(*, FMT = fmt_wr)  " CHAN",                     &
!                  (no+(nc-1)*n_data_in_Line, no=1, N_OUT)
      fmt_wr = "(5x, "//char_n_out//'(" ",'//char_n_digit1//'("-") ) )'
!      WRITE(*, FMT=fmt_wr) 
      WRITE(unit, FMT=fmt_wr) 
      fmt_wr = '(1x, I3, ":", '//char_n_out//TRIM(fmt_out(2:))
      DO ns = NZR_Core_Beg, NZR_Core_End
!      WRITE(*, FMT = fmt_wr)  ns,                                    &
!            (dist(Numb_Reactor_Core(no+(nc-1)*n_data_in_line), ns)*  &
!            scale_factor, no=1, N_OUT)
      WRITE(unit, FMT = fmt_wr)   ns,                                &
            (dist(Numb_Reactor_Core(no+(nc-1)*n_data_in_line), ns)*  &
            scale_factor, no=1, N_OUT)
      END DO
   END DO

   RETURN
END SUBROUTINE Write_3D_Dist

SUBROUTINE Write_ST_Data(unit)

   USE IO_Files
   USE Formats
   USE Model_Parameters
   USE Model_Description, ONLY : NZR_Core_Beg, NZR_Core_End, index_core
   USE Dist_Time_Int
   USE ST_Output
   USE Distributions

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: unit

   INTEGER :: n, NN, i

   INTEGER      :: n_digit_out
   CHARACTER*15 :: fmt_out
   REAL         :: scale_factor_dst  

! Preparing Average and Maximum Data
      DO i = 1, N_ST_DIST
         CALL Compute_3D_Average(dist_st(0,0,i), NZR_Core_Beg, NZR_Core_End, &
               Index_Core, dist_st_mm(3,i), dist_st_mm(-3,i),                   &
                  k_dist_st_mm(1,3,i), k_dist_st_mm(1,-3,i) )      
         CALL Compute_2D_Average(dist_st(0,0,i), NZR_Core_Beg, NZR_Core_End, &
            Index_Core, dist_st_mm(2,i), dist_st_mm(-2,i),                   &
            k_dist_st_mm(1,2,i), k_dist_st_mm(1,-2,i) )      
         CALL Compute_1D_Average(dist_st(0,0,i), NZR_Core_Beg, NZR_Core_End, &
            Index_Core, dist_st_mm(1,i), dist_st_mm(-1,i),                   &
            k_dist_st_mm(1,1,i), k_dist_st_mm(1,-1,i) )      
      END DO

! Starting Output
! SKETCH PostProcessing Title

   CALL OUTput_Write_Header(unit)

   WRITE(unit, '(/, A, A)' )                               &
" Input  SKETCH   *.grf File : ", TRIM(File_grf)         
   WRITE(unit, '(A, A)' )                                  &
" Input  PostProc *.dat File : ", TRIM(File_inp)        

! SKETCH Title

   DO n = 1, N_LINE_HEADER
       WRITE(unit, '(A)' ) HEADER(n)
   END DO

! Problem Title

   DO n = 1, N_LINE_PROBLEM_TITLE
       WRITE(unit, '(A)' ) PROBLEM_TITLE(n)
   END DO

   WRITE(unit, '(A)') "Output of the steady-state Distributions :"

   DO n = 1, N_OUT_ST

      CALL OUTput_Write_Separator(unit)

      SELECT CASE(type_out_st(n))
         CASE( "D" ) 
            NN = (N_POLY+1)*(NZR+1)
            scale_factor_dst = scale_factor_st(n)
            IF(scale_factor_st(n).EQ.0) THEN
              scale_factor_dst = 1./dist_st(0,0, i_out_st(n))
            END IF
            CALL check_real_format(NN, dist_st(0,0, i_out_st(n)),             &
                 scale_factor_dst, n_digit_st_r, fmt_st_r,  n_digit_out,      &
                 fmt_out)
            CALL Write_Dist_TM(unit, Name_ST_Dist(i_out_st(n)),              &
                 dist_st(0,0, i_out_st(n)), dist_st_mm(-3, i_out_st(n)),     &
                 k_dist_st_mm(1, -3, i_out_st(n)), scale_factor_dst,         &
                 n_digit_out, fmt_out)

        END SELECT
   END DO

!   CALL OUTput_Write_Separator(unit)

RETURN
END SUBROUTINE Write_st_Data


