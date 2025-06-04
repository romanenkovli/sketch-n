
      subroutine INPut_Set_Param(i_sk_end, Time)
!=====================================================================*
!        Initialization of the Parameters                             *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Time, dt_save,  xme_ini, i_adjoint, Br_search, NG
!      xme_ini - estimate of the maximum eigenvalue of the 
!                  iteration matrix
! dt_save - value of the time step size from the Restart File
! dt - Value of the time step size from 'Neutron.Ini'
! Output: 
!      i_sk_end - Criterion for the Steady-State Calculations
!      npolin - number of the current Chebyshev polynomial
!      delnp - residual used to calculate the maximum eigenvalue 
!              of the iteration matrix
!      rc_cheb - quotient of the residuals at the current and 
!                 previous iterations
!      xme_ - estimate of the maximum eigenvalue of the 
!                  iteration matrix
!     d_kef_l - Local Error in Eigenvalue
!     d_flux_ l - Local Error of the Neutron Flux  
!     NG_BEG, NG_END, NG_Step - Neutro Group Numbering for the 
!               Source Iterations
! I_Bor_Start - Parameter of the Initialization of the Boron Critical
!               Search

      integer i_sk_end
      real time
!      real dt_sketch

! iteration counters
      i_sk_end = 0

! Chebyshev Polynomials
      npolin = -1
      delnp = 1.
      rc_cheb = 1.
      xme_ = xme_ini
      kin_k_ef = 0.99


! Convergence Criteria
      d_kef_l = 1.
      d_flux_l = 1.

! Neutron Group Order  of the Source Iterations

      if(Steady_State_Type.EQ."BoronSearch") then
            i_Bor_Start = 1
      end if

      i_view = 1

      dt_sketch = dt_input(1)


      if( (Time + dt_sketch).gt.TTV(NP_VIEW) ) then

         call MSC_ERR_Add_Message('WARNING',&
        'Please, Check value of the Time Transient: '//&
        'TTV(NP_VIEW) in "Neutron.INI" ' //&
        '(Time + dt_sketch).gt.TTV(NP_VIEW)' )
      end if

      if(Problem_Type.NE."Kinetics") then
            eigenv_shift = delta_shift/(1. + k_ef*delta_shift)
!           write(*,*) 'eigenv_shift =', eigenv_shift
!           pause
            eigenv = k_ef/(1. - k_ef*eigenv_shift)
      else
            eigenv_shift = 1.
      end if
! 
      if(Kinetics_Method.NE."IQS") then 
         Point_Deriv = 0.
         Point_Ampl_R = 1.
      end if

      IF(XS_Model.EQ. "TABLE") CALL XSR_Initialize_XS_Table

      IF( Xe_Sm_Model(1:1) == "0" ) THEN
        conc_isotope(1:NH,1:NZ,1) = 0.
        conc_isotope(1:NH,1:NZ,2) = 0.
      END IF

      IF( Xe_Sm_Model(2:2) == "0" ) THEN
        conc_isotope(1:NH,1:NZ,3) = 0.
        conc_isotope(1:NH,1:NZ,4) = 0.
      END IF



      return 
      end


      subroutine INPut_Read_Namelist
      USE core_history, ONLY : file_burnup_history 
      implicit none
      include 'sketch.fh'
!=====================================================================*
!        Input from the NAMELIST "NEUTRON.INI"                        *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

!=============================================================================!
! INI_PROBLEM  - the basic problem description & solution methods             !
! Problem_Type - type of the computed problem                                 !
!                  "Steady-State", Burnup", "Kinetics"   / "Steady-State"   / !
! Steady_State_Type - type of the steady-state problem                        !
!                  "Eigenvalue",  "BoronSearch"          / "Eigenvalue"   /   !
! Nodal_Method - Choice of the Nodal Method:                                  !
!               "PNM", "SANM", "ANM", "MCFD", "PNM1"               / "SANM" / !
! TH_Model - Thermal-hydraulics Model Used for the Calculations    / 'None' / !
!            (Possible Choices 'Internal', 'External', 'None', 'SKAZKA')      !
! Kinetics_Method - Method of the Solution of Kinetics Equation               !
!      "DRT"-Direct; "IQS"- Improved Quasi-Static; "PNT" - point   / "DRT"  / !  
! Iter_Solver - iterative solvers for kinertics :                             !
!                  CSA, CSI, CG, BCGSTAB, TFQMR, FOM, GMRES        / "CSA" /  !
! TRL_Approx - transverse leakage approximation                    / "QLA" /  !
!            (Possible Choices 'QLA', 'Flat')                                 !
! NonlinearIterations - nonlinear iteration procedure              / "Smith"/ !
!            ( Possible Choices "Smith", "Moon")                              !
! Xe_Sm_Model   - model of Xe and Sm for burnup                     / "nn" /  !
!            ( Possible Choices "nn", "ss", "st", "ts", "tt" )                !    
!=============================================================================!
      
      namelist /ini_problem/   Problem_Type, Nodal_Method, &
         TH_Model, Kinetics_Method, Iter_Solver, TRL_Approx,&
         NonlinearIterations, Steady_State_Type, Xe_Sm_Model,&
         iflag_divide_keff
!=============================================================================!
! INI_FILES - name of the input and restart files                             !
! File_INPUT - Input File Name of the Geometry Data               /   ""   /  !
! File_Reference - Name of the Reference Solution File            /   ""   /  !
! File_Dist - name of the DIST*.txt file with distributions for               !
!              Ringhals                                           /   ""  /   !
! File_CD - name of the CD*.txt file with XS Tables (ringhals)    /   ""  /   !
! File_Map   - Input File for the Mapping to the TRAC Geometry                !
!               Only for the Coupling with the TRAC               /   ""   /  !
! File_DMP_In - Input Restart File                                /   ""   /  !
! File_DMP_Out_st - Output Restart File for                                   !
!                                  Steady-State Calculations      /   ""   /  !
! File_DMP_Out_Kin - Output Restart File for                                  !
!                                  Kinetics     Calculations      /   ""   /  !
! File_dmp_Skazka_in   - Input  Restart File for SKAZKA  Module   /   ""   /  !  
! File_dmp_Skazka_Out  - Output Restart File for SKAZKA  Module   /   ""   /  !  
! File_BRN_In - Input Restart File with BURNUP Distriibution &                !
!                 Xe and Sm Isotopes                              /    ""  /  ! 
! File_BRN_Out - Output Restart File with BURNUP Distriibution &              !
!                 Xe and Sm Isotopes                              /    ""  /  !   
! File_burnup_history - Name of the input file with reactor power             !
!                and control rod history for burnup calculations  /    ""  /  !   
!=============================================================================!
      namelist /ini_files/&  
        file_dmp_in, file_dmp_out_st, file_dmp_out_kin,&
        File_dmp_Skazka_Out, File_dmp_Skazka_in,&
        FILE_INPUT,  File_Reference,&
        File_Map, File_Dist, File_CD, File_BRN_In,&
        File_BRN_Out, file_burnup_history 

!=============================================================================!
! INI_CONVERGENCE  - convergence criterion                                    !
! N_Out_Max  - Convergence Criterion for the Nonlinear Iterations             !
!                                                      Steady-State   /1 /    !
!                                                      Kinetics       /1000 / !
! N_Outer - Number of outer per nonlinear (thermalhydraulics)                 !
!                                                      Steady-State   /10  /  !
!                                                      Kinetics       /1000/  !
! E_Boron_Start - Accuracy of the Eigenvalue when critical boron              !
!            search starts                                            /1.E-2/ !
! E_Critical - Convergence Criterion of the Boron Critical Search     /1.E-5/ !
! E_Outer_L  - Converence Criterion for the Eigenvalue (Local)        /1.E-5/ !
! E_Flux_L - Convergence Criterion for the Flux                               !
!                                                      Steady-State   /1.E-5/ !
!                                                      Kinetics       /1.E-4/ !
! N_Inter     - Number of internal iterations per outer               /2    / !
! E_Inter - Converence Criterion for the Innner Iterations            /1.E-8/ !
!=============================================================================!

      namelist /ini_convergence/  e_inter, e_outer_l, &
        e_boron_start, n_inter, e_flux_l, e_critical, n_outer, &
        n_out_max

!=============================================================================!
! INI_CHEBYSHEV - Parameters for Chebyshev acceleration procedure             !
! Xbe - Estimate of Minimum Eigenvalue of the Iteration Matrix        /0.0  / !
! Xme_Ini - Estimate of the Dominance Ratio (Spectral Radius)         /0.8  / !
! F_Cheb - Adaptive Parameter of the Chebyshev Procedure                      !
!           (Non-Adaptive - 0, Adaptive - 0.65-0.85)                          !
!                                                      Steady-State   /0.0  / !
!                                                      Kinetics       /0.8  / !
! NPolins - Number of the Iteration when Chebyshev Starts             /5    / !
! Delta_Shift - Inverse Value of the Wieland shift                    /0.  /  !
!=============================================================================!
      namelist /ini_Chebyshev/   xbe, xme_ini,  f_cheb, npolins, &
                                 delta_shift

!=============================================================================!
! INI_TIME_STEP - Time Step Size Selection                                    !
! I_Auto -  Flag of the Automatic Time Step Selection             /    0    / !
! NP_View - Number of Time Step Interval                          /    1    / !
! Ttv - Time Moments of the Time Step Intervals                               !
!       Ttv(NP_View) - End of the Transien                        /    1.   / !
! Dt_Input(NP_View) - Time Step Size for the Time Interval                    !
!         [Ttv(i-1) - Ttv(i)]                                     /  0.01   / !
! St_Eps - Accuracy Criterion of the Automatic Time Step Control  /  5.E-03 / !
! Facmax - Maximum Increase of the Time Step Size                 /    2    / !
! Dt_Step_Max - Maximum Time Step Size (s)                        /    1.   / !
! N_Zap -  Output into GRF file at the  N_ZAPth Time Step         /    1    / !
!=============================================================================!
      namelist /ini_time_step/  st_eps, facmax, i_auto, np_view, ttv, &
                               dt_input, Kinetics_Method, &
                               dt_Step_Max, N_ZAP

      integer n,i

! Local
      Character*100 Message
   

          Problem_Type      = "Steady-State"
          Steady_State_Type = "Eigenvalue"

          Nodal_Method = "SANM"
          TH_Model= 'None'
          Kinetics_Method = "DRT"
          Iter_Solver = "CSA"
          TRL_Approx = "QLA"
          NonlinearIterations = "Smith"
          Xe_Sm_Model         = "nn"
          iflag_divide_keff   = 1
       
         open(io_unit,file='Input/SKETCH.INI',status='OLD')
          read(io_unit, NML = ini_problem)
         close(io_unit)

! Checking Input Problem_Type
         IF((Problem_Type.NE."Steady-State").AND.&
           (Problem_Type.NE."Burnup").AND.&
           (Problem_Type.NE."Kinetics") ) THEN
            WRITE(Message, '(A,A,A)') &
            "Problem_Type =", TRIM(Problem_Type), &
            "in SKETCH.INI, known types 'Steady-State', "//&
            " 'Burnup' and 'Kinetics'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF

! Checking Input Steady_State_Type
         IF((Steady_State_Type.NE."Eigenvalue").AND.&
           (Steady_State_Type.NE."BoronSearch") ) THEN
            WRITE(Message, '(A,A,A)') &
            "Problem_Type =", Problem_Type, &
            "in SKETCH.INI, known types 'BoronSearch', "//&
            " 'Eigenvalue'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF

! Checking Input Xe_Sm_Model
         DO i = 1, 2 
         IF((Xe_Sm_Model(i:i).NE."s")  .AND.&
           (Xe_Sm_Model(i:i).NE."t")  .AND.&
           (Xe_Sm_Model(i:i).NE."n")  .AND.&
           (Xe_Sm_Model(i:i).NE."0") ) THEN
            WRITE(Message, '(A,i2,A,A,A)') "i=", i,&
            "Xe_Sm_Model =", Xe_Sm_Model, &
            "in SKETCH.INI, known types 's', "//&
            " 't', 'n' or '0' "
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF
        END DO

! Checking Input TH_Model
         IF((TH_Model.NE."Internal")   .AND.&
           (TH_Model.NE."External")   .AND.&
           (TH_Model.NE."SKAZKA")     .AND.&
           (TH_Model.NE."Simulator")  .AND.&
           (TH_Model.NE."Athlet")  .AND.&
           (TH_Model.NE."None") ) THEN
            WRITE(Message, '(A,A,A)') &
            "TH_Model =", Problem_Type, &
            "in SKETCH.INI, known types 'Internal', "//&
            " 'External' and 'None'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF

! Checking Input Nodal Method
         IF((Nodal_Method.NE."PNM").AND.&
           (Nodal_Method.NE."SANM").AND.&
           (Nodal_Method.NE."ANM").AND.&
           (Nodal_Method.NE."PNM1").AND.&
           (Nodal_Method.NE."MCFD") ) THEN
            WRITE(Message, '(A,A,A)') &
            " Nodal_Method =", Nodal_Method,&
            "in SKETCH.INI, known types 'PNM', 'SANM' "//&
            " 'ANM', PNM1' and 'MCFD'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF
            
! Checking Input Kinetics Method  
         IF((Kinetics_Method.NE."DRT").AND.&
           (Kinetics_Method.NE."IQS").AND.&
           (Kinetics_Method.NE."PNT") ) THEN
            WRITE(Message, '(A,A,A)') &
            "Kinetics_Method =", Kinetics_Method, &
            "in SKETCH.INI, known types 'DRT', "//&
            " 'IQS' and 'PNT'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF


! Checking Kinetics Iterative Solver  CSA, CSI, CG, BCGSTAB, TFQMR, FOM, GMRES 
         IF((Iter_Solver.NE."CSA").AND.&
           (Iter_Solver.NE."CSI").AND.&
           (Iter_Solver.NE."CG").AND.&
           (Iter_Solver.NE."BCGSTAB").AND.&
           (Iter_Solver.NE."TFQMR").AND.&
           (Iter_Solver.NE."FOM").AND.&
           (Iter_Solver.NE."GMRES")      ) THEN
            WRITE(Message, '(A,A,A)') &
            "Iter_Solver =", Iter_Solver, &
            "in SKETCH.INI, known types: "//&
            "CSA, CSI, CG, BCGSTAB, TFQMR, FOM, GMRES "
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF

         FILE_INPUT = ""         
         FILE_DMP_IN =""
         FILE_DMP_OUT_ST =""
         FILE_DMP_OUT_KIN =""
         FILE_MAP = ""
         File_Dist = ""
         File_CD = ""
         File_Reference = ""
         File_dmp_Skazka_Out =""
         File_dmp_Skazka_in  =""
         File_BRN_in   =""
         File_BRN_out  =""
         file_burnup_history =""

         open(io_unit, file = 'Input/SKETCH.INI', status = 'OLD')
          read(io_unit, NML = ini_files)
         close(io_unit)

         IF(Problem_Type.EQ."Kinetics") THEN
            IF(FILE_DMP_IN.EQ.FILE_DMP_OUT_KIN) THEN
               WRITE(Message, '(A)') &
            " FILE_DMP_IN = FILE_DMP_OUT_KIN. "//&
            "They are usually different in Kinetics Calculations "
              CALL MSC_ERR_Add_Message('Warning',Message)
            END IF
          END IF

          if(Problem_Type.EQ."Kinetics") then
             N_OUT_MAX = 1000
             N_OUTER = 1000
             E_FLUX_L = 1.E-04
           else
              N_OUT_MAX = 1
              N_OUTER = 10
              E_FLUX_L = 1.E-05
               IF(Problem_Type.EQ."Burnup") THEN
                E_FLUX_L = 1.E-04
              END IF
          end if

          E_BORON_START = 1.E-02
          E_CRITICAL = 1.E-05
          E_OUTER_L = 1.E-5
          N_INTER = 2
          E_INTER = 1.E-8


         open(io_unit, file = 'Input/SKETCH.INI', status = 'OLD')
          read(io_unit, NML = ini_convergence)
         close(io_unit)
! for the KINETICS N_OUTER = 1000 ALWAYS
!         if(Problem_Type.EQ."Kinetics") N_OUTER = 1000
!          WRITE(*,*) 'N_OUTER =', N_OUTER
!          pause
    

!        parameter for chebyshev acceleration
          xbe = 0.
          xme_ini = 0.8
          if(Problem_Type.EQ."Kinetics") then
             F_Cheb = 0.8
           else
             F_Cheb = 0.0
          end if

          Npolins = 5   
          Delta_Shift = 0.

         open(io_unit, file = 'Input/SKETCH.INI', status = 'OLD')
          read(io_unit, NML = ini_Chebyshev)
         close(io_unit)


         I_AUTO = 0
         NP_VIEW = 1
         do n = 1, NP_VIEW
            DT_INPUT(n) = 0.01  
         end do
         TTV(NP_VIEW) = 1.E+30
         ST_EPS = 5.E-03
         FACMAX = 2.
         dt_Step_Max = 1.
         N_ZAP =  1 

         open(io_unit,file='Input/SKETCH.INI',status='OLD')
          read(io_unit, NML = ini_time_step)
         close(io_unit)

         return
         end


      subroutine INPut_Read_Restart_File(Time)
!**********************************************************************
!       Input of the Neutron Flux data from the Restart File or       *
!        Flat Neutrn Flux for the Steady-State Calculations           *    
!**********************************************************************
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'

! Input From The File:  Flux(NG, 0:N_TOT) - Neutron Flux
!                       k_ef - Eigenvalue
!                       trl_xyz - Tansverse Lekage in X - Y - Z directions
!  For Neutro Kinetics ONLY  
!                       Prec(m,k) - Delayed Neutron Precursor Concentration
!                       xrods(NN_CRod) - position of the Control Rods
!                       Source(N_TOT) - the source Terms ??     
!                       p(NH, NZ) - Power Density 
!                       p_average - Average Power Density, 
!                       p_total - Total Reactor Power
!                       time - current time
!                       dt_save - last value of the time step size
      real time
! Local Variables:
      integer n,  k,  nl

      if(Problem_Type.NE."Kinetics") then
! reading total restart file
          if(FILE_DMP_IN.NE."") then
            open(io_unit,file=FILE_DMP_IN,status='old',&
                                                 form='UNFORMATTED')
              CALL INPut_read_restart_steady_state_data(io_unit)
!              IF(TH_Model.EQ."Internal") THEN
              CALL THM_read_data_restart_file(io_unit)
!              END IF
              CALL INput_read_restart_burnup_data(io_unit)
            close(io_unit)
! reading restart file for burnup only
          else

! Flat Neutron Flux

             k_ef = 1.
             do k = 1, N_TOT
                do n = 1, NG
                   Flux(n, k) = 1.E+0
                end do
             end do
! adjoint flux is set to 1 (not computed yet)
           DO k = 1, N_TOT
             DO n = 1, NG
                Flux_A(n,k) = 1. ! Flux(n,k)
             END DO
           END DO

             do nl = 1, N_FACES
                do n = 1, NG
                   MAT_Nod(n, nl) = 0.
                end do
             end do

         end if ! if(FILE_DMP_IN.NE."")

          IF(FILE_BRN_IN.NE."") then
            open(io_unit,file=FILE_BRN_IN,status='old',&
                                                 form='UNFORMATTED')
              CALL INput_read_restart_burnup_data(io_unit)
            close(io_unit)
          END IF 


      else ! Input for Neutron Kinetics calculations
        open(io_unit,file = FILE_DMP_IN,STATUS='OLD',form='UNFORMATTED')

         CALL INPut_read_restart_steady_state_data(io_unit)
!         IF(TH_Model.EQ."Internal") THEN
         CALL THM_read_data_restart_file(io_unit)
!         END IF
         CALL INput_read_restart_burnup_data(io_unit)
         CALL INPut_read_restart_kinetics_data(io_unit, time)
         CALL INPut_read_restart_adjoint_flux(io_unit)
!         write(*,*) 'Flux_a(:,:)'
!           write(*,*) Flux_a(:,:)
!           read(*,*)

        close(io_unit)

         do n = 1, NG
           do k = 1, N_TOT
              Flux_k(n,k) = Flux(n,k)
           end do
         end do

      end if


! Empty Elements 


      return
      end


      subroutine INPut_Read_Restart_Steady_State_Data(unit)
!**********************************************************************
!       Input of the steady-state data from the Restart File or       *
!        Flat Neutrn Flux for the Steady-State Calculations           *    
!**********************************************************************
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'
! Input: 
      integer unit
! Local:
      INTEGER n, k, nd, n1, i, nl

      read(unit) (( Flux(n,k), n=1,NG), k = 1, N_TOT)
      read(unit) k_ef
      read(unit) (((trl_xyz(n,k,nd), n = 1, NG), &
            k = 1,N_TOT), nd = 1, NDD)
      read(unit) ((p(k,n1),k=1,NH),n1=1,NZ), p_col(0,0), p_total 
      read(unit) (((fdback(k,n1,i),k=1,NH),n1=1,NZ), &
                         i = 1, N_FEEDBACK)
      read(unit)  ((((D_Nod(n,k,i,nd), n=1, NG), &
        k = 1, N_TOT),i=1,2),nd=1,NDIR) 
      read(unit)  ((MAT_FD(n, nl), n=1, NG), nl = 1, N_FACES)
      read(unit)  ((MAT_Nod(n, nl), n=1, NG), nl = 1, N_FACES)

      RETURN
      END

      SUBROUTINE INput_read_restart_burnup_data(unit)
!**********************************************************************
!      Output of the steady-state DATA into  Restart File             *    
!**********************************************************************
      USE GeoMeTry_Faces
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      INTEGER unit
! Local
      INTEGER k, n1, i

      READ(unit) ((brn(k,n1),k=1,NH),n1=1,NZ)
      READ(unit) (((conc_isotope(k,n1,i),k=1,NH),n1=1,NZ),&
       i=1,N_ISOTOPE)

      RETURN
      END


      subroutine INPut_Read_Restart_Kinetics_Data(unit, time)
!**********************************************************************
!       Input of the steady-state data from the Restart File or       *
!        Flat Neutrn Flux for the Steady-State Calculations           *    
!**********************************************************************
      implicit none
      include 'sketch.fh'
! Input: 
      integer unit
! Output:
      real time
! Local 
      INTEGER m, k, i, n1

      read(unit) ((Prec(m,k),m=1,MD), k=1, N_TOT)
      read(unit) (zrods(i),i=1,NN_CRod), Time_Scram, &
                           Flag_Set_Time_Scram
      read(unit) (source(k), k = 1, N_TOT)
      read(unit) time, dt_save
      read(unit) Pow_point, (Prec_point(m), m = 1, MD), react
      read(unit) (((p_dh(m, k, n1), m =1, MH), k = 1, NH), n1=1, NZ)


!      write(*,*) 'SOURCE =, k = 205'
!      k = 205
!      DO n1=1, NZ/2
!        kt = k + (n1-1)*NH
!        write(*,*) 'n1, source =', n1, source(kt),&
!        flux(1,kt), flux(2,kt)
!      END DO 
!      Pause
!      write(*,*) 'SOURCE =, k = 206'
!      k = 206
!      DO n1=1, NZ/2
!        kt = k + (n1-1)*NH
!        write(*,*) 'n1, source =', n1, source(kt)
!      END DO 
!      pause
      
      RETURN
      END

      subroutine INPut_Read_Restart_adjoint_flux(unit)
!**********************************************************************
!       Input of the steady-state data from the Restart File or       *
!        Flat Neutrn Flux for the Steady-State Calculations           *    
!**********************************************************************
      implicit none
      include 'sketch.fh'
! Input: 
      integer unit
! Local:
      INTEGER n, k

      read(unit) (( Flux_a(n,k), n=1,NG), k = 1, N_TOT)

      RETURN
      END


      SUBROUTINE INPut_Set_Distr_Names
      implicit none
      include 'sketch.fh'
      INTEGER n,  i

      NAME_ST_DIST(1) = "BURNUP, [MWd/tU]"
      NAME_ST_DIST(2) = "XENON, [1/cm^3]"
      NAME_ST_DIST(3) = "VOID HISTORY, [-]"
      NAME_ST_DIST(4) = "CONTROL ROD HISTORY, [-]"
      NAME_ST_DIST(5) = "CONVERSION HISTORY, [-]"
      NAME_ST_DIST(6) = "VOID from *.DST file, [-] "
      NAME_ST_DIST(7) = "POWER from *.DST file, [-]"

      NAME_TR_DIST(1) = "POWER, [Wt/cm^3] "
      NAME_TR_DIST(2) = "BORON CONCENTRATION, [ppm]"
      NAME_TR_DIST(3) = "COOLANT TEMPERATURE, [C]"
      NAME_TR_DIST(4) = "COOLANT DENSITY, [g/cm^3]"
      NAME_TR_DIST(5) = "DOPPLER FUEL TEMPERATURE, [K]"

      n_tr_dist = 1 + N_FEEDBACK 

      IF(TH_Model.EQ."Internal".OR.TH_Model.EQ."SKAZKA") THEN
          n_tr_dist = n_tr_dist +  N_OUT_TH_DIST 
          NAME_TR_DIST(6) = "FUEL CENTERLINE TEMPERATURE, [K]"
          NAME_TR_DIST(7) = "CLADDING INNER SURFACE TEMPERATURE, [K]"
          NAME_TR_DIST(8) = "FUEL ENTHALPY, [J/Kg]"   
      END IF


!      DO j = 1, NN_FRD_FA
!        WRITE(NAME_TR_DIST(6+(j-1)*N_OUT_TH_DIST), '(A, I3, A, I3)') &
!         "FUEL CENTERLINE TEMPERATURE, [K], Fuel Rod  ", j, " of ",&
!        NN_FRD_FA   
!        WRITE(NAME_TR_DIST(7+(j-1)*N_OUT_TH_DIST), '(A, I3, A, I3)') &
!         "CLADDING INNER SURFACE TEMPERATURE, [K], Fuel Rod  ", &
!         j, " of ",   NN_FRD_FA   
!       WRITE(NAME_TR_DIST(8 +(j-1)*N_OUT_TH_DIST), '(A, I3, A, I3)')&
!       "FUEL ENTHALPY, [J/Kg], , Fuel Rod  ", &
!         j, " of ",   NN_FRD_FA   
!      END DO

! burnup
!      write(*,*) 'Before burnup, problem_type=', TRIM(Problem_Type)
!      PAUSE

      IF( INDEX(Problem_Type,"Burnup") /= 0 ) THEN 
!        write(*,*) 'inside'
!        PAUSE
        n_tr_dist = n_tr_dist + 1
!        write(*,*) 'n_tr_dist'
        NAME_TR_DIST(n_tr_dist)="Burnup, [MWt days/kgHM]"
! isotopes
        DO i = 1, N_ISOTOPE
        NAME_TR_DIST(n_tr_dist+i)=isotope_name(i)
        END DO
      n_tr_dist = n_tr_dist + N_ISOTOPE
      END IF

      DO n = 1, NG 
         WRITE(NAME_TR_DIST(n+n_tr_dist), '(A,I3)') &
           "Neutron Flux, group ", n
      END DO
      n_tr_dist = n_tr_dist + NG

!      write(*,*) 'n_tr_dist=', n_tr_dist
!      DO i = 1, n_tr_dist 
!      WRITE(*,*) 'Names =', TRIM(NAME_TR_DIST(i)) 
!      END DO 
!      pause

      
      NAME_TR_SCAL(1) = "Total Reactor Power, [MWt]"
      NAME_TR_SCAL(2) = "Reactivity, [$]"
      NAME_TR_SCAL(3) = "Time Step Size of the SKETCH code, [s]"
      NAME_TR_SCAL(4) = "Time Step Size proposed by the TRAC code, [s]"
      NAME_TR_SCAL(5) = "Computed pressure error (SKAZKA module)"
      NAME_TR_SCAL(6) = "Computed enthalpy error (SKAZKA module)"
      NAME_TR_SCAL(7) = "k_{eff}"
        NAME_TR_SCAL(8) = &
         "Axial Power Offset (p_down - p_up)*100./(p_up+p_down)"
        NAME_TR_SCAL(9) = &
         "Axial Power Offset I-135 (p_down - p_up)*100./(p_up+p_down)"
        NAME_TR_SCAL(10) = &
         "Axial Power Offset Xe-135 (p_down - p_up)*100./(p_up+p_down)"
        NAME_TR_SCAL(11) = &
         "Average Reactor Pressure, Pa"


        DIM_TR_VEC(1) = NN_CRoD
        NAME_TR_VEC(1) = "Control Rod Positions, [cm]"

      RETURN
      END
