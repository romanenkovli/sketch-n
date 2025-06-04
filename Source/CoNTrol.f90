      subroutine CTRL_INPut_neutron(i_sk_end,Time)
!=====================================================================*
!            input data and computing initial data                    *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE HOMOGENIZATION_XS
      implicit none
      include 'sketch.fh'
      integer i_sk_end
      real Time
      LOGICAL :: Error
!      real dt_sketch
!Local Variablse:
!      INTEGER n1, i

! Initialization of the Error Message File, overwritten if exists
      CALL MSC_ERR_Init_File


      if(DEBUG) write(*,*) 'Input par'
      call INPut_Read_Namelist

      CALL CNT_Input


      if(DEBUG) write(*,*) 'Input XS'
        call XS_Input
!!!        call XS_Input_Old

      if(DEBUG) write(*,*) 'Input Geometry'
        call GMT_Input

      if(DEBUG) write(*,*) 'Input Control Rod Data'
        call CRD_Input

      if(DEBUG) write(*,*) 'Input Data for XS Homogenization Data'
        call HOM_XS_Input  

      IF(Problem_Type.EQ."Burnup") THEN
      if(DEBUG) write(*,*) 'Input Data for burnup calculations'
        CALL BRN_Input 
      END IF

      if(DEBUG) write(*,*) ' XS_init'
        call XS_Init

      if(DEBUG) write(*,*) 'Initialization of the Geometry'
        call GMT_Set

      if(DEBUG) write(*,*) 'Initialization of the CRD Module'
        call CRD_Init


! set names for output
      CALL INPut_Set_Distr_Names


      if(TH_Model.EQ.'External') then

         if(DEBUG) write(*,*) 'Input Mapping Matrices for TRAC'
         call PVM_Input_Map_Matrices
  
!         if(DEBUG) write(*,*) 'Mapping Matrices Input Finished'
!         if(Debug) call PVM_DBG_Check_Map

      end if

      if(File_DIST.NE."") then 
            Call RDS_Input_Distrib
      end if

      if(File_CD.NE."".AND.XS_model.EQ."TABLE") then
        call XSR_Input_CD_FILE(File_CD)
        if(DEBUG) write(*,*) 'Reading XS Tables from DIST file'
      end if


      call INPut_Read_Restart_File(Time)


!      write(*,*) 'Set feddbacks  to 1000.'

!      do i = 1, N_FEEDBACK 
!         write(*,*) 'i=', i, (fdback(388,n1,i), n1=1, NZ)
!      end do
!      pause

      if(DEBUG) write(*,*) 'Init Parameters'

      call INPut_Set_Param(i_sk_end, Time)

      if(DEBUG) write(*,*) 'Init param '

! we do not make it anymore nuSF are normilized later
!      if(Problem_Type.EQ."Kinetics") call XS_Set_Criticality
!      if(Problem_Type.EQ."Kinetics") call XS_Set_Crit


      if(Problem_Type.EQ."Burnup")  THEN
           call BRN_Initialize 
           CALL set_init_core_history(Error)
      end if
      if(DEBUG) write(*,*) 'Input OUT'

      return
      end


      subroutine CTRL_Solver(i_source, i_ssor, i_nonl, i_sk_end, &
              dt_kin)       
!=====================================================================*
!        Iteration Procedure for the Steady-state & Kinetics          *
!                                                 Calculations        *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input, NG_BEG, NG_END, NG_Step, n_inter, e_inter      
      real dt_kin
      logical ADJOINT
! Output: i_source, i_ssor, i_nonl, i_sk_end
      integer i_source, i_ssor, i_nonl, i_sk_end
! Local Variables
      integer i_source_old, n
      logical Conv_Outers, First_Outer_Iteration, boron_completed
      real two_node_kef
      integer   i_out, i_int, n_int
      real  e_int 
 
!      real time_beg, time_end, tolerance

      integer NN
      real scaled_norm_res, dt_adj


      i_source_old = i_source

      i_ssor = 0
      i_source = 0
      i_nonl = 0
      ADJOINT = .False.

!        boron_completed = .False. 

           write(*,*) 'boron concentration =', fdback(1,1,1)
!           write(*,*) 'boron completed =', boron_completed  
           write(*,*)        'd_kef_l=', d_kef_l,  e_boron_start

!           pause


! Critical Boron Search
      if(Steady_State_Type.EQ."BoronSearch") then
        IF(d_kef_l .LE. e_boron_start) THEN
!           call EIS_Find_Critical_Boron_Conc
            CALL EIS_Find_Critical_Boron_Conc_NEW(boron_completed)
!            write(*,*) ' boron_completed =', boron_completed 
           write(*,*) 'boron concentration =', fdback(1,1,1)
           write(*,*) 'boron completed =', boron_completed  
          else
            boron_completed  = .false.
        END IF
      end if

! Xe steady-state
      IF( Xe_Sm_Model(1:1) == "s") THEN
!        write(*,*) 'Update Xe steady-state'
!        pause  
        CALL Xe_Update_Steady_State
!        write(*,*)
      END IF
! Sm steady-state
      IF( Xe_Sm_Model(2:2) == "s") THEN
        CALL Sm_Update_Steady_State
      END IF

! Compute new material composi
      call CRD_Compute_Mat_Comp
      
! Computing New Macro XS
      call XS_Update 

! Compute New Finite-Difference Matrix
      
! COMPUTE FINITE-DIFFERENCE COUPLING MATRIX ONLY FOR THE 1st OUTER
!     ITERATION or in the case of the finite-difference method  
!      First_Outer_Iteration=(i_source_old.eq.0).AND.&
!                  (file_dmp_in.eq."")
!       write(*,*) 'First_Outer_Iteration=', First_Outer_Iteration
!       pause
! IN THE CASE OF Kinetics calculations we take the matrix from the
! previous time step (neglecting the change of the diffusion coefficients)
!      IF(Problem_Type.EQ."Kinetics") First_Outer_Iteration = .False.
!      IF (Nodal_Method .EQ."MCFD".OR. First_Outer_Iteration ) THEN
!         write(*,*) 'Update FD Matrix'
!       pause
         call MAT_Update_FD_Matrix
!      END IF 

      if(Problem_Type.EQ."Kinetics".AND.Kinetics_Method.EQ."PNT") then
! computing total matrix and exit
        call MAT_Set_Total_Matrix
        go to 100
      end if

      if(Problem_Type.EQ."Kinetics") call KIN_Init_Time_Step(dt_kin)


   
! Compute New Nodal Coupling Coefficients and Total Matrix

   12 First_Outer_Iteration=(i_source_old.eq.0).AND.&
                  (file_dmp_in.eq."")
      
      First_Outer_Iteration = .False.
      
      if(Nodal_Method .NE."MCFD" .AND. (.NOT.First_Outer_Iteration) )&
          then

          i_nonl = i_nonl + 1
! recomputing TRL at the beginner of time step
!             call TRL_Compute_TRL

          if(Problem_Type.EQ."Kinetics") then
            two_node_kef = 1.
          else
            two_node_kef = k_ef
          end if

             if(Problem_Type.EQ."Kinetics") then
! computing the Source term to update KIN_TRL
               call EIS_Compute_Source(ADJOINT)
            end if

!                  n = 4
!                  write(*,*) 'f(85, 86) =', flux(n,85), flux(n, 86)
!                  write(*,*) 'f(99, 100, 101) =', &
!                       flux(n,99), flux(n,100), flux(n,101)
!                  write(*,*) 'f(114, 115) =', flux(n,114), flux(n,115)
!                  PAUSE


             call TRL_Update_TRL(dt_kin)

             call MAT_Solve_two_node(two_node_kef)
!             CALL MAT_Solve_Two_Node_ANM_OLD(two_node_kef)

             call CPU_Timer( time_nonl_cpu )

      end if ! nodal_method /= "MCFD"


      call MAT_Set_Total_Matrix


         if(Problem_Type.NE."Kinetics") then
! Changing Wieland Shift in the Case of STeady-State Calculations
             call EIS_Update_Wieland_Shift
         end if
 
         call MAT_Compute_Block_Diag(ADJOINT, dt_kin)

         if(Problem_Type.EQ."Kinetics".AND.Debug) then
!            write(*,*) 'Matrices are written in the unformatted file'
!            call MAT_KIN_Output_Dump
!            read(*,*) 
         end if


         if(NG.eq. 2) then
! analytical inversion in the case of 2 neutron energy groups
           call MAT_Inverse_2x2_Block_Diag
         else
           call MAT_Inverse_Block_Diag
         end if


!         call CPU_Time(time_beg)

! TMP             
!          CALL PNT_Compute_Reactivity(1.0)
!           CALL POWer_Compute
!           write(*,*) 'Power Total =', p_total
!           read(*,*)
! TMP

         if(Problem_Type.EQ."Kinetics") then 
           if( iter_solver.eq."CSA" .OR. iter_solver.eq."CSI" ) then 

! Compute Initial Flux  Approximation in the case of Neutron Kinetics
           call CHB_KIN_Init_Iterations(i_source)
! Computing the 1st approximation for the neutron kinetics problem
           call MAT_Set_Kin_Block_RHS

           e_int = e_inter
           n_int = n_inter
           
           call PCD_Block_SSOR(i_int, e_int, n_int)

           i_ssor = i_ssor + i_int

! calculation of the initial residual vector if f(0) = f(t)
           call CHB_Set_First_Residual

! zero Block RHS in the case of kinetics for the rest of calculations
           call MSC_SSET(NG*N_TOT, 0., MAT_Block_RHS)

           end if ! iter_solver == "CSA" or "CSI

         end if  ! Kinetics
 
         if( (Problem_Type.NE."Kinetics").OR. &
            (iter_solver.eq."CSA" .OR. iter_solver.eq."CSI") ) then 
             call CTRL_Block_Iterations(i_out, i_int, &
                                              Conv_Outers, ADJOINT)
            i_ssor = i_ssor + i_int

         else
             
             call CGSolver_RHS_set
             call CGSolver(e_flux_l, N_OUTER, i_out)
             call CGSolver_Flux_set
! 1 SSOR iterations 
             i_ssor = i_ssor + i_out

         end if

! TMP  CHECKING CONVERGENCE TESTS
         if(Problem_Type.EQ."Kinetics") THEN
            if( iter_solver.eq."CSA" .OR. iter_solver.eq."CSI" ) then            
              NN = NG*N_TOT
              call CGS_compute_scaled_residual(NN, flux, mat_rhs_k, &
               scaled_norm_res)              
              write(*,*) 'error_estimate =', d_flux_l, &
                'scaled res norm=', scaled_norm_res
            else
              NN = NG*N_TOT
!              call CGS_compute_scaled_residual(NN, flux, mat_rhs_k, &
!               scaled_norm_res)              
!              write(*,*)  'scaled norm of the residual =', &
!               scaled_norm_res
            end if
          end if
! END TMP  CHECKING CONVERGENCE TESTS

!         call CPU_Time(time_end)
!         write(*,*) 'CMFD Iterations take ', time_end - time_beg, ' s'     

         if(Problem_Type.EQ."Kinetics") then
            call EIS_Compute_Source(Adjoint) ! New Source term to update precursors
         end if


!      write(*,*) 'FSI OUT'      

      i_source = i_source + i_out

!      write(*,*) 'TRL IN'      

      call CPU_Timer( time_cmfd_cpu )

! Moving into TRL
!      call TRL_Compute_TRL

!      write(*,*) 'TRL OUT'      
      if(Problem_Type.EQ."Kinetics") then 
! performing several nonlinear iteration per time step
        if( (Nodal_Method .NE."MCFD" ) .and. (i_out.GT. N_OUT_MAX)) &
                  go to 12
      else
! setting termination criteria in the case of steady-state iterations
! in this case we alsways have 1 nonlinear/thermal-hydraulics iteration
       if(Conv_Outers.AND. (i_out.LE. N_OUT_MAX) ) THEN
        IF( Steady_State_Type.EQ."BoronSearch") THEN
          IF( boron_completed ) i_sk_end = 1
        ELSE ! Steady_State_Type /= "BoronSearch"
          i_sk_end = 1
          boron_completed = .false. 
          write(*,*) 'Convergence, i_sk_end, boron_completed' ,&
          i_sk_end, boron_completed, N_OUT_MAX
!          pause
         END IF ! Steady_State_Type.EQ."BoronSearch"
        END IF ! Conv_Outers.AND. (i_out.LE. N_OUT_MAX)

      end if ! Problem_Type.EQ."Kinetics"

      if(Problem_Type.NE."Kinetics") then
          write(*,'(" i_source = ",I8," k eff  = ",F8.5)')  &
                 i_source, k_ef
          write(*,'(" k_ef_max=  ",  F8.5, " k_ef_min=",F8.5)')&
                 k_ef_max,k_ef_min
          write(*,*) 'Convergence  of Local K_eff       ',&
                            d_kef_l,e_outer_l
          write(*,*) 'Convergence of Local Neutron Flux ',&
                            d_flux_l,e_flux_l

          flux_a(:,:) = 1.
          dt_adj = 0.1          
!          write(*,*) 'xs_al(2,1)=', xs_al(2,1)
!          pause
          IF( xs_al(2,1) == 0.) THEN
            xs_al(:,:) = 1.
          END IF
          call ADJ_Normalize_Flux(NG, N_TOT,  xs_al, &
                 Flux_a(1,1), Flux(1,1), Volume(1))
          call PNT_Compute_Reactivity(dt_adj)
          write(*,*) 'reactivity =', (1.-1./k_ef)/bet
!          pause 


      else
          write(*,'(" Number of outer iterations = ",I8)' )  &
                 i_source
!         pause
      end if

! Output Matrix in COO Format
!      if(Problem_Type.EQ."Kinetics") then
!         write(*,*) 'In  MAT_KIN_Output_CSR_Format(dt_kin) '
!         call MAT_KIN_Output_CSR_Format(dt_kin)
!         write(*,*) 'Output of the matrix in CSR format'
!      end if


  100 continue ! Exit for Point Kinetics Solution

      return
      end


      subroutine CTRL_Block_Iterations(i_source, i_ssor, &
                                        Conv_Outers, ADJOINT)
!=====================================================================*
!        Fission Source Iterations for the STEADY-STATE &             *
!           Neutron Kinetics Calculations  19 May  1998               *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: i_int, i_nonl, dt_kin
      integer i_int
      logical ADJOINT
! Output: i_source, i_ssor - Counters of the Inter and Outer Iterations
!  Conv_Outers - logical variables - Convergence of the Outer Iterations
      logical Conv_Outers
      integer i_ssor, i_source
! Local Variables
      integer n_int
      real e_int, a_norm

      i_ssor = 0.
      i_source = 0.

      Call CHB_Init_Iterations

      if(Problem_Type.NE."Kinetics") then 
        call EIS_Compute_Source(ADJOINT)
        call MSC_Get_Norm_1(N_TOT, Source, S_Norm)
      end if

      do while(i_source .LT. N_OUTER) ! IN a case of Kinetics  N_OUTER = 1000

         i_source = i_source + 1

         call CHB_Set_Iteration_Param(i_source)

         if(Problem_Type.NE."Kinetics") then
! Scale the Source 
             a_norm = 1. / eigenv
             call MSC_SSCALE(N_TOT, a_norm, Source)
! Initial Block RHS
             call MAT_Compute_Block_RHS(ADJOINT)
         end if


         e_int = e_inter
         n_int = n_inter

        call PCD_Block_SSOR(i_int,e_int,n_int)

        i_ssor = i_ssor + i_int

!          
        if(Problem_Type.EQ."Kinetics") then
           call CHB_Iterate_Transient 
        else
           call CHB_Iterate_Steady_State
           call EIS_Compute_Eigenvalue(Adjoint)
        end if

        call CHB_Estimate_Max_Eigenv(i_source)

        if(Problem_Type.EQ."Kinetics") then
           d_kef_l = 0.
           if(iter_solver.eq."CSA") then
              d_flux_l = deln_dp/((1.- xme_)*(1.- kin_k_ef))
           else if(iter_solver.eq."CSI") then
             d_flux_l = deln_fp/(1.- xme_)
           end if ! CSA
         end if ! Kinetics

!          write(*,*) 'Conergence of Outer Iterations =',Conv_outers
!          write(*,'(" i_source = ",I8," k eff  = ",F8.5)')  &
!                 i_source, k_ef
!          write(*,'(" k_ef_max=  ",  F8.5, " k_ef_min=",F8.5)')&
!                 k_ef_max,k_ef_min
!          write(*,*) 'Convergence  of Local K_eff       ',&
!                            d_kef_l,e_outer_l
!          write(*,*) 'Convergence of Local Neutron Flux ',&
!                            d_flux_l,e_flux_l
!          write(*, *) 'xme_ =', xme_, 'npolin =', npolin
            
           if((d_flux_l.le.e_flux_l).and.(abs(d_kef_l).le.e_outer_l))&
                                                                then
               Conv_Outers = .TRUE.
               go to 11
           end if

      end do

   11 continue

      if(Problem_Type.EQ."Kinetics") then
         call CHB_Extrapolate_CSA
      end if


      if(DEBUG) then
!         write(*,*) 'Conergence of Outer Iterations =',Conv_outers
!         write(*,'(" i_source = ",I8," k eff  = ",F8.5)')  
!    &            i_source, k_ef
!         write(*,'(" k_ef_max=  ",  F8.5, " k_ef_min=",F8.5)')
!    &            k_ef_max,k_ef_min
!         write(*,*) 'Convergence  of Local K_eff       ',
!    &                       d_kef_l,e_outer_l
!         write(*,*) 'Convergence of Local Neutron Flux ',
!    &                       d_flux_l,e_flux_l

!        write(*, *) 'xme_ =', xme_, 'npolin =', npolin
!        write(*, *) 'kin_k_ef =', kin_k_ef
!        pause
      end if

      return
      end
            
      subroutine CNT_Input
!=====================================================================*
!        Input from the file FILE_INPUT                                *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Local Variables
      integer  line

      integer ios

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      character*100 Message
      logical error_find

!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'


      open(io_unit,file = FILE_INPUT,status='old', iostat=ios)

      call Iostat_Error_Check&
      (ios,"Could not find the input FILE_INPUT file "//FILE_INPUT)

! reading the FILE_INPUT Header
      do line = 1, N_LINE_PROBLEM_TITLE
         read(io_unit,'(A80)') PROBLEM_TITLE(line)
      end do

! reading CNT_RCT_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CNT_RCT_TYPE", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CNT_RCT_TYPE in the FILE_INPUT,  '//&
         'set reactor type to "PWR" ')
      else
         read(io_unit,  fmt=*, iostat=ios) REACTOR_TYPE
         call Iostat_Error_Check&
      (ios,"Error in reading the reactor type "//&
       "under identifier CNT_RCT_TYPE from the FILE_INPUT") 
        IF( (REACTOR_TYPE.NE."PWR").AND.(REACTOR_TYPE.NE."BWR") ) THEN
          WRITE(Message, '(A)') " Reactor Type ="//REACTOR_TYPE//&
           " is not known in the code, shopuld be PWR or BWR"
          CALL MSC_ERR_Add_Message('ERROR', Message)
        END IF
      end if


! reading CNT_RCT_POWR
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CNT_RCT_POWR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier CNT_RCT_POWR in the FILE_INPUT,  '//&
         'reactor power is not defined ')
      else
         read(io_unit,  fmt=*, iostat=ios) p_reactor
         call Iostat_Error_Check&
      (ios,"Error in reading reactor power "//&
       "under identifier CNT_CNT_RCT_POWR from the FILE_INPUT") 
      end if

      CLOSE(io_unit)

      RETURN
      END

      subroutine CNT_Output_Problem(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
      INTEGER nl

! PROBLEM_TITLE

      DO nl = 1, N_LINE_PROBLEM_TITLE
         WRITE(unit, '(A)') PROBLEM_TITLE(nl)
      END DO

!     PROBLEM DESCRIPTION
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     PROBLEM DESCRIPTION                                "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

! Reactor Type
      write(unit, *)

      write(unit, '(A, A5)')&
       " Reactor Type: ", REACTOR_TYPE

! Reactor Power [MWt]:
      write(unit, '(A, E13.6)')&
       " Reactor Power [MWt]:", P_Reactor

! " Type of the Computing Problem"
      write(unit, '(A, A)')&
       " Type of the Computing Problem: ", Problem_Type

! " Thermal-Hydraulics Model"
      write(unit, '(A, A)')&
       " Thermal-Hydraulics Model: ", TH_Model

      RETURN
      END


      subroutine CNT_Output_Code_Options(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
      CHARACTER*80 Message

!     PROBLEM DESCRIPTION
      write(unit, *)
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     SKETCH-N OPTIONS                                   "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

      write(unit, *)

! Nodal Method
      if(Nodal_Method.EQ."MCFD") then
             write(Message,'(A)') &
                       'Mesh-Centered Finite-Difference Method'
      else if(Nodal_Method.EQ."SANM") then
             write(Message,'(A)') &
                       'Semi-Analytic Method (SANM)'
      else if(Nodal_Method.EQ."PNM") then
             write(Message,'(A)') &
                       'Polynomial Nodal Method (PNM)'
      else if(Message.EQ."ANM") then
             write(Message,'(A)') &
                 'Analytic Nodal Method (ANM)'
      else if(Nodal_Method.EQ."PNM1") then
             write(Message,'(A)') &
          'Polynomial Nodal Method (PNM) using Matrix Functions'
      end if

      write(unit, '(A, A)')&
       "   Space Discretization            : ", Message

      write(unit, '(A, A)')&
       "   Transverse Leakage Approximation: ", TRL_Approx

      write(unit, '(A, A)')&
       "   Nonlinear Iteration Procedure   : ", NonlinearIterations


!     
      IF(Problem_Type.EQ."Kinetics") THEN
! Kinetics_Method
        IF(Kinetics_Method.EQ."DRT") then
             write(Message,'(A)') &
                       'Direct Method, Fully-Implicit Scheme'
        else if(Kinetics_Method.EQ."IQS") then
             write(Message,'(A)') &
                       'Improved Quasi-Static Method (NOT VERIFIED!)'
        else if(Kinetics_Method.EQ."IQS") then
             write(Message,'(A)') &
                       'Point Kinetics Model'
        end if

      write(unit, '(A, A)')&
       "   Time Discretization : ", Message
! Linear System Solver

        IF(Iter_Solver.EQ."CSA") then
             write(Message,'(A)') &
             'Chebyshev Semi-Analytical Method (CSA)'
        else if(Iter_Solver.EQ."CSI") then
             write(Message,'(A)') &
              'Chebyshev Semi-Iterative Method (CSI)'
        else if(Iter_Solver.EQ."CG") then
             write(Message,'(A)') &
              'Conjugate-Gradient Method (CG)'
        else if(Iter_Solver.EQ."BCGSTAB") then
             write(Message,'(A)') &
              'BiConjugate-Gradient Method Stabilized (BICGSTAB)'
        else if(Iter_Solver.EQ."TFQMR") then
             write(Message,'(A)') &
              'Transpose-Free Quasi-Minimal Residual Method (TFQMR)'
        else if(Iter_Solver.EQ."FOM") then
             write(Message,'(A)') &
              "Arnoldi's Method (FOM)"
        else if(Iter_Solver.EQ."GMRES") then
             write(Message,'(A)') &
              'Generalized Minimumal Residual Method (GMRES)'
        end if
        write(unit, '(A, A)')&
       "   Linear System Solver : ", Message
      END IF
 

      RETURN
      END

      subroutine CNT_Output_IO_Files(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
!      CHARACTER*80 Message

!     PROBLEM DESCRIPTION
      write(unit, *)
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     INPUT/OUTPUT FILES:                                "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

      write(unit, *)

! General Input Data
      write(unit, '(A, A)')&
       "   Input General Data: ", File_Input
! XS Table File
      IF(FILE_CD.NE."") THEN
      write(unit, '(A, A)')&
       "   Input File with XS Tables: ", File_CD
      END IF
! File with Distributions
      IF(FILE_DIST.NE."") THEN
      write(unit, '(A, A)')&
       "   Input File with Distributions : ", File_DIST
      END IF
! Input File with Mapping Matrices
      IF(FILE_Map.NE."") THEN
      write(unit, '(A, A)')&
       "   Input File with Mapping Matrices: ", File_Map
      END IF
! Input Restart File 
      IF(FILE_DMP_IN.NE."") THEN
      write(unit, '(A, A)')&
       "   Input Restart File : ", File_dmp_in
      END IF

! Output Restart File for steady-state calculation
      IF(FILE_DMP_OUT_ST.NE."") THEN
      write(unit, '(A, A)')&
       "   Output Steady-State Restart File: ", &
           File_dmp_out_st
      END IF

! Output Restart File for steady-state calculation
      IF(FILE_DMP_OUT_KIN.NE."") THEN
      write(unit, '(A, A)')&
       "   Output Kinetics Restart File: ", &
           File_dmp_out_KIN
      END IF
       
      write(unit,*)

      RETURN
      END
  
      subroutine CNT_Output_Conv_Criteria(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
!      CHARACTER*80 Message

!     PROBLEM DESCRIPTION
      write(unit, *)
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     ITERATIVE CONVERGENCE CRITERIA                     "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

      write(unit, *)

! "   Number of Outer Iterations per Nonlinear (convergence):"
      write(unit, '(A, I13)')&
       "   Number of Outer Iterations per Nonlinear (convergence) :", &
          N_Out_Max
! "   Number of Outer Iterations per Nonlinear (in calculation):"
      write(unit, '(A, I13)')&
       "   Number of Outer Iterations per Nonlinear (calculation) :", &
          N_Outer
! 
      IF(Steady_State_Type.EQ."BoronSearch") THEN
      write(unit, '(A, E13.5)')&
       "   Eigenvalue Accuracy when Critical Boron Search Starts  :", &
          E_Boron_Start
      write(unit, '(A, E13.5)')&
       "   Eigenvalue Accuracy when Critical Boron Search Finished:", &
          E_critical
      END IF

      IF(Problem_Type.NE."Kinetics") THEN
      write(unit, '(A, E13.5)')&
       "   Convergence Criteria for Eigenvalue                    :", &
          E_OUTER_L
      END IF

      write(unit, '(A, E13.5)')&
       "   Convergence Criteria on Neutron Flux                   :", &
          E_flux_l

      write(unit, '(A, I13)')&
       "   Number of Inner Iterations per Outer                   :", &
          N_Inter

      write(unit, '(A, E13.5)')&
       "   Convergence Criteria for Inner Iterations              :", &
          E_Inter

       
      write(unit,*)

      RETURN
      END

      subroutine CNT_Output_Iter_Parameters(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
!      CHARACTER*80 Message

!     PROBLEM DESCRIPTION
      write(unit, *)
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     ITERATIVE PARAMETERS                               "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

      write(unit, *)

! Estimate of Minimum Eigenvalue of the Iteration Matrix
      write(unit, '(A, E13.5)')&
       "   Estimate of Minimum Eigenvalue of the Iteration Matrix  :", &
          xbe
! Estimate of Maximum Eigenvalue of the Iteration Matrix
      write(unit, '(A, /, A,  E13.5)')&
       "   Estimate of Maximum Eigenvalue or Dominance Ratio",&
       "                                  of the Iteration Matrix  :", &
          Xme_ini
! Adaptivity Parameter of the Chebyshev Acceleration 
      write(unit, '(A,   E13.5)')&
       "   Adaptivity Parameter of the Chebyshev Acceleration      :",&
          F_Cheb
! Number of Outer Iteration when Chebyshev Acceleration Starts
      write(unit, '(A,   I13)')&
       "   Outer Iteration Number when Chebyshev Acceleration Start:",&
          NPolins
! 
      IF(Problem_Type.NE."Kinetics") THEN
      write(unit, '(A, E13.5)')&
       "   Inverse Value of the Wieland Shift                      :",&
        Delta_Shift 
      END IF

      write(unit,*)

      RETURN
      END

      subroutine CNT_Output_Time_Step_Control(unit)
      implicit none
      include 'sketch.fh'

! Input: 
      INTEGER unit
! Local
!      CHARACTER*80 Message
      INTEGER np

!     PROBLEM DESCRIPTION
      write(unit, *)
      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(A80)') &
        "!     TIME STEP CONTROL                                  "//&
        "                      !"
      CALL OUTput_Write_Separator(unit)

      write(unit, *)

      IF(I_AUTO.EQ.1) THEN
      write(unit, '(A)')&
       "   Automatic Time Step Size Control"
      ELSE
      write(unit, '(A)')&
       "   Constant Time Step Size in Time Intervals"
      END IF
! "   Maximum Number of the Time Step Size Intervals          :"
      write(unit, '(A, I13)')&
       "   Maximum Number of the Time Step Size Intervals          :",&
          NN_OUT_VIEW_MAX
! "   Number of the Time Step Size Intervals          :"
      write(unit, '(A, I13)')&
       "   Number of the Time Step Size Intervals                  :",&
          NP_View
!"   Time Intervals:                                         :"
      write(unit, '(A)')&
       "   Time Intervals:                                         :"
      write(unit, '(8x, 5E13.5)')&
          (TTV(np), np = 1, NP_View)

! "   Time Step Size inside Intervals:                        :"
      write(unit, '(A)' )&
       "   Time Step Size inside Intervals:                        :"
      write(unit, '(8x, 5E13.5)')&
          ( DT_INPUT(np), np = 1, NP_View)

! Accuracy Criterion of the Automatic Time Step Control
      IF(I_Auto.EQ.1) THEN
      write(unit, '(A, E13.5)')&
       "   Accuracy Criterion of the Automatic Time Step Control   :",&
          St_Eps

      write(unit, '(A, E13.5)')&
       "   Maximum Increase of the Time Step Size                  :",&
          Facmax
      END IF

      write(unit, '(A, E13.5)')&
       "   Maximum Time Step Size                                  :",&
          Dt_Step_Max
         
      write(unit, '(A, I13)')&
       "   Output into GRF file at the Nth Time Step               :",&
          N_Zap

      write(unit,*)

      RETURN
      END
