      SUBROUTINE sketch_init
      USE MAIN_VAR
      implicit none
      include 'sketch.fh'

! Initial 
      i_source = 0
      i_ssor = 0
      i_nonl_tot = 0
      i_therm = 0

! Initialization of Computing Time
      call CPU_Timer( time_sk_cpu)
      time_sk_cpu = 0.
      time_th_cpu = 0.
      time_XS_cpu = 0.
      time_nonl_cpu = 0.
      time_cmfd_cpu = 0.
! End of Initialization of Computing Time

! initialization of neutronics data

      call CTRL_Input_Neutron(i_sk_end, time_sk)


!     THERMAL HYDRAULICS MODEL 

      if (TH_Model.EQ."Internal") then

         call CPU_Timer(time_sk_cpu)

         call THM_Init(FILE_INPUT, file_dmp_in,&
                NZR_Core_Beg, NZR_Core_End, cool_heating, hz )

         call CPU_Timer(time_th_cpu)

      else if (TH_Model.EQ."SKAZKA") then

!        write(*,*) 'SKAZKA Input IN'

        call CPU_Timer(time_sk_cpu)

        call YH_Input_Data(FILE_INPUT)
        IF(file_dmp_skazka_in.ne."") THEN
           CALL YH_read_dynvar(file_dmp_skazka_in)
        ELSE
            call yh_initial_dynvar
        END IF                      
        CALL YH_Get_Core_Feedbacks 


        call CPU_Timer(time_th_cpu)

!        write(*,*) 'SKAZKA Input OUT'

      else if (TH_Model.EQ."Athlet") then

        write(*,*) 'Athlet'

        call CPU_Timer(time_sk_cpu)

        CALL ATHlet_Get_Core_Feedbacks

        call CPU_Timer(time_th_cpu)

!        write(*,*) 'SKAZKA Input OUT'


      else if (TH_Model.EQ."External") then
           call PVM_SKETCH_Enroll(tid_parent)
           call PVM_SKETCH_Receive(tid_parent, i_sk_end, i_trac_end)
           if(dt_trac.lt.dt_sketch) dt_sketch = dt_trac ! 1st time_sk step determined by SKETCH

      else if (TH_Model.EQ."Simulator") then

      end if

! END OF    THERMAL HYDRAULICS MODEL 

        i_trac_end = i_sk_end

!     END INPUT      
!

!  OUTPUT of the ALL DATA for all models except SIMULATOR      

      IF( TH_Model /= "Simulator" ) THEN
       
        CALL OUTput_Problem_Description

        CALL OUTput_NumResults(i_source,i_ssor, i_nonl_tot,&
                    i_therm, time_sk, 1 )
! Output into SKETCH.GRF file
         CALL GRaPhics_Write_INI_Data
         IF(Problem_Type.EQ."Kinetics".OR.Problem_Type.EQ."Burnup") &
            THEN
             CALL GRaPhics_Write_Distr_Data(time_sk)
          END IF

      END IF ! IF( TH_Model /= "Simulator" )

      call CPU_Timer(time_output_cpu)

! END  OUTPUT of the ALL DATA for all models except SIMULATOR      


      RETURN
      END  SUBROUTINE sketch_init
    
      SUBROUTINE sketch_compute_time_step(dt_sim)
      USE MAIN_VAR
      implicit none
      include 'sketch.fh'
      
      REAL :: dt_sim  


      IF( TH_Model == "Simulator") THEN
        dt_sketch = dt_sim
      END IF 

      if(TH_Model.NE."None") then
            i_therm = i_therm + 1
      end if


!      IF( TH_Model /= "Simulator") THEN
        WRITE(*,*) 'TEMPERATURE ITERATION NUMBER, time_sk = ',&
                      i_therm,  time_sk
!      END IF  
!         write(*,*)
!         write(*,*) 'NEUTRON CALCULATION'
! solution of the neutron diffusion equations
      IF (Problem_Type.EQ."Kinetics") then

!             CALL CPU_time_sk ( time_sk_begin )
              call Time_Step_Control(time_sk,i_int,i_out,i_nonl, &
               i_therm, i_trac_end, i_sk_end, tid_parent)

      ELSE IF(Problem_Type.EQ."Burnup") THEN

              call Burnup_Step_Control(time_sk,i_int,i_out,i_nonl, &
               i_therm, i_trac_end, i_sk_end, tid_parent)
      
      ELSE  
!            write(*,*) 'SOLVER IN'
             dt_sketch = 0.
             call CTRL_solver(i_out,i_int,i_nonl,i_sk_end, dt_sketch)
!            write(*,*) 'SOLVER OUT'

!            write(*,*) 'POWER IN'
             call POWer_Compute
             a_dt_core = 1.E-20
!            write(*,*) 'POWER OUT'

      END IF

          i_ssor  = i_ssor + i_int
          i_source = i_source + i_out
          i_nonl_tot = i_nonl_tot + i_nonl

          write(*,*) 'Power Total =', p_total

      IF(Problem_Type.EQ."Kinetics".OR.Problem_Type.EQ."Burnup") &
            then

! Computing average feedbacks for output
            CALL THM_Compute_Average_Feedbacks

            CALL CPU_Timer( time_sk_cpu )

            if(mod(i_therm,n_zap).eq.0) then
             CALL GRaPhics_Write_Distr_Data(time_sk)
            end if
            if(i_view_out.eq.1) then
             CALL OUTput_NumResults(i_source,i_ssor,&
                    i_nonl_tot, i_therm, time_sk, 0)
            end if
            CALL CPU_Timer( time_output_cpu )
      end if ! IF(Problem_Type.EQ."Kinetics".OR.Problem_Type.EQ."Burnup")


      IF(Problem_Type.EQ."Steady-State") then
           IF(TH_Model.EQ."External") then
!             write(*,*) 'Steady-State '
!            write(*,*) 'SEND IN'
             dt_sketch_send = DBLE(dt_sketch)
             call PVM_SKETCH_Send(tid_parent,i_sk_end, dt_sketch_send)
!            write(*,*) 'SEND OUT'
             call PVM_SKETCH_Receive(tid_parent, i_sk_end, i_trac_end)
!            write(*,*) 'RECEIVE OUT'
           ELSE IF(TH_Model.EQ."Internal") then

            call CPU_Timer( time_sk_cpu )

            call THM_Set_Core_Power
            call THM_Compute_Time_Step( a_dt_core)
            call THM_Get_Core_Feedbacks
            call CPU_Timer( time_th_cpu )
            i_trac_end = i_sk_end

           ELSE IF(TH_Model.EQ."Athlet") then

            call CPU_Timer( time_sk_cpu )
            call Athlet_Set_Core_Power
            call CPU_Timer( time_th_cpu )
            i_trac_end = i_sk_end


           ELSE IF(TH_Model.EQ."SKAZKA") then

            call CPU_Timer( time_sk_cpu )
!            write(*,*) 'SKAZKA IN'
            call YH_Set_Core_Power
!            write(*,*) 'time_sk =', time_sk
!            pause 
            call YH_Set_Inlet ( time_sk )
            write(*,*) 'a_dt_core =', a_dt_core,  data_scalar_skazka(1),&
            data_scalar_skazka(2)
            call YH_one_step(a_dt_core,  data_scalar_skazka(1), &
                  data_scalar_skazka(2) )
            CALL YH_Get_Core_Feedbacks 
            CALL THM_Prepare_Output_Dist 
!            write(*,*) 'SKAZKA OUT, ', data_scalar_skazka(1), &
!                     data_scalar_skazka(2)  

            call CPU_Timer( time_th_cpu )
            i_trac_end = i_sk_end

!           pause
           ELSE
            i_trac_end = i_sk_end
           END IF ! IF(TH_Model.EQ."External")
 
      END IF  ! IF(Problem_Type.EQ."Steady-State")

      RETURN 
      END SUBROUTINE sketch_compute_time_step


      SUBROUTINE sketch_end_remarks
      USE MAIN_VAR
      USE HOMOGENIZATION_XS 
      implicit none
      include 'sketch.fh'
      REAL :: time_st

      time_st = 0.0   

      IF(Problem_Type.NE."Kinetics".AND.Problem_Type.NE."Burnup") THEN
      write(*,*) 'Compute SVRK QL'
!      CALL POW_Compute_SVRK_QL(time_st)
!      pause
      END IF

! Neutron Diffusion Calculations
      call CPU_Timer( time_sk_cpu )

      IF(Problem_Type.NE."Kinetics")  THEN
!      .AND.Problem_Type.NE."Burnup") THEN
    write(*,*) 'KIN_Normalize_Flux '
           CALL KIN_Normalize_Flux 
         CALL OUTput_Compute_Average_Flux
      ! Computing the Surface Flux for Output
          IF( NonlinearIterations.EQ."Moon") THEN
            CALL HOM_Set_Surface_Flux
            CALL HOM_Average_Surface_Flux
          END IF  
      END IF

    write(*,*) 'THM_Compute_Average_Feedbacks'
      CALL THM_Compute_Average_Feedbacks
      CALL OUTput_NumResults(i_source,i_ssor, i_nonl_tot,&
                    i_therm, time_sk, 0)

      IF(Problem_Type.NE."Kinetics") THEN
        IF( FLAG_XS_HOPMOGENIZATION ) THEN
          CALL HOM_XS_Compute
          CALL HOM_XS_OUTPUT        
        END IF
      END IF ! (Problem_Type.NE."Kinetics") 


      IF(Problem_Type.NE."Kinetics".AND.Problem_Type.NE."Burnup") THEN
! Output into SKETCH.GRF file
         time_out = 0.
    write(*,*) 'GRaPhics_Write_Distr_Data'
         CALL GRaPhics_Write_Distr_Data(time_out)
      END IF

      IF(FILE_DMP_OUT_KIN.NE."") then
        if(Problem_Type.NE."Kinetics") then
    write(*,*) 'KIN_Initialize'
                call KIN_Initialize(time_sk)
         else
                dt_save = dt_sketch
         end if      
      end if

    write(*,*) 'ADJ_Normalize_Flux'
      call ADJ_Normalize_Flux(NG, N_TOT,  xs_al, Flux(1,1), &
                 Flux_A(1,1), Volume (1) )

    write(*,*) 'OUTput_Write_Restart_File(time_sk)'
      call OUTput_Write_Restart_File(time_sk)

    write(*,*) 'yh_write_dynvar(file_dmp_skazka_out) '
      IF(TH_Model.EQ."SKAZKA".AND.file_dmp_skazka_out.ne."") &
            CALL yh_write_dynvar(file_dmp_skazka_out) 
    write(*,*) 'yh_write_dynvar(file_dmp_skazka_out) end '

      if(Problem_Type.NE."Kinetics".AND. File_Reference.NE."" ) then
! A comparison of the SKETCH-N results with the reference solution
             write(*,*) 'reference comparison'
    write(*,*) 'OUTput_Reference_Comparison'
             call OUTput_Reference_Comparison(i_nonl_tot, i_source)
      end if

      call CPU_Timer( time_output_cpu )
      Call CPU_Timer_Output


!         call MAT_Output_COO_Format
!         write(*,*) 'Matrices are written in COO format'


!         if(Problem_Type.NE."Kinetics".AND. File_DMP_Out_Kin.NE."") &
!             call ADJ_Compute

 
         write(*,*) ' SKETCH-N Finished Calculations '
          
         if(TH_Model.EQ."External") then
               call PVM_SKETCH_Exit(tid_parent)
         end if

         write(*,*) 'end remarks end '


      RETURN
      END SUBROUTINE sketch_end_remarks

