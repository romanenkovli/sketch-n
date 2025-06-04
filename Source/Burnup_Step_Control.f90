      subroutine Burnup_Step_Control(time_sk, i_ssor, i_source, &
      i_nonl_tot,  i_therm, i_trac_end, i_sk_end, tid_parent)
!=====================================================================*
!        main subroutine of Automatic Time Step Control of the        *
!     the SKETCH-N Code ( Neutron Kinetics  Calculations )            *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
!        USE core_history
!        USE CR_POSITION

      implicit none
      include 'sketch.fh'
! Input: time, i_int, i_out, i_nonl 
!        dt_sketch - proposed value of the time step size (really 
!                     dt = 0.5*dt_sketch
      integer i_int, i_out, i_nonl, i_therm, i_therm_step
! Output: Time_sk - Current Time
!         i_int_a - Number of SSOR iterations
!         i_out_a -  Number of the Source iterations
!         i_nonl_a - Number of Nonlinear Iterations
!         dt_sketch  - Used value of the coarse mesh time step 
!         dt_sketch_new - Estimate of the next time step size
      real time_sk
!      real dt_sketch
      integer i_ssor, i_source, i_nonl_tot
! Local Variables
      integer  step_accept
! step_accept = 1 if time step was accepted
! i_time_fine = 1 if fine temporal mesh (dt = 0.5*dt_step)
      integer j, N_Step_Fine, i_trac_end, tid_parent, i_sk_end
      real dt_step, time_diff, dt_xesm
!      real*8 dt_sketch_send 

!   For the Coupling without Restart      
      Logical No_Restart
!      parameter (NO_Restart = .True.)


! we can not perform restart of the TRAC code
!                   as a result thre is no rejected time steps
!                   in this case 
      LOGICAL :: error_read
      REAL    :: dt_hist_input

! reading reactor history data
      CALL read_core_history(error_read, dt_hist_input)
      IF( error_read ) GO TO 1000
      dt_sketch = dt_hist_input

! End reading burnup history      
      if(TH_Model.EQ."External".OR.TH_Model.EQ."SKAZKA") then
         No_Restart = .True.
      else
         No_Restart = .False.
      end if

!      write(*,*) ' No_Restart =', No_Restart 
!      pause

      N_Step_Fine = 1

      i_ssor = 0       
      i_source = 0
      i_nonl_tot = 0
      step_accept = 0

      dt_step = dt_sketch

      if(i_auto.eq.1) then
          call TSC_Save_Data(time_sk)
          if(TH_Model.EQ."Internal") then
             call THM_Save_Data
          end if
      end if

! first two fine time steps 
!             i_time_fine = 1

! Updating Burnup and Isotope concentration for the next time step
      CALL BRN_Compute_Burnup(dt_step)
      CALL BRN_Compute_Burnup_Distribution

      IF( Xe_Sm_Model(1:1) == "t".OR.Xe_Sm_Model(2:2) == "t") THEN
! convert days to seconds
          dt_xesm = 24*3600*dt_step
        write(*,*) 'Update Xe transient'
          write(*,*) 'dt_xesm =', dt_xesm
!        pause  
        CALL XE_Sm_Transient_Concentrations(dt_xesm)
      END IF


             do j = 1, N_Step_Fine
               
               time_sk = time_sk + dt_step
               i_time = i_time + 1
               write(*,*) 'Time =',time_sk,' dt=',dt_step,&
                            'TIME STEP NUMBER =',i_time
!               write(*,*) 'Time =',time_sk24*3600,' dt=',dt_step24*3600,&
!                            'TIME STEP NUMBER =',i_time
!               read(*,*)
 
               CALL Burnup_Compute_Time_Step(time_sk,dt_step,i_int,&
                i_out,i_nonl, i_therm_step, tid_parent,&
                i_sk_end, i_trac_end)
               i_ssor = i_ssor + i_int
               i_source = i_source + i_out
               i_nonl_tot = i_nonl_tot + i_nonl
               i_therm = i_therm + i_therm_step

               if(abs(TTV(NP_VIEW)-time_sk) .LT. eps_round_off) &
                                                    i_sk_end = 1

            end do

!           write(*,*) 'PVM_SKETCH_Receive',  'i_trac_end =', i_trac_end

         if(i_auto.eq.1) then        
           call TSC_Estimate_Temporal_Error
           call TSC_Estimate_Time_Step_Size(time_sk, dt_step, &
            step_accept, No_Restart)
! if time step was rejected

             if(NO_Restart) then
               STEP_ACCEPT = 1
             end if
             
             if(step_accept.ne.1) dt_sketch = dt_step

         else
! there is no rejected time step in the case of the constant time step size
          step_accept = 1
!         dt_step = dt_sketch ! simply to prevent round_off errors
         end if


! Time Step Selection for the output moments
   
      if(time_sk.LT. (TTV(i_view) - eps_round_off)) then      
        i_view_out = 0
      else if(time_sk.lt.TTV(NP_VIEW) - eps_round_off) then
        i_view = i_view + 1 
        i_view_out = 1
      end if

      if(time_sk.LT. (TTV(NP_VIEW)-eps_round_off)) then
        time_diff = TTV(i_view) - time_sk
      else
        time_diff = BIG_VALUE
      end if

! In the case of automatic time step size control we can adjust
! the value of time step size if we want to get an output in 
! the certain time moment

      if(i_auto.ne.1) dt_step = dt_input(i_view)
! 1st minimum between SKETCH & TRAC            
        
        
      if(TH_Model.EQ."External") dt_step = amin1(dt_step, dt_trac)
! adjusting time step for output
!         write(*,*) 'dt_sketch (OLD),time_diff =',dt_sketch,time_diff 

         dt_sketch= time_diff/&
                (aint(time_diff/dt_step - eps_round_off)+1.)
!         write(*,*) 'dt_sketch (NEW) =', dt_sketch 
!         pause

!   allow to decrease the time step in the certain moments           
         if(i_view_out.eq.1) &
            dt_sketch = amin1(dt_input(i_view),dt_sketch)

!      else
!         dt_sketch = dt_input(i_view)
!         if(TH_Model.EQ."External") &
!                            dt_sketch = amin1(dt_sketch, dt_trac)
!      end if

      if(TH_Model.NE."External") then
         i_trac_end  = i_sk_end
      end if

!       write(*,*) 'time_diff, time_diff/dt_step =', time_diff, &
!                                                  time_diff/dt_step
        write(*,*)'TIME STEP SIZE =', dt_sketch, dt_step
!          'TE_Flux =', TE_FLUX
!       write(*,*) 'dt_step =', dt_step 
!     pause


!      write(*,*) 'i_ssor, i_source, i_nonl_tot =', &
!        i_ssor, i_source, i_nonl_tot

      GO TO 1001

 1000 i_sk_end=1 
      i_trac_end = 1
 1001 CONTINUE 


      return
      end SUBROUTINE Burnup_Step_Control


      subroutine BURnup_Compute_Time_Step(time,dt_step, i_ssor, &
         i_source, i_nonl_tot, i_therm, tid_parent,&
         i_sk_end, i_trac_end)
!=====================================================================*
!         Estimation of the  Time Step Size Using                     * 
! the Time Error Estimation (Time_error) & Given Tolerance (ST_EPS)   *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! input: Time, dt_step
!         i_time_fine =  1 if fine temporal mesh (real solution)
      real dt_step, Time
      integer i_flag_convergence
      integer i_int, i_out, i_nonl, i_therm, tid_parent, &
       i_sk_end, i_trac_end
!  i_int - number of internal iterations
!         i_out - number of the source iterations
!         i_nonl_tot - number of the nonlinear iterations
! Output: 
      integer i_ssor, i_source, i_nonl_tot
! Local Variables: 
      real dt_rods, time_rods, dt_kin, a_dt_core
      real*8 dt_sketch_send 
      real   index_not_used
!      real shape_norm, pow_scale
      index_not_used = 0

      dt_rods = dt_step
      time_rods = time

      if(Reactor_Type.eq."PWR") then        
          call CRD_Move_Rods(time, dt_rods)
!      else
!          call CRD_Move_Rods_BWR(time, dt_rods) ! NOT READY YET
      end if

      IF(CRD_Rod_Type.EQ."FUEL") THEN
! Move the delayed neutron precursors for VVER-440 control rod types
        CALL CRD_Fuel_Control_Rods
      END IF


      dt_kin = dt_step
      a_dt_core = 1.E-11

      i_flag_convergence = 0
      i_therm  = 0. 
      i_ssor = 0 
      i_source = 0
      i_nonl_tot = 0 
      DO WHILE ( i_flag_convergence /= 1 ) 

      call CTRL_Solver(i_out, i_int, i_nonl, i_flag_convergence, dt_kin)

      call POWer_Compute

!      WRITE(*,*) 'i_flag_convergence=', i_flag_convergence
!      pause

      i_ssor     = i_ssor     + i_int
      i_source   = i_source   + i_out
      i_nonl_tot = i_nonl_tot + i_nonl

      IF(TH_Model.EQ."External") then
!             write(*,*) 'Steady-State '
!            write(*,*) 'SEND IN'
             dt_sketch_send = DBLE(dt_kin)
             call PVM_SKETCH_Send(tid_parent,i_sk_end, dt_sketch_send)
!            write(*,*) 'SEND OUT'
             call PVM_SKETCH_Receive(tid_parent, i_sk_end, i_trac_end)
!            write(*,*) 'RECEIVE OUT'
      else if(TH_Model.EQ."Internal") then

            call CPU_Timer( time_sk_cpu )

            call THM_Set_Core_Power
            call THM_Compute_Time_Step( a_dt_core)
            call THM_Get_Core_Feedbacks
            call CPU_Timer( time_th_cpu )
            i_trac_end = i_sk_end

      else if(TH_Model.EQ."SKAZKA") then

            call CPU_Timer( time_sk_cpu )
            write(*,*) 'SKAZKA IN'
            call YH_Set_Core_Power
            call YH_Set_Inlet (time)
            call YH_one_step(a_dt_core,  data_scalar_skazka(1), &
                  data_scalar_skazka(2) )
            CALL YH_Get_Core_Feedbacks 
            CALL THM_Prepare_Output_Dist 
            write(*,*) 'SKAZKA OUT, ', data_scalar_skazka(1), &
                     data_scalar_skazka(2)  

            call CPU_Timer( time_th_cpu )
            i_trac_end = i_sk_end

!           pause
      else
            i_trac_end = i_sk_end
      end if

       IF(TH_Model.NE."None") then
        i_therm = i_therm + 1
       END IF

      END DO ! i_flag_convergence /= 1 

      write(*,*) 'Compute SVRK QL'
!        pause
!      CALL POW_Compute_SVRK_QL(time)

!      call PNT_Compute_Reactivity(dt_kin)

! Sm steady-state


      RETURN
      END


