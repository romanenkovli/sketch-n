      subroutine Time_Step_Control(time_sk, i_ssor, i_source, &
      i_nonl_tot,  i_therm, i_trac_end, i_sk_end, tid_parent)
!=====================================================================*
!        main subroutine of Automatic Time Step Control of the        *
!     the SKETCH-N Code ( Neutron Kinetics  Calculations )            *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

      implicit none
      include 'sketch.fh'
! Input: time, i_int, i_out, i_nonl 
!        dt_sketch - proposed value of the time step size (really 
!                     dt = 0.5*dt_sketch
      integer i_int, i_out, i_nonl, i_therm
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
      integer  step_accept, i_time_fine
! step_accept = 1 if time step was accepted
! i_time_fine = 1 if fine temporal mesh (dt = 0.5*dt_step)
      integer j, N_Step_Fine, i_trac_end, tid_parent, i_sk_end
      real dt_step, time_diff, a_dt_core
      real*8 dt_sketch_send 

!   For the Coupling without Restart      
      Logical No_Restart
!      parameter (NO_Restart = .True.)

! we can not perform restart of the TRAC code
!                   as a result thre is no rejected time steps
!                   in this case 

      
      if(TH_Model.EQ."External".OR.TH_Model.EQ."SKAZKA") then
         No_Restart = .True.
      else
         No_Restart = .False.
      end if

!      write(*,*) ' No_Restart =', No_Restart 
!      pause

      IF ( i_auto.eq.1 ) THEN
        N_Step_Fine = 2
      ELSE IF ( i_auto.eq.0 ) THEN
        N_Step_Fine = 1
      END IF      

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


      do while(step_accept.ne.1) 

         if(i_auto.eq.1) then
! 1st coarse mesh step size dt_step
             dt_step = dt_step*2.
             time_sk = time_sk + dt_step
             i_time = i_time + 1

             write(*,*)'TIME=',time_sk,'dt =',dt_step,&
                            'TIME STEP NUMBER=', i_time

             i_time_fine = 0
             call KIN_Compute_Time_Step(time_sk, dt_step, i_int, &
               i_out, i_nonl, i_time_fine)
 
             i_ssor = i_ssor + i_int
             i_source = i_source + i_out
             i_nonl_tot = i_nonl_tot + i_nonl

! save neutron flux for the error estimation
             call TSC_Save_Solution
! restart initial parameters in time = time_0
             call TSC_Recover_Data(time_sk)

! time step size for the next 2 time step
             dt_step = dt_step*0.5
           end if ! if(i_auto.eq.1) 

! first two fine time steps 
             i_time_fine = 1

             do j = 1, N_Step_Fine
               
               if( (j.EQ.N_Step_Fine) .AND. (i_auto==1) ) then 
! already did for the  first time step
                  if(TH_Model.EQ."External") then
! 
                     call PVM_SKETCH_RECEIVE&
                             (tid_parent, i_sk_end, i_trac_end)

! we should not change the time step in the middle 
!c                            if(dt_trac .LT. dt_sketch) dt_sketch = dt_trac
                  else if(TH_Model.EQ."Internal") then
                      call THM_Get_Core_Feedbacks
                  else if(TH_Model.EQ."SKAZKA") then
                      CALL YH_Get_Core_Feedbacks 
                      CALL THM_Prepare_Output_Dist
                  end if
               end if

               time_sk = time_sk + dt_step
               i_time = i_time + 1
               write(*,*) 'Time =',time_sk,' dt=',dt_step,&
                            'TIME STEP NUMBER =',i_time
               call KIN_Compute_Time_Step(time_sk,dt_step,i_int,i_out,&
                                        i_nonl, i_time_fine)
               i_ssor = i_ssor + i_int
               i_source = i_source + i_out
               i_nonl_tot = i_nonl_tot + i_nonl


              if(TH_Model.NE."None") then
                 i_therm = i_therm + 1
              end if

               if(abs(TTV(NP_VIEW)-time_sk) .LT. eps_round_off) &
                                                    i_sk_end = 1

               if(TH_Model.EQ."External") then

                 dt_sketch_send = DBLE(dt_step)
                 call PVM_SKETCH_Send(tid_parent,i_sk_end, &
                                              dt_sketch_send)
               
               else if(TH_Model.EQ."Internal") then
               
                 a_dt_core = 1./dt_step
                 call THM_Set_Core_Power

                 call CPU_Timer( time_sk_cpu )
                    call THM_Compute_Time_Step( a_dt_core)
                 call CPU_Timer( time_th_cpu )

              else if(TH_Model.EQ."SKAZKA") then

                call CPU_Timer( time_sk_cpu )
!                write(*,*) 'SKAZKA IN'
                 a_dt_core = 1./dt_step
                 call YH_Set_Core_Power
                 call YH_Set_Inlet(time_sk)
                 call YH_one_step(a_dt_core,  data_scalar_skazka(1), &
                  data_scalar_skazka(2) )
            write(*,*) 'SKAZKA OUT1, ', data_scalar_skazka(1), &
                     data_scalar_skazka(2)  
                call CPU_Timer( time_th_cpu )

!           pause


            end if

            end do ! do j = 1, N_Step_Fine

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
! end of the time stepping procedure
      end do

! receiving the message from TRAC
! already did it for the  first time step

      if(TH_Model.EQ."External") then

        call PVM_SKETCH_Receive(tid_parent, i_sk_end, i_trac_end)

        if(dt_trac .LT. dt_step) dt_step = dt_trac
      else if(TH_Model.EQ."Internal") then
        call THM_Get_Core_Feedbacks
      else if(TH_Model.EQ."SKAZKA") then
        CALL YH_Get_Core_Feedbacks 
        CALL THM_Prepare_Output_Dist 
      end if

! Time Step Selection for the output moments
   
      if(time_sk.LT. (TTV(i_view) - eps_round_off)) then      
        i_view_out = 0
      else if(time_sk.lt.TTV(NP_VIEW) - eps_round_off) then
        i_view = i_view + 1 
        i_view_out = 1
      end if

!      if(time_sk.LT. (TTV(NP_VIEW)- SMALL_VALUE)) then
        time_diff = TTV(i_view) - time_sk
!      else
!        time_diff = BIG_VALUE
!      end if

! In the case of automatic time step size control we can adjust
! the value of time step size if we want to get an output in 
! the certain time moment

      if(i_auto.ne.1) dt_step = dt_input(i_view)
! 1st minimum between SKETCH & TRAC            
        
      if(TH_Model.EQ."External") dt_step = amin1(dt_step, dt_trac)
! adjusting time step for output
!          dt_sketch = dt_step
!         dt_sketch= 0.5*time_diff/&
!                (aint(0.5*time_diff/dt_step - eps_round_off)+1.)

         IF( time_diff <= SMALL_VALUE ) THEN
           dt_sketch = dt_step
         ELSE
           IF( i_auto == 1 ) THEN
           dt_sketch= 0.5*time_diff/&
                (CEILING(time_diff/(2*dt_step) - SMALL_VALUE))
           ELSE
           dt_sketch= time_diff/&
                (CEILING(time_diff/(dt_step) - SMALL_VALUE))
           END IF 
         END IF
  
!         write(*,*)  '(time_diff', time_diff
!         write(*,*)  '(0.5*time_diff/dt_step - SMALL_VALUE)'
!         write(*,*)  (0.5*time_diff/dt_step - SMALL_VALUE)
!         write(*,*)  '(CEILING(0.5*time_diff/dt_step - SMALL_VALUE))'
!         write(*,*)  (CEILING(0.5*time_diff/dt_step - SMALL_VALUE))
!         write(*,*) 'dt_sketch, dt_step =', dt_sketch, dt_step
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

      return
      end


      subroutine TSC_Save_Data(time)
!=====================================================================*
!     Saving the Model Parameters before the Coarse Time Step         *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Time, Source(N_TOT), Flux(NG, N_TOT), TRL_XYZ(NG, N_TOT, NDIR),
!         Prec(MD, N_TOT), zrods(IRODS)
      real Time
! Output:   Old_Time - Time
!           Old_Flux(NG,N_TOT) - Neutro Flux
!           Old_TRL(NG, N_TOT, NDIR) -  Transverse Leakages
!           Old_prec(MD, N_TOT) - Precursor Concentration
!           Old_Source(N_TOT) - Source Term 
!           Old_Zrods(IRODS) - Control Rod Position
!           Old_Flux_Mom(MOMENTS, NG, N_TOT, NDIR) - Flux Moments
!Local Variables:
      integer k, n, m,  nd, n1

      do k = 1, N_TOT
         Old_Source(k) = Source(k)
         do n = 1, NG
            Old_Flux(n, k) = Flux(n, k)
         end do
         do m = 1, MD
            Old_Prec(m, k) =  Prec(m, k)
         end do
      end do

      Old_Time = Time

      do nd = 1, NDD
         do k = 1, N_TOT
            do n = 1, NG
               Old_TRL(n, k, nd) = TRL_XYZ(n, k, nd)
            end do
         end do
      end do

      call CRD_Save_Data


       Old_React = React 
       Old_Pow_Point = Pow_point
       do m = 1, MD
          Old_Prec_Point(m) = Prec_point(m)
       end do

       DO n1 = 1, NZ
          DO k = 1, NH
            Old_p(k,n1) = p(k,n1)
          END DO
      END DO  

       DO n1 = 1, NZ
          DO k = 1, NH
          DO  m = 1, MH
            OLD_p_dh(m,k,n1) = p_dh(m,k,n1)
          END DO
          END DO
      END DO  


!NOT USED      do nd = 1, NDIR
!NOT USED        do k = 1, N_TOT
!NOT USED          do m = 1, MOMENTS
!NOT USED            do n = 1, NG
!NOT USED               Old_Flux_Mom(m, n,  k, nd)  = Flux_Mom(m, n,  k, nd)
!NOT USED           end do
!NOT USED          end do
!NOT USED        end do
!NOT USED      end do

      return
      end

      subroutine TSC_Save_Solution
!=====================================================================*
!     Saving The Neutron Flux after the Coarse Time Step              *
!          (Used only for Time Error Estimator)                       *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Flux(NG, N_TOT)
! Output:   Flux_Coarse(NG, N_TOT)
! Local Variables:
      integer n, k, n1
  
      do k = 1, N_TOT 
         do n = 1, NG
            Flux_Coarse(n, k) = Flux(n, k)
         end do
      end do

      do n1 = 1, NZ
         do k = 1, NH
! change p(k,n1) to p_tot(k, n1) for decay heat calculations
            Pow_Coarse(k,n1) = p_tot(k,n1)
         end do
      end do

      return
      end 

      subroutine TSC_Recover_Data(time)
!=====================================================================*
!     Restoring the Model Parameters after the Coarse Time Step       *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Input: Old_Time, Old_Source(N_TOT), Old_Flux(NG, N_TOT), 
!       Old_TRL(NG, N_TOT, NDIR), Old_Prec(MD, N_TOT), Old_Zrods(IRODS)
! Output:   Time - Time
!           Flux(NG,N_TOT) - Neutro Flux
!           TRL_XYZ(NG, N_TOT, NDIR) -  Transverse Leakages
!           Prec(MD, N_TOT) - Precursor Concentration
!           Source(N_TOT) - Source Term 
!           Zrods(IRODS) - Control Rod Position
!           Flux_Mom(MOMENTS, NG, N_TOT, NDIR) - Flux Moments
      real Time
!Local Variables:
      integer k, n, m, nd, n1

      do k = 1, N_TOT
         Source(k) = Old_Source(k)
         do n = 1, NG
            Flux(n, k) = Old_Flux(n, k)
         end do
         do m = 1, MD
            Prec(m, k) =  Old_Prec(m, k)
         end do
      end do

      Time = Old_Time

      do nd = 1, NDD
         do k = 1, N_TOT
            do n = 1, NG
               TRL_XYZ(n, k, nd) = Old_TRL(n, k, nd)
            end do
         end do
      end do

      call CRD_Recover_data

       React = Old_React
       Pow_Point = Old_Pow_Point 
       do m = 1, MD
          Prec_Point(m) = Old_Prec_Point(m) 
       end do


       DO n1 = 1, NZ
          DO k = 1, NH
            p(k,n1) = Old_p(k,n1)
          END DO
      END DO  

       DO n1 = 1, NZ
          DO k = 1, NH
          DO  m = 1, MH
             p_dh(m,k,n1)= OLD_p_dh(m,k,n1)
          END DO
          END DO
      END DO  


!NOT USED      do nd = 1, NDIR
!NOT USED        do k = 1, N_TOT
!NOT USED          do m = 1, MOMENTS
!NOT USED            do n = 1, NG
!NOT USED                Flux_Mom(m, n,  k, nd)  = Old_Flux_Mom(m, n,  k, nd)
!NOT USED            end do
!NOT USED          end do
!NOT USED        end do
!NOT USED      end do

      return
      end


      subroutine TSC_Estimate_Temporal_Error
!=====================================================================*
!    Estimation of the Time Discretization Error Comparing the        *
! Solutions of the Coarse (Flux_Coarse) & Fine (Flux) Temporal Meshes *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Flux(NG, N_TOT), Flux_Coarse(NG, N_TOT)
! Output:   

! Local Variables
! i_order - an order of the time discretization method (equal to 1 for
! fully-implicit scheme)
      integer i_order
      parameter(i_order = 1)

      integer  k, n1
      real Flux_Norm, Eps

      TE_Flux = 0.
      Flux_Norm = 0.

!      do k = 1, N_TOT
!         do n = 1, NG
!            Eps = Flux(n, k) - Flux_Coarse(n, k)
!            Flux_Norm = Flux_Norm + Flux(n, k)*Flux(n, k)
!            TE_Flux = TE_Flux + Eps*Eps
!
!         end do
!      end do
!            TE_Flux = sqrt(TE_Flux)/sqrt(Flux_Norm)

      do n1 = 1, NZ      
         do k = 1, NH
! change p(k, n1) to p_tot for decay heat calculations
            Eps  = p_tot(k,n1) - Pow_Coarse(k,n1) 
            Flux_Norm = Flux_Norm + p_tot(k,n1)*p_tot(k,n1)
            TE_Flux = TE_Flux + Eps*Eps
         end do
      end do


!     write(*,*) 'TE_Flux = ', TE_Flux
!     write(*,*) 'Flux_Norm =', Flux_Norm
!     pause

            TE_Flux = sqrt(TE_Flux)/sqrt(Flux_Norm)

            Time_Error = TE_Flux/(2.**i_order - 1.)

      return
      end


      subroutine TSC_Estimate_Time_Step_Size&
                    (time, dt_step, step_accept,No_Restart)
!=====================================================================*
!         Estimation of the  Time Step Size Using                     * 
! the Time Error Estimation (Time_error) & Given Tolerance (ST_EPS)   *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: ST_EPS - A given Tolerance for the Error in Neutron Flux
!        Time_error - Computed Time Discretization Error
!        Time - Current Time
!        Dt_Step - Current Value of Time Step Size
!        facmax - Current value of the maximum increase of the time step
!                 Size (1 if the previous time step was rejected)
!        No_Restart - No reastrt for the calculations with TRAC
!                        No rejected time steps in this case
      real dt_step, Time
      Logical No_Restart
! Output: step_accept - if Time Step accepted = 1 , else Input Value (0)
!          dt_step - An estimate of the next time step size
!        facmax - Current value of the maximum increase of the time step
!                 Size (1 if the current time step was rejected)
      integer step_accept
! Local Variables:
! order - an order of the time discretization method (1 for the 
!           fully-imlicit scheme)
! safety_time - a safety parameter for the time step size
! facmax_in - maximum increase of the time step size
! facmin - maximum decrease of the time step size
      real  order, safety_time, facmax_in, facmin
      parameter(order = 1., safety_time = 0.8, facmax_in = 2., &
                facmin = 0.1)

      real dt_new, dt_factor

!      write(*,*) 'No_Restart=', No_Restart
!      pause

      dt_factor = (ST_EPS/Time_Error)**(1./(order+1.))
      dt_new = dt_step*amin1(facmax, amax1(facmin, safety_time*&
                     dt_factor))

!      write(*,*) 'dt_new =', dt_new, 'dt_step =', dt_step
!      write(*,*) 'dt_factor =', dt_factor
!      pause
      
      if(dt_new.gt.(dt_Step_Max)) dt_new = dt_Step_Max
  
      if(dt_factor.ge.1.) then
! time step was accepted 
         facmax = facmax_in
         step_accept = 1
      else
!* refuse from this time step and recalculation from time_0
         if(.NOT. No_Restart) then
            write(*,*) 'REFUSE from this time step'
!            pause 
            n_step_ref = n_step_ref + 1
            call TSC_Recover_Data(time)
            if(TH_Model.EQ."Internal") then
               call THM_Restore_Data
            end if
         end if
         facmax = 1.
      end if

      if((ST_EPS/Time_Error).LE.1) then
        write(*,*) 'ST_EPS/Time_Error,  dt_factor =',&
                    ST_EPS/TE_FLUX, dt_factor

!       write(*,*) 'dt_new =', dt_new, 'dt_step =', dt_step
!       pause

      end if

      dt_step = dt_new

      return
      end




