      subroutine KIN_Initialize(time)
!**********************************************************************
!       Computing of the initial data for the kinetics calculations
!**********************************************************************
      implicit none
      include 'sketch.fh'

! INPUT: N_TOT, NG, Flux(NG,N_TOT), XS_SF(NG_N_TOT), alfa(MD), beta(MD),
!        Pow_Norm, p_average (for 1st detector), p_col(NPOLY, NZ/NZC), 
!       NH_DAT(NDAT-1) - position of the detectors (Number of the Assembly)
!       NZ_DAT(NDAT-1) - position of the detectors (Axial Layer)
! Output: Flux(NG, N_TOT) -Neutron Flux normalized to the Reactor Power
!         trl_xyz(NG, N_TOT, NDD) - Transverse Leakage normalized to 
!                                          the Reactor Power
!         Source(N_TOT) - Fission Source Term
!         Prec(MD,N_TOT) - Steady-State Concentration of the Delayed Neutron
!                            Precursors
!        time - current time
!        dt_save - time step size

      real time
! Local Variables
      integer  k, m
!      integer  n_in_iter, n_out_iter, n_nonl_iter 
      Logical Adjoint
      real k_ef_norm

      time = 0.
      dt_save = dt_input(1)

      adjoint = .False.
      call EIS_Compute_Source(ADJOINT)
      k_ef_norm = 1./k_ef
      call MSC_SSCALE(N_TOT, k_ef_norm, Source)

!      write(*,*) 'SOURCE =, k = 205'
!      k = 205
!      DO n1=1, NZ/2
!        kt = k + (n1-1)*NH
!        write(*,*) 'n1, source =', n1, source(kt),&
!        flux(1,kt), flux(2,kt), xs_sf(1,kt), xs_sf(2,kt)
!      END DO 
!      Pause
!      write(*,*) 'SOURCE =, k = 206'
!      k = 206
!      DO n1=1, NZ/2
!        kt = k + (n1-1)*NH
!        write(*,*) 'n1, source =', n1, source(kt),      &
!        flux(1,kt), flux(2,kt), xs_sf(1,kt), xs_sf(2,kt)
!      END DO 
!      pause


      do k = 1, N_TOT
         do m = 1,MD
            Prec(m, k) = beta(m)*Source(k)/(alfa(m)*volume(k))
         end do
      end do


! Initialization of the Point Kinetics Equations
      Pow_Point = 1. !  P_Reactor
      do m = 1, MD
        Prec_Point(m) = beta(m)*Pow_Point/alfa(m)
      end do
      react = 0.
! Zero Omega
!OMEGA        do k = 1, N_TOT
!OMEGA          do n = 1, NG
!OMEGA              Omega_Flux(n,k)=0.
!OMEGA          end do
!OMEGA          do m = 1, MD
!OMEGA               Omega_Prec(m,k)=0.
!OMEGA           end do
!OMEGA        end do

! Setting data for CONTROL ROD SCRAM
        Time_scram = 1.E+30
        Flag_Set_Time_Scram = .False.

      return
      end


      subroutine KIN_Compute_Time_Step(time,dt_step, i_ssor, i_source, &
                        i_nonl_tot, i_time_fine)
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
      integer i_time_fine
      integer i_int, i_out, i_nonl
!  i_int - number of internal iterations
!         i_out - number of the source iterations
!         i_nonl_tot - number of the nonlinear iterations
! Output: 
      integer i_ssor, i_source, i_nonl_tot
! Local Variables: 
      real dt_rods, time_rods, dt_kin
      real Pow_Scale, shape_norm, index_not_used
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

      call CTRL_Solver(i_out, i_int, i_nonl, index_not_used, dt_kin)

      i_ssor  = i_int
      i_source = i_out
      i_nonl_tot = i_nonl


! Point Kinetics Equation Solution
      if(Kinetics_Method.NE."DRT") then
          shape_norm = 1. 
          call ADJ_Normalize_Flux(NG, N_TOT,  al, &
                 Flux_a(1,1), Flux(1,1), Volume(1))
      end if

      call PNT_Compute_Reactivity(dt_kin)


      if(Kinetics_Method .NE."DRT") then
           Pow_Scale = Pow_point
           call MSC_SSCALE(NG*N_TOT, Pow_Scale, Flux(1,1))
      end if 


      if(i_time_fine.eq.1) then
! Updating Delayed Neutron Precursors (only for the fine time steps)
         call KIN_Update_Precursors
         dt_dh=dt_kin ! передача шага по времени в
         call DEC_Heat_Update_Precursors_Euler()
      end if

      call POWer_Compute

      return
      end

      subroutine KIN_Update_Precursors
!=====================================================================*
! Calculation of delayed neutron precursors and Fission Source        *
!    at the End of Time Step                                          *
! (c) Slava 15.IV.1998                                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'         

! Input: Source_dt(N_TOT), Flux(NG, N_TOT), XS_SF(NG< N_TOT), 
!        Prec(MD, N_TOT), ecn(MD), eco(MD), volume(N_TOT), ec(MD)
! Output: Source(N_TOT) -Fission Source Term at the Current Time Step
!         Prec(MD, N_TOT) - Concentration of te Delayed Neutron at
!                             the Current Time Step
! Local Variables:
      real a_vol
      integer n, m, k

      do k = 1, N_TOT
         a_vol = 1./volume(k)
         do n = 1,NG
         end do

         do m = 1, MD
              Prec(m, k) = Prec(m, k)*ec(m) + (Source(k)*ecn(m) + &
                   Source_dt(k)*eco(m))*a_vol
         end do

      end do

      return
      end

      subroutine KIN_Init_Time_Step(dt_kin)
!=====================================================================*
! Computing the RHS for Neutron Kinetics Calculations                 *
!                   Vyachreslav Zimin (c) April 15 1998               *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'         
! Input: 
! dt_kin - Time Step Size
! alfa(MD) - Delayed neutron decay constant
! beta(MD) - delayed neutron yield fraction
! bet - total yield of delayed neutron per fission
! xp(NG) - prompt neutron fision spectrum
! xm(NG) - delayed neutron fission spectrum
! Source(N_TOT) - fission Source Term at the Previous Time Step
! Prec(MD, N_TOT) - Concentration of the delayed Neutron Precursors
! al(NG) = 1./v(NG)
! volume(N_TOT) - volume of the node
      real dt_kin

! Output: 
! ec(MD) - constant before delayed neutron precursor to compute
!           the delayed neutron precursors
! eco(MD) - constant before old fission source term to compute
!           the delayed neutron precursors
! ecn(MD) - constant before the new fission source term to compute
!           the delayed neutron precursors
! xpn(NG) - effective prompt fission spectrum
! RHS_K(NG, N_TOT) - Right Hand Side of the Neutron Flux Equations  
! Source_dt(NG, N_TOT) - the fission source term at the previous 
!                          time step
! Local Variables:
      real ecc(MD), ecsn, ecso, a1, sumc
      integer k, n, m

      ecsn=0.
      ecso=0.
      do m = 1,MD
         ec(m)=exp(-1.*alfa(m)*dt_kin)
         ecc(m) = ec(m)*alfa(m)
         ecn(m) = beta(m)*(1.-(1.-ec(m))/(alfa(m)*dt_kin))/alfa(m)
         eco(m) = beta(m)*((1.-ec(m))/(alfa(m)*dt_kin)-ec(m))/alfa(m)
         ecsn = ecsn + ecn(m)*alfa(m)
         ecso = ecso + eco(m)*alfa(m)
      end do

      do n = 1, NG
         xpn(n) = (1.-bet)*xp(n) + xm(n)*ecsn
      end do

! delayed neutron fission fission source calculation
      CALL EIS_Compute_Source(.False.)

      call MSC_SCOPY(N_TOT, Source, Source_dt)

      do k = 1, N_TOT
         sumc = 0.
         do m = 1, MD
            sumc = sumc + ecc(m)*Prec(m, k)
         end do
         sumc = sumc*volume(k) + ecso*Source_dt(k)
   
         a1 = volume(k) / dt_kin
         do n = 1,NG
!         the right side of the neutron diffusion equation
            MAT_RHS_K(n, k) = xs_al(n,k)*Flux(n, k)*a1 + &
                                          xm(n)*sumc
            ! MAT_RHS_K(n, k) = al(n)*Flux(n, k)*a1 + &
            !                               xm(n)*sumc
! saving Flux to compute the OMEGA
!OMEGA            Flux_Dt(n,k) = Flux(n,k)
        end do

      end do


      return
      end

      SUBROUTINE KIN_Normalize_Flux 
!**********************************************************************
!       Normalization of the Neutron Flux and Transverse Leakage
!**********************************************************************
      implicit none
      include 'sketch.fh'

      INTEGER n, k, nd

      do k = 1, N_TOT
         do n = 1, NG
            Flux(n,k)  = Flux(n, k)*Pow_Norm
          end do
      end do

      do nd = 1, NDIR
         do k = 1, N_TOT
           do n = 1, NG
            trl_xyz(n,k,nd)  = trl_xyz(n, k, nd)*Pow_Norm
           end do
         end do
      end do

!NOT USED      do nd = 1, NDIR
!NOT USED        do k = 1, N_TOT
!NOT USED        do m = 1, MOMENTS
!NOT USED           do n = 1, NG
!NOT USED            Flux_Mom(m, n,  k, nd)  = Flux_Mom(m, n,  k, nd)*Pow_Norm
!NOT USED           end do
!NOT USED          end do
!NOT USED        end do
!NOT USED      end do

      RETURN
      END
