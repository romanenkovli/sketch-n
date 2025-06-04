      SUBROUTINE THM_Compute_Average_Feedbacks
!=====================================================================*
! Computing Average & Maximum Values of the Feedbacks                 *
! 19.VII.2000 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE

      INCLUDE 'sketch.fh'

      REAL  Rel_Flow_Area(NN_FA_TYPE), &
       S_Fuel_per_Assembly(NN_FRD_FA, NN_FA_TYPE)
      
      COMMON /therm016/ Rel_Flow_Area, &
       S_Fuel_per_Assembly
  
      INTEGER  TH_TYPE_CORE(NP_REACTOR_CORE), N_FRD_FA(NN_FA_TYPE)  
      COMMON /therm017/   TH_TYPE_CORE, N_FRD_FA

! Local:
      INTEGER i, j, nc, np, i_type_fa
      REAL weighting_factor(N_POLY,2)
! weighting_factor(N_POLY,1) coolant volume
! weighting_factor(N_POLY,2) fuel volume 
! Setting up the weighting factor

      IF(TH_MODEL.NE."Internal") THEN

        CALL MSC_SSET(N_POLY, 1.0, weighting_factor(1,1) )
        CALL MSC_SSET(N_POLY, 1.0, weighting_factor(1,2) )

      ELSE
    
        CALL MSC_SSET(N_POLY, 0.0, weighting_factor(1,1))
        CALL MSC_SSET(N_POLY, 0.0, weighting_factor(1,2))

        DO nc = 1, NP_Reactor_Core
           np = Numb_Reactor_Core(nc)
           i_type_fa = TH_TYPE_CORE(nc) 
           weighting_factor(np,1) = Rel_Flow_Area(i_type_fa)
           weighting_factor(np,2) = 0.
           DO j = 1, N_FRD_FA(i_type_fa) 
           weighting_factor(np,2) = weighting_factor(np,2) + &
               S_Fuel_per_Assembly(j, i_type_fa)
           END DO  
         END DO
      END IF

      DO i = 1, N_FEEDBACK
      
      IF(i.eq.N_FEEDBACK) THEN
        j = 2 ! fuel 
      ELSE
        j = 1 ! coolant 
      END IF                 

      CALL OUTput_Compute_3D_Average(fdback(1,1,i), &
        NZ_Core_Beg, NZ_Core_End, Index_Core, fdback_col(0,0,i), &
        fdback_mm(-3,i), k_fdback_mm(1,-3,i), weighting_factor(1,j))

      CALL OUTput_Compute_2D_Average(fdback_col(0,0,i), NZR_Core_Beg, &
        NZR_Core_End, Index_Core, fdback_mm(-3,i), &
        k_fdback_mm(1,-3,i))        

      CALL OUTput_Compute_1D_Average(fdback_col(0,0,i), NZR_Core_Beg, &
        NZR_Core_End, Index_Core, fdback_mm(-3,i), &
        k_fdback_mm(1,-3,i), weighting_factor(1,j))        

      END DO

      RETURN
      END

      SUBROUTINE THM_Prepare_Output_Dist
!=====================================================================*
! Computing Average & Maximum Values of the Feedbacks                 *
! 19.VII.2000 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
! 
      INCLUDE 'sketch.fh'

      REAL  Temp_FR(NN_FRD_TOTAL+2, NN_FRD_FA, NZR, NP_Reactor_Core) 
      COMMON /therm003/ Temp_FR
      REAL Temp_Fuel_Av(NN_FRD_FA, NZR, NP_Reactor_Core),&
          Fuel_Enthalpy(NN_FRD_FA, NZR, NP_Reactor_Core) 
      COMMON /therm011/ Fuel_Enthalpy, Temp_Fuel_Av

      REAL  Rel_Flow_Area(NN_FA_TYPE), &
       S_Fuel_per_Assembly(NN_FRD_FA, NN_FA_TYPE)
      
      COMMON /therm016/ Rel_Flow_Area, &
       S_Fuel_per_Assembly
  
      INTEGER  TH_TYPE_CORE(NP_REACTOR_CORE), N_FRD_FA(NN_FA_TYPE)  
      COMMON /therm017/   TH_TYPE_CORE, N_FRD_FA

! Local:
      REAL dist_TH_Model(NP_Reactor_Core, NZR)
      INTEGER i, ns, nc, j, np, i_type_fa  
      REAL weighting_factor(N_POLY), weighting_rod_j

! Setting up the weighting factor
      CALL MSC_SSET(N_POLY, 0.0, weighting_factor)

      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         i_type_fa = TH_TYPE_CORE(nc) 
         weighting_factor(np) = 0.
         DO j = 1, N_FRD_FA(i_type_fa) 
           weighting_factor(np) = weighting_factor(np) + &
               S_Fuel_per_Assembly(j, i_type_fa)
         END DO  
      END DO

 
! Fuel Centerline Temperature
      i = 1
!        write(*,*) 'weighting_factor(np)=', weighting_factor(:)
!        pause
      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         i_type_fa = TH_TYPE_CORE(nc) 

            DO ns = 1, NZR
                  dist_TH_Model(nc, ns) = 0.
                DO j = 1, N_FRD_FA(i_type_fa) 
                weighting_rod_j = S_Fuel_per_Assembly(j, i_type_fa)/&
                  weighting_factor(np)
                dist_TH_Model(nc, ns) = dist_TH_Model(nc, ns)+&
                  Temp_FR(1, j, ns, nc)*weighting_rod_j
                END DO ! j = 1, N_FRD_FA(i_type_fa) 
            END DO ! ns = 1, NZR
      END DO

!      write(*,*) 'Temp_FR(1:, j, ns, nc)=', Temp_FR(:,1,2,1)
!      pause

         CALL OUTput_Convert_Dist_Core_Reactor(dist_TH_Model, &
        NZR_Core_Beg, NZR_Core_End,  dist_th_col(0,0,i), &
        dist_th_mm(-3,i), k_dist_th_mm(1,-3,i), &
        weighting_factor )


! Cladding Inner Surface Temperature
      i = 2
      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         i_type_fa = TH_TYPE_CORE(nc) 
            DO ns = 1, NZR
            dist_TH_Model(nc, ns) = 0.
            DO j = 1, N_FRD_FA(i_type_fa) 
               weighting_rod_j = S_Fuel_per_Assembly(j, i_type_fa)/&
               weighting_factor(np)
               dist_TH_Model(nc, ns) = dist_TH_Model(nc, ns)+&
              Temp_FR(NN_FRD_FUEL+2, j, ns, nc)*weighting_rod_j
            END DO ! j = 1, N_FRD_FA(i_type_fa) 
            END DO ! ns = 1, NZR
      END DO

      CALL OUTput_Convert_Dist_Core_Reactor(dist_TH_Model, &
        NZR_Core_Beg, NZR_Core_End,  dist_th_col(0,0,i), &
        dist_th_mm(-3,i), k_dist_th_mm(1, -3, i),&
        weighting_factor )

! Fuel Enthalpy
      i = 3
      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         i_type_fa = TH_TYPE_CORE(nc) 
            DO ns = 1, NZR
            dist_TH_Model(nc, ns) = 0.
            DO j = 1, N_FRD_FA(i_type_fa) 
             weighting_rod_j = S_Fuel_per_Assembly(j, i_type_fa)/&
               weighting_factor(np)
            dist_TH_Model(nc, ns) = dist_TH_Model(nc, ns)+&
              Fuel_Enthalpy(j, ns, nc)*weighting_rod_j
            END DO ! j = 1, N_FRD_FA(i_type_fa) 
            END DO ! ns = 1, NZR
      END DO

      CALL OUTput_Convert_Dist_Core_Reactor(dist_TH_Model, &
        NZR_Core_Beg, NZR_Core_End,  dist_th_col(0,0,i), &
        dist_th_mm(-3,i), k_dist_th_mm(1, -3, i),&
        weighting_factor )

!      write(*,*) 'dist_th_model(1,NZR_Core_Beg) =', &
!           dist_th_model(1,NZR_Core_Beg),&
!        Fuel_Enthalpy(1, NZR_Core_Beg, 1)

!      write(*,*) 'dist_th_col(0,0,i) =', dist_th_col(0,0,i)
!      pause  
!      write(*,*) 'dist_th_col(0,0,i) =', dist_th_col(0,0,i)
!      pause  

! 2D and 1D Distributions
      DO i = 1, N_OUT_TH_DIST

      CALL OUTput_Compute_2D_Average(dist_th_col(0,0,i), &
        NZR_Core_Beg, NZR_Core_End, Index_Core, dist_th_mm(-3,i), &
        k_dist_th_mm(1,-3, i) )        

      CALL OUTput_Compute_1D_Average(dist_th_col(0,0,i), &
        NZR_Core_Beg, NZR_Core_End, Index_Core, dist_th_mm(-3,i), &
        k_dist_th_mm(1,-3,i), weighting_factor )        

      END DO ! i = 1, N_OUT_TH_DIST

      CALL TH_COMPUTE_AVR_FROD(vol_ass, Numb_Reactor_Core, hz_core)


      RETURN
      END


      SUBROUTINE THM_Get_Core_Feedbacks
      implicit  none
      INCLUDE 'sketch.fh'

! common blocks for the couling
      real Temp_Cool(NZR, NP_Reactor_Core), &
             Temp_Doppl(NZR, NP_Reactor_Core),&
             Ro_Cool(NZR, NP_Reactor_Core)

      common /therm002/ Temp_Cool, Temp_Doppl,&
             Ro_Cool

      real  Temp_FR(NN_FRD_TOTAL+2, NN_FRD_FA, NZR, NP_Reactor_Core) 
      common /therm003/ Temp_FR
! End Commons blocks for data exchgange

! Local Variables
      integer n_c, np,  ns, j, k, n1, nz_in
!      real Convert_To_C, Convert_To_Gram
!      parameter (Convert_to_C = 273.15) ! Kelvin into Celcius
!      parameter (Convert_to_Gram = 1.E-03) ! Kg/m^3 into g/cm^3
            
      do n_c = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(n_c)
          n1 = 0
          do ns = 1, NZR
             do nz_in = 1, NPZ(ns)
                n1 = n1 + 1
                do j = 1, NCHM
                   k = poly_out(np,j)
                   if(k.ne.0)  then
                    fdback(k,n1,2) = temp_cool(ns, n_c)+Conv_Cool_Temp
                    fdback(k,n1,3) = ro_cool(ns, n_c)*Conv_Cool_Dens
                    fdback(k,n1,4) = temp_doppl(ns, n_c)+Conv_Fuel_Temp
!                    IF(n_c.eq.1) THEN
!                         WRITE(*,*) 'ns =', ns, 'Doppler Temp =', &
!                        temp_doppl(ns, n_c)
!                    PAUSE
!                    END IF
                   end if 
                end do ! NCHM
             end do ! NPZ(ns)
          end do ! NZR
       end do ! NP_Reactor_Core

!      call THM_Compute_Average_Feedbacks
      CALL THM_Prepare_Output_Dist

      return
      end 


      SUBROUTINE THM_Set_Core_Power
      implicit  none
      INCLUDE 'sketch.fh'

! common for the couling
      real  Pow_Core(NZR, NP_Reactor_Core )
      common /therm005/ Pow_Core

      integer ns, np, k
      real Convert_to_CI
      parameter (Convert_to_CI = 1.E+06) ! watts/cm^3 into watts/M^3


!      Power_Thermal = 0.
      do ns = 1, NZR
         do k = 1, NP_Reactor_Core
            np = Numb_Reactor_Core(k)
            Pow_Core(ns,k) = p_col(np,ns)*Convert_to_CI
!          Power_Thermal = Power_thermal + Pow_Core(n1,k)
         end do
      end do

!      WRITE(*,*) 'Total Thermal Power = ', Power_Thermal

      return
      end


      SUBROUTINE THM_Init(file_data, File_Dmp_In,&
        NZR_Core_Beg, NZR_Core_End, cool_heating_neutron, hz_neutron)
      implicit  none
      INCLUDE 'TH_Model.fh'
! Input:
      character*(*) file_data, file_dmp_in
      integer  NZR_Core_Beg, NZR_Core_End
      real hz_neutron(NZR)
      real cool_heating_neutron

! Local
      integer ns, k, i, NIN, i_type_fa, i_type_frd, j
      parameter (NIN = 0)

      real*8 T_PROP, P_PROP, RO_PROP, CP_PROP, Viscos_PROP,&
             Conduct_PROP, diffus_PROP
      real Ro_Cool_Inlet
! External function defined in FRD
      real FRD_fuel_enthalpy


!         call THM_Input

      cool_heating = cool_heating_neutron
      CALL THM_Input_Data(file_data, &
             hz_neutron, NZR_Core_Beg, NZR_Core_End) 

      CALL FRD_Input_Data(file_data)


      IF(File_DMP_In .EQ. "") THEN
! Inlet Temperature Distribution
           IF( TH_COOL_TYPE.EQ."WATER") THEN
            
            T_PROP = Temp_Cool_Inlet
            P_PROP = Pressure

            if(Steam_Table .EQ. 'JAERI') then
                call JAERI_STEAM(NIN, T_PROP, P_PROP, RO_PROP, &
                   CP_PROP, Viscos_PROP, Conduct_PROP)
            else if(Steam_Table.eq. 'TRAC') then
                 P_PROP = Pressure*1.E+06
            CALL PF1STEAM( P_PROP, T_PROP, RO_PROP, CP_PROP, &
              Viscos_PROP, Conduct_PROP)
            END IF ! (Steam_Table .EQ. 'JAERI')

           ELSE IF( TH_COOL_TYPE.EQ."LEAD") THEN

             CALL lead_properties(T_PROP, RO_PROP, Conduct_PROP,&
                  diffus_PROP, CP_PROP)
           END IF ! ( TH_COOL_TYPE.EQ."WATER")
  
            Ro_Cool_Inlet = RO_PROP

           do k = 1, NP_Reactor_Core
              i_type_fa = th_type_core(k)
              do ns = 1, NZR
                 DO j = 1, N_FRD_FA(i_type_fa) 
                    i_type_frd = frd_fa_type(j, i_type_fa)  
                    do i = 1, NN_FRD_TOTAL + 2
                       Temp_FR(i, j, ns, k) = Temp_Cool_Inlet
                     end do
                    Temp_Fuel_Av(j,ns,k) =  Temp_Cool_Inlet
                    Fuel_Enthalpy(j,ns,k) = &
                   FRD_fuel_enthalpy(Temp_Cool_Inlet,i_type_frd)
                 END DO
                 Temp_Doppl(ns, k) = Temp_Cool_Inlet
                 Temp_Cool(ns,k) = Temp_Cool_Inlet
                 Ro_Cool(ns,k) = RO_Cool_Inlet    
              end do
           end do
      END IF ! IF(File_DMP_In .EQ. "")

      CALL THM_Get_Core_Feedbacks

      RETURN
      END


      SUBROUTINE THM_Compute_Time_Step( a_dt)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

!Input:  
      real a_dt ! 1./Time_Step Size
! DEBUG
!      real vol_ass(N_POLY, NZR)

! Parameter for the steam table calculations ! NIN = 0 Subcooled Water
      integer  NIN
      parameter (NIN = 0)

      integer ns, k, i, i_type_frd, i_type_fa, j



      real Pow_Dens_Fuel, Pow_Dens_Cool, Pow_dens_FR
      real tc_int(NZR), tc_int_old
      real Time_Der, Part_Steady
      real*8 T_PROP, P_PROP, RO_PROP, CP_PROP, Viscos_PROP,&
             Conduct_PROP, diffus_PROP
      real Prandtl, Reynolds, Heat_Trans_Coeff, Res_Film_Clad
      real Temp_Old, Temp_FR_Old(5)
      real HTC_Constant
      parameter (HTC_Constant = 30000.) ! (W/m^2/K)
      Logical Constant_Heat_Trans_Coeff
      parameter (Constant_Heat_Trans_Coeff=.False.)
      REAL Part_Power_frd, Temp_Doppl_frd
      REAL Pekle, Nusselt, veloc_cool, Nusselt_pl, Nusselt_kon,&
          kappa
      REAL S_Fuel_per_Assembly_Total



!      WRITE(*,*) 'Internal TH Model'
!      WRITE(*,*) 'NZ_Reactor_Core(ns)', NZ_Reactor_Core
!      pause


      if(Constant_Heat_Trans_Coeff) then
        Heat_Trans_Coeff = HTC_Constant
      end if   

      do i = 1, 7
         e_temp_max(i) = 0.
      end do


      do k = 1, NP_Reactor_Core

         tc_int(1) = Temp_Cool_inlet
         i_type_fa = th_type_core(k)

         do ns = 1, NZR

            tc_int_old = tc_int(ns)
! compute the new value of the interface temperature
            if(ns.ne.1) then
              tc_int(ns) = 2.*temp_cool(ns-1,k) - tc_int(ns-1)
            end if        
! extrapolate the coolant temperature for the properties
            T_PROP = temp_cool(ns,k) 
            P_PROP = Pressure

        IF( TH_COOL_TYPE.EQ."WATER") THEN

        if(Steam_Table .EQ. 'JAERI') then


          call JAERI_STEAM(NIN, T_PROP, P_PROP, RO_PROP, CP_PROP,&
              Viscos_PROP, Conduct_PROP)

         else if(Steam_Table.eq. 'TRAC') then
            
             P_PROP = Pressure*1.E+06

            CALL PF1STEAM( P_PROP, T_PROP, RO_PROP, CP_PROP, &
              Viscos_PROP, Conduct_PROP)

         else
          WRITE(*,*) 'STEAM_TABLE is not specified'
          WRITE(*,*) 'Check the file "TH_Model.fh" '
          stop
         end if

         ELSE IF( TH_COOL_TYPE.EQ."LEAD") THEN

             CALL lead_properties(T_PROP, RO_PROP, Conduct_PROP,&
                  diffus_PROP, CP_PROP)

         END IF ! TH_COOL_TYPE

!           if(iteration.eq.1) then !DEBUG
!                Pow_Dens_Fuel = Pow_Core(n1,k)*(1. - Cool_Heating)&
!                                                   /Rel_Flow_Area
!           else  
            if(NZ_Reactor_Core(ns).EQ.1) then


              if(.NOT. Constant_Heat_Trans_Coeff) then
                  IF( TH_COOL_TYPE.EQ."WATER" ) THEN
                  Prandtl = CP_PROP*Viscos_PROP/Conduct_PROP
                  Reynolds = Diam_Eq(i_type_fa)*Flow_Rate_Core(k)/&
                 Viscos_PROP
                  Heat_Trans_Coeff = Const_DB*Conduct_PROP*&
                  (Reynolds**Pow_Reynolds)*(Prandtl**Pow_Prandtl)&
                                                  /Diam_Eq(i_type_fa)
                  ELSE IF( TH_COOL_TYPE.EQ."LEAD") THEN
                    veloc_cool = Flow_Rate_Core(k)/RO_PROP
                    Pekle = veloc_cool*Diam_Eq(i_type_fa)/diffus_PROP
                    Nusselt_kon = 7.55*ratio_s_to_d(i_type_fa)-&
                               14.*(ratio_s_to_d(i_type_fa)**(-5))+&
                 0.009*(Pekle**(0.64+0.246*ratio_s_to_d(i_type_fa))) 
                  kappa = (5. + 0.025*Pekle**0.8)/&
                              (7. + 0.025*Pekle**0.8)  
                    Nusselt_pl = Nusselt_kon*kappa/(kappa-1.)
! OLD                    Nusselt = 0.6 + 0.013*Pekle**0.7
!                    write(*,*) 'Nusselt OLD =', Nusselt
                    Nusselt = 1/( 1./Nusselt_pl+1./Nusselt_kon)
!                    write(*,*) 'Nusselt NEW =', Nusselt
!                    pause
!                    write(*,*) 'veloc_cool =', veloc_cool
!                    write(*,*) 'Pekle =', Pekle
!                    write(*,*) 'Nusselt =', Nusselt
!                    PAUSE
                    Heat_Trans_Coeff = Nusselt*Conduct_PROP/&
                         Diam_Eq(i_type_fa)
                  END IF 
              end if

              Res_Film_Clad = 1./Heat_Trans_Coeff 
              Pow_Dens_Fuel = 0.
              DO j = 1, N_FRD_FA(i_type_fa)
              Pow_Dens_Fuel = Pow_Dens_Fuel + gamma_HF(j,i_type_fa)*&
                (Temp_FR(NN_FRD_TOTAL+2, j, ns, k)-Temp_Cool(ns,k))/&
                 Res_Film_Clad
              END DO
! DEBUG
! All Power Generated in the Fuel Goes into the Coolant
!                IF (k.eq.1.and.ns.eq.17 ) THEN
!               write(*,*) 'Pow_Dens_Fuel (Fuel) =', Pow_Dens_Fuel 
!              END IF

!               Pow_Dens_Fuel = Pow_Core(ns,k)*(1.-Cool_Heating)/&
!                  Rel_Flow_Area(i_type_fa)

!                IF (k.eq.1.and.ns.eq.17 ) THEN
!               write(*,*) 'Pow_Dens_Fuel (Direct) =', Pow_Dens_Fuel 
!              END IF

!                IF (k.eq.1.and.ns.eq.17 ) THEN
!                  write(*,*) 'k=1, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel*&
!                           Rel_Flow_Area(i_type_fa)
!                  write(*,*) 'Total Generated Power=', Pow_Core(ns,k)*&
!                             (1.-Cool_Heating)*vol_ass(k, ns)*1.E-06
!                 write(*,*) 'Power From the Fuel=', Pow_Dens_Fuel*&
!                    Rel_Flow_Area(i_type_fa)*vol_ass(k, ns)*1.E-06 
!               END IF
!                IF (k.eq.5.and.ns.eq.17 ) THEN
!                  write(*,*) 'k=1, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel*&
!                           Rel_Flow_Area(i_type_fa)
!                  write(*,*) 'Total Generated Power=', Pow_Core(ns,k)*&
!                             (1.-Cool_Heating)*vol_ass(k, ns)*1.E-06
!                 write(*,*) 'Power From the Fuel=', Pow_Dens_Fuel*&
!                    Rel_Flow_Area(i_type_fa)*vol_ass(k, ns)*1.E-06 
!               END IF
!                IF (k.eq.7.and.ns.eq.17 ) THEN
!                  write(*,*) 'k=1, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel*&
!                           Rel_Flow_Area(i_type_fa)
!                  write(*,*) 'Total Generated Power=', Pow_Core(ns,k)*&
!                             (1.-Cool_Heating)*vol_ass(k, ns)*1.E-06
!                 write(*,*) 'Power From the Fuel=', Pow_Dens_Fuel*&
!                    Rel_Flow_Area(i_type_fa)*vol_ass(k, ns)*1.E-06 
!               END IF
!                IF (k.eq.7.and.ns.eq.17 ) THEN
!                  write(*,*) 'k=7, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel
!                  write(*,*) 'Pow_Core(ns,k)=', Pow_Core(ns,k)/&
!                   Rel_Flow_Area(i_type_fa)
!               END IF
           else
! Zero Heat Flux from the Fuel
              Pow_Dens_Fuel = 0. 
            end if
! Direct Heating of the Coolant

            Pow_Dens_Cool = Pow_Core(ns,k)*Cool_Heating/&
               Rel_Flow_Area(i_type_fa)


! computing node-average coolant temperature
            temp_old = temp_cool(ns,k)
            Time_Der = a_dt*hz(ns)*RO_PROP*CP_PROP
            Part_Steady = 2.*Flow_Rate_Core(k)*CP_PROP
            Temp_Cool(ns,k) = (Time_Der*temp_Cool(ns,k) + &
            Part_Steady*tc_int(ns) + &
             (Pow_Dens_Fuel + Pow_Dens_Cool)*hz(ns))/&
               (Time_Der + Part_Steady)

            e_temp_max(1) = amax1(e_temp_max(1), abs(temp_cool(ns,k) -&
                                 temp_old)) 
! coolant density
            
            T_PROP = temp_cool(ns,k)
            P_PROP = Pressure

        IF( TH_COOL_TYPE.EQ."WATER") THEN

        if(Steam_Table .EQ. 'JAERI') then
          call JAERI_STEAM(NIN, T_PROP, P_PROP, RO_PROP, CP_PROP,&
              Viscos_PROP, Conduct_PROP)
         else if(Steam_Table.eq. 'TRAC') then

            P_PROP = Pressure*1.E+06

            CALL PF1STEAM( P_PROP, T_PROP, RO_PROP, CP_PROP, &
              Viscos_PROP, Conduct_PROP)


         else
          WRITE(*,*) 'STEAM_TABLE is not specified'
          WRITE(*,*) 'Check file THERMAL.INI'
          stop
         end if

         ELSE IF( TH_COOL_TYPE.EQ."LEAD") THEN

             CALL lead_properties(T_PROP, RO_PROP, Conduct_PROP,&
                  diffus_PROP, CP_PROP)

         END IF ! TH_COOL_TYPE

          temp_old = ro_cool(ns,k)
          ro_cool(ns,k) = RO_PROP
          e_temp_max(2) = amax1(e_temp_max(2), abs(ro_cool(ns,k)/&
                                 temp_old-1.)) 

          temp_fr_old(1) = temp_FR(1, 1, ns,k )
          temp_fr_old(2) = temp_FR(NN_FRD_FUEL+1, 1, ns,k)
          temp_fr_old(3) = temp_FR(NN_FRD_FUEL+2, 1, ns,k)


          temp_fr_old(4) = temp_doppl(ns,k)
          Temp_fr_old(5) = Temp_FR(NN_FRD_TOTAL+2, 1, ns, k)

! computing the new heat transfer coefficient & Heat Conduction Equations
! for the Reactor Core ONLY
          if(NZ_Reactor_Core(ns).eq.1) then
!
          if(.NOT. Constant_Heat_Trans_Coeff) then
          IF( TH_COOL_TYPE.EQ."WATER" ) THEN
           Prandtl = CP_PROP*Viscos_PROP/Conduct_PROP
           Reynolds = Diam_Eq(i_type_fa)*Flow_Rate_Core(k)/Viscos_PROP
           Heat_Trans_Coeff = Const_DB*Conduct_PROP*&
                    (Reynolds**Pow_Reynolds)*&
              (Prandtl**Pow_Prandtl)/Diam_Eq(i_type_fa)
          ELSE IF( TH_COOL_TYPE.EQ."LEAD") THEN
                    veloc_cool = Flow_Rate_Core(k)/RO_PROP
                    Pekle = veloc_cool*Diam_Eq(i_type_fa)/diffus_PROP
                    Nusselt_kon = 7.55*ratio_s_to_d(i_type_fa)-&
                               14.*(ratio_s_to_d(i_type_fa)**(-5))+&
                 0.009*(Pekle**(0.64+0.246*ratio_s_to_d(i_type_fa))) 
                  kappa = (5. + 0.025*Pekle**0.8)/&
                              (7. + 0.025*Pekle**0.8)  
                    Nusselt_pl = Nusselt_kon*kappa/(kappa-1.)
! OLD                    Nusselt = 0.6 + 0.013*Pekle**0.7
!                    write(*,*) 'Nusselt OLD =', Nusselt
                    Nusselt = 1/( 1./Nusselt_pl+1./Nusselt_kon)
!                    write(*,*) 'Nusselt NEW =', Nusselt
!                    pause
!                    write(*,*) 'veloc_cool =', veloc_cool
!                    write(*,*) 'Pekle =', Pekle
!                    write(*,*) 'Nusselt =', Nusselt
!                    PAUSE
                    Heat_Trans_Coeff = Nusselt*Conduct_PROP/&
                         Diam_Eq(i_type_fa)
          END IF 
          end if

! Doppler Fuel Temperature is computed by volume-weighting
         Temp_Doppl(ns,k) = 0.

         S_Fuel_per_Assembly_Total = 0.
         DO j = 1, N_FRD_FA(i_type_fa) 
           S_Fuel_per_Assembly_Total = S_Fuel_per_Assembly_Total + &
               S_Fuel_per_Assembly(j, i_type_fa)
         END DO 

         DO j = 1, N_FRD_FA(i_type_fa) 
            i_type_frd = frd_fa_type(j, i_type_fa) 
            Pow_Dens_FR = Pow_Core(ns,k)*(1. - Cool_Heating)
            Part_Power_frd = &
           Pow_FRD_Assm(j, i_type_fa)/&
           S_Fuel_per_Assembly_Total

       
         call FRD_Compute_Time_Step( a_dt, Pow_Dens_FR, &
                Temp_Cool(ns,k), Heat_Trans_Coeff,&
                Temp_FR(1,j,ns,k), Temp_Doppl_frd,&
                Temp_Fuel_Av(j, ns,k), Fuel_Enthalpy(j,ns,k),&
                Part_Power_frd, i_type_frd, k, ns )

!                IF (k.eq.1.and.ns.eq.17 ) THEN
!                 Pow_Dens_Fuel = gamma_HF(j,i_type_fa)*
!     *                   (Temp_FR(NN_FRD_TOTAL+2, j, ns, k)- &
!                    Temp_Cool(ns,k))*Heat_Trans_Coeff
!                  write(*,*) 'k=1, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel*
!!                  write(*,*) 'Pow_Core(ns,k)=', Pow_Core(ns,k)
!                pause
!               END IF
!                IF (k.eq.7.and.ns.eq.17 ) THEN
!           Pow_Dens_Fuel = gamma_HF(j,i_type_fa)*
!     *                   (Temp_FR(NN_FRD_TOTAL+2,j,ns,k)- &
!                    Temp_Cool(ns,k))*Heat_Trans_Coeff
!                  write(*,*) 'k=7, ns=17', k, ns
!                  write(*,*) 'Pow_Dens_Fuel=', Pow_Dens_Fuel*
!!                  write(*,*) 'Pow_Core(ns,k)=', Pow_Core(ns,k)*&
!                       (1.-Cool_Heating)
!                  pause
!               END IF

         Temp_Doppl(ns,k) = Temp_Doppl(ns,k) + Temp_Doppl_frd*&
              TH_Dopl_Mix(j, i_type_fa)   

!         IF ( ns.EQ. 5.AND. k.eq.8) THEN
!           WRITE(*,*) 'i_type_frd, Temp_Doppl_frd', i_type_frd, &
!           Temp_Doppl_frd
!         END IF 

         END DO 


         end if

!          if(ns.eq.1) then
!           WRITE(*,*) 'temp_doppl(ns,k) =', temp_doppl(ns,k)
!           pause
!          end if


          e_temp_max(4) = amax1(e_temp_max(4), abs(temp_FR(1,1,ns,k)-&
                                 temp_fr_old(1))) 
          e_temp_max(5) = amax1(e_temp_max(5), &
              abs(temp_FR(NN_FRD_FUEL+1,1,ns,k)-temp_fr_old(2))) 
          e_temp_max(6) = amax1(e_temp_max(6), &
              abs(temp_FR(NN_FRD_FUEL+2,1,ns,k)-temp_fr_old(3))) 
          e_temp_max(7) = amax1(e_temp_max(7), abs(temp_Doppl(ns,k)-&
                                 temp_fr_old(4))) 
          e_temp_max(3) = amax1(e_temp_max(3), &
               abs(temp_FR(NN_FRD_TOTAL+2,1,ns,k)-temp_fr_old(5)))

! Maximum Fuel Temperature


         end do

      end do

      e_temp_max(2) = e_temp_max(2)*100.

      return
      end

      SUBROUTINE THM_Save_Data
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

      integer n, k, i, j, i_type_fa

       do n = 1, NZR
          do k = 1, NP_Reactor_Core
            i_type_fa = TH_TYPE_CORE(k)
            Old_Temp_Cool(n, k) = Temp_Cool(n, k)
            Old_Temp_Doppl(n,k) = Temp_Doppl(n,k)
            Old_Ro_Cool(n,k) = Ro_Cool(n,k)
                DO j = 1, N_FRD_FA(i_type_fa)
                    do i = 1, NN_FRD_TOTAL + 2
                     Old_Temp_FR(i, j, n, k) = Temp_FR(i, j, n, k)
                    end do
                END DO
          end do
       end do

      return
      end


      SUBROUTINE THM_Restore_Data
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

      integer n, k, i, i_type_fa, j

       do n = 1, NZR
          do k = 1, NP_Reactor_Core
            i_type_fa = TH_TYPE_CORE(k)
            Temp_Cool(n, k) = Old_Temp_Cool(n, k)
            Temp_Doppl(n,k) = Old_Temp_Doppl(n,k)
            Ro_Cool(n,k) = Old_Ro_Cool(n,k)
                DO j = 1, N_FRD_FA(i_type_fa)
                    do i = 1, NN_FRD_TOTAL + 2
                     Temp_FR(i, j, n, k) = Old_Temp_FR(i, j, n, k)
                    end do
                END DO
          end do
       end do

      return
      end

      SUBROUTINE THM_Get_Hot_Pellet_Data(k, ns, Fuel_Enth_Hot_Pellet,&
        Cool_Out_Temp_Hot_Pellet, Fuel_Temp_CL_Hot_Pellet,&
        Clad_Out_Temp_Hot_Pellet, Fuel_Temp_Av_Hot_Pellet)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'
! Input:
      integer k, ns
! Output:
      real Fuel_Enth_Hot_Pellet, Cool_Out_Temp_Hot_Pellet, &
        Fuel_Temp_CL_Hot_Pellet, Clad_Out_Temp_Hot_Pellet,&
        Fuel_Temp_Av_Hot_Pellet

! Local Variables:
      real Convert_K_to_C
      parameter(Convert_K_to_C = - 273.15)

      Fuel_Enth_Hot_Pellet = Fuel_Enthalpy(ns, k, 1) 
      Cool_Out_Temp_Hot_Pellet = Temp_Cool(NZR, k) + Convert_K_to_C
      Fuel_Temp_CL_Hot_Pellet = Temp_FR(1, 1, ns, k) + Convert_K_to_C
      Clad_Out_Temp_Hot_Pellet = Temp_FR(NN_FRD_TOTAL+2, 1, ns, k) +&
        Convert_K_to_C
      Fuel_Temp_Av_Hot_Pellet = Temp_Fuel_Av(1, ns, k) + Convert_K_to_C

      return
      end

      SUBROUTINE THM_read_data_restart_file(unit)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

! Input:
      INTEGER unit

      integer n, k, i, j


      read(unit) ((Temp_Cool(n, k),n=1,NZR), k = 1,&
           NP_Reactor_Core)
      read(unit) ((Temp_Doppl(n,k), n=1,NZR), k = 1,&
           NP_Reactor_Core)
      read(unit) ((Ro_Cool(n,k), n = 1, NZR), &
           k = 1, NP_Reactor_Core)
      read(unit) (((Temp_Fuel_Av(j,n,k), &
           j=1, NN_FRD_FA ),&
           n=1,NZR), k = 1,&
           NP_Reactor_Core)
      read(unit) (((Fuel_Enthalpy(j,n,k),&
           j=1, NN_FRD_FA ),&
           n=1,NZR), k = 1, NP_Reactor_Core)
      read(unit) ((((Temp_FR(i, j, n, k), i=1, NN_FRD_TOTAL + 2), &
           j=1, NN_FRD_FA),&
                          n = 1, NZR), k = 1, NP_Reactor_Core)

      return
      end


      SUBROUTINE THM_WRITE_data_restart_file(unit)
      implicit  none
      INCLUDE 'TH_Model.fh'

! Input:
      INTEGER unit

      integer k, n, i, j

      WRITE(unit) ((Temp_Cool(n, k),n=1,NZR), &
                k = 1, NP_Reactor_Core)
      WRITE(unit) ((Temp_Doppl(n,k), n=1,NZR), &
                   k = 1, NP_Reactor_Core)
      WRITE(unit) ((Ro_Cool(n,k), n = 1, NZR), &
            k = 1, NP_Reactor_Core)
      WRITE(io_unit) (((Temp_Fuel_Av(j, n,k), &
           j=1, NN_FRD_FA ),&
           n=1,NZR), k = 1,&
           NP_Reactor_Core)
      WRITE(io_unit) (((Fuel_Enthalpy(j, n,k), &
           j=1, NN_FRD_FA ),&
      n=1,NZR), k = 1,&
           NP_Reactor_Core)
      WRITE(unit) ((((Temp_FR(i, j, n, k), i=1, NN_FRD_TOTAL + 2), &
           j=1, NN_FRD_FA ),&
                          n = 1, NZR), k = 1, NP_Reactor_Core)

      return
      end


      SUBROUTINE THM_Input_Data(file_data, &
       hz_neutron, NZR_Core_Beg, NZR_Core_End) 
      implicit  none
      INCLUDE 'TH_Model.fh'
!=====================================================================*
! Input T/H Parameters from the file  "TH_Core.DAT"                   *
!=====================================================================*
! Input:
      character*(*) file_data
      real hz_neutron(NZR)
      integer NZR_Core_Beg, NZR_Core_End
! Local
      integer ns, ios, np, i, nt, j
      real Conv_cm_to_m
      parameter (Conv_CM_to_M = 1.E-02)
      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      logical error_find

! Marking Reactor core by 1
      DO ns = 1, NZR
         NZ_Reactor_Core(ns) = 0
      END DO
      DO ns = NZR_Core_Beg, NZR_Core_End
         NZ_Reactor_Core(ns) = 1
      END DO
! Convert neutronics axila mesh into SI units
      DO ns = 1, NZR
         hz(ns) = Conv_CM_to_M*hz_neutron(ns)
      END DO

!initialization of the identifiers
      WRITE(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'

      open(io_unit,file = file_data ,status='old', iostat=ios)
!TH_COOL_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_COOL_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier TH_COOL_TYPE in the FILE_INPUT' //&
         "set coolant to 'WATER' " )

         TH_COOL_TYPE="WATER"

      else
   
       read(io_unit,  fmt=*, iostat=ios) TH_COOL_TYPE
        call Iostat_Error_Check(ios,"Error Reading "// &
      "TH_COOL_TYPE "//&
      "from the FILE_INPUT file")

        IF(TH_COOL_TYPE.NE."WATER".AND.  TH_COOL_TYPE.NE."LEAD")&
       THEN
          ios = -1  
           call Iostat_Error_Check(ios,"Coolant type should "// &
         "be either 'WATER' or 'LEAD' ")
        END IF

      end if


!TH_HTCF_CNST
      IF(TH_COOL_TYPE.EQ."WATER") THEN
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_HTCF_CNST", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier TH_HTCF_CNST in the FILE_INPUT' //&
         "set the constants to 0.023, 0.8, 0.4" )

            Const_DB = 0.023
            Pow_Reynolds = 0.8
            Pow_Prandtl = 0.4
      else
   
       read(io_unit,  fmt=*, iostat=ios) Const_DB, Pow_Reynolds, &
         Pow_Prandtl
        call Iostat_Error_Check(ios,"Error Reading "// &
      "Heat Transfer Coefficient Constants under identifier" //&
      " TH_HTCF_CNST from the FILE_INPUT file")

      end if
      END  IF

!TH_TEMC_INLT
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_TEMC_INLT", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_TEMC_INLT in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) Temp_Cool_Inlet

       call Iostat_Error_Check(ios,"Error Reading "// &
      "Inlet Coolant Temperature from the FILE_INPUT file")

      end if

!TH_PRES_COOL
      IF(TH_COOL_TYPE.EQ."WATER") THEN
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_PRES_COOL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_PRES_COOL in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) Pressure

       call Iostat_Error_Check(ios,"Error Reading "// &
      "Coolant Pressure, Mpa under identifier TH_PRES_COOL"//&
      " from the FILE_INPUT file")

      end if
      END IF
!TH_N_FA_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_N_FA_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_N_FA_TYPE in the FILE_INPUT'//&
          'set TH_N_FA_TYPE = 1')
         
         N_FA_TYPE = 1

      else
   
       read(io_unit,  fmt=*, iostat=ios) N_FA_TYPE

       call Iostat_Error_Check(ios,"Error Reading "// &
      "N_FA_TYPE from the FILE_INPUT file")

      end if

!TH_TYPE_CORE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_TYPE_CORE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_Type_Core in the FILE_INPUT' //&
         ' set TH_TYPE_CORE =1' )
         DO np = 1, NP_Reactor_Core
            TH_TYPE_CORE(np) = 1
         END DO
      else
   
       read(io_unit,  fmt=*, iostat=ios) (TH_TYPE_CORE(i),i=1, &
          NP_Reactor_Core)

       call Iostat_Error_Check(ios,"Error Reading "// &
      " Fuel Assembly types under identifier TH_TYPE_CORE"//&
      " from the FILE_INPUT file")

      end if

!TH_MFRT_ASSM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_MFRT_ASSM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_MFRT_ASSM in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) (Flow_Rate(i),i=1, N_FA_TYPE)

        call Iostat_Error_Check(ios,"Error Reading "// &
      "Mass Flow Rate under identifier TH_MFRT_ASSM"//&
      " from the FILE_INPUT file")

      end if


!TH_MFRT_CORE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_MFRT_CORE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_MFRT_CORE in the FILE_INPUT')
         DO np = 1, NP_Reactor_Core
           nt = TH_TYPE_CORE(np)
           Flow_Rate_Core(np) = Flow_Rate(nt)
         END DO             
         Flag_Flow_Rate_Core = .False.

      else
   
       read(io_unit,  fmt=*, iostat=ios) (Flow_Rate_Core(i),&
            i=1, NP_Reactor_Core)

        call Iostat_Error_Check(ios,"Error Reading "// &
      "Mass Flow Rate from the FILE_INPUT file")

        Flag_flow_Rate_Core = .True.

      end if

!TH_HYDR_DIAM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_HYDR_DIAM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_HYDR_DIAM in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) (Diam_Eq(i),i=1, N_FA_TYPE)

        call Iostat_Error_Check(ios,"Error Reading "// &
      "Hydraulic Diameter under identifier TH_HYDR_DIAM"//&
      " from the FILE_INPUT file")

      end if

!TH_FLOW_AREA
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_FLOW_AREA", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_FLOW_AREA in the FILE_INPUT')

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (Rel_Flow_Area(i),i=1, N_FA_TYPE) 
        call Iostat_Error_Check(ios,"Error Reading "// &
      "Realtive Flow Area under identifier TH_FLOW_AREA " //&
      "from the FILE_INPUT file")

      end if

!TH_NFRD_ASSM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_NFRD_ASSM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_NFRD_ASSM in the FILE_INPUT'//&
         ' set N_FRD_FA = 1')
       DO i = 1, N_FA_TYPE
            N_FRD_FA(i) = 1
         END DO              

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (N_FRD_FA(i),i=1, N_FA_TYPE) 
        call Iostat_Error_Check(ios,"Error Reading "// &
      "number of the equivalent fuel rods in FA" //&
      "under idfentifer TH_NFRD_ASSM from the FILE_INPUT file")

      end if

!TH_FRD_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_FRD_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_FRD_TYPE in the FILE_INPUT'//&
         ' set FRD_TYPE_FA = 1')
         DO i = 1, N_FA_TYPE
            DO j = 1, N_FRD_FA(i)
               FRD_FA_TYPE(j, i) = 1
            END DO 
         END DO              

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (( FRD_FA_TYPE(j,i),j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

        call Iostat_Error_Check(ios,"Error Reading "// &
      " fuel rods types in fuel assembly under identifier" //&
      " TH_FRD_TYPE from the FILE_INPUT file")

      end if

!TH_SURF_VOLM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_SURF_VOLM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_SURF_VOLM in the FILE_INPUT')

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
         ((Gamma_HF(j,i) ,j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 
        call Iostat_Error_Check(ios,"Error Reading "// &
      "Ratio of the Surface Area of Cladding to the Coolant Volume "//&
      " under identifier  TH_SURF_VOLM"//&
      "from the FILE_INPUT file")

      end if

!TH_RAT_ASFL

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_RAT_FLAS", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_RAT_FLAS in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       (( S_Fuel_per_assembly(j,i),j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Ratio of the Fuel Volume to Assembly Volume "//&
       "under identifier TH_RAT_FLAS from the FILE_INPUT file")

      end if
!TH_POW_ASSM

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_POW_ASSM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_POW_ASSM in the FILE_INPUT'//&
         ' set TH_POW_ASSM = 1.')

         DO i=1, N_FA_TYPE
            DO j = 1, N_FRD_FA(i)
                Pow_FRD_Assm(j,i) = 1.
            END DO
         END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       ((Pow_FRD_Assm(j,i) ,j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Ratio of the Power Density in the Fuel Rod Type "//&
       "under identifier TH_POW_ASSM from the FILE_INPUT file")

      end if
! TH_DOPL_MIX

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_DOPL_MIX", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_DOPL_MIX in the FILE_INPUT'//&
         ' set TH_DOPL_MIX = 1./N_FRD_FA(i)')

         DO i=1, N_FA_TYPE
            DO j = 1, N_FRD_FA(i)
                TH_DOPL_MIX(j,i) = 1./( N_FRD_FA(i) )
            END DO
         END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       ((TH_DOPL_MIX(j,i) ,j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Coefficient to Compute Effective Doppler Temperature "//&
       "under identifier TH_DOPL_MIX from the FILE_INPUT file")

      end if

!TH_CONV_UNIT 
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_CONV_UNIT", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier TH_CONV_UNIT in the FILE_MAP,  '//&
         'constants of the unit conversion are not given'//&
         'set them to 1.E-03, -273.15, 0.' )
          Conv_Cool_Dens = 1.E-03
          Conv_Cool_Temp  = -273.15
          Conv_Fuel_Temp = 0.
      else
        read(io_unit,fmt=*,iostat=ios) &
           Conv_Cool_Dens, Conv_Cool_Temp, Conv_Fuel_Temp

        call Iostat_Error_Check&
      (ios,"Error in Reading Conv_Cool_Dens, Conv_Cool_Temp, &
                  Conv_Fuel_Temp from the FILE_Input")

      end if

! TH_MESH_DGR
      IF(TH_COOL_TYPE.EQ."LEAD") THEN

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_MESH_DGR", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Error',' Could not find '//&
         'identifier TH_MESH_DGR in the FILE_INPUT')

      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       (ratio_s_to_d(i), i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Ratio of the fuel rod mesh size to hydraulic diameter"//&
       "under identifier TH_MESH_DGR from the FILE_INPUT file")

      end if
      END IF


      close(io_unit)

      return
      end



      SUBROUTINE THM_Output_Parameters(unit)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

      integer unit

      WRITE(unit, *)

      WRITE(unit, '(A)')&
       "   Thermal-Hydraulics Parameters:"
      WRITE(unit, '(A, I8)') &
       "       Number of Heat Conduction Zones in Fuel Rod:", &
       NN_FRD_FUEL
      WRITE(unit, '(A, I8)') &
       "       Number of Heat Conduction Zones in Cladding:", &
       NN_FRD_CLAD
      WRITE(unit, '(A, I8)') &
       "       Number of Axial Layers                     :", &
       NZR
      WRITE(unit, '(A, I8)') &
       "       Number of Thermal-Hydraulics Channels      :", &
       NP_Reactor_Core

      IF(TH_COOL_TYPE.EQ."WATER") THEN
      WRITE(unit, '(A, 4x, A8)')&
       "       Steam Table                                :", &
       Steam_Table
      END IF  

      call OUTput_Write_Separator(unit)

      RETURN 
      END

      SUBROUTINE THM_Output_Data(unit)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

      integer unit, j, i
      CHARACTER*80 Header_Map

      WRITE(unit, '(A)') " Data of Internal Thermal-Hydraulics Model:"

! "     Inlet Coolant Temperature, [K]: "
      WRITE(unit, '(A)') &
       "           Coolant                                    : "
      WRITE(unit, '(8x,A)') TH_COOL_TYPE

! "     Inlet Coolant Temperature, [K]: "
      WRITE(unit, '(A)') &
       "     Inlet Coolant Temperature, [K]                   : "
      WRITE(unit, '(8x,F10.4)') Temp_Cool_Inlet

! "     Coolant Pressure, [Mpa]       : "
      WRITE(unit, '(A)') &
       "     Coolant Pressure, [Mpa]                          : "
      WRITE(unit, '(8x,F10.4)')  Pressure

! "     Number of FA Types
      WRITE(unit, '(A)') &
       "     Number of FA Types                               : "
      WRITE(unit, '(8x,I10)')  N_FA_TYPE

! "     Number of FA Types
      WRITE(unit, '(A)') 
      HEADER_MAP =&
       "     Core Loading by FA Types                         : "

      CALL TH_OUT_Map_INTEGER(unit, TH_TYPE_CORE, Header_Map)

! "     Coolant Flow Rate per Assembly , [kg/(m^2*s)]
!      IF ( Flag_flow_Rate_Core ) THEN
      WRITE(unit, '(A)') &
       "     Coolant Flow Rate In Core     , [kg/(m^2*s)]     : "

      Call OUT_Map_Flow_Rate_Core(unit, Flow_Rate_Core)  

!      END IF

      WRITE(unit, '(A)') &
       "     Coolant Flow Rate per Assembly, [kg/(m^2*s)]     : "
      WRITE(unit, '(8x,10E12.5)')  (Flow_Rate(i), i=1, N_FA_TYPE) 

!"     Equivalent Hydraulic Diameter, [m]    : "
      WRITE(unit, '(A)') &
       "     Equivalent Hydraulic Diameter, [m]                : "
      WRITE(unit, '(8x,10E12.5)')  (Diam_Eq(i), i=1, N_FA_TYPE)

!"     Heat Transfer Coefficient Constants : "
      WRITE(unit, '(A)') &
       "     Heat Transfer Coefficient Constants               : "
      WRITE(unit, '(8x,3E12.5)')  Const_DB, Pow_Reynolds, &
         Pow_Prandtl

!"     Ratio of Surface Area to Volume for Coolant, [1/m] :"
      WRITE(unit, '(A)') &
       "     Ratio of Surface Area to Volume for Coolant, [1/m]:"
      WRITE(unit, '(8x,10E12.5)') ((Gamma_HF(j,i), j=1, N_FRD_FA(i)),&
        i=1, N_FA_TYPE)  

!"     Relative Coolant Flow Area :"
      WRITE(unit, '(A)') &
       "     Relative Coolant Flow Area                        :"
      WRITE(unit, '(8x,10E12.5)')  (Rel_Flow_Area(i), i=1, N_FA_TYPE) 

! "     Ratio of the Fuel Volume to the Assembly Volume        : "
      WRITE(unit, '(A)') &
       "     Ratio of the Fuel Volume to the Assembly Volume        : "
      WRITE(unit, '(8x,10E12.5)') ((S_Fuel_per_Assembly(j, i), &
        j=1, N_FRD_FA(i)), i=1, N_FA_TYPE)  

! "Mixing Coefficients to Compute Effective Doppler Temperature"
      WRITE(unit, '(A)') &
       "     Mixing Coefficients to Compute Doppler Temperature     : "
      WRITE(unit, '(8x,10E12.5)') ((TH_Dopl_Mix(j,i), &
        j=1, N_FRD_FA(i)), i=1, N_FA_TYPE)  
! "Coefficient to compute the power density in the fuel rod type i"
      WRITE(unit, '(A)') &
       "      Coefficient for power density in fuel rod type i      : "
      WRITE(unit, '(8x,10E12.5)') ((Pow_FRD_Assm(j,i), &
        j=1, N_FRD_FA(i)), i=1, N_FA_TYPE)  

! "     Conversion Constants from T/H model to Feedbacks   :"
      WRITE(unit, '(A)') &
       "     Conversion Constants from T/H model to Feedbacks   :"
      WRITE(unit, '(8x,3E14.6)')  Conv_Cool_Dens, Conv_Cool_Temp, &
       Conv_Fuel_temp

      IF(TH_COOL_TYPE.EQ."LEAD") THEN
! "Ratio of the fuel rod mesh size to the hydraulic diameter"
      WRITE(unit, '(A)') &
       "     Ratio of the fuel rod mesh size to the "//&
       "hydraulic diameter: "
      WRITE(unit, '(8x,10E12.5)') (Ratio_s_to_d(i),i=1, N_FA_TYPE)  
      END IF


      call OUTput_Write_Separator(unit)


      RETURN
      END 

      subroutine THM_Feedbacks_Output(unit)
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      character*85 Header_Map
      character*4 val_fmt
      character*6 val_char(0:N_POLY)
      integer ind, ns, i
!      integer nlx_core, nly_core
! external function (in "OUTput.f")
      integer nlz_core
!      CHARACTER*30 feedb_name(N_FEEDBACK)

      real fdback_scaling_factor(N_FEEDBACK)
      data fdback_scaling_factor / 1.E-3, 1.E-3, 1., 1.E-3 /

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            FEEDBACK DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)


      DO i = 1, N_FEEDBACK
! 1st - Power, then FEEDBACKS
         WRITE(Header_Map, '(5x, A)') &
           NAME_TR_DIST(i+1)
!        write(*,*) 'i = ', i


      CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
           fdback_scaling_factor(i), N_POLY, NZR, &
           fdback_col(0,0,i), fdback_mm(-3,i), k_fdback_mm(1,-3,i) )

      END DO


!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            2D FEEDBACK DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)


      DO i = 1, N_FEEDBACK
      
      
! 1st - Power, then FEEDBACKS
      WRITE(Header_Map, '(5x, A)') &
           NAME_TR_DIST(1+i)

      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         fdback_col(ind, 0, i)*fdback_scaling_factor(i)
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL FEEDBACK DISTRIBUTIONS"
!      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(1x, 5A, /)') "  N  ", &
     " Br CONC",&
     " T  COOL",  &
     " RO COOL",&
     " T  DPLR"

       do ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End
          WRITE(unit,'(1x, I3,": ", 4F8.4)') &
         nlz_core( ns ) , &
         (fdback_col(0, ns, i)*fdback_scaling_factor(i), &
          i =1, N_FEEDBACK)
       end do

      CALL OUTput_Write_Separator(unit)

 
      RETURN
      END   

      subroutine THM_Internal_Distr_Output(unit)
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      character*85 Header_Map
      character*4 val_fmt
      character*6 val_char(0:N_POLY)
      integer ind, ns, i
!      integer nlx_core, nly_core
! external function (in "OUTput.f")
      integer nlz_core

      real dist_th_scaling_factor(N_OUT_TH_DIST)
      data dist_th_scaling_factor / 1.E-3, 1.E-3, 1./

!  Scaling Factor for Fuel Enthalpy = 1./dist_th_av!
!      write(*,*) 'dist_th_col(0,0,N_OUT_TH_DIST)=', &
!      dist_th_col(0,0,N_OUT_TH_DIST)
      
      IF(  dist_th_col(0,0,N_OUT_TH_DIST).LT.SMALL_VALUE) THEN
      dist_th_scaling_factor(N_OUT_TH_DIST) = 1.
        dist_th_mm(-3,N_OUT_TH_DIST) = 0.
      ELSE
      dist_th_scaling_factor(N_OUT_TH_DIST) = &
        1./dist_th_col(0,0,N_OUT_TH_DIST)
      END IF  
!      write(*,*) 'dist_th_scaling_factor(N_OUT_TH_DIST)=',&
!       dist_th_scaling_factor(N_OUT_TH_DIST)
!      pause  

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            T/H PARAMETER DISTRIBUTION"
      WRITE(unit,'(A, E12.6, A)') &
     "            Average Value of Fuel Enthalpy =", &
         dist_th_col(0,0, N_OUT_TH_DIST), " [J/Kg]"
      CALL OUTput_Write_Separator(unit)



      DO i = 1, N_OUT_TH_DIST
!         WRITE(*,*) 'I =', I
!         write(*,*) NAME_TR_DIST(1 + N_FEEDBACK + i )
!           IF ( i == 3) THEN
!             write(*,*)  dist_th_mm(-3,i)              
!             pause
!           END IF
 
! Flux(NG) + Power + FDBacks + T/H Distributions
         WRITE(Header_Map, '(5x, A)') &
           NAME_TR_DIST(1 + N_FEEDBACK + i )

      CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
          dist_th_scaling_factor(i), &
          N_POLY, NZR, dist_th_col(0,0,i), dist_th_mm(-3,i), &
          k_dist_th_mm(1,-3,i) )

      END DO


!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            2D T/H PARAMETER DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)


      DO i = 1, N_OUT_TH_DIST
      
      
      WRITE(Header_Map, '(5x, A)') &
           NAME_TR_DIST(1+N_FEEDBACK + i)

      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         dist_th_col(ind,0, i)*dist_th_scaling_factor(i)
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL T/H PARAMETER DISTRIBUTIONS "
!      CALL OUTput_Write_Separator(unit)

      WRITE(unit, '(1x, 4A, /)') "  N  ", &
     " TF CENT",&
     " TCL SRF",  &
     " FL ENTH"

       do ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End

          WRITE(unit,'(1x, I3,": ", 4F8.4)') &
         nlz_core( ns ) , &
         (dist_th_col(0, ns, i)*dist_th_scaling_factor(i), &
            i=1, N_OUT_TH_DIST)
       end do

      CALL OUTput_Write_Separator(unit)

      RETURN
      END   

      SUBROUTINE lead_properties( temp, ro, alfa, a, cp)
!======================================================================!
! Lead properties taken from Kirillov, Yuriev, Bobkov 1984 p. 227-228  !
! formula for conductivity below 550 C is changed to provide the       !
!  continuous function                                                 !
! (c) Slava  20 June 2001 MEPhI                                        ! 
!======================================================================!
      REAL*8, INTENT(IN)  :: temp
      REAL*8, INTENT(OUT) :: ro , alfa, a, cp
! temp - temperature (input in Kelvin, used in Celsius
!      ro - coolant density [Kg/m^3]
!      alfa - thermal conductivity  [ W/m*K ]
!      a -  - thermal diffusivity  [m^2/s]
!      cp - Specific Heat [J/(Kg K)]
      REAL*8 :: temp_c

      temp_c = temp - 273.15 

      ro   = 11072 -1.2*temp_c
      a = 1.E-06*(14.125 + temp_c*(- 2.11E-02 +2.5E-05*temp_c) )
      IF(temp_c.GT. 550.) THEN
        alfa = 34.89 + temp_c*( -0.073 + 6.9E-05*temp_c)
      ELSE
        alfa = 14.588888+temp_c*( 0.000822222 +1.88888888888E-6*temp_c)
      END IF                 
        cp = 147.3
      
      RETURN
      END       

      SUBROUTINE  OUT_Map_Flow_Rate_Core(unit, Flow_Rate_Core)
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:  
      INTEGER unit
      REAL  Flow_Rate_Core(NP_Reactor_Core)
! Local :
      INTEGER      np, n_c
      CHARACTER*80 Header_Map
      CHARACTER*4 val_fmt
      CHARACTER*6 val_char(0:N_POLY)
      REAL Flow_Rate_Core_Average

      Flow_Rate_Core_Average = 0. 
      DO np=1, NP_Reactor_Core
         Flow_Rate_Core_Average = Flow_Rate_Core_Average +&
           Flow_Rate_Core(np)
      END DO
         Flow_Rate_Core_Average = Flow_Rate_Core_Average/&
           NP_Reactor_Core
      WRITE(unit, '(A,  E12.5, A)') &
     " Average Coolant Mass Flow Rate per Assembly :", &
       Flow_Rate_Core_Average,  " [kg/(m^2*s)]"

      val_fmt = "A6"
      val_char(:) = ""

      DO n_c = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(n_c)
          WRITE(val_char(np), '(F6.3)')  Flow_Rate_Core(n_c)/&
                          Flow_Rate_Core_Average
      END DO

      Header_Map = " Relative Coolant Mass Flow Rate per Assembly :" 
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      RETURN
      END 

      SUBROUTINE  TH_OUT_Map_INTEGER(unit, TH_TYPE_CORE, Header_Map)
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:  
      INTEGER unit
      INTEGER  TH_TYPE_CORE(NP_Reactor_Core)
      CHARACTER*(*) Header_Map
! Local :
      INTEGER      np, n_c
      CHARACTER*4 val_fmt
      CHARACTER*6 val_char(0:N_POLY)
      
      val_fmt = "A6"
      val_char(:) = ""

      DO n_c = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(n_c)
          WRITE(val_char(np), '(I6)')  TH_TYPE_CORE(n_c)
      END DO

      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      RETURN
      END 


      SUBROUTINE TH_COMPUTE_AVR_FROD(vol_ass, Numb_Reactor_Core, &
            hz_core)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'
! Input 
      INTEGER Numb_Reactor_Core(NP_Reactor_Core)
      REAL vol_ass(N_POLY, NZR), hz_core

      REAL Convert_CM3_M3, Convert_CM_M
      PARAMETER ( Convert_CM3_M3 = 1.E-06, Convert_CM_M=1.E-02)

      REAL weighting_factor(NN_FRD_FA, NN_FA_TYPE), vol_fuel,&
        Part_Power_frd, S_Fuel_per_Assembly_Total, vol_fuel_3d

      INTEGER k, i, j, ns, np, i_type_fa


      Temp_FR_AVR(:,:,:,:) = 0.
      weighting_factor(:,:)  = 0.
      POWER_FRD_AVR(:, :, :) = 0.
      POWER_FRD_TOTAL = 0.

!      write(*,*) 'Cool_Heating=', Cool_Heating
!      write(*,*) 'vol_ass=', vol_ass(1,:)
!      write(*,*) 'Numb_Reactor_Core=', Numb_Reactor_Core(:)
!      write(*,*) 'hz_core=', hz_core
!      pause

! Temp_FR_AVR(NN_FRD_TOTAL+2, NN_FRD_FA, NZR, NN_FA_TYPE) 
      do k = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(k) 
         i_type_fa = th_type_core(k)


! Computing the total volume of the fuel
         S_Fuel_per_Assembly_Total = 0.
         DO j = 1, N_FRD_FA(i_type_fa) 
           S_Fuel_per_Assembly_Total = S_Fuel_per_Assembly_Total + &
               S_Fuel_per_Assembly(j, i_type_fa)
         END DO 


         DO j = 1, N_FRD_FA(i_type_fa) 

! Axial Layer is not important (FA area is constant)             
               vol_fuel=&
                        S_Fuel_per_Assembly(j, i_type_fa)*vol_ass(np,1)
               weighting_factor(j, i_type_fa) = &
                    weighting_factor(j, i_type_fa)+vol_fuel

         DO ns = 1, NZR 

            Part_Power_frd = Pow_Core(ns,k)*(1. - Cool_Heating)*&
               Pow_FRD_Assm(j, i_type_fa)/S_Fuel_per_Assembly_Total

            vol_fuel_3d =&
                    S_Fuel_per_Assembly(j, i_type_fa)*vol_ass(np,ns)*&
           Convert_CM3_M3

             POWER_FRD_AVR(j, ns, i_type_fa) = &
                   POWER_FRD_AVR(j, ns, i_type_fa)+ &
               Part_Power_frd*vol_fuel_3d 

               DO i = 1, NN_FRD_TOTAL+2
                 Temp_FR_AVR(i, j, ns, i_type_fa) =&
                    Temp_FR_AVR(i, j, ns, i_type_fa)+&
                      Temp_FR(i, j, ns, k)*vol_fuel                    
               END DO
! Coolant 
               i=NN_FRD_TOTAL+3
                Temp_FR_AVR(i, j, ns, i_type_fa)=&
               Temp_FR_AVR(i, j, ns, i_type_fa)+vol_fuel*&
                      Temp_Cool(ns, k)                                          
         END DO
         END DO
      END DO                              

      do i_type_fa = 1, NN_FA_TYPE
         DO j = 1, N_FRD_FA(i_type_fa) 
         POWER_FRD_AVR(j, 0, i_type_fa) = 0.
         DO ns = 1, NZR 
            POWER_FRD_AVR(j, 0, i_type_fa)=&
                       POWER_FRD_AVR(j, 0, i_type_fa)+&
                    POWER_FRD_AVR(j, ns, i_type_fa)
               DO i = 1, NN_FRD_TOTAL+3
                Temp_FR_AVR(i, j, ns, i_type_fa)= &
                  Temp_FR_AVR(i, j, ns, i_type_fa)/&
                  weighting_factor(j,i_type_fa) 
               END DO
         END DO ! NZR
         POWER_FRD_TOTAL = POWER_FRD_TOTAL+&
         POWER_FRD_AVR(j, 0, i_type_fa)
         DO ns=1, NZR  
            IF( POWER_FRD_AVR(j, 0, i_type_fa).ne.0) THEN
            POWER_FRD_AVR(j, ns, i_type_fa)=&
           POWER_FRD_AVR(j, ns, i_type_fa)*hz_core*Convert_CM_M/&
                    (POWER_FRD_AVR(j, 0, i_type_fa)*hz(ns))
            ELSE
!              write(*,*) 'POWER_FRD_AVR(j, 0, i_type_fa)==0'
!              write(*,*) 'j  =', j, 'i_type_fa =', i_type_fa
!              pause
            END IF
         END DO ! NZR
         END DO
      END DO

      WRITE(*,*) 'POWER_FRD_TOTAL=', POWER_FRD_TOTAL*1.E-06
!      pause
!      write(*,*) 'POWER_FRD_AVR(j, ns, i_type_fa)=', &
!           POWER_FRD_AVR(1, :, 1)
!      write(*,*) 'hz(ns)=', hz(:) 
!      pause


      RETURN
      END

      SUBROUTINE TH_OUT_Average_Fuel_Rods(unit, NZR_Beg, NZR_End)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'

! Input
      INTEGER unit, NZR_Beg, NZR_End

! Local
      INTEGER i, j, ns, i_type_fa

      WRITE(unit,'(A,I4)') &
     "     Temperature Distribution in the Average Fuel Rods, [K]"

      WRITE(unit,'(A,ES12.5)') &
     "     Total Power Generated in the Fuel Rods,  [MWt] :",&
         POWER_FRD_TOTAL*1.E-06

      CALL OUTput_Write_Separator(unit)

      do i_type_fa = 1, NN_FA_TYPE
         WRITE(unit,'(A,I4)') "     Fuel Assembly Type: ", i_type_fa
         DO j = 1, N_FRD_FA(i_type_fa) 
         WRITE(unit,'(A,I4, A, I4)') " Fuel Rod Number ", j, " of ", &
          N_FRD_FA(i_type_fa)           

!         WRITE(*,*)    " i_type_fa, j = ", i_type_fa, j
!         WRITE(*,*)    " Power of the FRD [MWt] = ", &
!          POWER_FRD_AVR(j, 0, i_type_fa)*1.E-06
         WRITE(unit,'(A,ES12.5)')    " Power of the FRD [MWt] = ", &
          POWER_FRD_AVR(j, 0, i_type_fa)*1.E-06
         CALL OUTput_Write_Separator(unit)

         WRITE(unit, '(4x, 15I7)')  (i,  i = 1, NN_FRD_TOTAL+3) 
         WRITE(unit, '(4x, 15A7)') "  POWER", &
           ("   FUEL",  i = 1, NN_FRD_FUEL+1),&
         ("   CLAD",  i = 1, NN_FRD_CLAD+1),      "   COOL"   
         DO ns = NZR_Beg, NZR_End
!      write(*,*) 'j, ns, i_type_fa =', j, ns, i_type_fa 
!      write(*,*) POWER_FRD_AVR(j, ns, i_type_fa)
!      pause
         WRITE(unit, '(I3,A,F7.4,15F7.1)' ) ns, ":", &
         POWER_FRD_AVR(j, ns, i_type_fa),&
        (Temp_FR_AVR(i, j, ns, i_type_fa), i = 1, NN_FRD_TOTAL+3) 
         END DO
         CALL OUTput_Write_Separator(unit)
         END DO
      END DO

      CALL OUTput_Write_Separator(unit)

      RETURN
      END



      SUBROUTINE THM_SKAZKA_Input_Data(file_data, cool_heating_skazka,&
       hz_skazka  ) 
      implicit  none
      INCLUDE 'TH_Model.fh'
!=====================================================================*
! Input T/H Parameters from the file  "TH_Core.DAT"                   *
!=====================================================================*
! Input:
      character*(*) file_data
      REAL cool_heating_skazka, hz_skazka(NZR)

! Local:
      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      logical error_find

!       REAL  Rel_Flow_Area(NN_FA_TYPE), &
!       S_Fuel_per_Assembly(NN_FRD_FA, NN_FA_TYPE)
!      
!      INTEGER  TH_TYPE_CORE(NP_REACTOR_CORE), N_FRD_FA(NN_FA_TYPE)  

! Local:
      INTEGER i, j, np, ios     



!initialization of the identifiers
      WRITE(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'

      open(io_unit,file = file_data ,status='old', iostat=ios)
      
!TH_N_FA_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_N_FA_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_N_FA_TYPE in the FILE_INPUT'//&
          'set TH_N_FA_TYPE = 1')
         
         N_FA_TYPE = 1

      else
   
       read(io_unit,  fmt=*, iostat=ios) N_FA_TYPE

       call Iostat_Error_Check(ios,"Error Reading "// &
      "N_FA_TYPE from the FILE_INPUT file")

      end if

!TH_TYPE_CORE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_TYPE_CORE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_Type_Core in the FILE_INPUT' //&
         ' set TH_TYPE_CORE =1' )
         DO np = 1, NP_Reactor_Core
            TH_TYPE_CORE(np) = 1
         END DO
      else
   
       read(io_unit,  fmt=*, iostat=ios) (TH_TYPE_CORE(i),i=1, &
          NP_Reactor_Core)

       call Iostat_Error_Check(ios,"Error Reading "// &
      " Fuel Assembly types under identifier TH_TYPE_CORE"//&
      " from the FILE_INPUT file")

      end if

!TH_NFRD_ASSM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_NFRD_ASSM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_NFRD_ASSM in the FILE_INPUT'//&
         ' set N_FRD_FA = 1')
         DO i = 1, N_FA_TYPE
            N_FRD_FA(i) = 1
         END DO              

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (N_FRD_FA(i),i=1, N_FA_TYPE) 
        call Iostat_Error_Check(ios,"Error Reading "// &
      "number of the equivalent fuel rods in FA" //&
      "under idfentifer TH_NFRD_ASSM from the FILE_INPUT file")

      end if

!TH_FRD_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_FRD_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_FRD_TYPE in the FILE_INPUT'//&
         ' set FRD_TYPE_FA = 1')
         DO i = 1, N_FA_TYPE
            DO j = 1, N_FRD_FA(i)
               FRD_FA_TYPE(j, i) = 1
            END DO 
         END DO              

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (( FRD_FA_TYPE(j,i),j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

        call Iostat_Error_Check(ios,"Error Reading "// &
      " fuel rods types in fuel assembly under identifier" //&
      " TH_FRD_TYPE from the FILE_INPUT file")

      end if


!TH_FLOW_AREA
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_FLOW_AREA", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_FLOW_AREA in the FILE_INPUT')

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
            (Rel_Flow_Area(i),i=1, N_FA_TYPE) 
        call Iostat_Error_Check(ios,"Error Reading "// &
      "Realtive Flow Area under identifier TH_FLOW_AREA " //&
      "from the FILE_INPUT file")

      end if


!TH_RAT_ASFL

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_RAT_FLAS", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier TH_RAT_FLAS in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       (( S_Fuel_per_assembly(j,i),j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Ratio of the Fuel Volume to Assembly Volume "//&
       "under identifier TH_RAT_FLAS from the FILE_INPUT file")

      end if

!TH_POW_ASSM

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "TH_POW_ASSM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier TH_POW_ASSM in the FILE_INPUT'//&
         ' set TH_POW_ASSM = 1.')

         DO i=1, N_FA_TYPE
            DO j = 1, N_FRD_FA(i)
                Pow_FRD_Assm(j,i) = 1.
            END DO
         END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios)  &
       ((Pow_FRD_Assm(j,i) ,j=1, N_FRD_FA(i) ),i=1, N_FA_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Ratio of the Power Density in the Fuel Rod Type "//&
       "under identifier TH_POW_ASSM from the FILE_INPUT file")

      end if

! Cool_Heating

      close(io_unit)

      Cool_Heating = Cool_Heating_skazka

      hz(:) = hz_skazka(:)

      return
      end

      SUBROUTINE THM_Set_Core_Inlet(temperature)
      IMPLICIT NONE
      INCLUDE 'TH_Model.fh'
!  
      REAL, INTENT(IN) :: temperature
      REAL, parameter :: Convert_to_K = 273.15 ! Celcius into Kelvin  

        
      Temp_Cool_inlet = temperature + Convert_to_K
      
      RETURN
      END SUBROUTINE       THM_Set_Core_Inlet     
