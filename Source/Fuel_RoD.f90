      subroutine FRD_Compute_Time_Step( a_dt, Power_Dens, &
                   Temp_Cool, Heat_Trans_Coeff, Temp_FR,  Temp_Doppl,&
               Temp_Fuel_Av, Fuel_Enthalpy, Part_Power,i_type_rod, &
               k, ns )

      implicit none
      include 'Fuel_RoD.fh'
! Input: 
      real Power_Dens ! Kg/M^3 
      real Temp_Cool ! k
      real Heat_Trans_Coeff 
      real Temp_FR(NN_FRD_TOTAL+2) ! Temperature from the previous iteration
      real a_dt ! 1./Time_Step_Size  (0 for the Steady-State Calculations)
      real Part_Power
      integer i_type_rod
! DEBUG k
      INTEGER k, ns
! Output:
!     real Temp_FR(NN_FRD_TOTAL) ! New Temperature 
      real Temp_Doppl, Fuel_Enthalpy, Temp_Fuel_Av
! Local Variables
      real  Power(NN_FRD_FUEL+1)
      integer n
      real FRD_fuel_enthalpy
      external FRD_fuel_enthalpy

      do n = 1, NN_FRD_FUEL
          Power(n) = Power_Dens*Part_Power*Power_Pin(n, i_type_rod)
      end do

!      call FR_Analytic(Temp_Cool, Heat_Trans_Coeff, Power, Temp_FR, &
!           i_type_rod, k, ns)

       


      call FRD_Form_Matrix(Temp_Cool, Heat_Trans_Coeff, Power, &
                                        Temp_FR, a_dt, i_type_rod) 

!DEBUG      open(io_unit, File = 'Output_Debug/Matrix_FR.dat', status ='unknown')
!DEBUG          do n = 1, NN_FRD_TOTAL + 2
!DEBUG            write(io_unit,*) n, (Matr_FR(i,n),i=1,3), RHS_FR(n)
!DEBUG          end do
!DEBUG      close(io_unit)

      call MSC_progonka(NN_FRD_TOTAL+2, Matr_FR, RHS_FR, Temp_FR)

! Computing Average Fuel Rod Temperature (equal areas of heat conduction nodes)
      Temp_Fuel_Av = 0.5*( Temp_FR(1) + Temp_FR(NN_FRD_FUEL+1) )
      do n = 2, NN_FRD_FUEL
         Temp_Fuel_Av = Temp_Fuel_Av + Temp_FR(n)
      end do
      Temp_Fuel_Av = Temp_Fuel_Av/real(NN_FRD_FUEL)

      IF(FRD_DOPL_TYPE.eq."LIN") THEN
          Temp_Doppl = (1. - Alfa_Doppl)*Temp_FR(1) + &
                          Alfa_Doppl*Temp_FR(NN_FRD_FUEL+1)
      ELSE
          Temp_Doppl = Temp_Fuel_Av
      END IF

!      IF(k.eq.1.and.ns.eq.17) THEN
!      write(*,*) 'Computed Temperature Cladding, IN, OUT'
!      write(*,*)  Temp_FR(NN_FRD_FUEL+2), Temp_FR(NN_FRD_TOTAL + 2)
!      write(*,*) 'Computed Temperature Fuel  IN, OUT'
!      write(*,*)   Temp_FR(1), Temp_FR(NN_FRD_FUEL+1)
!      END IF


!DEBUG call FR_Output(Temp_Cool, Heat_Trans_Coeff, Temp_FR, i_type_rod)


!Case3      Fuel_Enthalpy = FRD_fuel_enthalpy(Temp_Doppl, i_type_rod)

!Case1
      Fuel_Enthalpy = FRD_fuel_enthalpy(Temp_Fuel_Av, i_type_rod)

!Case2      Fuel_Enthalpy = 0.5*(  FRD_fuel_enthalpy( Temp_FR(1), i_type_rod ) + 
!Case2     &             FRD_fuel_enthalpy( Temp_FR(NN_FRD_FUEL+1), i_type_rod ) )

!Case2      do n = 2, NN_FRD_FUEL
!Case2         Fuel_Enthalpy = Fuel_Enthalpy + FRD_fuel_enthalpy(Temp_FR(n), i_type_rod)
!Case2      end do

!Case2      Fuel_Enthalpy = Fuel_Enthalpy/ real(NN_FRD_FUEL, i_type_rod)


      return
      end


      subroutine FRD_Read_Namelist
      implicit none
      include 'Fuel_RoD.fh'
!=====================================================================*
!        Input from the NAMELIST "THERMAL.INI"                          *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      namelist /ini_fuel_rod/ File_Fuel_Rod
        
      open(io_unit,file='Input/SKETCH.INI',status='OLD')
        read(io_unit, NML = ini_fuel_rod)
      close(io_unit)

      return
      end         


      Subroutine FRD_Form_Matrix(Temp_Cool, Heat_Trans_Coeff, Power, &
                                        Temp_FR, a_dt, i_type_rod) 

      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp_Cool ! Temperature of the Coolant
      real Heat_Trans_Coeff
      real Power(NN_FRD_FUEL)
      real Temp_FR(NN_FRD_TOTAL + 2)
      real a_dt ! 1. / Time Step Size
      INTEGER i_type_rod

!     real Conduct_gap
! Output:
!     RHS_FR(NN_FRD_TOTAL) 
!     Matr_FR(NN_FRD_TOTAL)            
! Local Variables
      real alfa1, alfa2,  Coupl, vol, Time_Der
      integer n, nc
      REAL GapConductance, Temp_Lead, Power_Density
! Functions
      real FRD_Cond_Clad, FRD_Cond_Fuel, FRD_Sp_Heat_Fuel, &
             FRD_Sp_Heat_Clad, FRD_Gap_Conductance

! Starting with the 1st Node of thew Fuel

      alfa1 = FRD_Cond_Fuel(Temp_FR(1), i_type_rod)
      alfa2 = FRD_Cond_Fuel(Temp_FR(2), i_type_rod)
      Matr_FR(3,1) = - 2.*alfa1*alfa2*Radius_Int(1, i_type_rod)/&
             ( Delta_RF(1, i_type_rod)*(alfa1 + alfa2) )

      Vol = 0.5*(Radius_Int(1, i_type_rod)**2 - &
                   Radius_Fuel(1, i_type_rod)**2 )
!      write(*,*) 'Radius_Fuel(1, i_type_rod) = ', &
!                     Radius_Fuel(1, i_type_rod)
       
      Time_Der = a_dt*Density_Fuel(i_type_rod)*&
           FRD_Sp_Heat_Fuel(Temp_FR(1), i_type_rod)

      Matr_FR(2,1) = - Matr_FR(3,1) + Time_Der*Vol
      Power_Density = Power(1)
      RHS_FR(1) = ( Power_Density + Time_Der*Temp_FR(1))*Vol

!     write(*,*) '1st node vol =', vol
!     pause

!  Internal Nodes
      do n = 2, NN_FRD_FUEL

         alfa1 = alfa2
         alfa2 = FRD_Cond_Fuel(Temp_FR(n+1), i_type_rod)

         Matr_FR(3,n) = - 2.*alfa1*alfa2*Radius_Int(n, i_type_rod)/&
             (DElta_RF(n, i_type_rod)* (alfa1 + alfa2))
         Matr_FR(1,n) = Matr_FR(3,n-1)

         Vol = 0.5*(Radius_Int(n, i_type_rod)**2 - &
              Radius_Int(n-1, i_type_rod)**2)
         Time_Der = a_dt*Density_Fuel(i_type_rod)*&
             FRD_Sp_Heat_Fuel(Temp_FR(n), i_type_rod)

         Matr_FR(2,n) = - Matr_FR(1,n) - Matr_FR(3,n) + Time_Der*Vol
         Power_Density = 0.5*( Power(n-1) + Power(n) ) 
         RHS_FR(n) = ( Power_Density  +&
              Time_Der*Temp_FR(n) )*Vol 
!        write(*,*) 'Fuel node vol =', vol, Radius_Int(n), Radius_Int(n-1)
!        pause

      end do
! Boundary Node of the Fuel
      n = NN_FRD_FUEL + 1

      IF( FlagComputeGapConductance ) THEN
        Temp_Lead = 0.5*( Temp_FR(n)+Temp_FR(n+1) )
        GapConductance =  FRD_Gap_Conductance&
          ( Temp_Lead, Radius_Fuel(n, i_type_rod), &
               Radius_Clad(1, i_type_rod))
!        write(*,*) 'i_type_rod =', i_type_rod
!        write(*,*) 'Temp_Lead =', Temp_Lead
!        write(*,*) 'R_Inner, R_Outer =', Radius_Fuel(n, i_type_rod),&
!       Radius_Clad(1, i_type_rod)
!        write(*,*) 'GapConductance =', GapConductance
!        PAUSE

      ELSE
        GapConductance =  Conduct_Gap(i_type_rod)
      END IF

      Matr_FR(3, n) = - Radius_Clad(1, i_type_rod)*&
              GapConductance
!          Conduct_Gap(i_type_rod)
      Matr_FR(1, n) = Matr_FR(3, n - 1)

      vol = 0.5*(Radius_Fuel(n, i_type_rod)**2 - &
                               Radius_Int(n-1, i_type_rod)**2)

!     write(*,*) 'last Fuel node vol =', vol
!     pause

      Time_Der = a_dt*Density_Fuel(i_type_rod)*&
          FRD_Sp_Heat_Fuel(Temp_FR(n), i_type_rod)

      Matr_FR(2, n) = - Matr_FR(1, n) - Matr_FR(3, n) + &
                                                 Time_Der*Vol
! Power generation in the last node of the fuel 
      Power_Density = Power(n-1)
      RHS_FR(n) = ( Power_Density + Time_Der*Temp_FR(n) )*vol
        
! 1st Node of the Cladding
      n = NN_FRD_FUEL + 1 + 1

      alfa1 = FRD_Cond_Clad(Temp_FR(n), i_type_rod) ! Cladding
      alfa2 = FRD_Cond_Clad(Temp_FR(n+1), i_type_rod) ! Cladding

      Matr_FR(1, n ) = Matr_FR(3, n - 1)
      Matr_FR(3, n ) = - 2.*alfa1*alfa2*&
              Radius_Int(NN_FRD_FUEL+1, i_type_rod)/&
             (Delta_RC(1, i_type_rod)* ( alfa1 + alfa2 ) )

      vol = 0.5*(Radius_Int(NN_FRD_FUEL+1, i_type_rod)**2-&
          Radius_Clad(1, i_type_rod)**2)

      Time_Der = a_dt*Density_Clad(i_type_rod)*&
         FRD_Sp_Heat_Clad(Temp_FR(n), i_type_rod)
      Matr_FR(2,n) = - Matr_FR(1, n) - Matr_FR(3, n) +&
                                                         Time_Der*Vol
      RHS_FR(n) = Time_Der*Temp_FR(n)*vol

!     write(*,*) '1st Clad node vol =', vol
!     pause

! Internal Cladding Nodes if any
      do nc = 2, NN_FRD_CLAD 
         n = nc + NN_FRD_FUEL + 1
         alfa1 = alfa2
         alfa2 = FRD_Cond_Clad(Temp_FR(n+1), i_type_rod ) ! Cladding
         Matr_FR(3, n) = - 2.*alfa1*alfa2*&
             Radius_Int(NN_FRD_FUEL + nc, i_type_rod)/&
             (Delta_RC(nc, i_type_rod)* (alfa1 + alfa2) )
         Matr_FR(1, n) = Matr_FR(3, n-1)
         vol = &
         0.5*(Radius_Int(NN_FRD_FUEL + nc, i_type_rod)**2-&
         Radius_Int(NN_FRD_FUEL+nc-1, i_type_rod)**2)
         Time_Der = a_dt*Density_Clad(i_type_rod)*&
            FRD_Sp_Heat_Clad(Temp_FR(n), i_type_rod)
         Matr_FR(2, n) = - Matr_FR(1, n) - &
            Matr_FR(3, n) + Time_Der*Vol
         RHS_FR(n) = Time_Der*Temp_FR(n)*vol
!       write(*,*) 'Clad nodes vol =', vol
!       pause

      end do
! Boundary Node of the Cladding

      n = NN_FRD_TOTAL + 2
      Matr_FR(1, n) = Matr_FR(3, n-1)
      Coupl = Radius_Clad(NN_FRD_CLAD+1, i_type_rod)*Heat_Trans_Coeff

      vol =  0.5*(Radius_Clad(NN_FRD_CLAD+1, i_type_rod)**2 - &
               Radius_Int(NN_FRD_FUEL+NN_FRD_CLAD, i_type_rod)**2)

!     write(*,*) 'Last Clad nodes vol =', vol
!     pause

      Time_Der = a_dt*Density_Clad(i_type_rod)*&
                    FRD_Sp_Heat_Clad(Temp_FR(n), i_type_rod)

      Matr_FR(2, n) = Coupl - Matr_FR(1, n) + &
                       Time_Der*Vol
      RHS_FR(n) = Coupl*Temp_Cool + Time_Der*Vol*Temp_FR(n)

      return
      end

      real function FRD_Cond_Fuel(Temp, i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp ! K
      INTEGER i_type_rod
! Output Fuel Conductiwity (watts/cm*K)
!     Real FRD_Cond_Fuel

      FRD_Cond_Fuel = const_alfa_fuel(1, i_type_rod) + &
        const_alfa_fuel(2, i_type_rod)/&
        (Temp + const_alfa_fuel(3, i_type_rod) )

      return
      end

      real function FRD_Sp_Heat_Fuel(Temp, i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp ! K
      INTEGER i_type_rod
! Output Fuel Specific Heat Capacity (J/KG*K)
!     Real 
! Local 

      FRD_Sp_Heat_Fuel = (( const_cp_fuel(4, i_type_rod)*Temp + &
        const_cp_fuel(3, i_type_rod))*&
        temp + const_cp_fuel(2, i_type_rod) )*temp + &
        const_cp_fuel(1, i_type_rod)

      return
      end

      real function FRD_Cond_Clad(Temp, i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp ! K
      INTEGER i_type_rod
! Output Cladding Conductiwity (watts/cm*K)
!     real FRD_Cond_Clad 
! Local variables 

      FRD_Cond_Clad = ((const_alfa_clad(4, i_type_rod)*Temp + &
        const_alfa_clad(3, i_type_rod))*&
        temp + const_alfa_clad(2, i_type_rod))*temp + &
        const_alfa_clad(1, i_type_rod)

      return
      end

      real function FRD_Sp_Heat_Clad(Temp, i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp ! K
      INTEGER i_type_rod
! Output Cladding Specific Heat Capacity (J/KG*K)
!     Real 
! Local 

      FRD_Sp_Heat_Clad = const_cp_clad(1, i_type_rod) + &
         const_cp_clad(2, i_type_rod)*Temp 

      return
      end


      Subroutine FR_Analytic(Temp_Cool, Heat_Trans_Coeff, Power, &
                    Temp_FR, i_type_rod, k, ns)
      implicit none
      include 'Fuel_RoD.fh'
! Input:
      real Temp_Cool ! Temperature of the Coolant
      real Heat_Trans_Coeff
      real Power(NN_FRD_FUEL + 1)
      real Temp_FR(NN_FRD_TOTAL + 2)
      INTEGER i_type_rod
!     real Conduct_gap
! DEBUG
      INTEGER k, ns
! Output:
      real Temp_Clad_In ! Temperature at the Inner Surface of the Caldding
      real Temp_Clad_Out
      real Temp_Fuel_Out ! Temperature at the Surface of the Fuel
      real Temp_Fuel_CL ! Centerline Temperature of the Fuel
! Local Variables
      real alfa,  Heat_Flux, Resist_Clad, Resist_Fuel, Resist_Gap,&
             Resist_Film
! Functions
      real FRD_Cond_Clad, FRD_Cond_Fuel

      IF(k.eq.1.and.ns.eq.17) THEN  
! resistance fuel
      alfa = FRD_Cond_Fuel(Temp_FR(1), i_type_rod)
      Resist_Fuel = Radius_Fuel(NN_FRD_FUEL+1, i_type_rod) /(2.*alfa)

! resistance clad
      alfa = FRD_Cond_Clad(Temp_FR(NN_FRD_FUEL+1), i_type_rod)
      Resist_Clad=log(Radius_Clad(NN_FRD_CLAD+1, i_type_rod)/&
                  Radius_Clad(1, i_type_rod))*&
                  Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)/(alfa)
! resistance gap
      Resist_Gap = Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)/&
          (Conduct_Gap(i_type_rod)*Radius_Clad(1, i_type_rod))
! resistance film
      Resist_Film = Radius_Fuel(NN_FRD_FUEL+1, i_type_rod) / &
         ( Heat_Trans_Coeff*Radius_Clad(NN_FRD_CLAD+1, i_type_rod) )

! heat generated in the fuel 
      Heat_Flux = Power(1)*Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)/2.

! Inner Cladding Temp
      Temp_Clad_Out = Temp_Cool + Heat_Flux*Resist_film
      Temp_Clad_In = Temp_Clad_Out + Heat_Flux*Resist_Clad
      Temp_Fuel_Out = Temp_Clad_Out + &
             Heat_Flux*(Resist_Clad+Resist_Gap)
      Temp_Fuel_CL = Temp_Clad_Out + &
                    Heat_Flux*(Resist_Clad+Resist_Gap+Resist_Fuel)

    
      write(*,*) 'Coolant Temperature =', Temp_Cool
      write(*,*) 'ANALYTICAL Temperature Cladding, IN, OUT'
      write(*,*)  Temp_Clad_In, Temp_Clad_Out
      write(*,*) 'ANALYTICAL Temperature Fuel  IN, OUT'
      write(*,*)   Temp_Fuel_CL, Temp_Fuel_Out
      write(*,*) 'Analytic Heat Flux =', Heat_Flux
      write(*,*) 'Total Generated Heat =',  Heat_Flux&
                      *2.*pi*Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)
      write(*,*) 'Heat Flux at the Outer Surface ', Heat_Flux*&
      Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)/&
      Radius_Clad(NN_FRD_CLAD+1, i_type_rod)

      END IF

      return
      end



      subroutine FR_Output(Temp_Cool, Heat_Trans_Coeff, Temp_FR, &
                             i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'

      real Temp_Cool, Heat_Trans_Coeff, Temp_FR(NN_FRD_TOTAL +2)
      INTEGER i_type_rod
      integer n

      real  Heat_Flux
         

      open(io_unit, File = 'Output/Temp_FR.dat', Status = 'Unknown')
         do n = 1, NN_FRD_FUEL + 1
            write(io_unit,*) Radius_Fuel(n, i_type_rod), &
                Temp_FR(n)
         end do
         do n = 1, NN_FRD_CLAD + 1
            write(io_unit,*) Radius_Clad(n, i_type_rod), &
                            Temp_FR(n + NN_FRD_FUEL+1)
         end do
         do n = 2, NN_FRD_FUEL
            write(io_unit,*) 'Radius_Int(n), Radius_Fuel(n+1)'
            write(io_unit,*) &
                 n, 0.5*(Radius_Int(n,i_type_rod)**2 - &
                               Radius_Int(n-1, i_type_rod)**2),&
                          Radius_Int(n, i_type_rod),  &
                          Radius_Int(n-1, i_type_rod)
          end do

          do n = 1, NN_FRD_CLAD
            write(io_unit,*) &
               'Radius_Clad(n+1), Radius_int(NN_FRD_FUEL+n)'
            write(io_unit,*) n,  Radius_int(NN_FRD_FUEL+n, i_type_rod), &
               Delta_RC(n, i_type_rod), Radius_Clad(n+1, i_type_rod)
          end do

      close(io_unit)

! checking comparing with the analytic solution

         Heat_Flux=(Temp_FR(NN_FRD_FUEL+1) - Temp_FR(NN_FRD_FUEL+2))*&
         Conduct_Gap(i_type_rod)*Radius_Clad(1, i_type_rod)/&
          Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)
         write(*,*) 'Heat Fux at the Surface of the Fuel', Heat_Flux
         write(*,*) 'Total Power Generated in the Fuel =',&
            Heat_Flux*2.*pi*Radius_Fuel(NN_FRD_FUEL+1, i_type_rod)

         Heat_Flux = (Temp_FR(NN_FRD_FUEL+NN_FRD_CLAD+2) - Temp_Cool)*&
                        Heat_Trans_Coeff
         write(*,*) 'Heat_Flux at the Outer Surface FD =', Heat_Flux
         write(*,*) 'Total Heat Out of the Pin =', Heat_Flux*&
               2.*pi*Radius_Clad(NN_FRD_CLAD+1, i_type_rod)

          write(*,*) 'Coolant Temperature =', Temp_Cool
          write(*,*) ' Finite-Difference Temperature Cladding IN, OUT'
          write(*,*)  Temp_FR(NN_FRD_FUEL+2), Temp_FR(NN_FRD_TOTAL+2)
          write(*,*) ' Finite-Difference  Temperature Fuel, IN, OUT '
          write(*,*)  Temp_FR(1), Temp_FR(NN_FRD_FUEL+1)
!         pause



      return
      end

      real function FRD_fuel_enthalpy(temp_fuel, i_type_rod)
      implicit none
      include 'Fuel_RoD.fh'
      real temp_fuel
      INTEGER i_type_rod
!     
      FRD_fuel_enthalpy = ((( 0.25*const_cp_fuel(4, i_type_rod)*&
          temp_fuel + &
          0.3333333*const_cp_fuel(3, i_type_rod))*temp_fuel + &
          0.5*const_cp_fuel(2, i_type_rod))*&
           temp_fuel + const_cp_fuel(1, i_type_rod))*temp_fuel
!      write(*,*) 'i_type_rod =', i_type_rod 
!      write(*,*) 'const_cp_fuel(:, i_type_rod)=', 
!      write(*,*) 'temp_fuel =' , temp_fuel
!      write(*,*) 'FRD_fuel_enthalpy=', FRD_fuel_enthalpy
!      pause
       
      return
      end


      subroutine FRD_Input_data(file_data) 
      implicit none
      include 'Fuel_RoD.fh'
! Input: 
!  real Delta_RF(NN_FRD_FUEL), DElta_RC(NN_FRD_CLAD), Delta_Gap from the file 
!                   "FR_Geometry.DAT"
! Output:
!  real radius_Fuel(NN_FRD_FUEL+1), Radius_CL(NN_FRD_CLAD+1)
! Input:
      character*(*) file_data

! Local Variables
      integer n, i
      real area_1, s_clad_inner, s_clad_outer, &
                  s_fuel_inner, s_fuel_outer

!      call FRD_Read_Namelist
! Local
      integer ios

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      logical error_find


!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'

      open(io_unit,file = file_data ,status='old', iostat=ios)

!FRD_DOPL_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_DOPL_TYP", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier FRD_DOPL_TYPE in the FILE_INPUT'//&
         ' set FRD_DOPL_TYPE = "LIN" ')
          FRD_DOPL_TYPE = "LIN"
      else
   
       read(io_unit,  fmt=*, iostat=ios) FRD_DOPL_TYPE

       call Iostat_Error_Check(ios,"Error Reading "// &
      "FRD_DOPL_TYPE under identifier FRD_DOPL_TYPE "//&
      "from the FILE_INPUT file")

      end if

      IF( FRD_DOPL_TYPE.eq."LIN" ) THEN
!FRD_ALF_DOPL

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_ALF_DOPL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier FRD_ALF_DOPL in the FILE_INPUT'//&
         ' set FRD_ALF_DOPL = 0.7 ')
         Alfa_Doppl = 0.7 
      else
   
       read(io_unit,  fmt=*, iostat=ios) Alfa_Doppl 

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Parameter to compute Doppler Temperature"//&
       "from the FILE_INPUT file")

      end if


      END IF ! FRD_DOPL_TYPE.eq."LIN"
!FRD_NROD_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_NROD_TYP", input_line, fmt_inp_ident, error_find)  


      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier FRD_NROD_TYPE in the FILE_INPUT'//&
         ' set N_FRD_TYPE=1')
          N_FRD_TYPE = 1
      else
   
       read(io_unit,  fmt=*, iostat=ios) N_FRD_TYPE

       call Iostat_Error_Check(ios,"Error Reading "// &
      "N_FRD_TYPE under identifier FRD_NROD_TYPE "//&
      "from the FILE_INPUT file")

      end if


!FRD_RAD_PELT
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_RAD_PELT", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier FRD_RAD_PELT in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) &
          ( Radius_Fuel(NN_FRD_FUEL+1,i), i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
      "Radius of the Fuel Pellet from the FILE_INPUT file")

      end if

!FRD_RAD_HOLE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_RAD_HOLE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier FRD_RAD_HOLE in the FILE_INPUT'// &
         ' Set Inner Radius of the fuel pellet to 0')
           DO i=1, N_FRD_TYPE
              Radius_Fuel(1,i) = 0.
           END DO
      else
   
       read(io_unit,  fmt=*, iostat=ios) (Radius_Fuel(1,i),&
          i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
      "Inner Radius of the Fuel Pellet from the FILE_INPUT file")

      end if

!FRD_RAD_CLAD
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_RAD_CLAD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier FRD_RAD_CLAD in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) (Radius_Clad(1,i), &
          Radius_Clad(NN_FRD_CLAD+1,i),i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
     "Radius of the Cladding (Inner, Outer) the FILE_INPUT file ")

      end if

!FRD_DNS_F&CL

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_DNS_F&CL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier FRD_DNS_F&CL in the FILE_INPUT')
      else
   
       read(io_unit,  fmt=*, iostat=ios) &
         (Density_Fuel(i), Density_Clad(i), i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
     "Densities of Fuel and Cladding from the FILE_INPUT file")

      end if
!FRD_GAP_COND

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_GAP_COND", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         FlagComputeGapConductance = .True.

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier FRD_GAP_COND in the FILE_INPUT'//&
         'Calculating using LEAD properties ')
      else
   
       read(io_unit,  fmt=*, iostat=ios) &
         (Conduct_Gap(i), i=1, N_FRD_TYPE) 

       call Iostat_Error_Check(ios,"Error Reading "// &
     "Gap Conductance from the FILE_INPUT file")

         FlagComputeGapConductance = .False.

      end if

!FRD_POW_FROD

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_POW_FROD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier FRD_POW_FROD in the FILE_INPUT'//&
         'set constant power distribution in fuel pellet')

        DO i=1, N_FRD_TYPE
         do n = 1, NN_FRD_FUEL
            Power_Pin(n,i) = 1.
         end do
        END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
          ((Power_Pin(n,i), n = 1, NN_FRD_FUEL),i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Power Distribution in Fuel Pellet "//&
       "from the FILE_INPUT file")

      end if

! constants for functions alfa and cp on temperature
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "FRD_ALF_FUEL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier FRD_ALF_FUEL in the FILE_INPUT'//&
         'set constants to typical values in PWR')

        DO i = 1, N_FRD_TYPE
         const_alfa_fuel(1,i) = 1.05 
         const_alfa_fuel(2,i) = 2150. 
         const_alfa_fuel(3,i) = -73.15 
        END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
          ((const_alfa_fuel(n,i), n = 1, 3), i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Constants ALFA fuel "  //&
       "from the FILE_INPUT file")

      end if

      rewind(io_unit)

      call MSC_Search_Header_In_File(io_unit, &
        "FRD_ALF_CLAD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier FRD_ALF_CLAD in the FILE_INPUT'//&
         'set constants to typical values in PWR')

       DO i=1, N_FRD_TYPE
         const_alfa_clad(1,i) =  7.51
         const_alfa_clad(2,i) =  2.09E-02 
         const_alfa_clad(3,i) = -1.45E-05
         const_alfa_clad(4,i) =  7.67E-09
       END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
       ((const_alfa_clad(n, i), n = 1, 4), i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Constants ALFA Cladding "  //&
       "from the FILE_INPUT file")

      end if

      rewind(io_unit)

      call MSC_Search_Header_In_File(io_unit, &
        "FRD_CP_FUEL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier FRD_CP_FUEL in the FILE_INPUT'//&
         'set constants to typical values in PWR')

        DO i=1, N_FRD_TYPE
         const_cp_fuel(1,i) =   162.3
         const_cp_fuel(2,i) =   0.3038
         const_cp_fuel(3,i) =  -2.391E-04
         const_cp_fuel(4,i) =   6.404E-08
        END DO


      else
   
       read(io_unit,  fmt=*, iostat=ios) &
          ((const_cp_fuel(n,i), n = 1, 4), i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Constants CP Fuel "  //&
       "from the FILE_INPUT file")

      end if

      rewind(io_unit)

      call MSC_Search_Header_In_File(io_unit, &
        "FRD_CP_CLAD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier FRD_CP_CLAD in the FILE_INPUT'//&
         'set constants to typical values in PWR')

        DO i=1, N_FRD_TYPE
         const_cp_clad(1,i) = 252.54
         const_cp_clad(2,i) =   0.11474
        END DO

      else
   
       read(io_unit,  fmt=*, iostat=ios) &
          ((const_cp_clad(n, i), n = 1, 2),i=1, N_FRD_TYPE)

       call Iostat_Error_Check(ios,"Error Reading "// &
       "Constants ALFA Cladding " //&
       "from the FILE_INPUT file")

      end if


      close(io_unit)

! fuel
      DO i = 1, N_FRD_TYPE
      s_fuel_inner = pi*Radius_Fuel(1,i)**2
      s_fuel_outer = pi*Radius_Fuel(NN_FRD_FUEL+1,i)**2
      area_1 = ( s_fuel_outer - s_fuel_inner )/NN_FRD_FUEL
      do n = 2, NN_FRD_FUEL
         Radius_Fuel(n,i) = sqrt((s_fuel_inner + area_1*(n-1))/pi)
      end do
! cladding
      s_clad_inner = pi*radius_Clad(1, i)**2
      s_clad_outer = pi*radius_Clad(NN_FRD_CLAD+1, i)**2
      area_1 = (s_clad_outer - s_clad_inner)/NN_FRD_CLAD
      do n = 2, NN_FRD_CLAD
         Radius_Clad(n,i) = sqrt((s_clad_inner + area_1*(n-1))/pi)
      end do

! Radius of the fuel and Cladding Mesh      
      do n = 1, NN_FRD_FUEL
         Delta_RF(n,i) = Radius_Fuel(n+1,i) - Radius_Fuel(n,i) 
      end do

      do n = 1, NN_FRD_CLAD
         Delta_RC(n,i) = Radius_Clad(n+1,i) - Radius_Clad(n,i) 
      end do
! Internal Radius between the nodes
      do n = 1, NN_FRD_FUEL
         Radius_Int(n,i) = 0.5*(Radius_Fuel(n+1,i) + Radius_Fuel(n,i))
      end do
      do n = 1, NN_FRD_CLAD
         Radius_Int(NN_FRD_FUEL+n,i) = &
                    0.5*(Radius_Clad(n+1,i) + Radius_Clad(n,i))
      end do
      END DO

      do n = 1, NN_FRD_TOTAL+2
         RHS_FR(n) = 0.
         do i = 1, 3
            Matr_FR(i, n) = 0.
         end do
      end do

!      Delta_RF(1) = Radius_Fuel(NN_FRD_FUEL+1)/NN_FRD_FUEL
!      Radius_Fuel(1) = 0.
!      do n = 2, NN_FRD_FUEL
!         Delta_RF(n) = Delta_RF(1)
!         Radius_Fuel(n) = Radius_Fuel(n-1) + Delta_RF(n-1)
!      end do
!      Delta_RC(1) = (Radius_Clad(NN_FRD_CLAD+1) - Radius_Clad(1))/NN_FRD_CLAD
!      do n = 2, NN_FRD_CLAD
!         Delta_RC(n) = Delta_RC(1)
!         Radius_Clad(n) = Radius_Clad(n-1) + Delta_RC(n-1)
!      end do
 
      return 
      end

      SUBROUTINE FRD_Output_Data(unit)
      IMPLICIT NONE
      INCLUDE 'Fuel_RoD.fh'

! input
      integer unit
! locals
      INTEGER n, i

      WRITE(unit, '(A)') &
       "     Method to Compute the Doppler Fuel Temperature       : "
      WRITE(unit, '(8x,A)')  FRD_DOPL_TYPE

! "     Interpolation Constant to Compute Doppler Temperature:"
      IF( FRD_DOPL_TYPE.eq."LIN" ) THEN
      WRITE(unit, '(A)') &
       "     Interpolation Constant to Compute Doppler Temperature : "
      WRITE(unit, '(8x,2E12.5)')  Alfa_Doppl
      END IF

      WRITE(unit, '(A)') " Data of Fuel Rod Model:"
! " Number of the Fuel Rod Types"
      WRITE(unit, '(A)') &
       "     Number of the Fuel Rod Types                         : "
      WRITE(unit, '(8x,I5)') N_FRD_TYPE

! "Radius of the Fuel Pellet, [m]"
      WRITE(unit, '(A)') &
       "     Inner Radius of the Fuel Pellet, [m]                  : "
      WRITE(unit, '(8x,E12.5)') (Radius_Fuel(1,i),i=1, N_FRD_TYPE)

! "Radius of the Fuel Pellet, [m]"
      WRITE(unit, '(A)') &
       "     Outer Radius of the Fuel Pellet, [m]                  : "
      WRITE(unit, '(8x,E12.5)') (Radius_Fuel(NN_FRD_FUEL+1,i)&
        ,i=1, N_FRD_TYPE)

! "Inner and Outer Radius of the Cladding, [m]"
      WRITE(unit, '(A)') &
       "     Inner and Outer Radius of the Cladding, [m]           : "
      WRITE(unit, '(8x,2E12.5)') (Radius_Clad(1, i), &
          Radius_Clad(NN_FRD_CLAD+1,i) ,i=1, N_FRD_TYPE)

! "Density of the Fuel and Cladding      , [kg/m^3]"
      WRITE(unit, '(A)') &
       "     Density of the Fuel and Cladding      , [kg/m^3]      : "
      WRITE(unit, '(8x,2E12.5)') (Density_Fuel(i), Density_Clad(i),&
       i=1, N_FRD_TYPE)

! "Fuel Rod Gap Conductance               , [W/(m^2*K)]"
      WRITE(unit, '(A)') &
       "     Fuel Rod Gap Conductance               , [W/(m^2*K)]  : "
      WRITE(unit, '(8x,E12.5)') (Conduct_Gap(i),i=1, N_FRD_TYPE)

! "     Form Function of the Power Distribution in  Fuel Pin  : "
      WRITE(unit, '(A)') &
       "     Form Function of the Power Distribution in  Fuel Pin  : "
      DO i = 1, N_FRD_TYPE
      WRITE(unit, '(8x,20F6.3)') (Power_Pin(n,i), n=1, NN_FRD_FUEL)
      END DO

! "     Constants of the fuel heat conductivity as function of   "
      WRITE(unit, '(A)') &
       "     Constants of the fuel heat conductivity as function of  "
      WRITE(unit, '(A)') &
       "           temperature alfa(T)  = C(1) + C(2)/(T+C(3))     : "
      DO i = 1, N_FRD_TYPE
      WRITE(unit, '(8x,5E14.5)') (const_alfa_fuel(n,i), n = 1, 3)
      END DO

! "     Constants of cladding heat conductivity as function of  "
      WRITE(unit, '(A)') &
       "     Constants of cladding heat conductivity as function of  "
      WRITE(unit, '(A)') &
       "      temperature alfa(T)  = C(1)+C(2)*T+C(3)*T^2+C(4)*T^3 : "
      DO i = 1, N_FRD_TYPE
      WRITE(unit, '(8x,5E14.5)') (const_alfa_clad(n,i), n = 1, 4)
      END DO

! "     Constants of fuel specific heat capacity as function of  "
      WRITE(unit, '(A)') &
       "     Constants of fuel specific heat capacity as function of "
      WRITE(unit, '(A)') &
       "      temperature cp(T)  = C(1)+C(2)*T+C(3)*T^2+C(4)*T^3   : "
      DO i = 1, N_FRD_TYPE
      WRITE(unit, '(8x,5E14.5)') (const_cp_fuel(n,i), n = 1, 4)
      END DO

! "     Constants of cladding specific heat capacity as function of  "
      WRITE(unit, '(A)') &
       "     Constants of cladding specific heat capacity as         "
      WRITE(unit, '(A)') &
       "       function of temperature cp(T)  = C(1)+C(2)*T        : "
      DO i = 1, N_FRD_TYPE
      WRITE(unit, '(8x,5E14.5)') (const_cp_clad(n,i), n = 1, 2)
      END DO

      call OUTput_Write_Separator(unit)

      RETURN
      END

      REAL FUNCTION FRD_Gap_Conductance&
                                ( Temp_Lead, R_Inner, R_Outer )
!======================================================================!
! Calculation of the Conductivity of the Lead Sublayer for BREST-300   !
!======================================================================!
      IMPLICIT NONE
      REAL, INTENT(IN) :: Temp_Lead, R_Inner, R_Outer
      REAL  :: Temp_c, alfa
        
      Temp_c=Temp_Lead-273.15 

      IF(temp_c.GT. 550.) THEN
        alfa = 34.89 + temp_c*( -0.073 + 6.9E-05*temp_c)
      ELSE
        alfa = 14.588888+temp_c*( 0.000822222 +1.88888888888E-6*temp_c)
      END IF
!        write(*,*) 'alfa =', alfa, 'temp_c =', temp_c
        FRD_Gap_Conductance=(alfa/R_Outer)/LOG(R_Outer/R_Inner)
!        write(*,*) 'FRD_Gap_Conductance=', FRD_Gap_Conductance

      RETURN
      END                        
