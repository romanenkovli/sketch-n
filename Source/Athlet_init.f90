      SUBROUTINE ATHlet_Get_Core_Feedbacks
      implicit  none
      INCLUDE 'sketch.fh'

! common blocks for the couling
      real    Temp_Cool(NZR, NP_Reactor_Core), &
             Temp_Doppl(NZR, NP_Reactor_Core),&
             Ro_Cool(NZR, NP_Reactor_Core)

      common /therm002/ Temp_Cool, Temp_Doppl,&
             Ro_Cool

      real  Temp_FR(NN_FRD_TOTAL+2, NN_FRD_FA, NZR, NP_Reactor_Core) 
      common /therm003/ Temp_FR
! End Commons blocks for data exchgange

      INTEGER, PARAMETER :: N_SKIP_LINE = 3      

! Local Variables
      integer n_c, np,  ns, j, k, n1, nz_in, i
      real    a   
      real Convert_To_C, Convert_To_Gram
      parameter (Convert_to_C = 273.15) ! Kelvin into Celcius
!      parameter (Convert_to_Gram = 1.E-03) ! Kg/m^3 into g/cm^3

      OPEN(io_unit, file ="Input/Athlet.dat", status ='old')
        DO i =1, N_SKIP_LINE
          read(io_unit,*) 
        END DO
! No coord,m TF,K TC,K TM,K roM,g/cm3

      Temp_FR(:, :, :, : ) = 1.

      do n_c = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(n_c)
          do ns = 1, NZR
                do j = 1, NCHM
                   k = poly_out(np,j)
                   if(k.ne.0)  then
                    read(io_unit,*) i, a, temp_doppl(ns, n_c), a, &
                     temp_cool(ns, n_c), ro_cool(ns, n_c)
                    write(*,'(I3,5F10.3)') i, a, temp_doppl(ns, n_c), a, &
                     temp_cool(ns, n_c), ro_cool(ns, n_c)
!                    IF(n_c.eq.1) THEN
!                         WRITE(*,*) 'ns =', ns, 'Doppler Temp =', &
!                        temp_doppl(ns, n_c)
!                    PAUSE
!                    END IF
  
                   end if 
                end do ! NCHM
          end do ! NZR
       end do ! NP_Reactor_Core

      close(io_unit)

            
      do n_c = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(n_c)
          n1 = 0
          do ns = 1, NZR
             do nz_in = 1, NPZ(ns)
                n1 = n1 + 1
                do j = 1, NCHM
                   k = poly_out(np,j)
                   if(k.ne.0)  then
                    fdback(k,n1,2) = temp_cool(ns, n_c) - Convert_to_C 
                    fdback(k,n1,3) = ro_cool(ns, n_c)
                    fdback(k,n1,4) = temp_doppl(ns, n_c)

!          write(*,*) 'Conv_Cool_Temp  =', Conv_Cool_Temp
!          write(*,*) 'fdback(k,n1,2)  =', fdback(k,n1,2)
!            pause

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

      call THM_Compute_Average_Feedbacks
!      CALL THM_Prepare_Output_Dist

      return
      end 

      SUBROUTINE Athlet_Set_Core_Power
      implicit  none
      INCLUDE 'sketch.fh'

! common for the couling
      real  Pow_Core(NZR, NP_Reactor_Core )
      common /therm005/ Pow_Core

      integer ns, np, k
      real Convert_to_CI, convert_to_M3, convert_WT_kWT
      parameter (Convert_to_CI = 1.E+06) ! watts/cm^3 into watts/M^3
      parameter (Convert_to_M3 = 1.E-06) ! watts/cm^3 into watts/M^3
      parameter (convert_WT_kWT = 1.E-03) ! watts into kilowatts
            
      real Power_thermal


      Power_Thermal = 0.
      do ns = 1, NZR
         do k = 1, NP_Reactor_Core
            np = Numb_Reactor_Core(k)
            Pow_Core(ns,k) = p_col(np,ns)*vol_ass(np, ns)
          Power_Thermal = Power_thermal + Pow_Core(ns,k)
         end do
      end do

!      WRITE(*,*) 'vol_ass(np, ns) = ', vol_ass(1, 2)
!      WRITE(*,*) 'Total Thermal Power = ', Power_Thermal*1.E-6
!       pause

      OPEN(io_unit, file ="Output/SKETCH_POW.dat", status ='unknown')
     
      do ns = NZR_CORE_BEG, NZR_CORE_END
         do k = 1, NP_Reactor_Core
            WRITE(io_unit,'(I3,ES14.5)') ns-NZR_CORE_BEG+1, &
                 Pow_Core(ns,k)*convert_WT_kWT
         end do
      end do

      close(io_unit)

      return
      end
