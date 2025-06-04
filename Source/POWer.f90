      SUBROUTINE POWer_Compute()
!=====================================================================*
! calculation of power distribution for the Steady-State Calculatione *
!    & Neutron Kinetics Calculations                                  *
! p_reactor - Total Reactor Power (Wt)                                *
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: conv_Wt_to_MWt - Coefficient to conert Wt int MWt
!        XS_SF_P(NG, N_TOT), Flux(NG, N_TOT), P_Reactor (WT),
!        NZ_Core_BEG, NZ_Core_End, NH_Core, np_Core(NH_Core), 
!        v_core
!      logical KINETICS

! Output: p(NH,NZ) - power density in the nodes (Wt/cm**3)
!         p_col(N_Poly, NZR) - Assembly Power Density (Wt/cm**3)
!         p_2D(N_POLY) - Axially-Averaged Assembly Power Density (Wt/cm**3)
!         P_Total - Total Reactor Power (MWt)
!         P_Average - Average over the Reactor Core Power Density (Wt/cm**3)
!         p3d_max - maximum power density (from p_col)
!         p2d_max - maximum power density (from pch)
!         k_p3d_max(2) - X-Y Coordinate of the Node
!                       with maximum power density (3D)
!         np_max - Axial Level of the Node with maximum power density
!         k_p2d_max(2) - X-Y Coordinate of the Assembly
!                       with maximum power density (2D)
! Pow_Norm - Constant Used for the Neutron Flux Normalization
      real conv_Wt_to_MWt
      parameter (conv_Wt_to_MWt = 1.E-06)

!     Local Variables:     
      INTEGER k, n1, n, kc, nn,  kt
      REAL  Pow_Core, a1, Pow_norm_min

      if(Problem_Type.NE."Kinetics") then

        Pow_Norm = 0.
!        v_core = 0.
!        write(*,*) 'v_core =', v_core
!        pause

        DO n1 =  NZ_Core_BEG, NZ_Core_End 
           nn = (n1 - 1)*NH
           DO kc = 1, NH_Core

              k = np_core(kc) + nn 
!              IF(XS_SF_P(NG,k).ne.0) THEN
                  DO n = 1,NG  
                  Pow_Norm = Pow_Norm + XS_SF_P(n,k)*Flux(n, k)&
                             *pow_conv(n)
                  END DO
!                  v_core = v_core + volume(k)
!              END IF
           END DO
        END DO

        Pow_Norm = P_Reactor/(Pow_Norm*conv_Wt_to_MWt)

      ELSE
        Pow_Norm = 1. ! - e_dc ! decay heat
      END IF


      if(Problem_Type.EQ."Kinetics") then
       DO n1 = 1, NZ
       DO k = 1, NH
        p_dt(k, n1) = p(k, n1)
       END DO
       END DO        


      end if ! if(Problem_Type.EQ."Kinetics") then

      Pow_Core = 0.

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
            kt = k + nn 
            p(k,n1) = 0.
            DO n = 1,NG  
               a1 = pow_conv(n)*XS_SF_P(n,kt)*&
                            Flux(n,kt)*Pow_Norm
              ! DO m = 1, MH
              !   a1 = a1 + lam_dc(m)*p_dh(m,k,n1)*volume(kt)
              ! END DO
               p(k,n1) = p(k,n1) + a1         
            END DO ! NG
            Pow_Core = Pow_Core + p(k,n1)
            p(k,n1) = p(k,n1)/volume(kt)
         END DO
      END DO

      p_total = Pow_Core*conv_Wt_to_MWt
      
      p_col(0,0) = Pow_Core/v_core


      write(*,*) 'p_total (prompt) = ', p_total
      write(*,*) 'Pow_Norm = ', pow_norm
      if(Problem_Type.NE."Kinetics") then
      call DEC_heat_Init()
      endif
      CALL POWer_Add_dh
      CALL POWer_Compute_Average_dh
      CALL OUTput_Compute_Average_Flux

      write(*,*) 'after p_total (prompt) = ', p_total


!        pause    
         

      RETURN
      END 

!       SUBROUTINE POW_Compute_SVRK_QL(time)
! !=====================================================================*
! ! computing the average values
! ! (c) Slava 4.III.1998 JAERI                                          *
! !=====================================================================*
!       IMPLICIT NONE
!       REAL :: time


!             call POW_SVRK_QL_CALC      
!         call POW_SVRK_KV7_CALC
!         call POW_compare_p_ql_with_limits(time)

!       RETURN
!       END SUBROUTINE POW_Compute_SVRK_QL



      SUBROUTINE POWer_Compute_Average
!=====================================================================*
! computing the average values
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
      REAL weighting_factor(N_POLY)
!
      CALL MSC_SSET(N_POLY, 1.0, weighting_factor) 

      CALL OUTput_Compute_3D_Average( p, NZ_Core_Beg, NZ_Core_End, &
        Index_Core, p_col,  p_mm, k_p_mm, weighting_factor )

      CALL OUTput_Compute_2D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm )

      CALL OUTput_Compute_1D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm, weighting_factor )

      CALL POWer_Compute_OFFSET


      RETURN
      END 

      SUBROUTINE POWer_Compute_OFFSET
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
     
! LOcal
      INTEGER  ns
!      REAL  pow_reactor_up, pow_reactor_down, pow_ax_offset

      IF ( NZ /= 1 ) THEN
! even
       pow_reactor_down =0.
!       write(*,*) 'n_down =', &
!        NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1  
       DO ns =  NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1
         pow_reactor_down = pow_reactor_down + p_col(0,ns)/p_col(0,0)
       END DO
       pow_reactor_up =0.
!       write(*,*) 'up =', NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
!       write(*,*) NZR_Core, (NZR_Core/2)*2
       DO ns =  NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
         pow_reactor_up = pow_reactor_up + p_col(0,ns)/p_col(0,0)
       END DO
       IF( NZR_Core .NE. (NZR_Core/2)*2 ) THEN
         ns = NZR_Core_Beg + NZR_Core/2
         pow_reactor_up = pow_reactor_up -  0.5*p_col(0,ns)/p_col(0,0)
         pow_reactor_down = pow_reactor_down + &
          0.5*p_col(0,ns)/p_col(0,0)
       END IF

       pow_ax_offset = (pow_reactor_down - pow_reactor_up)*100./&
         (pow_reactor_up+pow_reactor_down)

        ELSE
         pow_reactor_up   = 0.5
         pow_reactor_down = 0.5
         pow_ax_offset    = 1. 
      END IF !( NZ /= 1 )

      RETURN
      END SUBROUTINE POWer_Compute_OFFSET

      SUBROUTINE POWer_Output(unit)
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      CHARACTER*80 Header_Map
      CHARACTER*4 val_fmt
      CHARACTER*6 val_char(0:N_POLY)
      INTEGER ind, ns
      REAL scale_factor
! external function
      INTEGER nlz_core

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "     POWER DISTRIBUTION            "
      CALL OUTput_Write_Separator(unit)


      WRITE(unit, '(A, E12.5, A)') &
     "     Total Reactor Power                   : ", p_total, " MWt"
      WRITE(unit, '(A,  E12.5, A)') &
     "     Average Power Density                 : ", &
       p_col(0,0),  " Wt/cm^3"

      WRITE(Header_Map, '(5x, A)') &
       "NODAL POWER PEAKING FACTORS"
      scale_factor = 1./p_col(0,0)

      CALL OUTput_Distrb_Summary( unit, Header_Map, scale_factor, &
        N_POLY, NZR, p_col,  p_mm, k_p_mm)

      Header_Map = &
     "     2D RADIAL ASSEMBLY-AVERAGED POWER DENSITY"
      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  p_col(ind,0)/p_col(0,0)
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL AVERAGE POWER DENSITY"
!      CALL OUTput_Write_Separator(unit)

       DO ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End
          WRITE(unit,'(1x, I3,": ", F8.4)') &
         nlz_core( ns ) , p_col(0,ns)/p_col(0,0)
       END DO

      CALL OUTput_Write_Separator(unit)


      IF ( NZ /= 1 ) THEN
! Axial offset

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A,F8.2)') &
     "                    pow_reactor_down + pow_reactor_up =", &
           pow_reactor_down + pow_reactor_up
      WRITE(unit,'(A,2F8.2)') &
     "                    pow_reactor_down , pow_reactor_up =", &
           pow_reactor_down , pow_reactor_up
      WRITE(unit,'(A,F8.2)') &
     "POWER AXIAL OFFSET (p_down - p_up)*100./(p_up+p_down) =", &
                       pow_ax_offset
      CALL OUTput_Write_Separator(unit)

      END IF !( NZ /= 1 )

 
      RETURN
      END   


      SUBROUTINE OUTput_Compute_Average_Flux
!=====================================================================*
! computing the average values
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
!
      REAL FLUX_TMP(N_TOT)
      INTEGER n, k
      REAL weighting_factor(N_POLY)
!
      CALL MSC_SSET(N_POLY, 1.0, weighting_factor) 

      DO n = 1, NG

         DO k = 1, N_TOT
            flux_tmp(k) = flux(n, k)
         END DO

         CALL OUTput_Compute_3D_Average(flux_tmp, 1, NZ, &
        npoly, dist_flux(0,0,n), dist_flux_mm(-3, n),&
        k_dist_flux(1, -3, n), weighting_factor )

         CALL OUTput_Compute_2D_Average(dist_flux(0,0,n), &
        1,  NZR, npoly, dist_flux_mm(-3, n),&
        k_dist_flux(1, -3, n) )

         CALL OUTput_Compute_1D_Average(dist_flux(0,0,n), &
        1, NZR, npoly, dist_flux_mm(-3, n),&
        k_dist_flux(1, -3, n), weighting_factor )

      END DO

      RETURN
      END 


      subroutine POW_Neutron_Flux_Output(unit)
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      character*80 Header_Map
      character*4 val_fmt
      character*6 val_char(0:N_POLY)
      integer ind, ns, i, n
!      integer nlx_core, nly_core
! external function (in "OUTput.f")
!      integer nlz_core
!      CHARACTER*30 feedb_name(N_FEEDBACK)

      real scale_factor(NG)

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            NEUTRON FLUX"
      CALL OUTput_Write_Separator(unit)

      WRITE(unit, '(A,/,4x,30E12.5)') &
       "     Average Neutron Flux (1/cm^3)         : ",&
       (dist_flux(0,0,n), n = 1, NG)

      DO i = 1, NG
         WRITE(Header_Map, '(5x, A, I3)') &
            "Neutron Flux [1/(s cm^2)], group ", i
         scale_factor(i) = 1./dist_flux(0,0,i)
         CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
           scale_factor(i), N_POLY, NZR, &
           dist_flux(0,0,i), dist_flux_mm(-3,i), &
           k_dist_flux(1,-3,i) )

      END DO

      WRITE(unit,'(A)') &
     "            2D NEUTRON FLUX DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)

      DO i = 1, NG
         WRITE(Header_Map, '(5x, A, I3)') &
            "Neutron Flux [1/(s cm^2)], group ", i
            val_fmt = "A6"
         DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         dist_flux(ind, 0, i)*scale_factor(i)
         END DO
      
         call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
         NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
         NYR_E_Reactor, npoly, Header_Map,io_unit, val_char, val_fmt)

      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL NEUTRON FLUX DISTRIBUTIONS"
!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(1x, A3, 2x, 30I6 )') "NZ", (i, i=1, NG)

       do ns = 1, NZR ! NZ_Core_BEG, NZ_Core_End
          WRITE(unit,'(1x, I3,": ", 30F6.3)') &
         ns, &
         (dist_flux(0, ns, i)*scale_factor(i), &
          i =1, NG)
       end do

      CALL OUTput_Write_Separator(unit)

 
      RETURN
      END   
