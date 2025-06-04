      subroutine DEC_Heat_Update_Precursors
!=====================================================================*
! Calculation of delayed neutron precursors and Fission Source        *
!    at the End of Time Step                                          *
! (c) Slava 15.IV.1998                                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'         

! Input: Source_dt(N_TOT), Flux(NG, N_TOT), XS_SF(NG< N_TOT), 
!        P_DH(MH,NH,NZ), ecn(MD), eco(MD), volume(N_TOT), ec(MD)
! Output: Source(N_TOT) -Fission Source Term at the Current Time Step
!         P_DH(MH,NH,NZ) - Concentration of te Delayed Neutron at
!                             the Current Time Step
! Local Variables:
      integer n1, m, k

      DO n1 = 1, NZ  
       do k = 1, NH
          do m = 1, MH
              p_dh(m, k, n1) = p_dh(m, k,n1)*edh(m) + (p(k,n1)*edhn(m) + &
                   p_dt(k,n1)*edho(m))
         end do
      end do
      END DO

      return
      end


      subroutine DEC_Heat_Update_Precursors_Euler_
!=====================================================================*
! Calculation of delayed neutron precursors and Fission Source        *
!    at the End of Time Step                                          *
! (c) Slava 15.IV.1998                                                *
!  неявная схема эйлера для эмиттеров остаточного энерговыделения     *
!=====================================================================*
      implicit none
      include 'sketch.fh'         

! Input: Source_dt(N_TOT), Flux(NG, N_TOT), XS_SF(NG< N_TOT), 
!        P_DH(MH,NH,NZ), ecn(MD), eco(MD), volume(N_TOT), ec(MD)
! Output: Source(N_TOT) -Fission Source Term at the Current Time Step
!         P_DH(MH,NH,NZ) - Concentration of te Delayed Neutron at
!                             the Current Time Step
! Local Variables:
      integer n1, m, k

      write(*,*) 'dt_dh =', dt_dh,  'edho(m)=', edho(1)

      DO n1 = 1, NZ  
       do k = 1, NH
          do m = 1, MH
              p_dh(m, k, n1)=(p_dh(m, k,n1)/dt_dh + p(k,n1)*ej_dc(m))/&
               ( 1./dt_dh +lam_dc(m) )
         end do
      end do
      END DO

      return
      end


      subroutine DEC_Heat_Update_Precursors_Euler
!=====================================================================*
! Calculation of delayed neutron precursors and Fission Source        *
!    at the End of Time Step                                          *
! (c) Slava 15.IV.1998                                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'         

! Input: Source_dt(N_TOT), Flux(NG, N_TOT), XS_SF(NG< N_TOT), 
!        P_DH(MH,NH,NZ), ecn(MD), eco(MD), volume(N_TOT), ec(MD)
! Output: Source(N_TOT) -Fission Source Term at the Current Time Step
!         P_DH(MH,NH,NZ) - Concentration of te Delayed Neutron at
!                             the Current Time Step
! Local Variables:
      integer n1, m, k
      real*8 expv(1:MH)

      ! write(*,*) 'dt_dh =', dt_dh,  'p_dh(m,n1,k)=', p_dh(1,100,10)
      ! write(*,*) 'p(k,n1)', p(100,10)

      do m = 1, MH
        expv(m) = DEXP(-lam_dc(m)*dt_dh)
      enddo

      DO n1 = 1, NZ  
       do k = 1, NH
          do m = 1, MH
              p_dh(m, k, n1)=ej_dc(m)*(1.d0-expv(m))*p(k,n1) + &
                             p_dh(m, k, n1)*expv(m)
         end do
      end do
      END DO

      return
      end

      subroutine DEC_heat_Init_
!**********************************************************************
!       Computing of the initial data for the decay heat calculations
!**********************************************************************
      implicit none
      include 'sketch.fh'         

      integer n1, m, k

      DO n1 = 1, NZ  
       do k = 1, NH
          do m = 1, MH
              p_dh(m, k, n1) = ej_dc(m)*p(k, n1)/lam_dc(m)
         end do
      end do
      END DO

      RETURN 
      END


      subroutine DEC_heat_Init
!**********************************************************************
!       Computing of the initial data for the decay heat calculations
!**********************************************************************
      implicit none
      include 'sketch.fh'         

      integer n1, m, k
      ! write(*,*) 'dt_dh =', dt_dh,  'p_dh(m,n1,k)=', p_dh(1,100,10)
      ! write(*,*) 'p(k,n1)', p(100,10)

      DO n1 = 1, NZ  
       do k = 1, NH
          do m = 1, MH
              p_dh(m, k, n1) = ej_dc(m)*p(k, n1)
         end do
      end do
      END DO

      RETURN 
      END

      SUBROUTINE POWer_Add_dh_
!=====================================================================*
!  Add decay heat to the total power                                  *
!    for Neutron Kinetics Calculations                                *
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
      INTEGER k, n1, n, kc, nn,  kt, m
      REAL  Pow_Core, a1

      Pow_Core = 0.

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
            p_tot(k,n1) = (1.-e_dc)*p(k,n1)*volume(kt)
               a1 = 0.
               DO m = 1, MH
                 a1 = a1 + lam_dc(m)*p_dh(m,k,n1)*volume(kt)
               END DO
               p_tot(k,n1) = p_tot(k,n1) + a1         
               Pow_Core = Pow_Core + p_tot(k,n1)
               p_tot(k,n1) = p_tot(k,n1)/volume(kt)
         END DO
      END DO

      p_total = Pow_Core*conv_Wt_to_MWt
      
      p_col(0,0) = Pow_Core/v_core


      RETURN
      END 

      SUBROUTINE POWer_Add_dh
!=====================================================================*
!  Add decay heat to the total power                                  *
!    for Neutron Kinetics Calculations                                *
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
      INTEGER k, n1, n, kc, nn,  kt, m
      REAL  Pow_Core, a1

      Pow_Core = 0.

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
            p_tot(k,n1) = (1.-e_dc)*p(k,n1)*volume(kt)
               a1 = 0.
               DO m = 1, MH
                 a1 = a1 + p_dh(m,k,n1)*volume(kt)
               END DO
               p_tot(k,n1) = p_tot(k,n1) + a1         
               Pow_Core = Pow_Core + p_tot(k,n1)
               p_tot(k,n1) = p_tot(k,n1)/volume(kt)
         END DO
      END DO

      p_total = Pow_Core*conv_Wt_to_MWt
      
      p_col(0,0) = Pow_Core/v_core


      RETURN
      END 


      SUBROUTINE POWer_Compute_Average_dh
!=====================================================================*
! computing the average values
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
      REAL weighting_factor(N_POLY)
!
      CALL MSC_SSET(N_POLY, 1.0, weighting_factor) 

      CALL OUTput_Compute_3D_Average( REAL(p_tot), NZ_Core_Beg, &
        NZ_Core_End, &
        Index_Core, p_col,  p_mm, k_p_mm, weighting_factor )

      CALL OUTput_Compute_2D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm )

      CALL OUTput_Compute_1D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm, weighting_factor )

      CALL POWer_Compute_OFFSET


      RETURN
      END 
