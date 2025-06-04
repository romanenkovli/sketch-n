      SUBROUTINE Sm_Update_Steady_State
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
      REAL :: sa_micro(NG), SF_MACRO(NG), vol, flux_norm(1:NG)
      INTEGER :: n1, nn, kc, k, kt

!      conc_isotope_reactor(1:N_ISOTOPE)=0.
      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
              vol = volume(kt)
              SF_MACRO(1:NG)   = XS_SF_P(1:NG, kt)/vol
              Flux_norm(1:NG)  = Flux(1:NG, kt)*Pow_Norm
! Sm                  
                 sa_micro(1:NG) = sa_isotope(4, 1:NG, kt)
                 CALL Sm_Equilibrium( NG, SF_MACRO, Flux_norm, &
                 sa_micro, lambda_isotope(3), yields_isotope(3,kt),&
                 conc_isotope(k,n1,4), conc_isotope(k,n1,3))

                   
         END DO
      END DO
      

      RETURN
      END SUBROUTINE Sm_Update_Steady_State

      SUBROUTINE Xe_Update_Steady_State
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
      REAL :: sa_micro(NG), SF_MACRO(NG), vol, flux_norm(1:NG)
      INTEGER :: n1, nn, kc, k, kt
      REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24


      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
              vol = volume(kt)
              SF_MACRO(1:NG)   = XS_SF_P(1:NG, kt)/vol
              Flux_norm(1:NG)  = Flux(1:NG, kt)*Pow_Norm
! Xe
                sa_micro(1:NG) = sa_isotope(2, 1:NG, kt)
                CALL Xe_Equilibrium( NG, SF_MACRO, Flux_norm, &
                    sa_micro, lambda_isotope(2), lambda_isotope(1),&
                yields_isotope(2,kt), yields_isotope(1,kT),&
                conc_isotope(k,n1,2), conc_isotope(k,n1,1))


!              if ( k == 10 .and. n1 == 8 ) THEN
!                write(*,*) 'flux =', Flux_norm(1:NG)
!                write(*,*) 'conc_isotope(k,n1,2) =', &
!                 conc_isotope(k,n1,2)/Convert_Barns_to_cm2 
!                write(*,*) 'number_of_fission =',&
!                  SUM( SF_MACRO(1:NG)*FLUX_norm(1:NG) )    
!            sa_micro(1:NG) = sa_isotope(2,1:NG,kt)*Convert_Barns_to_cm2 
!               write(*,*) 'absorption_Xe=', &
!                   SUM( sa_micro(1:NG)*FLUX_NORM(1:NG) )
!               write(*,*) 'steady-state Xe'

!              end if
              
               


         END DO
      END DO

!      write(*,*) 'Conc I, Xe=', conc_isotope(1,1,2), &
!       conc_isotope(1,1,1), pow_norm
       
      
      RETURN
      END SUBROUTINE Xe_Update_Steady_State

      SUBROUTINE XE_Sm_Transient_Concentrations(dt)
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
      REAL, INTENT(IN) :: dt
      REAL :: sa_micro(NG), SF_MACRO(NG), vol, flux_norm(1:NG)
      INTEGER :: n1, nn, kc, k, kt

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
              vol = volume(kt)
              SF_MACRO(1:NG)   = XS_SF_P(1:NG, kt)/vol
              Flux_norm(1:NG)  = Flux(1:NG, kt)*Pow_Norm
! Xe
              sa_micro(1:NG) = sa_isotope(2, 1:NG, kt)
              IF(Xe_Sm_Model(1:1) == "t") THEN
                    CALL Xe_Time_Step( dt, NG, SF_MACRO, Flux_norm, sa_micro, &
                lambda_isotope(2), lambda_isotope(1),&
                yields_isotope(2,kt), yields_isotope(1,kt),&
                conc_isotope(k,n1,2), conc_isotope(k,n1,1))     
              END IF

! Sm                  
              sa_micro(1:NG) = sa_isotope(4, 1:NG, kt)
              IF(Xe_Sm_Model(2:2) == "t") THEN
!                 write(*,*) 'Sm transient' 
                 CALL Sm_Time_Step( dt, NG, SF_MACRO, Flux_norm, &
                     sa_micro, lambda_isotope(3), yields_isotope(3,kt),&
                conc_isotope(k,n1,4), conc_isotope(k,n1,3))
              END IF

                   
         END DO
      END DO
      
!                 write(*,*) 'conc Sm, Pm=', &
!                   conc_isotope(1,1,4), conc_isotope(1,1,3)
!                 pause


      RETURN
      END SUBROUTINE XE_Sm_Transient_Concentrations

      SUBROUTINE Xe_Equilibrium( NG, SF, FLUX, sa_XE, lambda_Xe,&
       lambda_I, y_Xe, y_I, N_Xe, N_I)
      IMPLICIT NONE

        INTEGER, INTENT(IN) :: NG 
        REAL, INTENT(IN)    :: SF(NG), FLUX(NG),  lambda_Xe,&
       lambda_I, y_Xe, y_I
        REAL                :: sa_Xe(NG)
        REAL, INTENT(OUT)   :: N_Xe, N_I
        REAL :: number_of_fission, absorption_Xe 
        REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24

! convert micr
        sa_Xe(1:NG) = sa_Xe(1:NG)*Convert_Barns_to_cm2 
        number_of_fission = SUM( SF(1:NG)*FLUX(1:NG) )
        absorption_Xe     = SUM( sa_Xe(1:NG)*FLUX(1:NG) )
        N_xe = (y_Xe+y_i)*number_of_fission/ &
           ( absorption_Xe +  lambda_Xe )
        N_I = y_I*number_of_fission / lambda_I 

        N_xe= N_xe*Convert_Barns_to_cm2 
        N_I = N_I*Convert_Barns_to_cm2 
      RETURN 
      END

      SUBROUTINE Xe_Time_Step( dt, NG, SF, FLUX, sa_XE, lambda_Xe,&
       lambda_I, y_Xe, y_I, N_Xe, N_I)
      IMPLICIT NONE
        INTEGER, INTENT(IN) :: NG 
        REAL, INTENT(IN)    :: SF(NG), FLUX(NG),  lambda_Xe,&
       lambda_I, y_Xe, y_I, dt
        REAL                :: sa_Xe(NG)
        REAL, INTENT(INOUT)   :: N_Xe, N_I
        REAL :: number_of_fission, absorption_Xe
        DOUBLE PRECISION :: exp_I, exp_Xe
        REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24

        sa_Xe(1:NG) = sa_Xe(1:NG)*Convert_Barns_to_cm2 
        N_xe= N_xe/Convert_Barns_to_cm2 
        N_I = N_I/Convert_Barns_to_cm2 

        exp_I = DEXP( DBLE(-lambda_I*dt) )
        absorption_Xe     = SUM( sa_Xe(1:NG)*FLUX(1:NG) )
!        write(*,*)  'absorption_Xe + lambda_Xe=', &
!        (absorption_Xe + lambda_Xe), 'dt=', dt
        exp_Xe = DEXP( DBLE(-(absorption_Xe + lambda_Xe)*dt) )

        number_of_fission = SUM( SF(1:NG)*FLUX(1:NG) )

        N_Xe = N_Xe*exp_Xe + &
        (y_Xe+y_i)*number_of_fission*(1.-exp_Xe)/&
          ( absorption_Xe +  lambda_Xe ) -&
         (y_I*number_of_fission - lambda_I*N_I)*(exp_Xe-exp_I)/&
         (lambda_I - lambda_Xe - absorption_Xe)  

        N_I = y_I*number_of_fission*( 1. - EXP_I)/ lambda_I +&
           N_I*EXP_I

        N_xe= N_xe*Convert_Barns_to_cm2 
        N_I = N_I*Convert_Barns_to_cm2 

      RETURN 
      END
        
      SUBROUTINE Sm_Equilibrium( NG, SF, FLUX, sa_Sm, lambda_Pm,&
        y_Nd, N_Sm, N_Pm)
      IMPLICIT NONE
        INTEGER, INTENT(IN) :: NG 
        REAL, INTENT(IN)    :: SF(NG), FLUX(NG),  lambda_Pm,&
        y_Nd
        REAL :: sa_Sm(NG)
        REAL, INTENT(OUT)   :: N_Sm, N_Pm
        REAL :: number_of_fission, absorption_Sm
        REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24

        sa_Sm(1:NG) = sa_Sm(1:NG)*Convert_Barns_to_cm2 

        number_of_fission = SUM( SF(1:NG)*FLUX(1:NG) )
        absorption_Sm     = SUM( sa_Sm(1:NG)*FLUX(1:NG) )
        N_Sm = y_Nd*number_of_fission/absorption_Sm 
        N_Pm = y_Nd*number_of_fission/lambda_Pm

        N_Sm= N_Sm*Convert_Barns_to_cm2 
        N_Pm = N_Pm*Convert_Barns_to_cm2 


      RETURN 
      END
       
      SUBROUTINE Sm_Time_Step( dt, NG, SF, FLUX, sa_Sm, lambda_Pm,&
       y_Nd, N_Sm, N_Pm)
      IMPLICIT NONE
        INTEGER, INTENT(IN) :: NG 
        REAL, INTENT(IN)    :: SF(NG), FLUX(NG), lambda_Pm,&
       y_Nd, dt
        REAL :: sa_Sm(NG)
        REAL, INTENT(INOUT)   :: N_Sm, N_Pm
        REAL :: number_of_fission, absorption_Sm
        DOUBLE PRECISION :: exp_Pm, exp_Sm  
        REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24

        sa_Sm(1:NG) = sa_Sm(1:NG)*Convert_Barns_to_cm2 
        N_Pm= N_Pm/Convert_Barns_to_cm2 
        N_Sm = N_Sm/Convert_Barns_to_cm2 

        exp_Pm = DEXP( DBLE(- lambda_Pm*dt) )
        absorption_Sm     = SUM( sa_Sm(1:NG)*FLUX(1:NG) )
        exp_Sm = DEXP( DBLE(- absorption_Sm*dt) )
        number_of_fission = SUM( SF(1:NG)*FLUX(1:NG) )

        N_Sm = N_Sm*exp_Sm + &
        y_Nd*number_of_fission*(1.-exp_Sm)/absorption_Sm -&
        (y_Nd*number_of_fission - lambda_Pm*N_Pm)*(exp_Sm-exp_Pm)/&
         (lambda_Pm - absorption_Sm)  

        N_Pm = y_Nd*number_of_fission*( 1. - EXP_Pm)/ lambda_Pm +&
           N_Pm*EXP_Pm

        N_Sm= N_Sm*Convert_Barns_to_cm2 
        N_Pm = N_Pm*Convert_Barns_to_cm2 


      RETURN 
      END

      REAL FUNCTION k_inf_compute( NG, SA, nuSF, S12)
      IMPLICIT NONE
      INTEGER :: NG
      REAL    :: SA(1:NG), nuSF(1:NG), S12
      REAL    :: rinf
        
      IF( NG == 2) THEN
        rinf =  S12/SA(2)
        k_inf_compute =  ( nuSF(1)+ rinf*nuSF(2))/( S12 + SA(1) )
      ELSE 
        WRITE(*,*) 'K_INF function is not implemented for NG>2'
        STOP
      END IF          

      RETURN
      END
 
      subroutine Xe_Sm_Distr_Output(unit)
!=====================================================================*
! Output of the Xe and Sm distribution into "SKETCH.lst"              *
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

      real dist_scaling_factor(N_ISOTOPE)
!      data dist_th_scaling_factor / 1.E-3, 1.E-3, 1./

!  Scaling Factor for Fuel Enthalpy = 1./dist_th_av

      dist_scaling_factor(1:N_ISOTOPE)=1.
      DO i = 1, N_ISOTOPE
        IF ( conc_isotope_col(0,0,i) > SMALL_VALUE ) THEN 
        dist_scaling_factor(i) = &
        1./conc_isotope_col(0,0,i)
        END IF
      END DO

      CALL OUTput_Write_Separator(unit)

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            ISOTOPE DISTRIBUTION"

      DO i = 1, N_ISOTOPE

      WRITE(unit,'(A,A, ES14.6, A)') &
     isotope_name(i)," Average Value of Isotope concentration =", &
         conc_isotope_col(0,0,i), " [nuclei/cm^3 x 10^24]"

! Flux(NG) + Power + FDBacks + T/H Distributions
         WRITE(Header_Map, '(5x, A)') &
           isotope_name(i)


      CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
          dist_scaling_factor(i), &
          N_POLY, NZR, conc_isotope_col(0,0,i), conc_isotope_mm(-3,i), &
          k_conc_isotope_mm(1,-3,i) )

      END DO


!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            2D ISOTOPE DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)


      DO i = 1, N_ISOTOPE
      
      WRITE(Header_Map, '(5x, A)') &
           isotope_name(i)

      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         conc_isotope_col(ind,0,i)*dist_scaling_factor(i)
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL ISOTOPE DISTRIBUTIONS "
!      CALL OUTput_Write_Separator(unit)

      WRITE(unit, '((1x, 5A7), /)') "  N  ", &
       (isotope_name(i), i=1, N_ISOTOPE)

       do ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End

          WRITE(unit,'(1x, I3,": ", 4F8.4)') &
         nlz_core( ns ) , &
         (conc_isotope_col(0,ns,i)*dist_scaling_factor(i), &
            i=1, N_ISOTOPE)
       end do

      CALL OUTput_Write_Separator(unit)


      IF ( NZ /= 1 ) THEN
! Axial offset

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A,F8.2)') 
      WRITE(unit,'(A,F8.2)') &
     "POWER AXIAL OFFSET Iod (i_down - i_up)*100./(i_up+i_down) =", &
                       conc_ax_offset_i
      WRITE(unit,'(A,F8.2)') &
     "POWER AXIAL OFFSET Xe (xe_down - xe_up)*100./(xe_up+xe_down) =", &
                  conc_ax_offset_xe
      WRITE(unit,'(A,F8.2)') 
      CALL OUTput_Write_Separator(unit)

      END IF !( NZ /= 1 )


      RETURN
      END   

      SUBROUTINE Xe_Sm_Compute_Distribution
      implicit none
      include 'sketch.fh'
      INTEGER :: i

      DO i = 1, N_ISOTOPE

       CALL Compute_Average_Core_Distr(NH, NZ, N_POLY, NZR,&
       NYR, NXR, NZ_Core_Beg, NZ_Core_End, NZR_Core_Beg, NZR_Core_End, &
       conc_isotope(1,1,i), conc_isotope_col(0,0,i), &
       conc_isotope_mm(-3,i),k_conc_isotope_mm(1,-3,i),Index_Core)

      END DO 

!      CALL Xe_Sm_Compute_OFFSET 


      RETURN
      END 


      SUBROUTINE Xe_Sm_Compute_OFFSET
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
     
! LOcal
      INTEGER  ns, kt, k, n1, nn, kc, n
!      REAL  conc_reactor_up_i, conc_reactor_down_i, conc_ax_offset_i
        REAL, PARAMETER :: Convert_Barns_to_cm2 = 1.E-24
        REAL :: SF_MACRO(1:NG), Flux_norm(1:NG), sa_micro(1:NG), vol,&
         N_xe, N_Xe_up, N_Xe_Down, number_of_Fissions, number_xe_abs
        REAL ::        conc_down_xe_nom, conc_down_xe_denom, &
                      conc_up_xe_nom, conc_up_xe_denom

      real conv_Wt_to_MWt
      parameter (conv_Wt_to_MWt = 1.E-06)

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

      write(*,*) 'Pow_norm =', pow_norm
!      pause

      IF ( NZ /= 1 ) THEN
! even

       conc_reactor_down_i =0.
!       write(*,*) 'n_down =', &
!        NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1  
       DO ns =  NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1
         conc_reactor_down_i = conc_reactor_down_i + &
               conc_isotope_col(0,ns,1)/conc_isotope_col(0,0,1)
       END DO
       conc_reactor_up_i =0.
!       write(*,*) 'up =', NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
!       write(*,*) NZR_Core, (NZR_Core/2)*2
       DO ns =  NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
         conc_reactor_up_i = conc_reactor_up_i +&
               conc_isotope_col(0,ns,1)/conc_isotope_col(0,0,1)
       END DO
       IF( NZR_Core .NE. (NZR_Core/2)*2 ) THEN
         ns = NZR_Core_Beg + NZR_Core/2
         conc_reactor_up_i = conc_reactor_up_i -&
               0.5*conc_isotope_col(0,ns,1)/conc_isotope_col(0,0,1)
         conc_reactor_down_i = conc_reactor_down_i + &
               0.5*conc_isotope_col(0,ns,1)/conc_isotope_col(0,0,1)
       END IF

       conc_ax_offset_i = (conc_reactor_down_i-conc_reactor_up_i)*100./&
         (conc_reactor_up_i+conc_reactor_down_i)

        ELSE
         conc_reactor_up_i   = 0.5
         conc_reactor_down_i = 0.5
         conc_ax_offset_i    = 1. 
      END IF !( NZ /= 1 )


      IF ( NZ /= 1 ) THEN
! even

       conc_reactor_down_xe =0.
!       write(*,*) 'n_down =', &
!        NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1  
       DO ns =  NZR_Core_Beg, NZR_Core_Beg + NZR_Core/2 - 1
         conc_reactor_down_xe = conc_reactor_down_xe + &
               conc_isotope_col(0,ns,2)/conc_isotope_col(0,0,2)
       END DO
       conc_reactor_up_xe =0.
!       write(*,*) 'up =', NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
!       write(*,*) NZR_Core, (NZR_Core/2)*2
       DO ns =  NZR_Core_Beg + NZR_Core/2,  NZR_Core_End
         conc_reactor_up_xe = conc_reactor_up_xe +&
               conc_isotope_col(0,ns,2)/conc_isotope_col(0,0,2)
       END DO
       IF( NZR_Core .NE. (NZR_Core/2)*2 ) THEN
         ns = NZR_Core_Beg + NZR_Core/2
         conc_reactor_up_xe = conc_reactor_up_xe -&
               0.5*conc_isotope_col(0,ns,2)/conc_isotope_col(0,0,2)
         conc_reactor_down_xe = conc_reactor_down_xe + &
               0.5*conc_isotope_col(0,ns,2)/conc_isotope_col(0,0,2)
       END IF

      conc_ax_offset_xe=(conc_reactor_down_xe-conc_reactor_up_xe)*100./&
         (conc_reactor_up_xe+conc_reactor_down_xe)

        ELSE
         conc_reactor_up_xe   = 0.5
         conc_reactor_down_xe = 0.5
         conc_ax_offset_xe    = 1. 
      END IF !( NZ /= 1 )


!             write(*,*)  'conc_reactor_down_xe, conc_reactor_up_xe'
!              write(*,*)  conc_reactor_down_xe, conc_reactor_up_xe
              write(*,*)  'conc_ax_offset_xe=', conc_ax_offset_xe
!              write(*,*)  conc_ax_offset_xe
!              pause 

      RETURN
      END SUBROUTINE Xe_Sm_Compute_OFFSET

