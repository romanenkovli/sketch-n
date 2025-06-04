      subroutine XSP_Compute
!=====================================================================*
! Current Value of Macro Cross-Sections due to Feedbacks              *
! All Nodes Without Control Rods                                      *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!=====================================================================*

      implicit none
      include 'sketch.fh'
! Input: NH, NZ, volume(N_TOT), np_out(NH), ns_out(NZ), l(N_POLY, NZR),
!        N_FEEDBACK, fdback(NH,NZ, N_FEEDACK), fdback00(N_FEEDBACK, NNODE),
!        NG, d(NNODE,NG), sf(NNODE,NG), sf_p(NNODE,NG), sa(NNODE,NG),
!        sik(NNODE,NG,NG), d_fb(N_FEEDBACK,NG,NNODE), 
!        sf_fb(N_FEEBACK,NG,NNODE), sik_fb(N_FEEDBACK, NG, NG,NNODE), 
!        sf_p_fb(N_FEEDBACK, NG, NNODE,), sa_fb(N_FEEDBACK, NG,NNODE) 
!        
! Output XS_D(NG, N_TOT), XS_SA(NG, N_TOT), XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),
!        XS_SF_P(NG,N_TOT) - Macro-Cross Sections* Volume (Except XS_D)
           
! Local Variables
      real delta_fb(N_FEEDBACK), vol, sa2, sf2, sf_p2,&
           sik2(NG), str2, d2
      integer np,ns,nl, i,k,n1,n, m, kt, nn, nd

!      write(*,*) 'NO SQRT OF DOPPLER TEMPERASTURE FOR WIGL PROBLEM'

       do n1 = 1, NZ 
        nn = (n1-1)*NH
        do k = 1, NH

          kt = k + nn
          vol = volume(kt)

          np = np_out(k)
          ns = ns_out(n1)
          nl = l(np,ns)

! Doppler

          delta_fb(N_FEEDBACK) = sqrt(fdback(k,n1,N_FEEDBACK)) - &
                         sqrt(fdback00(N_FEEDBACK,nl))

! Doppler Feddback for the WIGL Model
!cc          delta_fb(N_FEEDBACK) = (fdback(k,n1,N_FEEDBACK)) - 
!cc     &                   (fdback00(N_FEEDBACK,nl))

! everything else
          do i = 1, N_FEEDBACK - 1
            delta_fb(i) = fdback(k,n1,i) - fdback00(i,nl)
          end do

! ALL assembly without control rods

       do n = 1, NG
! initial values
         if(I_Diff_Coeff.ne.1) then
            str2 = 1./(3.*d(nl,n))
         ELSE
            d2 = d(nl,n)
         END IF 
         sf2 = sf(nl,n) 
         sf_p2 = sf_p(nl,n)
         do m = 1, NG
           sik2(m) = sik(nl,m,n) 
         end do          
         sa2 = sa(nl,n) 
! feedbacks

         do i = 1, N_FEEDBACK

         if(I_Diff_Coeff.ne.1) then
           str2 = str2 + d_fb(i,n,nl)*delta_fb(i)
         ELSE
           d2 = d2 + d_fb(i,n,nl)*delta_fb(i)
         END IF  
           sf2 =  sf2 + sf_fb(i,n,nl)*delta_fb(i) 
           sf_p2 = sf_p2 + sf_p_fb(i,n,nl)*delta_fb(i) 
           do m = 1, NG
             sik2(m)= sik2(m) + sik_fb(i,m,n,nl)*delta_fb(i)
           end do          
             sa2 = sa2 + sa_fb(i,n,nl)*delta_fb(i)
!          if(k.eq.6.and.n1.eq.10.and.n.eq.2) then
!             write(*,*) 'i=', i, 'sa,=', sa2,&
!                     'sa_fb(i,n,nl)*delta_fb(i)=',&
!                        sa_fb(i,n,nl)*delta_fb(i)
!          end if
!          write(*,*) 'n = ', n , 'sa_sik =', sa_sik
          end do ! FEEDBACKS
! multiplication by volume 

         if(I_Diff_Coeff.ne.1) then
            XS_D(n, kt) = 1./(3.*str2)
         ELSE
            XS_D(n, kt) = d2
         END IF
         XS_SF(n, kt) = sf2*vol
         XS_SF_P(n, kt) = sf_p2*vol
         XS_SA(n, kt) = sa2*vol
         do m = 1, NG
           XS_SIK(m,n, kt) = sik2(m)*vol
         end do    

! Assembly discontinuity factors
         DO nd = 1, NDD
            DO i = 1, 2  
            xs_adf(n, i, nd, kt) = adf(nl, n, i, nd)
            END DO
         END DO

       end do ! NG

       end do ! NH
      end do ! NZ

!      if(Debug) then
!        call XS_DBG_output
!      end if

      return
      end

      subroutine XSP_Compute_CRD
!=====================================================================*
! Current Value of Macro Cross-Sections due to Control Rods           *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!  27.IV.1998 EXCLUDING TOP AXIAL REFLECTOR FOR BWR PROBLEM           *
!  2.IX.1998 Flux-Weithing Procedure for the Bottom Part of Absorber  *
! 24.VI.1999 small revision  of the line 169 and cleaning garbage     *
!=====================================================================*
      implicit none

      include 'sketch.fh'


! Input: NN_CRod_El, NN_CRod, Mat_Com_rod(NN_CRod, NN_CRod_El), NCHM, 
!        poly_out(N_POLY, NCHM), nrods(NN_CRod), NZ_Core_Beg, NZ, 
!        NZ_Core_End, volume(N_TOT), rod_node(NN_CRod, NZ, NN_CRod_El),
!        xs_ca(NUMBER_XS, N_ROD_COMP), XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) 
! Output: XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) -
!        Corrected Value of the Macro XS Due to Control Rods 
!                                     (XS*Volume, Except XS_D)
!Local Variables
      integer ie, ir, nl, np, nch, k, n1, kt, n, m,&
        n1_rod, n1_nrod, ib
      real h_rod, h_nrod, sect_nrod, sect_rod, vol, str
      real xappa_rod(NG)
! real function
      real XS_Sect_CR

      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 

      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
         nl = Mat_Com_Rod(ir,ie)
         do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              vol = volume(kt)
              h_rod = rod_node(ir,n1,ie)
              h_nrod = 1. - h_rod
              n1_rod = n1 + 1
              n1_nrod = n1 - 1

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node

           call XS_Compute_Xappa_Rod(h_rod,e_rod_round_off, &
                 n1, n1_rod, n1_nrod, k, ir, ie, &
                 Xappa_Rod)

           do n = 1, NG 
! Diffusion Coefficients
              if(I_Diff_Coeff.ne.1) then
! differential cross sections are given for sigma transport 
              sect_nrod = 1./(3.*XS_D(n,kt))
              sect_rod  = sect_nrod + d_ca(n,nl)
              str = XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
              XS_D(n, kt) = 1./(3.*str)
              ELSE
! differential cross sections are given for diffusion coefficient
                 sect_nrod = XS_D(n,kt)
                 sect_rod  = sect_nrod + d_ca(n,nl)
                 XS_D(n, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
              END IF

! all other XS are multiplied by the node volume 
! scattering
             do m = 1, NG
              sect_nrod = XS_SIK(m, n, kt)
              sect_rod = sect_nrod  +  sik_ca(m,n,nl)*vol
                   XS_SIK(m, n, kt)=XS_Sect_CR(sect_nrod, sect_rod,&
              xappa_rod(n))
             end do ! scattering
!absorption 
              sect_nrod = XS_SA(n, kt)
              sect_rod = sect_nrod + sa_ca(n,nl)*vol
              XS_SA(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(n))

              if(n1.le.NZ_Core_END) then
! no need for nuSF of the UPPER REFLECTOR (PWR NEACRP PROBLEM)
              sect_nrod = XS_SF(n, kt)
              sect_rod = sect_nrod + sf_ca(n,nl)*vol
              XS_SF(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(n))
! Eliminating the round of the the nonfuel absorbers
!              IF(XS_SF(n,kt).LT.0) XS_SF(n,kt)=0.

              sect_nrod = XS_SF_P(n, kt)
              sect_rod = sect_nrod + sf_p_ca(n,nl)*vol
              XS_SF_P(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(n))
! Eliminating the round of the the nonfuel absorbers
!              IF(XS_SF_P(n,kt).LT.0) XS_SF_P(n,kt)=0.

             end if ! nuSF in the core 

           end do ! NG

             end if ! if there is anything in the node
           end do ! NZ
          end if ! k.ie.0
         end do ! NCHM

         end do
        end do ! NN_CRod
       end do ! NN_CRod_El

!       if(Debug) CALL XS_DBG_output

      return
      end

      subroutine XSP_Compute_Subset
!=====================================================================*
! Current Value of Macro Cross-Sections due to Feedbacks              *
! All Nodes Without Control Rods                                      *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!=====================================================================*
      USE PRECISION, ONLY : dp 
      USE termsort_lib, ONLY: ops_lib, get_lib_nx_ny, get_lib_titles
      implicit none
      include 'sketch.fh'
! Input: NH, NZ, volume(N_TOT), np_out(NH), ns_out(NZ), l(N_POLY, NZR),
!        N_FEEDBACK, fdback(NH,NZ, N_FEEDACK), fdback00(N_FEEDBACK, NNODE),
!        NG, d(NNODE,NG), sf(NNODE,NG), sf_p(NNODE,NG), sa(NNODE,NG),
!        sik(NNODE,NG,NG), d_fb(N_FEEDBACK,NG,NNODE), 
!        sf_fb(N_FEEBACK,NG,NNODE), sik_fb(N_FEEDBACK, NG, NG,NNODE), 
!        sf_p_fb(N_FEEDBACK, NG, NNODE,), sa_fb(N_FEEDBACK, NG,NNODE) 
!        
! Output XS_D(NG, N_TOT), XS_SA(NG, N_TOT), XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),
!        XS_SF_P(NG,N_TOT) - Macro-Cross Sections* Volume (Except XS_D)
           
! Local Variables
!      real delta_fb(N_FEEDBACK), vol, sa2, sf2, sf_p2,&
!           sik2(NG), str2, d2
      real vol
      integer np,ns,nl, i,k,n1,n, kt, nn, nd, i_xs_x, i_xs_y, i_xs_v,&
             i_xs_d, m
! Local Variables for XS Model SUBSET

      INTEGER, PARAMETER :: NNX=10, NNY =1000
      INTEGER            :: N_XS_X, N_XS_Y
      REAL(dp)  ::  X(NNX), Y(NNY)
      REAL(dp), PARAMETER :: Convert_Celcius_to_Kelvin = 273.15
      CHARACTER*18 :: title_x(NNX), title_y(NNY)
      CHARACTER*1000 :: Message
      REAL :: delta_cb
!      LOGICAL, PARAMETER :: flag_cb = .False.
      INTEGER :: i_isotope

      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST( YMNCNFNUM ), &
!           CONST2( YMNCNFNUM ),&
            DCONST( YMNCNFNUM * YMNCNARG )
      INTEGER KSDCA(YMNCNFNUM), nl_getera
      LOGICAL, PARAMETER ::  Flag_D2_from_Getera = .False.
      INTEGER, PARAMETER :: N_DELAYED_NEUTRON_DATA  = 13
      REAL               :: norm_beta


!      write(*,*) 'NO SQRT OF DOPPLER TEMPERASTURE FOR WIGL PROBLEM'
!       N_XS_X =  6
!       N_XS_Y = 30

       do n1 = 1, NZ 
        nn = (n1-1)*NH
        do k = 1, NH

          kt = k + nn
          vol = volume(kt)

          np = np_out(k)
          ns = ns_out(n1)
          nl = l(np,ns)

          CALL get_lib_nx_ny(nl, N_XS_X, N_XS_Y)

!          write(*,*) 'nl =', nl

          CALL get_lib_titles(nl, title_x, title_y)

!          write(*,*) 'nl =', nl, title_x

        
         DO i_xs_x = 1, N_XS_X 
          IF( INDEX( TRIM(title_x(i_xs_x)), "burnup" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = brn(k,n1) ! burnup
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "\rho" ) &
                  .GT. 0  ) THEN 
!\rho(g/cm^3)            T_c(K)            T_f(K)          c_b(ppm)
            X(i_xs_x) = fdback(k,n1,3) ! coolant density
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_c" ) &
                  .GT. 0  ) THEN 
!          write(*,*) 'fdback(k,n1,2) =', fdback(k,n1,2)
!            pause
            X(i_xs_x) = fdback(k,n1,2) + Convert_Celcius_to_Kelvin !   T_c
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_f" ) &
                  .GT. 0  ) THEN 
          X(i_xs_x) = fdback(k,n1,4) !+ Convert_Celcius_to_Kelvin !   T_f
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "c_b" ) &
                  .GT. 0  ) THEN 
!            IF( nl == 111 ) THEN
!             write(*,*) 'nl = 111'
!             write(*,*) 'c_b =', fdback(k,n1,1)
!             pause
!            END IF 
             
          X(i_xs_x) = fdback(k,n1,1) ! C_Bor
          ELSE
            WRITE(Message, '(A,I3,A)') &
            " Error in the neutron XS library, "//&
            "Fuel Type =", nl, "FA state variable "//&
              TRIM(title_x(i_xs_x))//" is unknown"//&
            " possible choice "//&
            "burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
!           write(*,*) TRIM(title_x(i_xs_x))
           write(*,*) " burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
            CALL MSC_ERR_Add_Message( 'ERROR',TRIM(Message) )
         END IF
         END DO ! i_xs_x = 1, N_XS_X 
         if(N_XS_X<=0 .or. N_XS_Y<=0) stop 11111

         CALL ops_lib(N_XS_X, N_XS_Y, x, nl, y)

! ALL assembly without control rods
          i_xs_y = 0
          IF(i_Diff_Coeff == 1) THEN
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_D(n, kt)  = Y(i_xs_y)
           END DO
          ELSE
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_D(n, kt)  = 1./Y(i_xs_y)
           END DO
          END IF 

          
!            XS_D(1, kt) = D(nl, 1) 
!            XS_D(2, kt) = D(nl, 2) 
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_SA(n, kt) = Y(i_xs_y)*vol
           END DO

!            XS_SA(2, kt) = Y(4)*vol
!            XS_SA(1, kt) = sa(nl,1)*vol
!            XS_SA(2, kt) = sa(nl,2)*vol

           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_SF(n, kt) = Y(i_xs_y)*vol
           END DO


           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_SF_P(n, kt) = Y(i_xs_y)*vol
           END DO

!            XS_SF_P(1, kt) = sf_p(nl,1)*vol
!            XS_SF_P(2, kt) = sf_p(nl,2)*vol
           IF( NG == 2 ) THEN
            i_xs_y = i_xs_y + 1
            XS_SIK(2, 1, kt) = Y(i_xs_y)*vol
!            XS_SIK(2, 1, kt) =  sik(nl,2,1)*vol
           ELSE
            WRITE(*,*) 'Scattering matrix for multigroup case' 
            STOP  
           END IF
! 
           i_xs_y = i_xs_y + 1
           k_inf_lib(kt) = Y(i_xs_y)

! neutron velocity 26, 27
!           write(*,*)  'I_flag_vel_Library =', I_flag_vel_Library
           IF( I_flag_vel_Library == 1 ) THEN
           IF ( N_XS_Y == 30 ) THEN 
             i_xs_v = 25
           ELSE IF ( N_XS_Y == 12) THEN
             i_xs_v = 10
           ELSE
             WRITE(*,*) 'LIbrary does not contains neutron velocity'
             WRITE(*,*) 'for fuel type nl = ', nl
             STOP
           END IF
           DO n = 1, NG
            i_xs_v = i_xs_v + 1
            xs_al(n, kt) = 1./(Y(i_xs_v)*100.)
!            IF( (n == 1).AND.(kt ==1) ) THEN
!                write(*,*) 'n, kt =', n, kt
!                write(*,*) 'xs_al(n, kt) = ', xs_al(n, kt)
!                pause
!            END IF
           END DO

           END IF ! IF( I_flag_vel_Library == 1 ) THEN

           IF( I_flag_deln_Library == 1 ) THEN
           IF ( N_XS_Y == 30 ) THEN 
             i_xs_d = 12   

             i_xs_d = i_xs_d + 1
             xs_bet(kt) = Y(i_xs_d)
             norm_beta = 0.
             DO m = 1, MD
               i_xs_d = i_xs_d + 1
               xs_beta(m, kt) = Y(i_xs_d)
               norm_beta = norm_beta + Y(i_xs_d)
             END DO

!           WRITE(*,*) 'norm_beta =', norm_beta
!             pause

             DO m = 1, MD
                xs_beta(m,kt) = xs_bet(kt)*xs_beta(m, kt)/norm_beta
             END DO
 
             DO m = 1, MD
               i_xs_d = i_xs_d + 1
               xs_alfa(m, kt) = Y(i_xs_d)
!               IF(  Y(i_xs_d) > 3 ) THEN
!                 write(*,*) 'kt =',kt,'m=',m,'Y(i_xs_d)=',Y(i_xs_d)
!                 pause
!               END IF
             END DO

!           write(*,*) 'alfa =', xs_alfa(1:MD,kt)


           END IF ! IF ( N_XS_Y == 30 ) THEN 
           END IF ! IF( I_flag_deln_Library == 1 ) THEN


 
! xe + sm
! BURNUP NODE 
           IF( N_XS_Y .GT. 12 ) THEN ! tmp for the old library
           IF ( k_inf_lib(kt) .GT. eps_round_off ) THEN
           i_isotope = 2
           i_xs_y = i_xs_y + 1
           sa_isotope(i_isotope, NG, kt)=Y(i_xs_y)
!           write(*,*) 'i_xs_y=', i_xs_y, 'sa_isotope(2=', Y(i_xs_y)

           i_isotope = 4
           i_xs_y = i_xs_y + 1
           sa_isotope(i_isotope, NG, kt)=Y(i_xs_y)
!           write(*,*) 'i_xs_y=', i_xs_y, 'sa_isotope(4=', Y(i_xs_y)

          IF(File_DMP_Out_Kin=="".AND.Problem_Type/="Kinetics") THEN
              i_xs_y = i_xs_y + N_DELAYED_NEUTRON_DATA + NG ! velocity(NG)
          ELSE
              i_xs_y = i_xs_y + N_DELAYED_NEUTRON_DATA + NG ! velocity(NG)
!            WRITE(*,*) 'Reading Delayed Neutron Data not DONE' 
!            STOP  
         END IF


! yields xe + sm               
! xe + sm
           i_isotope = 1
           i_xs_y = i_xs_y + 1
           yields_isotope(i_isotope,  kt)=Y(i_xs_y)
!           write(*,*) 'i_xs_y=',i_xs_y,'yields_isotope(1=', Y(i_xs_y)

           i_isotope = 2
           i_xs_y = i_xs_y + 1
           yields_isotope(i_isotope,  kt)=Y(i_xs_y)
!           write(*,*) 'i_xs_y=',i_xs_y,'yields_isotope(2=', Y(i_xs_y)

           i_isotope = 3
           i_xs_y = i_xs_y + 1
           yields_isotope(i_isotope,  kt)=Y(i_xs_y)
!           write(*,*) 'i_xs_y=',i_xs_y,'yields_isotope(3=', Y(i_xs_y)
!           pause

          END IF ! k_inf_lib(kt) .GT. eps_round_off 
          END IF ! ( N_XS_Y .GT. 12 ) THEN ! tmp for the old library

          IF (XS_Model.EQ."SUBSED" ) THEN  
 
               ARGUM(1)=fdback(k,n1,3)*1000.
                ARGUM(2)=fdback(k,n1,4)
             ARGUM(3)=fdback(k,n1,2)+ Convert_Celcius_to_Kelvin  !   T_c
             ARGUM(4)=0.
             ARGUM(5)=brn(k,n1) ! burnup
             ARGUM(6)=fdback(k,n1,1)              
             SELECT CASE (nl)
               CASE (1)
                nl_getera = -1
               CASE (2)
                nl_getera = -2
               CASE (3)
                nl_getera = -3
               CASE (9)
                nl_getera = -4
               CASE (16) 
                nl_getera = 11
               CASE (17) 
                nl_getera = 12
               CASE (30) 
                nl_getera = 1
               CASE (31) 
                nl_getera = 2
               CASE (42) 
                nl_getera = 7
               CASE (44) 
                nl_getera = 3
               CASE DEFAULT
                 WRITE(*,*) 'unknown fuel type, nl=', nl
                 STOP
               END SELECT      

           IF( nl_getera .GT. 0 ) THEN
!              write(*,*) 'nl_getera=', nl_getera
              CALL YMNCCALC(ARGUM, nl_getera, CONST, DCONST, KSDCA)
           ELSE 
              CALL YNREFsection ( ARGUM, nl_getera, CONST)
           END IF 

! ALL assembly without control rods

            XS_D(1, kt)  = const(1)
            IF (Flag_D2_from_Getera) THEN
                XS_D(2, kt)  = const(9)
            END IF

       END IF ! (XS_Model.EQ."SUBSED" )


!            write(*,*) 'nl =', nl
!            write(*,*) 'Diff (old) = ', D(nl, 1:2) 
!            write(*,*) 'Diff (new) = ', Y(1:2)
!            write(*,*) 'sa (old) = ', sa(nl, 1:2) 
!            write(*,*) 'sa (new) = ', Y(3:4)
!            write(*,*) 'sf (old) = ', sf(nl, 1:2) 
!            write(*,*) 'sf (new) = ', Y(5:6)
!            write(*,*) 'sf_p (old) = ', sf_p(nl, 1:2) 
!            write(*,*) 'sf_p (new) = ', Y(7:8)
!            write(*,*) 's1->2 (old) = ', sik(nl,2,1)
!            write(*,*) 's1->2 (new) = ', Y(9)
!            pause


! Assembly discontinuity factors
        DO n = 1, NG
         DO nd = 1, NDD
            DO i = 1, 2  
            xs_adf(n, i, nd, kt) = 1.
            END DO
         END DO
       END DO

!       end do ! NG

       end do ! NH
      end do ! NZ

!      if(Debug) then
!        call XS_DBG_output
!      end if

      return
      end


      subroutine XSP_Compute_CRD_Subset
!=====================================================================*
! Current Value of Macro Cross-Sections due to Control Rods           *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!  27.IV.1998 EXCLUDING TOP AXIAL REFLECTOR FOR BWR PROBLEM           *
!  2.IX.1998 Flux-Weithing Procedure for the Bottom Part of Absorber  *
! 24.VI.1999 small revision  of the line 169 and cleaning garbage     *
!=====================================================================*
      USE PRECISION, ONLY : dp 
      USE termsort_lib, ONLY: ops_lib, get_lib_nx_ny, get_lib_titles

      implicit none
      include 'sketch.fh'

! Input: NN_CRod_El, NN_CRod, Mat_Com_rod(NN_CRod, NN_CRod_El), NCHM, 
!        poly_out(N_POLY, NCHM), nrods(NN_CRod), NZ_Core_Beg, NZ, 
!        NZ_Core_End, volume(N_TOT), rod_node(NN_CRod, NZ, NN_CRod_El),
!        xs_ca(NUMBER_XS, N_ROD_COMP), XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) 
! Output: XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) -
!        Corrected Value of the Macro XS Due to Control Rods 
!                                     (XS*Volume, Except XS_D)
!Local Variables
      integer ie, ir, nl, np, nch, k, n1, kt,&
        n1_rod, n1_nrod, ib, i_xs_x, i_xs_y, n, i_xs_v,&
        i_xs_d, m
      real h_rod, h_nrod, sect_nrod, sect_rod, vol
      real xappa_rod(NG), a
! real function
      real XS_Sect_CR
      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 

      INTEGER, PARAMETER :: NNX=20, NNY =1000
      INTEGER            :: N_XS_X, N_XS_Y
      REAL(dp)  ::  X(NNX), Y(NNY)
      REAL(dp), PARAMETER :: Convert_Celcius_to_Kelvin = 273.15
      CHARACTER*18 :: title_x(NNX), title_y(NNY)
      CHARACTER*1000 :: Message
      INTEGER :: i_isotope

      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST( YMNCNFNUM ),&
!            CONST2( YMNCNFNUM ),&
            DCONST( YMNCNFNUM * YMNCNARG )
      INTEGER KSDCA(YMNCNFNUM), nl_getera
      LOGICAL, PARAMETER ::  Flag_D2_from_Getera = .False.
      INTEGER, PARAMETER :: N_DELAYED_NEUTRON_DATA  = 13
   
!      REAL  ::  s_12,sr2,sr1,nuSF2,nuSF1,rinf,KINF

!      write(*,*) 'NO SQRT OF DOPPLER TEMPERASTURE FOR WIGL PROBLEM'
!       N_XS_X =  6
!       N_XS_Y = 30

      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
         nl = Mat_Com_Rod(ir,ie)
         do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              vol = volume(kt)
              h_rod = rod_node(ir,n1,ie)
              h_nrod = 1. - h_rod
              n1_rod = n1 + 1
              n1_nrod = n1 - 1

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node

           call XS_Compute_Xappa_Rod(h_rod,e_rod_round_off, &
                 n1, n1_rod, n1_nrod, k, ir, ie, &
                 Xappa_Rod)

          CALL get_lib_nx_ny(nl, N_XS_X, N_XS_Y)

          CALL get_lib_titles(nl, title_x, title_y)
        
         DO i_xs_x = 1, N_XS_X 
          IF( INDEX( TRIM(title_x(i_xs_x)), "burnup" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = brn(k,n1) ! burnup
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "\rho" ) &
                  .GT. 0  ) THEN 
!\rho(g/cm^3)            T_c(K)            T_f(K)          c_b(ppm)
            X(i_xs_x) = fdback(k,n1,3) ! coolant density
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_c" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback(k,n1,2) + Convert_Celcius_to_Kelvin !   T_c
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_f" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback(k,n1,4) !+ Convert_Celcius_to_Kelvin !   T_f
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "c_b" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback(k,n1,1) ! C_Bor
          ELSE
            WRITE(Message, '(A,I3,A)') &
            " Error in the neutron XS library, "//&
            "Fuel Type =", nl, "FA state variable "//&
              TRIM(title_x(i_xs_x))//" is unknown"//&
            " possible choice "//&
            " burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
           write(*,*) TRIM(title_x(i_xs_x))
           write(*,*) "burnup,\rho(g/cm^3), T_c(K), T_f(K), c_b(ppm)"
            CALL MSC_ERR_Add_Message( 'ERROR',TRIM(Message) )
         END IF
         END DO ! i_xs_x = 1, N_XS_X 

         CALL ops_lib(N_XS_X, N_XS_Y, x, nl, y)

!           do n = 1, NG 
! Diffusion Coefficients
!              if(I_Diff_Coeff.ne.1) then
! differential cross sections are given for sigma transport 
!              sect_nrod = 1./(3.*XS_D(n,kt))
!              sect_rod  = sect_nrod + d_ca(n,nl)
!              str = XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
!              XS_D(n, kt) = 1./(3.*str)
!              ELSE
! differential cross sections are given for diffusion coefficient

          i_xs_y = 0
          IF(i_Diff_Coeff == 1) THEN
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            sect_nrod = XS_D(n,kt)
            sect_rod  = Y(i_xs_y)
            XS_D(n, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
           END DO
          ELSE
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            XS_D(n, kt)  = 1./Y(i_xs_y)
            sect_nrod = XS_D(n,kt)
            sect_rod  = 1./Y(n)
            XS_D(n, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
           END DO
          END IF 

          IF (XS_Model.EQ."SUBSED" ) THEN  

             ARGUM(1)=fdback(k,n1,3)*1000.
                ARGUM(2)=fdback(k,n1,4)
             ARGUM(3)=fdback(k,n1,2)+ Convert_Celcius_to_Kelvin  !   T_c
             ARGUM(4)=0.
             ARGUM(5)=brn(k,n1) ! burnup
             ARGUM(6)=fdback(k,n1,1)              

!           write(*,*) 'xappa_rod =', xappa_rod

             SELECT CASE (nl)
               CASE (1)
                nl_getera = -1
               CASE (2)
                nl_getera = -2
               CASE (3)
                nl_getera = -3
               CASE (9)
                nl_getera = -4
               CASE (16) 
                nl_getera = 11
               CASE (17) 
                nl_getera = 12
               CASE (30) 
                nl_getera = 1
               CASE (31) 
                nl_getera = 2
               CASE (42) 
                nl_getera = 7
               CASE (44) 
                nl_getera = 3
               CASE DEFAULT
                 WRITE(*,*) 'unknown fuel type, nl=', nl
                 STOP
                END SELECT      


           IF(nl_getera .GT. 0 ) THEN
              CALL YMNCCALC(ARGUM, nl_getera, CONST, DCONST, KSDCA)
           ELSE 
              CALL YNREFsection ( ARGUM, nl_getera, CONST)
           END IF 

! ALL assembly without control rods

! differential cross sections are given for diffusion coefficient
                 sect_nrod = XS_D(1,kt)
                 sect_rod  = const(1)
                 XS_D(1, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(1))

            IF (Flag_D2_from_Getera) THEN
                 sect_nrod = XS_D(2,kt)
                 sect_rod  = const(9)
                 XS_D(2, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
            END IF
           END IF ! xs_mode == SUBSED


! all other XS are multiplied by the node volume 
! scattering
!             end do ! scattering
!absorption 
           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            sect_nrod = XS_SA(n, kt)
            sect_rod = y(i_xs_y)*vol
            XS_SA(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(n))
           END DO

           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            IF( n1.le.NZ_Core_END ) THEN  
! no need for nuSF of the UPPER REFLECTOR (PWR NEACRP PROBLEM)
            sect_nrod = XS_SF(n, kt)
            sect_rod = y(i_xs_y)*vol
            XS_SF(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(n))
            END IF 
           END DO

           DO n = 1, NG
            i_xs_y = i_xs_y + 1
            IF( n1.le.NZ_Core_END ) THEN  
! no need for SF of the UPPER REFLECTOR (PWR NEACRP PROBLEM)
            sect_nrod = XS_SF_P(n, kt)
            sect_rod = y(i_xs_y)*vol
            XS_SF_P(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(n))
            END IF 
           END DO

           IF( NG == 2 ) THEN
            i_xs_y = i_xs_y + 1
              sect_nrod = XS_SIK(2, 1, kt)
              sect_rod =  y(i_xs_y)*vol
              XS_SIK(2, 1, kt)=XS_Sect_CR(sect_nrod, sect_rod,&
              xappa_rod(1))
           ELSE
            WRITE(*,*) 'Scattering matrix for multigroup case' 
            STOP  
           END IF

! 
           i_xs_y = i_xs_y + 1
           sect_nrod = k_inf_lib(kt) 
             sect_rod  = Y(i_xs_y)
           k_inf_lib(kt)=XS_Sect_CR(sect_nrod, sect_rod,&
              xappa_rod(NG))



! neutron velocity 26, 27
           IF( I_flag_vel_Library == 1 ) THEN
           IF ( N_XS_Y == 30 ) THEN 
             i_xs_v = 25
           ELSE IF ( N_XS_Y == 12) THEN
             i_xs_v = 10
           ELSE
             WRITE(*,*) 'LIbrary does not contains neutron velocity'
             WRITE(*,*) 'for fuel type nl = ', nl
             STOP
           END IF
           DO n = 1, NG
            i_xs_v = i_xs_v + 1
            sect_nrod = xs_al(n, kt)
            sect_rod = 1./(y(i_xs_v)*100.)
            xs_al(n, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(n))
           END DO

           END IF


           IF( I_flag_deln_Library == 1 ) THEN
           IF ( N_XS_Y == 30 ) THEN 
             i_xs_d = 12   

             i_xs_d = i_xs_d + 1

             sect_nrod = xs_bet(kt)
             sect_rod =  y(i_xs_d)
             xs_bet(kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(NG))

             DO m = 1, MD
               i_xs_d = i_xs_d + 1

               sect_nrod = xs_beta(m,kt)
               sect_rod = y(i_xs_d)*xs_bet(kt)
               xs_beta(m,kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(NG))
             END DO
          
             xs_bet(kt) = sum(xs_beta(1:MD,kt))
!             write(*,*) 'bet = ', xs_bet(kt)
!             write(*,*) 'SUM(beta) = ', sum(xs_beta(1:MD,kt))
!             pause

             DO m = 1, MD
               i_xs_d = i_xs_d + 1
               sect_nrod = 1./xs_alfa(m,kt)
               sect_rod = 1./y(i_xs_d)
               a = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(NG))
               xs_alfa(m,kt) = 1./a
             END DO

           END IF ! IF ( N_XS_Y == 30 ) THEN 
           END IF ! IF( I_flag_deln_Library == 1 ) THEN



           IF( N_XS_Y .GT. 12 ) THEN ! tmp for the old library
! xe + sm
! BURNUP NODE 
           IF ( k_inf_lib(kt) .GT. eps_round_off ) THEN
            i_isotope = 2
            i_xs_y = i_xs_y + 1
            sect_nrod = sa_isotope(i_isotope, NG, kt)
            sect_rod=Y(i_xs_y)
            sa_isotope(i_isotope, NG, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(NG))

            i_isotope = 4
            i_xs_y = i_xs_y + 1
            sect_nrod = sa_isotope(i_isotope, NG, kt)
            sect_rod=Y(i_xs_y)
            sa_isotope(i_isotope, NG, kt)=&
            XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(NG))
             
          IF(File_DMP_Out_Kin=="".AND.Problem_Type/="Kinetics") THEN
              i_xs_y = i_xs_y + N_DELAYED_NEUTRON_DATA + NG ! velocity(NG)
          ELSE
!            WRITE(*,*) 'Reading Delayed Neutron Data not DONE' 
!            STOP  
          END IF


! yields xe + sm               
! xe + sm
            i_isotope = 1
            i_xs_y = i_xs_y + 1
            sect_nrod = yields_isotope(i_isotope, kt)
            sect_rod=Y(i_xs_y)
            yields_isotope(i_isotope, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod,xappa_rod(NG))

            i_isotope = 2
            i_xs_y = i_xs_y + 1
            sect_nrod = yields_isotope(i_isotope, kt)
            sect_rod=Y(i_xs_y)
            yields_isotope(i_isotope, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod,xappa_rod(NG))

            i_isotope = 3
            i_xs_y = i_xs_y + 1
            sect_nrod = yields_isotope(i_isotope, kt)
            sect_rod=Y(i_xs_y)
            yields_isotope(i_isotope, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod,xappa_rod(NG))

          END IF ! k_inf_lib(kt) .GT. eps_round_off 
      END IF! ( N_XS_Y .GT. 12 ) THEN ! tmp for the old library 
!          s_12=XS_SIK(2,1,kt)
!          sr2=XS_SA(2,kt)
!           sr1=XS_SA(1,kt)+XS_SIK(2,1,kt)
!           nuSF2=XS_SF(2,kt)
!           nuSF1=XS_SF(1,kt)
!           rinf =  s_12/sr2
!           KINF =  (nuSF1 + rinf*nuSF2)/sr1
!           write(*,*) 'kinf from XS=', kinf
!           write(*,*) 'kinf from LIBRARY=', y(10)

!           s_12=Y(9)
!           sr1=Y(3)+Y(9)
!           sr2=Y(4)
!           nuSF1=Y(5)
!           nuSF2=Y(6)
!           rinf =  s_12/sr2
!           KINF =  (nuSF1 + rinf*nuSF2)/sr1
!           write(*,*) 'Diffusion coeff 1 ', XS_D(1,kt), 1/Y(1)
!           write(*,*) 'Diffusion coeff 2 ', XS_D(2,kt), 1/Y(2)
!           write(*,*) 'Absorption  1', XS_SA(1,kt)/vol, Y(3)  
!           write(*,*) 'Absorption  2', XS_SA(2,kt)/vol, Y(4)  
!
!           write(*,*) 'nuSF  1', XS_SF(1,kt)/vol, Y(5)  
!           write(*,*) 'nuSF  2', XS_SF(2,kt)/vol, Y(6)  
!           write(*,*) 'SF  1', XS_SF_P(1,kt)/vol, Y(7)  
!           write(*,*) 'SF  2', XS_SF_P(2,kt)/vol, Y(8)  
!           write(*,*) 'S 1->2', XS_SIK(2,1,kt)/vol, Y(9)  
!           write(*,*) 'KINF from XS library=', KINF
!           pause 

!           end do ! NG

             end if ! if there is anything in the node
           end do ! NZ
          end if ! k.ie.0
         end do ! NCHM

         end do
        end do ! NN_CRod
       end do ! NN_CRod_El

!       if(Debug) CALL XS_DBG_output

      return
      end

      subroutine XSP_Compute_LINTABLE
!=====================================================================*
! Current Value of Macro Cross-Sections due to Feedbacks              *
! All Nodes Without Control Rods                                      *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!=====================================================================*
!      USE PRECISION, ONLY : dp 
!      USE termsort_lib, ONLY: ops_lib, get_lib_nx_ny
      implicit none
      include 'sketch.fh'

      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST( YMNCNFNUM ),&
!            CONST2( YMNCNFNUM ),&
            DCONST( YMNCNFNUM * YMNCNARG )
      INTEGER KSDCA(YMNCNFNUM)

      INTEGER :: i_isotope


! Input: NH, NZ, volume(N_TOT), np_out(NH), ns_out(NZ), l(N_POLY, NZR),
!        N_FEEDBACK, fdback(NH,NZ, N_FEEDACK), fdback00(N_FEEDBACK, NNODE),
!        NG, d(NNODE,NG), sf(NNODE,NG), sf_p(NNODE,NG), sa(NNODE,NG),
!        sik(NNODE,NG,NG), d_fb(N_FEEDBACK,NG,NNODE), 
!        sf_fb(N_FEEBACK,NG,NNODE), sik_fb(N_FEEDBACK, NG, NG,NNODE), 
!        sf_p_fb(N_FEEDBACK, NG, NNODE,), sa_fb(N_FEEDBACK, NG,NNODE) 
!        
! Output XS_D(NG, N_TOT), XS_SA(NG, N_TOT), XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),
!        XS_SF_P(NG,N_TOT) - Macro-Cross Sections* Volume (Except XS_D)
           
! Local Variables
      real  vol
!          delta_fb(N_FEEDBACK), sa2, sf2, sf_p2,&
!           sik2(NG), str2, d2
      integer np,ns,nl, i,k,n1,n, kt, nn, nd
! Local Variables for XS Model SUBSET

      REAL, PARAMETER :: Convert_Celcius_to_Kelvin = 273.15


!      write(*,*) 'NO SQRT OF DOPPLER TEMPERASTURE FOR WIGL PROBLEM'
!       N_XS_X =  6
!       N_XS_Y = 30

       do n1 = 1, NZ 
        nn = (n1-1)*NH
        do k = 1, NH

          kt = k + nn
          vol = volume(kt)

          np = np_out(k)
          ns = ns_out(n1)
          nl = l(np,ns)

             ARGUM(1)=fdback(k,n1,3)*1000.
              ARGUM(2)=fdback(k,n1,4)
           ARGUM(3)=fdback(k,n1,2)+ Convert_Celcius_to_Kelvin  !   T_c
             ARGUM(4)=0.
             ARGUM(5)=brn(k,n1) ! burnup
             ARGUM(6)=fdback(k,n1,1)              

           IF( nl .GT. 0 ) THEN
              CALL YMNCCALC(ARGUM, nl, CONST, DCONST, KSDCA)
           ELSE 
              CALL YNREFsection ( ARGUM, nl, CONST)
           END IF 

! ALL assembly without control rods

          IF(i_Diff_Coeff == 1) THEN
            XS_D(1, kt)  = const(1)
            XS_D(2, kt)  = const(9)
          ELSE
            XS_D(1, kt)  = 1./const(1)
            XS_D(2, kt)  = 1./const(9)
          END IF 

!       difk(k,1)=const(1)
!       difk(k,2)=const(9)
!       signf(k,1)=const(3)
!       signf(k,2)=const(11)
!       sigp(k,1)=const(4)
!       sigp(k,2)=const(12)
!       sigr(k,1)=const(2)
!       sigr(k,2)=const(10)
!       sig12(k)=const(5)

            XS_SA(1, kt) = const(2)*vol
            XS_SA(2, kt) = const(10)*vol

            XS_SF(1, kt) = const(3)*vol
            XS_SF(2, kt) = const(11)*vol

            XS_SF_P(1, kt) = const(4)*vol
            XS_SF_P(2, kt) = const(12)*vol

            XS_SIK(2, 1, kt) = const(5)*vol

!  yields xe + sm               
!  xe + sm

! i-135
            i_isotope = 1
            yields_isotope(i_isotope,  kt)=const(15)
! xe-135
            i_isotope = 2
            yields_isotope(i_isotope,  kt)=const(16)
            sa_isotope(i_isotope, NG, kt)=const(7)
! pm-149
            i_isotope = 3
            yields_isotope(i_isotope,  kt)=const(17)
! sm-149
            i_isotope = 4
            sa_isotope(i_isotope, NG, kt)=const(8)

! kinf
!          write(*,*)
!          rinf =  XS_2G_MACRO%SCAT_SUM(1,2)/XS_2G_MACRO%SA(2)
!      XS_2G_MACRO%KINF = &
!          ( XS_2G_MACRO%nuSF(1)+ rinf*XS_2G_MACRO%nuSF(2))/&
!          ( XS_2G_MACRO%SCAT_SUM(1,2) + XS_2G_MACRO%SA(1) )


! Assembly discontinuity factors
        DO n = 1, NG
         DO nd = 1, NDD
            DO i = 1, 2  
            xs_adf(n, i, nd, kt) = 1.
            END DO
         END DO
       END DO

!       end do ! NG

       end do ! NH
      end do ! NZ

!      if(Debug) then
!        call XS_DBG_output
!      end if

      return
      end


      subroutine XSP_Compute_CRD_LINTABLE
!=====================================================================*
! Current Value of Macro Cross-Sections due to Control Rods           *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!  27.IV.1998 EXCLUDING TOP AXIAL REFLECTOR FOR BWR PROBLEM           *
!  2.IX.1998 Flux-Weithing Procedure for the Bottom Part of Absorber  *
! 24.VI.1999 small revision  of the line 169 and cleaning garbage     *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST( YMNCNFNUM ),&
!       CONST2( YMNCNFNUM ),&
            DCONST( YMNCNFNUM * YMNCNARG )
      INTEGER KSDCA(YMNCNFNUM)


! Input: NN_CRod_El, NN_CRod, Mat_Com_rod(NN_CRod, NN_CRod_El), NCHM, 
!        poly_out(N_POLY, NCHM), nrods(NN_CRod), NZ_Core_Beg, NZ, 
!        NZ_Core_End, volume(N_TOT), rod_node(NN_CRod, NZ, NN_CRod_El),
!        xs_ca(NUMBER_XS, N_ROD_COMP), XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) 
! Output: XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) -
!        Corrected Value of the Macro XS Due to Control Rods 
!                                     (XS*Volume, Except XS_D)
!Local Variables
      integer ie, ir, nl, np, nch, k, n1, kt,&
        n1_rod, n1_nrod, ib
      real h_rod, h_nrod, sect_nrod, sect_rod, vol
      real xappa_rod(NG)
! real function
      real XS_Sect_CR
      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 

      REAL, PARAMETER :: Convert_Celcius_to_Kelvin = 273.15
      INTEGER :: i_isotope

!     REAL :: rinf, kinf


!      write(*,*) 'NO SQRT OF DOPPLER TEMPERASTURE FOR WIGL PROBLEM'
!       N_XS_X =  6
!       N_XS_Y = 30

      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
         nl = Mat_Com_Rod(ir,ie)
         do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              vol = volume(kt)
              h_rod = rod_node(ir,n1,ie)
              h_nrod = 1. - h_rod
              n1_rod = n1 + 1
              n1_nrod = n1 - 1

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node

           call XS_Compute_Xappa_Rod(h_rod,e_rod_round_off, &
                 n1, n1_rod, n1_nrod, k, ir, ie, &
                 Xappa_Rod)


             ARGUM(1)=fdback(k,n1,3)*1000.
              ARGUM(2)=fdback(k,n1,4)
           ARGUM(3)=fdback(k,n1,2)+ Convert_Celcius_to_Kelvin  !   T_c
           ARGUM(4)=0.
           ARGUM(5)=brn(k,n1) ! burnup
           ARGUM(6)=fdback(k,n1,1)              

!           write(*,*) 'xappa_rod =', xappa_rod

           IF( nl .GT. 0 ) THEN
              CALL YMNCCALC(ARGUM, nl, CONST, DCONST, KSDCA)
           ELSE 
              CALL YNREFsection ( ARGUM, nl, CONST)
           END IF 

! ALL assembly without control rods

!       difk(k,1)=const(1)
!       difk(k,2)=const(9)
!       signf(k,1)=const(3)
!       signf(k,2)=const(11)
!       sigp(k,1)=const(4)
!       sigp(k,2)=const(12)
!       sigr(k,1)=const(2)
!       sigr(k,2)=const(10)
!       sig12(k)=const(5)


!           do n = 1, NG 
! Diffusion Coefficients
!              if(I_Diff_Coeff.ne.1) then
! differential cross sections are given for sigma transport 
!              sect_nrod = 1./(3.*XS_D(n,kt))
!              sect_rod  = sect_nrod + d_ca(n,nl)
!              str = XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(n))
!              XS_D(n, kt) = 1./(3.*str)
!              ELSE
! differential cross sections are given for diffusion coefficient
              IF( i_Diff_Coeff == 1) THEN
                 sect_nrod = XS_D(1,kt)
                 sect_rod  = const(1)
                 XS_D(1, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(1))

                 sect_nrod = XS_D(2,kt)
                 sect_rod  = const(9)
                 XS_D(2, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
              ELSE
                 sect_nrod = XS_D(1,kt)
                 sect_rod  = 1./const(1)
                 XS_D(1, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(1))

                 sect_nrod = XS_D(2,kt)
                 sect_rod  = 1./const(9)
                 XS_D(2, kt) =  &
                   XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
              END IF
!              END IF
! all other XS are multiplied by the node volume 
! scattering
              sect_nrod = XS_SIK(2, 1, kt)
              sect_rod =  const(5)*vol
               XS_SIK(2, 1, kt)=XS_Sect_CR(sect_nrod, sect_rod,&
              xappa_rod(1))
!             end do ! scattering
!absorption 
              sect_nrod = XS_SA(1, kt)
              sect_rod = const(2)*vol
              XS_SA(1, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(1))

              sect_nrod = XS_SA(2, kt)
              sect_rod = const(10)*vol
              XS_SA(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(2))

              if(n1.le.NZ_Core_END) then
! no need for nuSF of the UPPER REFLECTOR (PWR NEACRP PROBLEM)

              sect_nrod = XS_SF(1, kt)
              sect_rod = const(3)*vol
              XS_SF(1, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(1))

              sect_nrod = XS_SF(2, kt)
              sect_rod = const(11)*vol
              XS_SF(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(2))

! Eliminating the round of the the nonfuel absorbers
!              IF(XS_SF(n,kt).LT.0) XS_SF(n,kt)=0.

              sect_nrod = XS_SF_P(1, kt)
              sect_rod = const(4)*vol
              XS_SF_P(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(1))

              sect_nrod = XS_SF_P(2, kt)
              sect_rod = const(12)*vol
              XS_SF_P(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(2))

! i-135
            i_isotope = 1
            sect_nrod =  yields_isotope(i_isotope,  kt)
            sect_rod  = const(15)
            yields_isotope(i_isotope,  kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))

! xe-135
            i_isotope = 2
            sect_nrod =  yields_isotope(i_isotope,  kt)
            sect_rod  = const(16)
            yields_isotope(i_isotope,  kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
            sect_nrod = sa_isotope(i_isotope, NG, kt)
            sect_rod  = const(7)
            sa_isotope(i_isotope, NG, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))

! pm-149
            i_isotope = 3
            sect_nrod =  yields_isotope(i_isotope,  kt)
            sect_rod  = const(17)
            yields_isotope(i_isotope,  kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
! sm-149
            i_isotope = 4
            sect_nrod = sa_isotope(i_isotope, NG, kt)
            sect_rod  = const(8)
            sa_isotope(i_isotope, NG, kt)=&
              XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))


!              write(*,*) 'CR type =', nl
!         rinf(k) =  s_12(k)/sr2(k)
!          KINF(k) =  (nuSF1(k) + rinf(k)*nuSF2(k))/sr1(k)

!            rinf =  const(5)/const(10)
!              KINF =  (const(3) + rinf*const(11))/(const(2)+const(5))
!              write(*,*) 'KINF =', KINF
!              rinf = XS_SIK(2, 1, kt)/ XS_SA(2, kt)
!              kinf = (XS_SF(1, kt) + rinf*XS_SF(2, kt))/&
!               ( XS_SA(1, kt) +   XS_SIK(2, 1, kt) )
!              write(*,*) 'KINF in XS =', KINF
!      write(*,*) 's12 (lib, XS)=',  const(5), XS_SIK(2, 1, kt)/vol
!      write(*,*) 'sa1          =',  const(2), XS_SA(1, kt)/vol
!      write(*,*) 'sa2          =',  const(10), XS_SA(2, kt)/vol
!      write(*,*) 'nuSF1        =',  const(3), XS_SF(1, kt)/vol
!      write(*,*) 'nuSF2        =',  const(11), XS_SF(2, kt)/vol
!      read(*,*) 
! Eliminating the round of the the nonfuel absorbers
!              IF(XS_SF_P(n,kt).LT.0) XS_SF_P(n,kt)=0.

             end if ! nuSF in the core 

!           end do ! NG

             end if ! if there is anything in the node
           end do ! NZ
          end if ! k.ie.0
         end do ! NCHM

         end do
        end do ! NN_CRod
       end do ! NN_CRod_El

!       if(Debug) CALL XS_DBG_output

      return
      end

