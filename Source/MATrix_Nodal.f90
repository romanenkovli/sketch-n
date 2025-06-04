      subroutine MAT_Solve_Two_Node(k_eff_2)
!=====================================================================*
! Solution of the Two-Node Problems and Computing New Coupling        *
!     Coefficients (Polynomial & Semi_Analytical Nodal Method)        *
! Last Update              Slava (c) April 1998                       *
!=====================================================================*
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'

! Input: trl_0_xyz, trl_qla, all XS_, Flux, Matr_FD, Diag_Tot
      real k_eff_2
! Output: Matr_Tot(NG, NE_T, N_TOT), Diag_Tot(NG, N_TOT)

! Local Variables:
      integer k, n, kt, nd, nd_right, m, k_left, k_right, NN, nl

      real x1(NG), x2(NG), x3(NG),x4(NG) ! Neutron Flux Expansion Coeff
      real b8(2*NG), b8_l(2*NG) ! RHS of the 1st Expansion Coeffi
      real sa2(NG), sik2(NG,NG), sf2(NG), d2(NG), sf0(NG,NG) ! XS
      real adf_k(NG), adf_k_plus_1(NG)
      real a0(NG), a01(NG)
      real a_rl 
      real am2(NG,NG), am4(2*NG,2*NG), am4_l(2*NG,NG)
      real fm0(NG)
      real c3j(NG), c4j(NG), c2j

      REAL current(NG), flux_surface(NG)
      real a_vol, hk, sq,  d2_left(NG), hk_left,&
          flux_left(NG)  !, MAT_Nod

      real Coeff_2(NG), xappa(NG,NG), trans0(NG), transq(2*NG),&
          Coeff_42(NG), Coeff_31(NG), Coeff_1c(NG), Coeff_22

!      real Flux_Rodded(NG), Flux_Unrodded(NG)
      
!OMEGA      real xmm, xp_omega, xp_prompt(NG)

! Initial data for Polynomial Nodal Method
      real coeff_diag
      data coeff_diag /0.0/  ! diagonal coefficient for 3rd moment
!      = 0.0 in the both HEX and XYZ geometry
! Local data for conformal mapping of the hexagon into rectangle
      real const_hk ! a/H (x-side of rectangle / pitch of hex)
      parameter (const_hk = 1.0525557)
!      parameter (const_hk = 1.0) 
      real const_ju ! b/R (y-side of rectangle / side of hex)
!cc   real const_sq
!cc      parameter (const_sq = 0.822782)
      logical flag_hex

      real cg_nod(2, NDIR), hk_fd, fg_nod(NG, 2, NDIR)  ! 
!      modified constant for boundary conditions 
! cg_nod is modified for HEX geometry due to conformal mapping
! fg_nod is modified due to Assembly Discontinuity Factors 

! ANM Method
      real Matr_Funct(NG,NG,5)
      real am2_k(NG,NG), RHS_Flux_K(NG)
      Logical ANM
      integer i, k_out


      nl = 0

      do nd = 1, NDD

       nd_right = Index_Neib(nd,2)

       IF( (GMT_CRD_TYPE(1:4) .EQ. "HEXZ") .AND. (nd.NE.NDIR) ) THEN
            c2j = 2.58834 ! PNM 3.
            const_ju = 1.42510
            Coeff_22 =  3.03919  ! 0.0 in PNM
            DO n = 1, NG
               c3j(n) = 5.29438 ! PNM 6.
               c4j(n) = 8.90412 ! PNM 10.
               Coeff_42(n) = 0.0186293 ! PNM 0.028571429/ ! 1./35.
               Coeff_31(n) = 0.0529302 ! PNM 0.066666667/ ! 1./15.
            END DO
            Flag_hex = .True.
       ELSE ! XYZ geometry or Z ditrection in HEXZ
            c2j =  3.
            const_ju = 1.
            Coeff_22 = 0.0 
            DO n = 1, NG
               c3j(n) = 6.
               c4j(n) = 10.
               Coeff_42(n) = 0.028571429
               Coeff_31(n) = 0.066666667
            END DO
            Flag_hex = .False.
       END IF


       DO n = 1, NG
         Coeff_2(n)  = Coeff_42(n)*c4j(n)  ! 10/35=2/7
         Coeff_1c(n) = Coeff_31(n)*c3j(n)  ! 6/15=2/5
       END DO

! modifying the logarithmic derivative boundary condition for HEX
        cg_nod(1, nd) = cg(1, nd)*const_ju
        cg_nod(2, nd) = cg(2, nd)*const_ju

!        write(*,*) 'cg_nod(1:2, nd)' 
!        write(*,*) cg_nod(1:2, nd)
!        pause


       do kt = 1, N_TOT_DIR(nd)

          nl = nl + 1 ! Left Interface

          k = Numb(kt, nd)
          a_vol = 1./volume(k)
          sq = s_xyz(nd, k) 
          hk  = h_xyz(nd, k)
          hk_fd = hk

          IF ( Flag_Hex ) THEN
              hk = hk*const_hk
          END IF



          do n = 1, NG

             fm0(n) = Flux(n,k)

             sa2(n) = XS_SA(n,k)*a_vol
             sf2(n) = XS_SF(n,k)*a_vol
             d2(n) = XS_D(n,k)
             a0(n) = 2.*d2(n)/hk
             a01(n) = (hk*hk)/(4.*d2(n))

             do m = 1, NG
                sik2(n,m) = XS_SIK(n,m,k)*a_vol
             end do
          end do

             do n = 1, NG
                do m = 1, NG
                   sf0(n,m) = xp(n)*sf2(m)/k_eff_2 + sik2(n,m)
                end do
             end do

          do n = 1, NG
             do m = 1, NG
                   xappa(n,m) = - sf0(n,m)*a01(n)
             end do
                   xappa(n,n) = xappa(n,n) + sa2(n)* a01(n)
          end do

!          write(*,*) 'k_ef=', k_eff_2
!          IF ( k == 4 ) THEN
!        write(*,*) 'k = ', k 
!          DO n = 1, NG
!            write(*,*) 'xappa=', (xappa(n,m),m=1,NG)
!          END DO
!          DO n = 1, NG
!            write(*,*) 'sf0=', (sf0(n,m),m=1,NG)
!          END DO

!         write(*,*) 'fm0=', fm0(1:NG) 
!         write(*,*) 'flux_surface=', flux_surface(1:NG) 

!         pause
!          END IF 

! Cartesian geometry SANM
         if(Nodal_Method.EQ."SANM") then 
! Semi-Analytical Method

           call MAT_Set_SANM(NG,hk,d2,sa2, c3j, c4j, &
            Coeff_2, Coeff_42, Coeff_31, Coeff_1c, c2j)

         else if(Nodal_Method.EQ."ANM") then
! Analytic Nodal Method
            if(NG.eq.2) then
               call MAT_Set_MatFunc_ANM(xappa, Matr_Funct)
            else
               ANM = .True.
!               write(*,*) 'Multigroup ANM requires IMSL libraries'
!               write(*,*) 'Uncomment the line below'
!               write(*,*) 'Include file "ANM_Ng.for" into Makefile'
!               write(*,*) 'Link with ISML library'
              call MAT_Set_MatFunc_NG&
                                (xappa, Matr_Funct, NG, ANM)
            end if
! Polynomial Nodal method Based on Matrix Functions
         else if(Nodal_Method.EQ."PNM1") then
            if(NG.eq.2) then
               call MAT_Set_MatFunc_PNM(xappa, Matr_Funct)
            else
               ANM = .False.
!               write(*,*) 'Multigroup Matrix Functions PNM &
!                       requires IMSL libraries'
!               write(*,*) 'Uncomment the line below'
!               write(*,*) 'Include file "ANM_NG.for" into Makefile'
!               write(*,*) 'Link with ISML libraries'
               call MAT_Set_MatFunc_NG&
                                 (xappa, Matr_Funct, NG, ANM)
            end if
         end if ! Nodal Method

         k_left = Neib(nd,k)
         k_right = Neib(nd_right, k)
         IF (k_left.ne.I_BOUND_NODE) THEN
            hk_left = h_xyz(nd, k_left)
            DO n = 1, NG
               d2_left(n) = xs_d(n,k_left)
               flux_left(n) = flux(n, k_left)
            END DO
         END IF               

            DO n = 1, NG
             adf_k_plus_1(n) = xs_adf(n,1,nd,k)
             IF(k_left.ne.I_BOUND_NODE) adf_k(n) = &
                  xs_adf(n,2,nd,k_left)
            END DO 

          do n = 1, NG
             trans0(n)    = TRL_0_xyz(n, k, nd)*a01(n)
             transq(n)    = trl_qla(1, n, k, nd)*a01(n)
             transq(n+NG) = trl_qla(2, n, k, nd)*a01(n)
          end do

         IF(Nodal_Method.EQ."SANM".OR.Nodal_Method.EQ."PNM") then
! 2nd expansion coefficient
           call MAT_Compute_Mat2(NG, &
                fm0, transq(NG+1), xappa, Coeff_2, trans0, x2, am2, &
                c2j, coeff_22 )
            NN = NG
            if(NG.eq.2) then
               call MSC_Solve2x2(am2, x2)
            else
               call MSC_LU_Solve(am2, NN, x2)
            end if
! 4th expansion coefficient
            call MAT_Compute_3_4_Moments(&
                NG, x2, xappa, transq(NG+1), Coeff_42, coeff_22, x4)

!        IF (  nd ==  1) THEN
!        IF ( k == 4 ) THEN
!      write(*,'(4A14)')  'trans0'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  trans0(n)/fm0(n)
!        END DO 
!        END IF 
!        END IF  


         END IF

         if(k_left.eq.I_BOUND_NODE.AND. i_dr(1,nd) /= 2) then 
! Left Boundary Node ( not computed in the case of the REPEATED STRUCTURE
! Modifying the boundary conditions due to ADF
           DO n = 1, NG
              fg_nod(n,1,nd) = fg(n,1,nd)*adf_k_plus_1(n)
           END DO                
                
           a_rl = -1.         


         IF(Nodal_Method.EQ."SANM".OR.Nodal_Method.EQ."PNM") then
           call MAT_Compute_Mat1_Bound(NG, xappa, a0, &
                   c4j, Coeff_1c, Coeff_31, transq(1), fm0, x2, x4,&
                   a_rl, cg_nod(1, nd), fg_nod(1,1,nd), am4_l, am2, &
                   x1, b8_l, c2j)
           if(NG.eq.2) then
             call MSC_Solve2x2(am2, x1)
           else
             call MSC_LU_Solve(am2,NG,x1)
           end if
           call MAT_Compute_3_4_Moments&
             (NG, x1, xappa, transq(1), Coeff_31, coeff_diag, x3)
           CALL MAT_Compute_Current(NG, a0, x1, x2, x3, x4, c2j, &
                 c3j, c4j, a_rl, const_ju, current)
           IF( NonlinearIterations.EQ."Moon") THEN
           CALL MAT_Compute_Flux_Surface(NG, fm0, x1, x2, x3, x4, &
                  a_rl, flux_surface)
           END IF
          ELSE
           call MAT_Set_MatrFlux_Bound(NG, a0, fm0, trans0, &
                 transq(1), transq(1+NG), am2, &
                 am2_k, flux_surface, RHS_Flux_K, Matr_Funct, &
                 a_rl, fg_nod(1,1,nd), cg(1,nd) )
           if(NG.eq.2) then
             call MSC_Solve2x2(am2, flux_surface)
           else
             call MSC_LU_Solve(am2, NG, flux_surface)
           end if
           call MAT_Compute_Current_ANM(&
                 NG, a_rl, a0, fm0, flux_surface, trans0, transq(1),&
                 transq(NG+1), Matr_Funct, current)
          END IF


              IF( NonlinearIterations.EQ."Smith") THEN
                CALL MAT_Compute_Coupling_Coeff_Bound(NG, &
                 mat_fd(1,nl), current, fm0, sq, a_rl, mat_nod(1,nl))
              ELSE IF( NonlinearIterations.EQ."Moon") THEN
! we are working in cartesian coordinates
                DO n = 1, NG 
                   flux_surface(n) = flux_surface(n)*adf_k_plus_1(n)
                END DO                     

!                IF(nd.eq.1.AND.k.eq.1) THEN 
!                  WRITE(*,*) 'k=1'
!                  WRITE(*,*) 'Surface Flux =', Flux_Surface
!                END IF

                CALL MAT_Compute_Coupling_Coeff_Bound_Moon(NG, current,&
                  fm0, hk_fd, d2,  sq, flux_surface, a_rl, fg(1,1,nd), &
                     cg(1,nd) , mat_fd(1,nl), mat_nod(1,nl), &
                     d_nod(1,k,1,nd) )

              END IF


                                

!              MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_Nod(n,nl)

!              MAT_Tot(n, nd, k) = MAT_FD(n, nl) - MAT_Nod(n,nl)

         else if(k_left.eq.I_BOUND_NODE.AND. i_dr(1,nd) == 2) THEN
! in the case of the REPEATED STRUCTURE we computec only matrices
         a_rl = -1.         
         IF( i_dr(1,nd) == 2 )    THEN
           IF(Nodal_Method.EQ."SANM".OR.Nodal_Method.EQ."PNM") then
              call MAT_Compute_Mat1_Bound(NG, xappa, a0, &
                   c4j, Coeff_1c, Coeff_31, transq(1), fm0, x2, x4,&
                   a_rl, cg_nod(1, nd), fg_nod(1,1,nd), am4_l, am2, &
                   x1, b8_l, c2j)
           ELSE ! ANM
           call MAT_Set_MatrFlux_Bound(NG, a0, fm0, trans0, &
                 transq(1), transq(1+NG), am2, &
                 am2_k, flux_surface, RHS_Flux_K, Matr_Funct, &
                 a_rl, fg_nod(1,1,nd), cg(1,nd) )
           END IF ! NODAL METHODS
         END IF        

      ELSE 
! the normal case of the two-node problem
!         k_out = 2
!         IF( k == k_out.AND. nd==1)  THEN
!           write(*,*) 'k=', k_out, ' nd=1, AM4_L(4,2), B8_L(4)'
!          DO n = 1, 2*NG
!           write(*,'(4E14.5,5x,E14.5)') (AM4_L(n,i),i=1,NG), B8_L(n)
!            END DO   
!           write(*,*) 'trans0(NG)'
!           write(*,'(4E14.5)') trans0(1:NG)
!             write(*,*) 'transq(2*NG)'
!             write(*,'(4E14.5)') transq(1:2*NG)
!          pause
!         END IF
          
         IF(Nodal_Method.EQ."SANM".OR.Nodal_Method.EQ."PNM") then
          call MAT_Compute_Mat1&
                  (NG, xappa, a0, c4j, Coeff_1c, Coeff_31,&
                   transq(1), fm0, x2, x4, am4_l, am4, b8, b8_l, &
                   c2j, adf_k, adf_k_plus_1)

           NN = NG*2
           call MSC_LU_Solve(am4,NN,B8)
           do n = 1, NG
              x1(n) = b8(NG + n)
           end do
! calculation of the 3rd moment
         call MAT_Compute_3_4_Moments(NG, x1, xappa,&
                          transq(1),  Coeff_31, coeff_diag, x3)
         CALL MAT_Compute_Current(NG, a0, x1, x2, x3, x4, c2j, &
                 c3j, c4j, a_rl, const_ju, current)
         IF( NonlinearIterations.EQ."Moon") THEN
            CALL MAT_Compute_Flux_Surface(NG, fm0, x1, x2, x3, x4, &
                 a_rl, flux_surface)
         END IF

!        IF( k == k_out.AND. nd==1)  THEN
!           write(*,*) 'k=', k_out, ' nd=1, CURRENT'
!           write(*,'(4E14.5)') current(1:NG)
!           pause
!        END IF 
 
        ELSE ! ANM 
         call MAT_Set_MatrFlux(NG, a0, fm0, trans0(1), &
                transq(1),  transq(1+NG), am2, &
                am2_k, flux_surface, RHS_Flux_K, Matr_Funct,&
                adf_k, adf_k_plus_1)
           if(NG.eq.2) then           
               call MSC_Solve2x2(am2, flux_surface)
           else
               NN = NG
               call MSC_LU_Solve(am2, NN, flux_surface)
           end if
         call MAT_Compute_Current_ANM(NG, a_rl, a0, fm0,&
                flux_surface, trans0(1),  transq(1), transq(NG+1), &
                Matr_Funct, current)
!         write(*,*) 'k=', k, 'a_rl=', a_rl  
!         write(*,*) 'a0 =', a0(1:NG)
!         write(*,*) 'fm0=', fm0(1:NG)
!         write(*,*) 'flux_surface=', flux_surface(1:NG)
!         write(*,*) 'current=', current(1:NG)
!           write(*,*) 'Matrices(1)='
!           do n = 1, NG
!             write(*,*) (Matr_Funct(n,m,1),m=1,NG)
!           end do
!           write(*,*) 'Matrices(2)='
!           do n = 1, NG
!             write(*,*) (Matr_Funct(n,m,2),m=1,NG)
!           end do
!           write(*,*) 'trans0=', trans0(1:NG)
!           write(*,*) 'transq=', transq(1:2*NG)
!         pause
        END IF   

! calculation of the  nodal coupling coefficients
              IF( NonlinearIterations.EQ."Smith") THEN
                CALL MAT_Compute_Coupling_Coeff(NG, mat_fd(1,nl), &
                  current, fm0, flux(1,k_left), sq, mat_nod(1,nl))
              ELSE IF( NonlinearIterations.EQ."Moon") THEN
                DO n = 1, NG 
                   flux_surface(n) = flux_surface(n)*adf_k_plus_1(n)
                END DO                     


               CALL MAT_Compute_Coupling_Coeff_Moon(NG, current,&
                     fm0, flux_left, hk_fd, hk_left, d2, d2_left, sq, &
                     flux_surface, mat_fd(1,nl), mat_nod(1,nl), &
                     d_nod(1,k_left,2,nd), d_nod(1,k,1,nd) )

!        IF (  nd ==  1) THEN
!        IF ( k == 4 ) THEN
!      write(*,'(4A20)')  'left interface'

!      write(*,'(4A14)')  'x1(n)', 'x2(n)', 'x3(n)', 'x4(n)'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  x1(n)/fm0(n), x2(n)/fm0(n), &
!                    x3(n)/fm0(n), x4(n)/fm0(n)
!        END DO 
!      write(*,'(4A14)')  'fm0(n)', 'flux_surface', 'flux_surface',&
!                    'current'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  fm0(n), flux_surface(n)/fm0(n),&
!           (1.+(-x1(n)+x2(n)-&
!                    x3(n)+x4(n))/fm0(n)), current(n)/fm0(n)
!        END DO 
!        END IF 
!        END IF  


!                IF(nd.eq.1.AND.k.eq.2) THEN 
!                  WRITE(*,*) 'k=2'
!                  WRITE(*,*) 'Surface Flux =', Flux_Surface
!                  pause
!                END IF

              END IF

!               MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_Nod(n,nl)

!               MAT_Tot(n, nd, k) = MAT_FD(n, nl) - MAT_Nod(n,nl)

!               MAT_Diag_Tot(n, k_left) = MAT_Diag_Tot(n, k_left) - &
!                                           MAT_Nod(n,nl)

!               MAT_Tot(n, nd_right, k_left) = MAT_FD(n, nl) + &
!                    MAT_Nod(n,nl)

                 END IF


         if(k_right.eq.I_BOUND_NODE) then        
! Right Boundary Node
           nl = nl + 1 ! Right Interface

         IF( i_dr(2,nd) /= 2 ) THEN 

           a_rl = 1.         

        DO n = 1, NG
             adf_k_plus_1(n) = xs_adf(n,2,nd,k)
             fg_nod(n,2,nd) = adf_k_plus_1(n)*fg(n,2,nd)
          END DO

         IF(Nodal_Method.EQ."SANM".OR.Nodal_Method.EQ."PNM") then

          call MAT_Compute_Mat1_Bound(NG, xappa, a0, &
             c4j, Coeff_1c, Coeff_31, transq(1), fm0, x2, x4, a_rl, &
             cg_nod(2,nd), fg_nod(1,2,nd), am4_l, am2, x1, b8_l, c2j)


! LU decomposition and solution
          NN = NG
          if(NG.eq.2) then
             call MSC_Solve2x2(am2, x1)
          else
             call MSC_LU_Solve(am2,NN, x1)
          end if
! 3rd moment, current and diffusion coefficient 
         call MAT_Compute_3_4_Moments&
            (NG, x1, xappa, transq(1), Coeff_31, coeff_diag, x3)
         CALL MAT_Compute_Current(NG, a0, x1, x2, x3, x4, c2j, &
                 c3j, c4j, a_rl, const_ju, current)
         IF( NonlinearIterations.EQ."Moon") THEN
            CALL MAT_Compute_Flux_Surface(NG, fm0, x1, x2, x3, x4, &
                 a_rl, flux_surface)
         END IF

!        IF (  nd ==  1) THEN
!        IF ( k == 4 ) THEN

!      write(*,'(4A20)')  'right interface'

!      write(*,'(4A14)')  'x1(n)', 'x2(n)', 'x3(n)', 'x4(n)'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  x1(n)/fm0(n), x2(n)/fm0(n), &
!                    x3(n)/fm0(n), x4(n)/fm0(n)
!        END DO 
!      write(*,'(4A14)')  'fm0(n)', 'flux_surface', 'flux_surface',&
!                    'current'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  fm0(n), flux_surface(n)/fm0(n),&
!           (1.+(x1(n)+x2(n)+&
!                    x3(n)+x4(n))/fm0(n)), current(n)/fm0(n)
!        END DO 
!        END IF 
!        END IF  

        ELSE
           call MAT_Set_MatrFlux_Bound(NG, a0, fm0, trans0(1),&
                 transq(1), transq(1+NG), am2, &
                 am2_k, flux_surface, RHS_Flux_K, Matr_Funct, &
                 a_rl, fg_nod(1,2,nd), cg(2,nd) )
           if(NG.eq.2) then           
               call MSC_Solve2x2(am2, flux_surface)
           else
               NN = NG
               call MSC_LU_Solve(am2, NN, flux_surface)
           end if
           call MAT_Compute_Current_ANM(NG, a_rl, a0, fm0, &
                flux_surface, trans0(1), transq(1), transq(NG+1), &
                Matr_Funct, current)
        END IF


! calculation of the nodal coupling coefficients
             IF( NonlinearIterations.EQ."Smith") THEN
                CALL MAT_Compute_Coupling_Coeff_Bound(NG, &
                 mat_fd(1,nl), current, fm0, sq, a_rl, mat_nod(1,nl))
             ELSE IF( NonlinearIterations.EQ."Moon") THEN

              DO n = 1, NG 
                   flux_surface(n) = flux_surface(n)*adf_k_plus_1(n)
              END DO                     

              CALL MAT_Compute_Coupling_Coeff_Bound_Moon(NG, current,&
                fm0, hk_fd, d2,  sq, flux_surface, a_rl, fg(1,2,nd), &
                     cg(2,nd) , mat_fd(1,nl), mat_nod(1,nl), &
                     d_nod(1,k,2,nd) )



             END IF 
             
!             MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) - MAT_Nod(n,nl)
                     
!             MAT_Tot(n,nd_right,k) = MAT_FD(n, nl) + MAT_Nod(n,nl)


!         write(*,*) 'Neutron flux on the right'
!         DO n = 1, NG
!           write(*,*) 'n=', n, &
!                 (fm0(n) + x1(n) + x2(n) + x3(n) + x4(n))/fm0(n)
!         END DO
         
         END IF ! ( i_dr(2,nd /= 2) 
         end if ! k_right = 0                 

       end do !k
      end do ! NDD


!      if(Debug) then
!         call MAT_DBG_Check_Matr_tot
!   call Eigenv_Write
!   pause 'Matrix Written'
!      end if

       


      return
      end

      subroutine  MAT_Compute_3_4_Moments&
            (NG, x2, xappa, transq, Coeff_Mom, coeff_diag, x4)
!=====================================================================*
!  Computing 3rd and 4th Expansion Coefficients (Moments)             *
! Last Update              Slava (c) April 1998                        *
!=====================================================================*

      implicit none

! Input: x2(NG),sa2(NG),sf0(NG,NG),trlq(NG),a01(NG),x4(NG)
      integer NG
      real x2(NG), xappa(NG,NG), transq(NG), Coeff_Mom(NG)
      real coeff_diag
! Output: x4(NG) - 3rd or 4th Expansion Coefficient 
      real x4(NG) 
! Local Variables: 
      integer n, m

      do n = 1, NG

          x4(n) = transq(n)  

          do m = 1, NG
               x4(n) = x4(n) + xappa(n,m)*x2(m)
          end do
! additional term for hexagonal geometry              
          x4(n) = x4(n) - coeff_diag*x2(n)

          x4(n) = Coeff_Mom(n)*x4(n)

      end do


      return
      end

      subroutine MAT_Compute_Mat2(&
       NG, fm0, transq, xappa, Coeff_2, trans0,  x2, am2, c2j, &
       coeff_22 )
!=====================================================================*
! Compute Matrix & RHS of the 2nd Expansion Coefficient Equations     * 
! Last Update 7.V.1998        Slava (c) April 1998                    *
!=====================================================================*
      implicit none
! Input:
      integer NG
      real fm0(NG), transq(NG), xappa(NG,NG), Coeff_2(NG), trans0(NG)
      real c2j, coeff_22
! Output: am2 - Matrix of the 2nd Moment Equations
!         x2 - Right Side
      real am2(NG,NG), x2(NG)
! Local Variables: 
      integer n, m

      do n = 1, NG

         x2(n) = trans0(n) - Coeff_2(n)*transq(n)

         do m = 1, NG
             am2(n,m) = Coeff_2(n)*xappa(n,m)
             x2(n) =  x2(n) + xappa(n,m)*fm0(m)
         end do

         am2(n,n) = am2(n,n) + c2j - Coeff_2(n)*Coeff_22

      end do

      return
      end


      subroutine MAT_Compute_Mat1_Bound(NG, xappa, a0, Coeff_14,&
             Coeff_1c, Coeff_1f, &
             transq, fm0, x2, x4, a_rl, cg, fg, am4_l, am2, x1, b8_l,&
             c2j )
!=====================================================================*
! Compute Matrix & RHS of the 1st Expansion Coefficient Equations     * 
! Last Update Slava (c) April 1998 (PNM & SANM Method)                *
!               !!!! BOUNDARY NODES  !!!!                             *
!=====================================================================*
      implicit none
! Input: 
      integer NG
      real xappa(NG,NG), a0(NG), Coeff_14(NG), Coeff_1c(NG), &
                 Coeff_1f(NG), x2(NG), x4(NG), fm0(NG), transq(NG),&
                 a_rl, cg, fg(NG), c2j
! Output: AM2(NG,NG) - Matrix of the 1st Moment Equations 
!         x1(NG) - Right Side of the First Moment Equations
!         am4_l(2*NG,NG) - Matrix for the Next Two-Node Problem
!         b8_l(2*NG) - Right Side for the Next Two-Node Problem
      real  x1(NG), am2(NG,NG), b8_l(2*NG), am4_l(2*NG,NG)
! Local Variables: 
      integer n, m, ii
      real b8_e, b8_o, Coeff_Current, Coeff_Matr, Diag_El

!      real am2_check(NG,NG), x1_check(NG)


!     write(*,*) 'Coeff_14 =', Coeff_14
!     write(*,*) 'Coeff_1c =', Coeff_1c
!     pause

      do n = 1, NG
          
          Coeff_Current = a0(n)*Coeff_1c(n)

! current continuity
          b8_e =  a0(n) * (c2j*x2(n)+ Coeff_14(n)*x4(n) )
          b8_o = Coeff_Current*transq(n)
          b8_l(n) = b8_o + b8_e

! * flux continuity
          ii = n + NG

          b8_e = fm0(n) + x2(n) + x4(n) 
          b8_o = Coeff_1f(n)*transq(n)
          b8_l(ii) = b8_e + b8_o

          Coeff_Matr = a_rl*(fg(n)*Coeff_1f(n) + &
                                   cg*Coeff_Current)
          Diag_El = a_rl * (fg(n) + cg*a0(n))

          x1(n) = - fg(n)*(fm0(n) + x2(n) + x4(n)) - &
                        cg*a0(n)*(c2j*x2(n) + Coeff_14(n)*x4(n)) -&
                        Coeff_Matr*transq(n)

          do m = 1, NG
! right side of the matrix for the NEXT Two-Node problem
! current continuity
                am4_l(n,m) = Coeff_Current*xappa(n,m)
! flux continuity          
                am4_l(ii,m) = Coeff_1f(n)*xappa(n,m)
 
                am2(n,m) = Coeff_Matr*xappa(n,m)
          end do
! right half of the matrix for two-node problem (Diagonal Terms)
! current continuity

          am4_l(n,n) = am4_l(n,n) + a0(n)
! flux continuity          
          am4_l(ii,n) = am4_l(ii,n) + 1.

          am2(n,n) = am2(n,n) + Diag_El

      end do

!      write(*,*) 'am2_check =', am2_check
!      write(*,*) 'x1_check =', x1_check
!      write(*,*) 'b8_l =', b8_l

!      pause

      return
      end


      subroutine MAT_Compute_Mat1&
                  (NG, xappa, a0, Coeff_14, Coeff_1c, Coeff_1f,&
                     transq, fm0, x2, x4, am4_l, am4, b8, b8_l,&
                     c2j, adf_k, adf_k_plus_1 )
!=====================================================================*
! Compute Matrix & RHS of the1st Expansion Coefficient Equations      * 
! Last Update Slava (c) April 1998 (PNM & SANM)                       *
!               !!!! Normal Two-Node Problem  !!!!                    *
!=====================================================================*
      implicit none
! Input: 
      integer NG
      real xappa(NG,NG), a0(NG), Coeff_14(NG), Coeff_1c(NG), &
                 Coeff_1f(NG), x2(NG), x4(NG), fm0(NG), transq(NG),&
          c2j
! adf_k - adf at the RIGHT face of the node k 
! adf_k_plus_1 - adf at the LEFT face of the node k+1 
      real adf_k(NG), adf_k_plus_1(NG)  
! Output: 
!  AM4(2*NG,2*NG) - Matrix of the 1st Moment Equations
!         b8(2*NG) - Right Side of the First Moment Equations
! RHS(NG) - Right Hand Side for the 1st expansion coefficient node k+1
! AM2(NG,NG) - Matrix for the 1st expansion coefficient node k+1
!         am4_l(2*NG,NG) - Matrix for the Next Two-Node Problem
!         b8_l(2*NG) - Right Side for the Next Two-Node Problem
!     
      real  b8(2*NG),b8_l(2*NG), am4_l(2*NG,NG), AM4(2*NG,2*NG)
! Local Variables: 
      integer n, m, ii, jj
      real b8_e, b8_o, Coeff_Current

      do n = 1, NG
          
          Coeff_Current = a0(n)*Coeff_1c(n)

! current continuity
          b8_e =  a0(n) * (c2j*x2(n)+ Coeff_14(n)*x4(n) )
          b8_o = Coeff_Current*transq(n)
          b8(n) = b8_e - b8_o + b8_l(n)
          b8_l(n) = b8_o + b8_e

! * flux continuity
          ii = n + NG

          b8_e = fm0(n) + x2(n) + x4(n) 
          b8_o = Coeff_1f(n)*transq(n)
! incorporated adf
          b8(ii) = (b8_e - b8_o)*adf_k_plus_1(n) - b8_l(ii)*adf_k(n)
          b8_l(ii) = b8_e + b8_o


          do m = 1, NG
! left half of the matrix
! current 
                am4(n,m) = - am4_l(n,m)
! flux
                am4(ii,m) = am4_l(ii,m)
! right side of the matrix
! current continuity
                jj = NG + m
                am4(n,jj) = Coeff_Current*xappa(n,m)
                am4_l(n,m) = am4(n,jj)
! flux continuity          
                am4(ii,jj)= Coeff_1f(n)*xappa(n,m)
                am4_l(ii,m) = am4(ii,jj)
          end do
! right half of the matrix for two-node problem (Diagonal Terms)
! current continuity

          am4(n,ii) = am4(n,ii) + a0(n)
          am4_l(n,n) = am4(n,ii)
! flux continuity          
          am4(ii,ii) = am4(ii,ii) + 1.
          am4_l(ii,n) = am4(ii,ii)

! incorporated adf into the flux continuity equations
         DO m = 1, NG
! left half of the matrix
            am4(ii,m) = am4(ii,m)*adf_k(m)
            jj = NG + m
! right half of the matrix
            am4(ii,jj)=am4(ii,jj)*adf_k_plus_1(m) 
         END DO 

      end do
! 

      return
      end



      subroutine MAT_Set_SANM(NG, hk, d2, sa2, c3j, c4j, &
            Coeff_2, Coeff_42, Coeff_31, Coeff_1c, c2j)
!=====================================================================*
! Compute Hyperbolic Functions and their Moments for the SANM         * 
! Last Update Slava (c) April 1998 (PNM & SANM)                       *
!               !!!! Normal Two-Node Problem  !!!!                    *
!=====================================================================*
      implicit none
! Input: 
      integer NG
      real hk, d2(NG), sa2(NG) 
! Output: 
!      c3j(NG), c4j(NG) - coefficients to calculate neutron current
!     Coeff_2(NG) - Coefficients used for the Neutron Balance Equation
!     Coeff_42(NG) - Coefficient to compute the 4th moment from the 2nd
!     Coeff_31(NG) - Coefficient to compute the 3rd moment from the 1st

      real  c2j, c3j(NG), c4j(NG),&
       Coeff_2(NG), Coeff_42(NG), Coeff_31(NG), Coeff_1c(NG)
! Local Variables: 
      integer n
      real alf, ch, sh, m0cosh, m2cosh, m1sinh, alf2, ak, a1
      real a_alf, a_alf2, a_alf2_m1sinh

       do n = 1, NG 

        ak = sqrt(sa2(n)/d2(n))
        alf = ak * hk * 0.5
        alf2 = alf*alf

        ch = cosh(alf)
        sh = sinh(alf)

        a_alf = 1./alf
        a_alf2 = a_alf*a_alf

        m0cosh = sh / alf
        m2cosh = 5.*((3. + alf2)*sh - 3.*alf*ch)*(a_alf2*a_alf)
        m1sinh = 3.*(ch - sh*a_alf)*a_alf

        a_alf2_m1sinh = a_alf2/m1sinh

        c3j(n) = (alf * ch - m1sinh )/(sh - m1sinh)

        a1 =  ch - m0cosh - m2cosh
        c4j(n) = (alf * sh - c2j*m2cosh)/a1

        Coeff_42(n) = a1*a_alf2/m2cosh
        Coeff_31(n) = (sh - m1sinh) * a_alf2_m1sinh

        Coeff_2(n) = Coeff_42(n)*c4j(n)  ! 10/35=2/7
        Coeff_1c(n) = Coeff_31(n)*c3j(n) ! 6/15=2/5

      end do

      return 
      end

      SUBROUTINE MAT_Compute_Current(NG, a0, x1, x2, x3, x4, c2j, &
        c3j, c4j, a_rl, const_ju, current)
!=====================================================================!
! Surface current for PNM and SANM                                    !
!=====================================================================!

      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: a0(NG), x1(NG), x2(NG), x3(NG), x4(NG),&
                   c2j, c3j(NG), c4j(NG), a_rl, const_ju
      REAL, INTENT(OUT)   :: current(NG)  

      INTEGER :: n
      
      DO n = 1, NG
         current(n) = - const_ju*a0(n)*&
          (x1(n) + a_rl*c2j*x2(n) + c3j(n)*x3(n) +a_rl*c4j(n)*x4(n))
      END DO

      RETURN
      END

      SUBROUTINE MAT_Compute_Flux_Surface(NG, fm0, x1, x2, x3, x4, &
       a_rl, flux_surface)
!=====================================================================!
! Surface flux for PNM and SANM                                       !
!=====================================================================!

      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: fm0(NG), x1(NG), x2(NG), x3(NG), x4(NG),&
                    a_rl
      REAL, INTENT(OUT)   :: flux_surface(NG)  

      INTEGER :: n
      
      DO n = 1, NG
         flux_surface(n) =&
          (fm0(n) + a_rl*x1(n) + x2(n) +a_rl*x3(n) + x4(n))
      END DO

      RETURN
      END


      SUBROUTINE MAT_Compute_Coupling_Coeff_Bound(NG, mat_fd, current,&
       fm0, sq, a_rl, mat_nod)
!=====================================================================!
! Nodal coupling coefficients for a boundary node (Smith formulation) !
!=====================================================================!
      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: fm0(NG), mat_fd(NG), current(NG),&
                  sq, a_rl
      REAL, INTENT(OUT)   :: mat_nod(NG)  

      INTEGER :: n

!left
!        MAT_Nod(n,nl)  = -MAT_FD(n, nl) - Cur_ANM(n)*sq/fm0(n))
!right
!          MAT_Nod(n,nl) = MAT_FD(n, nl) - Cur_ANM(n)*sq/fm0(n)

        DO n = 1, NG  
             MAT_Nod(n) = a_rl*MAT_FD(n) - current(n)*sq/fm0(n)
        END DO 

      RETURN
      END

      SUBROUTINE MAT_Compute_Coupling_Coeff(NG, mat_fd, current,&
       fm0, flux_left, sq, mat_nod)
!=====================================================================!
! Nodal coupling coefficients for a two-node problem                  ! 
!                              (Smith formulation)                    !
!=====================================================================!
      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: fm0(NG), mat_fd(NG), current(NG),&
                  sq, flux_left(NG)
      REAL, INTENT(OUT)   :: mat_nod(NG)  


      INTEGER :: n

!               MAT_Nod(n,nl) = - (MAT_FD(n, nl)*(fm0(n)-Flux(n,k_left))&
!                    + current(n)*sq) /(fm0(n) + Flux(n, k_left))

        DO n = 1, NG  
             MAT_Nod(n) = - (MAT_FD(n)*(fm0(n)-flux_left(n))&
                    + current(n)*sq) /(fm0(n) + flux_left(n))
        END DO 

      RETURN
      END

      SUBROUTINE MAT_Compute_Coupling_Coeff_Moon(NG, current,&
       fm0, flux_left, hk, hk_left, d2, d2_left, sq, &
       flux_surface, mat_fd, mat_nod, d_nod_left, d_nod )
!=====================================================================!
! Nodal coupling coefficients for a two-node problem                  ! 
!                              (Moon formulation)                     !
!=====================================================================!
      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: fm0(NG), current(NG),&
                  sq, flux_left(NG), hk, hk_left, d2(NG), d2_left(NG),&
                  flux_surface(NG)
      REAL, INTENT(OUT)   :: mat_nod(NG), mat_fd(NG), &
          d_nod_left(NG), d_nod(NG)  

      INTEGER :: n
      REAL    :: denom

        DO n = 1, NG
                D_Nod_left(n) = (-hk_left*current(n)+&
                     2.*d2_left(n)*(flux_left(n)-flux_surface(n) ))/&
                     (2.*(flux_left(n)+flux_surface(n)))   
                D_Nod(n) = &
               (hk*current(n)+2.*d2(n)*(fm0(n)-flux_surface(n)))/&
               (2.*(fm0(n)+flux_surface(n)))   
                denom = ( d2(n) + D_Nod(n) )*hk_left +&
                   (d2_left(n) + D_Nod_left(n))*hk
                MAT_FD(n) = 2.*(d2(n)*d2_left(n) - &
                     D_Nod(n)*D_Nod_left(n))*sq/denom
                MAT_NOD(n) = 2.*(d2(n)*D_Nod_left(n)-&
                     d2_left(n)*D_Nod(n))*sq/denom
        END DO
! checking     
!                write(*,*) 'Internal Flux Surface =', flux_surface(n) 
!                
!                flux_surface(n)=( (d_left-D_Nod(n,k_left,2,nd))*&
!                        hk_fd*flux_left +&
!                       (d2(n)-D_Nod(n,k,1,nd))*hk_left*fm0(n) )/denom
!                write(*,*) 'Internal Flux Surface =', flux_surface(n) 
!                write(*,*) 'hk, hk_left=', hk, hk_left
!                pause 
!                write(*,*) 'Internal current =', current(n)
!               current(n) =  (- MAT_FD(n, nl)*(fm0(n)-flux_left) &
!                    - MAT_NOD(n, nl)*(fm0(n)+flux_left) )/sq
!               write(*,*) 'Internal current =', current(n)

      RETURN
      END         

      SUBROUTINE MAT_Compute_Coupling_Coeff_Bound_Moon(NG, current,&
       fm0, hk, d2,  sq, flux_surface, a_rl, fg, cg,&
       mat_fd, mat_nod, d_nod )
!=====================================================================!
! Nodal coupling coefficients for a boundary node                     ! 
!                              (Moon formulation)                     !
!=====================================================================!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NG
      REAL, INTENT(IN)    :: fm0(NG), current(NG),&
                  sq,  hk,  d2(NG), &
                  flux_surface(NG), fg(NG), cg, a_rl
      REAL, INTENT(OUT)   :: mat_nod(NG), mat_fd(NG), d_nod(NG)  
        REAL ::             flux_sf(NG), cur(NG) 


      INTEGER :: n
      REAL    :: denom

         DO n = 1, NG
                D_Nod(n) = &
             (-a_rl*hk*current(n)+2.*d2(n)*(fm0(n)-flux_surface(n)))/&
               (2.*(fm0(n)+flux_surface(n)))   
                denom = fg(n)*hk + &
                    2.*cg*(d2(n) + D_Nod(n))
                MAT_FD(n)  = 2.*fg(n)*d2(n)*sq/denom
                MAT_NOD(n) = a_rl*2.*fg(n)*D_Nod(n)*sq/denom
         END DO

! LEFT
!                D_Nod(n,k,1,nd) = &
!                 ( hk*cur_ANM(n) + 2.*d2(n)*(fm0(n)- flux_surface) )/&
!               (2.*( fm0(n)+flux_surface ))   
!                denom = fg(n,1,nd)*hk + &
!                   2.*cg(1, nd)*(d2(n) + D_Nod(n,k,1,nd))
!                MAT_FD(n, nl) = 2.*fg(n,1,nd)*d2(n)*sq/denom
!                MAT_NOD(n, nl)= -2.*fg(n,1,nd)*D_Nod(n,k,1,nd)*sq/denom
          
!                IF( a_rl == 1) THEN
! checking     
!               DO n = 1, NG
!                write(*,*) 'RIGHT Flux Surface =', flux_surface(n) 
!               flux_sf(n)=(2*cg*(d2(n)-D_Nod(n))*fm0(n))/&
!                ( 2*cg*(D_Nod(n)+d2(n))+fg(n)*hk )
!                 write(*,*) 'hk =', hk
!                 write(*,*) 'cg =', cg
!                 write(*,*) 'fg =', fg(n)
!                 write(*,*) 'd2 =', 2*d2(n)/hk
!                 write(*,*) 'd_nod =', 2*d_nod(n)/hk
!                 pause
!                write(*,*) 'RIGHT Flux Surface =', flux_sf(n) 
!                write(*,*) 'RIGHT current =', current(n)
!                 cur(n)=a_rl*(MAT_FD(n)-a_rl*MAT_NOD(n))*fm0(n) 
!                 write(*,*) 'RIGHT current =', cur(n)/sq
!                END DO  
!                 pause 
!                END IF

! checking     
!                IF( a_rl == -1 ) THEN
!                DO n = 1, NG
!                write(*,*) 'LEFT Flux Surface =', flux_surface(n) 
!                 flux_sf(n)=(2*cg*(d2(n) - D_Nod(n))*fm0(n))/&
!                   (2*cg*(D_Nod(n)+d2(n))+fg(n)*hk)
!                 write(*,*) 'LEFT Flux Surface =', flux_sf(n) 
!                pause 
!                write(*,*) 'LEFT current =', current(n)
!               cur(n)=a_rl*(MAT_FD(n)-a_rl*MAT_NOD(n))*fm0(n) 
!               write(*,*) 'LEFT current =', cur(n)/sq
!               END DO 
!                pause 
!               END IF            
      RETURN
      END


      SUBROUTINE test_rhs_symmatric_node(NG, fm0, x2, x4, a0, &
        c2j, c4j, k)
! Input
      INTEGER NG, k
      REAL fm0(NG), x2(NG), x4(NG), a0(NG), c2j, c4j(NG)
! Local
      REAL rhs_flux(NG), rhs_current(NG) 
      INTEGER sign_flux, n

      IF(k.eq.2) THEN
        sign_flux = -1
        DO n = 1, NG
          rhs_flux(n) = 0.
          rhs_current(n) = 0.
        END DO
      ELSE IF(k.eq.3) THEN
        sign_flux = +1 
      END IF
      
      DO n = 1, NG
        rhs_flux(n) = rhs_flux(n)+sign_flux*( fm0(n) + x2(n) + x4(n) )
        rhs_current(n) = rhs_current(n)+&
                             a0(n)*( c2j*x2(n) + c4j(n)*x2(n) )
      END DO       

      IF(k.eq.3) THEN
         WRITE(*,*) 'rhs_flux ='
         WRITE(*,'(4E12.5)') rhs_flux
         WRITE(*,*) 'rhs_current ='
         WRITE(*,'(4E12.5)') rhs_current
      END IF

      RETURN
      END


      SUBROUTINE test_nodal_balance(NG, fm0, x2, x4, xappa, &
        c2j, c4j, trans0)
! Input
      INTEGER NG
      REAL fm0(NG), x2(NG), x4(NG), xappa(NG, NG), c2j, c4j(NG),&
       trans0(NG)
! Local
      REAL res(NG)
      INTEGER n, m

        DO n = 1, NG
          res(n) = -( c2j*x2(n) + c4j(n)*x4(n) ) + trans0(n)
          DO m = 1, NG
             res(n) = res(n) + xappa(n,m)*fm0(m)
          END DO
! relative value
          res(n) = res(n)/fm0(n)
        END DO

         WRITE(*,*) 'Residual of the Nodal Balance'
         WRITE(*,'(4E12.5)') res

      RETURN
      END

      SUBROUTINE test_second_order_moment_eq(NG,  x2, x4, xappa, &
        coeff_22, coeff_42, transq)
! Input
      INTEGER NG
      REAL  x2(NG), x4(NG), xappa(NG, NG), coeff_22, coeff_42(NG),&
       transq(NG)
! Local
      REAL res(NG)
      INTEGER n, m

        DO n = 1, NG
          res(n) = -( x4(n)/coeff_42(n) + coeff_22*x2(n) ) + transq(n)
          DO m = 1, NG
             res(n) = res(n) + xappa(n,m)*x2(m)
          END DO
! relative value
          res(n) = res(n)/x2(n)
        END DO

         WRITE(*,*) 'Residual of the second order moment weighting eq.'
         WRITE(*,'(4E12.5)') res

      RETURN
      END

      

