      subroutine MAT_Compute_Eigenv_2x2(xappa, eigenv_xappa)
!=====================================================================*
! Compute Eigenvalue of 2x2 matrix                                    *
! is used to compute eigenvalue of buckling matrix in the case        *
!  of 2 group ANM                                                     *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
      integer NG
      parameter(NG=2)
! Input: xappa(2,2) - matrix
      real xappa(NG,NG)
! Output: eigenv_xappa(2) - eigenvalues
      real eigenv_xappa(NG)
! Local Variables:
!      determ - Determinant of the matrix
      real determ, d

      determ = xappa(1,1)*xappa(2,2) - xappa(1,2)*xappa(2,1)
      d = (xappa(1,1)+xappa(2,2))**2 - 4.*determ
      eigenv_xappa(2)  = 0.5*(xappa(1,1) + xappa(2,2) + sqrt(d))
      eigenv_xappa(1)  = 0.5*(xappa(1,1) + xappa(2,2) - sqrt(d))

      return
      end

      real function MAT_F1_ANM(gamma,eigenv)
!=====================================================================*
!     Response function for the average neutron flux                  *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
      real gamma, eigenv

      MAT_F1_ANM = -1./gamma**2 + eigenv

      return
      end

      real function MAT_F2_ANM(gamma,eigenv)
!=====================================================================*
! Response function for the neutron flux on the right boundary        *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
      real gamma, eigenv

      MAT_F2_ANM = 1./gamma

      return
      end

      real function MAT_F3_ANM(gamma,eigenv)
!=====================================================================*
! Response function for the average transverse leakage                *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
      real gamma, eigenv

      MAT_F3_ANM = 1. - (1.-gamma)/(eigenv*gamma**2)

      return
      end

      real function MAT_F4_ANM(gamma,eigenv)
!=====================================================================*
! Response function for the 1st TRL expansion coefficient             * 
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
      real gamma, eigenv

      MAT_F4_ANM = (1. - gamma)/(gamma*eigenv)

      return
      end

      real function MAT_F5_ANM(gamma,eigenv)
!=====================================================================*
! Response function for the 2nd TRL expansion coefficient             * 
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
      real gamma, eigenv

      MAT_F5_ANM = -(3. - gamma*(3.+eigenv))/((gamma*eigenv)**2)

      return
      end

      subroutine MAT_Set_MatFunc_ANM(xappa, Matr_Func)
!=====================================================================*
! Computing Response Matrix Functions for the node                    *
!       ANM (2 neutron energy groups)                                 * 
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
! Input: Eigenvalues of the Matrix matrix B^2
      integer NG ! NG = 2
      parameter (NG=2)
! eigenv(2) - eigenvalues of the buckling matrix
! xappa(2,2) -buckling matrix
      real eigenv(NG), xappa(NG,NG)
! Output:
      real Matr_Func(NG,NG,5)
! Response Matrix Functions with respect:
! 1) Average Flux 2) Right Boundary Flux 3) Average TRL 4) 1st TRL
! expansion coefficient 5) 2nd TRL expansion coefficient

! Local Variables:
! gamma - generic function: tanh(sqrt(x))/sqrt(x) if x > 0;
!                           tan(sqrt(-x))/sqrt(-x) if x <0.
! lambda -  generic function sqrt(x) if x > 0; sqrt(-x) if x < 0;
      real gamma(NG), lambda(NG) ! 

      integer n, m, i

      real a1(5), a2(5)
! Unity - Unity Marix
      real Unity(NG,NG)
      data Unity / 1., 0., 0., 1./ ! NG = 2
      real MAT_F1_ANM, MAT_F2_ANM, MAT_F3_ANM, &
              MAT_F4_ANM, MAT_F5_ANM
      real denom


      call MAT_Compute_Eigenv_2x2(xappa, eigenv)

      do n = 1, NG
         if(Eigenv(n).LT. 0) then
            lambda(n) = SQRT(-eigenv(n))
            gamma(n) = TAN(lambda(n))/lambda(n)
         else
            lambda(n) = SQRT(eigenv(n))
            gamma(n) = TANH(lambda(n))/lambda(n)
         end if
      end do

      denom = Eigenv(1) - Eigenv(2)
      a1(1) = (MAT_F1_ANM(gamma(1),eigenv(1))-&
               MAT_F1_ANM(gamma(2),eigenv(2)))/denom
      a1(2) = (MAT_F2_ANM(gamma(1),eigenv(1))-&
               MAT_F2_ANM(gamma(2),eigenv(2)))/denom
      a1(3) = (MAT_F3_ANM(gamma(1),eigenv(1))-&
               MAT_F3_ANM(gamma(2),eigenv(2)))/denom
      a1(4) = (MAT_F4_ANM(gamma(1),eigenv(1))-&
               MAT_F4_ANM(gamma(2),eigenv(2)))/denom
      a1(5) = (MAT_F5_ANM(gamma(1),eigenv(1))-&
               MAT_F5_ANM(gamma(2),eigenv(2)))/denom
                           
      a2(1) = (-Eigenv(2)*MAT_F1_ANM(gamma(1),eigenv(1))+&
                Eigenv(1)*MAT_F1_ANM(gamma(2),eigenv(2)))/denom
      a2(2) = (-Eigenv(2)*MAT_F2_ANM(gamma(1),eigenv(1))+&
                Eigenv(1)*MAT_F2_ANM(gamma(2),eigenv(2)))/denom
      a2(3) = (-Eigenv(2)*MAT_F3_ANM(gamma(1),eigenv(1))+&
                Eigenv(1)*MAT_F3_ANM(gamma(2),eigenv(2)))/denom
      a2(4) = (-Eigenv(2)*MAT_F4_ANM(gamma(1),eigenv(1))+&
                Eigenv(1)*MAT_F4_ANM(gamma(2),eigenv(2)))/denom
      a2(5) = (-Eigenv(2)*MAT_F5_ANM(gamma(1),eigenv(1))+&
                Eigenv(1)*MAT_F5_ANM(gamma(2),eigenv(2)))/denom

      do i = 1, 5
         do n = 1, NG
            do m = 1, NG
               Matr_Func(n,m,i) = xappa(n,m)*a1(i) + &
                                    a2(i)*Unity(n,m)
            end do
         end do
      end do



      return
      end
              

      subroutine MAT_Set_MatrFlux&
         (NG, a0, F0, S0, S1, S2, Matr_Flux, &
          Matr_Flux_k, RHS, RHS_K, Matr_Funct, adf_k, adf_k_plus_1)
!=====================================================================*
! Computing Matrix and RHS of the Equation for the Neutron Flux       * 
!   on Left  Boundary of Node  (Normal Tewo-node Problem)             *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
! Input:
      integer NG ! 
! a0(NG) - nondimensional diffusion coefficient
! f0(NG) - avergae neutron flux
! s0(NG) - average TRL (nondimensional)
! s1(NG) - 1st TRL Expansion Coefficient
! s2(NG) - 2nd TRL expansion Coefficient
! Matr_Funct(NG,NG,5) - Response matrix functions
! Matr_Flux_k(NG,NG) - a part of the matrix for the previous node
! RHS_K(NG) - a part of the RHS for the previous node
      real a0(NG), f0(NG), s0(NG), s1(NG), s2(NG), Matr_Funct(NG,NG,5), &
             Matr_Flux_k(NG,NG), RHS_K(NG)
      REAL adf_k(NG), adf_k_plus_1(NG)  
! Output:
! Matr_Flux(NG,NG) - Matrix for the Left Boundary Flux 
! RHS(NG) - Right Hand Side for the Left Boundary Flux
      real Matr_Flux(NG,NG), RHS(NG)

!     Matr_k(NG,NG) - a part of the matrix for the next two-node problem
!     RHS_K(NG) - a part of the RHS for the next two-node problem
! Local Variables:
      real tmp, RHS_Odd, RHS_Even
      integer n,m

      do n = 1, NG

         RHS_Odd = 0.
         RHS_Even = 0.
         do m = 1, NG
            Matr_Flux(n,m) = -a0(n)*Matr_Funct(n,m,2)
            tmp = Matr_Flux_K(n,m)
            Matr_Flux_K(n,m) = Matr_Flux(n,m)
            Matr_Flux(n,m) = Matr_Flux(n,m) + &
               tmp*adf_k_plus_1(m)/adf_k(m)
       
            RHS_Odd = RHS_Odd + Matr_Funct(n,m,4)*s1(m)
            RHS_Even = RHS_Even + Matr_Funct(n,m,1)*f0(m) + &
                Matr_Funct(n,m,3)*s0(m) + Matr_Funct(n,m,5)*s2(m)
         end do
         RHS_Even = a0(n)*RHS_Even
         RHS_Odd = a0(n)*RHS_Odd
         tmp = RHS_K(n)
         RHS_K(n) = RHS_Even + RHS_Odd 
         RHS(n) = ( RHS_Even - RHS_Odd ) + &
             tmp
      end do

      return
      end


      subroutine MAT_Set_MatrFlux_Bound(NG, a0, F0, S0, S1, S2, &
      Matr_Flux, Matr_Flux_k, RHS, RHS_K, Matr_Funct, a_rl, ag, bg)
!=====================================================================*
! Computing Matrix and RHS of the Equation for the Neutron Flux       * 
!  on the Left or Right Boundary of Node  (Boundary Nodes)            *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
! Input:
      integer NG ! NG = 2
! a0(NG) - nondimensional diffusion coefficient
! f0(NG) - avergae neutron flux
! s0(NG) - average TRL (nondimensional)
! s1(NG) - 1st TRL Expansion Coefficient
! s2(NG) - 2nd TRL expansion Coefficient
! Matr_Funct(NG,NG,5) - Response matrix functions
! a_rl = +1 for the right node; a_rl = -1 for the left node;
! bg, ag(NG) - coefficients of the mixed type boudary condition
      real a0(NG), f0(NG), s0(NG), s1(NG), s2(NG), Matr_Funct(NG,NG,5), &
                 a_rl, bg, ag(NG)
! Output:
! Matr_Flux(NG,NG) - Matrix for the Left Boundary Flux 
! RHS(NG) - Right Hand Side for the Left Boundary Flux
! Matr_k(NG,NG) - a part of the matrix for the next two-node problem
! RHS_K(NG) - a part of the RHS for the next two-node problem
      real Matr_Flux(NG,NG), RHS(NG)
      real Matr_Flux_k(NG,NG), RHS_K(NG)
! Local Variables:
      real RHS_Odd, RHS_Even
      integer n,m

      do n = 1, NG

         RHS_Odd = 0.
         RHS_Even = 0.

         do m = 1, NG
            Matr_Flux(n,m) = -a0(n)*Matr_Funct(n,m,2)
            Matr_Flux_K(n,m) = Matr_Flux(n,m)
            Matr_Flux(n,m) = Matr_Flux(n,m)*bg
            RHS_Odd = RHS_Odd + Matr_Funct(n,m,4)*s1(m)
            RHS_Even = RHS_Even + Matr_Funct(n,m,1)*f0(m) + &
                Matr_Funct(n,m,3)*s0(m) + Matr_Funct(n,m,5)*s2(m)
         end do

         Matr_Flux(n,n) = Matr_Flux(n,n) - ag(n)
         RHS_Even = a0(n)*RHS_Even
         RHS_Odd = a0(n)*RHS_Odd
         RHS_K(n) = RHS_Even + RHS_Odd 
         RHS(n) = (RHS_Even + a_rl*RHS_Odd)*bg

      end do

      return
      end


      subroutine MAT_Compute_Current_ANM(NG, a_rl, a0, f0, &
                   Flux, s0, s1, s2, Matr_Funct, current )
!=====================================================================*
! Computing Neutron Curreent                                          * 
!  on the Left or Right Boundary of Node                              *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
! Input:
      integer NG ! NG = 2
! a0(NG) - nondimensional diffusion coefficient
! f0(NG) - avergae neutron flux
! s0(NG) - average TRL (nondimensional)
! s1(NG) - 1st TRL Expansion Coefficient
! s2(NG) - 2nd TRL expansion Coefficient
! Matr_Funct(NG,NG,5) - Response matrix functions
! a_rl = +1 for the right node; a_rl = -1 for the left node;
      real a0(NG), f0(NG), Flux(NG), s0(NG), s1(NG), s2(NG), &
             Matr_Funct(NG,NG,5), a_rl
! Output:
! Cur_ANM(n)(NG) - Neutron Cur_ANM(n) on the Left or Right Boundary of the Node
      real current(NG)
! Local Variables:
      integer n,m

      do n = 1, NG
         current(n) = 0.
         do m = 1, NG
            current(n) = current(n)  + Matr_Funct(n,m,1)*f0(m) +&
               Matr_Funct(n,m,2)*Flux(m) + Matr_Funct(n,m,3)*s0(m) +&
               a_rl*Matr_Funct(n,m,4)*s1(m) + Matr_Funct(n,m,5)*s2(m)
         end do
            current(n) = - a_rl*a0(n)*current(n)
      end do

      return
      end

      real function MAT_F1_PNM(t)
      real t
        MAT_F1_PNM = &
           (-66. + 2.*t + 2875./(15. + t) - 2744./(21. + 2.*t))/5.
      return
      end

      real function MAT_F2_PNM(t)
      real t

        MAT_F2_PNM = (3.*(5. + 2.*t))/(15. + t)

      return
      end

      real function MAT_F3_PNM(t)
! f3 for zero moment og TRL
      real t

        MAT_F3_PNM = &
           (2.*(525. + t*(15. + 2.*t)))/(5.*(15. + t)*(21. + 2*t))
      return
      end

      real function MAT_F4_PNM(t)
! f4 for 1-st moment of TRL
      real t

        MAT_F4_PNM = 5./(15. + t)

      return
      end

      real function MAT_F5_PNM(t)
! f5 for 2-nd moment of TRL
      real t

       MAT_F5_PNM = &
         (21.*(5. + 2.*t))/(5.*(15. + t)*(21. + 2.*t))

      return
      end



      subroutine MAT_Set_MatFunc_PNM(xappa, Matr_Func)
!=====================================================================*
! Computing Response Matrix Functions for the node                    *
!       PNM (2 neutron energy groups)                                 * 
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x
! Computing Matricx Functions
! Input: Eigenvalues of the Matrix matrix B^2
!        local variables 
      implicit none
! Input:
      integer NG ! NG = 2
      parameter (NG=2)
! eigenv(2) - eigenvalues of the buckling matrix
! xappa(2,2) -buckling matrix
      real eigenv(NG), xappa(NG,NG)
! Output:
      real Matr_Func(NG,NG,5)
! Response Matrix Functions with respect:
! 1) Average Flux 2) Right Boundary Flux 3) Average TRL 4) 1st TRL
! expansion coefficient 5) 2nd TRL expansion coefficient

! Local Variables:

      integer n, m, i

      real a1(5), a2(5)
! Unity - Identity Marix
      real Unity(NG,NG)
      data Unity / 1., 0., 0., 1./ ! NG = 2
      real MAT_F1_PNM, MAT_F2_PNM, MAT_F3_PNM, &
               MAT_F4_PNM, MAT_F5_PNM
      real denom


      call MAT_Compute_Eigenv_2x2(xappa, eigenv)

      denom = Eigenv(1) - Eigenv(2)
      a1(1) = (MAT_F1_PNM(eigenv(1))-&
                MAT_F1_PNM(eigenv(2)))/ denom
      a1(2) = (MAT_F2_PNM(eigenv(1))-&
                MAT_F2_PNM(eigenv(2)))/denom
      a1(3) = (MAT_F3_PNM(eigenv(1))-&
                MAT_F3_PNM(eigenv(2)))/denom
      a1(4) = (MAT_F4_PNM(eigenv(1))-&
                MAT_F4_PNM(eigenv(2)))/denom
      a1(5) = (MAT_F5_PNM(eigenv(1))-&
                MAT_F5_PNM(eigenv(2)))/denom
                           
      a2(1) = (-Eigenv(2)*MAT_F1_PNM(eigenv(1))+&
            Eigenv(1)*MAT_F1_PNM(eigenv(2)))/denom
      a2(2) = (-Eigenv(2)*MAT_F2_PNM(eigenv(1))+&
            Eigenv(1)*MAT_F2_PNM(eigenv(2)))/denom
      a2(3) = (-Eigenv(2)*MAT_F3_PNM(eigenv(1))+&
            Eigenv(1)*MAT_F3_PNM(eigenv(2)))/denom
      a2(4) = (-Eigenv(2)*MAT_F4_PNM(eigenv(1))+&
            Eigenv(1)*MAT_F4_PNM(eigenv(2)))/denom
      a2(5) = (-Eigenv(2)*MAT_F5_PNM(eigenv(1))+&
            Eigenv(1)*MAT_F5_PNM(eigenv(2)))/denom

      do i = 1, 5
         do n = 1, NG
            do m = 1, NG
               Matr_Func(n,m,i) = xappa(n,m)*a1(i) + &
                                    a2(i)*Unity(n,m)
            end do
         end do
      end do


      return
      end

      subroutine MAT_Solve_Two_Node_ANM_OLD(k_eff_2)
!=====================================================================*
! Solution of the Two-Node Problems and Computing New Coupling        *
!     Coefficients (Polynomial & Semi_Analytical Nodal Method)        *
! Last Update              Slava (c) April 1998                       *
!=====================================================================*
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'

! Input: TRL_0_xyz, trl_qla, all XS_, Flux, Matr_FD, Diag_Tot
      real k_eff_2
! Output: MAT_Tot(NG, NE_T, N_TOT), MAT_Diag_Tot(NG, N_TOT)

! Local Variables:
      real sa2(NG), sik2(NG,NG), sf2(NG), d2(NG), sf0(NG,NG) ! XS
      real a0(NG), a01(NG)
      real a_rl 
      real fm0(NG)

      integer k, n, kt, nd, nd_right, m, k_left, k_right,  nl, NN
      real a_vol, hk, sq !,  MAT_Nod

      real  xappa(NG,NG), trans0(NG), transq(2*NG)

! ANM Method
      real Matr_Funct(NG,NG,5)
      real Matr_Flux_k(NG,NG), Matr_Flux(NG,NG), RHS_Flux(NG), &
                 RHS_Flux_K(NG), Cur_ANM(NG)
      Logical ANM
! Nonlinear Iterations by Moon
      REAL denom, flux_left, d_left, hk_left
      REAL adf_k(NG), adf_k_plus_1(NG)
!tmp
!      real Flux_Rodded(NG), Flux_Unrodded(NG)

      adf_k(:)=1. 
      adf_k_plus_1(:)=1.

      nl = 0

      do nd = 1, NDD

       nd_right = Index_Neib(nd,2)

       do kt = 1, N_TOT_DIR(nd)

          nl = nl + 1 ! Left Interface

          k = Numb(kt, nd)
          a_vol = 1./volume(k)
          sq = s_xyz(nd, k) 
          hk  = h_xyz(nd, k)

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
!A         end if

          do n = 1, NG
             trans0(n) = TRL_0_xyz(n, k, nd)*a01(n)
             transq(n) = trl_qla(1, n, k, nd)*a01(n)
             transq(n+NG) = trl_qla(2, n, k, nd)*a01(n)
          end do

          do n = 1, NG
             do m = 1, NG
                   xappa(n,m) = - sf0(n,m)*a01(n)
             end do
                   xappa(n,n) = xappa(n,n) + sa2(n)*a01(n)
          end do



       write(*,*) ' k= ', k, 'nd =', nd
       write(*,*) xappa(1,1),xappa(1,2)
       write(*,*) xappa(2,1), xappa(2,2)
!      write(*,*)  'Eigenv =', eigen22
!      write(*,*) 'sqrt(Eigenv) =', sqrt(abs(eigen22))
!      write(*,*) 'alfa=', sqrt(sa2(1)*a01(1)), sqrt(sa2(2)*a01(2))

         if(Nodal_Method.EQ."ANM") then

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

         end if


         k_left = Neib(nd,k)

         if(k_left.eq.I_BOUND_NODE.AND. i_dr(1,nd) /=2 ) then        
! Left Boundary Node NOT ccomputed in the case of the REPEATED structure
           a_rl = -1.         

           call MAT_Set_MatrFlux_Bound(NG, a0, fm0, trans0(1), &
                 transq(1), transq(1+NG), Matr_Flux, &
                 Matr_Flux_k, RHS_Flux, RHS_Flux_K, Matr_Funct, &
                 a_rl, fg(1,1,nd), cg(1,nd))

           if(NG.eq.2) then           
               call MSC_Solve2x2(Matr_Flux, RHS_Flux)
           else
               NN = NG
               call MSC_LU_Solve(Matr_Flux, NN, RHS_Flux)
           end if
          

           call MAT_Compute_Current_ANM(&
                 NG, a_rl, a0, fm0, RHS_Flux, trans0(1), transq(1),&
                 transq(NG+1), Matr_Funct, Cur_ANM)
      
            do n = 1, NG

            IF( NonlinearIterations.EQ."Smith") THEN
               MAT_Nod(n,nl)  =-(MAT_FD(n, nl) + Cur_ANM(n)*sq/fm0(n))
            ELSE IF( NonlinearIterations.EQ."Moon") THEN
! we are working in cartesian coordinates
            END IF

!              MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_Nod(n,nl)

!              MAT_Tot(n, nd, k) = MAT_FD(n, nl) - MAT_Nod(n,nl)

           end do

         else


! the normal case of the two-node problem

! calculation of the  nodal coupling coefficients
         a_rl = -1. ! Left interface

         call MAT_Set_MatrFlux(NG, a0, fm0, trans0(1), &
                transq(1),  transq(1+NG), Matr_Flux, &
                Matr_Flux_k, RHS_Flux, RHS_Flux_K, Matr_Funct,&
                adf_k, adf_k_plus_1)

           if(NG.eq.2) then           
               call MSC_Solve2x2(Matr_Flux, RHS_Flux)
           else
               NN = NG
               call MSC_LU_Solve(Matr_Flux, NN, RHS_Flux)
           end if

         call MAT_Compute_Current_ANM(NG, a_rl, a0, fm0,&
                RHS_Flux, trans0(1),  transq(1), transq(NG+1), &
                Matr_Funct, Cur_ANM)

!        if(nd.eq.3.and.k.eq.3) then
!         call Compute_Neutron_Flux(0, k, &
!         sf2, sa2, sik2, d2, &
!         sf2, sa2, sik2, d2, &
!         0.5, k_eff_2, xp, &
!         Flux_Rodded, Flux_Unrodded )
!         end if

!          if(nd.eq.3.and.k.eq.323) then 
!             open(2, file = 'Output_Debug/LENYA.dat', &
!                         status = 'unknown')
!             write(2,*) 'k_eff =', k_eff_2
!             write(2,*) 'Unrodded Node'
!             write(2,*) 'hk =', hk, 'vol =', volume(k)
!             write(2,*) "SA =", sa2
!             write(2,*) "SF =", sf2
!             write(2,*) "D =", d2
!             write(2,*) "SIK =", SIK2 
!             write(2,*) 'Average Flux = ', fm0
!             write(2,*) 'Transverse Leakage 0 =', TRL_0_xyz(:, k, nd)
!             write(2,*) 'Transverse Leakage 1 =', trl_qla(1, :, k, nd)
!             write(2,*) 'Transverse Leakage 2 =', trl_qla(2, :, k, nd)
!             write(2,*) 'Left Interface Flux =', RHS_Flux
!             write(2,*) 'Left Interface Cur_ANM(n) =', Cur_ANM
             
!             close(2)

!          else if(nd.eq.3.and.k.eq.387) then 
!             open(2, file = 'Output_Debug/LENYA.dat', &
!                         status = 'unknown', access ='append')
!             write(2,*) 'Rodded Node'
!             write(2,*) 'hk =', hk, 'vol =', volume(k)
!             write(2,*) "SA =", sa2
!             write(2,*) "SF =", sf2
!             write(2,*) "D =", d2
!             write(2,*) "SIK =", SIK2 
!             write(2,*) 'Average Flux = ', fm0
!             write(2,*) 'Transverse Leakage 0 =', TRL_0_xyz(:, k, nd)
!             write(2,*) 'Transverse Leakage 1 =', trl_qla(1, :, k, nd)
!             write(2,*) 'Transverse Leakage 2 =', trl_qla(2, :, k, nd)
!             write(2,*) 'Left Interface Flux =', RHS_Flux
!             write(2,*) 'Left Interface Cur_ANM(n) =', Cur_ANM
!             close(2)

!          else if(nd.eq.3.and.k.eq.451) then 
!             open(2, file = 'Output_Debug/LENYA.dat', &
!                         status = 'unknown', access ='append')
!             write(2,*) 
!             write(2,*) 'Right Interface Flux =', RHS_Flux
!             write(2,*) 'Right Interface Cur_ANM(n) =', Cur_ANM

!             close(2)

!          end if 
             


            do n = 1, NG

!               if(nl.eq.42.OR.nl.eq.43) then
!                  write(*,*) 'nl=', nl, 'Cur_ANM(n) =', Cur_ANM
!                  pause
!               end if

              IF( NonlinearIterations.EQ."Smith") THEN
               MAT_Nod(n,nl) =-( MAT_FD(n, nl)*(fm0(n)-Flux(n,k_left))&
                    + Cur_ANM(n)*sq) /(fm0(n) + Flux(n, k_left))
              ELSE IF( NonlinearIterations.EQ."Moon") THEN
              END IF
              

!               MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_Nod(n,nl)

!               MAT_Tot(n, nd, k) = MAT_FD(n, nl) - MAT_Nod(n,nl)

!               MAT_Diag_Tot(n, k_left) = MAT_Diag_Tot(n, k_left) - &
!                                                        MAT_Nod(n,nl)

!               MAT_Tot(n, nd_right, k_left) = MAT_FD(n, nl) + &
!                                                        MAT_Nod(n,nl)

             end do

         end if ! k_left = 0 & k_left /= 0

         k_right = Neib(nd_right, k)

         if(k_right.eq.I_BOUND_NODE) then        
! Right Boundary Node
          nl = nl + 1 ! Right Interface

         IF( i_dr(2,nd) /= 2 ) THEN 

          a_rl = 1.         

!ANM
           call MAT_Set_MatrFlux_Bound(NG, a0, fm0, trans0(1),&
                 transq(1), transq(1+NG), Matr_Flux, &
                 Matr_Flux_k, RHS_Flux, RHS_Flux_K, Matr_Funct, &
                 a_rl, fg(1,2,nd), cg(2,nd))

           if(NG.eq.2) then           
               call MSC_Solve2x2(Matr_Flux, RHS_Flux)
           else
               NN = NG
               call MSC_LU_Solve(Matr_Flux, NN, RHS_Flux)
           end if

           call MAT_Compute_Current_ANM(NG, a_rl, a0, fm0, &
                RHS_Flux, trans0(1), transq(1), transq(NG+1), &
                Matr_Funct, Cur_ANM)
       
! calculation of the nodal coupling coefficients
          do n = 1, NG

             IF( NonlinearIterations.EQ."Smith") THEN
                MAT_Nod(n,nl) = MAT_FD(n, nl) - Cur_ANM(n)*sq/fm0(n)
             ELSE IF( NonlinearIterations.EQ."Moon") THEN
             END IF 
             
!             MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) - MAT_Nod(n,nl)
                     
!             MAT_Tot(n,nd_right,k) = MAT_FD(n, nl) + MAT_Nod(n,nl)


          end do ! NG
          END IF ! i_dr(2,nd) /= 2
         end if ! k_right = 0                 

       end do !k
      end do ! NDD


      if(Debug) then
         call MAT_DBG_Check_Matr_tot
!        call Moment_Write
!        pause 'Matrix Written'
      end if

      return
      end


