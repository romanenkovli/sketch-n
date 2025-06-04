      subroutine CHB_Set_First_Residual
!=====================================================================*
! Initial Residual for the Neutron Kinetics Calculatiuons             *
! (c) Slava 15.IV.1998                                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'      
! Input: Flux(NG, N_TOT) - Computed Neutron Flux 
!        Flux_k(NG, N_TOT) - Initial Approximation of the Neutron Flux
! Output: Flux(NG, N_TOT) - Initial Residual
!         Flux(NG, N_TOT) - Neutron Flux (k-1) for the Chebyshev 
!          Iterative Procedure
      integer n, k

      do k = 1, N_TOT
         do n = 1, NG
            Flux(n, k) = Flux(n, k) - Flux_k(n, k)
            Flux_k1(n, k) = 0. ! Flux_k(n, k)
         end do
      end do

      return 
      end

      subroutine CHB_KIN_Init_Iterations(i_out)
!=====================================================================*
! Initial Iterative Approximation for the Neutron Flux During Neutron *
!    Kinetics Calculatiuons                                           *
! (c) Slava 15.IV.1998                                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'         

! Input: I_CSA - if the CSA method = 1, if CSI 0
!        i_out = 0 for the 1st outer iteration at the time step 
!        Flux(NG, N_TOT) - Neutron Flux at the Previous Nonlinear
!        Iteration 
!        Flux_k(NG, N_TOT) - Added Extrapolation Error
      integer i_out
! Output: Flux(NG, N_TOT) - Initial Neutron Flux Approximation
! FLUX_K(NG, N_TOT)  = Flux(NG, N_TOT) - used to compute the residual
! Local Variables: 
      integer n, k
      real harm_norm, flux_norm, alfa_res
      logical New_Initial_Flux
      parameter (New_Initial_Flux = .True.)       

      if(iter_solver.eq. "CSA" ) then

      harm_norm = 0. 
      flux_norm = 0.

      if(New_Initial_Flux) then
         do k = 1, N_TOT
            do n = 1, NG
            flux_norm = flux_norm + Flux(n,k)*Flux_k(n,k)*&
                   Volume(k)
            harm_norm = harm_norm + Flux_k(n,k)*Flux_k(n,k)*&
                   Volume(k)
            end do
         end do
      else     
        do k = 1, N_TOT
         do n = 1, NG
          flux_norm = flux_norm + xs_al(n,k)*Flux(n,k)*Flux_a(n,k)*&
                   Volume(k)
          harm_norm = harm_norm + xs_al(n,k)*Flux_k(n,k)*Flux_a(n,k)*&
                   Volume(k)
         end do
        end do
      end if

      alfa_res = Flux_Norm/Harm_Norm
!cc   write(*,*) 'alfa_res = ', alfa_res
!     pause

        if(i_out.eq.0) then
!           write(*,*) 'Initial Neutron Flux: 0'
           do k = 1, N_TOT
              do n = 1, NG
!TMP                 Flux(n,k) = 0.
                  Flux(n, k) = Flux(n, k) - alfa_res*Flux_k(n,k)
              end do
          end do
        else

! We Continue from the last iterative approximation

          do k = 1, N_TOT
             do n = 1, NG
                Flux(n, k) =  Flux(n, k) - Flux_k(n, k)
             end do
          end do

        end if
      end if
! In the case of the CSI method Flux(0) = Flux

      do k = 1, N_TOT
         do n = 1, NG
            Flux_k(n, k) = Flux(n, k)
         end do
      end do

      return
      end

      subroutine CHB_Init_Iterations
!=====================================================================*
! Initial Chebyshev Iteration Parameters                              *
!      Slava (c) 23.II.1998                                           *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: xme_ini, kin_k_ef

! Output:
! npolin - degree of the Chebyshev polinomial
! delnp - norm of the residula vector used in estimation of the eigenvalue
! is_polin - number of the Chebyshev polynomial sequences
! xme_ - estimate of the eigenvalue used to compute the Chebyshev iteration
! parameters

! Local variables:
      integer n, k


      npolin = -1
      delnp = 1
      is_polin = 0
      T_Cheb = 0

      if(Problem_Type.EQ."Kinetics") xme_ = xme_ini*kin_k_ef

      do k = 1, N_TOT
        do n = 1, NG
            Flux_n(n,k) = Flux(n,k)
            Flux_n1(n,k) = 0. ! Flux(n,k)
        end do
       end do

      return
      end


      subroutine CHB_Iterate_CSI_Residual(N, x, xn, xn1, &
         ro_cheb, gamma_cheb)  

!=====================================================================*
!        A Chyebyshev Iteration for PseudoResidual                    *
! Hageman & Young "Applied Iterative Methods", Academic Press, 1981   *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x(N), xn(N), xn1(N)
      real ro_cheb, gamma_cheb
! Output:
!      real x(N), xn(N), xn1(N)
! Local Variables
      real aa, bb, cc
      integer k

      aa = ro_cheb * gamma_cheb
      bb = ro_cheb * (1. - gamma_cheb)
      cc = 1. - ro_cheb

      do k = 1, N

       x(k) = aa*x(k) + bb*xn(k) + cc*xn1(k)

       xn1(k) = xn(k)
       xn(k) = x(k)

      end do

      return
      end
       


      subroutine CHB_Iterate_CSI(N, xn, res, xn1, &
         ro_cheb, gamma_cheb)  

!=====================================================================*
!        A Chyebyshev Iteration for Solution                          *
! Hageman & Young "Applied Iterative Methods", Academic Press, 1981   *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
! Input:
      integer N
      real xn(N), res(N), xn1(N)
      real ro_cheb, gamma_cheb
! Output:
!      real x(N), xn(N), xn1(N)
! Local Variables
      integer k
      real x_tmp

!O            Flux_k1(n, k) = ro_cheb*
!O     &        (gamma*Flux_n(n, k) + Flux_k(n, k))+cc*Flux_k1(n, k)

      do k = 1, N

       xn1(k)=ro_cheb*( gamma_cheb*res(k) + xn(k) - xn1(k) ) + xn1(k)

       x_tmp = xn(k)
       xn(k) = xn1(k)
       xn1(k) = x_tmp

      end do

      return
      end
       
      
      subroutine CHB_Estimate_Res_CSA&
          (N, Res_old, Res_new, Eigenvalue, Res_Norm_2_CSA)
      implicit none
!Input:
      integer N
      real Res_old(N), Res_new(N), Eigenvalue
! Output:
      real Res_Norm_2_CSA
! Local
      real Res_CSA
      integer k

      Res_Norm_2_CSA = 0.

      do k = 1, N
         Res_CSA = Res_old(k)*Eigenvalue - Res_New(k)
         Res_Norm_2_CSA = Res_Norm_2_CSA + Res_CSA*Res_CSA
      end do

      Res_Norm_2_CSA = sqrt(Res_Norm_2_CSA)

      return
      end




      subroutine CHB_Iterate_Steady_State
!=====================================================================*
! Chebyshev Iterations for the Steady-State Calculations              *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! INPut: Flux(NG, N_TOT), Flux_N(NG, N_TOT), Flux_N1(NG, N_TOT),
!         ro_cheb, gamma
!      logical ADJOINT
! Output:  Flux(NG, N_TOT), Flux_N(NG, N_TOT), Flux_N1(NG, N_TOT),
!          delnp, d_flux_l
! Local Variables
      real Error_Flux_Norm_Inf, Error_Flux_Norm_2, Flux_Norm_Inf
      integer NN

! Chebyshev Extrapolation of the Neuron Flux

      NN = NG*N_TOT

      call MSC_Get_Abs_Error_Norm_Inf(&
             NN, Flux, Flux_N, Error_Flux_Norm_Inf)

      call MSC_Get_Abs_Error_Norm_2(&
             NN, Flux, Flux_N, Error_Flux_Norm_2)

      call CHB_Iterate_CSI_Residual(&
             NN, Flux, Flux_N, Flux_N1, ro_cheb, gamma )  

      call MSC_Get_Norm_Inf(&
             NN, Flux, Flux_Norm_Inf)

      delnp = Error_Flux_Norm_2/real(NN)
      d_flux_l = Error_Flux_Norm_Inf/Flux_Norm_Inf

      return
      end


      subroutine CHB_Iterate_Transient
!=====================================================================*
!  Computing Eigenvalue and Chebyshev Extrapolation                   *
!                 of Neutron Flux (Neutron Kinetics Calculations)     *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! INPUT: 
!  Flux(NG, N_TOT) - the residual (k+1)
!  Flux_N(NG, N_TOT) - Residual k
!  Flux_N1(NG, N_TOT) - Residual k-1
!  Flux_K(NG, N_TOT) - Neutron Flux k
!  Flux_K1(NG, N_TOT) - Neutron Flux k-1
!  ro_cheb, gamma - Chebyshev Iterative parameters
! Output: *  Flux(NG, N_TOT) - the residual (k+1)
!  Flux_N(NG, N_TOT) - Residual k
!  Flux_N1(NG, N_TOT) - Residual k-1
!  Flux_K(NG, N_TOT) - Neutron Flux k
!  Flux_K1(NG, N_TOT) - Neutron Flux k-1
! kin_k_ef - maximum eigenvalue of the iteration matrix (computed
!    using Rayleigh quotient)
! k_ef_min, k_ef_max - lower and upper bound for the maximum eigenvalue 
! d_kef_l - local error in eigenvalue = 0.5*(k_ef_max - k_ef_min)
!  deln_fp - the relative norm of the residual for the CSI
!  deln_dp - relative norm of the residual   for the CSA
!  delnp - relative norm of the residual which are used for
!      the Eigenvalue Estimate in the Chebyshev Iterative Proicedure
! Local Variables
      real  Flux_Norm, k_ef_new,  Rayleigh
!      real a_k, aa, bb, cc,&
!                 Flux_Tmp&
!                                      k_min_max
!      integer n, k       

      integer NN
      real Res_Norm_2, Flux_Norm_2, Res_Norm_2_CSA

      NN = NG*N_TOT

! Iteration of the Chebyshev Semi-Iterative Method
      call CHB_Iterate_CSI(NN, Flux_k, Flux_n, &
        Flux_k1, ro_cheb, gamma) 

! Estimate the 2-norm of the residual
      call MSC_DOT_PRODUCT(NN, Flux_n, Flux_n, Res_Norm_2)
     
      deln_fp  =       sqrt(Res_Norm_2)

! Get 2-norm of the flux
      call MSC_Get_Norm_2(NN, Flux_k, Flux_Norm_2)

      Flux_norm = Flux_Norm_2

! Get Residual dot product for Rayleigh eigenvalue estimate
      call MSC_DOT_PRODUCT(NN, Flux, Flux_n, Rayleigh)

! Get Minimum and Maximum Estimate of the Eigenvalues
!      call MSC_MAX_MIN_RATIO(NN, Flux, Flux_n, &
!                           k_ef_max, k_ef_min ) 

! local error of eigenvalue estimate        
!      d_kef_l = (k_ef_max - k_ef_min)*0.5

      k_ef_new  = Rayleigh / Res_Norm_2

      if(k_ef_new.LT.1) then ! It can happens sometimes
              kin_k_ef  = k_ef_new
      else
         if(Debug) then
             write(*,*) 'Maximum  Eigenvalue > 1, Time Step is Large'
             write(*,*) 'k_ef_max, k_ef_min, k_ef_new :' 
!                     0.5*(k_ef_max + k_ef_min):'
             write(*,*)  k_ef_max, k_ef_min, k_ef_new 
!                              0.5*(k_ef_max+k_ef_min)
          end if
      end if

      deln_fp = deln_fp/real(NN)
      Flux_Norm = Flux_Norm/real(NN)
      deln_fp = deln_fp / Flux_Norm

! calculation of the new residual vector  vector 
! estimate of the error if analytical summation can be done 


      call CHB_Estimate_Res_CSA(NN, Flux_n, Flux, kin_k_ef,&
          Res_Norm_2_CSA)

      call CHB_Iterate_CSI_Residual(NN, Flux, Flux_n, Flux_n1, &
         ro_cheb, gamma)  


      deln_dp = Res_Norm_2_CSA / real(NN)
      deln_dp = deln_dp / Flux_Norm

      if(iter_solver.eq."CSA") then
         delnp = deln_dp ! CSA
      else if(iter_solver.eq."CSI") then
         delnp = deln_fp ! CSI
      else 
       write(*,*) 'Iterative solver is not CHEBYSHEV type'
       stop     
      end if

 
      return 
      end


      subroutine CHB_Set_Iteration_Param(ip)
!=====================================================================*
! Comuting New Chebyshev Acceleration Parameters (For details         *
! see Algorithm 6.4.1 from Hageman and Young " Applied Iterative      *
!  Methods" Academic Press, 1981.                                     *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: ip - counter of the current iteration
! NPOLINS - Number of Iteration when Chebyshev Acceleration Starts
! Tau_Cheb - Upper Bound of Certain Eigenvalue of the Iteration Matrix
! Tau_CSA - Upper Bound of 2nd Eigenvalue (or dominance ratio)
!                                     of the Iteration Matrix
! Tau_CSI - Upper Bound of Maximum Eigenvalue of the Iteration Matrix
! D_Cheb -  strategy parameter used todetermine the minimum degree 
!            required for generated Chebyshev sequence 
!            (recommended in the range  [0.03, 0.15]
!  xme_ - current estimate of the maximum eigenvalue
!  xme - the estimate of the Maximum eigenvalue used for generated 
!        Chebyshev polynomilas
! delnp -  norm of the residual vector used in addaptive Chebyshev procedure
      integer ip

! Output: 
!  delnp0 - norm of the reidual vector at the previous iteration
!  npolin - counter of the Chebyshev polynomial currently 
!                  being used
! is_polin - counter of the chebyshev polynomial sequence
! gamma, ro_cheb - parameters of the Chebyshev Acceleration procedure
! r_cheb - theoretical  convergence rate
! npolin_min - the degree of the Chebyshev polynomial when the
!               maximum eigenvalue estimate can be changed

! Local Variables: 
      real Tau_Cheb, Tau_CSA, Tau_CSI
      parameter (Tau_CSA = 0.98, Tau_CSI = 0.995)

      real D_Cheb
      parameter(D_Cheb = 0.1)
      real sigme, a
      save sigme

      if(Problem_Type.NE."Kinetics") then 
         Tau_Cheb = Tau_CSA
      else ! KINETICS
         if(iter_solver .eq. "CSA") then
            Tau_Cheb = Tau_CSA
         else if(iter_solver .eq. "CSI") then
            Tau_Cheb = Tau_CSI
         else
! not used for CG type methods
         end if
      end if

      delnp0  = delnp

!      write(*,*) 'ip =', ip, 'npolins  =', npolins

      if(ip.lt.npolins) then
! simple iterations
         npolin = -1
         gamma = 1.
         ro_cheb = 1.
         xme = xme_
      else
          npolin = npolin + 1
!          write(*,*) 'ip =', ip, 'npolin  =', npolin
          if(npolin.eq.0) then
! initialization of parameters for a new polinomial

             is_polin = is_polin + 1
             if(xme_.gt.Tau_Cheb) xme_ = Tau_Cheb
             xme = xme_
             ro_cheb = 1.0
             gamma = 2. / (2. - xme - xbe)
             sigme = (xme - xbe)/(2. - xme - xbe)
             a = sqrt(1. - sigme*sigme)
             r_cheb = (1. - a)/(1. + a)
             npolin_min = alog(d_cheb)/alog(r_cheb)
             if(npolin_min.LT.4) npolin_min = 4

         else if(npolin.eq.1) then
              ro_cheb = 1./(1. - 0.5*sigme*sigme)
         else
              ro_cheb = 1./(1. - 0.25*sigme*sigme*ro_cheb)
         end if

      end if


      return
      end


      subroutine CHB_Estimate_Max_Eigenv(i_source)
!=====================================================================*
!     Estimation of the Maximum Eigenvalue (Or Dominance Ratio)       *
! of the Iteration Matrix see Algorithm 6.4.1 from Hageman and Young  *
! " Applied Iterative  Methods" Academic Press, 1981.                 *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: * F_Cheb - damping factor used in the parameter change test 
!        for addaptivity in te range [0.65-0.8] if = 0, Non-adaptive
!        procedure
!        rc_cheb - residual quotient at the previoius iteration
!        delnp - norm of the residual vector
!        delnp0 - norm of the residula vector at the previous iteration
!        npolin - degree of the Chebyshev polynomial
!        npolin_min - minimum degree of the Chebyshev polynomial when
!                     Maximum Eigenvalue Estimate can be changed
!        delnpi - norm of the residual vector at the 1st Chebyshev iteration
!        xbe - estimate of the minimum eigenvalue of the iteration matrix
!               in our case xbe = 0
!       r_cheb - theoretical  convergence rate
!       i_source - Number of Iterations (Used only for Output)
      integer i_source
! Output: Alsp_Cheb - quotient of the practical to theoretical Convergence 
!         rate
!         Csp_Cheb - quotient of log(practical)/log(theoretical) 
!         convergence rate
!         npolin = -1 if new eigenvalue estimate is accepted
!         xme_ - maximum eigenvalue estimate

! Local Variables:
! q_cheb - theoretical Convergence rate
! b_chen - practical convergence rate
! x_cheb 
      real q_cheb, b_cheb, x_cheb
! Local Variables: 
      real Tau_Cheb, Tau_CSA, Tau_CSI
      parameter (Tau_CSA = 0.98, Tau_CSI = 0.995)

      if(Problem_Type.NE."Kinetics") then 
         Tau_Cheb = Tau_CSA
      else ! KINETICS
         if(iter_solver.eq."CSA") then
            Tau_Cheb = Tau_CSA
         else if(iter_solver.eq."CSI") then
            Tau_Cheb = Tau_CSI
         end if
      end if

      rc_cheb0 = rc_cheb
      rc_cheb = delnp/delnp0          

      if(npolin .LE. 2) then
           
        if(npolin.eq.0) delnpi = delnp

      else

            Q_cheb = (2.*r_cheb**(npolin*0.5))/&
           (1. + r_cheb**(npolin))
            b_cheb = delnp/delnpi

            Alsp_Cheb = b_cheb/q_cheb
            Csp_Cheb = alog(b_cheb)/alog(q_cheb)

            if( b_cheb .LT. 1. ) then
! if the iterations converges (the decrease of the residuals)
                if(B_Cheb .GT. Q_cheb) then
! if theoretical convergence rate > the practical one we
! compute new eigenvalue estimate xme_
                   x_cheb = (0.5*(1.+r_cheb**(npolin))*&
                    (b_cheb + sqrt(b_cheb*b_cheb - q_cheb*q_cheb)))&
                     **(1./npolin)
                   xme_ = 0.5*(xme + xbe + (2. - xme - xbe)*&
                     (x_cheb*x_cheb + r_cheb)/((1.+r_cheb)*x_cheb))

!                 write(*,*) 'i_source =', i_source , 'XME_ =', xme_
                 
               else
! teoretical convegence rate <= practical => the same eigenvalue estimate
                  xme_ = xme
               end if !! B_Cheb > Q_Cheb
            else ! B_Cheb >= 1
! No convergence t_cheb is indicator that the iterations diverge
               T_Cheb = T_Cheb + 1
            end if !! B_CHEB < 1
      end if 

! Test of the new estimate change
      if(npolin.gt.npolin_min) then
! there are at least Npolin_Min (4) Chebyshev iterations shpould be performed
! before we change the eigenvalue estimate
            if(b_cheb.gt.q_cheb**F_Cheb) then
! theoretical convergence rate < than the practical one (F_Cheb = 0.65 - 
!  Safety parameter, if 0 Non Adaptive procedure) 
               if(T_Cheb.eq.0) then
! if there is a convegence
                  if(xme.eq.Tau_Cheb) then
! when the eigenvalue estimate is equal to the Maximum it is better to start
! from the beginner (it is not important really)
                     xme_ = xme_ini
                  end if
!npolin = -1 indicator that we start new Chebyshev sequence
                  npolin = -1
!                 write(*,*) 'NEW CHEBYSHEV SEQUENCE XME =',&
!                                    i_source, XME
              else  ! (T_CHEB /=0
!  it was no nonvergence
                if(rc_cheb.GT.1) then
! at the last iteration the residula is also increasing (New Chebyshev sequence
! with initial eigenvalue estimate)
                   npolin = -1
                   xme_ = xme_ini
!                    write(*,*) 'Convergence is not optimal',&
!                               'XME_ =', XME_
!                    pause
                else 
! if the last residual is decreasing we simply continue 
                   T_Cheb = 0
                end if
             end if ! T_Cheb = 0
            end if ! B_Cheb > Q_Cheb
      end if ! Npolin > NPolin_Min

!     write(*,*) 'xme_ =', xme_, 'xme =', xme

      if(DEBUG) call view_par(i_source)

      return
      end 


      subroutine view_par(i_source)
      implicit none
      include 'sketch.fh'
      integer i_source
!***********************************************************************
! output Chebyshev parameter into the files 'Output_Debug/domin.dat'   *
!  and 'Output_Debug/res_cheb.dat' (when Debug = .True.)               *
!***********************************************************************

      open(io_unit,file = 'Output_Debug/domin.dat',access='append',&
         status='unknown')
        write(io_unit,3) i_source, npolin,  kin_k_ef, xme, xme_
      close(io_unit)

      open(io_unit,file = 'Output_Debug/res_cheb.dat',access='append',&
         status='unknown')
        write(io_unit,1) i_source, npolin, delnp, d_flux_l, d_kef_l,&
                   1./((1-xme_)*(1.-kin_k_ef))

      close(io_unit)

 1    format(1x,2I4,4E12.5)
 3    format(1x,2I4,5E12.5)

      return
      end


! IT IS PART OF MY 'sketch.fh' where Chebyshev parameters are 
! saved 

! Module Chebyshev Parametrs
!         real xme, xme_, xbe, xme_ini, ro_cheb, gamma, r_cheb,&
!                        f_cheb, delnp,delnpi,delnp0,&
!                        rc_cheb, rc_cheb0, alsp_cheb,csp_cheb
!         integer         npolin, npolins, npolin_min, is_polin, t_cheb
!         common/chebysh/ xme, xme_, xbe, xme_ini, ro_cheb, gamma, &
!                        r_cheb,f_cheb, delnp,delnpi,delnp0,&
!                        rc_cheb, rc_cheb0, alsp_cheb,csp_cheb,&
!                        npolin, npolins, npolin_min, is_polin,&
!                        t_cheb
! End Module Chebyshev Parameters

      subroutine CHB_Extrapolate_CSA
!=====================================================================*
!  Lyusternik Extrapolation of the Neutron Flux for Neutron Kinetics  *
!                   Calculations                                      *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Flux_k(NG, N_TOT) - Neutron Flux (k)
!        Flux(NG_N_TOT) - Residual (teta in Zimin & Ninokata paper)
!        kin_k_ef - maximum eigenvalue of the iteration matrix
!        i_CSA - Iteration method (Extrapolation only for CSA)
!          if CSI we simply add the last residual

! Output: Flux(NG, N_TOT) - Neutron Flux 
!         Flux_k(NG, N_TOT) - Added Error 
! Local Variables:
      integer n, k
      real extrap , Flux_Tmp

      if(iter_solver.eq."CSA") then
        extrap = 1. / (1. - kin_k_ef)
      else if(iter_solver.eq."CSI") then
        extrap = 1.
      end if


         do k = 1, N_TOT
            do n = 1, NG
               Flux_Tmp = extrap*Flux(n, k)
               Flux(n, k) = Flux_k(n, k) + Flux_Tmp
               Flux_k(n, k) = Flux_tmp
            end do
         end do

       return
       end
