      subroutine  XSP_Compute_Deln_Par
!=====================================================================*
!                Computing Delayed Neutron Parameters                 *
!                      Debug Procedure for Kinetics                   *
! (c) Slava 15.IV.1998  vgzimin@mail.ru                 *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      real fiss, fiss_d(MD)
      INTEGER :: n, m, kt
 

      IF( I_flag_deln_Library == 1 ) THEN

         fiss = 0.
         fiss_d(1:MD) = 0.
         DO kt = 1, N_TOT
          DO n = 1, NG
            fiss = fiss + XS_SF(n,kt)*Flux(n, kt)
          DO m = 1, MD
            fiss_d(m) = fiss_d(m) + XS_SF(n,kt)*Flux(n, kt)*&
                                          xs_beta(m,kt)
          END DO 

          END DO
         END DO  
         
        DO m = 1, MD
         beta(m) = 0. 
         alfa(m) = 0.
         DO kt = 1, N_TOT
           DO  n = 1, NG
           beta(m) =beta(m)+xs_beta(m,kt)*XS_SF(n,kt)*Flux(n,kt)

           IF ( xs_alfa(m,kt) .GT. SMALL_VALUE ) THEN 
           alfa(m) =alfa(m)+xs_beta(m,kt)*XS_SF(n,kt)*Flux(n,kt)/&
                                                       xs_alfa(m,kt)
!           write(*,*)  'xs_alfa(m,kt) =', m, xs_alfa(m,kt)
           END IF


           END DO
         END DO
!         PAUSE 
       END DO 

         DO m = 1, MD 
         beta(m) = beta(m)/fiss
         alfa(m) = alfa(m)/fiss_d(m)
         alfa(m) = 1./alfa(m)
         END DO

         bet = SUM( beta(1:MD) )



!      write(*,*) 'alfa =', alfa(1:MD)

!      write(*,*) 'beta =', beta(1:MD)
         
      write(*,*) 'bet =', bet
        
!        pause

       

      END IF !       IF( I_flag_deln_Library == 1 ) THEN

      RETURN 
      END 
 



      subroutine PNT_Compute_Reactivity(dt_kin)
!=====================================================================*
!                Computing Reactivity  from the Neutron Balance       *
!                          weighting with ADJOINT Flux                *
!                      Debug Procedure for Kinetics                   *
! (c) Slava 15.IV.1998  vgzimin@mail.ru                 *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: Diag_Tot(NG, N_TOT), al(NG), volume(N_TOT), dt_kin, Flux,
!        XS_SIK, XS_SF, Matr_Tot(n, i, k)

! Module MATrix: 
!   real  MAT_Diag_Tot(NG, N_TOT) 
!   real  MAT_Tot(NG, NE_T, N_TOT) 
      real dt_kin
! Output: alp - Prompt Neutron Life Time
! ro - reactivity (absolute value)
! local Variables
      real vol, ff, Sum_Non_Diag , fiss, scat, k_ef_norm
      integer k, n , i,  m , next

      double precision Point_Tol, beta_p(MD), alfa_p(MD), dt_point, &
                           bet_p
      parameter (Point_Tol = 1.E-05)
      logical Extrapolation 
      integer N_Point_Step 

      react_dt = react

      if(Problem_Type .EQ. "Kinetics" ) then
         k_ef_norm = 1.
      else
         k_ef_norm = 1.
!         k_ef_norm =  1./k_ef
      end if

         react = 0.
         al_prompt = 0.

         ff = 0.

       do k = 1, N_TOT

          vol = volume(k)

! absorption + leakagec out
          do n = 1, NG

!            write(*,*) 'Flux_A=', Flux_A(:,:)

           al_prompt=al_prompt+Flux(n, k)*xs_al(n,k)*Flux_A(n, k)*vol

           react = react - MAT_Diag_tot(n, k)*&
               Flux(n, k)*Flux_A(n, k)

! scattering + multiplication
            scat = 0.
            fiss = 0.
            do m = 1,NG
             fiss = fiss + XS_SF(m,k)*Flux(m, k)
             scat  = scat + XS_SIK(n, m, k)*Flux(m, k)
            end do

            fiss = fiss*k_ef_norm

            react = react + (scat + xp(n)*fiss)*Flux_A(n, k)

            ff = ff + xp(n)*fiss*Flux_A(n, k)

! leakage in
             Sum_Non_Diag = 0.
             do i = 1, NE_T
                 next  = Neib(i,k)
                 if(next.ne.I_BOUND_NODE) then
                 Sum_Non_Diag = Sum_Non_Diag + &
                    Flux(n,next) * MAT_Tot(n, i, k)
                 end if
             end do

             react = react + Sum_Non_Diag*Flux_A(n, k)
             end do
          end do

      react = react/ff
      al_prompt = al_prompt/ff

!      write(*,*) 'react =', react, 'bet = ', bet
!      read(*,*)

      write(*,'(" reactivity =",F9.5," ($)"," al=", E12.5," (s)")')&
             react/bet, al_prompt
!      write(*,*) 'flux_a =', flux_a(1,1:10)
!     pause

      do m = 1, MD
         alfa_p(m) = alfa(m)
         beta_p(m) = beta(m)
      end do

      dt_point = dt_kin
      bet_p = bet

      Pow_Point_dt = Pow_Point

      if(Kinetics_Method.EQ."PNT") then
! Point Kinetics Calculations
! Fully-implicit scheme with external time step control 
          extrapolation = .False. ! 
          N_Point_Step = 1        ! calculation of the one time step
!          write(*,*) 'Pow_Point =', Pow_Point, 'al_prompt=', &
!           'al_prompt=', al_prompt
!          write(*,*) 'bet =', bet_p,  'beta=', beta_p(:)
!          write(*,*)  'alfa_p=', alfa_p(:)
!          write(*,*)  'Prec_point=', prec_point(:)
          
          call PNT_Solver(MD,  dt_point, Pow_Point, al_prompt, &
                beta_p, bet_p, alfa_p, Prec_Point, react, react_dt, &
                Point_Tol, Point_Deriv, extrapolation, N_Point_Step)
!          write(*,*) 'Pow_Point =', Pow_Point, 'al_prompt=', &
!           'al_prompt=', al_prompt
!      read(*,*)

      else if (Kinetics_Method .EQ. "IQS") then
! Point Kinetics Calculations
! Fully-implicit scheme with internal Stoer-Burlich extrapolation to get 
! accurate solution of the point kinetics equations
          extrapolation = .True.
          N_Point_Step = 1 ! NOT USED
          call PNT_Solver(MD,  dt_point, Pow_Point, al_prompt, &
                beta_p, bet_p, alfa_p, Prec_Point, react, react_dt, &
                Point_Tol, Point_Deriv, extrapolation, N_Point_Step)
      else if (Kinetics_Method .EQ. "DRT") then
           Point_Ampl_R = 1.
           Point_Deriv = 0.
      else
          write(*,'(A)') ' Neutron Kinetics Method is undefined'
          write(*,'(A)') ' It should be DRT, IQS or PNT        '
          write(*,'(A)') ' Please Check Namelist File SKETCH.INI'
          stop
      end if

      return
      end

      subroutine PNT_Solver(MD,  dt_sketch, power, gtime, beta, bet, &
        lambda, c, react_new, react_old, Tolerance, Point_Deriv, &
        extrapolation, N_Point_Step)
!=====================================================================*
!                   Point Kinetics Solver                             *
!=====================================================================*
      implicit none
! Input:
      integer MD, MDMAX
      parameter (MDMAX = 6)
      double precision dt_sketch,  beta(MD), bet, lambda(MD), Tolerance
      double precision  power, gtime, react_new, &
             c(MD),  react_old
      logical Extrapolation
      integer N_Point_Step
!      parameter (Extrapolation = .True.)
! Local Input
      double precision ex1_point(MDMAX), ex_point(MDMAX),  C_Old(MDMAX), &
                   C_Coarse(MDMAX)
! Output:  
! double precision Power, C(MD), Power_Der
!        dt_sketch - proposed value of the time step size (double precisionly 
      integer j, m
      double precision Point_Deriv
! Local Variables
      Logical  Step_Accept
      integer i_time, j_step
      integer NMAX, KMAXX
      parameter (NMAX = 7, KMAXX = 8)
      integer  N_Step_Fine(KMAXX), N_STEP_Burilsch(KMAXX)
! Stoer-Burlich Sequence
      data N_STEP_Burilsch /  2, 4,  6, 8,  12, 16, 26, 32/
!      data N_Point_Step / 20 /

! Saving for Restart
      double precision Pow_Old
! Saving for Error Control
      double precision Pow_Coarse
!     
      double precision dt_step
      double precision TE_Power, TE_Prec, Err_Max

      double precision react, delta_react
      double precision y_est(NMAX), y_extr(NMAX), &
               y_err(NMAX)

      Delta_React = React_New - React_Old

      i_time = 0
      j_step = 0
      step_accept = .False.
      dt_step = dt_sketch
     
      if(Extrapolation) then
        do j = 1, KMAXX
           N_STEP_Fine(j) = N_STEP_Burilsch(j)
        end do
      else 
! Constant Time Step for Point Kinetics
          N_STEP_Fine(1) = N_Point_STEP
      end if

! Saving for Restart
      Pow_old = Power
      call MSC_DReplace_Array(MD, C, C_Old )

!     write(*,*) 'react_old, react_new =', react_old, react_new

      do while(.NOT. Step_Accept .AND. i_time.lt.KMAXX ) 

           i_time = i_time + 1

           dt_step = dt_sketch/N_Step_Fine(i_time)
!           write(*,*) i_time, 'dt = ', dt_step
!           write(*,*) 'N_Step_Fine(i_time) = ', N_Step_Fine(i_time)
!           pause

           do j = 1, N_STEP_FINE(i_time)

! Linera Interpolation of the Reactivity
!                write(*,*) 'Step Begins, j/ N_Step_Fine(j)',&
!                      dble(j)/dble(N_Step_Fine(i_time))
                react = react_old + delta_react*&
                     dble(j)/dble(N_Step_Fine(i_time))
!               write(*,*) 'react =', react, 'j =', j
 
                call PNT_Compute_Time_Step(power, gtime, beta, lambda,&
                         c,dt_step, MD, react, ex1_point, ex_point)
                j_step = j_step + 1
 

          end do

!         write(*,*) 'Steps Finished '
!         pause

          if(Extrapolation) then

            y_est(1) = power
            do m = 1, MD
               y_est(m+1) = c(m)
            end do

            call PNT_Rational_Extrapolation&
                   (i_time, dt_step, y_est, y_extr, y_err, MD + 1)

            if(i_time.ne.1) then
               
               TE_Power = dabs(y_err(1)/y_extr(1))
               TE_Prec = 0.
               do m = 2, MD + 1
                  Te_Prec = dmax1(Te_Prec, dabs(y_err(m)/y_extr(m)))
               end do

               Err_max =dmax1(Te_Power, TE_Prec)

               Err_Max = Err_max / Tolerance

               if(Err_max.LT. 1 .OR.  i_time.eq.KMAXX ) then
                 Step_Accept = .True.
! gives the Extrapolated Values to the Power
                 Power = y_extr(1)
                 do m = 1, MD
                   C(m) = y_extr(m+1)
                 end do
               end if
           end if

        else
! Constant Time Step 
             Step_Accept = .True. 
             TE_Power = 1.
             TE_Prec = 1.
        end if

       if(Step_Accept.OR.i_time.EQ.KMAXX) then
!          write(*,*) 'Point Solved', 'Step_accept=', Step_Accept,&
!                  'i_time =', i_time
! Time Step is Accepted
              call Comp_Point_Deriv(MD, Power, &
                   react, bet, gtime, lambda, c,&
                              Point_Deriv)         
       else
! Saving Current Values for the Error Control
              Pow_Coarse = Power
              call MSC_DReplace_Array(MD, C, C_Coarse )

! Restart with the Old Values
              Power = Pow_Old
              call MSC_DReplace_Array(MD, C_Old, C)

          end if


        end do

      write(*,*) 'Time Step of Point Kinetics  = ', &
                dt_step, N_Step_fine(i_time)
!      write(*,*) 'Error of the Power, Precursors =', TE_Power, TE_Prec
!     pause

      return
      end


      subroutine PNT_Compute_Time_Step(power, gtime, beta, lambda, c, &
               dt, MD, react, ex1_point, ex_point )
!=====================================================================*
! point neutron kinetics time step                                    *
! input: power, gtime, beta, lambda,c(MD),MD,void,void_react,ro_ext,  *
!             dt,MD                                                   *
! output: power, c(MD), react                                         *
!=====================================================================*
      implicit none

! Input:
      integer MD
      double precision dt, beta(MD), lambda(MD)
      double precision  power, gtime, react,  c(MD),&
                   ex1_point(MD), ex_point(MD)
! Output:
!     double precision power, c(MD)
! Local Variables
      integer m
      double precision  sum_c, beta_up, beta_down, &
             power_old, argument
      double precision dt_relative

!     write(*,*) 'Solver In '


!     write(*, *) 'C =', c

      dt_relative = dt/gtime

      sum_c = 0.
      beta_up = 0.
      beta_down = 0.
      do m = 1, MD
          argument = - lambda(m)*dt
          ex_point(m) = dexp(argument)
          sum_c = sum_c + lambda(m)*c(m)*ex_point(m)
          ex1_point(m) = (1. - ex_point(m))/(-argument)
          beta_up = beta_up + beta(m)*(ex1_point(m)-ex_point(m))
          beta_down = beta_down + beta(m)*ex1_point(m)
      end do

!      write(*,*) 'react =', react
!      write(*,*) 'dt_relative*(react - beta_down)=', dt_relative*&
!                (react - beta_down)
!      write(*,*) 'react - beta_down =', react - beta_down


!     write(*,*) 'beta_up =', beta_up, 'beta_down =', beta_down
!     write(*,*) 'sumc_up =', sum_c
!     write(*,*) 'ex1 =', ex1_point
!     write(*,*) 'ex_point =', ex_point
!     write(*,*) '(1.  - dt_relative * (react - beta_down))  =',&
!               (1.  - dt_relative * (react - beta_down))  
!     pause
 
      power_old = power
      power = (power_old*(1. + dt_relative * beta_up ) + &
                                     dt_relative * sum_c)/&
                       (1.  - dt_relative * (react - beta_down))  

      do m  = 1, MD
         c(m) = c(m)*ex_point(m) + beta(m)*&
            (power_old*(ex1_point(m)-ex_point(m))+&
            power*(1. - ex1_point(m)))/lambda(m)
      end do

!     write(*,*) 'Inside Solver power, react =', power, react
!     pause
!     write(*,*) 'Solver Out '

      return
      end


      SUBROUTINE PNT_Rational_Extrapolation(iest,xest,yest,yz,dy,nv)
!=====================================================================*
!         Rational Extrapolation in the Burlich-Stoer Method          *
! Taken from Numerical Recipes (c) 1992 Cambridge University Press    *
!         Section 16.4 page 725                                       *
!=====================================================================*

! Input: iest - number of calls
!       NV - dimension of the extrapolated vector
!       xest - h/nseq(k)**order 
!       yseq(NV) - vector to extrapolat
! Output:
!        yz(NV) - extraplated function value
!        dy(NV) their estimate error        
      INTEGER iest,nv,IMAX,NMAX
      double precision xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=7)
!     Exact substitute for pzextr, but uses diagonal rational function 
!     extrapolation instead of polynomial extrapolation.
      INTEGER j,k
      double precision b,b1,c,ddy,v,yy,d(NMAX,IMAX),fx(IMAX),x(IMAX)
      SAVE d,x

      x(iest)=xest ! Save current independent variable.
      if(iest.eq.1) then
        do 11 j=1,nv
           yz(j)=yest(j)
           d(j,1)=yest(j)
           dy(j)=yest(j)
 11     end do
      else
        do 12 k=1,iest-1
           fx(k+1)=x(iest-k)/xest
 12     end do 
        do 14 j=1,nv ! Evaluate next diagonal in tableau.
         yy=yest(j)
         v=d(j,1)
         c=yy
         d(j,1)=yy
         do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.) then
               b=(c-v)/b
               ddy=c*b
               c=b1*b
            else  ! Care needed to avoid division by 0.
               ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
               d(j,k)=ddy
               yy=yy+ddy
 13      enddo 
         dy(j)=ddy
         yz(j)=yy
 14    end do
      end if

      return
      END

      subroutine Comp_Point_Deriv(MD, Power, &
                   react, bet, gtime, alfa, c,&
                              Point_Deriv)
      implicit none
! input
      integer MD
      double precision Power, react, bet, gtime, alfa(MD), c(MD)
! output 
      double precision Point_Deriv
! local variables
      integer m
      double precision sumc

      sumc = 0.
      do m = 1, MD
        sumc = sumc + alfa(m)*c(m)
      end do

      Point_Deriv = (react - bet)/gtime +  sumc/(Power*gtime)

!     write(*,*) 'Point Deriv = ', Point_Deriv
!     pause

      return
      end
