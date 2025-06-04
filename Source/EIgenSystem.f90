      subroutine EIS_Update_Wieland_Shift
!=====================================================================*
! Updating the Wieland Shift (Steady-State Calculations)              *
!                   Vyachreslav Zimin (c) May 20 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      eigenv_shift = delta_shift/(1. + k_ef*delta_shift)

      return
      end

      subroutine EIS_Compute_Eigenvalue(ADJOINT)
!=====================================================================*
! calculation new source term, norm of the source,  eigenvalue, k_eff *
!   and error in k_eff   *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! INPut: Flux(NG, N_TOT), Source(N_TOT)
      logical ADJOINT
! Output: Source(N_TOT),  eigenv_min, eigenv_max
!          d_kef_l, S_Norm, eigenv
! Local Variables
      real S_Norm_Old,  eigenv_max, eigenv_min

! Computing New Fission Source & Eigenvalue

      S_Norm_Old  = S_Norm

      call MSC_SCOPY(N_TOT, Source, Source_dt)

      call EIS_Compute_Source(ADJOINT)

      call MSC_Get_Norm_1(N_TOT, Source, S_Norm)

      call MSC_MAX_MIN_RATIO(N_TOT, Source, Source_dt, &
                           eigenv_max, eigenv_min )
       
      eigenv = eigenv*S_Norm/S_Norm_Old

      k_ef = eigenv / (1. + eigenv*eigenv_shift)
      k_ef_min = eigenv_min / (1. + eigenv_min*eigenv_shift)
      k_ef_max = eigenv_max / (1. + eigenv_max*eigenv_shift)

      d_kef_l = 0.5*(k_ef_max - k_ef_min)

      return
      end


      subroutine EIS_Compute_Source(ADJOINT)
!=====================================================================*
!        Computing the Source Term                                    *
!         Slava (c) 23.II.1998                                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input:  XS_SF(NG<N_TOT), FLUX(NG,N_TOT), XP(NG)
      Logical ADJOINT
! Output: Source(N_TOT)

! Local Variables
      
      integer n, k

      if(ADJOINT) then

           do k = 1, N_TOT
              Source(k) = 0.
              do n = 1,NG
                Source(k) = Source(k) + xp(n)*Flux(n,k)
              end do
           end do

      else

           do k = 1, N_TOT
              Source(k) = 0.
              do n = 1,NG
                Source(k) = Source(k) + XS_SF(n,k)*Flux(n,k)
              end do
           end do

      end if

      return
      end

      subroutine EIS_Find_Critical_Boron_Conc
!=====================================================================*
!    Computing New Boron Concentration for the Criticality Search     *
!         Slava (c) 23.II.1998                                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: 
      integer k , n1 
      real Minimum_Difference
      parameter (Minimum_Difference = 1.E-06)

      real x, x_old, f, f_old, df, df_dx, dx
      save x_old, f_old, x, f
      real dx_1
      parameter (dx_1 = 200.)

      f = k_ef - 1



      IF( abs(f) .GT. e_critical ) THEN

        df = f - f_old

        IF( abs(df) .GT. Minimum_Difference ) then

          IF(I_Bor_Start.eq.1) THEN

             x_old = fdback(1,1,1)
             IF( f .GT. 0) then
                x = x_old + dx_1
             ELSE
                x = AMAX1(x_old - dx_1, 0.)
             END IF
             do n1 = 1, NZ
                do k = 1, NH
!                boron_prev(k,n1)  = fdback(k,n1,1)
                fdback(k,n1,1) = x
!                  fdback(k,n1,1) + delta_boron_1*(k_ef-1.)
                end do
             end do
             I_Bor_Start = 0
          ELSE
            x = fdback(1,1,1)
            dx = x - x_old
            IF( ABS(dx) .GT. Minimum_Difference) THEN
               df_dx = (f - f_old)/dx
               x_old = x
               f_old = f
                 x = AMAX1( x - f/df_dx, 0.)
!                 WRITE(*,*) 'Derivative =', df_dx, 'f =', f,&
!                'x =', x
!                 PAUSE
               do n1 = 1, NZ
                  do k = 1, NH
!               boron_tmp  = fdback(k,n1,1)
                     fdback(k,n1,1)= x
!               (fdback(k,n1,1)-boron_prev(k,n1))*&
!                 (1.-k_ef)/(k_ef - k_ef_prev) + fdback(k,n1,1)
!               boron_prev(k,n1) = boron_tmp
                  end do
               end do
            END IF ! dx > e_round_off
          END IF ! I_Bor_Start /= 1

!          k_ef_prev = k_ef

        END IF ! df > e_round_off

      END IF ! f /= 0


      return 
      end

      subroutine EIS_Find_Critical_Boron_Conc_NEW(completed)
!=====================================================================*
!    Computing New Boron Concentration for the Criticality Search     *
!         Slava (c) 23.II.1998                                        *
!=====================================================================*
      implicit none
        include 'sketch.fh'
! Input: 
!        CHARACTER*2  name_crit_parameter
      integer k , n1 
      real*8 Minimum_Difference
      parameter (Minimum_Difference = 5.E-06)
!        REAL, PARAMETER :: DX_MIN_FORMAT = 1.1E-04

      real*8 x, x_old, f, f_old, df, df_dx, dx_rel, dx, dx_new
        save x_old, f_old, x, f, df_dx
      real*8 dx_1
      parameter (dx_1 = 200.)
      LOGICAL  completed, first_iteration
!        real     e_critical
        REAL, PARAMETER :: xmin =0.0, xmax=2000.

      completed = .False.
! reading the previous iteration values
      IF( I_Bor_Start == 1) THEN
        first_iteration = .True.
          I_Bor_Start = 0
      ELSE
        first_iteration = .False.
      END IF         
! reading BUCKL
      f = k_ef - 1
      x = fdback(1,1,1)
!
      IF( x <= xmin .AND. f < 0 )  THEN
          completed = .True.
          x_old = x
          f_old = f
          GO TO 10
      END IF

      IF( x >= xmax .AND. f > 0 )  THEN
          completed = .True.
          x_old = x
          f_old = f
          GO TO 10
      END IF

!      write(*,*) 'f=', f, 'e_critical =', e_critical
!        write(*,*) 'completed =' , completed

      IF( abs(f) .GT. e_critical ) THEN

          df = f - f_old
          dx = x - x_old
          dx_rel = 1. - x_old/x

!          write(*,*) 'x=', x, 'x_old=', x_old
!          write(*,*) 'df=', df, 'dx_rel=', dx_rel

!        write(*,*) 'first_iteration =', first_iteration

          IF(first_iteration) THEN

              x_old = x
              f_old = f
            IF( f .GT. 0) then
                x = x_old + dx_1
              x = MAX( x, xmax )
              ELSE
               x = MAX(x_old - dx_1, xmin)
              x = MIN( x, xmax )
              END IF

          ELSE ! /= first iteration

            x_old = x
            f_old = f

!            write(*,*) ' ABS(dx_rel), ABS(df), TOL =',&
!                ABS(dx_rel), ABS(df), Minimum_Difference
            IF( ABS(dx_rel) .GT. Minimum_Difference.AND.&
            ABS(df) .GT. Minimum_Difference ) THEN
! update df_dx
              df_dx = df/dx
              write(*,*) 'updated df_dx=', df_dx
            ELSE
                  
              write(*,*) 'used old df_dx=', df_dx
 
            END IF ! dx > e_round_off
                IF ( abs(df_dx) > Minimum_Difference ) THEN
                   dx_new =  -f/df_dx
                ELSE
                IF( f .GT. 0) then
                   dx_new =  dx_1
                 ELSE
                   dx_new =  -dx_1
                 END IF
                END IF
            x = MAX( x + dx_new, xmin)
          x = MIN( x, xmax )
            write(*,*) 'x_old, x_new =', x_old, x
         END IF ! first_iteration
      ELSE   ! 
         x_old = x
         f_old = f
         completed = .True.
      END IF ! f /= 0

   10 CONTINUE

!              pause 


               do n1 = 1, NZ
                  do k = 1, NH
                     fdback(k,n1,1)= x
                  end do
               end do
              
      RETURN
      END subroutine EIS_Find_Critical_Boron_Conc_NEW
