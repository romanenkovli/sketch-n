      subroutine PCD_Block_SSOR(i_ssor,E_INT,N_INT)
!=====================================================================*
!        Block Symmetric Gauss-Seidel Iterative Solver                *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Input: Flux(NG,N_TOT),  
! real MAT_Block_Diag_Inv(NG,NG,N_TOT),
!        Neib_REP_STR(NE_T,N_TOT), Block_RHS(NG,N_TOT) 
! Module MATrix:
!     real MAT_Tot(NG,NE_T,N_TOT)
      integer n,  N_INT
      real E_INT
! Output: i_int, Flux(NG, N_TOT), 
      integer i_ssor
! Local Variables:
      real RHS_SSOR(NG), Flux_Old(NG,N_TOT), Delta_Flux, Delta
      integer k, i, next, m

      i_ssor = 0
      Delta_Flux = 1.

! Symmetric Successive OverRelaxation method  (start)
      do while((Delta_Flux .GT. E_INT).AND. (i_ssor.LT. N_INT))

       i_ssor = i_ssor + 1
           
       Delta_Flux = 0.

! Direct Passage
       do k = 1, N_TOT

         do n = 1, NG
            RHS_SSOR(n) = 0.
         end do

         do i = 1, NE_T
           next  = Neib_REP_STR(i,k)
!           if(next.ne.I_BOUND_NODE) then
           if(next.GT.0) then
           do n = 1, NG
             RHS_SSOR(n) = RHS_SSOR(n) + &
                    Flux(n,next) * MAT_Tot(n, i, k)
           end do ! NG
           end if
         end do

         do n = 1, NG
            Flux_Old(n, k) = Flux(n, k)
            RHS_SSOR(n) = RHS_SSOR(n) + MAT_Block_RHS(n,k)
         end do

         do n = 1, NG   
            Flux(n,k) = 0.
            do m = 1, NG
              Flux(n, k)= Flux(n,k) + MAT_Block_Diag_Inv(n,m,k)*&
                                   RHS_SSOR(m)
            end do ! m
         end do ! NG

       end do
 
! Reverse Passage

      do k =  N_TOT, 1, -1 

         do n = 1, NG
            RHS_SSOR(n) = 0.
         end do

         do i = 1, NE_T
            next  = Neib_REP_STR(i,k)
!            if(next.ne.I_BOUND_NODE) then
            if(next.GT.0) then
            do n = 1, NG
               RHS_SSOR(n) = RHS_SSOR(n) + &
                    Flux(n,next) * MAT_Tot(n, i, k)
            end do ! NG
            end if
         end do

         do n = 1, NG
            RHS_SSOR(n) = RHS_SSOR(n) + MAT_Block_RHS(n,k)
         end do

         do n = 1, NG   
            Flux(n,k) = 0.
            do m = 1, NG
              Flux(n, k)= Flux(n,k) + MAT_Block_Diag_Inv(n,m,k)&
                                                      *RHS_SSOR(m)
            end do ! m
            Delta = Flux(n,k) - Flux_Old(n,k)
            Delta_Flux = AMAX1(Delta_Flux,abs(Delta/Flux(n,k)))
         end do ! NG

      end do ! Reverse Passage

      end do ! SSOR Iterations

!     call PCD_DBG_Output_Block_SSOR(Delta_Flux, i_ssor)
!     pause 'Debug SSOR'

      return
      end

      subroutine PCD_DBG_Output_Block_SSOR_Output(Delta_Flux, i_int)
!=====================================================================*
!        Debug Output of Neutron Flux into the file                   *
!          'Output_Debug/SSOR.dat'                                    * 
!                   Vyachreslav Zimin (c) 19 May 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'

      integer n, k, i_int
      real Delta_Flux

      n = 2

      open(io_unit,file ='Output_Debug/Block_SSOR.dat', &
        status ='unknown', access = 'Append')

      write(io_unit,*) 'Neutron Energy Group =', n, 'N_Iter =', i_int
      write(io_unit,*) 'Delta Flux =', Delta_flux
      write(io_unit,1) (Flux(n,k), k = 1, NH)

    1 Format(5(1x,6E12.5/),(1x,5E12.5/))

      close(io_unit)

      return
      end


      subroutine PCD_Block_SSOR_NEW(i_ssor,E_INT,N_INT)
!=====================================================================*
!        Block Symmetric Gauss-Seidel Iterative Solver                *
!  implemented as a Preconditioner (residual & SSOR solver )          *         
!                   Vyachreslav Zimin (c) May 8 2000                  *
! FINISHED WORK SLOWLY THAN THE PREVIOUS VERSION                      *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Input: Flux(NG,N_TOT),  
! real MAT_Block_Diag_Inv(NG,NG,N_TOT),
!        Neib_REP_STR(NE_T,N_TOT), Block_RHS(NG,N_TOT) 
! Module MATrix:
!     real MAT_Tot(NG,NE_T,N_TOT)
      integer n,  N_INT
      real E_INT
! Output: i_int, Flux(NG, N_TOT), 
      integer i_ssor
! Local Variables:
      real Delta_Flux
      real res(NG, N_TOT), delta(NG, N_TOT)
      integer k

      i_ssor = 0
      Delta_Flux = 1.

! Symmetric Successive OverRelaxation method  (start)
      do while((Delta_Flux .GT. E_INT).AND. (i_ssor.LT. N_INT))

       i_ssor = i_ssor + 1
       Delta_Flux = 0.

! computed residual r = b - A*x

       call AMUX_SKETCH(flux, res) ! res = A*x
!       write(*,*) 'flux =', flux(:,1), 'res(:,1) =', res(:,1)
!       pause

       do k = 1, N_TOT
          do n = 1, NG
             res(n,k) = MAT_Block_RHS(n,k) - res(n,k)
          end do
       end do


! solving M*d = r                                    
       call block_ssor_solve(res, delta)

! new solution x = x + d

       do k = 1, N_TOT
          do n = 1, NG
             flux(n,k) = flux(n,k) + delta(n,k)
          end do
       end do

     
      end do

      return
      end 

