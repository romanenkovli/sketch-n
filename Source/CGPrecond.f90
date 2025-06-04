      subroutine block_ssor_solve(y, x)
!=====================================================================!
!   Precponditioner solver of the system M*x=y,                       !
!               with M = (D+L)D^-1(D+U)                               !
!        matrix is stored in the internal sketch format               !
!   NG*N_TOT - size of the matrix                                     !
!   MAT_Block_Diag_Inv(NG,NG,N_TOT) - inverse block diagonal D^-1     !
!   MAT_Tot(NG,NE_T,N_TOT):                                           !
!     MAT_Tot(NG,1:3,N_TOT) - (-L)                                    !
!     MAT_Tot(NG,1:3,N_TOT) - (-U)                                    !
!     Neib(1:6,N_TOT) - indexing array                                !
! NG, N_TOT, Mat_Block_Diag_Inv, Mat_Tot, Neib are passed through     !
!  common block                                                       !
! the subroutine can be in place                                      !
!=====================================================================!
      implicit none    
      include "sketch.fh"
      
      real x(NG, N_TOT), y(NG, N_TOT)

!     working array 
      real wrk(NG)

      integer n, i, next, m, k

! forward solve
      do k = 1, N_TOT
         do n = 1, NG
            wrk(n) = y(n,k)
            do i = 1, 3
               next = Neib_REP_STR(i,k)
               if(next.ne.I_BOUND_NODE.and.next.GT.0) then
               wrk(n) =  wrk(n) + MAT_Tot(n,i,k)*x(n,next)
               end if
            end do
        end do ! NG
        do n = 1, NG
           x(n,k) = 0.
           do m = 1, NG
              x(n,k) = x(n,k) + MAT_Block_Diag_Inv(n, m, k)*wrk(m)
           end do 
        end do ! n
      end do ! k

! backward solve

      do k = N_TOT, 1, -1
         do n = 1, NG
            wrk(n) = 0.
            do i = 4, NE_T
               next = Neib_REP_STR(i,k)
               if(next.ne.I_BOUND_NODE.and.next.GT.0) then
               wrk(n) =  wrk(n) + MAT_Tot(n,i,k)*x(n,next)
               end if
            end do
        end do ! NG
        do n = 1, NG
           do m = 1, NG
              x(n,k) = x(n,k) + MAT_Block_Diag_Inv(n, m, k)*wrk(m)
           end do 
        end do ! n
      end do ! k

      return
      end

      subroutine AMUX_SKETCH(x, y)
!=====================================================================*
!  matrix vector product in the sketch  format   y = A x              !
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include "sketch.fh"

      real x(NG, N_TOT), y(NG, N_TOT)

      integer k, n, m, next, i

      do k = 1, N_TOT
        do n = 1, NG
! NONDIAGONAL ELEMENTS

         y(n,k) = 0.0

! BLOCK DIAGONAL ELEMENTS
         do m = 1, NG
            y(n,k) = y(n,k) + MAT_Block_Diag(n,m,k)*x(m,k)
         end do ! m

         do i = 1, NE_T

           next  = Neib_REP_STR(i,k)
           if(next.ne.I_BOUND_NODE.and.next.GT.0) then
             y(n,k)  = y(n,k) - x(n,next)*MAT_Tot(n, i, k)
           end if

         end do


       end do
      end do



      return
      end


      subroutine CGSolver_RHS_set

      implicit none    
      include "sketch.fh"

      integer i, k, n

      i = 0
      
      do k = 1, N_TOT
         do n = 1, NG
            i = i + 1
            rhs(i) = Mat_rhs_k(n,k)
            xran(i) = flux(n,k)
         end do
      end do

      return
      end


      subroutine CGSolver_Flux_Set

      implicit none    
      include "sketch.fh"

      integer i, k, n
      
      i = 0
 
      do k = 1, N_TOT
         do n = 1, NG
            i = i + 1
            flux(n,k) = sol(i)
         end do
      end do

      return
      end

      subroutine CGS_compute_scaled_residual(NN, x, f, res_norm)
! Input:
      integer NN
      real x(NN), f(NN)
! external function : dnrm2
      real dnrm2
! Output:
      real res_norm
! Local Variables
      real res(NN)
      integer n

! new convergence criteria
!       res = A*x
      call AMUX_SKETCH( x, res )
      do n = 1, NN
         res(n) = f(n) - res(n) 
      end do
      res_norm = dnrm2(NN, res, 1)/ dnrm2(n, f , 1)

      return
      end
