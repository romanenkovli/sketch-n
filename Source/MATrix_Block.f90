      subroutine MAT_Set_Kin_Block_RHS
!=====================================================================*
!   Set Block RHS to the RHS of the the Neutron Kinetics Equations    *
!=====================================================================*         
      implicit none    
      include 'sketch.fh'
! Input: RHS_K(NG,N_TOT)
! Output: MAT_Block_RHS(NG,N_TOT) 
! Local Variables:
      integer k, n 

      do k = 1, N_TOT
         do n = 1, NG
            MAT_Block_RHS(n,k) = MAT_RHS_K(n,k)
         end do
      end do

      return
      end 


      subroutine MAT_Compute_Block_RHS(ADJOINT)
!=====================================================================*
!        Computing Block RHS for the Steady-State Calculations        *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Input: Source(N_TOT), xp(n), XS_SF(NG,N_TOT)
      logical Adjoint
! Output: Block_RHS(NG,N_TOT) 

! Local Variables:
      integer k, n 

      if(ADJOINT) then
! Adjoint Calculations
         do k = 1, N_TOT
            do n = 1, NG
               MAT_Block_RHS(n,k) = XS_SF(n,k)*Source(k)
            end do
         end do
      else 
! Steady-State Diffusion Calculations
        do k = 1, N_TOT
           do n = 1, NG
              MAT_Block_RHS(n,k) = xp(n)*Source(k)
           end do
        end do
      end if


      return
      end


      subroutine MAT_Compute_Block_Diag(ADJOINT, dt_kin)
!=====================================================================*
!     Computing Block Diagonal for the Steady-State Calculations      *
!       & Kinetics Calculations                                       *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Input: 
!  XS_SIK(NG,NG,N_TOT)
!Module MATrix:
! real Diag_Tot(NG,N_TOT)
      real dt_kin
      logical ADJOINT
! Output: MAT_Block_Diag(NG,NG,N_TOT) 
! Local Variables:
      integer k, n , m
      real Time_Der(NG,N_TOT)

      if(ADJOINT) then
! adjoint calculations
       do k = 1, N_TOT
          do n = 1, NG
             do m = 1, NG
                MAT_Block_Diag(n,m,k) = - XS_SIK(m,n,k) -&
                            eigenv_shift*xp(m)*XS_SF(n,k)
             end do
                MAT_Block_Diag(n,n,k) = MAT_Block_Diag(n,n,k) + &
                                                   MAT_Diag_TOT(n,k)
          end do
        end do

      else
! Steady-State and Kinetics Diffusion Calculations

      if(Problem_Type.EQ."Kinetics") then
       DO k = 1, N_TOT
        do n = 1, NG
          Time_Der(n,k) = xs_al(n,k)/dt_kin
!             Time_Der(n) = al(n)*(1./dt_kin + Point_Deriv)
        end do
      END DO  
      else
            Time_Der(:,:) = 0.
      end if

        do k = 1, N_TOT
           do n = 1, NG
              do m = 1, NG
              MAT_Block_Diag(n,m,k) = - XS_SIK(n,m,k) - &
                            eigenv_shift*xpn(n)*XS_SF(m,k)
              end do
              MAT_Block_Diag(n,n,k) = MAT_Block_Diag(n,n,k) + &
                   MAT_Diag_TOT(n,k) + Time_Der(n,k)*volume(k)
           end do
        end do

      end if


      return
      end

            
      subroutine MAT_Inverse_2x2_Block_Diag
!=====================================================================*
!        Computing Inverse of the Block Diagonal Matrix               *
!             2 Neutron  Energy Groups !!!!!!!!!!!!!!                 *                               
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Input: MAT_Block_Diag(NG,NG,N_TOT) 
! Output: MAT_Block_Diag_Inv(NG,NG,N_TOT) 
! Local Variables:
      real tmp_diag, determ_inv
      integer k
!     integer  n , m


      do k = 1, N_TOT

         determ_inv = 1./&
          (MAT_Block_Diag(1,1,k)*MAT_Block_Diag(2,2,k) - &
           MAT_Block_Diag(1,2,k)*MAT_Block_Diag(2,1,k))
         tmp_diag = MAT_Block_Diag(1,1,k)
         MAT_Block_Diag_Inv(1,1,k) = determ_inv*MAT_Block_Diag(2,2,k)
         MAT_Block_Diag_Inv(1,2,k) = - determ_inv*MAT_Block_Diag(1,2,k)
         MAT_Block_Diag_Inv(2,1,k) = - determ_inv*MAT_Block_Diag(2,1,k)
         MAT_Block_Diag_Inv(2,2,k) = determ_inv*tmp_diag

      end do ! N_TOT

      return
      end

      Subroutine MAT_Inverse_Block_Diag
!=====================================================================*
!        Computing Inverse of the Block Diagonal Matrix               *
!                    Using LU Decomosition                            *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'

! Input: MAT_Block_Diag(NG,NG,N_TOT) 
! Output: MAT_Block_Diag_Inv(NG,NG,N_TOT) 
! Local Variables:

      integer k, n , m
!      real Matr_Inv(NG,NG)

      do k = 1, N_TOT
         call MSC_LU_Decomp(MAT_Block_Diag(1,1,k),NG)
         do n = 1, NG
            do m = 1, NG
               MAT_Block_Diag_Inv(n,m,k) = 0.
            end do
            MAT_Block_Diag_Inv(n,n,k) = 1.
         end do
         do m = 1, NG
            call MSC_LUBKSB&
                (MAT_Block_Diag(1,1,k), NG, MAT_Block_Diag_Inv(1,m,k))
         end do
          
!         do n = 1, NG
!            do m = 1, NG
!               MAT_Block_Diag(n,m,k) = Matr_Inv(n,m)
!            end do
!         end do

      end do

      return
      end


      
      subroutine MAT_KIN_Output_COO_Format(dt_kin)
!=====================================================================*
!  Output of the Matrices of the Neutron Kinetics  Problem            *
!    in the Coordinate Format                                         *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Preparing the Matrices
      integer k, n, m, next, i, n_nonzero_AA
      real dt_kin
      real Time_Der(NG,N_TOT)

      real res_norm2, rhs_norm2, RHS_SSOR(NG), lhs_tmp


      DO k = 1, N_TOT
      do n = 1, NG
         Time_Der(n,k) = xs_al(n,k)/dt_kin
      end do
      END DO


      open(2, file ='Output_Debug/KIN_COO.dat', status ='unknown') 
      open(3, file ='Output_Debug/KIN_JR_COO.dat', status ='unknown') 
      open(4, file ='Output_Debug/KIN_JC_COO.dat', status ='unknown') 

      n_nonzero_AA = 0
      
      do k = 1, N_TOT
         do n = 1, NG
            do m = 1, NG
              MAT_Block_Diag(n,m,k) = - XS_SIK(n,m,k) -&
                            eigenv_shift*xpn(n)*XS_SF(m,k)

            end do
              MAT_Block_Diag(n,n,k) = MAT_Block_Diag(n,n,k) + &
                   MAT_Diag_TOT(n,k)  + Time_Der(n,k)*volume(k)
         end do


         do i = 1, NE_T
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
! None-zero elements
              write(2,*) -MAT_Tot(1, i, k), -MAT_Tot(2, i, k)
! Row Index          
                 write(3,*) 1 + (k-1)*NG, 2 + (k-1)*NG
! Column Index
              write(4,*) 1 + (next-1)*NG, 2 + (next-1)*NG
              n_nonzero_AA = n_nonzero_AA + 2
           end if

        end do


! None-zero elements
!         write(*,'(2E12.5)') ((MAT_Block_Diag(n,m,k), m=1, NG), n=1, NG)
!         pause

         write(2,*) MAT_Block_Diag(1,1,k), MAT_Block_Diag(1,2,k), &
                   MAT_Block_Diag(2,1,k), MAT_Block_Diag(2,2,k)
! Row Index          
         write(3,*) 1+(k-1)*NG,1+(k-1)*NG,2+(k-1)*NG,2+(k-1)*NG 

! Column Index
         write(4,*) 1+(k-1)*NG,2+(k-1)*NG,1+(k-1)*NG,2+(k-1)*NG

             n_nonzero_AA=n_nonzero_AA + 4

         do i = 1, NE_T
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
! None-zero elements
              write(2,*) -MAT_Tot(1, i, k), -MAT_Tot(2, i, k)
! Row Index          
                 write(3,*) 1 + (k-1)*NG, 2 + (k-1)*NG
! Column Index
              write(4,*) 1 + (next-1)*NG, 2 + (next-1)*NG
              n_nonzero_AA = n_nonzero_AA + 2
           end if

        end do

      end do

      write(2,*) 'Total number of nonezero elements =', n_nonzero_AA

      close(2) 
      close(3) 
      close(4) 

      open(2, file ='Output_Debug/KIN_RHS.dat', status ='unknown') 
      
      do k = 1, N_TOT
         write(2,*) MAT_RHS_K(1,k), MAT_RHS_K(2,k)
      end do

! Row Index          
      close(2) 

      open(2, file ='Output_Debug/KIN_Sol.dat', status ='unknown') 
      
      do k = 1, N_TOT
         write(2,*) Flux(1,k), Flux(2,k)
      end do

! Row Index          
      close(2) 

! Checking convergence of the solution

      rhs_norm2 = 0.
      res_norm2 = 0.

       do k = 1, N_TOT

         do n = 1, NG
            RHS_SSOR(n) = 0.
         end do

         do i = 1, NE_T
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
           do n = 1, NG
             RHS_SSOR(n) = RHS_SSOR(n) + &
                    Flux(n,next) * MAT_Tot(n, i, k)
           end do ! NG
           end if
         end do

         do n = 1, NG
            RHS_SSOR(n) = RHS_SSOR(n) + MAT_RHS_K(n,k)
            rhs_norm2 = rhs_norm2 + MAT_RHS_K(n,k)*MAT_RHS_K(n,k)
         end do

         do n = 1, NG   
              lhs_tmp = 0.
            do m = 1, NG
              lhs_tmp = lhs_tmp + MAT_Block_Diag(n,m,k)*flux(m,k)
            end do ! m
              res_norm2 = res_norm2 + (lhs_tmp - RHS_SSOR(n))**2
         end do ! NG

      end do
         
      res_norm2 = sqrt(res_norm2)
      rhs_norm2 = sqrt(rhs_norm2)

      write(*,*) 'rhs_norm2 = ', rhs_norm2
      write(*,*) 'res_norm2 = ', res_norm2
      write(*,*) 'Ratio res_norm2/rhs_norm2 =', res_norm2/rhs_norm2


      return 
      end

      subroutine MAT_KIN_Output_CSR_Format(dt_kin)
!=====================================================================*
!  Output of the Matrices of the Neutron Kinetics  Problem            *
!    in the Compressed Sparse Row Format                              *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'
! Preparing the Matrices
      integer k, n, m, next, i, n_nonzero_AA, j
      real dt_kin
      real Time_Der(NG,N_TOT), Mat_Block_Diag_Output(NG, NG)

      real res_norm2, rhs_norm2, lhs_tmp
! Sorting Matrix in Column order
      integer i_order_left(NDIR)
      integer i_order_right(NDIR)
!      data i_order_left  /3, 2, 1/
!      data i_order_right /4, 5, 6/

      character*10 frmt

        IF ( NDIR == 3 )  THEN
          i_order_left(1) =  3
          i_order_left(2) =  2
          i_order_left(3) =  1
          i_order_right(1) = 4
          i_order_right(2) = 5
          i_order_right(3) = 6
        ELSE    
          write(*,*) 'in subroutine MAT_KIN_Output_CSR_Format(dt_kin)'
          write(*,*) ' i_order_left and i_order_right should be '
          write(*,*) ' defined'
          STOP 
        END IF

      DO k = 1, N_TOT  
      do n = 1, NG
         Time_Der(n,k) = xs_al(n,k)/dt_kin
      end do
      END DO

      frmt = "(E15.7)"


      open(2, file ='Output_Debug/KIN_CSR.dat', status ='unknown') 
      open(3, file ='Output_Debug/KIN_IA_CSR.dat', status ='unknown') 
      open(4, file ='Output_Debug/KIN_JA_CSR.dat', status ='unknown') 

      n_nonzero_AA = 0

! Checking convergence of the solution

      rhs_norm2 = 0.
      res_norm2 = 0.
      
      do k = 1, N_TOT

         do n = 1, NG
            write(3,*) n_nonzero_AA + 1
            do m = 1, NG
              MAT_Block_Diag_Output(n,m) = - XS_SIK(n,m,k) -&
                            eigenv_shift*xpn(n)*XS_SF(m,k)
            end do
            MAT_Block_Diag_Output(n,n) = MAT_Block_Diag_Output(n,n) + &
              MAT_Diag_TOT(n,k)  + Time_Der(n,k)*volume(k)

!        write(*,*) 'Neibouring nodes =', (Neib(i,k), i=1, NE_T)
!        pause

! LEFT NEIGHBOURS
         do j = 1, NDIR
           i = i_order_left(j)
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
! None-zero elements
              write(2, frmt) -MAT_Tot(n, i, k)
! Column Index
              write(4,*) n + (next-1)*NG
              n_nonzero_AA = n_nonzero_AA + 1
           end if
        end do

! BLOCK DIAGONAL NONZERO ELEMENTS

         do m = 1, NG
!            if(abs(MAT_Block_Diag_Output(n,m)).GE. eps_round_off) then
               write(2, frmt) MAT_Block_Diag_Output(n,m)
! Column Index
               write(4,*) m+(k-1)*NG
                   n_nonzero_AA=n_nonzero_AA + 1
!             end if
          end do

! RIGHT NEIGHBOURS
         do j = 1, NDIR
           i = i_order_right(j) ! NDIR+1, NE_T
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
! None-zero elements
              write(2, frmt) -MAT_Tot(n, i, k)
! Column Index
              write(4,*) n + (next-1)*NG
              n_nonzero_AA = n_nonzero_AA + 1
           end if
        end do

!        write(*,*) 'Neibouring nodes =', (Neib(i,k), i=1, NE_T)
!        pause
      end do ! NG

! CHECKING CONVERGENCE

      do n = 1, NG
      ! NONDIAGONAL ELEMENTS

         lhs_tmp = 0.
         do i = 1, NE_T
           next  = Neib(i,k)
           if(next.ne.I_BOUND_NODE) then
             lhs_tmp  = lhs_tmp - &
                    Flux(n,next)*MAT_Tot(n, i, k)
           end if
         end do

! BLOCK DIAGONAL ELEMENTS
         do m = 1, NG
            lhs_tmp = lhs_tmp + MAT_Block_Diag_Output(n,m)*flux(m,k)
         end do ! m

         rhs_norm2 = rhs_norm2 + MAT_RHS_K(n,k)*MAT_RHS_K(n,k)
         
         res_norm2 = res_norm2 + (lhs_tmp - MAT_RHS_K(n,k))**2

      if((k.eq.1).OR.K.eq.N_TOT) then
          write(*,'(2I10,A,E12.5)') k, n, "left hand side = ", lhs_tmp
      end if

      end do ! NG

      end do !! k

      write(2,*) 'Total number of nonezero elements =', n_nonzero_AA

      write(3,*) n_nonzero_AA + 1

      close(2) 
      close(3) 
      close(4) 

      open(2, file ='Output_Debug/KIN_RHS.dat', status ='unknown') 
      
      do k = 1, N_TOT
         write(2, frmt) (MAT_RHS_K(n,k), n = 1, NG)
      end do

! Row Index          
      close(2) 

      open(2, file ='Output_Debug/KIN_Sol.dat', status ='unknown') 
      
      do k = 1, N_TOT
         write(2, frmt) (Flux(n,k), n = 1, NG)
      end do

! Row Index          
      close(2) 

      res_norm2 = sqrt(res_norm2)
      rhs_norm2 = sqrt(rhs_norm2)

      write(*,*) 'rhs_norm2 = ', rhs_norm2
      write(*,*) 'res_norm2 = ', res_norm2
      write(*,*) 'Ratio res_norm2/rhs_norm2 =', res_norm2/rhs_norm2


      return 
      end

      subroutine MAT_KIN_Output_Dump
!=====================================================================*
!  Output of the SKETCH Matrices  into the dump unformatted file      *
!    "Output_Debug/Mat_Kin.dmp"                                       *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none    
      include 'sketch.fh'

      open(io_unit, file = "Output_Debug/Mat_Kin.dmp", &
                              form ='unformatted', status ='unknown')

      write(io_unit) NG, NX, NY, NZ, NH, N_TOT, NE_T
      write(io_unit) Neib
      write(io_unit) Mat_Tot
      write(io_unit) MAT_Block_Diag
      write(io_unit) Mat_rhs_k
      write(io_unit) flux
      close(io_unit)

      return
      end

         
