      subroutine MAT_Update_FD_Matrix
!=====================================================================*
!      Normalization of the Initial Neutron Flux Distribution         *
!      Slava (c) 23.II.1998                                           *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      call MAT_Compute_FD_Matrix

      return
      end

      subroutine MAT_Transpose_Matrix
!=====================================================================*
! Transpose Matrix of the Couling Coefficients                        *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input:
!   Index_Oppos(NE_T) - couling coefficient index in the neibouring node
!   Neib(NE_T, N_TOT) - neibouring Node number

! Module MATrix
!     real  MAT_Tot(NG, NE_T, N_TOT) - matrix of the coupling coefficients
! Output:

! real  MAT_TOT(NG, NE_T, N_TOT) - transposed matrix of the coupling 
!                                    coefficients
! Local variables:
      integer nd, i_right, next, n, k
      real tmp

      do k = 1, N_TOT
         do nd = 1, NDD 
            i_right = Index_Neib(nd, 2) ! right neibour
            next = Neib(i_right, k)
            if(next.ne.I_BOUND_NODE.and.next.GT.0) THEN !
            do n = 1, NG
               tmp = MAT_Tot(n, i_right, k)
               MAT_Tot(n, i_right, k) = MAT_Tot(n, nd, next)
               MAT_Tot(n, nd, next) = tmp
            end do
            end if
         end do
      end do

      return
      end

      subroutine MAT_Compute_FD_Matrix
!=====================================================================*
!             Calculation of Coupling Matrix                          *
!  of the  Finite-Difference Mesh-Centered Method                     *
!             Cartesian X-Y-Z Geometry                                *
!=====================================================================*
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'

! Input: N_TOT NDIR, N_BOUND, k_bound(N_BOUND), nl_bound(N_BOUND),
!         nd_bound(N_BOUND), Index_Neib(NDIR,2), Neib(NE_T, N_TOT),
!         s_xyz(NDIR,N_TOT), h_xyz(NDIR,N_TOT), XS_D(NG,N_TOT),
!         i_dr(2,3), dr(n,2,3), MAT_Tot(NG, NE_T, N_TOT)

! Output: MAT_FD(NG, NTSJ) - Coupling Matrix of the Mesh-Centered
!                       Finite-Difference Method
!         MAT_Diag_FD(NG, N_TOT) - Diagonal of the Coupling Matrix

! Internal variables
      integer n, k, nl, nd,  k_left, k_right, &
             i_right, kt
      real hk, hk_neib, dif, dif_neib, denom
      LOGICAL flag_rep_str_left,  flag_rep_str_right

!      do k=1,N_TOT
!         do n = 1, NG
!            MAT_Diag_FD(n,k)=0.
!         end do
!      end do
        
! Coupling Coefficients between the Nodes (Internal Interface)
!      write(*,*) 'NTSJ =', NTSJ
!      pause

      nl = 0

      do nd = 1, NDD

        i_right = Index_Neib(nd,2)

        flag_rep_str_left=.False.
        flag_rep_str_right=.False.

        IF( i_dr(1, nd) == 2 ) flag_rep_str_left=.True.
        IF( i_dr(2, nd) == 2 ) flag_rep_str_right=.True.

        do kt = 1, N_TOT_DIR(nd)

         k = Numb(kt, nd)

!         write(*,*) 'nd =', nd, 'k=', k

         hk = h_xyz(nd,k)


         k_left = Neib(nd, k)

         nl = nl + 1 ! Left Interface

!        write(*,*) 'k,nl, nd, k_left =', k, nl, nd, k_left
!        pause


         if(k_left.eq.I_BOUND_NODE) then

!          write(*,*) 'hk, s_xyz=', hk, s_xyz(nd,k)
!          write(*,*) 'nd, k(left boundary) =', nd, k
!          pause

           IF( .NOT. flag_rep_str_left ) THEN
            DO n = 1, NG
               dif = XS_D(n,k)
               denom = fg(n,1,nd)*hk + &
                  2.*cg(1, nd)*(dif + D_Nod(n,k,1,nd) )
               MAT_FD(n, nl) = 2.*fg(n,1,nd)*dif*s_xyz(nd,k)/denom
            END DO
           END IF

! Left Boundary Node
!           if(i_dr(1,nd).eq.0) then
! Extrapolation Distance
!              do n = 1,NG
!                 dif = XS_D(n,k)
!                 coupl_coef = dif*s_xyz(nd,k)/(0.5*hk + dr(n,1,nd))
!                 MAT_FD(n,nl) = coupl_coef 
!                 MAT_Diag_FD(n,k) = MAT_Diag_FD(n,k) + coupl_coef
!              end do
!             else
! Logarithmic Derivative
!              do n = 1,NG
!                 dif = XS_D(n,k)
!
!                 coupl_coef = dr(n,1,nd)*dif*s_xyz(nd,k)/&
!                                  (0.5*hk*dr(n,1,nd) + dif)
!                 MAT_FD(n,nl) = coupl_coef 
!                 MAT_Diag_FD(n,k) = MAT_Diag_FD(n,k) + coupl_coef
!                write(*,*) 'nl =', nl, 'k =', k, 'coupl_coef=', &
!                             coupl_coef
!                 pause
!              end do
!            end if 

          else if(k_left.ne.0) THEN ! Left /= 0
! Normal Interface

            hk_neib = h_xyz(nd,k_left)

!          write(*,*) 'nd, k, k_left =', nd, k, k_left
!          write(*,*) 'hk, hk_neib, s_xyz=', hk,  hk_neib, s_xyz(nd,k)
!          pause

            do n = 1, NG
                dif =  XS_D(n,k)
                dif_neib = XS_D(n,k_left)
                denom = ( dif + D_Nod(n,k,1,nd) )*hk_neib +&
                   (dif_neib + D_Nod(n,k_left,2,nd))*hk
                MAT_FD(n, nl) = 2.*(dif*dif_neib - D_Nod(n,k,1,nd)&
               *D_Nod(n,k_left,2,nd))*s_xyz(nd,k)/denom

!                coupl_coef = 2.*s_xyz(nd,k)*dif*dif_neib/&
!                                 (dif_neib*hk + dif*hk_neib)
!                MAT_FD(n,nl) = coupl_coef

!                MAT_Diag_FD(n,k) = MAT_Diag_FD(n,k) + coupl_coef
!                MAT_Diag_FD(n,k_left) = MAT_Diag_FD(n,k_left) + &
!                                                        coupl_coef
             end do ! NG

         end if ! Left / = 0

         k_right = Neib(i_right, k)
         if(k_right.eq.I_BOUND_NODE) then

            nl = nl + 1 ! Right Boundary Interface

            IF ( .NOT. flag_rep_str_right ) THEN
            DO n = 1, NG
                 dif = XS_D(n,k)
                 denom = fg(n,2,nd)*hk + &
                    2.*cg(2, nd)*(dif + D_Nod(n,k,2,nd))
                 MAT_FD(n, nl) = 2.*fg(n,2,nd)*dif*s_xyz(nd,k)/denom
            END DO 
            END IF

!            if(i_dr(2,nd).eq.0) then
! Extrapolation Distance
!              do n = 1,NG
!                 dif = XS_D(n,k)
!                 coupl_coef = dif*s_xyz(nd,k)/(0.5*hk + dr(n,2,nd))
!                 MAT_FD(n,nl) = coupl_coef 
!                 MAT_Diag_FD(n,k) = MAT_Diag_FD(n,k) + coupl_coef
!              end do
!             else
! Logarithmic Derivative
!              do n = 1,NG
!                 dif = XS_D(n,k)
!                 coupl_coef = dr(n,2,nd)*dif*s_xyz(nd,k)/&
!                                  (0.5*hk*dr(n,2,nd) + dif)
!                 MAT_FD(n,nl) = coupl_coef 
!                 MAT_Diag_FD(n,k) = MAT_Diag_FD(n,k) + coupl_coef
!              end do
!            end if ! Boudary Condition

         end if ! k_right.eq.0)
        end do ! N_TOT
       end do ! NDD

      if(Debug) CALL MAT_DBG_Write_Matr_FD

      return
      end

      subroutine MAT_Set_Total_Matrix
!=====================================================================!
! Set Total Matrix to the Sum of Finite-Difference and Nodal Matrices !
!             Cartesian X-Y-Z Geometry                                !
! (c) Slava 1 Dec. 1999                                               !
!=====================================================================!
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'

! Input:  MAT_FD(NG, NE_T, N_TOT) 
! Output:  MAT_Tot(NG, NE_T, N_TOT)
! Local Variables

      integer k, n, nd, k_right, k_left, nd_right, nl , kt
!      real coupl_coef

      LOGICAL flag_rep_str_left,  flag_rep_str_right


      do k=1,N_TOT
         do n = 1, NG
            MAT_Diag_Tot(n,k)=XS_SA(n,k)
         end do
      end do


      nl = 0
      DO nd = 1, NDD
        nd_right = Index_Neib(nd,2)
        flag_rep_str_left=.False.
        flag_rep_str_right=.False.
        IF( i_dr(1, nd) == 2 ) flag_rep_str_left=.True.
        IF( i_dr(2, nd) == 2 ) flag_rep_str_right=.True.

        DO kt = 1, N_TOT_DIR(nd)

          k = Numb(kt, nd)
          nl = nl + 1 ! Left Interface

          k_left = Neib(nd, k)
          k_right = Neib(nd_right, k)

! LEFT Interface          

!           Coupl_Coef = MAT_FD(n, nl) 
!           MAT_Tot(n, nd, k) = Coupl_Coef 

            IF(k_left.eq.I_BOUND_NODE.AND.flag_rep_str_left) THEN 
! THE CASE OF THE REPETEATED STRUCTURE at left interface of the boundary node
            DO n = 1, NG
             MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) - MAT_Nod(n,nl+1) &
                    + MAT_FD(n,nl+1) 
             MAT_Tot(n,nd, k) = MAT_FD(n, nl+1) + MAT_Nod(n,nl+1)
            END DO
            ELSE 
! any other interface
            DO n = 1, NG
              MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_FD(n, nl) + &
                MAT_Nod(n,nl)
              MAT_Tot(n, nd, k) = MAT_FD(n, nl) - MAT_Nod(n,nl)
            END DO 
            END IF ! k_left.eq.I_BOUND_NODE.AND.flag_rep_str_left

           if(k_left.ne.I_BOUND_NODE.AND.k_left.GT.0) then
              DO n = 1, NG
!              MAT_Tot(n, i_right, k_left) = Coupl_Coef 
               MAT_Diag_Tot(n, k_left) = MAT_Diag_Tot(n, k_left) - &
                                 MAT_Nod(n,nl) + MAT_FD(n,nl)

               MAT_Tot(n, nd_right, k_left) = MAT_FD(n, nl) + &
                    MAT_Nod(n,nl)
              END DO 
           end if

! RIGHT Interface

          if(k_right.eq.I_BOUND_NODE) then
            nl = nl + 1 ! Right Interface
            IF ( .NOT. flag_rep_str_right ) THEN
            do n = 1, NG
!              Coupl_Coef = MAT_FD(n, nl) 
!              MAT_Tot(n, i_right, k) = Coupl_Coef 
             MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) - MAT_Nod(n,nl) &
                    + MAT_FD(n,nl) 
             MAT_Tot(n,nd_right,k) = MAT_FD(n, nl) + MAT_Nod(n,nl)
            end do
            ELSE
! THE CASE of the REPETEAD STRUCTURE
            do n = 1, NG
            MAT_Diag_Tot(n, k) = MAT_Diag_Tot(n, k) + MAT_FD(n, nl-1)+ &
                MAT_Nod(n,nl-1)
            MAT_Tot(n,nd_right, k) = MAT_FD(n, nl-1) - MAT_Nod(n,nl-1)
            end do
            END IF ! .NOT. flag_rep_str_right 

          end if
        end do
      end do

      if(Debug) then
        call MAT_DBG_Check_Matr_Tot
      end if

      if(Debug) call MAT_DBG_Write_Matr_Tot

      return
      end

      subroutine MAT_DBG_Check_Matr_Tot
!=====================================================================*
! Check Diagonal Dominance  of the Iteration Matrix                   *
! And output if there is Non-Diagonal Dominant Elements               *        
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer n, k, nn, i, j, n1
      real Diag_Dom

   
      open(io_unit, file = 'Output_Debug/Matrix.dat', status ='unknown')
      
      do n = 1, NG
        do n1 = 1, NZ
           nn = (n1 - 1)*NH
           do k = 1, NH
              Diag_Dom = 0.
              do i = 1, NE_T
                Diag_Dom = Diag_Dom + abs(MAT_Tot(n, i, k + nn))
              end do
              IF(  MAT_Diag_Tot(n, k + nn) == 0. ) THEN
                WRITE(*,*) 'n, k, nn =', n, k, nn, &
                        MAT_Diag_Tot(n, k + nn) 
                READ(*,*) 
              END IF !
              Diag_Dom = ( MAT_Diag_Tot(n, k + nn) - Diag_Dom ) / &
                         MAT_Diag_Tot(n, k + nn)
              if(Diag_Dom.LE.0) then
                write(*,*) 'Matrix is not Diagonally Dominant'
                write(io_unit,*) 'ngroup, nz=', n, n1, ' Diag_Dom = ',&
                           Diag_Dom
                write(io_unit,*) 'k,   (MAT_Tot(n, i, k + nn),i=1,6),&
                                   MAT_Diag_Tot(n,k)'
                write(io_unit,1) k,   (MAT_Tot(n, i, k + nn),i=1,6),&
                                   MAT_Diag_Tot(n,k)

                write(io_unit,*)  'Assembly Number = ', np_out(k)
                write(io_unit,*) ' Numbers inside Assembly' 
                write(io_unit,9) (poly_out(np_out(k),j),j=1,NCHM)
               end if
          end do
        end do
      end do

      close(io_unit)

    9 format(20x,10I5)  
    1 format(1x, i4, 7E12.5)

      return
      end

      subroutine MAT_DBG_Write_Matr_Tot
      implicit none
      include 'sketch.fh'
!=====================================================================*
!                  DEBUG Output of the matrix                         *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

! Local Variables
      integer k, i, kt, n1

      open(io_unit,file ='Output_Debug/Matrix_TOT.dat', &
           status ='unknown')
        write(io_unit,*) 'Thermal Energy Group NG'
        do n1 = 1, 1
           write(io_unit,*) 'n1=', n1
           write(io_unit,*)
           do kt = 1, NH
              k = kt + (n1-1)*NH
              write(io_unit,1) k, (MAT_Tot(NG,i,k),i=1,NE_T),&
                          MAT_Diag_Tot(NG,k)
           end do
        end do

      close(io_unit)

    1 format(1x, I4, 9E12.5)

      return
      end


      subroutine MAT_DBG_Write_Matr_FD
      USE GeoMeTry_Faces
      implicit none
      include 'sketch.fh'
!=====================================================================*
!                  DEBUG Output of the matrix                         *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

! Local Variables
      integer k,  n1

      open(io_unit,file ='Output_Debug/Matrix_FD.dat', &
           status ='unknown')
        write(io_unit,*) 'Thermal Energy Group NG'
        do n1 = 1, 1
           write(io_unit,*) 'n1=', n1, 'N_FACES=', N_FACES
           write(io_unit,*)
           do k = 1, N_FACES              
              write(io_unit,1) k, MAT_FD(NG,k)
           end do
        end do

      close(io_unit)

    1 format(1x, I10, 9E12.5)

      return
      end
