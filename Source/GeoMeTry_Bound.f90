      MODULE GeoMeTry_Boundary

      INTEGER :: N_BOUND_M, N_BOUND

      INTEGER, DIMENSION(:), ALLOCATABLE ::&
       k_bound(:), nl_bound(:), nd_bound(:) 

      CONTAINS

      SUBROUTINE Allocate_Boundary
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

            N_BOUND_M = 2*(NH + NZ*(NX + NY + NV))

!      WRITE(*,*) 'N_BOUND_M =', N_BOUND_M


            ALLOCATE( k_bound(N_BOUND_M), nl_bound(N_BOUND_M), &
           nd_bound(N_BOUND_M) ) 

      RETURN
      END SUBROUTINE Allocate_Boundary

      Subroutine GMT_Set_Boundary
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
       
      INTEGER nd, kb, k, kt

      CALL Allocate_Boundary

      kb = 0
      do nd = 1, NDIR
!      write(*,*) 'nd, N_TOT_DIR(nd)  = ', nd, N_TOT_DIR(nd) 
!      write(*, '(10i3)') ( Numb(kt, nd), kt =1, N_TOT_DIR(nd) )
!      pause
       DO k = 1, N_TOT_DIR(nd) 
          kt = Numb(k, nd) 
!          write(*,*) 'kt =', kt
! Left Boundary Nodes
          IF(Neib(nd, kt).eq.I_BOUND_NODE) THEN
                  kb = kb + 1
!                  write(*,*) 'kb =', kb
!                  write(*,*) 'k =', kt
                  k_bound(kb) = kt
                  nl_bound(kb) = 1
                  nd_bound(kb) = nd            
          END IF
! Right Boundary Nodes
          IF(Neib(nd+NDIR, kt).eq.I_BOUND_NODE) THEN
                  kb = kb + 1
                  k_bound(kb) = kt
                  nl_bound(kb) = 2
                  nd_bound(kb) = nd            
          END IF
        END DO
      END DO

      N_BOUND = kb

      RETURN
      END SUBROUTINE GMT_Set_Boundary

      END MODULE GeoMeTry_Boundary


      Subroutine GMT_Set_Boundary_Conditions
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

      INTEGER nd, n

      DO nd = 1, NDD

         if(i_dr(1,nd).eq.0) then
               cg(1, nd) = 0.
         else
               cg(1, nd) = 1.
         end if

         if(i_dr(2,nd).eq.0) then
               cg(2, nd) = 0.
         else
               cg(2, nd) = 1.
         end if

         do n = 1, NG

         if(i_dr(1,nd).eq.0) then
            fg(n, 1, nd) = 1.
         else
            fg(n, 1, nd) =  dr(n,1,nd)
         end if

          if(i_dr(2,nd).eq.0) then
           fg(n, 2, nd) = 1.
          else
           fg(n, 2, nd) =  dr(n,2,nd)
         end if

         END DO ! NG

      END DO ! NDD


      RETURN
      END