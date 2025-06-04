      MODULE GeoMeTry_Triangle
      
      IMPLICIT NONE

      INTEGER, SAVE :: NP_SIDE
      INTEGER, ALLOCATABLE, SAVE :: n_fine_wv(:, :, :) 
!                               2*NP_SIDE, 2*NP_SIDE, 2
!      INTEGER, ALLOCATABLE, SAVE :: la_hex(:, :) ! N_POLY, NE+1 
      INTEGER, ALLOCATABLE, SAVE :: la_hex(:,:)  

      CONTAINS

       SUBROUTINE Allocate_Triangle_LA_HEX(NE, N_POLY)
       
       INTEGER, INTENT(IN) :: NE, N_POLY

        ALLOCATE ( la_hex(N_POLY, NE+1) )

       RETURN 
       END SUBROUTINE Allocate_Triangle_LA_hex          

       SUBROUTINE Set_Triangle_Local_Mesh(NCHM) 
!======================================================================!
! Setting up a triangular spatial mesh (v,w)                           !
!                                  inside the hexagonal region         !
!======================================================================!

         INTEGER, INTENT(IN):: NCHM

         INTEGER :: nlv, nlw, n_f

         NP_SIDE = NINT( SQRT(REAL(NCHM/6))) 
                                
         ALLOCATE ( n_fine_wv(2*NP_SIDE, 2*NP_SIDE, 2 ) )

         n_fine_wv(:,:,:) = 0

         n_f = 0
         DO nlv = 1, NP_SIDE
            DO nlw = 1, nlv + NP_SIDE - 1
              n_f = n_f + 1
              n_fine_wv(nlw, nlv,1) = n_f
              n_f = n_f + 1 
              n_fine_wv(nlw, nlv,2) = n_f
            END DO
            nlw = nlv + NP_SIDE
            n_f = n_f + 1
            n_fine_wv(nlw, nlv,1) = n_f
         END DO

         DO nlv = NP_SIDE+1, 2*NP_SIDE
            nlw  = nlv - NP_SIDE 
            n_f = n_f + 1
            n_fine_wv(nlw, nlv,2) = n_f
              DO nlw = nlv - NP_SIDE + 1, 2*NP_SIDE 
              n_f = n_f + 1
              n_fine_wv(nlw, nlv,1) = n_f
              n_f = n_f + 1 
              n_fine_wv(nlw, nlv,2) = n_f
            END DO
         END DO
      RETURN
      END SUBROUTINE Set_Triangle_Local_Mesh      

      SUBROUTINE Output_Triangle_Local_Mesh 
      INTEGER :: nlv, nlw

      DO nlv = 1, 2*NP_SIDE
            WRITE(*, '(10(3x,2I3))' ) ( n_fine_wv(nlw, nlv,1),&
               n_fine_wv(nlw, nlv,2), nlw=1, 2*NP_SIDE)
      END DO

      RETURN
      END SUBROUTINE Output_Triangle_Local_Mesh 

      SUBROUTINE Set_LA_Triangle(NCHM, NH, NE, NXR, NYR, npoly, la, &
        I_BOUND_NODE)
!======================================================================!
! Setting up LA(NH, NE+1) for trinagular mesh inside hexagons          !
!  new variant LA(k, nd) = I_BOUND_NODE if there is no neighbours                !
!======================================================================!
      INTEGER, INTENT(IN)  :: NE, NXR, NYR, NH, NCHM, I_BOUND_NODE
      INTEGER, INTENT(IN)  :: npoly(NYR, NXR)
      INTEGER, INTENT(OUT) :: la(NH, NE+1)
! locals
      INTEGER :: nlx, nly, nlv, nlw, nf, nf_right, np, np_shift, nd
!
      nf=0
      nd = 1
! x direction        
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            IF(np /= 0) THEN
               DO nlv = 1, NP_SIDE
                  DO nlw = 1, nlv + NP_SIDE - 1
                     nf = nf + 1
                     la(nf, nd + NE/2) = nf+1
                     la(nf+1, nd ) = nf
                     la(nf, NE+1) = nf
                     nf = nf +1 
                     la(nf, NE+1) = nf
                  END DO ! nlw
                  nf = nf +1 
                  la(nf, NE+1) = nf
                  IF(nlx /= NXR .AND. npoly(nly, nlx+1) /= 0) THEN
                     nf_right = np*NCHM + n_fine_wv(nlv, nlv+NP_SIDE, 2)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                  ELSE
                     la(nf, nd + NE/2) = I_BOUND_NODE
                  END IF
               END DO ! nlv
!               write(*,*) 'np, nf =', np, nf
!               pause
               DO nlv = 1, NP_SIDE
                  nf = nf + 1
                  la(nf, NE+1) = nf
                  IF(nlx == 1 .OR. npoly(nly, nlx-1) == 0) THEN
                     la(nf, nd) = I_BOUND_NODE
                  END IF                    
                  DO nlw = nlv + 1, 2*NP_SIDE
                      nf= nf + 1
                     la(nf, nd + NE/2) = nf+1
                     la(nf+1, nd ) = nf
                     la(nf, NE+1) = nf
                     nf = nf +1 
                     la(nf, NE+1) = nf
                  END DO ! nlw
               END DO ! nlv
            END IF ! np /= 0
!               write(*,*) 'np, nf =', np, nf
!               pause
         END DO ! nlx
      END DO ! nly

!     WRITE(*,*) 'nf = ', nf
! y - direction
      nd = 2
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            np_shift = (np-1)*NCHM
            IF(np /= 0) THEN
               DO nlw = 1, NP_SIDE
                  IF( nly==1 .OR. npoly(nly-1, nlx) ==0) THEN
                     nlv = 1
                     nf = np_shift + n_fine_wv(nlw, nlv, 2)
                     la(nf, nd) = I_BOUND_NODE
                  END IF
                  DO nlv = 1, nlw + NP_SIDE - 1
                     nf = np_shift + n_fine_wv(nlw, nlv, 1)
                     nf_right = np_shift + n_fine_wv(nlw, nlv+1, 2)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                  END DO
               END DO
               DO nlw = 1, NP_SIDE
                  DO nlv = nlw, 2*NP_SIDE - 1
                     nf = np_shift + n_fine_wv(nlw+NP_SIDE, nlv, 1)
                     nf_right = np_shift + &
                          n_fine_wv(nlw+NP_SIDE, nlv+1, 2)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                  END DO
                     nf=np_shift+n_fine_wv(nlw+NP_SIDE, 2*NP_SIDE, 1)
                     IF( nly/=NYR .AND. &
                          npoly(nly+1, nlx) /=0) THEN
                     nf_right = (npoly(nly+1, nlx)-1)*NCHM + &
                          n_fine_wv(nlw, 1, 2)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                     ELSE
                     la(nf, nd + NE/2) = I_BOUND_NODE
                     END IF
              END DO
            END IF ! np /= 0
         END DO ! nlx
      END DO ! nly


! w direction        
      nd = 3
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            np_shift = (np-1)*NCHM
            IF(np /= 0) THEN
               DO nlv = 1, NP_SIDE
                  IF (nlx == 1 .OR. nly==1 .OR. &
                          npoly(nly-1, nlx-1) ==0) THEN
                     nlw = 1
                     nf = np_shift + n_fine_wv(nlw, nlv, 1)
                     la(nf, nd) = I_BOUND_NODE
                  END IF                    
                  DO nlw = 1, nlv + NP_SIDE - 1
                     nf = np_shift + n_fine_wv(nlw, nlv, 2)
                     nf_right = np_shift + n_fine_wv(nlw+1, nlv, 1)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                  END DO
               END DO
               DO nlv = 1, NP_SIDE
                  DO nlw = nlv,  2*NP_SIDE - 1
                     nf = np_shift + n_fine_wv(nlw, nlv+NP_SIDE, 2)
                     nf_right = np_shift + &
                          n_fine_wv(nlw+1, nlv+NP_SIDE, 1)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                  END DO
                     nf=np_shift+n_fine_wv(2*NP_SIDE, nlv+NP_SIDE, 2)
                     IF( nlx /= NXR .AND. nly/=NYR .AND. &
                          npoly(nly+1, nlx+1) /=0) THEN
                     nf_right = (npoly(nly+1, nlx+1)-1)*NCHM + &
                          n_fine_wv(1, nlv, 1)
                     la(nf_right, nd) = nf
                     la(nf, nd + NE/2) = nf_right
                     ELSE
                     la(nf, nd + NE/2) = I_BOUND_NODE
                     END IF
              END DO
            END IF ! np /= 0
         END DO ! nlx
      END DO ! nly
      

      RETURN
      END SUBROUTINE Set_LA_Triangle


      SUBROUTINE Set_Numb_Triangle(N_TOT, NCHM, NDIR, NH, NXR, NYR, NZ,&
       npoly, Numb, N_TOT_DIR)
!======================================================================!
! Setting up Numb(N_TOT, NDIR) for trinagular mesh inside hexagons 
!======================================================================!
      include 'units.fh'
      INTEGER, INTENT(IN)  :: N_TOT, NDIR,  NXR, NYR, NH, NCHM, NZ

      INTEGER, INTENT(IN)  :: npoly(NYR, NXR)
      INTEGER, INTENT(OUT) :: Numb(N_TOT, NDIR), N_TOT_DIR(NDIR)
! locals
      INTEGER :: nlx, nly, nlv, nlw, nf,  np, np_shift, nd, nn, n1, k, &
       nd_tot 

      nd = 1
      nd_tot = 0
      DO n1 = 1, NZ
         nn = (n1 - 1)*NH
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            np_shift = (np-1)*NCHM
            IF(np /= 0) THEN
               DO nlv = 1, NP_SIDE
                  DO nlw = 1, nlv + NP_SIDE - 1
                     nd_tot = nd_tot + 1
                     nf = nn +  np_shift + n_fine_wv(nlw, nlv, 2)
                     Numb(nd_tot, nd) = nf
                  END DO ! nlw
! boundary right face
                  IF( nlx == NXR .OR. npoly(nly, nlx+1) == 0) THEN
                      nlw = nlv + NP_SIDE
                      nd_tot = nd_tot + 1
                      nf = nn +  np_shift + n_fine_wv(nlw, nlv, 1)
                      Numb(nd_tot, nd) = nf
                  END IF 
               END DO ! nlv
!               write(*,*) 'np, nf =', np, nf
!               pause
               DO nlv = 1, NP_SIDE
! left boundary face (always numbering)
!                  IF(nlx == 1 .OR. npoly(nly, nlx-1) == 0) THEN
!                   nf = nn+np_shift+n_fine_wv(nlv,nlv+NP_SIDE,2)
!                   nd_tot = nd_tot + 1
!                   Numb(nd_tot, nd) = nf
!                  END IF                    
                  DO nlw = nlv, 2*NP_SIDE
                     nd_tot = nd_tot + 1
                     nf = nn +  np_shift + &
                          n_fine_wv(nlw, nlv+NP_SIDE, 2)
                     Numb(nd_tot, nd) = nf
                  END DO ! nlw
               END DO ! nlv
            END IF ! np /= 0
!               write(*,*) 'np, nf =', np, nf
!               pause
         END DO ! nlx
      END DO ! nly
      END DO ! NZ

      N_TOT_DIR(nd) = nd_tot

      nd = 2
      nd_tot = 0
      DO n1 = 1, NZ
         nn = (n1 - 1)*NH
         DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            np_shift = (np-1)*NCHM
            IF(np /= 0) THEN
               DO nlw = 1, NP_SIDE
! left boundary node (always numbering)
!                     IF(nly==1 .OR. npoly(nly-1,nlx)==0) THEN 
!                        nd_tot = nd_tot + 1
!                        nf = np_shift + n_fine_wv(nlw, 1, 2) + nn 
!                        Numb(nd_tot, nd) = nf
!                     END IF
                  DO nlv = 1, nlw + NP_SIDE 
                     nf = nn +  np_shift + n_fine_wv(nlw, nlv,2)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf 
                  END DO
               END DO
               DO nlw = 1, NP_SIDE
                  DO nlv = nlw, 2*NP_SIDE - 1
                     nf = nn  + np_shift + &
                           n_fine_wv(nlw+NP_SIDE, nlv+1, 2)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf
                  END DO
! boundary right face
                     IF(nly == NYR.OR.npoly(nly+1,nlx)==0) THEN
                     nf=nn + np_shift+&
                          n_fine_wv(nlw+NP_SIDE, 2*NP_SIDE, 1)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf
                     END IF
              END DO
            END IF ! np /= 0
         END DO ! nlx
      END DO ! nly
      END DO ! NZ

      N_TOT_DIR(nd) = nd_tot

      nd = 3
      nd_tot = 0
      DO n1 = 1, NZ
         nn = (n1 - 1)*NH
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            np_shift = (np-1)*NCHM
            IF(np /= 0) THEN
               DO nlv = 1, NP_SIDE
! left boundary node (always numbering)
!                  nk = np_shift + n_fine_wv(1, nlv, 1)
!                  IF(la(nk, nd).eq.0) THEN
!                   nf = nn  + nk
!                   nd_tot = nd_tot + 1
!                   Numb(nd_tot, nd) = nf
!                  END IF                    
                 DO nlw = 1, nlv + NP_SIDE
                     nf = nn + np_shift + n_fine_wv(nlw, nlv, 1)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf
                  END DO
               END DO
               DO nlv = 1, NP_SIDE
                  DO nlw = nlv,  2*NP_SIDE - 1
                     nf = nn + np_shift+n_fine_wv(nlw+1, nlv+NP_SIDE,1)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf
                  END DO
! boundary right face
                     IF(nly == NYR.OR.nlx==NXR.OR. &
                            npoly(nly+1,nlx+1)==0) THEN
                     nf=np_shift+n_fine_wv(2*NP_SIDE, nlv+NP_SIDE, 2)
                     nd_tot = nd_tot + 1
                     Numb(nd_tot, nd) = nf
                     END IF
              END DO
            END IF ! np /= 0
         END DO ! nlx
      END DO ! nly
      END DO ! 

      N_TOT_DIR(nd) = nd_tot
     

      IF(DEBUG) THEN
        OPEN(io_unit, file='Output_Debug\Numb.dat', status='unknown')
           write(io_unit, *) ' ND =1', N_TOT_DIR(1) 
           write(io_unit, '(10I5)') (Numb(k,1), k=1, N_TOT)
           write(io_unit, *) ' ND =2', N_TOT_DIR(2) 
           write(io_unit, '(10I5)') (Numb(k,2), k=1, N_TOT)
           write(io_unit, *) ' ND =3', N_TOT_DIR(3) 
           write(io_unit, '(10I5)') (Numb(k,3), k=1, N_TOT)
        CLOSE(io_unit)
      END IF           

      RETURN
      END SUBROUTINE Set_Numb_Triangle



      END MODULE GeoMeTry_Triangle
      
!      PROGRAM MAIN

!      USE GeoMeTry_Triangle

!      INTEGER :: NCHM

!      WRITE(*,*) 'Input NCHM'
!      READ(*,*) NCHM 

!      CALL Set_Triangle_Local_Mesh(NCHM) 
!      CALL Output_Triangle_Local_Mesh 

!      STOP
!      END PROGRAM MAIN
