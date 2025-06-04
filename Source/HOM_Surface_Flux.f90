      SUBROUTINE HOM_Set_Surface_Flux
!=====================================================================!        
! Calculation of the Surface Neutron Flux for the Fine Spatial Mesh   !
!  (c) Slava 5 June 2001 MEPhI                                        !
!=====================================================================!        
      USE GeoMeTry_Faces, ONLY : FluxSurface, CurrentSurface 
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

      INTEGER :: nd, k, k_right, k_left, nl, nd_right, kt, n
      REAL    :: d0(NG), d0_neib(NG), d0_nod(NG), d0_nod_neib(NG), &
                hk, hk_neib  
      INTEGER :: i_right
!      real const_hk ! a/H (x-side of rectangle / pitch of hex)
!      parameter (const_hk = 1.0525557)
!      real const_ju ! b/R (y-side of rectangle / side of hex)
!      real cg_nod(2, NDIR), hk_nod
!      logical Flag_Hex 


!      Flag_Hex  = .False. 
      nl = 0

      do nd = 1, NDD

!       IF( (GMT_CRD_TYPE(1:4) .EQ. "HEXZ") .AND. (nd.NE.NDIR) ) THEN
!             const_ju = 1.42510
!             cg_nod(1, nd) = cg(1, nd)*const_ju
!             cg_nod(2, nd) = cg(2, nd)*const_ju
!             Flag_Hex  = .True. 
!         ELSE 
!             cg_nod(1, nd) = cg(1, nd)
!             cg_nod(2, nd) = cg(2, nd)
!         END IF  



      nd_right = Index_Neib(nd,2)

      do kt = 1, N_TOT   

         nl = nl + 1 ! Left Interface

         k = Numb(kt, nd)
         k_left = Neib(nd,k)
         k_right = Neib(nd_right, k)
      
         hk  = h_xyz(nd, k)
          

         IF(k_left.eq.I_BOUND_NODE) THEN

!              IF (Flag_HEX) THEN
!                  hk_nod = hk*const_hk
!              ELSE 
!                  hk_nod = hk
!              END IF 

              DO n = 1, NG
                   d0(n) = 2.*XS_D(n,k)/hk
                 d0_nod(n) = 2.*d_nod(n, k, 1, nd)/hk
              END DO

              i_right = -1
            CALL HOM_Compute_Surface_Flux_Bound(i_right, NG, &
              cg(1, nd), fg(1,1,nd), d0, Flux(1,k),&
              d0_nod, FluxSurface(1,nl), CurrentSurface(1,nl) )

!          IF(nd.eq.1.AND.k.eq.1) THEN 
!                 WRITE(*,*) 'k=1'
!                 WRITE(*,*) 'HOM Surface Flux =', FluxSurface(:,nl)
!               END IF

         ELSE

         DO n = 1, NG
              d0(n) = 2.*XS_D(n,k)/hk
            d0_nod(n) = 2.*d_nod(n, k, 1, nd)/hk
         END DO


         hk_neib  = h_xyz(nd, k_left)
         DO n = 1, NG
            d0_neib(n)= 2.*XS_D(n,k_left)/hk_neib
            d0_nod_neib(n)= 2.*d_nod(n, k_left, 2, nd)/hk_neib
!            d0_neib(n)= XS_D(n,k_left)
!            d0_nod_neib(n)= d_nod(n, k_left, 2, nd)
         END DO


            CALL HOM_Compute_Surface_Flux(NG, &
              d0, d0_neib, Flux(1,k), Flux(1, k_left),&
              d0_nod, d0_nod_neib,&
              FluxSurface(1,nl), CurrentSurface(1,nl) )

!        IF (  nd ==  1) THEN
!        IF ( k == 4 ) THEN
!      write(*,'(4A20)')  'homog left interface'

!      write(*,'(4A14)')  'fm0(n)', 'flux_surface', 'flux_surface',&
!                    'current'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  Flux(n,k), FluxSurface(n,nl)/Flux(n,k),&
!                   CurrentSurface(n,nl)/Flux(n,k)
!        END DO 
!        END IF 
!        END IF  


         END IF

! Right Face of the Boundary Node

         IF(k_right .EQ. I_BOUND_NODE) THEN

         nl = nl + 1

!              IF (Flag_HEX) THEN
!                  hk_nod = hk*const_hk
!                  hk_nod = hk
!              ELSE 
!                  hk_nod = hk
!              END IF 

              DO n = 1, NG
                   d0(n) = 2.*XS_D(n,k)/hk
                 d0_nod(n) = 2.*d_nod(n, k, 2, nd)/hk
              END DO

         i_right = 1
         CALL HOM_Compute_Surface_Flux_Bound(i_right, NG, &
              cg(2, nd), fg(1,2,nd), d0, Flux(1,k),&
              d0_nod, FluxSurface(1,nl), CurrentSurface(1,nl) )

!        IF (  nd ==  1) THEN
!        IF ( k == 4 ) THEN
!      write(*,'(4A20)')  'homog right interface'
!
!      write(*,'(4A14)')  'fm0(n)', 'flux_surface', 'flux_surface',&
!                    'current'
!      DO n = 1, NG
!          write(*,'(4ES14.5)')  Flux(n,k), FluxSurface(n,nl)/Flux(n,k),&
!                   CurrentSurface(n,nl)/Flux(n,k)
!        END DO 
!        END IF 
!        END IF  


         END IF

      END DO  ! kt
      END DO ! nd                          

      RETURN
      END 

      SUBROUTINE HOM_Compute_Surface_Flux_Bound(i_right, NG,&
              cg, fg, d0, f0, d_nod, Flux_Surface, Current_Surface)
      IMPLICIT NONE
! Input
      INTEGER :: NG, i_right  
! i_right = 1 for right face, i_right = -1 for left face
      REAL    :: cg, fg(NG), d0(NG), f0(NG), d_nod(NG)       
!      REAL    :: hk 
! Output
      REAL    :: flux_surface(NG), current_surface(NG)
! Local
      INTEGER :: n
            
      DO n = 1, NG       
       Flux_Surface(n) =   cg*( d0(n)-D_Nod(n) )*f0(n) /&
                ( fg(n) + cg*( D_Nod(n)+d0(n) ) )
       current_surface(n) =  i_right*fg(n)*( d0(n)-D_Nod(n) )*f0(n) /&
                ( fg(n) + cg*( D_Nod(n)+d0(n) ) )

!       flux_surface(n)=&
!                (2.*cg*(d0(n)-D_Nod(n))*f0(n))/&
!                ( 2*cg*(D_Nod(n)+d0(n))+fg(n)*hk)
      END DO

      RETURN
      END

      SUBROUTINE HOM_Compute_Surface_Flux(NG, &
              d0, d0_neib, f0, f0_neib,&
              d_nod, d_nod_neib, Flux_Surface, Current_Surface)
!     , hk_neib, hk )
      IMPLICIT NONE
! Input
      INTEGER :: NG
      REAL    :: d0(NG), d0_neib(NG), f0(NG), f0_neib(NG), &
         d_nod(NG), d_nod_neib(NG)
!      REAL    :: hk, hk_neib  
! Output
      REAL    :: flux_surface(NG), current_surface(NG)
! Local
      INTEGER :: n
      REAL    :: denom

      DO n = 1, NG
         denom = d0(n) + d0_neib(n) + d_nod(n) + d_nod_neib(n)
         flux_surface(n) = ( ( d0(n) - d_nod(n) )*f0(n) +  &
              ( d0_neib(n) - d_nod_neib(n) )*f0_neib(n) )/denom
! minus sign before expression for the current 
! because we consider the two-node
! problem (k-1,k), not (k, k+1) as in the manual
         current_surface(n) = -&
         ( -( d0(n)*d0_neib(n) - d_nod_neib(n)*d_nod(n) )*&
           ( f0_neib(n) - f0(n) ) &
           -( d0_neib(n)*d_nod(n) - d0(n)*d_nod_neib(n) )*&
           ( f0_neib(n) + f0(n) )  )  /denom

!                denom = ( d0(n) + D_Nod(n) )*hk_neib +&
!                   (d0_neib(n) + D_Nod_neib(n))*hk
! checking     
!                flux_surface(n)=( (d0_neib(n)-D_Nod_neib(n))*&
!                        hk*f0_neib(n) +&
!                       (d0(n)-D_Nod(n))*hk_neib*f0(n) )/denom

      END DO

      RETURN
      END

      SUBROUTINE HOM_Average_Surface_Flux
!=====================================================================!        
! Calculation of the Surface Neutron Flux for the Fine Spatial Mesh   !
!  (c) Slava 5 June 2001 MEPhI                                        !
!=====================================================================!        
      USE GeoMeTry_Faces_Poly, ONLY : FluxSurfacePoly, N_FACES_POLY,&
        CurrentSurfacePoly
      USE GeoMeTry_Faces, ONLY : FluxSurface, CurrentSurface, N_FACES
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Local
      INTEGER :: ns, nd, nl, ix, iz, iy, nlx, nly, nlz, n_in, k, kh, &
        nfz, iz_shift, np, n
      LOGICAL   flag_boundary_face

      iz_shift = 0
      DO nlz = 1, NZR
         DO nly = 1, NYR
            DO  nlx = 1, NXR
               np = npoly(nly,nlx)
               IF ( np /= 0 ) THEN
! X-Direction left Interface
                  nd = 1   
                  ns = Index_Left_Interface_Poly(np, nlz, nd)

                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO

                  DO iz = 1, NPZ(nlz)
                     DO iy = 1, NPY(nly)
                        n_in = 1 + (iy -1)*NPX(nlx)
                        kh = poly_out(np, n_in)
                        nfz  = iz + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd)
                        DO n = 1, NG
                      FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
                        CurrentSurfacePoly(n, ns)=&
                      CurrentSurfacePoly(n, ns)+CurrentSurface(n, nl)
                        END DO
                     END DO ! iy
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPY(nly) )  
                  CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPY(nly) )  
                  END DO
! Checking that the right face at the boundary
                  flag_boundary_face = .False.
                  IF( nlx.EQ.NXR ) THEN
                            flag_boundary_face = .True.   
                  ELSE IF(  npoly(nly, nlx+1).EQ.0 ) THEN
                            flag_boundary_face = .True.   
                  END IF
!        
                  IF(flag_boundary_face) THEN
                     ns = ns + 1                    
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO

                  DO iz = 1, NPZ(nlz)
                     DO iy = 1, NPY(nly)
                        n_in = iy*NPX(nlx)
                        kh = poly_out(np, n_in)
                        nfz  = iz + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd) + 1
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
                  CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                        END DO
                     END DO ! iy
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPY(nly) )  
                CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPY(nly) )  
                  END DO
                  END IF
! END Checking that the right face at the boundary
! END X-Direction left Interface
                  nd = nd + 1  
! Y-Direction
                  ns = Index_Left_Interface_Poly(np, nlz, nd)
                  
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO
                  DO iz = 1, NPZ(nlz)
                     DO ix = 1, NPX(nlx)
! Left interface n_in = ix
                        n_in = ix 
                        kh = poly_out(np, n_in)
                        nfz  = iz + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd)
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
              CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                        END DO
                     END DO ! ix
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPX(nlx) )  
               CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPX(nlx) )  
                  END DO
! Checking that the right face at the boundary

                  flag_boundary_face = .False.
                  IF( nly.EQ.NYR  ) THEN
                            flag_boundary_face = .True.   
                  ELSE IF(  npoly(nly+1, nlx).EQ.0 ) THEN
                            flag_boundary_face = .True.   
                  END IF

                  IF(flag_boundary_face) THEN
                     ns = ns + 1                    
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO

                  DO iz = 1, NPZ(nlz)
                     DO ix = 1, NPX(nlx)
                        n_in = ix + (NPY(nly)-1)*NPX(nlx)
                        kh = poly_out(np, n_in)
                        nfz  = iz + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd) + 1
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
               CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                        END DO
                     END DO ! ix
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPX(nlx) )  
                CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz)*NPX(nlx) )  
                  END DO
!                  IF( ns == 13 ) THEN
!                  write(*,*) 'ns=12', 'FluxSurfacePoly, FluxSurface=',&
!                    FluxSurfacePoly(1, ns), FluxSurface(1, 12)
!                   PAUSE
                   
!                  END IF ! ns == 12
                  END IF ! nly.EQ.NYR .OR. npoly(nly+1, nlx).EQ.0 
! END Checking that the right face at the boundary
       IF (gmt_crd_type(1:4).EQ."HEXZ") THEN
          nd = nd + 1
          ns = Index_Left_Interface_Poly(np, nlz, nd)

                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO
                  DO iz = 1, NPZ(nlz)
                     nfz  = iz + iz_shift
                     k = np + (nfz -1)*NH
                     nl = Index_Left_Interface(k, nd)
                     DO n = 1, NG
                       FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
             CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                     END DO
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz))  
           CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz))  
                  END DO
! Checking that the right face at the boundary
                  IF( nly.EQ.NYR  ) THEN
                            flag_boundary_face = .True.   
                  ELSE IF(  nlx.EQ.NXR ) THEN
                            flag_boundary_face = .True.   
                  ELSE IF(  npoly(nly+1, nlx+1).EQ.0 ) THEN
                            flag_boundary_face = .True.   
                  END IF

                  IF(flag_boundary_face) THEN
                     ns = ns + 1                    
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO

                  DO iz = 1, NPZ(nlz)
                        nfz  = iz + iz_shift
                        k = np + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd) + 1
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
                CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                        END DO
                  END DO ! iz
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPZ(nlz) )  
                 CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPZ(nlz) )  
                  END DO
                  END IF ! right face
      END IF ! (gmt_crd_type(1:4).EQ."HEXZ")

      nd = nd + 1       
      IF( nd.LE. NDD) THEN     

                  ns = Index_Left_Interface_Poly(np, nlz, nd)
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO
                  DO ix = 1, NPX(nlx)
                     DO iy = 1, NPY(nly)
! Left interface n_in = ix
                        n_in = ix + (iy-1)*NPX(nlx) 
                        kh = poly_out(np, n_in)
                        nfz  = 1 + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd)
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
               CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                          CurrentSurface(n, nl)
                        END DO
                     END DO ! iy
                  END DO ! ix
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPY(nly)*NPX(nlx) )  
              CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPY(nly)*NPX(nlx) )  
                  END DO
                  IF( nlz.eq.NZR) THEN
! Checking that the right face at the boundary
                  ns = ns + 1                    
                  DO n = 1, NG
                     FluxSurfacePoly(n, ns) = 0.
                     CurrentSurfacePoly(n, ns) = 0.
                  END DO
                  DO ix = 1, NPX(nlx)
                     DO iy = 1, NPY(nly)
                        n_in = ix + (iy-1)*NPX(nlx)
                        kh = poly_out(np, n_in)
                        nfz  = NPZ(nlz) + iz_shift
                        k = kh + (nfz -1)*NH
                        nl = Index_Left_Interface(k, nd) + 1
                        DO n = 1, NG
                        FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)+&
                           FluxSurface(n, nl)
                CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)+&
                           CurrentSurface(n, nl)
                        END DO
                     END DO ! ix
                  END DO ! iy
                  DO n = 1, NG
                    FluxSurfacePoly(n, ns)=FluxSurfacePoly(n, ns)/&
                       ( NPY(nly)*NPX(nlx) )  
              CurrentSurfacePoly(n, ns)=CurrentSurfacePoly(n, ns)/&
                       ( NPY(nly)*NPX(nlx) )  
                  END DO
                  END IF 

             END IF 
 
                  END IF ! ( np /= 0 )
             END DO ! nlx
           END DO ! nly
           iz_shift = iz_shift + npz(nlz)
          END DO ! nlz

      RETURN
      END 


      SUBROUTINE OUTput_Average_Surface_Flux
         USE GeoMeTry_Faces_Poly, ONLY : FluxSurfacePoly, N_FACES_POLY,&
        GMT_Faces_Poly_Output_Debug
         USE GeoMeTry_Faces,      ONLY : FluxSurface
         USE Index_Map_Faces_Class, ONLY : Index_Map_Faces, n_out_faces
         USE Core_Map_Class, ONLY : Faces_Map, output

         IMPLICIT NONE
         INCLUDE 'sketch.fh'

           CHARACTER*80 Header_Map
         CHARACTER*4 val_fmt
         CHARACTER*8, DIMENSION(:), ALLOCATABLE :: val_char
         INTEGER ns, nly, nlx, nlz, n_max_faces, nd, n, n_shift, ND_RAD
         REAL Scale_Factor 

         IF( DEBUG ) THEN
           CALL GMT_Faces_Poly_Output_Debug(io_unit+1)
         END IF

         n_max_faces = MAXVAL(n_out_faces)
         ALLOCATE ( val_char(0:n_max_faces) ) 
         val_fmt = "A8"
         val_char(0) =""

      WRITE(io_unit,'(A)') &
     "            SURFACE-AVERAGED NEUTRON FLUX"
      CALL OUTput_Write_Separator(io_unit)

      WRITE(io_unit, '(A,/,4x,30E12.5)') &
       "     Average Neutron Flux (1/cm^3)         : ",&
       (dist_flux(0,0,n), n = 1, NG)

      WRITE(io_unit, '(A)') 


      n = 1  
      Scale_Factor =  1./dist_flux(0,0,n)

      IF( gmt_crd_type(1:4).EQ."HEXZ") THEN
         ND_RAD = 3
      ELSE
         ND_RAD = 2
      END IF                    

      DO n = 1, NG

         Scale_Factor =  1./dist_flux(0,0,n)
         n_shift = 0

         DO nd = 1, ND_RAD
         DO nlz = 1, NZR
         WRITE(Header_Map, '(5x, 3(A, I4) )') " Group= ", n,&
            " Axial Layer = ", nlz, "Direction (x-y-(v)-z) =", nd
          DO nly = 1, Index_Map_Faces(nd)%NYR
             DO nlx = 1, Index_Map_Faces(nd)%NXR
             ns =  Index_Map_Faces(nd)%npoly(nly, nlx)
             IF( ns.ne.0) THEN
             WRITE(&
             val_char(ns), '(F8.4)') FluxSurfacePoly(n, ns+n_shift)&
                *Scale_Factor 

             END IF
             END DO ! nlx
          END DO ! nly
          n_shift = n_shift + n_out_faces(nd)


         CALL output(Faces_Map(nd), n_out_faces(nd), val_char, &
             header_map, io_unit, val_fmt, gmt_crd_type )

          END DO ! nlz 
         END DO ! nd

! Axial Direction
         nd = ND_RAD+1

         IF(nd <= NDD) THEN

             DO nlz = 1, NZR+1

         WRITE(Header_Map, '(5x, 3(A, I4) )') " Group= ", n,&
            " Axial Layer = ", nlz, "Direction (x-y-(v)-z) =", nd
          DO nly = 1, Index_Map_Faces(nd)%NYR
             DO nlx = 1, Index_Map_Faces(nd)%NXR
              ns =  Index_Map_Faces(nd)%npoly(nly, nlx)
             IF( ns.ne.0) THEN

             WRITE(&
             val_char(ns), '(F8.4)') &
              FluxSurfacePoly(n, nlz + (ns-1)*(NZR+1)+n_shift)*&
                  Scale_Factor 
             END IF
             END DO ! nlx
          END DO ! nly

         CALL output(Faces_Map(nd), n_out_faces(nd), val_char, &
             header_map, io_unit, val_fmt, gmt_crd_type )

          END DO ! nlz 

         END IF


! END Axial Direction

       END DO ! NG

      DEALLOCATE ( val_char ) 

!      write(*,*) 

      RETURN 
      END

      SUBROUTINE OUTput_Average_Surface_Current
         USE GeoMeTry_Faces_Poly, ONLY : CurrentSurfacePoly, &
         N_FACES_POLY,  GMT_Faces_Poly_Output_Debug
         USE GeoMeTry_Faces,      ONLY : CurrentSurface
         USE Index_Map_Faces_Class, ONLY : Index_Map_Faces, n_out_faces
         USE Core_Map_Class, ONLY : Faces_Map, output

         IMPLICIT NONE
         INCLUDE 'sketch.fh'

           CHARACTER*80 Header_Map
         CHARACTER*4 val_fmt
         CHARACTER*8, DIMENSION(:), ALLOCATABLE :: val_char
         INTEGER ns, nly, nlx, nlz, n_max_faces, nd, n, n_shift, ND_RAD
         REAL Scale_Factor 

         IF( DEBUG ) THEN
           CALL GMT_Faces_Poly_Output_Debug(io_unit+1)
         END IF

         n_max_faces = MAXVAL(n_out_faces)
         ALLOCATE ( val_char(0:n_max_faces) ) 
         val_fmt = "A8"
         val_char(0) =""

      WRITE(io_unit,'(A)') &
     "            SURFACE-AVERAGED NEUTRON CURRENT * 100"
      CALL OUTput_Write_Separator(io_unit)

      WRITE(io_unit, '(A,/,4x,30E12.5)') &
       "     Average Neutron Flux (1/cm^3)         : ",&
       (dist_flux(0,0,n), n = 1, NG)

      WRITE(io_unit, '(A)') 


      n = 1  
      Scale_Factor =  1./dist_flux(0,0,n)

      IF( gmt_crd_type(1:4).EQ."HEXZ") THEN
         ND_RAD = 3
      ELSE
         ND_RAD = 2
      END IF                    

      DO n = 1, NG

         Scale_Factor =  1./dist_flux(0,0,n)
         n_shift = 0

         DO nd = 1, ND_RAD
         DO nlz = 1, NZR
         WRITE(Header_Map, '(5x, 3(A, I4) )') " Group= ", n,&
            " Axial Layer = ", nlz, "Direction (x-y-(v)-z) =", nd
          DO nly = 1, Index_Map_Faces(nd)%NYR
             DO nlx = 1, Index_Map_Faces(nd)%NXR
             ns =  Index_Map_Faces(nd)%npoly(nly, nlx)
             IF( ns.ne.0) THEN
             WRITE(&
             val_char(ns), '(F8.4)') CurrentSurfacePoly(n, ns+n_shift)&
                *Scale_Factor*100 

             END IF
             END DO ! nlx
          END DO ! nly
          n_shift = n_shift + n_out_faces(nd)


         CALL output(Faces_Map(nd), n_out_faces(nd), val_char, &
             header_map, io_unit, val_fmt, gmt_crd_type )

          END DO ! nlz 
         END DO ! nd

! Axial Direction
         nd = ND_RAD+1

         IF(nd <= NDD) THEN

             DO nlz = 1, NZR+1

         WRITE(Header_Map, '(5x, 3(A, I4) )') " Group= ", n,&
            " Axial Layer = ", nlz, "Direction (x-y-(v)-z) =", nd
          DO nly = 1, Index_Map_Faces(nd)%NYR
             DO nlx = 1, Index_Map_Faces(nd)%NXR
              ns =  Index_Map_Faces(nd)%npoly(nly, nlx)
             IF( ns.ne.0) THEN

             WRITE(&
             val_char(ns), '(F8.4)') &
              CurrentSurfacePoly(n, nlz + (ns-1)*(NZR+1)+n_shift)*&
                  Scale_Factor 
             END IF
             END DO ! nlx
          END DO ! nly

         CALL output(Faces_Map(nd), n_out_faces(nd), val_char, &
             header_map, io_unit, val_fmt, gmt_crd_type )

          END DO ! nlz 

         END IF


! END Axial Direction

       END DO ! NG

      DEALLOCATE ( val_char ) 

!      write(*,*) 

      RETURN 
      END