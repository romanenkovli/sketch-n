      MODULE GeoMeTry_Faces

      INTEGER :: N_FACES

      REAL, DIMENSION(:,:), ALLOCATABLE :: MAT_FD, MAT_Nod 

      REAL, DIMENSION(:,:), ALLOCATABLE :: FluxSurface, &
              CurrentSurface   


      CONTAINS

      SUBROUTINE GMT_Allocate_Faces(NX, NY, NZ, N_POLY, NH,&
        NCHM, Flag_HEX, NG, NXR, NYR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, NY, NZ, N_POLY, NCHM, NH, NG, &
                            NXR, NYR
      LOGICAL, INTENT(IN) :: Flag_HEX

      INTEGER :: N_FACES_HEX_XYV, N_FACES_XYV, NP_SIDE, NVR

      IF( .NOT.Flag_HEX ) THEN 

             N_FACES_XYV = NH + NX + NH + NY

      ELSE 

             NVR = NXR + NYR - 1

              N_FACES_HEX_XYV = N_POLY + NXR+ N_POLY + NYR +&
              N_POLY + NVR
              
             IF(NCHM.EQ.1) THEN
                N_FACES_XYV = N_FACES_HEX_XYV
             ELSE
                NP_SIDE = NINT( SQRT(REAL(NCHM/6)))
                N_FACES_XYV  =  N_FACES_HEX_XYV*NP_SIDE +&
               6*N_POLY*NP_SIDE + 6*3*NP_SIDE*(NP_SIDE-1)/2
              END IF
      END IF
!      WRITE(*,*) 'N_FACES_HEX_XYV=' , N_FACES_HEX_XYV
!      WRITE(*,*) 'NP_SIDE =',NP_SIDE, 'N_FACES_XYV =', N_FACES_XYV
!      PAUSE
      N_FACES = N_FACES_XYV*NZ + NH*(NZ+1) 

!         WRITE(*,*) 'N_FACES =', N_FACES


            ALLOCATE( MAT_FD(NG, N_FACES), MAT_Nod(NG, N_FACES) ) 

            MAT_FD(:, :) = 0.
            MAT_Nod(:, :) =0.

            ALLOCATE( FluxSurface(NG, N_FACES) ) 
            ALLOCATE( CurrentSurface(NG, N_FACES) ) 

            FluxSurface(:, :) =0.
            CurrentSurface(:, :) =0.

      RETURN
      END SUBROUTINE GMT_Allocate_Faces

      END MODULE GeoMeTry_Faces


      MODULE GeoMeTry_Faces_Poly

      INTEGER :: N_FACES_POLY

      REAL, DIMENSION(:,:), ALLOCATABLE :: FluxSurfacePoly, &
              CurrentSurfacePoly 


      CONTAINS

      SUBROUTINE GMT_Allocate_Faces_Poly(NXR, NYR, NZR, N_POLY, &
         Flag_HEX, NG )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NXR, NYR, NZR, N_POLY,  NG
      LOGICAL, INTENT(IN) :: Flag_HEX

      INTEGER ::  N_FACES_XYV, NVR


      IF( .NOT.Flag_HEX ) THEN 

             N_FACES_XYV = 2*N_POLY + NXR +  NYR

      ELSE 

             NVR = NXR + NYR - 1

              N_FACES_XYV = N_POLY + NXR+ N_POLY + NYR +&
              N_POLY + NVR
              
      END IF
!      WRITE(*,*) 'N_FACES_HEX_XYV=' , N_FACES_HEX_XYV
!      WRITE(*,*) 'NP_SIDE =',NP_SIDE, 'N_FACES_XYV =', N_FACES_XYV
!      PAUSE
      N_FACES_POLY = N_FACES_XYV*NZR + N_POLY*(NZR+1) 

!         WRITE(*,*) 'N_FACES_POLY =', N_FACES_POLY

            ALLOCATE( FluxSurfacePoly(NG, N_FACES_POLY) ) 
            ALLOCATE( CurrentSurfacePoly(NG, N_FACES_POLY) ) 

            FluxSurfacePoly(:, :) =0.
            CurrentSurfacePoly(:, :) =0.

      RETURN
      END SUBROUTINE GMT_Allocate_Faces_Poly

      SUBROUTINE GMT_Faces_Poly_Output_Debug(io)
      USE Index_Map_Faces_Class, ONLY : Index_Map_Faces, n_out_faces
      INCLUDE 'sketch.fh' ! USE : gmt_crd_type(1:4), io_unit
      INTEGER :: io, ND_RAD, nlx, nly, nd
        
      IF( gmt_crd_type(1:4).EQ."HEXZ") THEN
         ND_RAD = 3
      ELSE
         ND_RAD = 2
      END IF                    

       OPEN(io, file ='Output_Debug/Faces_Poly.DAT')
          DO nd = 1, ND_RAD
           write(io, *) 'nd=', nd
          DO nly = 1, Index_Map_Faces(nd)%NYR
           write(io, '(A,I5)') 'nl =', nly
           write(io, '(10I5)')&
           (Index_Map_Faces(nd)%npoly(nly, nlx), &
            nlx=1, Index_Map_Faces(nd)%NXR )
           END DO
          END DO
       CLOSE(io)

      END SUBROUTINE GMT_Faces_Poly_Output_Debug

      END MODULE GeoMeTry_Faces_Poly