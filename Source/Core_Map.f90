      MODULE Index_Map_Class
            IMPLICIT NONE
         PRIVATE :: new_index_map
         TYPE Index_Map
             INTEGER :: NXR, NYR
             INTEGER, DIMENSION(:,:), POINTER :: npoly
         END TYPE Index_Map
         INTERFACE new
             MODULE PROCEDURE new_index_map
         END INTERFACE new
      CONTAINS
      
      SUBROUTINE new_index_map(NXR, NYR, this)
      INTEGER, INTENT(IN)  :: NXR, NYR
      TYPE (Index_Map), INTENT(OUT) :: this
      INTEGER nly, nlx

      this%nxr = NXR
      this%nyr = NYR
      ALLOCATE ( this%npoly(NYR, NXR) ) 
        DO nly = 1, NYR
           DO nlx = 1, NXR
              this%npoly(nly, nlx) = 0
           END DO
      END DO 


      RETURN
      END SUBROUTINE new_index_map

      END MODULE Index_Map_Class
!=====================================================================!
      MODULE Index_Map_Faces_Class
         USE Index_Map_Class
            IMPLICIT NONE
         include 'sketch.fh'
         
         TYPE(Index_Map), DIMENSION(NDIR) :: Index_Map_Faces
         INTEGER, DIMENSION(NDIR)         :: n_out_faces                      

      CONTAINS

      SUBROUTINE Set_Index_Map_Faces(la_poly)
      IMPLICIT NONE
      include 'sketch.fh'
      INTEGER, INTENT(IN) :: la_poly(N_POLY, NE+1) 


      INTEGER nlx, nly, np, nd, k, nom

      nd = 1
      n_out_faces(nd) = 0

      CALL new(NXR+1, NYR, Index_Map_Faces(nd) )

         DO nly = 1, NYR
            DO nlx = 1, NXR
               np = npoly(nly, nlx)
               IF ( np /= 0 ) THEN
                 n_out_faces(nd) = n_out_faces(nd) + 1
                 Index_Map_Faces(nd)%npoly(nly, nlx) = n_out_faces(nd)
                 IF( nlx == NXR ) THEN
                  n_out_faces(nd)=n_out_faces(nd)+1
                  Index_Map_Faces(nd)%npoly(nly, nlx+1)=n_out_faces(nd)
                 ELSE IF ( npoly(nly, nlx+1) == 0) THEN
                  n_out_faces(nd)=n_out_faces(nd)+1
                  Index_Map_Faces(nd)%npoly(nly, nlx+1)=n_out_faces(nd)
                 END IF
                END IF
            END DO
          END DO

! Y Direction
      nd = nd + 1   
      n_out_faces(nd)       = 0

      CALL new(NXR, NYR+1, Index_Map_Faces(nd) )

        DO nlx = 1, NXR
            DO nly = 1, NYR
               np = npoly(nly, nlx)
               IF ( np /= 0 ) THEN
                   n_out_faces(nd) = n_out_faces(nd)+1
                 Index_Map_Faces(nd)%npoly(nly, nlx) = n_out_faces(nd) 
                 IF( nly == NYR ) THEN
                   n_out_faces(nd) = n_out_faces(nd)+1
                  Index_Map_Faces(nd)%npoly(nly+1,nlx)=n_out_faces(nd) 
                 ELSE IF ( npoly(nly+1, nlx) == 0) THEN
                   n_out_faces(nd) = n_out_faces(nd)+1
                 Index_Map_Faces(nd)%npoly(nly+1,nlx)=n_out_faces(nd) 
               END IF
               END IF
            END DO
          END DO
! v-direction for Hexagonal Mesh
      IF (gmt_crd_type(1:4).EQ."HEXZ") THEN
         nd = nd + 1
         n_out_faces(nd) = 0

           CALL new(NXR+1, NYR+1, Index_Map_Faces(nd) )

          DO k = 1, N_POLY
             nom = k
! we start the line if there is no left neighbour
            IF(la_poly(nom,3).eq.I_BOUND_NODE) THEN
              n_out_faces(nd) = n_out_faces(nd)+1
              nlx = N_Coord(nom, 1)
              nly = N_Coord(nom, 2)
              Index_Map_Faces(nd)%npoly(nly, nlx) = n_out_faces(nd) 
! do while there is the right neighbour
              DO while( la_poly(nom, 3 + NE/2) .NE. I_BOUND_NODE ) 
                n_out_faces(nd) = n_out_faces(nd)+1
                nom = la_poly(nom, 3 + NE/2)
                  nlx = N_Coord(nom, 1)
                nly = N_Coord(nom, 2)
                Index_Map_Faces(nd)%npoly(nly, nlx) = n_out_faces(nd) 
!!              write(*,*) 'nom, ns =', nom, ns
              END DO ! while
              n_out_faces(nd) = n_out_faces(nd)+1
              Index_Map_Faces(nd)%npoly(nly+1, nlx+1)=n_out_faces(nd) 
            END IF ! la_poly(k,3) 
          END DO ! N_POLY
      END IF ! gmt_crd_type == hexz 

! Z-direction
        nd = NDIR

        n_out_faces(nd)       = N_POLY

         CALL new(NXR, NYR, Index_Map_Faces(nd) )

        DO nlx = 1, NXR
            DO nly = 1, NYR
             Index_Map_Faces(nd)%npoly(nly,nlx)=npoly(nly, nlx)
            END DO
         END DO

!      write(*,*) 'n_out_faces =', n_out_faces
!            N_FACES = NZ*( (NDIR-1)*NH + NX + NY + NV )  + (NZ+1)*NH
!        write(*,*) N_POLY + NY, N_POLY + NX, N_POLY + NX, N_POLY
!   pause


      RETURN
      END SUBROUTINE Set_Index_Map_Faces

      END MODULE Index_Map_Faces_Class

!======================================================================!
! 
      MODULE Core_Map_Class

      USE Index_Map_Class
      IMPLICIT NONE
      include 'sketch.fh'

      PRIVATE :: init_core_map, output_core_map


      TYPE Core_Map
         PRIVATE
         TYPE(Index_Map) :: Core_Index 
         INTEGER :: NXR_Beg_Min, NXR_Max
         INTEGER, DIMENSION(:), POINTER :: NXR_Beg, NXR_End
         INTEGER :: NYR_Beg, NYR_End
      END TYPE Core_Map

      TYPE(Core_Map), DIMENSION(NDIR)  :: Faces_Map

      
      INTERFACE new
         MODULE PROCEDURE init_core_Map
      END INTERFACE new
                          
      INTERFACE output
         MODULE PROCEDURE output_core_Map
      END INTERFACE output

      CONTAINS

      SUBROUTINE set_faces_map
      USE Index_Map_Faces_Class, ONLY : Index_Map_Faces
      INTEGER nd
      
             DO nd = 1, NDIR

         CALL new( Index_Map_Faces(nd), Faces_Map(nd) )

        END DO

      RETURN
      END SUBROUTINE set_faces_map
      
      SUBROUTINE init_core_map(Core_Index, this)
         TYPE(Index_Map), INTENT(IN) :: Core_Index
         TYPE(Core_Map),  INTENT(OUT) :: this
! local
         INTEGER :: nlx, nly

         CALL new( core_index%NXR, core_index%NYR, this%core_index )
         this%core_index%nxr = core_index%NXR 
         this%core_index%nyr = core_index%NYR
         ALLOCATE ( this%nxr_beg(core_index%NYR), &
                 this%nxr_end(core_index%NYR) )
         DO nly = 1, core_index%NYR
            DO nlx = 1, core_index%NXR
               this%core_index%npoly(nly, nlx) = &
                   core_index%npoly(nly, nlx)
            END DO
         END DO
         CALL   GMT_Reactor_Numbering(core_index%NXR, core_index%NYR, &
          core_index%npoly, &
          this%Nxr_Beg, this%Nxr_End, this%Nyr_Beg, this%Nyr_End, &
              this%Nxr_Max, this%Nxr_Beg_Min)
      RETURN
      END SUBROUTINE init_core_map

      SUBROUTINE output_core_map( this, ndata, value, header_map,&
        io_out, val_fmt, gmt_crd_type )

      INTEGER, INTENT(IN)        :: ndata, io_out
      TYPE(Core_Map), INTENT(IN) :: this
      CHARACTER*(*), INTENT(IN)  :: value(0:ndata), val_fmt, header_map,&
                                   gmt_crd_type  

! Local Variables      
      INTEGER  N_Position, N_Blank, nly, nlx, N_Digits, N_Blank_HEX
      CHARACTER*5 C_Blank, C_Position
      CHARACTER Char_Digits
      CHARACTER*100 fmt
      
      integer n_hex_0 , nx_blank

!      CALL OUTput_Write_Separator(io_out)
      WRITE(io_out,'(/,A)') Header_Map
!      CALL OUTput_Write_Separator(io_out)

      READ(val_fmt, '(1X, I1)') N_Digits
      N_Position = this%NXR_Max - this%NXR_Beg_Min + 1

         WRITE(C_Position, '(I5)') N_Position
         WRITE(Char_Digits, '(I1)') N_Digits

          IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
           N_Blank_HEX = (this%Nyr_End - this%Nyr_Beg)*N_Digits/2
           n_hex_0 = 1000
           DO nly = this%NYR_BEG, this%NYR_END
              nx_blank = (this%Nyr_End - nly)*N_Digits/2 +&
                  (this%NXR_Beg(nly)-this%NXR_Beg_Min)*N_Digits
              n_hex_0 = min(n_hex_0, nx_blank)
           END DO
           N_Blank_HEX = N_Blank_HEX - n_hex_0
        ELSE
           N_Blank_HEX = 0
        END IF 
        WRITE(C_Blank, '(I5)') N_Blank_HEX

         IF (  N_Blank_HEX /= 0 ) THEN
         fmt = '(/,5x,'//C_Blank//'x,'//C_Position//&
         'I'//Char_Digits//')'
         ELSE

         fmt = '(/,5x,'//C_Position//&
         'I'//Char_Digits//')'
         END IF 
      

!         WRITE(*,'(A)') fmt&
!                (nlx,  nlx = this%Nxr_Beg_Min, this%Nxr_Max)
!         pause


         WRITE(io_out,fmt) &
                (nlx,  nlx = this%Nxr_Beg_Min, this%Nxr_Max)

          

         WRITE(Char_Digits, '(I1)') N_Digits - 1

         IF (  N_Blank_HEX /= 0 ) THEN
         fmt = '(5x,'//C_Blank//'x,'//C_Position//'(1x,'//&
          Char_Digits//'("-")))'
         ELSE
         fmt = '(5x,'//C_Position//'(1x,'//&
          Char_Digits//'("-")))'
         END IF
  
!         WRITE(*,'(A)') fmt&
!                (nlx,  nlx = this%Nxr_Beg_Min, this%Nxr_Max)
!         pause

         WRITE(io_out, fmt) 
 
         DO nly = this%Nyr_Beg, this%Nyr_End
            N_Blank = this%Nxr_Beg(nly) - this%Nxr_Beg_Min
            N_Position = this%Nxr_End(nly) - this%Nxr_Beg(nly) + 1
              IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
              N_Blank_HEX = (this%Nyr_End - nly)*N_Digits/2 - n_hex_0
           ELSE
              N_Blank_HEX = 0
           END IF 

           WRITE(C_Position, '(I5)') N_Position
           WRITE(C_Blank, '(I5)') N_Blank*N_Digits+1+N_Blank_HEX

            fmt ='(I3,":",'//C_Blank//'x,'//C_Position//val_fmt//')'
            WRITE(io_out, fmt)&
                nly,(value(this%core_index%npoly(nly,nlx)),&
                     nlx=this%NXR_Beg(nly), this%NXR_End(nly))

         END DO ! nly

      WRITE(io_out,*)
!      CALL OUTput_Write_Separator(io_out)

      RETURN
      END SUBROUTINE output_core_map
      
      END MODULE Core_Map_Class                    

