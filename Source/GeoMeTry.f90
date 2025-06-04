      subroutine GMT_Input
!=====================================================================*
!        Input from the file FILE_INPUT                                *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Local Variables

      integer i, j, n,  np, ios, nd,  nly

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      character*200 Message
      logical error_find

!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'


      open(io_unit,file = FILE_INPUT,status='old', iostat=ios)

!      call Iostat_Error_Check&
!      (ios,"Could not find the input FILE_INPUT file "//FILE_INPUT)

! reading the FILE_INPUT Header
!      do line = 1, N_LINE_PROBLEM_TITLE
!         read(io_unit,'(A80)') PROBLEM_TITLE(line)
!      end do

!     Reading Geometry Type (XYZ or HEXZ
      REWIND(io_unit)
      CALL MSC_Search_Header_In_File(io_unit, &
        "GMT_CRD_TYPE", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',&
           " identifier GMT_CRD_TYPE is not found " //&
           "in the FILE_INPUT, Cartesian XYZ geomertry "//&
           "is assumed")

           GMT_CRD_TYPE(1:3) = "XYZ"

      else

      read(io_unit,fmt=*,iostat=ios) GMT_CRD_TYPE
      call Iostat_Error_Check&
      (ios,"Error in reading the geometry coordinate type,&
       under identifier GMT_CRD_TYPE from the FILE_INPUT")

      end if

! reading GMT_NUM_BNDL
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_NUM_BNDL", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier GMT_NUM_BNDL in the FILE_INPUT,  '//&
         'bundle (assembly) numbers are not specified')
      else
!        IF(gmt_crd_type(1:3).EQ."XYZ") THEN
!         read(io_unit,  fmt=*, iostat=ios) &
!                   ((npoly(i,j),j=1,NXR),i=1,NYR)
!        ELSE IF (gmt_crd_type(1:4).EQ."HEXZ") THEN
           DO nly = 1, NYR
              read(io_unit,  fmt=*, iostat=ios) &
                   n_b(nly), n_e(nly) 
!             write(*,*) 'nly, n_b, n_e, NYR =', nly, n_b(nly), &
!              n_e(nly) , NYR
           END DO
!        END IF ! gmt_crd_type
! setting up npoly
         call Iostat_Error_Check&
      (ios,"Error in reading the geometry array npoly(NXR,NYR),&
       under identifier GMT_NUM_BNDL from the FILE_INPUT") 
      end if

!      read(io_unit,fmt=*,iostat=ios) ((npoly(i,j),j=1,NXR),i=1,NYR)

      IF(gmt_crd_type(1:3).EQ."XYZ") THEN
! reading GMT_MSH_XDIR
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_XDIR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier GMT_MESH_XDIR in the FILE_INPUT,  '//&
         'spatial mesh in X direction is not given')
      else
      read(io_unit,fmt=*,iostat=ios) (npx(i),hx(i),i=1,NXR)

      call Iostat_Error_Check&
      (ios,"Error in reading the geometry arrays npx(NXR),hx(NXR),&
       under identifier GMT_MESH_XDIR from the FILE_INPUT")
      end if

! reading GMT_MSH_YDIR
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_YDIR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier GMT_MSH_YDIR in the FILE_INPUT,  '//&
         'spatial mesh in Y direction is not given')
      else
      read(io_unit,fmt=*,iostat=ios) (npy(i),hy(i),i=1,NYR)

      call Iostat_Error_Check&
      (ios,"Error in reading the geometry arrays npy(NYR),hy(NYR),&
       under identifier GMT_MSH_YDIR from the FILE_INPUT")
      end if

      NV = 0 ! Third Direction used only for HEXagonal geometry

! END OF Cartesian coordinates
      ELSE IF (gmt_crd_type(1:4).EQ."HEXZ") THEN
! Hexagonal geometry
      rewind(io_unit)
! GMT_MSH_RDIR
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_RDIR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier GMT_MSH_RDIR in the FILE_INPUT,  '//&
         'HEXagonal spatial mesh is not given')
      else
      read(io_unit,fmt=*,iostat=ios) hr_hex
      call Iostat_Error_Check&
      (ios,"Error in reading the spatial mesh HR_HEX&
       under identifier MSH_RDIR from the FILE_INPUT")
      end if
! "GMT_MSH_NV"
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_NV", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier GMT_MSH_NV in the FILE_INPUT,  '//&
         'HEXagonal spatial mesh is not given, '//&
         ' Set NV to NX (OK for Symmetric Core')
          NV = NX
      else
      read(io_unit,fmt=*,iostat=ios) NV
      call Iostat_Error_Check&
      (ios,"Error in reading the dimension NV" //&
       "under identifier GMT_MSH_NV from the FILE_INPUT")
      end if 
      END IF ! (gmt_crd_type(1:4).EQ."HEXZ")
! reading GMT_MSH_ZDIR
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_ZDIR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         
! in the case of 2D calculations we do not need the axial mesh,
! no warning
         IF( (NDD.EQ.2.AND. gmt_crd_type(1:3).EQ."XYZ").OR.&
            (NDD.EQ.3.AND. gmt_crd_type(1:4).EQ."HEXZ") ) THEN
             DO i = 1, NZR
               npz(i) = 1
               hz(i) = 1.0
             END DO
         ELSE 
         WRITE(*,*)  'gmt_crd_type(1:3).EQ."XYZ")', &
                gmt_crd_type(1:3).EQ."XYZ" 
         WRITE(*,*) 'NDD =', NDD
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier GMT_MSH_ZDIR in the FILE_INPUT,  '//&
         'spatial mesh in Z direction is not given')
         END IF

      else
      read(io_unit,fmt=*,iostat=ios) (npz(i),hz(i),i=1,NZR)
      call Iostat_Error_Check&
      (ios,"Error in reading the geometry arrays npz(NZR),hz(NZR),&
       under identifier GMT_MSH_ZDIR from the FILE_INPUT")
      end if
     

! reading GMT_MSH_ZNEW
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_MSH_ZNEW", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
!         call MSC_ERR_Add_Message('WARNING',&
!           " identifier GMT_MSH_ZNEW is not found " //&
!           "in the FILE_INPUT, quasy-uniform spatial mesh in Z "//&
!           "direction is used")

           iz_non_equal = 0

      else
      read(io_unit,fmt=*,iostat=ios) (hzt(n),n=1,NZ)
      call Iostat_Error_Check&
      (ios,"Error in reading the geometry arrays hzt(NZ),&
       under identifier GMT_MSH_ZNEW from the FILE_INPUT")

        iz_non_equal = 1

      end if

! reading GMT_NUM_CORE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_NUM_CORE", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier GMT_NUM_CORE in the FILE_INPUT,  '//&
         'core geometry is defined by the fission materials')
          flag_core_boundary = .False.
      else
         read(io_unit,  fmt=*, iostat=ios) &
                   (n_core(np), np = 1, N_POLY)
! setting up npoly
         call Iostat_Error_Check&
      (ios,"Error in reading the core boundarr n_core(N_POLY),"//&
       "under identifier GMT_NUM_CORE from the FILE_INPUT") 
         flag_core_boundary = .True.
      end if

! reading GMT_COR_LOAD
!      write(*,*) 'GMT_COR_LOAD'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_COR_LOAD", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',&
           "identifier GMT_COR_LOAD is not found " //&
           "in the FILE_INPUT, Core Loading is not specified")
      else
       read(io_unit,  fmt=*, iostat=ios) &
         ( Core_Load(np), np = 1, N_POLY )

       call Iostat_Error_Check(ios,"Error Reading Core Bundle&
      Loading Core_Load(1:N_POLY) under identifier GMT_COR_LOAD "&
       //" from the FILE_INPUT file")

       do i = 1, N_Bundle_Type
         read(io_unit,  fmt=*, iostat=ios) &
                             (Bundle_Compos(n, i), n = 1, NZR)
         if(ios .NE. 0) then
              write(*,*) 'i=', i
              write(Message, '(A,I5,A)')&
       "Error Reading Bundle Compositions, Bundle type  =", i,&
       " under identifier GMT_COR_LOAD "
              call Iostat_Error_Check(ios, Message)
         end if
       end do

       do n = 1, NZR
        do np = 1, N_POLY
         l(np, n) = Bundle_Compos(n, Core_Load(np) )
        end do
       end do
      
      end if

!GMT_BND_COND
!      write(*,*) 'GMT_BND_COND'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "GMT_BND_COND", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',&
           "identifier GMT_BND_COND is not found " //&
           "in the FILE_INPUT")
      else

       read(io_unit,fmt=*, iostat=ios) ((i_dr(i,nd),i=1,2),nd=1,NDIR-1)

       call Iostat_Error_Check(ios,"Error Reading Boundrary Condition"&
       // "type I_DR(1:2,1:2) from  the FILE_INPUT file " )

! 
       read(io_unit,  fmt=*, iostat=ios) &
                 (((dr(n,i,nd),n=1,NG),i=1,2),nd=1,NDIR-1)

       call Iostat_Error_Check(ios,"Error Reading Boundrary Condition "&
       // "values DR(1:NG, 1:2,2) from  the FILE_INPUT file " )


       if (NZ.eq.1) then

        read(io_unit,  fmt=*, iostat=ios) (b2(n),n=1,NG)


        call Iostat_Error_Check(ios,"Error Reading Buckling b2(NG)&
       from  the FILE_INPUT file " )

           DO i = 1, 2
              i_dr(i, NDIR) = 1
              DO n = 1, NG
                 dr(n,i, NDIR) = 0.
              END DO
            END DO

       else
         read(io_unit,  fmt=*, iostat=ios) (i_dr(i,NDIR),i=1,2)
         call Iostat_Error_Check(ios,"Error Reading Axial " //&
     "Boundrary Condition type I_DR(1:2,3) from  the FILE_SECT file ")
         read(io_unit,fmt=*, iostat=ios) ((dr(n,i,NDIR),n=1,NG),i=1,2)
         call Iostat_Error_Check(ios,"Error Reading Axial " // &
     "Boundrary Condition Values I_DR(1:2,3) from the FILE_SECT file")
       end if
      end if

      close(io_unit)

!      write(*,*) 'Reading Done'
      

      return 
      end

      subroutine GMT_Set
!=====================================================================*
!        Computing Gemetry Parameters of the Reactor Model            *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Boundary  
      USE GeoMeTry_Faces
      USE GeoMeTry_Faces_Poly
      USE GeoMeTry_Triangle
      implicit none
      include 'sketch.fh'
      LOGICAL Flag_HEX

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Numbering_XY"
      IF(gmt_crd_type(1:3).EQ."XYZ") THEN
         CALL GMT_Set_Numbering_XY
      ELSE IF(gmt_crd_type(1:4).EQ."HEXZ") THEN
         CALL GMT_Set_Numbering_HEX
      END IF

      if(DEBUG) write(*,'(A)') "Geometry: GMT_DBG_Check_Neiboughrs"
      if(DEBUG) call GMT_DBG_Check_Neiboughrs

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Spatial_Mesh"
      call GMT_Set_Spatial_Mesh

! NOT NEEDED    
!      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Red_Black_Numbering"
!      call GMT_Set_Red_Black_Numbering

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Numbering_3D"
      call GMT_Set_Numbering_3D

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Assembly_Volume"
      call GMT_Set_Assembly_Volume

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_OUT_Core_Numbering"
      call GMT_Set_OUT_Core_Numbering

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_3D_Neighbours "
      call GMT_Set_3D_Neighbours
      if(DEBUG) write(*,'(A)') "Geometry: GMT_Allocate_Faces"
      
      IF( gmt_crd_type(1:3).EQ."XYZ") THEN
         Flag_HEX=.False.
      ELSE
         Flag_HEX=.True. 
      END IF  

      CALL GMT_Allocate_Faces(NX, NY, NZ, N_POLY, NH,&
        NCHM, Flag_HEX, NG, NXR, NYR)
      CALL GMT_Allocate_Faces_Poly(NXR, NYR, NZR, N_POLY, &
         Flag_HEX, NG )

      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Boundary"
      CALL GMT_Set_Boundary
      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Boundary Conditions"
      CALL GMT_Set_Boundary_Conditions


      if(DEBUG) write(*,'(A)') "Geometry: GMT_Set_Index_Left_Interface"
      IF(gmt_crd_type(1:4).EQ."HEXZ".AND.NCHM.NE.1) THEN 
        CALL GMT_Set_Index_Left_Interface(la_hex)
      ELSE
        CALL GMT_Set_Index_Left_Interface(la)
      END IF

! renumbering of the boundary nodes in the case of the repeated structures
      CALL GMT_Set_3D_Neighbours_Repeated_Structure

      if(DEBUG) write(*,'(A)') "Geometry: GMT_DBG_Output"
      if(Debug) call GMT_DBG_Output



      return
      end


      subroutine GMT_Set_Numbering_HEX
!=====================================================================*
!        Numbering of the Nodes in Hexagonal Plane                    *
!                   Vyachreslav Zimin (c) 21 November 2000            *
!               e-mail: slava@ets.mephi.ru                            *
!=====================================================================*
      USE GeoMeTry_Triangle

      implicit none
      include 'sketch.fh'

! Input: n_b(NY), n_e(NY)

! Output NPOLY(NY, NX)       - numbering of the nodes
!          N_Coord(np, 2) 
!          n_fine(NY, NX)
!          poly_out(N_POLY, NCHM)   
! local variables 
        INTEGER nlx, nly, n_in, j, np, k, nf
      CHARACTER*100 Message

      NV = MAX(NX, NY)

      DO nly = 1, NYR   
         DO nlx = 1, NXR
           npoly(nly, nlx) = 0
         END DO
      END DO  

      np = 0
      DO nly = 1, NYR
!         write(*,*) 'nly =', nly, 'n_b(nly), n_e(nly)=', &
!          n_b(nly), n_e(nly)
         DO nlx = n_b(nly), n_e(nly)
             np = np + 1
             npoly(nly, nlx) = np
         END DO
      END DO

      IF(np.NE. N_POLY) then 
         WRITE(Message, '(A, I5, I5)') &
        "total number of the assemblies Computed /= N_POLY ", &
        np, N_POLY
         call MSC_ERR_Add_Message('ERROR',Message)
      END IF

!
      nf=0
      np_out(0) = 0
      DO nly = 1, NYR
         DO nlx = 1, NXR
            np = npoly(nly, nlx)
            IF(np /= 0) THEN
               N_Coord(npoly(nly,nlx),1) = nlx
               N_Coord(npoly(nly,nlx),2) = nly
! n_fine used only in hexagonal geometry with 1 point per hex
               n_fine(nly, nlx) =  npoly(nly, nlx)
            DO n_in = 1, NCHM
               nf = nf + 1
               poly_out(np,n_in) = np
               np_out(nf) = np
            END DO ! n_in
            END IF
         END DO ! nlx
      END DO !nly


! hexagonal geometry  neighbour's number

      IF (NCHM.EQ.1) THEN 
         CALL GMT_Set_LA(NH, NE, NX, NY, n_fine, la, I_BOUND_NODE)
      ELSE 
        CALL Set_Triangle_Local_Mesh(NCHM) 
        CALL Set_LA_Triangle(NCHM, NH, NE, NXR, NYR, npoly, la, &
            I_BOUND_NODE)

         CALL Allocate_Triangle_LA_HEX(NE, N_POLY) 
         CALL GMT_Set_LA(N_POLY, NE, NXR, NYR, npoly, la_hex, &
         I_BOUND_NODE)
!          CALL GMT_Set_Rhomb_Neighbours
      END IF 

      IF( DEBUG ) THEN
           OPEN(io_unit, file ='Output_Debug/LA.dat', &
                          status = 'unknown')
           DO k = 1, NH
                 WRITE(io_unit, '(10I8)' ) (la(k, j), j=1, NE+1) 
           END DO 
           CLOSE(io_unit)

      END IF  

      return
      end

      subroutine GMT_Set_Numbering_XY
!=====================================================================*
!        Numbering of the Nodes in X - Y Plane                        *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Input: NXR, NPX(NXR), NYR, NPY(NYR), NH, NX, NY, N_POLY, NCHM, 
!        NZR, 
!       types l(N_POLY, 1 or 2), sf(NNODE, NG)

! Output NP_OUT(N_POLY,NCHM) - Node Number in assembly (2D)
!        N_FINE(NLY, NLX, 1) - Node Number in X-Y Coordinate
!        NP_OUT(NH) - Assembly Number for the node k
!        LA(NH,5) - Neighbours in X-Y plane { West, North, East, South, k}
!*        NP_AZ -  Number of assemblies in the reactor core ! Not Used
!*        KP_AZ(NP_AZ) - Node Number of the Reactor Core ! Not Used
! local variables 
      integer k, i, j, n_f, nly, nlx, i_out, j_out, n_in, np
      CHARACTER*100 Message

      DO nly = 1, NYR   
         DO nlx = 1, NXR
           npoly(nly, nlx) = 0
         END DO
      END DO  

      np = 0
      DO nly = 1, NYR
!         write(*,*) 'nly =', nly, 'n_b(nly), n_e(nly)=', &
!          n_b(nly), n_e(nly)
         DO nlx = n_b(nly), n_e(nly)
             np = np + 1
             npoly(nly, nlx) = np
         END DO
      END DO

      IF(np.NE. N_POLY) then 
         WRITE(Message, '(A, I5, I5)') &
        "total number of the assemblies Computed /= N_POLY ", &
        np, N_POLY
         call MSC_ERR_Add_Message('ERROR',Message)
      END IF

      np_out(0) = 0

      do k = 1, N_POLY
         do j = 1,NCHM
            poly_out(k,j) = 0.
         end do
      end do

      do i = 1, NYR
         do j = 1, NXR
            if(npoly(i,j).ne.0) then
               N_Coord(npoly(i,j),1) = j
               N_Coord(npoly(i,j),2) = i
            end if
        end do
      end do

! n_f - number of fine-mesh point
           n_f = 0
! nly number of fine-mesh line in y direction 
           nly = 0
! i - line's number of the region in y direction
          do i = 1, NYR
! i_out - internal number of line in region
           do i_out = 1, NPY(i)
           nly = nly + 1
! nlx - number of fine-mesh line in x direction
            nlx = 0
! j - line's number of the region in x direction
            do j = 1, NXR
!              write(*,*) 'i, j =', i, j

               do j_out = 1,NPX(j)
!                  write(*,*) 'i_out, j_out =', i_out, j_out
                   nlx = nlx + 1
                   if(npoly(i,j).ne.0) then
                     np = npoly(i,j)
                     n_f = n_f + 1
!                    write(*,*) 'nlx, nly, n_f', nlx, nly, n_f
!                    pause
                     n_fine(nly, nlx) = n_f
! n_in - number of internal point into the  region
                     n_in = j_out + (i_out-1)*NPX(j)
!                    write(*,*) 'np, n_in', np, n_in

                     poly_out(np,n_in) = n_f

                     np_out(n_f) = np
                 end if
               end do
             end do
          end do
         end do

         CALL GMT_Set_LA(NH, NE, NX, NY, n_fine, la, I_BOUND_NODE)           

      return
      end

      subroutine GMT_DBG_Check_Neiboughrs
!=====================================================================*
!        Check of the Neighbours in LA Array                          *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Output: Checked LA
! Input:  LA(NH,5), NH, NE/2, 
! Local Variables 
      integer k, j, next, i

! cartesian geometry
      do k = 1, NH
          do j = 1, NE/2
          next=la(k,j)
          if(next.ne.I_BOUND_NODE.and.next.GT.0) then
           i = la(next,j+NE/2)
           if(k.ne.i) then
            write(*,*) ' k =',k,'next =',next
            write(*,*) 'i =',i
            write(*,*) '                  ATTENTION !!           '
            write(*,*) '          Error in neigbours calculation !'
            write(*,*) ' Please Check Geometry Input Data - NH, NX, NY'
            stop
           end if
          end if
         end do
         end do

         return
         end



      subroutine GMT_Set_Spatial_Mesh
!=====================================================================*
!  Spatial Mesh Size in X-Y-Z Directions                              *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Triangle
      implicit none
      include 'sketch.fh'

! Input:  NZR, NPZ(NZR), NXR, NPX(NXR), NYR, NPY(NYR), NZ, 
!         NH, iz_non_equal 
! Output: 
!       hzp(NZR) - mesh size inside axial layer in Z direction
!       hzt(NZ) - axial mesh size 
!       ns_out(NZ) - Axial Layer Number
!       hxp(NXR), hyp(NYR)  -  mesh size for X-Y directions
!       hxy(NH,2) - mesh size for the node k in X-Y direction 
!       sxy(NH) - area of the node k
!       volume(N_TOT) - volume of the node k
!       h_xyz(NDIR, N_TOT) - spatial mesh size of the node k 
!       s_xyz(NDIR, N_TOT) - area of the face orthogonal nd direction

! Local Variables
      integer n_in, ns, n, n1, nlx, nly, k, kt, nn
      REAL hex_area, hex_side, constant_triangle_mesh

! Z - direction

      if(iz_non_equal.eq.1) then
! Non-uniform mesh for all Layers
           do ns = 1,NZR
            do n = 1, NPZ(ns)
               n_in = n_in + 1
               ns_out(n_in) = ns
            end do
           end do
       else
! Uniform Axial Mesh inside Axial Layer
           n_in = 0
           do ns = 1,NZR
            hzp(ns) = hz(ns)/npz(ns)
            do  n = 1, NPZ(ns)
               n_in = n_in + 1
               hzt(n_in) = hzp(ns)
               ns_out(n_in) = ns
            end do
           end do
       end if

!       write(*,*) 'ns_out =', ns_out, 'npz=', npz
!       pause

! X - direction
         IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN

            IF(NCHM.eq.1) THEN
               NP_SIDE = 1
               constant_triangle_mesh = 1.
            ELSE
               constant_triangle_mesh = 0.3333333
            END IF

!         WRITE(*,*) NP_SIDE, constant_triangle_mesh
!         pause

            DO ns = 1, NXR
               npx(ns) = NP_SIDE
               hx(ns) = constant_triangle_mesh*HR_HEX
            END DO
            DO ns = 1, NYR
               npy(ns) = NP_SIDE
               hy(ns) = constant_triangle_mesh*HR_HEX
            END DO
              hex_area = 0.5*sqrt(3.)*HR_HEX*HR_HEX
            hex_side = HR_HEX/sqrt(3.)
         END IF

      IF( gmt_crd_type(1:4).EQ."HEXZ" .AND.(NCHM.ne.1) ) THEN
      n_in = 0

      do ns = 1,NXR
         hxp(ns) = hx(ns)/npx(ns)
         do  n = 1, 2*NPX(ns)
            n_in = n_in + 1
            hxy(n_in,1) = hxp(ns)
         end do
      end do

!       write(*,*) 'hxy(n_in,1)=', (hxy(n_in,1), n_in = 1, NX)
!      pause
! Y - Direction
      n_in = 0
      do ns = 1,NYR
         hyp(ns) = hy(ns)/npy(ns)
         do  n = 1, 2*NPY(ns)
            n_in = n_in + 1
            hxy(n_in,2) = hyp(ns)
         end do
      end do
      ELSE
      n_in = 0
      do ns = 1,NXR
         hxp(ns) = hx(ns)/npx(ns)
         do  n = 1, NPX(ns)
            n_in = n_in + 1
            hxy(n_in,1) = hxp(ns)
         end do
      end do

! Y - Direction

      n_in = 0
      do ns = 1,NYR
         hyp(ns) = hy(ns)/npy(ns)
         do  n = 1, NPY(ns)
            n_in = n_in + 1
            hxy(n_in,2) = hyp(ns)
         end do
      end do
      END IF

! setyting up the 3rd direction in radial plane for HEX-Z geometry 
      IF ( gmt_crd_type(1:3).EQ."HEXZ" ) THEN
        DO n = 1, NV
           hxy(n,3) = constant_triangle_mesh*HR_HEX/NP_SIDE
        END DO
      END IF                    

      IF( gmt_crd_type(1:3).EQ."XYZ" ) THEN
         k = 0
         DO nly = 1, NY
            DO nlx = 1, NX
               IF(n_fine(nly, nlx).NE.0) THEN
                  k = k + 1
                  sxy(k) = hxy(nly,2)*hxy(nlx,1)
               END IF
            END DO
         END DO
         
      ELSE IF ( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
         do k = 1, NH
            sxy(k) = hex_area/NCHM
         end do        
      END IF

      do n1 = 1, NZ
        nn = (n1-1)*NH
        do k = 1, NH
           kt = k + nn
           volume(kt) = hzt(n1)*sxy(k)
        end do
      end do

      IF( gmt_crd_type(1:3).EQ."XYZ" ) THEN
         do n1 = 1, NZ
         nn = (n1 - 1)*NH
         DO nly = 1, NY
            DO nlx = 1, NX
               k = n_fine(nly, nlx) 
               IF(k.NE.0) THEN
                  kt = k + nn
                  h_xyz(1, kt) = hxy(nlx, 1)
                  h_xyz(2, kt) = hxy(nly, 2)
                  s_xyz(1, kt) = hxy(nly,2)*hzt(n1)
                  s_xyz(2, kt) = hxy(nlx,1)*hzt(n1)
               END IF
           END DO ! NX
         END DO ! NY
         END DO !NZ
      ELSE IF ( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
      DO n1 = 1, NZ    
        nn = (n1 - 1)*NH
        DO k = 1, NH
           kt = k + nn
           h_xyz(1, kt) = constant_triangle_mesh*HR_HEX/NP_SIDE
           h_xyz(2, kt) = constant_triangle_mesh*HR_HEX/NP_SIDE
           h_xyz(3, kt) = constant_triangle_mesh*HR_HEX/NP_SIDE
           s_xyz(1, kt) = hex_side*hzt(n1)/NP_SIDE
           s_xyz(2, kt) = hex_side*hzt(n1)/NP_SIDE
           s_xyz(3, kt) = hex_side*hzt(n1)/NP_SIDE
        END DO
      END DO
      END IF

         do n1 = 1, NZ
            nn = (n1 - 1)*NH
            do k = 1, NH
               kt = k + nn
               s_xyz(NDIR, kt) = sxy(k)
               h_xyz(NDIR, kt) = hzt(n1)
            end do
         end do

      return
      end




      subroutine GMT_Set_Red_Black_Numbering
!=====================================================================*
!  Red-Black Numbering of Nodes                                       *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

      implicit none
      include 'sketch.fh'

! Input: NX, NY, n_fine(NY, NX, 1)
! Output: 1 red points, 2 black points
!  NRB(2) Number of tbe Red and Black Points in X-Y plane
!  NPRB(NRB(2), 2) - Node Number for the Red-Black Numbering 

! Local Variables
      integer nly, nlx, n_f
      
      nrb(1) = 0
      nrb(2) = 0

      do nly = 1, NY
         do nlx = 1, NX
            n_f =  n_fine(nly, nlx)
            if(n_f.ne.0) then
                if( ((nlx+nly)/2)*2.eq.(nlx+nly)) then
! Red Points
                  nrb(1) = nrb(1) + 1
                  nprb(nrb(1),1) = n_f
                else
! Black Points
                  nrb(2) = nrb(2) + 1
                  nprb(nrb(2),2) = n_f
                end if  
             end if ! n_f /= 0
          end do ! NX
      end do ! NY


      return
      end

      subroutine GMT_Set_Numbering_3D
!=====================================================================*
!  Numbering the Nodes in 3D                                          *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Input: NX, NY, N_FINE(NLY,NLX)
! Output: 
!            Index_NEIB(NDIR,1) - Number of the Left (-) Interface in Neib
!            Index_NEIB(NDIR,2) - Number of the Right (+) Interface in Neib

! Local Variables
      integer nd !, k, nk, i, j, i_flag, nly, nlx, ntr


! ir_la number of right interface for nd direction in la() massive 
      do nd = 1, NDIR ! X-Y-Z Direction 
        Index_Neib(nd,1) = nd  
        Index_Neib(nd,2) = nd + NDIR
      end do

      DO nd = 1, NDIR-1
         IF ( gmt_crd_type(1:4).EQ."HEXZ" .AND. (NCHM.NE.1) ) THEN
            N_TOT_DIR(nd) = N_TOT !2*N_TOT/3
         ELSE 
            N_TOT_DIR(nd) = N_TOT
         END IF
      END DO           
         
         N_TOT_DIR(NDIR) = N_TOT




      return 
      end

      subroutine GMT_Set_Assembly_Volume
!=====================================================================*
!  computing the volume of the NODES  for the Couling under PVM       *                          
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none      
      include 'sketch.fh'
! Input: nk_sketch(N_POLY), nz_sketch(N_POLY), NH, NZ, np_out(NH),
!        sxy(NH), hzt(NZ)
! Output: vol_ass(N_POLY, NZR) - volume of the nodes 
! Local Variables
      integer k, ns, np

      do ns = 1, NZR
         do np = 1, N_POLY
             vol_ass(np,ns) = 0.
         end do
      end do

      v_reactor = 0.
      do ns = 1, NZR
         do k = 1, NH
            np = np_out(k)
            vol_ass(np,ns) = vol_ass(np,ns) + sxy(k)*hz(ns)
            v_reactor = v_reactor + vol_ass(np,ns)
         end do
      end do


      return
      end


      subroutine GMT_Set_OUT_Core_Numbering
!=====================================================================*
!           Numbering of the Reactor Core  (without Reflector)        *                            
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE PRECISION, ONLY : dp 
      USE termsort_lib, ONLY: ops_lib, get_lib_nx_ny, get_lib_titles

      implicit none
      include 'sketch.fh'
! Input: npz(NZR), NZ, ns_out(NZR), NH, np_out(NH), sf(NNODE, NG), 
!        l(N_POLY, NZR)
! Output: NZ_Core_BEG - 1st Axial Layer of the Reactor Core
!         NZ_Core_END - Last Axial Layer of ther Reactor Core
!         NH_Core - Number of the Assemblies in the Core
!         N_Core_Tot - Total Numer of the Nodes in the Core
!         np_core(NH_Core) - Reactor Assemly Numbers of the Core
!         v_core - Volume of the Ractor Core
!         Numb_Poly_Core(N_POLY) - numbering of the Assembly in the Core
! Reactor Core
!      integer Nxr_B_Core(NYR), Nxr_E_Core(NYR), Nyr_B_Core,&
!              Nyr_E_Core, Nxr_Max_Core, Nxr_B_Min_Core
! Reactor
!      integer Nxr_B_Reactor(NYR), Nxr_E_Reactor(NYR), Nyr_B_Reactor, &
!              Nyr_E_Reactor, Nxr_Max_Reactor, Nxr_B_Min_Reactor

! Local Variables
      integer k,np,n1,ns,n_type, kt, nn, np_c, nc
      integer nly, nlx !, nlz
      LOGICAL first_x_line
!      real s_col

      integer k_reactor_core, flag_core

      INTEGER, PARAMETER :: NNX=20, NNY =200
      INTEGER            :: N_XS_X, N_XS_Y
      REAL(dp)  ::  X(NNX), Y(NNY)
      real nu_sigma_fission
      CHARACTER*18 :: title_x(NNX), title_y(NNY)
      CHARACTER*1000 :: Message


      INCLUDE 'YMNC_DAT.fh'
      REAL*4 ARGUM ( YMNCNARG  )
      REAL CONST( YMNCNFNUM ), &
            DCONST( YMNCNFNUM * YMNCNARG )
      INTEGER KSDCA(YMNCNFNUM)
      REAL, PARAMETER :: Convert_C_to_Kelvin=273.15   
      INTEGER :: i_xs_x 

      X(1:NNX) = 0.
      write(*,*) 'XS_MODEL=', TRIM(XS_MODEL)

! Setting the axial boundary of the reactor core
      do ns = 1, NZR 
         do k = 1, NH
           np  = np_out(k)
           n_type = l(np,ns)

         IF( XS_MODEL.EQ."SUBSET".OR.XS_MODEL.EQ."SUBSED") THEN

!          np = np_out(k)
!          ns = ns_out(n1)
!          nl = l(np,ns)

          CALL get_lib_nx_ny(n_type, N_XS_X, N_XS_Y)

          CALL get_lib_titles(n_type, title_x, title_y)

         DO i_xs_x = 1, N_XS_X 
          IF( INDEX( TRIM(title_x(i_xs_x)), "burnup" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = 0.0 ! burnup
            
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "\rho" ) &
                  .GT. 0  ) THEN 
!\rho(g/cm^3)            T_c(K)            T_f(K)          c_b(ppm)
            X(i_xs_x) = fdback0(3) ! coolant density
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_c" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback0(2) + Convert_C_to_Kelvin !   T_c
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_f" ) &
                  .GT. 0  ) THEN 
          X(i_xs_x) = fdback0(4) !+ Convert_C_to_Kelvin !   T_f
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "c_b" ) &
                  .GT. 0  ) THEN 
!            IF( nl == 111 ) THEN
!             write(*,*) 'nl = 111'
!             write(*,*) 'c_b =', fdback(k,n1,1)
!             pause
!            END IF 
          X(i_xs_x) = fdback0(1) ! C_Bor
          ELSE
            WRITE(Message, '(A,I3,A)') &
            " Error in the neutron XS library, "//&
            "Fuel Type =", n_type, "FA state variable "//&
              TRIM(title_x(i_xs_x))//" is unknown"//&
            " possible choice "//&
            "burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
!           write(*,*) TRIM(title_x(i_xs_x))
           write(*,*) " burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
            CALL MSC_ERR_Add_Message( 'ERROR',TRIM(Message) )
         END IF
         END DO ! i_xs_x = 1, N_XS_X 

!             write(*,*) 'n_type=', n_type
!           write(*,*) 'x =', x(1:N_xs_x)
!           pause

             CALL ops_lib(N_XS_X, N_XS_Y, x, n_type, y)  

           nu_sigma_fission = Y(10) ! k_inf

         ELSE IF( XS_MODEL.EQ."LINTAB") THEN 

           ARGUM(1)=fdback0(3)*1000.
            ARGUM(2)=fdback0(4)
         ARGUM(3)=fdback0(2)+Convert_C_to_Kelvin
           ARGUM(4)=0.
           ARGUM(5)=0. ! burnup
           ARGUM(6)=fdback0(1) 
           IF( n_type .GT. 0 ) THEN               
              CALL YMNCCALC  (ARGUM, n_type, CONST, DCONST, KSDCA)
           ELSE
              CALL YNREFsection (  ARGUM, n_type, CONST)
           END IF 

           nu_sigma_fission = const(12) 

         ELSE
!               write(*,*) 'XS_MODEL = POLYNOM'
!               write(*,*) 'ns =', ns, 'np =', np, 'n_type =', n_type
!               pause
             nu_sigma_fission = sf(n_type,NG)
!               IF(nzr == 3) THEN
!                 write(*,*) 'np=', np, 'nu_sigma_fission =',&
!                     nu_sigma_fission
!               END IF    
!               if(nu_sigma_fission.ne.0) then
!                  write(*,*) 'ns=', ns, 'np=', np
!               end if 

         END IF    

           if(nu_sigma_fission.ne.0) then
               k_reactor_core = np
               go to 10
! exit from the cycle
           end if
        end do
      end do
 
 10   continue

!      write(*,*) 'k_reactor_core=', k_reactor_core

      flag_core = 0
      do n1 = 1, NZ
        ns = ns_out(n1) 
        n_type = l(k_reactor_core,ns)

       IF( XS_MODEL.EQ."SUBSET".OR.XS_MODEL.EQ."SUBSED") THEN

          CALL get_lib_nx_ny(n_type, N_XS_X, N_XS_Y)
          CALL get_lib_titles(n_type, title_x, title_y)

         DO i_xs_x = 1, N_XS_X 
          IF( INDEX( TRIM(title_x(i_xs_x)), "burnup" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = 0.0 ! burnup
            
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "\rho" ) &
                  .GT. 0  ) THEN 
!\rho(g/cm^3)            T_c(K)            T_f(K)          c_b(ppm)
            X(i_xs_x) = fdback0(3) ! coolant density
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_c" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback0(2) + Convert_C_to_Kelvin !   T_c
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_f" ) &
                  .GT. 0  ) THEN 
          X(i_xs_x) = fdback0(4) !+ Convert_C_to_Kelvin !   T_f
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "c_b" ) &
                  .GT. 0  ) THEN 
!            IF( nl == 111 ) THEN
!             write(*,*) 'nl = 111'
!             write(*,*) 'c_b =', fdback(k,n1,1)
!             pause
!            END IF 
          X(i_xs_x) = fdback0(1) ! C_Bor
          ELSE
            WRITE(Message, '(A,I3,A)') &
            " Error in the neutron XS library, "//&
            "Fuel Type =", n_type, "FA state variable "//&
              TRIM(title_x(i_xs_x))//" is unknown"//&
            " possible choice "//&
            "burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
           write(*,*) TRIM(title_x(i_xs_x))
           write(*,*) " burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
            CALL MSC_ERR_Add_Message( 'ERROR',TRIM(Message) )
         END IF
         END DO ! i_xs_x = 1, N_XS_X 

!             write(*,*) 'n_type=', n_type
!           write(*,*) 'x =', x(1:N_xs_x)

             CALL ops_lib(N_XS_X, N_XS_Y, x, n_type, y)  
             nu_sigma_fission = Y(10) ! k_inf


         ELSE IF( XS_MODEL.EQ."LINTAB") THEN 

           ARGUM(1)=fdback0(3)*1000.
            ARGUM(2)=fdback0(4)
         ARGUM(3)=fdback0(2)+Convert_C_to_Kelvin
           ARGUM(4)=0.
           ARGUM(5)=0. ! burnup
           ARGUM(6)=fdback0(1)              
           IF( n_type .GT. 0 ) THEN               
              CALL YMNCCALC  (ARGUM, n_type, CONST, DCONST, KSDCA)
           ELSE
              CALL YNREFsection (  ARGUM, n_type, CONST)
           END IF 
           nu_sigma_fission = const(12) 
         ELSE
             nu_sigma_fission = sf(n_type,NG)
         END IF    

        if(nu_sigma_fission.ne.0 .AND. Flag_Core.eq. 0 ) then 
            NZ_Core_BEG = n1
            flag_core = 1
        end if
        if(nu_sigma_fission.eq.0 .AND. Flag_core.eq.1 ) then 
           NZ_Core_END = n1 - 1
           flag_core = 0
        end if
      end do


      if(NZ_Core_END.eq.0) NZ_Core_END = NZ

! set up the axial core height and axial reflectors heights
      NZR_Core_Beg = ns_out(NZ_Core_Beg)
      NZR_Core_End = ns_out(NZ_Core_End)


      hz_axial_reflector(1:2)=0.
      hz_core = 0.

      DO ns = 1, NZR_Core_Beg - 1
         hz_axial_reflector(1) = hz_axial_reflector(1) + hz(ns)
      END DO
      
      DO ns = NZR_Core_Beg, NZR_Core_End
         hz_core = hz_core + hz(ns)
      END DO

      DO ns = NZR_Core_End+1, NZR
         hz_axial_reflector(2) = hz_axial_reflector(2) + hz(ns)
      END DO

!        write(*,*) 'hz_core =', hz_core
!      PAUSE
! END Setting the axial boundary of the reactor core

! Reactor core position in XY plane

! 1st Axial Layer of the reactor core
      n1 = NZ_Core_Beg
      ns = ns_out(n1)


      IF(.NOT.flag_core_boundary) THEN
        DO np = 1, N_POLY
           n_core(np) = 0
        END DO
! the reactor core is defined by the fission materials 
        do nly = 1, NYR
         first_x_line =  .True.
         do nlx = 1, NXR
            np = npoly(nly,nlx)
            if(np.ne.0)  then
                n_type = l(np,ns)

           IF( XS_MODEL.EQ."SUBSET".OR.XS_MODEL.EQ."SUBSED") THEN
          CALL get_lib_nx_ny(n_type, N_XS_X, N_XS_Y)
          CALL get_lib_titles(n_type, title_x, title_y)

         DO i_xs_x = 1, N_XS_X 
          IF( INDEX( TRIM(title_x(i_xs_x)), "burnup" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = 0.0 ! burnup
            
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "\rho" ) &
                  .GT. 0  ) THEN 
!\rho(g/cm^3)            T_c(K)            T_f(K)          c_b(ppm)
            X(i_xs_x) = fdback0(3) ! coolant density
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_c" ) &
                  .GT. 0  ) THEN 
            X(i_xs_x) = fdback0(2) + Convert_C_to_Kelvin !   T_c
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "T_f" ) &
                  .GT. 0  ) THEN 
          X(i_xs_x) = fdback0(4) !+ Convert_C_to_Kelvin !   T_f
          ELSE IF( INDEX( TRIM(title_x(i_xs_x)), "c_b" ) &
                  .GT. 0  ) THEN 
!            IF( nl == 111 ) THEN
!             write(*,*) 'nl = 111'
!             write(*,*) 'c_b =', fdback(k,n1,1)
!             pause
!            END IF 
          X(i_xs_x) = fdback0(1) ! C_Bor
          ELSE
            WRITE(Message, '(A,I3,A)') &
            " Error in the neutron XS library, "//&
            "Fuel Type =", n_type, "FA state variable "//&
              TRIM(title_x(i_xs_x))//" is unknown"//&
            " possible choice "//&
            "burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
           write(*,*) TRIM(title_x(i_xs_x))
           write(*,*) " burnup,\rho(g/cm^3),T_c(K),T_f(K),c_b(ppm)"
            CALL MSC_ERR_Add_Message( 'ERROR',TRIM(Message) )
         END IF
         END DO ! i_xs_x = 1, N_XS_X 

!             write(*,*) 'n_type=', n_type
!           write(*,*) 'x =', x(1:N_xs_x)

             CALL ops_lib(N_XS_X, N_XS_Y, x, n_type, y)  
             nu_sigma_fission = Y(10) ! k_inf

           ELSE IF( XS_MODEL.EQ."LINTAB") THEN 

           ARGUM(1)=fdback0(3)*1000.
            ARGUM(2)=fdback0(4)
         ARGUM(3)=fdback0(2)+Convert_C_to_Kelvin
           ARGUM(4)=0.
           ARGUM(5)=0. ! burnup
           ARGUM(6)=fdback0(1)              
           IF( n_type .GT. 0 ) THEN               
              CALL YMNCCALC  (ARGUM, n_type, CONST, DCONST, KSDCA)
           ELSE
              CALL YNREFsection (  ARGUM, n_type, CONST)
           END IF 
           nu_sigma_fission = const(12) 
           ELSE
             nu_sigma_fission = sf(n_type,NG)
           END IF    


                if(nu_sigma_fission.ne.0) then
                      IF(first_x_line) THEN
                     first_x_line = .false.
                   END IF
                   n_core(np) = 1
                   Index_Core(nly,nlx) = np
                end if
            end if
         end do
        end do
      ELSE
!  the reactor core is defined by the code user
! New Vatria
        do nly = 1, NYR
         do nlx = 1, NXR
            np = npoly(nly, nlx)
            IF(np.NE.0) THEN
              IF(n_core(np).NE.0) Index_Core(nly, nlx)=npoly(nly,nlx)
            END IF
         end do
        end do
      END IF


      nc = 0
      do nly = 1, NYR
         do nlx = 1, NXR
            np = npoly(nly, nlx)
            np_c = Index_Core(nly, nlx)
            IF(np_c.NE.0) THEN
              nc = nc + 1
              Numb_Reactor_Core(nc) = np
              Numb_Poly_Core(np) = nc
            END IF
         end do
      end do

      NH_Core = 0
      do k = 1, NH
       np  = np_out(k)
       nlx = N_Coord(np,1)
       nly = N_Coord(np,2)
       if(Index_Core(nly, nlx).ne.0) THEN
         NH_Core = NH_Core + 1
         np_core(NH_Core) = k
       end if
      end do

! Attention vcore include not only fission materials but also 
! holes in the core
      v_core = 0.
      do n1 = NZ_Core_Beg, NZ_Core_End
      nn = (n1 - 1)*NH
         do k = 1, NH_Core
            kt = np_core(k) + nn
            v_core = v_core + volume(kt)
         end do
      end do

      N_Core_Tot = NH_Core*(NZ_Core_END - NZ_Core_BEG + 1)
 
      call GMT_Reactor_Numbering(NXR, NYR, Index_Core, &
          Nxr_B_Core, Nxr_E_Core, Nyr_B_Core, Nyr_E_Core, &
              Nxr_Max_Core, Nxr_B_Min_Core)

      call GMT_Reactor_Numbering(NXR, NYR, npoly, Nxr_B_Reactor, &
          Nxr_E_Reactor, Nyr_B_Reactor, Nyr_E_Reactor, &
              Nxr_Max_Reactor, Nxr_B_Min_Reactor)

      return
      end

      Subroutine GMT_Reactor_Numbering(NXR, NYR, npoly, &
          Nxr_B_Core, Nxr_E_Core, Nyr_B_Core, Nyr_E_Core, &
              Nxr_Max_Core, Nxr_B_Min_Core)
!=====================================================================*
! Compute Numbering of the Reactor Core for Output Procedure          *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
! Input: 
      integer NXR, NYR
      integer  npoly(NYR,NXR)
! Output:
      integer Nxr_B_Core(NYR), Nxr_E_Core(NYR), Nyr_B_Core, Nyr_E_Core, &
              Nxr_Max_Core, Nxr_B_Min_Core
! Local Variables:
      integer Flag_X_B, Flag_Y_B
      integer nlx, nly, np

! Computing the Reactor Position (for Output)
       
      Nyr_B_Core = 0
      Flag_Y_B = 0

      do nly = 1, NYR
          Nxr_B_Core(nly) = 0
          Flag_X_B = 0
          do nlx = 1, NXR

          np = npoly(nly,nlx)

          if(np.ne.0) then
               if(Flag_X_B.eq.0) then
                 Nxr_B_Core(nly) = nlx
                 Flag_X_B = 1
               end if
               Nxr_E_Core(nly) = nlx
          end if ! Assembly /= 0
          end do ! NXR

          if(Flag_X_B.eq.1) then
            if(Flag_Y_B.eq.0) then
              Nyr_B_Core = nly
              Flag_Y_B = 1
            end if
            Nyr_E_Core = nly
          end if
      end do ! NLY


! Compting the line with minimum position in the left
       Nxr_B_Min_Core = 10000
       Nxr_Max_Core = 0
       do nly = Nyr_B_Core, Nyr_E_Core
         Nxr_B_Min_Core = min(Nxr_B_Min_Core, Nxr_B_Core(nly))
         Nxr_Max_Core = max(Nxr_Max_Core, Nxr_E_Core(nly))
       end do

      return
      end


      subroutine GMT_Set_3D_Neighbours
!=====================================================================*
!  Numbering of the Boundary Nodes & Interface in 3D                  *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Triangle 
      implicit none
      include 'sketch.fh'

! Input: NH, NZ, LA(NE+1,NH), N_FINE(NY,NX)
! Output: Neib(NE_T, N_TOT) - Neibours of the node k in 3D
!         k_bound(N_BOUND_M) - boundary nodes
!         nd_bound(N_BOUND_M) - direction of the boundary interface
!         nl_bound(N_BOUND_M) - left or right boundary interface
!         Numb(N_TOT, NDIR) - Node Numbering for X-Y-Z directions

! Local Variables
      
      integer n1, nn, nn_left, nn_right, kt, nd, nly, nlx, kt_neib,&
              k, np, nt, nom


      IF ( gmt_crd_type(1:4).EQ."HEXZ" .AND. NCHM.NE.1) THEN
        CALL Set_Numb_Triangle(N_TOT, NCHM, NDIR, NH, NXR, NYR, NZ,&
       npoly, Numb, N_TOT_DIR)
      ELSE 


! Numbering of the node in X-Y-Z directions
      do n1 = 1, NZ
       nn = (n1-1)*NH
! Numbering of the nodes in X direction
       np = 0
       DO nly = 1, NY
          DO nlx = 1, NX
             k = n_fine(nly, nlx)
             IF( k.NE.0) THEN
               np = np+1
                  nt = np + nn
               kt = k  + nn
               Numb(nt, 1) = kt
             END IF
          END DO
       END DO

! Numbering of the nodes in Y direction
       np = 0
       DO nlx = 1, NX
             DO nly = 1, NY
             k = n_fine(nly, nlx)
             IF( k.NE.0) THEN
               np = np+1
               nt = np + nn
               kt = k  + nn
               Numb(nt, 2) = kt
             END IF
          END DO
       END DO    
       
       IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
          np = 0
          DO k = 1, NH
             nom = k
! we start the line if there is no left neighbour
            IF(la(nom,3).eq.I_BOUND_NODE) THEN
              np = np + 1
              nt = np + nn
!              write(*,*) 'k = ', k, 'np =', np
              Numb(nt, 3) = nom + nn
! do while there is the right neighbour
              DO while( la(nom, 3 + NE/2) .NE. I_BOUND_NODE ) 
                 np = np + 1
                 nt = np + nn
                 Numb(nt, 3) =  la(nom, 3 + NE/2) + nn
!                  write(*,*) 'np =', np, 'la(nom, 3 + NE/2)', &
!                                            la(nom, 3 + NE/2)
                 nom = la(nom, 3 + NE/2)
              END DO ! while
            END IF ! la(k,3) 
          END DO 
       END IF ! gmt_crd_type == hexz 
               
      END DO ! NZ

      END IF ! IF ( gmt_crd_type(1:4).EQ."HEXZ" .AND. NCHM.EQ.1) THEN

!      write(*,*) '(Numb(nt,3) ='
!      write(*, '(10I3)') (Numb(nt,3), nt = 1, N_TOT)
!      pause

! Numbering in Z direction 
      nn = 0
      do k = 1, NH
         do n1 = 1, NZ
            nn = nn + 1
            kt = k + (n1 - 1)*NH
            Numb(nn, NDIR) = kt
         end do
      end do
   

! Neighbors in 3D 
      do n1 = 1, NZ
       nn = (n1 - 1)*NH
       nn_left = (n1 - 2)*NH
       nn_right = n1*NH
       do k = 1, NH
          kt = k + nn
          do nd = 1, NDIR-1 !
! Left
            kt_neib = la(k,nd)
            if(kt_neib.ne.I_BOUND_NODE) THEN
               kt_neib = kt_neib + nn
            else
               kt_neib = I_BOUND_NODE
            end if
            Neib(nd, kt) = kt_neib
! Right
            kt_neib = la(k, nd + NE/2)
            if(kt_neib.ne.I_BOUND_NODE) THEN
              kt_neib = kt_neib + nn
            else
             kt_neib = I_BOUND_NODE
            end if
            Neib(nd + NDIR, kt) = kt_neib

         end do
! Neighbours in axial directions
        if(n1.eq.1) then
           kt_neib = I_BOUND_NODE
        else
           kt_neib = k + nn_left
        end if
        Neib(NDIR,kt) = kt_neib

        if(n1.eq.NZ) then
           kt_neib = I_BOUND_NODE
        else
           kt_neib = k + nn_right
        end if
        Neib(2*NDIR,kt) = kt_neib

       end do
      end do


      return
      end


      subroutine GMT_Set_Index_Left_Interface(la_poly)
!=====================================================================*
! Set index of the left interface number of the node in X-Y-Z direct  *
!                   Vyachreslav Zimin (c) 18 August 1999              *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE Index_Map_Faces_Class, ONLY : Set_Index_Map_Faces
      USE Core_Map_Class, ONLY : set_faces_map

      implicit none
      include 'sketch.fh'


      INTEGER la_poly(N_POLY, NE+1)
      integer nl, nd, nd_right, kt, k, k_right
      INTEGER ns, np, nom, nlx, nly , n1
      LOGICAL flag_boundary_face


      nl = 0

      do nd = 1, NDIR

       nd_right = Index_Neib(nd,2)

       do kt = 1, N_TOT_DIR(nd)

          k = Numb(kt, nd) 

          nl = nl + 1 ! Left Interface

          Index_Left_Interface(k, nd) = nl

          k_right = Neib(nd_right, k)

! Right Boundary Node
          if(k_right.eq.I_BOUND_NODE) nl = nl + 1 ! Right Interface

        end do
      
      end do

! Setting Array of the Left Faces for Polycells
      ns =  0
! X-direction 
      nd = 1
      DO n1 = 1, NZR
         DO nly = 1, NYR
            DO nlx = 1, NXR
               np = npoly(nly, nlx)
               IF ( np /= 0 ) THEN
                 ns = ns + 1
                 Index_Left_Interface_Poly(np, n1, nd) = ns
                 flag_boundary_face = .FALSE.
                     IF ( nlx == NXR ) THEN
                  flag_boundary_face = .TRUE.
                     ELSE IF ( npoly(nly, nlx+1) == 0 ) THEN
                    flag_boundary_face = .TRUE.
                     END IF ! nlx == NXR       
                 IF(  flag_boundary_face ) ns=ns+1
               END IF
            END DO
          END DO
      END DO                          

! Y Direction
      nd = nd + 1   
      DO n1 = 1, NZR
        DO nlx = 1, NXR
            DO nly = 1, NYR
               np = npoly(nly, nlx)
               IF ( np /= 0 ) THEN
                 ns = ns + 1
!                 write(*,*) 'nd=2', 'np=', np, 'ns =', ns
!                 pause 
                 Index_Left_Interface_Poly(np, n1, nd) = ns
                 flag_boundary_face = .False.
                     IF ( nly == NYR ) THEN
                  flag_boundary_face = .TRUE.
                     ELSE IF ( npoly(nly+1, nlx) == 0 ) THEN
                    flag_boundary_face = .TRUE.
                     END IF ! nly == NYR  
                 IF(  flag_boundary_face ) ns=ns+1

               END IF ! np /= 0 
            END DO
          END DO
      END DO                          
! v-direction for Hexagonal Mesh
      IF (gmt_crd_type(1:4).EQ."HEXZ") THEN
         nd = nd + 1
         DO n1 = 1, NZR
          DO k = 1, N_POLY
             nom = k
! we start the line if there is no left neighbour
            IF(la_poly(nom,3).eq.I_BOUND_NODE) THEN
              ns = ns + 1
              Index_Left_Interface_Poly(nom, n1, nd) = ns
!              write(*,*) 'nom, ns =', nom, ns
!              pause
! do while there is the right neighbour
              DO while( la_poly(nom, 3 + NE/2) .NE. I_BOUND_NODE ) 
                 ns  = ns + 1
                nom = la_poly(nom, 3 + NE/2)
                Index_Left_Interface_Poly(nom, n1, nd) = ns
!!              write(*,*) 'nom, ns =', nom, ns
              END DO ! while
              ns = ns +1
            END IF ! la_poly(k,3) 
          END DO ! N_POLY
         END DO ! NZR
      END IF ! gmt_crd_type == hexz 

! Z-Direction
      nd = nd + 1   
            DO nly = 1, NYR
             DO nlx = 1, NXR
               np = npoly(nly, nlx)
               IF ( np /= 0 ) THEN
                   DO n1 = 1, NZR
                    ns = ns + 1
                    Index_Left_Interface_Poly(np, n1, nd) = ns
                 END DO 
                 ns=ns+1
                END IF 
            END DO
        END DO

        CALL Set_Index_Map_Faces(la_poly)
        CALL set_faces_map

! Setting Arrays of the Faces for Printing the Core Map
! X-direction 

      return
      end


      subroutine GMT_DBG_Output
!=====================================================================*
!               Output Gemetry Data                                   *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Boundary
      implicit none
      include 'sketch.fh'

           character*80 Header_Map
      character*4 val_fmt
      character*4 val_char(0:N_POLY)

! Input: All Geometry Data
! Output: Files: 
!             Output/Lines.dat, Output/Points.dat, Output/Red_Black.dat
!      Local Variables
      integer nly, nlx, i, j, k, nd, ind
!      character*100 Message
      
      open(io_unit,file = 'Output_Debug/GEOMETRY.dat',status='UNKNOWN')
         write(io_unit,'("Index_Neib(nd,1): ", 3I4)')&
                            (Index_Neib(nd,1),nd=1,NDIR)

         write(io_unit,'("Index_Neib(nd,2): ", 3I4)') &
                            (Index_Neib(nd,2),nd=1,NDIR)

         write(io_unit,'(A)') "(n_b(nly),nly=1,NY):"
         write(io_unit,1) (n_b(nly),nly=1,NYR)

         write(io_unit,'(A)') "(n_e(nly),nly=1,NYR):"
         write(io_unit,1) (n_e(nly),nly=1,NYR)

         write(io_unit,'(A)') '2D Neighboring Nodes:'
         do k = 1, NH
            write(io_unit,2) (la(k,i),i=1,NE+1)
         end do

         write(io_unit,*) '3D Neighbors Neib'
!         write(*,*) 'NE_T =', NE_T
!         pause
         do k = 1, N_TOT
            write(io_unit,'("k=",I8, " Neib=", 10I8)') &
             k, (Neib(i,k),i=1,NE_T)
         end do

         write(io_unit,*) '3D Neighbors Neib_REP_STR'
!         write(*,*) 'NE_T =', NE_T
!         pause
         do k = 1, N_TOT
            write(io_unit,'("k=",I8, " Neib_REP_STR=", 10I8)') &
             k, (Neib_REP_STR(i,k),i=1,NE_T)
         end do

         write(io_unit,*) '3D Numbering'
         do nd = 1, NDIR
            write(io_unit,*) 'nd =', nd
            write(io_unit,111) (Numb(k,nd), k = 1, N_TOT_DIR(nd))
         end do
!      close(io_unit)

!      open(io_unit,file = 'Output_Debug/Points.dat',status='UNKNOWN') 
 
          write(io_unit,'(A)') '(np_out(i),i=1,NH)'
          write(io_unit,1) (np_out(i),i=1,NH)

          write(io_unit,'(A)') '(ns_out(i),i=1,NZ)'
          write(io_unit,1) (ns_out(i),i=1,NZ)

          write(io_unit,'(A)') '((hxy(i,j),j=1,2),i=1,NXYM)'
          write(io_unit,3) ((hxy(i,j),j=1,2),i=1,NXYM)


          write(io_unit,'(A)') '( hzt(i),i=1,NZ)'
          write(io_unit,3) ( hzt(i),i=1,NZ)

          write(io_unit,'(A)') '(hzp(i),i=1,NZR)'
          write(io_unit,3) (hzp(i),i=1,NZR)

          write(io_unit,'(A)') '(hxp(i),i=1,NXR)'
          write(io_unit,3) (hxp(i),i=1,NXR)

          write(io_unit,'(A)') '(hyp(i),i=1,NYR)'
          write(io_unit,3) (hyp(i),i=1,NYR)

          DO nd = 1, NDIR
          write(io_unit,'(A,I2,A)') 'nd =',nd,'(h_xyz(nd, i),i=n_tot)'
          write(io_unit,3) (h_xyz(nd,i),i=1, N_TOT)
          write(io_unit,'(A,I2,A)') 'nd =',nd,'(s_xyz(nd, i),i=n_tot)'
          write(io_unit,3) (s_xyz(nd, i),i=1, N_TOT)
          END DO

          write(io_unit,'(A)') '(n_fine(nly,nlx),nlx=1,NX),nly=1,NY)'
          write(io_unit,1) ((n_fine(nly,nlx),nlx=1,NX),nly=1,NY)

          write(io_unit,'(A)') &
                           'Boundary Nodes (k_bound(k),k=1,N_BOUND)'
          write(io_unit,4) (k_bound(k),k=1,N_BOUND)
          write(io_unit,'(A)') '(nl_bound(k),k=1,N_BOUND)'
          write(io_unit,4) (nl_bound(k),k=1,N_BOUND)
          write(io_unit,'(A)') '(nd_bound(k),k=1,N_BOUND)'
          write(io_unit,4) (nd_bound(k),k=1,N_BOUND)

          write(io_unit,'(A)') 'Reactor Core Numbering'

          write(io_unit,1)  (Numb_Reactor_Core(k), &
                    k = 1, NP_Reactor_Core)

         Header_Map = "Core Numbers from N_POLY to NP_Reactor_Core"
         val_fmt = "A4"
         DO ind = 1, N_POLY
             write(*,*) 'written in val_char, ind, NUmb_poly_core', &
                        Numb_Poly_Core(ind)
             write(val_char(ind), '(I4)')  Numb_Poly_Core(ind)
         END DO
         write(*,*) 'written in val_char'
         call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, io_unit, val_char, val_fmt)

! N_COORD(:,1)
         Header_Map = "N_COORD(N_POLY, 1)"
         val_fmt = "A4"
         DO ind = 1, N_POLY
             write(val_char(ind), '(I4)')  N_Coord(ind,1)
         END DO
         call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, io_unit, val_char, val_fmt)

! N_COORD(:,2)
         Header_Map = "N_COORD(N_POLY, 2)"
         val_fmt = "A4"
         DO ind = 1, N_POLY
             write(val_char(ind), '(I4)')  N_Coord(ind,2)
         END DO
         call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, io_unit, val_char, val_fmt)

          write(io_unit,'(A)') 'RED-Blasck Numbering'
          write(io_unit,'(A)') '(nrb(i),i=1,2)'
          write(io_unit,1) (nrb(i),i=1,2)
          write(io_unit,'(A)') 'red points'
          write(io_unit,1) (nprb(i,1),i=1,nrb(1))
          write(io_unit,'(A)') 'black points'
          write(io_unit,1) (nprb(i,2),i=1,nrb(2))

          DO i = 1, NZR
          write(Header_Map, '("Core Material Properties, Axial Layer ", &
               I3)') i
          DO ind = 1, N_POLY
              write(val_char(ind), '(I4)')  l(ind,i)
          END DO
          call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
         NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
         NYR_E_Reactor, npoly, Header_Map,io_unit, val_char, val_fmt)
         END DO



        DO nd = 1, 1 ! NDD 
          DO i = 1, 1 ! NZR
          write(Header_Map, '("Index_Left_Interface_Poly (nz, nd)", &
               I3, I3)') i, nd
          DO ind = 1, N_POLY
              write(val_char(ind), '(I4)')  &
                       Index_Left_Interface_Poly(ind, i, nd)
          END DO
          call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
         NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
         NYR_E_Reactor, npoly, Header_Map,io_unit, val_char, val_fmt)
         END DO
        END DO

      close(io_unit)

    1 FORMAT(20I4)
    2 FORMAT(10I5)
    3 FORMAT(10F8.4)
    4 FORMAT(10I8)
  111 FORMAT(7I8)
      return
      end

      subroutine GMT_Output_Parameter(unit)
!=====================================================================*
!               Output Geometry Parameters                            *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer unit

!      call OUTput_Write_Separator(io_unit)
      write(unit, *)
      write(unit, '(A)') " Geometry Parameters"
      write(unit, '(A, I8)') &
        "    Number of Dimensions:                               ", &
        NDD
      write(unit, '(A, I2,"x", I2,"x",I2)') &
        "    Initial Spatial Mesh (x,y,z):                       ",&
        NXR, NYR, NZR
      write(unit, '(A, I2,"x", I2,"x",I2)') &
        "    Neuronics Spatial Mesh (x,y,z):                     ", &
        NX, NY, NZ
      write(unit, '(A, I8)') &
        "    Number of Bundles (Assemblies) in x-y plane:        ", &
        N_POLY
      write(unit, '(A, I8)') &
        "    Number of Nodes per Bundle (Assembly) in x-y plane: ", &
        NCHM
      write(unit, '(A, I8)') &
        "    Number of Neutronics Nodes in x-y plane:            ", &
        NH
      write(unit, '(A, I8)') &
        "    Number of Bundle (Assembly) Types:                  ", &
        N_BUNDLE_TYPE
!      call OUTput_Write_Separator(io_unit)

      return
      end
      
      subroutine GMT_Output_Data(unit)
!=====================================================================*
!               Output Geometry Parameters                            *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer unit
      integer i, j, ind, n, nd

           character*80 Header_Map
      character*4 val_fmt
      character*4 val_char(0:N_POLY)

!      call OUTput_Write_Separator(io_unit)
! " Geometry Data :"
      write(unit, '(A)') " Geometry Data :"

! "    Numbering the Reactor Bundles (Assemblies):"       
      Header_Map = "    Numbering the Reactor Bundles (Assemblies):"
      val_fmt = "A4"
      ind = 0
      DO i=1,NYR
        DO j=1,NXR
           IF(npoly(i,j) .NE. 0) THEN
             ind = ind + 1
             write(val_char(ind), '(I4)')  npoly(i,j)
           END IF
         END DO
      END DO
      call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, unit, val_char, val_fmt)
      
! "    Initial Spatial Mesh X, Y, Z"
      WRITE(unit, '(A)') &
       "    Initial Spatial Mesh X:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hx(ind), ind = 1, NXR)
      WRITE(unit, '(A)') &
       "    Initial Spatial Mesh Y:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hy(ind), ind = 1, NYR)
      WRITE(unit, '(A)') &
       "    Initial Spatial Mesh Z:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hz(ind), ind = 1, NZR)

!  "    Neutronics Spatial Mesh X, Y, Z"
      WRITE(unit, '(A)') &
       "    Neutronics Spatial Mesh X:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hxy(ind,1), ind = 1, NX)
      WRITE(unit, '(A)') &
       "    Neutronics Spatial Mesh Y:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hxy(ind,2), ind = 1, NY)
      WRITE(unit, '(A)') &
       "    Neutronics Spatial Mesh Z:"
      WRITE(unit, '(8x, 10F7.3)') &
       (hzt(ind), ind = 1, NZ)

! "    Reactor Loading by the Bundle Types        "
      Header_Map = "    Reactor Loading by the Bundle Types        "
      val_fmt = "A4"
      DO ind = 1, N_POLY
          write(val_char(ind), '(I4)')  Core_Load(ind)
      END DO
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, unit, val_char, val_fmt)

! "    Bundle Material Compositions:"      
      write(unit, '(A)')&
       "    Bundle Material Compositions:"
      DO ind = 1, N_Bundle_Type
         write(unit, ' (8x, I3,":", 10I4, /, 12x, 10I4,/, 12x, 10I4)')&
           ind, (Bundle_Compos(n, ind), n = 1, NZR)
      END DO

      write(unit,*)
      
!"    Boundary Conditions Types "
      write(unit, '(A)')&
       "    Boundary Conditions Types "//&
      "(left, right, up, down, top, bottom) :"
      write(unit, '(A)')&
       '      "0" is zero flux, "1" is logarithmic derivative'
      write(unit, '(8x, 6I2)') ((i_dr(i,nd),i=1,2),nd=1,3)

! "    Boundary Conditions Constants:"
      write(unit, '(A)')&
       "    Boundary Conditions Constants:"
      write(unit, '("       left  x:  ", 10F7.3)') (dr(n,1,1),n=1,NG)
      write(unit, '("       right x:  ", 10F7.3)') (dr(n,2,1),n=1,NG)
      write(unit, '("          up y:  ", 10F7.3)') (dr(n,1,2),n=1,NG)
      write(unit, '("        down y:  ", 10F7.3)') (dr(n,2,2),n=1,NG)
      write(unit, '("         top z:  ", 10F7.3)') (dr(n,1,3),n=1,NG)
      write(unit, '("      bottom z:  ", 10F7.3)') (dr(n,2,3),n=1,NG)

! "    Axial Buckling :"
      IF(NZ.eq.1) THEN
         WRITE(unit, '(A)')&
          "    Axial Buckling :"
          WRITE(unit, '(8x, 10F7.3)') &
         (b2(n), n = 1, NG)
      END IF

      call OUTput_Write_Separator(io_unit)

      return
      end


      subroutine GMT_Output_Core_Parameter(unit)
      implicit none
      include 'sketch.fh'
      integer unit

      write(unit, *)
!      call OUTput_Write_Separator(unit)

      write(unit, '(A)')&
       "    Reactor Core Parameters:"
      write(unit, '(A, I8)') &
       "       Number of Bundles (Assemblies) in x-y plane:        ", &
       NP_Reactor_Core
      write(unit, '(A, I8)') &
       "       Number of Axial Layers                     :        ", &
       NZR_Core

!      call OUTput_Write_Separator(unit)

      RETURN 
      END


      subroutine GMT_Output_Core_Data(unit)
      implicit none
      include 'sketch.fh'
      integer unit

      integer ind
      character*80 Header_Map
      character*4 val_fmt
      character*4 val_char(0:N_POLY)

!      call OUTput_Write_Separator(unit)
      write(unit, '(A, F8.3)') &
       "       Reactor Core Height, [cm]                  :        ", &
               hz_core

      Header_Map = "    Reactor Core Geometry:"
      val_fmt = "A4"
      DO ind = 1, N_POLY
          write(val_char(ind), '(I4)')  Numb_Poly_Core(ind)
      END DO

      call OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      call OUTput_Write_Separator(unit)

      RETURN 
      END

      SUBROUTINE GMT_Set_LA(NH, NE, NX, NY, n_fine, la, I_BOUND_NODE)
      IMPLICIT NONE
! Input:
      INTEGER NX, NY, NH, NE, I_BOUND_NODE
      INTEGER n_fine(NY, NX)
! Output: 
      INTEGER la(NH, NE+1)
! Locals:
      INTEGER nlx, nly, n_f, k, j        
      LOGICAL Flag_HEX 

         do k = 1, NH
          do j = 1, NE+1
           la(k,j) = 0.
          end do
         end do

        IF(NE.eq.4) THEN
          Flag_HEX = .False.
        ELSE IF(NE.eq.6) THEN
           Flag_HEX = .True.
        END IF 
                          
              

! hexagonal geometry  neighbour's number

          do nly = 1, NY
            do nlx = 1, NX
               n_f =  n_fine(nly, nlx)
               if(n_f.ne.0) then

                    la(n_f,NE+1) = n_f

!left
                    IF( nlx==1 ) THEN
                       la(n_f,1) = I_BOUND_NODE
                    ELSE IF ( n_fine(nly,nlx-1)==0 ) THEN
                       la(n_f,1) = I_BOUND_NODE
                    ELSE
                       la(n_f,1) = n_fine(nly,nlx-1)
                    END IF 
! up
                    IF(nly==1 ) THEN
                      la(n_f,2) = I_BOUND_NODE
                    ELSE IF ( n_fine(nly-1,nlx)==0) THEN
                      la(n_f,2) = I_BOUND_NODE
                    ELSE
                      la(n_f,2) = n_fine(nly-1,nlx)
                    END IF
! hexagonal         
                    IF( Flag_Hex ) THEN
                    IF( nlx==1 ) THEN
                      la(n_f,3) = I_BOUND_NODE
                    ELSE IF (nly==1) THEN
                      la(n_f,3) = I_BOUND_NODE
                    ELSE IF ( n_fine(nly-1, nlx-1)==0 ) THEN
                      la(n_f,3) = I_BOUND_NODE
                    ELSE 
                      la(n_f,3) = n_fine(nly-1,nlx-1)
                    END IF
                    END IF 
! right
                    IF(nlx==NX) THEN
                        la(n_f,NE/2+1) = I_BOUND_NODE
                    ELSE IF (n_fine(nly,nlx+1)==0) THEN
                        la(n_f,NE/2+1) = I_BOUND_NODE
                    ELSE                          
                        la(n_f,NE/2+1) = n_fine(nly,nlx+1)
                    END IF 
! down
                              
                    IF(nly==NY) THEN
                        la(n_f,NE/2+2) = I_BOUND_NODE
                    ELSE IF (n_fine(nly+1,nlx)==0) THEN
                        la(n_f,NE/2+2) = I_BOUND_NODE
                    ELSE
                       la(n_f,NE/2+2) = n_fine(nly+1,nlx)
                    END IF

! hexagonal         
                    IF( Flag_Hex ) THEN
                    IF( nlx==NX) THEN
                       la(n_f,NE/2+3) = I_BOUND_NODE
                    ELSE IF(nly==NY) THEN
                       la(n_f,NE/2+3) = I_BOUND_NODE
                    ELSE IF( n_fine(nly+1,nlx+1)==0) THEN
                       la(n_f,NE/2+3) = I_BOUND_NODE
                    ELSE
                       la(n_f,NE/2+3) = n_fine(nly+1,nlx+1) 
                    END IF 
                    END IF

                  end if ! n_f /= 0
            end do
           end do

      RETURN
      END 


      subroutine GMT_Set_3D_Neighbours_Repeated_Structure
!=====================================================================*
! Renumbering of the Boundary Nodes in the Case of the                *
!            REPEATED STRUCTURES                                      *  
!                   Vyachreslav Zimin (c) September 6 2005            *
!               e-mail: slava@ets.mephi.ru                            *
!=====================================================================*
!     USE GeoMeTry_Triangle 
      implicit none
      include 'sketch.fh'

! Input: NH, NZ, LA(NE+1,NH), 
!        Neib(NE_T, N_TOT) - Neibours of the node k in 3D
! Output: Neib_REP_STR(NE_T, N_TOT) - Neibours of the node k in 3D
!          for the repeated structure
! Local Variables
      
      integer n1, nn, nn_left, nn_right, kt, nd,  kt_neib,&
              k,  kt_neib_right, kt_neib_left

      LOGICAL flag_repeated_structure 

! first of all set up Neib_REP_STR
      Neib_REP_STR(:, :)=Neib(:, :)

      flag_repeated_structure = .False.

! checking that there are repeated structure
      DO nd=1, NDIR-1
       IF( i_dr(1,nd)==2 .OR. i_dr(2,nd) == 2 ) THEN
         flag_repeated_structure  = .True.
         EXIT
       END IF
      END DO !nd=1, NDIR-1

      IF ( flag_repeated_structure ) THEN   
! Renumbering of the BOUNDARY Neighbors in 3D for the repaeted structure
      do n1 = 1, NZ
       nn = (n1 - 1)*NH
       nn_left = (n1 - 2)*NH
       nn_right = n1*NH
       do k = 1, NH
          kt = k + nn
          do nd = 1, NDIR-1 !
! Left
          IF( i_dr(1,nd)==2 ) THEN  
            kt_neib = la(k,nd)
            if(kt_neib == I_BOUND_NODE) THEN
               kt_neib_right = la(k, nd + NE/2)+nn
               Neib_Rep_STR(nd, kt) = kt_neib_right
            end if
          END IF !( i_dr(1,nd)==2 ) 
! Right
          IF( i_dr(2,nd)==2 ) THEN  
            kt_neib = la(k, nd + NE/2)
            if(kt_neib == I_BOUND_NODE) THEN
              kt_neib_left = la(k,nd)+nn
              Neib_REP_STR(nd + NDIR, kt) = kt_neib_left
            end if
          END IF !( i_dr(2,nd)==2 ) 

         end do ! nd = 1, NDIR-1 !

       end do ! k = 1, NH
      end do ! n1 = 1, NZ

      END IF  ! ( flag_repeated_structure )

      return
      end

