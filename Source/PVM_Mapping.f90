      SUBROUTINE PVM_Input_Map_Matrices
!=====================================================================*
!        Input Mapping Matrices from the file FILE_MAP                *
!          for TRAC/SKETCH coupled calculations                       * 
!                   Vyachreslav Zimin (c) 18 October 1999             *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
! Input Mapping Arrays
      implicit none
      include 'sketch.fh'

! Local Variables
      integer L_HEADER, line
      parameter (L_HEADER = 5)

      integer ios
      character*132 Message

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      logical error_find

!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'

      if(TH_Model.EQ."External") then

        open(io_unit,file = FILE_MAP,status='old', iostat=ios)
        call Iostat_Error_Check&
      (ios,"Could not find the input FILE_MAP file "//FILE_MAP)


! reading the FILE Header
        do line = 1, L_HEADER
         read(io_unit,*) 
        end do

! reading PVM_EXT_CODE
!      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_EXT_CODE", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier PVM_EXT_CODE in the FILE_MAP,  '//&
         'exteranl T/H code name is not specified')
      else
        read(io_unit,fmt=*,iostat=ios) &
           TRAC_Version
        call Iostat_Error_Check&
     (ios,"Error in Reading External T/H Code Name from the FILE_MAP")

        IF( (TRAC_Version .NE. "TRAC-BF1").AND.&
           (TRAC_Version .NE. "TRAC-PF1") ) THEN
            WRITE(Message, '(A,A,A)') &
            "External T/H Code Name =", TRAC_Version, &
            " in FILE_MAP,  known types 'TRAC-BF1', "//&
            " and 'TRAC-PF1'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
         END IF
      end if
          
! PVM_CNV_UNIT: Conv_Cool_Dens, Conv_Cool_Temp, Conv_Fuel_Temp
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_CNV_UNIT", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier PVM_CNV_UNIT in the FILE_MAP,  '//&
         'constants of the unit conversion are not given')
      else
        read(io_unit,fmt=*,iostat=ios) &
           Conv_Cool_Dens, Conv_Cool_Temp, Conv_Fuel_Temp

        call Iostat_Error_Check&
      (ios,"Error in Reading Conv_Cool_Dens, Conv_Cool_Temp, &
                  Conv_Fuel_Temp from the FILE_MAP")

      end if

! PVM_MAP_NTHC
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_MAP_NTHC", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier PVM_MAP_NTHC in the FILE_MAP,  '//&
         'mapping matrices from neutronics spatial mesh into '//&
         'heat conduction spatial mesh')
      else

! input of 2D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NN_RT_HC_TRAC, &
          2*N_POLY, ia_2D_NTHC, ja_2D_NTHC, Map_2D_NTHC, &
         " 2D: Neutronics => Heat Conduction" )
! input of 1D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NN_Z_HC_TRAC, &
          2*NZ, ia_1D_NTHC, ja_1D_NTHC, Map_1D_NTHC, &
         " 1D: Neutronics => Heat Conduction" )

      end if

! PVM_MAP_NTFD
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_MAP_NTFD", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier PVM_MAP_NTFD in the FILE_MAP,  '//&
         'mapping matrices from neutronics spatial mesh into '//&
         'fluid dynamics spatial mesh')
      else

! input of 2D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NN_RT_FD_TRAC, &
          2*N_POLY, ia_2D_NTFD, ja_2D_NTFD, Map_2D_NTFD, &
         " 2D: Neutronics => Fluid Dynamics" )
! input of 1D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NN_Z_FD_TRAC, &
          2*NZ, ia_1D_NTFD, ja_1D_NTFD, Map_1D_NTFD, &
         " 1D: Neutronics => Fluid Dynamics" )

      end if

! PVM_MAP_HCNT
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_MAP_HCNT", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier PVM_MAP_HCNT in the FILE_MAP,  '//&
         'mapping matrices from heat conduction spatial mesh '//&
         'into neutronics spatial mesh; Taking TRANSPOSE')

         call csrcsc2 (NN_RT_HC_TRAC, N_POLY,1,1, &
             Map_2D_NTHC,ja_2D_NTHC,ia_2D_NTHC,&
             Map_2D_HCNT,ja_2D_HCNT,ia_2D_HCNT)

         call csrcsc2 (NN_Z_HC_TRAC, NZR, 1,1, &
             Map_1D_NTHC,ja_1D_NTHC,ia_1D_NTHC,&
             Map_1D_HCNT,ja_1D_HCNT,ia_1D_HCNT)
      else

! input of 2D Matrix
        CALL PVM_read_mapping_matrix(io_unit, N_POLY, &
          2*N_POLY, ia_2D_HCNT, ja_2D_HCNT, Map_2D_HCNT, &
         " 2D: Heat Conduction => Neutronics" )
! input of 1D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NZR, &
          2*NZ, ia_1D_HCNT, ja_1D_HCNT, Map_1D_HCNT, &
         " 1D: Heat Conduction => Neutronics " )

      end if


! PVM_MAP_FDNT
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "PVM_MAP_FDNT", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier PVM_MAP_HCNT in the FILE_MAP,  '//&
         'mapping matrices from fluid dynamics spatial mesh '//&
         'into neutronics spatial mesh; Taking TRANSPOSE')

         call csrcsc2 (NN_RT_FD_TRAC, N_POLY, 1,1, &
             Map_2D_NTFD,ja_2D_NTFD,ia_2D_NTFD,&
             Map_2D_FDNT,ja_2D_FDNT,ia_2D_FDNT)

         call csrcsc2 (NN_Z_FD_TRAC, NZR,1,1, &
             Map_1D_NTFD, ja_1D_NTFD, ia_1D_NTFD,&
             Map_1D_FDNT, ja_1D_FDNT, ia_1D_FDNT)
      else

! input of 2D Matrix
        CALL PVM_read_mapping_matrix(io_unit, N_POLY, &
          2*N_POLY, ia_2D_FDNT, ja_2D_FDNT, Map_2D_FDNT, &
         " 2D: Neutronics => Heat Conduction" )
! input of 1D Matrix
        CALL PVM_read_mapping_matrix(io_unit, NZR, &
          2*NZ, ia_1D_FDNT, ja_1D_FDNT, Map_1D_FDNT, &
         " 2D: Neutronics => Heat Conduction" )

      end if

      close(io_unit)

      end if

      return
      end



      SUBROUTINE PVM_Output_Channels(time)    
!=====================================================================*
!      Output of the TRAC Channel Data  into the file                 *
!              "Output/Chan_TRAC.dat"                                 *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Input: time
      real time

      character*100 title, file_name


! output ALL TRAC data
      file_name = "Output/Chan_TRAC.dat"
      title = "Channel Power Distribution"
      CALL OUT_Write_Channel_Real8(file_name, &
         title , NN_RT_HC_TRAC, NN_Z_HC_TRAC, p_hc_trac, io_unit)

      title = "Channel Doppler Temperature"
      CALL OUT_Write_Channel_Real8(file_name, &
        title, NN_RT_HC_TRAC, NN_Z_HC_TRAC, tf_trac, io_unit)

      title = "Channel Coolant Temperature"
      CALL OUT_Write_Channel_Real8(file_name, &
        title , NN_RT_FD_TRAC, NN_Z_FD_TRAC, tc_trac, io_unit)

      title =  "Channel Coolant Density" 
      CALL OUT_Write_Channel_Real8(file_name, &
        title ,NN_RT_FD_TRAC, NN_Z_FD_TRAC,  ro_c_trac, io_unit)


      title =  "Channel Fuel Centerline Temperature" 
      CALL OUT_Write_Channel_Real8(file_name, &
        title ,NN_RT_HC_TRAC, NN_Z_HC_TRAC,  tf_trac_cl, io_unit)

      title =  "Channel Fuel Surface Temperature" 
      CALL OUT_Write_Channel_Real8(file_name, &
        title ,NN_RT_HC_TRAC, NN_Z_HC_TRAC,  tf_trac_sf, io_unit)

      RETURN
      END



      SUBROUTINE PVM_Output_Mapping_Parameters(unit)
      implicit none
      include 'sketch.fh'

      integer unit

      write(unit, *)

      write(unit, '(A, A)')&
       "   Version of the TRAC code: ", TRAC_Version
      write(unit, '(A)')&
       "   Parameters of the Mapping to/from the TRAC code:"
      write(unit, '(A, I8)') &
       "       Number of the TRAC heat conduction channels    :",&
       NN_RT_HC_TRAC
      write(unit, '(A, I8)') &
       "       Number of the TRAC fluid dynamics channels     :", &
       NN_RT_FD_TRAC
      write(unit, '(A, I8)') &
       "       Number of the TRAC heat conduction axial layers:", &
       NN_Z_HC_TRAC
      write(unit, '(A, I8)') &
       "       Number of the TRAC fluid dynamics axial layers :", &
       NN_Z_FD_TRAC

      RETURN 
      END

      SUBROUTINE PVM_Output_Mapping_Data(unit)
      implicit none
      include 'sketch.fh'
! Input:
      integer unit
! Local
! * Checking the Mapping Matrix
      INTEGER NXY_MAX_TRAC, NXY_MAX, NZ_MAX_TRAC, NZ_MAX, N_MAX
      PARAMETER ( NXY_MAX_TRAC = MAX(NN_RT_FD_TRAC, NN_RT_HC_TRAC) )
      PARAMETER ( NXY_MAX = MAX(NXY_MAX_TRAC, N_POLY) )
      PARAMETER ( NZ_MAX_TRAC = MAX(NN_Z_FD_TRAC, NN_Z_HC_TRAC) )
      PARAMETER ( NZ_MAX = MAX(NZ_MAX_TRAC, NZR) )
      PARAMETER ( N_MAX = MAX(NXY_MAX, NZ_MAX) )
      real  map_unity_vector(N_MAX)
!      , map_HC_pow_z(NZR),&
!            map_FD_pow_xy(N_POLY), map_FD_pow_z(NZR),&
!            map_HC_xy(N_POLY), map_HC_z(NZR),&
!            map_FD_xy(N_POLY), map_FD_z(NZR),&
!            sum_HC_pow, sum_HC_pow_z, sum_HC, sum_HC_Z,&
!            sum_FD, sum_FD_Z, sum_FD_pow, sum_FD_pow_z

      real Identity(N_MAX)
      data Identity /N_MAX*1./

      CHARACTER*80 HEADER


      write(unit, '(A)') " Data of the Mapping o/from the TRAC code:"

!"     Conversion Coefficients for Coolant Temperature, "
      WRITE(unit, '(A)') &
       "     Conversion Coefficients for Coolant Temperature, "
      WRITE(unit, '(A)') &
       "        Coolant Temperature and Fuel Temperature:"
     
      WRITE(unit, '(8x,3E13.5)') Conv_Cool_Dens, Conv_Cool_Temp, &
        Conv_Fuel_Temp

!"     Mapping Power into TRAC Heat Conduction Mesh (x-y)"
        HEADER = &
      "     Mapping Power into TRAC Heat Conduction Mesh (x-y) :"
        CALL PVM_Output_Mapping_Matrix(unit, header, NN_RT_HC_TRAC,&
       ia_2D_NTHC, ja_2D_NTHC, Map_2D_NTHC, &
       identity,  map_unity_vector, .False.)

!"     Mapping Power into TRAC Heat Conduction Mesh (z)"
        HEADER = &
      "     Mapping Power into TRAC Heat Conduction Mesh (z) :"
        CALL PVM_Output_Mapping_Matrix(unit, header, NN_Z_HC_TRAC,&
       ia_1D_NTHC, ja_1D_NTHC, Map_1D_NTHC, &
       identity,  map_unity_vector, .False.)

!"     Mapping Power into TRAC Fluid Dynamics Mesh (x-y)"
        HEADER = &
      "     Mapping Power into TRAC Fluid Dynamics Mesh (x-y)"
        CALL PVM_Output_Mapping_Matrix(unit, header, NN_RT_FD_TRAC,&
       ia_2D_NTFD, ja_2D_NTFD, Map_2D_NTFD, &
       identity,  map_unity_vector, .False.)

!"     Mapping Power into TRAC Fluid Dynamics Mesh (Z)"
        HEADER = &
      "     Mapping Power into TRAC Fluid Dynamics Mesh (x-y)"
        CALL PVM_Output_Mapping_Matrix(unit, header, NN_Z_FD_TRAC,&
       ia_1D_NTFD, ja_1D_NTFD, Map_1D_NTFD, &
       identity,  map_unity_vector, .False.)



!"     Mapping Fuel Temperature into SKETCH Spatial Mesh (x-y)"
        HEADER = &
      "     Mapping Fuel Temperature into SKETCH Spatial Mesh (x-y)"
        CALL PVM_Output_Mapping_Matrix(unit, header, N_POLY,&
       ia_2D_HCNT, ja_2D_HCNT, Map_2D_HCNT, &
       identity,  map_unity_vector, .True.)

!"     Mapping Fuel Temperature into SKETCH Spatial Mesh (z)"
        HEADER = &
      "     Mapping Fuel Temperature into SKETCH Spatial Mesh (z)"
        CALL PVM_Output_Mapping_Matrix(unit, header, NZR,&
       ia_1D_HCNT, ja_1D_HCNT, Map_1D_HCNT, &
       identity,  map_unity_vector, .False.)

!"     Mapping Coolant Density into SKETCH Spatial Mesh (x-y)"
        HEADER = &
      "     Mapping Coolant Density into SKETCH Spatial Mesh (x-y):"
        CALL PVM_Output_Mapping_Matrix(unit, header, N_POLY,&
       ia_2D_FDNT, ja_2D_FDNT, Map_2D_FDNT, &
       identity,  map_unity_vector, .True.)

!"     Mapping Coolant Density into SKETCH Spatial Mesh (z):"
        HEADER = &
      "     Mapping TRACer into TRAC Heat Conduction Mesh (z) :"
        CALL PVM_Output_Mapping_Matrix(unit, header, NZR, &
       ia_1D_FDNT, ja_1D_FDNT, Map_1D_FDNT, &
       identity,  map_unity_vector, .False.)


      RETURN 
      END


      SUBROUTINE PVM_Output_Mapping_Matrix(unit, header, NN_RT_HC_TRAC,&
       ia_2D_NTHC, ja_2D_NTHC, Map_2D_NTHC, &
       identity,  map_unity_vector, OUT_MAP)

      IMPLICIT NONE
! Input:
      INTEGER unit,  NN_RT_HC_TRAC
      INTEGER ia_2D_NTHC(NN_RT_HC_TRAC+1), ja_2D_NTHC(*)
      REAL identity(*), map_unity_vector(*),  Map_2D_NTHC(*)
      CHARACTER*(*) header
      LOGICAL OUT_MAP
! Local:
      INTEGER none_zero, np
      CHARACTER*80 header_map
      CHARACTER*7 fmt_real, fmt_char
          
! "        Mapping Matrix Z:"
      WRITE(unit,'(A)') &
       header
      WRITE(unit, '(A)') &
       "        Mapping Matrix :"
! "           Number of nonzero elements:", none_zero 
      none_zero =  ia_2D_NTHC(NN_RT_HC_TRAC+1) - 1
      WRITE(unit, '(A, I8)') &
       "           Number of nonzero elements:", none_zero 

! "        Row Pointers IA:"
      header_map = "        Row Pointers IA:"
      IF(out_map) THEN
        fmt_real = "(I5)"
        fmt_char = "A5"
        CALL OUT_Write_Map_Integer(unit, Header_Map, ia_2D_NTHC, &
        fmt_real, fmt_char)
      ELSE
         WRITE(unit, '(A)') &
          header_map
         WRITE(unit, '(8x, 10I7)') &
        (ia_2D_NTHC(np),np = 1,NN_RT_HC_TRAC+1)
      END IF

! "        Column Indices JA:"
      header_map = "        Column Indices JA:"
!      IF(out_map) THEN
!        fmt_real = "(I5)"
!        fmt_char = "A5"
!        CALL OUT_Write_Map_Integer(unit, Header_Map, ja_2D_NTHC, &
!        fmt_real, fmt_char)
!      ELSE
        WRITE(unit, '(A)') header_map
        WRITE(unit, '(8x, 10I7)') &
        (ja_2D_NTHC(np), np = 1, none_ZERO)
!      END IF

! "        Nonezero Matrix Elements AA :"
      header_map = "        Nonezero Matrix Elements AA :"
!      IF(out_map) THEN
!        fmt_real = "(F5.2)"
!        fmt_char = "A5"
!        CALL OUT_Write_Map_Real(unit, Header_Map, Map_2D_NTHC, &
!        fmt_real, fmt_char)
!      ELSE
        WRITE(unit, '(A)') header_map
        WRITE(unit, '(8x, 10F7.4)') &
        (Map_2D_NTHC(np), np = 1, none_ZERO)
!      END IF

      CALL amux(NN_RT_HC_TRAC, Identity, map_unity_vector, &
             Map_2D_NTHC, ja_2D_NTHC, ia_2D_NTHC)

! "        Check, Mapping of the unity distribution:"
      header_map = "        Check, Mapping of the unity distribution:"
      IF(out_map) THEN
        fmt_real = "(F5.2)"
        fmt_char = "A5"
        CALL OUT_Write_Map_Real(unit, Header_Map, map_unity_vector, &
        fmt_real, fmt_char)
      ELSE

        WRITE(unit, '(A)') header_map
        WRITE(unit, '(8x, 10F7.4)') &
        (map_unity_vector(np), np = 1, NN_RT_HC_TRAC)
      END IF

      RETURN
      END

      SUBROUTINE PVM_read_mapping_matrix(io_unit, M, nnz_max, &
      ia, ja, val, mapping_matrix_name)
!=====================================================================!
!     Input of the mapping matrix in CRS format {IA, JA, VAL}         !
! 16.VIII.1999 (c) Slava                                              ! 
!=====================================================================!

       
      IMPLICIT NONE
! input:
      INTEGER  M, NNZ_MAX, io_unit
      CHARACTER*(*) mapping_matrix_name
! output:
      INTEGER ia(M+1), ja(*)
      REAL val(*)
! local:
      INTEGER np, none_zero, ios


! reading pointers IA
        read(io_unit,fmt=*,iostat=ios) &
           (ia(np),np = 1, M+1)
        call Iostat_Error_Check&
      (ios,"Error in reading IA of "//mapping_matrix_name//&
           "from the FILE_MAP " )

        none_zero =  ia(M+1) - 1
        if(none_zero.GT.(NNZ_MAX) ) then
           write(*,*) &
      " Please, change the dimension of the arrays JA and VAL of "//&
      " the matrix "//mapping_matrix_name
           write(*,*) &
          'Number of nonzero elements ', none_zero
           write(*,*) &
          'Dimension of Arrays JA and Val = ', NNZ_MAX
           stop
        end if

! reading JA 
        read(io_unit,fmt=*,iostat=ios) &
         (ja(np), np = 1, none_zero)
        call Iostat_Error_Check&
      (ios,"Error in reading JA of "//mapping_matrix_name//&
           "from the FILE_MAP " )

! reading VAL
        read(io_unit,fmt=*,iostat=ios) &
         (val(np), np = 1, none_zero)
        call Iostat_Error_Check&
      (ios,"Error in reading VAL of "//mapping_matrix_name//&
           "from the FILE_MAP " )

      RETURN
      END


