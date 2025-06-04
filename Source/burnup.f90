      subroutine BRN_Input
!=====================================================================*
!        Input from the file FILE_INPUT                                *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Local Variables
      integer  line

      integer ios

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      character*100 Message
      logical error_find

!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'


      open(io_unit,file = FILE_INPUT,status='old', iostat=ios)

      call Iostat_Error_Check&
      (ios,"Could not find the input FILE_INPUT file "//FILE_INPUT)

! reading BRN_REF_POWR
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "BRN_REF_POWR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier BRN_REF_POWR in the FILE_INPUT, set '//&
        'set REF_POWER to 3000 MWt' )
            REF_POWER = 3000.
      else
         read(io_unit,  fmt=*, iostat=ios) ref_power
         call Iostat_Error_Check&
      (ios,"Error in reading reference reactor power "//&
       "under identifier BRN_REF_POWR from the FILE_INPUT") 
      end if


! reading BRN_FRM_TIME
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "BRN_FRM_TIME", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier BRN_FRM_TIME in the FILE_INPUT, set '//&
        'set form of the input time table to DATE' )
            brn_form_time_table = "DATE"
      else

        read (io_unit,fmt=*,iostat=ios) brn_form_time_table

        call Iostat_Error_Check&
        (ios,"error in reading BRN Model "//&
        " from the FILE_INP")

        IF( (brn_form_time_table .NE. "DATE").AND.&
           (brn_form_time_table .NE. "TIME") ) THEN
            WRITE(Message, '(A,A,A)') &
            "brn_form_time_table =",&
              xs_model, &
            " in FILE_INP,  known types , "//&
            " 'DATE' and 'TIME' "
            CALL MSC_ERR_Add_Message('ERROR',Message)
        END IF

      end if

      CLOSE(io_unit)

      RETURN
      END SUBROUTINE BRN_Input

      SUBROUTINE BRN_Initialize
      USE PRECISION   , ONLY : dp
      USE termsort_lib, ONLY : get_lib_hm_mass_of_fuel
      IMPLICIT NONE
      include 'sketch.fh'
      INTEGER  :: np, ns, n_type
      REAL(dp) :: hm_mass 

      IF( XS_MODEL.EQ."SUBSET".OR.XS_MODEL.EQ."SUBSED") THEN

      DO ns = 1, NZR
         DO np = 1, N_POLY
                 n_type = l(np,ns)
                 CALL get_lib_hm_mass_of_fuel(n_type, hm_mass)  
!                 IF( n_type > 10 ) THEN
!                 write(*,*) 'np, ns, n_type =', np, ns, n_type
!                 write(*,*) 'hm_mass =', hm_mass
!                 PAUSE  
!                 END IF 
                 hm_mass_of_fuel(np,ns) = hm_mass
     
           END DO 
        END DO
      ELSE
! temporary I do no know where they get heavy metal mass of fuel
! in their LINTAB
       hm_mass_of_fuel(1:N_POLY,1:NZR)  = 2.41163E-3
      END IF
   
      RETURN 
      END SUBROUTINE BRN_Initialize
      
      SUBROUTINE BRN_Compute_Burnup(days)
      implicit none
      include 'sketch.fh'

!      REAL, PARAMETER :: hm_mass_of_Fuel = 2.41163E-3
      REAL, PARAMETER :: convert_Wt_to_MWt = 1.E-06
! Heavy Metal mass [kg] of 1 cm^3 of FA
      REAL, INTENT(IN) :: days
      REAL ::  a1
!      , vol
      INTEGER :: n1, nn, kc, k, kt, ns, np
       

!      days = 11.16495
!        burnup = 0.
!        burnup = burnup + pow_ref*days*1.E-6/hm_mass_of_fuel

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         ns = ns_out(n1)
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
              kt = k + nn 
!              vol = volume(kt) 
              np = np_out(k)
              a1 = p(k,n1)*convert_Wt_to_MWt*days/&
                   hm_mass_of_fuel(np, ns)
!            write(*,*)  'np, ns =', np, ns
!              write(*,*)  'hm_mass_of_fuel(np, ns)=',&
!              hm_mass_of_fuel(np, ns)
!              PAUSE                
              brn(k,n1) = brn(k,n1) + a1         
!              brn_reactor = brn_reactor + brn(k,n1)*vol
         END DO
      END DO

!        = brn_reactor/v_core
!      brn_col(0,0) = brn_reactor

      RETURN
      END 

      SUBROUTINE BRN_Compute_Burnup_Distribution
      implicit none
      include 'sketch.fh'

      CALL Compute_Average_Core_Distr(NH, NZ, N_POLY, NZR,&
       NYR, NXR, NZ_Core_Beg, NZ_Core_End, NZR_Core_Beg, NZR_Core_End, &
       brn, brn_col, brn_mm, k_brn_mm, Index_Core )

      write(*,*) 'Average Core burnup  =', &
       brn_col(0,0), ' [MWt*day/kgHM]'

      RETURN
      END 


      subroutine BURnup_Distr_Output(unit)
!=====================================================================*
! Output of the burnup distribution into "SKETCH.lst"              *
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      character*85 Header_Map
      character*4 val_fmt
      character*6 val_char(0:N_POLY)
      integer ind, ns
!      integer nlx_core, nly_core
! external function (in "OUTput.f")
      integer nlz_core

      real dist_scaling_factor
!      data dist_th_scaling_factor / 1.E-3, 1.E-3, 1./

!  Scaling Factor for Fuel Enthalpy = 1./dist_th_av

      dist_scaling_factor=1.
      IF ( brn_col(0,0) > SMALL_VALUE ) THEN 
        dist_scaling_factor = &
        1./brn_col(0,0)
        END IF

      CALL OUTput_Write_Separator(unit)

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            BURNUP DISTRIBUTION"

!      DO i = 1, N_ISOTOPE

      WRITE(unit,'(A, E12.6, A)') &
        " Average Value of burnup =", &
         brn_col(0,0), " [MWt*day/kgHM]"

      CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
          dist_scaling_factor, &
          N_POLY, NZR, brn_col(0,0), brn_mm(-3), &
          k_brn_mm(1,-3) )

!brn_col, brn_mm, k_brn_mm
!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            2D BURNUP DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)
      
      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         brn_col(ind,0)*dist_scaling_factor
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

!      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL BURNUP DISTRIBUTIONS "
!      CALL OUTput_Write_Separator(unit)

      WRITE(unit, '(1x, 5A7, /)') "  N  ", &
       "BURNUP"

       do ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End

          WRITE(unit,'(1x, I3,": ", 4F8.4)') &
             nlz_core(ns), brn_col(0,ns)*dist_scaling_factor 
       end do

      CALL OUTput_Write_Separator(unit)

      RETURN
      END SUBROUTINE BURnup_Distr_Output 



      SUBROUTINE Compute_Average_Core_Distr(NH, NZ, N_POLY, NZR,&
       NYR, NXR, NZ_Core_Beg, NZ_Core_End, NZR_Core_Beg, NZR_Core_End, &
       p, p_col, p_mm, k_p_mm, Index_Core )

!=====================================================================*
! computing the average values over the reactor core                  *
! (c) Slava 4.III.1998 JAERI                                           *
!=====================================================================*
      IMPLICIT NONE

      INTEGER::NH, NZ, NYR, NXR,N_POLY, NZR, NZ_Core_Beg, NZ_Core_End, &
                NZR_Core_Beg, NZR_Core_End
      integer k_p_mm(3, -3:3), Index_Core(NYR, NXR)
      real :: p(NH,NZ), p_col(0:N_POLY,0:NZR), p_mm(-3:3)

      REAL weighting_factor(N_POLY)
!
      CALL MSC_SSET(N_POLY, 1.0, weighting_factor) 

      CALL OUTput_Compute_3D_Average( p, NZ_Core_Beg, NZ_Core_End, &
        Index_Core, p_col,  p_mm, k_p_mm, weighting_factor )

      CALL OUTput_Compute_2D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm )

      CALL OUTput_Compute_1D_Average(p_col, NZR_Core_Beg, NZR_Core_End, &
       Index_Core, p_mm, k_p_mm, weighting_factor )

      RETURN
      END