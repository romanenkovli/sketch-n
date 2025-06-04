      SUBROUTINE OUTput_NumResults(i_source,i_ssor,i_nonl,&
             i_therm, time, flag_initial_data)
!=====================================================================*
!      output of the power distribution  into the files               *
!              "Output/SKETCH.lst"  , "Output/Power3D.dat"            *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER flag_initial_data, i_source, i_ssor, i_nonl, i_therm
      REAL  time
! Local:
      REAL, EXTERNAL :: CONVERT_PPM_TO_GKG ! external function


!      real scaling_factor
      INTEGER iter_out

!      real peak_factor
!      character*72 Nodal_Method_Title


! Computing Power Distribution 
      IF(flag_initial_data.EQ.1) THEN
           CALL CRD_Compute_Mat_Comp
         CALL XS_Update 
         CALL POWer_Compute
         CALL THM_Set_Core_Power      
           CALL THM_Compute_Average_Feedbacks
           IF (TH_Model.EQ."Internal".OR.TH_Model.EQ."SKAZKA") THEN
                     CALL THM_Prepare_Output_Dist
           END IF

      END IF


! Starting output
      OPEN(io_unit,file='Output/SKETCH.lst', status='unknown',&
               access = 'Append')


      WRITE(*, '(/," Time =", E11.4, " sec", /)') time 

!      DO io = 6, io_unit, io_unit-6 ! Output into Screan and into the File pch_n
              
      CALL OUTput_Write_Header(io_unit)

      IF(flag_initial_data.EQ.1) THEN
         WRITE(io_unit,'(A)') "                 INITIAL DATA SET"
      ELSE
         WRITE(io_unit,'(A)') "                 NUMERICAL RESULTS"
      END IF

      IF(Problem_Type.NE."Kinetics") THEN
         WRITE(io_unit,'(A)')&
     "                 STEADY-STATE CALCULATIONS"
      ELSE
         WRITE(io_unit,'(A, A, E11.4, A)')&
     "                 TRANSIENT CALCULATION:"," Time =",time," sec"
      END IF

      CALL OUTput_Write_Separator(io_unit)

      IF(Problem_Type.NE."Kinetics") THEN
        WRITE(io_unit,'(A, F9.6)')&
     "     k_ef                                  : ", k_ef
        IF(Steady_State_Type.EQ."BoronSearch") THEN
           WRITE(io_unit,'(A,ES12.5,A,ES12.5,A)')&
     "     Boron Concentration                   : ",&
       fdback(1,1,1), " ppm", CONVERT_PPM_TO_GKG(fdback(1,1,1)),&
        " g/kg"
        END IF
      END IF


!     Control Rod Positions 
      CALL OUTput_Write_Separator(io_unit)

      CALL CRD_Output_CR_Positions(io_unit)

      CALL POWer_Output(io_unit)

      CALL THM_Feedbacks_Output(io_unit)
      IF (TH_Model.EQ."Internal".OR.TH_Model.EQ."SKAZKA") THEN
            CALL THM_Internal_Distr_Output(io_unit)
          CALL TH_OUT_Average_Fuel_Rods(io_unit,&
                NZR_Core_Beg, NZR_Core_End) 
      END IF



      CALL POW_Neutron_Flux_Output(io_unit)

      IF(flag_initial_data.NE.1) THEN
        IF( NonlinearIterations.EQ."Moon") THEN 
           CALL OUTput_Average_Surface_Flux
           CALL OUTput_Average_Surface_Current
        END IF
      END IF

! burnup 
         CALL BRN_Compute_Burnup_Distribution
        CALL BURnup_Distr_Output(io_unit)

!      IF( Xe_Sm_Model /= "nn") THEN 
       IF(flag_initial_data.NE.1) THEN
       CALL Xe_Sm_Compute_Distribution 
       CALL Xe_Sm_Distr_Output(io_unit)
       END IF       
!      END IF

      IF(flag_initial_data.NE.1) THEN

         WRITE(io_unit,'(" Number of Thermal-Hydraulics Calculations =",&
             I8)') i_therm 
         WRITE(io_unit,'(" Number of Nonlinear Iterations            =",&
             I8)') i_nonl

         IF(i_nonl.NE.0) THEN
             iter_out = i_source/i_nonl
         ELSE
            iter_out = i_source
         END IF                                    
         WRITE(io_unit,'(" Number of Source Iteration per Nonlinear  =",&
              I8)') iter_out

         IF(i_source.NE.0) THEN
          iter_out = i_ssor/i_source
         ELSE
          iter_out = i_ssor
         END IF              
         WRITE(io_unit,'(" Number of SSOR iterations  per Outer      ="&
              I8)') iter_out

         WRITE(io_unit,'(" n_inter =", I3, " e_inter =", E10.5)')&
             n_inter, e_inter
!           write(io,'(" xme_ini =", F8.5, "npolins =", I3)')&
!             xme_ini, npolins
         WRITE(io_unit,'(" Convergence  of k_eff (actual, tolerance) = "&
            , 2E12.5)') d_kef_l, e_outer_l
         WRITE(io_unit,'(" Convergence  of flux  (actual, tolerance) = "&
            , 2E12.5)') d_flux_l, e_flux_l
      END IF

      CLOSE(io_unit)

!      if(NZ.ne.1) then
! Output of 3D Power Distribution

!      file_name = "Output/Power3D.dat"
!      title = "3D Power Distribiution in Core (pow/p_average) "
!      scaling_factor = 1./p_average

!      CALL OUT_Write_Channel_Real&
!         (file_name, time, title, N_POLY, NZR, p_col, scaling_factor,&
!         io_unit)

!      open(UNIT = io_unit, FILE = file_name, ACCESS ='Append')

!      write(io_unit, '(/, " Average Power Density = ", E12.5,&
!             " Wt/cm^3")') p_average
!      write(io_unit, '(   " Total Reactor Power   = ", E12.5,&
!             " MWt")') p_total
!      write(io_unit,'(" 3D Nodal Power Peaking Factor = ", F10.5)')&
!        p3d_max/p_average
!      write(io_unit,'(" Location of 3D Maximum (Channel Number; nz): "&
!       2I3)') npoly(k_p3d_max(2),k_p3d_max(1)) , k_p3d_max(3) 
!      close(io_unit)
!      end if

!      CALL OUTput_Feedbacks_Dist(time)

!      IF(TH_Model.EQ."External") CALL PVM_Output_Channels(time)

      RETURN
      END

      SUBROUTINE OUT_Write_Map_Errors_Real(Header_Map, io_out,&
          Dist_av_rin, Dist_1d_rin, Dist_2d_rin, Scale_Factor_Rin,&
          Dist_av, Dist_1d, Dist_2d, Scale_Factor,&
          Error_2d, Error_1D, Error_2d_Max, Error_2d_Av,&
          Error_1d_Max, Error_1d_Av, error_av, k_error_2D_max,&
          k_error_1D_max)

!=====================================================================*
! Output of a compariosn 2D and 1D Distributions                      *
!     25.VI.1999 (c) Slava                                            * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      CHARACTER*(*) Header_Map
      INTEGER io_out
      REAL Dist_av_rin, Dist_1d_rin(NZR_Core),&
          Dist_2d_rin(NP_Reactor_Core), Scale_Factor_Rin
      REAL Dist_av, Dist_1d(NZR),&
          Dist_2d(N_POLY), Scale_Factor
      REAL Error_2d(NP_reactor_Core), Error_1D(NZR_Core),&
          Error_2d_Max, Error_2d_Av, Error_1d_Max, Error_1d_Av
      REAL error_av
      INTEGER k_error_2D_max(2), k_error_1D_max

! Local Variables:
      INTEGER N_Blank, N_Position
      CHARACTER*5 C_Blank, C_Position
      INTEGER nly, nlx,  nc
      CHARACTER*100 fmt


      CALL OUTput_Write_Separator(io_unit)
        WRITE(io_out, '(A)') Header_Map
      CALL OUTput_Write_Separator(io_unit)

!2D Power distributio_unitn
      WRITE(io_out,*)
      CALL OUTput_Write_Separator(io_out)
      WRITE(io_out,*) 'Normalized 2D Distribution Value/Value_Av'
      CALL OUTput_Write_Separator(io_out)

      N_Position = Nxr_Max_Core - Nxr_B_Min_Core + 1

      WRITE(C_Position, '(I5)') N_Position
      fmt = '(/, 5x, '//C_Position//'I5)'

      WRITE(io_out, fmt )&
      (nlx-Nxr_B_Min_Core+1,  nlx = Nxr_B_Min_Core, Nxr_Max_Core)

      WRITE(io_out, '(5x,'//C_Position//'("  ---"))')

      DO nly = Nyr_B_Core, Nyr_E_Core
        N_Blank = Nxr_B_Core(nly) - Nxr_B_Min_Core
        N_Position = Nxr_E_Core(nly) - Nxr_B_Core(nly) + 1

      WRITE(C_Position, '(I5)') N_Position
      WRITE(C_Blank, '(I5)') N_Blank*5 + 1
      

      WRITE(io_out, '(I3,":", '//C_Blank//'x,'//C_Position//'F5.1)')&
       nly - Nyr_B_Core + 1,&
        (Dist_2D(npoly(nly,nlx))*Scale_Factor,&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))

      WRITE(io_out, '(I3,":", '//C_Blank//'x,'//C_Position//'F5.1)')&
       nly - Nyr_B_Core + 1,&
       (Dist_2d_rin(Numb_poly_core(npoly(nly,nlx)))*Scale_Factor_Rin,&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))

      WRITE&
       (io_out,'(I3,":",'//C_Blank//'x,'//C_Position//'F5.1)')&
         nly - Nyr_B_Core + 1,&
         (Error_2d(Numb_poly_core(npoly(nly,nlx))),&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))
        WRITE(io_out,*)

      END DO

      WRITE(io_out,*)&
        'Errors in 2d Power Distribution: AVERAGE, MAX(X,Y)'
      WRITE(io_out, '(15x, 2F8.2, 2I3)' ) Error_2d_Av, Error_2d_max,&
         k_error_2D_max

      WRITE(io_out,'(/,A,/)') ' Axial 1D Distribution: Ringhals,&
        SKETCH, diff(%)'
       DO nc = 1, NZR_CORE
          WRITE(io_out,'(I3,":   ", 2F8.4, F8.2)') nc,&
            Dist_1d_rin(nc)*Scale_Factor_Rin,&
            Dist_1D(nc + NZR_Core_Beg - 1)*Scale_Factor,&
            Error_1d(nc)
       END DO

      WRITE(io_out,*)&
        ' Errors in 1d Power Distribution: AVERAGE, MAX(Z)'
      WRITE(io_out, '(15x, 2F8.2, I3)' ) Error_1d_Av, Error_1d_max,&
         k_error_1d_max

      WRITE(io_out,&
      '(" Average Distributions (RINGHALS, SKETCH, diff(%) ):",&
        /,1x, 3E12.5)' )&
        Dist_av_rin, Dist_av, error_av

      RETURN
      END


      SUBROUTINE OUT_Write_Channel_Real8&
         (file_name, header, N_Poly, NZ, p_trac, io_out)
      IMPLICIT NONE
!=====================================================================*
!  WRITE TRAC Channel DATA (REAL*8) into the file                     *
!          (c) Slava 25.II.1999                                       *
!=====================================================================*
! Input:
      CHARACTER*(*) file_name, header
      INTEGER N_Poly, NZ
      INTEGER io_out
      REAL*8 p_trac(N_POLY, NZ)     
! Local Variable
      INTEGER N_In_Line ! Number of data in 1 line 
!       (when changing change FORMAT)
      PARAMETER (N_In_Line = 8)
      INTEGER N_Out, N_Left, N_Cycle
      INTEGER nc, no, n1

      CHARACTER*100 fmt
      CHARACTER*2 char_n_out
      
      OPEN(io_out, file = File_Name, status = 'unknown', access=&
                  'append')
                     
      CALL OUTput_Write_Header(io_out)
               
      WRITE(io_out, '(A)') header
        
      N_Cycle = N_poly/N_In_Line + 1
      N_Left = MOD(N_poly,N_In_Line)
      IF(N_Left.EQ.0) THEN
          N_Cycle = N_Cycle - 1
          N_Left = N_In_Line
      END IF
      N_Out = N_In_Line
        
      DO nc = 1, N_Cycle

           IF(nc.EQ.N_Cycle) N_out = N_Left

           WRITE(char_n_out, '(I2)' )  n_out

           WRITE(io_out,*)           
           WRITE(io_out,'(" CHAN", 8I15)') (no+(nc-1)*N_In_Line,&
                no=1, N_OUT)
           fmt = '(5x,'//char_n_out//'("  _____________"))'
           WRITE(io_out, FMT = fmt)

           DO n1=1, NZ
             WRITE(io_out, '(1x,I3,":",8E15.7)') n1,&
            (REAL(p_trac(no+(nc-1)*N_In_Line, n1)),no=1, N_OUT)
           END DO
      END DO
      CLOSE(io_out)

      RETURN
      END                       

      SUBROUTINE OUT_Write_Map(N_POLY_L,&
        NXR_Beg_Min, NXR_Max,&
        NXR_Beg, NXR_End, NYR_Beg, NYR_End, npoly_l,&
        Header_Map, io_out, Value, val_fmt)
!=====================================================================*
!             WRITE CHARACTER Array Map (Core or Reactor)             *
!                   Vyachreslav Zimin (c) 31 May 2000                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      IMPLICIT NONE
      include 'sketch.fh'
! Input:
      INTEGER N_POLY_L !, NYR_L, NXR_L
      INTEGER NXR_Beg_Min, NXR_Max, NXR_Beg(NYR), NXR_End(NYR),&
          NYR_Beg, NYR_End, npoly_l(NYR, NXR)
      CHARACTER*(*) value(0:N_POLY_L) ! Output Value
      CHARACTER*(*) val_fmt
      INTEGER io_out
      CHARACTER*(*) Header_Map

! Local Varia      
      INTEGER  N_Position, N_Blank, nly, nlx, N_Digits, N_Blank_HEX
      CHARACTER*5 C_Blank, C_Position
      CHARACTER Char_Digits
      CHARACTER*100 fmt
      
      integer n_hex_0 , nx_blank

      value(0) = ""

!      CALL OUTput_Write_Separator(io_out)
      WRITE(io_out,'(/,A)') Header_Map
!      CALL OUTput_Write_Separator(io_out)

      READ(val_fmt, '(1X, I1)') N_Digits
      N_Position = NXR_Max - NXR_Beg_Min + 1

         WRITE(C_Position, '(I5)') N_Position
         WRITE(Char_Digits, '(I1)') N_Digits

          IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
           N_Blank_HEX = (Nyr_End - Nyr_Beg)*N_Digits/2
           n_hex_0 = 1000
           DO nly = NYR_BEG, NYR_END
              nx_blank = (Nyr_End - nly)*N_Digits/2 +&
                  (NXR_Beg(nly)-NXR_Beg_Min)*N_Digits
              n_hex_0 = min(n_hex_0, nx_blank)
           END DO
           N_Blank_HEX = N_Blank_HEX - n_hex_0
        ELSE
           N_Blank_HEX = 0
        END IF 
        WRITE(C_Blank, '(I5)') N_Blank_HEX

        IF( N_Blank_HEX /= 0) THEN
         fmt = '(/,5x,'//C_Blank//'x,'//C_Position//&
         'I'//Char_Digits//')'
         ELSE
 
         fmt = '(/,5x,'//C_Position//&
         'I'//Char_Digits//')'
         END IF 

!         write(*,'(A)') 'fmt =', fmt
!         write(*,*) 'Nxr_Beg_Min, Nxr_Max', Nxr_Beg_Min, Nxr_Max
!         read(*,*) 

         WRITE(io_out,fmt)&
                (nlx,  nlx = Nxr_Beg_Min, Nxr_Max)

         WRITE(Char_Digits, '(I1)') N_Digits - 1

        IF( N_Blank_HEX /= 0) THEN
         fmt = '(5x,'//C_Blank//'x,'//C_Position//'(1x,'//&
          Char_Digits//'("-")))'
        ELSE 
         fmt = '(5x,'//C_Position//'(1x,'//&
          Char_Digits//'("-")))'
        END IF 

!         write(*,'(A)') 'fmt =', fmt
         WRITE(io_out, fmt) 
 
         DO nly = Nyr_Beg, Nyr_End
            N_Blank = Nxr_Beg(nly) - Nxr_Beg_Min
            N_Position = Nxr_End(nly) - Nxr_Beg(nly) + 1
              IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
              N_Blank_HEX = (Nyr_End - nly)*N_Digits/2 - n_hex_0
           ELSE
              N_Blank_HEX = 0
           END IF 

           WRITE(C_Position, '(I5)') N_Position
           WRITE(C_Blank, '(I5)') N_Blank*N_Digits+1+N_Blank_HEX

            fmt ='(I3,":",'//C_Blank//'x,'//C_Position//val_fmt//')'
            WRITE(io_out, fmt)&
                nly,(value(npoly_l(nly,nlx)),&
                     nlx=NXR_Beg(nly), NXR_End(nly))


         END DO ! nly

      WRITE(io_out,*)
!      CALL OUTput_Write_Separator(io_out)

      RETURN
      END

      SUBROUTINE OUT_Write_Map_Real(unit, Header_Map, value,&
        fmt_real, fmt_char)

      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      CHARACTER*(*) Header_Map, fmt_real, fmt_char
      INTEGER unit
      REAL value(N_POLY)
! Local:
      INTEGER ind
      CHARACTER*5 val_char(0:N_POLY)


      DO ind=1, N_POLY
         WRITE(val_char(ind), FMT=fmt_real ) value(ind)
      END DO

      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core,&
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, fmt_char)

      RETURN
      END

      SUBROUTINE OUT_Write_Map_Integer(unit, Header_Map, value,&
        fmt_integer, fmt_char)

      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      CHARACTER*(*) Header_Map, fmt_integer, fmt_char
      INTEGER unit
      INTEGER value(N_POLY)
! Local:
      INTEGER ind
      CHARACTER*5 val_char(0:N_POLY)


      DO ind=1, N_POLY
         WRITE(val_char(ind), FMT=fmt_integer ) value(ind)
      END DO

      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core,&
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, fmt_char)

      RETURN
      END

      SUBROUTINE OUTput_Write_Restart_File(Time)
!**********************************************************************
!       Output of the Neutron Flux DATA into the Restart File         *
!**********************************************************************
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Output into the File FILE_OUT:  Flux(NG, 0:N_TOT) - Neutron Flux
!                       k_ef - Eigenvalue
!                       trl_xyz - Tansverse Lekage in X - Y - Z directions
! + Output into The File FILE_D:  Flux(NG, 0:N_TOT) - Neutron Flux
!                       k_ef - Eigenvalue
!                       trl_xyz - Tansverse Lekage in X - Y - Z directions
!                       Prec(m,k) - Delayed Neutron Precursor Concentration
!                       xrods(NN_CRod) - position of the Control Rods
!                       Source(N_TOT) - the source Terms ??     
!                       p(NH, NZ) - Power Density 
!                       p_average - Average Power Density, 
!                       p_total - Total Reactor Power
!                       fd(NDAT) - Detectors DATA
!                       c_re(NDAT, MD) - Delayed Neutron Precursors for Detrector 
!                       time - current time
!                       dt_save - last value of the time step size
      REAL time
! Local Variables 
      INTEGER n,  k

      IF(FILE_DMP_OUT_ST.NE."") THEN
         OPEN(io_unit,file=FILE_DMP_OUT_ST,form='UNFORMATTED',&
          status='UNKNOWN')
         CALL OUTput_write_restart_steady_state_data(io_unit)
!         IF(TH_Model.EQ."Internal") THEN
         CALL THM_write_Data_Restart_File(io_unit)
!         END IF
         CALL OUTput_write_restart_burnup_data(io_unit)
      END IF

      IF(FILE_BRN_OUT.NE."") THEN
         OPEN(io_unit,file=FILE_BRN_OUT,form='UNFORMATTED',&
          status='UNKNOWN')
         CALL OUTput_write_restart_burnup_data(io_unit)
      END IF


      IF(FILE_DMP_OUT_KIN.NE."") THEN
      OPEN(io_unit,file = FILE_DMP_OUT_KIN, STATUS='Unknown',&
        form='UNFORMATTED')
           CALL OUTput_write_restart_steady_state_data(io_unit)
!         IF(TH_Model.EQ."Internal") THEN
           CALL THM_write_Data_Restart_File(io_unit)
!         END IF
           CALL OUTput_write_restart_burnup_data(io_unit)
         CALL OUTput_write_restart_kinetics_data(io_unit, time)
         CALL OUTput_write_restart_adjoint_flux(io_unit)
      CLOSE(io_unit)
      END IF

      RETURN
      END

      SUBROUTINE OUTput_Compute_3D_Average(var, NZ_Beg, NZ_End,&
        Index_Map, var_col, var_mm, k_mm, weighting_factor)
!=====================================================================*
! Computing Average & Maximum Values of the Power                     *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      REAL var(NH, NZ)
      INTEGER NZ_Beg, NZ_End, Index_Map(NYR, NXR)
      REAL weighting_factor(N_POLY)
! Output:
      REAL var_col(0:N_POLY, 0:NZR), var_mm(-3:3)
      INTEGER k_mm(3,-3:3)
! Local Variables
      INTEGER ns, k, np, kt, n1,  n_shift, nlx, nly
      INTEGER NZR_Beg, NZR_End
      REAL vol, vol_total


      NZR_Beg = ns_out(NZ_Beg)
      NZR_End = ns_out(NZ_End)


      DO np = 1, N_POLY
         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
         IF(Index_Map(nly, nlx) .NE. 0) THEN
            DO ns = NZR_Beg, NZR_End
                 var_col(np, ns) = 0.
            END DO ! NZR
         END IF
      END DO

      vol_total = 0.
      DO k = 1, NH
         np = np_out(k)
         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
         IF(Index_Map(nly, nlx) .NE. 0) THEN
            DO n1 = NZ_Beg, NZ_End
               ns = ns_out(n1)
               n_shift = (n1 - 1)*NH
               kt = k + n_shift
               vol = volume(kt)              
! Computing the average for all materials in the core              
               var_col(np, ns) = var_col(np, ns) + var(k,n1)*vol
               vol_total = vol_total + vol*weighting_factor(np)
           END DO ! n1
         END IF 
      END DO 

      var_mm(3) = 0.
      var_mm(-3) = Big_Value ! see in "units.fh"

      DO np = 1, N_POLY
         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
         IF(Index_Map(nly, nlx) .NE. 0) THEN
              DO ns = NZR_Beg, NZR_End
               var_col(0,0) = var_col(0,0) + var_col(np, ns)*&
                 weighting_factor(np)
                  var_col(np, ns) = var_col(np, ns)/&
                 ( vol_ass(np,ns))
                  IF( var_col(np, ns) .GT. var_mm(3) ) THEN
                  var_mm(3) = var_col(np, ns)
                  k_mm(1,3) = N_Coord(np, 1)
                  k_mm(2,3) = N_Coord(np, 2)
                  k_mm(3,3) = ns
                END IF
                IF(var_col(np, ns) .LT. var_mm(-3)) THEN
                   var_mm(-3) = var_col(np, ns)
                   k_mm(1,-3) = N_Coord(np, 1)
                   k_mm(2,-3) = N_Coord(np, 2)
                   k_mm(3,-3) = ns
                END IF
              END DO
          END IF
      END DO

      var_col(0,0) = var_col(0,0) / vol_total


      RETURN
      END 


      SUBROUTINE OUTput_Compute_2D_Average(var_col, NZR_Beg, NZR_End,&
       Index_Map,  var_mm, k_mm )
!=====================================================================*
! Computing 2D Average & Maximum Values of the Power                  *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      REAL var_col(0:N_POLY, 0:NZR) ! INOUT
      INTEGER NZR_Beg, NZR_End, Index_Map(NYR, NXR)
! Output:
      REAL  var_mm(-3:3) 
      INTEGER k_mm(3,-3:3)
! Local Variables
      INTEGER n1, np, nlx, nly
      REAL vol, vol_chan

      var_mm(2) = 0. 
      var_mm(-2) = 1.E+30

      DO np = 1, N_POLY
         var_col(np,0) = 0.
         nlx = N_Coord(np, 1)
         nly = N_Coord(np, 2)
! Check IF we are inside the Reactor Core 
         IF(Index_Map(nly, nlx) .NE. 0) THEN

            vol_chan = 0.
            DO n1 = NZR_Beg, NZR_End
               vol = vol_ass(np, n1)              
               var_col(np,0) = var_col(np,0) + var_col(np,n1)*vol
               vol_chan = vol_chan + vol
            END DO ! NZ
            IF(vol_chan .GT. SMALL_VALUE ) THEN 
               var_col(np,0) = var_col(np,0)/vol_chan
            END IF
               IF( var_col(np,0) .GT. var_mm(2) ) THEN
               var_mm(2) = var_col(np,0)
               k_mm(1,2) = nlx
               k_mm(2,2) = nly
            END IF
            IF(var_col(np,0) .LT. var_mm(-2) ) THEN
               var_mm(-2) = var_col(np,0)
               k_mm(1,-2) = nlx
               k_mm(2,-2) = nly
            END IF
         END IF ! Index_Map
      END DO ! NH

      RETURN
      END 


      SUBROUTINE OUTput_Compute_1D_Average(var_col, NZR_Beg, NZR_End,&
       Index_Map, var_mm, k_mm, weighting_factor )
!=====================================================================*
! Computing 1D Average & Maximum Values of the Power                  *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      REAL var_col(0:N_POLY, 0:NZR), weighting_factor(N_POLY)
      INTEGER NZR_Beg, NZR_End, Index_Map(NYR, NXR)
! Output:
      REAL var_mm(-3:3)
      INTEGER k_mm(3,-3:3)
! Local Variables
      INTEGER n1, np, nlx, nly
      REAL vol, vol_plane

        var_mm(1) = 0. 
        var_mm(-1) = 1.E+30
      

      DO n1 = NZR_Beg, NZR_End
         var_col(0,n1) = 0.
         vol_plane = 0.
         DO np = 1, N_POLY
            nlx = N_Coord(np, 1)
            nly = N_Coord(np, 2)
! Check IF we are inside the Reactor Core 
           IF(Index_Map(nly, nlx) .NE. 0) THEN
              vol = vol_ass(np, n1)              
! Computing the average of the fission materials
              var_col(0,n1) = var_col(0,n1) + var_col(np,n1)*vol*&
             weighting_factor(np) 
              vol_plane = vol_plane + vol*weighting_factor(np)
           END IF ! Index_Map
         END DO ! NH
         IF(vol_plane .GT. SMALL_VALUE ) THEN
            var_col(0,n1) = var_col(0,n1)/vol_plane
          END IF
         IF(var_col(0,n1) .GT. var_mm(1) ) THEN
            var_mm(1) = var_col(0,n1)
            k_mm(1,1) = n1
         END IF

         IF(var_col(0,n1) .LT. var_mm(-1) ) THEN
            var_mm(-1) = var_col(0,n1)
            k_mm(1,-1) = n1
         END IF

      END DO ! NZ

      RETURN
      END 


      SUBROUTINE OUTPut_Write_Separator(io)
!=====================================================================*
! Writes a separator line  for the SKETCH output files                *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INTEGER io

      WRITE(io, '( "!",78("-"), "!" )' )

      RETURN
      END

      SUBROUTINE OUTPut_Write_Separator_Bold(io)
!=====================================================================*
! Writes a separator line  for the SKETCH output files                *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INTEGER io

      CHARACTER*80 string

      WRITE(string, '( "!",78("="), "!" )' )

      WRITE(io, '(A)') string

      RETURN
      END

      SUBROUTINE OUTPut_Write_Separator_Bold_unformatted(io)
!=====================================================================*
! Writes a separator line  for the SKETCH output files                *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INTEGER io

      CHARACTER*80 string

      WRITE(string, '( "!",78("="), "!" )' )

      WRITE(io) string

      RETURN
      END



      SUBROUTINE OUTput_Write_Header(Unit)
!=====================================================================*
! Writes Header for the SKETCH output files                           *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INTEGER Unit

        CALL OUTPut_Write_Separator_Bold(Unit)
        WRITE(Unit,'(A)')&
      "!                                   SKETCH-N"//&
      "                                   !"
        WRITE(Unit, '(A)')&
       "!  Nodal Neutron Diffusion Code for Solving"//&
      " Steady-State and Kinetics Problems !"
        WRITE(Unit,'(A)')&
       "!              Version 1.0 (c)"//&
       " Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !" 
!        write(Unit, '(/," Time =", E11.4, " sec", /)') time 
        CALL CPU_write_date_and_time(Unit)
        CALL OUTPut_Write_Separator_Bold(Unit)

      RETURN
      END

      SUBROUTINE OUTput_Write_Header_unformatted(unit)
!=====================================================================*
! Writes Header for the SKETCH *.GRF files                           *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INTEGER unit

      INTEGER N_LINE_HEADER_TITLE
      PARAMETER ( N_LINE_HEADER_TITLE = 6 )

        WRITE(unit) N_LINE_HEADER_TITLE

        CALL OUTPut_Write_Separator_Bold_Unformatted(unit)
        WRITE(unit)&
      "!                                   SKETCH-N"//&
      "                                   !"
        WRITE(unit)&
       "!  Nodal Neutron Diffusion Code for Solving"//&
      " Steady-State and Kinetics Problems !"
        WRITE(unit)&
       "!              Version 1.0 (c)"//&
       " Slava 2000 e-mail: na.vzimin@na-net.ornl.gov    !" 
!        write(Unit, '(/," Time =", E11.4, " sec", /)') time 
        CALL CPU_write_date_and_time_unformatted(unit)
        CALL OUTPut_Write_Separator_Bold_unformatted(unit)

      RETURN
      END


      SUBROUTINE OUTput_Reference_Comparison(i_nonl, i_source)
!=====================================================================*
! A comparison of the SKETCH results WITH the reference solution      *
! taken from the file FILE_REFERENCE specified in the NAMELIST        *
!  'Input/SKETCH.INI' (ONLY IF File_Reference/="");                   *
! 22.II.1999 (c) Slava                                                * 
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input:
      INTEGER i_nonl, i_source
!      INTEGER i_nodal, i_nonl, i_source
!      INTEGER N_Poly, NXR, NYR, NYR_B_Core,&
!       NYR_E_Core, NXR_B_Core(NYR), NXR_E_Core(NYR),&
!       NXR_B_Min, NXR_Max_Core, npoly(NYR,NXR)

      REAL Pow_Sketch(N_Poly), pow_sketch_max
! Internal Files
      CHARACTER*5 C_Blank, C_Position

      CHARACTER*100 Comments
! from the file:
      REAL k_ef_ref, Pow_Ref(NYR, NXR)
! Output:

      INTEGER n_pow_error_max, N_FA(NYR,NXR)
! Error_Pow (%), Error_k_Ef(pcm)
      REAL Error_Pow(NYR,NXR), Error_k_ef, Error_Pow_Max,&
          Error_Pow_Av, Error_Pow_RMS, err_abs
      REAL  pow_ref_max
      INTEGER  :: n_pow_ref_max, n_pow_sketch_max
! Local variables 
      INTEGER nly, nlx, N_Blank, N_Position,  np
      REAL Vol_Core, vol
      REAL PERCENT, PCM
      PARAMETER (PERCENT = 1.E+02, PCM = 1.E+05)
      INTEGER N_Digits, n_hex_0 , nx_blank, N_Blank_HEX, i, iy
      CHARACTER*100 fmt
      REAL, EXTERNAL :: CONVERT_PPM_TO_GKG ! external function

       

      DO np = 1, N_Poly
         Pow_SKetch(np) = p_col(np, 0)/p_col(0,0)
      END DO


      OPEN(io_unit, File = File_Reference, Status = 'Old')

        READ(io_unit,'(A)') Comments ! Reading Header
        READ(io_unit,*) k_ef_ref
        DO nly = NYR_B_Core, NYR_E_Core
           READ(io_unit,*) (Pow_Ref(nly,nlx), nlx = NXR_B_Core(nly),&
                      NXR_E_Core(nly))
!           write(*,'(8F6.3)') (Pow_Ref(nly,nlx), nlx = NXR_B_Core(nly),&
!                      NXR_E_Core(nly))
!           write(*,*) 'nly=', nly
!           pause
        END DO
      CLOSE(io_unit)

      i = 0
      N_FA(1:NYR,1:NXR)  = 0
      DO  nly = NYR_B_Core, NYR_E_Core
         DO nlx = NXR_B_Core(nly), NXR_E_Core(nly)
         i = i + 1
         N_FA(nly,nlx)  = i
         END DO
      END DO 
           


      Error_k_ef = (k_ef/k_ef_ref - 1.)*PCM

      pow_ref_max = 0.
      pow_sketch_max = 0.
 
      Error_Pow_Max = 0.
      Error_Pow_Av  = 0.
      Error_Pow_RMS = 0. 
      Vol_Core = 0.
      DO nly = NYR_B_Core, NYR_E_Core
         DO nlx = NXR_B_Core(nly), NXR_E_Core(nly)
          vol = hx(nlx)*hy(nly)
          Vol_Core = Vol_Core + vol 
          Error_pow(nly,nlx) = (Pow_Sketch(npoly(nly,nlx))-&
                       Pow_Ref(nly,nlx))*PERCENT
          err_abs = (Pow_Sketch(npoly(nly,nlx))-&
                       Pow_Ref(nly,nlx))
          Error_Pow_Av  = Error_Pow_Av + ABS(err_abs)*vol
          Error_Pow_RMS = Error_Pow_RMS + (err_abs*err_abs)*vol

          IF(ABS(Error_pow(nly,nlx)) .GT. ABS(Error_pow_Max)) THEN
              Error_Pow_Max = Error_pow(nly,nlx)
              n_pow_error_max = N_FA(nly,nlx)
          END IF

          IF( Pow_Sketch(npoly(nly,nlx)) .GT.  pow_sketch_max) THEN
              pow_sketch_max = Pow_Sketch(npoly(nly,nlx))
              n_pow_sketch_max = N_FA(nly,nlx)
          END IF

          IF( Pow_Ref(nly,nlx) .GT.  pow_ref_max) THEN
              pow_ref_max = Pow_Ref(nly,nlx)
              n_pow_ref_max = N_FA(nly,nlx)
          END IF

         END DO
      END DO

      Error_Pow_Av = Error_Pow_Av*PERCENT/Vol_Core
      Error_Pow_RMS = SQRT( Error_Pow_RMS/(Vol_Core) )*PERCENT

      OPEN(io_unit, File = 'Output/Comparison.dat',status ='Unknown')

      CALL OUTput_Write_Header(io_unit)

      WRITE(io_unit,*) 'A Comparison of the SKETCH Results with  ',&
                  File_Reference
      WRITE(io_unit,*) Comments
        
      IF(Nodal_Method.EQ."MCFD") THEN
             WRITE(io_unit,*)&
             'MESH-CENTERED FINITE-DIFFERENCE METHOD'
      ELSE IF(Nodal_Method.EQ."SANM") THEN
             WRITE(io_unit,*) 'SEMI-ANALYTICAL NODAL METHOD'
      ELSE IF(Nodal_Method.EQ."PNM") THEN
             WRITE(io_unit,*) 'POLYNOMIAL NODAL METHOD, TRADITIONAL' 
      ELSE IF(Nodal_Method.EQ."ANM") THEN
             WRITE(io_unit,*) 'ANALYTICAL NODAL METHOD, '
      ELSE IF(Nodal_Method.EQ."PNM1") THEN
             WRITE(io_unit,*)&
           'POLYNOMIAL NODAL METHOD based on MATRIX FUNCTIONS'
      END IF

      WRITE(io_unit,*)

      IF( Steady_State_Type =="BoronSearch"  ) THEN  
       WRITE(io_unit,'(A,F10.0,F10.2)')&
      "critical boron concentration"//&
        " experiment [ppm, g/kg]  :",&
         k_ef_ref, CONVERT_PPM_TO_GKG(k_ef_ref)
       WRITE(io_unit,'(A,F10.0,F10.2)')&
      "critical boron concentration"//&
        " SKETCH-N   [ppm, g/kg]  :",&
        fdback(1,1,1), CONVERT_PPM_TO_GKG( fdback(1,1,1) ) 
       WRITE(io_unit,'(A,F10.0,F10.2)')&
        "Error                       "//&
        "            [ppm, g/kg]  :",&
        fdback(1,1,1)-k_ef_ref,&
        CONVERT_PPM_TO_GKG( fdback(1,1,1) ) -&
        CONVERT_PPM_TO_GKG(k_ef_ref)
      ELSE IF( Steady_State_Type == "Eigenvalue") THEN
       WRITE(io_unit,'("k_eff Reference        :",F10.5)')  k_ef_ref
       WRITE(io_unit,'("k_eff SKETCH-N         :",F10.5)')  k_ef
       WRITE(io_unit,'("Error, [pcm]           :",F7.0)' ) Error_k_ef
      END IF

      WRITE(io_unit, '(/,A,/)')&
        'Comparison of the 2D Power Distribution'

      N_Position = Nxr_Max_Core - Nxr_B_Min_Core + 1
      N_Digits = 6

          IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
           N_Blank_HEX = (Nyr_E_Core - Nyr_B_Core)*N_Digits/2
           n_hex_0 = 1000
           DO nly = Nyr_B_Core, Nyr_E_Core
              nx_blank = (Nyr_E_Core - nly)*N_Digits/2 +&
                  (Nxr_B_Core(nly) - Nxr_B_Min_Core)*N_Digits
              n_hex_0 = min(n_hex_0, nx_blank)
           END DO
           N_Blank_HEX = N_Blank_HEX - n_hex_0
        ELSE
           N_Blank_HEX = 0
        END IF 
        WRITE(C_Blank, '(I5)') N_Blank_HEX
        WRITE(C_Position, '(I5)') N_Position
      IF( N_Blank_HEX == 0) THEN
         fmt = '(/,5x,'//C_Position//&
         'I6)'
      ELSE
         fmt = '(/,5x,'//C_Blank//'x,'//C_Position//&
         'I6)'
      END IF
     

!      write(*,*) 'fmt =', fmt
!      write(*,*) (nlx,  nlx = Nxr_B_Min_Core, Nxr_Max_Core)
!      pause

      WRITE(io_unit,fmt)&
                          (nlx,  nlx = Nxr_B_Min_Core, Nxr_Max_Core)
      IF( N_Blank_HEX == 0) THEN
      WRITE(io_unit,&
        '(5x,'//C_Position//'("  ----"))')
      ELSE
      WRITE(io_unit,&
        '(5x,'//C_Blank//'x,'//C_Position//'("  ----"))')
      END IF 

      iy = 0
      DO nly = Nyr_B_Core, Nyr_E_Core
        N_Blank = Nxr_B_Core(nly) - Nxr_B_Min_Core
        N_Position = Nxr_E_Core(nly) - Nxr_B_Core(nly) + 1

              IF( gmt_crd_type(1:4).EQ."HEXZ" ) THEN
              N_Blank_HEX = (Nyr_E_Core - nly)*N_Digits/2 - n_hex_0
           ELSE
              N_Blank_HEX = 0
           END IF 

        WRITE(C_Position, '(I5)') N_Position
        WRITE(C_Blank, '(I5)') N_Blank*6+1+N_Blank_HEX 


        WRITE&
       (io_unit,'(I3,":",'//C_Blank//'x,'//C_Position//'I6)')&
                nly, (iy + 1 + nlx - Nxr_B_Core(nly),&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))
        iy = iy +  Nxr_E_Core(nly) - Nxr_B_Core(nly) + 1

        WRITE&
       (io_unit,'(I3,":",'//C_Blank//'x,'//C_Position//'F6.2)')&
                nly, (Pow_Ref(nly,nlx),&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))
        WRITE&
       (io_unit,'(I3,":",'//C_Blank//'x,'//C_Position//'F6.2)')&
             nly,(Pow_Sketch(npoly(nly,nlx)),&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))
        WRITE&
       (io_unit,'(I3,":",'//C_Blank//'x,'//C_Position//'F6.1)')&
           nly,(Error_Pow(nly,nlx),&
          nlx=Nxr_B_Core(nly), Nxr_E_Core(nly))
        WRITE(io_unit,*)
      END DO

      WRITE(io_unit,5)  Error_Pow_Av
      WRITE(io_unit,55) Error_Pow_RMS
      WRITE(io_unit,6) Error_Pow_Max
      WRITE(io_unit,7) n_pow_error_max

      write(io_unit,*) 
      WRITE(io_unit,16) Pow_ref_Max
      WRITE(io_unit,17) n_pow_ref_max
      WRITE(io_unit,26) Pow_sketch_Max
      WRITE(io_unit,27) n_pow_sketch_max
      WRITE(io_unit,36) (Pow_sketch_Max-Pow_ref_Max)*PERCENT

!      WRITE(io_unit,9) i_nonl
!      WRITE(io_unit,10) i_source

      CLOSE(io_unit)

    5 FORMAT( 'Average Error in Power Distribution, (abs*100): ', F6.1)
   55 FORMAT( '    RMS Error in Power Distribution, (abs*100): ', F6.1)
    6 FORMAT( 'Maximum Error in Power Distribution, (abs*100): ', F6.1)
    7 FORMAT( 'Number of FA with the Maximum Error           : ', I6)

   16 FORMAT( 'Maximum Reference Power Density               : ', F6.2)
   17 FORMAT( 'Number of FA with the Maximum Reference Power : ', I6)
   26 FORMAT( 'Maximum SKETCH-N  Power Density               : ', F6.2)
   27 FORMAT( 'Number of FA with the Maximum SKETCH-N Power  : ', I6)
   36 FORMAT( 'Error in Maximum Power Density       (abs*100): ', F6.1)

    9 FORMAT('Number of Nonlinear Iterations :', I4)
   10 FORMAT('Number of Source Iterations    :', i4,/)

      RETURN
      END

      SUBROUTINE OUTput_Distrb_Summary( unit, title, scale_factor,&
        N_POLY, NZR, dist_col, dist_mm, k_dist_mm)
!=====================================================================!
! Writing Distribution Summary into "SKETCH.lst"                      !
! (c) Slava 21.VII.2000 JAERI                                         !
!=====================================================================!
      IMPLICIT NONE

! Input:
      INTEGER unit, N_POLY, NZR 
      CHARACTER*(*) title
      REAL scale_factor, dist_col(0:N_POLY, 0:NZR), dist_mm(-3:3)
      INTEGER k_dist_mm(3, -3:3)

!LOCAL Valus
      INTEGER n, nlx_core, nly_core
! external function
      INTEGER nlz_core

      WRITE(unit, '(/,A)') title
!           feedb_name(i), "units in ", feedb_units(i)

      WRITE(unit, '(A,  E12.5)')&
     "     Scaling Factor                        : ",&
       scale_factor


      WRITE(unit, '(A,  F8.5)')&
     "     Average Value                         : ",&
       dist_col(0,0)*scale_factor

      WRITE(unit, '(/,A, F8.5)')&
     "     3D Maximum Value                      : ",&
       dist_mm(3)*scale_factor
      WRITE(unit,'(A, 3I3)')&
     "     Reactor Location of 3D Maximum (X,Y,Z):",&
           (k_dist_mm(n,3), n=1,3)

!      write(*,*) 'nlx =', k_dist_3D_max(1), 'nly=',&
!       k_dist_3D_max(2)

      CALL OUT_Convert_Coord_Core_xy&
          (k_dist_mm(1,3), k_dist_mm(2,3),&
           nlx_core, nly_core)
      WRITE(unit,'(A, 3I3)')&
     "     Core    Location of 3D Maximum (X,Y,Z):",&
           nlx_core, nly_core, nlz_core( k_dist_mm(3,3) )

      WRITE(unit, '(/,A, F8.5)')&
     "     3D Minimum Value                      : ",&
       dist_mm(-3)*scale_factor
      WRITE(unit,'(A, 3I3)')&
     "     Reactor Location of 3D minimum (X,Y,Z):",&
           (k_dist_mm(n,-3), n=1,3)
      CALL OUT_Convert_Coord_Core_xy&
          (k_dist_mm(1,-3), k_dist_mm(2,-3),&
           nlx_core, nly_core)
      WRITE(unit,'(A, 3I3)')&
     "     Core    Location of 3D minimum (X,Y,Z):",&
           nlx_core, nly_core, nlz_core( k_dist_mm(3,-3) )

      WRITE(unit,'(/,A, F8.5)')&
     "     2D Maximum Value                      : ",&
           dist_mm(2)*scale_factor
      WRITE(unit,'(A, 2I3)')&
     "     Reactor Location of 2D Maximum (X,Y)  :",&
           k_dist_mm(1,2), k_dist_mm(2,2)

      CALL OUT_Convert_Coord_Core_xy&
          (k_dist_mm(1,2), k_dist_mm(2,2),&
          nlx_core, nly_core)

      WRITE(unit,'(A, 2I3)')&
     "     Core    Location of 2D Maximum (X,Y)  :",&
           nlx_core, nly_core

      WRITE(unit,'(/,A, F8.5)')&
     "     2D Minimum Value                      : ",&
           dist_mm(-2)*scale_factor
      WRITE(unit,'(A, 2I3)')&
     "     Reactor Location of 2D minimum (X,Y)  :",&
           k_dist_mm(1,-2), k_dist_mm(2,-2)
      CALL OUT_Convert_Coord_Core_xy&
          (k_dist_mm(1,-2), k_dist_mm(2,-2),&
        nlx_core, nly_core)
      WRITE(unit,'(A, 2I3)')&
     "     Core    Location of 2D minimum (X,Y)  :",&
           nlx_core, nly_core

      WRITE(unit,'(/,A, F8.5)')&
     "     1D Maximum                            : ",&
           dist_mm(1)*scale_factor
      WRITE(unit,'(A,6x,I3)')&
     "     Reactor Location of the 1D Maximum (Z):",&
           k_dist_mm(1,1) 
      WRITE(unit,'(A,6x,I3)')&
     "     Core    Location of the 1D Maximum (Z):",&
           nlz_core ( k_dist_mm(1,1) )

      WRITE(unit,'(/,A, F8.5)')&
     "     1D Minimum                            : ",&
           dist_mm(-1)*scale_factor
      WRITE(unit,'(A,6x,I3)')&
     "     Reactor Location of the 1D minimum (Z):",&
           k_dist_mm(1,-1) 
      WRITE(unit,'(A,6x,I3)')&
     "     Core    Location of the 1D minimum (Z):",&
           nlz_core ( k_dist_mm(1,-1) )

      CALL OUTput_Write_Separator(unit)

      RETURN
      END

      SUBROUTINE OUTput_Convert_Dist_Core_Reactor(dist_core,&
        NZR_Beg, NZR_End, dist_col, dist_mm, k_dist_mm,&
        weighting_factor )

!=====================================================================!
! Converting Distribution Computed for the Reactor Core into Reactor  !
!     Distributions  (T/H Distributions)                              !
! (c) Slava 21.VII.2000 JAERI                                         !
!=====================================================================!

      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:  
      REAL Dist_core(NP_Reactor_Core, NZR)
      INTEGER NZR_Beg, NZR_End
      REAL weighting_Factor(N_POLY)  
! Output:
      REAL dist_col(0:N_POLY, 0:NZR)
      REAL dist_mm(-3:3)
      INTEGER k_dist_mm(3,-3:3)
! Local:
      INTEGER ns, nc, np
      REAL volume_total

      dist_mm(3) = 0.
      dist_mm(-3) = Big_Value ! see in "units.fh"
      dist_col(0,0) = 0.

      DO ns = NZR_Beg, NZR_End
        DO np = 1, N_POLY
         dist_col(np, ns) = 0.
        END DO
      END DO


      volume_total = 0. 
      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         DO ns = NZR_Beg, NZR_End
            IF( dist_core(nc, ns).GT.Small_Value) THEN 
            dist_col(np, ns) = dist_core(nc, ns)
            dist_col(0, 0) = dist_col(0,0) + dist_core(nc, ns)*&
              vol_ass(np, ns)*weighting_factor(np)
            volume_total = volume_total +&
              vol_ass(np, ns)*weighting_factor(np) 
               IF(dist_core(nc, ns) .GT. dist_mm(3) ) THEN
                  dist_mm(3) = dist_core(nc, ns)
                  k_dist_mm(1,3) = N_Coord(np, 1)
                  k_dist_mm(2,3) = N_Coord(np, 2)
                  k_dist_mm(3,3) = ns
            END IF
            IF(dist_core(nc, ns) .LT. dist_mm(-3)) THEN
                   dist_mm(-3) = dist_core(nc, ns)
                   k_dist_mm(1,-3) = N_Coord(np, 1)
                   k_dist_mm(2,-3) = N_Coord(np, 2)
                   k_dist_mm(3,-3) = ns
            END IF
            END IF ! IF( dist_core(nc, ns).GT.Small_Value) THEN 
         END DO
      END DO                  

!      write(*,*) 'volume_total =', volume_total
      IF ( volume_total /= 0. ) THEN 
        dist_col(0,0) = dist_col(0,0) / volume_total
      END IF 
           
      RETURN
      END 

      SUBROUTINE OUTput_Convert_Ring_Dist_Core_Reactor(dist_core,&
        dist_col, dist_mm, k_dist_mm )

!=====================================================================!
! Converting Distribution Computed for the Reactor Core into Reactor  !
!     Distributions  (Ringhals Distributions)                         !
! (c) Slava 21.VII.2000 JAERI                                         !
!=====================================================================!

      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:  
      REAL Dist_core(NZR_Core, NP_Reactor_Core)
! Output:
      REAL dist_col(0:N_POLY, 0:NZR)
      REAL dist_mm(-3:3)
      INTEGER k_dist_mm(3,-3:3)
! Local:
      INTEGER ns, nc, np, nsc

      dist_mm(3) = 0.
      dist_mm(-3) = Big_Value ! see in "units.fh"
      dist_col(0,0) = 0.

      DO nc = 1, NP_Reactor_Core
         np = Numb_Reactor_Core(nc)
         DO ns = NZR_Core_Beg, NZR_Core_End

            nsc = ns - NZR_Core_Beg + 1

            dist_col(np, ns) = dist_core(nsc, nc)
!            write(*,*) 'ns =', ns, "nsc=", nsc
!            pause

            dist_col(0,0) =  dist_col(0,0) + dist_core(nsc, nc)*&
              vol_ass(np, ns)

               IF(dist_core(nsc, nc) .GT. dist_mm(3) ) THEN
                  dist_mm(3) = dist_core(nsc, nc)
                  k_dist_mm(1,3) = N_Coord(np, 1)
                  k_dist_mm(2,3) = N_Coord(np, 2)
                  k_dist_mm(3,3) = ns
            END IF
            IF(dist_core(nsc, nc) .LT. dist_mm(-3)) THEN
                   dist_mm(-3) = dist_core(nsc, nc)
                   k_dist_mm(1,-3) = N_Coord(np, 1)
                   k_dist_mm(2,-3) = N_Coord(np, 2)
                   k_dist_mm(3,-3) = ns
            END IF
         END DO
      END DO                  

      dist_col(0,0) = dist_col(0,0) / v_core
           
      RETURN
      END 

      SUBROUTINE OUT_Convert_Coord_Core_xy&
          (nlx, nly, nlx_core, nly_core)
!=====================================================================*
! convert reactor coordinates into core coordinates (x, y)            !
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER nlx, nly
! Output:
      INTEGER nlx_core, nly_core

      if(nly.eq.0) nly=1
      nlx_core = nlx - NXR_B_Core(nly) + 1
      nly_core = nly - NYR_B_Core      + 1

!      write(*,*) 'nlx_core =', nlx_core, 'nly_core=', nly_core


      RETURN
      END
                        
      INTEGER FUNCTION  nlz_core(nlz)
!=====================================================================*
! convert reactor coordinates into core coordinates (z)               !
!=====================================================================*
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Input:
      INTEGER nlz
! Output:
!      integer nlz_core

      nlz_core = nlz - NZR_Core_Beg + 1

      RETURN
      END

      SUBROUTINE OUTput_Problem_Description
!=====================================================================*
!      output of the problem description into the file                *
!              "Output/SKETCH.lst"                                    *
!=====================================================================*
!      USE HOMOGENIZATION_XS
      IMPLICIT NONE
      INCLUDE 'sketch.fh'

! Starting output
      OPEN(io_unit,FILE='Output/SKETCH.lst',STATUS='UNKNOWN')

      CALL OUTput_Write_Header(io_unit)

      CALL CNT_Output_IO_Files(io_unit)

      CALL CNT_Output_Code_Options(io_unit)

      CALL CNT_Output_Conv_Criteria(io_unit)

      IF(Problem_Type.NE."Kinetics".OR.(Iter_Solver.EQ."CSA").OR.&
         (Iter_Solver.EQ."CSI") ) THEN
        CALL CNT_Output_Iter_Parameters(io_unit)
      END IF

      IF(Problem_Type.EQ."Kinetics") THEN
         CALL CNT_Output_Time_Step_Control(io_unit)
      END IF

      CALL CNT_Output_Problem(io_unit)

      CALL GMT_Output_Parameter(io_unit)
      

      CALL GMT_Output_Data(io_unit)

      CALL GMT_Output_Core_Parameter(io_unit)
      CALL GMT_Output_Core_Data(io_unit)

      CALL CRD_Output_Parameters(io_unit)
        CALL CRD_Output_Data(io_unit)

        CALL XS_Output_Parameters(io_unit)
        IF(XS_Model == "POLYNOM" ) THEN 
          CALL XS_Output_Data(io_unit)
        END IF



      IF(TH_Model.EQ."Internal") THEN
        WRITE(io_unit, '(A)')&
          " Internal Thermal-Hydraulics Model"
        CALL THM_Output_Parameters(io_unit)
        CALL FRD_Output_Data(io_unit)
        CALL THM_Output_Data(io_unit)
      ELSE IF(TH_Model.EQ."External") THEN
        WRITE(io_unit, '(A)')&
         "  External Thermal-Hydraulics Model (TRAC code)"
        CALL PVM_Output_Mapping_Parameters(io_unit)
        CALL PVM_Output_Mapping_Data(io_unit)
      ELSE IF(TH_Model.EQ."None") THEN
         WRITE(io_unit, '(A)')&
          " No Feedbacks, Pure Neutronics Calculations"
      END IF

      IF(File_Dist.NE."") CALL RDS_Distr_Output(io_unit)

      CLOSE(io_unit)

      RETURN
      END

      SUBROUTINE OUTput_write_restart_steady_state_data(unit)
!**********************************************************************
!      Output of the steady-state DATA into  Restart File             *    
!**********************************************************************
      USE GeoMeTry_Faces
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      INTEGER unit
! Local
      INTEGER k, n, n1, nl, i, nd

      WRITE(unit) (( Flux(n,k), n=1,NG), k = 1, N_TOT)
      WRITE(unit) k_ef
      WRITE(unit) (((trl_xyz(n,k,nd), n = 1, NG),&
            k = 1,N_TOT), nd = 1, NDD)
      WRITE(unit) ((p(k,n1),k=1,NH),n1=1,NZ), p_col(0,0), p_total 
      WRITE(unit) (((fdback(k,n1,i),k=1,NH),n1=1,NZ),&
                         i = 1, N_FEEDBACK)
      WRITE(unit)  ((((D_Nod(n,k,i,nd), n=1, NG),k=1,N_TOT),i=1,2),&
                      nd=1, NDIR)       
      WRITE(unit)  ((MAT_FD(n, nl), n=1, NG), nl = 1, N_FACES)
      WRITE(unit)  ((MAT_Nod(n, nl), n=1, NG), nl = 1, N_FACES)

      RETURN
      END

      SUBROUTINE OUTput_write_restart_burnup_data(unit)
!**********************************************************************
!      Output of the steady-state DATA into  Restart File             *    
!**********************************************************************
      USE GeoMeTry_Faces
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      INTEGER unit
! Local
      INTEGER k, n1, i

      WRITE(unit) ((brn(k,n1),k=1,NH),n1=1,NZ)
      WRITE(unit) (((conc_isotope(k,n1,i),k=1,NH),n1=1,NZ),&
       i=1,N_ISOTOPE)

      RETURN
      END


      SUBROUTINE OUTput_write_restart_kinetics_data(unit, time)
!**********************************************************************
!       OUTput of the neutron kinetics DATA into Restart File         !
!**********************************************************************
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      INTEGER unit
      REAL time
! Local
      INTEGER k, m, i, n1

      WRITE(unit) ((Prec(m,k),m=1,MD), k=1, N_TOT)
      WRITE(unit) (zrods(i),i=1,NN_CRod), Time_Scram,&
                           Flag_Set_Time_Scram
      WRITE(unit) (source(k), k = 1, N_TOT)
      WRITE(unit) time, dt_save
      WRITE(unit) Pow_point, (Prec_point(m), m = 1, MD), react
      WRITE(unit) (((p_dh(m, k, n1), m =1, MH), k = 1, NH), n1=1, NZ)


      RETURN
      END

      SUBROUTINE OUTput_write_restart_adjoint_flux(unit)
!**********************************************************************
!       Output of the adjont neutron flux into  Restart File          *
!**********************************************************************
      IMPLICIT NONE
      INCLUDE 'sketch.fh'
! Input: 
      INTEGER unit
! Output:
!      REAL time
! Local
      INTEGER k, n

      WRITE(unit) (( Flux_a(n,k), n=1,NG), k = 1, N_TOT)

      RETURN
      END

      