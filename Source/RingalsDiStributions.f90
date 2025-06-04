      subroutine RDS_Input_Distrib
!=====================================================================*
! Input Distributions  from the distribution file                     *
!     BWR Ringals Stability Benchmark                                 *
!                         Slava (c) 22.VI.1999                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Local
      INTEGER i, nc, ns
      REAL weighting_factor(N_POLY)
!
      CALL MSC_SSET(N_POLY, 1.0, weighting_factor) 

            call RDS_Input_Fuel_Types
            if(DEBUG) write(*,*) 'Reading Fuel Types from DIST file'

            call RDS_Read_Real_Dist(io_unit, "BURNUP", File_DIST, &
                NZR_Core, NP_Reactor_Core, burnup, Debug )
            if(DEBUG) write(*,*) 'Reading BURNUP from DIST file'
           i = 1
           CALL OUTput_Convert_Ring_Dist_Core_Reactor(burnup, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

           call RDS_Read_Real_Dist(io_unit, "XENON ", File_DIST, &
                NZR_Core, NP_Reactor_Core, xenon, Debug )
           if(DEBUG) write(*,*) 'Reading XENON from DIST file'
           i = 2
           CALL OUTput_Convert_Ring_Dist_Core_Reactor(xenon, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

!            XENON(:,:) = 0.
                        
            call RDS_Read_Real_Dist(io_unit, "VHIST ", File_DIST, &
                NZR_Core, NP_Reactor_Core, vhist, Debug )
            if(DEBUG) write(*,*) 'Reading VHIST from DIST file'
           i = 3
           CALL OUTput_Convert_Ring_Dist_Core_Reactor(vhist, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

            call RDS_Read_Real_Dist(io_unit, "SSHIST", File_DIST, &
                NZR_Core, NP_Reactor_Core, sshist, Debug )
            if(DEBUG) write(*,*) 'Reading SSHIST from DIST file'
 
           i = 4
           CALL OUTput_Convert_Ring_dist_Core_Reactor(sshist, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

!           write(*,'(A)') "maximum value of CR History ="
!           np=npoly( k_dist_rg_3D_max(2,i), k_dist_rg_3D_max(1,i) )
!           ns=k_dist_rg_3D_max(3,i)
!           write(*,*) dist_rg_col(np,ns,i)
!           write(*,'(A, 3I3)') "coordinates =", &
!             ( k_dist_rg_3D_max(nd,i), nd=1,3)
!           nc = Numb_Poly_Core(np)
!           nsc= ns - NZR_Core_Beg + 1
!           write(*,*) 'sshist =', sshist(nsc, nc)
!           write(*,*) 'nc =', nc, 'nsc=', nsc
!           write(*,'(A,/,3(8F9.5))' ) 'sshist(*,648) =', &
!          (sshist(nsc, nc), nsc=1, NZR_Core)
!           pause

!  Computing Conversion History
      DO nc = 1, NP_Reactor_Core
         do ns = 1, NZR_Core
         conv_hist(ns, nc) =  vhist(ns,nc) + sshist(ns, nc)
         end do
      END do
      i = 5
      CALL OUTput_Convert_Ring_Dist_Core_Reactor(conv_hist, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )


!  VOID are temporally, will not  be used after DEBUG
            call RDS_Read_Real_Dist(io_unit, "VOID  ", File_DIST, &
                NZR_Core, NP_Reactor_Core, void, Debug )
            if(DEBUG) write(*,*) 'Reading VOID from DIST file'
           i = 6
           CALL OUTput_Convert_Ring_Dist_Core_Reactor(void, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

            call RDS_Read_Real_Dist(io_unit, "POWER ", File_DIST, &
                NZR_Core, NP_Reactor_Core, pow_rin, Debug )
            if(DEBUG) write(*,*) 'Reading POWER from DIST file'
           i = 7
           CALL OUTput_Convert_Ring_Dist_Core_Reactor(pow_rin, &
                dist_rg_col(0,0,i), dist_rg_mm(-3, i), &
                k_dist_rg_mm(1,-3, i) )

!           write(*,'(A)') "minimum value of POWER ="

!           np=npoly( k_dist_rg_3D_min(2,i), k_dist_rg_3D_min(1,i) )
!           ns=k_dist_rg_3D_min(3,i)
!           write(*,*) dist_rg_col(np,ns,i), "dist_rg_3D_min(i)=",&
!             dist_rg_3D_min(i)
!           write(*,'(A, 3I3)') "coordinates =", &
!             ( k_dist_rg_3D_min(nd,i), nd=1,3)
!           nc = Numb_Poly_Core(np)
!           nsc= ns - NZR_Core_Beg + 1
!           write(*,*) 'POWER =', POW_RIN(nsc, nc)
!           write(*,*) 'nc =', nc, 'nsc=', nsc
!           write(*,'(A,/,3(8F9.5))' ) 'POW_RIN(*,MIN) =', &
!          (POW_RIN(nsc, nc), nsc=1, NZR_Core)
!           pause

! Preparing 2D & 1D data for OUTPUT

      DO i = 1, N_OUT_RG_DIST

      CALL OUTput_Compute_2D_Average(dist_RG_col(0,0,i), NZR_Core_Beg, &
        NZR_Core_End, Index_Core, dist_rg_mm(-3,i), &
        k_dist_RG_mm(1,-3, i))        

      CALL OUTput_Compute_1D_Average(dist_RG_col(0,0,i), NZR_Core_Beg, &
        NZR_Core_End, Index_Core, dist_rg_mm(-3,i), &
        k_dist_RG_mm(1,-3, i), weighting_factor)        

      END DO


      return
      end


      subroutine RDS_Distr_Output(unit)
!=====================================================================*
! Output of the power distribution into "SKETCH.lst"
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'sketch.fh'
     
!     Input:
      INTEGER unit
! LOcal

      character*80 Header_Map
      character*4 val_fmt
      character*6 val_char(0:N_POLY)
      integer ind, ns, i
!      integer nlx_core, nly_core
! external function (in "OUTput.f")
      integer nlz_core

!      CHARACTER*15 dist_RG_units(N_OUT_RG_DIST)

      real dist_RG_scaling_factor(N_OUT_RG_DIST)
      data dist_RG_scaling_factor / N_OUT_RG_DIST*1./

!  Scaling Factor for Burnup = 1./burnup_av
!  Scaling Factor for Xenon = 1./xenon_av
      dist_RG_scaling_factor(1) = &
        1./dist_RG_col(0,0,1)
      dist_RG_scaling_factor(2) = &
        1./dist_RG_col(0,0,2)

!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            RINGHALS DISTRIBUTIONS"
      WRITE(unit,'(A, E12.6, A)') &
     "            Average Value of Burnup =", &
         dist_rg_col(0,0,1), " [MWd/tU]"
      WRITE(unit,'(A, E12.6, A)') &
     "            Average Value of Xenon =", &
         dist_rg_col(0,0,2), " [1/cm^3]"
      CALL OUTput_Write_Separator(unit)



      DO i = 1, N_OUT_RG_DIST
         WRITE(Header_Map, '(5x, A)') &
           NAME_ST_DIST(i)

      CALL OUTput_Distrb_Summary( io_unit, Header_Map, &
        dist_RG_scaling_factor(i),  N_POLY, NZR, &
        dist_rg_col(0,0,i), dist_rg_mm(-3,i), k_dist_rg_mm(1,-3,i) )

      END DO


!      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A)') &
     "            2D PARAMETER DISTRIBUTION"
      CALL OUTput_Write_Separator(unit)


      DO i = 1, N_OUT_RG_DIST
      
      
      WRITE(Header_Map, '(5x, A)') &
           NAME_ST_DIST(i)

      val_fmt = "A6"
      DO ind = 1, N_POLY
          WRITE(val_char(ind), '(F6.3)')  &
         dist_rg_col(ind,0,i)*dist_RG_scaling_factor(i)
      END DO
      
      CALL OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

      END DO

      CALL OUTput_Write_Separator(unit)
      WRITE(unit,'(A, /)') &
     "     1D AXIAL PARAMETER DISTRIBUTIONS "
!      CALL OUTput_Write_Separator(unit)
      WRITE(unit, '(1x, 8(A), /)') "  N  ", &
     "  BURNUP",&
     "   XENON", &
     " VOID HS",&
     " CRD  HS",&
     " CONV HS",&
     "    VOID",&
     "   POWER"

       do ns = NZR_Core_Beg, NZR_Core_End ! NZ_Core_BEG, NZ_Core_End
          WRITE(unit,'(1x, I3,": ", 7F8.4)') &
         nlz_core( ns ) , &
         (dist_RG_col(0, ns, i)*dist_RG_scaling_factor(i), &
            i=1, N_OUT_RG_DIST)
       end do

      CALL OUTput_Write_Separator(unit)

 
      RETURN
      END   



      subroutine RDS_Compute_2d_1D_Ring_Dist(MIN_MAX,&
        Ring_Dist, Ring_Dist_2D,&
        Ring_Dist_1D, RING_Dist_av,&
        Ring_3D_Max, Ring_2D_Max, Ring_1D_max,&
        k_3d_max, k_2d_Max, k_1d_max  )
!=====================================================================*
!  Computing Average 1-d and Average 2d Ringhals distributions        * 
! Last update 1.11.1999
! 25.VI.1999 (c) Slava                                                * 
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: 
      real Ring_Dist(NZR_Core, NP_Reactor_Core) 
      CHARACTER*3 MIN_MAX
!      pow_rin(NZ_Core, NP_Reactor_Core)
! Output:
      real Ring_Dist_2D(NP_Reactor_Core)
      real Ring_Dist_1D(NZR_Core)
      real Ring_Dist_Av
      real Ring_3D_Max, Ring_2D_Max, Ring_1D_max
      integer k_3d_max(2), k_2d_Max, k_1d_max
!      pow_av_rin
!      pow_1d_rin(NZ_Core)
!      pow_2d_rin(NP_Reactor_Core)
! Local variables:
      integer nc, kc

      if(MIN_MAX.EQ."MAX") then
        Ring_3D_Max = 0. 
        Ring_2D_Max = 0. 
        Ring_1D_max = 0.
      else if(MIN_MAX.EQ."MIN") then
        Ring_3D_Max = 1.E+30 
        Ring_2D_Max = 1.E+30 
        Ring_1D_max = 1.E+30 
      end if

      Ring_Dist_Av = 0.
      do nc = 1, NZR_Core
         Ring_Dist_1D(nc) = 0.
      end do
      do kc = 1, NP_Reactor_Core
         Ring_Dist_2D(kc) = 0.
      end do

      do kc = 1, NP_Reactor_Core
         do nc = 1, NZR_Core
         Ring_Dist_Av = Ring_Dist_Av + Ring_Dist(nc, kc)
         Ring_Dist_1D(nc) = Ring_Dist_1D(nc) + Ring_Dist(nc, kc)
         Ring_Dist_2D(kc) = Ring_Dist_2D(kc) + Ring_Dist(nc, kc)

         if(MIN_MAX.EQ."MAX") then
            if(Ring_3D_Max .LT. Ring_Dist(nc, kc) ) then
               Ring_3d_Max = Ring_Dist(nc, kc)
               k_3d_max(1) = kc
               k_3d_max(2) = nc
            end if
         else if (MIN_MAX.EQ."MIN") then
            if(Ring_3D_Max .GT. Ring_Dist(nc, kc) ) then
               Ring_3d_Max = Ring_Dist(nc, kc)
               k_3d_max(1) = kc
               k_3d_max(2) = nc
            end if
         end if

         end do
      end do

      do kc = 1, NP_Reactor_Core

         Ring_Dist_2D(kc) = Ring_Dist_2D(kc)/real(NZR_Core)

         if(MIN_MAX.EQ."MAX") then
            if(Ring_2D_Max .LT. Ring_Dist_2D(kc) ) then
               Ring_2d_Max = Ring_Dist_2D(kc)
               k_2d_max = kc
            end if
         else if (MIN_MAX.EQ."MIN") then
            if(Ring_2D_Max .GT. Ring_Dist_2D(kc) ) then
               Ring_2d_Max = Ring_Dist_2D(kc)
               k_2d_max = kc
            end if
        end if

      end do

      do nc = 1, NZR_Core
         Ring_Dist_1D(nc) = Ring_Dist_1D(nc)/real(NP_Reactor_Core)
         if(MIN_MAX.EQ."MAX") then
            if(Ring_1D_Max .LT. Ring_Dist_1D(nc) ) then
               Ring_1d_Max = Ring_Dist_1D(nc)
               k_1d_max = nc
            end if
         else if (MIN_MAX.EQ."MIN") then
            if(Ring_1D_Max .GT. Ring_Dist_1D(nc) ) then
               Ring_1d_Max = Ring_Dist_1D(nc)
               k_1d_max = nc
            end if
        end if
      end do

      Ring_Dist_Av = Ring_Dist_Av/real(NZR_Core*NP_Reactor_Core)


      return
      end                                      

      subroutine RDS_Compute_2D_1D_Ring_Errors(Dist_2D_Rin,&
        Dist_1D_Rin, Dist_av_Rin, &
        Dist_2D, Dist_1D, Dist_av, &
        Error_2D, Error_2D_max, k_error_2d_max, &
        Error_2D_av, &
        Error_1D, Error_1D_max,&
        k_error_1D_max, error_1D_av, error_av)
!=====================================================================*
!  Computing Errors in Average 1-d and Average 2d Ringhals            *
!                    distributions                                    * 
! 1.11.1999  (c) Slava                                                * 
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input:
      real Dist_2D_Rin(NP_Reactor_Core), Dist_1D_Rin(NZR_Core),&
          Dist_av_Rin 
      real Dist_2D(N_POLY), Dist_1D(NZR), Dist_Av
! Ouput:

      real Error_2D(NP_Reactor_Core), Error_2D_max, &
        Error_2D_av
      integer k_error_2d_max(2) 
      real Error_1D(NZR_Core), Error_1D_max,&
         error_1D_av
      integer k_error_1D_max
      real Error_av
! Local:
      integer kc, nc

      Error_2d_Max = 0.
      Error_2d_Av = 0.

      do kc  = 1, NP_Reactor_Core
         Error_2d(kc) = &
            ( (Dist_2D(numb_reactor_core(kc))/Dist_av)/&
               (Dist_2D_Rin(kc)/Dist_Av_Rin )  - 1.)*100.
         Error_2d_Av = Error_2d_av + abs(Error_2d(kc))
         if(abs(Error_2d(kc)) .GT. abs(Error_2d_max) )  then
             Error_2d_Max =  Error_2d(kc)
             k_error_2d_max(1) = N_Coord(numb_reactor_core(kc), 1) - &
                Nxr_B_Min_Core + 1
             k_error_2d_max(2) = N_Coord(numb_reactor_core(kc), 2) -&
                Nyr_B_Core + 1
         end if 
      end do
      Error_2d_Av = Error_2d_av/NP_Reactor_Core

      Error_1d_Max = 0.
      Error_1d_Av = 0.
      do nc  = 1, NZR_Core
         Error_1d(nc) = &
            ( (Dist_1D(nc + NZR_Core_Beg - 1)/Dist_Av)/&
               (Dist_1D_Rin(nc)/Dist_Av_Rin )  - 1.)*100.
         Error_1d_Av = Error_1d_av + abs( Error_1d(nc) )
         if(abs(Error_1d(nc)).GT. abs(Error_1d_max) )  then
             Error_1d_Max =  Error_1d(nc)
             k_error_1d_max = nc
         end if 
      end do

      Error_1d_av = Error_1D_Av/NZR_Core

      Error_av = (Dist_av - Dist_Av_Rin)*100./Dist_Av_Rin

      return 
      end
           


      subroutine RDS_Input_Fuel_Types
!=====================================================================*
! Input Fuel Rod Types for the Reactor Core from the DIST file        *
!     BWR Ringals Stability Benchmark                                 *
!                         Slava (c) 22.VI.1999                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      integer LEN_HEADER 
      parameter (LEN_HEADER = 6)
      character*6 header, line
      data header /"FUEROD"/
      Logical Error
      integer N_LINE_COMMENTS
      parameter (N_LINE_COMMENTS = 3)
      integer FUEROD(NZR_CORE)

      integer n, k, n1

      open(io_unit, file = FILE_DIST, status = 'old', ERR = 100) 

      call MSC_Search_Header_In_File(io_unit,  header, line,  &
        '(A6)', error)  
      if(Error) go to 200

      do n = 1, N_LINE_COMMENTS
        read(io_unit, *)
      end do

      do k = 1, NP_Reactor_Core
         read(io_unit, '(I6,1x,8(1x,I8.8)/(7x,8(1x,I8.8)) )',ERR=250)&
            n, (FUEROD(n1), n1=1, NZR_CORE)
! Assigning values 
         do n1 = 1, NZR_Core
           CR_Type(n1, k) = FUEROD(n1) / 1000000
           CR_Presence(n1,k)=(FUEROD(n1)-1000000*CR_Type(n1, k))/1000
           Fuel_Type_Core(n1, k)=FUEROD(n1)-1000000*CR_Type(n1, k)-&
               CR_Presence(n1,k)*1000
         end do

      end do

      close(io_unit)  

      go to 300
  100 write(*,*) 'Distribution File ', FILE_DIST, ' does not exist'
      stop
  200 write(*,*) 'Could not find FUEROD in the distribution file ',&
       FILE_DIST
      stop
  250 write(*,*) 'Error during reading FUEROD in distribution file ',&
       FILE_DIST 
      write(*,*) 'The Last Input: Assembly Number ', n
      write(*, '( "FUEROD =", 25I4)' ) FUEROD
      stop

  300 continue

      if(DEBUG) call Output_Fuel_Type_Core

      return
      end     



      subroutine Output_Fuel_Type_Core
!=====================================================================*
! DEBUG Output of the core Fuel Types                                 *
!                         Slava (c) 22.VI.1999                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      integer LEN_HEADER 
      parameter (LEN_HEADER = 6)
      character*6 header
      data header /"FUEROD"/
!      Logical Error
      integer N_LINE_COMMENTS
      parameter (N_LINE_COMMENTS = 3)
      integer FUEROD(NZR_CORE)

      integer n, k, n1

      open(io_unit, file = 'Output_Debug/FILE_DIST.dat', &
         status = 'unknown') 

      write(io_unit,'(A6)') header
      do n = 1, N_LINE_COMMENTS
        write(io_unit, *)
      end do

      do k = 1, NP_Reactor_Core
        do n1 = 1, NZR_Core
        FUEROD(n1) = Fuel_Type_Core(n1, k) + 1000000*CR_Type(n1, k) +&
               CR_Presence(n1,k)*1000
        end do

         write(io_unit,'(I6,1x,8(1x,I8.8)/(7x,8(1x,I8.8)) )')&
            k, (FUEROD(n1), n1=1, NZR_CORE)
! Assigning values 
      end do

      close(io_unit)  

      return
      end

 
      
      subroutine RDS_Read_Real_Dist(io_unit, header, File_DIST, &
          NZR_Core, NP_Reactor_Core, Distrib, Debug )
!=====================================================================*
! reading Real distribution from the DIST file                        *
!     BWR Ringals Stability Benchmark                                 *
!                         Slava (c) 22.VI.1999                        *
!=====================================================================*
      implicit none
! Input:
      integer io_unit
      integer LEN_HEADER 
      parameter (LEN_HEADER = 6)
      character*6 header, line
      character*100 File_DIST
      integer NZR_Core, NP_Reactor_Core
      Logical Debug
! Output:
      real Distrib(NZR_Core, NP_Reactor_Core)
! Local Variables
      Logical Error
      integer N_LINE_COMMENTS
      parameter (N_LINE_COMMENTS = 2)

      real Dist_Scale
      integer n, k, n1

      open(io_unit, file = FILE_DIST, status = 'old', ERR = 100) 

      call MSC_Search_Header_In_File(io_unit,  header, &
        line, '(A6)', error)  
      if(Error) go to 200

      do n = 1, N_LINE_COMMENTS
        read(io_unit, *)
      end do

      read(io_unit, '(50x, 1PE7.1)', ERR = 225) Dist_Scale

      do k = 1, NP_Reactor_Core

         read(io_unit, '(I6,1x,8F9.5/(7X, 8F9.5) )',ERR=250)&
            n, (Distrib(n1, k), n1=1, NZR_CORE)
! Assigning values 
         do n1 = 1, NZR_Core
           Distrib(n1, k) = Distrib(n1, k)*Dist_Scale
         end do

      end do

      close(io_unit)  

      go to 300
  100      write(*,*) 'Distribution File ', FILE_DIST, ' does not exist'
      stop
  200 write(*,*) 'Could not find ', header,' in the distribution file ',&
       FILE_DIST
      stop
  225 write(*,*) 'Error reading scaling factor for ', header, &
        'in the distribution file', FILE_DIST
      stop
  250 write(*,*) 'Error during reading', header,&
       '  in distribution file ',FILE_DIST 
      write(*,*) 'The Last Input: Assembly Number ', n
      write(*, '( "Dist_Input =", I6,1x,8F9.5/(7X, 8F9.5))' ) &
          ( Distrib(n1,k), n1=1, NZR_Core)
      stop

  300 continue

      if(DEBUG) call RDS_Output_Distrib( io_unit, header, &
          NZR_Core, NP_Reactor_Core, Distrib, Dist_Scale )

      return
      end     


      subroutine RDS_Output_Distrib( io_unit, header, &
          NZR_Core, NP_Reactor_Core, Distrib, Dist_Scale )
!=====================================================================*
! DEBUG Output of DISTRIBUTIONS                                       *
!                         Slava (c) 22.VI.1999                        *
!=====================================================================*
      implicit none

      integer LEN_HEADER 
      parameter (LEN_HEADER = 6)
      character*6 header
      integer N_LINE_COMMENTS
      parameter (N_LINE_COMMENTS = 2)
      integer NZR_Core, NP_Reactor_Core, io_unit
      real Distrib(NZR_Core, NP_Reactor_Core), Dist_Scale

      integer n, k, n1

      open(io_unit, file = 'Output_Debug/FILE_DIST.dat', &
         status = 'unknown', access = 'append') 

      write(io_unit,'(A6)') header
      do n = 1, N_LINE_COMMENTS
        write(io_unit, *)
      end do
      write(io_unit, '(50x, 1PE7.1)') Dist_Scale

      do k = 1, NP_Reactor_Core
         write(io_unit,'(I6,1x,8F9.5/(7X, 8F9.5) )')&
            k, (Distrib(n1,k)/Dist_Scale, n1=1, NZR_CORE)
      end do


      close(io_unit)  

      return
      end

