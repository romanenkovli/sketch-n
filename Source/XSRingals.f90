      subroutine XSR_Compute_BWR_CR_XS
!=====================================================================*
! Current Value of Macro Cross-Sections due to Control Rods           *
!   Ringhals BWR   ,      Slava (c) 24.VI.1999                        *
!=====================================================================*
      IMPLICIT NONE

      include 'sketch.fh'

        
! Input: NN_CRod_El, NN_CRod, Mat_Com_rod(NN_CRod, NN_CRod_El), NCHM, 
!        poly_out(N_POLY, NCHM), nrods(NN_CRod, NN_Crod_Bundle),
!        NZ_Core_Beg, NZ, 
!        NZ_Core_End, volume(N_TOT), rod_node(NN_CRod, NZ, NN_CRod_El),
!        xs_ca(NUMBER_XS, N_ROD_COMP), XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) 
! Output: XS_D(NG, N_TOT), XS_SA(NG, N_TOT), 
!        XS_SIK(NG,NG,N_TOT), XS_SF(NG, N_TOT),XS_SF_P(NG,N_TOT) -
!        Corrected Value of the Macro XS Due to Control Rods 
!                                     (XS*Volume, Except XS_D)
!Local Variables
      integer ie, ir, icr, np, nch, k, n1, kt, &
        n1_rod, n1_nrod, kc, nc, ib
      real h_rod, h_nrod, sect_nrod, sect_rod, vol,  str
      real xappa_rod(NG)
!External Function
      real XS_sect_CR

      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 
      Logical Rodded
      parameter( Rodded = .True. )
      integer N_XS_TAB
      parameter (N_XS_TAB = 8)
      real XS_INTER(N_XS_TAB)
!      real Temp_Fuel_Const
!      real k_inf, r_inf


!      write(*,*) 'NZ_CORE_BEG, NZ_COE_END =', NZ_CORE_BEG, &
!               NZ_CORE_END 

      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod

         icr = Mat_Com_Rod(ir,ie)
         do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)

         kc =  Numb_Poly_Core(np)

         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then

           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
              kt = k + (n1-1)*NH
              vol = volume(kt)
              h_rod = rod_node(ir,n1,ie)
              h_nrod = 1. - h_rod
! in the BWR case, Control Rods starts from the bottom
              n1_rod = n1 - 1
              n1_nrod = n1 + 1

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node
!           write(*,*) 'Rodded Node = k, kc, ir, n1, icr =', &
!                k, kc, ir, n1, icr
!           pause

           call XS_Compute_Xappa_Rod(h_rod,e_rod_round_off, &
                 n1, n1_rod, n1_nrod, k, ir, ie, &
                 Xappa_Rod)

           nc = n1 - NZ_Core_Beg + 1


         if(TH_Model.EQ."External") then
         call XSR_Interpolate_XS_Table(&
           Fuel_Type_Core(nc, kc), icr, Rodded,  &
           n_burnup_int(nc, kc), n_hist_void_int(nc,kc),&
           burnup(nc, kc), conv_hist(nc, kc), &
           fdback(k,n1,3), fdback(k,n1,4), xenon(nc,kc), XS_INTER)
         else
!         Temp_Fuel_Const = 794.
         call XSR_Interpolate_XS_Table(&
           Fuel_Type_Core(nc, kc), icr, Rodded,  &
           n_burnup_int(nc, kc), n_hist_void_int(nc,kc),&
           burnup(nc, kc), conv_hist(nc, kc), &
           fdback(k,n1,3), fdback(k,n1,4), xenon(nc,kc) , XS_INTER) ! ZERO XENON
         end if
         if( (h_rod+e_rod_round_off).GE.1 ) then
! if the node is full rodded
           call XSR_Convert_Ringhals_to_SKETCH(XS_INTER, vol,  &
            XS_D(1,kt), XS_SIK(2,1,kt), XS_SA(1,kt), XS_SF(1,kt), &
            XS_SF_P(1,kt))
         else
! if the node is partially rodded

! Diffusion Coefficients
              sect_nrod = 1./(3.*XS_D(1,kt) )
              sect_rod  = 1./(3.*XS_INTER(1) )
              str = XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(1))
              XS_D(1, kt) = 1./(3.*str)

              sect_nrod = 1./(3.*XS_D(2,kt) )
              sect_rod  = 1./(3.*XS_INTER(2) )
              str = XS_Sect_CR(sect_nrod, sect_rod, xappa_rod(2))
              XS_D(2, kt) = 1./(3.*str)

! all other XS are multiplied by the node volume 
! scattering
              sect_nrod = XS_SIK(2, 1, kt)
              sect_rod = XS_INTER(3)*vol
              XS_SIK(2, 1, kt)  =  XS_Sect_CR(sect_nrod, sect_rod,&
              xappa_rod(1) )
     
! removal from the 1st group
              sect_nrod = XS_SA(1, kt)
              sect_rod = ( XS_INTER(4) + XS_INTER(3) )*vol
              XS_SA(1, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(1))
! thermal absorption
              sect_nrod = XS_SA(2, kt)
              sect_rod = XS_INTER(5)*vol
              XS_SA(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
                     xappa_rod(2))

              sect_nrod = XS_SF(1, kt)
              sect_rod =  XS_INTER(6)*vol
              XS_SF(1, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(1))

              sect_nrod = XS_SF(2, kt)
              sect_rod =  XS_INTER(7)*vol
              XS_SF(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(2))

              sect_nrod =  XS_SF_P(1, kt)
              sect_rod = XS_INTER(6)*vol/XS_INTER(8)
              XS_SF_P(1, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(1))

              sect_nrod =  XS_SF_P(2, kt)
              sect_rod = XS_INTER(7)*vol/XS_INTER(8)
              XS_SF_P(2, kt) = XS_Sect_CR(sect_nrod, sect_rod, &
              xappa_rod(2))

           end if ! selection for fully & partially rodded nodes

             end if ! if there is anything in the node
           end do ! NZ
          end if ! k.ie.0
         end do ! NCHM
         end do ! ib 
       end do ! NN_CRod
      end do ! NN_CRod_El

!      if(Debug) call XS_DBG_output

      return
      end

      
      subroutine XSR_Convert_Ringhals_to_SKETCH(XS_INTER, vol,  &
        XS_D, XS_SIK, XS_SA, XS_SF, XS_SF_P)
!=====================================================================*
!      Converting Ringhals XS into SKETCH internal XS set             *
!   (c)  Slava 24.VI.1999                                             *
!=====================================================================*
      implicit none
!Input:
      real XS_INTER(8), vol
!Output:
      real XS_D(2), XS_SIK, XS_SA(2), XS_SF(2), XS_SF_P(2)

!      "D1", "D2", "SR1", "SA1", "SA2", "nuSF1"  "nuSF2", "NY"
!        1    2      3      4      5      6        7      8   
     
      XS_D(1) = XS_INTER(1)
      XS_D(2) = XS_INTER(2)
      XS_SIK = XS_INTER(3)*vol
! Absorption + Scattering
      XS_SA(1) = ( XS_INTER(4) + XS_INTER(3) )*vol
      XS_SA(2) = XS_INTER(5)*vol
      XS_SF(1) = XS_INTER(6)*vol
      XS_SF(2) = XS_INTER(7)*vol
      XS_SF_P(1) = XS_INTER(6)*vol/XS_INTER(8)
      XS_SF_P(2) = XS_INTER(7)*vol/XS_INTER(8)

      return 
      end

      subroutine XSR_Compute_BWR_XS
!=====================================================================*
!  Computing XS for the Reactor Core                                  *
!    Table Interpolation                                              *
!     All cross-Sections without CR                                   * 
!    Ringhals BWR (c) Slava 23.VI.1999                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer n1, nc, nn, kc, k, kt, np, kh
      integer icr
      parameter(icr = 0)
      Logical Rodded
      parameter(Rodded = .False.)
      real vol
      integer N_XS_TAB
      parameter (N_XS_TAB = 8)
      real XS_INTER(N_XS_TAB)
!      real Temp_Fuel_Const


      do n1 = NZ_Core_Beg, NZ_Core_End
         nc = n1 - NZ_Core_Beg + 1
         nn = (n1 - 1)*NH
         do kh = 1, NH_Core
         k = np_core(kh)
         kt = k + nn
         np = np_out(k)
         kc =  Numb_Poly_Core(np)

         vol = volume(kt)

!         write(*,*) 'k, n1, kc, nc =', k, n1, kc, nc
!         write(*,*) 'NZ_Coe_Beg, NZ_Core_End =', NZ_Core_Beg, &
!         NZ_Core_End
!         pause

         if(TH_Model.EQ."External") then
         call XSR_Interpolate_XS_Table(&
           Fuel_Type_Core(nc, kc), icr, Rodded,  &
           n_burnup_int(nc, kc), n_hist_void_int(nc,kc),&
           burnup(nc, kc), conv_hist(nc, kc), &
           fdback(k,n1,3), fdback(k,n1,4), xenon(nc,kc), XS_INTER)
         else
!         Temp_Fuel_Const = 794.
         call XSR_Interpolate_XS_Table(&
           Fuel_Type_Core(nc, kc), icr, Rodded,  &
           n_burnup_int(nc, kc), n_hist_void_int(nc,kc),&
           burnup(nc, kc), conv_hist(nc, kc), &
           fdback(k,n1,3), fdback(k,n1,4), xenon(nc,kc) , XS_INTER) ! ZERO XENON
         end if

! multiplication by volume 
         call XSR_Convert_Ringhals_to_SKETCH(XS_INTER, vol,  &
            XS_D(1,kt), XS_SIK(2,1,kt), XS_SA(1,kt), XS_SF(1,kt), &
            XS_SF_P(1,kt))


         end do
      end do


      return
      end

      subroutine XSR_Initialize_XS_Table
!=====================================================================*
!  Computing two pointers for burnup & conversion history             *
! and values of the Cross-Sections for Reflectors  (once and for ALL) *
!    Ringhals BWR (c) Slava 23.VI.1999                                *
!=====================================================================*
      implicit none
      include 'sketch.fh'

      integer n1, k, nn, kt, np, ns, nl, n, m

      real vol, sa_sik

      do k = 1, NP_Reactor_Core
         do n1 = 1, NZR_Core

         call XSR_Find_Pointers(Fuel_Type_Core(n1, k), &
         burnup(n1, k), CONV_HIST(n1, k), n_burnup_int(n1,k),&
         n_hist_void_int(n1,k) )

         end do
      end do

! Initial values of the Cross-Sections for Reflectors ONLY
      do n1 = 1, NZ 
        nn = (n1-1)*NH
        do k = 1, NH

          kt = k + nn
          vol = volume(kt)
          np = np_out(k)
          ns = ns_out(n1)
          nl = l(np,ns)

! ALL assembly without control rods and without FEEDBACKS

         do n = 1, NG
! multiplication by volume 
         XS_D(n, kt) = d(nl,n)
         XS_SF(n, kt) = sf(nl,n)*vol
         XS_SF_P(n, kt) = sf_p(nl,n)*vol
         XS_SA(n, kt) = sa(nl,n)*vol
         do m = 1, NG
           XS_SIK(m, n, kt) = sik(nl,m,n)*vol
         end do    
         end do ! NG

! Adding scattering
         do n = 1, NG
           sa_sik = 0.
            do m = 1, NG
               sa_sik = sa_sik + XS_SIK(m,n,kt)
            end do
            XS_SA(n, kt) = XS_SA(n, kt) + sa_sik
         end do ! NG

       end do ! NH
      end do ! NZ


      return
      end 


      subroutine XSR_Interpolate_XS_Table(n, icr, Rodded,  &
        n_burnup_int, n_hist_void_int,&
        burnup, history_void, void, fuel_temp, xe, XS_INTER)
!=====================================================================*
!  Interpolating cross-section for in XS Tables                       *
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none 
      include 'cd_file.fh'  
! Input:
      integer n, icr
      Logical Rodded
      real  burnup, history_void, void, fuel_temp, xe
      integer n_burnup_int, n_hist_void_int
! Output:
      real XS_INTER(N_XS)
! Local Values:
      integer ift ! Internal number of fuel types
      integer n_void_int
      logical Comp_Ref_Xenon
      parameter (Comp_Ref_Xenon = .True.)
!      real k_inf, r_inf
!      real Power_Table
!      real XS_SA1_OLD, Delta_SA1

      ift = IND_FTYP(n)

      call XSR_Find_1D_Pointer( NV(ift), &
               void, VENTRY(1, ift),  n_void_int)


      call XSR_Interpolate_XS(ift, n_burnup_int, &
            n_void_int, n_hist_void_int,  burnup, void,&
            history_void, Rodded, icr,  XS_INTER )

!      if(Comp_Ref_Xenon) then
!        Power_Table = Powden(ift)*U_Weight(ift)
!         Power_Table = Powden(ift)*U_Weight(ift)*25.
!      end if

!      write(*,*) 'Power_Table =', Power_Table
!      pause

      call XSR_Compute_XS_XE(N_XS, burnup, void, Xe, &
            Sigma_Xe(1, ift), XS_INTER, Comp_Ref_Xenon)

!      XS_SA1_OLD = XS_INTER(4)
      call XSR_Compute_XS_Doppler( Fuel_temp, ift, Burnup, &
                        Void, XS_INTER )
!      Delta_SA1 = (XS_INTER(4)/XS_SA1_OLD - 1.)*100
!      write(*,*) 'Doppler Effect in SA1, SA1_OLD, SA1, Change (%)'
!      write(*,*)  XS_SA1_OLD, XS_INTER(4), Delta_SA1
!      read(*,*)       

      return
      end

      subroutine XSR_Find_Pointers(n, burnup, history_void, &
         n_burnup_int, n_hist_void_int)
!=====================================================================*
!  Computing Pointers to the Burnup Interval and History Void Interval*
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none 
      include 'cd_file.fh'  
! Input:
      integer n
      real  burnup, history_void
! Output: 
      integer n_burnup_int, n_hist_void_int
! Local Values:
      integer ift

      ift = IND_FTYP(n)

      call XSR_Find_1D_Pointer( NE(ift), &
               burnup, BENTRY( 1, ift),  n_burnup_int)


      call XSR_Find_1D_Pointer( NC(ift), &
                history_void, CENTRY(1, ift), n_hist_void_int)

      return
      end

      subroutine XSR_Input_CD_File(File_CD)
!=====================================================================*
!  Input all data from FILE_CD and scaling the input data             *
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none 
      include 'cd_file.fh'  
      character*100 File_CD                 

      call XSR_Read_XS(File_CD)
      call XSR_Set_Fuel_Types
      call XSR_Read_DXSEC(File_CD)
      call XSR_Read_DOXESH(File_CD)
      if(DEBUG) call XSR_Write_cd_file

      call XSR_Scale_Input

      return
      end      

      subroutine XSR_Read_XS(File_CD)
!=====================================================================*
!  Reading Cross Sections Tables from the File FILE_CD                *
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none
      include 'cd_file.fh'  
! Input:
!      Logical DEBUG
      character*100 FILE_CD    
      character*6 line
! Local Variables
      Logical Error
      integer ios 
      integer n   

      open(io_unit, file = FILE_CD, status ='old', iostat=ios)
 
      if(ios .NE. 0) then
        write(*,*) "XS File ", FILE_CD, " does not exists"
        stop
      end if
  
      n = 0

      do while(.True.) 

        call MSC_Search_Header_In_File(io_unit,  "D1    ", line, &
        '(A6)' , error) 

        if(ERROR) go to 100 
        BACKSPACE (io_unit)

        n = n + 1

        call XSR_Read_fmt_a11(io_unit, N_XS,  NE_M, NV_M, NC_M, &
                   NE(n), NV(n), NC(n), RSCALE(1,n),&
                   FTYPE(n,1), NAMN, &
                   BENTRY(1,n), CENTRY(1,n), VENTRY(1,n), &
                   XS_R(1,1,1,1,n)  )


      end do

  100 close(io_unit) 

      NR_FUEL_TYPE = n 

      if(DEBUG) call Write_Ringhals_Fuel_Types

      return
      end        

      subroutine Write_Ringhals_Fuel_Types
!=====================================================================*
!  DEBUG Output of the Ringhals Fuel Types                            *
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none
      include 'cd_file.fh'  
      integer n, ic, iv

      open(io_unit, file = 'Output_Debug/R_Fuel_Types.dat', &
             status ='unknown')
        write(io_unit,'("Input finished with ", I3 , " Fuel Types")' )&
         NR_FUEL_TYPE

        write(io_unit, '(10(I4, "[",I2, "]" ) )') &
          ( FTYPE( n,1 ), n, n=1, NR_FUEL_TYPE)
!       
        n = 1
        write(io_unit,*) ' 1st Fuel Type, 1st Cross Section '
        write(io_unit,*) 'Void History Entries :', NC(n)
        write(io_unit, '(10E12.5)') &
        (CENTRY(ic, n) , ic = 1, NC(n))

        write(io_unit,*) 'Void Entries :', NV(n)
        write(io_unit, '(10E12.5)') &
        (VENTRY(iv, n) , iv = 1, NV(n))

      close(io_unit)


      return
      end        

      subroutine XSR_Set_Fuel_Types
!=====================================================================*
!    Set array IND_FTYP(m) converting from  Ringhals fuel types into  *
!      the SKETCH internal fuel types                                 *
!    Ringhals BWR (c) Slava 23.IV.1999                                *
!=====================================================================*
      implicit none
      include 'cd_file.fh'      
! Input: FTYPE(1,n,1)
! Output: IND_FTYPE(m)
      integer m, n
 
      do n = 1, NR_FUEL_TYPE
         m = FTYPE(n,1)
         IND_FTYP(m) = n
      end do

      return
      end


      
      
      subroutine XSR_Scale_Input
      implicit none
!=====================================================================*
!  Scaling the macro cross-sections XS_R and differential  control    *
!         rod cross sections  D_FXS                                   *
!   BWR (c) Slava 23.IV.1999                                          *
!=====================================================================*

      include 'cd_file.fh'      
! Local Variables
      integer m, n, ie, ic, iv, icr
      
      do m= 1, N_XS
        do n = 1, NR_FUEL_TYPE
           do ie = 1, NE(n)
              do iv = 1, NV(n)
                 do ic = 1, NC(n)
                     XS_R(m, ie, iv, ic,  n) = RSCALE(m,n)*&
                      XS_R(m, ie, iv, ic,  n)
                 end do
              end do
           end do
        end do
      end do

! Scaling CR differential cross-section

      do icr = 1, N_CR_TYPE
        do m= 1, N_FDCR
           do n = 1, NR_FUEL_TYPE
             do ie = 1, NE( n)
               do iv = 1, NV( n)
                 do ic = 1, NC( n)

                     D_FXS(m, ie, iv, ic, n, icr) =&
                      RSCALE_DFSX(m, n, icr)*&
                      D_FXS(m, ie, iv, ic, n, icr)

                 end do
              end do
           end do

          end do
        end do
      end do


      return
      end 


      subroutine XSR_Compute_K_Inf(XS, r_inf, k_inf )
!=====================================================================!
!   Computation of the infinite multiplication factor                 !
! Stammler & Abbate "Methods of Steady-State Reactor Physics in       !
!                       Nuclear Design"             p.352             !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
      implicit none
!  input:
      real  XS(8)   ! - macro cross section
!      "D1", "D2", "SR1", "SA1", "SA2", "nuSF1"  "nuSF2", "NY"
!        1    2      3      4      5      6        7      8   
!      integer n

!   output
      real k_inf ! infinite multiplication factor
      real r_inf ! spectrum index of the infinite system
! r_inf = S12/SA2
! k_inf = (nuSF1 + r_inf*nuSF2) / (SA1 + S12)

      r_inf =  XS(3) / XS(5)
      k_inf = ( XS(6) + r_inf*XS(7) ) /  ( XS(3) + XS(4) )

!      open(2,file ='XS_Ringhals', status = 'unknown')
!      write(2,*) (XS(n), n=1,8)
!      write(2,*) 'K-Inf, r_inf =', k_inf, r_inf
!      close(2)
!      pause 'XS Done'

      return
      end

      subroutine XSR_Write_cd_file
!=====================================================================!
! Write a copy of the "CD-file.txt" as "Output_Debug/Cd-file_Copy.txt'!
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
      implicit none

      include 'cd_file.fh'      
! Local Variables
      integer n, icr
      

      open(io_unit, file ='Output_Debug/Cd-file_Copy.txt', &
             status ='unknown')

      do n = 1, NR_FUEL_TYPE

        call XSR_write_fmt_a11(io_unit, N_XS,  NE_M, NV_M, NC_M, &
                   NE(n), NV(n), NC(n), RSCALE(1,n),&
                   FTYPE(n,1), NAMN, &
                   BENTRY(1,n), CENTRY(1,n), VENTRY(1,n), &
                   XS_R(1,1,1,1,n)  )

! DXSEC
        do icr = 1, N_CR_TYPE

              call XSR_write_fmt_a11(io_unit, N_FDCR, NE_M, NV_M, NC_M, &
                   NE(n), NV(n), NC(n), rscale_dfsx(1, n, icr),&
                   FTYPE(n,1), NAME_FDCR, &
                   BENTRY(1,n), CENTRY(1,n), VENTRY(1,n), &
                   D_FXS(1, 1, 1, 1, n, icr) )

              call XSR_write_fmt_A12(io_unit, "DXSEC ", &
                icr, FTYPE(n,1), N_DCR, NAME_DCR, D_XS(1, n, icr) )


        end do ! N_CR_TYPE


! End DXSEC

! DOXESH
         call XSR_write_fmt_a13(io_unit, "DOXESH",&
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE(n,1), Powden(n), &
          Temp_Fuel_Nom(n), Temp_Fuel(1,n), Doppl_Coeff(1,n),&
          Sigma_Xe(1,n), SSHP(1,n)  )

! END DOSESH

      end do ! N

      close(io_unit)

      return
      end        


      Subroutine XSR_Find_1D_Pointer( NP, x, End_Int, int_pointer)
!=====================================================================!
!      Finding Pointers  in a one-dimensional table                   !
!   (c) Slava 27.V.1999                                               !
!=====================================================================!
      implicit none
! Input : burnup, history_void,
      integer NP ! number of points (number of intervals = NP - 1)
      real End_Int(NP)
      real x ! Value
!     real BENTRY(ie,  n) - burnup entries
!     real CENTRY(ic,  n) - history void entries
! Output: 
      integer int_pointer
! Local Variables
      Logical Inside_Int

      int_pointer = 0
      Inside_Int  =  .False.

      do while ( .NOT. (Inside_int) )
  
         int_pointer = int_pointer + 1
         if( (int_pointer + 1) .EQ. NP ) then
           Inside_Int = .TRUE. ! Extrapolation
         else
           Inside_Int  =  &
           (x .LE. End_Int(int_pointer + 1) )
         end if         

      end do

      return
      end       


      real function XSR_Interpolate_3D( xsi, ysi, zsi, f111, f112,&
       f121, f122, f211, f212, f221, f222)
!=====================================================================*
!      3D Tensor Product Linear Interpolation                         *
!   C. Ueberhuber "Numerical Computations, Methods, Software and      *
!                          Analysis", Springer, 1997                  *
!    Sect. 9.9.1 page 431                                             *
!  (c) Slava 27.V.1999                                                *
!=====================================================================*
      implicit none
! Input:
      real xsi, ysi, zsi
      real f111, f112,&
       f121, f122, f211, f212, f221, f222

      real XSR_Interpolate_1D
      external XSR_Interpolate_1D

      
      XSR_Interpolate_3D = (1.- xsi)*( &
       (1. - ysi)*XSR_Interpolate_1D( zsi, f111, f112) + &
         ysi*XSR_Interpolate_1D( zsi, f121, f122) ) +&
        xsi*(&
       (1. - ysi)*XSR_Interpolate_1D( zsi, f211, f212) + &
             ysi*XSR_Interpolate_1D( zsi, f221, f222) )

      return
      end



      real function XSR_Interpolate_1D( xsi, y1, y2)
!=====================================================================*
!        Lagrange one-dimensional Linear Intrerpolation               *
!  (c) Slava 27.V.1999                                                *
!=====================================================================*
      implicit none
! Input:
      real y1, y2, xsi
! Output:
      XSR_Interpolate_1D = (1. - xsi)*y1 + xsi*y2

      return
      end

      

      subroutine XSR_Read_DOXESH(File_CD)
!=====================================================================*
! read Doppler, Xenon and SSHIS coefficients                          *
!             Appendix 1.3 of the benchmark specification             *
!  (c) Slava 27.V.1999                                                *
!=====================================================================*
      implicit none
      include 'cd_file.fh'      
! Local Variables
      character*100 File_CD
      Logical Error
      integer n
      character*6 header, line
      data header / "DOXESH"/

      open(io_unit, file = File_CD, status ='old')

      n = 0
  
      do while(.True.)

        call MSC_Search_Header_In_File(io_unit,  header, line, &
        '(A6)' , error) 
        if(ERROR) go to 100

        n = n + 1

        BACKSPACE( io_unit )            

        call XSR_Read_fmt_a13(io_unit, header,&
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE(n,1), Powden(n), &
          Temp_Fuel_Nom(n), Temp_Fuel(1,n), Doppl_Coeff(1,n),&
          Sigma_Xe(1,n), SSHP(1,n)  )

      end do


  100 continue

      write(*,*) 'DOXESH, last input for Ringhals fuel type ', &
          FTYPE(n,1)
!      pause

      return
      end        




      subroutine XSR_Read_DXSEC(FILE_CD)
!=====================================================================*
! read differential C.R. cross section representation                 *
!             Appendix 1.2 of the benchmark specification             *
!  (c) Slava 27.V.1999                                                *
!=====================================================================*
      implicit none

      include 'cd_file.fh'      
! Local Variables
      character*100 FILE_CD
      Logical Error
      integer n, icr
      character*6 header, line 

      open(io_unit, file = FILE_CD, status ='old')

      n = 0

       do while(.True.)

        call MSC_Search_Header_In_File(io_unit,  "DSA2   ", line, &
        '(A6)' , error) 
        if(ERROR) go to 100
        BACKSPACE( io_unit)          
        n = n + 1

! reading DSA2 & DNSF2
        do icr = 1, N_CR_TYPE

              call XSR_Read_fmt_a11(io_unit, N_FDCR, NE_M, NV_M, NC_M, &
                   NE(n), NV(n), NC(n), rscale_dfsx(1, n, icr),&
                   FTYPE(n,1), NAME_FDCR, &
                   BENTRY(1,n), CENTRY(1,n), VENTRY(1,n), &
                   D_FXS(1, 1, 1, 1, n, icr) )

              call XSR_Read_fmt_A12(io_unit, header, &
                icr, FTYPE(n,1), N_DCR, NAME_DCR, D_XS(1, n, icr) )


        end do ! N_CR_TYPE


      end do

  100 continue

      close(io_unit) 

      return
      end        

          
      subroutine XSR_Interpolate_XS(ift, ix, iy, iz,  &
       burnup, void, history_void, Rodded, icr,  XS_INTER )
!=====================================================================*
! Computing cross sections by tensor product linear interpolation     *
!             Appendix 1.2 of the benchmark specification             *
!  (c) Slava 27.V.1999                                                *
!=====================================================================*

      implicit none 
      include 'cd_file.fh'      
! Input: 
      integer ift, ix, iy, iz, icr
      real  burnup,  void, history_void 
      Logical Rodded
!
! Used functions:
      real  XSR_Interpolate_3D !( xsi, ysi, zsi, f111, f112,&
!       f121, f122, f211, f212, f221, f222)
      external XSR_Interpolate_3D
! Used from the module   
!   XS_R(m, ie, iv, ic,  n) - cross sections

! Output:
      real XS_INTER(N_XS)
! Local Variables:
      real xsi, ysi, zsi
!     computing xsi, ysi, zsi
      integer i_xs
      real f111, f112, f121, f122, f211, f212, f221, f222
      real D_FXS_INT(N_FDCR)
! 
      xsi =  (burnup - BENTRY(ix , ift)) / &
        ( BENTRY(ix + 1 , ift) - BENTRY(ix , ift) )

      ysi =  (void - VENTRY(iy , ift)) / &
        ( VENTRY(iy + 1 , ift) - VENTRY(iy , ift) )

      zsi =  (history_void - CENTRY(iz , ift)) / &
        ( CENTRY(iz + 1 , ift) - CENTRY(iz , ift) )

      do i_xs = 1, N_XS

        f111 = XS_R(i_xs, ix, iy, iz,  ift)

        f112 = XS_R(i_xs, ix, iy, iz + 1,  ift)
        f121 = XS_R(i_xs, ix, iy + 1, iz,  ift)
        f211 = XS_R(i_xs, ix + 1, iy, iz,  ift)

        f122 = XS_R(i_xs, ix, iy + 1, iz + 1,  ift)
        f212 = XS_R(i_xs, ix + 1, iy, iz + 1,  ift)
        f221 = XS_R(i_xs, ix + 1, iy + 1, iz,  ift)

        f222 = XS_R(i_xs, ix + 1, iy + 1, iz + 1,  ift)

        XS_INTER(i_xs) = XSR_Interpolate_3D&
           ( xsi, ysi, zsi, f111, f112, &
                      f121, f122, f211, f212, f221, f222)

      end do

!  Treatment of the Rodded Nodes
      if(Rodded) then 

       do i_xs = 1,  N_FDCR
     
        f111 = D_FXS(i_xs, ix, iy, iz,  ift, icr)

        f112 = D_FXS(i_xs, ix, iy, iz + 1,  ift, icr)
        f121 = D_FXS(i_xs, ix, iy + 1, iz,  ift, icr)
        f211 = D_FXS(i_xs, ix + 1, iy, iz,  ift, icr)

        f122 = D_FXS(i_xs, ix, iy + 1, iz + 1,  ift, icr)
        f212 = D_FXS(i_xs, ix + 1, iy, iz + 1,  ift, icr)
        f221 = D_FXS(i_xs, ix + 1, iy + 1, iz,  ift, icr)

        f222 = D_FXS(i_xs, ix + 1, iy + 1, iz + 1,  ift, icr)

        D_FXS_INT(i_xs) = XSR_Interpolate_3D( &
       xsi, ysi, zsi, f111, f112, f121, f122, f211, f212, f221, f222)

       end do

       XS_INTER(1) = XS_INTER(1) + D_XS(1, ift, icr)
       XS_INTER(2) = XS_INTER(2) + D_XS(2, ift, icr)
       XS_INTER(3) = XS_INTER(3) + D_XS(3, ift, icr)
       XS_INTER(4) = XS_INTER(4) + D_XS(4, ift, icr)
      
       XS_INTER(5) = XS_INTER(5) + D_FXS_INT(1)
      
       XS_INTER(6) = XS_INTER(6) + D_XS(5, ift, icr)

       XS_INTER(7) = XS_INTER(7) + D_FXS_INT(2)

       XS_INTER(8) = XS_INTER(8) + D_XS(6, ift, icr)

      end if
          
      return
      end


      subroutine XSR_write_fmt_a11(io_unit, N_XS, NE_M, NV_M, NC_M,&
       NE, NV, NC, RSCALE,&
       FTYPE, NAMN, BENTRY, CENTRY, VENTRY, XS)
      implicit none 
!=====================================================================!
!  Writing XS data                                                    !
!     format described in Appendix  1.1                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
! Input :
      integer io_unit, N_XS, NE_M, NV_M, NC_M
! Output into the  file:
      integer  NE, NV, NC, FTYPE
      character*6 NAMN(N_XS)
      real RSCALE(N_XS)
      real BENTRY(NE_M), CENTRY(NC_M), VENTRY(NV_M)
      real XS(N_XS, NE_M, NV_M, NC_M)

! Local Variables:
      integer m, ie, iv, ic, i 

      do m = 1, N_XS
         write(io_unit, 10 ) NAMN(m),&
          NE, NV, NC, RSCALE(m), FTYPE, (0, i=1,17)

         write(io_unit, 20) (BENTRY(ie), ie = 1, NE )

           do ic = 1, NC
              write(io_unit, 30) CENTRY(ic)
              do iv = 1, NV
                write(io_unit, 40 )  VENTRY(iv), &
                   (XS(m, ie, iv, ic), ie = 1, NE )
              end do ! NV
           end do ! NC
      end do


   10 FORMAT(A6, 1x, I3, I3, I3, E9.1, 1x, 18I3)
   20 FORMAT(6x,8F9.1)
   30 FORMAT(9x, F5.3)
   40 FORMAT(F5.3, 1X, 8F9.6/(6x,8F9.6) ) 


      return
      end

      subroutine XSR_write_fmt_A12(io_unit, header,&
           icr, FTYPE, N_DCR, NAME_DCR, D_XS )
      implicit none
!=====================================================================!
!  Writing DXSEC data                                                 !
!     format described in Appendix  1.2                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!

      integer io_unit, icr, FTYPE, N_DCR
      character*6 NAME_DCR(N_DCR), header
      real D_XS(N_DCR)
! Local Variables
      integer i, m


        write(io_unit,10) header,  icr,  &
                              FTYPE, (0,i=1,17)
          do m = 1, N_DCR          
             write(io_unit,20) NAME_DCR(m), D_XS(m)
          end do

         write(io_unit, *) 


   10 FORMAT(A6, I1, 19x, 18I3)
   20 FORMAT(A6, F9.6)

      return
      end


      subroutine XSR_write_fmt_A13(io_unit, header,&
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE, Powden, &
          Temp_Fuel_Nom, Temp_Fuel, Doppl_Coeff,&
          Sigma_Xe, SSHP  )
!=====================================================================!
!  Writing DOXESH data                                                ! 
!     format described in Appendix  1.3                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
       
      integer io_unit, &
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE
      character*6 header
      real Powden, Temp_Fuel_Nom, Temp_Fuel(N_FT), &
          Doppl_Coeff(NC_Doppler), Sigma_Xe(NC_XE), SSHP(NC_SSHP)
! Local Variables:
      character*6 Sub_Titl(6)      
      data Sub_Titl / "FTYP  ", "POWDEN", "TEMP  ",&
        "DOPPLE", "SIGAXE",        "SSHPOL" / 
      integer i

      write(io_unit, 10) header
      write(io_unit,  20) Sub_Titl(1), &
             FTYPE, (0,i=1,19)
      write(io_unit, 30 ) Sub_Titl(2), Powden
      write(io_unit, 40 ) Sub_Titl(3), Temp_Fuel_Nom, &
                (Temp_Fuel(i) , i = 1, N_FT)

      write(io_unit, 50 )  Sub_Titl(4),&
                 (Doppl_Coeff(i), i=1, NC_Doppler) 

      write(io_unit, 50)  Sub_Titl(5),&
                  (Sigma_Xe(i), i=1, NC_XE)

      write(io_unit, 50) Sub_Titl(6),&
                  (SSHP(i),i=1, NC_SSHP)

   10 FORMAT( A6 )
   20 FORMAT( A6, 20I3 )   
   30 FORMAT( A6, E12.4 )
   40 FORMAT( A6, 7f8.1 )
   50 FORMAT( A6, 5E12.4 ) 

      return
      end



      subroutine XSR_Read_fmt_A13(io_unit, header,&
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE, Powden, &
          Temp_Fuel_Nom, Temp_Fuel, Doppl_Coeff,&
          Sigma_Xe, SSHP  )
!=====================================================================!
!  Reading DOXESH data from "CD-file.txt",                            !
!     format described in Appendix  1.3                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
       
      integer io_unit, &
          N_FT, NC_Doppler, NC_XE, NC_SSHP, FTYPE
      character*6 header
      real Powden, Temp_Fuel_Nom, Temp_Fuel(N_FT), &
          Doppl_Coeff(NC_Doppler), Sigma_Xe(NC_XE), SSHP(NC_SSHP)
! Local Variables:
      character*6 Sub_Titl(6)      
!      data Sub_Titl / "FTYP  ", "POWDEN", "TEMP  ",&
!        "DOPPLE", "SIGAXE",        "SSHPOL" / 
      integer i

      read(io_unit, 10) header
      read(io_unit,  20) Sub_Titl(1), &
             FTYPE
      read(io_unit, 30 ) Sub_Titl(2), Powden
      read(io_unit, 40 ) Sub_Titl(3), Temp_Fuel_Nom, &
                (Temp_Fuel(i) , i = 1, N_FT)

      read(io_unit, 50 )  Sub_Titl(4),&
                 (Doppl_Coeff(i), i=1, NC_Doppler) 

      read(io_unit, 50)  Sub_Titl(5),&
                  (Sigma_Xe(i), i=1, NC_XE)

      read(io_unit, 50) Sub_Titl(6),&
                  (SSHP(i),i=1, NC_SSHP)

   10 FORMAT( A6 )
   20 FORMAT( A6, 20I3 )   
   30 FORMAT( A6, E12.4 )
   40 FORMAT( A6, 7f8.1 )
   50 FORMAT( A6, 5E12.4 ) 

      return
      end


       

      subroutine XSR_Read_fmt_A12(io_unit, header,&
           icr, FTYPE, N_DCR, NAME_DCR, D_XS )
      implicit none
!=====================================================================!
!  Reading DXSEC data from "CD-file.txt",                             !
!     format described in Appendix  1.2                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!

      integer io_unit, icr, FTYPE, N_DCR
      character*6 NAME_DCR(N_DCR), header
      real D_XS(N_DCR)
! Local Variables
      integer  m


       read(io_unit,10) header,  icr,  &
                              FTYPE
          do m = 1, N_DCR          
             read(io_unit,20) NAME_DCR(m), D_XS(m)
          end do

         read(io_unit, *) 


   10 FORMAT(A6, I1, 19x, 18I3)
   20 FORMAT(A6, F9.6)

      return
      end



      subroutine XSR_Read_fmt_a11(io_unit, N_XS, NE_M, NV_M, NC_M,&
       NE, NV, NC, RSCALE,&
       FTYPE, NAMN, BENTRY, CENTRY, VENTRY, XS)
!=====================================================================!
!  Reading XS data from "CD-file.txt",                                !
!     format described in Appendix A1.1                               !
!     of the Ringhals Benchmark Specification                         !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
      implicit none 
! Input :
      integer io_unit, N_XS, NE_M, NV_M, NC_M
! Output into the  file:
      integer  NE, NV, NC, FTYPE
      character*6 NAMN(N_XS)
      real RSCALE(N_XS)
      real BENTRY(NE_M), CENTRY(NC_M), VENTRY(NV_M)
      real XS(N_XS, NE_M, NV_M, NC_M)

! Local Variables:
      integer m, ie, iv, ic

      do m = 1, N_XS
         read(io_unit, 10 ) NAMN(m),&
          NE, NV, NC, RSCALE(m), FTYPE

         read(io_unit, 20) (BENTRY(ie), ie = 1, NE )

           do ic = 1, NC
              read(io_unit, 30) CENTRY(ic)
              do iv = 1, NV
                read(io_unit, 40 )  VENTRY(iv), &
                   (XS(m, ie, iv, ic), ie = 1, NE )
              end do ! NV
           end do ! NC
      end do


   10 FORMAT(A6, 1x, I3, I3, I3, E9.1, 1x, 18I3)
   20 FORMAT(6x,8F9.1)
   30 FORMAT(9x, F5.3)
   40 FORMAT(F5.3, 1X, 8F9.6/(6x,8F9.6) ) 


      return
      end


      subroutine XSR_Compute_XS_Doppler&
                        (Fuel_temp, ift, Burnup, Void, XS)
!=====================================================================!
!   Computation of the corrected thermal group absorption             !
!       cross-section due to Doppler feedback from thge condition to  !
!       preserve infinite multiplication factor                       !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
      implicit none
      include 'cd_file.fh'      
! Input: Macro Cross Sections
      integer ift ! Fuel Type
      real Burnup, Void
      real  XS(N_XS)   ! - macro cross section
!      "D1", "D2", "SR1", "SA1", "SA2", "nuSF1"  "nuSF2", "NY"
!        1    2      3      4      5      6        7      8   
      real Fuel_Temp
! Output:
!     XS(8) - corrected for Doppler Feedback
! Local Variables:
      real k_inf, r_inf

      call XSR_Compute_K_Inf(  XS, r_inf, k_inf )

      k_inf = k_inf + &
       ( sqrt(Fuel_Temp) - sqrt(Temp_Fuel_Nom(ift)) )*&
       (1. + Burnup*(Doppl_Coeff(1, ift) + Burnup*&
       (Doppl_Coeff(2, ift) + Burnup*Doppl_Coeff(3, ift) )))*&
       ( Doppl_Coeff(4, ift) + void*Doppl_Coeff(5, ift) )

!      write(*,*) Fuel_Temp, Temp_Fuel_Nom(ift)
!      read(*,*) 
! change of the thermal group absorption cross section
!      XS(5) = XS(3)*XS(7) / ( k_inf*(XS(3) + XS(4)) - XS(6) )
! change of the fast group absorption cross section
       XS(4) = ( XS(6) + XS(3)*(XS(7)/XS(5) - k_inf) )/k_inf 

      return
      end     



      subroutine XSR_Compute_XS_XE&
         (N_XS, Burnup, Void, Xe, Const_Xe, XS , Comp_Ref_Xenon )
!          Power_Table)
!=====================================================================!
!      Compute the thermal group absorption cross-section             !
!          taking into account the Xenon poisoning;                   !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
      implicit none
!      "D1", "D2", "SR1", "SA1", "SA2", "nuSF1"  "nuSF2", "NY"
!        1    2      3      4      5      6        7      8   
!      data  NAME_FDCR /'DSA2  ','DNSF2 '/
!      data NAME_DCR / "DD1   ", "DD2   ", "DSIGR1 ", "DSA1  ",&
!       "DNSF1 ", "DNY   " /
! Input:
      integer N_XS
      real XS(N_XS), Burnup, Void, Const_XE(5), XE !, Power_Table
      logical Comp_Ref_Xenon
! Added the parameters to compute the Reference Xenon concentration:
      real Vcore ! Core Volume [cm^3]
      parameter (Vcore = 5.564E+07)
      real Power_Nom   ! Wt        
      parameter (Power_Nom = 2.270E+09)
      real Energy_Conv_Factor
      parameter (Energy_Conv_Factor = 3.237E-11)
      real yield_iodine
      parameter (yield_iodine = 0.064) !!! CHECK
      real yield_xenon 
      parameter (yield_xenon = 0.003)
      real const_decay_xenon 
      real V_Fuel_Bundle
      parameter (V_Fuel_Bundle = 85863.83)
      parameter (const_decay_xenon = 2.095E-05)
      real Sigma_Fission, Flux_Ref, Xenon_Ref
      real Sigma_Micro_Xenon
! Output:
!      real XS(N_XS) ! XS corrected for Xe
! functions:
      real XSR_Get_Micro_XE_XS
      external XSR_Get_Micro_XE_XS


      Sigma_Micro_Xenon = XSR_Get_Micro_XE_XS&
                             ( Const_Xe , Burnup, Void )

!      write(*,*) 'Sigma_Micro_Xenon = ', Sigma_Micro_Xenon
!      write(*,*) 'Xe concentration = ', Xe
!      pause
!      write(*,*) 'Power Table =', power_table*648
!      pause

      if(Comp_Ref_Xenon) then

        Sigma_Fission = XS(7)/XS(8) + ( XS(6)/XS(8) )*( XS(5)/XS(3) )

     
       Flux_Ref = Power_Nom/(Vcore*Energy_Conv_Factor*Sigma_Fission)
!        Flux_Ref = Power_Table/&
!         (V_Fuel_Bundle*Energy_Conv_Factor*Sigma_Fission)

        Xenon_Ref=(yield_iodine + yield_xenon)*Sigma_Fission*Flux_Ref/&
        (const_decay_xenon + Sigma_Micro_Xenon*Flux_Ref)

      else

        Xenon_Ref = 0.

      end if

!      write(*,*) 'Sigma_Micro_Xenon = ', Sigma_Micro_Xenon
!      write(*,*) 'Xe concentration = ', Xe
!      write(*,*) 'Reference Xe concentration = ', Xenon_Ref

!      pause

         XS(5) = XS(5) + &
              Sigma_Micro_Xenon*(Xe - Xenon_Ref)

      return
      end

      real function XSR_Get_Micro_XE_XS(C_XE, Burnup, Void)
!=====================================================================!
! real function to compute Xenon microscopic cross-sections           !
!   Ringhals-1 Problem    Slava (c) 21.V.1999                         !
!=====================================================================!
! real function to compute Xenon microscopic cross-sections
!   Input!
      real C_XE(5)
      real burnup, void
     
      XSR_Get_Micro_XE_XS = ( C_Xe(1) + C_XE(2)*Burnup)*(C_Xe(3) + &
         Void*( C_XE(4) + Void*C_Xe(5) ) )

      return 
      end  


      block data XSR_init_names_Diff_CR_XS
      implicit none
      include 'cd_file.fh'      
      
      data NAME_FDCR / &
       "DSA2  ",&
       "DNSF2 " /

      data NAME_DCR /&
       "DD1   ", &
       "DD2   ", &
       "DSIGR1",&
       "DSA1  ",&
       "DNSF1 ", &
       "DNY   " /
     
      end 
           

      block data XSR_init_U_Weight
      implicit none
      include 'cd_file.fh'      

! weight of uranium (grams)

      data U_Weight / &
      7150., 7283., 7238., 7239., 7234., 6993., 7281., 7285.,&
      7127., 7113., 7118., 7112., 7098., 7103., 7275., 7284.,&
      7256., 7245., 7240., 7257., 7252., 7372., 7374., 7366.,&
      7365., 7352., 7365., 7365., 7365., 7353., 7347.,&
      7390., 7370., 7376., 7241., 7216., 7241., 7370.,&
      7127., 7113., 7118., 7112., 7098., 7103. /                     

      end


      
