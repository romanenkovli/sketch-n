      MODULE HOMOGENIZATION_XS
      IMPLICIT NONE

      LOGICAL FLAG_XS_HOPMOGENIZATION

      INTEGER :: N_HOM_ZONE, N_HOM_NZR, N_COND_NG

      INTEGER, DIMENSION(:), ALLOCATABLE ::&
       ind_cond_ng

      INTEGER, DIMENSION(:), ALLOCATABLE ::&
       ind_hom_zone, ind_hom_nzr

      REAL, DIMENSION(:,:,:), ALLOCATABLE ::&
       XS_HOM_D, XS_HOM_D_STR, XS_HOM_SA, &
       XS_HOM_SF, XS_HOM_SF_P, XS_HOM_FLUX,&
       XS_HOM_DB2 

      REAL, DIMENSION(:,:,:,:), ALLOCATABLE ::&
           XS_HOM_SIK
      
      REAL,DIMENSION(:,:), ALLOCATABLE ::  &
       XS_HOM_K_INF, XS_HOM_POW

      REAL,DIMENSION(:), ALLOCATABLE ::  &
       XS_HOM_XP

      REAL ::  vol_hom_core, pow_hom_total
 

      CONTAINS

      SUBROUTINE HOM_XS_Input  
      INCLUDE 'sketch.fh'
! Local Variables

      integer   n,  np, ios

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
!      character*200 Message
      logical error_find

      FLAG_XS_HOPMOGENIZATION = .False.
       
!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'

!     
      CALL HOM_XS_Allocate_Index   

      OPEN(io_unit,FILE=file_input,STATUS='old',IOSTAT=ios)

!  HOM_IND_ZONE
      REWIND(io_unit)
      CALL MSC_Search_Header_In_File(io_unit, &
        "HOM_IND_ZONE", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
! No need top prepare homogenious XS exit
       GO TO 100
      else
      
      FLAG_XS_HOPMOGENIZATION = .TRUE.  

      read(io_unit,fmt=*,iostat=ios) &
       (ind_hom_zone(np),np=1, N_POLY)
      call Iostat_Error_Check&
      (ios,"Error in reading indexes of the homogenization zones,"&
      //"under identifier HOM_IND_ZONE from the FILE_INPUT")
      end if

!  HOM_IND_NZR
      REWIND(io_unit)
      CALL MSC_Search_Header_In_File(io_unit, &
        "HOM_IND_NZR", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',&
           " identifier HOM_IND_NZR is not found " //&
           "axial zone structure is preserved")
           DO n = 1, NZR
              ind_hom_nzr(n) = n
           END DO
      else
      read(io_unit,fmt=*,iostat=ios) ind_hom_nzr(1:NZR)
      call Iostat_Error_Check&
     (ios,"Error in reading indexes of the axial homogenization "&
      //"zones, under identifier HOM_IND_NZR from the FILE_INPUT")
      end if


!  HOM_IND_NG
      REWIND(io_unit)
      CALL MSC_Search_Header_In_File(io_unit, &
        "HOM_IND_NG", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',&
           " identifier HOM_IND_NG is not found " //&
           "neutron group structure is preserved")
           DO n = 1, NG
              ind_cond_ng(n) = n
           END DO
      else
      read(io_unit,fmt=*,iostat=ios) ind_cond_ng(1:NG)
      call Iostat_Error_Check&
     (ios,"Error in reading indexes of the condensed neutron groups,"&
      //"under identifier HOM_IND_NG from the FILE_INPUT")
      end if

! defining number of homogenization zones and neutron energy groups
      N_HOM_ZONE=MAXVAL( ind_hom_zone(1:N_POLY) )
      N_HOM_NZR =MAXVAL( ind_hom_nzr(1:NZR) )
      N_COND_NG =MAXVAL( ind_cond_ng(1:NG) )

  100 CONTINUE

      CLOSE( io_unit )

      END SUBROUTINE HOM_XS_Input  

      SUBROUTINE HOM_XS_Compute
      INCLUDE 'sketch.fh'

      INTEGER :: n1, nn, k, kt, np, ns, np_hom, ns_hom, n, m,&
                 n_hom, m_hom
      REAL    :: flux_volume_hom(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
                volume_hom(N_HOM_ZONE, N_HOM_NZR),&
                a1, vol, sa_sik,&
                A(N_COND_NG,N_COND_NG), B(N_COND_NG),&
                 pow_hom_norm, fission_source
      REAL conv_Wt_to_MWt
      PARAMETER (conv_Wt_to_MWt = 1.E-06)


        CALL HOM_XS_Allocate_XS

        volume_hom(:,:)          = 0.

        flux_volume_hom(:,:,:)   = 0.
        XS_HOM_D       (:,:,:)   = 0.
        XS_HOM_D_STR   (:,:,:)   = 0. 
        XS_HOM_SA      (:,:,:)   = 0.
        XS_HOM_SF      (:,:,:)   = 0.
        XS_HOM_SF_P    (:,:,:)   = 0.
        XS_HOM_FLUX    (:,:,:)   = 0.
        XS_HOM_SIK     (:,:,:,:) = 0.
        XS_HOM_K_INF   (:,:) = 0.

! Condensed prompt neutron spectrum
        XS_HOM_XP(:) = 0.
        DO n = 1, NG 
            n_hom = ind_cond_ng(n)
            XS_HOM_XP(n_hom)=XS_HOM_XP(n_hom) + xp(n)
        END DO


      DO n1 = 1, NZ 
        nn = (n1-1)*NH
        DO k = 1, NH

          kt = k + nn
          vol = volume(kt)
          np = np_out(k)
          ns = ns_out(n1)

          np_hom = ind_hom_zone(np)
          ns_hom = ind_hom_nzr(ns) 

          volume_hom( np_hom,ns_hom )=&
             volume_hom( np_hom,ns_hom )+vol

          DO n = 1, NG 
            n_hom = ind_cond_ng(n)
            a1 = Flux(n,kt)*vol
! flux - volume   
           flux_volume_hom(n_hom,np_hom,ns_hom)=&
             flux_volume_hom(n_hom,np_hom,ns_hom)+a1

! 1st variant collapsing D
           XS_HOM_D(n_hom, np_hom, ns_hom)= &
      XS_HOM_D(n_hom, np_hom, ns_hom)+XS_D(n,kt)*a1

! 2nd variant collapsing 1/3str
           XS_HOM_D_STR(n_hom, np_hom, ns_hom)= &
      XS_HOM_D_STR(n_hom, np_hom, ns_hom)+a1*(1/XS_D(n,kt))

! nuSF ! volume is included already
           XS_HOM_SF(n_hom, np_hom, ns_hom)= &
      XS_HOM_SF(n_hom, np_hom, ns_hom)+XS_SF(n,kt)*Flux(n,kt)

! SF
           XS_HOM_SF_P(n_hom, np_hom, ns_hom)= &
      XS_HOM_SF_P(n_hom, np_hom, ns_hom)+XS_SF_P(n,kt)*Flux(n,kt)

! SA  ( substracting scattering, which was included ) 
          sa_sik = 0.
          DO m = 1, NG
             sa_sik = sa_sik + XS_SIK(m,n,kt)
           end do
           XS_HOM_SA(n_hom, np_hom, ns_hom)= &
            XS_HOM_SA(n_hom, np_hom, ns_hom)+&
             (XS_SA(n,kt) - sa_sik)*Flux(n,kt)

! SIK
           DO m = 1, NG
              m_hom = ind_cond_ng(m)
            XS_HOM_SIK(m_hom, n_hom, np_hom, ns_hom)= &
                XS_HOM_SIK(m_hom, n_hom, np_hom, ns_hom) + &
                XS_SIK(m,n,kt)*Flux(n,kt)
           END DO 

        END DO ! n = 1, NG 
       END DO ! k2
      END DO ! n1
              
! Normalize the flux
      DO ns_hom =  1, N_HOM_NZR
         DO np_hom = 1, N_HOM_ZONE


           DO n_hom = 1, N_COND_NG

! XS hom flux
            XS_HOM_FLUX( n_hom, np_hom, ns_hom )= &
              flux_volume_hom( n_hom,np_hom,ns_hom )/&
              volume_hom( np_hom,ns_hom )  

            XS_HOM_D(n_hom, np_hom, ns_hom)= &
              XS_HOM_D(n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom,np_hom,ns_hom)

            XS_HOM_D_STR(n_hom, np_hom, ns_hom)= &
              XS_HOM_D_STR(n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom,np_hom,ns_hom)
            XS_HOM_D_STR(n_hom, np_hom, ns_hom)=1./ &
              XS_HOM_D_STR(n_hom, np_hom, ns_hom)

            XS_HOM_SF(n_hom, np_hom, ns_hom)= &
              XS_HOM_SF(n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom,np_hom,ns_hom)

            XS_HOM_SF_P(n_hom, np_hom, ns_hom)= &
              XS_HOM_SF_P(n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom,np_hom,ns_hom)

            XS_HOM_SA(n_hom, np_hom, ns_hom)= &
              XS_HOM_SA(n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom,np_hom,ns_hom)


            DO m_hom = 1, N_COND_NG
              XS_HOM_SIK(m_hom, n_hom, np_hom, ns_hom)= &
              XS_HOM_SIK(m_hom, n_hom, np_hom, ns_hom)/&
              flux_volume_hom(n_hom, np_hom,ns_hom)
            END DO ! m
         END DO ! n

!        write(*,*) 'XS_HOM_SIK(2,1)=', XS_HOM_SIK(2,1,1, 1)
!        pause

! computing K_INF
! form the RHS 
          B(1:N_COND_NG)=XS_HOM_XP(1:N_COND_NG)
          A(1:N_COND_NG,1:N_COND_NG) = 0.
                
         DO n_hom=1,N_COND_NG
! set diagonal to zero due to diffusion approximation
            XS_HOM_SIK(n_hom, n_hom, np_hom, ns_hom)=0.
!            write(*,*) n_hom, &
!             (XS_HOM_SIK(n_hom, m_hom, np_hom, ns_hom)&
!              ,m_hom=1,N_COND_NG)

! add removal due to scattering
           sa_sik = 0.
           do m_hom = 1,N_COND_NG
             sa_sik = sa_sik + &
                XS_HOM_SIK(m_hom, n_hom, np_hom, ns_hom)
!             XS_SIK(m,n,k)
           end do
            A(n_hom, n_hom) = XS_HOM_SA(n_hom, np_hom, ns_hom)+&
                                   sa_sik
         END DO
!         PAUSE

            A(:,:) =A(:,:) - &
             XS_HOM_SIK(:,:, np_hom, ns_hom)

!           write(*,*) 

           CALL  MSC_LU_Solve(A, N_COND_NG, B)
! now RHS contains the solution - neutron flux spectrum
!           write(*,*) 'neutron flux spectrum=', B(:)
!           pause
          
           XS_HOM_K_INF(np_hom, ns_hom) = DOT_PRODUCT( B(1:N_COND_NG),&
              XS_HOM_SF(1:N_COND_NG, np_hom, ns_hom)        )                

! computing DB2
!  the fission source
!            B(1:N_COND_NG) = 0.
            fission_source = DOT_PRODUCT( &
             XS_HOM_SF(1:N_COND_NG, np_hom, ns_hom),&
             XS_HOM_FLUX(1:N_COND_NG, np_hom, ns_hom) ) / k_ef  
            B(1:N_COND_NG) = XS_HOM_XP(1:N_COND_NG)*fission_source
! in scattering
            B(1:N_COND_NG)=B(1:N_COND_NG) + MATMUL(&
           XS_HOM_SIK(1:N_COND_NG, 1:N_COND_NG, np_hom, ns_hom),&
           XS_HOM_FLUX(1:N_COND_NG, np_hom, ns_hom) )

!            write(*,*) 'In scattering =', MATMUL(&
!           XS_HOM_SIK(1:N_COND_NG, 1:N_COND_NG, np_hom, ns_hom),&
!           XS_HOM_FLUX(1:N_COND_NG, np_hom, ns_hom) )
!            pause

! absorption + out-scattering
           DO n_hom = 1,N_COND_NG 
            B(n_hom) = B(n_hom)-&
             (XS_HOM_SA(n_hom, np_hom, ns_hom)+&
              SUM( XS_HOM_SIK(1:N_COND_NG, n_hom, np_hom, ns_hom) ))*&
                XS_HOM_FLUX(n_hom, np_hom, ns_hom)                

            XS_HOM_DB2(n_hom, np_hom, ns_hom) = B(n_hom)/&
               XS_HOM_FLUX(n_hom, np_hom, ns_hom) 
         END DO
!         PAUSE

! computing DB2

        END DO
      END DO    

! COMPUTING POWER
! in the following assume that all fission have the same power
! pow_conv(1:NG) = CONSTANT
      pow_hom_norm = 0. 
      vol_hom_core = 0.
      XS_HOM_POW(:,:) = 0.

      DO ns_hom =  1, N_HOM_NZR
         DO np_hom = 1, N_HOM_ZONE
          DO n_hom = 1, N_COND_NG
           pow_hom_norm = pow_hom_norm + &
            XS_HOM_SF_P(n_hom, np_hom, ns_hom)*&
            XS_HOM_FLUX(n_hom, np_hom, ns_hom)*&
            volume_hom(np_hom,ns_hom)
          END DO
          vol_hom_core = vol_hom_core + volume_hom( np_hom,ns_hom )
         END DO
      END DO          

!      write(*,*) 'pow_hom_norm=', pow_hom_norm
!      pause

        pow_hom_norm = P_Reactor/(pow_hom_norm*conv_Wt_to_MWt)

!      write(*,*) 'pow_hom_norm=', pow_hom_norm
!      pause
         
      pow_hom_total = 0.

      DO ns_hom =  1, N_HOM_NZR
         DO np_hom = 1, N_HOM_ZONE
          DO n_hom = 1, N_COND_NG
           XS_HOM_POW(np_hom, ns_hom) = XS_HOM_POW(np_hom, ns_hom) + &
           XS_HOM_SF_P(n_hom, np_hom, ns_hom)*&
           XS_HOM_FLUX(n_hom, np_hom, ns_hom)*&
           volume_hom(np_hom,ns_hom)*pow_hom_norm
          END DO
            pow_hom_total = pow_hom_total + XS_HOM_POW(np_hom, ns_hom)
            XS_HOM_POW(np_hom, ns_hom)=XS_HOM_POW(np_hom, ns_hom)/&
               volume_hom( np_hom,ns_hom )
         END DO
      END DO          

      
      XS_HOM_POW(0,0) = pow_hom_total/vol_hom_core
      pow_hom_total = pow_hom_total*conv_Wt_to_MWt

                  
      RETURN
      END SUBROUTINE HOM_XS_Compute

      SUBROUTINE HOM_XS_Allocate_Index
      INCLUDE 'sketch.fh'

      ALLOCATE( ind_hom_zone(N_POLY), ind_hom_nzr(NZR), &
        ind_cond_ng(NG) )
         
      END SUBROUTINE HOM_XS_Allocate_Index  
        
      SUBROUTINE HOM_XS_DEallocate_Index
      INCLUDE 'sketch.fh'

      DEALLOCATE( ind_hom_zone, ind_hom_nzr, ind_cond_ng )
         
      END SUBROUTINE HOM_XS_Deallocate_Index  

      SUBROUTINE HOM_XS_Allocate_XS
      INCLUDE 'sketch.fh'

      ALLOCATE( &
       XS_HOM_XP(N_COND_NG),&
       XS_HOM_D (N_COND_NG, N_HOM_ZONE, N_HOM_NZR), &
       XS_HOM_D_STR (N_COND_NG, N_HOM_ZONE, N_HOM_NZR), &
       XS_HOM_SA(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),  &
       XS_HOM_SF(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_SF_P(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_FLUX(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_DB2(N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_SIK(N_COND_NG, N_COND_NG, N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_K_INF(N_HOM_ZONE, N_HOM_NZR),&
       XS_HOM_POW(0:N_HOM_ZONE, 0:N_HOM_NZR)  )  

      END SUBROUTINE HOM_XS_Allocate_XS
        
      SUBROUTINE HOM_XS_DEallocate_XS
      INCLUDE 'sketch.fh'

      DEALLOCATE( &
       XS_HOM_XP,&
       XS_HOM_D , &
       XS_HOM_D_STR , &
       XS_HOM_SA,  &
       XS_HOM_SF,&
       XS_HOM_SF_P,&
       XS_HOM_FLUX,&
       XS_HOM_DB2,&
       XS_HOM_SIK,&
       XS_HOM_K_INF,&
       XS_HOM_POW )  
         
      END SUBROUTINE HOM_XS_DEallocate_XS

      SUBROUTINE HOM_XS_Output_Parameters(unit)
      INTEGER, INTENT(IN) :: unit

      write(unit, *)

      write(unit, '(A)')&
       "   Homogenization of Neutron Cross Section parameters:"
      write(unit, '(A, 2I8)') &
       "       Number of Homogenization Zones (X-Y, Z)         :", &
       N_HOM_ZONE, N_HOM_NZR
      write(unit, '(A, I8)') &
       "       Number of Condensed Neutron Energy Groups       :", &
       N_COND_NG

!      call OUTput_Write_Separator(unit)

      RETURN 
      END SUBROUTINE HOM_XS_Output_Parameters


      SUBROUTINE HOM_XS_Output_Data(unit)
!=====================================================================*
!               Output XS   Parameters                                *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      include 'sketch.fh'
      INTEGER, INTENT(IN) :: unit

           CHARACTER*80 Header_Map
      CHARACTER*4 val_fmt
      CHARACTER*4 val_char(0:N_POLY)

      INTEGER :: np, n, m, np_hom, ns_hom

!      call OUTput_Write_Separator(io_unit)
! " Geometry Data :"
      write(unit, '(A)') " Results of the XS homogenization  :"

! " Geometry of Homogenization Zones :"
      Header_Map = " Geometry of Homogenization Zones :"
      val_fmt = "A4"
      DO np=1, N_POLY
             write(val_char(np), '(I4)')  ind_hom_zone(np)
      END DO
      call OUT_Write_Map(N_POLY,  NXR_B_Min_Reactor, &
        NXR_Max_Reactor, NXR_B_Reactor, NXR_E_Reactor, NYR_B_Reactor,&
        NYR_E_Reactor, npoly, Header_Map, unit, val_char, val_fmt)
!   
      WRITE(unit, '(A)') &
       "    Homogenization of the Axial Zones:"
      WRITE(unit, '(5x,10I3)') (ind_hom_nzr(n), n= 1, NZR)

      WRITE(unit, '(A,E14.5)') &
       "    Total power of the homogenized zone (MWt)", &
        pow_hom_total

      WRITE(unit, '(A,E14.5)') &
       "    Total volume of the homogenized zones (cm^3)", &
        vol_hom_core

      WRITE(unit, '(A,E14.5)') &
       "    Average power density of the homogenized zones(Wt/cm^3)", &
         XS_HOM_POW(0,0) 


      WRITE(unit, '(A,10E14.5)') &
       "    Group condensed prompt neutron spectrum :", &
         XS_HOM_XP(1:N_COND_NG) 


      WRITE(unit, '(A)') &
       "    Homogenized Macro Cross Sections, neutron flux etc. :"
      IF(N_COND_NG.EQ.2) then
!         input data for 2 group energy case
            WRITE(unit, '(9A14,/,9A14,/)')&
             "D(1)", "D(1/3STR)(1)", "DB2(1)",  &
             "SA(1)", "nuSF(1)", "SF(1)",  &
             "FLUX(1)", " S(1->2)", "REL POWER",&
             "D(2)", "D(1/3STR)(2)", "DB2(2)",&
             "SA(2)", "nuSF(2)", "SF(2)",  &
             "FLUX(2)", "K_INF"," NUM_HOM_ZONE"
      ELSE ! N_COND_NG.EQ.2  

            WRITE(unit, '(7A14)')&
             "D(n)", "D(1/3STR)(n)", "DB2(n)", &
             "SA(n)", "nuSF(n)", "SF(n)",  &
             "FLUX(n)"

            WRITE(unit, '(A14)')&
             "Sg->g'(n,m)"
            WRITE(unit, '(2A14)')&
             "K_INF", "REL POWER"
      END IF ! N_COND_NG.EQ.2 

      DO ns_hom = 1, N_HOM_NZR
      WRITE(unit, '(A, I4)') &
       "    Axial Layer :", ns_hom
      DO np_hom = 1, N_HOM_ZONE  
!        DO n_hom = 1, N_COND_NG            
        IF(N_COND_NG.EQ.2) then
!         input data for 2 group energy case
            WRITE(unit, '(9E14.6,/,8E14.6,I14,/)')&
                 XS_HOM_D(1,np_hom, ns_hom),&
                 XS_HOM_D_STR(1,np_hom, ns_hom),&
                 XS_HOM_DB2(1,np_hom, ns_hom),&
                 XS_HOM_SA(1,np_hom, ns_hom),&
                 XS_HOM_SF(1,np_hom, ns_hom),&
                 XS_HOM_SF_P(1,np_hom, ns_hom),&
                 XS_HOM_FLUX(1,np_hom, ns_hom),&
                 XS_HOM_SIK(2,1,np_hom, ns_hom),&
                 XS_HOM_POW(np_hom, ns_hom)/XS_HOM_POW(0,0),&
                 XS_HOM_D(2,np_hom, ns_hom),&
                 XS_HOM_D_STR(2,np_hom, ns_hom),&
                 XS_HOM_DB2(2,np_hom, ns_hom),&
                 XS_HOM_SA(2,np_hom, ns_hom),&
                 XS_HOM_SF(2,np_hom, ns_hom),&
                 XS_HOM_SF_P(2,np_hom, ns_hom),&
                 XS_HOM_FLUX(2,np_hom, ns_hom),&
                 XS_HOM_K_INF(np_hom, ns_hom), np_hom

!            WRITE(unit, '(8E14.6,/,7E14.6,I14,/)')&
!                 XS_HOM_FLUX(1,np_hom, ns_hom)/&
!                 XS_HOM_FLUX(2,np_hom, ns_hom),
!&
!                 XS_HOM_SA(2,np_hom, ns_hom)/&
!                 XS_HOM_SIK(2,1,np_hom, ns_hom)

         ELSE ! N_COND_NG /= 2
      WRITE(unit, '(A, I4)') &
       "    X-Y homogenization Zone :", np_hom

            WRITE(unit, '(7E14.6)' )&
               ( XS_HOM_D(n,np_hom, ns_hom),&
                 XS_HOM_D_STR(n,np_hom, ns_hom),&
                 XS_HOM_DB2(n,np_hom, ns_hom),&
                 XS_HOM_SA(n,np_hom, ns_hom),&
                 XS_HOM_SF(n,np_hom, ns_hom),&
                 XS_HOM_SF_P(n,np_hom, ns_hom),&
                 XS_HOM_FLUX(n,np_hom, ns_hom), n=1,N_COND_NG)

            DO n = 1, N_COND_NG
            WRITE(unit, '(10E14.6)' )&
                (XS_HOM_SIK(n,m,np_hom,ns_hom), m=1,N_COND_NG)
            END DO
            WRITE(unit,'(2E14.6)') XS_HOM_K_INF(np_hom, ns_hom),&
                 XS_HOM_POW(np_hom, ns_hom)/XS_HOM_POW(0,0)

         END IF ! N_COND_NG==2

       END DO !  np_hom
      END DO !  ns_hom

!      CALL OUTput_Write_Separator(unit)
! Deallocate all the variables
      CALL HOM_XS_Deallocate_Index  
      CALL HOM_XS_DEallocate_XS

      RETURN 
      END SUBROUTINE HOM_XS_Output_Data

      SUBROUTINE HOM_XS_OUTPUT
      INCLUDE 'units.fh' ! io_unit

!      IF(  FLAG_XS_HOPMOGENIZATION ) THEN
       OPEN(io_unit,file='Output/SKETCH.lst', status='unknown',&
               access = 'Append')

       call OUTput_Write_Separator(io_unit)


       WRITE(io_unit,'(A)') "  HOMOGENIZATION of the Neutron XS"

       call OUTput_Write_Separator(io_unit)

       CALL HOM_XS_Output_Parameters(io_unit)  

       call OUTput_Write_Separator(io_unit)

       CALL HOM_XS_Output_Data(io_unit)        

!       call OUTput_Write_Separator(io_unit)

       CLOSE(io_unit) 

!      END IF

      END SUBROUTINE HOM_XS_OUTPUT

      END MODULE HOMOGENIZATION_XS
