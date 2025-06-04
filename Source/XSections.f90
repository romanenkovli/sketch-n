      subroutine XS_Input
!=====================================================================*
!        Input from the file FILE_INPUT                               *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
! Output: xpn(NG) - Effective Prompt Fission Spectrum 
      USE termsort_lib, ONLY : init_lib

      implicit none
      include 'sketch.fh'

      integer L_HEADER, line, k, n,  m, i,  nrc, nd, i_xs_unit
      integer  ios
      parameter (L_HEADER = 5)

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      character*100 Message
      logical error_find

      write(file_inp_ident, '(I2)') LEN_INP_IDENT
      fmt_inp_ident = '(A'//file_inp_ident//')'


      open(io_unit,file = FILE_INPUT, status='OLD', action = 'READ',&
        iostat = ios)

      call Iostat_Error_Check&
      (ios,"Could not find the input FILE_INPUT file "//FILE_INPUT)

! reading the FILE_INPUT Header
       do line = 1, L_HEADER
         read(io_unit,*) 
       end do

! reading XS_MODL_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_MODL_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier XS_MODL_TYPE in the FILE_INPUT, set '//&
        'set macro cross section model to POLYNOM' )
            XS_Model = "POLYNOM"
      else

        read (io_unit,fmt=*,iostat=ios) XS_Model

        call Iostat_Error_Check&
        (ios,"error in reading XS Model "//&
        " from the FILE_INP")

        IF( (xs_model .NE. "POLYNOM").AND.&
           (xs_model .NE. "SUBSET").AND.         &
           (xs_model .NE. "SUBSED").AND.         &
           (xs_model .NE. "LINTAB").AND.         &
           (xs_model.NE. "TABLE") ) THEN
            WRITE(Message, '(A,A,A)') &
            "XS model =",&
              xs_model, &
            " in FILE_INP,  known types 'POLYNOM', "//&
            " SUBSET, and 'TABLE' 'LINTAB' "
            CALL MSC_ERR_Add_Message('ERROR',Message)
        END IF

      end if

      IF( XS_MODEL.eq."SUBSET".OR.XS_MODEL.eq."SUBSED" ) THEN
! No need for XS data in the input file
         i_xs_unit = io_unit+1
         OPEN(i_xs_unit, file = file_CD, STATUS ='OLD') 
           call init_lib( i_xs_unit )
         CLOSE(i_xs_unit)

         DO k = 1, NNODE
             DO n = 1, NG
                DO i = 1, 2
                   DO nd = 1, NDIR
                      adf(k, n, i, nd) = 1.
                   END DO
                END DO
             END DO
          END DO                                                                         
          flag_adf = .False.

        IF( XS_MODEL.eq."SUBSED" ) THEN
           CALL YMNCFALI( "Input\ymconst.lib" )
           write(*,*) 'Initialization of YMNCFALI done'
        END IF


      ELSE IF( XS_MODEL.eq."LINTAB" ) THEN
! No need for XS data in the input file
!        CAKLLinit_lib( i_xs_unit)
           CALL YMNCFALI( file_CD )

         DO k = 1, NNODE
             DO n = 1, NG
                DO i = 1, 2
                   DO nd = 1, NDIR
                      adf(k, n, i, nd) = 1.
                   END DO
                END DO
             END DO
          END DO                                                                         
          flag_adf = .False.

!      END IF

      ELSE ! XS_MODEL /= SUBSET AND XS_MODEL /= LINTAB 

! XS_BASE_DATA
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_BASE_DATA", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier XS_BASE_DATA in the FILE_INPUT' )
      else

      do k = 1, NNODE
        write(*,*) 'ng=========',NG
        if(NG.EQ.2) then
! input data for 2 group energy case
          read(io_unit,fmt=*, iostat=ios) &
                d(k,1),sik(k,2,1),sa(k,1),sf(k,1),sf_p(k,1),&
                d(k,2), sa(k,2),    sf(k,2),   sf_p(k,2), i
         else
           read(io_unit,fmt=*, iostat=ios) &
                (d(k,n), sa(k,n), sf(k,n), sf_p(k,n), &
                               n = 1, NG)
            n=1
            write(*,*) 'aa', d(k,n), sa(k,n), sf(k,n)
           read(io_unit,fmt=*, iostat=ios) &
                ((sik(k,n,m), m = 1, NG), n = 1, NG), i
         end if

        if(ios .NE. 0) then
        write(Message, '("Error Reading basic XS data, type =", I4, &
          " from the FILE_INPUT file ")' ) k
!        write(*,*) &
!                d(k,1),sik(k,2,1),sa(k,1),sf(k,1),sf_p(k,1),&
!                d(k,2), sa(k,2),    sf(k,2),   sf_p(k,2), i
        call Iostat_Error_Check(ios, Message)
        end if

      end do 

      end if


! XS_ADF_DATA (Base Assembly Discontinuity Factors)
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_ADF_DATA", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier XS_ADF_DATA in the FILE_INPUT' // &
         ' assembly  discontinuity factors set to 1')

          DO k = 1, NNODE
             DO n = 1, NG
                DO i = 1, 2
                   DO nd = 1, NDIR
                      adf(k, n, i, nd) = 1.
                   END DO
                END DO
             END DO
          END DO                                                                         
          flag_adf = .False.

      else

      DO k = 1, NNODE
! input data for 2 group energy case
         DO n = 1, NG
          read(io_unit,fmt=*, iostat=ios) &
                (( adf(k, n, i, nd), i = 1,2), nd = 1, NDD)
         END DO
      END DO           

        if(ios .NE. 0) then
        write(Message, '("Error Reading base ADF data, type =", I4, &
          " from the FILE_INPUT file ")' ) k
        call Iostat_Error_Check(ios, Message)
        end if

        flag_adf = .True.

      end if

!XSP_POL_COEF
      IF(XS_Model.EQ."POLYNOM") THEN
         rewind(io_unit)
         call MSC_Search_Header_In_File(io_unit, &
        "XSP_POL_COEF", input_line, fmt_inp_ident, error_find)  

         if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XSP_POL_COEF in the FILE_INPUT'//&
         'set XS polynomial coefficients to 0' )


            do k = 1, NNODE
              do i = 1, N_FEEDBACK
                 do n = 1, NG
                 d_fb(i, n, k) = 0.
                 sa_fb(i, n, k) = 0.
                 sf_fb(i, n, k) = 0.
                 sf_p_fb(i, n, k) = 0.
                 do m = 1, NG
                    sik_fb(i, n, m, k) = 0.
                 end do ! m
              end do ! n
              fdback00(i,k) = 1.
              end do ! i
            end do ! k nnode

         else

         do k = 1, NNODE
          if(NG.eq.2) then
             do i = 1, N_FEEDBACK
             read(io_unit, fmt=*, iostat=ios) &
                          d_fb(i,1,k),sik_fb(i,2,1,k),sa_fb(i,1,k),&
                                       sf_fb(i,1,k),sf_p_fb(i,1,k),&
                          d_fb(i,2,k),sa_fb(i,2,k),sf_fb(i,2,k),&
                                       sf_p_fb(i,2,k), fdback00(i,k)
              end do ! Feddback
           else
              do i = 1, N_FEEDBACK
              read(io_unit, fmt=*, iostat=ios) &
                 (d_fb(i, n, k),sa_fb(i, n, k),sf_fb(i, n, k), &
                   sf_p_fb(i, n, k), n = 1, NG)
              read(io_unit,*) ((sik_fb(i,n,m,k),m=1,NG) ,n=1,NG),  &
                               fdback00(i,k)
              end do ! Feddback
           end if ! NG.eq.2
           if(ios .NE. 0) then
           write(Message, '("Error Reading XS Polynomial Coefficients,"&
            //"type =", I4, " from the FILE_INPUT file ")' ) k
              call Iostat_Error_Check(ios, Message)
           end if
        end do ! NNODE

        end if
      END IF ! XS_Model.EQ.POLYNOM
      END IF ! XS_MODEL.ne.SUBSET   


!XS_DIFF_FLAG
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_DIFF_FLAG", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_DIFF_FLAG in the FILE_INPUT'//'&
         set I_Diff_Coeff to 1' )
         
              i_Diff_Coeff = 1

      else

       read(io_unit, fmt=*, iostat=ios) I_Diff_Coeff

       call Iostat_Error_Check(ios,"Error Reading I_Diff_Coeff from &
          the FILE_INPUT file " )
      
      end if      


! initial feedback distribution
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_FEED_INIT", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_FEED_INIT in the FILE_INPUT'//'&
         set fdback0(N_FEEDBACK) = 1' )

          do i = 1, N_FEEDBACK 
             fdback0(i) = 1.
          end do

      else

         read(io_unit, fmt=*, iostat=ios) (fdback0(i),i=1,N_FEEDBACK)
         call Iostat_Error_Check(ios,&
          "Error Reading initial feedbacks &
          FDBACK0 from the FILE_INPUT file " )

      end if 


!XSP_POL_COEF
      IF(XS_Model.EQ."POLYNOM") THEN

!XS_CROD_COEF
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_CROD_COEF", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_CROD_COEF in the FILE_INPUT'//&
         'set control rod differential XS  to 0.' )

        do nrc = 1, NN_CRod_Comp
           do n = 1, NG
              d_ca(n,nrc) = 0.
              sa_ca(n,nrc) = 0.
              sf_ca(n,nrc) = 0.
              sf_p_ca(n,nrc) = 0.
              do m = 1, NG
                sik_ca(n,m,nrc) = 0.
              end do ! m
           end do  ! NG
        end do ! nrc
      else 
      if(NG.eq.2) then
        do nrc = 1, NN_CRod_Comp
        read(io_unit,  fmt=*, iostat=ios) &
                      d_ca(1,nrc),sik_ca(2,1,nrc),&
                      sa_ca(1,nrc), sf_ca(1,nrc),sf_p_ca(1,nrc),&
                      d_ca(2,nrc),sa_ca(2,nrc), sf_ca(2,nrc), &
                      sf_p_ca(2,nrc)
         if(ios .NE. 0) then
              write(Message, '("Error Reading CR Differential XS, &
      CR type =", I4, " from the FILE_INPUT file ")' ) nrc
              call Iostat_Error_Check(ios, Message)
         end if
        end do
      else
        do nrc = 1, NN_CRod_Comp
         read(io_unit, fmt=*, iostat=ios) &
                  (d_ca(n,nrc),  sa_ca(n,nrc), sf_ca(n,nrc), &
                   sf_p_ca(n,nrc),  n = 1, NG)
         read(io_unit,  fmt=*, iostat=ios) &
             (( sik_ca(n,m,nrc), m = 1, NG),n=1, NG )
         if(ios .NE. 0) then
              write(Message, '("Error Reading CR Differential XS, &
      CR type =", I4, " from the FILE_INPUT file ")' ) nrc
              call Iostat_Error_Check(ios, Message)
         end if
        end do
      end if ! NG
      end if

      END IF ! XS_Model.EQ.POLYNOM

! reading XS_CROD_CUSP
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_CROD_CUSP", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier XS_CROD_CUSP in the FILE_INPUT, set '//&
        'homogenization of XS for the partial rodded nodes'//&
        ' to the flux-weighting method' )
            xs_crod_cusp = "FLUX-WEIGHTING"
      else

        read (io_unit,fmt=*,iostat=ios) xs_crod_cusp

        call Iostat_Error_Check&
        (ios,"error in reading XS_CROD_CUSP "//&
        " from the FILE_INP")

        IF( (xs_crod_cusp.NE. "FLUX-WEIGHTING").AND.&
           (xs_crod_cusp.NE. "VOLUME-WEIGHTING") ) THEN
            WRITE(Message, '(A,A,A)') &
            "homogenization of XS for the partial rodded nodes",&
              xs_crod_cusp, &
            " in FILE_INP,  known types 'FLUX-WEIGHTING', "//&
            " and 'VOLUME-WEIGHTING'"
            CALL MSC_ERR_Add_Message('ERROR',Message)
        END IF


      end if

! XS_NEUT_SPEC
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_NEUT_SPEC", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier XS_NEUT_SPEC in the FILE_INPUT')

      else

       write(*,*) 'NG = ', NG

       read(io_unit, fmt=*, iostat=ios) (xp(n), n = 1, NG)

       call Iostat_Error_Check(ios,"Error Reading Prompt Neutron &
       Spectrum xp(NG) from  the FILE_INPUT file " )

      end if

! XS_POWR_CONV

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_POWR_CONV", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('ERROR',' Could not find '//&
         'identifier XS_POWR_CONV in the FILE_INPUT')

      else

       read(io_unit, fmt=*, iostat=ios) (pow_conv(n), n = 1, NG)

       call Iostat_Error_Check(ios,"Error Reading Power Release&
       pow_conv(NG) from  the FILE_INPUT file ")
      
      end if

! XS_POWR_COOL

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "XS_POWR_COOL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_POWR_COOL in the FILE_INPUT'//&
         ' set direct coolant heating 0' )
         cool_heating = 0.

      else

       read(io_unit, fmt=*, iostat=ios) cool_heating

       call Iostat_Error_Check(ios,"Error Reading Part of the Power&
      Directly Deposited into the Coolant from FILE_INPUT file ")
      
      end if

!XS_NEUT_VELC


      IF( FILE_CD /=  "")  I_flag_vel_Library  = 1

      IF(File_DMP_Out_Kin .NE."".OR.Problem_Type.EQ."Kinetics") THEN
! Reading Neutron KInetics Data
         IF( I_flag_vel_Library == 0 ) THEN
         rewind(io_unit)
         call MSC_Search_Header_In_File(io_unit, &
        "XS_NEUT_VELC", input_line, fmt_inp_ident, error_find)  

         if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_NEUT_VELC in the FILE_INPUT' //&
         ' set neutron velocities to 1' )

             do n = 1, NG
                v(n) = 1.
            end do

            else
   
         read(io_unit,  fmt=*, iostat=ios) (v(n), n = 1,NG)

         call Iostat_Error_Check(ios,"Error Reading Neutron Velocities"&
        //"v(NG) from the FILE_INPUT file")

         end if
         END IF


         IF( FILE_CD /=  "") I_flag_deln_Library  = 1


         IF( I_flag_deln_Library == 0 ) THEN


!XS_PREC_ALFA
         rewind(io_unit)
         call MSC_Search_Header_In_File(io_unit, &
        "XS_PREC_ALFA", input_line, fmt_inp_ident, error_find)  

         if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_PREC_BETA in the FILE_INPUT' //&
         ' set delayed neutron constants to 0' )

          do m = 1, MD
             alfa(m) = 0.
          end do

         else
             
         read(io_unit,fmt=*, iostat=ios) (alfa(m),m=1,MD)
         call Iostat_Error_Check(ios,"Error Reading Delayed  "//&
        "Neutron Constant alfa(MD) from the FILE_INPUT file")

         end if

!XS_PREC_BETA
         rewind(io_unit)
         call MSC_Search_Header_In_File(io_unit, &
        "XS_PREC_BETA", input_line, fmt_inp_ident, error_find)  

         if(error_find) then

         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier XS_PREC_BETA in the FILE_INPUT' //&
         ' set delayed neutron yields to 0' )
          
          bet = 0.
          do m = 1, MD
             beta(m) = 0.
          end do

         else

         read(io_unit,fmt=*, iostat=ios) i_beta

         call Iostat_Error_Check(ios,"Error Reading i_beta "//&
        " from the FILE_INPUT file")

         if(i_beta.eq.1)  read(io_unit,*) bet

!         write(*,*) 'bet =', bet
!         read(*,*)    

         call Iostat_Error_Check(ios,"Error Reading total "//&
     "yield of delayed neutrons bet (i_bet = 1) the FILE_INPUT file")
       
         read(io_unit,fmt=*, iostat=ios) (beta(m),m=1,MD)

         call Iostat_Error_Check(ios,"Error Reading Delayed "//&
           "Neutron Constant beta(MD) from the FILE_INPUT file")

         end if

         END IF ! IF( I_flag_deln_Library == 0 ) THEN

!XS_PREC_SPEC
         rewind(io_unit)
         call MSC_Search_Header_In_File(io_unit, &
        "XS_PREC_SPEC", input_line, fmt_inp_ident, error_find)  

         if(error_find) then

         call MSC_ERR_Add_Message('Warning',' Could not find '//&
         'identifier XS_PREC_SPEC in the FILE_INPUT' //&
         ' set delayed neutron spectrum to 1. 0.' )
          
          xm(1) = 1.
          do n = 2, NG
             xm(n) = 0.
          end do

         else

         read(io_unit,fmt=*, iostat=ios) (xm(n),n=1,NG)

         call Iostat_Error_Check(ios,"Error Reading Delayed Neutron "&
        //" Spectrum xm(NG) from  the FILE_INPUT file " )

         end if
      END IF ! IF(File_DMP_Out_Kin /= "") THEN

! Decay Heat Data
        i_flag_dh = 0 
        !XS_PREC_ALFA
                 rewind(io_unit)
                 call MSC_Search_Header_In_File(io_unit, &
                "XS_DECH_ALFA", input_line, fmt_inp_ident, error_find)  
        
                 if(error_find) then
        
                 call MSC_ERR_Add_Message('WARNING',' Could not find '//&
                 'identifier XS_DECH_ALFA in the FILE_INPUT' //&
                 ' set delayed neutron constants to 0' )
        
                  do m = 1, MH
                     lam_dc(m) = 1.
                  end do
        
                 else
                     
                 read(io_unit,fmt=*, iostat=ios) (lam_dc(m),m=1,MH)
                 call Iostat_Error_Check(ios,"Error Reading Decay  "//&
                "Heat Constant lam_dc(MH) from the FILE_INPUT file")
        
                 end if
        
        !XS_PREC_BETA
                 rewind(io_unit)
                 call MSC_Search_Header_In_File(io_unit, &
                "XS_DECH_BETA", input_line, fmt_inp_ident, error_find)  
        
                 if(error_find) then
        
                 call MSC_ERR_Add_Message('WARNING',' Could not find '//&
                 'identifier XS_DECH_BETA in the FILE_INPUT' //&
                 ' set delayed neutron yields to 0' )
                  
                  e_dc = 0.
                  do m = 1, MH
                     ej_dc(m) = 0.
                  end do
        
                 else
        
                 read(io_unit,fmt=*, iostat=ios) (ej_dc(m),m=1,MH)
        
                 call Iostat_Error_Check(ios,"Error Reading Delayed "//&
                   "Neutron Constant ej_dc(MH) from the FILE_INPUT file")
        
                  e_dc = SUM ( ej_dc(1:MH) )
                  i_flag_dh = 1
                  
        
                 end if
        
      close(io_unit)

! if(DEBUG) call Copy_Input_XS
          
      return
      end
             

      subroutine XS_Init
!=====================================================================*
!        Initialize Macro_Cross Section Data                          *
! for the Steady-State & Neutron Kinetics Calculations                * 
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
! Output: xpn(NG) - Effective Prompt Fission Spectrum 
      implicit none
      include 'sketch.fh'

      integer  k, n, i, n1, m

      isotope_name(1) = "i-135"
      isotope_name(2) = "Xe135"
      isotope_name(3) = "Pm149"
      isotope_name(4) = "Sm149"

      IF( XS_MODEL.ne."SUBSET".AND.XS_MODEL.ne."SUBSED" ) THEN

       if(I_Diff_Coeff.ne.1) then
         do k = 1, NNODE
          do n = 1,NG
             d(k,n) = 1./(d(k,n)*3.)
          end do
         end do
       end if

!adding buckling in the case of 2D calculations
      IF(NZ.eq.1) THEN
          do k = 1,NNODE
             do n=1,NG
               sa(k,n) = sa(k,n) + b2(n)*d(k,n)
             end do
           end do
      END IF
      END IF

      IF( XS_MODEL.eq."SUBSET".or.XS_MODEL.eq."SUBSED" ) THEN
! i-135
        lambda_isotope(1) = ALOG(2.0)/2.3652E+04 
! xe - 135
        lambda_isotope(2) = ALOG(2.0)/3.2904E+04
! pm-149
        lambda_isotope(3) = ALOG(2.0)/1.9109E+05

      ELSE IF (XS_MODEL.eq."LINTAB") THEN
       lambda_isotope(1) = 2.86e-5
       lambda_isotope(2) = 2.07e-5
         lambda_isotope(3) = 3.51288e-6
      END IF ! IF( XS_MODEL.ne."SUBSET" ) THEN

!Set initial feedback distribution 
      if(FILE_DMP_IN .EQ. "" ) then ! else feedbacks from the restart file
         do i = 1, N_FEEDBACK
            do n1 = 1, NZ
               do k = 1, NH
                  fdback(k,n1,i) = fdback0(i)
               end do
             end do 
         end do
      end if

       if(Problem_Type.NE."Kinetics") then
          do n = 1, NG
           xpn(n) = xp(n)
          end do
       end if


      IF(File_DMP_Out_Kin .NE."".OR.Problem_Type.EQ."Kinetics") THEN
       IF ( XS_MODEL /= "SUBSET" ) THEN
       do n=1,NG
          al(n) = 1./ v(n)
       end do
       ELSE
          al(:) = 0.  
       END IF

       IF( I_flag_vel_Library == 0 ) THEN
          DO n = 1, NG
            xs_al(n, 1:N_TOT) = al(n)
          END DO
       END IF

!        write(*,*) ' xs_al(1,1)=', xs_al(1,1), xs_al(2,1)
!        pause
 
       if(i_beta.eq.0) then
         bet = 0.
         do m = 1, MD
             bet = bet + beta(m)
         end do
       else if(i_beta.eq.1) then
         do m = 1, MD
            beta(m) = beta(m)*bet
       end do
       end if
      ELSE
        bet = 0.0065 ! for the output of the reactivity only  
      END IF ! File_DMP_Out_Kin .NE."".OR.Problem_Type.EQ."Kinetics" 

!      write(*,*) 'bet = ', bet
!      write(*,*) 'al = ', al
!      pause



      return
      end

      subroutine XS_Step_Perturbation
!=====================================================================*
!        Step Perturbation of the Multiplication Cross Section        * 
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'         
! Input: k_ef - k_eff
!        sf(NNODE, NG) - Multiplication Cross sections
!        sf_fb(N_FEEBACK,NG,NNODE) - Feedback Coefficients of the 
!                                     multiplication cross ections
!       sf_ca(n, NN_CRod_Comp)  - Control Rod Coefficients of
!                                    the multiplication Cross Section
! Output: 
!        sf(NNODE, NG)/k_ef  ; 
!        sf_fb(N_FEEBACK,NG,NNODE) / k_ef 
!        sf_ca(n, NN_CRod_Comp) / k_ef
! Loca Variables: 
      real Step_Pert
      data Step_Pert / 0.25 / ! 0.25 beta
      integer n, k, nr, i 
      
      real k_eff_pert, reakt_pert

      reakt_pert = step_pert*bet

      k_eff_pert = 1./(1. - reakt_pert)

      do n=1,NG
        do k=1,NNODE
          sf(k,n)=sf(k,n)*k_eff_pert
        end do
      end do   

      do i = 1, N_FEEDBACK
        do n = 1, NG
           do k  = 1, NNODE
             sf_fb(i,n,k) = sf_fb(i,n,k)*k_eff_pert
           end do
        end do
      end do

      do nr = 1, NN_CRod_Comp 
         do  n = 1, NG
           sf_ca(n, nr)  = sf_ca(n, nr)*k_eff_pert
         end do
      end do


      return 
      end

      subroutine XS_Set_Criticality
!=====================================================================*
!        Normalization of the Fission Cross Section for Neutron       * 
!              Kinetics Calculations                                  *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'         
! Input: k_ef - k_eff
!        sf(NNODE, NG) - Multiplication Cross sections
!        sf_fb(N_FEEBACK,NG,NNODE) - Feedback Coefficients of the 
!                                     multiplication cross ections
!       sf_ca(n, NN_CRod_Comp)  - Control Rod Coefficients of
!                                    the multiplication Cross Section
! Output: 
!        sf(NNODE, NG)/k_ef  ; 
!        sf_fb(N_FEEBACK,NG,NNODE) / k_ef 
!        sf_ca(n, NN_CRod_Comp) / k_ef
! Loca Variables: 
      integer n, k, nr, i 
      
      do n=1,NG
        do k=1,NNODE
          sf(k,n)=sf(k,n)/k_ef
        end do
      end do   

      do i = 1, N_FEEDBACK
        do n = 1, NG
           do k  = 1, NNODE
             sf_fb(i,n,k) = sf_fb(i,n,k)/k_ef
           end do
        end do
      end do

      do nr = 1, NN_CRod_Comp 
         do  n = 1, NG
           sf_ca(n, nr)  = sf_ca(n, nr)/ k_ef
         end do
      end do

      write(*,*) 'Set criticality', 'k_eff = ', k_ef
      read(*,*)


      return 
      end

      subroutine XS_Set_Crit
!=====================================================================*
!        Normalization of the Fission Cross Section for Neutron       * 
!              Kinetics Calculations  ALL XS MODELS                   *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'         
! Input: k_ef - k_eff
!        XS_SF(NG, N_TOT) - Multiplication Cross sections
!                                    the multiplication Cross Section
! Output: 
!        XS_SF(NG, N_TOT)/k_ef  ; 
      Integer :: n, kt      
      
      do n=1,NG
        do kt = 1, N_TOT
          XS_SF(n, kt) = XS_SF(n, kt)/k_ef 
        end do
      end do   

      return 
      end


      subroutine XS_Update 
!=====================================================================*
!        Update of the Macro Cross Sections                           *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'

!      REAL, PARAMETER :: D1_multiplier = 0.9
!        LOGICAL, PARAMETER :: Flag_change_D1 = .True.  
        
      call CPU_Timer(time_sk_cpu)

      if(XS_Model .EQ. "POLYNOM") then
        call XSP_Compute
        call XSP_Compute_CRD
          call XS_Add_Removal
      else if(XS_Model .EQ. "SUBSET".OR.XS_MODEL.eq."SUBSED") then  
        call XSP_Compute_Subset
        call XSP_Compute_CRD_Subset
          call XSP_Compute_Deln_Par
        CALL XS_Add_Xe_Sm
          call XS_Add_Removal
!          IF ( Flag_change_D1 ) THEN
!            CALL XS_Sensitivity_CR            
!            CALL XS_Sensitivity_Reflector
!          END IF         
      else if(XS_Model .EQ. "LINTAB") then  
        call XSP_Compute_LINTABLE
        call XSP_Compute_CRD_LINTABLE
        CALL XS_Add_Xe_Sm
        call XS_Add_Removal
      else if(XS_Model .EQ. "TABLE") then
        if(FILE_CD .NE. "" ) then
           call XSR_Compute_BWR_XS
           call XSR_Compute_BWR_CR_XS
        end if
      else 

        write(*,*) 'XS Mode =', TRIM(XS_MODEL)
        write(*,*) 'This model is not defined'
        stop

      end if

      if(Problem_Type.EQ."Kinetics".AND. &
       iflag_divide_keff == 1) call XS_Set_Crit

! TMP
!      CALL XS_OUTPUT_DEBUG

      call CPU_Timer(time_XS_cpu)


!      TMP 

      return
      end

      Subroutine  XS_Compute_Xappa_Rod(h_rod,e_rod_round_off, &
                 n1, n1_rod, n1_nrod, k, ir, ie, &
                 Xappa_Rod)
!=====================================================================*
! Computing Flux-Weighting Factor for XS homogenization for CR        *
!  Slava 24.VI.1999                                                   *
!=====================================================================*
      IMPLICIT NONE

      INCLUDE 'sketch.fh'


!Input:
      real h_rod, e_rod_round_off
      integer n1, n1_rod, n1_nrod, kt, k, ir, ie
!Output:      
      real Xappa_rod(NG)
! Local: 
      integer kt_rod, kt_nrod, n
      real hz_rod, hz_nrod, alfa_rod, beta_rod
      real flux_rod(NG), flux_nrod(NG)

      real CRD_Xappa_OUtput(NG, NN_CRod), CRD_h_rod_output(NN_CRod)
      common /CRD_Xappa_save/ crd_h_rod_output, crd_xappa_output

! flat neutron flux for volume weghting
      do n = 1, NG
! Volume-Weighting Homogenization Procedure
         xappa_rod(n) = h_rod
      end do

      if(xs_crod_cusp.EQ."FLUX-WEIGHTING") then 
 
         if((h_rod+e_rod_round_off).LT.1) then
! if there is a partially rodded node 

            if(ie.eq.1.AND. &
                     rod_node( ir, n1_nrod , 1 ).eq.0) then
! if it is an absorber and it is the bottom part of it 
!                     absorber is the 1st element of the node
              if(n1.ne.1 .and. n1.ne.NZ) then
! if it is not the 1st and the last axial layer
                     kt      = k + (n1 - 1)*NH
                     kt_rod  = k + (n1_rod - 1)*NH
                     kt_nrod = k + (n1_nrod - 1)*NH
                     hz_rod = h_rod*hzt(n1)
                     hz_nrod = (1. - h_rod)*hzt(n1)
                     alfa_rod = 1. / (1. + hzt(n1_rod)/hz_rod)
                     beta_rod = 1. / (1. + hzt(n1_nrod)/hz_nrod)
                    do n = 1, NG
! Neibouring Spectrum Index Method
!                     IF (n == 1) THEN
                       flux_rod(n) = alfa_rod*flux(n,kt) + &
                          (1. - alfa_rod)*flux(n,kt_rod)
                       flux_nrod(n) = beta_rod*flux(n,kt) + &
                          (1. - beta_rod)*flux(n,kt_nrod)
!                     ELSE
!                       flux_rod(n) = &
!                       (flux(n,kt_rod)/flux(1,kt_rod))*flux_rod(1) 
!                       flux_nrod(n) = &
!                       (flux(n,kt_nrod)/flux(1,kt_nrod))*flux_nrod(1) 
!                     END IF 
! IF negative flux in the control rod absorber we set him to 0.
!                     if(flux_rod(n).LT.0) flux_rod(n) = 0.
                      xappa_rod(n) = flux_rod(n)*h_rod/ &
                    ( flux_rod(n)*h_rod + flux_nrod(n)*(1.-h_rod) )
                     IF( xappa_rod(n).LT.0) xappa_rod(n) = h_rod
                     IF( xappa_rod(n).GT.1) xappa_rod(n) = h_rod
!                     write(*,*) 'All'
!                     select case (ir)
!                       case(1)
!                        xappa_rod(1)=0.60482E+00
!                        xappa_rod(2)=0.52575E+00
!                       case(5)
!                        xappa_rod(1)=0.60482E+00
!                        xappa_rod(2)= 0.52575E+00
!                       case(7)
!                        xappa_rod(1)=0.60021E+00 
!                        xappa_rod(2)=0.52064E+00
!                       case(10)
!                        xappa_rod(1)=0.60021E+00
!                        xappa_rod(2)= 0.52064E+00
!                       case(11)
!                        xappa_rod(1)=0.58691E+00
!                        xappa_rod(2)=0.50658E+00
!                       case default
!                         write(*,*) 'something wrong '
!                         stop
!                     end select

                       
                    end do

                    do n = 1, NG
                       CRD_Xappa_Output(n, ir) = xappa_rod(n)
                    end do

                       CRD_h_rod_output(ir) = h_rod

             end if ! not the 1st and the last axial zone
          end if ! if it is an absorber and it is the bottom part of it
      end if ! if there is a partially rodded node 
      end if ! if we want a flux-weighting procedure

!           IF(k.eq.169.and.n1.eq.3) THEN
!                     kt = k + (n1 - 1)*NH
!                     kt_rod = k + (n1_rod - 1)*NH
!                     kt_nrod = k + (n1_nrod - 1)*NH
!              write(*,*) 'k =', 'n1=', n1
!              write(*,*) 'n1_rod, n1_nrod =', n1_rod, n1_nrod
!              write(*,*) 'h_rod =', h_rod, hz_rod, hz_nrod
!              write(*,*) 'flux(n,kt_nrod)=', flux(:, kt_nrod)
!              write(*,*) 'flux(n,kt)=', flux(:, kt)
!              write(*,*) 'flux(n,kt_rod)=', flux(:, kt_rod)
!              write(*,*) 'XS_SF(n,kt)=', XS_SF(:, kt)
!              write(*,*) 'XS_SF_P(n,kt)=', XS_SF_P(:, kt)
!             write(*,*) 'flux_rod =', flux_rod
!             write(*,*) 'flux_nrod =', flux_nrod
!              write(*,*) 'xappa_rod =', xappa_rod
!              pause
!           END IF


      return
      end


      subroutine XS_DBG_output
!=====================================================================*
!   DEBUG SUBOUTINE - OUTPUT OF THE CROSS-SECTIONS                    *
!    (c) Slava 2.III.1998  JAERI                                      *
!=====================================================================*
      include 'sketch.fh'
      integer n, k, nn
      REAL  s_12(NH), sr2(NH), sr1(NH), nuSF2(NH), nuSF1(NH), rinf(NH), &
           KINF(NH) 
            
      real CRD_Xappa_OUtput(NG, NN_CRod), CRD_h_rod_output(NN_CRod)
      common /CRD_Xappa_save/ crd_h_rod_output, crd_xappa_output
      
      open(io_unit,file = 'Output_Debug/xs.dat', status = 'unknown')
          do n= 1, NZ
           nn = (n-1)*NH
           write(io_unit,*) 'AXIAL ZONE NUMBER =', n                 
           write(io_unit,*) 'Volume' 
           write(io_unit,2) (volume(k+nn),k=1,NH)
           write(io_unit,*) 'Fast Group Diffusion Coefficient' 
           write(io_unit,1) (XS_D(1,k+nn),k=1,NH)
           write(io_unit,*) 'Thermal Group Diffusion Coefficient' 
           write(io_unit,1) (XS_D(2,k+nn),k=1,NH)
           write(io_unit,*) 'Removal Cross Section n =1'
           write(io_unit,1) (XS_SA(1,k+nn)/volume(k+nn),k=1,NH)
           write(io_unit,*) 'Removal Cross Section n = 2'
           write(io_unit,1) (XS_SA(2,k+nn)/volume(k+nn),k=1,NH)
           write(io_unit,*) 'Scattering Cross Section 1 -> 2'
           write(io_unit,1) &
            (XS_SIK(2,1,k+nn)/volume(k+nn),k=1,NH)

          write(io_unit,*) 'Multiplication Cross Section n = 2'
          write(io_unit,1) (XS_SF(2,k+nn)/volume(k+nn),k=1,NH)

          write(io_unit,*) 'KINF'
          DO k = 1, NH
           s_12(k)=XS_SIK(2,1,k+nn)/volume(k+nn)
           sr2(k)=XS_SA(2,k+nn)/volume(k+nn)
           sr1(k)=XS_SA(1,k+nn)/volume(k+nn)
           nuSF2(k)=XS_SF(2,k+nn)/volume(k+nn)
           nuSF1(k)=XS_SF(1,k+nn)/volume(k+nn)
         rinf(k) =  s_12(k)/sr2(k)
           KINF(k) =  (nuSF1(k) + rinf(k)*nuSF2(k))/sr1(k)
          END DO
          write(io_unit,1) (kinf(k),k=1,NH)

           write(io_unit,*) 'Coolant Temperature'
           write(io_unit,2) (fdback(k,n,2),k=1,NH)
           write(io_unit,*) 'Fuel Temparature'
           write(io_unit,2) (fdback(k,n,4),k=1,NH)
           write(io_unit,*) 'Coolant density'
           write(io_unit,1) (fdback(k,n,3),k=1,NH)
           write(io_unit,*) 'Boron concentration'
           write(io_unit,2) (fdback(k,n,1),k=1,NH)
           end do
       close(io_unit)

    1 format(10es14.5)
!    1 format(10f8.5)
    2 format(10f8.1)

!      open(io_unit,file = 'Output_Debug/CRD_Xappa.dat', &
!           status = 'unknown')
!        write(io_unit,'(A)') 'Weighting factors for Partially'// &
!                  ' Rodded Nodes'
!        write(io_unit, '(A)') ' rod number, rodded part, xappa(NG)'
!        do ir = 1, NN_CROd
!          write(io_unit, '(I4, 3E12.5)') ir, CRD_h_rod_output(ir),&
!          (CRD_Xappa_OUtput(n, ir), n = 1, NG)
!        end do

!      close(io_unit)

      return
      end


      real function XS_Sect_CR(sect_nrod, sect_rod, xappa_rod)
      real sect_nrod, sect_rod, xappa_rod          

      XS_Sect_CR= sect_nrod*(1.-xappa_rod) + sect_rod*xappa_rod

      return
      end

      subroutine XS_Add_Xe_Sm
!=====================================================================*
! Add absorption due to Xe and Sm to macro absorption Xs              *
!                         Slava (c) 04.X.2004                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      INTEGER :: n1, nn, kc, k, kt, n, i
      REAL :: a, vol
!      REAL :: conve

      DO n1 =  NZ_Core_BEG, NZ_Core_End 
         nn = (n1 - 1)*NH
         DO kc = 1, NH_Core
            k = np_core(kc)
            kt = k + nn 
            vol = volume(kt)
            a = 0.
              DO n = 1, NG
               a = 0. 
                 DO  i = 1, N_ISOTOPE
                 a = a + sa_isotope(i, n, kt)*conc_isotope(k,n1,i)*&
                     vol
               END DO
               XS_SA(n, kt)= XS_SA(n, kt)+a
            END DO
         END DO
      END DO

      RETURN
      END subroutine XS_Add_Xe_Sm


      subroutine XS_Add_Removal
!=====================================================================*
! Current Value of Removal Cros Section (Absorption + Scattering)     *
!                         Slava (c) 19.III.1998                       *
!=====================================================================*
      implicit none
       include 'sketch.fh'
! Input: XS_SIK(NG,NG,N_TOT), XS_SA(NG, N_TOT)
! Output: XS_SA(NG, N_TOT) - Removal XS * Volume
! Local Variables
      real sa_sik
      integer k, n, m


      do k = 1, N_TOT
        do n = 1, NG
           sa_sik = 0.
           do m = 1, NG
             sa_sik = sa_sik + XS_SIK(m,n,k)
           end do
           XS_SA(n, k) = XS_SA(n, k) + sa_sik
        end do
       end do

       if(Debug) call XS_DBG_output

      return
      end

      SUBROUTINE XS_Output_Parameters(unit)
      implicit none
      include 'sketch.fh'

      integer unit

      write(unit, *)

      write(unit, '(A)')&
       "   Cross Section parameters:"
      write(unit, '(A, I8)') &
       "       Number of Material Compositions                 :", &
       NNODE
      write(unit, '(A, I8)') &
       "       Number of Neutron Energy Groups                 :", &
       NG
      write(unit, '(A, I8)') &
       "       Number of Delayed Neutron Precursor Groups      :", &
       MD

      write(unit, '(A, I8)') &
       "       Number of Feedbacks                             :", &
       N_FEEDBACK

      write(unit, '(A, 1x, A10)') &
       "       Macro Cross Section Model                       :", &
       XS_Model

      IF(XS_Model.EQ."TABLE" .AND. FILE_CD.NE."" ) THEN
      write(unit, '(A, 1x, A10)') &
       "       File Name of the XS Table Data                  :", &
       FILE_CD
      END IF

      write(unit, '(A, 1x, A)') &
       "       Homogenization Method for Partially Roded Nodes :", &
       xs_crod_cusp

!      call OUTput_Write_Separator(unit)

      RETURN 
      END

      subroutine XS_Output_Data(unit)
!=====================================================================*
!               Output XS   Parameters                                *
!                   Vyachreslav Zimin (c) April 1 1998                *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer unit
      integer k, nrc, m, n, i, nd
      real d_out(NNODE,NG)

      if(I_Diff_Coeff.ne.1) then
         do k = 1, NNODE
          do n = 1,NG
             d_out(k,n) = 1./(d(k,n)*3.)
          end do
         end do
      ELSE
         do k = 1, NNODE
          do n = 1,NG
             d_out(k,n) = d(k,n)
          end do
         end do
      end if


      write(unit, '(A)') " Macro Cross Section Data :"

      WRITE(unit, '(A, I2)') &
       "    Diffusion Coefficient Flag in the Input Data : ", &
      I_Diff_Coeff


! "    Basic Macro Cross Section Set:"
      WRITE(unit, '(A)') &
       "    Basic Macro Cross Section Set:"

      DO k = 1, NNODE
        IF(NG.EQ.2) then
!         input data for 2 group energy case
            WRITE(unit, '(5E14.6,/,4E14.6,I14,/)')&
                d_out(k,1),sik(k,2,1),sa(k,1),sf(k,1),sf_p(k,1),&
                d_out(k,2), sa(k,2),    sf(k,2),   sf_p(k,2), k

         ELSE
               WRITE(unit, '(A,I4)') &
       "           Material Composition number :", k

            WRITE(unit, '(4E14.6)' )&
                (d_out(k,n), sa(k,n), sf(k,n), sf_p(k,n), &
                               n = 1, NG)
            DO n = 1, NG
            WRITE(unit, '(10E14.6)' )&
                (sik(k,n,m), m = 1, NG)
            END DO
         END IF
      END DO
! "    Assembly Discontinuity Factors"
      IF(Flag_ADF) THEN  
      WRITE(unit, '(A)') &
       "    Assembly Discontinuity Factors: Left, Right, Up, Down,"//&
            "Bottom, Top)"
      DO k = 1, NNODE
               WRITE(unit, '(A,I4)') &
       "           Material Composition number :", k
! input data for 2 group energy case
         DO n = 1, NG
          WRITE(io_unit,'(10E14.6)') &
                (( adf(k, n, i, nd), i = 1,2), nd = 1, NDD)
         END DO
      END DO           
      END IF


! "    Polynomial Coefficients for Macro Cross Sections"
      IF(XS_Model.EQ."POLYNOM") THEN
            WRITE(unit, '(A)') &
       "    Polynomial Coefficients for Macro Cross Sections"

        DO k = 1, NNODE
               WRITE(unit, '(A,I4)') &
       "           Material Composition number :", k
          IF(NG.eq.2) THEN
             DO i = 1, N_FEEDBACK
             WRITE(unit, '(5E14.6)') &
                          d_fb(i,1,k),sik_fb(i,2,1,k),sa_fb(i,1,k),&
                                       sf_fb(i,1,k),sf_p_fb(i,1,k),&
                          d_fb(i,2,k),sa_fb(i,2,k),sf_fb(i,2,k),&
                                       sf_p_fb(i,2,k), fdback00(i,k)
             END DO ! Feddback
          ELSE
             DO i = 1, N_FEEDBACK
               WRITE(unit, '(4E14.6)') &
                 (d_fb(i,n,k),sa_fb(i,n,k),sf_fb(i,n,k), &
                   sf_p_fb(i,n,k), n = 1, NG)
              DO n = 1, NG
               WRITE(unit, '(10E14.6)') &
             (sik_fb(i, n, m, k ),m=1,NG) 
              END DO
              WRITE(unit, '(10E14.6)')  fdback00(i,k)
             END DO ! N_FEEDBACK
          END IF ! (NG.EQ.2)
         END DO ! NNODE
      END IF ! (XS_Model.EQ."POLYNOM") 

! "    Initial Feedback Values                         "
      WRITE(unit, '(A)') &
       "    Initial Feedback Values                         "
      WRITE(unit, '(5E14.6)') (fdback0(i),i = 1, N_FEEDBACK)

! "    Differential Cross Sections for Control Rods    "
      WRITE(unit, '(A)') &
       "    Differential Cross Sections for Control Rods    "
        IF(NG.eq.2) then
           DO nrc = 1, NN_CRod_Comp
              WRITE(unit, '(5E14.6,/,4E14.6,I14,/)')&
                      d_ca(1,nrc),sik_ca(2,1,nrc),&
                      sa_ca(1,nrc), sf_ca(1,nrc),sf_p_ca(1,nrc),&
                      d_ca(2,nrc),sa_ca(2,nrc), sf_ca(2,nrc), &
                      sf_p_ca(2,nrc), nrc
           END DO
        ELSE
           DO nrc = 1, NN_CRod_Comp
               WRITE(unit, '(A,I4)') &
       "           CR Material Composition number :", nrc
            WRITE(unit, '(4E14.6)') &
                  (d_ca(n,nrc),  sa_ca(n,nrc), sf_ca(n,nrc), &
                   sf_p_ca(n,nrc),  n = 1, NG)
            WRITE(unit, '(4E14.6)')&
                  (( sik_ca(n,m,nrc), m = 1, NG),n=1, NG )
           END DO
        END IF

! "    Prompt Neutron Spectrum:   "
      WRITE(unit, '(A)') &
       "    Prompt Neutron Spectrum:   "
      WRITE(unit, '(5F10.5)') (xp(i),i = 1, NG)

! "    Energy Release per Fission:   "
      WRITE(unit, '(A)') &
       "    Energy Release per Fission:   "
      WRITE(unit, '(5E14.6)') (pow_conv(i),i = 1, NG)

! "    Part of the Power Directly Deposited into Coolant:   "
      WRITE(unit, '(A)') &
       "    Part of the Power Directly Deposited into Coolant:   "
      WRITE(unit, '(5E14.6)') cool_heating

      IF(File_DMP_Out_Kin .NE. "".OR.Problem_Type.EQ."Kinetics") THEN

! "    Neutron Velocities:   "
      WRITE(unit, '(A)') &
       "    Neutron Velocities:   "
      WRITE(unit, '(5E14.6)') (v(i),i = 1, NG)

! "    Delayed Neutron Decay Constants:   "
      WRITE(unit, '(A)') &
       "    Delayed Neutron Decay Constants:   "
      WRITE(unit, '(6E13.6)') (alfa(i),i = 1, MD)

! "    Delayed Neutron Yields:   "
      WRITE(unit, '(A)') &
       "    Delayed Neutron Yields:   "
      WRITE(unit, '(6E13.6)') (beta(i),i = 1, MD)

! "    Total Yield of Delayed Neutrons:   "
      WRITE(unit, '(A)') &
       "    Total Yield of Delayed Neutrons:   "
      WRITE(unit, '(E13.6)') bet

! "    Delayed Neutron Spectrum:   "
      WRITE(unit, '(A)') &
       "    Delayed Neutron Spectrum:   "
      WRITE(unit, '(5F10.5)') (xm(i),i = 1, NG)

      END IF !       IF(File_DMP_Out_Kin /= "") THEN

      call OUTput_Write_Separator(io_unit)

      return
      end


      SUBROUTINE XS_Multiply_D1 ( D1_multiplier )
      implicit none
      include 'sketch.fh'
      REAL, INTENT(IN) :: D1_multiplier
      INTEGER :: k

          DO k= 1, N_TOT
! multipliying only in the reactor core
            IF ( XS_SF(2,k).LT.1.E-10) THEN 
              XS_D(1,k) = XS_D(1,k)*D1_multiplier
            END IF
          END DO             
      RETURN
      END SUBROUTINE XS_Multiply_D1 

      SUBROUTINE XS_Multiply_D1_CR ( D1_multiplier )
      implicit none
      include 'sketch.fh'
      REAL, INTENT(IN) :: D1_multiplier
      INTEGER :: ie, ir, ib, np , nch, n1, kt, k
      REAL    ::  h_rod

      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 


      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
       do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              h_rod = rod_node(ir,n1,ie)

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node
              XS_D(1, kt) =  XS_D(1, kt)*D1_multiplier
           end if
          end do ! n1
         end if ! (k.NE.0)
        end do ! nch
      end do ! ie 
       end do ! ir
      end  do ! ib       
      RETURN
      END SUBROUTINE XS_Multiply_D1_CR 

      SUBROUTINE XS_Multiply_ADF_CR ( D1_multiplier )
      implicit none
      include 'sketch.fh'
      REAL, INTENT(IN) :: D1_multiplier
      INTEGER :: ie, ir, ib, np , nch, n1, kt, k
      REAL    ::  h_rod

      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 


      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
       do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              h_rod = rod_node(ir,n1,ie)

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node

              XS_D(2, kt) =  XS_D(2, kt)/D1_multiplier
              XS_SA(2, kt) =  XS_SA(2, kt)/D1_multiplier


           end if
          end do ! n1
         end if ! (k.NE.0)
        end do ! nch
      end do ! ie 
       end do ! ir
      end  do ! ib       
      RETURN
      END SUBROUTINE XS_Multiply_ADF_CR

      SUBROUTINE XS_Multiply_SF1_CR ( D1_multiplier )
      implicit none
      include 'sketch.fh'
      REAL, INTENT(IN) :: D1_multiplier
      INTEGER :: ie, ir, ib, np , nch, n1, kt, k
      REAL    ::  h_rod

      real e_rod_round_off
      parameter(e_rod_round_off = 1.E-03) 


      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
       do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              h_rod = rod_node(ir,n1,ie)

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node
              XS_SF(1, kt) =  XS_SF(1, kt)*D1_multiplier
              XS_SF_P(1, kt) =  XS_SF_P(1, kt)*D1_multiplier
           end if
          end do ! n1
         end if ! (k.NE.0)
        end do ! nch
      end do ! ie 
       end do ! ir
      end  do ! ib       
      RETURN
      END SUBROUTINE XS_Multiply_SF1_CR 


      SUBROUTINE XS_Sensitivity_CR
      implicit none
      include 'sketch.fh'
!      REAL, INTENT(IN) :: D1_multiplier
      INTEGER :: ie, ir, ib, np , nch, n1, kt, k, n
      REAL    ::  h_rod

      real e_rod_round_off, delta_adf(NG), delta_diff(NG)
      parameter(e_rod_round_off = 1.E-03)
      CHARACTER*100 fname
      LOGICAL exists
 
      fname = "Input/Sensitivity_CR_ADF.dat"       

!     Get the name of a file:

!     INQUIRE about file's existence:
      INQUIRE (FILE = fname, EXIST = exists)


      IF (exists) THEN

      write(*,*) 'CHANGING DIFF & ADF of FA with CR' 
      
      OPEN(io_unit, file = fname, action='read')
      READ(io_unit, *)  delta_adf(1:NG), delta_diff(1:NG)
      CLOSE(io_unit)

      WRITE(*,*) 'delta adf, d=',  delta_adf(1:NG), delta_diff(1:NG)       

      do ie = 1, NN_CRod_El
        do ir = 1, NN_CRod
       do ib = 1, NN_CRod_Bundle
         np = nrods(ir, ib)
         do nch = 1, NCHM
           k = poly_out(np,nch)
           if(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
! 
!           do n1 = NZ_Core_BEG, NZ ! 
              kt = k + (n1-1)*NH
              h_rod = rod_node(ir,n1,ie)

           if((h_rod - e_rod_round_off).GT.0) then 
! if there is control rod in this node in the node

           DO n = 1, NG
              XS_D(n, kt) =  XS_D(n, kt) + XS_D(n, kt)*delta_diff(n)
              xs_adf(n, 1:2, 1:(NDIR-1), kt) = &
              xs_adf(n, 1:2, 1:(NDIR-1), kt) + &
              xs_adf(n, 1:2, 1:(NDIR-1), kt)*delta_adf(n)
          END DO

           end if
          end do ! n1
         end if ! (k.NE.0)
        end do ! nch
      end do ! ie 
       end do ! ir
      end  do ! ib       

      END IF ! exists

      RETURN
      END SUBROUTINE XS_Sensitivity_CR

      SUBROUTINE XS_Sensitivity_Reflector
      implicit none
      include 'sketch.fh'
      INTEGER ::  n1, kt, k, n, nl, np, ns, nn

      real e_rod_round_off, delta_adf(NG), delta_diff(NG)
      parameter(e_rod_round_off = 1.E-03)
      CHARACTER*100 fname
      LOGICAL exists
 
      fname = "Input/Sensitivity_REFLECTOR.dat"       

!     Get the name of a file:

!     INQUIRE about file's existence:
      INQUIRE (FILE = fname, EXIST = exists)


      IF (exists) THEN

      write(*,*) 'CHANGING DIFF & ADF of REFLECTOR' 
      
      OPEN(io_unit, file = fname, action='read')
      READ(io_unit, *)  delta_adf(1:NG), delta_diff(1:NG)
      CLOSE(io_unit)

      WRITE(*,*) 'delta adf, d=',  delta_adf(1:NG), delta_diff(1:NG)       

       do n1 = 1, NZ 
        nn = (n1-1)*NH
        do k = 1, NH

          kt = k + nn

          np = np_out(k)
          ns = ns_out(n1)
          nl = l(np,ns)

! multipliying only in the reflector (only radial)
           IF ( ( nl == 1).OR.( nl == 2).OR.( nl == 3)) THEN 
            DO n = 1, NG
              XS_D(n, kt) =  XS_D(n, kt) + XS_D(n, kt)*delta_diff(n)
              xs_adf(n, 1:2, 1:(NDIR-1), kt) = &
              xs_adf(n, 1:2, 1:(NDIR-1), kt) + &
              xs_adf(n, 1:2, 1:(NDIR-1), kt)*delta_adf(n)
           END DO
          END IF
       END DO    
       END DO             

      END IF ! exists

      RETURN
      END SUBROUTINE XS_Sensitivity_REFLECTOR

      SUBROUTINE   XS_OUTPUT_DEBUG
      implicit none
      include 'sketch.fh'
      INTEGER ::  n1,  k, n

      OPEN(io_unit, File='Output_Debug/XS.dat', STATUS='UNKNOWN')

      WRITE(io_unit,'(A)') '   SA'
        DO k = 1, NH
          DO n = 1, NG          
          WRITE(io_unit,'(I10,I5, 100ES14.4)') &
                 k, n,(XS_SA(n,k + (n1-1)*NH),n1=1, NZ)
        END DO
        END DO

      WRITE(io_unit,'(A)') '   SF'
        DO k = 1, NH
          DO n = 1, NG          
          WRITE(io_unit,'(I10,I5, 100ES14.4)') &
                 k, n,(XS_SF(n,k + (n1-1)*NH),n1=1, NZ)
        END DO
        END DO

      WRITE(io_unit,'(A)') '   D'
        DO k = 1, NH
          DO n = 1, NG          
          WRITE(io_unit,'(I10, I5, 100ES14.4)') &
                 k, n,(XS_D(n,k + (n1-1)*NH),n1=1, NZ)
        END DO
        END DO

      WRITE(io_unit,'(A)') '   SF_P'
        DO k = 1, NH
          DO n = 1, NG          
          WRITE(io_unit,'(I10, I5, 100ES14.4)') &
                 k, n,(XS_SF_P(n,k + (n1-1)*NH),n1=1, NZ)
        END DO
        END DO

       CLOSE(Io_unit)
       write(*,*) 'XS are written'
       read(*,*)      

       RETURN
       END            
      
      
