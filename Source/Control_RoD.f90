      subroutine CRD_Init
      USE CR_POSITION
      USE CR_LINEAR

      implicit none
      INTEGER, PARAMETER :: FLAG_CR_FROM_FILE=0

      CALL SET_INDEX_CR_in_BANK

      IF( Input_Bank_Position ) THEN

         IF( Bank_Position_in_Percents ) &
            CALL SET_ALL_CR_POSITION_PERC
                         
         IF( Bank_Position_Absolute ) &
            CALL SET_ALL_CR_POSITION_ABS

      END IF

         write(*,*) 'FLAG_CR_FROM_FILE=', FLAG_CR_FROM_FILE

         IF(FLAG_CR_FROM_FILE==1) &
          CALL  CR_LINEAR_READ


      RETURN
      END

      subroutine CRD_Input
!=====================================================================*
!   Input Control Rod Data from the FILE_INPUT                         *
!                         Slava (c) 17.IX.1999                        *
!=====================================================================*
      USE CR_POSITION
      implicit none
      include 'sketch.fh'

!Input: from the FILE_INPUT
!      nrods(NN_CRod, NN_Crod_Bundle)
!      zrods(NN_CRod)
!      h_rod_el(NN_CRod_El)
!      v_rods
!      NROD_MOVING
!      NROD_MOVING
!      Flag_SCRAM
! Settimng the following values
!        Flag_Rod_Top = .False.
!        Flag_Rod_Bottom = .False.
!        Time_scram = 1.E+30
!        if(Flag_Scram) then
! Setting Scram Time During Calvculations
!           Flag_Set_Time_Scram = .False.
!         else
! NO Scram
!           Flag_Set_Time_Scram = .True.
!         end if
!      end if

      character*12 input_line
      character*5 fmt_inp_ident
      character*2 file_inp_ident
      character*100 Message
      logical error_find
      integer i, ie, ios, ir, n,  nr , n1, ib, j, i_cr



!initialization of the identifiers
      write(file_inp_ident, '(I2)') LEN_INP_IDENT
!      fmt_inp_ident = '(A'//trim(file_inp_ident)//')'
      fmt_inp_ident = '(A'//file_inp_ident//')'

!Starting input
      open(io_unit,file = FILE_INPUT,status='old', iostat=ios)
      call Iostat_Error_Check&
      (ios,"Could not find the input FILE_INPUT file "//FILE_INPUT)


! reading CRD_COR_LOAD
!      write(*,*) 'CRD_COR_LOAD'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_COR_LOAD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_COR_LOAD in the FILE_INPUT, types '//&
        'of the control rods set to 1')
         do ir = 1, NN_CRod
            CR_Load(ir) = 1
         end do
      else
 
        read(io_unit, fmt=*, iostat=ios) &
           ( CR_Load( ir ), ir = 1, NN_CRod)
        call Iostat_Error_Check(ios,"Error Reading CR "//&
            "Loading CR_Load(NN_CRod) from the FILE_SECT file")
 
      end if

! reading CRD_MAT_COMP
!      write(*,*) 'CRD_MAT_COMP'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_MAT_COMP", input_line, fmt_inp_ident, error_find)  
      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_MAT_COMP in the FILE_INPUT, material '//&
         'compositions of the control rods set to 1')

       do n = 1, NN_CRod_Type
          do ie = 1, NN_CRod_El
            CR_Compos(ie,n)=1
          end do
       end do

      else
       do n = 1, NN_CRod_Type
         read(io_unit,  fmt=*, iostat=ios) &
                   (CR_Compos(ie, n), ie = 1, NN_CRod_El)
         if(ios .NE. 0) then
              write(Message, '("Error Reading CR Material "//&
                    "Compositions, CR Type =", I4, &
                    " from the FILE_SECT file ")' ) n
              call Iostat_Error_Check(ios, Message)
         end if
       end do
      end if

! reading CRD_LOCATION   
!      write(*,*) 'CRD_LOCATION'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_LOCATION", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_LOCATION in the FILE_INPUT, set locations '//&
         'of the control rods  to 0')
         do i = 1, NN_CRod
            do ib=1,NN_CRod_Bundle
               nrods(i, ib) = 1
            end do 
         end do
      else
         read (io_unit,fmt=*,iostat=ios) ( (nrods(i, ib), &
            ib=1,NN_CRod_Bundle), i=1,NN_CRod)
         call Iostat_Error_Check&
        (ios,"Error in Reading Control Rod Locations nrods(NN_CRod),"//&
         "from the FILE_INPUT")
      end if

! reading CRD_POSITION   
!      write(*,*) 'CRD_POSITION'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_POSITION", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_POSITION in the FILE_INPUT, set positions '//&
        'of the control rods  to 0')
         do i = 1, NN_CRod
           zrods(i) = 0
         end do
      else

        read (io_unit,fmt=*,iostat=ios) (zrods(i),i=1,NN_CRod)

        call Iostat_Error_Check&
        (ios,"Error in Reading Control Rod Positions zrods(NN_CRod),"//&
        " from the FILE_INPUT")

      end if
! reading CRD_H_ROD_EL   
      write(*,*) 'CRD_H_ROD_EL'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_H_ROD_EL", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message("WARNING"," Could not find "//&
         "identifier CRD_H_ROD_EL in the FILE_INPUT, set heights"//&
         " of the control rod elements  to 0")
         do ie=1,NN_CRod_El
           h_rod_el(ie) = 0
         end do
      else
        read (io_unit,fmt=*,iostat=ios) (h_rod_el(ie),ie=1,NN_CRod_El)
        call Iostat_Error_Check&
       (ios,"Error in Reading CR elements lengths " //&
        "h_rod_el(NN_CRod_El) from the FILE_INPUT")
      end if

! Control Rod Bank Information
! reading 
!      write(*,*) 'CRD_POSITION'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_BNK_INDX", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_BNK_INDX the FILE_INPUT, no information '//&
        'on CR bank is available')
         GO TO 100
      else

        READ(io_unit,fmt=*,iostat=ios) N_CR_BANK
        call Iostat_Error_Check&
        (ios,"Error in Reading CR Bank Information N_CR_BANK,"//&
        "under identifier CRD_BNK_INDX")

        READ(io_unit,fmt=*,iostat=ios) N_CR_in_Bank(1:N_CR_BANK) 

        call Iostat_Error_Check&
        (ios,"Error in Reading CR Bank N_CR_in_Bank(1:N_CR_BANK) ,"//&
        "under identifier CRD_BNK_INDX")

! Checking that the number of CR is equal to the number of CR in BANKS

        N_CR_TOTAL = 0
        DO n = 1, N_CR_BANK
           N_CR_TOTAL = N_CR_TOTAL + N_CR_in_Bank(n) 
        END DO

        IF( N_CR_TOTAL /= NN_CRod ) THEN
           ios = -100
           call Iostat_Error_Check&
        (ios,"Number of the Control Rods in Banks is not equal ,"//&
        "to the total number of Control Rods NN_CRod")
        END IF

      end if


! * reading CRD_BNK_LOAD
!      write(*,*) 'CRD_BNK_LOAD'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_BNK_LOAD", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_BNK_LOAD in the FILE_INPUT, types '//&
        'of the control rods set to 1')

        N_CR_TYPE_in_BANK(1:N_CR_BANK)=1  

      else
 
        read(io_unit, fmt=*, iostat=ios) &
           N_CR_TYPE_in_BANK(1:N_CR_BANK)
        call Iostat_Error_Check(ios,"Error Reading CR "//&
            "Loading CR_Load(NN_CRod) from the FILE_SECT file")
! setting up the CR loads from the BANK load types
         call MSC_ERR_Add_Message('WARNING',' SEt UP CR types '//&
        'from the the BANK types ')

        i_cr = 0
        DO n = 1, N_CR_BANK
           DO j = 1, N_CR_in_Bank(n) 
              i_cr = i_cr + 1
              CR_Load(i_cr) = N_CR_TYPE_in_BANK(n)
          END DO
        END DO

!      write(*,*) 'CR_load(:)=', CR_Load(1:10)
!      pause
                      
      end if

! reading CRD_BNK_POSA
!      write(*,*) 'CRD_POSITION'

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_BNK_POSA", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_BANK_POSA in the FILE_INPUT')

      else

        read (io_unit,fmt=*,iostat=ios) BANK_POS(1:N_CR_BANK)

        call Iostat_Error_Check&
        (ios,"Error in Reading Control Rod BANK Positions ,"//&
        " under identifier CRD_BNK_POSA")

        Input_Bank_Position = .True.
        Bank_Position_Absolute = .True.

      end if

      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_BNK_POS%", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
        'identifier CRD_BANK_POS% in the FILE_INPUT')

      else

        read (io_unit,fmt=*,iostat=ios) BANK_POS(1:N_CR_BANK)

        call Iostat_Error_Check&
        (ios,"Error in Reading Control Rod BANK Positions ,"//&
        " under identifier CRD_BNK_POSA")

        Input_Bank_Position = .True.
        Bank_Position_in_Percents = .True.        

      end if

  100 CONTINUE 

      IF(Problem_Type.EQ."Kinetics") THEN


! reading CRD_ROD_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_ROD_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_ROD_TYPE in the FILE_INPUT, set    '//&
         'CRD_ROD_TYPE to Normal (No moving Fuel) ')
          CRD_ROD_Type = "NORMAL"
      else
        read(io_unit,fmt=*,iostat=ios) &
         CRD_ROD_Type
         call Iostat_Error_Check&
        (ios,"Error in Reading Type of the CR movement "//&
         "CRD_ROD_Type, from the FILE_INPUT")
      end if

!      WRITE(*,*)    "CRD_ROD_Type =",     CRD_ROD_Type
!      PAUSE



! reading CRD_MOV_TYPE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_MOV_TYPE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_MOV_TYPE in the FILE_INPUT, set    '//&
         'CRD Moving Type to LINear ')
          CRD_Mov_Type = "LIN"
      else
        read(io_unit,fmt=*,iostat=ios) &
         CRD_Mov_Type
         call Iostat_Error_Check&
        (ios,"Error in Reading Type of the CR movement "//&
         "CRD_Mov_Type, from the FILE_INPUT")
      end if


! reading CRD_NUM_MOVE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_NUM_MOVE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_NUM_MOVE in the FILE_INPUT, set number'//&
         'of the moving control rods to 0')
         NROD_MOVING = 0
      else
         read(io_unit,fmt=*,iostat=ios) NROD_MOVING
         call Iostat_Error_Check&
        (ios,"Error in Reading Number of Moving CR NROD_MOVING, "//&
         "from the FILE_INPUT")


      end if

! reading CRD_TIM_STRT
!      write(*,*) 'CRD_TIM_STRT'
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
          "CRD_TIM_STRT", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_TIM_STRT in the FILE_INPUT, set start '//&
         'time of the control rod movement  to 0')
          CRD_Time_Start(:) = 0.
      else
        read(io_unit,fmt=*,iostat=ios) &
         (CRD_Time_start(i), i=1, NROD_MOVING)
         call Iostat_Error_Check&
        (ios,"Error in Reading Start Time of the CR movement "//&
         "CRD_Time_Start, from the FILE_INPUT")
      end if

! reading CRD_TIM_END
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_TIM_END", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_TIM_END in the FILE_INPUT, set end   '//&
         'time of the control rod movement  to 1.E+30')
          CRD_Time_End(:) = 1.E+30
      else
        read(io_unit,fmt=*,iostat=ios) &
         (CRD_Time_End(i), i=1, NROD_MOVING)
         call Iostat_Error_Check&
        (ios,"Error in Reading End Time of the CR movement "//&
         "CRD_Time_End, from the FILE_INPUT")
      end if

! reading CRD_IND_MOVE
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_IND_MOVE", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_IND_MOVE in the FILE_INPUT, set indexes'//&
         'of the moving control rods to 0')
         do i= 1, NROD_MOVING
           Num_Mov_Rod(i) = 0
         end do
      else

         read(io_unit,fmt=*,iostat=ios) &
            (Num_Mov_Rod(i), i= 1, NROD_MOVING)

         call Iostat_Error_Check&
        (ios,"Error in Reading Indexes of Moving CR "//&
         "Num_Mov_Rod(NROD_MOVING), from the FILE_INPUT")

      end if


      if(CRD_Mov_Type.EQ."LIN") then 
! reading CRD_VEL_RODS ! only for LINEAR PERTURBATION 
        rewind(io_unit)
        call MSC_Search_Header_In_File(io_unit, &
        "CRD_VEL_RODS", input_line, fmt_inp_ident, error_find)  

        if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_VEL_RODS in the FILE_INPUT, set velocity '//&
         'of the moving control rods to 0')
          v_rods = 0.
        else
         read(io_unit,fmt=*,iostat=ios) (v_rods(i), i=1, NROD_MOVING)
         call Iostat_Error_Check&
        (ios,"Error in Reading CR velocity v_rods, from FILE_INPUT")
        end if

      else if(CRD_Mov_Type.EQ."SIN") then 

! CRD_SIN_AMPL
        rewind(io_unit)
        call MSC_Search_Header_In_File(io_unit, &
        "CRD_SIN_AMPL", input_line, fmt_inp_ident, error_find)  

        if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_SIN_AMPL in the FILE_INPUT, set amplitude '//&
         'of the moving (Sinus) control rods to 0')
         do i= 1, NROD_MOVING
           CRD_Sin_Ampl(i) = 0.
         end do

         else
           read(io_unit,fmt=*,iostat=ios) &
                (CRD_Sin_Ampl(i),i=1,NROD_MOVING)
           call Iostat_Error_Check&
        (ios,"Error in Reading CR Sinus Amplitude CRD_Sin_Freq, &
        from FILE_INPUT")
        end if


! CRD_SIN_FREQ
        rewind(io_unit)
        call MSC_Search_Header_In_File(io_unit, &
        "CRD_SIN_FREQ", input_line, fmt_inp_ident, error_find)  

        if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_SIN_FREQ in the FILE_INPUT, set frequency '//&
         'of the moving (Sinus) control rods to 1')
         do i= 1, NROD_MOVING
           CRD_Sin_Freq(i) = 1.
         end do
        else

        read(io_unit,fmt=*,iostat=ios) &
                (CRD_Sin_Freq(i),i=1,NROD_MOVING)
         call Iostat_Error_Check&
        (ios,"Error in Reading CR Sinus Frequency (CRD_SIN_FREQ),"//&
         " from FILE_INPUT")
        end if

      end if ! "SIN".OR."LIN"         

! reading CRD_TOP_POSI
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_TOP_POSI", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_TOP_POSI in the FILE_INPUT, set start '//&
         'set TOP Position to the Reactor Top Position ')

          Z_Rod_Top = 0.
          do n1 = 1, NZR
             Z_RoD_Top = Z_Rod_Top + hz(n1)
          end do

      else
        read(io_unit,fmt=*,iostat=ios) &
         Z_Rod_Top
         call Iostat_Error_Check&
        (ios,"Error in Reading Maximum Top Position of CR "//&
         "Z_Rod_Top, from the FILE_INPUT")
      end if

! reading CRD_BTM_POSI
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_BTM_POSI", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_BTM_POSI in the FILE_INPUT, set start '//&
         'set Minimum Bottom Position to 0. ')

          Z_Rod_Bottom = 0.

      else
        read(io_unit,fmt=*,iostat=ios) &
         Z_Rod_Bottom
         call Iostat_Error_Check&
        (ios,"Error in Reading Minimum Bottom Position of CR "//&
         "Z_Rod_Bottom, from the FILE_INPUT")
      end if

! reading CRD_FLG_SCRM
      rewind(io_unit)
      call MSC_Search_Header_In_File(io_unit, &
        "CRD_FLG_SCRM", input_line, fmt_inp_ident, error_find)  

      if(error_find) then
         call MSC_ERR_Add_Message('WARNING',' Could not find '//&
         'identifier CRD_FLG_SCRM in the FILE_INPUT, set flag '//&
         'of the flag of the scram signal to .False.')

           Flag_Scram = .False.

      else

           Flag_Scram = .True.        

! Input of the data for scram
         read(io_unit,fmt=*,iostat=ios) v_rods_scram
         call Iostat_Error_Check&
      (ios,"SCRAM, Error in Reading Number of CR scram velocity,"//&
       "v_rods_scram from the FILE_INPUT")

         read(io_unit,fmt=*,iostat=ios) Nrod_Scram
         call Iostat_Error_Check&
        (ios,"SCRAM, Error in Reading Number of scram CR, "//&
        "Nrod_Scram from the FILE_INPUT")

         read(io_unit,fmt=*,iostat=ios) &
          (Num_Rod_Scram(i), i= 1, Nrod_Scram)
         call Iostat_Error_Check&
      (ios,"SCRAM, Error in Reading Indexes of scram CR, "//&
      "Num_Rod_Scram(Nrod_Scram) from the FILE_INPUT")

         read(io_unit,fmt=*,iostat=ios) Power_Scram
         call Iostat_Error_Check&
      (ios,"SCRAM, Error in Reading Power Level for Scram, "//&
       "Power_Scram from the FILE_INPUT")

         read(io_unit,fmt=*,iostat=ios) Time_Scram_Delay
         call Iostat_Error_Check&
      (ios,"SCRAM, Error in Reading Time Delay of Scram "//&
       "Time_Scram_Delay from the FILE_INPUT")

      end if

      END IF ! IF(Problem_Type.EQ."Kinetics") THEN

      close(io_unit)


! Input: NN_CRod, nrods(NN_CRod)
! Output: NRS(N_POLY) - Control Rod Number in Assembly np 
!                            (Looks Like No need)
!         Mov_Rods(NN_CRod) - Indicator the Control Rod has moved to 
!                 recalculate the material properties
! Local Variables


! Setting up their locations

!      do nr = 1, NN_CRod
!         np = nrods(nr)
!         nrs(np) = nr
!      end do

! initialization control rods routine
      do nr = 1 , NN_CRod
         Mov_Rods(nr) = .TRUE.
      end do

! Initial Control Rod Position
      do nr = 1 , NN_CRod
         zrods_ini(nr) = zrods(nr)
         zrods_dt(nr) = zrods(nr)
      end do

! Initialization of the material compositions for CR
       do ie = 1, NN_CRod_El
          do ir = 1, NN_CRod
            Mat_Com_Rod(ir,ie) = CR_Compos(ie, CR_Load(ir) )
          end do
       end do

      return
      end 

      subroutine CRD_Compute_Mat_Comp
!=====================================================================*
!     Calculation of the control rod elements in the node             *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
! The program modification to compute BWR control rods 21.VI.1999     *
!=====================================================================*
      implicit none
      include 'sketch.fh'
! Input: NN_CRod, Mov_Rods(NN_CRod), zrods(NN_CRod), NN_CRod_El, H_Rod_El(NN_CRod_El),
!         NZ, HZT(NZ), REACTOR_TYPE
! Output: Rod_Node(NN_CRod, NZ, NN_CRod_El) - Control Rod Element in the Node
!                      for the Rod NN_CRod in Axial Layer NZ 

! Local Variables
      integer ir, n1, ie
      integer NZ_CR_START, NZ_CR_END, NZ_CR_STEP
      real h_0, h_1, z0_pos(NN_CRod_El), z1_pos(NN_CRod_El),&
        point_up, point_down

      if(REACTOR_TYPE.eq."PWR") then

         NZ_CR_START = 1
         NZ_CR_END = NZ
         NZ_CR_STEP = 1
      
      else if(REACTOR_TYPE.eq."BWR") then
    
         NZ_CR_START = NZ
         NZ_CR_END = 1
         NZ_CR_STEP = -1

      end if


          do ir = 1,NN_CRod
! recalculate the material composition if control rod has moved
          if(Mov_Rods(ir)) then
             z0_pos(1) = zrods(ir)
             z1_pos(1) =  z0_pos(1) +  h_rod_el(1)
             do ie = 2, NN_CRod_El
               z0_pos(ie) =  z0_pos(ie-1) + h_rod_el(ie-1)
               z1_pos(ie) =  z0_pos(ie) + h_rod_el(ie)
             end do                   

             do ie = 1, NN_CRod_El
               h_0 = 0.
               h_1 = 0.
               do n1 = NZ_CR_START, NZ_CR_END, NZ_CR_STEP
                  h_1 = h_1 + hzt(n1)
                  point_up = amin1(z1_pos(ie), h_1)
                  point_down = amax1(z0_pos(ie), h_0)
                  if(point_up .GT. point_down) then 
                    rod_node(ir,n1,ie)=(point_up-point_down)/hzt(n1)
                  else
                    rod_node(ir,n1,ie) = 0.
                  end if
                  h_0 = h_0 + hzt(n1)
               end do ! NZ
            end do ! NN_CRod_El
          end if ! Moving
        end do ! NN_CRod


       if(Debug) then
          CALL CRD_Debug_Output_Rod_Comp
       end if

      return
      end 


!      subroutine CRD_Move_Rods_PWR_LMW(time, dt_rods)
!====================================================================*
! Simulation of the Control Rods Movement for the LMW Benchmark      *
!              rod group 1 moves from t=0   to t=26.6 s              *
!              rod group 2 moves from t=7.5 to t=47.5 s              *
!      speed  - 3 cm/sec                                             *
!====================================================================*      
!      implicit none
!      include 'sketch.fh'
! Input: zrods(NN_CRod) 
!      real time, dt_rods
! Output: zrods(NN_CRod) - Position of the Control Rods from the Bottom
! Mov_Rods() - Indicator that Control Rods Moved
! Local Variables:      
!     real v_rods ! Control Rod Speed NOW IN THE GEOMETRY FILE
!      parameter( v_rods = 3. ) 
!      integer NROD1, NROD2
!      parameter(NROD1 = 2, NROD2 =  2) ! 2 Groups of Control Rods
!      integer nrods1(NROD1),nrods2(NROD2)
!      data nrods1 /2, 4/
!      data nrods2 /1, 3/
!      real TM_R1, TM_R2 ! Time Moment when the Group of CR Starts Moving
!      parameter(TM_R1=0.0, TM_R2=7.5)
!      real z1_end, z2_end ! Final Position of the CR Group
!      parameter(z1_end = 180., z2_end = 60.) ! cm from the Bottom 

!      real time_rods ! Previous Time Moment when we computed CR Position 
!      real dt_r1, dt_r2
!      integer n, nn
 
!      time_rods = time - dt_rods

! 1st group of control rods
!      if(time.gt.TM_R1) then

!         if(time_rods.lt.TM_R1) then
!             dt_r1 = time - TM_R1
!         else
!             dt_r1 = dt_rods
!         end if

!         do nn = 1,NROD1
!            n = nrods1(nn)
!            zrods(n)  = zrods(n) + v_rods*dt_r1
!            if(zrods(n).GT.z1_end) then
!              zrods(n) = z1_end
!            end if
!            Mov_Rods(n) = .True.
!         end do
!      end if

! 2nd group of control rods

!      if(time.gt.TM_R2) then

!         if(time_rods.lt.TM_R2) then
!             dt_r2 = time - TM_R2
!         else
!             dt_r2 = dt_rods
!         end if

!         do nn = 1,NROD2
!            n = nrods2(nn)
!            zrods(n)  = zrods(n) - v_rods*dt_r2

!            if(zrods(n).LT.z2_end) then
!               zrods(n) = z2_end
!            end if
!            Mov_Rods(n) = .True.
!         end do
!      end if

!      write(*,*) 'time = ',time,'zrods group 1 =', zrods(nrods1(1))
!      write(*,*) 'time = ',time,'zrods group 2 =', zrods(nrods2(1))
!      write(*,*) 'time_rods = ', time_rods, dt_rods


!      return
!      end 

!      subroutine CRD_Move_Rods_BWR_LRA(time, dt_rods)
!====================================================================*
! Simulation of the Control Rods Movement for the BWR Benchmark      *
!              rods  move from t=0   to t = 2 s                      *
!      speed  - 150 cm/sec                                           *
! NOT USED
!====================================================================*      
!      implicit none
!      include 'sketch.fh'
! Input: zrods(NN_CRod) 
!      real time, dt_rods
! Output: zrods(NN_CRod) - Position of the Control Rods from the Bottom
! Mov_Rods() - Indicator that Control Rods Moved
! Local Variables:      
!      real v_rods
!      parameter( v_rods= 150. ) ! NOW IN THE GEOMETRY FILE 
!!      real Time_Start 
!!      parameter(Time_Start=0.0)
!      real Z_Rod_End
!      parameter(Z_Rod_End = 330.) !Top Reflector

!      real time_rods, dt_mov
!      integer n
 
!      time_rods = time - dt_rods

!      if(time.gt.CRD_Time_Start) then

!          if(time_rods.lt.CRD_Time_Start) then
!             dt_mov = time - CRD_Time_Start
!          else
!             dt_mov = dt_rods
!          end if

!          do n = 1, NN_CRod
!              zrods(n)  = zrods(n) + v_rods*dt_mov
!              if(zrods(n).gt.Z_Rod_End) then
!                 zrods(n) = Z_Rod_End
!              end if
!              Mov_Rods(n) = .TRUE.
!          end do
!      end if

!      write(*,*) 'time = ',time,'zrods  =', zrods(1)
!      write(*,*) 'time_rods = ', time_rods, dt_rods

!      return
!      end 

      subroutine CRD_Move_Rods(time, dt_rods)
!====================================================================!
! Simulation of the Control Rods Movement                            !
! Lastr update December 21 1999                                      !
!====================================================================!      
      implicit none
      include 'sketch.fh'
! Input: zrods(NN_CRod) 
      real time, dt_rods
! Output: zrods(NN_CRod) - Position of the Control Rods from the Bottom
! Mov_Rods() - Indicator that Control Rods Moved
! Local Variables:      
!  v_rods - Speed of the Conrol Rod (It is different for
!                                       different Cases)
!      real v_rods ! 
! CRD_Time_Start -Time when Control Rod Starts Moving 
! NROD - Number of the Moving Control Rods
! time_Rods - time passed from the previous calculation of the CA Positions
! dt_mov - time of the control rod moving during current tme step
      real time_rods, dt_mov

      integer n, i

      time_rods = time - dt_rods

!      write(*,*) 'CRD_Time_Start, CRD_Time_End = ', &
!             CRD_Time_Start, CRD_Time_End
!      write(*,*) 'NROD_MOVING =', NROD_MOVING
!      write(*,*) 'CRD_MOV_TYPE =', CRD_MOV_TYPE
!      pause

!      write(*,*) 'NROD_MOVING =', NROD_MOVING
!      write(*,*) 'CRD_Time_Start =', (CRD_Time_Start(i), &
!      i = 1, NROD_MOVING )
!      write(*,*) 'CRD_Time_End =', (CRD_Time_End(i),&
!      i = 1, NROD_MOVING )
!      write(*,*) 'vrods =', (v_rods(i),&
!      i = 1, NROD_MOVING )
!      pause


      if(time.lt.Time_Scram) then ! CR STOP Moving When Scram happens

         do i = 1, NROD_MOVING

         n = Num_Mov_Rod(i)

        if(time.gt.CRD_Time_Start(i).and.time_rods.LE.CRD_Time_End(i))&
                                                                  then

           if(CRD_Mov_Type.eq."LIN") then 

! Linear Perturbation 
            if(time_rods.lt.CRD_Time_Start(i)) then
               dt_mov = time - CRD_Time_Start(i)
            else
               dt_mov = dt_rods
            end if

              zrods_dt(n) = zrods(n)
              zrods(n)  = zrods(n) + v_rods(i)*dt_mov
              if(zrods(n) .GT. Z_Rod_Top) then
                 zrods(n) = Z_Rod_Top
              end if
              if(zrods(n) .LT. Z_Rod_Bottom) then
                 zrods(n) = Z_Rod_Bottom
              end if
              Mov_Rods(n) = .TRUE.
           else if(CRD_Mov_Type.eq."SIN") then 
! Sinusoidal Perturbation 
!            write(*,*) 'Sinus Inside =', Time

              if (Time.LE. CRD_Time_End(i)) then

!                     write(*,*) 'n =', n, 'zrods(n) = ', zrods(n)
                     zrods_dt(n) = zrods(n)
                     zrods(n)  = zrods_ini(n) + CRD_Sin_Ampl(i)*&
              (COS( 2.*PI*CRD_Sin_Freq(i)* (time-CRD_Time_Start(i))/&
              (CRD_Time_End(i)-CRD_Time_Start(i)) ) - 1. )
              else
                zrods(n)  = zrods_ini(n)
              end if
              Mov_Rods(n) = .TRUE.
          
          end if ! CRD_Mov_Type.eq."SIN" 

        end if ! Time > CRD_Time_Start & Time_Rods < CRD_Time_End
        end do ! NROD_MOVING
      end if ! Time.LT. Time_Scram

      if(Flag_Scram) then
          if((.NOT.Flag_Set_Time_Scram) .AND. &
                                 (P_Total.GT.Power_Scram)) then
            Time_Scram = Time_Rods + Time_Scram_Delay
            Flag_Set_Time_Scram = .True.
          end if

          if(time.gt.Time_Scram) then

            if(time_rods.lt.Time_Scram) then
               dt_mov = time - Time_Scram
            else
               dt_mov = dt_rods
            end if

            do i = 1, NROD_SCRAM
              n = Num_Rod_Scram(i)
              zrods_dt(n) = zrods(n)
              zrods(n)  = zrods(n) - v_rods_scram*dt_mov
              if(zrods(n).lt.Z_Rod_Bottom) then
                 zrods(n) = Z_Rod_Bottom
              end if
              if(zrods(n) .GT. Z_Rod_Top) then
                 zrods(n) = Z_Rod_Top
              end if
              Mov_Rods(n) = .TRUE.
            end do
        end if ! Time > Time_Scram
      end if ! Flag_Scram

!      write(*,*) 'time = ',time,'zrods  =', zrods(1)
!      write(*,*) 'time_rods = ', time_rods, dt_rods
!      write(*,*) 'Time_Scram =', Time_Scram

      return
      end 

      subroutine CRD_Save_Data
!====================================================================*
! Saving CRD data for restart                                        *
!====================================================================*      
      implicit none
      include 'sketch.fh'

      integer i

      do i = 1, NN_CRod
         Old_Zrods(i) = Zrods(i)
      end do

      Old_Time_Scram = Time_Scram
      Old_Flag_Set_Time_scram = Flag_Set_Time_scram

      return
      end ! subroutine CRD_Save_Data


      subroutine CRD_Recover_Data
!====================================================================*
! Recovering CRD data for restart                                    *
!====================================================================*      
      implicit none
      include 'sketch.fh'

      integer i

      do i = 1, NN_CRod
         Zrods(i) = Old_Zrods(i)
         zrods_dt(i) = zrods(i)
         Mov_Rods(i) = .TRUE.
      end do

      Time_Scram = Old_Time_Scram
      Flag_Set_Time_scram = Old_Flag_Set_Time_scram

      return
      end 


      subroutine CRD_Debug_Output_Rod_Comp
!=====================================================================*
! Output of the Control Rod positions into the file "Output/rods.dat" *
!   NEACRP PWR PEA ,      Slava (c) 19.II.1998                        *
!=====================================================================*
      implicit none
      include 'sketch.fh'

! Input: Rod_Node(NN_CRod, NZ, NN_CRod_El), NN_CRod_El, NN_CRod, NZ
! Output: File "Output/rods.dat"
! Local Variables
      integer ie, ir, n1

      open(io_unit,file = 'Output_Debug/rods.dat', status ='unknown')
        do ir = 1, NN_CRod
           write(io_unit,*) 'rods number =', ir, 'zrods=', zrods(ir)
           do ie = 1, NN_CRod_El
            write(io_unit,*) 'Rod element number =', ie
            write(io_unit,*) 'Rod_Node=',(rod_node(ir,n1,ie),n1=1,NZ)
           end do
        end do
      close(io_unit)

      return
      end 


      subroutine CRD_Compute_CRD_Position_in_Percents&
        (CR_Pos_2D, CR_Number_1D)
!=====================================================================*
! 2D collapsed CR positions  only for DEBUG                           *
! 25.VI.1999 (c) Slava                                                * 
!=====================================================================*
      implicit none
      include 'sketch.fh'

      real CR_Pos_2D(N_POLY), CR_Number_1D(NZR)

      integer   ns,  ir, ie, np,  n1, icr,  ib,&
       n1_Start, n1_End
      real  h_rod

      do np = 1, N_POLY
         CR_Pos_2D(np) = 0.
      end do

      do ns = 1, NZR
         CR_Number_1D(ns) = 0.
      end do

      do ie = 1, NN_CRod_El
         do ir = 1, NN_CRod
            icr = Mat_Com_Rod(ir,ie)
            do ib = 1, NN_Crod_Bundle         
               np = nrods(ir, ib)
               n1_End = NZ_Core_Beg -1
               do ns = NZR_Core_BEG, NZR_CORE_END ! NZ_Core_END (adding the upper reflector)
                  n1_Start = n1_End + 1
                  n1_End = n1_Start + NPZ(ns) - 1
                  do n1 = n1_Start, n1_End  
                    h_rod = rod_node(ir,n1,ie)*hzt(n1)
                    CR_Pos_2D(np) = CR_Pos_2D(np) + 100.*h_rod/hz_core
                  end do ! N1
               end do ! ns            
            end do ! ib        
         end do ! ir
      end do ! ie

      do ie = 1, NN_CRod_El
         do ir = 1, NN_CRod
            icr = Mat_Com_Rod(ir,ie)
            do ib = 1, NN_Crod_Bundle         
               np = nrods(ir, ib)
               n1_End = 0
               do ns = 1, NZR
                  n1_Start = n1_End + 1
                  n1_End = n1_Start + NPZ(ns) - 1
                  do n1 = n1_Start, n1_End  
                     h_rod = rod_node(ir,n1,ie)*hzt(n1)
                     CR_Number_1D(ns) = CR_Number_1D(ns) + h_rod / &
                                   ( real(NN_Crod_Bundle)*hz(ns) )
                  end do ! N1
               end do ! ns            
            end do ! ib        
         end do ! ir
      end do ! ie

      return
      end 


      SUBROUTINE CRD_Output_Parameters(unit)
      implicit none
      include 'sketch.fh'
      integer unit


      write(unit, *)
      write(unit, '(A)') " Control Rod Parameters:"
      write(unit, '(A, I8)') &
     "    Number of Control Rods                                 :", &
        NN_CRod
      write(unit, '(A, I8)') &
     "    Number of Control Rod Types                            :", &
        NN_CRod_Type
      write(unit, '(A, I8)') &
     "    Number of Control Rod Elements (Absorber, Driver, etc.):", &
        NN_CRod_El
      write(unit, '(A, I8)') &
     "    Number of Control Rod Material Compositions            :", &
        NN_CRod_Comp
      write(unit, '(A, I8)') &
     "    Number of Bundles Covered by Control Rod               :", &
        NN_Crod_Bundle
      

      return
      end

      SUBROUTINE CRD_Output_Data(unit)
      implicit none
      include 'sketch.fh'
      integer unit

      integer ind, i, nb, n, ie
      character*80 Header_Map
      character*4 val_fmt
      character*5 val_char(0:N_POLY)

! " Control Rod Data :"
      write(unit, '(A)') " Control Rod Data :"

! "    Location of the Control Rods in the Core  :"       
      Header_Map = "    Location of the Control Rods in the Core  :"
      val_fmt = "A5"
      ind = 0
      DO ind = 1, N_POLY 
         val_char(ind) = "    -"
      END DO
      DO i=1,NN_CRoD
        DO nb = 1, NN_CRod_Bundle
           write(val_char(nrods(i,nb)), '(I5)' ) i
        END DO
      END DO

      call OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

!"    Control Rod Types                         :"       
      Header_Map = "    Control Rod Types                         :"
      val_fmt = "A5"
      ind = 0
      DO ind = 1, N_POLY 
         val_char(ind) = "    -"
      END DO
      DO i=1,NN_CRoD
        DO nb = 1, NN_CRod_Bundle
           write(val_char(nrods(i,nb)), '(I5)' ) CR_Load(i)
        END DO
      END DO
      call OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)

!"    Control Rod Material Compositions:"
      write(unit, '(A)')&
       "    Control Rod Material Compositions:"
      DO ind = 1, NN_CRoD_Type
         write(unit, ' (8x, I3,":", 10I4, /, 12x, 10I4,/, 12x, 10I4)')&
           ind, (CR_Compos(n, ind), n = 1, NN_CRod_El)
      END DO


! "    Length of the Control Rod Elements (Absorber, Driver): "
      write(unit, '(A)')&
       "    Length of the Control Rod Elements (Absorber, Driver): "
      write(unit, ' (8x,  5F10.3)') &
         (h_rod_el(ie),ie=1,NN_CRod_El)

! "    Control Rod Positions "

!      CALL CRD_Output_CR_Positions(unit)
      IF(Problem_Type.EQ."Kinetics") THEN

!"    Type of the Control Rods (with Fuel, Normal)       : "
            write(unit, '(A)')&
       "    Type of the Control Rods (with Fuel, Normal)       : "
      write(unit, ' (8x,  A)') &
         CRD_Rod_Type

! "    Type of the Control Rod Movement                   : "
            write(unit, '(A)')&
       "    Type of the Control Rod Movement                   : "
      write(unit, ' (8x,  A)') &
         CRD_Mov_Type

! "    Number of the Moving Control Rods                  : "
      write(unit, '(A)')&
       "    Number of the Moving Control Rods                  : "
      write(unit, ' (8x, I3)') &
         NROD_Moving

! "    Indexes of the Moving Control Rods                  : "
      write(unit, '(A)')&
       "    Indexes of the Moving Control Rods                 : "
      write(unit, ' (8x, 10I4)') &
         (Num_Mov_Rod(i), i=1, NROD_Moving)

! "    Start and End Time of the Control Rod Movement (s) : "
            write(unit, '(A)')&
       "    Start Time of the Control Rod Movement (s) : "
            write(unit, ' (8x,  10E10.2)') &
         (CRD_Time_Start(i), i=1, NROD_Moving)

            write(unit, '(A)')&
       "    End Time of the Control Rod Movement (s) : "
            write(unit, ' (8x,  10E10.2)') &
         (CRD_Time_End(i), i=1, NROD_Moving)

! "    Velocity of the Moving Control Rods (cm/s)         : "
      if(CRD_Mov_Type.EQ."LIN") then 
      write(unit, '(A)')&
       "    Velocity of the Moving Control Rods (cm/s)         : "
      write(unit, ' (8x, 10E12.5)') &
         (v_rods(i), i=1, NROD_Moving)
! "    Amplitude of the Moving Control Rods (cm)         : "
      else if(CRD_Mov_Type.EQ."SIN") then 
      write(unit, '(A)')&
       "    Amplitude of the Moving Control Rods (cm)          : "
      write(unit, ' (8x, E12.5)') &
         (CRD_Sin_Ampl(i), i= 1, NROD_MOVING)
! "    Frequency of the Moving Control Rods (cm)         : "
      write(unit, '(A)')&
       "    Frequency of the Moving Control Rods (cm)          : "
      write(unit, ' (8x, E12.5)') &
         (CRD_Sin_Freq(i), i= 1, NROD_MOVING)
      end if ! "SIN".OR."LIN"         

!"    Top and Bottom Positions the Control Rods (cm)    : "
      write(unit, '(A)')&
       "    Top and Bottom Positions the Control Rods (cm)     : "
      write(unit, ' (8x, 2E12.5)') &
         Z_Rod_Top, Z_Rod_Bottom

!  SCRAM DATA   
      IF( Flag_Scram ) THEN
      WRITE(unit, '(A)')&
       "    Scram is  ON "
      write(unit, '(A)')&
       "    Velocity of the Scram Control Rod Movement (cm/s)  : "
      write(unit, ' (8x, E12.5)') &
         v_rods_scram

      write(unit, '(A)')&
       "    Number of Control Rods in Scram                    : "
      write(unit, ' (8x, I4)') &
         Nrod_Scram

      write(unit, '(A)')&
       "    Indexes of the Scram Control Rods                  : "
      write(unit, ' (8x, 10I4)') &
         (Num_Rod_Scram(i), i= 1, Nrod_Scram)

      write(unit, '(A)')&
       "    Power Level of the Scram (MWt)                     : "
      write(unit, ' (8x, E12.5)') &
         Power_Scram

      write(unit, '(A)')&
       "    Time Delay in Scram (s)                            : "
      write(unit, ' (8x, E12.5)') &
         Time_Scram_Delay

      ELSE
      WRITE(unit, '(A)')&
       "    Scram is  OFF "
      END IF

      END IF ! IF(Problem_Type.EQ."Kinetics") THEN

      call OUTput_Write_Separator(io_unit)
 
      RETURN
      END


      SUBROUTINE CRD_Output_CR_Positions(unit)
      implicit none
      include 'sketch.fh'
      integer unit

      integer ind, i, nb
      real CR_Pos_2D(N_POLY), CR_Number_1D(NZR)

      character*80 Header_Map
      character*4 val_fmt
      character*5 val_char(0:N_POLY)

! "    Control Rod Positions (cm, PWR from the bottom, BWR from the top):"
      Header_Map = &
       "    Control " //&
       "Rod Positions (cm, PWR from the bottom, BWR from the top):"
      val_fmt = "A5"
      ind = 0
      DO ind = 1, N_POLY 
         val_char(ind) = "   -"
      END DO
      DO i=1,NN_CRoD
        DO nb = 1, NN_CRod_Bundle
           write(val_char(nrods(i,nb)), '(F5.0)' ) zrods(i)
        END DO
      END DO
      call OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)


! "    Control Rod Positions (% of the Covered Reactor Core)           :"

      CALL CRD_Compute_CRD_Position_in_Percents(CR_Pos_2D, CR_Number_1D)

      Header_Map = &
       "    Control " //&
       "Rod Positions (% of the Covered Reactor Core)           :"

      val_fmt = "A5"
      ind = 0
      DO ind = 1, N_POLY 
         val_char(ind) = "   -"
      END DO
      DO i=1,NN_CRoD
        DO nb = 1, NN_CRod_Bundle
           write(val_char(nrods(i,nb)), '(F5.0)' ) &
                                           CR_Pos_2D(nrods(i,nb))
        END DO
      END DO
      call OUT_Write_Map(N_POLY,  NXR_B_Min_Core, &
        NXR_Max_Core, NXR_B_Core, NXR_E_Core, NYR_B_Core, NYR_E_Core,&
        index_core, Header_Map, unit, val_char, val_fmt)
      

! " Number of the CR in the Axial Layer:"
      write(io_unit,*) " Number of the CR in the Axial Layer:"
      WRITE(unit, '(8x, 10F7.3)') &
       (CR_Number_1D(ind), ind = 1, NZR)

      return
      end

      SUBROUTINE CRD_Fuel_Control_Rods
      implicit none
      include 'sketch.fh'

      INTEGER ir, ib, nch, n1, k, kt, np, m
      REAL    value_dt(NZ), value(NZ), z_new, z_old

      DO ir = 1,NN_CRod
! recalculate the material composition if control rod has moved
!         if(Mov_Rods(ir)) then
         z_new = zrods(ir)
         z_old = zrods_dt(ir)
         IF( ABS( z_new - z_old).GT. EPS_ROUND_OFF ) THEN
!         write(*,*) 'ir =', ir, 'z_new =', z_new, 'z_old =', z_old
!         pause 
           DO ib = 1, NN_CRod_Bundle
            np = nrods(ir, ib)
            do nch = 1, NCHM
               k = poly_out(np,nch)
               If(k.NE.0) then
! ATTENTION FOR BWR PROBLEM NO TREATMENT OF THE CRDS IN THE ALL AXIAL REFLECTORS
!           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
!              kt = k + (n1-1)*NH
! moving the values of the source term
!              value_dt(n1) = source(kt) 
!           END DO
!           CALL CRD_move_precursors(NZ, hzt, value_dt, value, z_new, &
!                z_old) 
!           do n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
!              kt = k + (n1-1)*NH
!              source(kt) = value(n1)  
!           END DO
           DO m = 1, MD
              DO n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
                kt = k + (n1-1)*NH
                value_dt(n1) = prec(m, kt) 
!                write(*,*) 'sf =', XS_SF(1, kt), XS_SF(2, kt)
!                write(*,*) 'flux =', flux(1, kt), flux(2, kt)
!                write(*,*) 'source =', source(kt)
              END DO  
!             IF(m.eq.1) THEN
!               write(*,*) 'np, k =', np, k
!             WRITE(*, *) 'Precursors old'
!             WRITE(*, '(5E13.5)') (value_dt(n1), n1=1,NZ/2)
!             END IF

               CALL CRD_move_precursors(NZ, hzt, value_dt, value, z_new, &
                z_old) 
!             IF(m.eq.1) THEN
!             WRITE(*, *) 'Precursors new'
!             WRITE(*, '(5E12.5)') (value(n1), n1=1,NZ/2)
!             pause
!             END IF
              DO n1 = NZ_Core_BEG, NZ_CORE_END ! NZ_Core_END (adding the upper reflector)
                 kt = k + (n1-1)*NH
                 prec(m, kt) = value(n1)  
              END DO ! n1
           END DO ! MD
          END IF ! k/=0
          END DO ! NCHM
        END DO !ib = 1, NN_CRod_Bundle
        END IF ! Mov_Rods 
      END DO ! NN_CRRODS

      RETURN
      END             

      SUBROUTINE CRD_find_axial_layer(NZ, hzt, z_pos, h1, n1)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NZ
      REAL,    INTENT(IN)    :: z_pos, hzt(NZ)

      REAL,    INTENT(OUT)   :: h1
      INTEGER, INTENT(OUT)   :: n1

      h1 = 0

      DO n1 = 1, NZ
         h1 = h1 + hzt(n1)
         IF( h1.GE.z_pos) EXIT
      END DO
      
      RETURN
      END           

      SUBROUTINE CRD_compute_average(NZ, hzt, z1, h1, n1, z0, n0,&
       c, c_av )
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NZ, n1, n0
      REAL,    INTENT(IN)    :: hzt(NZ), z1,  z0,  c(NZ)

      REAL,    INTENT(OUT)   :: c_av

      REAL, INTENT(INOUT) :: h1

      REAL :: dz
      INTEGER :: n

      IF(n1.NE.n0) THEN
        dz = z1 - ( h1-hzt(n1) )        
      ELSE
        dz = z1 - z0
      END IF               

       c_av = c(n1)*dz/hzt(n1)
!      WRITE(*,*) 'C(N1) =', C(n1)
!      WRITE(*,*) 'DZ =', DZ
!      WRITE(*,*) 'HZT(N1) =', HZT(N1)
!      WRITE(*,*) 'C_AV =', C_AV
!      PAUSE


      IF( n1.ne.n0 ) THEN  
         h1 = h1 - hzt(n1)
         DO n = n1 - 1, n0+1, -1
            dz = hzt(n)
            c_av = c_av + c(n)*dz/HZT(N)
            h1 = h1 - hzt(n)
         END DO
         dz = h1 - z0
         c_av = c_av + c(n0)*dz/HZT(N0)
      END IF

      RETURN
      END

      SUBROUTINE CRD_move_precursors(NZ, hzt, c, c_new, z, z_dt)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NZ
      REAL, INTENT(IN)    :: hzt(NZ), c(NZ), z, z_dt
             
      REAL, INTENT(OUT)    :: c_new(NZ)

      REAL :: z1_new, z1_old, z0_old,  h1_new, h1_old, &
         h0_old, dz, c_av
      INTEGER :: n1_new, n1_old, n0_old, n
      REAL, PARAMETER  :: EPS_ROUND_OFF = 1.E-05

      z1_new = z

      CALL CRD_find_axial_layer(NZ, hzt, z1_new, h1_new, n1_new)

! Setting to zero the precursors at the absorber (top of the core)
      DO n = NZ, n1_new+1, -1
         c_new(n) = 0.
      END DO
      
      z1_old = z_dt

      CALL CRD_find_axial_layer(NZ, hzt, z1_old, h1_old, n1_old)

      dz = z1_new - ( h1_new - hzt(n1_new) )
      z0_old = amax1( z1_old - dz, 0.)

      CALL CRD_find_axial_layer(NZ, hzt, z0_old, h0_old, n0_old)

      CALL CRD_compute_average(NZ, hzt, z1_old, h1_old, n1_old, &
        z0_old, n0_old, c, c_av )

!      dz = ( z1_new - (h1_new - hzt(n1_new) ) )

      IF( dz .GT. EPS_ROUND_OFF ) THEN  
        c_new(n1_new) = c_av ! hzt(n1_new)
      ELSE
        c_new(n1_new) = 0.
      END IF
     
      z1_old = z0_old 
      n1_old = n0_old

      DO n = n1_new - 1, 1, -1

         dz = hzt(n)
            z0_old = amax1( z1_old - dz, 0.)
         CALL CRD_find_axial_layer(NZ, hzt, z0_old, h0_old, n0_old)

         CALL CRD_compute_average(NZ, hzt, z1_old, h1_old, n1_old, &
            z0_old,  n0_old, c, c_av )
         c_new(n) = c_av !dz

         z1_old = z0_old 
         n1_old = n0_old

      END DO
     
      RETURN
      END 
    
