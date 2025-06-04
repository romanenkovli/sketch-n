      MODULE CR_POSITION


        INTEGER, PARAMETER   :: NN_CR_BANK=100, NN_CR_in_Bank=100
        INTEGER              :: N_CR_BANK, N_CR_TOTAL
        INTEGER              :: N_CR_in_Bank(NN_CR_BANK) 
          INTEGER              :: N_CR_POS(NN_CR_in_Bank,NN_CR_BANK) 
        INTEGER              :: N_CR_TYPE_in_BANK(NN_CR_BANK)  
        REAL                 :: BANK_POS(NN_CR_BANK), &
                               CR_POS(NN_CR_in_Bank,NN_CR_BANK)
        REAL                 :: ONE_PERCENT_CORE
!                               H_AXIAL_REFLECTOR=35.5            
        INTEGER              :: Index_CR_IN_Bank( NN_CR_in_Bank, &
                                                       NN_CR_BANK )
        LOGICAL              :: Input_Bank_Position, &
        Bank_Position_in_Percents, Bank_Position_Absolute
          REAL,PARAMETER  :: HZ_NK = 0.0! 17.25

      END MODULE CR_POSITION
         
      SUBROUTINE SET_CR_POSITION(i_bank, position )
      USE CR_POSITION
      IMPLICIT NONE
      include 'sketch.fh'

      INTEGER, INTENT(in) :: i_bank
      REAL,    INTENT(in) :: position

      INTEGER :: n, i_cr

!         HZ_NK = 0.0 ! 17.25   
!        DO n = 1,  
        CR_POS(1:N_CR_in_Bank(i_bank) ,i_bank)=HZ_AXIAL_REFLECTOR(1)+&
           position + HZ_NK

       DO n = 1, N_CR_in_Bank(i_bank)
          i_cr = Index_CR_IN_Bank( n,  i_bank )
          Mov_Rods(i_cr) = .True.
          zrods(i_cr) = CR_POS(n,i_bank)
       END DO
          
      RETURN
      END SUBROUTINE SET_CR_POSITION

      SUBROUTINE SET_ALL_CR_POSITION_ABS
      USE CR_POSITION
      IMPLICIT NONE
      include 'sketch.fh'

      INTEGER :: i_bank
!      REAL,    INTENT(in) :: position
!      REAL                :: HZ_NK 


      INTEGER :: n, i_cr

!        DO n = 1,  
!         HZ_NK = 0.0 !17.25   


       DO i_bank = 1, N_CR_BANK
        CR_POS(1:N_CR_in_Bank(i_bank) ,i_bank)=HZ_AXIAL_REFLECTOR(1)+&
           BANK_POS(i_bank) + HZ_NK

        DO n = 1, N_CR_in_Bank(i_bank)
          i_cr = Index_CR_IN_Bank( n,  i_bank )
          Mov_Rods(i_cr) = .True.
          zrods(i_cr) = CR_POS(n,i_bank)
        END DO
       END DO
          
      RETURN
      END SUBROUTINE SET_ALL_CR_POSITION_ABS

      SUBROUTINE SET_ALL_CR_POSITION_PERC
      USE CR_POSITION
      IMPLICIT NONE
      include 'sketch.fh'
!      REAL                :: HZ_NK 


      INTEGER :: i_bank
!      REAL,    INTENT(in) :: position

      INTEGER :: n, i_cr

       ONE_PERCENT_CORE =  hz_core / 100.

!       HZ_NK = 0.0 !  17.25   

       DO i_bank = 1, N_CR_BANK
        CR_POS(1:N_CR_in_Bank(i_bank) ,i_bank)=HZ_AXIAL_REFLECTOR(1)+&
           (100.-BANK_POS(i_bank))*ONE_PERCENT_CORE + HZ_NK

        DO n = 1, N_CR_in_Bank(i_bank)
          i_cr = Index_CR_IN_Bank( n,  i_bank )
          Mov_Rods(i_cr) = .True.
          zrods(i_cr) = CR_POS(n,i_bank)
        END DO
       END DO
          
      RETURN
      END SUBROUTINE SET_ALL_CR_POSITION_PERC

      SUBROUTINE SET_INDEX_CR_in_BANK
      USE CR_POSITION
      IMPLICIT NONE

      INTEGER :: n, i_cr, nn

        i_cr = 0
        DO  n = 1, N_CR_BANK
           DO  nn = 1, N_CR_in_Bank(n)
             i_cr = i_cr + 1
               Index_CR_IN_Bank( nn, n ) = i_cr
           END DO                 
        END DO

      RETURN
      END SUBROUTINE SET_INDEX_CR_in_BANK   


      subroutine set_cr_bank_pos_kln3(time, dt)
        USE CR_POSITION
      REAL time, dt
      real :: time_start_10, time_end_10
      real, parameter :: hz_cr_core = 376.
        real ::  v_suz      

        
!          time = 0.
!          dt =   0.1
!          time_end = 100.
        time_end_10 = 71.3
        time_start_10 = 1.4
!          h_10_0 = 302.
!          h_9_0  = 362.
        v_suz = -2.
!        h_9 = h_9_0
!        h_10 = h_10_0

           
!      DO WHILE ( time < TIME_END)        
!        time = time + dt
        IF( time .GE. time_start_10 .and. time .LE. time_end_10 ) then       
        CALL MOVE_SUZ_POS( BANK_POS(10), dt, v_suz) 
          END IF
        IF( BANK_POS(10) < hz_cr_core/2 .and. time .LE. time_end_10  ) then
         CALL MOVE_SUZ_POS( BANK_POS(9), dt, v_suz) 
        END IF

        write(*,'(A, 3F8.2)') 'time, h_9, h_10 =', &
         time, BANK_POS(9), BANK_POS(10)

!       END DO
      
      RETURN
      END

      SUBROUTINE MOVE_SUZ_POS( h, dt, v_suz) 
      real, intent(in):: dt, v_suz
      real, intent(inout) :: h

      h = h + v_suz*dt
      
      return
      end