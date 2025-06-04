      SUBROUTINE CPU_write_date_and_time(UNIT)
!      USE DFPORT ! DEC NT
      INTEGER unit
      CHARACTER*1 mer
      CHARACTER*9 month(12)
      data month /"January", "February", "March",    "April",&
                 "May ",    "June",     "July",     "August", &
                 "September", "October","November", "December"/

      INTEGER*4 DMY(3), HMS(3)    

      CHARACTER*80 output_string


      CALL IDATE(DMY)
      CALL ITIME (HMS)

      IF ( hms(1) .GT. 12) THEN
         mer = 'p'
         hms(1) = hms(1) - 12
      ELSE
         mer = 'a'
      END IF

      WRITE(output_string, &
     '( A, I3, A10, I4, I4,":", I2.2, ":", I2.2, " ", A, "m", &
        14(" "), "!"  ) ' )&
         "!      Date & Time of Calculation :", &
         DMY(2), month(DMY(2)), DMY(3), &
         HMS(1), HMS(2), HMS(3),  mer

      WRITE(unit, '(A)' ) output_string
      

      RETURN
      END

      SUBROUTINE CPU_write_date_and_time_unformatted(UNIT)
!      USE DFPORT ! DEC NT
      INTEGER unit
      CHARACTER*1 mer
      CHARACTER*9 month(12)
      data month /"January", "February", "March",    "April",&
                 "May ",    "June",     "July",     "August", &
                 "September", "October","November", "December"/

      INTEGER*4 DMY(3), HMS(3)    

      CHARACTER*80 output_string


      CALL IDATE(DMY)
      CALL ITIME (HMS)

      IF ( hms(1) .GT. 12) THEN
         mer = 'p'
         hms(1) = hms(1) - 12
      ELSE
         mer = 'a'
      END IF

      WRITE(output_string, &
     '( A, I3, A10, I4, I4,":", I2.2, ":", I2.2, " ", A, "m", &
        14(" "), "!"  ) ' )&
         "!      Date & Time of Calculation :", &
         DMY(2), month(DMY(2)), DMY(3), &
         HMS(1), HMS(2), HMS(3),  mer

      WRITE(unit) output_string
      

      RETURN
      END



      subroutine CPU_timer( Time_cpu )
!      USE DFPORT ! DEC
      implicit none
! Input Variables:
      Real Time_CPU
! Output Variables        
!      Real  Time_CPU
! Local Variables
      real Time_comp, Time_Array(2)     
      real   DTIME ! UNIX

      Time_Comp =  DTIME( Time_Array )

      Time_cpu = Time_cpu + Time_Comp

      return
      end

      subroutine CPU_total_timer(Time_Tot, Time_cpu)
!      USE DFPORT ! DEC
      implicit none

! Output Variables        
      Real Time_tot(2), Time_CPU
      real   ETIME !UNIX

      Time_CPU =  ETIME( Time_Tot )

      return
      end

      subroutine CPU_Timer_Output
      implicit none
      include 'sketch.fh'

      call CPU_Timer( time_sk_cpu )


      open(io_unit, file = 'Output/SKETCH.lst', &
         status = 'unknown', access = 'append')
                     
      call OUTput_Write_Separator(io_unit)


      time_sk_cpu = time_sk_cpu + time_th_cpu + time_nonl_cpu +&
         time_cmfd_cpu + time_XS_CPU+time_output_CPU
!     + time_tr_cpu 

      write(io_unit,*) 'SKETCH has used', time_sk_cpu, &
          'seconds of CPU time'

      write(io_unit,*) '   This Include: '
    
      write(io_unit,*) '      Internal  T/H Module     ', &
          time_th_cpu, ' s ',&
          time_th_cpu*100./time_sk_cpu, ' %'

      write(io_unit,*) '      Nodal Coefficient Update:',&
          time_nonl_cpu, ' s ',&
          time_nonl_cpu*100./time_sk_cpu, ' %'

      write(io_unit,*) '      CMFD Iterations          ',&
        time_cmfd_cpu, ' s ',&
        time_cmfd_cpu*100./time_sk_cpu, ' %'

      write(io_unit,*) '      XS Updates               ',&
        time_XS_cpu, ' s ',&
        time_XS_cpu*100./time_sk_cpu, ' %'

      write(io_unit,*) '      Output Module     ', &
          time_output_cpu, ' s ',&
          time_output_cpu*100./time_sk_cpu, ' %'

!      write(io_unit,*) '      TRAC                     ', &
!        time_tr_cpu,  ' s ', &
!        time_tr_cpu*100./time_sk_cpu, ' %'

      call CPU_total_timer( Time_ch_tot, Time_ch_cpu)

      write(io_unit,*) 'CHECKING SUM', time_ch_cpu,&
         ' seconds of CPU time'
      write(io_unit,*) 'This include user time, system time = ', &
          time_ch_tot(1),&
          time_ch_tot(2)

      call OUTput_Write_Separator(io_unit)

      close(io_unit)
      
      return
      end




