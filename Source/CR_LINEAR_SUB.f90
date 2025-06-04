      MODULE CR_LINEAR
      IMPLICIT NONE
      INTEGER, PARAMETER :: NP=401, NQ=2
      REAL :: time(NP),x_grp(NP,NQ)
      
      
      CONTAINS

        SUBROUTINE CR_LINEAR_READ
        INTEGER :: i,j 

!       write(*,*) 'CR_LINEAR_READ ='
!      pause
         

          open (unit=1, file='Input\CR_POS_ABS.dat', action='read')
           do i=1,NP
            read (unit=1, FMT=*) time(i), &
                     (x_grp(i,j),j=1,NQ)
!            write(*,*) time(i)
           END DO
           close(unit=1)
                               

         RETURN
           END SUBROUTINE CR_LINEAR_READ

           SUBROUTINE LINEAR_COMPUTE(x_time, CR_POS)

         REAL :: x_time, CR_POS(1:NQ)

           CALL linear_interpolation&
            (time, x_grp, NP, x_time, &
                 CR_POS, NQ)
           
           RETURN
         END SUBROUTINE LINEAR_COMPUTE

           END MODULE CR_LINEAR      
      

      SUBROUTINE linear_interpolation(xa,ya,n,x,y,NQ)
      INTEGER n
      REAL x,y(NQ),xa(n),ya(n,NQ)
      INTEGER k,khi,klo
      REAL a,b,h
!      write(*,*) (xa(i), i=1,3)
!      pause
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      write(*,*) klo, khi, h, xa(klo), xa(khi)
      if (h.eq.0.) then
        write(*,*) 'bad xa input in splint'
        read(*,*)
      endif
!      a=(xa(khi)-x)/h
!      b=(x-xa(klo))/h
!      y(1:NQ)=a*ya(klo,1:NQ)+b*ya(khi,1:NQ)+((a**3-a)*y2a(klo,1:NQ)+ &
!               (b**3-b)*y2a(khi,1:NQ))*(h**2)/6.
      y(1:NQ)=ya(klo,1:NQ)*(1-(x-xa(klo))/h)+ya(khi,1:NQ)*&
               ((x-xa(klo))/h)
      return
      END


      SUBROUTINE CR_LINEAR_KLN3(x_time)
      USE CR_LINEAR
        USE CR_POSITION

      INTEGER, PARAMETER :: NQ1 = 2
        REAL :: x_time, CR_POS_LINEAR(1:NQ1)  
      
        CALL CR_LINEAR_READ

!        x_time = 1.5      

      CALL  LINEAR_COMPUTE(x_time, CR_POS_LINEAR(1:NQ1)  )
             BANK_POS(9) = CR_POS_LINEAR(1)
             BANK_POS(10) = CR_POS_LINEAR(2)

!      OPEN(2, FILE ='Result_1s.dat', action = 'WRITE')
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,p_enter(1:Nq1)
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,p_exit(1:Nq1)
!      WRITE(2,'(ES12.5,163ES14.5)') x_time, CR_POS(1:Nq1)
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,g_enter(1:Nq1)
!        CLOSE(2)

        RETURN
      END


!      PROGRAM MAIN_CR_LINEAR
!      USE CR_LINEAR

!      INTEGER, PARAMETER :: NQ1 = 2
!        REAL :: x_time, CR_POS(1:NQ1)  
      
!        CALL CR_LINEAR_READ

!        x_time = 1.5      
!      CALL  LINEAR_COMPUTE(x_time, CR_POS)

!      OPEN(2, FILE ='Result_1s.dat', action = 'WRITE')
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,p_enter(1:Nq1)
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,p_exit(1:Nq1)
!      WRITE(2,'(ES12.5,163ES14.5)') x_time, CR_POS(1:Nq1)
!       WRITE(2,'(ES12.5,163ES14.5)') x_time,g_enter(1:Nq1)
!        CLOSE(2)

!        STOP
!      END
      



