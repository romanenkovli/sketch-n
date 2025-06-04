	MODULE TH_SPLINES
	IMPLICIT NONE
	INTEGER, PARAMETER :: NP=278, NQ=163
	REAL :: x_p_enter(NP),y_p_enter(NP,NQ),y2_p_enter(NP,NQ)
	REAL :: x_p_exit(NP),y_p_exit(NP,NQ),y2_p_exit(NP,NQ)
	REAL :: x_t_enter(NP),y_t_enter(NP,NQ),y2_t_enter(NP,NQ)
	REAL :: x_g_enter(NP),y_g_enter(NP,NQ),y2_g_enter(NP,NQ)


	
	CONTAINS

	  SUBROUTINE TH_SPLINE_READ
	  INTEGER :: i,j 
         

 	   open (unit=10, file='Input\ASS_P_ENTER.spl', action='read')
           do i=1,NP
            read (unit=10, FMT=*) x_p_enter(i), 
     &                (y_p_enter(i,j), y2_p_enter(i,j),j=1,NQ)
           END DO
           close(unit=10)
	 

 	   open (unit=10, file='Input\ASS_T_ENTER.spl', action='read')
           do i=1,NP
            read (unit=10, FMT=*) x_t_enter(i),    
     &                 (y_t_enter(i,j), y2_t_enter(i,j),j=1,NQ)
           END DO
           close(unit=10)

 	   open (unit=10, file='Input\ASS_G_ENTER.spl', action='read')
           do i=1,NP
            read (unit=10, FMT=*) x_g_enter(i),   
     &                 (y_g_enter(i,j), y2_g_enter(i,j),j=1,NQ)
           END DO
           close(unit=10)

 	   open (unit=10, file='Input\ASS_P_EXIT.spl', action='read')
           do i=1,NP
            read (unit=10, FMT=*) x_p_exit(i),  
     &                 (y_p_exit(i,j), y2_p_exit(i,j),j=1,NQ)
           END DO
           close(unit=10)

	   RETURN
           END SUBROUTINE TH_SPLINE_READ

           SUBROUTINE SPLINE_COMPUTE(x_time,p_enter, p_exit, t_enter, 
     &             g_enter)

	   REAL :: x_time, p_enter(1:NQ), p_exit(1:NQ), t_enter(1:NQ),  
     &             g_enter(1:NQ)

           CALL splint_nq
     &       (x_p_enter,y_p_enter,y2_p_enter,NP,x_time,p_enter,NQ)
           CALL splint_nq
     &       (x_g_enter,y_g_enter,y2_g_enter,NP,x_time,g_enter,NQ)
           CALL splint_nq
     &       (x_t_enter,y_t_enter,y2_t_enter,NP,x_time,t_enter,NQ)
           CALL splint_nq
     &       (x_p_exit,y_p_exit,y2_p_exit,NP,x_time,p_exit,NQ)

           RETURN
	   END SUBROUTINE SPLINE_COMPUTE

           END MODULE TH_SPLINES	
	

      SUBROUTINE splint_nq(xa,ya,y2a,n,x,y,NQ)
      INTEGER n
      REAL x,y(NQ),xa(n),y2a(n,NQ),ya(n,NQ)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y(1:NQ)=a*ya(klo,1:NQ)+b*ya(khi,1:NQ)+((a**3-a)*y2a(klo,1:NQ)+ 
     &          (b**3-b)*y2a(khi,1:NQ))*(h**2)/6.
      return
      END


	PROGRAM MAIN_TH_SPLINE
	USE TH_SPLINES
        IMPLICIT NONE

	INTEGER, PARAMETER :: NQ1 = 163
        REAL :: x_time, p_enter(1:NQ1), p_exit(1:NQ1), t_enter(1:NQ1),  
     &             g_enter(1:NQ1)

        real*8 T_PROP, P_PROP, RO_PROP, CP_PROP, Viscos_PROP,
     &        Conduct_PROP, diffus_PROP
        real y_g_enter_total, y_flow_rate

        INTEGER, PARAMETER :: N_GRAPH=14
        INTEGER N_OUT(N_GRAPH), n, k
	DATA N_OUT /12, 13, 28, 45, 54, 71, 77, 80, 89, 98, 106, 112,
     &   131, 154/ 

        external PF1STEAM

	
        CALL TH_SPLINE_READ

!        x_time = 9.29333E+00
!	CALL  SPLINE_COMPUTE(x_time,p_enter, p_exit, t_enter, g_enter)


	OPEN(2, FILE ='p_enter_graph.dat', action = 'WRITE')
         DO n = 1, NP
	 WRITE(2,'(ES12.5,163ES14.5)') x_p_enter(n), 
     & (y_p_enter(n,n_out(k)),k=1,N_GRAPH)
        END DO 
        CLOSE(2)

	OPEN(2, FILE ='t_enter_graph.dat', action = 'WRITE')
         DO n = 1, NP
	 WRITE(2,'(ES12.5,163ES14.5)') x_t_enter(n), 
     & (y_t_enter(n,n_out(k)),k=1,N_GRAPH)
        END DO 
        CLOSE(2)

	OPEN(2, FILE ='g_enter_graph.dat', action = 'WRITE')
         DO n = 1, NP
	 WRITE(2,'(ES12.5,163ES14.5)') x_g_enter(n), 
     & (y_g_enter(n,n_out(k)),k=1,N_GRAPH)
        END DO 
        CLOSE(2)

	OPEN(2, FILE ='p_exit_graph.dat', action = 'WRITE')
         DO n = 1, NP
	 WRITE(2,'(ES12.5,163ES14.5)') x_p_exit(n), 
     & (y_p_exit(n,n_out(k)),k=1,N_GRAPH)
        END DO 
        CLOSE(2)

	OPEN(2, FILE ='g_total_graph.dat', action = 'WRITE')
         DO n = 1, NP
	 WRITE(2,'(ES12.5,163ES14.5)') x_g_enter(n), 
     & (SUM(y_g_enter(n,1:NQ1))*3600./1000.)
        END DO 
        CLOSE(2)



	OPEN(2, FILE ='flow_rate_graph.dat', action = 'WRITE')
         DO n = 1, NP
          T_PROP = SUM(y_t_enter(n,1:NQ1))/NQ1
          P_PROP = SUM(y_P_enter(n,1:NQ1))/NQ1 + 273.15
          CALL PF1STEAM( P_PROP, T_PROP, RO_PROP, CP_PROP, 
     &         Viscos_PROP, Conduct_PROP)
          y_g_enter_total =  SUM(y_g_enter(n,1:NQ1))
          y_flow_rate = y_g_enter_total*3600./ RO_PROP
	 WRITE(2,'(ES12.5,163ES14.5)') x_g_enter(n), 
     &  y_flow_rate
        END DO 
        CLOSE(2)



        STOP
	END
	



