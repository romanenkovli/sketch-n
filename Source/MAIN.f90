      PROGRAM SKETCH_N
!======================================================================c
!          SKETCH-N version 0.95: Nodal Neutron Diffusion Code for     c
!            Solving Steady-State & Kinetics Problems                  c
!                                                                      c
!             Moscow Engineering Physics Institute                     c
!                Tokyo Institute of Technology                         c
!              Japan Atomic Energy Research Institute                  c 
!                                                                      c
!       Author:  Vyacheslav G. Zimin                                   c
!                                                                      c
!                    (C) 1999 All Rights Reserved                      c
!                                                                      c
!                               NOTICE                                 c
!                                                                      c
!  Permission to use, copy, modify, and distribute this software and   c
!  its documentation for any purpose and without fee is hereby granted c
!  provided that the above copyright notice appear in all copies and   c
!  that both the copyright notice and this permission notice appear in c
!  supporting documentation.                                           c
!                                                                      c
!  Neither the Institutions  nor the Authors make any                  c
!  representations about the suitability of this software for any      c
!  purpose.  This software is provided ``as is'' without express or    c
!  implied warranty.                                                   c
!======================================================================c
      USE MAIN_VAR 
      IMPLICIT NONE
      include 'sketch.fh'

      REAL :: dt_sim

! Initial 

      CALL sketch_init


!     MAIN COMPUTING PROCEDURE

      do while(i_trac_end.ne.1)

            dt_sim = dt_sketch

            CALL sketch_compute_time_step(dt_sim)

      END DO


!     END MAIN COMPUTING PROCEDURE
      write(*,*) 'CALL sketch_end_remarks'
!      read(*,*)

      CALL sketch_end_remarks


     

!      dt_sim = 0.1
!      CALL  PNT_Compute_Reactivity(dt_sim)


      STOP
      END

