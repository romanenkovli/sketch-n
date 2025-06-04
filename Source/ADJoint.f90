      subroutine ADJ_Compute
!=====================================================================*
! Iteration Procedure for the Steady_State Adjojnt Calculations       *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*

      implicit none
      include 'sketch.fh'
! Input, NG_BEG, NG_END, NG_Step, n_inter(NG), e_inter      
      integer i_out, i_int
! Output: i_source, i_ssor, i_nonl, i_sk_end
      integer i_source, i_ssor, i_nonl, i_sk_end
! Local Variables
      real dt_kin
      logical Conv_Outers, ADJOINT
!
      i_ssor = 0
      i_source = 0
      i_nonl = 0
      ADJOINT = .TRUE.


! Compute Block Diagonal

         call EIS_Update_Wieland_Shift

         call MAT_Compute_Block_Diag(ADJOINT, dt_kin)
           if(NG.eq. 2) then
! analytical inversion in the case of 2 neutron energy groups
              call MAT_Inverse_2x2_Block_Diag
           else
              call MAT_Inverse_Block_Diag
           end if

          call MAT_Transpose_Matrix


          N_OUTER = 1000


       call CTRL_Block_Iterations(i_out, i_int, &
                                            Conv_Outers, ADJOINT)


      i_ssor = i_ssor + i_int
      i_source = i_source + i_out


      if(Conv_Outers) i_sk_end = 1

       write(*,'(" Adjoint Calculations,  Number of Iterations = ", I3 &
         )') i_source
       write(*, '(" k eff =",&
         F8.5, " k_ef_max =", F8.5, " k_ef_min =", F8.5)')&
          k_ef, k_ef_max, k_ef_min 

      IF(i_sk_end.ne.1) THEN
         write(*,*) 'COULD NOT COMPUTE ADJOINT'
      ELSE
            call ADJoint_Normalize
      END IF

      return
      end



      subroutine ADJoint_Normalize
!**********************************************************************
!       Input of the Neutron Flux data from the Restart File, 
!        Normalize the Adjoint Flux and Output into the File FILE_A   
!**********************************************************************
      implicit none
      include 'sketch.fh'


! Input into the File FILE_D:  Flux(NG, 0:N_TOT) - Neutron Flux
      real time
! + Output into The File FILE_A:  Flux(NG, 0:N_TOT) - Normalized Adjoint
      integer n, k

      do k = 1, N_TOT
         do n = 1, NG
            Flux_A(n, k) = Flux(n,k)
         end do
      end do

      if (FILE_DMP_OUT_KIN.NE."") then
        open(io_unit,file=FILE_DMP_OUT_KIN,&
         form = 'UNFORMATTED',status='OLD')
         CALL INPut_read_restart_steady_state_data(io_unit)
         CALL THM_read_data_restart_file(io_unit)
         CALL INput_read_restart_burnup_data(io_unit)
         CALL INPut_read_restart_kinetics_data(io_unit, time)
       close(io_unit)
      else
         write(*,*) 'Adjoint Neutron Flux is not written in the FILE'
         write(*,*) 'Please,  specify the FILE_DMP_OUT_KIN &
           in the NEUTRON.INI'
         RETURN
      end if

! Setting computed adjoint

      call ADJ_Normalize_Flux(NG, N_TOT,  xs_al, Flux(1,1), &
                 Flux_A(1,1), Volume (1) )

      open(io_unit,file = FILE_DMP_OUT_KIN,form = 'UNFORMATTED',&
         status='OLD')

         CALL OUTput_write_restart_steady_state_data(io_unit)
         CALL THM_write_data_restart_file(io_unit)
         CALL OUTput_write_restart_burnup_data(io_unit)
         CALL OUTput_write_restart_kinetics_data(io_unit, time)
         CALL OUTput_write_restart_adjoint_flux(io_unit)

         write(*,*) 'Flux_a(:,:)'
           write(*,*) Flux_a(:,:)
           read(*,*)

      close(io_unit)
 
      return
      end

      subroutine ADJ_Normalize_Flux(NG, N_TOT, xs_al, Flux, &
                  Flux_a, Volume)
      implicit none 
!=====================================================================*
!  Normalization of the Flux or Adjoint from the Condition            *
!                 P_Reactor = <Flux, al*Flux_a>                       *
!=====================================================================*


!Input:  
      integer NG, N_TOT
      real xs_al(NG, N_TOT), volume(N_TOT)
      real Flux(NG, N_TOT), Flux_a(NG, N_TOT)
!      real P_Reactor
! Output:
!      real FLux_a(NG, N_TOT)
! Local Variables
      real adjoint_norm
      integer k, n

!     write(*,*) 'Flux_a(1,1), Flux_a(1,2) =', Flux_a(1,1),Flux_a(1,2)


      adjoint_norm = 0. 
      do k = 1, N_TOT
         do n = 1, NG
            adjoint_norm=adjoint_norm+xs_al(n,k)*Flux(n,k)*Flux_a(n,k)*&
                   Volume(k)
         end do
      end do

      adjoint_norm = 1. / Adjoint_Norm ! P_reactor

 
      call MSC_SSCALE(N_TOT*NG, adjoint_norm, Flux_a)

      return
      end


