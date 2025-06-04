      subroutine PVM_SKETCH_Enroll(tid_parent)
!=====================================================================*
!       Enroll SKETCH into PVM and Get Parents ID                     *
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'fpvm3.h'        
! Input: -
! Output: tid_parent - TID of TRAC Process
      integer tid_parent
! Local Variables: mytid - TID of SKETCH

!PVM      integer mytid

! get SKETCH TID
!PVM  call pvmfmytid(mytid)
! get TRAC TID
!PVM  call pvmfparent(tid_parent)

      return
      end


      subroutine PVM_SKETCH_SEND(tid_parent,i_sk_end, dt_sketch_send )
!=====================================================================*
!      Converting the power density distribution for TRAC nodes       *
!            Volume-weighting of the power densities                  *
!       &  Sending Message to the TRAC                                *
! (c) Slava 4.III.1998 JAERI                                          *
! Last modification 2.06.98 Sending Power in Wt                       *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      include 'fpvm3.h'
! Input: tid_parent - TID of TRAC Process
!        i_sk_end - SKETCH Convergence (1 if convergent)
!        dt_sketch -  time step size which was used 
!        p_col(N_POLY, NZR) - power density
!        vol_ass(N_POLY, NZR)
!        ia_1D_NTHC(N_Z_HC_POW), ja_1D_NTHC(NON_ZERO_Z), 
!        a_Z_HC_POW(none_ZERO_Z) - Mapping Array for Power - Z Direction
!        ia_2D_NTHC(N_RT_HC_POW), ja_2D_NTHC(NON_RTERO_RT), 
!        a_RT_HC_POW(none_RTERO_RT) - Mapping Array for Power - RT Direction
      integer i_sk_end, tid_parent
      real*8 dt_sketch_send 
!      real dt_sketch
!     msgtype - Meesage identifier
      integer msgtype
      parameter  (msgtype = 2)
! Output: p_trac(N_RT_HC_POW, N_RT_HC_POW) - power density distribution 
!                                                  for TRAC
!*NEW       real*8 pow_trac_new(N_RT_HC_POW, N_Z_HC_POW+1) 
!   pow_trac_new    - Power on the Staggered Heat Conduction Mesh
          
! Local Variables:
      integer  k, n1, n1_trac, k_trac, kp, nzp

!PVM  integer info         

      real  p_sketch
!*NEW      , p_trac_new    ! Total Power CHECKING
        
 

      do n1_trac = 1, NN_Z_HC_TRAC
            do k_trac = 1, NN_RT_HC_TRAC
            p_hc_trac(k_trac, n1_trac) = 0.

             do nzp = ia_1D_NTHC(n1_trac), ia_1D_NTHC(n1_trac+1)-1
                n1 = ja_1D_NTHC(nzp)
                do kp=ia_2D_NTHC(k_trac), ia_2D_NTHC(k_trac+1)-1
                   k = ja_2D_NTHC(kp)
                   p_hc_trac(k_trac, n1_trac) = &
                        p_hc_trac(k_trac, n1_trac) +&
                    (1. - Cool_Heating)*vol_ass(k,n1)*p_col(k,n1)*&
                    Map_2D_NTHC(kp)*Map_1D_NTHC(nzp)
               end do
            end do
         end do          
       end do

!      write(*,*)  'Power HC Compute'


       do n1_trac = 1, NN_Z_FD_TRAC
            do k_trac = 1, NN_RT_FD_TRAC
!            write(*,*) 'k_trac, n1_trac =', k_trac, n1_trac
             p_FD_trac(k_trac, n1_trac) = 0.

             do nzp = ia_1D_NTFD(n1_trac), ia_1D_NTFD(n1_trac+1)-1
!                write(*,*) 'nzp =', nzp
                 n1 = ja_1D_NTFD(nzp)
!                write(*,*) 'n1 =', n1
                 do kp=ia_2D_NTFD(k_trac), ia_2D_NTFD(k_trac+1)-1
                    k = ja_2D_NTFD(kp)
!                   write(*,*)  'k,   kp', k,  kp
                    p_FD_trac(k_trac, n1_trac) = &
                        p_FD_trac(k_trac, n1_trac) +&
                    Cool_Heating*vol_ass(k,n1)*p_col(k,n1)*&
                    Map_2D_NTFD(kp)*Map_1D_NTFD(nzp)
               end do
            end do
         end do          
       end do

!      write(*,*)  'Power FD Compute'

! CHECKING POWER
         
          p_sketch = 0.
         do k = 1, NN_RT_HC_TRAC
         do n1 = 1, NN_Z_HC_TRAC
           p_sketch = p_sketch + p_hc_trac(k,n1)
          end do ! N_Z_HC_POW
         end do

        do k = 1, NN_RT_FD_TRAC
         do n1 = 1, NN_Z_FD_TRAC
          p_sketch = p_sketch + p_FD_trac(k,n1)
          end do ! N_Z_FD_POW
        end do

          write(*,*) 'Sketch Sending Power: ',&
               p_sketch/1.E+06

!          WRITE(*,'(A,12ES12.3)') "POW SKETCH :", p_col(1,1:NZR)
!        WRITE(*,'(A,12ES12.3)') "POW HC TRAC:",  p_hc_trac(1,:)
!        WRITE(*,'(A,12ES12.3)') "POW FD TRAC:",  p_FD_trac(1,:)

!          PAUSE 
! END CHECKING POWER

! OLD Variant
!         dt_sketch_send = dble(dt_sketch)

! initializing of the send buffer    
!PVM      call pvmfinitsend(PVMDEFAULT, info)
! packing data into the buffer      
!PVM      call pvmfpack(INTEGER4, i_sk_end, 1, 1, info)
!PVM      call pvmfpack(REAL8, p_hc_trac, NN_RT_HC_TRAC*NN_Z_HC_TRAC,
!PVM     &                                            1, info)
          if(TRAC_Version.eq.'TRAC-PF1') then
!PVM      call pvmfpack(REAL8, p_fd_trac, NN_RT_FD_TRAC*NN_Z_FD_TRAC,
!PVM     &                                            1, info)
          end if
!PVM      call pvmfpack(REAL8, dt_sketch_send, 1, 1, info)
           
! sending message to PARENT      
!PVM      call pvmfsend(tid_parent, msgtype, info)

      return
      end

      
      subroutine PVM_SKETCH_RECEIVE(tid_parent, i_sk_end, i_trac_end)
!=====================================================================*
!       Reseiving Message from the TRAC                               *
!      & Converting Feedbacks into SKETCH nodes                       *
! (c) Slava 4.III.1998 JAERI                                          *
! + Converting Units Water Density from kg/m^3 into g/cm^3            *
!                    Water Temperature from K into C                  *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      include 'fpvm3.h'
! INPUT: NN_RT_FD_TRAC = NTSX*NRSX , NN_Z_FD_TRAC = NASX 
!   NN_RT_HC_TRAC = NRODS,  NN_Z_HC_TRAC = NCRZ  - dimension of the TRAC 
!                                           variables
!    tf_trac(NN_RT_HC_TRAC, NN_Z_HC_TRAC) - fuel temperature
!    tc_trac(NN_RT_FD_TRAC, NN_Z_FD_TRAC) - coolant temperature
!    ro_c_trac(NN_RT_FD_TRAC, NN_Z_FD_TRAC) - coolant density
! msgtype - message identifier for the TRAC
      integer  msgtype
      parameter (msgtype = 1)
! tid_parent - TID of the TRAC process
      integer  i_sk_end, i_trac_end, tid_parent
! Output: fdback(NH,NZ,N_FEEDACK) - current values of the feedbacks
! 1 - Boron (Not Used) 2 - Coolant Tempearture, 3 - Coolant Density,
! 4 - Fuel Temperature
! dt_trac - time step size required for the TRAC code
!      real dt_trac
!PVM      real*8 dt_trac_receive      
! i_trac_end - Convergence of the TRAC Calculations for the
!              Steady-State case
! Local Variables:
      real Doppler_Factor
      parameter (Doppler_Factor = 0.7)

      integer  k, n1, n1_trac, k_trac, np,  kp, nzp, ns

!PVM   integer info

      Logical SHIFT
      parameter(SHIFT=.True.)
      
!      real Conv_Dens, Conv_Temp
!      parameter (Conv_Dens = 1.E-03, Conv_Temp = 273.15)

      real ro_tmp, ro2_tmp, temp2_tmp, temp_tmp

      Logical Zero_Map ! used to Reflector which is treated separately

!PVM       call pvmfrecv(tid_parent, msgtype, info)

!PVM       call pvmfunpack(INTEGER4, i_trac_end, 1, 1, info)
      if(TRAC_Version.eq.'TRAC-PF1') then
!PVM       call pvmfunpack(REAL8, tf_trac_cl, NN_RT_HC_TRAC*NN_Z_HC_TRAC,
!PVM     &                                                 1, info)
!PVM       call pvmfunpack(REAL8, tf_trac_sf, NN_RT_HC_TRAC*NN_Z_HC_TRAC,
!PVM     &                                                 1, info)
      else
!PVM       call pvmfunpack(REAL8, tf_trac, NN_RT_HC_TRAC*NN_Z_HC_TRAC,
!PVM     &                                                 1, info)
      end if
            
!PVM       call pvmfunpack(REAL8, tc_trac, NN_RT_FD_TRAC*NN_Z_FD_TRAC, 
!PVM     &                                                 1, info)
!PVM       call pvmfunpack(REAL8, ro_c_trac, NN_RT_FD_TRAC*NN_Z_FD_TRAC,
!PVM     &                                                 1, info)
!PVM       call pvmfunpack(REAL8, dt_trac_receive, 1, 1, info)

!converting the TRAC variables into the SKETCH feedbacks
! and real*8 type to real type

        dt_trac = 1.E+34
      i_trac_end = i_sk_end

      tf_trac(1:NN_RT_HC_TRAC,1:NN_Z_HC_TRAC) = 553.15
            
      tc_trac(1:NN_RT_FD_TRAC,1:NN_Z_FD_TRAC) = 553.15

      ro_c_trac(1:NN_RT_FD_TRAC,1:NN_Z_FD_TRAC) = 0.76538*1.E+3


! Maximum Fuel Centerline Temperature

      do n1 = 1, NZ
         do k = 1, NH

               np=np_out(k)
             ns=ns_out(n1)

               zero_map = (ia_1D_HCNT(ns).GT.(ia_1D_HCNT(ns+1)-1))&
              .OR. (ia_2D_HCNT(np).GT.(ia_2D_HCNT(np+1)-1))

               if (.NOT. Zero_Map) then

               fdback(k, n1, 4) = 0.

               do nzp = ia_1D_HCNT(ns), ia_1D_HCNT(ns+1)-1
                n1_trac = ja_1D_HCNT(nzp)
                do kp=ia_2D_HCNT(np), ia_2D_HCNT(np+1)-1
!                  write(*,*) 'n1=', n1, 'ns= ', ns, 'n1_HC=', n1_trac
                   k_trac = ja_2D_HCNT(kp)
                   fdback(k, n1, 4) = fdback(k, n1, 4) +&
                    (sngl(tf_trac(k_trac,n1_trac)) + Conv_Fuel_Temp)*&
                         Map_2D_HCNT(kp)*Map_1D_HCNT(nzp)
                end do
               end do
               end if

         end do    
       end do

!      pause 

! Converting the Coolant Temperature and Coolant Density
         if(TRAC_Version.eq."TRAC-PF1") then
! linear interpolation of the coolant temperature and coolant density to get the 
!   node-avergaed values
         do k = 1, NN_RT_FD_TRAC
              ro_tmp =  ro_c_trac(k,1)   
              temp_tmp =  tc_trac(k,1)   
              do n1 = 2, NN_Z_FD_TRAC

               ro2_tmp = ro_c_trac(k,n1)
               ro_c_trac(k,n1) = ro_tmp + &
                            0.5*(ro_c_trac(k,n1)-ro_tmp)
               ro_tmp = ro2_tmp     

               temp2_tmp = tc_trac(k,n1)
               tc_trac(k,n1) = temp_tmp + &
                            0.5*(tc_trac(k,n1)-temp_tmp)
               temp_tmp = temp2_tmp     

               end do
            end do
        end if  


       do n1 = 1, NZ
            do k = 1, NH
              
             ns=ns_out(n1) 
               np=np_out(k)


               zero_map = (ia_1D_FDNT(ns).GT.(ia_1D_FDNT(ns+1)-1))&
              .OR. (ia_2D_FDNT(np).GT.(ia_2D_FDNT(np+1)-1))
               
!               write(*,*) 'zero_map=', zero_map 

               if(.NOT. Zero_Map) then

               fdback(k, n1, 2) = 0.
               fdback(k, n1, 3) = 0.
               do nzp = ia_1D_FDNT(ns), ia_1D_FDNT(ns+1)-1
                  n1_trac = ja_1D_FDNT(nzp)
                  do kp=ia_2D_FDNT(np), ia_2D_FDNT(np+1)-1
!                  write(*,*) 'n1=', n1, 'ns= ', ns, 'n1_trac=', n1_trac
                   k_trac = ja_2D_FDNT(kp)
! Converting Kelvin into Celcius             
                   fdback(k, n1, 2) = fdback(k, n1, 2) +&
                    (sngl(tc_trac(k_trac,n1_trac)) + Conv_Cool_Temp)*&
                         Map_2D_FDNT(kp)*Map_1D_FDNT(nzp)
! Converting Kg/m^3 into g/cm^3                
                   fdback(k, n1, 3) = fdback(k, n1, 3) +&
                    (sngl(ro_c_trac(k_trac,n1_trac))*Conv_Cool_Dens)*&
                         Map_2D_FDNT(kp)*Map_1D_FDNT(nzp)
                  end do
                end do
               end if
         end do    
       end do

!          WRITE(*,'(A,12ES12.3,12ES12.3)') &
!                        "fdback(k, n1, 2)", fdback(1, :, 2)
!        WRITE(*,'(A,12ES12.3,12ES12.3)') &
!                       "fdback(k, n1, 3)", fdback(1, :, 3)
!        WRITE(*,'(A,12ES12.3,12ES12.3)') &
!                       "fdback(k, n1, 4)", fdback(1, :, 4)

       
      return
      end



      subroutine PVM_SKETCH_Exit(tid_parent)
!=====================================================================*
!              Exit from PVM                                          *
! (c) Slava 4.III.1998 JAERI                                          *
!=====================================================================*
      implicit none
      include 'fpvm3.h'
! Local Variables 
      integer info, tid_parent, msgtype
      parameter (msgtype = 3)
      character*15 Sketch_Status
      parameter (Sketch_Status = 'SKETCH FINISHED')
      
! PVM      call pvmfinitsend(PVMDEFAULT, info)
! packing data into the buffer      
! PVM      call pvmfpack(STRING,Sketch_Status, 15, 1, info)
! sending message to PARENT      
! PVM      call pvmfsend(tid_parent, msgtype, info)

! exit from PVM
! PVM      call pvmfexit(info)
      write(*,*) 'PVM EXIT DONE, info ='!, info

      return
      end

