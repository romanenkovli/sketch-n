      subroutine TRL_Update_TRL(dt_kin)
!=====================================================================*
! Update Transverse Leakage Terms                                     *
!         Slava (c) 7.07.1999                                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      real dt_kin

! new variant everything inside
      call TRL_Compute_TRL
 
      call TRL_Compute_TRL0(dt_kin)

      IF(Debug) CALL TRL_DBG_Write_TRL0 

      call TRL_Compute_QLA

      IF(Debug) CALL TRL_DBG_Write_QLA  

      return
      end

      subroutine TRL_Compute_TRL
!=====================================================================*
!           Average over the Node Leakages in X-Y-Z Direction         *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Faces
      implicit none 
      include 'sketch.fh'
! 
! Input: Numb(N_TOT,NDIR), NDD, N_TOT,  Flux(NG,N_TOT),Neib(NDIR,N_TOT),
!        Index_Neib(NDIR,2)
!        Numb(N_TOT, NDIR)
! Module MATrix
!        MAT_FD(NG,NTSJ) 
!        MAT_Tot(NG,NDIR,N_TOT) 
! Output:  TRL_XYZ(NG, N_TOT, NDIR) - Transverse Leakage 
! Local Variables
!      real MAT_Nod 
      real Cur_Left(NG), Cur_Right(NG) 
      integer k, k_left, i_right, nd, n, k_right, kt, nl
      real Flux_Right, Flux_left           

! hexagonal geometry const_leak = 1.0 for XYZ, 2/3 for HEX
      real const_leak, const_trl, const_trl_1 

      IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
         DO nd=1, NDIR-1
            DO k = 1, N_TOT
               DO n = 1, NG
                  trl_1_uvw(n, k, nd, 1) = 0.
                  trl_1_uvw(n, k, nd, 2) = 0.
                END DO
            END DO
         END DO
      END IF

      nl = 0

      do nd = 1, NDD
         i_right = Index_Neib(nd,2)
        IF ( (gmt_crd_type(1:4).EQ."HEXZ") .AND. (nd.ne.NDIR) ) THEN
              const_leak =  0.6666667
        ELSE
              const_leak = 1.0
        END IF
        do kt = 1, N_TOT_DIR(nd)

           k = Numb(kt, nd)
           k_left = Neib(nd, k)
           k_right = Neib(i_right, k)

! Neutron Current on the Left Interface
           if(k_left.eq.I_BOUND_NODE) then
! Boundary Node on the Left, Computing Current
              nl = nl + 1 ! Left Interface on the Boundary
              IF ( i_dr(1, nd) /= 2 ) THEN
              do n = 1,NG 
!                  MAT_Nod = MAT_FD(n, nl) - MAT_Tot(n,nd,k)
                  Cur_Left(n)= - (MAT_Nod(n,nl) + MAT_FD(n,nl))&
                                              * Flux(n,k) 
              end do
              ELSE ! the case of the repeated structure
              DO n = 1, NG
                 IF(k_right.eq.I_BOUND_NODE) then 
                   Flux_Right = 0.
                 ELSE
                   Flux_Right = Flux(n,k_right) 
                 END IF
                 Cur_Left(n) = MAT_FD(n, nl+1)*&
                              ( Flux_Right  - Flux(n,k)) +&
                          MAT_Nod(n,nl+1)*( Flux(n,k) + Flux_Right)
              END DO
              END IF ! i_dr(1, nd) /= 2
            else ! normal internal case
              do n = 1, NG
                 Cur_Left(n) = Cur_Right(n)
               end do
            end if
! Neutron Current on the Right Interface (k_right can be boundary node)
            nl = nl + 1 ! Right Interface 

            const_trl = const_leak/(h_xyz(nd,k)*s_xyz(nd,k))
            const_trl_1 = const_trl

            do n = 1, NG

! Normal case and the bounday node (not repeated structure
            if(k_right.eq.I_BOUND_NODE) then 
               Flux_Right = 0.
            else
               Flux_Right = Flux(n,k_right) 
            end if            

!              MAT_Nod =  MAT_Tot(n, i_right, k) - MAT_FD(n, nl) 
              Cur_Right(n) = - MAT_FD(n, nl)*&
                              ( Flux_Right  - Flux(n,k)) - &
                          MAT_Nod(n,nl)*( Flux(n,k) + Flux_Right) 


! the case of the REPEATED STRUCTURE
              IF(k_right.eq.I_BOUND_NODE.and.i_dr(2,nd)==2) then 
              Flux_left = Flux( n,k_left)
              Cur_Right(n) = MAT_FD(n, nl-1)*&
                              ( Flux(n,k)  - Flux_left)  +&
                          MAT_Nod(n,nl-1)*( Flux_Left + Flux(n,k) ) 
                  END IF

              TRL_XYZ(n,k,nd)=(Cur_Right(n) - Cur_Left(n))*const_trl

! HEXAGONAL GEOMETRY
               IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
              IF ( nd.EQ.1 ) THEN
                  trl_1_uvw(n,k,2,1)=trl_1_uvw(n,k,2,1)+&
                           Cur_Right(n)*const_trl_1
                  trl_1_uvw(n,k,2,2)=trl_1_uvw(n,k,2,2)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,1)=trl_1_uvw(n,k,3,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,2)=trl_1_uvw(n,k,3,2)+&
                           Cur_right(n)*const_trl_1
              ELSE IF(nd.eq.2) THEN
                  trl_1_uvw(n,k,1,1)=trl_1_uvw(n,k,1,1)+&
                           Cur_Right(n)*const_trl_1
                  trl_1_uvw(n,k,1,2)=trl_1_uvw(n,k,1,2)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,1)=trl_1_uvw(n,k,3,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,2)=trl_1_uvw(n,k,3,2)+&
                           Cur_right(n)*const_trl_1
              ELSE IF(nd.eq.3) THEN
                  trl_1_uvw(n,k,1,1)=trl_1_uvw(n,k,1,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,1,2)=trl_1_uvw(n,k,1,2)+&
                           Cur_right(n)*const_trl_1
                  trl_1_uvw(n,k,2,1)=trl_1_uvw(n,k,2,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,2,2)=trl_1_uvw(n,k,2,2)+&
                           Cur_right(n)*const_trl_1
              END IF ! nd
              END IF ! HEX
             end do ! NG
        end do ! N_TOT
      end do ! nd


      return
      end


      subroutine TRL_Compute_TRL0(dt_kin)
!=====================================================================*
!   Average over the Node Transverse Leakages for X-Y-Z Direction     *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
         implicit none 
         include 'sketch.fh'
! 
! Input:  NDD, N_TOT, NG, n_op_xyz, 
!         TRL_XYZ(NG, N_TOT, NDIR),
! for Neutron Kinetics Calculations
!  XS_SF(NG, N_TOT), Flux(NG, N_TOT), al(NG), xp(NG), xpn(NG), 
!               MAT_RHS_K(NG, N_TOT), volume(N_TOT)
      real dt_kin

! Output: TRL_0_XYZ(NG, N_TOT, N_DIR) - Average over the node TRL

! Local Variables:
      real Kin_Trl(NG, N_TOT) ! - additional Terms for Neutron Kinetics
      integer N_DIM
      parameter (N_DIM = NG*N_TOT)
      data Kin_Trl /N_DIM*0./
      integer k, n, nd, ntr
      real a_vol

      Logical Adjoint

!      
      if(Problem_Type.EQ."Kinetics") then

       Adjoint = .False.

        do k = 1, N_TOT
          a_vol = 1./volume(k)
          do n = 1, NG
          kin_trl(n, k) = xs_al(n,k)*Flux(n,k)/dt_kin&
                + ((xp(n)-xpn(n))*Source(k) - MAT_RHS_K(n, k))*a_vol
          end do
        end do
      end if

      do nd = 1, NDD
            do k = 1, N_TOT
               do n = 1, NG
                  TRL_0_xyz(n,k,nd) = kin_trl(n, k)
               end do ! NG
            end do ! N_TOT
      end do ! NDD

!  XYZ Geometry   
      IF ( GMT_CRD_TYPE(1:3).EQ."XYZ") THEN

      do nd = 1, NDD
         do ntr = 1, NDD 
! ntr - opposite direction
            IF(ntr.NE.nd) THEN
            do k = 1, N_TOT
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
             end do ! N_TOT
            END IF
          end do ! ntr
      end do ! nd
! HEXAGONAL GEOMETRY 

      ELSE IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
! uvw directions
         DO nd = 1, NDIR-1
! axial leakage
            ntr = NDIR
            do k = 1, N_TOT
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
            end do ! N_TOT 
       END DO ! uvw
! Axial direction
       nd = NDIR
       do ntr = 1, NDIR-1
! ntr - opposite direction
            do k = 1, N_TOT
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
             end do ! N_TOT
       end do ! ntr
      END IF ! geometry
            
!      if(Debug) call TRL_DBG_Write_TRL0

      return
      end

      subroutine TRL_Compute_QLA
!=====================================================================*
!   Quadratic Leakage Approximation for  X-Y-Z Direction              *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
         implicit none 
         include 'sketch.fh'
! 
         INTEGER nd
      
         DO nd = 1, NDD
          IF( (GMT_CRD_TYPE(1:4) .EQ. "HEXZ") .AND. (nd.NE.NDIR) ) THEN
! Hexagonal geometry 
             CALL TRL_Compute_QLA_HEX(nd)
         ELSE !   
! Cartesian Geometry or axial direction in Hexagonal Geometry
            IF( TRL_Approx(1:3).eq."QLA") THEN
              CALL TRL_Compute_QLA_XYZ(nd)
            END IF
         END IF
         END DO              

! ONLY Hexagonal geometry (radial faces)

!         DO nd = 1, NDD
!           IF( (GMT_CRD_TYPE(1:4) .EQ. "HEXZ") .AND. (nd.NE.NDIR) ) THEN
!              CALL TRL_Compute_QLA_HEX(nd)
!          END IF
!         END DO              
                             
      RETURN
      END 

      subroutine TRL_Compute_QLA_HEX(nd)
!=====================================================================*
!   Quadratic Leakage Approximation for  Hexagonal Geometry           *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none 
      include 'sketch.fh'

      INTEGER i_right, nd, n, k, k_left, k_right, kt
      REAL trl_0_uvw, hk, a0(NG), c_trl(NG)

      real u0_trl_hex, const_trl_vw_11, const_trl_vw_12
      parameter (u0_trl_hex = 0.5)
      parameter (const_trl_vw_11=0.346862)
      parameter (const_trl_vw_12 =0.701622)

!      real const_leak1_z, const_leak2_z 
      REAL  const_leak2_vw, const_leak1_uvw
!      data const_leak1_z /0.0/ !/1.0/
!      data const_leak2_z /0.0/ !/1.15905/
      parameter (const_leak2_vw  =-0.275028)
      parameter (const_leak1_uvw  =1.557158)

      real const_hk ! a/H (x-side of rectangle / pitch of hex)
      parameter (const_hk = 1.0525557)

      real adf_k(NG)
      REAL  cg_nod(2, NDIR), fg_nod(NG, 2, NDIR)  ! 

      real  const_ju 
      parameter (const_ju = 1.42510)

!      REAL a01(NG), transq(2*NG) ! CHECKING

        i_right = Index_Neib(nd,2)

        cg_nod(1, nd) = cg(1, nd)*const_ju
        cg_nod(2, nd) = cg(2, nd)*const_ju


         do kt = 1, N_TOT_DIR(nd)

           k = Numb(kt, nd)
           k_left = Neib_REP_STR(nd, k)
           k_right = Neib_REP_STR(i_right, k)
           hk  = h_xyz(nd, k)
           hk = hk*const_hk


! Step function evaluation
           DO n = 1, NG

             trl_0_uvw = TRL_1_uvw(n, k, nd, 1)+&
                    TRL_1_uvw(n, k, nd, 2)

             TRL_0_xyz(n, k, nd)=( TRL_0_xyz(n, k, nd)+TRL_0_uvw )

             trl_qla(1, n, k, nd)= &
! const_leak1_z=0
!             const_leak1_z*trl_qla(1, n, k, nd) +
                   const_leak1_uvw*&
              ( trl_1_uvw(n, k, nd, 2)-trl_1_uvw(n, k, nd, 1) )

             trl_qla(2, n, k, nd) =  const_leak2_vw*trl_0_uvw
! const_leak2_z=0
!    &              +const_leak2_z*trl_qla(2, n, k, nd) )
                                   
           END DO
! Linear Terms
          IF( TRL_Approx(1:3).EQ."QLA" ) THEN

! left half of the hexagonal node
          IF ( k_left.eq.I_BOUND_NODE ) THEN 
!    Adding left moments into transverse leakage
           DO n = 1, NG
              a0(n) = 2.*XS_D(n,k)/hk
                adf_k(n) = xs_adf(n,1,nd,k)
              fg_nod(n,1,nd) = fg(n,1,nd)*adf_k(n)
              c_trl(n) =  fg_nod(n,1,nd)*TRL_1_uvw(n, k, nd, 1)/&
          ( fg_nod(n,1,nd)*(1. - u0_trl_hex)+cg_nod(1,nd)*a0(n)*&
            sqrt(3.)*0.5 ) 
          END DO
          DO n = 1, NG
             trl_qla(1, n, k, nd) = trl_qla(1, n, k, nd) + &
                  c_trl(n)*const_trl_vw_11
             trl_qla(2, n, k, nd) = trl_qla(2, n, k, nd) - c_trl(n)*&
                  const_trl_vw_12
           END DO
!           IF(k.eq.1.and.nd.eq.1) THEN
!           DO n =1, NG
!             a01(n) = (hk*hk)/(4.*XS_D(n,k))
!             transq(n) = trl_qla(1, n, k, nd)*a01(n)
!             transq(n+NG) = trl_qla(2, n, k, nd)*a01(n)
!           END DO
!            write(*,*) 'c_trl(n)=', c_trl(:)
!            write(*,*) 'trans(1) NEW=', transq(1:NG)
!            write(*,*) 'trans(2) NEW=', transq(NG+1:2*NG)
!              PAUSE
!           END IF 
         END IF ! IF ( k_left.eq.0 )
!  right half of the hexagonal node           
          IF ( k_right.eq.I_BOUND_NODE ) THEN 
           DO  n = 1, NG
              a0(n) = 2.*XS_D(n,k)/hk
              adf_k(n) = xs_adf(n,2,nd,k) 
              fg_nod(n,2,nd) = adf_k(n)*fg(n,2,nd)
              c_trl(n) = - fg_nod(n,2,nd)*TRL_1_uvw(n, k, nd, 2)/&
               (fg_nod(n,2,nd)*(1. - u0_trl_hex)+cg_nod(2, nd)*a0(n)*&
                  sqrt(3.)*0.5 ) 
            END DO
          DO n = 1, NG
             trl_qla(1, n, k, nd) = trl_qla(1, n, k, nd) + &
                  c_trl(n)*const_trl_vw_11
             trl_qla(2, n, k, nd) = trl_qla(2, n, k, nd) + c_trl(n)*&
                  const_trl_vw_12
           END DO
          END IF ! ( k_right.eq.0 ) 
! Adding moments due to the linear term
          END IF ! IF( TRL_Approx(1:3).EQ."QLA") 
        END DO ! kt
      RETURN
      END           

      subroutine TRL_Compute_QLA_XYZ(nd)
!=====================================================================*
!   Quadratic Leakage Approximation for  X-Y-Z Direction              *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
         implicit none 
         include 'sketch.fh'

! 
! Input: Numb(N_TOT,NDIR), NDD, N_TOT, ,Neib_REP_STR(NDIR,N_TOT),
!        Index_Neib(NDIR,2), h_xyz(N_TOT, NDIR) 
!        Numb(N_TOT, NDIR), TRL_0_XYZ(NG, N_TOT, NDIR) 
      INTEGER nd 
! Output: TRL_QLA(2, NG, N_TOT, N_DIR) - expansion coefficients
!         of the transverse leakage
! Local Variables
      integer k, k_left, k_right, i_right, n, k_neib, i_neib, &
             kt
      real     &
       a_left, ah_r, ah_l, g1, a2_r, a1_r, a2_l, a1_l, del_trl_r,&
        del_trl_l

      REAL fg_nod(NG), cg_nod, ah, s0, d1, d2(NG), hk

!      REAL res(NG), c1, c2

! Internal 

         i_right = Index_Neib(nd,2)
         do kt = 1, N_TOT_DIR(nd)

           k = Numb(kt, nd)
           k_left = Neib_REP_STR(nd, k)
           k_right = Neib_REP_STR(i_right, k)

           IF(k_left.ne.I_BOUND_NODE.AND.k_right.ne.I_BOUND_NODE) then
! Usual Case of the QLA (Two Neighbors on the Left and on the Right)
              ah_r = h_xyz(nd,k_right)/h_xyz(nd,k)
              ah_l = h_xyz(nd,k_left)/h_xyz(nd,k)
              g1 = (1.+ah_l)*(1.+ah_r)*(1.+ah_l+ah_r)
              g1 = 0.5/g1
              a2_r = (1. + ah_r)*g1
              a1_r = (1. + 2.*ah_r)*a2_r
              a2_l = (1. + ah_l)*g1
              a1_l = (1. + 2.*ah_l)*a2_l
              do n = 1, NG
                Del_TRL_L = TRL_0_xyz(n,k_left,nd) - TRL_0_xyz(n,k,nd)
                Del_TRL_R = TRL_0_xyz(n,k_right,nd) - TRL_0_xyz(n,k,nd)
                trl_qla(1,n,k,nd) = Del_TRL_R*a1_l - Del_TRL_L*a1_r
                trl_qla(2,n,k,nd) = Del_TRL_R*a2_l + Del_TRL_L*a2_r
              end do ! NG
          ELSE IF(k_left.eq.I_BOUND_NODE.AND.k_right.eq.I_BOUND_NODE) &
                                                                 THEN
! only for 1D cases
           do n = 1, NG
              trl_qla(2,n,k,nd) = 0.
              trl_qla(1,n,k,nd) = 0.
           end do
          ELSE ! k_left.eq.0 or k_right.eq.0
! left boundary node
            IF(k_left.eq.I_BOUND_NODE) THEN
             a_left = -1.
             i_neib = Index_Neib(nd,2)
             k_neib = Neib(i_neib, k)
             cg_nod = cg(1, nd)
             DO n = 1, NG
                fg_nod(n) = fg(n, 1, nd)*xs_adf(n,1,nd,k)
             END DO                  
           ELSE ! k_right.eq.I_BOUND_NODE
             a_left = 1.
             k_neib = Neib(nd, k)
             cg_nod = cg(2, nd)
             DO n = 1, NG
                fg_nod(n) = fg(n, 2, nd)*xs_adf(n,2,nd,k)
             END DO                  
           END IF  

            ah = h_xyz(nd,k_neib)/h_xyz(nd,k)
            hk = h_xyz(nd,k)

            DO n = 1, NG
             d2(n) = 2.*XS_D(n,k)/hk
             g1 = 2.*(ah+1.)*( (ah+2.)*cg_nod*d2(n)+fg_nod(n)*(ah+1.))
             g1 = 1./g1
             d1 = TRL_0_xyz(n, k, nd) - TRL_0_xyz(n, k_neib, nd) 
             s0 = TRL_0_xyz(n, k, nd)
             trl_qla(1,n,k,nd) = g1*a_left*&
                 (d1*(3.*cg_nod*d2(n)+fg_nod(n))-&
                       (ah+1.)*(2.*ah+1.)*fg_nod(n)*s0) 
             trl_qla(2,n,k,nd) = -g1*&
                 (d1*(cg_nod*d2(n)+fg_nod(n))+&
                       (ah+1.)*fg_nod(n)*s0) 


            END DO

          end if ! k_left & k_right /= I_BOUND_NODE
         end do ! N_TOT

!      if(Debug) call TRL_DBG_Write_QLA
!     pause


      return
      end


      subroutine TRL_DBG_Write_TRL0
!=====================================================================*
!   Debug Output of the Neutron Flux and Transverse Leakages  into    *
!  the files 'Output_Debug/Flux.dat' and 'Output_Debug/TRL_XYZ.dat'   *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
      include 'sketch.fh'
      integer nd, n, k

      n = 2
!      open(io_unit,file='Output_Debug/Flux.dat', status ='unknown')
!           write(io_unit,1) (Flux(n, k),k=1,N_TOT)
!      close(io_unit) 

      open(io_unit,file='Output_Debug/TRL_XYZ.dat',status='unknown')
        do nd = 1, NDD
           write(io_unit,*) 'ND = ', nd
           write(io_unit,*) 'nd,n,n1 = ',nd,n,1
           write(io_unit,1) &
            (TRL_0_XYZ(n, Numb(k,nd), nd),k=1,N_TOT_DIR(nd))
        end do ! NDD
      close(io_unit) 

    1 Format(4(1x,9E12.5/),2(1x,8E12.5/),1x,7E12.5/1x,6E12.5/&
                1x,4E12.5)


      return
      end

      subroutine TRL_DBG_Write_QLA
!=====================================================================*
! Debug output 1st and 2nd order transverse leakages coefficients     *
!  into the file 'Output_Debug/TRL_QLA.dat' (when Debug =.True.)         
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      include 'sketch.fh'
      integer nd, n, k

      n = 2
      open(io_unit, file = 'Output_Debug/TRL_QLA.dat', &
                    status ='unknown')
        do nd = 1, 1
           write(io_unit,*) 'ND = ', nd
           write(io_unit,*) 'nd,n, n1 = ',nd,n, 1
           write(io_unit,*) '1st Expansion Coefficient'
           write(io_unit,1) (TRL_QLA(1, n, Numb(k,nd), nd),k=1,NH)
           write(io_unit,*) '2nd Expansion Coefficient'
           write(io_unit,1) (TRL_QLA(2, n, Numb(k,nd), nd),k=1,NH)
       end do ! NDD
      close(io_unit) 

    1 Format(4(1x,9E12.5/),2(1x,8E12.5/),1x,7E12.5/1x,6E12.5/&
                1x,4E12.5)


      return
      end

      subroutine TRL_Compute_QLA_GL
!=====================================================================*
!   Quadratic Leakage Approximation for for X-Y-Z Direction           *
!         German Approach, Linearized Equations for the TRL           *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit NONE 
      include 'sketch.fh'

! Input: Numb(N_TOT,NDIR), NDD, N_TOT, Neib(NDIR,N_TOT),
!        Index_Neib(NDIR,2), h_xyz(N_TOT, NDIR) 
!         TRL_0_XYZ(NG, N_TOT, NDIR), XS_D(NG,N_TOT)
! Output: TRL_QLA(2, NG, N_TOT, N_DIR) - expansion coefficients
!         of the transverse leakage
! Local Variables
      integer k, k_left, k_right, i_right, nd, n,  &
            kt, n_interface , k1, kt1, i
      real hk, d2, d2_p, trl2, trl2_p
      real d2_dimens, trl_left, trl_right
      real TRL_S(NXYZ_MAX+1)

      do nd = 1, NDD
         i_right = Index_Neib(nd,2)

          do n = 1, NG

! End of the Boundary Conditions

           do kt = 1, N_TOT_DIR(nd)

              k = Numb(kt, nd)
              k_left = Neib(nd, k)
              k_right = Neib(i_right, k)
              hk = h_xyz(nd,k)
              d2_dimens = XS_D(n,k)
              d2 = 2.*d2_dimens/hk
              trl2 = trl_0_xyz(n,k,nd)/d2_dimens

!              if(nd.eq.3.and.kt.le.NH) then        
!               write(*,*) 'k =', k, 'kt =', kt
!               write(*,*) 'k_left =', k_left, 'k_right=', k_right
!               write(*,*) 'd2 =', d2
!               pause
!              end if  

              if(k_left.eq.I_BOUND_NODE) then

                 n_interface = 1
                 trl_s(n_interface) = cg(1,nd)*d2*trl2/&
                           (cg(1,nd)*d2 + fg(n,1,nd)) 
                 d2_p = d2
                 trl2_p = trl2
              
              else 
                   
                n_interface = n_interface + 1
                trl_s(n_interface) = (d2*trl2 +d2_p*trl2_p)/&
                                                 (d2+d2_p)
                d2_p = d2
                trl2_p = trl2
       
             end if

             if(k_right.eq.I_BOUND_NODE) then

                 n_interface = n_interface + 1
                 trl_s(n_interface) = cg(2,nd)*d2*trl2/&
                           (cg(2,nd)*d2 + fg(n,2,nd)) 

!                 write(*,*) 'trl_s(i) =', (trl_s(i),i=1, n_interface)
!           pause

! Computing the Expansion Coefficients Using the Surface Values
                do i = 1, n_interface - 1 
                   kt1 = kt - (n_interface - 1 - i)

                   k1 = Numb(kt1,nd)
                   trl_left = d2_dimens*trl_s(i)
                   trl_right = d2_dimens*trl_s(i+1)
!                   trl_qla(1,n,k1,nd) = 0.5*(trl_s(i+1) - trl_s(i))
!                   trl_qla(2,n,k1,nd) = 0.5*(trl_s(i+1) + trl_s(i)) -&
!                               trl0_xyz(n,k1,nd)
                   trl_qla(1,n,k1,nd) = 0.5*(trl_right - trl_left)
                   trl_qla(2,n,k1,nd) = 0.5*(trl_right + trl_left) -&
                               trl_0_xyz(n,k1,nd)
!                 if(nd.eq.3.and.kt1.le.NH) then        
!                   write(*,*) 'i =', i
!                   write(*,*) 'kt1=', kt1
!                   write(*,*) 'nd =', nd, 'k1 =', k1, 'n_interface=', &
!                       n_interface
!                   write(*,*) 'trl_s(i), trl_s(i+1) =', trl_s(i), &
!                                      trl_s(i+1)
!                   write(*,*) 'trl0 =', trl_0_xyz(n,k1,nd)
!                   trl_tl = trl_0_xyz(n,k1,nd) - trl_qla(1,n,k1,nd) +&
!                    trl_qla(2,n,k1,nd)
!                   trl_tr = trl_0_xyz(n,k1,nd) + trl_qla(1,n,k1,nd) +&
!                    trl_qla(2,n,k1,nd)
!                   write(*,*) 'trl test: left, right: =',trl_tl,trl_tr
!                   pause
!                end if

                end do

            end if                          
 
         end do ! kt
        end do ! NG
      end do ! NDD

!      if(Debug) call QLA_Write
!     pause

      return 
      end

      subroutine TRL_Compute_QLA_G
!=====================================================================*
!   Quadratic Leakage Approximation for for X-Y-Z Direction           *
!         German Approach, We are solving Continuity Equations        *
!          for the transverse leakage at the interface                *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit NONE 
      include 'sketch.fh'

! Input: Numb(N_TOT,NDIR), NDD, N_TOT, Neib(NDIR,N_TOT),
!        Index_Neib(NDIR,2), h_xyz(N_TOT, NDIR) 
!         TRL0_XYZ(NG, N_TOT, NDIR), XS_D(NG,N_TOT)
! Output: TRL_QLA(2, NG, N_TOT, N_DIR) - expansion coefficients
!         of the transverse leakage
! Local Variables
      integer k, k_left, k_right, i_right, nd, n,  &
            kt, n_interface , k1, kt1, i
      real hk, d2, d2_p, trl2, trl2_p
      real d2_dimens, trl_left, trl_right
      real AA(3,NXYZ_MAX+1), bb(NXYZ_MAX+1), TRL_S(NXYZ_MAX+1)
!        AAA(3,NXYZ_MAX+1), res(NXYZ_MAX+1)
! Used Subroutines: 
!                  progonka

      do nd = 1, NDD
         i_right = Index_Neib(nd,2)

         do n = 1, NG

          do kt = 1, N_TOT_DIR(nd)

              k = Numb(kt, nd)
              k_left = Neib(nd, k)
              k_right = Neib(i_right, k)
              hk = h_xyz(nd,k)
              d2_dimens = XS_D(n,k)
              d2 = 2.*d2_dimens/hk
              trl2 = trl_0_xyz(n,k,nd)/d2_dimens

              if(k_left.eq.I_BOUND_NODE) then

                 n_interface = 1
                 AA(2,n_interface) = fg(n,1,nd) + 2.*cg(1,nd)*d2
                 AA(3,n_interface) = cg(1,nd)*d2
                 bb(n_interface) = 3.*trl2*cg(1,nd)*d2
                 d2_p = d2
                 trl2_p = trl2
               else 
                 n_interface = n_interface + 1
                 AA(1,n_interface) = d2_p
                 AA(2,n_interface) = 2.*(d2 + d2_p)
                 AA(3,n_interface) = d2
                 bb(n_interface) = 3.*(d2*trl2 + d2_p*trl2_p)
                 d2_p = d2
                trl2_p = trl2
              end if

              if(k_right.eq.I_BOUND_NODE) then

                 n_interface = n_interface + 1
                 AA(1,n_interface) = cg(2,nd)*d2
                 AA(2,n_interface) = fg(n,2,nd) + 2.*cg(2,nd)*d2
                 bb(n_interface) = 3.*trl2*cg(2,nd)*d2
! Solving the Equations for the TRL at the interfaces
!               do i = 1, n_interface
!                   AAA(1,i) = AA(1,i)
!                   AAA(2,i) = AA(2,i)
!                   AAA(3,i) = AA(3,i)
!               end do

               call MSC_progonka(n_interface, AA, BB, TRL_S)
! Computing the Expansion Coefficients Using the Surface Values
                 
!                 write(*,*) 'trl_s(i) =', (trl_s(i),i=1, n_interface)
!           pause

                do i = 1, n_interface - 1 
                   kt1 = kt - (n_interface - 1 - i)
                   k1 = Numb(kt1,nd)
                   trl_left = d2_dimens*trl_s(i)
                   trl_right = d2_dimens*trl_s(i+1)
                   trl_qla(1,n,k1,nd) = 0.5*(trl_right - trl_left)
                   trl_qla(2,n,k1,nd) = 0.5*(trl_right + trl_left) -&
                               trl_0_xyz(n,k1,nd)
!                   trl_qla(1,n,k1,nd) = 0.5*(trl_s(i+1) - trl_s(i))
!                   trl_qla(2,n,k1,nd) = 0.5*(trl_s(i+1) + trl_s(i)) -&
!                               trl0_xyz(n,k1,nd)
!                   write(*,*) 'nd =', nd, 'k1 =', k1
!                   write(*,*) 'trl_s(i), trl_s(i+1) =', trl_s(i), &
!                                      trl_s(i+1)
!                   write(*,*) 'trl0 =', trl0_xyz(n,k1,nd)
!                   pause
!                 if(nd.eq.3.and.kt1.le.NH) then        
!                   write(*,*) 'i =', i
!                   write(*,*) 'kt1=', kt1
!                   write(*,*) 'nd =', nd, 'k1 =', k1, 'n_interface=', &
!                       n_interface
!                   write(*,*) 'trl_s(i), trl_s(i+1) =', trl_s(i), &
!                                      trl_s(i+1)
!                   write(*,*) 'trl0 =', trl0_xyz(n,k1,nd)
!                   trl_tl = trl0_xyz(n,k1,nd) - trl_qla(1,n,k1,nd) +&
!                    trl_qla(2,n,k1,nd)
!                   trl_tr = trl0_xyz(n,k1,nd) + trl_qla(1,n,k1,nd) +&
!                    trl_qla(2,n,k1,nd)
!                   write(*,*) 'trl test: left, right: =',trl_tl,trl_tr
!                   pause
!                 end if

                end do

            end if                          
 
         end do ! kt
       end do ! NG
      end do ! NDD

!      if(Debug) call QLA_Write
!     pause

      return 
      end


      subroutine TRL_Compute_TRL_OLD
!=====================================================================*
!           Average over the Node Leakages in X-Y-Z Direction         *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      USE GeoMeTry_Faces
      implicit none 
      include 'sketch.fh'
! 
! Input: Numb(N_TOT,NDIR), NDD, N_TOT,  Flux(NG,N_TOT),Neib(NDIR,N_TOT),
!        Index_Neib(NDIR,2)
!        Numb(N_TOT, NDIR)
! Module MATrix
!        MAT_FD(NG,NTSJ) 
!        MAT_Tot(NG,NDIR,N_TOT) 
! Output:  TRL_XYZ(NG, N_TOT, NDIR) - Transverse Leakage 
! Local Variables
!      real MAT_Nod 
      real Cur_Left(NG), Cur_Right(NG) 
      integer k, k_left, i_right, nd, n, k_right, kt, nl
      real Flux_Right           

! hexagonal geometry const_leak = 1.0 for XYZ, 2/3 for HEX
      real const_leak, const_trl, const_trl_1 

      IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
         DO nd=1, NDIR-1
            DO k = 1, N_TOT
               DO n = 1, NG
                  trl_1_uvw(n, k, nd, 1) = 0.
                  trl_1_uvw(n, k, nd, 2) = 0.
                END DO
            END DO
         END DO
      END IF

      nl = 0

      do nd = 1, NDD
         i_right = Index_Neib(nd,2)
        IF ( (gmt_crd_type(1:4).EQ."HEXZ") .AND. (nd.ne.NDIR) ) THEN
              const_leak =  0.6666667
        ELSE
              const_leak = 1.0
        END IF
        do kt = 1, N_TOT_DIR(nd)

           k = Numb(kt, nd)
           k_left = Neib(nd, k)

! Neutron Current on the Left Interface
           if(k_left.eq.I_BOUND_NODE) then
! Boundary Node on the Left, Computing Current
              nl = nl + 1 ! Left Interface on the Boundary
              do n = 1,NG 
!                  MAT_Nod = MAT_FD(n, nl) - MAT_Tot(n,nd,k)
                  Cur_Left(n)= - (MAT_Nod(n,nl) + MAT_FD(n,nl))&
                                              * Flux(n,k) 
              end do
            else
              do n = 1, NG
                 Cur_Left(n) = Cur_Right(n)
               end do
            end if
! Neutron Current on the Right Interface (k_right can be boundary node)
            nl = nl + 1 ! Right Interface 
            k_right = Neib(i_right, k)

            const_trl = const_leak/(h_xyz(nd,k)*s_xyz(nd,k))
            const_trl_1 = const_trl

            do n = 1, NG

            if(k_right.eq.I_BOUND_NODE) then 
               Flux_Right = 0.
            else
               Flux_Right = Flux(n,k_right) 
            end if            

!              MAT_Nod =  MAT_Tot(n, i_right, k) - MAT_FD(n, nl) 
              Cur_Right(n) = - MAT_FD(n, nl)*&
                              ( Flux_Right  - Flux(n,k)) - &
                          MAT_Nod(n,nl)*( Flux(n,k) + Flux_Right) 

              TRL_XYZ(n,k,nd)=(Cur_Right(n) - Cur_Left(n))*const_trl
! HEXAGONAL GEOMETRY
               IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
              IF ( nd.EQ.1 ) THEN
                  trl_1_uvw(n,k,2,1)=trl_1_uvw(n,k,2,1)+&
                           Cur_Right(n)*const_trl_1
                  trl_1_uvw(n,k,2,2)=trl_1_uvw(n,k,2,2)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,1)=trl_1_uvw(n,k,3,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,2)=trl_1_uvw(n,k,3,2)+&
                           Cur_right(n)*const_trl_1
              ELSE IF(nd.eq.2) THEN
                  trl_1_uvw(n,k,1,1)=trl_1_uvw(n,k,1,1)+&
                           Cur_Right(n)*const_trl_1
                  trl_1_uvw(n,k,1,2)=trl_1_uvw(n,k,1,2)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,1)=trl_1_uvw(n,k,3,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,3,2)=trl_1_uvw(n,k,3,2)+&
                           Cur_right(n)*const_trl_1
              ELSE IF(nd.eq.3) THEN
                  trl_1_uvw(n,k,1,1)=trl_1_uvw(n,k,1,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,1,2)=trl_1_uvw(n,k,1,2)+&
                           Cur_right(n)*const_trl_1
                  trl_1_uvw(n,k,2,1)=trl_1_uvw(n,k,2,1)-&
                           Cur_left(n)*const_trl_1
                  trl_1_uvw(n,k,2,2)=trl_1_uvw(n,k,2,2)+&
                           Cur_right(n)*const_trl_1
              END IF ! nd
              END IF ! HEX
             end do ! NG
        end do ! N_TOT
      end do ! nd


      return
      end


      subroutine TRL_Compute_TRL0_OLD(dt_kin)
!=====================================================================*
!   Average over the Node Transverse Leakages for X-Y-Z Direction     *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
         implicit none 
         include 'sketch.fh'
! 
! Input:  NDD, N_TOT, NG, n_op_xyz, 
!         TRL_XYZ(NG, N_TOT, NDIR),
! for Neutron Kinetics Calculations
!  XS_SF(NG, N_TOT), Flux(NG, N_TOT), al(NG), xp(NG), xpn(NG), 
!               MAT_RHS_K(NG, N_TOT), volume(N_TOT)
      real dt_kin

! Output: TRL_0_XYZ(NG, N_TOT, N_DIR) - Average over the node TRL

! Local Variables:
      real Kin_Trl(NG, N_TOT) ! - additional Terms for Neutron Kinetics
      integer N_DIM
      parameter (N_DIM = NG*N_TOT)
      data Kin_Trl /N_DIM*0./
      integer k, n, nd, ntr
      real a_vol

      Logical Adjoint

!      
      if(Problem_Type.EQ."Kinetics") then

       Adjoint = .False.

        do k = 1, N_TOT
          a_vol = 1./volume(k)
          do n = 1, NG
          kin_trl(n, k) = xs_al(n,k)*Flux(n,k)/dt_kin&
                + ((xp(n)-xpn(n))*Source(k) - MAT_RHS_K(n, k))*a_vol
          end do
        end do
      end if

      do nd = 1, NDD
            do k = 1, N_TOT
               do n = 1, NG
                  TRL_0_xyz(n,k,nd) = kin_trl(n, k)
               end do ! NG
            end do ! N_TOT
      end do ! NDD

!  XYZ Geometry   
      IF ( GMT_CRD_TYPE(1:3).EQ."XYZ") THEN

      do nd = 1, NDD
         do ntr = 1, NDD 
! ntr - opposite direction
            IF(ntr.NE.nd) THEN
            do k = 1, N_TOT_DIR(nd)
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
             end do ! N_TOT
            END IF
          end do ! ntr
      end do ! nd
! HEXAGONAL GEOMETRY 

      ELSE IF ( GMT_CRD_TYPE(1:4).EQ."HEXZ") THEN
! uvw directions
         DO nd = 1, NDIR-1
! axial leakage
            ntr = NDIR
            do k = 1, N_TOT_DIR(nd)
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
            end do ! N_TOT 
       END DO ! uvw
! Axial direction
       nd = NDIR
       do ntr = 1, NDIR-1
! ntr - opposite direction
            do k = 1, N_TOT_DIR(nd)
                  do n = 1, NG
                     TRL_0_xyz(n,k,nd) = TRL_0_xyz(n,k,nd) + &
                                   TRL_xyz(n,k,ntr)
                  end do ! NG
             end do ! N_TOT
       end do ! ntr
      END IF ! geometry
            
!      if(Debug) call TRL_DBG_Write_TRL0

      return
      end


      subroutine TRL_Compute_QLA_OLD
!=====================================================================*
!   Quadratic Leakage Approximation for for X-Y-Z Direction           *
!                   Vyachreslav Zimin (c) 1988-1998                   *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
           USE GeoMeTry_Boundary
         implicit none 
         include 'sketch.fh'

! 
! Input: Numb(N_TOT,NDIR), NDD, N_TOT, ,Neib(NDIR,N_TOT),
!        Index_Neib(NDIR,2), h_xyz(N_TOT, NDIR) 
!        Numb(N_TOT, NDIR), TRL_0_XYZ(NG, N_TOT, NDIR) 
! Output: TRL_QLA(2, NG, N_TOT, N_DIR) - expansion coefficients
!         of the transverse leakage
! Local Variables
      integer k, k_left, k_right, i_right, nd, n, kb, k_neib, i_neib, &
            nl, kt, nd_right
      real alfa_1(NG), alfa_2(NG), t_extr, pol_1, pol_2, alfa_0, &
       a_left, ah_r, ah_l, g1, a2_r, a1_r, a2_l, a1_l, del_trl_r,&
        del_trl_l, a22, a11, ah_hk, hk

!      REAL res(NG), c1, c2, cg_nod, fg_nod(NG), d2(NG), s0  

! Internal 
      do nd = 1, NDD
         i_right = Index_Neib(nd,2)
         do kt = 1, N_TOT_DIR(nd)

           k = Numb(kt, nd)
           k_left = Neib(nd, k)
           k_right = Neib(i_right, k)

           if(k_left.ne.I_BOUND_NODE.AND.k_right.ne.I_BOUND_NODE) then
! Usual Case of the QLA (Two Neighbors on the Left and on the Right)
              ah_r = h_xyz(nd,k_right)/h_xyz(nd,k)
              ah_l = h_xyz(nd,k_left)/h_xyz(nd,k)
              g1 = (1.+ah_l)*(1.+ah_r)*(1.+ah_l+ah_r)
              g1 = 0.5/g1
              a2_r = (1. + ah_r)*g1
              a1_r = (1. + 2.*ah_r)*a2_r
              a2_l = (1. + ah_l)*g1
              a1_l = (1. + 2.*ah_l)*a2_l
              do n = 1, NG
                Del_TRL_L = TRL_0_xyz(n,k_left,nd) - TRL_0_xyz(n,k,nd)
                Del_TRL_R = TRL_0_xyz(n,k_right,nd) - TRL_0_xyz(n,k,nd)
                trl_qla(1,n,k,nd) = Del_TRL_R*a1_l - Del_TRL_L*a1_r
                trl_qla(2,n,k,nd) = Del_TRL_R*a2_l + Del_TRL_L*a2_r
              end do ! NG

          end if ! k_left & k_right /= 0
         end do ! N_TOT
       end do ! NDD

! Boundary nodes

      do kb = 1, N_BOUND


         k = k_bound(kb)
         nd = nd_bound(kb)
         nl = nl_bound(kb)
         hk = h_xyz(nd,k)

! treatment of 1D case
         k_left = Neib(nd,k)
         nd_right = Index_Neib(nd,2)
         k_right = Neib(nd_right, k)

         if(k_left.eq.I_BOUND_NODE.AND.k_right.eq.I_BOUND_NODE) then
           do n = 1, NG
              trl_qla(2,n,k,nd) = 0.
              trl_qla(1,n,k,nd) = 0.
           end do
         else
            if(nl.eq.1) then
! Left Boundary 
             a_left = 1.
             i_neib = Index_Neib(nd,2)
             k_neib = Neib(i_neib, k)
            else
! Right Boundary (nl = 2)
             a_left = -1.
             k_neib = Neib(nd, k)
            end if

            ah_hk = h_xyz(nd,k_neib)/hk
            a11 = 1. + ah_hk
            a22 = 1. + 2.*ah_hk

            if(i_dr(nl,nd).eq.0) then
! Extrapolation Distance
               do n = 1, NG
                  t_extr = 1. + 2.*dr(n,nl,nd)/hk
                  pol_1 = t_extr
                  pol_2 = 0.5*(3.*t_extr**2 - 1.)
                  alfa_1(n) = 1. / pol_1
                  alfa_2(n) = pol_2/pol_1
               end do
            else 
! Logarithmic Derivative
               do n = 1, NG
                  alfa_0 = dr(n,nl,nd)*hk/(2.*XS_D(n,k))
                  alfa_1(n) = alfa_0 / (alfa_0 + 1)
                  alfa_2(n) = (alfa_0 + 3.)/(alfa_0 + 1)
               end do
            end if

            do n = 1, NG
               trl_qla(2,n,k,nd)=(TRL_0_xyz(n, k_neib, nd) - &
                 TRL_0_xyz(n, k, nd)*(1. + alfa_1(n)*a11))/&
                 (a11*(a22 + alfa_2(n)))
               trl_qla(1,n,k,nd)=a_left*(alfa_2(n)*trl_qla(2,n,k,nd)+ &
                      alfa_1(n)*TRL_0_xyz(n,k,nd))
            end do !

         end if ! (k_left == I_BOUND_NODE) .AND. (k_right == I_BOUND_NODE)

      end do ! N_BOUND

!      if(Debug) call TRL_DBG_Write_QLA
!     pause


      return
      end

