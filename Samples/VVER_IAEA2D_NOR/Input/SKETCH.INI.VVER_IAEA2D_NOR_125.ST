c============================================================================c
c          SKETCH-N version 1.0: Nodal Neutron Diffusion Code for            c
c            Solving Steady-State & Kinetics Problems                        c
c                                                                            c
c             Moscow Engineering Physics Institute                           c
c                Tokyo Institute of Technology                               c
c              Japan Atomic Energy Research Institute                        c 
c                                                                            c
c       Author:  Vyacheslav G. Zimin                                         c
c                                                                            c
c                    (c) 1999 All Rights Reserved                            c
c                                                                            c
c                               NOTICE                                       c
c                                                                            c
c  Permission to use, copy, modify, and distribute this software and         c
c  its documentation for any purpose and without fee is hereby granted       c
c  provided that the above copyright notice appear in all copies and         c
c  that both the copyright notice and this permission notice appear in       c
c  supporting documentation.                                                 c
c                                                                            c
c  Neither the Institutions  nor the Author make any                         c
c  representations about the suitability of this software for any            c
c  purpose.  This software is provided ``as is'' without express or          c
c  implied warranty.                                                         c
c============================================================================c

!=============================================================================!
! INI_PROBLEM  - the basic problem description & solution methods             !
! Problem_Type - type of the computed problem                                 !
!                  "Eigenvalue", BoronSearch", "Kinetics"          /   ""   / !
! Nodal_Method - Choice of the Nodal Method:                                  !
!               "PNM", "SANM", "ANM", "MCFD", "PNM1"               / "SANM" / !
! TH_Model - Thermal-hydraulics Model Used for the Calculations    / 'None' / !
!            (Possible Choices 'Internal', 'External', 'None')                !
! Kinetics_Method - Method of the Solution of Kinetics Equation               !
!      "DRT"-Direct; "IQS"- Improved Quasi-Static; "PNT" - point   / "DRT"  / !  
! Iter_Solver - iterative solvers for kinertics :                             !
!                  CSA, CSI, CG, BCGSTAB, TFQMR, FOM, GMRES        / "CSA" /  !
!=============================================================================!

!=============================================================================!
! INI_FILES - name of the input and restart files                             !
! File_Input- Input File Name of the Geometry Data                /   ""   /  !
! File_Reference - Name of the Reference Solution File            /   ""   /  !
! File_Dist - name of the DIST*.txt file with distributions for               !
!              Ringhals                                           /   ""  /   !
! File_CD - name of the CD*.txt file with XS Tables (Ringhals)    /   ""  /   !
! File_Map   - Input File for the Mapping to the TRAC Geometry                !
!               Only for the Coupling with the TRAC               /   ""   /  !
! File_DMP_In - Input Restart File                                /   ""   /  !
! File_DMP_Out_st - Output Restart File for                                   !
!                                  Steady-State Calculations      /   ""   /  !
! File_DMP_Out_Kin - Output Restart File for                                  !
!                                  Kinetics     Calculations      /   ""   /  !
!=============================================================================!

!=============================================================================!
! INI_CONVERGENCE  - convergence criterion                                    !
! N_Out_Max  - Convergence Criterion for the Nonlinear Iterations             !
!                                                      Steady-State   /1 /    !
!                                                      Kinetics       /1000 / !
! N_Outer - Number of outer per nonlinear (thermalhydraulics)                 !
!                                                      Steady-State   /10  /  !
!                                                      Kinetics       /1000/  !
! E_Boron_Start - Accuracy of the Eigenvalue when critical boron              !
!            search starts                                            /1.E-2/ !
! E_Critical - Convergence Criterion of the Boron Critical Search     /1.E-5/ !
! E_Outer_L  - Converence Criterion for the Eigenvalue (Local)        /1.E-5/ !
! E_Flux_L - Convergence Criterion for the Flux                               !
!                                                      Steady-State   /1.E-5/ !
!                                                      Kinetics       /1.E-4/ !
! N_Inter     - Number of internal iterations per outer               /2    / !
! E_Inter - Converence Criterion for the Innner Iterations            /1.E-8/ !
!=============================================================================!

!=============================================================================!
! INI_CHEBYSHEV - Parameters for Chebyshev acceleration procedure             !
! Xbe - Estimate of Minimum Eigenvalue of the Iteration Matrix        /0.0  / !
! Xme_Ini - Estimate of the Dominance Ratio (Spectral Radius)         /0.8  / !
! F_Cheb - Adaptive Parameter of the Chebyshev Procedure                      !
!           (Non-Adaptive - 0, Adaptive - 0.65-0.85)                          !
!                                                      Steady-State   /0.0  / !
!                                                      Kinetics       /0.8  / !
! NPolins - Number of the Iteration when Chebyshev Starts             / 5   / !
! Delta_Shift - Inverse Value of the Wieland shift                    / 0.  / !
!=============================================================================!

!=============================================================================!
! INI_TIME_STEP - Time Step Size Selection                                    !
! I_Auto -  Flag of the Automatic Time Step Selection             /    0    / !
! NP_View - Number of Time Step Interval                          /    1    / !
! Ttv - Time Moments of the Time Step Intervals                               !
!       Ttv(NP_View) - End of the Transien                        /    1.   / !
! Dt_Input(NP_View) - Time Step Size for the Time Interval                    !
!         [Ttv(i-1) - Ttv(i)]                                     /  0.01   / !
! St_Eps - Accuracy Criterion of the Automatic Time Step Control  /  5.E-03 / !
! Facmax - Maximum Increase of the Time Step Size                 /    2    / !
! Dt_Step_Max - Maximum Time Step Size (s)                        /    1.   / !
! N_Zap -  Output into GRF file at the  N_ZAPth Time Step         /    1    / !
!=============================================================================!

 &INI_PROBLEM
   Problem_Type = "Steady-State"
   Steady_State_Type = "Eigenvalue"
   TH_Model = "None"
   Nodal_Method = "PNM"
 &END 

 &INI_FILES 
  FILE_INPUT      = "Input/VVER_IAEA2D_NOR_125.DAT"
  FILE_REFERENCE  = "Input/VVER_IAEA2D_NOR_125.ref"
  FILE_DMP_OUT_ST = "Restart/VVER_IAEA2D_NOR_125.dat"
 &END

 &INI_CONVERGENCE
 &END

 &INI_CHEBYSHEV
  Delta_Shift = 0.
 &END

 &INI_TIME_STEP
 &END
