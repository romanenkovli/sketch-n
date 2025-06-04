      complex function MAT_F_COMPL_ANM(gamma,eigenv, m)
!=====================================================================*
!     Response function for Multigroup ANM                            *
! 1 Node Average Flux                                                 *
! 2 neutron flux on the right boundary                                *
! 3-5 Transverse Leakage Expansion Coefficients                       *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x, 
!   m - number of the matrix function
      complex eigenv
      complex gamma
      integer m

      SELECT CASE (m)
             CASE (1)
                  MAT_F_COMPL_ANM = -1./gamma**2 + eigenv
             CASE (2)
                  MAT_F_COMPL_ANM = 1./gamma
             CASE (3)
                  MAT_F_COMPL_ANM = 1. - &
                                        (1.-gamma)/(eigenv*gamma**2)
             CASE (4)
                  MAT_F_COMPL_ANM = (1. - gamma)/(gamma*eigenv)
             CASE (5)
                  MAT_F_COMPL_ANM = &
                  -(3. - gamma*(3.+eigenv))/((gamma*eigenv)**2)
      END SELECT

      return
      end


      complex function MAT_F_COMPL_PNM(t, m)
!=====================================================================*
!     Response function for Multigroup PNM                            *
! 1 Node Average Flux                                                 *
! 2 neutron flux on the right boundary                                *
! 3-5 Transverse Leakage Expansion Coefficients                       *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
! Input: gamma=tanh(sqrt(x))/sqrt(x), eigenv = x, 
!   m - number of the matrix function
      complex t
      integer m

      SELECT CASE (m)
             CASE (1)
                MAT_F_COMPL_PNM = &
             (-66. + 2.*t + 2875./(15. + t) - 2744./(21. + 2.*t))/5.
             CASE (2)
                MAT_F_COMPL_PNM = (3.*(5. + 2.*t))/(15. + t)
             CASE (3)
                MAT_F_COMPL_PNM = (2.*(525. + t*(15. + 2.*t)))/&
                                     (5.*(15. + t)*(21. + 2*t)) 
             CASE (4)
                MAT_F_COMPL_PNM = 5./(15. + t)
             CASE (5)
                MAT_F_COMPL_PNM = &
                   (21.*(5. + 2.*t))/(5.*(15. + t)*(21. + 2.*t))
      END SELECT

      return
      end


      subroutine MAT_Set_MatFunc_NG(&
                                         xappa, Matr_Func, NG, ANM)
!=====================================================================*
! Computing Response Matrix Functions for the node                    *
!       ANM & PNM (> 2 neutron energy groups)                         * 
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
!      INCLUDE 'MATHS.FI'
! Input: Eigenvalues of the Matrix matrix B^2
      integer NG ! NG > 2
      real xappa(NG,NG) ! -buckling matrix
      LOGICAL ANM ! True if ANM, else if PNM
! Output:
      integer N_Matr_Func
      parameter (N_Matr_Func = 5)
      real Matr_Func(NG, NG, N_Matr_Func)
! Response Matrix Functions with respect:
! 1) Average Flux 2) Right Boundary Flux 3) Average TRL 4) 1st TRL
! expansion coefficient 5) 2nd TRL expansion coefficient
! Local Variables:
      complex eigenv(NG), Pol_Coeff(NG), Funct_Val(NG)
      complex lambda(NG), gamm(NG)
      real Real_Pol_Coeff(NG)

      external EVLRG, S_POLRG , CTAN ! ISML Library subroutines 
!     EVLRG - compute eigenvalues of the real matrix
!     POLGR  - compute matrix Polynomial Using Horner's Rule
      external  MSC_POLCOE
      external MAT_F_COMPL_ANM, MAT_F_COMPL_PNM
      complex MAT_F_COMPL_ANM, MAT_F_COMPL_PNM, CTAN
!     MAT_F_COMPL_ANM - definition of the 5th matrix functions
!     MSC_POLCOE - find coefficients of the matrix polynomial

      integer m, n

! Computing Eigenvalues of the Buckling Matrix
      call EVLRG(NG, Xappa, NG, Eigenv)

      if(ANM) then
         do n = 1, NG
            lambda(n) = CSQRT( - eigenv(n) )
            gamm(n) = CTAN( lambda(n))/ lambda(n) 
         end do
      end if

        do m = 1, N_Matr_Func

            if(ANM) then
               do n = 1, NG
                  Funct_Val(n)=&
                    MAT_F_COMPL_ANM( gamm(n), eigenv(n), m )
              end do
            else
               do n = 1, NG
                  Funct_Val(n) = &
                    MAT_F_COMPL_PNM( eigenv(n), m )
               end do
            end if

            call MSC_POLCOE(Eigenv, Funct_Val, NG, Pol_Coeff)

            Real_Pol_Coeff(:) = REAL(Pol_Coeff(:))

            call S_POLRG(NG, Xappa, NG, NG, Real_Pol_Coeff, &
                                  Matr_Func(1,1,m), NG)

      end do

      return
      end

