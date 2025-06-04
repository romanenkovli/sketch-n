       MODULE XS_POL_FTYPES
         
       USE PRECISION, ONLY : dp
       IMPLICIT NONE
       INTEGER, PARAMETER :: NN_FTYPE = 200
       INTEGER            :: N_FTYPE
       INTEGER            :: index_ftype(NN_FTYPE)
       INTEGER            :: index_internal(NN_FTYPE)
       CHARACTER*300      :: title_ftype(NN_FTYPE)
       REAL(dp)           :: hm_mass_of_fuel(NN_FTYPE) 

       END MODULE XS_POL_FTYPES

       module spline_lib

       USE PRECISION, ONLY : dp
       USE XS_POL_FTYPES

       PRIVATE

       INTEGER, PARAMETER :: NNKNOT=200
       INTEGER, PARAMETER :: N_SPLINE_ORDER=3
       INTEGER, PARAMETER :: max_ny=100

       INTEGER  :: NKNOT(NN_FTYPE), NY_SP(NN_FTYPE)

       REAL(dp) :: X_SP(1:NNKNOT, NN_FTYPE)

       REAL(dp) ::  &
        SP_COEFF(1:max_ny, 1:NNKNOT-1, 0:N_SPLINE_ORDER, NN_FTYPE)

       PUBLIC :: spline_eval, SPLINE_INPUT, spline_eval_der

      CONTAINS

      SUBROUTINE SPLINE_INPUT(id, i_ftype)
      INTEGER, INTENT(IN) :: id, i_ftype
      INTEGER :: i
      CHARACTER           :: end_of_data

         READ(id, *) NKNOT(i_ftype), NY_SP(i_ftype)
         READ(id,*) (X_SP(i, i_ftype), i=1,NKNOT(i_ftype)) 
         READ(id,*) (((SP_COEFF( i, j, k, i_ftype), &
                                      i=1,NY_SP(i_ftype)), &
                                      j=1,NKNOT(i_ftype)-1), &
                                      k=0,N_SPLINE_ORDER)
         READ(id,*) end_of_data 
      RETURN
      END SUBROUTINE SPLINE_INPUT 

      function spline_eval(m, u, i_ftype) result(value)
!     SP_COEFF goes through MODULE
      INTEGER, INTENT(IN) :: m, i_ftype
      INTEGER             :: N
! N   -  number of knots
! N-1 -  number of intervals
! M   -  number of y 
      REAL(dp), INTENT(IN) ::  u
!      REAL(dp), INTENT(IN) ::  y(m,n-1), b(m,n-1), c(m,n-1), d(m,n-1)
      REAL(dp) :: value(M), X(NNKNOT)
      integer i, im
      REAL(dp) ::  dx
!      data i/1/

      i = 1

      IF( M .GT. NY_SP(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Error! M .GT. NY_SP(i_ftype) :',M,NY_SP(i_ftype)
        stop
      ELSE IF ( M .LT. NY_SP(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Warning! M .LT. NY_SP(i_ftype) :',M,NY_SP(i_ftype)
      END IF

      N = NKNOT(i_ftype)

      X(1:N) = X_SP(1:N, i_ftype)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!

! Binary search should be changed in the fu 
!      SAVE :: i
!      if ( i .ge. n ) i = 1
!      if ( u .lt. x(i) ) go to 10
!      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
!   10 i = 1
!      j = n+1
!   20 k = (i+j)/2
!      if ( u .lt. x(k) ) j = k
!      if ( u .ge. x(k) ) i = k
!      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline

      CALL XSR_Find_1D_Pointer( N, u, x, i)

!
   30  dx = u - x(i)

!   IF(u.eq.0) THEN 
!       write(*,*) 'i=', i
!       write(*,*) 'x(i) =', x(i) , 'dx=', u - x(i)
!       write(*,*) 'sp_coeff(1,i,0:3,i_ftype)=', sp_coeff(1,i,0:3,i_ftype)
!       pause
!   END IF 

!   IF(u.eq.60) THEN 
!       write(*,*) 'i=', i
!       write(*,*) 'x(i) =', x(i) , 'dx=', u - x(i)
!       write(*,*) 'sp_coeff(1,i,0:3,i_ftype)=', sp_coeff(1,i,0:3,i_ftype)
!       pause
!   END IF 

!      write(*,*) 'spline_eval i=', i
!      pause

      DO im = 1, M
!        value(im) = y(im,i) + dx*(b(im,i) + dx*(c(im,i) + dx*d(im,i)))
! SP_COEFF(NNY,NNROW,0:3)
        value(im) = sp_coeff(im,i,0,i_ftype) + &
           dx*(sp_coeff(im,i,1,i_ftype) + &
           dx*(sp_coeff(im,i,2,i_ftype) + &
           dx*sp_coeff(im,i,3,i_ftype)))
      END DO
      return

      END function spline_eval 

      function spline_eval_der(m, u, i_ftype) result(value_der)
!     SP_COEFF goes through MODULE
      INTEGER, INTENT(IN) :: m, i_ftype
      INTEGER             :: N
! N   -  number of knots
! N-1 -  number of intervals
! M   -  number of y 
      REAL(dp), INTENT(IN) ::  u
!      REAL(dp), INTENT(IN) ::  y(m,n-1), b(m,n-1), c(m,n-1), d(m,n-1)
      REAL(dp) :: value_der(M), X(NNKNOT)
      integer i, im
      REAL(dp) ::  dx
!      data i/1/

     i = 1

      IF( M .GT. NY_SP(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Error! M .GT. NY_SP(i_ftype) :', M,NY_SP(i_ftype)
        stop
      ELSE IF ( M .LT. NY_SP(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Warning! M .LT. NY_SP(i_ftype) :', M,NY_SP(i_ftype)
      END IF

      N = NKNOT(i_ftype)

      X(1:N) = X_SP(1:N, i_ftype)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!

! Binary search should be changed in the fu 
!      SAVE :: i
!      if ( i .ge. n ) i = 1
!      if ( u .lt. x(i) ) go to 10
!      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
!   10 i = 1
!      j = n+1
!   20 k = (i+j)/2
!      if ( u .lt. x(k) ) j = k
!      if ( u .ge. x(k) ) i = k
!      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline

      CALL XSR_Find_1D_Pointer( N, u, x, i)

!
   30 dx = u - x(i)

!       write(*,*) 'dx, u, x = ', dx, u, x(i)
!       pause

!   IF(u.eq.0) THEN 
!       write(*,*) 'i=', i
!       write(*,*) 'x(i) =', x(i) , 'dx=', u - x(i)
!       write(*,*) 'sp_coeff(1,i,0:3,i_ftype)=', sp_coeff(1,i,0:3,i_ftype)
!       pause
!   END IF 

!   IF(u.eq.60) THEN 
!       write(*,*) 'i=', i
!       write(*,*) 'x(i) =', x(i) , 'dx=', u - x(i)
!       write(*,*) 'sp_coeff(1,i,0:3,i_ftype)=', sp_coeff(1,i,0:3,i_ftype)
!       pause
!   END IF 

!      write(*,*) 'spline_eval i=', i
!      pause

      DO im = 1, M
!        value(im) = y(im,i) + dx*(b(im,i) + dx*(c(im,i) + dx*d(im,i)))
! SP_COEFF(NNY,NNROW,0:3)
        value_der(im) =  sp_coeff(im,i,1,i_ftype) + dx*(             &
          2.*sp_coeff(im,i,2,i_ftype)+3.*dx*sp_coeff(im,i,3,i_ftype))
      END DO
      return

      END function spline_eval_der


      Subroutine XSR_Find_1D_Pointer( NP, x, End_Int, int_pointer)
!=====================================================================!
!      Finding Pointers  in a one-dimensional table                   !
!   (c) Slava 27.V.1999                                               !
!=====================================================================!
      implicit none
! Input : burnup, history_void,
      integer, INTENT(IN)  :: NP ! number of points (number of intervals = NP - 1)
      REAL(dp), INTENT(IN) :: End_Int(NP)
      REAL(dp), INTENT(IN) :: x ! Value
!     REAL(dp) BENTRY(ie,  n) - burnup entries
!     REAL(dp) CENTRY(ic,  n) - history void entries
! Output: 
      INTEGER, INTENT(OUT) :: int_pointer
! Local Variables
      LOGICAL Inside_Int

      int_pointer = 0
      Inside_Int  =  .False.

      do while ( .NOT. (Inside_int) )
  
         int_pointer = int_pointer + 1
         if( (int_pointer + 1) .EQ. NP ) then
           Inside_Int = .True. ! Extrapolation
         else
           Inside_Int  =  (x .LE. End_Int(int_pointer + 1) )
         end if         

      end do

      return
      end    Subroutine XSR_Find_1D_Pointer   


      END MODULE spline_lib   

      MODULE POLYNOMIAL_ARRAY

      CONTAINS 

      FUNCTION LEGENDRE_ARRAY_old(n,x) RESULT (pol)
      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n) :: pol

      REAL(dp) :: pol0, pol1
      INTEGER  :: m

      IF( n .LT. 0) THEN
        WRITE(*,*) 'Computing Legendre Polynomials with negative power'
        STOP
      END IF

      pol0 = 1.
      pol(0) = 1.
      IF( n == 0) RETURN
         pol1 = x
         pol(1) = x
      IF( n == 1) RETURN

      DO m = 2, n
        pol(m) = ((2*m+1)*x*pol1 - m*pol0)/(m+1)
        pol0=pol1
        pol1=pol(m)
      END DO

      RETURN
      END FUNCTION LEGENDRE_ARRAY_old     


      FUNCTION LEGENDRE_ARRAY_NEW(n,x) RESULT (pol)

      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n) :: pol

      REAL(dp) :: pol0, pol1
      INTEGER  :: m

      IF( n .LT. 0) THEN
        WRITE(*,*) 'Computing Legendre Polynomials with negative power'
      STOP
      END IF
 
 
      pol0 = 1.
      pol(0) = 1.
      IF( n == 0) RETURN
        pol1 = x
        pol(1) = x
      IF( n == 1) RETURN

      DO m = 1, n-1
        pol(m+1) = ((2*m+1)*x*pol1 - m*pol0)/(m+1)
        pol0=pol1
        pol1=pol(m+1)
      END DO

      RETURN

      END FUNCTION LEGENDRE_ARRAY_NEW     

      FUNCTION LEGENDRE_ARRAY_DER_NEW(n,x,pol) RESULT (pol_der)

      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n), INTENT(IN) :: pol
      REAL(dp), DIMENSION(0:n) :: pol_der


      REAL(dp) :: pold0, pold1
      INTEGER  :: m

      IF( n .LT. 0) THEN
        WRITE(*,*) 'Computing Legendre Polynomials with negative power'
        STOP
      END IF
 
      pold0 = 0.
      pol_der(0) = 0.
      IF( n == 0) RETURN
        pold1 = 1.
        pol_der(1) = 1.

      IF( n == 1) RETURN
 
      DO m = 1, n-1
         pol_der(m+1) = x * pold1 + (m+1) * pol(m)
         pold0=pold1
         pold1=pol_der(m+1)
      END DO

      RETURN

      END FUNCTION LEGENDRE_ARRAY_DER_NEW


      FUNCTION LEGENDRE_ARRAY_DER_OLD(n,x,pol) RESULT (pol_der)

      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n), INTENT(IN) :: pol
      REAL(dp), DIMENSION(0:n) :: pol_der


      REAL(dp) :: pold0, pold1
      INTEGER  :: m


      IF( n .LT. 0) THEN
       WRITE(*,*) 'Computing Legendre Polynomials with negative power'
       STOP
      END IF

      IF( n .GT. 8) THEN
         WRITE(*,*) 'We can compute the derivative only for n <=8'
         STOP
      END IF

 
      pol_der(0) = 0.
      IF( n == 0) RETURN

      pol_der(1) = 1.
      IF( n == 1) RETURN
 
      DO m = 2, n
        IF ( m == 2 ) pol_der(m) = 10.*x/3.
        IF ( m == 3 ) pol_der(m) = - 23./12. + 35.*(X**2)/4.
        IF ( m == 4 ) pol_der(m) = (7.*x/30.)*(-41.+90*x**2)
        IF ( m == 5 ) pol_der(m) = (103 - 1344*x**2+1925*x**4)/40.
        IF ( m == 6 ) pol_der(m) = 2487*x/140.-506*x**3/5 + 429*x**5/4
        IF ( m == 7 ) pol_der(m) = (-6967+177771*x**2-625625*x**4+ &
                                   525525*x**6)/2240.
        IF ( m == 8 ) pol_der(m) = 11*x*(-25237+265590*x**2        &
                                 -667485*x**4 + 464100*x**6)/10080.
      END DO

      RETURN

      END FUNCTION LEGENDRE_ARRAY_DER_OLD

      FUNCTION CHEBYSHEV_ARRAY(n,x) RESULT (pol)
      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n) :: pol

      REAL(dp) :: pol0, pol1
      INTEGER  :: m

      pol0 = 1.
      pol(0) = 1.
      IF(n==0) RETURN
       pol1 = x
       pol(1) = x
      IF(n==1) RETURN

      DO m = 2, n
         pol(m) = 2.*x*pol1 - pol0
         pol0=pol1
         pol1=pol(m)
      END DO

      RETURN
      END FUNCTION CHEBYSHEV_ARRAY     


      FUNCTION CHEBYSHEV_ARRAY_DER(n,x,pol) RESULT (pol_der)

      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), INTENT(IN), DIMENSION(0:n) :: pol
      REAL(dp), DIMENSION(0:n) :: pol_der

      REAL(dp) :: u0, u1
      INTEGER  :: m

      u0 = 0
      pol_der(0) = 0.
      IF(n==0) RETURN

      DO m = 0, n-1
         u1 = x * u0 + pol(m)
         pol_der(m+1) = (m+1) * u1
         u0 = u1
      END DO

      RETURN

      END FUNCTION CHEBYSHEV_ARRAY_DER


      FUNCTION MONOMIAL_ARRAY(n,x) RESULT (pol)
      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x 
      REAL(dp), DIMENSION(0:n) :: pol

      INTEGER  :: m

      pol(0) = 1.
      DO m = 1, n
          pol(m) = x*pol(m-1)
      END DO

      RETURN
      END FUNCTION MONOMIAL_ARRAY     


      FUNCTION MONOMIAL_ARRAY_DER(n,pol) RESULT (pol_der)
      USE PRECISION, ONLY : dp
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN), DIMENSION(0:n) :: pol
      REAL(dp), DIMENSION(0:n) :: pol_der

      INTEGER  :: m

      pol_der(0) = 0.
      IF( n == 0) RETURN
        pol_der(1) = 1.

      DO m = 2, n
        pol_der(m) = m * pol(m-1) ! /(m-1)*x*pol_der(m-1)
      END DO

      RETURN
      END FUNCTION MONOMIAL_ARRAY_DER

      END MODULE POLYNOMIAL_ARRAY

      module termsort_lib

      USE PRECISION, ONLY : dp, eps_round_off, real_min_value
      USE XS_POL_FTYPES

      implicit none

! private
      public

      integer, parameter :: max_nx=10, max_ny=100

      integer, parameter :: max_total_term=1000,max_term=200

      integer, parameter :: max_1d_power=10

      integer n_1d_power(max_nx, NN_FTYPE)
      real(dp) :: pol_xxx( 0:max_1d_power) 
      real(dp) :: pol_1d_value( max_nx, 0:max_1d_power) 
      real(dp) :: pol_1d_dvalue( max_nx, 0:max_1d_power) 

      CHARACTER     :: pol_type(NN_FTYPE)
      CHARACTER*6   :: flag_spline(NN_FTYPE)

 
      integer n_y(NN_FTYPE) ! количество полиномов
      integer n_x(NN_FTYPE) ! количество независимых переменных
      integer n_total_term(NN_FTYPE) ! количество слагаемых во всех полиномах

! integer count
      integer pol_power( max_nx, max_total_term, NN_FTYPE)
      real(dp) :: xmin_scale( max_nx, NN_FTYPE ), &
                        xmax_scale( max_nx, NN_FTYPE )
      real(dp) pol_term_value(max_total_term)
      real(dp) pol_term_dvalue(max_total_term, max_nx)
      integer n_term(max_ny, NN_FTYPE)
      real(dp) pol_coeff(max_ny,max_term, NN_FTYPE)
      integer index_term(max_ny,max_term, NN_FTYPE)
      CHARACTER*18  x_title(max_nx, NN_FTYPE),y_title(max_ny, NN_FTYPE)

      real(dp) :: x_scale(max_nx)

      PUBLIC:: init_lib,ops_lib, write_n_term_lib, write_termref_lib, &
            scale_back, get_lib_hm_mass_of_fuel

      integer,  public :: n_x_der(NN_FTYPE, max_nx) ! производные, которые нужно взять в зависимости от типа кассеты и индекса xi
      integer,  public :: nxder_xxx(max_nx) ! вспомогательный массив
      REAL(dp), public :: y_der(MAX_NX, MAX_NY)

! Debug
      PUBLIC:: lib_flag_spline, get_lib_titles, get_lib_nx_ny
!  , y_title

      contains

      function array_sum(v,dim) result(summa)

      implicit none
      integer :: i, dim
      real(dp) :: v(dim), summa

      summa = 0
      do i=1,dim
        summa = summa + v(i)
      enddo

      return
      end function


      function term_eq(a1,a2, n_x) result(eqv)
      logical eqv
      integer n_x 
      integer, intent(in):: a1(n_x),a2(n_x)
      integer i
      do i=1,n_x
         if( a1(i).ne.a2(i) )then
          eqv=.false.
          return
         endif
      enddo
      eqv=.true.
      end function term_eq

! subroutine term_rd(id,a1)
! integer id,i
! type(term) :: a1
! return
! end subroutine term_rd


      function term_val(n_pol_power, n_x, nder) result( value )
      real(dp) value
      integer :: n_x, nder
      integer :: n_pol_power(1:n_x)
      integer i
 
      value=1.

      do i=1,n_x
        if (i == nder) then
           value=value*pol_1d_dvalue( i, n_pol_power(i) )
!      write(*,*) 'i, pow, pol = ', i, n_pol_power(i), pol_1d_dvalue( i, n_pol_power(i) )
       else
           value=value*pol_1d_value( i, n_pol_power(i) )     
       endif
      enddo

      return
      end function term_val


      function array_norm(v,dim) result (norm)

      implicit none

      integer :: i, dim
      real(dp) :: v(dim)
      real(dp) :: norm

      norm = 0
      do i=1,dim
         norm = norm + v(i) ** 2
      enddo
      norm = norm ** 0.5

      return
      end function array_norm


      subroutine ops_lib(NX_IN, NY_IN, x, i_ext_ftype, y, y_der_)

      USE POLYNOMIAL_ARRAY !  functions chebyshev_array, legendre_array, monomial_array
      USE SPLINE_LIB, ONLY : spline_eval, spline_eval_der

      logical :: ifout ! выходит ли точка за границу допустимой области
!      logical :: ifder ! нужно ли считать производные внутри допустимой области
      INTEGER, INTENT(IN)  :: NX_IN, NY_IN
      integer :: ort(nx_in) ! направление, по которому входная точка выходит за границу допустимой области

      real(dp), INTENT(IN) :: x(NX_IN)
      integer,  INTENT(IN) :: i_ext_ftype
! character,INTENT(IN) :: polynomial_type
      real(dp), INTENT(OUT):: y(NY_IN)            ! macrosections
      real(dp), INTENT(OUT), OPTIONAL:: y_der_(NX_IN, NY_IN) ! derivatives of macrosections

      real(dp) :: Xb(NX_IN), Xb_scale(NX_IN), alpha(NX_IN)
      real(dp) :: dx_scale(nx_in), dx_scale_norm
      real(dp) :: scalar, dx_norm, dX(NX_IN)

! Local
      REAL(dp) :: sum, ywrap
      integer i,j, n, num, i_ftype, dim

! EXTERNAL :: chebyshev_array, legendre_array, monomial_array
!  chebyshev_array, legendre_array, monomial_array
!    write(*,*) '========================= dX(NX_IN)',NX_IN
! Changing external fuel type numbering to internal
      i_ftype = index_internal(i_ext_ftype)
!------------------------------------------------------------------------------!
! First of all cheking NX_IN <= n_x(i_ftype) , NY_IN <= n_y(i_ftype)
!------------------------------------------------------------------------------!
      IF( NX_IN .GT. N_X(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Error! NX_IN .GT. N_X(i_ftype):',NX_IN,N_X(i_ftype)
        stop
      ELSE IF ( NX_IN .LT. N_X(i_ftype) ) THEN
        write(*,*) 'i_ftype =', i_ftype
        write(*,*) 'Warning! NX_IN .LT. N_X(i_ftype) :', &
                                                   NX_IN,  N_X(i_ftype)
      END IF

      IF( NY_IN .GT. N_Y(i_ftype) ) THEN
         write(*,*) 'i_ftype =', i_ftype
         write(*,*) 'Error! NY_IN .GT. N_Y(i_ftype) :', &
                                                   NY_IN,  N_Y(i_ftype)
         stop
      ELSE IF ( NY_IN .LT. N_Y(i_ftype) ) THEN
         write(*,*) 'i_ftype =', i_ftype
         write(*,*) 'Warning! NY_IN .LT. N_Y(i_ftype) :', &
                                                   NY_IN,  N_Y(i_ftype)
      END IF

!------------------------------------------------------------------------------!
! First of all scaling X
!------------------------------------------------------------------------------!
       n_x_der(i_ftype, :) = 0
       CALL SCALE_X(NX_IN, i_ftype, x, x_scale)
       CALL REGIONOUT(NX_IN, X_scale, Xb_scale, ort, ifout)

!       if( IFOUT )  then
!      WRITE(*,*) 'i_ext_ftype =', i_ext_ftype
!        WRITE(*,*) 'NX_IN =',  NX_IN
!      WRITE(*,*) 'X_scale= ', X_scale(1:nx_in)
!      WRITE(*,*) 'XB_scale= ', XB_scale(1:nx_in)
!        write(*,*) 'X =', x(1:nx_in)
!        pause 
!       END IF       


       n_x_der(i_ftype,1:nx_in) = ort(:)

!  write(*,*) 'x_scale =', x_scale(1:n_x(i_ftype))
!  pause
!------------------------------------------------------------------------------!
! Computing 1D polynomial terms
!------------------------------------------------------------------------------!
! write(*,*) 'nx, xb_scale: ', NX_IN, xb_scale

      IF( pol_type(i_ftype).eq. "M") THEN
       DO n = 1, NX_IN
          pol_1d_value(n, 0:n_1d_power(n,i_ftype) ) =  &
           monomial_array( n_1d_power(n,i_ftype), xb_scale(n) )
       END DO
      ELSE IF (pol_type(i_ftype) .eq. "C") THEN
       DO n = 1, NX_IN
       pol_1d_value(n, 0:n_1d_power(n,i_ftype) ) =     &
           chebyshev_array( n_1d_power(n,i_ftype), xb_scale(n) )
       END DO
      ELSE IF (pol_type(i_ftype).eq. "L") THEN
       DO n = 1, NX_IN
         pol_1d_value(n, 0:n_1d_power(n,i_ftype) ) =    &
           legendre_array_old( n_1d_power(n,i_ftype), xb_scale(n) )
       END DO
      ELSE IF (pol_type(i_ftype).eq. "P") THEN
       DO n = 1, NX_IN
         pol_1d_value(n, 0:n_1d_power(n,i_ftype) ) =  &
           legendre_array_new( n_1d_power(n,i_ftype), xb_scale(n) )
       END DO
      ELSE
       write(*,*) 'FUEL TYPE =', i_ftype
       write(*,*) 'polynomial_type=', pol_type(i_ftype)
       write(*,*) 'polynomial_type can be M, C or L, P'
       stop
      END IF

! computing multidimensional polynomial terms
       do i=1,n_total_term(i_ftype)
         pol_term_value(i)=term_val(pol_power(1,i,i_ftype), NX_IN, 0)
       enddo


! write(*,*) 'pol_term_value(i)=', pol_term_value(1:5)
! pause

       do i=1, NY_IN
          sum=0.
          do j=1,n_term(i,i_ftype)
          sum=sum+pol_coeff(i,j,i_ftype)*   &
                    pol_term_value(index_term(i,j,i_ftype))

!   if (i==2) write(*,*) 'coef, term: ', pol_coeff(i,j,i_ftype), pol_term_value(index_term(i,j,i_ftype))
          enddo
          y(i)=sum
       enddo

!  write(*,*) 'y = ',y


        IF( INDEX(flag_spline(i_ftype), "SPLINE") /=0 ) THEN
           y(1:NY_IN) = y(1:NY_IN) + spline_eval(NY_IN, x(1), i_ftype) 
!    WRITE(*,*) 'Spline =', spline_eval(NY_IN, x(1), i_ftype) 
!    pause
        END IF 


!              write(*,*)  'n_x_der(i_ftype, 1:NX_IN)', n_x_der(i_ftype, 1:NX_IN)

!         pause


        DO n = 1, NX_IN
           num = n_x_der(i_ftype,n) !номер индекса производной
           if ( (num /= 0).or.( PRESENT(y_der_) )) then 
             dim = n_1d_power(n,i_ftype)
             pol_xxx(0:dim) = pol_1d_value(n, 0:dim)

            SELECT CASE (pol_type(i_ftype))  
                 CASE("M") 
                  pol_1d_dvalue(n, 0:dim ) = monomial_array_der &
                                               ( dim, pol_xxx )
                 CASE("C") 
                  pol_1d_dvalue(n, 0:dim ) = chebyshev_array_der &
                                        ( dim, xb_scale(n), pol_xxx )
                 CASE("L") 
                  pol_1d_dvalue(n, 0:dim ) = legendre_array_der_old &
                                        ( dim, xb_scale(n), pol_xxx )
                 CASE("P") 
                  pol_1d_dvalue(n, 0:dim ) = legendre_array_der_new &
                                        ( dim, xb_scale(n), pol_xxx )
                 CASE DEFAULT 
                  write(*,*) 'FUEL TYPE =', i_ftype
                  write(*,*) 'polynomial_type=', pol_type(i_ftype)
                  write(*,*) 'polynomial_type can be M, C or L, P'
                  stop
            END SELECT

          END IF
        END DO

    ! computing derivative of multidimensional polynomial terms
        do i=1,n_total_term(i_ftype)
           do n = 1, NX_IN
              num = n_x_der(i_ftype,n) !номер индекса производной
              if ( (num /= 0).or.( PRESENT(y_der_) ) ) then 
                pol_term_dvalue(i, n) = term_val  &
                              (pol_power(1,i,i_ftype), NX_IN, n)
!             write(*,*) 'dterm(i,n): ',i,n,pol_term_dvalue(i,n)
              endif
           enddo
       enddo

    ! sum of terms
       do n = 1, NX_IN ! по направлениям
        do i=1, NY_IN ! по элементам выходного вектора y
            sum=0.
            do j=1,n_term(i,i_ftype) ! по членам полинома
                num = n_x_der(i_ftype,n) !номер индекса производной
                if ( (num /= 0).or.( PRESENT(y_der_) ) ) then 
                    sum=sum+pol_coeff(i,j,i_ftype)* &
                          pol_term_dvalue(index_term(i,j,i_ftype), n)
                end if
            enddo
            y_der(n,i)=sum 
        enddo
       enddo


       
 
       IF ( ifout ) THEN

!       write(*,*) 'ifout =', ifout
!       PAUSE
!      WRITE(*,*) 'i_ext_ftype=', i_ext_ftype


!       write(*,*) 'we are inside extrapolation =', ifout

         CALL SCALE_BACK(NX_IN, i_ext_ftype, Xb_scale, Xb)

      ! ВЫЧИСЛЯЕМ ВЕКТОР dX И ЕГО НОРМУ
         dX_norm = 0
         dX(1:NX_IN) = 0.

         DO n=1,NX_IN

            num = n_x_der(i_ftype,n) !номер индекса производной

            if (num /= 0) then 

                dX(n) = X(n) - Xb(n)
                dX_norm = dx_norm + dX(n) ** 2

            end if

         END DO

         dx_norm = sqrt(dx_norm)

!         write(*,*) 'xb_scale(1:NX_IN)', xb_scale(1:NX_IN)
!         write(*,*) 'x_scale(1:NX_IN) ', x_scale(1:NX_IN)


!        dx_scale(:) = x_scale(:) - xb_scale(:)
 
        dx_scale(1:NX_IN) = x_scale(1:NX_IN) - xb_scale(1:NX_IN)
         dx_scale_norm = array_norm(dx_scale,nx_in)
         do i =1,NX_IN      
           alpha(i) = 1 - 2*acos(abs(dx_scale(i))/dx_scale_norm)/  &
                                                         3.14159265
         enddo

         do j=1,ny_in
          ywrap = 0
          do i=1,nx_in
            ywrap = ywrap + y_der(i,j) * alpha(i)
          enddo
          y_der(1:nx_in,j) = ywrap/array_sum(alpha,nx_in) * &
                           abs(dx_scale(1:nx_in))/dx_scale_norm
         enddo

      ! ВЫЧИСЛЯЕМ ЭКСТРАПОЛИРОВАННОЕ ЗНАЧЕНИЕ Y
         DO i=1,NY_IN
          DO n=1,NX_IN
                y(i) = y(i) + y_der(n,i) * dx_scale(n)
          END DO
         END DO

!        write(*,*) 'we are inside extrapolation' 
        END IF

        nxder_xxx(1:NX_IN) = n_x_der(i_ftype,1:NX_IN)  
        CALL SCALE_DER_X(nxder_xxx, i_ftype, NX_IN, NY_IN, &
                                           PRESENT(y_der_) ) !y_der(1:NX_IN, 1:NY_IN)

!  write(*,*)  'y_der = ', y_der(2,10)
!  write(*,*)  'y, y_der, dx = ', x(2), xb(2), y(10), y_der(2,10), dX(2)

!        IF ( PRESENT(y_der_) ) THEN
!          DO n=1,NX_IN
!               write(*,*) 'n =', n,  'Y_DER(n,10)=', Y_DER(n,10)
!          end do
!        END IF


        IF ( PRESENT(y_der_) ) &
                     y_der_(1:NX_IN, 1:NY_IN) = y_der(1:NX_IN, 1:NY_IN)

        IF ( PRESENT(y_der_) ) THEN
         IF( INDEX(flag_spline(i_ftype), "SPLINE") /=0 ) THEN
                   y_der_(1, 1:NY_IN) = y_der(1, 1:NY_IN) + &
                                spline_eval_der(NY_IN, x(1), i_ftype) 
         END IF
        END IF           

!        IF ( PRESENT(y_der_) ) THEN
!          DO n=1,NX_IN
!               write(*,*) 'n =', n,  'Y_DER(n,10)=', Y_DER_(n,10)
!          end do
!        END IF


!  write(*,*) 'xdata: ', xb_scale(1:NX_IN)!, y_der(3,5)
!  write(*,'(5F20.6,E15.7)') x(1:NX_IN), y_der(3,2)
!  pause

        return
        end subroutine ops_lib
 
        subroutine add_poly(id, iy)
        integer id, iy
        integer i,j
        integer :: a1(max_nx)
        logical flag_found

! write(*,*) 'iy =', iy, 'n_ftype  =', n_ftype
         read(id,*) n_term(iy, n_ftype)
! write(*,*) 'n_term(iy, n_ftype)=', n_term(iy, n_ftype)
         read(id,*) (pol_coeff(iy,j, n_ftype),j=1,n_term(iy,n_ftype))
         do i=1,n_term(iy, n_ftype)

           read(id,*) a1(1:n_x(n_ftype))
   
           index_term(iy,i, n_ftype) = 0
           flag_found = .False.
           do j=1,n_total_term(n_ftype)
              if( term_eq(a1, pol_power(1, j, n_ftype), n_x(n_ftype))) &
                                                                   then 
                 index_term(iy, i, n_ftype)=j
                 flag_found =.True.
                 exit
              endif
           enddo

!  if(j>n_total_term)then
!  if( index_term(n_y,i) == 0 ) THEN
          if( flag_found .eqv. .False. ) THEN
            n_total_term(n_ftype) = n_total_term(n_ftype)+1
            index_term(iy,i,n_ftype)=n_total_term(n_ftype)
            pol_power(1:n_x(n_ftype), n_total_term(n_ftype),n_ftype)= &
                                                    a1(1:n_x(n_ftype))
         endif

         enddo

! write(*,*) 'iy =', iy, 'n_ftype=', n_ftype  
! write(*,*) 'n_term(iy, n_ftype)=', n_term(iy, n_ftype)
! write(*,*) 'n_total_term(n_ftype)=', n_total_term(n_ftype)
! write(*,*) 'index_term(iy,i,n_ftype)=',  &
!       index_term(iy,1:n_term(iy, n_ftype),n_ftype)
! pause

         end subroutine add_poly

         subroutine init_lib(id)
         USE SPLINE_LIB
         integer id
         INTEGER iy, i, n, ios, i_ftype

         N_FTYPE = 1

         DO 
!
           read (id, FMT=*, IOSTAT=ios) index_ftype(n_ftype), &
                                        hm_mass_of_fuel(n_ftype) 
!   pause
           IF( ios /= 0 ) THEN 
               n_ftype = n_ftype - 1
               EXIT
           END IF
!           write(*,*) 'index_ftype(n_ftype)=', index_ftype(n_ftype)

           index_internal(index_ftype(n_ftype))=n_ftype
           read(id,*) ! header
           read(id,'(A)') title_ftype(n_ftype)
           write(*,'(A)') TRIM(title_ftype(n_ftype))
           read(id,*) ! header

           read(id,*) pol_type(n_ftype)
           write(*,*) 'pol_type(n_ftype)=', pol_type(n_ftype)

           read (id,*) n_x(n_ftype), n_y(n_ftype)! количество x! количество сечений

!           write(*,*)  'x_title in'
           read (id,'(5A18)') (x_title(i,n_ftype), i=1,n_x(n_ftype))
           read (id,'(5A18)') (y_title(i,n_ftype), i=1,n_y(n_ftype))

!           write(*,*)  'y_title out'


           read(id,*) xmin_scale(1:n_x(n_ftype), n_ftype)
           read(id,*) xmax_scale(1:n_x(n_ftype), n_ftype)

!           write(*,*) 'xmin_scale(1:n_x(n_ftype), n_ftype)', &
!                     xmin_scale(1:n_x(n_ftype), n_ftype)
!           write(*,*) 'xmax_scale(1:n_x(n_ftype), n_ftype)', &
!                     xmax_scale(1:n_x(n_ftype), n_ftype)
!   pause
!
           n_total_term(n_ftype)=0 !
           do iy=1, n_y(n_ftype)
!              write(*,*) 'add_poly in'
              call add_poly(id, iy)
!              write(*,*) 'add_poly out'
           enddo

           read(id,'(A)') flag_spline(n_ftype)

           IF( INDEX(flag_spline(n_ftype),"SPLINE") /= 0 ) THEN
!              write(*,*) 'spline in'
               CALL SPLINE_INPUT(id, n_ftype)
!              write(*,*) 'spline out'
           END IF

           n_ftype = n_ftype + 1

       END DO ! input

!       write(*,*) 'input of ', n_ftype, ' XS types'
!       pause
  
!input finished

         DO i_ftype=1, n_ftype
            DO iy = 1, n_y(i_ftype)
!   write(*,*) 'iy, n_total_term', iy, n_term(iy, i_ftype)
!   write(*,'(20I4)') index_term(iy,1:n_term(iy,i_ftype),i_ftype)
            END DO

            n_1d_power(:,i_ftype) = 0
            DO i = 1, n_total_term(i_ftype)
               DO n = 1, n_x(i_ftype)
               n_1d_power(n,i_ftype) = MAX( pol_power(n, i,i_ftype),  &
                                               n_1d_power(n,i_ftype) )
               END DO
            END DO
 
! WRITE(*,*) 'maximum polynomial power ',  n_1d_power(1:n_x(i_ftype), i_ftype)     
       end do ! n_ftype

       RETURN
       end subroutine init_lib

       subroutine write_n_term_lib(io_unit, i_ext_ftype)
       INTEGER, INTENT(IN):: io_unit, i_ext_ftype
       integer n, i_ftype

       i_ftype = index_internal(i_ext_ftype)
 
       DO n = 1, n_y(i_ftype)
           write(io_unit, '(1x,2I5)') n, n_term(n,i_ftype)
       END DO 

       return 
       end subroutine write_n_term_lib

       subroutine write_termref_lib(io_unit, i_ext_ftype)
! integer index_term(max_ny,max_term)
       INTEGER, INTENT(IN):: io_unit, i_ext_ftype
       integer n, j, i_ftype

       i_ftype = index_internal(i_ext_ftype)
       DO n = 1, n_y(i_ftype)
           do j = 1, n_term(n, i_ftype)
           write(io_unit, '(1x,2I5)') n, index_term(n,j,i_ftype)
           end do
       END DO 
       end subroutine write_termref_lib

       SUBROUTINE SCALE_X(NX, i_ftype, X, x_scale)
!=========================================================================!
! Scale of the original data, XMAX, XMIN are given                        !
!=========================================================================!

! xmin_scale( max_nx, NN_FTYPE ), xmax_scale( max_nx, NN_FTYPE ) 
       INTEGER, INTENT(IN) :: NX, i_ftype
       REAL (dp), INTENT(IN) :: X(1:NX)
       REAL (dp), INTENT(OUT) :: x_scale(1:NX)

       INTEGER :: i
       REAL(dp)    :: bb, aa, delta_xmax_xmin

       DO i = 1, NX
!      write(*,*) 'i,i_ftype=', i,i_ftype
!      write(*,*) 'xmin_scale, xmax_scale=', xmin_scale(i,i_ftype), &
!                    xmax_scale(i,i_ftype)
        IF( ABS(xmax_scale(i,i_ftype)) .GT. REAL_MIN_VALUE ) THEN
            delta_xmax_xmin =   &
                ABS((xmax_scale(i,i_ftype)-xmin_scale(i,i_ftype)) &
                                         /xmax_scale(i,i_ftype) )
        ELSE
            delta_xmax_xmin = ABS( xmax_scale(i,i_ftype) )
        END IF    
        IF( delta_xmax_xmin .GT. EPS_ROUND_OFF ) THEN
           aa = 2./(xmax_scale(i,i_ftype) - xmin_scale(i,i_ftype))
           bb = -  (xmax_scale(i,i_ftype) + xmin_scale(i,i_ftype) )/ &
               ( xmax_scale(i,i_ftype) - xmin_scale(i,i_ftype))
           X_SCALE(i)= aa*X(i)+bb
        ELSE
           X_SCALE(i)= X(i)
        END IF
       END DO

       END SUBROUTINE SCALE_X

       SUBROUTINE REGIONOUT(NX, X_scale, Xb_scale, ort, ifout)
 ! проверяем, лежит ли точка внутри допустимой области

        logical :: ifout ! =True если точка выходит за хотя бы одну границу
        INTEGER , INTENT(IN) :: NX
        REAL (dp),INTENT(IN) :: X_scale(1:NX)  ! входная точка
        REAL (dp),INTENT(OUT):: Xb_scale(1:NX) ! точка, лежащая внутри или на границе области
        INTEGER, INTENT(OUT) :: ort(NX)      ! направление, по которому точка выходит за пределы области
!    REAL (dp) :: xmax, xmin
        REAL (dp) :: xo
        INTEGER :: i = 0

        Xb_scale(1:NX) = X_scale(1:NX)
        ort(:) = 0
        xo = 1.
        ifout = .False.

        DO i=1,NX 

        IF (abs(X_scale(i)) .gt. 1.) THEN
            xo = abs(X_scale(i))
            Xb_scale(i) = X_scale(i) / xo            
            ort(i) = 1
            ifout = .True.
        END IF 

        END DO


        END SUBROUTINE REGIONOUT

!CALL SCALE_DER_X(y_der, nxder_xxx, i_ftype, NX_IN, NY_IN)
        SUBROUTINE SCALE_DER_X(nxder_xxx, i_ftype, NX, NY, ifder)
!=========================================================================!
! Scale of the original data for derives, XMAX, XMIN are given                        !
!=========================================================================!

        logical :: ifder
        INTEGER, INTENT(IN) :: NX, NY, i_ftype, nxder_xxx(1:NX)
!  REAL (dp) :: Y_DER(1:NX, 1:NY)

        INTEGER :: i, k
        REAL(dp)    :: aa, delta_xmax_xmin

        DO k = 1, NX
         IF( ABS(xmax_scale(k,i_ftype)) .GT. REAL_MIN_VALUE ) THEN
          delta_xmax_xmin = &
          ABS( (xmax_scale(k,i_ftype) - xmin_scale(k,i_ftype)) &
           /xmax_scale(k,i_ftype) )
         ELSE
           delta_xmax_xmin = ABS( xmax_scale(k,i_ftype) )
         END IF    

        IF( delta_xmax_xmin .GT. EPS_ROUND_OFF ) THEN
           aa = 2./(xmax_scale(k,i_ftype) - xmin_scale(k,i_ftype))
        ENDIF

        IF ( (nxder_xxx(k) /= 0).or.(ifder)) THEN
            DO i = 1, NY
!               IF( i == 1 ) THEN 
!               WRITE(*,*) 'k=', k, 'i_ftype=', i_ftype
!               WRITE(*,*) 'delta_xmax_xmin=', delta_xmax_xmin
!               WRITE(*,*) 'ABS(xmax_scale(k,i_ftype)) =', &
!                         ABS(xmax_scale(k,i_ftype))
!               WRITE(*,*) 'ABS(xmin_scale(k,i_ftype)) =', &
!                         ABS(xmin_scale(k,i_ftype))
!               write(*,*) 'k =', k, 'aa=', aa
!             WRITE(*,*) 'Y_DER(k,i)=', Y_DER(k,i)
!               end if
                Y_DER(k,i)= aa * Y_DER(k,i)
!                if((i==1 .or. i==26).and.(k==2)) &
!                     write(*,*) 'i, y_der, aa: ',i,y_der(k,i), aa
            END DO
        END IF

       END DO


!      write(*,*) 'we are inside scale back'

      END SUBROUTINE SCALE_DER_X


      SUBROUTINE SCALE_BACK(NX, i_ext_ftype, ksi, x)
!=========================================================================!
! External procedure i_ext_ftype  
! recover the original data from the scaled, XMAX, XMIN are given                        !
!=========================================================================!

! xmin_scale( max_nx, NN_FTYPE ), xmax_scale( max_nx, NN_FTYPE ) 
      INTEGER, INTENT(IN) :: NX, i_ext_ftype
      REAL (dp), INTENT(IN) :: ksi(1:NX)
      REAL (dp), INTENT(OUT) :: x(1:NX)

      INTEGER :: i, i_ftype
      REAL(dp)    :: bb, aa, delta_xmax_xmin

      i_ftype = index_internal(i_ext_ftype)
      DO i = 1, NX
        IF( ABS(xmax_scale(i,i_ftype)) .GT. REAL_MIN_VALUE ) THEN
           delta_xmax_xmin = &
                ABS( (xmax_scale(i,i_ftype) - xmin_scale(i,i_ftype)) &
                 /xmax_scale(i,i_ftype) )
        ELSE
           delta_xmax_xmin = ABS( xmax_scale(i,i_ftype) )
        END IF    
        IF( delta_xmax_xmin .GT. EPS_ROUND_OFF ) THEN
!      write(*,*) 'i,i_ftype=', i,i_ftype
!      write(*,*) 'xmin_scale, xmax_scale=', xmin_scale(i,i_ftype), &
!                    xmax_scale(i,i_ftype)
           aa =  xmax_scale(i,i_ftype) - xmin_scale(i,i_ftype)
           bb =  xmax_scale(i,i_ftype) + xmin_scale(i,i_ftype) 
           X(i)= 0.5*(aa*ksi(i)+bb)
        ELSE
           X(i)= ksi(i)
        END IF
       END DO

      END SUBROUTINE SCALE_BACK

      LOGICAL FUNCTION  lib_flag_spline(i_ext_ftype)
      INTEGER, INTENT (IN) :: i_ext_ftype
      INTEGER              :: i_ftype 

      i_ftype = index_internal(i_ext_ftype)
 
      IF( flag_spline(i_ftype)=="SPLINE" ) THEN
          lib_flag_spline = .True.
      ELSE
          lib_flag_spline = .False.
      END IF     

      RETURN
      END FUNCTION lib_flag_spline 


      SUBROUTINE get_lib_nx_ny(i_ext_ftype, NX, NY)
      INTEGER, INTENT (IN) :: i_ext_ftype
      INTEGER, INTENT(OUT) :: NX, NY
      INTEGER :: i_ftype

!      write(*,*) 'i_ext_ftype=', i_ext_ftype
!      write(*,*) 'index_internal(8)=', &&
!               index_internal(i_ext_ftype)
! pause
        i_ftype = index_internal(i_ext_ftype)
!       write(*,*) 'i_ftype =', i_ftype, 'i_ext_ftype=', i_ext_ftype
!       write(*,*) 'n_x =', n_x(i_ftype), 'n_y=', n_y(i_ftype)
        NX = n_x(i_ftype)
        NY = n_y(i_ftype)

      RETURN
      END SUBROUTINE get_lib_nx_ny   

      SUBROUTINE get_lib_titles(i_ext_ftype, title_x, title_y)
!=========================================================================!
! External procedure i_ext_ftype  
! recover the original data from the scaled, XMAX, XMIN are given                        !
!=========================================================================!

! xmin_scale( max_nx, NN_FTYPE ), xmax_scale( max_nx, NN_FTYPE ) 
        INTEGER, INTENT(IN)       :: i_ext_ftype
        CHARACTER*18, INTENT(OUT) :: title_x(*), title_y(*)

        INTEGER :: i_ftype

        i_ftype = index_internal(i_ext_ftype)

!        write(*,*) 'i_ext_ftype,i_ftype=', i_ext_ftype,i_ftype

!        write(*,*) 'x_title(1:n_x(i_ftype), i_ftype)=', &
!               x_title(1:n_x(i_ftype), i_ftype)
          
!        write(*,*) 'n_x(i_ftype)', n_x(i_ftype)
          
        title_x(1:n_x(i_ftype)) = x_title(1:n_x(i_ftype), i_ftype)
        title_y(1:n_y(i_ftype)) = y_title(1:n_y(i_ftype), i_ftype)

         


      END SUBROUTINE get_lib_titles


      SUBROUTINE get_lib_hm_mass_of_fuel(i_ext_ftype, hm_mass)
!=========================================================================!
! External procedure i_ext_ftype  
! recover the original data from the scaled, XMAX, XMIN are given                        !
!=========================================================================!

! xmin_scale( max_nx, NN_FTYPE ), xmax_scale( max_nx, NN_FTYPE ) 
        INTEGER, INTENT(IN)       :: i_ext_ftype
        REAL(dp), INTENT(OUT)     :: hm_mass

        INTEGER :: i_ftype

        i_ftype = index_internal(i_ext_ftype)

        hm_mass = hm_mass_of_fuel(i_ftype)

!  IF( i_ext_ftype > 10 ) THEN
!  write(*,*) 'i_ext_ftype, i_ftype =', i_ext_ftype, i_ftype
!  write(*,*) 'hm_mass_of_fuel(i_ftype)', hm_mass_of_fuel(i_ftype)
!  read(*,*) 
!  END IF

      END SUBROUTINE get_lib_hm_mass_of_fuel

      end module termsort_lib
