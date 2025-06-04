      subroutine MSC_progonka(NN, A, b, x)
      implicit none
!=====================================================================*
! Solution of the triagonal system of the liner equations: A x = b    *
!            Diagonal of the Matrix A is destroyed                    *
! Sect. 1.2 (page 33) V. P. Il'in "Incomplete Factorization Methods"  *
!=====================================================================*
! Input:
      integer NN
      real A(3,NN), b(NN) ! matrix and the LHS
! Output:
      real x(NN) ! solution
! Local Variables:
      integer n
      real denom

      denom = 1. / A(2,1)
      A(2,1) = -A(3,1)*denom
      x(1) = b(1)*denom

      do n = 2, NN
         denom = 1./(a(2,n) + a(1,n)*a(2,n-1))
         a(2,n) = -a(3,n)*denom
         x(n) = (-a(1,n)*x(n-1) + b(n))*denom
      end do

      do n = NN-1, 1, -1
           x(n) = x(n+1)*a(2,n) + x(n)
      end do

      return
      end

      Subroutine MSC_Solve2x2(A, B)
!=====================================================================*
! Solving 2 Linear Equations A*x=B                                    *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      implicit none
      integer NG
      parameter (NG=2)
! Input: 
      real A(NG,NG), B(NG)
! Output: RHS - Solution
!Local Variables 
      real tmp_b, determ

      determ = a(2,2)*a(1,1) - a(1,2)*a(2,1)

      tmp_b = b(1)
      b(1) = (a(2,2)*b(1) - a(1,2)*b(2))/determ
      b(2) = (a(1,1)*b(2) - a(2,1)*tmp_b)/determ

      return
      end



      subroutine MSC_LU_Solve(A1, N, B1)
      implicit none 
!=====================================================================*
! Solution of the equation A1 * x = B1 Using LU Decomposition         * 
! Last Update              Slava (c) April 1998                       *
!=====================================================================*
! Input: A1 - Matrix, N - Dimension, B1 - Right Side
      integer N
      real A1(N,N),B1(N)
! Output: A1 - LU Decomposition of the MAtrix A1, B - Solution
! Local Variables: 
      integer k, i, j
      real sum 

! LU Decomposition
      do j = 1,N
! first part beta_ij expression (2.3.12)
         do i= 1, j
            sum = a1(i,j)
            do k = 1, i-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
            a1(i,j) = sum
         end do
! second part alfa_ij expression (2.3.13)
       do i = j+1,N
            sum = a1(i,j)
            do k = 1, j-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
         a1(i,j) = sum/a1(j,j)
       end do
      end do

! Solution LU x = b

! forward substitution
      do i = 1,N
         sum = b1(i)
         do j = 1, i-1
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i) = sum
       end do

! backsubstitution

      do i = N, 1,-1
         sum = b1(i)
         do j = i+1,N
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i)=sum/a1(i,i)
      end do

      return
      end
     


      SUBROUTINE MSC_Polcoe(x,y,n,cof)
!=====================================================================*
!    Coefficients of the interpolating polynomial                     *
!      "Numerical Recipes", Section 3.5 Pages 114-115                 *
!  Given arrays x(1:n) and y(1:n) containing a tabulated function     *
!   y(i) = f(xi), this routine returns an array of coefficients       *
!                   cof(1:n)                                          *
! Last Update              Slava (c) 5.06.1998                        *
!=====================================================================*
      INTEGER n,NMAX
      complex cof(n),x(n),y(n)
      PARAMETER (NMAX=15) ! Largest anticipated value of n.

      INTEGER i,j,k
      complex b,ff,phi,s(NMAX)
      do 11 i=1,n
         s(i)=0.
         cof(i)=0.
   11 enddo ! 11
      s(n)=-x(1)
      do 13 i=2,n 
! Coeffcients si of the master polynomial P(x) are found by recurrence. 
         do 12 j=n+1-i,n-1
            s(j)=s(j)-x(i)*s(j+1)
   12    enddo ! 12
         s(n)=s(n)-x(i)
   13 enddo ! 13
      do 16 j=1,n
         phi=n
         do 14 k=n-1,1,-1 !The quantity 
            phi=k*s(k+1)+x(j)*phi
   14    enddo !14
         ff=y(j)/phi
         b=1. 
         do 15 k=n,1,-1
            cof(k)=cof(k)+b*ff
            b=s(k)+x(j)*b
   15    enddo ! 15
   16 enddo ! 16
      return
      END


      subroutine MSC_LU_Decomp(A1, N)
      implicit none 
!=====================================================================*
!  LU Decomposition of the Matrix A                                   * 
! Last Update              Slava (c) April 1998                       *
!=====================================================================*
! Input: A1 - Matrix, N - Dimension
      integer N
      real A1(N,N)
! Output: A1 - Matrix Inverse
! Local Variables:
      integer k, i, j
      real sum 

! LU Decomposition
      do j = 1,N
! first part beta_ij expression (2.3.12)
         do i= 1, j
            sum = a1(i,j)
            do k = 1, i-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
            a1(i,j) = sum
         end do
! second part alfa_ij expression (2.3.13)
       do i = j+1,N
            sum = a1(i,j)
            do k = 1, j-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
         a1(i,j) = sum/a1(j,j)
       end do
      end do

      return
      end



      subroutine MSC_LUBKSB(A1,N,B1)
      implicit none
      integer N
      real A1(N,N),B1(N)
      integer i,j
      real sum

! forward substitution
      do i = 1,N
         sum = b1(i)
         do j = 1, i-1
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i) = sum
       end do

! backsubstitution

      do i = N, 1,-1
         sum = b1(i)
         do j = i+1,N
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i)=sum/a1(i,i)
      end do

      return
      end


      subroutine MSC_SSET(N, a, x)
!=====================================================================*
!  BLAS 1 Like Subroutine:  Set a Vector to Scalar                    *
!                       x(N) = a                                      *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
!input
      integer N
      real a
! output 
      real x(N)
! Local 
      integer i

      do i = 1, N
         x(i) = a
      end do

      return
      end


      subroutine MSC_SCOPY(N, x1, x2)
!=====================================================================*
!  BLAS 1 Like Subroutine:  Copy a vector x1 onto  x2                 *
!                       x2(:) = x1(:)                                 *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
!input
      integer N
      real x1(N)
! output 
      real x2(N)
! Local 
      integer i

      do i = 1, N
         x2(i) = x1(i)
      end do

      return
      end



      subroutine MSC_SSCALE(N, a, x1)
!=====================================================================*
!  BLAS 1 Like Subroutine: scale  a vector x1 by scalar a             *
!                       x(:) = a*X(:)                                 *
!                   Vyachreslav Zimin (c) May 18 1998                 *
!               e-mail:       vgzimin@mail.ru                         *
!=====================================================================*
      implicit none
!input
      integer N
      real x1(N), a
! output 
!      real x1(N)
! Local 
      integer i

      do i = 1, N
         x1(i) = a*x1(i)
      end do

      return
      end

      subroutine MSC_DReplace_Array(N, A1, A2)
!=====================================================================*
! Replacing Array A2 by Array A1
!=====================================================================*
      implicit none
      integer N
      double precision A1(N), A2(N)
! Input: A1(N), A2(N)
! Output: A2(N) = A1(N(
! Local Variables:
      integer k
  
      do k = 1, N 
         A2(k) = A1(k)
      end do

      return
      end 




      subroutine MSC_Search_Header_In_File&
                (io_unit, header, line, fmt, error)
      implicit none
! Input:
      integer io_unit !, ios
      character*(*) header, line
      character*(*) fmt
!Local Variables
      Logical Found, Error
      Found = .False.
      do while (.NOT. Found )
         read(io_unit, FMT, ERR=100, END=100 ) LINE
!        read(io_unit, FMT, iostat=ios) LINE
!        IF(ios /=0) GO TO 100
!       write(*,*) 'LINE =', LINE
!       pause
        IF( INDEX( TRIM(Line), TRIM(header) ) .GT. 0  ) THEN 
!        if(Line .EQ. header) then
          Found =.True.
          Error = .False.
        end if
      end do
      go to 200
  100 Error = .True.
  200 continue
      return
      end

      subroutine MSC_ERR_Add_MESSAGE(Message_Type,Message)
      implicit none
      include 'units.fh'
! Input:
      character*(*) Message, Message_type
! Local
      CHARACTER*10 format_output
      CHARACTER*3 format_char
      INTEGER n_line_output, n_beg, n_end, n_line, length_message,&
               n_left

! Preparing the output
      
      length_message = LEN(Message)
      n_line_output = length_message/length_output
      n_left = MOD(length_message, length_output)

      WRITE(format_char, '("A", I2)') LENGTH_OUTPUT
      format_output = '('//format_char//')'
!      WRITE(format_output, FMT='( "'(A", I2, ")'" )' ) LENGTH_OUTPUT

      OPEN(out_unit_err, FILE ='Output/Errors.msg', STATUS = &
        'UNKNOWN', ACCESS ='APPEND')

        WRITE(out_unit_err, '(A7,":")') Message_Type
        IF(Error_Message_on_Screan) then
            WRITE(*, '(A7,":")') Message_Type
        END IF

        n_end = 0
        DO n_line = 1, n_line_output
           n_beg = n_end + 1
           n_end = n_beg + length_output - 1
           WRITE(out_unit_err, FMT = format_output) &
             Message(n_beg:n_end)
           IF(Error_Message_on_Screan) then
              WRITE(*, FMT = format_output) &
                 Message(n_beg:n_end)
           END IF
        END DO
        IF(n_left .NE. 0) THEN
           WRITE(format_char, '("A", I2)') n_left+1
           format_output = '('//format_char//')'
           n_beg = n_end + 1
           n_end = n_beg + n_left - 1
           WRITE(out_unit_err, FMT = format_output) &
             Message(n_beg:n_end)
           IF(Error_Message_on_Screan) then
              WRITE(*, FMT = format_output) &
                 Message(n_beg:n_end)
           END IF
        END IF
      CLOSE(out_unit_err)

      if(Message_Type.EQ.'ERROR') then
        stop
      end if

      return
      end


      subroutine MSC_ERR_Init_File
      implicit none
      include 'units.fh'
! Input:

      OPEN(out_unit_err, file ='Output/Errors.msg', status = &
        'unknown')

       CALL OUTput_Write_Header(out_unit_err)

      CLOSE(out_unit_err)

     
      return
      end

      subroutine Iostat_Error_Check(ios, Message)
      implicit none

      integer ios
      character*(*) Message

      if(ios .NE. 0) then
        write(*,*) Message
        stop
      end if

      return
      end

      REAL FUNCTION CONVERT_PPM_TO_GKG( ppm )
        CONVERT_PPM_TO_GKG = ppm*5.7222/1000.
      RETURN
      END                   
 