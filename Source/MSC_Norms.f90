      subroutine  MSC_Get_Norm_1(N, x, x_norm)
!=====================================================================*
!              Infinite Norm of the Vector                            *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x(N)
! Output:
      real x_norm
! Local:
      integer k

      x_norm = 0.
      do k = 1, N
         x_norm = x_norm + abs( x(k) )
      end do

      return
      end

      subroutine  MSC_Get_Norm_Inf(N, x, x_norm)
!=====================================================================*
!              Infinite Norm of the Vector                            *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x(N)
! Output:
      real x_norm
! Local:
      integer k

      x_norm = 0.
      do k = 1, N
         x_norm = amax1(abs(x(k)), x_norm)
      end do

      return
      end

      subroutine  MSC_Get_Norm_2(N, x, x_norm)
!=====================================================================*
!              2 Norm of the Vector                                   *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                 *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x(N)
! Output:
      real x_norm
! Local:
      integer k

      x_norm = 0.
      do k = 1, N
         x_norm = x_norm + x(k)*x(k)
      end do

      x_norm = sqrt(x_norm)

      return
      end


      subroutine MSC_Get_Abs_Error_Norm_Inf&
                  (N, x_new, x_old, x_error)
!=====================================================================*
!    Absolute Difference between 2 Vectors  ( Infinite Norm )         *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x_new(N), x_old(N)
! Output:
      real x_error
! Local:
      integer k

      x_error = 0.
      do k = 1, N
         x_error  =  amax1( abs(x_new(k)-x_old(k)), x_error )
      end do

      return
      end

      subroutine MSC_Get_Abs_Error_Norm_2(N, x_new, x_old, x_error)
!=====================================================================*
!    Absolute Difference between 2 Vectors  ( 1 Norm )                *                                   *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x_new(N), x_old(N)
! Output:
      real x_error
! Local:
      integer k
      double precision diff, x_error_double

      x_error_double = 0.
      do k = 1, N
         diff = x_new(k)-x_old(k)
         x_error_double  = x_error_double + diff*diff
      end do

!      WRITE(*,*) 'MAX_NEW', MAXVAL( x_new, N )
!      WRITE(*,*) 'MIN_NEW', MINVAL( x_new, N )
!      WRITE(*,*) 'AV_NEW', SUM( x_new, N )/N
!      WRITE(*,*) 'MAX_OLD', MAXVAL( x_old, N )
!      WRITE(*,*) 'MIN_OLD', MINVAL( x_old, N )
!      WRITE(*,*) 'AV_OLD', SUM( x_old, N )/N
!      WRITE(*,*) 'MSC_Get_Abs_Error_Norm_2', x_error
!      PAUSE

      x_error = REAL( sqrt(x_error_double) )

      return
      end

      subroutine  MSC_DOT_PRODUCT(N, x1, x2, dot_product)
!=====================================================================*
!       Dot Product of 2 Vectors                                      *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x1(N), x2(N)
! Output:
      real dot_product
! Local:
      integer k

      dot_product = 0.
      do k = 1, N
         dot_product = dot_product + x1(k)*x2(k)
      end do

      return
      end


      subroutine  MSC_MAX_MIN_RATIO(N, x1, x2, &
                           Ratio_Max, Ratio_Min)
!=====================================================================*
! Maximum and minimum ratio of the vector components  x1(i)/x2(i)     *
! (c) Slava 15.VII.1999  slava@lstf3.tokai.jaeri.go.jp                *
!=====================================================================*
      implicit none
! Input:
      integer N
      real x1(N), x2(N)
! Output:
      real ratio_max, ratio_min
! Local:
      integer k
      real ratio
      real BIG_VALUE
      parameter(BIG_VALUE = 1.E+30)

      ratio_max = 0.
      ratio_min = BIG_VALUE
 
      do k = 1, N

         if( x2(k) .NE. 0 ) then 
           ratio = abs( x1(k) / x2(k) )
           ratio_min = amin1(ratio_min, ratio)
           ratio_max = amax1(ratio_max, ratio)
         end if

      end do

      return
      end
 