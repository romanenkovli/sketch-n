      subroutine amux(n, x, y, a, ja, ia) 
      real  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      real t
      integer i, k
!-----------------------------------------------------------------------
      do 100 i = 1,n
!
!     compute the inner product of row i with vector x
! 
         t = 0.0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
!
!     store result in y(i) 
!
         y(i) = t
 100  continue
!
      return
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine atmux(n,n2, x, y, a, ja, ia)
      real x(*), y(*), a(*) 
      integer n2, n, ia(*), ja(*)
!-----------------------------------------------------------------------
!         transp( A ) times a vector
! Modified 25.01.99 RECTANGULAR VERSION
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!----------------------------------------------------------------------- 
! multiplies the transpose of a matrix by a vector when the original
! matrix is stored in compressed sparse row storage. Can also be
! viewed as the product of a matrix by a vector when the original
! matrix is stored in the compressed sparse column format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!     local variables 
!
      integer i, k 
!-----------------------------------------------------------------------
!
!     zero out output vector
! 
      do 1 i=1,n2 
         y(i) = 0.0
 1    continue
!
! loop over the rows
!
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
!
      return
      end

      subroutine csrcsc2(n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer n, n2
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real  a(*),ao(*)
      integer next, k, i, j, ipos, job
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n = number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!       for any other normal usage, enter ipos=1.
! a= real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja = integer array of length nnz containing the column positions
!        of the corresponding elements in a.
! ia = integer of size n+1. ia(k) contains the position in a, ja of
!       the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao = real array of size nzz containing the "a" part of the transpose
! jao = integer array of size nnz containing the column indices.
! iao = integer array of size n+1 containing the "ia" index array of
!       the transpose. 
!
!----------------------------------------------------------------------- 
!----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
!---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
!--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
!-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
!--------------- end of csrcsc2 ---------------------------------------- 
!-----------------------------------------------------------------------
      end
