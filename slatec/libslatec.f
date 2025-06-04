! Functions and subroutins for SKETCH based on SLATEC library
!
! SUBROUTINE EVLRG(N, A, LDA, EV)
! SUBROUTINE POLRG(N, A, LDA, NCOEF, COEF, B, LDB) 
! SUBROUTINE MULTIPLY(N1, M1, A, N2, M2, B, C) 
! SUBROUTINE MLINE(N, M, A, B, c1, c2) 
! SUBROUTINE JUST(Z,Y)
! COMPLEX FUNCTION JUST2(Z)
!
!*******************************************************************

      SUBROUTINE EVLRG(N, A, LDA, EV)
    ! subroutin for eigenvalues calculation
    ! N - order of matix
    ! A - matrix
    ! LDA (= N) - leading dimension of matrix
    ! EV - real part of eigen values
      INTEGER N
      REAL WI(N) ! imaginary part of eigen values 
      REAL Z(N,N) ! real and imaginary part of eigen vectors
      REAL FV1(N) !tenporary storage
      INTEGER NM ! row dimension
      INTEGER MATZ ! = 0 if only eigenvalues desired
      INTEGER IV1(N) ! tenporary storage
      INTEGER IERR ! = 0 - for normal return
      MATZ = 1.0
      IERR = 0
      NM = LDA
      CALL RG (NM, N, A, EV, WI, MATZ, Z, IV1, FV1, IERR)
      END 
      

      SUBROUTINE S_POLRG(N, A, LDA, NCOEF, COEF, B, LDB) 
    ! The routine POLRG computes the matrix polynomial
    ! using Horner's scheme
      INTEGER N
      REAL A(LDA, N), B(LDA, N) ! input matrices
      REAL COEF(NCOEF) ! coefficients
      REAL B2(LDA, N) ! temporary matrix
      CALL MLINE(LDA, N, A, B, COEF(NCOEF), COEF(NCOEF-1))
      DO k = 2, NCOEF-1
        CALL MULTIPLY(LDA, N, A, N, LDA, B, B2)
        CALL MLINE(LDA, N, B2, B, 1.0, COEF(NCOEF-k))
      ENDDO
      END

      SUBROUTINE MULTIPLY(N1, M1, A, N2, M2, B, C) 
    ! myltiply matrices C = A * B
    ! A - (N1, M1) matrix 
    ! B - (N2, M2) matrix
    ! C - (N1, M2) matrix
    ! Must be M1 = N2
      REAL A(N1, M1), B(N2, M2), C(N1, M2) ! matrices
      REAL VAL ! temporary value
      IF (M1 .eq. N2) THEN
          DO n = 1, N1
            DO m = 1, M2
              VAL = 0.0
              DO i = 1, N1
                  VAL = VAL + A(n,i) * B(i,m)
              END DO
              C(n,m) = VAL
            END DO
          END DO
      ELSE
        WRITE(*,*) 'CAN NOT MULTYPLY MATRICES WITH dim[A(1,:)] 
     & NONEQUAL dim[B(:,1)]'
      END IF
      END

      SUBROUTINE MLINE(N, M, A, B, c1, c2) 
    ! calculate B = c1*A + c2*I matrix
    ! A, B - (N, M) matrices
    ! a, b - real coefficients
      REAL A(N, M), B(N, M) ! matrices
      DO nn = 1, N
        DO mm = 1, M
            IF (nn .eq.mm) THEN
                B(nn,mm) = c1 * A(nn,mm) + c2
            ELSE
                B(nn,mm) = c1 * A(nn,mm)
            ENDIF
        ENDDO
      ENDDO
      END


      SUBROUTINE JUST(Z,Y)
      COMPLEX :: Z, Y
      Y = Z
      END

      COMPLEX FUNCTION JUST2(Z)
      COMPLEX :: Z
      JUST2 = Z
      RETURN
      END

      SUBROUTINE ZBREN(F, ERRABS, ERRREL, A, B, MAXFN)
    ! Finds a zero of a real function that changes sign in a given interval
    ! F  User-supplied FUNCTION to compute the value of the function of which a zero will be found
    ! ERRABS - absolute error
    ! ERRREL - relative error
    ! (A, B) - interval for x
    ! MAXFN - upper bound on the number of function evaluations required for convergence


      INTEGER    NOUT, MAXFN, IFLAG
      REAL       A, B, R, ERRABS, ERRREL, F
      EXTERNAL F

      R = 0.5*(A+B)
      R = F(R)

      call FZERO (F, A, B, R, ERREL, ERRABS, IFLAG)
      end SUBROUTINE ZBREN


      SUBROUTINE DZBREN(F, ERRABS, ERRREL, A, B, MAXFN)
    ! Finds a zero of a real function that changes sign in a given interval
    ! F  User-supplied FUNCTION to compute the value of the function of which a zero will be found
    ! ERRABS - absolute error
    ! ERRREL - relative error
    ! (A, B) - interval for x
    ! MAXFN - upper bound on the number of function evaluations required for convergence


      INTEGER  :: NOUT, MAXFN, IFLAG
      REAL*8   :: A, B, R, ERRABS, ERRREL, F
      EXTERNAL F


      R = 0.5*(A+B)
      R = F(R)

      call DFZERO (F, A, B, R, ERRREL, ERRABS, IFLAG)
      if(IFLAG.ne.1)then
      stop IFLAG
      endif
      end SUBROUTINE DZBREN
