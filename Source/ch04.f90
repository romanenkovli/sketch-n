
      SUBROUTINE DPCHEZ(N,X,F,D,SPLINE,WK,LWK,IERR)
!***BEGIN PROLOGUE  DPCHEZ
!***DATE WRITTEN   870821   (YYMMDD)
!***REVISION DATE  870908   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
!             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!***PURPOSE  Easy to use spline or cubic Hermite interpolation.
!***DESCRIPTION
!
!          DPCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C. Moler, S. Nash
!               Prentice Hall 1988
!
!     Sets derivatives for spline (two continuous derivatives) or
!     Hermite cubic (one continuous derivative) interpolation.
!     Spline interpolation is smoother, but may not "look" right if the
!     data contains both "steep" and "flat" sections.  Hermite cubics
!     can produce a "visually pleasing" and monotone interpolant to
!     monotone data. This is an easy to use driver for the routines
!     by F. N. Fritsch in reference (4) below. Various boundary
!     conditions are set to default values by DPCHEZ. Many other choices
!     are available in the subroutines PCHIC, DPCHIM and DPCHSP.
!
!     Use PCHEV to evaluate the resulting function and its derivative.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:   CALL  DPCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
!
!     INTEGER  N, IERR,  LWK
!     DOUBLE PRECISION  X(N), F(N), D(N), WK(*)
!     LOGICAL SPLINE
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(I) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!
!     SPLINE -- (input) logical variable to specify if the interpolant
!           is to be a spline with two continuous derivaties
!           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
!           continuous derivative (set SPLINE=.FALSE.).
!        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
!           default "not-a-knot" boundary condition, with a continuous
!           third derivative at X(2) and X(N-1). See reference (3).
!              If SPLINE=.FALSE. the interpolating Hermite cubic will be
!           monotone if the input data is monotone. Boundary conditions
!           computed from the derivative of a local quadratic unless thi
!           alters monotonicity.
!
!     WK -- (scratch) real work array, which must be declared by the cal
!           program to be at least 2*N if SPLINE is .TRUE. and not used
!           otherwise.
!
!     LWK -- (input) length of work array WK. (Error return if
!           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means tha
!                 IERR switches in the direction of monotonicity were de
!                 When SPLINE=.FALSE.,  DPCHEZ guarantees that if the inp
!                 data is monotone, the interpolant will be too. This wa
!                 is to alert you to the fact that the input data was no
!                 monotone.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
! ----------------------------------------------------------------------
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
!                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
!                 1980), 238-246.
!               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
!                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
!                 PREPRINT UCRL-87559 (APRIL 1982).
!               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
!               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
!                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
!                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
!                 AUGUST 1982.
!***ROUTINES CALLED  DPCHIM,DPCHSP
!***END PROLOGUE  DPCHEZ
      INTEGER  N, LWK, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), WK(LWK)
      LOGICAL SPLINE
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER IC(2), INCFD
      DOUBLE PRECISION  VC(2)
      DATA IC(1) /0/
      DATA IC(2) /0/
      DATA INCFD /1/
!
!
!***FIRST EXECUTABLE STATEMENT  DPCHEZ
!
      IF ( SPLINE ) THEN
        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
!
!  ERROR CONDITIONS ALREADY CHECKED IN DPCHSP OR DPCHIM

      RETURN
!------------- LAST LINE OF DPCHEZ FOLLOWS ------------------------------
      END
      SUBROUTINE DPCHIM(N,X,F,D,INCFD,IERR)
!***BEGIN PROLOGUE  DPCHIM
!***DATE WRITTEN   811103   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(DPCHIM-S DPCHIM-D),
!             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine a monotone piecewise
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See DPCHIC if user control
!            is desired over boundary or switch conditions.)
!***DESCRIPTION
!
!       **** Double Precision version of DPCHIM ****
!
!          DPCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubic
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See DPCHIC if user control of boundary con-
!     ditions is desired.)
!
!     If the data are only piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See DPCHIC if user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).  DPCHIM is designed for monotonic data, but it will
!           work for any F-array.  It will force extrema at points where
!           monotonicity switches direction.  If some other treatment of
!           switch points is desired, DPCHIC should be used instead.
!                                     -----
!     D -- (output) real*8 array of derivative values at the data
!           points.  If the data are monotonic, these values will
!           determine a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
!                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
!                 1980), 238-246.
!               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
!                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
!                 PREPRINT UCRL-87559 (APRIL 1982).
!***ROUTINES CALLED  DPCHST,XERROR
!***END PROLOGUE  DPCHIM
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-02-01   1. Introduced  DPCHST  to reduce possible over/under-
!                   flow problems.
!                2. Rearranged derivative formula for same reason.
!     82-06-02   1. Modified end conditions to be continuous functions
!                   of data when monotonicity switches in next interval.
!                2. Modified formulas so end conditions are less prone
!                   of over/underflow problems.
!     82-08-03   Minor cosmetic changes for release 1.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. To produce a single precision version, simply:
!        a. Change DPCHIM to PCHIM wherever it occurs,
!        b. Change DPCHST to PCHST wherever it occurs,
!        c. Change all references to the Fortran intrinsics to their
!           single precision equivalents,
!        d. Change the double precision declarations to real, and
!        e. Change the constants ZERO and THREE to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE,&
           H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      DOUBLE PRECISION  DPCHST
      DATA  ZERO /0.D0/, THREE/3.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
!
!  NORMAL CASE  (N .GE. 3).
!
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (DABS(D(1,1)) .GT. DABS(DMAX))  D(1,1) = DMAX
      ENDIF
!
!  LOOP THROUGH INTERIOR POINTS.
!
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
!
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) )  42, 41, 45
!
!        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( DPCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
!
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = DMAX1( DABS(DEL1), DABS(DEL2) )
         DMIN = DMIN1( DABS(DEL1), DABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!
   50 CONTINUE
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (DABS(D(1,N)) .GT. DABS(DMAX))  D(1,N) = DMAX
      ENDIF
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'&
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIM -- INCREMENT LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIM -- X-ARRAY NOT STRICTLY INCREASING'&
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHIM FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHST(ARG1,ARG2)
!***BEGIN PROLOGUE  DPCHST
!***REFER TO  DPCHCE,DPCHCI,DPCHCS,DPCHIM
!***ROUTINES CALLED  (NONE)
!***DESCRIPTION
!
!         DPCHST:  DPCHIP Sign-Testing Routine.
!
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
!     The object is to do this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!  Fortran intrinsics used:  SIGN.
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHST to PCHST wherever it occurs,
!        b. Change all references to the Fortran intrinsics to their
!           single presision equivalents,
!        c. Change the double precision declarations to real, and
!        d. Change the constants  ZERO  and  ONE  to single precision.
!***END PROLOGUE  DPCHST
!
!  DECLARE ARGUMENTS.
!
      DOUBLE PRECISION  ARG1, ARG2
!
!  DECLARE LOCAL VARIABLES.
!
      DOUBLE PRECISION  ONE, ZERO
      DATA  ZERO /0.D0/,  ONE/1.D0/
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = DSIGN(ONE,ARG1) * DSIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  DPCHST = ZERO
!
      RETURN
!------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
      END
      SUBROUTINE DPCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)
!***BEGIN PROLOGUE  DPCHSP
!***DATE WRITTEN   820503   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(PCHSP-S DPCHSP-D),
!             CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
!             SPLINE INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine the Hermite represen-
!            tation of the cubic spline interpolant to given data, with
!            specified boundary conditions.
!***DESCRIPTION
!
!       **** Double Precision version of PCHSP ****
!
!          DPCHSP:   Piecewise Cubic Hermite Spline
!
!     Computes the Hermite representation of the cubic spline inter-
!     polant to the data given in X and F satisfying the boundary
!     conditions specified by IC and VC.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
!     NOTE:  This is a modified version of C. de Boor'S cubic spline
!            routine CUBSPL.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor'S cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine the cubic spline
!           interpolant with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real*8 array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK.LT.2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
!              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978), PP. 53-59.
!***ROUTINES CALLED  DPCHDF,XERROR
!***END PROLOGUE  DPCHSP
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-04   Converted to SLATEC library version.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHSP to PCHSP wherever it occurs,
!        b. Change the double precision declarations to real, and
!        c. Change the constants ZERO, HALF, ... to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER IBEG, IEND, INDEX, J, NM1
      DOUBLE PRECISION  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4),&
       ZERO
      DOUBLE PRECISION  DPCHDF
!
      DATA  ZERO /0.D0/, HALF/.5D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHSP
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
!
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
      IF ( NWK .LT. 2*N )  GO TO 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
!        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
!           INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
!                 --------------------------------
         D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)
!                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
!
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
!        PICK UP LAST IEND POINTS.
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
!           INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
!                 --------------------------------
         D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)
!                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
!
! --------------------( BEGIN CODING FROM CUBSPL )--------------------
!
!  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
!           NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
!           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)&
                             + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
!        SLOPE PRESCRIBED AT LEFT END.
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
!        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1)&
                       + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
      IF (IEND .EQ. 1)  GO TO 30
!
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
!           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
!           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
!           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
!           KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
!           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)&
                       + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
!        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
! --------------------(  END  CODING FROM CUBSPL )--------------------
!
!  NORMAL RETURN.
!
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'&
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHSP -- INCREMENT LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHSP -- X-ARRAY NOT STRICTLY INCREASING'&
                , 0, IERR, 1)
      RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERROR ('DPCHSP -- IC OUT OF RANGE'&
                , 0, IERR, 1)
      RETURN
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
      IERR = -7
      CALL XERROR ('DPCHSP -- WORK ARRAY TOO SMALL'&
                , 0, IERR, 1)
      RETURN
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
!   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      CALL XERROR ('DPCHSP -- SINGULAR LINEAR SYSTEM'&
                , 0, IERR, 1)
      RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM DPCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERROR ('DPCHSP -- ERROR RETURN FROM DPCHDF'&
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHSP FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHDF(K,X,S,IERR)
!***BEGIN PROLOGUE  DPCHDF
!***REFER TO  DPCHCE,DPCHSP
!***ROUTINES CALLED  XERROR
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!          DPCHDF:   DPCHIP Finite Difference Formula
!
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.
!
!     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
!     derivative approximations.
!
! ----------------------------------------------------------------------
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K.LT.2 .
!        DPCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
! ----------------------------------------------------------------------
!
!  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
!              Verlag (New York, 1978), pp. 10-16.
!
!***END PROLOGUE  DPCHDF
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHDF to PCHDF wherever it occurs,
!        b. Change the double precision declarations to real, and
!        c. Change the constant Zero to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER K, IERR
      DOUBLE PRECISION  X(K), S(K)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, J
      DOUBLE PRECISION  VALUE, ZERO
      DATA  ZERO /0.D0/
!
!  CHECK FOR LEGAL VALUE OF K.
!
!***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
!
!  EVALUATE DERIVATIVE AT X(K).
!
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
!
!  NORMAL RETURN.
!
      IERR = 0
      DPCHDF = VALUE
      RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!     K.LT.3 RETURN.
      IERR = -1
      CALL XERROR ('DPCHDF -- K LESS THAN THREE'&
                , 0, IERR, 1)
      DPCHDF = ZERO
      RETURN
!------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
      END
      SUBROUTINE DPCHEV(N,X,F,D,NVAL,XVAL,FVAL,DVAL,IERR)
!***BEGIN PROLOGUE  DPCHEV
!***DATE WRITTEN   870828   (YYMMDD)
!***REVISION DATE  870828   (YYMMDD)
!***CATEGORY NO.  E3,H1
!***KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE
!             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
!***AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             ROOM A161, TECHNOLOGY BUILDING
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!***PURPOSE  Evaluates the function and first derivative of a piecewise
!            cubic Hermite or spline function at an array of points
!            XVAL.  It is easy to use.
!***DESCRIPTION
!
!          DPCHEV:  Piecewise Cubic Hermite or Spline Derivative
!                   Evaluator. Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C. Moler, S. Nash
!                 Prentice Hall 1988
!
!     Evaluates the function and first derivative of the cubic Hermite
!     or spline function defined by  N, X, F, D, at the array of points
!     XVAL.
!
!
!     This is an easy to use driver for the routines by F.N. Fritsch
!     described in reference (2) below. Those also have other 
!     capabilities.
!
! ----------------------------------------------------------------------
!
!  Calling sequence: CALL  DPCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
!
!     INTEGER  N, NVAL, IERR
!     DOUBLE PRECISION  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) double precision array of independent variable 
!           values.  The elements of X must be strictly increasing:
!             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
!
!     F -- (input) double precision array of function values.  F(I) is
!           the value corresponding to X(I).
!
!     D -- (input) double precision array of derivative values.  
!          D(I) is the value corresponding to X(I).
!
!  NVAL -- (input) number of points at which the functions are to be
!           evaluated. ( Error return if NVAL.LT.1 )
!
!  XVAL -- (input) double precision array of points at which the 
!          functions are to be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XVAL are increasing relative to X;
!              that is,   XVAL(J) .GE. X(I)
!              implies    XVAL(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XVAL are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  FVAL -- (output) double precision array of values of the cubic 
!          Hermite function defined by  N, X, F, D  at the points  XVAL.
!
!  DVAL -- (output) double precision array of values of the 
!          first derivative of the same function at the points  XVAL.
!
!  IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NVAL.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine DCHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
! ----------------------------------------------------------------------
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
!                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
!                 1980), 238-246.
!               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
!                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
!                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
!                 AUGUST 1982.
!***ROUTINES CALLED  DPCHFD
!***END PROLOGUE  DPCHEV
      INTEGER  N, NVAL, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL),&
     DVAL(NVAL)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER INCFD
      LOGICAL SKIP
      DATA SKIP /.TRUE./
      DATA INCFD /1/

!
!
!***FIRST EXECUTABLE STATEMENT  DPCHEV
!
      CALL DPCHFD(N,X,F,D,INCFD,SKIP,NVAL,XVAL,FVAL,DVAL,IERR)
!
!
 5000 CONTINUE
      RETURN
!
!------------- LAST LINE OF DPCHEV FOLLOWS ------------------------------
      END
      SUBROUTINE DPCHFD(N,X,F,D,INCFD,SKIP,NE,XE,FE,DE,IERR)
!***BEGIN PROLOGUE  DPCHFD
!***DATE WRITTEN   811020   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E3,H1
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(PCHFD-S DPCHFD-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             HERMITE INTERPOLATION,PIECEWISE CUBIC EVALUATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate a piecewise cubic hermite function and its first
!            derivative at an array of points.  May be used by itself
!            for Hermite interpolation, or as an evaluator for DPCHIM
!            or DPCHIC. If only function values are required, use
!            DPCHFE instead.
!***DESCRIPTION
!
!       **** Double Precision version of PCHFD ****
!
!          DPCHFD:  Piecewise Cubic Hermite Function and Derivative
!                  evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!     gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use DPCHFE, instead.
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
!                          DE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine DCHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFDV,XERROR
!***END PROLOGUE  DPCHFD
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. To produce a single precision version, simply:
!        a. Change DPCHFD to PCHFD, and DCHFDV to CHFDV, wherever they
!           occur,
!        b. Change the double precision declaration to real,
!
!     2. Most of the coding between the call to DCHFDV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. DCHFDV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of DPCHFD that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points .LT.
!        X(N-1), followed by points .GT.X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD, NE, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),&
      DE(NE)
      LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHFD
      IF (SKIP)  GO TO 5
!
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
         IF (JFIRST .GT. NE)  GO TO 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
!
   40    CONTINUE
         NJ = J - JFIRST
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
         IF (NJ .EQ. 0)  GO TO 50
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
        CALL DCHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR)&
                   ,NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
!       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
!
         IF (NEXT(2) .EQ. 0)  GO TO 42
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
            IF (IR .LT. N)  GO TO 41
!           IF (IR .EQ. N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
!           ENDIF
!        ENDIF
   42    CONTINUE
!
         IF (NEXT(1) .EQ. 0)  GO TO 49
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
            IF (IR .GT. 2)  GO TO 43
!           IF (IR .EQ. 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN DCHFDV.
               GO TO 5005
!
   45          CONTINUE
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!              AT THIS POINT, EITHER  XE(J) .LT. X(1)
!                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
               IR = MAX0(1, I-1)
!           ENDIF
!        ENDIF
   49    CONTINUE
!
         JFIRST = J
!
!     END OF IR-LOOP.
!
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'&
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHFD -- INCREMENT LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHFD -- X-ARRAY NOT STRICTLY INCREASING'&
                , 0, IERR, 1)
      RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -4
      CALL XERROR ('DPCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM DCHFDV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      CALL XERROR ('DPCHFD -- ERROR RETURN FROM DCHFDV -- FATAL'&
                , 0, IERR, 2)
      RETURN
!------------- LAST LINE OF DPCHFD FOLLOWS -----------------------------
      END
      SUBROUTINE DCHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
!***BEGIN PROLOGUE  DCHFDV
!***DATE WRITTEN   811019   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E3,H1
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(CHFDV-S DCHFDV-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
!            first derivative at an array of points.  While designed for
!            use by DPCHFD, it may be useful directly as an evaluator
!            for a piecewise cubic Hermite function in applications,
!            such as graphing, where the interval is known in advance.
!            If only function values are required, use DCHFEV instead.
!***DESCRIPTION
!
!       **** Double Precision version of CHFDV ****
!
!        DCHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use DCHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE),
!                          DE(NE)
!
!        CALL  DCHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1.EQ.X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1.EQ.X2 .
!                (Output arrays have not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  DCHFDV
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!     87-07-07   Corrected XERROR calls for d.p. names(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DCHFDV to CHFDV wherever it occurs,
!        b. Change the double precision declaration to real,
!        c. Change the constant Zero to single precision, and
!        d. Change the names of the Fortran functions:  AMAX1, AMIN1.
!
!  DECLARE ARGUMENTS.
!
      INTEGER NE, NEXT(2), IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I
      DOUBLE PRECISION  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X,&
       XMI, XMA, ZERO
      DATA  ZERO /0.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DCHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
!
!  INITIALIZE.
!
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = DMIN1(ZERO, H)
      XMA = DMAX1(ZERO, H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
!                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
!
!  EVALUATION LOOP.
!
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
!          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
!
!  NORMAL RETURN.
!
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -1
      CALL XERROR ('DCHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'&
                , 51, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERROR ('DCHFDV -- INTERVAL ENDPOINTS EQUAL'&
                , 34, IERR, 1)
      RETURN
!------------- LAST LINE OF DCHFDV FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHQA(N,X,F,D,A,B,IERR)
!***BEGIN PROLOGUE  DPCHQA
!***DATE WRITTEN   870829   (YYMMDD)
!***REVISION DATE  870829   (YYMMDD)
!***CATEGORY NO.  E3,H2A2
!***KEYWORDS  EASY TO USE CUBIC HERMITE OR SPLINE INTEGRATION
!             NUMERICAL INTEGRATION, QUADRATURE
!***AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             ROOM A161, TECHNOLOGY BUILDING
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!***PURPOSE  Evaluates the definite integral of a piecewise cubic Hermit
!            or spline function over an arbitrary interval, easy to use.
!***DESCRIPTION
!
!          DPCHQA:  Piecewise Cubic Hermite or Spline Integrator,
!                  Arbitrary limits, Easy to Use.
!
!          From the book "Numerical Methods and Software"
!                  by  D. Kahaner, C. Moler, S. Nash
!                          Prentice Hall 1988
!
!     Evaluates the definite integral of the cubic Hermite or spline
!     function defined by  N, X, F, D  over the interval [A, B].  This
!     is an easy to use driver for the routine DPCHIA by F.N. Fritsch
!     described in reference (2) below. That routine also has other
!     capabilities.
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!           VALUE = DPCHQA (N, X, F, D, A, B, IERR)
!
!     INTEGER  N, IERR
!     DOUBLE PRECISION  X(N), F(N), D(N), A, B
!
!   Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) double precision array of independent variable
!           values.  The elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) double precision array of function values.
!           F(I) is the value corresponding to X(I).
!
!     D -- (input) double precision array of derivative values.  D(I) is
!           the value corresponding to X(I).
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (Value has not been computed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
!                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
!                 1980), 238-246.
!               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
!                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
!                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
!                 AUGUST 1982.
!***ROUTINES CALLED  DPCHIA
!***END PROLOGUE  DPCHQA
      INTEGER  N, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), A, B
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER  INCFD
      DOUBLE PRECISION  DPCHIA
      LOGICAL SKIP
!
!  INITIALIZE.
!
      DATA  INCFD /1/
      DATA  SKIP /.TRUE./
!
!
!***FIRST EXECUTABLE STATEMENT  DPCHQA

      DPCHQA  =  DPCHIA( N, X, F, D, INCFD, SKIP, A, B, IERR )
!
! ERROR MESSAGES ARE FROM LOWER LEVEL ROUTINES
      RETURN
!
!------------- LAST LINE OF DPCHQA FOLLOWS ------------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHIA(N,X,F,D,INCFD,SKIP,A,B,IERR)
!***BEGIN PROLOGUE  DPCHIA
!***DATE WRITTEN   820730   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E3,H2A2
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(PCHIA-S DPCHIA-D),
!             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
!             QUADRATURE
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate the definite integral of a piecewise cubic
!            Hermite function over an arbitrary interval.
!***DESCRIPTION
!
!       **** Double Precision version of PCHIA ****
!
!          DPCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [A, B].
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
!        DOUBLE PRECISION  VALUE, DPCHIA
!        LOGICAL  SKIP
!
!        VALUE = DPCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!   Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR.GE.0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (Value has not been computed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFIV,DPCHID,XERROR
!***END PROLOGUE  DPCHIA
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-04   Converted to SLATEC library version.
!     87-07-07   Corrected conversion to double precision.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHIA to PCHIA wherever it occurs,
!        b. Change DPCHID to PCHID wherever it occurs,
!        c. Change DPCHIV to PCHIV wherever it occurs,
!        d. Change the double precision declarations to real,  and
!        e. Change the constant  ZERO  to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
      LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, IA, IB, IERD, IERV, IL, IR
      DOUBLE PRECISION  VALUE, XA, XB, ZERO
      DOUBLE PRECISION  DCHFIV, DPCHID
!
!  INITIALIZE.
!
      DATA  ZERO /0.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHIA
      IF (SKIP)  GO TO 5
!
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
    5 CONTINUE
      SKIP = .TRUE.
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
!
!  COMPUTE INTEGRAL VALUE.
!
      IF (A .EQ. B)  THEN
         VALUE = ZERO
      ELSE
         XA = DMIN1 (A, B)
         XB = DMAX1 (A, B)
         IF (XB .LE. X(2))  THEN
!           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!                   --------------------------------------------
            VALUE = DCHFIV (X(1),X(2), F(1,1),F(1,2),&
                                     D(1,1),D(1,2), A, B, IERV)
!                   --------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE IF (XA .GE. X(N-1))  THEN
!           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!                   -----------------------------------------------
            VALUE = DCHFIV(X(N-1),X(N), F(1,N-1),F(1,N),&
                                      D(1,N-1),D(1,N), A, B, IERV)
!                   -----------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE
!           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
!      ......LOCATE IA AND IB SUCH THAT
!               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
!             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
!             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
!
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
!             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
!             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
!
!     ......COMPUTE THE INTEGRAL.
            IERV = 0
            IF (IB .LT. IA)  THEN
!              THIS MEANS IB = IA-1 AND
!                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
!                      ------------------------------------------------
               VALUE = DCHFIV (X(IB),X(IA), F(1,IB),F(1,IA),&
                                          D(1,IB),D(1,IA), A, B, IERV)
!                      ------------------------------------------------
               IF (IERV .LT. 0)  GO TO 5004
            ELSE
!
!              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
               IF (IB .EQ. IA)  THEN
                  VALUE = ZERO
               ELSE
!                         ---------------------------------------------
                  VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
!                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5005
               ENDIF
!
!              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX0 (1, IA-1)
                  IR = IL + 1
!                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),&
                                     D(1,IL),D(1,IR), XA, X(IA), IERV)
!                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
!
!              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN0 (IB+1, N)
                  IL = IR - 1
!                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),&
                                     D(1,IL),D(1,IR), X(IB), XB, IERV)
!                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
!
!              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
!
!  NORMAL RETURN.
!
      DPCHIA = VALUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIA -- NUMBER OF DATA POINTS LESS THAN TWO'&
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIA -- INCREMENT LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIA -- X-ARRAY NOT STRICTLY INCREASING'&
                , 0, IERR, 1)
      RETURN
!
 5004 CONTINUE
!     TROUBLE IN DCHFIV.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERROR ('DPCHIA -- TROUBLE IN DCHFIV'&
                , 0, IERR, 1)
      RETURN
!
 5005 CONTINUE
!     TROUBLE IN DPCHID.  (SHOULD NEVER OCCUR.)
      IERR = -5
      CALL XERROR ('DPCHIA -- TROUBLE IN DPCHID'&
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHIA FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHID(N,X,F,D,INCFD,SKIP,IA,IB,IERR)
!***BEGIN PROLOGUE  DPCHID
!***DATE WRITTEN   820723   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B,H2A2
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=DOUBLE PRECISION(PCHID-S DPCHID-D),
!             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
!             QUADRATURE
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate the definite integral of a piecewise cubic
!            Hermite function over an interval whose endpoints are data
!            points.
!***DESCRIPTION
!
!       **** Double Precision version of PCHID ****
!
!          DPCHID:  Piecewise Cubic Hermite Integrator, Data limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IA, IB, IERR
!        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
!
!   Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
!
!     IA,IB -- (input) indices in X-array for the limits of integration.
!           both must be in the range [1,N].  (Error return if not.)
!           No restrictions on their relative values.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IA or IB is out of range.
!                (Value has not been computed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  DPCHID
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-04   Converted to SLATEC library version.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHID to PCHID wherever it occurs,
!        b. Change the double precision declarations to real,  and
!        c. Change the constants ZERO, HALF, SIX to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD, IA, IB, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, IUP, LOW
      DOUBLE PRECISION  H, HALF, SIX, SUM, VALUE, ZERO
!
!  INITIALIZE.
!
      DATA  ZERO /0.D0/,  HALF/.5D0/, SIX/6.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHID
      IF (SKIP)  GO TO 5
!
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
    5 CONTINUE
      SKIP = .TRUE.
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
!
!  COMPUTE INTEGRAL VALUE.
!
      IF (IA .EQ. IB)  THEN
         VALUE = ZERO
      ELSE
         LOW = MIN0(IA, IB)
         IUP = MAX0(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +&
                           (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
!
!  NORMAL RETURN.
!
      DPCHID = VALUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHID -- NUMBER OF DATA POINTS LESS THAN TWO'&
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHID -- INCREMENT LESS THAN ONE'&
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHID -- X-ARRAY NOT STRICTLY INCREASING'&
                , 0, IERR, 1)
      RETURN
!
 5004 CONTINUE
!     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERROR ('DPCHID -- IA OR IB OUT OF RANGE'&
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHID FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DCHFIV(X1,X2,F1,F2,D1,D2,A,B,IERR)
!***BEGIN PROLOGUE  DCHFIV
!***REFER TO  DPCHIA
!***ROUTINES CALLED  XERROR
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!          DCHFIV:  Cubic Hermite Function Integral Evaluator.
!
!     Called by  DPCHIA  to evaluate the integral of a single cubic (in
!     Hermite form) over an arbitrary interval (A,B).
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  IERR
!        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
!        DOUBLE PRECISION  VALUE, DCHFIV
!
!        VALUE = DCHFIV (X1, X2, F1, F2, D1, D2, A, B, IERR)
!
!   Parameters:
!
!     VALUE -- (output) VALUE of the requested integral.
!
!     X1,X2 -- (input) endpoints if interval of definition of cubic.
!           (Must be distinct.  Error return if not.)
!
!     F1,F2 -- (input) function values at the ends of the interval.
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!
!     A,B -- (input) endpoints of interval of integration.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0 (no errors).
!           "Recoverable errors":
!              IERR = -1  if X1.EQ.X2 .
!                (VALUE has not been set in this case.)
!
!***END PROLOGUE  DCHFIV
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!     87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DCHFIV to CHFIV wherever it occurs,
!        b. Change the double precision declarations to real, and
!        c. Change the constants HALF, TWO, ... to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
!
!  DECLARE LOCAL VARIABLES.
!
      DOUBLE PRECISION  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2,&
           PHIB1, PHIB2, PSIA1, PSIA2, PSIB1, PSIB2, TA1, TA2, TB1,&
           TB2, THREE, TWO, UA1, UA2, UB1, UB2
!
!  INITIALIZE.
!
      DATA  HALF/.5D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/, SIX/6.D0/
!
!  VALIDITY CHECK INPUT.
!
!***FIRST EXECUTABLE STATEMENT  DCHFIV
      IF (X1 .EQ. X2)  GO TO 5001
      IERR = 0
!
!  COMPUTE INTEGRAL.
!
      H = X2 - X1
      TA1 = (A - X1) / H
      TA2 = (X2 - A) / H
      TB1 = (B - X1) / H
      TB2 = (X2 - B) / H
!
      UA1 = TA1**3
      PHIA1 = UA1 * (TWO - TA1)
      PSIA1 = UA1 * (THREE*TA1 - FOUR)
      UA2 = TA2**3
      PHIA2 =  UA2 * (TWO - TA2)
      PSIA2 = -UA2 * (THREE*TA2 - FOUR)
!
      UB1 = TB1**3
      PHIB1 = UB1 * (TWO - TB1)
      PSIB1 = UB1 * (THREE*TB1 - FOUR)
      UB2 = TB2**3
      PHIB2 =  UB2 * (TWO - TB2)
      PSIB2 = -UB2 * (THREE*TB2 - FOUR)
!
      FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
      DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
!
!  RETURN VALUE.
!
      DCHFIV = (HALF*H) * (FTERM + DTERM)
      RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('DCHFIV -- X1 EQUAL TO X2'&
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DCHFIV FOLLOWS -----------------------------
      END

      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Symbolic dump (should be locally written).
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  23 May 1979
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
!***BEGIN PROLOGUE  J4SAVE
!***REFER TO  XERROR
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                 = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                 = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                 = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                 = 6 Refers to the 2nd unit for error messages
!                 = 7 Refers to the 3rd unit for error messages
!                 = 8 Refers to the 4th unit for error messages
!                 = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!    Adapted from Bell Laboratories PORT Library Error Handler
!     Latest revision ---  23 MAY 1979
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      FUNCTION NUMXER(NERR)
!***BEGIN PROLOGUE  NUMXER
!***REFER TO  XERROR
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        NUMXER returns the most recent error number,
!        in both NUMXER and the parameter NERR.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 JUNE 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  NUMXER
!***FIRST EXECUTABLE STATEMENT  NUMXER
      NERR = J4SAVE(1,0,.FALSE.)
      NUMXER = NERR
      RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERABT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Aborts program execution and prints error message.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        ***Note*** machine dependent routine
!        XERABT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG and NMESSG are as in XERROR, except that NMESSG may
!        be zero, in which case no message is being supplied.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCLR
!***BEGIN PROLOGUE  XERCLR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Resets current error number to zero.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        This routine simply resets the current error number to zero.
!        This may be necessary to do in order to determine that
!        a certain error has occurred again since the last time
!        NUMXER was referenced.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 June 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XERCLR
!***FIRST EXECUTABLE STATEMENT  XERCLR
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
!***BEGIN PROLOGUE  XERCTL
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Allows user control over handling of individual errors.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCTL.
!        If the user has provided his own version of XERCTL, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        MESSG1 - the first word (only) of the error message.
!        NMESSG - same as in the call to XERROR or XERRWV.
!        NERR   - same as in the call to XERROR or XERRWV.
!        LEVEL  - same as in the call to XERROR or XERRWV.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
!***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERDMP
!***BEGIN PROLOGUE  XERDMP
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Prints the error tables and then clears them.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XERDMP prints the error tables, then clears them.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 June 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  XERSAV
!***END PROLOGUE  XERDMP
!***FIRST EXECUTABLE STATEMENT  XERDMP
      CALL XERSAV(' ',0,0,0,KOUNT)
      RETURN
      END
      SUBROUTINE XERMAX(MAX)
!***BEGIN PROLOGUE  XERMAX
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Sets maximum number of times any error message is to be
!            printed.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XERMAX sets the maximum number of times any message
!        is to be printed.  That is, non-fatal messages are
!        not to be printed after they have occured MAX times.
!        Such non-fatal messages may be printed less than
!        MAX times even if they occur MAX times, if error
!        suppression mode (KONTRL=0) is ever in effect.
!
!     Description of Parameter
!      --Input--
!        MAX - the maximum number of times any one message
!              is to be printed.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 June 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XERMAX
!***FIRST EXECUTABLE STATEMENT  XERMAX
      JUNK = J4SAVE(4,MAX,.TRUE.)
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERPRT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  870916   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Prints error messages.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        Print the Hollerith message in MESSG, of length NMESSG,
!        on each file indicated by XGETUA.
!     Latest revision ---  16 SEPT 1987
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
!***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
!     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
!***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
!         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            IF(IUNIT.EQ.0)THEN
              WRITE (*,'(1X,A)') MESSG(ICHAR:LAST)
            ELSE
              WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
!***BEGIN PROLOGUE  XERROR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes an error (diagnostic) message.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XERROR processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed, containing
!                no more than 72 characters.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!
!     Examples
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
!                    43,2,1)
!        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
!    1ULLY COLLAPSED.',65,3,0)
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  XERRWV
!***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
!***BEGIN PROLOGUE  XERRWV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  870916   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes error message allowing 2 integer and two real
!            values to be included in the message.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XERRWV processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!        In addition, up to two integer values and two real
!        values may be printed along with the message.
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!        NI    - number of integer values to be printed. (0 to 2)
!        I1    - first integer value.
!        I2    - second integer value.
!        NR    - number of real values to be printed. (0 to 2)
!        R1    - first real value.
!        R2    - second real value.
!
!     Examples
!        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
!    1   1,NUM,0,0,0.,0.)
!        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
!    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
!
!     Latest revision ---  16 SEPT 1987
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
!                    XGETUA
!***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
!     GET FLAGS
!***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
!     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.&
         (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',&
       29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
!     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
!     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
!     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
!     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))&
     .OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))&
     .OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))&
     .OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
!           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT&
     ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT&
           ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
!        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
!            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) I1
                 IF (I.EQ.2) WRITE (*,FORM) I2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) I1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) I2
               ENDIF
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',&
              I2,'.',I2,')')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) R1
                 IF (I.EQ.2) WRITE (*,FORM) R2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) R1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) R2
               ENDIF
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
!              ERROR NUMBER
               IF(IUINT.EQ.0)THEN
                 WRITE(*,30) LERR
               ELSE
                 WRITE (IUNIT,30) LERR
               ENDIF
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
!        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))&
     IFATAL = 1
!     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
!        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT&
        ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT&
        ('JOB ABORT DUE TO FATAL ERROR.',29)
!        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
!     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
!***BEGIN PROLOGUE  XERSAV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Records that an error occurred.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        Record that this error occurred.
!
!     Description of Parameters
!     --Input--
!       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
!       except that when NMESSG=0 the tables will be
!       dumped and cleared, and when NMESSG is less than zero the
!       tables will be dumped and not cleared.
!     --Output--
!       ICOUNT will be the number of times this message has
!       been seen, or zero if the table has overflowed and
!       does not contain this message specifically.
!       When NMESSG=0, ICOUNT will not be altered.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 Mar 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
!***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
!     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
!     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),&
          KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)&
          /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
!***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
!     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
!        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/&
           51H MESSAGE START             NERR     LEVEL     COUNT)
!           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
!           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
!        CLEAR THE ERROR TABLES
         DO  I=1,10
          KOUNT(I) = 0
         end do
         KOUNTX = 0
         RETURN
   80 CONTINUE
!     PROCESS A MESSAGE...
!     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
!     THREE POSSIBLE CASES...
!     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
!     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
!     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETF(KONTRL)
!***BEGIN PROLOGUE  XGETF
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Returns current value of error control flag.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!   Abstract
!        XGETF returns the current value of the error control flag
!        in KONTRL.  See subroutine XSETF for flag value meanings.
!        (KONTRL is an output parameter only.)
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 June 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XGETF
!***FIRST EXECUTABLE STATEMENT  XGETF
      KONTRL = J4SAVE(2,0,.FALSE.)
      RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
!***BEGIN PROLOGUE  XGETUA
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Returns unit number(s) to which error messages are being
!            sent.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE XGETUN(IUNIT)
!***BEGIN PROLOGUE  XGETUN
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Returns the (first) output file to which messages are being
!            sent.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XGETUN gets the (first) output file to which error messages
!        are being sent.  To find out if more than one file is being
!        used, one must use the XGETUA routine.
!
!     Description of Parameter
!      --Output--
!        IUNIT - the logical unit number of the  (first) unit to
!                which error messages are being sent.
!                A value of zero means that the default file, as
!                defined by the I1MACH routine, is being used.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision --- 23 May 1979
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XGETUN
!***FIRST EXECUTABLE STATEMENT  XGETUN
      IUNIT = J4SAVE(3,0,.FALSE.)
      RETURN
      END
      SUBROUTINE XSETF(KONTRL)
!***BEGIN PROLOGUE  XSETF
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3A
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Sets the error control flag.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XSETF sets the error control flag value to KONTRL.
!        (KONTRL is an input parameter only.)
!        The following table shows how each message is treated,
!        depending on the values of KONTRL and LEVEL.  (See XERROR
!        for description of LEVEL.)
!
!        If KONTRL is zero or negative, no information other than the
!        message itself (including numeric values, if any) will be
!        printed.  If KONTRL is positive, introductory messages,
!        trace-backs, etc., will be printed in addition to the message.
!
!              IABS(KONTRL)
!        LEVEL        0              1              2
!        value
!          2        fatal          fatal          fatal
!
!          1     not printed      printed         fatal
!
!          0     not printed      printed        printed
!
!         -1     not printed      printed        printed
!                                  only           only
!                                  once           once
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE,XERRWV
!***END PROLOGUE  XSETF
!***FIRST EXECUTABLE STATEMENT  XSETF
      IF ((KONTRL.GE.(-2)).AND.(KONTRL.LE.2)) GO TO 10
         CALL XERRWV('XSETF  -- INVALID VALUE OF KONTRL (I1).',33,1,2,&
       1,KONTRL,0,0,0.,0.)
         RETURN
   10 JUNK = J4SAVE(2,KONTRL,.TRUE.)
      RETURN
      END
      SUBROUTINE XSETUA(IUNITA,N)
!***BEGIN PROLOGUE  XSETUA
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3B
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Sets up to 5 unit numbers to which messages are to be sent.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XSETUA may be called to declare a list of up to five
!        logical units, each of which is to receive a copy of
!        each error message processed by this package.
!        The purpose of XSETUA is to allow simultaneous printing
!        of each error message on, say, a main output file,
!        an interactive terminal, and other files such as graphics
!        communication files.
!
!     Description of Parameters
!      --Input--
!        IUNIT - an array of up to five unit numbers.
!                Normally these numbers should all be different
!                (but duplicates are not prohibited.)
!        N     - the number of unit numbers provided in IUNIT
!                must have 1 .LE. N .LE. 5.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE,XERRWV
!***END PROLOGUE  XSETUA
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XSETUA
      IF ((N.GE.1).AND.(N.LE.5)) GO TO 10
         CALL XERRWV('XSETUA -- INVALID VALUE OF N (I1).',34,1,2,&
       1,N,0,0,0.,0.)
         RETURN
   10 CONTINUE
      DO 20 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         JUNK = J4SAVE(INDEX,IUNITA(I),.TRUE.)
   20 CONTINUE
      JUNK = J4SAVE(5,N,.TRUE.)
      RETURN
      END
      SUBROUTINE XSETUN(IUNIT)
!***BEGIN PROLOGUE  XSETUN
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3B
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Sets output file to which error messages are to be sent.
!***DESCRIPTION
!    From the book "Numericl Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!           Prentice Hall 1988
!     Abstract
!        XSETUN sets the output file to which error messages are to
!        be sent.  Only one file will be used.  See XSETUA for
!        how to declare more than one file.
!
!     Description of Parameter
!      --Input--
!        IUNIT - an input parameter giving the logical unit number
!                to which error messages are to be sent.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  7 June 1978
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XSETUN
!***FIRST EXECUTABLE STATEMENT  XSETUN
      JUNK = J4SAVE(3,IUNIT,.TRUE.)
      JUNK = J4SAVE(5,1,.TRUE.)
      RETURN
      END

      REAL FUNCTION R1MACH(I) 
!***BEGIN PROLOGUE  R1MACH
!***DATE WRITTEN   790101   (YYMMDD)
!***REVISION DATE  831014   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS) 
!***PURPOSE  Returns single precision machine dependent constants
!***DESCRIPTION
!     From the book, "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988
!
!
!     R1MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function 
!     subroutine with one (input) argument, and can be called
!     as follows, for example
!
!          A = R1MACH(I)
!
!     where I=1,...,5.  The (output) value of A above is
!     determined by the (input) value of I.  The results for 
!     various values of I are discussed below.
!
!  Single-Precision Machine Constants
!  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. 
!  R1MACH(3) = B**(-T), the smallest relative spacing.
!  R1MACH(4) = B**(1-T), the largest relative spacing. 
!  R1MACH(5) = LOG10(B)
!***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
!                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
!                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
!                 PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  R1MACH
!
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
!
      REAL RMACH(5) 
!
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
!
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
!
!      DATA RMACH(1) / O"00014000000000000000" /
!      DATA RMACH(2) / O"37767777777777777777" /
!      DATA RMACH(3) / O"16404000000000000000" /
!      DATA RMACH(4) / O"16414000000000000000" /
!      DATA RMACH(5) / O"17164642023241175720" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA RMACH(1) / X'9000400000000000' /
!     DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /
!     DATA RMACH(3) / X'FFA3400000000000' /
!     DATA RMACH(4) / X'FFA4400000000000' /
!     DATA RMACH(5) / X'FFD04D104D427DE8' /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR THE CRAY 1
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
!
      DATA RMACH/1.18E-38,3.40E+38,0.595E-07,1.19E-07,0.30102999566/
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!
!     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
!       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
!       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
!
!     DATA SMALL(1) / X'00800000' /
!     DATA LARGE(1) / X'7F7FFFFF' /
!     DATA RIGHT(1) / X'33800000' /
!     DATA DIVER(1) / X'34000000' /
!     DATA LOG10(1) / X'3E9A209B' /
!
!
!     MACHINE CONSTANTS FOR THE VAX 11/780
!    (EXPRESSED IN INTEGER AND HEXADECIMAL)
!  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)&
        CALL XERROR ( 'R1MACH -- I OUT OF BOUNDS',25,1,2)
!
      R1MACH = RMACH(I)
      RETURN
!
      END 
      DOUBLE PRECISION FUNCTION D1MACH(I)
!***BEGIN PROLOGUE  D1MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  831014   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS)
!***PURPOSE  Returns double precision machine dependent constants
!***DESCRIPTION
!     From the book, "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988
!
!
!     D1MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function
!     subprogram with one (input) argument, and can be called
!     as follows, for example
!
!          D = D1MACH(I)
!
!     where I=1,...,5.  The (output) value of D above is
!     determined by the (input) value of I.  The results for
!     various values of I are discussed below.
!
!  Double-precision machine constants
!  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!  D1MACH( 3) = B**(-T), the smallest relative spacing.
!  D1MACH( 4) = B**(1-T), the largest relative spacing.
!  D1MACH( 5) = LOG10(B)
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  D1MACH
!
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
!
      DOUBLE PRECISION DMACH(5)
!
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
!
!      DATA SMALL(1) / O"00604000000000000000" /
!      DATA SMALL(2) / O"00000000000000000000" /
!
!      DATA LARGE(1) / O"37767777777777777777" /
!      DATA LARGE(2) / O"37167777777777777777" /
!
!      DATA RIGHT(1) / O"15604000000000000000" /
!      DATA RIGHT(2) / O"15000000000000000000" /
!
!      DATA DIVER(1) / O"15614000000000000000" /
!      DATA DIVER(2) / O"15010000000000000000" /
!
!      DATA LOG10(1) / O"17164642023241175717" /
!      DATA LOG10(2) / O"16367571421742254654" /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA SMALL(1) / X'9000400000000000' /
!     DATA SMALL(2) / X'8FD1000000000000' /
!
!     DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
!     DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
!
!     DATA RIGHT(1) / X'FF74400000000000' /
!     DATA RIGHT(2) / X'FF45000000000000' /
!
!     DATA DIVER(1) / X'FF75400000000000' /
!     DATA DIVER(2) / X'FF46000000000000' /
!
!     DATA LOG10(1) / X'FFD04D104D427DE7' /
!     DATA LOG10(2) / X'FFA17DE623E2566A' /
!
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /
!
!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /
!
!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /
!
!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /
!
!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /
!
!     MACHINE CONSTANTS FOR THE CRAY 1
!
!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /
!
!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777774B /
!
!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /
!
!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /
!
!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /
!
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
!     DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!     DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
!     DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
!     DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
!
!     MACHINE CONSTATNS FOR THE IBM PC FAMILY (D. KAHANER NBS)
!
      DATA DMACH/2.23D-308,1.79D+308,1.11D-16,2.22D-16,&
       0.301029995663981195D0/
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!     DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!     DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1),LOG10(2) / "177464202324, "476747767461 /
!
!
!     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
!       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
!       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
!
!      DATA SMALL(1),SMALL(2) / X'00100000', X'00000000' /
!      DATA LARGE(1),LARGE(2) / X'7FEFFFFF', X'FFFFFFFF' /
!      DATA RIGHT(1),RIGHT(2) / X'3CA00000', X'00000000' /
!      DATA DIVER(1),DIVER(2) / X'3CB00000', X'00000000' /
!      DATA LOG10(1),LOG10(2) / X'3FD34413', X'509F79FF' /
!
!
!     MACHINE CONSTANTS FOR VAX 11/780
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!
!    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
!
!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!
!    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)&
        CALL XERROR( 'D1MACH -- I OUT OF BOUNDS',25,1,2)
!
      D1MACH = DMACH(I)
      RETURN
!
      END
      INTEGER FUNCTION I1MACH(I)
!***BEGIN PROLOGUE  I1MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  840405   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS) 
!***PURPOSE  Returns integer machine dependent constants
!***DESCRIPTION
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!   These machine constant routines must be activated for
!   a particular environment.
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
!     I1MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function 
!     subroutine with one (input) argument, and can be called
!     as follows, for example
!
!          K = I1MACH(I)
!
!     where I=1,...,16.  The (output) value of K above is
!     determined by the (input) value of I.  The results for 
!     various values of I are discussed below.
!
!  I/O unit numbers.
!    I1MACH( 1) = the standard input unit.
!    I1MACH( 2) = the standard output unit.
!    I1MACH( 3) = the standard punch unit.
!    I1MACH( 4) = the standard error message unit.
!
!  Words.
!    I1MACH( 5) = the number of bits per integer storage unit.
!    I1MACH( 6) = the number of characters per integer storage unit.
!
!  Integers. 
!    assume integers are represented in the S-digit, base-A form
!
!               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!    I1MACH( 7) = A, the base.
!    I1MACH( 8) = S, the number of base-A digits.
!    I1MACH( 9) = A**S - 1, the largest magnitude. 
!
!  Floating-Point Numbers.
!    Assume floating-point numbers are represented in the T-digit,
!    base-B form
!               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               where 0 .LE. X(I) .LT. B for I=1,...,T,
!               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, the base.
!
!  Single-Precision
!    I1MACH(11) = T, the number of base-B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!  Double-Precision
!    I1MACH(14) = T, the number of base-B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment,
!  the desired set of DATA statements should be activated by
!  removing the C from column 1.  Also, the values of
!  I1MACH(1) - I1MACH(4) should be checked for consistency
!  with the local operating system.
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL 
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  I1MACH
!
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
!
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
!
!      DATA IMACH( 1) /    5 / 
!      DATA IMACH( 2) /    6 / 
!      DATA IMACH( 3) /    7 / 
!      DATA IMACH( 4) /    6 / 
!      DATA IMACH( 5) /   60 / 
!      DATA IMACH( 6) /   10 / 
!      DATA IMACH( 7) /    2 / 
!      DATA IMACH( 8) /   48 / 
!      DATA IMACH( 9) / O"00007777777777777777" /
!      DATA IMACH(10) /    2 / 
!      DATA IMACH(11) /   48 / 
!      DATA IMACH(12) / -974 / 
!      DATA IMACH(13) / 1070 / 
!      DATA IMACH(14) /   96 / 
!      DATA IMACH(15) / -927 / 
!      DATA IMACH(16) / 1070 / 
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
!
!     DATA IMACH( 1) /      5 /
!     DATA IMACH( 2) /      6 /
!     DATA IMACH( 3) /      7 /
!     DATA IMACH( 4) /      6 /
!     DATA IMACH( 5) /     64 /
!     DATA IMACH( 6) /      8 /
!     DATA IMACH( 7) /      2 /
!     DATA IMACH( 8) /     47 /
!     DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
!     DATA IMACH(10) /      2 /
!     DATA IMACH(11) /     47 /
!     DATA IMACH(12) / -28625 /
!     DATA IMACH(13) /  28718 /
!     DATA IMACH(14) /     94 /
!     DATA IMACH(15) / -28625 /
!     DATA IMACH(16) /  28718 /
!
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!     DATA IMACH( 1) /    5 / 
!     DATA IMACH( 2) /    6 / 
!     DATA IMACH( 3) /    7 / 
!     DATA IMACH( 4) /6LOUTPUT/
!     DATA IMACH( 5) /   60 / 
!     DATA IMACH( 6) /   10 / 
!     DATA IMACH( 7) /    2 / 
!     DATA IMACH( 8) /   48 / 
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /    2 / 
!     DATA IMACH(11) /   47 / 
!     DATA IMACH(12) / -929 / 
!     DATA IMACH(13) / 1070 / 
!     DATA IMACH(14) /   94 / 
!     DATA IMACH(15) / -929 / 
!     DATA IMACH(16) / 1069 / 
!
!     MACHINE CONSTANTS FOR THE CRAY 1
!
!     DATA IMACH( 1) /   100 /
!     DATA IMACH( 2) /   101 /
!     DATA IMACH( 3) /   102 /
!     DATA IMACH( 4) /   101 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    63 /
!     DATA IMACH( 9) /  777777777777777777777B /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -8189 /
!     DATA IMACH(13) /  8190 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -8099 /
!     DATA IMACH(16) /  8190 /
!
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  32 /
!     DATA IMACH( 6) /   4 /
!     DATA IMACH( 7) /  16 /
!     DATA IMACH( 8) /  31 /
!     DATA IMACH( 9) / Z7FFFFFFF /
!     DATA IMACH(10) /  16 /
!     DATA IMACH(11) /   6 /
!     DATA IMACH(12) / -64 /
!     DATA IMACH(13) /  63 /
!     DATA IMACH(14) /  14 /
!     DATA IMACH(15) / -64 /
!     DATA IMACH(16) /  63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
!
      DATA IMACH/5,6,0,6,32,4,2,31,2147483647,2,24,&
      -125,127,53,-1021,1023/
!               NOTE! I1MACH(3) IS NOT WELL DEFINED AND IS SET TO ZERO.
!
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!     DATA IMACH( 1) /    5 / 
!     DATA IMACH( 2) /    6 / 
!     DATA IMACH( 3) /    5 / 
!     DATA IMACH( 4) /    6 / 
!     DATA IMACH( 5) /   36 / 
!     DATA IMACH( 6) /    5 / 
!     DATA IMACH( 7) /    2 / 
!     DATA IMACH( 8) /   35 / 
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 / 
!     DATA IMACH(11) /   27 / 
!     DATA IMACH(12) / -128 / 
!     DATA IMACH(13) /  127 / 
!     DATA IMACH(14) /   54 / 
!     DATA IMACH(15) / -101 / 
!     DATA IMACH(16) /  127 / 
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!     DATA IMACH( 1) /    5 / 
!     DATA IMACH( 2) /    6 / 
!     DATA IMACH( 3) /    5 / 
!     DATA IMACH( 4) /    6 / 
!     DATA IMACH( 5) /   36 / 
!     DATA IMACH( 6) /    5 / 
!     DATA IMACH( 7) /    2 / 
!     DATA IMACH( 8) /   35 / 
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 / 
!     DATA IMACH(11) /   27 / 
!     DATA IMACH(12) / -128 / 
!     DATA IMACH(13) /  127 / 
!     DATA IMACH(14) /   62 / 
!     DATA IMACH(15) / -128 / 
!     DATA IMACH(16) /  127 / 
!
!
!     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
!       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
!       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    6 /
!      DATA IMACH( 4) /    0 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -125 /
!      DATA IMACH(13) /  128 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!
!     MACHINE CONSTANTS FOR THE VAX 11/780
!
!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    5 /
!     DATA IMACH(4) /    6 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) /2147483647 /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -127 /
!     DATA IMACH(13)/  127 /
!     DATA IMACH(14)/   56 /
!     DATA IMACH(15)/ -127 /
!     DATA IMACH(16)/  127 /
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) &
        CALL XERROR ( 'I1MACH -- I OUT OF BOUNDS',25,1,2)
!
      I1MACH=IMACH(I)
      RETURN
!
      END
