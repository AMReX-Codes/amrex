*DECK DVODE
      SUBROUTINE DVODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF,
     2            RPAR, IPAR)
      EXTERNAL F, JAC
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK, RPAR
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW,
     1        MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C DVODE: Variable-coefficient Ordinary Differential Equation solver,
C with fixed-leading-coefficient implementation.
C This version is in double precision.
C
C DVODE solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C DVODE is a package based on the EPISODE and EPISODEB packages, and
C on the ODEPACK user interface standard, with minor modifications.
C-----------------------------------------------------------------------
C Authors:
C               Peter N. Brown and Alan C. Hindmarsh
C               Center for Applied Scientific Computing, L-561
C               Lawrence Livermore National Laboratory
C               Livermore, CA 94551
C and
C               George D. Byrne
C               Illinois Institute of Technology
C               Chicago, IL 60616
C-----------------------------------------------------------------------
C References:
C
C 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
C    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
C    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.
C 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
C    Numerical Solution of Ordinary Differential Equations,"
C    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
C 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
C    for the Integration of Systems of Ordinary Differential
C    Equations," LLNL Report UCID-30112, Rev. 1, April 1977.
C 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
C    Package for the Integration of Systems of Ordinary Differential
C    Equations with Banded Jacobians," LLNL Report UCID-30132, April
C    1976.
C 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
C    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
C    North-Holland, Amsterdam, 1983, pp. 55-64.
C 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
C    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
C    Trans. Math. Software, 6 (1980), pp. 295-318.
C-----------------------------------------------------------------------
C Summary of usage.
C
C Communication between the user and the DVODE package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
C           DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use a method flag MF = 10.  If it is stiff, there are four standard
C choices for MF (21, 22, 24, 25), and DVODE requires the Jacobian
C matrix in some form.  In these cases (MF .gt. 0), DVODE will use a
C saved copy of the Jacobian matrix.  If this is undesirable because of
C storage limitations, set MF to the corresponding negative value
C (-21, -22, -24, -25).  (See full description of MF below.)
C The Jacobian matrix is regarded either as full (MF = 21 or 22),
C or banded (MF = 24 or 25).  In the banded case, DVODE requires two
C half-bandwidth parameters ML and MU.  These are, respectively, the
C widths of the lower and upper parts of the band, excluding the main
C diagonal.  Thus the band consists of the locations (i,j) with
C i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
C
C C. If the problem is stiff, you are encouraged to supply the Jacobian
C directly (MF = 21 or 24), but if this is not feasible, DVODE will
C compute it internally by difference quotients (MF = 22 or 25).
C If you are supplying the Jacobian, provide a subroutine of the form:
C           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
C           DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
C which supplies df/dy by loading PD as follows:
C     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
C the partial derivative of f(i) with respect to y(j).  (Ignore the
C ML and MU arguments in this case.)
C     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
C df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
C PD from the top down.
C     In either case, only nonzero elements need be loaded.
C
C D. Write a main program which calls subroutine DVODE once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by DVODE.  On the first call to DVODE, supply arguments as follows:
C F      = Name of subroutine for right-hand side vector f.
C          This name must be declared external in calling program.
C NEQ    = Number of first order ODEs.
C Y      = Array of initial values, of length NEQ.
C T      = The initial value of the independent variable.
C TOUT   = First point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = Relative tolerance parameter (scalar).
C ATOL   = Absolute tolerance parameter (scalar or array).
C          The estimated local error in Y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of Y at t = TOUT.
C ISTATE = Integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional input used.
C RWORK  = Real work array of length at least:
C             20 + 16*NEQ                      for MF = 10,
C             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
C             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
C LRW    = Declared length of RWORK (in user's DIMENSION statement).
C IWORK  = Integer work array of length at least:
C             30        for MF = 10,
C             30 + NEQ  for MF = 21, 22, 24, or 25.
C          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
C          and upper half-bandwidths ML,MU.
C LIW    = Declared length of IWORK (in user's DIMENSION statement).
C JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
C          If used, this name must be declared external in calling
C          program.  If not used, pass a dummy name.
C MF     = Method flag.  Standard values are:
C          10 for nonstiff (Adams) method, no Jacobian used.
C          21 for stiff (BDF) method, user-supplied full Jacobian.
C          22 for stiff method, internally generated full Jacobian.
C          24 for stiff method, user-supplied banded Jacobian.
C          25 for stiff method, internally generated banded Jacobian.
C RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL, RPAR, and IPAR.
C
C E. The output from the first call (or any call) is:
C      Y = Array of computed values of y(t) vector.
C      T = Corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if DVODE was successful, negative otherwise.
C          -1 means excess work done on this call. (Perhaps wrong MF.)
C          -2 means excess accuracy requested. (Tolerances too small.)
C          -3 means illegal input detected. (See printed message.)
C          -4 means repeated error test failures. (Check all input.)
C          -5 means repeated convergence failures. (Perhaps bad
C             Jacobian supplied or wrong choice of MF or tolerances.)
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C
C F. To continue the integration after a successful return, simply
C reset TOUT and call DVODE again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C EXAMPLE PROBLEM
C
C The following is a simple example problem, with the coding
C needed for its solution by DVODE.  The problem is from chemical
C kinetics, and consists of the following three rate equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C
C The following coding solves this problem with DVODE, using MF = 21
C and printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.
C At the end of the run, statistical quantities of interest are
C printed. (See optional output in the full description below.)
C To generate Fortran source code, replace C in column 1 with a blank
C in the coding below.
C
C     EXTERNAL FEX, JEX
C     DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
C     NEQ = 3
C     Y(1) = 1.0D0
C     Y(2) = 0.0D0
C     Y(3) = 0.0D0
C     T = 0.0D0
C     TOUT = 0.4D0
C     ITOL = 2
C     RTOL = 1.D-4
C     ATOL(1) = 1.D-8
C     ATOL(2) = 1.D-14
C     ATOL(3) = 1.D-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 67
C     LIW = 33
C     MF = 21
C     DO 40 IOUT = 1,12
C       CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
C    1            IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
C       WRITE(6,20)T,Y(1),Y(2),Y(3)
C 20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
C       IF (ISTATE .LT. 0) GO TO 80
C 40    TOUT = TOUT*10.
C     WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),
C    1            IWORK(20),IWORK(21),IWORK(22)
C 60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
C    1       '   No. J-s =',I4,'   No. LU-s =',I4/
C    2       '  No. nonlinear iterations =',I4/
C    3       '  No. nonlinear convergence failures =',I4/
C    4       '  No. error test failures =',I4/)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt: ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
C     DOUBLE PRECISION RPAR, T, Y, YDOT
C     DIMENSION Y(NEQ), YDOT(NEQ)
C     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
C     YDOT(3) = 3.D7*Y(2)*Y(2)
C     YDOT(2) = -YDOT(1) - YDOT(3)
C     RETURN
C     END
C
C     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
C     DOUBLE PRECISION PD, RPAR, T, Y
C     DIMENSION Y(NEQ), PD(NRPD,NEQ)
C     PD(1,1) = -.04D0
C     PD(1,2) = 1.D4*Y(3)
C     PD(1,3) = 1.D4*Y(2)
C     PD(2,1) = .04D0
C     PD(2,3) = -PD(1,3)
C     PD(3,2) = 6.D7*Y(2)
C     PD(2,2) = -PD(1,2) - PD(3,2)
C     RETURN
C     END
C
C The following output was obtained from the above program on a
C Cray-1 computer with the CFT compiler.
C
C At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
C At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
C At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
C At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
C At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
C At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
C At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
C At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
C At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
C At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
C At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
C At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
C
C No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
C  No. nonlinear iterations = 831
C  No. nonlinear convergence failures =   0
C  No. error test failures =  22
C-----------------------------------------------------------------------
C Full description of user interface to DVODE.
C
C The user interface to DVODE consists of the following parts.
C
C i.   The call sequence to subroutine DVODE, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is
C        * a description of optional input available through the
C          call sequence,
C        * a description of optional output (in the work arrays), and
C        * instructions for interrupting and restarting a solution.
C
C ii.  Descriptions of other routines in the DVODE package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      COMMON, and obtain specified derivatives of the solution y(t).
C
C iii. Descriptions of COMMON blocks to be declared in overlay
C      or similar environments.
C
C iv.  Description of two routines in the DVODE package, either of
C      which the user may replace with his own version, if desired.
C      these relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part i.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional input and optional output.  (The term output here refers
C to the return from subroutine DVODE to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 in the input.
C
C The descriptions of the call arguments are as follows.
C
C F      = The name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
C               DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared EXTERNAL in the calling program.
C
C          Subroutine F may access user-defined real and integer
C          work arrays RPAR and IPAR, which are to be dimensioned
C          in the main program.
C
C          If quantities computed in the F routine are needed
C          externally to DVODE, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use DVINDY instead.
C
C NEQ    = The size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may not be increased during the problem, but
C          can be decreased (with ISTATE = 3 in the input).
C
C Y      = A real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  In the output, Y contains the computed solution
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F and JAC.
C
C T      = The independent variable.  In the input, T is used only on
C          the first call, as the initial point of the integration.
C          In the output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = The next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal t interval, whose endpoints are
C          TCUR - HU and TCUR.  (See optional output, below, for
C          TCUR and HU.)
C
C ITOL   = An indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = A relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = An absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C          The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector e = (e(i)) of estimated local errors
C          in Y, according to an inequality of the form
C                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
C          and the rms-norm (root-mean-square norm) here is
C          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part iv below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = An index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at T = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          In the input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
C             and any of the optional input except H0.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful to include
C          the initial conditions in the output.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 in the input.
C
C          In the output, ISTATE has the following values and meanings.
C           1  means nothing was done, as TOUT was equal to T with
C              ISTATE = 1 in the input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again.
C              (The excess work step counter will be reset to 0.)
C              In addition, the user may increase MXSTEP to avoid
C              this error return.  (See optional input below.)
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other input, before
C          calling the solver again.
C
C IOPT   = An integer flag to specify whether or not any optional
C          input is being used on this call.  Input only.
C          The optional input is listed separately below.
C          IOPT = 0 means no optional input is being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means optional input is being used.
C
C RWORK  = A real working array (double precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LWM = length of work space for matrix-related data:
C          LWM = 0             if MITER = 0,
C          LWM = 2*NEQ**2 + 2  if MITER = 1 or 2, and MF.gt.0,
C          LWM = NEQ**2 + 2    if MITER = 1 or 2, and MF.lt.0,
C          LWM = NEQ + 2       if MITER = 3,
C          LWM = (3*ML+2*MU+2)*NEQ + 2 if MITER = 4 or 5, and MF.gt.0,
C          LWM = (2*ML+MU+1)*NEQ + 2   if MITER = 4 or 5, and MF.lt.0.
C          (See the MF description for METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is:
C             20 + 16*NEQ                    for MF = 10,
C             22 + 16*NEQ + 2*NEQ**2         for MF = 11 or 12,
C             22 + 16*NEQ + NEQ**2           for MF = -11 or -12,
C             22 + 17*NEQ                    for MF = 13,
C             22 + 18*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
C             22 + 17*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
C             20 +  9*NEQ                    for MF = 20,
C             22 +  9*NEQ + 2*NEQ**2         for MF = 21 or 22,
C             22 +  9*NEQ + NEQ**2           for MF = -21 or -22,
C             22 + 10*NEQ                    for MF = 23,
C             22 + 11*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
C             22 + 10*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
C          The first 20 words of RWORK are reserved for conditional
C          and optional input and optional output.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = The length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = An integer work array.  The length of IWORK must be at least
C             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
C             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
C          The first 30 words of IWORK are reserved for conditional and
C          optional input and optional output.
C
C          The following 2 words in IWORK are conditional input:
C            IWORK(1) = ML     These are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       banded Jacobian, excluding the main diagonal.
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if MITER is 4 or 5, and
C                       ignored otherwise.  ML and MU may in fact be
C                       the band parameters for a matrix to which
C                       df/dy is only approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to DVODE
C for the same problem, except possibly for the conditional and
C optional input, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside DVODE between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
C          compute the Jacobian matrix, df/dy, as a function of
C          the scalar t and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
C                               RPAR, IPAR)
C               DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
C          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
C          PD is to be loaded with partial derivatives (elements of the
C          Jacobian matrix) in the output.  PD must be given a first
C          dimension of NROWPD.  T and Y have the same meaning as in
C          Subroutine F.
C               In the full matrix case (MITER = 1), ML and MU are
C          ignored, and the Jacobian is to be loaded into PD in
C          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C               In the band matrix case (MITER = 4), the elements
C          within the band are to be loaded into PD in columnwise
C          manner, with diagonal lines of df/dy loaded into the rows
C          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
C          ML and MU are the half-bandwidth parameters. (See IWORK).
C          The locations in PD in the two triangular areas which
C          correspond to nonexistent matrix elements can be ignored
C          or loaded arbitrarily, as they are overwritten by DVODE.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               In either case, PD is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to F with the same
C          arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user COMMON block by F and not recomputed by JAC,
C          if desired.  Also, JAC may alter the Y array, if desired.
C          JAC must be declared external in the calling program.
C               Subroutine JAC may access user-defined real and integer
C          work arrays, RPAR and IPAR, whose dimensions are set by the
C          user in the main program.
C
C MF     = The method flag.  Used only for input.  The legal values of
C          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
C          -11, -12, -14, -15, -21, -22, -24, -25.
C          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
C          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
C            JSV =  1 means a copy of the Jacobian is saved for reuse
C                     in the corrector iteration algorithm.
C            JSV = -1 means a copy of the Jacobian is not saved
C                     (valid only for MITER = 1, 2, 4, or 5).
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on backward
C                     differentiation formulas (BDF-s).
C          MITER indicates the corrector iteration method:
C            MITER = 0 means functional iteration (no Jacobian matrix
C                      is involved).
C            MITER = 1 means chord iteration with a user-supplied
C                      full (NEQ by NEQ) Jacobian.
C            MITER = 2 means chord iteration with an internally
C                      generated (difference quotient) full Jacobian
C                      (using NEQ extra calls to F per df/dy value).
C            MITER = 3 means chord iteration with an internally
C                      generated diagonal Jacobian approximation
C                      (using 1 extra call to F per df/dy evaluation).
C            MITER = 4 means chord iteration with a user-supplied
C                      banded Jacobian.
C            MITER = 5 means chord iteration with an internally
C                      generated banded Jacobian (using ML+MU+1 extra
C                      calls to F per df/dy evaluation).
C          If MITER = 1 or 4, the user must supply a subroutine JAC
C          (the name is arbitrary) as described above under JAC.
C          For other values of MITER, a dummy argument can be used.
C
C RPAR     User-specified array used to communicate real parameters
C          to user-supplied subroutines.  If RPAR is a vector, then
C          it must be dimensioned in the user's main program.  If it
C          is unused or it is a scalar, then it need not be
C          dimensioned.
C
C IPAR     User-specified array used to communicate integer parameter
C          to user-supplied subroutines.  The comments on dimensioning
C          RPAR apply to IPAR.
C-----------------------------------------------------------------------
C Optional Input.
C
C The following is a list of the optional input provided for in the
C call sequence.  (See also Part ii.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of this input requires IOPT = 1, and in that
C case all of this input is examined.  A value of zero for any
C of these optional input variables will cause the default value to be
C used.  Thus to use a subset of the optional input, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C NAME    LOCATION      MEANING AND DEFAULT VALUE
C
C H0      RWORK(5)  The step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  The maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  The minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  The maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C-----------------------------------------------------------------------
C Optional Output.
C
C As optional additional output from DVODE, the variables listed
C below are quantities related to the performance of DVODE
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of this output is defined
C on any successful return from DVODE, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, output relevant to the error will be defined,
C as noted below.
C
C NAME    LOCATION      MEANING
C
C HU      RWORK(11) The step size in t last used (successfully).
C
C HCUR    RWORK(12) The step size to be attempted on the next step.
C
C TCUR    RWORK(13) The current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  In the output,
C                   TCUR will always be at least as far from the
C                   initial value of t as the current argument T,
C                   but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) The number of steps taken for the problem so far.
C
C NFE     IWORK(12) The number of f evaluations for the problem so far.
C
C NJE     IWORK(13) The number of Jacobian evaluations so far.
C
C NQU     IWORK(14) The method order last used (successfully).
C
C NQCUR   IWORK(15) The order to be attempted on the next step.
C
C IMXER   IWORK(16) The index of the component of largest magnitude in
C                   the weighted local error vector ( e(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) The length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) The length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NLU     IWORK(19) The number of matrix LU decompositions so far.
C
C NNI     IWORK(20) The number of nonlinear (Newton) iterations so far.
C
C NCFN    IWORK(21) The number of convergence failures of the nonlinear
C                   solver so far.
C
C NETF    IWORK(22) The number of error test failures of the integrator
C                   so far.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional output.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C NAME    BASE ADDRESS      DESCRIPTION
C
C YH      21             The Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the
C                        solution, evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
C                        corrections on each step, scaled in the output
C                        to represent the estimated local error in Y
C                        on the last step.  This is the vector e in
C                        the description of the error control.  It is
C                        defined only on a successful return from DVODE.
C
C-----------------------------------------------------------------------
C Interrupting and Restarting
C
C If the integration of a given problem by DVODE is to be
C interrrupted and then later continued, such as when restarting
C an interrupted run or alternating between two or more ODE problems,
C the user should save, following the return from the last DVODE call
C prior to the interruption, the contents of the call sequence
C variables and internal COMMON blocks, and later restore these
C values before the next DVODE call for that problem.  To save
C and restore the COMMON blocks, use subroutine DVSRCO, as
C described below in part ii.
C
C In addition, if non-default values for either LUN or MFLAG are
C desired, an extra call to XSETUN and/or XSETF should be made just
C before continuing the integration.  See Part ii below for details.
C
C-----------------------------------------------------------------------
C Part ii.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with DVODE.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     FORM OF CALL                  FUNCTION
C  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
C                             output of messages from DVODE, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C  CALL XSETF(MFLAG)          Set a flag to control the printing of
C                             messages by DVODE.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C  CALL DVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
C                             the internal COMMON blocks used by
C                             DVODE. (See Part iii below.)
C                             RSAV must be a real array of length 49
C                             or more, and ISAV must be an integer
C                             array of length 40 or more.
C                             JOB=1 means save COMMON into RSAV/ISAV.
C                             JOB=2 means restore COMMON from RSAV/ISAV.
C                                DVSRCO is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with DVODE.
C
C  CALL DVINDY(,,,,,)         Provide derivatives of y, of various
C        (See below.)         orders, at a specified point T, if
C                             desired.  It may be called only after
C                             a successful return from DVODE.
C
C The detailed instructions for using DVINDY are as follows.
C The form of the call is:
C
C  CALL DVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = Value of independent variable where answers are desired
C             (normally the same as the T last returned by DVODE).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional output for TCUR and HU.)
C K         = Integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional output).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by DVODE directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with DVINDY.
C RWORK(21) = The base address of the history array YH.
C NYH       = Column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = A real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = Integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part iii.  COMMON Blocks.
C If DVODE is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to DVODE,
C   (2) the two internal COMMON blocks
C         /DVOD01/  of length  81  (48 double precision words
C                         followed by 33 integer words),
C         /DVOD02/  of length  9  (1 double precision word
C                         followed by 8 integer words),
C
C If DVODE is used on a system in which the contents of internal
C COMMON blocks are not preserved between calls, the user should
C declare the above two COMMON blocks in his main program to insure
C that their contents are preserved.
C
C-----------------------------------------------------------------------
C Part iv.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the DVODE package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) DEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the DVODE call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by DEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparison with
C errors in Y(i).  The EWT array returned by DEWSET is passed to the
C DVNORM routine (See below.), and also used by DVODE in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of DEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C Optional Output.  In DEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of h**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in DEWSET the statements:
C     DOUBLE PRECISION RVOD, H, HU
C     COMMON /DVOD01/ RVOD(48), IVOD(33)
C     COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C     NQ = IVOD(28)
C     H = RVOD(21)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) DVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = DVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
C DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by subroutine DEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of DVNORM suitable for use in the error control in DVODE.
C None of the arguments should be altered by DVNORM.
C For example, a user-supplied DVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of Y.
C-----------------------------------------------------------------------
C REVISION HISTORY (YYYYMMDD)
C  19890615  Date Written.  Initial release.
C  19890922  Added interrupt/restart ability, minor changes throughout.
C  19910228  Minor revisions in line format,  prologue, etc.
C  19920227  Modifications by D. Pang:
C            (1) Applied subgennam to get generic intrinsic names.
C            (2) Changed intrinsic names to generic in comments.
C            (3) Added *DECK lines before each routine.
C  19920721  Names of routines and labeled Common blocks changed, so as
C            to be unique in combined single/double precision code (ACH).
C  19920722  Minor revisions to prologue (ACH).
C  19920831  Conversion to double precision done (ACH).
C  19921106  Fixed minor bug: ETAQ,ETAQM1 in DVSTEP SAVE statement (ACH).
C  19921118  Changed LUNSAV/MFLGSV to IXSAV (ACH).
C  19941222  Removed MF overwrite; attached sign to H in estimated second 
C            deriv. in DVHIN; misc. comment changes throughout (ACH).
C  19970515  Minor corrections to comments in prologue, DVJAC (ACH).
C  19981111  Corrected Block B by adding final line, GO TO 200 (ACH).
C  20020430  Various upgrades (ACH): Use ODEPACK error handler package.
C            Replaced D1MACH by DUMACH.  Various changes to main
C            prologue and other routine prologues.
C-----------------------------------------------------------------------
C Other Routines in the DVODE Package.
C
C In addition to subroutine DVODE, the DVODE package includes the
C following subroutines and function routines:
C  DVHIN     computes an approximate step size for the initial step.
C  DVINDY    computes an interpolated value of the y vector at t = TOUT.
C  DVSTEP    is the core integrator, which does one step of the
C            integration and the associated error control.
C  DVSET     sets all method coefficients and test constants.
C  DVNLSD    solves the underlying nonlinear system -- the corrector.
C  DVJAC     computes and preprocesses the Jacobian matrix J = df/dy
C            and the Newton iteration matrix P = I - (h/l1)*J.
C  DVSOL     manages solution of linear system in chord iteration.
C  DVJUST    adjusts the history array on a change of order.
C  DEWSET    sets the error weight vector EWT before each step.
C  DVNORM    computes the weighted r.m.s. norm of a vector.
C  DVSRCO    is a user-callable routine to save and restore
C            the contents of the internal COMMON blocks.
C  DACOPY    is a routine to copy one two-dimensional array to another.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C            systems of linear algebraic equations.
C  DGBFA and DGBSL   are routines from LINPACK for solving banded
C            linear systems.
C  DAXPY, DSCAL, and DCOPY are basic linear algebra modules (BLAS).
C  DUMACH    sets the unit roundoff of the machine.
C  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all
C            error messages and warnings.  XERRWD is machine-dependent.
C Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
C
      include "vode.H"
C
C Type declarations for local variables --------------------------------
C
      EXTERNAL DVNLSD
      LOGICAL IHIT
      DOUBLE PRECISION ATOLI, BIG, EWTI, FOUR, H0, HMAX, HMX, HUN, ONE,
     1   PT2, RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW,
     1   LENWM, LF0, MBAND, MFA, ML, MORD, MU, MXHNL0, MXSTP0, NITER, 
     2   NSLAST
      CHARACTER*80 MSG
C
C Type declaration for function subroutines called ---------------------
C
      DOUBLE PRECISION DUMACH, DVNORM
C
      DIMENSION MORD(2)

      DATA  MORD(1) /12/, MORD(2) /5/
      PARAMETER(MXSTP0 = 500, MXHNL0 = 10)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0)
      PARAMETER(PT2 = 0.2D0, HUN = 100.0D0)
C-----------------------------------------------------------------------
C The following internal COMMON blocks contain variables which are
C communicated between subroutines in the DVODE package, or which are
C to be saved between calls to DVODE.
C In each block, real variables precede integers.
C The block /DVOD01/ appears in subroutines DVODE, DVINDY, DVSTEP,
C DVSET, DVNLSD, DVJAC, DVSOL, DVJUST and DVSRCO.
C The block /DVOD02/ appears in subroutines DVODE, DVINDY, DVSTEP,
C DVNLSD, DVJAC, and DVSRCO.
C
C The variables stored in the internal COMMON blocks are as follows:
C
C ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
C CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
C CONP   = The saved value of TQ(5).
C CRATE  = Estimated corrector convergence rate constant.
C DRC    = Relative change in H*RL1 since last DVJAC call.
C EL     = Real array of integration coefficients.  See DVSET.
C ETA    = Saved tentative ratio of new to old H.
C ETAMAX = Saved maximum value of ETA to be allowed.
C H      = The step size.
C HMIN   = The minimum absolute value of the step size H to be used.
C HMXI   = Inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C HNEW   = The step size to be attempted on the next step.
C HSCAL  = Stepsize in scaling of YH array.
C PRL1   = The saved value of RL1.
C RC     = Ratio of current H*RL1 to value on last DVJAC call.
C RL1    = The reciprocal of the coefficient EL(1).
C TAU    = Real vector of past NQ step sizes, length 13.
C TQ     = A real vector of length 5 in which DVSET stores constants
C          used for the convergence test, the error test, and the
C          selection of H at a new order.
C TN     = The independent variable, updated on each step taken.
C UROUND = The machine unit roundoff.  The smallest positive real number
C          such that  1.0 + UROUND .ne. 1.0
C ICF    = Integer flag for convergence failure in DVNLSD:
C            0 means no failures.
C            1 means convergence failure with out of date Jacobian
C                   (recoverable error).
C            2 means convergence failure with current Jacobian or
C                   singular matrix (unrecoverable error).
C INIT   = Saved integer flag indicating whether initialization of the
C          problem has been done (INIT = 1) or not.
C IPUP   = Saved flag to signal updating of Newton matrix.
C JCUR   = Output flag from DVJAC showing Jacobian status:
C            JCUR = 0 means J is not current.
C            JCUR = 1 means J is current.
C JSTART = Integer flag used as input to DVSTEP:
C            0  means perform the first step.
C            1  means take a new step continuing from the last.
C            -1 means take the next step with a new value of MAXORD,
C                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
C          On return, DVSTEP sets JSTART = 1.
C JSV    = Integer flag for Jacobian saving, = sign(MF).
C KFLAG  = A completion code from DVSTEP with the following meanings:
C               0      the step was succesful.
C              -1      the requested error could not be achieved.
C              -2      corrector convergence could not be achieved.
C              -3, -4  fatal error in VNLS (can not occur here).
C KUTH   = Input flag to DVSTEP showing whether H was reduced by the
C          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
C L      = Integer variable, NQ + 1, current order plus one.
C LMAX   = MAXORD + 1 (used for dimensioning).
C LOCJS  = A pointer to the saved Jacobian, whose storage starts at
C          WM(LOCJS), if JSV = 1.
C LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
C          to segments of RWORK and IWORK.
C MAXORD = The maximum order of integration method to be allowed.
C METH/MITER = The method flags.  See MF.
C MSBJ   = The maximum number of steps between J evaluations, = 50.
C MXHNIL = Saved value of optional input MXHNIL.
C MXSTEP = Saved value of optional input MXSTEP.
C N      = The number of first-order ODEs, = NEQ.
C NEWH   = Saved integer to flag change of H.
C NEWQ   = The method order to be used on the next step.
C NHNIL  = Saved counter for occurrences of T + H = T.
C NQ     = Integer variable, the current integration method order.
C NQNYH  = Saved value of NQ*NYH.
C NQWAIT = A counter controlling the frequency of order changes.
C          An order change is about to be considered if NQWAIT = 1.
C NSLJ   = The number of steps taken as of the last Jacobian update.
C NSLP   = Saved value of NST as of last Newton matrix update.
C NYH    = Saved value of the initial value of NEQ.
C HU     = The step size in t last used.
C NCFN   = Number of nonlinear convergence failures so far.
C NETF   = The number of error test failures of the integrator so far.
C NFE    = The number of f evaluations for the problem so far.
C NJE    = The number of Jacobian evaluations so far.
C NLU    = The number of matrix LU decompositions so far.
C NNI    = Number of nonlinear iterations so far.
C NQU    = The method order last used.
C NST    = The number of steps taken for the problem so far.
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .NE. 1) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all input and various initializations.
C
C First check legality of the non-optional input NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = SIGN(1,MF)
      MFA = ABS(MF)
      METH = MFA/10
      MITER = MFA - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional input. ---------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      JCO = MAX(0,JSV)
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        LENWM = 2 + (1 + JCO)*N*N
        LOCJS = N*N + 3
      ENDIF
      IF (MITER .EQ. 3) LENWM = 2 + N
      IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MBAND = ML + MU + 1
        LENP = (MBAND + ML)*N
        LENJ = MBAND*N
        LENWM = 2 + LENP + JCO*LENJ
        LOCJS = LENP + 3
        ENDIF
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 30 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL DCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NJE = 0
      NNI = 0
      NCFN = 0
      NETF = 0
      NLU = 0
      NSLJ = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR)
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      CALL DCOPY (N, Y, 1, RWORK(LYH), 1)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
C Call DVHIN to set initial step size H0 to be attempted. --------------
      CALL DVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,
     1   UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,
     2   NITER, IER)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL DSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DVSTEP.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DVODE--  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
C-----------------------------------------------------------------------
C CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
C              WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
C-----------------------------------------------------------------------
      CALL DVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), Y, RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, F, DVNLSD, RPAR, IPAR)
      KGO = 1 - KFLAG
C Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
C  KFLAG .eq. 0,   -1,  -2
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DVODE.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional output is loaded into the work
C arrays before returning.
C-----------------------------------------------------------------------
 400  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = HNEW
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NEWQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C if there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH, and T is set to TN.
C The optional output is loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine:   see TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional output. --------------------------------
 580  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.   If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DVODE--  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DVODE--  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
      CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'DVODE--  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'DVODE--  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'DVODE--  MF (=I1) illegal     '
      CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 609  MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
      GO TO 700
 610  MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DVODE--  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine:   see TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
      CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
C
 700  CONTINUE
      ISTATE = -3
      RETURN
C
 800  MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      RETURN
C----------------------- End of Subroutine DVODE -----------------------
      END
