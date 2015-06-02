# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( void );
void assemble ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int n, int indx[], int node[], int nu, int nl, int nquad, 
  int nmax, double ul, double ur, double wquad[], double xn[], double xquad[] );
double ff ( double x );
void geometry ( double h[], int ibc, int indx[], int n, int nl, int nmax, 
  int node[], int nquad, int *nu, double wquad[], double xn[], double xquad[] );
double get_alpha ( void );
double get_beta ( void );
int get_problem ( void );
void init ( int *ibc, int n, double *tol, double *ul, double *ur, double *xl, 
  double xn[], double *xr );
void output ( double f[], int ibc, int indx[], int n, int nu, double ul, 
  double ur, double xn[] );
void phi ( int il, double x, double *phii, double *phiix, double xleft, 
  double xrite );
double pp ( double x );
void prsys ( double adiag[], double aleft[], double arite[], double f[], 
  int nu );
double qq ( double x );
double r8_max ( double x, double y );
void solve ( double adiag[], double aleft[], double arite[], double f[], 
  int nu );
void solvex ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int ibc, int indx[], int kount, int n, int nl, int nmax, 
  int node[], int nquad, int *nu, double ul, double ur, double wquad[], 
  double xn[], double xquad[] );
void solvey ( double eta[], double f[], double h[], int n, int nu, double ul, 
  double ur, double xn[] );
int subdiv ( double eta[], int kount, int *n, int nmax, double tol, 
  double xn[] );
void timestamp ( void );
double uexact ( double x );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_ADAPTIVE.
//
//  Discussion:
//
//    FEM1D_ADAPTIVE solves a 1D problem using an adaptive finite element method.
//
//    The equation to be treated is:
//
//      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
//
//    by the finite-element method using piecewise linear basis
//    functions.
//
//    An adaptive method is used to try to reduce the maximum
//    error by refining the mesh in certain places.
//
//    Here U is an unknown scalar function of X defined on the
//    interval [XL,XR], and P, Q and F are given functions of X.
//
//    The values of U at XL and XR are also specified.
//
//    The interval [XL,XR] is "meshed" with N+1 points,
//
//      XN(0) = XL, XN(1) = XL+H, XN(2) = XL+2*H, ..., XN(N) = XR.
//
//    This creates N subintervals, with interval I having endpoints 
//    XN(I-1) and XN(I).
//
//
//    The algorithm tries to guarantee a certain amount
//    of accuracy by examining the current solution, estimating the error
//    in each subinterval, and, if necessary, subdividing one or more
//    subintervals and repeating the calculation.
//
//    We can think of the adaptive part of the algorithm as a refined
//    problem.  The program re-solves the problem on the pair of
//    intervals J and J+1, which extend from node J-1 to node J+1.
//    The values of U that were just computed at nodes J-1 and J+1
//    will be used as the boundary values for this refined problem.
//    The intervals J and J+1 will each be evenly divided into NY
//    smaller subintervals.  This boundary value problem is solved,
//    and the derivatives of the original and refined solutions are
//    then compared to get an estimate of the error.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    double ADIAG(NU).
//    ADIAG(I) is the "diagonal" coefficient of the I-th
//    equation in the linear system.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    double ALEFT(NU).
//    ALEFT(I) is the "left hand" coefficient of the I-th
//    equation in the linear system.  That is, ALEFT(I) is the
//    coefficient of the (I-1)-th unknown in the I-th equation.
//    There is no value in ALEFT(1), since the first equation
//    does not refer to a "0-th" unknown.
//
//    double ARITE(NU).
//    ARITE(I) is the "right hand" coefficient of the I-th
//    equation in the linear system.  ARITE(I) is the coefficient
//    of the (I+1)-th unknown in the I-th equation.  There is
//    no value in ARITE(NU) because the NU-th equation does not
//    refer to an "NU+1"-th unknown.
//
//    double ETA(N).
//    ETA(I) is the error estimate for interval I.  It is computed
//    as the sum of two quantities, one associated with the left
//    and one with the right node of the interval.
//
//    double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    double FY(M).
//    FY is the right hand side of the linear system of the refined
//    problem.
//
//    double H(N)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    double HY(M).
//    HY(I) is the length of subinterval I in the refined problem.
//
//    int IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    int IBCY.
//    IBCY declares the boundary conditions for the refined problem
//    which should always be that the value of U is specified at
//    both the left and right endpoints.  This corresponds to a
//    value of IBCY = 3.
//
//    int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    int INDY(0:M).
//    INDY(I) records the index of the unknown associated with
//    node I for the refined problem.
//
//    int JADD(N).
//    JADD(I) is 1 if the error estimates show that interval I
//    should be subdivided.
//
//    int KOUNT, the number of adaptive steps that have been taken.
//
//    int M.
//    M is the number of subintervals used in the refined problem.
//    M is equal to NY for computations centered at node 0 or node N,
//    and otherwise, M is equal to 2*NY.
//
//    int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    int NMAX, the maximum number of unknowns that can be handled.
//
//    int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    int NODEY(NL,M).
//    NODEY performs the same function for the refined problem that
//    NODE performs for the full problem, recording the node numbers
//    associated with a particular subinterval.
//
//    int NQUAD
//    The number of quadrature points used in a subinterval.
//
//    int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    int NUY.
//    The number of unknowns in the refined problem.
//
//    int NY.
//    NY is the number of subintervals into which a given interval
//    will be subdivided, before solving the refined probelm.
//
//    int PROBLEM, chooses the problem to be solved.
//    The user must choose this value by setting it in routine GETPRB.
//    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
//    The program should find the solution exactly, and the
//    adaptive code should find that there is no reason to
//    subdivide any interval.
//    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
//    This problem should find the solution exactly, and
//    the adaptive code should again find there is nothing
//    to do.
//    *3, u = sin(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*sin(pi*x/2), 
//    ul = 0, ur = 1.
//    *4, u = cos(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*cos(pi*x/2), 
//    ul = 1, ur = 0.
//    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1, ibc = 3, 
//    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
//    ul = 0, ur = 1/((beta+2)*(beta+1))
//    (beta must be greater than -2, and not equal to -1)
//    *6: u = atan((x-0.5)/alpha), p = 1, q = 0, ibc = 3, 
//    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
//    ul = u(0), ur = u(1)
//
//    int STATUS, reports status of subdivision.
//    0, a new subdivision was carried out.
//    1, no more subdivisions are needed.
//    -1, no more subdivisions can be carried out.
//
//    double TOL.
//    A tolerance that is used to determine whether the estimated
//    error in an interval is so large that it should be subdivided
//    and the problem solved again.
//
//    double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    double WQUAD(NQUAD).
//    WQUAD(I) is the weight associated with the I-th point
//    of an NQUAD point Gaussian quadrature rule.
//
//    double XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    double XQUAD(NQUAD,NMAX), the I-th quadrature point
//    in interval J.
//
//    double XQUADY(NQUAD,NMAY ), the I-th quadrature point
//    in subinterval J of the refined problem.
//
//    double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
//    Workspace, double precision XT(0:NMAX), used to compute a new
//    set of nodes.
//
//    double YN(0:M).
//    YN(I) is the location of the I-th node in the refined
//    problem.
//
{
# define NL 2
# define NMAX 30
# define NQUAD 2

  double alpha;
  double adiag[NMAX+1];
  double aleft[NMAX+1];
  double arite[NMAX+1];
  double beta;
  double eta[NMAX];
  double f[NMAX+1];
  double h[NMAX];
  int ibc;
  int indx[NMAX+1];
  int jadd[NMAX];
  int kount;
  int n;
  int node[NL*NMAX];
  int nu;
  int problem;
  int status;
  double tol;
  double ul;
  double ur;
  double wquad[NQUAD];
  double xl;
  double xn[NMAX+1];
  double xquad[NQUAD*NMAX];
  double xr;
  double xt[NMAX+1];

  timestamp ( );

  cout << "\n";
  cout << "FEM1D_ADAPTIVE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "Solve the two-point boundary value problem:\n";
  cout << "\n";
  cout << "  -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)\n";
  cout << "\n";
  cout << "on the interval [0,1], specifying the value\n";
  cout << "of U at each endpoint.\n";
  cout << "\n";
  cout << "  The number of basis functions per element is " << NL << "\n";
  cout << "\n";
  cout << "  The number of quadrature points per element is " << NQUAD << "\n";

  problem = get_problem ( );

  cout << "\n";
  cout << "  Problem index = " << problem << "\n";
  cout << "\n";

  if ( problem == 1 )
  {
    cout << "\n";
    cout << "  \"Linear\" problem:\n";
    cout << "  (No refinement needed)\n";
    cout << "\n";
    cout << "  U(X) =  X\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  0.0\n";
    cout << "  F(X) =  0.0\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  0.0\n";
    cout << "  UR   =  1.0\n";
  }
  else if ( problem == 2 )
  {
    cout << "\n";
    cout << "  \"Quadratic\" problem:\n";
    cout << "  (No refinement needed)\n";
    cout << "\n";
    cout << "  U(X) =  X*X\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  0.0\n";
    cout << "  F(X) = -2.0\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  0.0\n";
    cout << "  UR   =  1.\n";
  }
  else if ( problem == 3 )
  {
    cout << "\n";
    cout << "  \"SINE\" problem:\n";
    cout << "\n";
    cout << "  U(X) =  SIN(PI*X/2)\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  0.0\n";
    cout << "  F(X) =  PI*PI*SIN(PI*X/2)/4\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  0.0\n";
    cout << "  UR   =  1.0\n";
  }
  else if ( problem == 4 )
  {
    cout << "\n";
    cout << "  \"COSINE\" problem:\n";
    cout << "\n";
    cout << "  U(X) =  COS(PI*X/2)\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  0.0\n";
    cout << "  F(X) =  PI*PI*COS(PI*X/2)/4\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  0.0\n";
    cout << "  UR   =  1.0\n";
  }
  else if ( problem == 5 )
  {
    beta = get_beta ( );

    cout << "\n";
    cout << "  \"RHEINBOLDT\" problem:\n";
    cout << "\n";
    cout << "  U(X) =  X**(B+2)/((B+2)*(B+1))\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  1.0\n";
    cout << "  F(X) =  -X**B+(X**B+2))/((B+2)*(B+1))\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  0.0\n";
    cout << "  UR   =  1/((B+2)*(B+1))\n";
    cout << "  B    = " << beta << "\n";
  }
  else if ( problem == 6 )
  {
    alpha = get_alpha ( );

    cout << "\n";
    cout << "  \"ARCTAN\" problem:\n";
    cout << "\n";
    cout << "  U(X) =  ATAN((X-0.5)/A)\n";
    cout << "  P(X) =  1.0\n";
    cout << "  Q(X) =  0.0\n";
    cout << "  F(X) =  2*A*(X-0.5)/(A**2+(X-0.5)**2)**2\n";
    cout << "  IBC  =  3\n";
    cout << "  UL   =  ATAN(-0.5/A)\n";
    cout << "  UR   =  ATAN( 0.5/A)\n";
    cout << "  A    = " << alpha << "\n";
  }
//
//  Start out with just 4 subintervals.
//
  n = 4;
//
//  Initialize values that define the problem.
//
  init ( &ibc, n, &tol, &ul, &ur, &xl, xn, &xr );
//
//  Start the iteration counter off at 0.
//
  kount = 0;
//
//  Begin the next iteration.
//
  for ( ; ; )
  {
    kount = kount + 1;

    cout << "\n";
    cout << "  Begin new iteration with " << n << " nodes.\n";
    cout << "\n";
//
//  Solve the regular problem.
//
    solvex ( adiag, aleft, arite, f, h, ibc, indx, kount, n, NL, NMAX, 
      node, NQUAD, &nu, ul, ur, wquad, xn, xquad );
//
//  Solve N subproblems to get the error estimators.
//
    solvey ( eta, f, h, n, nu, ul, ur, xn );
//
//  Examine the error estimators, and see how many intervals should
//  be subdivided.
//
    status = subdiv ( eta, kount, &n, NMAX, tol, xn );

    if ( status != 0 )
    {
      break;
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM1D_ADAPTIVE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
# undef NL
# undef NMAX
# undef NQUAD
}
//****************************************************************************80

void assemble ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int n, int indx[], int node[], int nu, int nl, int nquad, 
  int nmax, double ul, double ur, double wquad[], double xn[], double xquad[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE assembles the global matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double ADIAG(NU).
//    ADIAG(I) is the "diagonal" coefficient of the I-th
//    equation in the linear system.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    Output, double ALEFT(NU).
//    ALEFT(I) is the "left hand" coefficient of the I-th
//    equation in the linear system.  That is, ALEFT(I) is the
//    coefficient of the (I-1)-th unknown in the I-th equation.
//    There is no value in ALEFT(1), since the first equation
//    does not refer to a "0-th" unknown.
//
//    Output, double ARITE(NU).
//    ARITE(I) is the "right hand" coefficient of the I-th
//    equation in the linear system.  ARITE(I) is the coefficient
//    of the (I+1)-th unknown in the I-th equation.  There is
//    no value in ARITE(NU) because the NU-th equation does not
//    refer to an "NU+1"-th unknown.
//
//    Output, double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    Input, double H(N)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    Input, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NQUAD
//    The number of quadrature points used in a subinterval.
//
//    Input, int NMAX, the maximum number of unknowns that can be handled.
//
//    Input, double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    Input, double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    Input, double WQUAD(NQUAD).
//    WQUAD(I) is the weight associated with the I-th point
//    of an NQUAD point Gaussian quadrature rule.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Input, double XQUAD(NQUAD,NMAX), the I-th quadrature point
//    in interval J.
//
{
  double aij;
  double he;
  int i;
  int ie;
  int ig;
  int il;
  int iq;
  int iu;
  int jg;
  int jl;
  int ju;
  double phii;
  double phiix;
  double phij;
  double phijx;
  double wquade;
  double x;
  double xleft;
  double xquade;
  double xrite;
//
//  Zero out the entries.
//
  for ( i = 0; i < nu; i++ )
  {
    f[i] = 0.0;
  }
  for ( i = 0; i < nu; i++ )
  {
    aleft[i] = 0.0;
  }
  for ( i = 0; i < nu; i++ )
  {
    arite[i] = 0.0;
  }
  for ( i = 0; i < nu; i++ )
  {
    adiag[i] = 0.0;
  }
//
//  For each interval,
//
  for ( ie = 0; ie < n; ie++ )
  {
    he = h[ie];
    xleft = xn[node[0+ie*2]];
    xrite = xn[node[1+ie*2]];
//
//  For each quadrature point in the interval,
//
    for ( iq = 0; iq < nquad; iq++ )
    {
      xquade = xquad[iq+ie*nquad];
      wquade = wquad[iq];
//
//  Pick a basis function which defines the equation,
//
      for ( il = 1; il <= nl; il++ )
      {
        ig = node[il-1+ie*nl];
        iu = indx[ig];

        if ( 0 < iu )
        {
          phi ( il, xquade, &phii, &phiix, xleft, xrite );
          f[iu-1] = f[iu-1] + he * wquade * ff ( xquade ) * phii;
//
//  Take care of boundary conditions specifying the value of U'.
//
          if ( ig == 0 )
          {
            x = 0.0;
            f[iu-1] = f[iu-1] - pp ( x ) * ul;
          }
          else if ( ig == n )
          {
            x = 1.0;
            f[iu-1] = f[iu-1] + pp ( x ) * ur;
          }
//
//  Pick a basis function which defines the coefficient
//  being computed.
//
          for ( jl = 1; jl <= nl; jl++ )
          {
            jg = node[jl-1+ie*nl];
            ju = indx[jg];
            phi ( jl, xquade, &phij, &phijx, xleft, xrite );

            aij = he * wquade * 
               ( pp ( xquade ) * phiix * phijx 
               + qq ( xquade ) * phii * phij );
//
//  Decide where the coefficient is to be added.
//
            if ( ju <= 0 )
            { 
              if ( jg == 0 )
              {
                f[iu-1] = f[iu-1] - aij * ul;
              }
              else if ( jg == n )
              {
                f[iu-1] = f[iu-1] - aij * ur;
              }
            }
            else if ( iu == ju )
            {
              adiag[iu-1] = adiag[iu-1] + aij;
            }
            else if ( ju < iu )
            {
              aleft[iu-1] = aleft[iu-1] + aij;
            }
            else
            {
              arite[iu-1] = arite[iu-1] + aij;
            }
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

double ff ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FF evaluates the function F in the differential equation.
//
//  Discussion:
//
//    This is the function F(X) that appears on the right hand
//    side of the equation:
//
//      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double FF, the value of F(X).
//
{
  double alpha;
  double beta;
  int problem;
  double pi = 3.141592653589793;
  double value;
//
//  Find out which problem we're working on.
//
  problem = get_problem ( );

  if ( problem == 1 )
  {
    value = 0.0;
  }
  else if ( problem == 2 )
  {
    value = -2.0 * x;
  }
  else if ( problem == 3 )
  {
    value = 0.25 * pi * pi * sin ( 0.5 * pi * x );
  }
  else if ( problem == 4 )
  {
    value = 0.25 * pi * pi * cos ( 0.5 * pi * x );
  }
  else if ( problem == 5 )
  {
    beta = get_beta ( );

    value = - pow ( x, beta ) + ( pow ( x, beta + 2.0 ) ) 
      / ( ( beta + 2.0 ) * ( beta + 1.0 ) );
  }
  else if ( problem == 6 )
  {
    alpha = get_alpha ( );

    value = 2.0 * alpha * ( x - 0.5 ) 
      / pow ( alpha * alpha + ( x - 0.5 ) * ( x - 0.5 ), 2 );
  }

  return value;
}
//****************************************************************************80

void geometry ( double h[], int ibc, int indx[], int n, int nl, int nmax, 
  int node[], int nquad, int *nu, double wquad[], double xn[], double xquad[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRY sets up some of the geometric information for the problem.  
//
//  Discussion:
//
//    Note, however, that the location of the nodes
//    is done outside of this routine, and, in fact, before this
//    routine is called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double H(N)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    Input, int IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    Output, int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NMAX, the maximum number of unknowns that can be handled.
//
//    Output, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NQUAD
//    The number of quadrature points used in a subinterval.
//
//    Output, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Output, double WQUAD(NQUAD).
//    WQUAD(I) is the weight associated with the I-th point
//    of an NQUAD point Gaussian quadrature rule.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Output, double XQUAD(NQUAD,NMAX), the I-th quadrature point
//    in interval J.
//
{
  double alfa;
  int i;
  int igl;
  int igr;
  double xl;
  double xr;
//
//  Store in NODE the fact that interval I has node I-1
//  as its left endpoint, and node I as its right endpoint.
//
  for ( i = 0; i < n; i++ )
  {
    node[0+i*2] = i;
    node[1+i*2] = i + 1;
  }
//
//  For every node that is associated with an unknown, we
//  record the number of the unknown in INDX.
//
  *nu = 0;
  for ( i = 0; i <= n; i++ )
  {
    if ( i == 0 && ( ibc == 1 || ibc == 3 ) )
    {
      indx[i] = -1;
    }
    else if ( i == n && ( ibc == 2 || ibc == 3 ) )
    {
      indx[i] = -1;
    }
    else
    {
      *nu = *nu + 1;
      indx[i] = *nu;
    }
  }
//
//  We compute the width of each interval.
//
  for ( i = 0; i < n; i++ )
  {
    igl = node[0+i*2];
    igr = node[1+i*2];
    h[i] = xn[igr] - xn[igl];
  }
//
//  We compute the location of the quadrature points in each
//  interval.
//
  for ( i = 0; i < n; i++ )
  {
    xl = xn[node[0+i*2]];
    xr = xn[node[1+i*2]];

    if ( nquad == 1 )
    {
      xquad[0+i*nquad] = 0.5 * ( xl + xr );
    }
    else if ( nquad == 2 )
    {
      alfa = -0.577350;
      xquad[0+i*nquad] = ( ( 1.0 - alfa ) * xl   
                       + ( 1.0 + alfa ) * xr ) 
                       /   2.0;
      alfa = +0.577350;
      xquad[1+i*nquad] = ( ( 1.0 - alfa ) * xl   
                       + ( 1.0 + alfa ) * xr ) 
                       /   2.0;
    }
    else if ( nquad == 3 )
    {
      alfa = -0.774597;
      xquad[0+i*nquad] = ( ( 1.0 - alfa ) * xl   
                       + ( 1.0 + alfa ) * xr ) 
                       /   2.0;
      xquad[1+i*nquad] = 0.5 * ( xl + xr );
      alfa = +0.774597;
      xquad[2+i*nquad] = ( ( 1.0 - alfa ) * xl   
                       + ( 1.0 + alfa ) * xr ) 
                       /   2.0;
    }
  }
//
//  Store the weights for the quadrature rule.
//
  if ( nquad == 1 )
  {
    wquad[0] = 1.0;
  }
  else if ( nquad == 2 )
  {
    wquad[0] = 0.5;
    wquad[1] = 0.5;
  }
  else if ( nquad == 3 )
  {
    wquad[0] = 4.0 / 9.0;
    wquad[1] = 5.0 / 18.0;
    wquad[2] = 4.0 / 9.0;
  }

  return;
}
//****************************************************************************80

double get_alpha ( void )

//****************************************************************************80
//
//  Purpose:
//
//    GET_ALPHA returns the value of ALPHA, for use by problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double GET_ALPHA, the value of ALPHA.
//
{
  double value;

  value = 0.01;

  return value;
}
//****************************************************************************80

double get_beta ( void )

//****************************************************************************80
//
//  Purpose:
//
//    GET_BETA returns the value of BETA, for use by problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double VALUE, the value of BETA.
//
{
  double value;

  value = -0.9;

  return value;
}
//****************************************************************************80

int get_problem ( void )

//****************************************************************************80
//
//  Purpose:
//
//    GETPRB returns the value of the current problem number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, int GET_PROBLEM, the index of the problem.
//
{
  int value;

  value = 6;

  return value;
}
//****************************************************************************80

void init ( int *ibc, int n, double *tol, double *ul, double *ur, double *xl, 
  double xn[], double *xr )

//****************************************************************************80
//
//  Purpose:
//
//    INIT initializes some parameters that define the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, int *IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Output, double *TOL.
//    A tolerance that is used to determine whether the estimated
//    error in an interval is so large that it should be subdivided
//    and the problem solved again.
//
//    Output, double *UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    Output, double *UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    Output, double *XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    Output, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Output, double *XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
{
  double alpha;
  double beta;
  int i;
  int problem;

  *tol = 0.01;
//
//  Find out which problem we're working on.
//
  problem = get_problem ( );
//
//  Set the boundary conditions for the problem, and
//  print out its title.
//
  if ( problem == 1 )
  {
    *ibc = 3;
    *ul = 0.0;
    *ur = 1.0;
    *xl = 0.0;
    *xr = 1.0;
    cout << "\n";
    cout << "Exact solution is U = X\n";
  }
  else if ( problem == 2 )
  {
    *ibc = 3;
    *ul = 0.0;
    *ur = 1.0;
    *xl = 0.0;
    *xr = 1.0;
    cout << "\n";
    cout << "Exact solution is U = X*X\n";
  }
  else if ( problem == 3 )
  {
    *ibc = 3;
    *ul = 0.0;
    *ur = 1.0;
    *xl = 0.0;
    *xr = 1.0;
    cout << "\n";
    cout << "Exact solution is U = SIN(PI*X/2)\n";
  }
  else if ( problem == 4 )
  {
    *ibc = 3;
    *ul = 1.0;
    *ur = 0.0;
    *xl = 0.0;
    *xr = 1.0;
    cout << "\n";
    cout << "Exact solution is U = COS(PI*X/2)\n";
  }
  else if ( problem == 5 )
  {
    *ibc = 3;
    beta = get_beta ( );
    *ul = 0.0;
    *ur = 1.0 / ( ( beta + 2.0 ) * ( beta + 1.0 ) );
    *xl = 0.0;
    *xr = 1.0;
    cout << "\n";
    cout << "Rheinboldt problem\n";
  }
  else if ( problem == 6 )
  {
    *ibc = 3;
    alpha = get_alpha ( );
    *xl = 0.0;
    *xr = 1.0;
    *ul = uexact ( *xl );
    *ur = uexact ( *xr );
    cout << "\n";
    cout << "Arctangent problem\n";
  }
//
//  The nodes are defined here, and not in the geometry routine.
//  This is because each new iteration chooses the location
//  of the new nodes in a special way.
//
  for ( i = 0; i <= n; i++ )
  {
    xn[i] = ( ( double ) ( n - i ) * ( *xl )   
            + ( double ) (     i ) * ( *xr ) )
            / ( double ) ( n     );
  }

  cout << "The equation is to be solved for \n";
  cout << "X greater than " << *xl << "\n";
  cout << " and less than " << *xr << "\n";
  cout << "\n";
  cout << "The boundary conditions are:\n";
  cout << "\n";

  if ( *ibc == 1 || *ibc == 3 )
  {
    cout << "  At X = XL, U = " << *ul << "\n";
  }
  else
  {
    cout << "  At X = XL, U' = " << *ul << "\n";
  }

  if ( *ibc == 2 || *ibc == 3 )
  {
    cout << "  At X = XR, U= " << *ur << "\n";
  }
  else
  {
    cout << "  At X = XR, U' = " << *ur << "\n";
  }

  return;
}
//****************************************************************************80

void output ( double f[], int ibc, int indx[], int n, int nu, double ul, 
  double ur, double xn[] )

//****************************************************************************80
//
//  Purpose:
//
//    OUTPUT prints out the computed solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    Input, int IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    Input, int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    Input, double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
{
  double error;
  int i;
  double u;
  double uex;

  cout << "\n";
  cout << "Node    X(I)        U(X(I))        Uexact        Error\n";
  cout << "\n";

  for ( i = 0; i <= n; i++ )
  {
    if ( i == 0 )
    {
      if ( ibc == 1 || ibc == 3 )
      {
        u = ul;
      }
      else
      {
        u = f[indx[i]-1];
      }
    }
    else if ( i == n )
    {
      if ( ibc == 2 || ibc == 3 )
      {
        u = ur;
      }
      else
      {
        u = f[indx[i]-1];
      }
    }
    else
    {
      u = f[indx[i]-1];
    }

    uex = uexact ( xn[i] );
    error = u - uex;

    cout << "  " << setw(4) << i
         << "  " << setw(12) << xn[i]
         << "  " << setw(12) << u
         << "  " << setw(12) << uex
         << "  " << setw(12) << error << "\n";
  }
  return;
}
//****************************************************************************80

void phi ( int il, double x, double *phii, double *phiix, double xleft, 
  double xrite )

//****************************************************************************80
//
//  Purpose:
//
//    PHI evaluates a linear basis function and its derivative.
//
//  Discussion:
//
//    The functions are evaluated at a point X in an interval.  In any
//    interval, there are just two basis functions.  The first
//    basis function is a line which is 1 at the left endpoint
//    and 0 at the right.  The second basis function is 0 at
//    the left endpoint and 1 at the right.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int IL, the local index of the basis function.
//
//    Input, double X, the evaluation point.
//
//    Output, double *PHII, *PHIIX, the value of the basis function
//    and its derivative.
//
//    Input, double XLEFT, XRITE, the endpoints of the interval.
//
{
  if ( xleft <= x && x <= xrite )
  {
    if ( il == 1 )
    {
      *phii = ( xrite - x ) / ( xrite - xleft );
      *phiix = -1.0 / ( xrite - xleft );
    }
    else
    {
      *phii = ( x - xleft ) / ( xrite - xleft );
      *phiix = 1.0 / ( xrite - xleft );
    }
  }
//
//  If X is outside of the interval, then the basis function
//  is always zero.
//
  else
  {
    *phii = 0.0;
    *phiix = 0.0;
  }
  return;
}
//****************************************************************************80

double pp ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    PP evaluates the function P in the differential equation.
//
//  Discussion:
//
//    The function P(X) occurs in the differential equation:
//
//      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double PP, the value of P(X).
//
{
  int problem;
  double value;
//
//  Find out which problem we're working on.
//
  problem = get_problem ( );

  if ( problem == 1 )
  {
    value = 1.0;
  }
  else if ( problem == 2 )
  {
    value = 1.0;
  }
  else if ( problem == 3 )
  {
    value = 1.0;
  }
  else if ( problem == 4 )
  {
    value = 1.0;
  }
  else if ( problem == 5 )
  {
    value = 1.0;
  }
  else if ( problem == 6 )
  {
    value = 1.0;
  }

  return value;
}
//****************************************************************************80

void prsys ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )

//****************************************************************************80
//
//  Purpose:
//
//    PRSYS prints out the tridiagonal linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double ADIAG(NU).
//    ADIAG(I) is the "diagonal" coefficient of the I-th
//    equation in the linear system.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    Input, double ALEFT(NU).
//    ALEFT(I) is the "left hand" coefficient of the I-th
//    equation in the linear system.  That is, ALEFT(I) is the
//    coefficient of the (I-1)-th unknown in the I-th equation.
//    There is no value in ALEFT(1), since the first equation
//    does not refer to a "0-th" unknown.
//
//    Input, double ARITE(NU).
//    ARITE(I) is the "right hand" coefficient of the I-th
//    equation in the linear system.  ARITE(I) is the coefficient
//    of the (I+1)-th unknown in the I-th equation.  There is
//    no value in ARITE(NU) because the NU-th equation does not
//    refer to an "NU+1"-th unknown.
//
//    Input, double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
{
  int i;

  cout << "\n";
  cout << "Printout of tridiagonal linear system:\n";
  cout << "\n";
  cout << "Equation  A-Left  A-Diag  A-Rite  RHS\n";
  cout << "\n";

  for ( i = 0; i < nu; i++ )
  {
    cout << setw(4) << i+1;
    if ( i == 0 )
    {
      cout << "              ";
    }
    else
    {
      cout << "  " << setw(12) << aleft[i];
    }
    cout << "  " << setw(12) << adiag[i];
    if ( i < nu - 1 )
    {
      cout << "  " << setw(12) << arite[i];
    }
    else
    {
      cout << "              ";
    }
    cout << "  " << setw(12) << f[i] << "\n";
  }

  return;
}
//****************************************************************************80

double qq ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    QQ evaluates the function Q in the differential equation.
//
//  Discussion:
//
//    The function Q(X) occurs in the differential equation:
//
//      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double QQ, the value of Q(X).
//
{
  int problem;
  double value;
//
//  Find out which problem we're working on.
//
  problem = get_problem ( );

  if ( problem == 1 )
  {
    value = 0.0;
  }
  else if ( problem == 2 )
  {
    value = 0.0;
  }
  else if ( problem == 3 )
  {
    value = 0.0;
  }
  else if ( problem == 4 )
  {
    value = 0.0;
  }
  else if ( problem == 5 )
  {
    value = 1.0;
  }
  else if ( problem == 6 )
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

void solve ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )

//****************************************************************************80
//
//  Purpose:
//
//    SOLVE solves a tridiagonal matrix system of the form A*x = b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input/output, double ADIAG(NU), ALEFT(NU), ARITE(NU).
//    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
//    left and right entries of the equations.
//    On output, ADIAG and ARITE have been changed in order
//    to compute the solution.
//    Note that for the first equation, there is no ALEFT
//    coefficient, and for the last, there is no ARITE.
//    So there is no need to store a value in ALEFT(1), nor
//    in ARITE(NU).
//
//    Input/output, double F(NU).
//    On input, F contains the right hand side of the linear
//    system to be solve
//    On output, F contains the solution of the linear system.
//
//    Input, int NU, the number of equations to be solved.
//
{
  int i;
//
//  Handle the special case of a single equation.
//
  if ( nu == 1 )
  {
    f[0] = f[0] / adiag[0];
  }
//
//  The general case, when NU is greater than 1.
//
  else
  {
    arite[0] = arite[0] / adiag[0];
    for ( i = 2; i <= nu - 1; i++ )
    {
      adiag[i-1] = adiag[i-1] - aleft[i-1] * arite[i-2];
      arite[i-1] = arite[i-1] / adiag[i-1];
    }
    adiag[nu-1] = adiag[nu-1] - aleft[nu-1] * arite[nu-2];

    f[0] = f[0] / adiag[0];
    for ( i = 2; i <= nu; i++ )
    {
      f[i-1] = ( f[i-1] - aleft[i-1] * f[i-2] ) / adiag[i-1];
    }

    for ( i = nu-1; 1 <= i; i-- )
    {
      f[i-1] = f[i-1] - arite[i-1] * f[i];
    }
  }

  return;
}
//****************************************************************************80

void solvex ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int ibc, int indx[], int kount, int n, int nl, int nmax, 
  int node[], int nquad, int *nu, double ul, double ur, double wquad[], 
  double xn[], double xquad[] )

//****************************************************************************80
//
//  Purpose:
//
//    SOLVEX discretizes and solves a differential equation given the nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Workspace, double ADIAG(NU).
//    ADIAG(I) is the "diagonal" coefficient of the I-th
//    equation in the linear system.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    Workspace, double ALEFT(NU).
//    ALEFT(I) is the "left hand" coefficient of the I-th
//    equation in the linear system.  That is, ALEFT(I) is the
//    coefficient of the (I-1)-th unknown in the I-th equation.
//    There is no value in ALEFT(1), since the first equation
//    does not refer to a "0-th" unknown.
//
//    Workspace, double ARITE(NU).
//    ARITE(I) is the "right hand" coefficient of the I-th
//    equation in the linear system.  ARITE(I) is the coefficient
//    of the (I+1)-th unknown in the I-th equation.  There is
//    no value in ARITE(NU) because the NU-th equation does not
//    refer to an "NU+1"-th unknown.
//
//    Output, double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    Output, double H(N)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    Input, int IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    Workspace, int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    Input, int KOUNT, the number of adaptive steps that have been taken.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NMAX, the maximum number of unknowns that can be handled.
//
//    Workspace, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Workspace, int NQUAD
//    The number of quadrature points used in a subinterval.
//
//    Output, int *NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    Input, double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    Workspace, double WQUAD(NQUAD).
//    WQUAD(I) is the weight associated with the I-th point
//    of an NQUAD point Gaussian quadrature rule.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Workspace, double XQUAD(NQUAD,NMAX), the I-th quadrature point
//    in interval J.
//
{
//
//  Given a set of N nodes (where N increases on each iteration),
//  compute the other geometric information.
//
  geometry ( h, ibc, indx, n, nl, nmax, node, nquad, nu, wquad, xn, xquad );
//
//  Assemble the linear system.
//
  assemble ( adiag, aleft, arite, f, h, n, indx, node, *nu, nl, 
    nquad, nmax, ul, ur, wquad, xn, xquad );
//
//  Print out the linear system, just once.
//
  if ( kount == 1 )
  {
    prsys ( adiag, aleft, arite, f, *nu );
  }
//
//  Solve the linear system.
//
  solve ( adiag, aleft, arite, f, *nu );
//
//  Print out the solution.
//
  cout << "\n";
  cout << "Basic solution\n";

  output ( f, ibc, indx, n, *nu, ul, ur, xn );

  return;
}
//****************************************************************************80

void solvey ( double eta[], double f[], double h[], int n, int nu, double ul, 
  double ur, double xn[] )

//****************************************************************************80
//
//  Purpose:
//
//    SOLVEY computes error estimators for a finite element solution.
//
//  Discussion:
//
//    SOLVEY accepts information about the solution of a finite element
//    problem on a grid of nodes with coordinates XN.  It then starts
//    at node 0, and for each node, computes two "error estimators",
//    one for the left, and one for the right interval associated with the
//    node.  These estimators are found by solving a finite element problem
//    over the two intervals, using the known values of the original
//    solution as boundary data, and using a mesh that is "slightly"
//    refined over the original one.
//
//    Note that the computations at the 0-th and N-th nodes only involve
//    a single interval.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double ETA(N).
//    ETA(I) is the error estimate for interval I.  It is computed
//    as the sum of two quantities, one associated with the left
//    and one with the right node of the interval.
//
//    Input, double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    Input, double H(N)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    Input, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    Input, double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
{
# define NL 2
# define NY 2
# define NQUAD 2
# define NMAY ( 2 * NY )

  double adiag[NMAY];
  double aleft[NMAY];
  double arite[NMAY];
  double fy[NMAY];
  double hy[NMAY];
  int i;
  int ibcy;
  int indy[NMAY+1];
  int j;
  int jhi;
  int jlo;
  int jmid;
  int k;
  int m;
  int nodey[NL*NMAY];
  int nuy;
  double total;
  double uleft;
  double ulval;
  double uly;
  double uprime;
  double urite;
  double urval;
  double ury;
  double uval;
  double vlval;
  double vprime;
  double vrval;
  double vval;
  double wquad[NQUAD];
  double xquady[NQUAD*NMAY];
  double y;
  double yl;
  double ym;
  double yn[NMAY+1];
  double yr;
//
//  Initialize the error estimators to zero.
//
  for ( j = 0; j < n; j++ )
  {
    eta[j] = 0.0;
  }
//
//  Set the boundary conditions for each subproblem to be
//  known values of U at the left and right.
//
//
//  For each node, subdivide its left and right hand intervals
//  into NY subintervals.
//
//  Set up and solve the differential equation again on this
//  smaller region.
//
//  The 0-th and N-th nodes are special cases.
//
  ibcy = 3;

  for ( j = 0; j <= n; j++ )
  {
    if ( j == 0 )
    {
      m = NY;
      jlo = j;
      jmid = j + 1;
      jhi = j + 1;
    }
    else if ( j == n )
    {
      m = NY;
      jlo = j - 1;
      jmid = j;
      jhi = j;
    }
    else
    {
      m = 2 * NY;
      jlo = j - 1;
      jmid = j;
      jhi = j + 1;
    }
//
//  Set the location of the nodes in the subintervals.
//
    yl = xn[jlo];
    ym = xn[jmid];
    yr = xn[jhi];

    for ( i = 0; i <= NY; i++ )
    {
      yn[i] = ( ( double ) ( NY - i ) * yl   
              + ( double ) (      i ) * ym ) 
              / ( double ) ( NY     );
    }

    for ( i = NY+1; i <= m; i++ )
    {
      yn[i] = ( ( double ) ( m - i      ) * ym   
              + ( double ) (     i - NY ) * yr ) 
              / ( double ) ( m -     NY );
    }
//
//  Set up the geometry of the sub-problem.
//
    geometry ( hy, ibcy, indy, m, NL, NMAY, nodey, NQUAD, &nuy, 
      wquad, yn, xquady );
//
//  Set the boundary values for the sub-problem.
//
    if ( j <= 1 )
    {
      uly = ul;
    }
    else
    {
      uly = f[j-2];
    }

    if ( n - 1 <= j )
    {
      ury = ur;
    }
    else
    {
      ury = f[j];
    }
//
//  Assemble the matrix for the sub-problem.
//
    assemble ( adiag, aleft, arite, fy, hy, m, indy, nodey, nuy, NL, 
      NQUAD, NMAY, uly, ury, wquad, yn, xquady );
//
//  Solve the system.
//
    solve ( adiag, aleft, arite, fy, nuy );
//
//  Compute the weighted sum of the squares of the differences
//  of the original computed slope and the refined computed slopes.
//
//  Calculation for left interval.
//
    if ( 1 <= j )
    {
      if ( j <= 1 )
      {
        uleft = ul;
        urite = f[0];
      }
      else if ( j == n )
      {
        uleft = f[j-2];
        urite = ur;
      }
      else
      {
        uleft = f[j-2];
        urite = f[j-1];
      }

      uprime = ( urite - uleft ) / h[j-1];

      total = 0.0;

      for ( i = 1; i <= NY; i++ )
      {
        yl = yn[i-1];
        yr = yn[i];

        if ( i == 1 )
        {
          vlval = uly;
          vrval = fy[i-1];
        }
        else if ( i == m )
        {
          vlval = fy[i-2];
          vrval = ury;
        }
        else
        {
          vlval = fy[i-2];
          vrval = fy[i-1];
        }

        vprime = ( vrval - vlval ) / hy[i-1];

        ulval = ( ( double ) ( NY - i + 1 ) * uleft   
                + ( double ) (      i - 1 ) * urite ) 
                / ( double ) ( NY         );

        urval = ( ( double ) ( NY - i ) * uleft   
                + ( double ) (      i ) * urite ) 
                / ( double ) ( NY     );
//
//  Compute the integral of
//
//    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
//
        for ( k = 0; k < NQUAD; k++ )
        {
          y  =  xquady[k+(i-1)*NQUAD];

          uval = ( ( yl - y      ) * urval   
                 + (      y - yr ) * ulval ) 
                 / ( yl     - yr );

          vval = ( ( yl - y      ) * vrval   
                 + (      y - yr ) * vlval ) 
                 / ( yl     - yr );

          total = total + 0.5 * wquad[k] * hy[i-1] * 
            ( pp ( y ) * pow ( uprime - vprime, 2 ) 
            + qq ( y ) * pow ( uval - vval, 2 ) );
        }
      }
      eta[j-1] = eta[j-1] + 0.5 * sqrt ( total );
    }
//
//  Calculation for right interval.
//
    if ( j <= n - 1 )
    {
      if ( j == 0 )
      {
        uleft = ul;
        urite = f[j];
      }
      else if ( n - 1 <= j )
      {
        uleft = f[j-1];
        urite = ur;
      }
      else
      {
        uleft = f[j-1];
        urite = f[j];
      }

      uprime = ( urite - uleft ) / h[j];

      total = 0.0;

      for ( i = m+1-NY; i <= m; i++ )
      {
        yl = yn[i-1];
        yr = yn[i];

        if ( i == 1 )
        {
          vlval = uly;
          vrval = fy[i-1];
        }
        else if ( i == m )
        {
          vlval = fy[i-2];
          vrval = ury;
        }
        else
        {
          vlval = fy[i-2];
          vrval = fy[i-1];
        }

        vprime = ( vrval - vlval ) / hy[i-1];

        ulval = ( ( double ) (      m - i + 1 ) * uleft   
                + ( double ) ( NY - m + i - 1 ) * urite ) 
                / ( double ) ( NY             );

        urval = ( ( double ) (      m - i ) * uleft   
                + ( double ) ( NY - m + i ) * urite ) 
                / ( double ) ( NY         );
//
//  Compute the integral of
//
//    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
//
        for ( k = 0; k < NQUAD; k++ )
        {
          y  =  xquady[k+(i-1)*NQUAD];

          uval = ( ( yl - y      ) * urval   
                 + (      y - yr ) * ulval ) 
                 / ( yl     - yr );

          vval = ( ( yl - y      ) * vrval   
                 + (      y - yr ) * vlval ) 
                 / ( yl     - yr );
 
          total = total + 0.5 * wquad[k] * hy[i-1] * 
            ( pp ( y ) * pow ( uprime - vprime, 2 )
            + qq ( y ) * pow ( uval - vval, 2 ) );
        }
      }

      eta[j] = eta[j] + 0.5 * sqrt ( total );
    }
  }
//
//  Print out the error estimators.
//
  cout << "\n";
  cout << "ETA\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout << setw(12) << eta[j] << "\n";
  }
  return;
# undef NL
# undef NMAY
# undef NQUAD
# undef NY
}
//****************************************************************************80

int subdiv ( double eta[], int kount, int *n, int nmax, double tol, 
  double xn[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUBDIV decides which intervals should be subdivided.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double ETA(N).
//    ETA(I) is the error estimate for interval I.  It is computed
//    as the sum of two quantities, one associated with the left
//    and one with the right node of the interval.
//
//    Input, int KOUNT, the number of adaptive steps that have been taken.
//
//    Input/output, int N
//    The number of subintervals into which the interval
//    [XL,XR] is broken.  
//
//    Input, int NMAX, the maximum number of unknowns that can be handled.
//
//    Input, double TOL.
//    A tolerance that is used to determine whether the estimated
//    error in an interval is so large that it should be subdivided
//    and the problem solved again.
//
//    Input/output, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Output, int SUBDIV, reports status of subdivision.
//    0, a new subdivision was carried out.
//    1, no more subdivisions are needed.
//    -1, no more subdivisions can be carried out.
//
//  Local Parameters:
//
//    Local, int JADD(N).
//    JADD(I) is 1 if the error estimates show that interval I
//    should be subdivided.
//
{
  double ave;
  int j;
  int *jadd;
  int k;
  int status;
  double temp;
  double *xt;

  status = 0;
//
//  Add up the ETA's, and get their average.
//
  ave = 0.0;
  for ( j = 0; j < *n; j++ )
  {
    ave = ave + eta[j];
  }
  ave = ave / ( double ) ( *n );
//
//  Look for intervals whose ETA value is relatively large,
//  and note in JADD that these intervals should be subdivided.
//
  k = 0;
  temp = r8_max ( 1.2 * ave + 0.00001, tol * tol / ( double ) ( *n ) );

  cout << "\n";
  cout << "Tolerance = " << temp << "\n";
  cout << "\n";

  jadd = new int[*n];

  for ( j = 0; j < *n; j++ )
  {
    if ( temp < eta[j] )
    {
      k = k + 1;
      jadd[j] = 1;
      cout << "Subdivide interval " << j + 1 << "\n";
    }
    else
    {
      jadd[j] = 0;
    }
  }
//
//  If no subdivisions needed, we're done.
//
  if ( k <= 0 )
  {
    cout << "Success on step " << kount << "\n";
    status = 1;
    return status;
  }
//
//  See if we're about to go over our limit.
//
  if ( nmax < *n + k )
  {
    cout << "\n";
    cout << "The iterations did not reach their goal.\n";
    cout << "The next value of N is " << *n + k << "\n";
    cout << "which exceeds NMAX = " << nmax << "\n";
    status = -1;
    return status;
  }
//
//  Insert new nodes where needed.
//
  xt = new double[nmax+1];

  k = 0;
  xt[0] = xn[0];
  for ( j = 0; j < *n; j++ )
  {
    if ( 0 < jadd[j] )
    {
      xt[j+1+k] = 0.5 * ( xn[j+1] + xn[j] );
      k = k + 1;
    }
    xt[j+1+k] = xn[j+1];
  }
//
//  Update the value of N, and copy the new nodes into XN.
//
  *n = *n + k;

  for ( j = 0; j <= *n; j++ )
  {
    xn[j] = xt[j];
  }

  delete [] jadd;
  delete [] xt;

  return status;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double uexact ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    UEXACT returns the value of the exact solution at any point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double UEXACT, the value of the exact solution at X.
//
{
  double alpha;
  double beta;
  int problem;
  double pi = 3.141592653589793;
  double value;
//
//  Find out which problem we're working on.
//
  problem = get_problem ( );

  if ( problem == 1 )
  {
    value = x;
  }
  else if ( problem == 2 )
  {
    value = x * x;
  }
  else if ( problem == 3 )
  {
    value = sin ( pi * x / 2.0 );
  }
  else if ( problem == 4 )
  {
    value = cos ( pi * x / 2.0 );
  }
  else if ( problem == 5 )
  {
    beta = get_beta ( );

    value = ( pow ( x, beta + 2.0 ) ) / ( ( beta + 2.0 ) * ( beta + 1.0 ) );
  }
  else if ( problem == 6 )
  {
    alpha = get_alpha ( );
    value = atan ( ( x - 0.5 ) / alpha );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
