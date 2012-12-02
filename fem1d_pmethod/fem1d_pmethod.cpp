# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( void );
void alpbet ( double a[], double alpha[], double beta[], int np, int problem, 
  int quad_num, double quad_w[], double quad_x[] );
void exact ( double alpha[], double beta[], double f[], int np, int nprint, 
  int problem, int quad_num, double quad_w[], double quad_x[] );
double ff ( double x, int problem );
void ortho ( double a[], double alpha[], double beta[], int np, int problem, 
  int quad_num, double quad_w[], double quad_x[] );
void out ( double alpha[], double beta[], double f[], int np, int nprint );
void phi ( double alpha[], double beta[], int i, int np, double *phii, 
  double *phiix, double x );
double pp ( double x, int problem );
void quad ( int quad_num, double quad_w[], double quad_x[] );
double qq ( double x, int problem );
void sol ( double a[], double alpha[], double beta[], double f[], int np, 
  int problem, int quad_num, double quad_w[], double quad_x[] );
void timestamp ( void );
double uex ( double x, int problem );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_PMETHOD.
//
//  Discussion:
//
//    FEM1D_PMETHOD implements the P-version of the finite element method.
//
//    Program to solve the one dimensional problem:
//
//      - d/dX (P dU/dX) + Q U  =  F
//
//    by the finite-element method using a sequence of polynomials
//    which satisfy the boundary conditions and are orthogonal
//    with respect to the inner product:
//
//      (U,V)  =  Integral (-1 to 1) P U' V' + Q U V dx
//
//    Here U is an unknown scalar function of X defined on the
//    interval [-1,1], and P, Q and F are given functions of X.
//
//    The boundary values are U(-1) = U(1)=0.
//
//    Sample problem #1:
//
//      U=1-x**4,        P=1, Q=1, F=1.0+12.0*x**2-x**4
//
//    Sample problem #2:
//
//      U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x)
//
//    The program should be able to get the exact solution for
//    the first problem, using NP = 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    int NP.
//    The highest degree polynomial to use.
//
//    int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
//    int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    int QUAD_NUM, the order of the quadrature rule.
//
//    double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
# define NP 2
# define QUAD_NUM 10

  double a[NP+1];
  double alpha[NP];
  double beta[NP];
  double f[NP+1];
  int nprint = 10;
  int problem = 2;
  double quad_w[QUAD_NUM];
  double quad_x[QUAD_NUM];

  timestamp ( );

  cout << "\n";
  cout << "FEM1D_PMETHOD\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Solve the two-point boundary value problem\n";
  cout << "\n";
  cout << "  - d/dX (P dU/dX) + Q U  =  F\n";
  cout << "\n";
  cout << "  on the interval [-1,1], with\n";
  cout << "  U(-1) = U(1) = 0.\n";
  cout << "\n";
  cout << "  The P method is used, which represents U as\n";
  cout << "  a weighted sum of orthogonal polynomials.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Highest degree polynomial to use is " << NP << "\n";
  cout << "  Number of points to be used for output = " << nprint << "\n";

  if ( problem == 1 )
  {
    cout << "\n";
    cout << "  Problem #1:\n";
    cout << "  U=1-x**4,\n";
    cout << "  P=1,\n";
    cout << "  Q=1,\n";
    cout << "  F=1 + 12 * x**2 - x**4\n";
  }
  else if ( problem == 2 )
  {
    cout << "\n";
    cout << "  Problem #2:\n";
    cout << "  U=cos(0.5*pi*x),\n";
    cout << "  P=1,\n";
    cout << "  Q=0,\n";
    cout << "  F=0.25*pi*pi*cos(0.5*pi*x)\n";
  }
//
//  Get quadrature abscissas and weights for interval [-1,1].
//
  quad ( QUAD_NUM, quad_w, quad_x );
//
//  Compute the constants for the recurrence relationship
//  that defines the basis functions.
//
  alpbet ( a, alpha, beta, NP, problem, QUAD_NUM, quad_w, quad_x );
//
//  Test the orthogonality of the basis functions.
//
  ortho ( a, alpha, beta, NP, problem, QUAD_NUM, quad_w, quad_x );
//
//  Solve for the solution of the problem, in terms of coefficients
//  of the basis functions.
//
  sol ( a, alpha, beta, f, NP, problem, QUAD_NUM, quad_w, quad_x );
//
//  Print out the solution, evaluated at each of the NPRINT points.
//
  out ( alpha, beta, f, NP, nprint );
//
//  Compare the computed and exact solutions.
//
  exact ( alpha, beta, f, NP, nprint, problem, QUAD_NUM, quad_w, quad_x );

  cout << "\n";
  cout << "PMETHOD\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
# undef NP
# undef QUAD_NUM
}
//****************************************************************************80

void alpbet ( double a[], double alpha[], double beta[], int np, int problem, 
  int quad_num, double quad_w[], double quad_x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ALPBET calculates the coefficients in the recurrence relationship.
//
//  Discussion:
//
//    ALPHA and BETA are the coefficients in the three
//    term recurrence relation for the orthogonal basis functions
//    on [-1,1].
//
//    The routine also calculates A, the square of the norm of each basis
//    function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Output, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Output, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
  int i;
  int iq;
  int k;
  double q;
  double qm1;
  double qm1x;
  double qm2;
  double qm2x;
  double qx;
  double s;
  double ss;
  double su;
  double sv;
  double t;
  double u;
  double v;
  double x;

  ss = 0.0;
  su = 0.0;

  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];

    s = 4.0 * pp ( x, problem ) * x * x 
      + qq ( x, problem ) * ( 1.0 - x * x ) * ( 1.0 - x * x );

    u = 2.0 * pp ( x, problem ) * x * ( 3.0 * x * x - 1.0 ) 
      + x * qq ( x, problem ) * ( 1.0 - x * x ) * ( 1.0 - x * x );

    ss = ss + s * quad_w[iq];
    su = su + u * quad_w[iq];

  }

  a[0] = ss;
  alpha[0] = su / ss;
  beta[0] = 0.0;

  for ( i = 1; i <= np; i++ )
  {
    ss = 0.0;
    su = 0.0;
    sv = 0.0;

    for ( iq = 0; iq < quad_num; iq++ )
    {
      x = quad_x[iq];
      q = 1.0;
      qm1 = 0.0;
      qx = 0.0;
      qm1x = 0.0;

      for ( k = 0; k <= i-1; k++ )
      {
        qm2 = qm1;
        qm1 = q;
        qm2x = qm1x;
        qm1x = qx;
        q = ( x - alpha[k] ) * qm1 - beta[k] * qm2;
        qx = qm1 + ( x - alpha[k] ) * qm1x - beta[k] * qm2x;
      }

      t = 1.0 - x * x;

      s = pp ( x, problem ) * pow ( t * qx - 2.0 * x * q, 2 )
        + qq ( x, problem ) * t * t * q * q;

      u = pp ( x, problem ) 
        * ( x * t * qx + ( 1.0 - 3.0 * x * x ) * q ) 
        * ( t * qx - 2.0 * x * q ) + x * qq ( x, problem ) 
        * t * t * q * q;

      v = pp ( x, problem ) 
        * ( x * t * qx + ( 1.0 - 3.0 * x * x ) * q ) 
        * ( t * qm1x - 2.0 * x * qm1 ) 
        + x * qq ( x, problem ) * t * t * q * qm1;

      ss = ss + s * quad_w[iq];
      su = su + u * quad_w[iq];
      sv = sv + v * quad_w[iq];
    }

    a[i] = ss;

    if ( i < np )
    {
      alpha[i] = su / ss;
      beta[i] = sv / a[i-1];
    }
  }

  return;
}
//****************************************************************************80

void exact ( double alpha[], double beta[], double f[], int np, int nprint, 
  int problem, int quad_num, double quad_w[], double quad_x[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT compares the computed and exact solutions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
  double big_l2;
  double error;
  int i;
  int ip;
  int j;
  int k;
  int nsub = 10;
  double phii;
  double phiix;
  double ue;
  double up;
  double x;
  double xl;
  double xr;

  cout << "\n";
  cout << "Comparison of computed and exact solutions:\n";
  cout << "\n";
  cout << "    X        U computed    U exact     Difference\n";
  cout << "\n";

  for ( i = 0; i <= nprint; i++ )
  {
    x = ( double ) ( 2 * i - nprint ) / ( double ) ( nprint );
    ue = uex ( x, problem );
    up = 0.0;
    for ( j = 0; j <= np; j++ )
    {
      phi ( alpha, beta, j, np, &phii, &phiix, x );
      up = up + phii * f[j];
    }
    cout << "  " << setw(8) << x
         << "  " << setw(12) << up
         << "  " << setw(12) << ue
         << "  " << setw(12) << ue - up << "\n";
  }
//
//  Compute the big L2 error.
//
  big_l2 = 0.0;

  for ( i = 1; i <= nsub; i++ )
  {
    xl = ( double ) ( 2 * i - nsub - 1 ) / ( double ) ( nsub );
    xr = ( double ) ( 2 * i - nsub     ) / ( double ) ( nsub );

    for ( j = 0; j < quad_num; j++ )
    {
      x = ( xl * ( 1.0 - quad_x[j] ) 
          + xr * ( 1.0 + quad_x[j] ) ) / 2.0;

      up = 0.0;
      for ( k = 0; k <= np; k++ )
      {
        phi ( alpha, beta, k, np, &phii, &phiix, x );
        up = up + phii * f[k];
      }

      big_l2 = big_l2 + pow ( up - uex ( x, problem ), 2 ) * quad_w[j]
        * ( xr - xl ) / 2.0;
    }
  }

  big_l2 = sqrt ( big_l2 );

  cout << "\n";
  cout << "Big L2 error = " << big_l2 << "\n";

  return;
}
//****************************************************************************80

double ff ( double x, int problem )

//****************************************************************************80
//
//  Purpose:
//
//    FF evaluates the right hand side function F(X) at any point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double FF, the value of F(X).
//
{
  double pi = 3.141592653589793;
  double value;
//
//  Test problem 1
//
  if ( problem == 1 )
  {
    value = 1.0 + 12.0 * x * x - x * x * x * x;
  }
//
//  Test problem 2
//
  else if ( problem == 2 )
  {
    value = 0.25 * pi * pi * cos ( 0.5 * pi * x );
  }

  return value;
}
//****************************************************************************80

void ortho ( double a[], double alpha[], double beta[], int np, int problem, 
  int quad_num, double quad_w[], double quad_x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ORTHO tests the basis functions for orthogonality.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Input, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
  double *b;
  double bij;
  int i;
  int iq;
  int j;
  double phii;
  double phiix;
  double phij;
  double phijx;
  double x;
//
//  Zero out the B array, so we can start summing up the dot products.
//
  b = new double[(np+1)*(np+1)];

  for ( j = 0; j <= np; j++ )
  {
    for ( i = 0; i <= np; i++ )
    {
       b[i+j*(np+1)] = 0.0;
    }
  }
//
//  Approximate the integral of the product of basis function
//  I and basis function J over the interval [-1,1].
//
//  We expect to get zero, except when I and J are equal,
//  when we should get A(I).
//
  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      for ( j = 0; j <= np; j++ )
      {
        phi ( alpha, beta, j, np, &phij, &phijx, x );

        bij = pp ( x, problem ) * phiix * phijx 
            + qq ( x, problem ) * phii * phij;

        b[i+j*(np+1)] = b[i+j*(np+1)] + bij * quad_w[iq];
      }
    }
  }
//
//  Print out the results of the test.
//
  cout << "\n";
  cout << "Basis function orthogonality test:\n";
  cout << "\n";
  cout << "   i   j     b(i,j)/a(i)\n";
  cout << "\n";
  for ( i = 0; i <= np; i++ )
  {
    cout << "\n";
    for ( j = 0; j <= np; j++ )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(6) << j
           << "  " << setw(12) << b[i+j*(np+1)] / a[i] << "\n";
    }
  }

  delete [] b;

  return;
}
//****************************************************************************80

void out ( double alpha[], double beta[], double f[], int np, int nprint )

//****************************************************************************80
//
//  Purpose:
//
//    OUT prints the computed solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
{
  int i;
  int ip;
  double phii;
  double phiix;
  double up;
  double x;

  cout << "\n";
  cout << "Representation of solution:\n";
  cout << "\n";
  cout << "Basis function coefficients:\n";
  cout << "\n";
  for ( i = 0; i <= np; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(12) << f[i] << "\n";
  }

  cout << "\n";
  cout << "\n";
  cout << "       X     Approximate Solution\n";
  cout << "\n";
  for ( ip = 0; ip <= nprint; ip++ )
  {
    x = ( double ) ( 2 * ip - nprint ) / ( double ) ( nprint );
    up = 0.0;
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      up = up + phii * f[i];
    }
    cout << "  " << setw(12) << x
         << "  " << setw(12) << up << "\n";
  }

  cout << "\n";

  return;
}
//****************************************************************************80

void phi ( double alpha[], double beta[], int i, int np, double *phii, 
  double *phiix, double x )

//****************************************************************************80
//
//  Purpose:
//
//    PHI evaluates the I-th basis function at the point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, int I, the index of the basis function.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Output, double PHII, PHIIX, the value of the basis
//    function and its derivative.
//
//    Input, double X, the evaluation point.
//
{
  int j;
  double q;
  double qm1;
  double qm1x;
  double qm2;
  double qm2x;
  double qx;
  double t;

  qm1 = 0.0;
  q = 1.0;
  qm1x = 0.0;
  qx = 0.0;

  for ( j = 1; j <= i; j++ )
  {
    qm2 = qm1;
    qm1 = q;
    qm2x = qm1x;
    qm1x = qx;
    t = x - alpha[j-1];
    q = t * qm1 - beta[j-1] * qm2;
    qx = qm1 + t * qm1x - beta[j-1] * qm2x;
  }

  t = 1.0 - x * x;
  *phii = t * q;
  *phiix = t * qx - 2.0 * x * q;

  return;
}
//****************************************************************************80

double pp ( double x, int problem )

//****************************************************************************80
//
//  Purpose:
//
//    PP returns the value of the coefficient function P(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double PP, the value of P(X).
//
{
  double value;
//
//  Test problem 1
//
  if ( problem == 1 )
  {
    value = 1.0;
  }
//
//  Test problem 2
//
  else if ( problem == 2 )
  {
    value = 1.0;
  }

  return value;
}
//****************************************************************************80

void quad ( int quad_num, double quad_w[], double quad_x[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD returns the abscissas and weights for gaussian quadrature on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Output, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Output, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
//
//  Quadrature points on [-1,1]
//
  quad_x[1-1] = -0.973906528517172;
  quad_x[2-1] = -0.865063366688985;
  quad_x[3-1] = -0.679409568299024;
  quad_x[4-1] = -0.433395394129247;
  quad_x[5-1] = -0.148874338981631;
  quad_x[6-1] =  0.148874338981631;
  quad_x[7-1] =  0.433395394129247;
  quad_x[8-1] =  0.679409568299024;
  quad_x[9-1] =  0.865063366688985;
  quad_x[10-1] = 0.973906528517172;
//
//  Weight factors
//
  quad_w[1-1] =  0.066671344308688;
  quad_w[2-1] =  0.149451349150581;
  quad_w[3-1] =  0.219086362515982;
  quad_w[4-1] =  0.269266719309996;
  quad_w[5-1] =  0.295524224714753;
  quad_w[6-1] =  0.295524224714753;
  quad_w[7-1] =  0.269266719309996;
  quad_w[8-1] =  0.219086362515982;
  quad_w[9-1] =  0.149451349150581;
  quad_w[10-1] = 0.066671344308688;

  return;
}
//****************************************************************************80

double qq ( double x, int problem )

//****************************************************************************80
//
//  Purpose:
//
//    QQ returns the value of the coefficient function Q(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double QQ, the value of Q(X).
//
{
  double value;
//
//  Test problem 1
//
  if ( problem == 1 )
  {
    value = 1.0;
  }
//
//  Test problem 2
//
  else if ( problem == 2 )
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void sol ( double a[], double alpha[], double beta[], double f[], int np, 
  int problem, int quad_num, double quad_w[], double quad_x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SOL solves a linear system for the finite element coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Input, double ALPHA(NP).
//    ALPHA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Input, double BETA(NP).
//    BETA(I) contains one of the coefficients of a recurrence
//    relationship that defines the basis functions.
//
//    Output, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
{
  int i;
  int iq;
  double phii;
  double phiix;
  double t;
  double x;

  for ( i = 0; i <= np; i++ )
  {
    f[i] = 0.0;
  }

  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];
    t = ff ( x, problem ) * quad_w[iq];
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      f[i] = f[i] + phii * t;
    }
  }

  for ( i = 0; i <= np; i++ )
  {
    f[i] = f[i] / a[i];
  }

  return;
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

double uex ( double x, int problem )

//****************************************************************************80
//
//  Purpose:
//
//    UEX returns the value of the exact solution at a point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double UEX, the exact value of U(X).
//
{
  double pi = 3.141592653589793;
  double value;
//
//  Test problem 1
//
  if ( problem == 1 )
  {
    value = 1.0 - pow ( x, 4 );
  }
//
//  Test problem 2
//
  else if ( problem == 2 )
  {
    value = cos ( 0.5 * pi * x );
  }

  return value;
}
