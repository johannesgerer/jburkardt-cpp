# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "toms655.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( int nt, int kind, double alpha, double beta );
void test11 ( int nt, int kind, double alpha, double beta, double a, double b );
double f ( double x, int i );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS655_PRB.
//
//  Discussion:
//
//    TOMS655_PRB tests the TOMS655 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int kind;
  int nt;

  timestamp ( );
  cout << "\n";
  cout << "TOMS655_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS655 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
//
//  Compute 15 points of an example of each rule.
//
  for ( kind = 1; kind <= 9; kind++ )
  {
    nt = 15;
    if ( kind == 8 )
    {
      alpha = 1.0;
      beta = - alpha - 2 * nt - 2;
    }
    else
    {
      alpha = 0.0;
      beta = 0.0;
    }
    test10 ( nt, kind, alpha, beta );
  }
//
//  Compute 15 points of an example of each rule using nondefault A, B.
//
  for ( kind = 1; kind <= 9; kind++ )
  {
    nt = 15;

    if ( kind == 1 )
    {
      alpha = 0.0;
      beta = 0.0;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 2 )
    {
      alpha = 0.0;
      beta = 0.0;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 3 )
    {
      alpha = 1.0;
      beta = 0.0;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 4 )
    {
      alpha = 1.5;
      beta = 0.5;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 5 )
    {
      alpha = 1.0;
      beta = 0.0;
      a = 1.0;
      b = 1.0;
    }
    else if ( kind == 6 )
    {
      alpha = 1.0;
      beta = 0.0;
      a = 0.0;
      b = 0.5;
    }
    else if ( kind == 7 )
    {
      alpha = 1.0;
      beta = 0.0;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 8 )
    {
      alpha = 1.0;
      beta = - alpha - 2 * nt - 2;
      a = 0.0;
      b = 1.0;
    }
    else if ( kind == 9 )
    {
      alpha = 0.0;
      beta = 0.0;
      a = 0.0;
      b = 1.0;
    }

    test11 ( nt, kind, alpha, beta, a, b );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS655_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double alpha;
  double beta;
  int i;
  int key;
  int kind;
  int lu;
  int *mlt;
  int *ndx;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double *t;
  double *wts;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test CIQFS.\n";
//
//  Number of knots.
//
  nt = 5;
//
//  Set the knots in the default interval [-1,+1].
//
  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
//
//  Set the knot multiplicities.
//
  mlt = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
//
//  Set the size of the weights array.
//
  nwts = 0;
  for ( i = 0; i < nt; i++ )
  {
    nwts = nwts + mlt[i];
  }
//
//  Because KEY = 1, NDX will be set up for us.
//
  ndx = new int[nt];
//
//  KEY = 1 indicates that the WTS array should hold the weights
//  in the usual order.
//
  key = 1;
//
//  Request Legendre weight function.
//
  kind = 1;
//
//  ALPHA, BETA not used in Legendre weight function but set anyway.
//
  alpha = 0.0;
  beta  = 0.0;
//
//  LU controls printing.
//  A positive value requests that we compute and print weights, and
//  conduct a moments check.
//
  lu = 6;
//
//  This call returns the WTS array.
//
  wts = ciqfs ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, lu );

  delete [] mlt;
  delete [] ndx;
  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int key;
  int kind;
  int lu;
  int *mlt;
  int *ndx;
  int nt;
  int nwts;
  double *t;
  double *wts;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test CIQF, CIQFS, CGQF and CGQFS\n";
  cout << "  with all classical weight functions.\n";
//
//  Try all weight functions.
//
  for ( kind = 1; kind <= 9; kind++ )
  {
//
//  Number of knots.
//
    nt = 5;
//
//  Set parameters ALPHA and BETA.
//
    alpha = 0.5;
    if ( kind != 8 )
    {
      beta  = 2.0;
    }
    else
    {
      beta = - 16.0;
    }
//
//  Set A and B.
//
    a = - 0.5;
    b = 2.0;
//
//  Have CGQF compute the knots and weights.
//
    lu = 6;
    t = new double[nt];
    wts = new double[nt];

    cout << "\n";
    cout << "  Knots and weights of Gauss quadrature formula\n";
    cout << "  computed by CGQF.\n";
    cgqf ( nt, kind, alpha, beta, a, b, lu, t, wts );
//
//  Now compute the weights for the same knots by CIQF.
//
//  Set the knot multiplicities.
//
    mlt = new int[nt];
    for ( i = 0; i < nt; i++ )
    {
      mlt[i] = 2;
    }
//
//  Set the size of the weights array.
//
    nwts = 0;
    for ( i = 0; i < nt; i++ )
    {
      nwts = nwts + mlt[i];
    }
//
//  We need to deallocate and reallocate WTS because it is now of
//  dimension NWTS rather than NT.
//
    delete [] wts;
    wts = new double[nwts];
//
//  Because KEY = 1, NDX will be set up for us.
//
    ndx = new int[nt];
//
//  KEY = 1 indicates that the WTS array should hold the weights
//  in the usual order.
//
    key = 1;
//
//  LU controls printing.
//  A positive value requests that we compute and print weights, and
//  conduct a moments check.
//
    lu = 6;

    cout << "\n";
    cout << "  Weights of Gauss quadrature formula computed from the\n";
    cout << "  knots by CIQF.\n";

    wts = ciqf ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, a, b, lu );

    delete [] mlt;
    delete [] ndx;
    delete [] t;
    delete [] wts;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CEIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double alpha;
  double beta;
  int i;
  int kind;
  int lu;
  int *mlt;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test CEIQFS.\n";
//
//  Number of knots.
//
  nt = 5;
//
//  Set the knots in the default interval [-1,+1].
//
  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
//
//  Set the knot multiplicities.
//
  mlt = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
//
//  Set KIND to the Legendre weight function.
//
  kind = 1;
//
//  ALPHA, BETA not used in Legendre weight function but set anyway.
//
  alpha = 0.0;
  beta  = 0.0;
//
//  Call CEIQFS to set up the quadrature formula and evaluate it on F.
//
  qfsum = ceiqfs ( nt, t, mlt, kind, alpha, beta, f );

  cout << "\n";
  cout << "  Integral of sin(x) on -1, 1 by Fejer type rule\n";
  cout << "  with " << nt << " points of multiplicity 2.\n";
  cout << "  Quadrature formula:" << setw(24) << setprecision(16) << qfsum << "\n";

  qfsx = cos ( - 1.0 ) - cos ( 1.0 );
  cout << "  Exact value       :" << setw(24) << setprecision(16) << qfsx << "\n";
  cout << "  Error             :" << r8_abs ( qfsum - qfsx ) << "\n";

  delete [] mlt;
  delete [] t;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CEIQF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int kind;
  int lu;
  int *mlt;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test CEIQF.\n";
//
//  Number of knots.
//
  nt = 5;
//
//  Set the knots in the default interval [-1,+1].
//
  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
//
//  Set the knot multiplicities.
//
  mlt = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
//
//  Set KIND to the Legendre weight function.
//
  kind = 1;
//
//  ALPHA, BETA not used in Legendre weight function but set anyway.
//
  alpha = 0.0;
  beta  = 0.0;
//
//  Set nonstandard interval A, B.
//
  a = -0.5;
  b = 2.0;
//
//  Shift knots from [-1,1] to [A,B].
//
  for ( i = 0; i < nt; i++ )
  {
    t[i] = ( ( b - a ) * t[i] + ( a + b ) ) / 2.0;
  }
//
//  Call CEIQF to set up the quadrature formula and evaluate it on F.
//
  qfsum = ceiqf ( nt, t, mlt, kind, alpha, beta, a, b, f );

  cout << "\n";
  cout << "  Integral of sin(x) from " << a << " to " << b << "\n";
  cout << "  by Fejer type rule with " << nt << " points\n";
  cout << "  of multiplicity 2.\n";
  cout << "  Quadrature formula:" << setw(24) << setprecision(16) << qfsum << "\n";

  qfsx = cos ( a ) - cos ( b );
  cout << "  Exact value       :" << setw(24) << setprecision(16) << qfsx << "\n";
  cout << "  Error             :" << r8_abs ( qfsum - qfsx ) << "\n";

  delete [] mlt;
  delete [] t;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CLIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double alpha;
  double beta;
  int i;
  int kind;
  int lu;
  int nt;
  double pi = 3.14159265358979323846264338327950;
  double *t;
  double *wts;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST05\n";
  cout << "  Test CLIQFS.\n";
//
//  Number of knots.
//
  nt = 5;
//
//  Set the knots in the default interval [-1,+1].
//
  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
//
//  Request Legendre weight function.
//
  kind = 1;
//
//  ALPHA, BETA not used in Legendre weight function but set anyway.
//
  alpha = 0.0;
  beta  = 0.0;
//
//  LU controls printing.
//  A positive value requests that we compute and print weights, and
//  conduct a moments check.
//
  lu = 6;
//
//  This call returns the WTS array.
//
  wts = cliqfs ( nt, t, kind, alpha, beta, lu );

  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests CLIQF and EIQFS..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int kind;
  int lu;
  int nt;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;
  double *wts;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST06\n";
  cout << "  Test CLIQF and EIQFS.\n";
//
//  Number of knots.
//
  nt = 5;
//
//  Set the knots in the default interval [-1,+1].
//
  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
//
//  Set KIND to the Legendre weight function.
//
  kind = 1;
//
//  ALPHA, BETA not used in Legendre weight function but set anyway.
//
  alpha = 0.0;
  beta  = 0.0;
//
//  Set nonstandard interval A, B.
//
  a = -0.5;
  b = 2.0;
//
//  Shift knots from [-1,1] to [A,B].
//
  for ( i = 0; i < nt; i++ )
  {
    t[i] = ( ( b - a ) * t[i] + ( a + b ) ) / 2.0;
  }
//
//  LU controls printout.
//
  lu = 6;
//
//  Call CLIQF to set up the quadrature formula.
//
  wts = cliqf ( nt, t, kind, alpha, beta, a, b, lu );
//
//  Call EIQFS to evaluate the quadrature formula.
//
  qfsum = eiqfs ( nt, t, wts, f );

  cout << "\n";
  cout << "  Integral of sin(x) from " << a << " to " << b << "\n";
  cout << "  by Fejer type rule with " << nt << " points\n";
  cout << "  of multiplicity 1.\n";
  cout << "  Quadrature formula:" << setw(24) << setprecision(16) << qfsum << "\n";

  qfsx = cos ( a ) - cos ( b );
  cout << "  Exact value       :" << setw(24) << setprecision(16) << qfsx << "\n";
  cout << "  Error             :" << r8_abs ( qfsum - qfsx ) << "\n";

  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests CEGQF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double a;
  double alpha;
  double b;
  double beta;
  int kind;
  int nt;
  double qfsum;
  double qfsx;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST07\n";
  cout << "  Test CEGQF.\n";
//
//  Number of knots.
//
  nt = 12;
//
//  Request exponential weight function.
//
  kind = 7;
//
//  Set ALPHA and BETA.
//
  alpha = 1.0;
  beta  = 0.0;
//
//  Set interval [A,B].
//
  a = -0.5;
  b = 2.0;
//
//  Call CEGQF to compute and evaluate the Gauss quadrature formula.
//
  qfsum = cegqf ( nt, kind, alpha, beta, a, b, f );

  cout << "\n";
  cout << "  Integral of x*sin(x) from " << a << " to " << b << "\n";
  cout << "  by Gauss-exponential rule with " << nt << " points\n";
  cout << "  Quadrature formula:" << setw(24) << setprecision(16) << qfsum << "\n";

  qfsx = ( b - a ) * 0.5 * ( cos ( a ) - cos ( b ) )
    + sin ( b ) + sin ( a ) - 2.0 * sin ( ( a + b ) / 2.0 );

  cout << "  Exact value       :" << setw(24) << setprecision(16) << qfsx << "\n";
  cout << "  Error             :" << r8_abs ( qfsum - qfsx ) << "\n";

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CEGQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
{
  double alpha;
  double beta;
  int kind;
  int nt;
  double qfsum;
  double qfsx;

  cout << "  ----------------------------------------\n";
  cout << "\n";
  cout << "TEST08\n";
  cout << "  Test CEGQFS.\n";
//
//  Number of knots.
//
  nt = 12;
//
//  Request exponential weight function.
//
  kind = 7;
//
//  Set ALPHA and BETA.
//
  alpha = 1.0;
  beta  = 0.0;
//
//  Call CEGQFS to compute and evaluate the Gauss quadrature formula.
//
  qfsum = cegqfs ( nt, kind, alpha, beta, f );

  cout << "\n";
  cout << "  Integral of x*sin(x) from -1 to +1\n";
  cout << "  by Gauss-exponential rule with " << nt << " points.\n";
  cout << "  Quadrature formula: " << setw(24) << setprecision(16) << qfsum << "\n";

  qfsx = cos ( -1.0 ) - cos ( +1.0 );

  cout << "  Exact value       :" << setw(24) << setprecision(16) << qfsx << "\n";
  cout << "  Error             :" << r8_abs ( qfsum - qfsx ) << "\n";

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 calls CGQFS to compute and print generalized Gauss-Hermite rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  int io;
  int kind;
  int nt;
  double *t;
  double *wts;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  Call CGQFS to compute generalized Hermite rules.\n";

  nt = 15;
  kind = 6;
  alpha = 1.0;
  beta = 0.0;
  io = - 6;
  t = new double[nt];
  wts = new double[nt];

  cout << "\n";
  cout << "  NT = " << nt << "\n";
  cout << "  ALPHA = " << alpha << "\n";

  cgqfs ( nt, kind, alpha, beta, io, t, wts );

  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

void test10 ( int nt, int kind, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 calls CDGQF to compute a quadrature formula.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double *t;
  double *wts;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  Call CDGQF to compute a quadrature formula.\n";
  cout << "\n";
  cout << "  KIND = " << kind << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta << "\n";

  t = new double[nt];
  wts = new double[nt];

  cdgqf ( nt, kind, alpha, beta, t, wts );

  cout << "\n";
  cout << " Index     Abscissas                 Weights\n";
  cout << "\n";
  for ( i = 0; i < nt; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(24) << setprecision(16) << t[i]
         << "  " << setw(24) << setprecision(16) << wts[i] << "\n";
  }

  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

void test11 ( int nt, int kind, double alpha, double beta, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 calls CGQF to compute a quadrature formula.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int lo;
  double *t;
  double *wts;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  Call CGQF to compute a quadrature formula with nondefault\n";
  cout << "  values of parameters A and B.\n";
  cout << "\n";
  cout << "  KIND =  " << kind << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta << "\n";
  cout << "  A =     " << a << "\n";
  cout << "  B  =    " << b << "\n";

  lo = 0;
  t = new double[nt];
  wts = new double[nt];

  cgqf ( nt, kind, alpha, beta, a, b, lo, t, wts );

  cout << "\n";
  cout << " Index     Abscissas                 Weights\n";
  cout << "\n";
  for ( i = 0; i < nt; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(24) << setprecision(16) << t[i]
         << "  " << setw(24) << setprecision(16) << wts[i] << "\n";
  }

  delete [] t;
  delete [] wts;

  return;
}
//****************************************************************************80

double f ( double x, int i )

//****************************************************************************80
//
//  Purpose:
//
//    F returns values of the integrand or its derivatives.
//
//  Discussion:
//
//    This function is an example of an integrand function.
//
//    The package can generate quadrature formulas that use derivative
//    information as well as function values.  Therefore, this routine is
//    set up to provide derivatives of any order as well as the function
//    value.  In an actual application, the highest derivative needed
//    is of order one less than the highest knot multiplicity.
//
//    In other words, in the usual case where knots are not repeated,
//    this routine only needs to return function values, not any derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int I, the order of the derivative of F to
//    be evaluated.
//
//    Output, double F, the value of the I-th derivative of F at X.
//
{
  int l;
  double value;

  l = ( i % 4 );

  if ( l == 0 )
  {
    value = sin ( x );
  }
  else if ( l == 1 )
  {
    value = cos ( x );
  }
  else if ( l == 2 )
  {
    value = - sin ( x );
  }
  else if ( l == 3 )
  {
    value = - cos ( x );
  }

  return value;
}
