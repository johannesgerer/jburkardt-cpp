# include <cstdlib>
# include <iostream>
# include <cmath>
# include <iomanip>
# include <ctime>
# include <string>

using namespace std;

# include "toms655.hpp"

//****************************************************************************80

double *cawiq ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int nst, double aj[], double bj[], int *jdf, double zemu )

//****************************************************************************80
//
//  Purpose:
//
//    CAWIQ computes quadrature weights for a given set of knots.
//
//  Discussion:
//
//    This routine is given a set of distinct knots, T, their multiplicities MLT, 
//    the Jacobi matrix associated with the polynomials orthogonal with respect 
//    to the weight function W(X), and the zero-th moment of W(X).
//
//    It computes the weights of the quadrature formula
//
//      sum ( 1 <= J <= NT ) sum ( 0 <= I <= MLT(J) - 1 ) wts(j) d^i/dx^i f(t(j))
//
//    which is to approximate
//
//      integral ( a < x < b ) f(t) w(t) dt
//
//    The routine makes various checks, as indicated below, sets up
//    various vectors and, if necessary, calls for the diagonalization
//    of the Jacobi matrix that is associated with the polynomials
//    orthogonal with respect to W(X) on the interval A, B. 
//
//    Then for each knot, the weights of which are required, it calls the 
//    routine CWIQD which to compute the weights.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int NWTS, the number of weights.
//
//    Input/output, int NDX[NT], associates with each distinct 
//    knot T(J), an integer NDX(J) which is such that the weight to the I-th 
//    derivative value of F at the J-th knot, is stored in
//      WTS(abs(NDX(J))+I) for J = 1,2,...,NT, and I = 0,1,2,...,MLT(J)-1.
//    The sign of NDX includes the following information:
//    > 0, weights are wanted for this knot
//    < 0, weights not wanted for this knot but it is included in the quadrature
//    = 0. means ignore this knot completely.
//
//    Input, int KEY, indicates structure of WTS and NDX.
//    KEY is an integer with absolute value between 1 and 4.
//    The sign of KEY choosed the form of WTS:
//    0 < KEY, WTS in standard form.
//    0 > KEY, J]WTS(J) required.
//    The absolute value has the following effect:
//    1, set up pointers in NDX for all knots in T array (routine CAWIQ does 
//    this).  the contents of NDX are not tested on input and weights are 
//    packed sequentially in WTS as indicated above.
//    2, set up pointers only for knots which have nonzero NDX on input.  All 
//    knots which have a non-zero flag are allocated space in WTS.
//    3, set up pointers only for knots which have NDX > 0 on input.  Space in 
//    WTS allocated only for knots with NDX > 0.
//    4, NDX assumed to be preset as pointer array on input.
//
//    Input, int NST, the dimension of the Jacobi matrix.  
//    NST should be between (N+1)/2 and N.  The usual choice will be (N+1)/2.
//
//    Input/output, double AJ[NST], BJ[NST].
//    If JDF = 0 then AJ contains the  diagonal of the Jacobi matrix and
//    BJ(1:NST-1) contains the subdiagonal.
//    If JDF = 1, AJ contains the eigenvalues of the Jacobi matrix and
//    BJ contains the squares of the elements of the first row of U, the
//    orthogonal matrix which diagonalized the Jacobi matrix as U*D*U'.
//
//    Input/output, int *JDF, indicates whether the Jacobi
//    matrix needs to be diagonalized.
//    0, diagonalization required;
//    1, diagonalization not required.
//
//    Input, double ZEMU, the zero-th moment of the weight 
//    function W(X).
//
//    Output, double CAWIQ[NWTS], the weights.
//
{
  int i;
  int ip;
  int j;
  int jj;
  int jp;
  int k;
  int l;
  int m;
  int mnm;
  int mtj;
  int n;
  double p;
  double prec;
  double *r;
  double tmp;
  double *xk;
  double *wtmp;
  double *wts;
  double *z;

  prec = r8_epsilon ( );

  if ( nt < 1 )
  {
    cerr << "\n";
    cerr << "CAWIQ - Fatal error!\n";
    cerr << "  NT < 1.\n";
    exit ( 1 );
  }
//
//  Check for indistinct knots.
//
  if ( 1 < nt )
  {
    k = nt - 1;
    for ( i = 1; i <= k; i++ )
    {
      tmp = t[i-1];
      l = i + 1;
      for ( j = l; j <= nt; j++ )
      {
        if ( r8_abs ( tmp - t[j-1] ) <= prec )
        {
          cerr << "\n";
          cerr << "CAWIQ - Fatal error!\n";
          cerr << "  Knots too close.\n";
          exit ( 1 );
        }
      }
    }
  }
//
//  Check multiplicities,
//  Set up various useful parameters and
//  set up or check pointers to WTS array.
//
  l = abs ( key );

  if ( l < 1 || 4 < l )
  {
    cerr << "\n";
    cerr << "CAWIQ - Fatal error!\n";
    cerr << "  Magnitude of KEY not between 1 and 4.\n";
    exit ( 1 );
  }

  k = 1;

  if ( l == 1 )
  {
    for ( i = 1; i <= nt; i++ )
    {
      ndx[i-1] = k;
      if ( mlt[i-1] < 1 )
      {
        cerr << "\n";
        cerr << "CAWIQ - Fatal error!\n";
        cerr << "  MLT(I) < 1.\n";
        exit ( 1 );
      }
      k = k + mlt[i-1];
    }
    n = k - 1;
  }
  else if ( l == 2 || l == 3 )
  {
    n = 0;

    for ( i = 1; i <= nt; i++ )
    {
      if ( ndx[i-1] == 0 )
      {
        continue;
      }

      if ( mlt[i-1] < 1 )
      {
        cerr << "\n";
        cerr << "CAWIQ - Fatal error!\n";
        cerr << "  MLT(I) < 1.\n";
        exit ( 1 );
      }

      n = n + mlt[i-1];

      if ( ndx[i-1] < 0 && l == 3 )
      {
        continue;
      }

      ndx[i-1] = abs ( k ) * i4_sign ( ndx[i-1] );
      k = k + mlt[i-1];
    }

    if ( nwts + 1 < k )
    {
      cerr << "\n";
      cerr << "CAWIQ - Fatal error!\n";
      cerr << "  NWTS + 1 < K.\n";
      exit ( 1 );
    }
  }
  else if ( l == 4 )
  {
    for ( i = 1; i <= nt; i++ )
    {
      ip = abs ( ndx[i-1] );

      if ( ip == 0 )
      {
        continue;
      }

      if ( nwts < ip + mlt[i-1] )
      {
        cerr << "\n";
        cerr << "CAWIQ - Fatal error!\n";
        cerr << "  NWTS < IPM.\n";
        exit ( 1 );
      } 

      if ( i == nt )
      {
        break;
      }

      l = i + 1;
      for ( j = l; j <= nt; j ++ )
      {
        jp = abs ( ndx[j-1] );
        if ( jp != 0 )
        {
          if ( jp <= ip + mlt[i-1] && ip <= jp + mlt[j-1] )
          {
            break;
          }
        }
      }
    }
  }
//
//  Test some parameters.
//
  if ( nst < ( n + 1 ) / 2 )
  {
    cerr << "\n";
    cerr << "CAWIQ - Fatal error!\n";
    cerr << "  NST < ( N + 1 ) / 2.\n";
    exit ( 1 );
  }

  if ( zemu <= 0.0 )
  {
    cerr << "\n";
    cerr << "CAWIQ - Fatal error!\n";
    cerr << "  ZEMU <= 0.\n";
    exit ( 1 );
  }

  wts = new double[nwts];
//
//  Treat a quadrature formula with 1 simple knot first.
//
  if ( n <= 1 )
  {
    for ( i = 0; i < nt; i++ )
    {
      if ( 0 < ndx[i] )
      {
        wts[ abs ( ndx[i] ) - 1 ] = zemu;
        return wts;
      }
    }
  }
//
//  Carry out diagonalization if not already done.
//
  if ( *jdf == 0 )
  {
//
//  Set unit vector in work field to get back first row of Q.
//
    z = new double[nst];

    for ( i = 0; i < nst; i++ )
    {
      z[i] = 0.0;
    }
    z[0] = 1.0;
//
//  Diagonalize the Jacobi matrix.
//
    imtqlx ( nst, aj, bj, z );
//
//  Signal Jacobi matrix now diagonalized successfully.
//
    *jdf = 1;
//
//  Save squares of first row of U in subdiagonal array.
//
    for ( i = 0; i < nst; i++ )
    {
      bj[i] = z[i] * z[i];
    }
    delete [] z;
  }
//
//  Find all the weights for each knot flagged.
//
  for ( i = 1; i <= nt; i++ )
  {
    if ( ndx[i-1] <= 0 )
    {
      continue;
    }
    m = mlt[i-1];
    mnm = i4_max ( n - m, 1 );
    l = i4_min ( m, n - m + 1 );
//
//  Set up K-hat matrix for CWIQD with knots according to their multiplicities.
//
    xk = new double[mnm];

    k = 1;
    for ( j = 1; j <= nt; j++ )
    {
      if ( ndx[j-1] != 0 )
      {
        if ( j != i )
        {
          for ( jj = 1; jj <= mlt[j-1]; jj++ )
          {
            xk[k-1] = t[j-1];
            k = k + 1;
          }
        }
      }
    }
//
//  Set up the right principal vector.
//
    r = new double[l];

    r[0] = 1.0 / zemu;
    for ( j = 1; j < l; j++ )
    {
      r[j] = 0.0;
    }
//
//  Pick up pointer for the location of the weights to be output.
//
    k = ndx[i-1];
//
//  Find all the weights for this knot.
//
    wtmp = cwiqd ( m, mnm, l, t[i-1], xk, nst, aj, bj, r );

    delete [] r;
    delete [] xk;

    for ( j = 0; j < m; j++ )
    {
      wts[k-1+j] = wtmp[j];
    }
    delete [] wtmp;

    if ( key < 0 )
    {
      continue;
    }
//
//  Divide by factorials for weights in standard form.
//
    tmp = 1.0;
    for ( j = 1; j < m - 1; j++ )
    {
      p = j;
      tmp = tmp * p;
      wts[k-1+j] = wts[k-1+j] / tmp;
    }
  }
  return wts;
}
//****************************************************************************80

void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with a classical weight function with default values for A and B,
//    and only simple knots.
//
//    There are no moments checks and no printing is done.
//
//    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
  double *aj;
  double *bj;
  double zemu;

  parchk ( kind, 2 * nt, alpha, beta );
//
//  Get the Jacobi matrix and zero-th moment.
//
  aj = new double[nt];
  bj = new double[nt];

  zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );
//
//  Compute the knots and weights.
//
  sgqf ( nt, aj, bj, zemu, t, wts );

  delete [] aj;
  delete [] bj;

  return;
}
//****************************************************************************80

double cegqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    CEGQF computes a quadrature formula and applies it to a function.
//
//  Discussion:
//
//    The user chooses the quadrature formula to be used, as well as the
//    interval (A,B) in which it is applied.
//
//    Note that the knots and weights of the quadrature formula are not
//    returned to the user.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The value I will always be 0.  The value X will always be a knot.
//
//    Output, double CEGQF, the value of the quadrature formula 
//    applied to F.
//
{
  double qfsum;
  double *t;
  double *wts;

  t = new double[nt];
  wts = new double[nt];

  cgqf ( nt, kind, alpha, beta, a, b, t, wts );
//
//  Evaluate the quadrature sum.
//
  qfsum = eiqfs ( nt, t, wts, f );

  delete [] t;
  delete [] wts;

  return qfsum;
}
//****************************************************************************80

double cegqfs ( int nt, int kind, double alpha, double beta, 
  double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    CEGQFS estimates an integral using a standard quadrature formula.
//
//  Discussion:
//
//    The user chooses one of the standard quadrature rules
//    with the default values of A and B.  This routine determines
//    the corresponding weights and evaluates the quadrature formula
//    on a given function.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The value  I will always be 0.  The value X will always be a knot.
//
//    Output, double CEGQFS, the value of the quadrature formula 
//    applied to F.
//
{
  int lu;
  double qfsum;
  double *t;
  double *wts;

  lu = 0;

  t = new double[nt];
  wts = new double[nt];

  cgqfs ( nt, kind, alpha, beta, lu, t, wts );
//
//  Evaluate the quadrature sum.
//
  qfsum = eiqfs ( nt, t, wts, f );

  delete [] t;
  delete [] wts;

  return qfsum;
}
//****************************************************************************80

double ceiqf ( int nt, double t[], int mlt[], int kind, double alpha, 
  double beta, double a, double b, double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    CEIQF constructs and applies a quadrature formula based on user knots.
//
//  Discussion:
//
//    The knots may have multiplicity.  The quadrature interval is over
//    any valid A, B.  A classical weight function is selected by the user.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The highest value of I will be the maximum value in MLT minus
//    one.  The value X will always be a knot.
//
//    Output, double CEIQF, the value of the quadrature formula 
//    applied to F.
//
{
  int i;
  int key;
  int lu;
  int m;
  int n;
  int *ndx;
  double qfsum;
  double *wts;

  lu = 0;
  n = 0;
  for ( i = 0; i < nt; i++ )
  {
    n = n + mlt[i];
  }

  key = 1;
  ndx = new int[nt];

  wts = ciqf ( nt, t, mlt, n, ndx, key, kind, alpha, beta, a, b, lu );

  qfsum = eiqf ( nt, t, mlt, wts, n, ndx, key, f );

  delete [] ndx;
  delete [] wts;

  return qfsum;
}
//****************************************************************************80

double ceiqfs ( int nt, double t[], int mlt[], int kind, double alpha, 
  double beta, double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    CEIQFS computes and applies a quadrature formula based on user knots.
//
//  Discussion:
//
//    The knots may have multiplicity.  The quadrature interval is over
//    the standard interval A, B for the classical weight function selected 
//    by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The highest value of I will be the maximum value in MLT minus
//    one.  The value X will always be a knot.
//
//    Output, double CEIQFS, the value of the quadrature formula 
//    applied to F.
//
{
  int i;
  int key;
  int lu;
  int n;
  int *ndx;
  double qfsum;
  double *wts;

  lu = 0;
  n = 0;
  for ( i = 0; i < nt; i++ )
  {
    n = n + mlt[i];
  }
  ndx = new int[nt];
  key = 1;

  wts = ciqfs ( nt, t, mlt, n, ndx, key, kind, alpha, beta, lu );

  qfsum = eiqf ( nt, t, mlt, wts, n, ndx, key, f );

  delete [] ndx;
  delete [] wts;

  return qfsum;
}
//****************************************************************************80

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    CGQF computes knots and weights of a Gauss quadrature formula.
//
//  Discussion:
//
//    The user may specify the interval (A,B).
//
//    Only simple knots are produced.
//
//    Use routine EIQFS to evaluate this quadrature formula.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints, or
//    other parameters.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
  int i;
  int *mlt;
  int *ndx;
//
//  Compute the Gauss quadrature formula for default values of A and B.
//
  cdgqf ( nt, kind, alpha, beta, t, wts );
//
//  Prepare to scale the quadrature formula to other weight function with 
//  valid A and B.
//
  mlt = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 1;
  }
  ndx = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    ndx[i] = i + 1;
  }
  scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b );

  delete [] mlt;
  delete [] ndx;

  return;
}
//****************************************************************************80

void cgqfs ( int nt, int kind, double alpha, double beta, int lo, double t[], 
  double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    CGQFS computes knots and weights of a Gauss quadrature formula.
//
//  Discussion:
//
//    This routine computes the knots and weights of a Gauss quadrature
//    formula with:
//
//    * a classical weight function with default values for A and B;
//    * only simple knots
//    * optionally print knots and weights and a check of the moments
//
//    Use routine EIQFS to evaluate a quadrature formula computed by 
//    this routine.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, int LO, selects the action.
//    > 0, compute and print knots and weights.  Print moments check.
//    = 0, compute knots and weights.
//    < 0, compute and print knots and weights.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
  int i;
  int key;
  int m;
  int mex;
  int *mlt;
  int mmex;
  int mop;
  int *ndx;
  double *w;
//
//  Check there is enough workfield and assign workfield
//
  key = 1;
  mop = 2 * nt;
  m = mop + 1;
  mex = m + 2;
  mmex = i4_max ( mex, 1 );
//
//  Compute the Gauss quadrature formula for default values of A and B.
//
  cdgqf ( nt, kind, alpha, beta, t, wts );
//
//  Exit if no print required.
//
  if ( lo != 0 )
  {
    mlt = new int[nt];
    for ( i = 0; i < nt; i++ )
    {
      mlt[i] = 1;
    }
    ndx = new int[nt];
    for ( i = 0; i < nt; i++ )
    {
      ndx[i] = i + 1;
    }
    w = new double[mmex];

    chkqfs ( t, wts, mlt, nt, nt, ndx, key, w, mop, mmex, kind, alpha, 
      beta, lo );

    delete [] mlt;
    delete [] ndx;
    delete [] w;
  }
  return;
}
//****************************************************************************80

void chkqf ( double t[], double wts[], int mlt[], int nt, int nwts, int ndx[], 
  int key, int mop, int mex, int kind, double alpha, double beta, int lo, 
  double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CHKQF computes and prints the moments of a quadrature formula.
//
//  Discussion:
//
//    The quadrature formula is based on a clasical weight function with 
//    any valid A, B.
//
//    No check can be made for non-classical weight functions.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, double T[NT], the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int NT, the number of knots.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.  
//    If KEY = 1, then NDX need not be preset.  For more details see the 
//    comments in CAWIQ.
//
//    Input, int KEY, indicates the structure of the WTS
//    array.  It will normally be set to 1.  This will cause the weights to be 
//    packed sequentially in array WTS.  For more details see the comments 
//    in CAWIQ.
//
//    Input, int MOP, the expected order of precision of the 
//    quadrature formula.
//
//    Input, int MEX, the number of moments required to be 
//    tested.  Set MEX = 1 and LO < 0 for no moments check.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, int LO, selects the action to carry out.
//     > 0, print weights and moment tests.
//     = 0, print nothing. compute moment test.
//     < 0, print weights only. don't compute moment tests.
//
//    Input, double A, B, the interval endpoints.
//
{
  int i;
  int izero;
  int neg;
  double *t2;
  double tmp;
  double *w;

  w = new double[mex];

  parchk ( kind, mex, alpha, beta );

  if ( lo != 0 )
  {
    izero = 0;

    cout << "\n";
    cout << "  Interpolatory quadrature formula\n";
    cout << "\n";
    cout << "  Type  Interval       Weight function               Name\n";
    cout << "\n";
    if ( kind == 1 )
    {
      cout << "    1    (a,b)              1.0                    Legendre\n";
    }
    else if ( kind == 2 )
    {
      cout << "    2    (a,b)      ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1\n";
    }
    else if ( kind == 3 )
    {
      cout << "    3    (a,b)      ((b-x)*(x-a))^alpha           Gegenbauer\n";
    }
    else if ( kind == 4 )
    {
      cout << "    4    (a,b)    (b-x)^alpha*(x-a)^beta          Jacobi\n";
    }
    else if ( kind == 5 )
    {
      cout << "    5   (a,+oo)  (x-a)^alpha*exp(-b*(x-a))      Gen Laguerre\n";
    }
    else if ( kind == 6 )
    {
      cout << "    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)  Gen Hermite\n";
    }
    else if ( kind == 7 )
    {
      cout << "    7    (a,b)      |x-(a+b)/2.0|^alpha        Exponential\n";
    }
    else if ( kind == 8 )
    {
      cout << "    8   (a,+oo)    (x-a)^alpha*(x+b)^beta         Rational\n";
    }
    else if ( kind == 9 )
    {
      cout << "    9   (a,b)     (b-x)*(x-a)^(+0.5)         Chebyshev Type 2\n";
    }

    cout << "\n";
    cout << "     Parameters   A          " << a << "\n";
    cout << "                  B          " << b << "\n";
    if ( 3 <= kind && kind <= 8 )
    {
      cout << "                  alpha      " << alpha << "\n";
    }

    if ( kind == 4 || kind == 8 )
    {
      cout << "                  beta       " << beta << "\n";
    }
    chkqfs ( t, wts, mlt, nt, nwts, ndx, key, w, mop, mex, izero, 
      alpha, beta, - abs ( lo ) );
  }

  if ( 0 <= lo )
  {
//
//  Compute the moments in W.
//
    w = scmm ( mex, kind, alpha, beta, a, b );

    if ( kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 7 || kind == 9 )
    {
      tmp = ( b + a ) / 2.0;
    }
    else if ( kind == 5 || kind == 6 || kind == 8 )
    {
      tmp = a;
    }

    t2 = new double[nt];

    for ( i = 0; i < nt; i++ )
    {
      t2[i] = t[i] - tmp;
    }

    neg = -1;
//
//  Check moments.
//
    chkqfs ( t2, wts, mlt, nt, nwts, ndx, key, w, mop, mex, neg, alpha, beta, 
      lo );

    delete [] t2;
  }

  delete [] w;

  return;
}
//****************************************************************************80

void chkqfs ( double t[], double wts[], int mlt[], int nt, int nwts, int ndx[], 
  int key, double w[], int mop, int mex, int kind, double alpha, double beta, 
  int lo )

//****************************************************************************80
//
//  Purpose:
//
//    CHKQFS checks the polynomial accuracy of a quadrature formula.
//
//  Discussion:
//
//    This routine will optionally print weights, and results of a moments test.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, double T[NT], the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int NT, the number of knots.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.  
//    If KEY = 1, then NDX need not be preset.  For more details see the 
//    comments in CAWIQ.
//
//    Input, int KEY, indicates the structure of the WTS
//    array.  It will normally be set to 1.  This will cause the weights to be 
//    packed sequentially in array WTS.  For more details see the comments 
//    in CAWIQ.
//
//    Input/output, double W[MEX], the moments array.
//    This is input only if KIND = 0.
//
//    Input, int MOP, the expected order of precision of the
//    quadrature formula.
//
//    Input, int MEX, the number of moments to be tested.
//    MEX must be at least 1.  Set MEX = 1 and LO < 0 for no moment check.
//
//    Input, int KIND, the rule.
//    0, unknown weight function (the user must set the first MEX moments in
//       array W in this case.)
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, int LO, selects the action to carry out.
//     > 0, print weights and moment tests.
//     = 0, print nothing.  Dompute moment test.
//     < 0, print weights only.  Do not compute moment tests.
//
//  Local Parameters:
//
//    Local, double E[MEX], ER[MEX], the absolute and relative 
//    errors of the quadrature formula applied to (X-DEL)^n.
//
//    Local, double QM[MEX], the values of the quadrature formula
//    applied to (X-DEL)^N.
//
{
  double *e;
  double ek;
  double emn;
  double emx;
  double erest;
  double ern;
  double erx;
  double *er;
  int i;
  int j;
  int jl;
  int k;
  int kindp;
  int kjl;
  int l;
  int m;
  int mx;
  double px;
  double tmp;
  double tmpx;
  double prec;
  double *qm;
//
//  KIND may be set to -1 to allow printing of moments only.
//
//  This feature is only used internally, by CHKQF.
//
  kindp = i4_max ( 0, kind );

  if ( lo != 0 && kind != -1 )
  {
    if ( kindp != 0 )
    {
      cout << "\n";
      cout << "  Interpolatory quadrature formula\n";
      cout << "\n";
      cout << "  Type  Interval       Weight function               Name\n";
      cout << "\n";
      if ( kindp == 1 )
      {
        cout << "    1    (-1,+1)            1.0                    Legendre\n";
      }
      else if ( kindp == 2 )
      {
        cout << "    2    (-1,+1)    ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1\n";
      }
      else if ( kindp == 3 )
      {
        cout << "    3    (-1,+1)    ((b-x)*(x-a))^alpha           Gegenbauer\n";
      }
      else if ( kindp == 4 )
      {
        cout << "    4    (-1,+1)  (b-x)^alpha*(x-a)^beta          Jacobi\n";
      }
      else if ( kindp == 5 )
      {
        cout << "    5   (a,+oo)   (x-a)^alpha*exp(-b*(x-a))      Gen Laguerre\n";
      }
      else if ( kindp == 6 )
      {
        cout << "    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)  Gen Hermite\n";
      }
      else if ( kindp == 7 )
      {
        cout << "    7    (-1,+1)    |x-(a+b)/2.0|^alpha        Exponential\n";
      }
      else if ( kindp == 8 )
      {
        cout << "    8   (0,+oo)    (x-a)^alpha*(x+b)^beta         Rational\n";
      }
      else if ( kindp == 9 )
      {
        cout << "    9    (-1,+1)    ((b-x)*(x-a))^(+0.5)          Chebyshev Type 2\n";
      }

      if ( 3 <= kindp && kindp <= 8 )
      {
        cout << "                  alpha      " << alpha << "\n";
      }

      if ( kindp == 4 || kindp == 8 )
      {
        cout << "                  beta       " << beta << "\n";
      }

    }

    if ( kind != -1 )
    {
      prec = r8_epsilon ( );
      cout << "\n";
      cout << "  Machine precision = " << prec << "\n";
    }

    cout << "\n";
    cout << "           Knots               Mult                Weights\n";
    cout << "\n";

    for ( i = 1; i <= nt; i++ )
    {
      k = abs ( ndx[i-1] );
      if ( k != 0 )
      {
        cout << setw(4) << i
             << setw(26) << setprecision(17) << t[i-1]
             << setw(4) << mlt[i-1]
             << setw(26) << setprecision(17) << wts[k-1] << "\n";
        for ( j = k + 1; j <= k + mlt[i-1] - 1; j++ )
        {
          cout << "                                  "
               << setw(26) << setprecision(17) << wts[j-1] << "\n";
        }
      }
    }
  }

  if ( lo < 0 )
  {
    return;
  }
//
//  Compute the moments in W.
//
  if ( kindp != 0 )
  {
    w = wm ( mex, kindp, alpha, beta );
  }

  e = new double[mex];
  er = new double[mex];
  qm = new double[mex];

  for ( j = 0; j < mex; j++ )
  {
    qm[j] = 0.0;
  }
  erest = 0.0;

  for ( k = 1; k <= nt; k++ )
  {
    tmp = 1.0;
    l = abs ( ndx[k-1] );
    if ( l == 0 )
    {
      continue;
    }

    erest = erest + r8_abs ( wts[l-1] );
    for ( j = 1; j <= mex; j++ )
    {
      qm[j-1] = qm[j-1] + tmp * wts[l-1];
      tmpx = tmp;
      px = 1.0;
      for( jl = 2; jl <= i4_min ( mlt[k-1], mex - j + 1 ); jl++ )
      {
        kjl = j + jl - 1;
        tmpx = tmpx * ( kjl - 1 );
        qm[kjl-1] = qm[kjl-1] + tmpx * wts[l+jl-2] / px;
        if ( key <= 0 )
        {
          px = px * jl;
        }
      }
      tmp = tmp * t[k-1];
    }

  }
  for ( j = 0; j < mex; j++ )
  {
    e[j] = w[j] - qm[j];
    er[j] = e[j] / ( r8_abs ( w[j] ) + 1.0 );
  }
//
//  For some strange weight functions W(1) may vanish.
//
  erest = erest / ( r8_abs ( w[0] ) + 1.0 );

  if ( 0 < lo ) 
  {
    m = mop + 1;
    mx = i4_min ( mop, mex );

    emx = r8_abs ( e[0] );
    emn = emx;
    erx = r8_abs ( er[0] );
    ern = erx;
    for ( k = 1; k < mx; k++ )
    {
      emx = r8_max ( r8_abs ( e[k] ), emx );
      emn = r8_min ( r8_abs ( e[k] ), emn );
      erx = r8_max ( r8_abs ( er[k] ), erx );
      ern = r8_min ( r8_abs ( er[k] ), ern );
    }

    cout << "\n";
    cout << "  Comparison of moments\n";
    cout << "\n";
    cout << "  Order of precision " << mop << "\n";
    cout << "  Errors :    Absolute    Relative\n";
    cout << "  ---------+-------------------------\n";
    cout << "  Minimum :" << setw(12) << setprecision(3) << emn 
         << "  " << setw(12) << setprecision(3) << ern << "\n";
    cout << "  Maximum :" << setw(12) << setprecision(3) << emx 
         << "  " << setw(12) << setprecision(3) << erx << "\n";
    cout << "\n";
    cout << "  Weights ratio       " 
         << setw(12) << setprecision(3) << erest << "\n";

    if ( m <= mex )
    {
      ek = e[m-1];
      for ( j = 1; j <= mop; j++ )
      {
        ek = ek / ( double ) ( j );
      }

      cout << "  Error in " << mop << "th power " 
           << setw(12) << setprecision(3) << e[m-1] << "\n";
      cout << "  Error constant      " 
           << setw(12) << setprecision(3) << ek << "\n";
    }

    cout << "\n";
    cout << "  Moments:\n";
    cout << "\n";
    cout << "            True             from QF            Error      Relative\n";
    cout << "\n";
    for ( j = 1; j <= mx; j++ )
    {
      cout << setw(4) << j
           << setw(19) << setprecision(10) << w[j-1]
           << setw(19) << setprecision(10) << qm[j-1]
           << setw(12) << setprecision(3) << e[j-1]
           << setw(12) << setprecision(3) << er[j-1] << "\n";
    }
    cout << "\n";
    for ( j = m; j <= mex; j++ )
    {
      cout << setw(4) << j
           << setw(19) << setprecision(10) << w[j-1]
           << setw(19) << setprecision(10) << qm[j-1]
           << setw(12) << setprecision(3) << e[j-1]
           << setw(12) << setprecision(3) << er[j-1] << "\n";
    }
  }

  delete [] e;
  delete [] er;
  delete [] qm;

  return;
}
//****************************************************************************80

double *ciqf ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int kind, double alpha, double beta, double a, double b, int lo )

//****************************************************************************80
//
//  Purpose:
//
//    CIQF computes weights for a classical weight function and any interval.
//
//  Discussion:
//
//    This routine compute somes or all the weights of a quadrature formula
//    for a classical weight function with any valid A, B and a given set of 
//    knots and multiplicities.  
//
//    The weights may be packed into the output array WTS according to a 
//    user-defined pattern or sequentially. 
//
//    The routine will also optionally print knots and weights and a check 
//    of the moments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int NWTS, the number of weights.
//
//    Input/output, int NDX[NT], used to index the output
//    array WTS.  If KEY = 1, then NDX need not be preset.  For more
//    details see the comments in CAWIQ.
//
//    Input, int KEY, indicates the structure of the WTS
//    array.  It will normally be set to 1.  This will cause the weights to be 
//    packed sequentially in array WTS.  For more details see the comments 
//    in CAWIQ.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, int LO, selects the actions to perform.
//     > 0, compute and print weights.  Print moments check.
//     = 0, compute weights.
//     < 0, compute and print weights.
//
//    Output, double CIQF[NWTS], the weights.
//
{
  int j;
  int k;
  int l;
  int lu;
  int m;
  int mex;
  int mop;
  double *st;
  double *wts;

  m = 1;
  l = abs ( key );

  for ( j = 0; j < nt; j++ )
  {
    if ( l == 1 || abs ( ndx[j] ) != 0 )
    {
      m = m + mlt[j];
    }
  }

  if ( nwts + 1 < m )
  {
    cerr << "\n";
    cerr << "CIQF - Fatal error!\n";
    cerr << "  NWTS + 1 < M.\n";
    exit ( 1 );
  }

  mex = 2 + m;
//
//  Scale the knots to default A, B.
//
  st = sct ( nt, t, kind, a, b );
//
//  Compute the weights.
//
  lu = 0;

  wts = ciqfs ( nt, st, mlt, nwts, ndx, key, kind, alpha, beta, lu );
//
//  Don't scale user's knots - only scale weights.
//
  scqf ( nt, st, mlt, wts, nwts, ndx, wts, st, kind, alpha, beta, a, b );

  if ( lo != 0 )
  {
    mop = m - 1;

    chkqf ( t, wts, mlt, nt, nwts, ndx, key, mop, mex, kind,
      alpha, beta, lo, a, b );
  }

  return wts;
}
//****************************************************************************80

double *ciqfs ( int nt, double t[], int mlt[], int nwts, int ndx[], int key, 
  int kind, double alpha, double beta, int lo )

//****************************************************************************80
//
//  Purpose:
//
//    CIQFS computes some weights of a quadrature formula in the default interval.
//
//  Discussion:
//
//    This routine computes some or all the weights of a quadrature formula 
//    for a classical weight function with default values of A and B,
//    and a given set of knots and multiplicities. 
//
//    The weights may be packed into the output array WTS according to a 
//    user-defined pattern or sequentially. 
//
//    The routine will also optionally print knots and weights and a check of 
//    the moments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, int NWTS, the number of weights.
//
//    Input/output, int NDX[NT],  used to index the output
//    array WTS.  If KEY = 1, then NDX need not be preset.  For more
//    details see the comments in CAWIQ.
//
//    Input, int KEY, indicates the structure of the WTS
//    array.  It will normally be set to 1.  For more details see
//    the comments in CAWIQ.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, int LO, selects the actions to perform.
//     > 0, compute and print weights.  Print moments check.
//     = 0, compute weights.
//     < 0, compute and print weights.
//
//    Output, double CIQFS[NWTS], the weights.
//
{
  double *aj;
  double *bj;
  int j;
  int jdf;
  int l;
  int m;
  int mex;
  int mmex;
  int mop;
  int n;
  int nst;
  double *w;
  double *wts;
  double zemu;

  jdf = 0;
  n = 0;
  l = abs ( key );

  for ( j = 0; j < nt; j++ )
  {
    if ( l == 1 || abs ( ndx[j] ) != 0 )
    {
      n = n + mlt[j];
    }
  }
//
//  N knots when counted according to multiplicity.
//
  if ( nwts < n )
  {
    cerr << "\n";
    cerr << "CIQFS - Fatal error!\n";
    cerr << "  NWTS < N.\n";
    exit ( 1 );
  }

  m = n + 1;
  mex = 2 + m;
  nst = m / 2;
//
//  Get the Jacobi matrix.
//
  aj = new double[nst];
  bj = new double[nst];

  zemu = class_matrix ( kind, nst, alpha, beta, aj, bj );
//
//  Call weights routine.
//
  wts = cawiq ( nt, t, mlt, n, ndx, key, nst, aj, bj, &jdf, zemu );

  delete [] aj;
  delete [] bj;
//
//
//  Call checking routine.
//
  if ( lo != 0 )
  {
    mop = m - 1;

    w = new double[mex];

    chkqfs ( t, wts, mlt, nt, n, ndx, key, w, mop, mex, kind, 
      alpha, beta, lo );

    delete [] w;
  }
  return wts;
}
//****************************************************************************80

double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] )

//****************************************************************************80
//
//  Purpose:
//
//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
//
//  Discussion:
//
//    This routine computes the diagonal AJ and sub-diagonal BJ
//    elements of the order M tridiagonal symmetric Jacobi matrix
//    associated with the polynomials orthogonal with respect to
//    the weight function specified by KIND.
//
//    For weight functions 1-7, M elements are defined in BJ even
//    though only M-1 are needed.  For weight function 8, BJ(M) is
//    set to zero.
//
//    The zero-th moment of the weight function is returned in ZEMU.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, int M, the order of the Jacobi matrix.
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
//    of the Jacobi matrix.
//
//    Output, double CLASS_MATRIX, the zero-th moment.
//
{
  double a2b2;
  double ab;
  double aba;
  double abi;
  double abj;
  double abti;
  double apone;
  int i;
  double pi = 3.14159265358979323846264338327950;
  double temp;
  double temp2;
  double zemu;

  temp = r8_epsilon ( );

  parchk ( kind, 2 * m - 1, alpha, beta );

  temp2 = 0.5;

  if ( 500.0 * temp < r8_abs ( pow ( gamma ( temp2 ), 2 ) - pi ) )
  {
    cerr << "\n";
    cerr << "CLASS_MATRIX - Fatal error!\n";
    cerr << "  Gamma function does not match machine parameters.\n";
    exit ( 1 );
  }

  if ( kind == 1 )
  {
    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 2 )
  {
    zemu = pi;

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] = sqrt ( 0.5 );
    for ( i = 1; i < m; i++ )
    {
      bj[i] = 0.5;
    }
  }
  else if ( kind == 3 )
  {
    ab = alpha * 2.0;
    zemu = pow ( 2.0, ab + 1.0 ) * pow ( gamma ( alpha + 1.0 ), 2 )
      / gamma ( ab + 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
    for ( i = 2; i <= m; i++ )
    {
      bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * pow ( i + alpha, 2 ) - 1.0 ) );
    }
  }
  else if ( kind == 4 )
  {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow ( 2.0, ab + 1.0 ) * gamma ( alpha + 1.0 ) 
      * gamma ( beta + 1.0 ) / gamma ( abi );
    aj[0] = ( beta - alpha ) / abi;
    bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
      / ( ( abi + 1.0 ) * abi * abi ) );
    a2b2 = beta * beta - alpha * alpha;

    for ( i = 2; i <= m; i++ )
    {
      abi = 2.0 * i + ab;
      aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi * abi;
      bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
        / ( ( abi - 1.0 ) * abi ) );
    }
  }
  else if ( kind == 5 )
  {
    zemu = gamma ( alpha + 1.0 );

    for ( i = 1; i <= m; i++ )
    {
      aj[i-1] = 2.0 * i - 1.0 + alpha;
      bj[i-1] = sqrt ( i * ( i + alpha ) );
    }
  }
  else if ( kind == 6 )
  {
    zemu = gamma ( ( alpha + 1.0 ) / 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
    }
  }
  else if ( kind == 7 )
  {
    ab = alpha;
    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 8 )
  {
    ab = alpha + beta;
    zemu = gamma ( alpha + 1.0 ) * gamma ( - ( ab + 1.0 ) ) 
      / gamma ( - beta );
    apone = alpha + 1.0;
    aba = ab * apone;
    aj[0] = - apone / ( ab + 2.0 );
    bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for ( i = 2; i <= m; i++ )
    {
      abti = ab + 2.0 * i;
      aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
    }

    for ( i = 2; i <= m - 1; i++ )
    {
      abti = ab + 2.0 * i;
      bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
        / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
    }
    bj[m-1] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      bj[i] = sqrt ( bj[i] );
    }
  }
  else if ( kind == 9 )
  {
    zemu = pi / 2.0;

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 0; i < m; i++ )
    {
      bj[i] = 0.5;
    }
  }

  return zemu;
}
//****************************************************************************80

double *cliqf ( int nt, double t[], int kind, double alpha, double beta, 
  double a, double b, int lo )

//****************************************************************************80
//
//  Purpose:
//
//    CLIQF computes a classical quadrature formula, with optional printing.
//
//  Discussion:
//
//    This routine computes all the weights of an interpolatory
//    quadrature formula with
//    1. only simple knots and
//    2. a classical weight function with any valid A and B, and
//    3. optionally prints the knots and weights and a check of the moments.
//
//    To evaluate this quadrature formula for a given function F,
//    call routine EIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
//    Input, int LO, indicates what is to be done.
//    > 0, compute and print weights and moments check.
//    = 0, compute weights.
//    < 0, compute and print weights.
//
//    Output, double CLIQF[NT], the weights.
//
{
  int i;
  int key;
  int *mlt;
  int *ndx;
  double *wts;

  key = 1;
  mlt = new int[nt];
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 1;
  }
  ndx = new int[nt];

  wts = ciqf ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, a, b, lo );

  delete [] mlt;
  delete [] ndx;

  return wts;
}
//****************************************************************************80

double *cliqfs ( int nt, double t[], int kind, double alpha, double beta, 
  int lo )

//****************************************************************************80
//
//  Purpose:
//
//    CLIQFS computes the weights of a quadrature formula in the default interval.
//
//  Discussion:
//
//    This routine computes the weights of an interpolatory quadrature formula
//    with a classical weight function, in the default interval A, B,
//    using only simple knots.
//
//    It can optionally print knots and weights and a check of the moments.
//
//    To evaluate a quadrature computed by CLIQFS, call EIQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, int LO, chooses the printing option.
//     > 0, compute weights, print them, print the moment check results.
//     0, compute weights.
//     < 0, compute weights and print them.
//
//    Output, double CLIQFS[NT], the weights.
//
{
  int i;
  int key;
  int *mlt;
  int *ndx;
  double *wts;

  key = 1;
  mlt = new int[nt];

  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 1;
  }
  ndx = new int[nt];

  wts = ciqfs ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, lo );

  delete [] mlt;
  delete [] ndx;
 
  return wts;
}
//****************************************************************************80

double *cwiqd ( int m, int nm, int l, double v, double xk[], int nstar, 
  double phi[], double a[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CWIQD computes all the weights for a given knot.
//
//  Discussion:
//
//    The variable names correspond to the 1982 reference, and explanations of
//    some of the terminology may be found there.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Jaroslav Kautsky, Sylvan Elhay,
//    Calculation of the Weights of Interpolatory Quadratures,
//    Numerische Mathematik,
//    Volume 40, 1982, pages 407-422.
//
//  Parameters:
//
//    Input, int M, the multiplicity of the knot in question.
//
//    Input, int NM, is equal to max ( N - M, 1 ), where N is 
//    the number of knots used, counted according to multiplicity.
//
//    Input, int L, min ( M, N - M + 1), where N is the number
//    of knots used, counted according to multiplicity.
//
//    Input, double V,  the knot in question.
//
//    Input, double XK[NM], all but the last M entries in the 
//    diagonal of K-hat.
//
//    Input, int NSTAR, the dimension of the Jacobi matrix.
//
//    Input, double PHI[NSTAR], the eigenvalues of the Jacobi matrix.
//
//    Input, double A[NSTAR], the square of the first row of the 
//    orthogonal matrix that diagonalizes the Jacobi matrix.
//
//    Input, double R[L], used to compute the right 
//    principal vectors.
//
//    Output, double CWIQD[M], the weights.
//
{
  double *d;
  int i;
  int j;
  int jr;
  int k;
  int last;
  int minil;
  double sum;
  double tmp;
  double *wf;
  double *y;
  double *z;

  d = new double[m];
  wf = new double[nstar];
  y = new double[m];
  z = new double[m];
//
//  Compute products required for Y-hat.
//
  for ( j = 0; j < nstar; j++ )
  {
    wf[j] = a[j];
    for (i = 0; i < nm; i++ )
    {
      wf[j] = wf[j] * ( phi[j] - xk[i] );
    }
  }
//
//  Compute Y-hat.
//
  for ( i = 0; i < m; i++ )
  {
    sum = 0.0;
    for ( j = 0; j < nstar; j++ )
    {
      sum = sum + wf[j];
      wf[j] = wf[j] * ( phi[j] - v );
    }
    y[i] = sum;
  }
//
//  If N = 1 the right principal vector is already in R.
//  Otherwise compute the R-principal vector of grade M-1.
//
  for ( i = 1; i <= nm; i++ )
  {
    tmp = v - xk[i-1];

    last = i4_min ( l, i + 1 );
    for ( jr = 2; jr <= last; jr++ )
    {
      j = last - jr + 2;
      r[j-1] = tmp * r[j-1] + r[j-2];
    }
    r[0] = tmp * r[0];
  }
//
//  Compute left principal vector(s) and weight for highest derivative.
//  The following statement contains the only division in this
//  routine.  Any test for overflow should be made after it.
//
  d[m-1] = y[m-1] / r[0];

  if ( m == 1 )
  {
    delete [] wf;
    delete [] y;
    delete [] z;

    return d;
  }
//
//  Compute left principal vector.
//
  z[0] = 1.0 / r[0];
  for ( i = 2; i <= m; i++ )
  {
    sum = 0.0;
    minil = i4_min ( i, l );
    for ( j = 2; j <= minil; j++ )
    {
      k = i - j + 1;
      sum = sum + r[j-1] * z[k-1];
    }
    z[i-1] = - sum * z[0];
  }
//
//  Accumulate weights.
//
  for ( i = 2; i <= m; i++ )
  {
    sum = 0.0;
    for ( j = 1; j <= i; j++ )
    {
      k = m - i + j;
      sum = sum + z[j-1] * y[k-1];
    }
    k = m - i + 1;
    d[k-1] = sum;
  }

  delete [] wf;
  delete [] y;
  delete [] z;

  return d;
}
//****************************************************************************80

double eiqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  int key, double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    EIQF evaluates an interpolatory quadrature formula.
//
//  Discussion:
//
//   The knots, weights and integrand are supplied.
//
//   All knots with nonzero NDX are used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.  
//    If KEY = 1, then NDX need not be preset.  For more details see the 
//    comments in CAWIQ.
//
//    Input, int KEY, indicates the structure of the WTS
//    array.  It will normally be set to 1.  This will cause the weights to be 
//    packed sequentially in array WTS.  For more details see the comments 
//    in CAWIQ.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The highest value of I will be the maximum value in MLT minus
//    one.  The value X will always be a knot.
//
//    Output, double EIQF, the value of the quadrature formula 
//    applied to F.
//
{
  int i;
  int j;
  int l;
  double p;
  double qfsum;

  l = abs ( key );

  if ( l < 1 || 4 < l )
  {
    cerr << "\n";
    cerr << "EIQF - Fatal error!\n";
    cerr << "  Magnitude of KEY must be between 1 and 4.\n";
    exit ( 1 );
  }

  qfsum = 0.0;
  for ( j = 0; j < nt; j++ )
  {
    l = abs ( ndx[j] );
    if ( l != 0 )
    {
      p = 1.0;
      for ( i = 0; i < mlt[j]; i++ )
      {
        qfsum = qfsum + wts[l+i-1] * f ( t[j], i ) / p;
        if ( key <= 0 )
        {
          p = p * ( i + 1 );
        }
      }
    }
  }
  return qfsum;
}
//****************************************************************************80

double eiqfs ( int nt, double t[], double wts[], double f ( double x, int i ) )

//****************************************************************************80
//
//  Purpose:
//
//    EIQFS evaluates a quadrature formula defined by CLIQF or CLIQFS.
//
//  Discussion:
//
//    This routine evaluates an interpolatory quadrature formula with all knots 
//    simple and all knots included in the quadrature.  This routine will be used
//    typically after CLIQF or CLIQFS has been called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the knots.
//
//    Input, double WTS[NT], the weights.
//
//    Input, double F ( double X, int I ), the name of a routine which
//    evaluates the function and some of its derivatives.  The routine
//    must return in F the value of the I-th derivative of the function
//    at X.  The value of I will always be 0.  The value X will always be a knot.
//
//    Output, double EIQFS, the value of the quadrature formula 
//    applied to F.
//
{
  int j;
  double qfsum;

  qfsum = 0.0;
  for ( j = 0; j < nt; j++ )
  {
    qfsum = qfsum + wts[j] * f ( t[j], 0 );
  }
  return qfsum;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  int value;

  if ( i < 0 ) 
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentially) of applying the orthogonal transformations
//    directly to Z as they are generated.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( r8_abs ( e[m-1] ) <= prec * ( r8_abs ( d[m-1] ) + r8_abs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cerr << "\n";
        cerr << "IMTQLX - Fatal error!\n";
        cerr << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r = sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + r8_abs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( r8_abs ( g ) <= r8_abs ( f ) )
        {
          c = g / f;
          r = sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r = sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************80

void parchk ( int kind, int m, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    PARCHK checks parameters ALPHA and BETA for classical weight functions. 
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, int M, the order of the highest moment to
//    be calculated.  This value is only needed when KIND = 8.
//
//    Input, double ALPHA, BETA, the parameters, if required
//    by the value of KIND.
//
{
  double tmp;

  if ( kind <= 0 )
  {
    cerr << "\n";
    cerr << "PARCHK - Fatal error!\n";
    cerr << "  KIND <= 0.\n";
    exit ( 1 );
  }
//
//  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
//
  if ( 3 <= kind && kind <= 8 && alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "PARCHK - Fatal error!\n";
    cerr << "  3 <= KIND and ALPHA <= -1.\n";
    exit ( 1 );
  }
//
//  Check BETA for Jacobi.
//
  if ( kind == 4 && beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "PARCHK - Fatal error!\n";
    cerr << "  KIND == 4 and BETA <= -1.0.\n";
    exit ( 1 );
  }
//
//  Check ALPHA and BETA for rational.
//
  if ( kind == 8 )
  {
    tmp = alpha + beta + m + 1.0;
    if ( 0.0 <= tmp || tmp <= beta )
    {
      cerr << "\n";
      cerr << "PARCHK - Fatal error!\n";
      cerr << "  KIND == 8 but condition on ALPHA and BETA fails.\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the 
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  double value;

  value = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + value )  )
  {
    value = value / 2.0;
  }

  value = 2.0 * value;

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

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double *scmm ( int m, int kind, double alpha, double beta, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SCMM computes moments of a classical weight function scaled to [A,B].
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int M, the number of moments.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
//    Output, double W(M), the scaled moments.
//
{
  double al;
  double be;
  int i;
  double p;
  double q;
  double temp;
  double tmp;
  double *w;

  temp = r8_epsilon ( );

  if ( kind == 1 )
  {
    al = 0.0;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = pow ( q, al + be + 1.0 );
  }
  else if ( kind == 2 )
  {
    al = -0.5;
    be = -0.5;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = pow ( q, al + be + 1.0 );
  }
  else if ( kind == 3 )
  {
    al = alpha;
    be = alpha;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = pow ( q, al + be + 1.0 );
  }
  else if ( kind == 4 )
  {
    al = alpha;
    be = beta;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = pow ( q, al + be + 1.0 );
  }
  else if ( kind == 5 )
  {
    if ( b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B <= 0!\n";
      exit ( 1 );
    }
    q = 1.0 / b;
    p = pow ( q, alpha + 1.0 );
  }
  else if ( kind == 6 )
  {
    if ( b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B <= 0!\n";
      exit ( 1 );
    }
    q = 1.0 / sqrt ( b );
    p = pow ( q, alpha + 1.0 );
  }
  else if ( kind == 7 )
  {
    al = alpha;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = pow ( q, al + be + 1.0 );
  }
  else if ( kind == 8 )
  {
    if ( a + b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  A + B <= 0\n";
      exit ( 1 );
    }
    q = a + b;
    p = pow ( q, alpha + beta + 1.0 );
  }
  else if ( kind == 9 )
  {
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCMM - Fatal error!\n";
      cerr << "  B - A too small!\n";
      exit ( 1 );
    }
    q = ( b - a ) / 2.0;
    p = q * q;
  }
//
//  Compute the moments in W.
//
  w = wm ( m, kind, alpha, beta );

  tmp = p;

  for ( i = 0; i < m; i++ )
  {
    w[i] = w[i] * tmp;
    tmp = tmp * q;
  }

  return w;
}
//****************************************************************************80

void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    SCQF scales a quadrature formula to a nonstandard interval.
//
//  Discussion:
//
//    The arrays WTS and SWTS may coincide.
//
//    The arrays T and ST may coincide.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the original knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.  
//    For more details see the comments in CAWIQ.
//
//    Output, double SWTS[NWTS], the scaled weights.
//
//    Output, double ST[NT], the scaled knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
{
  double al;
  double be;
  int i;
  int k;
  int l;
  double p;
  double shft;
  double slp;
  double temp;
  double tmp;

  temp = r8_epsilon ( );

  parchk ( kind, 1, alpha, beta );

  if ( kind == 1 )
  {
    al = 0.0;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 2 )
  {
    al = -0.5;
    be = -0.5;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 3 )
  {
    al = alpha;
    be = alpha;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 4 )
  {
    al = alpha;
    be = beta;

    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 5 )
  {
    if ( b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  B <= 0\n";
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 6 )
  {
    if ( b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  B <= 0.\n";
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / sqrt ( b );
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 7 )
  {
    al = alpha;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 8 )
  {
    if ( a + b <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  A + B <= 0.\n";
      exit ( 1 );
    }
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  }
  else if ( kind == 9 )
  {
    al = 0.5;
    be = 0.5;
    if ( r8_abs ( b - a ) <= temp )
    {
      cerr << "\n";
      cerr << "SCQF - Fatal error!\n";
      cerr << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }

  p = pow ( slp, al + be + 1.0 );

  for ( k = 0; k < nt; k++ )
  {
    st[k] = shft + slp * t[k];
    l = abs ( ndx[k] );

    if ( l != 0 )
    {
      tmp = p;
      for ( i = l - 1; i <= l - 1 + mlt[k] - 1; i++ )
      {
        swts[i] = wts[i] * tmp;
        tmp = tmp * slp;
      }
    }
  }
  return;
}
//****************************************************************************80

double *sct ( int nt, double t[], int kind, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SCT rescales distinct knots to an interval [A,B].
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the original knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double A, B, the interval endpoints for which the
//    knots ST should be scaled.
//
//    Output, double SCT[NT], the scaled knots.
//
{
  double bma;
  int i;
  double shft;
  double slp;
  double *st;
  double tmp;

  if ( kind < 1 || 9 < kind )
  {
    cerr << "\n";
    cerr << "SCT - Fatal error!\n";
    cerr << "  KIND falls outside range of 1 to 8.\n";
    exit ( 1 );
  }

  if ( kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 7 || kind == 9 )
  {
    tmp = r8_epsilon ( );
    bma = b - a;

    if ( bma <= tmp )
    {
      cerr << "\n";
      cerr << "SCT - Fatal error!\n";
      cerr << "  B - A too small.\n";
      exit ( 1 );
    }
    slp = 2.0 / bma;
    shft = - ( a + b ) / bma;
  }
  else if ( kind == 5 )
  {
    if ( b < 0.0 )
    {
      cerr << "\n";
      cerr << "SCT - Fatal error!\n";
      cerr << "  B < 0.\n";
      exit ( 1 );
    }
    slp = b;
    shft = - a * b;
  }
  else if ( kind == 6 )
  {
    if ( b < 0.0 )
    {
      cerr << "\n";
      cerr << "SCT - Fatal error!\n";
      cerr << "  B < 0.\n";
      exit ( 1 );
    }
    slp = sqrt ( b );
    shft = - a * slp;
  }
  else if ( kind == 8 )
  {
    slp = 1.0 / ( a + b );

    if ( slp <= 0.0 )
    {
      cerr << "\n";
      cerr << "SCT - Fatal error.\n";
      cerr << "  1 / ( A + B ) <= 0.\n";
      exit ( 1 );
    }
    shft = - a * slp;
  }

  st = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    st[i] = shft + slp * t[i];
  }
  return st;
}
//****************************************************************************80

void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    SGQF computes knots and weights of a Gauss Quadrature formula.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with simple knots from the Jacobi matrix and the zero-th
//    moment of the weight function, using the Golub-Welsch technique.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double AJ[NT], the diagonal of the Jacobi matrix.
//
//    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
//    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
//
//    Input, double ZEMU, the zero-th moment of the weight function.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
  int i;
//
//  Exit if the zero-th moment is not positive.
//
  if ( zemu <= 0.0 )
  {
    cerr << "\n";
    cerr << "SGQF - Fatal error!\n";
    cerr << "  ZEMU <= 0.\n";
    exit ( 1 );
  }
//
//  Set up vectors for IMTQLX.
//
  for ( i = 0; i < nt; i++ )
  {
    t[i] = aj[i];
  }
  wts[0] = sqrt ( zemu );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i];
  }

  return;
}
//****************************************************************************80

void timestamp ( )

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
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double *wm ( int m, int kind, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    WM evaluates the first M moments of classical weight functions.
//
//  Discussion:
//
//    W(K) = Integral ( A <= X <= B ) X**(K-1) * W(X) dx
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int M, the number of moments to evaluate.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double WM[M], the first M moments.
//
{
  double als;
  int i;
  int ja;
  int jb;
  int k;
  double pi = 3.14159265358979323846264338327950;
  double rk;
  double sum;
  double tmpa;
  double tmpb;
  double trm;
  double *w;

  parchk ( kind, m, alpha, beta );

  w = new double[m];

  for ( k = 2; k <= m; k = k + 2 )
  {
    w[k-1] = 0.0;
  }

  if ( kind == 1 )
  {
    for ( k = 1; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = 2.0 / rk;
    }
  }
  else if ( kind == 2 )
  {
    w[0] = pi;
    for ( k = 3; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = w[k-3] * ( rk - 2.0 ) / ( rk - 1.0 );
    }
  }
  else if ( kind == 3 )
  {
    w[0] = sqrt ( pi ) * gamma ( alpha + 1.0 ) 
      / gamma ( alpha + 3.0 / 2.0 );

    for ( k = 3; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = w[k-3] * ( rk - 2.0 ) / ( 2.0 * alpha + rk );
    }
  }
  else if ( kind == 4 )
  {
    als = alpha + beta + 1.0;
    w[0] = pow ( 2.0, als ) * gamma ( alpha + 1.0 ) 
      / gamma ( als + 1.0 ) * gamma ( beta + 1.0 );

    for ( k = 2; k <= m; k++ )
    {
      sum = 0.0;
      trm = 1.0;
      rk = ( double ) ( k );

      for ( i = 0; i <= ( k - 2 ) / 2; i++ )
      {
        tmpa = trm;
        for ( ja = 1; ja <= 2 * i; ja++ )
        {
          tmpa = tmpa * ( alpha + ja ) / ( als + ja );
        }

        for ( jb = 1; jb <= k - 2 * i - 1; jb++ )
        {
          tmpa = tmpa * ( beta + jb ) / ( als + 2 * i + jb );
        }

        tmpa = tmpa / ( 2 * i + 1.0 ) * 
          ( 2 * i * ( beta + alpha ) + beta - ( rk - 1.0 ) * alpha ) 
          / ( beta + rk - 2 * i - 1.0 );
        sum = sum + tmpa;

        trm = trm * ( rk - 2 * i - 1.0 ) 
          / ( 2 * i + 1.0 ) * ( rk - 2 * i - 2.0 ) / ( 2 * i + 2.0 );
      }

      if ( ( k % 2 ) != 0 )
      {
        tmpb = 1.0;
        for ( i = 1; i <= k - 1; i++ )
        {
          tmpb = tmpb * ( alpha + i ) / ( als + i );
        }
        sum = sum + tmpb;
      }
      w[k-1] = sum * w[0];
    }
  }
  else if ( kind == 5 )
  {
    w[0] = gamma ( alpha + 1.0 );

    for ( k = 2; k <= m; k++ )
    {
      rk = ( double ) ( k );
      w[k-1] = ( alpha + rk - 1.0 ) * w[k-2];
    }
  }
  else if ( kind == 6 )
  {
    w[0] = gamma ( ( alpha + 1.0 ) / 2.0 );

    for ( k = 3; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = w[k-3] * ( alpha + rk - 2.0 ) / 2.0;
    }
  }
  else if ( kind == 7 )
  {
    als = alpha;
    for ( k = 1; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = 2.0 / ( rk + als );
    }
  }
  else if ( kind == 8 )
  {
    w[0] = gamma ( alpha + 1.0 ) 
      * gamma ( - alpha - beta - 1.0 ) 
      / gamma ( - beta );

    for ( k = 2; k <= m; k++ )
    {
      rk = ( double ) ( k );
      w[k-1] = - w[k-2] * ( alpha + rk - 1.0 ) / ( alpha + beta + rk );
    }

  }
  else if ( kind == 9 )
  {
    w[0] = pi / 2.0;
    for ( k = 3; k <= m; k = k + 2 )
    {
      rk = ( double ) ( k );
      w[k-1] = w[k-3] * ( rk - 2.0 ) / ( rk + 1.0 );
    }
  }

  return w;
}
//****************************************************************************80

double *wtfn ( double t[], int nt, int kind, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    WTFN evaluates the classical weight functions at given points.
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
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, double T[NT], the points where the weight function
//    is to be evaluated.
//
//    Input, int NT, the number of evaluation points.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double WTFN[NT], the value of the weight function.
//
{
  int i;
  double *w;

  parchk ( kind, 1, alpha, beta );

  w = new double[nt];

  if ( kind == 1 )
  {
    for ( i = 0; i < nt; i++ )
    {
      w[i] = 1.0;
    }
  }
  else if ( kind == 2 )
  {
    for ( i = 0; i < nt; i++ )
    {
      w[i] = 1.0 / sqrt ( ( 1.0 - t[i] ) * ( 1.0 + t[i] ) );
    }
  }
  else if ( kind == 3 )
  {
    if ( alpha == 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = pow ( ( 1.0 - t[i] ) * ( 1.0 + t[i] ), alpha );
      }
    }
  }
  else if ( kind == 4 )
  {
    if ( alpha == 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = pow ( 1.0 - t[i], alpha );
      }
    }
    if ( beta != 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = w[i] * pow ( 1.0 + t[i], beta );
      }
    }
  }
  else if ( kind == 5 )
  {
    if ( alpha == 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = exp ( - t[i] );
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = exp ( - t[i] ) * pow ( t[i], alpha );
      }
    }
  }
  else if ( kind == 6 )
  {
    if ( alpha == 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = exp ( - t[i] * t[i] );
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = exp ( - t[i] * t[i] ) * pow ( r8_abs ( t[i] ), alpha );
      }
    }
  }
  else if ( kind == 7 )
  {
    if ( alpha != 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = pow ( r8_abs ( t[i] ), alpha );
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = 1.0;
      }
    }
  }
  else if ( kind == 8 )
  {
    if ( alpha == 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = pow ( t[i], alpha );
      }
    }
    if ( beta != 0.0 )
    {
      for ( i = 0; i < nt; i++ )
      {
        w[i] = w[i] * pow ( 1.0 + t[i], beta );
      }
    }
  }
  else if ( kind == 9 )
  {
    for ( i = 0; i < nt; i++ )
    {
      w[i] = sqrt ( ( 1.0 - t[i] ) * ( 1.0 + t[i] ) );
    }
  }
  return w;
}


