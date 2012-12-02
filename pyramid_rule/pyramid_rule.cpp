# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <ctime>
# include <string>
# include <fstream>

using namespace std;

int main ( int argc, char *argv[] );
void jacobi_compute ( int order, double alpha, double beta, double x[], 
  double w[] );
void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] );
void jacobi_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] );
void legendre_compute ( int order, double x[], double w[] );
void pyramid_handle ( int legendre_order, int jacobi_order, string filename );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_gamma ( double x );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PYRAMID_RULE.
//
//  Discussion:
//
//    This program computes a quadrature rule for a pyramid
//    and writes it to a file.
//
//    The user specifies:
//    * the LEGENDRE_ORDER (number of points in the X and Y dimensions)
//    * the JACOBI_ORDER (number of points in the Z dimension)
//    * FILENAME, the root name of the output files.
//
//    The integration region is:
//
//      - ( 1 - Z ) <= X <= 1 - Z
//      - ( 1 - Z ) <= Y <= 1 - Z
//                0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int jacobi_order;
  int legendre_order;

  timestamp ( );
  cout << "\n";
  cout << "PYRAMID_RULE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compute a quadrature rule for approximating\n";
  cout << "  the integral of a function over a pyramid.\n";
  cout << "\n";
  cout << "  The user specifies:\n";
  cout << "\n";
  cout << "  LEGENDRE_ORDER, the order of the Legendre rule for X and Y.\n";
  cout << "  JACOBI_ORDER, the order of the Jacobi rule for Z,\n";
  cout << "  FILENAME, the prefix of the three output files:\n";
  cout << "\n";
  cout << "    filename_w.txt - the weight file\n";
  cout << "    filename_x.txt - the abscissa file.\n";
  cout << "    filename_r.txt - the region file.\n";
//
//  Get the Legendre order.
//
  if ( 1 < argc )
  {
    legendre_order = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the Legendre rule order:\n";
    cin >> legendre_order;
  }

  cout << "\n";
  cout << "  The requested Legendre order of the rule is " << legendre_order << "\n";
//
//  Get the Jacobi order.
//
  if ( 2 < argc )
  { 
    jacobi_order = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the Jacobi rule order:\n";
    cin >> jacobi_order;
  }

  cout << "\n";
  cout << "  The requested Jacobi order of the rule is " << jacobi_order << "\n";
//
//  Get the output option or quadrature file root name:
//
  if ( 3 < argc )
  {
    filename = argv[3];
  }
  else
  {
    cin >> filename;
  }

  pyramid_handle ( legendre_order, jacobi_order, filename );

  cout << "\n";
  cout << "PYRAMID_RULE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void jacobi_compute ( int order, double alpha, double beta, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_COMPUTE computes a Jacobi quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//    Thanks to Xu Xiang of Fudan University for pointing out that
//    an earlier implementation of this routine was incorrect!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    1 <= ORDER.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.  For simple Legendre quadrature,
//    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  double an;
  double *b;
  double bn;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x0;

  if ( order < 1 )
  {
    std::cerr << "\n";
    std::cerr << "JACOBI_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of ORDER = " << order << "\n";
    std::exit ( 1 );
  }

  b = new double[order];
  c = new double[order];
//
//  Check ALPHA and BETA.
//
  if ( alpha <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "JACOBI_COMPUTE - Fatal error!\n";
    std::cerr << "  -1.0 < ALPHA is required.\n";
    std::exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    std::cerr << "\n";
    std::cerr << "JACOBI_COMPUTE - Fatal error!\n";
    std::cerr << "  -1.0 < BETA is required.\n";
    std::exit ( 1 );
  }
//
//  Set the recursion coefficients.
//
  for ( i = 1; i <= order; i++ )
  {
    if ( alpha + beta == 0.0 || beta - alpha == 0.0 )
    {
      b[i-1] = 0.0;
    }
    else
    {
      b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
             ( ( alpha + beta + ( double ) ( 2 * i ) ) 
             * ( alpha + beta + ( double ) ( 2 * i - 2 ) ) );
    }

    if ( i == 1 )
    {
      c[i-1] = 0.0;
    }
    else
    {
      c[i-1] = 4.0 * ( double ) ( i - 1 ) 
         * ( alpha + ( double ) ( i - 1 ) ) 
          * ( beta + ( double ) ( i - 1 ) ) 
            * ( alpha + beta + ( double ) ( i - 1 ) ) / 
            ( ( alpha + beta + ( double ) ( 2 * i - 1 ) ) 
            * std::pow ( alpha + beta + ( double ) ( 2 * i - 2 ), 2 ) 
            * ( alpha + beta + ( double ) ( 2 * i - 3 ) ) );
    }
  }

  delta = r8_gamma ( alpha        + 1.0 ) 
        * r8_gamma (         beta + 1.0 ) 
        / r8_gamma ( alpha + beta + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * std::pow ( 2.0, alpha + beta + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );
      bn = beta / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 1.48 * an + 0.96 * bn 
        + 0.452 * an * an + 0.83 * an * bn;

      x0 = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * beta * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x0 = x0 - r1 * r2 * r3 * ( 1.0 - x0 );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * beta / 
        ( ( 6.28 + beta ) * ( double ) ( order * order ) );

      x0 = x0 - r1 * r2 * r3 * ( x[0] - x0 );
    }
    else if ( i < order - 1 )
    {
      x0 = 3.0 * x[i-2] - 3.0 * x[i-3] + x[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }

    jacobi_root ( &x0, order, alpha, beta, &dp2, &p1, b, c );

    x[i-1] = x0;
    w[i-1] = cc / ( dp2 * p1 );
  }
//
//  Reverse the order of the values.
//
  for ( i = 1; i <= order/2; i++ )
  {
    temp       = x[i-1];
    x[i-1]     = x[order-i];
    x[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp       = w[i-1];
    w[i-1]     = w[order-i];
    w[order-i] = temp;
  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_RECUR evaluates a Jacobi polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of J(ORDER)(X).
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, BETA, the exponents of (1+X) and
//    (1-X) in the quadrature rule.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0 );
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( x - b[i-1] ) *  ( *p1 ) - c[i-1] * p0;
    *dp2 = ( x - b[i-1] ) * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
//****************************************************************************80

void jacobi_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, BETA, the exponents of (1+X) and
//    (1-X) in the quadrature rule.
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    jacobi_recur ( &p2, dp2, p1, *x, order, alpha, beta, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
//****************************************************************************80

void legendre_compute ( int order, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE computes a Legendre quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = 1.0.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    1 <= ORDER.
//
//    Output, double X[ORDER], the abscissas of the rule.
//
//    Output, double W[ORDER], the weights of the rule.
//    The weights are positive, symmetric, and should sum to 2.
//
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    std::cerr << "\n";
    std::cerr << "LEGENDRE_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of ORDER = " << order << "\n";
    std::exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 =  std::cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
//
//  Initial approximation H:
//
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
//
//  Refine H using one step of Newton's method:
//
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    x[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    w[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    x[0] = 0.0;
  }
//
//  Shift the data up.
//
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    x[iback-1] = x[iback-ncopy-1];
    w[iback-1] = w[iback-ncopy-1];
  }
//
//  Reflect values for the negative abscissas.
//
  for ( i = 1; i <= order - nmove; i++ )
  {
    x[i-1] = - x[order-i];
    w[i-1] = w[order-i];
  }

  return;
}
//****************************************************************************80

void pyramid_handle ( int legendre_order, int jacobi_order, string filename )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_HANDLE computes the requested pyramid rule and outputs it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LEGENDRE_ORDER, JACOBI_ORDER, the orders 
//    of the component Legendre and Jacobi rules.
//
//    Input, string FILENAME, the rootname for the files,  
//    write files 'file_w.txt' and 'file_x.txt', and 'file_r.txt', weights,
//    abscissas, and region.
//
{
# define DIM_NUM 3

  string filename_r;
  string filename_w;
  string filename_x;
  int i;
  int j;
  double jacobi_alpha;
  double jacobi_beta;
  double *jacobi_w;
  double *jacobi_x;
  int k;
  int l;
  double *legendre_w;
  double *legendre_x;
  int pyramid_order;
  double pyramid_r[DIM_NUM*5] = {
    -1.0, -1.0, 0.0,
    +1.0, -1.0, 0.0,
    -1.0, +1.0, 0.0,
    +1.0, +1.0, 0.0,
     0.0,  0.0, 1.0 };
  double *pyramid_w;
  double *pyramid_x;
  double volume;
  double wi;
  double wj;
  double wk;
  double xi;
  double xj;
  double xk;
//
//  Compute the factor rules.
//
  legendre_w = new double[legendre_order];
  legendre_x = new double[legendre_order];

  legendre_compute ( legendre_order, legendre_x, legendre_w );

  jacobi_w = new double[jacobi_order];
  jacobi_x = new double[jacobi_order];

  jacobi_alpha = 2.0;
  jacobi_beta = 0.0;

  jacobi_compute ( jacobi_order, jacobi_alpha, jacobi_beta, jacobi_x, jacobi_w );
//
//  Compute the pyramid rule.
//
  pyramid_order = legendre_order * legendre_order * jacobi_order;

  pyramid_w = new double[pyramid_order];
  pyramid_x = new double[DIM_NUM*pyramid_order];

  volume = 4.0 / 3.0;

  l = 0;
  for ( k = 0; k < jacobi_order; k++ )
  {
    xk = ( jacobi_x[k] + 1.0 ) / 2.0;
    wk = jacobi_w[k] / 2.0;
    for ( j = 0; j < legendre_order; j++ )
    {
      xj = legendre_x[j];
      wj = legendre_w[j];
      for ( i = 0; i < legendre_order; i++ )
      {
        xi = legendre_x[i];
        wi = legendre_w[i];
        pyramid_w[l] = wi * wj * wk / 4.0 / volume;
        pyramid_x[0+l*3] = xi * ( 1.0 - xk );
        pyramid_x[1+l*3] = xj * ( 1.0 - xk );
        pyramid_x[2+l*3] =              xk;
        l = l + 1;
      }
    }
  }

  delete [] jacobi_w;
  delete [] jacobi_x;
  delete [] legendre_w;
  delete [] legendre_x;
//
//  Write the rule to files.
//
  filename_w = filename + "_w.txt";
  filename_x = filename + "_x.txt";
  filename_r = filename + "_r.txt";

  cout << "\n";
  cout << "  Creating quadrature files.\n";
  cout << "\n";
  cout << "  \"Root\" file name is   \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Weight file will be   \"" << filename_w << "\".\n";
  cout << "  Abscissa file will be \"" << filename_x << "\".\n";
  cout << "  Region file will be   \"" << filename_r << "\".\n";

  r8mat_write ( filename_w, 1,       pyramid_order, pyramid_w );
  r8mat_write ( filename_x, DIM_NUM, pyramid_order, pyramid_x );
  r8mat_write ( filename_r, DIM_NUM, 5,             pyramid_r );

  delete [] pyramid_w;
  delete [] pyramid_x;

# undef DIM_NUM
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
//    18 February 2008
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
    value = -x;
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
//    18 February 2008
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

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / std::sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * std::log ( y );
      res = std::exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
