# include <cmath>
# include <cstdlib>
# include <ctime>
# include <cstring>
# include <iomanip>
# include <iostream>

using namespace std;

# include "triangle_integrals.hpp"

//*****************************************************************************80

void i4_to_pascal ( int k, int &i, int &j )

//*****************************************************************************80
//
//  Purpose:
//
//    I4_TO_PASCAL converts a linear index to Pascal triangle coordinates.
//
//  Discussion:
//
//    We describe the grid points in Pascal's triangle in two ways:
//
//    As a linear index K:
//
//                     1
//                   2   3
//                 4   5   6
//               7   8   9   10
//
//    As elements (I,J) of Pascal's triangle:
//
//                     0,0
//                  1,0   0,1
//               2,0   1,1    0,2
//            3,0   2,1   1,2    0,3
//
//  Example:
//
//     K  I  J
//
//     1  0  0
//     2  1  0
//     3  0  1
//     4  2  0
//     5  1  1
//     6  0  2
//     7  3  0
//     8  2  1
//     9  1  2
//    10  0  3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the linear index of the (I,J) element.
//    1 <= K.
//
//    Output, int &I, &J, the Pascal indices.
//
{
  int d;

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "I4_TO_PASCAL - Fatal error!\n";
    cerr << "  K must be positive.\n";
    exit ( 1 );
  }

  d = i4_to_pascal_degree ( k );

  j = k - ( d * ( d + 1 ) ) / 2 - 1;
  i = d - j;

  return;
}
//****************************************************************************80

int i4_to_pascal_degree ( int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_PASCAL_DEGREE converts a linear index to a Pascal triangle degree.
//
//  Discussion:
//
//    We describe the grid points in Pascal's triangle in two ways:
//
//    As a linear index K:
//
//                     1
//                   2   3
//                 4   5   6
//               7   8   9   10
//
//    As elements (I,J) of Pascal's triangle:
//
//                     0,0
//                  1,0   0,1
//               2,0   1,1    0,2
//            3,0   2,1   1,2    0,3
//
//    The quantity D represents the "degree" of the corresponding monomial,
//    that is, D = I + J.
//
//    We can compute D directly from K using the quadratic formula.
//
//  Example:
//
//     K  I  J  D
//
//     1  0  0  0
//
//     2  1  0  1
//     3  0  1  1
//
//     4  2  0  2
//     5  1  1  2
//     6  0  2  2
//
//     7  3  0  3
//     8  2  1  3
//     9  1  2  3
//    10  0  3  3
//
//    11  4  0  4
//    12  3  1  4
//    13  2  2  4
//    14  1  3  4
//    15  0  4  4
//
//    16  5  0  5
//    17  4  1  5
//    18  3  2  5
//    19  2  3  5
//    20  1  4  5
//    21  0  5  5
//
//    22  6  0  6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the linear index of the (I,J) element.
//    1 <= K.
//
//     Output, int I4_TO_PASCAL_DEGREE, the degree (sum) of the corresponding 
//     Pascal indices.
//
{
  double arg;
  int d;

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "I4_TO_PASCAL_DEGREE - Fatal error!\n";
    cerr << "  K must be positive.\n";
    exit ( 1 );
  }

  arg = ( double ) ( 1 + 8 * ( k - 1 ) );

  d = ( int ) ( 0.5 * ( -1.0 + sqrt ( arg ) ) );

  return d;
}
//****************************************************************************80

int pascal_to_i4 ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    PASCAL_TO_I4 converts Pacal triangle coordinates to a linear index.
//
//  Discussion:
//
//    We describe the grid points in a Pascal triangle in two ways:
//
//    As a linear index K:
//
//                     1
//                   2   3
//                 4   5   6
//               7   8   9   10
//
//    As elements (I,J) of Pascal's triangle:
//
//                     0,0
//                  1,0   0,1
//               2,0   1,1    0,2
//            3,0   2,1   1,2    0,3
//
//  Example:
//
//     K  I  J
//
//     1  0  0
//     2  1  0
//     3  0  1
//     4  2  0
//     5  1  1
//     6  0  2
//     7  3  0
//     8  2  1
//     9  1  2
//    10  0  3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the row and column indices.  I and J 
//    must be nonnegative.
//
//    Output, int PASCAL_TO_I4, the linear index of the (I,J) element.
//
{
  int d;
  int k;

  if ( i < 0 )
  {
    cerr << "\n";
    cerr << "PASCAL_TO_I4 - Fatal error!\n";
    cerr << "  I < 0.\n";
    cerr << "  I = " << i << "\n";
    exit ( 1 );
  }
  else if ( j < 0 )
  {
    cerr << "\n";
    cerr << "PASCAL_TO_I4 - Fatal error!\n";
    cerr << "  J < 0.\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }

  d = i + j;

  k = ( d * ( d + 1 ) ) / 2 + j + 1;

  return k;
}
//*****************************************************************************80

double *poly_power ( int d1, double p1[], int n )

//*****************************************************************************80
//
//  Purpose:
//
//    POLY_POWER computes a power of a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D1, the degree of the polynomial.
//
//    Input, double P1(M1), the polynomial coefficients.
//    M1 = ((D1+1)*(D1+2))/2.
//
//    Input, int N, the nonnegative integer power.
//
//    Output, double POLY_POWER[M2], the polynomial power.
//    D2 = N * D1.
//    M2 = ((D2+1)*(D2+2))/2.
//
{
  int d2;
  int d3;
  int i;
  int m2;
  double *p2;
  double *p3;
//
//  Set D2 to 0, to indicate that P2 currently contains only
//  a constant term.
//
  d2 = 0;
  m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  p2 = new double[m2];
  p2[0] = 1.0;
//
//  Iterate N times:
//    P3 <= P1 * P2
//    P2 <= P3
//
  for ( i = 1; i <= n; i++ )
  {
    d3 = d1 + d2;
    p3 = poly_product ( d1, p1, d2, p2 );

    delete [] p2;

    d2 = d3;
    p2 = p3;
  }

  return p2;
}
//*****************************************************************************80

double *poly_power_linear ( int d1, double p1[], int n )

//*****************************************************************************80
//
//  Purpose:
//
//    POLY_POWER_LINEAR computes the polynomial ( a + b*x + c*y ) ^ n.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D1, the degree of the linear polynomial,
//    which should be 1 (or possibly 0).
//
//    Input, double P1(M1), the coefficients of the linear polynomial.
//    M1 = ( (D1+1)*(D1+2) ) / 2, which should be 3.
//
//    Input, int N, the power to which the polynomial is to be 
//    raised.  0 <= N.
//
//    Output, double P2(M2), the coefficients of the power polynomial.
//    D2 = N * D1;
//    M2 = ( (D2+1)*(D2+2) ) / 2.
//
{
  int d2;
  int i;
  int j;
  int k;
  int l;
  int m2;
  double *p2;

  if ( d1 < 0 )
  {
    cerr << "\n";
    cerr << "POLY_POWER_LINEAR - Fatal error!\n";
    cerr << "  D1 < 0.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "POLY_POWER_LINEAR - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  d2 = n * d1;
  m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  p2 = new double[m2];

  if ( d1 == 0 )
  {
    p2[0] = pow ( p1[0], n );
    return p2;
  }

  if ( n == 0 )
  {
    p2[0] = 1.0;
    return p2;
  }
//
//  Use the Trinomial formula.
//
  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n - i; j++ )
    {
      for ( k = 0; k <= n - i - j; k++ )
      {
//
//  We store X^J Y^K in location L.
//
        l = pascal_to_i4 ( j, k );
        p2[l-1] = ( double ) ( trinomial ( i, j, k ) ) 
          * pow ( p1[0], i ) * pow ( p1[1], j ) * pow ( p1[2], k );
      }
    }
  }

  return p2;
}
//*****************************************************************************80

void poly_print ( int d, double p[], string title )

//*****************************************************************************80
//
//  Purpose:
//
//    POLY_PRINT prints an XY polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the degree of the polynomial.
//
//    Output, double P[M], the coefficients of all monomials of 
//    degree 0 through D.  P must contain ((D+1)*(D+2))/2 entries.
//
//    Input, string TITLE, a title string.
//
{
  int i;
  int j;
  int k;
  int km1;
  int m;

  m = ( ( d + 1 ) * ( d + 2 ) ) / 2;

  for ( km1 = 0; km1 < m; km1++ )
  {
    if ( p[km1] != 0.0 )
    {
      break;
    }
    cout << title << " = 0\n";
    return;
  }

  cout << title << "\n";

  for ( km1 = 0; km1 < m; km1++ )
  {
    k = km1 + 1;
    i4_to_pascal ( k, i, j );

    if ( p[km1] != 0.0 )
    {

      if ( p[km1] < 0.0 )
      {
        cout << "  -" <<  fabs ( p[km1] );
      }
      else
      {
        cout << "  +" << p[km1];
      }

      if ( i + j != 0 )
      {
        cout << " ";
      }

      if ( i == 0 )
      {
      }
      else if ( i == 1 )
      {
        cout << "x";
      }
      else
      {
        cout << "x^" << i;
      }

      if ( j == 0 )
      {
      }
      else if ( j == 1 )
      {
        cout << "y";
      }
      else
      {
        cout << "y^" << j;
      }
      cout << "\n";
    }
  }

  return;
}
//*****************************************************************************80

double *poly_product ( int d1, double p1[], int d2, double p2[] )

//*****************************************************************************80
//
//  Purpose:
//
//    POLY_PRODUCT computes P3(x,y) = P1(x,y) * P2(x,y) for polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D1, the degree of factor 1.
//
//    Input, double P1[M1], the factor 1 coefficients.
//    M1 = ((D1+1)*(D1+2))/2.
//
//    Input, int D2, the degree of factor 2.
//
//    Input, double P2[M2], the factor2 coefficients.
//    M2 = ((D2+1)*(D2+2))/2.
//
//    Output, double POLY_PRODUCT[M3], the result coefficients.
//    D3 = D1 + D2;
//    M3 = ((D3+1)*(D3+2))/2.
//
{
  int d3;
  int i1;
  int i2;
  int i3;
  int j1;
  int j2;
  int j3;
  int k1;
  int k1m1;
  int k2;
  int k2m1;
  int k3;
  int k3m1;
  int m1;
  int m2;
  int m3;
  double *p3;

  m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
//
//  Consider each entry in P1:
//    P1(K1) * X^I1 * Y^J1
//  and multiply it by each entry in P2:
//    P2(K2) * X^I2 * Y^J2
//  getting 
//    P3(K3) = P3(K3) + P1(K1) * P2(X2) * X^(I1+I2) * Y(J1+J2)
//
  d3 = d1 + d2;
  m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
  p3 = new double[m3];

  for ( k3m1 = 0; k3m1 < m3; k3m1++ )
  {
    p3[k3m1] = 0.0;
  }

  for ( k1m1 = 0; k1m1 < m1; k1m1++ )
  {
    k1 = k1m1 + 1;
    i4_to_pascal ( k1, i1, j1 );
    for ( k2m1 = 0; k2m1 < m2; k2m1++ )
    {
      k2 = k2m1 + 1;
      i4_to_pascal ( k2, i2, j2 );
      i3 = i1 + i2;
      j3 = j1 + j2;
      k3 = pascal_to_i4 ( i3, j3 );
      k3m1 = k3 - 1;
      p3[k3m1] = p3[k3m1] + p1[k1m1] * p2[k2m1];
    }
  }

  return p3;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

void rs_to_xy_map ( double t[], double &a, double &b, double &c, double &d, 
  double &e, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    RS_TO_XY_MAP returns the linear map from reference to physical triangle.
//
//  Discussion:
//
//    This function returns the coefficients of the linear map that sends
//    the vertices of the reference triangle, (0,0), (1,0) and (0,1), to
//    the vertices of a physical triangle T, of the form:
//
//      X = A + B * R + C * S;
//      Y = D + E * R + F * S.
//
//  Reference Element:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2,3], the coordinates of the vertices.  The
//    vertices are assumed to be the images of (0,0), (1,0) and (0,1) 
//    respectively.
//
//    Output, double &A, &B, &C, &D, &E, &F, the mapping coefficients.
//
{
  a = t[0+0*2];
  b = t[0+1*2] - t[0+0*2];
  c = t[0+2*2] - t[0+0*2];

  d = t[1+0*2];
  e = t[1+1*2] - t[1+0*2];
  f = t[1+2*2] - t[1+0*2];

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

double triangle_area ( double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA returns the area of a triangle.
//
//  Discussion:
//
//    If the vertices are given in counter clockwise order, the area
//    will be positive.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double value;

  value = 0.5 * 
    ( 
        ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
      - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
    );

  return value;
}
//****************************************************************************80

double triangle_monomial_integral ( int i, int j, double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_MONOMIAL_INTEGRAL integrates a monomial over an arbitrary triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the exponents of X and Y in the monomial.
//    0 <= I, J.
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_MONOMIAL_INTEGRAL, the integral 
//    of X^I * Y^J over triangle T.
//
{
  double a;
  double b;
  double c;
  double d;
  int d1;
  int d2;
  int d3;
  int d4;
  int d5;
  double e;
  double f;
  int m1; 
  int m2;
  double *p1;
  double *p2;
  double *p3;
  double *p4;
  double *p5;
  double q;
//
//  Get map coefficients from reference RS triangle to general XY triangle.
//    R = a+b*X+c*Y
//    S = d+e*X+f*Y
//
  rs_to_xy_map ( t, a, b, c, d, e, f );
//
//  Set
//    P1(R,S) = a+b*R+c*S
//    P2(R,S) = d+e*R+f*S
//
  d1 = 1;
  m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  p1 = new double[m1];
  p1[0] = a;
  p1[1] = b;
  p1[2] = c;

  d2 = 1;
  m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  p2 = new double[m2];
  p2[0] = d;
  p2[1] = e;
  p2[2] = f;
//
//  Exponentiate:
//    P3(R,S) = P1(R,S)^i
//    P4(R,S) = P2(R,S)^j
//
  d3 = i * d1;
  p3 = poly_power_linear ( d1, p1, i );

  d4 = j * d2;
  p4 = poly_power_linear ( d2, p2, j );
//
//  Compute the product 
//    P5(R,S) = P3(R,S) * P4(R,S)
//
  d5 = d3 + d4;
  p5 = poly_product ( d3, p3, d4, p4 );
//
//  Compute the integral of P5(R,S) over the reference triangle.
//
  q = triangle01_poly_integral ( d5, p5 );
//
//  Multiply by the area of the physical triangle T(X,Y) divided by
//  the area of the reference triangle.
//
  q = q * triangle_area ( t ) / 0.5;

  delete [] p1;
  delete [] p2;
  delete [] p3;
  delete [] p4;
  delete [] p5;

  return q;
}
//****************************************************************************80

double triangle_poly_integral ( int d, double p[], double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_POLY_INTEGRAL: polynomial integral over a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the degree of the polynomial.
//
//    Input, double P[M], the polynomial coefficients.
//    M = ((D+1)*(D+2))/2.
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_POLY_INTEGRAL, the integral.
//
{
  int i;
  int j;
  int k;
  int km1;
  int m;
  double q;

  m = ( ( d + 1 ) * ( d + 2 ) ) / 2;

  q = 0.0;
  for ( km1 = 0; km1 <= m; km1++ )
  {
    k = km1 + 1;
    i4_to_pascal ( k, i, j );
    q = q + p[km1] * triangle_monomial_integral ( i, j, t );
  }

  return q;
}
//****************************************************************************80

double triangle_xy_integral ( double x1, double y1, double x2, double y2, 
  double x3, double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_XY_INTEGRAL computes the integral of XY over a triangle.
//
//  Discussion:
//
//    This function was written as a special test case for the general
//    problem of integrating a monomial x^alpha * y^beta over a general 
//    triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    triangle vertices.
//
//    Output, double TRIANGLE_XY_INTEGRAL, the integral of X*Y 
//    over the triangle.
//
{
  double det;
  double p00;
  double p01;
  double p02;
  double p10;
  double p11;
  double p20;
  double q;
//
//  x = x1 * ( 1.0 - xi - eta )
//    + x2 *         xi
//    + x3 *              eta;
//
//  y = y1 * ( 1.0 - xi - eta )
//    + y2 *         xi
//    + y3 *              eta;
//
//  Rewrite as linear polynomials in (xi,eta):
//
//  x = x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta
//  y = y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta
//
//  Jacobian:
//
//    J = [ ( x2 - x1 )  ( x3 - x1 ) ]
//        [ ( y2 - y1 )  ( y3 - y1 ) ]
//
//    det J = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
//
//  Integrand
//
//    x * y = ( x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta )
//          * ( y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta )
//
//  Rewrite as linear combination of monomials:
//
//    x * y = 1      * x1 * y1
//          + eta    * ( x1 * ( y3 - y1 ) + ( x3 - x1 ) * y1 )
//          + xi     * ( x1 * ( y2 - y1 ) + ( x2 - x1 ) * y1 )
//          + eta^2  * ( x3 - x1 ) * ( y3 - y1 )
//          + xi*eta * ( ( x2 - x1 ) * ( y3 - y1 ) + ( x3 - x1 ) * ( y2 - y1 ) )
//          + xi^2   * ( x2 - x1 ) * ( y2 - y1 )
//
  det = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 );

  p00 = x1 * y1;

  p01 = x1 * ( y3 - y1 ) + ( x3 - x1 ) * y1;
  p10 = x1 * ( y2 - y1 ) + ( x2 - x1 ) * y1;

  p02 = ( x3 - x1 ) * ( y3 - y1 );
  p11 = ( x2 - x1 ) * ( y3 - y1 ) + ( x3 - x1 ) * ( y2 - y1 );
  p20 = ( x2 - x1 ) * ( y2 - y1 );

  q = 0.0;
  q = q + p00 * triangle01_monomial_integral ( 0, 0 );
  q = q + p10 * triangle01_monomial_integral ( 1, 0 );
  q = q + p01 * triangle01_monomial_integral ( 0, 1 );
  q = q + p20 * triangle01_monomial_integral ( 2, 0 );
  q = q + p11 * triangle01_monomial_integral ( 1, 1 );
  q = q + p02 * triangle01_monomial_integral ( 0, 2 );

  q = q * det;

  return q;
}
//****************************************************************************80

double triangle01_monomial_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_MONOMIAL_INTEGRAL: monomial integrals in the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the exponents.  
//    Each exponent must be nonnegative.
//
//    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the integral.
//
{
  int k;
  int l;
  double q;

  k = 0;
  q = 1.0;

  for ( l = 1; l <= i; l++ )
  {
    k = k + 1;
    q = q * ( double ) ( l ) / ( double ) ( k );
  }

  for ( l = 1; l <= j; l++ )
  {
    k = k + 1;
    q = q * ( double ) ( l ) / ( double ) ( k );
  }

  for ( l = 1; l <= 2; l++ )
  {
    k = k + 1;
    q = q / ( double ) ( k );
  }

  return q;
}
//****************************************************************************80

double triangle01_poly_integral ( int d, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_POLY_INTEGRAL: polynomial integral over the unit triangle.
//
//  Discussion:
//
//    The unit triangle is T = ( (0,0), (1,0), (0,1) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer D, the degree of the polynomial.
//
//    Input, double P[M], the polynomial coefficients.
//    M = ((D+1)*(D+2))/2.
//
//    Output, double TRIANGLE01_POLY_INTEGRAL, the integral.
//
{
  int i;
  int j;
  int k;
  int km1;
  int m;
  double q;

  m = ( ( d + 1 ) * ( d + 2 ) ) / 2;

  q = 0.0;
  for ( km1 = 0; km1 < m; km1++ )
  {
    k = km1 + 1;
    i4_to_pascal ( k, i, j );
    q = q + p[km1] * triangle01_monomial_integral ( i, j );
  }

  return q;
}
//****************************************************************************80

int trinomial ( int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    TRINOMIAL computes a trinomial coefficient.
//
//  Discussion:
//
//    The trinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where I objects are of type 1, J of type 2, and K of type 3.
//    and N = I + J + K.
//
//    T(I,J,K) = N! / ( I! J! K! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, K, the factors.
//    All should be nonnegative.
//
//    Output, int TRINOMIAL, the trinomial coefficient.
//
{
  int l;
  int t;
  int value;
//
//  Each factor must be nonnegative.
//
  if ( i < 0 || j < 0 || k < 0 )
  {
    cerr << "\n";
    cerr << "TRINOMIAL - Fatal error!\n";
    cerr << "  Negative factor encountered.\n";
    exit ( 1 );
  }

  value = 1;

  t = 1;

  for ( l = 1; l <= i; l++ )
  {
//
//  value = value * t / l;
//
    t = t + 1;
  }

  for ( l = 1; l <= j; l++ )
  {
    value = value * t / l;
    t = t + 1;
  }

  for ( l = 1; l <= k; l++ )
  {
    value = value * t / l;
    t = t + 1;
  }
  
  return value;
}
//****************************************************************************80

void xy_to_rs_map ( double t[], double &a, double &b, double &c, double &d, 
  double &e, double &f )

//****************************************************************************80
//
//  Purpose:
//
//    XY_TO_RS_MAP returns the linear map from physical to reference triangle.
//
//  Discussion:
//
//    Given the vertices T of an arbitrary triangle in the (X,Y) coordinate
//    system, this function returns the coefficients of the linear map
//    that sends the vertices of T to (0,0), (1,0) and (0,1) respectively
//    in the reference triangle with coordinates (R,S):
//
//      R = A + B * X + C * Y;
//      S = D + E * X + F * Y.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the X and Y coordinates
//    of the vertices.  The vertices are assumed to be the images of
//    (0,0), (1,0) and (0,1) respectively.
//
//    Output, double &A, &B, &C, &D, &E, &F, the mapping coefficients.
//
{
  double g;

  g =    ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] )   
         - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) );

  a = ( - ( t[1+2*2] - t[1+0*2] ) * t[0+0*2]  
        + ( t[0+2*2] - t[0+0*2] ) * t[1+0*2] ) / g;

  b =     ( t[1+2*2] - t[1+0*2] ) / g;

  c =   - ( t[0+2*2] - t[0+0*2] ) / g;

  d = (   ( t[1+1*2] - t[1+0*2] ) * t[0+0*2] 
        - ( t[0+1*2] - t[0+0*2] ) * t[1+0*2] ) / g;

  e =   - ( t[1+1*2] - t[1+0*2] ) / g;

  f =     ( t[0+1*2] - t[0+0*2] ) / g;

  return;
}
