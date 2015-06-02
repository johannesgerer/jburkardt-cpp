# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "tetrahedron_arbq_rule.hpp"

//****************************************************************************80

double *kjacoypols3 ( double x, double y, double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    KJACOYPOLS3 evaluates modified Jacobi polynomials.
//
//  Discussion:
//
//    This procedure evaluates Jacobi polynomials multiplied by
//    specific polynomials given by the formula
//      P_n^{(a,b)} (x/y) y^n
//    at the user-provided point x/y.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, Y, define the evaluation point X/Y.
//
//    Input, double A, B, the parameters.
//
//    Input, int N, the highest degree to compute.
//
//    Output, double KJACOYPOLS3[N+1], the polynomial values.
//
{
  double alpha;
  double beta;
  int k;
  double k_r8;
  double pk;
  double pkm1;
  double pkp1;
  double *pols;

  pols = new double[n+1];

  pkp1 = 1.0;
  pols[0] = pkp1;

  if ( n == 0 )
  {
    return pols;
  }

  pk = pkp1;
  pkp1 = 0.5 * ( ( a - b ) * y + ( 2.0 + a + b ) * x );
  pols[1] = pkp1;

  if ( n == 1 )
  {
    return pols;
  }

  for ( k = 2; k <= n; k++ )
  {
    k_r8 = ( double ) ( k );

    alpha = ( 2.0 * k_r8 + a + b - 1.0 ) 
      * ( a - b ) * ( a + b ) * y 
      + ( 2.0 * k_r8 + a + b - 1.0 ) 
      * ( 2.0 * k_r8 + a + b - 2.0 ) 
      * ( 2.0 * k_r8 + a + b ) * x;

    beta = 2.0 * ( k_r8 + a - 1.0 ) 
      * ( k_r8 + b - 1.0 ) 
      * ( 2.0 * k_r8 + a + b ) * y * y;

    pkm1 = pk;
    pk = pkp1;
    pkp1 = ( alpha * pk - beta * pkm1 ) 
      / ( 2.0 * k_r8 * ( k_r8 + a + b ) 
      * ( 2.0 * k_r8 + a + b - 2.0 ) );

    pols[k] = pkp1;
  }

  return pols;
}
//****************************************************************************80

double *klegeypols ( double x, double y, int n )

//****************************************************************************80
//
//  Purpose:
//
//    KLEGEYPOLS evaluates scaled Legendre polynomials.
//
//  Discussion:
//
//    This routine evaluate a sequence of scaled Legendre polynomials
//    P_n(x/y) y^n, with the parameter y in [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double Y, the parameter.
//
//    Input, int N, the highest degree to be evaluated.
//
//    Output, double KLEGEYPOLS[N+1], the polynomial values.
//
{
  int k;
  double pk;
  double pkm1;
  double pkp1;
  double *pols;

  pols = new double[n+1];

  pkp1 = 1.0;
  pols[0] = pkp1;
  if ( n == 0 )
  {
    return pols;
  }

  pk = pkp1;
  pkp1 = x;
  pols[1] = pkp1;
  if ( n == 1 )
  {
    return pols;
  }

  for ( k = 1; k < n; k++ )
  {
    pkm1 = pk;
    pk = pkp1;
    pkp1 = ( ( 2.0 * k + 1.0 ) * x * pk 
      - k * pkm1 * y * y ) / ( k + 1.0 );
    pols[k+1] = pkp1;
  }

  return pols;
}
//****************************************************************************80

double *ortho3eva ( int degree, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    ORTHO3EVA evaluates polynomials orthogonal in the reference triangle.
//
//  Discussion:
//
//    This procedure evaluates the Koornwinder's orthogonal polynomial
//    up to order DEGREE on the reference tetrahedron with vertices
//      (-1, -1/Sqrt(3), -1/Sqrt(6)),
//      ( 0,  2/Sqrt(3), -1/Sqrt(6)),
//      ( 1, -1/Sqrt(3), -1/Sqrt(6)),
//      ( 0,  0,      3/Sqrt(6))
//
//    The polynomials are ordered by their order, and in each order,
//    they are ordered lexicographically in terms of their indices
//    (m,n,k).
//
//    The total number of polynomials should be
//      NVALS = ( ( DEGREE + 1 ) * ( DEGREE + 2 ) * ( DEGREE + 3 ) ) / 6.
//
//    The calculations are based on Koornwinder's representation
//    of the orthogonal polynomials on the right triangle
//      (-1,-1,-1), (-1,1,-1), (1,-1,-1),(-1,-1,1)
//    given by:
//      K_mnk(x,y,z) =
//        P_m ((2x+2+y+z)/(-y-z)) * ((-y-z)/2)^m *
//        P_n^{2m+1,0}((2y+1+z)/(1-z)) * ((1-z)/2)^{n}
//        P_k^{2m+2n+2,0} (z)
//    with the input (x,y,z) transformed as
//      x = -1/2 + xold -yold/s3 - zold/s6
//      y = -1/2 + 2/s3 * yold - zold/s6
//      z = -1/2 + s6/2 * zold
//    where
//      s3=sqrt(3)
//      s6=sqrt(6)
//    and
//      P_m is the Legendre polynomial of degree m,
//      P_n^{2m+1,0} the Jacobi polynomial of degree n, parameters 2*m+1 and 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the maximum degree.
//
//    Input, double XYZ[3], the evaluation point.
//
//    Output, double ORTHO3EVA[NVALS], the polynomial values.
//
{
  double *f;
  double *f1;
  double *f2s;
  double *f3s;
  double *fvals;
  int i;
  int j;
  int k;
  int m;
  int mmax;
  int n;
  int ncount;
  int nvals;
  double p1;
  double p2;
  double scale;
  double *uvw;
  double x;
  double x1;
  double y;
  double y1;
  double z;
//
//  Convert coordinates from reference to Koornwinder tetrahedron.
//
  uvw = ref_to_koorn ( xyz );
//
//  Compute F1.
//
  p1 = 0.5 * ( 2.0 * uvw[0] + 2.0 + uvw[1] + uvw[2] );
  p2 = - 0.5 * ( uvw[1] + uvw[2] );

  f1 = klegeypols ( p1, p2, degree );
//
//  Compute F2S.
//
  f2s = new double[(degree+1)*(degree+1)];

  for ( j = 1; j <= degree + 1; j++ )
  {
    for ( i = 1; i <= degree + 1; i++ )
    {
      f2s[i-1+(j-1)*(degree+1)] = 0.0;
    }
  }

  for ( m = 0; m <= degree; m++ )
  {
    x1 = 0.5 * ( 2.0 * uvw[1] + 1.0 + uvw[2] );
    y1 = 0.5 * ( 1.0 - uvw[2] );
    p1 = ( double ) ( 2 * m + 1 );
    p2 = 0.0;

    f = kjacoypols3 ( x1, y1, p1, p2, degree - m );

    for ( i = 1; i <= degree + 1 - m; i++ )
    {
      f2s[i-1+m*(degree+1)] = f[i-1];
    }
    delete [] f;
  }
//
//  Compute F3S.
//
  f3s = new double[(degree+1)*(degree+1)];

  for ( j = 1; j <= degree + 1; j++ )
  {
    for ( i = 1; i <= degree + 1; i++ )
    {
      f3s[i-1+(j-1)*(degree+1)] = 0.0;
    }
  }

  x1 = uvw[2];
  y1 = 1.0;
  p2 = 0.0;

  for ( j = 1; j <= degree + 1; j++ )
  {
    p1 = ( double ) ( 2 * j );

    f = kjacoypols3 ( x1, y1, p1, p2, degree + 1 - j );

    for ( i = 1; i <= degree + 2 - j; i++ )
    {
      f3s[i-1+(j-1)*(degree+1)] = f[i-1];
    }
    delete [] f;
  }
//
//  Construct FVALS.
//
  nvals = ( ( degree + 1 ) * ( degree + 2 ) * ( degree + 3 ) ) / 6;
  fvals = new double[nvals];

  ncount = 0;

  for ( mmax = 0; mmax <= degree; mmax++ )
  {
    for ( m = 0; m <= mmax; m++ )
    {
      for ( n = 0; n <= mmax - m; n++ )
      {
        k = mmax - m - n;

        scale = sqrt 
        ( 
          4.0 
          / ( double ) ( 2 * mmax + 3 ) 
          / ( double ) ( 2 * m + 1 ) 
          / ( double ) ( n + m + 1 ) 
          / sqrt ( 2.0 ) 
        );

        fvals[ncount] = 
          f1[m] * 
          f2s[n+m*(degree+1)] * 
          f3s[k+(m+n)*(degree+1)] / scale;

        ncount = ncount + 1;
      }
    }
  }

  delete [] f1;
  delete [] f2s;
  delete [] f3s;
  delete [] uvw;

  return fvals;
}
//****************************************************************************80

void r8mat_row_copy ( int m, int n, int i, double v[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ROW_COPY copies a vector into a row of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Input, int I, the index of the row.
//    0 <= I <= M-1.
//
//    Input, double V[N], the row to be copied.
//
//    Input/output, double A[M*N], the matrix into which
//    the row is to be copied.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    a[i+j*m] = v[j];
  }
  return;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

double *ref_to_koorn ( double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    REF_TO_KOORN maps points from the reference to Koornwinder's tetrahedron.
//
//  Discussion:
//
//    The reference tetrahedron has vertices:
//
//      (-1, -1/Sqrt(3), -1/Sqrt(6) )
//      ( 0,  2/Sqrt(3), -1/Sqrt(6) )
//      ( 1, -1/Sqrt(3), -1/Sqrt(6) )
//      ( 0,  0,      3/Sqrt(6) )
//
//    Koornwinder's tetrahedron has vertices:
//
//      ( -1, -1, -1 )
//      ( -1, +1, -1 )
//      ( +1, -1, -1 )
//      ( -1, -1, +1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R[3], the coordinates of a point in the
//    reference tetrahedron.
//
//    Output, double REF_TO_KOORN[3], the coordinates of the point in
//    the Koornwinder tetrahedron.
//
{
  double a10;
  double a11;
  double a12;
  double a13;
  double a20;
  double a21;
  double a22;
  double a23;
  double a30;
  double a31;
  double a32;
  double a33;
  double s3;
  double s6;
  double *u;

  s3 = sqrt ( 3.0 );
  s6 = sqrt ( 6.0 );

  a10 = - 0.5;
  a11 =   1.0;
  a12 = - 1.0 / s3;
  a13 = - 1.0 / s6;

  a20 = - 0.5;
  a21 =   0.0;
  a22 =   2.0 / s3;
  a23 = - 1.0 / s6;

  a30 = - 0.5;
  a31 =   0.0;
  a32 =   0.0;
  a33 =   0.5 * s6;

  u = new double[3];

  u[0] = a10 + a11 * r[0] + a12 * r[1] + a13 * r[2];
  u[1] = a20 +              a22 * r[1] + a23 * r[2];
  u[2] = a30 +                           a33 * r[2];

  return u;
}
//****************************************************************************80

void rule01 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE01 returns the rule of degree 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
       0.00000000000000000 };
  double ys[] = { 
       0.00000000000000000 };
  double zs[] = { 
       0.00000000000000000 };
  double ws[] = { 
       0.9709835434146467 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule02 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE02 returns the rule of degree 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.2677028585910073,0.1510931841533142, 
    -.1367699195237390,0.8067449309051964 };
  double ys[] = { 
    0.5939017006199488E-01,0.4354087732476309, 
    -.3280341115590410,-.3400977314285288 };
  double zs[] = { 
    -.3969426941142150,0.2151719254306919, 
    0.2950846519937133,-.3430378951002550 };
  double ws[] = { 
    0.2918008865477151,0.2706884392356724, 
    0.3098349601753470,0.9865925745591227E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule03 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE03 returns the rule of degree 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.1685037180276000,0.2783799427534418E-01, 
    -.3512216177343445,0.4308532463549043, 
    -.4676763747967377,0.1492831253848431 };
  double ys[] = { 
    0.1910914916271708,-.2304932838839657E-01, 
    0.1835144026339993,-.2474715823180446, 
    -.4235250827264375,0.6397847685164516 };
  double zs[] = { 
    -.3896267314585163,0.5481350663241830, 
    0.5147815330343534E-01,-.1683315641007033, 
    -.1586973077889307,-.1080219253055393 };
  double ws[] = { 
    0.1287213727402025,0.2179034339695993, 
    0.1243503792062836,0.2446917182410072, 
    0.1365439875826512,0.1187726516749031 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule04 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE04 returns the rule of degree 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.1459612280979987E-01,-.2937244896309993, 
    0.3230878337158274,-.8819654179396583E-01, 
    0.1270203663314710,-.4414307688091184, 
    0.3939418589633429,-.6034469714614210E-01, 
    -.8914059834601185E-01,0.6545213033603760, 
    -.7307642259677857 };
  double ys[] = { 
    0.2426564579199708E-02,0.1764396506764613, 
    -.1932070410956556,0.2317884270105980E-01, 
    0.5451410677215219,-.3848225631590180, 
    0.2068744670639530,-.4573739074927080, 
    0.7268092659599459,-.3970356243870571, 
    -.2669420088648982 };
  double zs[] = { 
    0.7093565780103633,0.1108860494941134, 
    0.1695426396704650,-.2011819391325586, 
    0.1309503990759315,0.9225429679162532E-01, 
    -.3111560426198242,-.3302215329322376, 
    -.3507931737363739,-.2970424833951137, 
    -.3861037570241846 };
  double ws[] = { 
    0.1033787090646894,0.1070356256090164, 
    0.1504792582740940,0.1877987156186583, 
    0.7395274312521298E-01,0.7712199925411270E-01, 
    0.6745419368701999E-01,0.5819413648173244E-01, 
    0.5378646711152148E-01,0.5366953183949744E-01, 
    0.3811216334909158E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule05 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE05 returns the rule of degree 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.1664180721577685E-15,0.4089925917487008, 
    0.6338467275429999E-16,-.4089925917487006, 
    0.3723900931035782E-16,-.2435436770532023, 
    0.7718657600705262E-17,0.2435436770532026, 
    0.6290589987564350,-.6290589987564351, 
    -.1407768953799327E-15,0.4089925917487006, 
    -.1047965211464895E-15,-.4089925917487007 };
  double ys[] = { 
    0.1094362324146226E-15,-.2361319829426750, 
    0.4722639658853500,-.2361319829426750, 
    0.3548001489301807E-16,0.1406100075060977, 
    -.2812200150121953,0.1406100075060977, 
    -.3631873822681842,-.3631873822681841, 
    0.7263747645363684,0.2361319829426750, 
    -.4722639658853500,0.2361319829426750 };
  double zs[] = { 
    0.7704367825296720,0.3339410527875835, 
    0.3339410527875834,0.3339410527875834, 
    -.2982788694307592,0.9942628981025299E-01, 
    0.9942628981025306E-01,0.9942628981025299E-01, 
    -.2568122608432240,-.2568122608432239, 
    -.2568122608432239,-.3339410527875835, 
    -.3339410527875834,-.3339410527875835 };
  double ws[] = { 
    0.7136053542145054E-01,0.4131148601232373E-01, 
    0.4131148601232375E-01,0.4131148601232376E-01, 
    0.1094181214137256,0.1094181214137255, 
    0.1094181214137255,0.1094181214137255, 
    0.7136053542145050E-01,0.7136053542145052E-01, 
    0.7136053542145050E-01,0.4131148601232370E-01, 
    0.4131148601232373E-01,0.4131148601232375E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule06 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE06 returns the rule of degree 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.4919941951523113E-02,0.3415058060219544E-01, 
    -.4561198328617235,0.2144499229769493, 
    -.1540693059995943,0.3064525403222947, 
    -.1305882959514773E-01,0.1535244689838072E-01, 
    -.1567784487775340,-.3973616940506757, 
    0.9764884616649552E-02,-.5228087984647453, 
    0.5565605702210328,0.8423081979023464E-01, 
    0.2174166091160236,0.3212158945225647, 
    -.2593648078545398,0.5636454767202584, 
    -.6410949644625729,-.3344151090582155, 
    -.1894588386892019E-01,-.8214655896199351, 
    0.8597415475615655 };
  double ys[] = { 
    -.1392238503430922E-01,0.2542073830543402, 
    -.3110385447445375,-.1234903873249722, 
    -.6865214472302168E-01,0.1818495501923178, 
    -.2488891169073396E-01,-.3261143100532925, 
    0.2522454301000923,-.1977159083378496, 
    0.6504184291133676,-.6635055304123420E-01, 
    -.3186690054220115,0.2138059822864521E-01, 
    0.4599507004997997,-.4260052421545051, 
    0.3851980211137577,-.1111606227821696, 
    -.1934318338507612,-.4491758699723651, 
    0.9010599459930961,-.5109860359355179, 
    -.5208085141225744 };
  double zs[] = { 
    0.1066228263967720E+01,0.6315886967039234, 
    0.2300267655106985,0.5632351636302827, 
    0.6206264896282621,0.1254685876442911, 
    -.3357705926282663,0.1147788985927003, 
    0.2140597800705522,-.7963361802882329E-01, 
    0.5868446281578688E-02,0.1934272688455952E-02, 
    0.1379740251424839E-01,0.3195756126912550E-01, 
    -.3027232303753162,-.2994855602604263, 
    -.2825197099253122,-.3068140495300409, 
    -.3660473605161587,-.3102423292746128, 
    -.3473952050436902,-.2635645240101362, 
    -.3603432585369335 };
  double ws[] = { 
    0.6889900028407776E-02,0.3059511845881479E-01, 
    0.2352818079375775E-01,0.4747145913756433E-01, 
    0.5149448476741100E-01,0.4198464552012798E-01, 
    0.6518777271610311E-01,0.5808554412413398E-01, 
    0.6079762507769643E-01,0.6246727661754427E-01, 
    0.4232555916776466E-01,0.2191089157324707E-01, 
    0.4526058614388953E-01,0.1093321687620292, 
    0.5237000773398501E-01,0.4832452079715138E-01, 
    0.6123639953387314E-01,0.4577996169297793E-01, 
    0.2476427914920903E-01,0.3856173525882120E-01, 
    0.1546337820963114E-01,0.1012088745519476E-01, 
    0.7031160695311595E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule07 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE07 returns the rule of degree 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.3082522913562920,0.3230068564765168, 
    0.1790030943099486E-02,-.1749858562850217, 
    0.3515818224512229E-01,-.1148912802505062, 
    0.1080242619580473,0.9960886019209524E-01, 
    0.2552589657173040,0.3018972093964544, 
    -.3637995251662868E-01,-.1979748432674351, 
    0.1705547608587194,0.5382714575357732, 
    -.5826255951559661,-.6638827484255140E-01, 
    0.1394865450483365,-.2525448211561483, 
    -.3364179479356783,0.5088174095091343, 
    -.3997676771299581,0.2233579208997628, 
    -.2887689832327979,0.2430056805988458, 
    -.2402459184668794,-.5737122495168941, 
    0.5530740586345179,0.7632160871540358, 
    0.1108323525316276E-02,-.7666589725519221, 
    0.3314379268437029 };
  double ys[] = { 
    -.1978432879678884,-.1166608144800289, 
    0.2760223041671345E-02,0.1492553071272583, 
    0.3222570009831960,0.3598150246997987E-01, 
    -.2123211082745265,0.3370600346165786, 
    0.2345931387652712E-01,-.2797440812728529, 
    0.6784158147792628,-.1099260120706853, 
    0.5854481149716676E-01,-.3326202438156060, 
    -.2866394288669600,0.3396643732745171, 
    -.3078502318368667,-.9687707852320762E-01, 
    -.4222809305623556,-.4607960014652947E-01, 
    0.1694117934001565,0.5070263922318456, 
    -.4320411612241725,0.5262155488082486, 
    0.4919239659783042,-.7453182334351910E-01, 
    -.3355026530466742E-02,-.4345394288375978, 
    0.8835815093231754,-.4540327424996255, 
    -.5419034033722774 };
  double zs[] = { 
    0.6553771405341420,0.5963217268212818, 
    0.9191670342435462,0.5162921387398244, 
    0.5492103599470217,-.3886917653966460, 
    0.5138214657068572,-.2732502226402171, 
    -.1289266426427488,-.3368155028155599, 
    0.4965154590499942E-01,0.3688101367310880, 
    0.3245398450336227,0.6530998463389483E-01, 
    0.4050073689194288E-01,0.4339137621154699E-01, 
    -.8554188228126980E-02,-.1212234343128167, 
    0.2244936992133664E-01,-.5389155149780182E-01, 
    -.1331769110761944E-01,-.3996173060097081E-01, 
    -.3304370757138854,-.3633697548034391, 
    -.3184670075538235,-.3348000741414298, 
    -.3507444035401421,-.3120479891575616, 
    -.3155336884007630,-.3148702231578014, 
    -.3104410268762152 };
  double ws[] = { 
    0.7484517057961263E-02,0.1165548467333077E-01, 
    0.2254313978794985E-01,0.2176106698437745E-01, 
    0.2261516504544270E-01,0.2653689704536702E-01, 
    0.2796071748247393E-01,0.4174371077929897E-01, 
    0.5125012874361469E-01,0.4062452567161350E-01, 
    0.2011109093086575E-01,0.5855177928202636E-01, 
    0.6477548431925077E-01,0.3049968783939012E-01, 
    0.2749125083938057E-01,0.6480294592357103E-01, 
    0.6232630627722486E-01,0.7994900688083093E-01, 
    0.2026705321864424E-01,0.2384971520555609E-01, 
    0.2170019222676944E-01,0.2153672490681497E-01, 
    0.3586406280991152E-01,0.1621736817489311E-01, 
    0.3281522904442925E-01,0.2724694235165297E-01, 
    0.2048072851756261E-01,0.1999288967511468E-01, 
    0.1864561424137827E-01,0.1753312507928042E-01, 
    0.1215099239866789E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule08 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE08 returns the rule of degree 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.1611772551670525E-01,-.7148827236840712E-02, 
    0.4589923372088558,-.4520990492494178, 
    0.1109560356430631,0.1825224293648890, 
    -.1507489805419971,-.2877324734878600, 
    0.3680747212378487,0.4788156424890513, 
    -.3478779367400819,-.3548285868734744E-01, 
    0.2591650442934043,-.1305037378946309E-02, 
    -.3493996070828597E-01,-.1660462859707297, 
    0.1010085522422429,-.3587468163605144, 
    -.3473067563985682,-.8635741121401184E-01, 
    0.2901168201778487,0.1042093355150105, 
    -.6419387451918364E-01,-.2016986972118968, 
    0.1846170677356399,0.3788893958160157, 
    -.4818716163165796,0.6821719050778170, 
    0.3978642293015137,-.6975889470159129, 
    0.2742964749336314E-02,0.1645687860050991E-01, 
    -.4618349214785903,0.1575090568397094, 
    0.5399175203215238,-.1403905262089725, 
    0.4453396570037857E-01,-.4234360302104348, 
    -.7214176376894034,0.7278118716190342, 
    0.4901890779701343,-.8758887352507779, 
    0.7632377205619759E-02,0.9317431147169187 };
  double ys[] = { 
    0.7998600119500696E-03,0.1723914595007140, 
    -.1744716351651229,-.2632941629071404, 
    -.6295720376005641E-01,0.1678382818593359, 
    -.7826941873249307E-01,-.1234686974783635, 
    -.2169876766339741,-.4629437300594242E-01, 
    -.1295630336164733,-.2268108031307481, 
    0.1330071443027294,0.4919409115894537, 
    0.3001995169122194,0.1407736457687509, 
    -.2259257575315161,-.1083199795573463, 
    -.3742304739030127,-.1478062482085008, 
    -.3449512373726749,0.2883671141680468E-01, 
    0.2558125570853043,0.4404368079041183, 
    0.3998568605285963,-.9956585162815586E-01, 
    -.1509747383909747E-01,-.3956518618817587, 
    0.2256616013248351,-.3899651023918449, 
    0.8134998230879619,-.4451036700294634, 
    -.4535647479937944,0.6650145460171435, 
    -.4674624469973033,0.7012653983759632, 
    -.5107241847309127,0.2263324997687613, 
    -.2727155170636352,-.2861627803070728, 
    0.1856306523332632,-.5566374678887420, 
    0.1056056946669113E+01,-.5664247124705913 };
  double zs[] = { 
    0.1140750120480647E+01,0.7860696072793454, 
    -.3991430912330994,0.3536092960240441, 
    0.8189404142301465,0.4764962716361752, 
    0.7628974295081231,-.3880058323817796, 
    0.3997553695274408,0.1144656464995089, 
    0.2509171968778639,0.4174615578343315, 
    -.2747778345118032,0.3031929324753210, 
    -.3442503524230918,0.3970537434180833, 
    -.3070893258601431,-.1989904396903917, 
    0.1526258958568332E-01,-.1542350189843949E-03, 
    0.4121026001314489E-01,0.3248913848279931, 
    -.6957900515388780E-01,-.3340987564879216E-01, 
    -.1303527431545757E-01,-.8390322512212595E-01, 
    -.4501404031339341E-01,-.7451523522744341E-01, 
    -.3993939022009145,-.1254210458205579, 
    -.1415815795205319,-.1858244299410800, 
    -.3245512424809017,-.3317121001508232, 
    -.3168159505701121,-.3372150239150091, 
    -.3950013670511707,-.3336444006035868, 
    -.3519392970001664,-.3226167582915808, 
    -.2780042325707565,-.3720577799104590, 
    -.3972751008974220,-.3821073616790632 };
  double ws[] = { 
    0.1694544493019982E-02,0.1462667227182682E-01, 
    0.9589935274093377E-02,0.9835015847019813E-02, 
    0.1710148279355958E-01,0.1516325692442914E-01, 
    0.1746691605511999E-01,0.1604974314640067E-01, 
    0.2185156320663785E-01,0.1242799062066869E-01, 
    0.2646409426919389E-01,0.3135597218974045E-01, 
    0.2736783716876112E-01,0.2345735439005822E-01, 
    0.3130572641601596E-01,0.3548953658742431E-01, 
    0.4036134913447077E-01,0.3792924852161072E-01, 
    0.2440963887102594E-01,0.5132903383025206E-01, 
    0.3177764520649842E-01,0.5417358749775593E-01, 
    0.5180566415546370E-01,0.2633068360731991E-01, 
    0.3893527990282211E-01,0.5027716272897634E-01, 
    0.2208377845212325E-01,0.1429278762793107E-01, 
    0.6697289387082959E-02,0.1648266825522313E-01, 
    0.1688462952130244E-01,0.2190464958161271E-01, 
    0.2281706804403004E-01,0.1956212468616595E-01, 
    0.1994889218942623E-01,0.1737944394368884E-01, 
    0.7302913219204132E-02,0.2128832734454982E-01, 
    0.1333720911018169E-01,0.1494192725418561E-01, 
    0.1152112823254964E-01,0.2472998490996246E-02, 
    0.2018756899942699E-02,0.1470016064285112E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule09 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE09 returns the rule of degree 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.2108052017103265,0.1283729355850052E-02, 
    -.2918734175018473E-01,-.4357949297689266, 
    -.1799341637870203,0.6942222186204158E-01, 
    -.1005642991135230,0.1904855696753355, 
    0.3231238851996552,0.6819987943340138E-01, 
    0.9708131623556143E-01,0.4246239788055203, 
    -.3569664256203109,-.2145752061141752, 
    -.4818674305781609E-02,-.4314954104621454, 
    -.1312453567417025,-.4083558222913988, 
    -.1245612605849638,0.1136348780330888, 
    0.3351668846328467,0.2159589333207287, 
    -.1228639493897767,0.1380276243456064, 
    0.3905319783823730,-.4580867844324543, 
    0.2310228564389770,0.3554068884610396, 
    -.6505241861040961E-01,-.1786506541288197, 
    0.1265104521534247E-01,-.4573827823867556, 
    -.5542731199877091E-01,-.6387630000068671E-01, 
    -.3413749732146179,0.5966077890976830, 
    0.6917314710518101,0.1078283022003015E-01, 
    0.1833983579629060,-.1566832100689738, 
    -.7215475195156269,0.5634418100399829, 
    -.6038175317729408,0.4424617155479634, 
    0.5376954508544742,0.5801674951413821, 
    0.1418634443888370,0.3802052805108170E-01, 
    -.6766546338276799,-.1483779428740810, 
    0.7588538984987838,-.4009332704109223, 
    0.4294022281773516,-.4911820218387012, 
    -.8230436696819400,0.8307385272442451, 
    0.2847756354434797E-02 };
  double ys[] = { 
    -.8955776048075967E-01,0.6121197754911472E-03, 
    0.2122758579659528,-.2839355480318341, 
    -.9582556152029258E-01,-.1394666680707507, 
    0.4836397528261835,0.3248843951738821, 
    -.3091097002134723,-.1799541017656358, 
    0.1346690341273634,-.1974037451364399, 
    0.5240665094885540E-01,0.2761387382865308, 
    0.4840041246351735,-.3972752286396523, 
    0.7569076187744307E-01,-.1239859742013766, 
    -.2344790300437711,-.3305476987034425E-01, 
    0.1770622072827525E-01,0.2385246764606510, 
    0.5313763678851683E-01,0.7344859267512056E-02, 
    -.2434286634053894,-.2589556174266634, 
    -.2578369917972710,0.1458366018963618, 
    -.3194064610603558,-.8361295659204318E-01, 
    0.2789606323941811,-.3620544381079852, 
    0.4583494650640864,-.4110512725485398, 
    0.1453586652323175,-.2817335623108552, 
    -.4244951895917183,0.7740686042096789, 
    0.5096501737807261,0.5871627777302516, 
    -.4069288134308984,-.8745287448038948E-01, 
    -.1605011166676197,-.4658949840050826, 
    -.4822919051612324,-.2630487459758545E-01, 
    0.7093134008050596,-.4981319573109797, 
    -.1608102554139131,0.7476578114817578, 
    -.2863864956290275,0.3184909341350718, 
    0.3342451266830045,-.5259906817121476, 
    -.4752298758892802,-.4894484556416516, 
    0.9800293546928647 };
  double zs[] = { 
    0.7735882074261043,0.1020940836140086E+01, 
    0.7799358723388078,0.3644377643783332, 
    0.7572044348058194,0.7055487604983020, 
    -.4071603611745789,0.2784337835342537, 
    0.3220003691351323,-.3955135683176328, 
    0.6632209417507691,0.3671897521183273, 
    -.3747481356463467,0.2325932263186099, 
    0.3487840306818015,-.4035704541972099, 
    0.5374900196303547,0.2654318104765763, 
    0.3698883029355781,0.4369873036895401, 
    0.1860393220233105,-.3479911186004724, 
    -.2731450419753432,-.8193864872695196E-01, 
    -.2812924318728223,-.2546313790752990, 
    0.8678030145840356E-01,-.1428438498378993, 
    -.2156622116049685,0.9674903232992661E-01, 
    0.1583036737187768,-.2371578880066309E-01, 
    -.2083640565350588,-.1327703130967625E-01, 
    -.8664488303113425E-01,-.7748272430107629E-01, 
    -.9669501947122619E-02,-.6316745396431407E-01, 
    -.1041997605600258,-.7558004547861685E-01, 
    -.7999199041063003E-01,-.3751069250486347, 
    -.1299227031705554,-.1892728310331580, 
    -.3910612496171095,-.1907765273821742, 
    -.3477134529420335,-.3474543531762018, 
    -.3606479588180986,-.3440499686491684, 
    -.3474563411684222,-.3276892436683218, 
    -.3542380917718721,-.3040019129455090, 
    -.3440705957046394,-.3298160634602457, 
    -.3378012113084769 };
  double ws[] = { 
    0.5680408778529046E-02,0.6887854441113948E-02, 
    0.6899595536569939E-02,0.6990873936739353E-02, 
    0.1142122940842178E-01,0.1220808590880294E-01, 
    0.7651487209643779E-02,0.9540263940869291E-02, 
    0.9495638963926394E-02,0.1385248437188017E-01, 
    0.1446333231867884E-01,0.1128763818550772E-01, 
    0.1290073048331596E-01,0.1468008559737106E-01, 
    0.1526131467815127E-01,0.7960414659472547E-02, 
    0.2551792590874883E-01,0.1931170779707098E-01, 
    0.2561026230172263E-01,0.3169262013385435E-01, 
    0.2455872440888530E-01,0.2465571070153197E-01, 
    0.3419839251013865E-01,0.4066027350436995E-01, 
    0.3092126671144214E-01,0.2731106642703576E-01, 
    0.3711204543388613E-01,0.2490755526192437E-01, 
    0.3425568258765496E-01,0.4767735489615401E-01, 
    0.4389950894493362E-01,0.2234183537671572E-01, 
    0.3468016547781724E-01,0.1599733331684861E-01, 
    0.3216001820424372E-01,0.2173130208768448E-01, 
    0.4172519829863862E-02,0.1261067900439839E-01, 
    0.2073467359947484E-01,0.1259876220642243E-01, 
    0.6842300320404830E-02,0.1001228308171065E-01, 
    0.1328037050111305E-01,0.1493735138014892E-01, 
    0.6932861063756816E-02,0.6827617690151055E-02, 
    0.1270616397760185E-01,0.1288875403806118E-01, 
    0.8899252930026689E-02,0.8881169613715286E-02, 
    0.5288619835392034E-02,0.1284718318260967E-01, 
    0.6409598847098365E-02,0.8055508398015583E-02, 
    0.7788463427997768E-02,0.6545788868966830E-02, 
    0.5341431206058898E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule10 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE10 returns the rule of degree 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.4135441926108436,0.6920344574148628E-02, 
    -.4204645371849727,0.1212827019050737, 
    0.2550746585402268,-.3763573604446380, 
    -.4793578815415788E-03,-.9545817212686222E-01, 
    0.9593753000785579E-01,0.3126301401176562, 
    -.6029042538437952,0.2902741137267026, 
    0.1793124962974703,-.3341633308670676, 
    0.1548508345696874,0.4006666604011494, 
    -.3875674230187136,-.1309923738259358E-01, 
    -.1159749898125322,-.6294076737424348, 
    0.7453826635553922,-.1465165887848432, 
    -.4643141868299223,0.6108307756147643, 
    0.1526773226986329,-.6092465369941313, 
    0.4565692142952937,-.2349836288652037, 
    0.5557206177695652,-.3207369889038079, 
    0.2617190792164227,-.2353754650695125, 
    -.2634361414712093E-01,-.1518707842414898, 
    -.2491539125596758,0.4010246968008866, 
    0.2108318943729519,0.7968056612120857E-02, 
    -.2187999509850842,-.6451845341149960, 
    0.5779406907458116,0.6724384336957467E-01, 
    -.1117027003111241,-.2717471055583957, 
    0.3834498058697393,0.2353864631064914, 
    0.1649339445101810,-.4003204076167762, 
    0.7813866474350551,0.1310349087707534E-02, 
    -.7826969965232643,-.1446486485942199, 
    -.1795841794748884E-02,0.1464444903888344, 
    0.5662897411515679,0.1839447497576211E-02, 
    -.5681291886496517,0.1959958820407995E-02, 
    -.3041697942355038,0.3022098354145521, 
    -.7138061627291874E-01,-.3075192230121023, 
    0.3788998392856530,-.3187848524903898E-01, 
    -.9081742825948123,0.9400527678433584, 
    0.9976676831438748E-01,-.8701982408315725, 
    0.7704314725177221,-.8207757778140335E-14, 
    -.3478417597234231E-12,0.4790952016774938E-12, 
    0.1374787712765062E-12,0.1420358457394610E-12 };
  double ys[] = { 
    -.2467507765313092,0.4815151646539115, 
    -.2347643881227681,-.3645574461360758, 
    0.2873126239566789,0.7724482217897545E-01, 
    0.1105023601492690,-.5566631617759855E-01, 
    -.5483604397160680E-01,0.5156767709512212, 
    0.1290725785546668E-01,-.5285840288063388, 
    0.2823324600522004,0.1412294698335291E-01, 
    -.2964554070359728,0.2161993077800489, 
    0.2388878524669672,-.4550871602470489, 
    0.7937355713187059,-.4973050730407310, 
    -.2964304982776649,0.6207352334853357, 
    -.4372547047061620,-.1834805287791172, 
    0.6153490109135903,-.1754520654179484, 
    -.4398969454955785,-.5060230351463111, 
    0.4950972550220373E-01,0.4565133096445042, 
    0.1206845953990400,0.1663130735568127, 
    -.2869976689562073,0.3753807951356096, 
    -.3192143548135662,-.5616644032198326E-01, 
    -.1309245702301312,0.2480480615696827, 
    -.1171234913398124,-.2948509623067830, 
    -.4113207155193597,0.7061716778256595, 
    0.3782781131850991,-.2858764327334718, 
    -.9240168045211911E-01,-.3263497523611192, 
    0.3670255329377183,-.4067578057659137E-01, 
    -.4526468520348789,0.9030241128745099, 
    -.4503772608395989,0.8558659569093126E-01, 
    -.1680627021511816,0.8247610646027727E-01, 
    -.3290715455089481,0.6549570744945620, 
    -.3258855289855604,0.3500934424096372, 
    -.1733493470760727,-.1767440953334651, 
    0.3963042303943790,-.2599695422276034, 
    -.1363346881672204,0.1067074385094486E+01, 
    -.5611447706066169,-.5059296144870963, 
    0.9472180066426262,-.3872084475068061, 
    -.5600095591352807,-.1466739111569232E-12, 
    0.2419854376387837E-13,-.1552488340937713E-12, 
    -.2306323605219641E-12,-.1706174368865759E-12 };
  double zs[] = { 
    -.2246678712539802E-01,-.2246678712556187E-01, 
    -.2246678712549978E-01,-.3674985430672734, 
    -.3674985430672704,-.3674985430673924, 
    0.9168884763080956,0.9168884763077307, 
    0.9168884763085783,-.3553336558198409, 
    -.3553336558198190,-.3553336558198354, 
    -.1983165808272536,-.1983165808275630, 
    -.1983165808272446,-.1447858405450569, 
    -.1447858405451866,-.1447858405451966, 
    -.3542393281108072,-.3542393281107521, 
    -.3542393281107588,-.1251511927365907, 
    -.1251511927364737,-.1251511927364298, 
    -.1387177797398434,-.1387177797398834, 
    -.1387177797396668,-.3549137865746493, 
    -.3549137865746083,-.3549137865746373, 
    0.8659130596210637E-01,0.8659130596171444E-01, 
    0.8659130596207916E-01,0.2597687530320102, 
    0.2597687530321350,0.2597687530322305, 
    0.2630632081142404,0.2630632081143893, 
    0.2630632081144756,-.3453833089184027, 
    -.3453833089183551,-.3453833089183527, 
    -.2327889534298980,-.2327889534297419, 
    -.2327889534296802,0.2415268120578684, 
    0.2415268120576211,0.2415268120576665, 
    -.2133068003220094,-.2133068003224423, 
    -.2133068003227821,0.6055816506735499, 
    0.6055816506736809,0.6055816506737702, 
    0.1414131684170811,0.1414131684165572, 
    0.1414131684162387,0.5593538346997974, 
    0.5593538346994973,0.5593538347003276, 
    -.3961007879024273,-.3961007879023746, 
    -.3961007879023505,-.3950741090480416, 
    -.3950741090482457,-.3950741090476022, 
    -.3779492208572051,-.3779492208572193, 
    -.3779492208571862,-.3464548867786856E-01, 
    0.1201009856985286E+01,0.4068123865912653, 
    -.3341407786078385,0.5699475205900819 };
  double ws[] = { 
    0.2007080596188233E-01,0.2007080596192337E-01, 
    0.2007080596192014E-01,0.1202850982705908E-01, 
    0.1202850982706234E-01,0.1202850982703673E-01, 
    0.6428246923028596E-02,0.6428246923032827E-02, 
    0.6428246923025550E-02,0.5961017737373873E-02, 
    0.5961017737376429E-02,0.5961017737354793E-02, 
    0.2466985253776778E-01,0.2466985253776116E-01, 
    0.2466985253775616E-01,0.1220343371669821E-01, 
    0.1220343371668563E-01,0.1220343371666686E-01, 
    0.6782282187245643E-02,0.6782282187250343E-02, 
    0.6782282187251558E-02,0.1276699476009991E-01, 
    0.1276699476009263E-01,0.1276699476009060E-01, 
    0.1330898468041284E-01,0.1330898468040929E-01, 
    0.1330898468040468E-01,0.8872828138103311E-02, 
    0.8872828138107657E-02,0.8872828138112581E-02, 
    0.2604852195750169E-01,0.2604852195750782E-01, 
    0.2604852195751094E-01,0.1091598279364878E-01, 
    0.1091598279363667E-01,0.1091598279364397E-01, 
    0.2727741583607240E-01,0.2727741583609616E-01, 
    0.2727741583610115E-01,0.1214067234472188E-01, 
    0.1214067234470734E-01,0.1214067234471387E-01, 
    0.2653047387010774E-01,0.2653047387011107E-01, 
    0.2653047387012045E-01,0.1147407678956960E-01, 
    0.1147407678960903E-01,0.1147407678957783E-01, 
    0.6941339438094368E-02,0.6941339438087489E-02, 
    0.6941339438089641E-02,0.1618088676127407E-01, 
    0.1618088676125095E-01,0.1618088676125473E-01, 
    0.9879744049167588E-02,0.9879744049174469E-02, 
    0.9879744049169482E-02,0.1106165676903961E-01, 
    0.1106165676903509E-01,0.1106165676904102E-01, 
    0.7973762810133187E-02,0.7973762810143904E-02, 
    0.7973762810145713E-02,0.7671732380845806E-03, 
    0.7671732380905287E-03,0.7671732380951896E-03, 
    0.1643453575966134E-02,0.1643453575957693E-02, 
    0.1643453575955291E-02,0.4522697832344839E-01, 
    0.3708344580879497E-03,0.1141401437229523E-01, 
    0.2710167289616961E-01,0.1108569325544364E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule11 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE11 returns the rule of degree 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.9651177776899700E-01,-.7093619453993043, 
    0.6128501676297941,-.4527836424396549E-02, 
    -.4612438490926584,0.4657716855189367, 
    0.2438655531438615E-01,0.3220723301923039, 
    -.3464588855081770,-.1267500664585161, 
    -.5368057828378728,0.6635558492970384, 
    -.9547972451775022E-01,-.2357077639614438, 
    0.3311874884792527,-.7908328086591814E-02, 
    -.7854922635641037,0.7934005916498822, 
    -.1494644432070149,-.3155198764172914, 
    0.4649843196230657,0.1165276860785583, 
    -.4038862538014747,0.2873585677227137, 
    -.2047755067090988E-02,-.9205200026937792E-01, 
    0.9409975533966558E-01,0.3468354574885807, 
    0.1682368211989309,-.5150722786890125, 
    -.1144319727418180,-.7020417246489494, 
    0.8164736973908203,-.3048702567593000, 
    -.1613129802547456,0.4661832370146523, 
    0.2101640308259509,0.3270815374994656, 
    -.5372455683250684,-.6404070599046534E-01, 
    -.5801912598472644,0.6442319658370544, 
    0.5470009647103604E-02,-.2517469475214085, 
    0.2462769378777770,0.1008461981380965E-01, 
    0.1731441764929287,-.1832287963083161, 
    0.1748982078445799,-.3794913490849509, 
    0.2045931412407059,0.1065391638487668E-01, 
    -.3452291078932529,0.3345751915058239, 
    -.7083686826885625E-03,0.2101573677291210, 
    -.2094489990466886,-.2851273048890965, 
    -.3267460563855006,0.6118733612746845, 
    -.2323209396920720E-02,-.2164769352136706, 
    0.2188001446096009,-.2051577904029304E-01, 
    -.5570838828147304,0.5775996618543481, 
    -.3026455196997997E-02,0.1457949533626576, 
    -.1427684981666913,-.1973538549140939, 
    -.1582348302473636,0.3555886851614904, 
    0.1068857744362264,-.4126086236419754, 
    0.3057228492051779,-.3003364082942412, 
    0.3997857603256970,-.9944935203203158E-01, 
    0.1221835246093986E-01,-.8799955089249785, 
    0.8677771564639284,0.1534568838650121, 
    -.7253511736488428,0.5718942897840599, 
    0.2235456216473980E-02,-.5785135622840735, 
    0.5762781060658508,0.9436416592156352E-01, 
    0.3995572440244000,-.4939214099453045, 
    0.1578951817317316E-11,-.7234508899977211E-12, 
    -.1195010503052243E-12,0.1273193864088492E-12, 
    -.1398313045262104E-12 };
  double ys[] = { 
    0.7633795193833395,-.2981081083795290, 
    -.4652714110037153,0.5352126684505152, 
    -.2715275555928822,-.2636851128574337, 
    -.3859766773459927,0.2141077150870582, 
    0.1718689622588539,0.6930291114382893, 
    -.4562833332036821,-.2367457782343889, 
    0.3272971265990525,-.2463364302786462, 
    -.8096069632058686E-01,0.9115742149797984, 
    -.4626359205144959,-.4489382944649009, 
    0.4506243076880464,-.3547521586231100, 
    -.9587214906505584E-01,0.3990903837832518, 
    -.9862925550286258E-01,-.3004611282798994, 
    0.1074747662091034,-.5551079101315480E-01, 
    -.5196397519770712E-01,-.3945086927598004, 
    0.4976226634994518,-.1031139707392109, 
    0.8767152876835527,-.5374586392410802, 
    -.3392566484421098,0.3622851099556234, 
    -.4451679421897679,0.8288283223322354E-01, 
    -.4990194872156122,0.4315171332646398, 
    0.6750235395125617E-01,0.7069210789517217, 
    -.4089214177391241,-.2979996612119608, 
    0.2875342242963981,-.1390299448345922, 
    -.1485042794638735,-.2057520317788134, 
    0.1116095528368255,0.9414247894240614E-01, 
    0.3372213377187678,-.1714437779017452E-01, 
    -.3200769599289869,0.3924851952555710, 
    -.1870160353880812,-.2054691598668353, 
    -.2422598488114232,0.1205164591315570, 
    0.1217433896809788,0.5419121734527502, 
    -.5178835760728716,-.2402859737978143E-01, 
    0.2513073392081394,-.1276656279597885, 
    -.1236417112471668,0.6551098499602489, 
    -.3453221108076170,-.3097877391525884, 
    -.1666021864185573,0.8068010612598300E-01, 
    0.8592208029345208E-01,0.2966561449375769, 
    -.3192415243590572,0.2258537942164062E-01, 
    0.4147288692159332,-.1147986386430089, 
    -.2999302305730746,-.2882335265091744, 
    -.1159821960086673,0.4042157225193962, 
    0.1009077018844530E+01,-.4939571057987255, 
    -.5151199130457724,0.7489650175178515, 
    -.2415849489462737,-.5073800685715635, 
    0.6667192805807411,-.3314236784178127, 
    -.3352956021620059,-.5158501413847826, 
    0.3396468355868865,0.1762033057983631, 
    -.8069526265224892E-12,0.4596851666550035E-12, 
    0.1156254386801148E-12,0.5765654131403360E-13, 
    -.1226198587930237E-12 };
  double zs[] = { 
    -.1868268863122697,-.1868268863123053, 
    -.1868268863112705,0.4194136456710284, 
    0.4194136456717518,0.4194136456711013, 
    0.5030086583798776E-01,0.5030086583729960E-01, 
    0.5030086583707338E-01,-.1505956145194099, 
    -.1505956145190243,-.1505956145192429, 
    -.1725275785423899,-.1725275785429845, 
    -.1725275785414913,-.1709810615420194, 
    -.1709810615421794,-.1709810615408042, 
    0.7674742839401762E-01,0.7674742839470597E-01, 
    0.7674742839681360E-01,-.3720134569923326, 
    -.3720134569922893,-.3720134569921418, 
    0.9342960665019454,0.9342960665021878, 
    0.9342960664953203,-.3611736289802844E-01, 
    -.3611736289860155E-01,-.3611736289902880E-01, 
    -.3535560367475518,-.3535560367476202, 
    -.3535560367476884,-.1926114352197085, 
    -.1926114352198034,-.1926114352197086, 
    -.2818523519641055,-.2818523519643090, 
    -.2818523519648298,-.3888552761464975, 
    -.3888552761464523,-.3888552761454053, 
    0.6654848063127317,0.6654848063130780, 
    0.6654848063078891,0.1794683181875785E-01, 
    0.1794683181975163E-01,0.1794683181866871E-01, 
    -.2272794453431112,-.2272794453432979, 
    -.2272794453429204,0.8078480103364902E-01, 
    0.8078480103319574E-01,0.8078480103608091E-01, 
    0.3406111942551650,0.3406111942551089, 
    0.3406111942537021,-.3549131236009490, 
    -.3549131236011042,-.3549131236012491, 
    0.4213057676587053,0.4213057676584742, 
    0.4213057676619541,-.2430250873757710, 
    -.2430250873755122,-.2430250873736784, 
    0.6939507005767215,0.6939507005768448, 
    0.6939507005757009,-.3547895907184612, 
    -.3547895907186004,-.3547895907184947, 
    0.3273634204692249,0.3273634204691543, 
    0.3273634204697847,0.3797479097987934, 
    0.3797479098019154,0.3797479097984070, 
    -.3626672020688875,-.3626672020689380, 
    -.3626672020686431,-.3656623860600263, 
    -.3656623860600038,-.3656623860598407, 
    0.1056953360682411,0.1056953360682777, 
    0.1056953360707847,-.4015122352990265, 
    -.4015122352990430,-.4015122352995484, 
    0.1118223609071020E+01,0.6664104583675564, 
    -.3809990858052303,-.2253945502682283, 
    0.2604940999755755 };
  double ws[] = { 
    0.6188005592967251E-02,0.6188005592972490E-02, 
    0.6188005592984317E-02,0.1421063031043561E-02, 
    0.1421063031047919E-02,0.1421063030985486E-02, 
    0.9481646089825704E-02,0.9481646089870772E-02, 
    0.9481646089801222E-02,0.6814905118771046E-02, 
    0.6814905118757536E-02,0.6814905118753749E-02, 
    0.1912041674313553E-01,0.1912041674309802E-01, 
    0.1912041674324984E-01,0.3076029506690844E-02, 
    0.3076029506690773E-02,0.3076029506710983E-02, 
    0.1135836724916717E-01,0.1135836724921729E-01, 
    0.1135836724932070E-01,0.1091603321391938E-01, 
    0.1091603321391454E-01,0.1091603321395138E-01, 
    0.4213691001392970E-02,0.4213691001384827E-02, 
    0.4213691001487200E-02,0.1206837849134008E-01, 
    0.1206837849134312E-01,0.1206837849132945E-01, 
    0.2578861292324762E-02,0.2578861292345221E-02, 
    0.2578861292334968E-02,0.1252456101442731E-01, 
    0.1252456101442267E-01,0.1252456101444618E-01, 
    0.8063913155156468E-02,0.8063913155154395E-02, 
    0.8063913155182456E-02,0.4993280626450522E-02, 
    0.4993280626459499E-02,0.4993280626580781E-02, 
    0.6150188384809636E-02,0.6150188384813883E-02, 
    0.6150188384677927E-02,0.2678388417077401E-01, 
    0.2678388417077927E-01,0.2678388417078649E-01, 
    0.1946856127752551E-01,0.1946856127748248E-01, 
    0.1946856127749281E-01,0.2102094189170934E-01, 
    0.2102094189169974E-01,0.2102094189174428E-01, 
    0.1685367818736497E-01,0.1685367818729912E-01, 
    0.1685367818732997E-01,0.5900869705091857E-02, 
    0.5900869705093673E-02,0.5900869705073673E-02, 
    0.1690563324458130E-01,0.1690563324459602E-01, 
    0.1690563324456019E-01,0.1476341187397782E-01, 
    0.1476341187398414E-01,0.1476341187397479E-01, 
    0.6992794755775152E-02,0.6992794755747075E-02, 
    0.6992794755746150E-02,0.1513056742691520E-01, 
    0.1513056742688911E-01,0.1513056742692685E-01, 
    0.7015883109811540E-02,0.7015883109809210E-02, 
    0.7015883109746114E-02,0.4965174122051565E-02, 
    0.4965174121903291E-02,0.4965174122041231E-02, 
    0.2712699778890213E-02,0.2712699778885687E-02, 
    0.2712699778899985E-02,0.5900718944404856E-02, 
    0.5900718944405426E-02,0.5900718944424988E-02, 
    0.8768888887154097E-02,0.8768888887160882E-02, 
    0.8768888887196669E-02,0.3175432879836244E-02, 
    0.3175432879833069E-02,0.3175432879794895E-02, 
    0.1078986495318477E-02,0.1747404076054822E-01, 
    0.1129982873905144E-01,0.2890093028298477E-01, 
    0.2624431483486055E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule12 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE12 returns the rule of degree 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.3958544098912666,-.4870511331351648, 
    0.9119672324389515E-01,-.5926431804946591, 
    0.5193516147881150,0.7329156570653904E-01, 
    -.4544141850245140E-01,-.2756993864806933, 
    0.3211408049831440,-.3413860630862161, 
    0.5248742662396945,-.1834882031534880, 
    -.1319722431139201,0.2771078804366659, 
    -.1451356373227354,0.2916719156941026, 
    -.6416905801467393,0.3500186644526196, 
    -.4070439485536645E-01,-.4815006651936328, 
    0.5222050600490048,-.9820856549819251E-01, 
    -.4276075601779065,0.5258161256761047, 
    0.2184089033681613,-.4656593625800790, 
    0.2472504592119225,-.6606777491144762E-02, 
    0.1271987824323769,-.1205920049412379, 
    -.1862346674942773,-.1097662119993451, 
    0.2960008794936306,-.5890402754221185, 
    0.6424505517737594,-.5341027635165645E-01, 
    -.1196400210307819E-01,-.7808240524677493, 
    0.7927880545708264,0.6235239461686299E-01, 
    -.4408318818916021,0.3784794872747496, 
    0.2235903863093564,-.5240787479715297E-01, 
    -.1711825115122079,0.1666477050331823, 
    -.2729955606600447,0.1063478556268736, 
    -.7415049139037595E-01,0.7389518464069887E-01, 
    0.2553067496801190E-03,0.1955908717228039, 
    0.1141488778114551,-.3097397495342445, 
    -.1802680258591541E-01,-.2161777731993448, 
    0.2342045757852541,-.1532613811702923, 
    -.1800985367097719,0.3333599178800782, 
    0.5776336665233524E-02,-.6762886678803212, 
    0.6705123312150877,-.8304590790913337E-01, 
    0.3554044991859694,-.2723585912768360, 
    0.5423845315967933,-.2526088517557438, 
    -.2897756798410528,0.5573431967310101, 
    -.3595768372350124,-.1977663594959926, 
    0.4731366238769863,-.1975148047452081E-01, 
    -.4533851434024622,0.3153219324427421, 
    -.2494039878683618E-01,-.2903815336558955, 
    -.1831653764008917E-01,-.9317223914837645, 
    0.9500389291238570,-.6699974553713883, 
    0.5322783600952914,0.1377190952760970, 
    0.2782783282780965,-.3217122219145910, 
    0.4343389363650290E-01,0.9746770878915403E-01, 
    -.3733686186921615,0.2759009099030095, 
    0.1089081465931238E-01,-.8195375330783097, 
    0.8086467184189974,-.1251938886030335, 
    0.3429675592810749E-01,0.9089713267491691E-01, 
    0.4868852723886817,0.1198290818690802, 
    -.6067143542577578,0.3382540465835279, 
    -.4603539675728300,0.1220999209893009, 
    -.9480210749913888E-01,-.6451285680033118, 
    0.7399306755024498,0.1140113067761596, 
    -.8552704011446526,0.7412590943684850, 
    -.1658421801556803,-.6245895743792733, 
    0.7904317545349575,0.2557413056493574, 
    0.3270954418922731E-01,-.2884508498385872, 
    -.1613157569571640,-.4583565335140027, 
    0.6196722904711670,0.2156998987651826E-14, 
    0.2157685773910294E-14,0.9329761593328095E-14, 
    -.1937826095544789E-14,0.5112100598082816E-14 };
  double ys[] = { 
    0.3338515555387722,0.1758941973965540, 
    -.5097457529353365,-.2575328894119814, 
    -.3844776049819921,0.6420104943939733, 
    0.3445858452048264,-.2116463454095397, 
    -.1329394997952951,-.4089732623879643, 
    -.9116237193661333E-01,0.5001356343245807, 
    -.2437824086419522,0.7599889189904093E-02, 
    0.2361825194520526,0.5725635993272842, 
    -.3368651110208424E-01,-.5388770882252032, 
    0.5794897706560004,-.3249959253184209, 
    -.2544938453375776,0.5504594216795844, 
    -.3602808234304559,-.1901785982491242, 
    0.4115986775195391,-.1665168003024589E-01, 
    -.3949469974892933,-.1430620777928657, 
    0.6580940175195596E-01,0.7725267604090996E-01, 
    0.2342697395017627,-.2784188228662802, 
    0.4414908336450658E-01,-.4017554364367290, 
    -.3092461241493650,0.7110015605860860, 
    0.9085253735987711,-.4646238165515825, 
    -.4439015570471948,0.4730296395382988, 
    -.1825160620441602,-.2905135774941487, 
    -.6857456848032077E-01,0.2279222388260412, 
    -.1593476703457152,0.2190140235085516, 
    0.3481413428683416E-01,-.2538281577953898, 
    -.4251600332346974E-01,-.4295820758543226E-01, 
    0.8547421090889233E-01,-.2447322131044543, 
    0.2917527702125314,-.4702055710805701E-01, 
    0.2600283704245143,-.1456258542006649, 
    -.1144025162238410,0.2964453769751457, 
    -.2809509380001344,-.1549443897500673E-01, 
    0.7775759193725890,-.3837855053933897, 
    -.3937904139792007,-.3624391892660122, 
    0.1092997287033503,0.2531394605626634, 
    -.2145827819997083E-01,0.4804479220825398, 
    -.4589896438825594,0.9342132288033096E-01, 
    0.4359627055553254,-.5293840284356405, 
    -.2503585120211285,0.5349275917488283, 
    -.2845690797277077,-.1532525106706588, 
    0.3497030592011267,-.1964505485304861, 
    0.1086435405003431E+01,-.5590802894074118, 
    -.5273551155960269,-.2277988977546228, 
    -.4663353679452422,0.6941342656998574, 
    0.2108172081069621,0.1355874975580061, 
    -.3464047056649888,0.3748559371110348, 
    -.1030184566954502,-.2718374804155939, 
    0.9400326158922729,-.4605845857832677, 
    -.4794480301090108,0.3267824275100505E-01, 
    -.1247602093042907,0.9208196655328471E-01, 
    -.4194700484257879,0.6313900388299887, 
    -.2119199904042083,0.3362799093518891, 
    0.1247966425982771,-.4610765519501666, 
    0.7996643270816305,-.4819331969673737, 
    -.3177311301142533,0.9217567340036824, 
    -.3621416790150311,-.5596150549886635, 
    0.8169629451576943,-.5521050136126617, 
    -.2648579315450346,-.1854220399449956, 
    0.3141894874618389,-.1287674475168422, 
    0.6224002317220075,-.4509036594166286, 
    -.1714965723053837,-.1747363709921250E-14, 
    -.7841339838758449E-14,-.1038742413216324E-13, 
    -.6434088705165255E-14,0.8511611593137641E-14 };
  double zs[] = { 
    -.3905577451659811,-.3905577451659810, 
    -.3905577451659885,-.4194761577597843E-02, 
    -.4194761577604291E-02,-.4194761577584145E-02, 
    0.1874027021856000,0.1874027021856146, 
    0.1874027021856214,-.3071102013390232, 
    -.3071102013390279,-.3071102013390283, 
    -.2701961597911837E-01,-.2701961597911512E-01, 
    -.2701961597910549E-01,-.3650223874563840, 
    -.3650223874563843,-.3650223874563830, 
    -.8612392000622994E-01,-.8612392000621880E-01, 
    -.8612392000621745E-01,0.1070791221730147, 
    0.1070791221730203,0.1070791221730225, 
    0.3230374822491166E-01,0.3230374822490626E-01, 
    0.3230374822491062E-01,0.3289827199375365, 
    0.3289827199375268,0.3289827199375392, 
    -.2002324469987915,-.2002324469988037, 
    -.2002324469988010,-.2771825970811163, 
    -.2771825970811178,-.2771825970811177, 
    -.3794145952130761,-.3794145952130777, 
    -.3794145952130777,0.2887592197883361, 
    0.2887592197883261,0.2887592197883259, 
    0.7307953050563157,0.7307953050563196, 
    0.7307953050563138,0.4124918246328610, 
    0.4124918246328517,0.4124918246328527, 
    0.9839020121236216,0.9839020121236166, 
    0.9839020121236378,0.6829751492607604E-01, 
    0.6829751492606395E-01,0.6829751492607447E-01, 
    0.5090267712108498,0.5090267712108671, 
    0.5090267712108663,0.3851617855258526, 
    0.3851617855258565,0.3851617855258507, 
    0.2348420776072121E-01,0.2348420776071207E-01, 
    0.2348420776071097E-01,0.7747281662419676E-01, 
    0.7747281662420324E-01,0.7747281662419568E-01, 
    -.1630338566084656,-.1630338566084718, 
    -.1630338566084695,-.3584562395667113, 
    -.3584562395667106,-.3584562395667084, 
    0.3698626219613211,0.3698626219613297, 
    0.3698626219613188,-.3667298228875690, 
    -.3667298228875713,-.3667298228875721, 
    -.3790081337166384,-.3790081337166408, 
    -.3790081337166418,-.2197227060814419, 
    -.2197227060814434,-.2197227060814370, 
    -.3654645869052863,-.3654645869052863, 
    -.3654645869052867,-.2014726326233556, 
    -.2014726326233535,-.2014726326233503, 
    -.2464405285751480,-.2464405285751552, 
    -.2464405285751559,0.7375395719855852, 
    0.7375395719855836,0.7375395719855902, 
    -.3689240261223760,-.3689240261223751, 
    -.3689240261223750,-.2202508600260902, 
    -.2202508600260935,-.2202508600260907, 
    -.1594720366091775,-.1594720366091723, 
    -.1594720366091740,-.3758886073617869, 
    -.3758886073617894,-.3758886073617892, 
    -.3542700621002714,-.3542700621002700, 
    -.3542700621002702,0.6895589091912253, 
    0.6895589091912366,0.6895589091912278, 
    -.4073916692470354,-.4073916692470370, 
    -.4073916692470400,0.6482117053575507, 
    0.1198603539538932E+01,-.4041509145471567, 
    -.2710584986072719,0.8757555073873696E-02 };
  double ws[] = { 
    0.2025764222341204E-02,0.2025764222341283E-02, 
    0.2025764222340314E-02,0.6116300216663094E-02, 
    0.6116300216663310E-02,0.6116300216663541E-02, 
    0.1275819159212241E-01,0.1275819159212196E-01, 
    0.1275819159212162E-01,0.8334684897884583E-02, 
    0.8334684897884422E-02,0.8334684897884272E-02, 
    0.1507614954107500E-01,0.1507614954107359E-01, 
    0.1507614954107202E-01,0.2734637742973762E-02, 
    0.2734637742973644E-02,0.2734637742974016E-02, 
    0.1154303798541227E-01,0.1154303798541274E-01, 
    0.1154303798541255E-01,0.7163952410599667E-02, 
    0.7163952410599083E-02,0.7163952410598778E-02, 
    0.7580503447649259E-02,0.7580503447650207E-02, 
    0.7580503447650126E-02,0.1894249152240793E-01, 
    0.1894249152240805E-01,0.1894249152240727E-01, 
    0.1641784697591419E-01,0.1641784697591236E-01, 
    0.1641784697591317E-01,0.8436640400263927E-02, 
    0.8436640400264633E-02,0.8436640400264418E-02, 
    0.3036004079076158E-02,0.3036004079076059E-02, 
    0.3036004079076094E-02,0.7492415688372185E-02, 
    0.7492415688372261E-02,0.7492415688372197E-02, 
    0.3062456464426824E-02,0.3062456464427586E-02, 
    0.3062456464427105E-02,0.8856133214568079E-02, 
    0.8856133214567973E-02,0.8856133214568157E-02, 
    0.3071402699546415E-02,0.3071402699546453E-02, 
    0.3071402699546178E-02,0.1823835135198804E-01, 
    0.1823835135198757E-01,0.1823835135198796E-01, 
    0.1121658861163283E-01,0.1121658861163258E-01, 
    0.1121658861163293E-01,0.5177567930581361E-02, 
    0.5177567930581481E-02,0.5177567930580846E-02, 
    0.3242723020974106E-02,0.3242723020974469E-02, 
    0.3242723020974430E-02,0.1136422978208590E-01, 
    0.1136422978208548E-01,0.1136422978208555E-01, 
    0.8174247103577847E-02,0.8174247103578128E-02, 
    0.8174247103578329E-02,0.4083194002127334E-02, 
    0.4083194002127503E-02,0.4083194002128160E-02, 
    0.3128101586450114E-02,0.3128101586450456E-02, 
    0.3128101586450350E-02,0.1081664090579396E-01, 
    0.1081664090579364E-01,0.1081664090579339E-01, 
    0.4910680133298682E-03,0.4910680133297917E-03, 
    0.4910680133297780E-03,0.8074513061447679E-02, 
    0.8074513061447566E-02,0.8074513061448123E-02, 
    0.1112708329974781E-01,0.1112708329974803E-01, 
    0.1112708329974812E-01,0.2131526487678964E-01, 
    0.2131526487679003E-01,0.2131526487678986E-01, 
    0.3571026646377801E-02,0.3571026646377656E-02, 
    0.3571026646377794E-02,0.7768447076866043E-02, 
    0.7768447076865915E-02,0.7768447076866019E-02, 
    0.7618830679555172E-02,0.7618830679555604E-02, 
    0.7618830679555573E-02,0.1180399216908969E-01, 
    0.1180399216908947E-01,0.1180399216908933E-01, 
    0.2753708273597594E-02,0.2753708273597404E-02, 
    0.2753708273597407E-02,0.1252040626033346E-02, 
    0.1252040626033167E-02,0.1252040626033129E-02, 
    0.2035505637330654E-02,0.2035505637330660E-02, 
    0.2035505637330673E-02,0.2152733016776995E-02, 
    0.2152733016776493E-02,0.2152733016776350E-02, 
    0.2538969965943544E-02,0.2538969965943405E-02, 
    0.2538969965943241E-02,0.1078372873958245E-01, 
    0.1943908493037603E-03,0.4799528389705479E-02, 
    0.2477673544485754E-01,0.2864883777302190E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule13 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE13 returns the rule of degree 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.4777385017466611E-01,-.7568808830760128, 
    0.8046547332507396,0.3005907944967990E-01, 
    -.3472424886466524,0.3171834091971876, 
    -.3529824413110864E-01,-.4209487812102883, 
    0.4562470253414500,-.7153537080888878, 
    0.5989558965388279,0.1163978115501301, 
    0.3710519456025519E-02,0.4652335033519287, 
    -.4689440228078527,-.3599160575353618E-01, 
    0.1939924024452817,-.1580007966914678, 
    -.2161537359848647E-01,0.2912599485905940, 
    -.2696445749921000,-.4913573074965884, 
    0.2288888183023164,0.2624684891944182, 
    -.1426281098822192,-.3118486240867259, 
    0.4544767339690260,0.7618754763810419, 
    -.1422306180040716,-.6196448583768703, 
    0.5776897880027803E-01,0.1164689180950719, 
    -.1742378968955057,-.1414719142296262, 
    -.3229372289650024,0.4644091431945710, 
    -.1398258170326137,-.5098850522745304, 
    0.6497108693071825,0.3610318921371553, 
    0.1054629929099188,-.4664948850470391, 
    0.2459013465212476E-01,-.7178146556427887, 
    0.6932245209907225,-.1065595141175732E-01, 
    -.9217121855780308,0.9323681369897492, 
    -.4115285528872360,0.2668526805669106, 
    0.1446758723202021,0.1571190412397715E-02, 
    -.6922217370423006E-01,0.6765098329171684E-01, 
    -.2991189861629327,-.1245621803681009, 
    0.4236811665311025,0.3871696430393148, 
    0.1500785041234294,-.5372481471625822, 
    -.8064886300342536E-01,-.4844598247542774, 
    0.5651086877577737,0.6705603853993255, 
    -.2613693382805031,-.4091910471186840, 
    -.2994914861160214,-.2674832138785850, 
    0.5669746999946459,0.4510147973428894E-01, 
    0.2172069032075606,-.2623083829418457, 
    -.2495785492724049,-.2085110500880744, 
    0.4580895993605156,0.2707504828489703, 
    -.4821600751443009,0.2114095922954477, 
    -.4554280384514680E-01,-.7996627976605564, 
    0.8452056015056444,0.2029674430134095E-01, 
    -.7997175993984930,0.7794208550971885, 
    0.1441756236446537,-.6915442312985614, 
    0.5473686076540003,0.3528170518753222E-01, 
    -.3311753761570134,0.2958936709696454, 
    -.7827835701762442E-02,0.2760083802741976, 
    -.2681805445724308,0.6806655606966093E-01, 
    -.5973241026022454,0.5292575465325308, 
    -.9685067665231929E-01,0.2259602937925122, 
    -.1291096171402099,-.1001708030503953E-01, 
    -.5630399697700922,0.5730570500751420, 
    0.3465920273834422E-01,-.1777529113458167, 
    0.1430937086076311,0.7285263734559304E-01, 
    -.8758889229618175,0.8030362856162614, 
    0.1343389663918027,0.3445210697535890, 
    -.4788600361453874,-.1348321313217423E-01, 
    0.1154680832801965,-.1019848701479677, 
    -.5176848990301737E-01,-.2375948027597131, 
    0.2893632926627010,0.1729207472857360E-01, 
    -.4389324777852445,0.4216404030566085, 
    -.3914791775164566E-02,-.2309684344749786, 
    0.2348832262501366,-.6881144138718593, 
    0.4266860230809687,0.2614283907909357, 
    0.4960973006288234E-02,-.2137441872471665, 
    0.2087832142407756,-.7325006065249269E-01, 
    -.2659488483143491,0.3391989089668666, 
    -.4678501489682662E-01,-.5410293892575722, 
    0.5878144041543488,-.2500405721723705E-01, 
    -.6313699376952513,0.6563739949124212, 
    0.1364199568193839,-.2996193314702689, 
    0.1631993746509763,0.1198743702984181E-12, 
    -.7662979863046558E-13,-.4829662024390088E-14, 
    -.2111050829994401E-13,-.1270671025297558E-13 };
  double ys[] = { 
    0.9015530084355994,-.4921498721057091, 
    -.4094031363297998,0.3836064709764307, 
    -.1657713090704851,-.2178351619061508, 
    0.5064492350445715,-.2837937936488191, 
    -.2226554413958301,-.2786050402678961, 
    -.4802119637624816,0.7588170040302994, 
    -.5393476461992264,0.2728872272096964, 
    0.2664604189895141,-.2032233682745297, 
    0.7044203923151297E-01,0.1327813290428652, 
    -.3238383776801540,0.1431997261914898, 
    0.1806386514886541,0.1938723202896774E-01, 
    -.4352215266416702,0.4158342946127344, 
    0.4424381517603146,-.3447386423319839, 
    -.9769950942847148E-01,-.2756352401942940, 
    0.7976211371634238,-.5219858969692130, 
    -.1678396578899558,0.1339492321368408, 
    0.3389042575316226E-01,0.4545746399118207, 
    -.3498055916007531,-.1047690483110467, 
    0.6694930174765074,-.4558392183934259, 
    -.2136537990830517,-.3302200348036231, 
    0.4777728075688429,-.1475527727653846, 
    0.8146638484662959,-.3860362429420576, 
    -.4286276055244945,0.1070453773333782E+01, 
    -.5444552112909320,-.5259985620427312, 
    -.7053881312997473E-01,-.3211247746179037, 
    0.3916635877479940,0.7902375403647920E-01, 
    -.3815118620693258E-01,-.4087256782951131E-01, 
    0.3165284439137526,-.4173088627282204, 
    0.1007804188144937,-.3968282271412536, 
    0.5337128600166933,-.1368846328756112, 
    0.6059686632316962,-.3728282957631877, 
    -.2331403674686730,-.8534490338988292E-01, 
    0.6233947802222222,-.5380498768323672, 
    0.4817745012020841,-.5002544857946896, 
    0.1847998459249884E-01,-.2768482795388788, 
    0.1774831669675317,0.9936511257127259E-01, 
    0.3848620644010382,-.4085723961101285, 
    0.2371033170906487E-01,0.4004326341980497, 
    0.3426047913492075E-01,-.4346931133331124, 
    0.9496652130401684,-.5142738316095435, 
    -.4353913814305450,0.9117160117907287, 
    -.4382805097163334,-.4734355020744341, 
    0.7152866610716716,-.2327835778531099, 
    -.4825030832186010,0.3620384831590592, 
    -.1504643885983839,-.2115740945607743, 
    -.3141876222502283,0.1503147065506677, 
    0.1638729156994478,0.6504322183920584, 
    -.2662687424915907,-.3841634759004985, 
    -.2049997086581880,0.1862470797446676E-01, 
    0.1863750006837488,0.6559259202331654, 
    -.3366380061325622,-.3192879141007369, 
    0.1852408823986869,-.6260469115301880E-01, 
    -.1226361912458057,0.9693279211217602, 
    -.4215717258868716,-.5477561952349157, 
    -.4753793031363919,0.3540306091815932, 
    0.1213486939547357,-.1255465211978510, 
    0.5109645550180968E-01,0.7445006569601835E-01, 
    0.3042393982437668,-.1969525264934776, 
    -.1072868717503699,0.4968519844113532, 
    -.2334506162065882,-.2634013682047193, 
    0.2689595817220580,-.1378700999889159, 
    -.1310894817332517,-.9541153848828572E-01, 
    -.5482187938791565,0.6436303323674273, 
    0.2439463089891027,-.1176768258436301, 
    -.1262694831454210,0.3493822205657847, 
    -.2381275236366834,-.1112546969291192, 
    0.6517382679995124,-.3663861454168269, 
    -.2853521225825763,0.7434793061384100, 
    -.3933938018169897,-.3500855043213671, 
    0.2672085045651081,-.1546110409379758E-01, 
    -.2517474004713509,-.6195967642546056E-13, 
    0.3002370388193685E-13,-.5457794979132388E-13, 
    0.1164494042115649E-12,0.3123283355285955E-13 };
  double zs[] = { 
    -.2327715323996341,-.2327715323996078, 
    -.2327715323995431,0.3496762696550786, 
    0.3496762696549099,0.3496762696549561, 
    0.2717988565409110,0.2717988565409942, 
    0.2717988565407555,-.3874128688211647, 
    -.3874128688212131,-.3874128688212028, 
    -.3783145937148669,-.3783145937148664, 
    -.3783145937148837,-.3908054363428844, 
    -.3908054363428517,-.3908054363428521, 
    0.2214621132658698,0.2214621132658886, 
    0.2214621132658825,-.5332156026748039E-01, 
    -.5332156026744841E-01,-.5332156026745260E-01, 
    -.1366126720097484,-.1366126720097481, 
    -.1366126720097814,-.3405718815923187, 
    -.3405718815923307,-.3405718815923774, 
    0.2482158558753812,0.2482158558753821, 
    0.2482158558753401,0.1493078550521198, 
    0.1493078550521607,0.1493078550521824, 
    -.1246861001425667,-.1246861001425311, 
    -.1246861001424509,-.8804090671666169E-01, 
    -.8804090671661184E-01,-.8804090671646460E-01, 
    -.2903904416529835,-.2903904416529482, 
    -.2903904416530703,-.3588716811515187, 
    -.3588716811515073,-.3588716811515079, 
    0.1748649644824907,0.1748649644825259, 
    0.1748649644824762,0.1005626022480115E+01, 
    0.1005626022480066E+01,0.1005626022480308E+01, 
    -.6085363845277859E-01,-.6085363845275459E-01, 
    -.6085363845275343E-01,-.2907416212662024, 
    -.2907416212662147,-.2907416212661059, 
    -.2928675049685535,-.2928675049684912, 
    -.2928675049686165,-.3986190384060752, 
    -.3986190384060841,-.3986190384061074, 
    -.2763005028815692,-.2763005028815834, 
    -.2763005028815310,-.3171132373617958E-01, 
    -.3171132373613690E-01,-.3171132373618723E-01, 
    -.3670072924650872,-.3670072924650821, 
    -.3670072924650642,-.3918966809387769, 
    -.3918966809387540,-.3918966809387686, 
    -.3945112928982144,-.3945112928982283, 
    -.3945112928982709,-.1689813827450514, 
    -.1689813827450869,-.1689813827450585, 
    -.1988456360856388,-.1988456360856537, 
    -.1988456360856610,-.3581788966687771, 
    -.3581788966687613,-.3581788966687278, 
    -.2652570624981906,-.2652570624981415, 
    -.2652570624981432,0.9127655651197962E-01, 
    0.9127655651192465E-01,0.9127655651200027E-01, 
    0.4712027028361149,0.4712027028361768, 
    0.4712027028361248,-.4789101434801453E-01, 
    -.4789101434802212E-01,-.4789101434808656E-01, 
    0.5562071598537035,0.5562071598536467, 
    0.5562071598533351,-.3666902411131864, 
    -.3666902411132085,-.3666902411132190, 
    -.2500081968302335,-.2500081968302558, 
    -.2500081968302762,0.7953910781041313, 
    0.7953910781041719,0.7953910781041406, 
    0.1246295654192337,0.1246295654192380, 
    0.1246295654192949,0.4430780439124138, 
    0.4430780439124452,0.4430780439125641, 
    -.1789190687389412,-.1789190687388659, 
    -.1789190687388622,-.3621753277358900, 
    -.3621753277358689,-.3621753277358712, 
    0.7590201625172059,0.7590201625171641, 
    0.7590201625173268,0.5385097619725802, 
    0.5385097619725214,0.5385097619724014, 
    -.4034889863680245,-.4034889863679839, 
    -.4034889863680722,0.9129480050841532E-01, 
    0.9129480050844316E-01,0.9129480050850451E-01, 
    0.5076476008019328,0.5076476008020036, 
    0.5076476008019484,0.7587686360919684, 
    0.1174064403165666E+01,-.2917588131252908, 
    0.4072850512806022,-.2737975651676979E-02 };
  double ws[] = { 
    0.2068909148962774E-02,0.2068909148960302E-02, 
    0.2068909148962221E-02,0.6862179883351744E-02, 
    0.6862179883357749E-02,0.6862179883341077E-02, 
    0.4745055930506758E-02,0.4745055930499521E-02, 
    0.4745055930498638E-02,0.2639271833496386E-02, 
    0.2639271833492975E-02,0.2639271833493415E-02, 
    0.2420284628312029E-02,0.2420284628310808E-02, 
    0.2420284628311964E-02,0.5428029730193044E-02, 
    0.5428029730198515E-02,0.5428029730197811E-02, 
    0.7363269664287002E-02,0.7363269664287796E-02, 
    0.7363269664288242E-02,0.4622589345403796E-02, 
    0.4622589345400746E-02,0.4622589345402929E-02, 
    0.1113014179512893E-01,0.1113014179512515E-01, 
    0.1113014179513268E-01,0.3019898782414975E-02, 
    0.3019898782415646E-02,0.3019898782415196E-02, 
    0.1497918593186746E-01,0.1497918593186317E-01, 
    0.1497918593186622E-01,0.6570279132551560E-02, 
    0.6570279132554112E-02,0.6570279132553927E-02, 
    0.4504033084972921E-02,0.4504033084968776E-02, 
    0.4504033084972839E-02,0.1239429415995445E-01, 
    0.1239429415996196E-01,0.1239429415995679E-01, 
    0.5367281142254301E-02,0.5367281142251953E-02, 
    0.5367281142250357E-02,0.6316169904420855E-03, 
    0.6316169904431759E-03,0.6316169904440820E-03, 
    0.8876101030265546E-02,0.8876101030282007E-02, 
    0.8876101030274161E-02,0.2162421606464688E-02, 
    0.2162421606465037E-02,0.2162421606462924E-02, 
    0.8054280667477870E-02,0.8054280667477506E-02, 
    0.8054280667478849E-02,0.9049885628396499E-02, 
    0.9049885628398672E-02,0.9049885628405025E-02, 
    0.9196546074381425E-02,0.9196546074380063E-02, 
    0.9196546074381355E-02,0.1389679560298547E-02, 
    0.1389679560298251E-02,0.1389679560296835E-02, 
    0.5469361348986402E-02,0.5469361348984248E-02, 
    0.5469361348984599E-02,0.1725450508839409E-01, 
    0.1725450508839278E-01,0.1725450508839148E-01, 
    0.6733913747887457E-02,0.6733913747885728E-02, 
    0.6733913747886398E-02,0.4022564578398640E-02, 
    0.4022564578402713E-02,0.4022564578399363E-02, 
    0.1301944787144861E-02,0.1301944787144253E-02, 
    0.1301944787143136E-02,0.2027792916616635E-02, 
    0.2027792916617886E-02,0.2027792916616735E-02, 
    0.4379548102164245E-02,0.4379548102164029E-02, 
    0.4379548102160293E-02,0.9788773323829804E-02, 
    0.9788773323830982E-02,0.9788773323832441E-02, 
    0.1491498412616913E-01,0.1491498412617151E-01, 
    0.1491498412617137E-01,0.3816494659634458E-02, 
    0.3816494659634902E-02,0.3816494659631925E-02, 
    0.1053753256419563E-01,0.1053753256419800E-01, 
    0.1053753256419685E-01,0.9378412737876001E-02, 
    0.9378412737873921E-02,0.9378412737876034E-02, 
    0.1034943926448387E-01,0.1034943926449078E-01, 
    0.1034943926449734E-01,0.1253515438165336E-02, 
    0.1253515438163976E-02,0.1253515438162651E-02, 
    0.7631115735488028E-02,0.7631115735488942E-02, 
    0.7631115735487678E-02,0.4662247664665583E-02, 
    0.4662247664663175E-02,0.4662247664669866E-02, 
    0.1924508745358327E-01,0.1924508745357723E-01, 
    0.1924508745358518E-01,0.2150599791391388E-02, 
    0.2150599791389261E-02,0.2150599791389503E-02, 
    0.1815794181007045E-01,0.1815794181007355E-01, 
    0.1815794181007060E-01,0.2074509937828351E-02, 
    0.2074509937826631E-02,0.2074509937827053E-02, 
    0.3856301126724407E-02,0.3856301126725213E-02, 
    0.3856301126727490E-02,0.2342776832115179E-02, 
    0.2342776832117610E-02,0.2342776832119853E-02, 
    0.2485903899501151E-02,0.2485903899504086E-02, 
    0.2485903899497665E-02,0.1958269886432824E-02, 
    0.1958269886433108E-02,0.1958269886434393E-02, 
    0.3112037442333638E-02,0.3112037442336383E-02, 
    0.3112037442329971E-02,0.8380938111968996E-02, 
    0.2408616607470909E-03,0.1651223797192308E-01, 
    0.1363274166717042E-01,0.2507433395640973E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule14 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE14 returns the rule of degree 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    -.5752422835841802,0.4314967658691802, 
    0.1437455177139820,-.3303076984155172, 
    0.4866309056595519,-.1563232072429510, 
    -.5770168798465761E-02,-.7016819527781105, 
    0.7074521215770434,-.3032416710412216, 
    -.2339888323701999,0.5372305034104569, 
    0.2995840525991514,-.4325354814520043, 
    0.1329514288534004,0.2313978508942549, 
    -.4060579401191659E-01,-.1907920568806182, 
    0.1198459969693036,-.5485013403177168, 
    0.4286553433493205,-.4223067018453286, 
    0.2051376117542425,0.2171690900904567, 
    -.1995552059487330,-.3336302067365579, 
    0.5331854126840413,-.3362367775957516E-02, 
    -.2082502914064486,0.2116126591830406, 
    0.3901789196843061,-.4123770748876975E-03, 
    -.3897665426036371,0.5440475521812093E-02, 
    -.5671976659757272,0.5617571904593656, 
    0.5109903009675820E-01,-.8812309166767075, 
    0.8301318865798517,-.3334033396451877, 
    -.4849473999181729E-01,0.3818980796337723, 
    0.2967209217147687,0.6320475373066934E-02, 
    -.3030413970906470,0.3947081540530285, 
    -.2978075731654519,-.9690058088341083E-01, 
    0.4602277425410929,-.2588695285519879, 
    -.2013582139898412,0.4320204013442107, 
    0.2158666376867163,-.6478870390304328, 
    0.3243502588043465E-02,-.8153377567221353, 
    0.8120942541315551,-.2929639748058975, 
    0.3937181324850053E-01,0.2535921615554882, 
    -.2911712364207900,0.4448737339859005, 
    -.1537024975672874,0.3770577628402459, 
    -.4566844608555435,0.7962669801609278E-01, 
    -.7343568438702067E-03,-.6012198115212343E-01, 
    0.6085633799849959E-01,-.7860612444751536, 
    0.7072341550237046,0.7882708945288715E-01, 
    -.1798079430257987,-.5716591672723068, 
    0.7514671102985679,0.2349121749284232E-01, 
    -.2977965988756665,0.2743053813811669, 
    -.3508185900315901,0.3561089208016359, 
    -.5290330768361618E-02,-.4799997070573949, 
    0.3284477153105494,0.1515519917450605, 
    0.5503749362395663,-.1236979845612246, 
    -.4266769516797990,-.1160913996018244, 
    0.1218716230766462,-.5780223478460763E-02, 
    -.1795966192456130,-.5071498846594161, 
    0.6867465039037227,-.4731311596169461E-02, 
    -.6889415863849245,0.6936728979813912, 
    0.1838019755826437,-.7414338401434005, 
    0.5576318645605018,-.3157482778904093, 
    0.1186260398816200E-01,0.3038856739015107, 
    -.3760156844932397,0.2533165859904101, 
    0.1226990985037814,0.3132847887437243E-01, 
    0.1567813343830147,-.1881098132555542, 
    -.5940398288585815,0.6767061397330466, 
    -.8266631087504343E-01,0.5799063749491890, 
    0.7755651500224905E-01,-.6574628899522570, 
    0.1959288531923632,0.8941019681549953E-01, 
    -.2853390500081148,0.2938476934905306, 
    0.2019773488401272E-02,-.2958674669799366, 
    0.3655020525321305,-.5615458414677339, 
    0.1960437889355228,-.5706595245720943E-01, 
    -.7715171223939032,0.8285830748522059, 
    0.1011528057210596,-.1051837271533038, 
    0.4030921431232027E-02,-.2124197258439737, 
    0.3605072524272380,-.1480875265845512, 
    -.1651782199629099,0.7878215733569789E-02, 
    0.1573000042284481,-.3617927525610758, 
    -.2009201308789826,0.5627128834404685, 
    -.5015820600281335,0.4144850071548899, 
    0.8709705287395679E-01,0.9434068372117834, 
    0.4463509082121593E-02,-.9478703462960307, 
    0.2005751720829750,-.1945620110873523, 
    -.6013160997523850E-02,-.1255658338458569E-01, 
    -.4674319437690188,0.4799885271519170, 
    -.1716908737168969E-01,0.1956591346765856, 
    -.1784900473056204,-.5462757915789873E-03, 
    -.5780834828017836,0.5786297585932533, 
    -.1833877013210589,0.2818210581349795, 
    -.9843335681437262E-01,0.1916617895075464, 
    -.7885465203364792E-03,-.1908732429835174, 
    -.4813891231174998E-01,-.8618855209039119, 
    0.9100244332179739,-.4424462871708599, 
    0.5232153544540827,-.8076906728405848E-01, 
    0.2330164669103972,-.5137198179676699, 
    0.2807033510564895,0.3967080654003286E-12, 
    0.1722425458471176E-11,0.3502410373140708E-12, 
    -.2947942461866533E-11,-.1032091026644403E-11, 
    -.4359393565820349E-12 };
  double ys[] = { 
    -.1661332605832194,-.4151078006247312, 
    0.5812410612053462,-.3712097301616572, 
    -.1004499928153783,0.4716597229753685, 
    0.8135639371543177,-.4117790813394923, 
    -.4017848558132318,0.4452636911188574, 
    -.4852468361674822,0.3998314505168596E-01, 
    0.3264840198869604,0.9620539017564665E-01, 
    -.4226894100630985,-.8671007929617915E-01, 
    0.2437514569045547,-.1570413776087273, 
    0.5641616743542007,-.1782911592597502, 
    -.3858705150939825,0.6946377256844247E-02, 
    -.3692015206135402,0.3622551433574067, 
    0.5004562312099742,-.4230479934143149, 
    -.7740823779440852E-01,0.2424079875448238, 
    -.1241158896833233,-.1182920978621629, 
    -.2247937322813834,0.4503017226062891, 
    -.2255079903239726,0.6518023902669200, 
    -.3211896051207328,-.3306127851443094, 
    0.9880557751410730,-.4497748293987996, 
    -.5382809457428603,0.2484874102685690, 
    -.4129794669742624,0.1644920567082871, 
    -.1786101603431211,0.3462729362116060, 
    -.1676627758691386,0.1159937060740129, 
    0.2838304354535469,-.3998241415280471, 
    0.3320417294495800E-01,0.3819668300942372, 
    -.4151710030395980,-.4986884177668357, 
    0.6234848514007479,-.1247964336323278, 
    0.9395983095546842,-.4669901991392865, 
    -.4726081104148907,0.1236801757609007, 
    -.3155543324553964,0.1918741566946638, 
    -.3455881484181401,-.7936761338079117E-01, 
    0.4249557617992578,0.3096393919455638, 
    0.1717219053442890,-.4813612972863747, 
    0.6984686512807867E-01,-.3555940424550629E-01, 
    -.3428746088226146E-01,-.3628109884665136, 
    -.4993435124130320,0.8621545008801837, 
    0.7639073125277335,-.5376719027263360, 
    -.2262354098013041,0.3303032323049931, 
    -.1448076250378262,-.1854956072669932, 
    -.2086539551781068,-.1994908334986071, 
    0.4081447886767858,-.1021307936198597, 
    -.3646265433107137,0.4667573369306311, 
    -.1749249882226239,0.5641011705032603, 
    -.3891761822782046,-.7369982797233314E-01, 
    -.6368818722810367E-01,0.1373880152000548, 
    0.6892964013212863,-.5001834353611165, 
    -.1891129659598641,0.7982528447344757, 
    -.4032238584029919,-.3950289863308680, 
    0.7500159343054114,-.2158307870331475, 
    -.5341851472730395,0.1685995980230182, 
    -.3577458288655779,0.1891462308433379, 
    -.7541204156110640E-01,-.2879331142125014, 
    0.3633451557721906,-.1991229969303482, 
    0.1266927570309605,0.7243023989959384E-01, 
    -.4384238887740483,-.2952416382641486, 
    0.7336655270377042,-.4243636513101600, 
    0.7143954781786902,-.2900318268675506, 
    -.2163615785319861,0.2778601534641544, 
    -.6149857493356037E-01,-.1719852784726787, 
    0.3404722066425572,-.1684869281697858, 
    0.4373945770482864,0.9783677410376379E-01, 
    -.5352313511523523,0.9238182796099402, 
    -.5113297043238586,-.4124885752868345, 
    0.6305510675913090E-01,0.5607334603928896E-01, 
    -.1191284527989022,-.2936373325714834, 
    -.3714221255853142E-01,0.3307795451305264, 
    0.8626870981050186E-01,-.1861828895445977, 
    0.9991417973422896E-01,0.4408837263795412, 
    -.5337635778128931,0.9287985143315733E-01, 
    -.1890175235323525,-.3398740442993634, 
    0.5288915678318241,-.5498302077599997, 
    0.1091929391010154E+01,-.5420991832506766, 
    0.1088587293550235,0.1192738297139006, 
    -.2281325590696736,0.5469934639217403, 
    -.2843710521572813,-.2626224117648735, 
    -.2160151309343958,0.9313869964408443E-01, 
    0.1228764312899015,0.6678287012945778, 
    -.3343874393606580,-.3334412619346900, 
    -.2195399888316988,-.4904841366906035E-01, 
    0.2685884025019004,-.1097454506740659, 
    0.2208567039838767,-.1111112533091888, 
    0.1023012688991847E+01,-.5531958654682623, 
    -.4698168235239601,-.3487105684781340, 
    -.2088144402600705,0.5575250087387288, 
    0.4586604304883936,-.2753203539748657E-01, 
    -.4311283950878784,0.1424849530503414E-12, 
    0.4163207520226179E-12,0.2564407634679153E-12, 
    -.9489253825677672E-12,0.5725707601904812E-12, 
    -.4307858177802721E-12 };
  double zs[] = { 
    -.3002408419822953,-.3002408419833086, 
    -.3002408419833821,-.2367243537500485, 
    -.2367243537456669,-.2367243537470976, 
    -.1124874139319764E-01,-.1124874139141120E-01, 
    -.1124874139183386E-01,-.2091200491939335, 
    -.2091200491951212,-.2091200491943568, 
    -.3871070169557922,-.3871070169557220, 
    -.3871070169561885,0.1576758205577684, 
    0.1576758205557184,0.1576758205553296, 
    -.3804994466636671,-.3804994466631279, 
    -.3804994466634047,-.2566785126781401, 
    -.2566785126792934,-.2566785126782982, 
    -.7550445623527601E-02,-.7550445625559176E-02, 
    -.7550445625114355E-02,0.3966685271521035, 
    0.3966685271469443,0.3966685271536675, 
    0.5436295368746583,0.5436295368772041, 
    0.5436295368824269,0.2108533972691848, 
    0.2108533972755834,0.2108533972664533, 
    -.3846227632566736,-.3846227632563715, 
    -.3846227632570768,0.4320199156927147E-01, 
    0.4320199156864359E-01,0.4320199157222126E-01, 
    -.3942552261036596,-.3942552261036111, 
    -.3942552261035123,-.2900823721656310, 
    -.2900823721654074,-.2900823721650085, 
    -.4001954700908103,-.4001954700902223, 
    -.4001954700893645,-.2514599966999538, 
    -.2514599967001752,-.2514599967005651, 
    -.2052993465249218,-.2052993465263826, 
    -.2052993465223245,0.2946594740406009, 
    0.2946594740416467,0.2946594740425612, 
    -.5795226862357377E-01,-.5795226862121653E-01, 
    -.5795226862189089E-01,-.2212157144490026, 
    -.2212157144492802,-.2212157144486162, 
    0.1031915041617615E+01,0.1031915041621368E+01, 
    0.1031915041615854E+01,-.2843664701680810, 
    -.2843664701691355,-.2843664701684534, 
    -.3916632709551901,-.3916632709550742, 
    -.3916632709550005,-.5853404503565347E-01, 
    -.5853404503336719E-01,-.5853404503471410E-01, 
    0.4349032593314631,0.4349032593318472, 
    0.4349032593321483,-.7166469320204438E-01, 
    -.7166469320319802E-01,-.7166469320166638E-01, 
    -.3651020078507813,-.3651020078509204, 
    -.3651020078510956,0.6666741608359100, 
    0.6666741608427517,0.6666741608415513, 
    -.2850607980246553,-.2850607980246767, 
    -.2850607980239678,-.3767139049041499, 
    -.3767139049041446,-.3767139049041464, 
    -.3847323113293710,-.3847323113294659, 
    -.3847323113294736,-.1430711791014896E-01, 
    -.1430711790970269E-01,-.1430711791041984E-01, 
    0.2226893434846138,0.2226893434844011, 
    0.2226893434848597,0.1861368718435910, 
    0.1861368718416253,0.1861368718431125, 
    -.1072795622749953,-.1072795622745489, 
    -.1072795622748554,-.6485516948547954E-01, 
    -.6485516948561629E-01,-.6485516948569679E-01, 
    0.5725622693459446,0.5725622693472699, 
    0.5725622693460108,-.2791087112264201, 
    -.2791087112266160,-.2791087112254542, 
    -.3714154371111505,-.3714154371112492, 
    -.3714154371110855,-.3161738417128371, 
    -.3161738417122791,-.3161738417134797, 
    0.8373277100424955,0.8373277100425663, 
    0.8373277100420310,0.2498873918545072, 
    0.2498873918552929,0.2498873918556660, 
    0.5108904140436563,0.5108904140435110, 
    0.5108904140435037,-.3706147756740159, 
    -.3706147756740524,-.3706147756740479, 
    0.2433004560275582,0.2433004560300772, 
    0.2433004560291551,-.3538702545014873, 
    -.3538702545021283,-.3538702545024522, 
    -.3582651744000784,-.3582651743999515, 
    -.3582651744000865,0.8755062971445479E-01, 
    0.8755062971455986E-01,0.8755062971379197E-01, 
    -.1397307998126365,-.1397307998123260, 
    -.1397307998124737,-.2183761124774179, 
    -.2183761124772129,-.2183761124776710, 
    0.5763507816362148,0.5763507816368744, 
    0.5763507816365103,0.8062208425433149, 
    0.8062208425454057,0.8062208425485021, 
    -.4049658784349283,-.4049658784337947, 
    -.4049658784361176,0.2225706813796632, 
    0.2225706813790471,0.2225706813801265, 
    -.4096725239466332E-02,-.4096725240921904E-02, 
    -.4096725238250503E-02,-.2592447810085974, 
    0.1191149121741652E+01,0.4120653127511454, 
    0.8435262601344188,-.4003639673386698, 
    0.3321002545060888E-01 };
  double ws[] = { 
    0.4194830168311239E-02,0.4194830168326290E-02, 
    0.4194830168354955E-02,0.6505251369757091E-02, 
    0.6505251369749837E-02,0.6505251369724089E-02, 
    0.1511857423469778E-02,0.1511857423531536E-02, 
    0.1511857423435278E-02,0.3576693638591392E-02, 
    0.3576693638610991E-02,0.3576693638618915E-02, 
    0.3373656058688358E-02,0.3373656058691433E-02, 
    0.3373656058671795E-02,0.1076784154760114E-01, 
    0.1076784154761919E-01,0.1076784154751974E-01, 
    0.3758938857122322E-02,0.3758938857165091E-02, 
    0.3758938857143189E-02,0.8302444588410703E-02, 
    0.8302444588357655E-02,0.8302444588429193E-02, 
    0.2942777031213149E-02,0.2942777031243674E-02, 
    0.2942777031245480E-02,0.8843589566618218E-02, 
    0.8843589566569863E-02,0.8843589566624911E-02, 
    0.9759028519288406E-03,0.9759028519397565E-03, 
    0.9759028519575427E-03,0.2004802238525724E-02, 
    0.2004802238518510E-02,0.2004802238530063E-02, 
    0.9695765259457370E-03,0.9695765259577475E-03, 
    0.9695765259399248E-03,0.1954897938090337E-02, 
    0.1954897938078232E-02,0.1954897938036924E-02, 
    0.3691846518220940E-02,0.3691846518218663E-02, 
    0.3691846518232707E-02,0.7970801887754822E-02, 
    0.7970801887704541E-02,0.7970801887626933E-02, 
    0.1979189926163023E-02,0.1979189926220139E-02, 
    0.1979189926268928E-02,0.3698742482337427E-02, 
    0.3698742482364018E-02,0.3698742482353030E-02, 
    0.1687317723859350E-02,0.1687317723863163E-02, 
    0.1687317723875827E-02,0.3835254087109209E-02, 
    0.3835254087144783E-02,0.3835254087114301E-02, 
    0.1016260451104831E-01,0.1016260451098046E-01, 
    0.1016260451102241E-01,0.5267012161743862E-02, 
    0.5267012161757533E-02,0.5267012161774349E-02, 
    0.1528430187195588E-02,0.1528430187157122E-02, 
    0.1528430187208298E-02,0.2982437800867588E-02, 
    0.2982437800865273E-02,0.2982437800863566E-02, 
    0.1466557323209961E-02,0.1466557323217223E-02, 
    0.1466557323219012E-02,0.1461057795662250E-01, 
    0.1461057795662684E-01,0.1461057795667047E-01, 
    0.5561076940885510E-02,0.5561076940819906E-02, 
    0.5561076940859880E-02,0.9751977499183117E-02, 
    0.9751977499261835E-02,0.9751977499259098E-02, 
    0.5302530732607095E-02,0.5302530732562353E-02, 
    0.5302530732501225E-02,0.7136244952723934E-02, 
    0.7136244952687442E-02,0.7136244952687213E-02, 
    0.3962851386077182E-02,0.3962851386075291E-02, 
    0.3962851386117417E-02,0.3518864309482585E-02, 
    0.3518864309488224E-02,0.3518864309487660E-02, 
    0.1910503466789297E-02,0.1910503466785413E-02, 
    0.1910503466785542E-02,0.1114372702758599E-01, 
    0.1114372702763622E-01,0.1114372702759230E-01, 
    0.9339073844350225E-02,0.9339073844312899E-02, 
    0.9339073844325385E-02,0.1427519172629088E-01, 
    0.1427519172633135E-01,0.1427519172637184E-01, 
    0.4705321020321325E-02,0.4705321020382762E-02, 
    0.4705321020369120E-02,0.4787930353379107E-02, 
    0.4787930353341169E-02,0.4787930353337981E-02, 
    0.3363836009342347E-02,0.3363836009325845E-02, 
    0.3363836009335837E-02,0.1284569887457350E-01, 
    0.1284569887456974E-01,0.1284569887463706E-01, 
    0.2752818417406538E-02,0.2752818417393949E-02, 
    0.2752818417406429E-02,0.2544892612698709E-02, 
    0.2544892612676544E-02,0.2544892612707782E-02, 
    0.3052619386398206E-02,0.3052619386307692E-02, 
    0.3052619386349196E-02,0.8865680919273325E-02, 
    0.8865680919243184E-02,0.8865680919267345E-02, 
    0.9030103498125034E-02,0.9030103498122236E-02, 
    0.9030103498131978E-02,0.2991177381103681E-02, 
    0.2991177381087045E-02,0.2991177381095963E-02, 
    0.2211164766946887E-02,0.2211164766969601E-02, 
    0.2211164766996503E-02,0.2973562062279324E-03, 
    0.2973562062207861E-03,0.2973562062075669E-03, 
    0.9172628789611696E-02,0.9172628789611484E-02, 
    0.9172628789627685E-02,0.9783236108387362E-02, 
    0.9783236108363432E-02,0.9783236108405431E-02, 
    0.1738125981964304E-01,0.1738125981965133E-01, 
    0.1738125981965202E-01,0.9766706161581359E-02, 
    0.9766706161581049E-02,0.9766706161590663E-02, 
    0.3005216389199320E-02,0.3005216389178637E-02, 
    0.3005216389160241E-02,0.2785222700150757E-02, 
    0.2785222700144989E-02,0.2785222700178524E-02, 
    0.3753359578155891E-03,0.3753359578293081E-03, 
    0.3753359577958143E-03,0.2191182863239439E-02, 
    0.2191182863316704E-02,0.2191182863281961E-02, 
    0.2531026444894308E-02,0.2531026444936393E-02, 
    0.2531026444945725E-02,0.1445748928029682E-01, 
    0.1287820052356487E-03,0.1574864463007665E-01, 
    0.4621613088145150E-02,0.2792502892512972E-02, 
    0.1851564470805342E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule15 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE15 returns the rule of degree 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[] = { 
    0.2581521194788427,-.3924996680742822, 
    0.1343475485954617,0.1747936627600043, 
    -.3848787335398570,0.2100850707798483, 
    -.2686158385486640,0.3159539810469303, 
    -.4733814249823749E-01,0.5356808415765576E-02, 
    -.8659172184635908E-01,0.8123491343054050E-01, 
    -.8144631134845116,0.8470041241206219, 
    -.3254101063610144E-01,0.2249358441935477, 
    -.7323416768865754,0.5074058326930275, 
    -.5473570384341295,0.4843265718105462, 
    0.6303046662358983E-01,-.2378791442311287E-01, 
    -.2016252113668169,0.2254131257898949, 
    -.7391871819240758E-01,0.5113049221289290, 
    -.4373862039365108,-.2968004530699073, 
    0.1205383163903803,0.1762621366794993, 
    -.4241098255347354E-02,-.3585193417940622, 
    0.3627604400494130,0.4202386364031931E-02, 
    0.1599370643332776E-01,-.2019609279738776E-01, 
    -.5388551273867861E-01,-.7018353952881240, 
    0.7557209080267955,0.3646027806298292E-01, 
    -.6343658219832248,0.5979055439202614, 
    0.1258503466822884E-01,0.2393555213511634, 
    -.2519405560194015,0.3851625317792869E-02, 
    0.1534194382150093,-.1572710635328328, 
    -.3485301409937485E-02,0.1619818108911392, 
    -.1584965094812263,0.5304635114647127E-01, 
    0.3545913748848099,-.4076377260313040, 
    0.6653928546486268E-01,-.2403180998360646, 
    0.1737788143712072,-.1667348931062531, 
    -.2713747582479675,0.4381096513542055, 
    0.9493532822319652E-01,0.1868738963948708E-02, 
    -.9680406718713280E-01,-.5105981722468006E-01, 
    0.3251121848754134,-.2740523676507155, 
    -.2080245286773480E-01,-.1497129535134266, 
    0.1705154063812837,0.2921781238709439, 
    0.1381062592802215,-.4302843831511968, 
    -.8780905828000338E-01,0.2673000828923510, 
    -.1794910246123719,-.1944587456799732E-01, 
    -.2823011707456111,0.3017470453136293, 
    -.7313431320411609E-01,-.5236633756518274, 
    0.5967976888559381,-.3762507052448836, 
    0.2871315908267298,0.8911911441816055E-01, 
    -.3797896271838643E-01,-.4650879543817582, 
    0.5030669171001346,0.2750555493489867E-02, 
    -.7277556508895471,0.7250050953960729, 
    0.6123006138565600,-.1192649298638286, 
    -.4930356839927218,-.2645727220089428, 
    -.1730736170401170,0.4376463390490505, 
    -.6630962451131286E-01,-.1763224174379790, 
    0.2426320419492632,-.4200384093190587E-02, 
    -.5418782850881636,0.5460786691813587, 
    0.5394663228522688E-01,0.1507657818092707, 
    -.2047124140944948,0.1018171426999433, 
    -.5394039692520657,0.4375868265521267, 
    0.1245513301303478,-.6455557549708775, 
    0.5210044248405197,-.9855929555799212E-01, 
    -.7473211662573033,0.8458804618152944, 
    -.1755388020283085,-.5318861672552897, 
    0.7074249692836033,-.6267629157220122E-02, 
    0.6632219831012681E-01,-.6005456915289890E-01, 
    -.2023381108267559,0.1094534186561336, 
    0.9288469217061772E-01,-.6473499743579236E-01, 
    -.3234505788317799,0.3881855762675665, 
    -.1233025545427212,-.2864321270324415, 
    0.4097346815751876,0.1989648914526389, 
    -.5364226464888502,0.3374577550361906, 
    0.2178991630195039E-01,-.5426229008583486, 
    0.5208329845563997,-.2014523059829409, 
    -.1081005526244278,0.3095528586073633, 
    -.1169567835631198E-01,0.1260374294457714, 
    -.1143417510894922,0.6829681010950600E-01, 
    -.2897229857831844,0.2214261756736523, 
    -.4298175003068387,0.2110526002984961, 
    0.2187649000083270,0.6578307149076966, 
    -.1808749039318607,-.4769558109758388, 
    -.3704630736636612,0.3211543940711805, 
    0.4930867959248819E-01,-.1193005517808560, 
    -.3998990653549201,0.5191996171357768, 
    -.1504520555101911,0.4460532087610402, 
    -.2956011532508667,0.1975830054753616, 
    -.6709593039696905,0.4733762984943273, 
    -.7463137388742009E-02,-.4211815306853685, 
    0.4286446680740814,0.2049473200886011E-02, 
    -.2026603403979144,0.2006108671970308, 
    0.7152086010197926E-01,-.7652447634530471, 
    0.6937239033510662,0.2166839105383212, 
    0.3458968024421906,-.5625807129804981, 
    0.1279853939984014,0.3081799665774158, 
    -.4361653605757960,-.1577812374505563E-01, 
    -.6857246456300986,0.7015027693751511, 
    -.7105996230949263E-03,-.9227565022751169, 
    0.9234671018982127,-.1882955447783474E-01, 
    -.7344571519218671,0.7532867063996859, 
    -.3375217299349303,-.2041740513278847, 
    0.5416957812628146,0.4996795038021838, 
    0.1342360979189166,-.6339156017211199, 
    0.7596408799506882,0.8346506611711336E-01, 
    -.8431059460678011,-.8820027427037123, 
    0.8657493528015160,0.1625338990220750E-01, 
    0.7067742334596752,-.2761780760336056, 
    -.4305961574260422,-.9186459446440826E-01, 
    -.1188760024130717,0.2107405968774760, 
    0.5405133012525005,-.2620825757773786, 
    -.2784307254751027,0.1458059413495244, 
    0.3829161242112978,-.5287220655607965, 
    -.1823158001392015E-13,-.1699106971810907E-13, 
    0.1569922552833598E-12,0.5359985656924711E-13 };
  double ys[] = { 
    0.3041753823660573,0.7147860232646702E-01, 
    -.3756539846925305,0.3435025125820099, 
    -.2037550392031693E-01,-.3231270086616979, 
    -.2097468053232621,-.1277547373803833, 
    0.3375015427036143,0.9689475305431051E-01, 
    -.4380824435589394E-01,-.5308650869840935E-01, 
    -.5078056203162217,-.4514429365648450, 
    0.9592485568810678,0.7157685583829504, 
    -.1630841238981528,-.5526844344847958, 
    -.2432354197382121,-.3524073903550577, 
    0.5956428100933059,0.2465506989117151, 
    -.1438762876493171,-.1026744112624022, 
    -.5477270770116898,0.2098480507360208, 
    0.3378790262756754,0.3217216264418476E-01, 
    -.2731228135353417,0.2409506508911702, 
    0.4164310762083722,-.2118884369332560, 
    -.2045426392751013,-.2089419032777848E-01, 
    0.1408646851164702E-01,0.6807721816119913E-02, 
    0.8415205240779101,-.4674264849665979, 
    -.3740940391113035,0.7114522048190210, 
    -.3241505753779200,-.3873016294411025, 
    -.2836499225217341,0.1527239209910285, 
    0.1309260015306806,-.1793772448187649, 
    0.9302422778043601E-01,0.8635301703831089E-01, 
    -.1850282445363872,0.8949576270733584E-01, 
    0.9553248182902167E-01,-.4400731765980834, 
    0.2659760759699537,0.1740971006281070, 
    0.2390789648881394,-.6191477088183470E-01, 
    -.1771641940062932,0.4096210148696401, 
    -.3492071605621237,-.6041385430752706E-01, 
    -.5696877119302565E-01,0.1107007915544150, 
    -.5373202036138137E-01,-.3459278156898313, 
    0.1287448090157230,0.2171830066740836, 
    0.1848839297874223,-.1104574175380845, 
    -.7442651224930893E-01,-.3281604904126481, 
    0.4171139229086880,-.8895343249598349E-01, 
    -.2579549661893832,0.5293260794183284E-01, 
    0.2050223582475361,0.3372003947615109, 
    -.1854408187554140,-.1517595760060760, 
    0.6468984972100610,-.3867854217281152, 
    -.2601130754819361,-.1143225565573846, 
    -.2686813906552311,0.3830039472125820, 
    0.5589644756672905,-.3123729843571679, 
    -.2465914913101748,0.8387518079361266, 
    -.4169938530361755,-.4217579548999488, 
    -.2157966455115084,0.6381662091083642, 
    -.4223695635968309,0.3525993310475737, 
    -.4054263639319154,0.5282703288436957E-01, 
    0.2418834699054176,-.1783675542949067, 
    -.6351591561047869E-01,0.6281322404142361, 
    -.3177037595374698,-.3104284808767618, 
    -.2052354320960713,0.1493368700556973, 
    0.5589856204040973E-01,0.5640658989533277, 
    -.1938567173577759,-.3702091815955630, 
    0.6735138338400053,-.2288923009519804, 
    -.4446215328880427,0.9198353888410614, 
    -.5452725481528549,-.3745628406882141, 
    0.7155166182904237,-.5097793710516186, 
    -.2057372472388190,-.7296366071409344E-01, 
    0.3105390428539375E-01,0.4190975642870679E-01, 
    -.9565958696542072E-02,-.1704469647814408, 
    0.1800129234779798,0.4108633257116824, 
    -.2614938151491584,-.1493695105625245, 
    0.4019320943505255,-.3077491917607414, 
    -.9418290258973967E-01,0.5045350850600017, 
    -.7995889207082835E-01,-.4245761929892107, 
    0.6139865417154907,-.2881226497939412, 
    -.3258638919215674,0.2411323094026556, 
    -.2950289693335045,0.5389665993088002E-01, 
    -.1387829845896145,0.5926273772376190E-01, 
    0.7952024686583803E-01,0.2951121059631410, 
    -.8840928042931671E-01,-.2067028255338328, 
    0.4452698313519278E-02,-.3744592234135989, 
    0.3700065251000597,-.1709423913837525, 
    0.6551693061916171,-.4842269148078915, 
    -.1569501964323080,-.2423553347406391, 
    0.3993055311729469,0.5306418717478397, 
    -.3686382444016378,-.1620036273461921, 
    -.4281943455532308,0.8380187065320079E-01, 
    0.3443924749000115,0.6606824681258688, 
    -.1592293319651774,-.5014531361606821, 
    0.4906473846181765,-.2517869588796617, 
    -.2388604257384930,0.2328287402613732, 
    -.1146394742743420,-.1181892659870265, 
    0.8423359525185951,-.3592290945104593, 
    -.4831068580081231,-.5245097380820163, 
    0.4499086401585516,0.7460109792348690E-01, 
    -.4297479750019432,0.3257125900169633, 
    0.1040353849849871,0.8009161214805088, 
    -.4141223167275365,-.3867938047529831, 
    0.1065917694853711E+01,-.5335742447523749, 
    -.5323434501013379,0.8589493170871557, 
    -.4457815310633328,-.4131677860238335, 
    0.4306281486266579,-.5076164667662489, 
    0.7698831813957928E-01,-.4434925905656182, 
    0.6544814393259342,-.2109888487602922, 
    -.5349560233082641,0.9253463114445929, 
    -.3903902881363264,-.4904567228554135, 
    -.5186084199612601,0.1009065142816689E+01, 
    -.8915332085968262E-01,0.6566611013461990, 
    -.5675077804865333,0.1903042323297749, 
    -.1747091886794188,-.1559504365035713E-01, 
    -.9438608628742183E-02,0.4728175542824327, 
    -.4633789456536745,-.5263345542684673, 
    0.3894389263656265,0.1368956279028533, 
    -.1738343488671503E-14,-.2490206331023665E-15, 
    0.1028727833546315E-12,0.7810443937532853E-13 };
  double zs[] = { 
    0.1149415369298274,0.1149415369297980, 
    0.1149415369297721,-.1425141964887960, 
    -.1425141964887534,-.1425141964887875, 
    -.3897391650202160,-.3897391650202120, 
    -.3897391650202178,0.7425283167186030, 
    0.7425283167186189,0.7425283167186136, 
    -.2856112242082660,-.2856112242082698, 
    -.2856112242082572,-.3891867546444293, 
    -.3891867546444389,-.3891867546444343, 
    0.1847325582324109,0.1847325582324027, 
    0.1847325582324080,0.5802704682522071, 
    0.5802704682522261,0.5802704682522338, 
    -.3914413911049179,-.3914413911049140, 
    -.3914413911049208,0.4061500673930710, 
    0.4061500673931097,0.4061500673930845, 
    0.5808632391400850,0.5808632391400872, 
    0.5808632391401026,0.1122300778479799E+01, 
    0.1122300778479790E+01,0.1122300778479803E+01, 
    -.1702382601410010,-.1702382601409925, 
    -.1702382601409868,-.2928517214460706E-01, 
    -.2928517214460719E-01,-.2928517214460035E-01, 
    -.2378639271913669E-01,-.2378639271913052E-01, 
    -.2378639271915441E-01,-.3876002482780669, 
    -.3876002482780664,-.3876002482780634, 
    -.1573276400467179,-.1573276400467495, 
    -.1573276400467503,-.9846677659403662E-01, 
    -.9846677659401415E-01,-.9846677659402990E-01, 
    0.4195749508887953,0.4195749508888120, 
    0.4195749508888156,0.1698706478280650, 
    0.1698706478280417,0.1698706478280619, 
    0.1001045633396894E+01,0.1001045633396896E+01, 
    0.1001045633396904E+01,0.1278252986614525, 
    0.1278252986614319,0.1278252986614156, 
    0.6699851292944298E-01,0.6699851292952046E-01, 
    0.6699851292943053E-01,0.1888603052005784, 
    0.1888603052005364,0.1888603052005152, 
    0.4162319117452538,0.4162319117452708, 
    0.4162319117452612,-.1429377419107346, 
    -.1429377419107034,-.1429377419107302, 
    0.7526757476604347E-01,0.7526757476604869E-01, 
    0.7526757476604784E-01,0.5502259745637508E-01, 
    0.5502259745639024E-01,0.5502259745637604E-01, 
    -.3783645730290734,-.3783645730290731, 
    -.3783645730290733,-.3557858186989302E-01, 
    -.3557858186987815E-01,-.3557858186989157E-01, 
    -.1180094770935530,-.1180094770935434, 
    -.1180094770935428,-.1637480027557202, 
    -.1637480027557221,-.1637480027557144, 
    0.2977454860475704,0.2977454860475816, 
    0.2977454860475674,0.2596712604598457, 
    0.2596712604598533,0.2596712604598482, 
    0.2271508926955859,0.2271508926955554, 
    0.2271508926955697,-.2779926178974129, 
    -.2779926178973848,-.2779926178974117, 
    -.3823432777224465,-.3823432777224400, 
    -.3823432777224473,-.3775386322444910, 
    -.3775386322444932,-.3775386322444943, 
    -.2829688627011475,-.2829688627011491, 
    -.2829688627011544,0.9266659902735654, 
    0.9266659902735608,0.9266659902735637, 
    0.7042408708164912,0.7042408708164883, 
    0.7042408708164785,0.4116884047492347, 
    0.4116884047492301,0.4116884047492345, 
    0.1190964787818936E-01,0.1190964787820830E-01, 
    0.1190964787819720E-01,-.9824886296290353E-01, 
    -.9824886296290938E-01,-.9824886296292371E-01, 
    -.1210076932793752,-.1210076932793519, 
    -.1210076932793500,-.2970323317075053, 
    -.2970323317075050,-.2970323317075011, 
    0.5828065214604544,0.5828065214604353, 
    0.5828065214604685,-.3017876570907457, 
    -.3017876570907429,-.3017876570907401, 
    -.3862244277990527,-.3862244277990511, 
    -.3862244277990521,-.3804346304483138, 
    -.3804346304483134,-.3804346304483132, 
    0.4614518667838061,0.4614518667838234, 
    0.4614518667838163,-.2664567796530529, 
    -.2664567796530461,-.2664567796530568, 
    -.3834440210666236,-.3834440210666228, 
    -.3834440210666243,-.2809660745834901, 
    -.2809660745835049,-.2809660745834982, 
    0.1843225423322873,0.1843225423322890, 
    0.1843225423322869,0.7949931745648482, 
    0.7949931745648681,0.7949931745648552, 
    -.2296019504049502,-.2296019504049503, 
    -.2296019504049410,-.3768873297184356, 
    -.3768873297184376,-.3768873297184345, 
    -.2779742390429553,-.2779742390429600, 
    -.2779742390429549,-.3088199065240608, 
    -.3088199065240601,-.3088199065240508, 
    -.3819799973245133,-.3819799973245119, 
    -.3819799973245152,-.4015361795222331, 
    -.4015361795222343,-.4015361795222284, 
    -.3045707129688655,-.3045707129688613, 
    -.3045707129688553,-.3905755491124832E-01, 
    -.3905755491127184E-01,-.3905755491126780E-01, 
    -.3726108784532423,-.3726108784532439, 
    -.3726108784532431,-.2575990026123318, 
    -.2575990026123464,-.2575990026123427, 
    -.3954674471205766,-.3954674471205773, 
    -.3954674471205850,0.7158508859841455, 
    0.7158508859841465,0.7158508859841323, 
    -.9143187100549784E-01,-.9143187100549642E-01, 
    -.9143187100549681E-01,-.2639610543972574, 
    -.2639610543972473,-.2639610543972535, 
    -.3003943802989339,0.4275416212259733, 
    0.1454645115620674,-.1072798570797851 };
  double ws[] = { 
    0.2565379947273549E-02,0.2565379947271849E-02, 
    0.2565379947272177E-02,0.6723785911414816E-02, 
    0.6723785911411640E-02,0.6723785911413937E-02, 
    0.3085562545744757E-02,0.3085562545745181E-02, 
    0.3085562545744572E-02,0.4446982777989968E-02, 
    0.4446982777989351E-02,0.4446982777990241E-02, 
    0.1185105459270273E-02,0.1185105459270278E-02, 
    0.1185105459270443E-02,0.8364165694309715E-03, 
    0.8364165694306848E-03,0.8364165694306367E-03, 
    0.1867492375962935E-02,0.1867492375963348E-02, 
    0.1867492375962387E-02,0.4264615556690722E-02, 
    0.4264615556689985E-02,0.4264615556689865E-02, 
    0.1036750961119152E-02,0.1036750961119384E-02, 
    0.1036750961118855E-02,0.2897017453521770E-02, 
    0.2897017453522957E-02,0.2897017453522235E-02, 
    0.9931774636746942E-03,0.9931774636749676E-03, 
    0.9931774636745359E-03,0.2968267531327685E-03, 
    0.2968267531331858E-03,0.2968267531325610E-03, 
    0.2082824878589059E-02,0.2082824878589280E-02, 
    0.2082824878589393E-02,0.3532360099744144E-02, 
    0.3532360099743258E-02,0.3532360099742647E-02, 
    0.9852850436647919E-02,0.9852850436649995E-02, 
    0.9852850436648092E-02,0.4465622015302662E-02, 
    0.4465622015302984E-02,0.4465622015303046E-02, 
    0.1064447126836857E-01,0.1064447126836189E-01, 
    0.1064447126836531E-01,0.4362160716691730E-02, 
    0.4362160716692162E-02,0.4362160716691511E-02, 
    0.7013986090707679E-02,0.7013986090706681E-02, 
    0.7013986090707788E-02,0.4089047890127443E-02, 
    0.4089047890127423E-02,0.4089047890126841E-02, 
    0.8331696920946515E-03,0.8331696920946096E-03, 
    0.8331696920947063E-03,0.5692536319085584E-02, 
    0.5692536319084854E-02,0.5692536319085869E-02, 
    0.1193994326736248E-01,0.1193994326737174E-01, 
    0.1193994326736250E-01,0.4894866890556946E-02, 
    0.4894866890557526E-02,0.4894866890558155E-02, 
    0.4559845100025953E-02,0.4559845100026132E-02, 
    0.4559845100026570E-02,0.1112520926249036E-01, 
    0.1112520926249222E-01,0.1112520926248886E-01, 
    0.2944181342559152E-02,0.2944181342559592E-02, 
    0.2944181342559789E-02,0.8953416346063897E-02, 
    0.8953416346061954E-02,0.8953416346063781E-02, 
    0.3952683058023190E-02,0.3952683058023330E-02, 
    0.3952683058022607E-02,0.1439142717612964E-02, 
    0.1439142717612978E-02,0.1439142717612582E-02, 
    0.5039224704625283E-02,0.5039224704625251E-02, 
    0.5039224704625310E-02,0.7617567889070472E-02, 
    0.7617567889069997E-02,0.7617567889070361E-02, 
    0.1054110867291372E-01,0.1054110867291354E-01, 
    0.1054110867291408E-01,0.1741617181694419E-02, 
    0.1741617181694472E-02,0.1741617181694386E-02, 
    0.1100640373498494E-01,0.1100640373498203E-01, 
    0.1100640373498438E-01,0.6926338144663090E-02, 
    0.6926338144663959E-02,0.6926338144663389E-02, 
    0.2939871895670393E-02,0.2939871895670890E-02, 
    0.2939871895670184E-02,0.1098371982889171E-02, 
    0.1098371982888983E-02,0.1098371982889088E-02, 
    0.2867324218775673E-02,0.2867324218775326E-02, 
    0.2867324218775453E-02,0.2627584056743849E-02, 
    0.2627584056743863E-02,0.2627584056743943E-02, 
    0.2450903116528004E-02,0.2450903116528749E-02, 
    0.2450903116529087E-02,0.3781264707127557E-02, 
    0.3781264707127677E-02,0.3781264707128319E-02, 
    0.9717660774153364E-02,0.9717660774155940E-02, 
    0.9717660774153394E-02,0.5844696436169811E-02, 
    0.5844696436168778E-02,0.5844696436168716E-02, 
    0.7798609566561363E-02,0.7798609566560966E-02, 
    0.7798609566562299E-02,0.9537993724461966E-02, 
    0.9537993724462559E-02,0.9537993724462396E-02, 
    0.6470868601806971E-02,0.6470868601806687E-02, 
    0.6470868601806209E-02,0.9788338589999005E-02, 
    0.9788338589999012E-02,0.9788338589999477E-02, 
    0.4149912310618085E-02,0.4149912310618563E-02, 
    0.4149912310618460E-02,0.2801823664802904E-02, 
    0.2801823664802910E-02,0.2801823664803036E-02, 
    0.3473047605746596E-02,0.3473047605746489E-02, 
    0.3473047605747130E-02,0.8146696772355114E-02, 
    0.8146696772356077E-02,0.8146696772355570E-02, 
    0.4104346629768971E-02,0.4104346629769071E-02, 
    0.4104346629768999E-02,0.3673293016416287E-02, 
    0.3673293016415687E-02,0.3673293016415978E-02, 
    0.7385814346097036E-02,0.7385814346096885E-02, 
    0.7385814346097676E-02,0.1952835342529582E-02, 
    0.1952835342529196E-02,0.1952835342529467E-02, 
    0.2956577367937866E-02,0.2956577367938132E-02, 
    0.2956577367937872E-02,0.2429727102722713E-02, 
    0.2429727102722471E-02,0.2429727102722800E-02, 
    0.7838569405690889E-02,0.7838569405690287E-02, 
    0.7838569405690948E-02,0.4786830078294293E-02, 
    0.4786830078294275E-02,0.4786830078294352E-02, 
    0.6464879373388097E-03,0.6464879373388156E-03, 
    0.6464879373387403E-03,0.1094759729336998E-02, 
    0.1094759729336943E-02,0.1094759729337179E-02, 
    0.4183486775030887E-02,0.4183486775030987E-02, 
    0.4183486775031043E-02,0.1388382329059297E-02, 
    0.1388382329058254E-02,0.1388382329057900E-02, 
    0.1547261679853963E-02,0.1547261679853909E-02, 
    0.1547261679853913E-02,0.5245559039887955E-03, 
    0.5245559039888078E-03,0.5245559039885548E-03, 
    0.4714073461717449E-03,0.4714073461717063E-03, 
    0.4714073461712849E-03,0.1734297963526248E-02, 
    0.1734297963526404E-02,0.1734297963527124E-02, 
    0.2194380264311005E-02,0.2194380264311055E-02, 
    0.2194380264311861E-02,0.1505920846848116E-02, 
    0.1505920846847692E-02,0.1505920846848549E-02, 
    0.1092191482369307E-01,0.1402579642206070E-01, 
    0.9137907285788647E-02,0.8807144129412933E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void tetrahedron_arbq ( int degree, int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_ARBQ returns a quadrature rule for the tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the degree of the quadrature rule.
//    1 <= DEGREE <= 15.
//
//    Input, int N, the number of nodes.
//    This can be determined by a call to TETRAHEDRON_ARBQ_SIZE(DEGREE).
//
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double volume;
  double w_sum;

  if ( degree == 1 )
  {
    rule01 ( n, x, w );
  }
  else if ( degree == 2 )
  {
    rule02 ( n, x, w );
  }
  else if ( degree == 3 )
  {
    rule03 ( n, x, w );
  }
  else if ( degree == 4 )
  {
    rule04 ( n, x, w );
  }
  else if ( degree == 5 )
  {
    rule05 ( n, x, w );
  }
  else if ( degree == 6 )
  {
    rule06 ( n, x, w );
  }
  else if ( degree == 7 )
  {
    rule07 ( n, x, w );
  }
  else if ( degree == 8 )
  {
    rule08 ( n, x, w );
  }
  else if ( degree == 9 )
  {
    rule09 ( n, x, w );
  }
  else if ( degree == 10 )
  {
    rule10 ( n, x, w );
  }
  else if ( degree == 11 )
  {
    rule11 ( n, x, w );
  }
  else if ( degree == 12 )
  {
    rule12 ( n, x, w );
  }
  else if ( degree == 13 )
  {
    rule13 ( n, x, w );
  }
  else if ( degree == 14 )
  {
    rule14 ( n, x, w );
  }
  else if ( degree == 15 )
  {
    rule15 ( n, x, w );
  }
  else
  {
    cerr << "\n";
    cerr << "TETRAHEDRON_ARBQ - Fatal error\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  w_sum = r8vec_sum ( n, w );
  volume = sqrt ( 8.0 ) / 3.0;

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * volume / w_sum;
  }

  return;
}
//****************************************************************************80

void tetrahedron_arbq_gnuplot ( int n, double x[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_ARBQ_GNUPLOT: plot of a quadrature rule for the tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double N, the number of nodes.
//
//    Input, double X[3*N], the coordinates of the nodes.
//
//    Input, string HEADER, a string to be used to identify
//    the files created.
//
{
  string command_filename;
  ofstream command_unit;
  int i;
  int j;
  int l;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  string vertex_filename;
  ofstream vertex_unit;
//
//  Create the vertex file.
//
  tetrahedron_ref ( v1, v2, v3, v4 );

  vertex_filename = header + "_vertices.txt";
  vertex_unit.open ( vertex_filename.c_str ( ) );
  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << "\n";
  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << "\n";
  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << "\n";
  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit.close ( );
  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Create node file.
//
  node_filename = header + "_nodes.txt";
  node_unit.open ( node_filename.c_str ( ) );
  for ( j = 0; j < n; j++ )
  {
    node_unit << x[0+j*3] << "  "
              << x[1+j*3] << "  "
              << x[2+j*3] << "\n";
  }
  node_unit.close ( );
  cout << "  Created node file '" << node_filename << "'\n";
//
//  Create graphics command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  plot_filename = header + ".png";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set view equal xyz\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "splot '" << vertex_filename << "' with lines lw 3, \\\n";
  command_unit << "      '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

int tetrahedron_arbq_size ( int degree )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_ARBQ_SIZE returns the size of quadrature rule for a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired degree of exactness.
//    1 <= DEGREE <= 15.
//
//    Output, int TETRAHEDRON_ARBQ_SIZE, the number of points in the
//    corresponding rule.
//
{
  int n;
  const int n_save[15] = {
      1,   4,   6,  11,  14, 
     23,  31,  44,  57,  74, 
     95, 122, 146, 177, 214 };

  if ( degree < 1 || 15 < degree ) 
  {
    cerr << "\n";
    cerr << "TETRAHEDRON_ARBQ_SIZE - Fatal error!\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  n = n_save[degree-1];

  return n;
}
//****************************************************************************80

void tetrahedron_ref ( double v1[], double v2[], double v3[], double v4[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_REF returns the vertices of the reference tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Output, double V1[3], V2[3], V3[3], V4[3], the vertices.
//
{

  v1[0] = - 1.0;
  v1[1] = - 1.0 / sqrt ( 3.0 );
  v1[2] = - 1.0 / sqrt ( 6.0 );

  v2[0] =   0.0;
  v2[1] =   2.0 / sqrt ( 3.0 );
  v2[2] = - 1.0 / sqrt ( 6.0 );

  v3[0] =   1.0;
  v3[1] = - 1.0 / sqrt ( 3.0 );
  v3[2] = - 1.0 / sqrt ( 6.0 );

  v4[0] =   0.0;
  v4[1] =   0.0;
  v4[2] =   3.0 / sqrt ( 6.0 );

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
