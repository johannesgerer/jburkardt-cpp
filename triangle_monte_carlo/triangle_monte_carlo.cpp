# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_monte_carlo.hpp"

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

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are halfway along the sides of
//    the physical triangle.
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
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                 + t[i+1*2] *       + ref[0+j*2]
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
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

double triangle_area ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the area of a triangle in 2D.
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
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.  AREA will
//    be nonnegative.
//
{
  double area;

  area = fabs ( 0.5 * (
    t[0+0*2] * ( t[1+2*2] - t[1+1*2] ) +
    t[0+1*2] * ( t[1+0*2] - t[1+2*2] ) +
    t[0+2*2] * ( t[1+1*2] - t[1+0*2] ) ) );

  return area;
}
//****************************************************************************80

double *triangle_integrand_01 ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_01 evaluates 1 integrand function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP[F_NUM*P_NUM], the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = 1.0;
  }
  return fp;
}
//****************************************************************************80

double *triangle_integrand_02 ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_02 evaluates 2 integrand functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP[F_NUM*P_NUM], the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = p[0+j*2];
    fp[1+j*f_num] = p[1+j*2];
  }

  return fp;
}
//****************************************************************************80

double *triangle_integrand_03 ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_03 evaluates 3 integrand functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP[F_NUM*P_NUM], the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = p[0+j*2] * p[0+j*2];
    fp[1+j*f_num] = p[0+j*2] * p[1+j*2];
    fp[2+j*f_num] = p[1+j*2] * p[1+j*2];
  }

  return fp;
}
//****************************************************************************80

double *triangle_integrand_04 ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_04 evaluates 4 integrand functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP[F_NUM*P_NUM], the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = p[0+j*2] * p[0+j*2] * p[0+j*2];
    fp[1+j*f_num] = p[0+j*2] * p[0+j*2] * p[1+j*2];
    fp[2+j*f_num] = p[0+j*2] * p[1+j*2] * p[1+j*2];
    fp[3+j*f_num] = p[1+j*2] * p[1+j*2] * p[1+j*2];
  }

  return fp;
}
//****************************************************************************80

double *triangle_integrand_05 ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_05 evaluates 5 integrand functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP[F_NUM*P_NUM], the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = p[0+j*2] * p[0+j*2] * p[0+j*2] * p[0+j*2];
    fp[1+j*f_num] = p[0+j*2] * p[0+j*2] * p[0+j*2] * p[1+j*2];
    fp[2+j*f_num] = p[0+j*2] * p[0+j*2] * p[1+j*2] * p[1+j*2];
    fp[3+j*f_num] = p[0+j*2] * p[1+j*2] * p[1+j*2] * p[1+j*2];
    fp[4+j*f_num] = p[1+j*2] * p[1+j*2] * p[1+j*2] * p[1+j*2];
  }
  return fp;
}
//****************************************************************************80

double *triangle_monte_carlo ( double t[], int p_num, int f_num,
  double *triangle_unit_sample ( int p_num, int *seed ),
  double *triangle_integrand ( int p_num, double p[], int f_num ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
//
//  Discussion:
//
//    The function f(x,y) is to be integrated over a triangle T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, int P_NUM, the number of sample points.
//
//    Input, int F_NUM, the number of functions to integrate.
//
//    Input, external TRIANGLE_UNIT_SAMPLE, the sampling routine.
//
//    Input, external TRIANGLE_INTEGRAND, the integrand routine.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, dobule TRIANGLE_MONTE_CARLO[F_NUM], the approximate integrals.
//
{
  double area;
  double *fp;
  double fp_sum;
  int i;
  int j;
  double *p;
  double *p2;
  double *result;

  area = triangle_area ( t );

  p = triangle_unit_sample ( p_num, seed );

  p2 = new double[2*p_num];

  reference_to_physical_t3 ( t, p_num, p, p2 );

  fp = triangle_integrand ( p_num, p2, f_num );

  result = new double[f_num];

  for ( i = 0; i < f_num; i++ )
  {
    fp_sum = 0.0;
    for ( j = 0; j < p_num; j++ )
    {
      fp_sum = fp_sum + fp[i+j*f_num];
    }
    result[i] = area * fp_sum / ( double ) ( p_num );
  }

  delete [] fp;
  delete [] p;
  delete [] p2;

  return result;
}
//****************************************************************************80

double *triangle_unit_sample_01 ( int p_num, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SAMPLE_01 selects points from the unit triangle.
//
//  Discussion:
//
//    The unit triangle has vertices (1,0), (0,1), (0,0).
//
//    Any point in the unit simplex CAN be chosen by this algorithm.
//
//    However, the points that are chosen tend to be clustered near
//    the centroid.
//
//    This routine is supplied as an example of "bad" sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, double TRIANGLE_UNIT_SAMPLE_01[2*P_NUM], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;

  x = new double[2*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    e = r8vec_uniform_01_new ( 3, seed );

    e_sum = r8vec_sum ( 3, e );

    for ( i = 0; i < 2; i++ )
    {
      x[i+j*2] = e[i] / e_sum;
    }

    delete [] e;
  }

  return x;
}
//****************************************************************************80

double *triangle_unit_sample_02 ( int p_num, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SAMPLE_02 selects points from the unit triangle.
//
//  Discussion:
//
//    The unit triangle has vertices (1,0), (0,1), (0,0).
//
//    The sampling is uniform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double TRIANGLE_UNIT_SAMPLE_02[2*P_NUM], the points.
//
{
  int j;
  double *r;
  double *x;

  x = new double[2*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }
    x[0+j*2] = r[0];
    x[1+j*2] = r[1];

    delete [] r;
  }

  return x;
}
//****************************************************************************80

double *triangle_unit_sample_03 ( int p_num, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SAMPLE_03 selects points from the unit triangle.
//
//  Discussion:
//
//    The unit triangle has vertices (1,0), (0,1), (0,0).
//
//    This routine uses Turk's rule 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double TRIANGLE_UNIT_SAMPLE_03[2*P_NUM], the points.
//
{
  double a;
  double b;
  double c;
  int j;
  double *r;
  double total;
  double *x;

  x = new double[2*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    r[1] = sqrt ( r[1] );

    a = 1.0 - r[1];
    b = ( 1.0 - r[0] ) * r[1];
    c = r[0] * r[1];

    x[0+j*2] = a;
    x[1+j*2] = b;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

double *triangle_unit_sample_04 ( int p_num, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SAMPLE_04 selects points from the unit triangle.
//
//  Discussion:
//
//    The unit triangle has vertices (1,0), (0,1), (0,0).
//
//    The sampling is uniform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double TRIANGLE_UNIT_SAMPLE_04[2*P_NUM], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;
//
//  The construction begins by sampling DIM_NUM+1 points from the
//  exponential distribution with parameter 1.
//
  x = new double[2*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    e = r8vec_uniform_01_new ( 3, seed );

    for ( i = 0; i <= 2; i++ )
    {
      e[i] = - log ( e[i] );
    }

    e_sum = r8vec_sum ( 3, e );

    for ( i = 0; i < 2; i++ )
    {
      x[i+2*j] = e[i] / e_sum;
    }
    delete [] e;
  }

  return x;
}
