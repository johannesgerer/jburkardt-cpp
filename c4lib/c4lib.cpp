# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <string>

using namespace std;

# include "c4lib.hpp"

//****************************************************************************80

float c4_abs ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ABS returns the absolute value of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_ABS, the magnitude of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

float c4_argument ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ARGUMENT returns the argument of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose argument is desired.
//
//    Output, float C4_ARGUMENT, the argument of X.
//
{
  float value;

  if ( imag ( x ) == 0.0 && real ( x ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = atan2 ( imag ( x ), real ( x ) );
  }

  return value;
}
//****************************************************************************80

complex <float> c4_cube_root ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_CUBE_ROOT returns the principal cube root of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the number whose cube root is desired.
//
//    Output, complex <float> C4_CUBE_ROOT, the cube root of X.
//
{
  float argument;
  float magnitude;
  complex <float> value;

  argument = c4_argument ( x );
  magnitude = c4_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = pow ( magnitude, ( float ) ( 1.0 / 3.0 ) ) 
      * complex <float> ( cos ( argument / 3.0 ), sin ( argument / 3.0 ) );
  }

  return value;
}
//****************************************************************************80

complex <float> c4_i ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_I returns the value of the imaginary unit, i as a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_I, the value of complex i.
//
{
  complex <float> value;

  value = complex <float> ( 0.0, 1.0 );

  return value;
}
//****************************************************************************80

bool c4_le_l1 ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_L1 := X <= Y for C4 values, and the L1 norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L1 norm can be defined here as:
//
//      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_L1, is TRUE if X <= Y.
//
{
  bool value;

  if ( r4_abs ( real ( x ) ) + r4_abs ( imag ( x ) ) <= 
       r4_abs ( real ( y ) ) + r4_abs ( imag ( y ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool c4_le_l2 ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_L2 := X <= Y for C4 values, and the L2 norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L2 norm can be defined here as:
//
//      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_L2, is TRUE if X <= Y.
//
{
  bool value;

  if ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) <= 
       pow ( real ( y ), 2 ) + pow ( imag ( y ), 2 ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool c4_le_li ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_LI := X <= Y for C4 values, and the L-oo norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L-oo norm can be defined here as:
//
//      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_LI, is TRUE if X <= Y.
//
{
  bool value;

  if ( r4_max ( r4_abs ( real ( x ) ), r4_abs ( imag ( x ) ) ) <= 
       r4_max ( r4_abs ( real ( y ) ), r4_abs ( imag ( y ) ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

float c4_magnitude ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_MAGNITUDE returns the magnitude of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_MAGNITUDE, the magnitude of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

float c4_norm_l1 ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_L1 evaluates the L1 norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie along diamonds centered at (0,0).
//
//    The L1 norm can be defined here as:
//
//      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_L1, the norm of X.
//
{
  float value;

  value = r4_abs ( real ( x ) ) + r4_abs ( imag ( x ) );

  return value;
}
//****************************************************************************80

float c4_norm_l2 ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_L2 evaluates the L2 norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie on circles centered at (0,0).
//
//    The L2 norm can be defined here as:
//
//      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_L2, the 2-norm of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 )
               + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

float c4_norm_li ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_LI evaluates the L-oo norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie along squares whose centers are at (0,0).
//
//    The L-oo norm can be defined here as:
//
//      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_LI, the L-oo norm of X.
//
{
  float value;

  value = r4_max ( r4_abs ( real ( x ) ), r4_abs ( imag ( x ) ) );

  return value;
}
//****************************************************************************80

complex <float> c4_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORMAL_01 returns a unit pseudonormal C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, complex <float> C4_NORMAL_01, a unit pseudonormal value.
//
{
  float pi = 3.141592653589793;
  float v1;
  float v2;
  complex <float> value;
  float x_c;
  float x_r;

  v1 = r4_uniform_01 ( seed );
  v2 = r4_uniform_01 ( seed );

  x_r = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * pi * v2 );
  x_c = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * pi * v2 );

  value = complex <float> ( x_r, x_c );

  return value;
}
//****************************************************************************80

complex <float> c4_one ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ONE returns the value of complex 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_ONE, the value of complex 1.
//
{
  complex <float> value;

  value = complex <float> ( 1.0, 0.0);

  return value;
}
//****************************************************************************80

void c4_print ( complex <float> a, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4_PRINT prints a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> A, the value to be printed.
//
//    Input, string TITLE, a title.
//
{
  cout << title
       << "  ( " << setw(14) << real ( a )
       << ", "   << setw(14) << imag ( a ) << " )\n";

  return;
}
//****************************************************************************80

complex <float> c4_sqrt ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SQRT returns the principal square root of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the number whose square root is desired.
//
//    Output, complex <float> C4_SQRT, the square root of X.
//
{
  float argument;
  float magnitude;
  complex <float> value;

  argument = c4_argument ( x );
  magnitude = c4_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = sqrt ( magnitude ) 
      * complex <float> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
  }

  return value;
}
//****************************************************************************80

void c4_swap ( complex <float> *x, complex <float> *y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SWAP swaps two C4's.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, complex <float> *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  complex <float> z;

   z = *x; 
  *x = *y;
  *y =  z;

  return;
}
//****************************************************************************80

complex <float> c4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01 returns a unit pseudorandom C4.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
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
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
//
{
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;
  complex <float> value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  theta = 2.0 * pi *
    ( ( float ) ( *seed ) * 4.656612875E-10 );

  value = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

complex <float> c4_zero ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ZERO returns the value of 0 as a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_ZERO, the value of complex 0.
//
{
  complex <float> value;

  value = complex <float> ( 0.0, 0.0 );

  return value;
}
//****************************************************************************80

void c4mat_copy ( int m, int n, complex <float> a1[], complex <float> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_COPY copies one C4MAT to another.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex <float> values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> A1[M*N], the matrix to be copied.
//
//    Output, complex <float> A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

complex <float> *c4mat_copy_new ( int m, int n, complex <float> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_COPY_NEW copies one C4MAT to a "new" C4MAT.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex <float> values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> A1[M*N], the matrix to be copied.
//
//    Output, complex <float> C4MAT_COPY_NEW[M*N], the copy of A1.
//
{
  complex <float> *a2;
  int i;
  int j;

  a2 = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

complex <float> *c4mat_identity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_IDENTITY sets a C4MAT to the identity.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_IDENTITY[N*N], the matrix.
//
{
  complex <float> *a;
  int i;
  int j;

  a = new complex <float> [n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = complex <float> ( 1.0, 0.0 );
      }
      else
      {
        a[i+j*n] = complex <float> ( 0.0, 0.0 );
      }
    }
  }
  return a;
}
//****************************************************************************80

complex <float> *c4mat_indicator_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_INDICATOR_NEW returns the C4MAT indicator matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, complex <float> C4MAT_INDICATOR_NEW[M*N], the matrix.
//
{
  complex <float> *a;
  int i;
  int j;

  a = new complex <float> [m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <float> ( i, j );
    }
  }
  return a;
}
//****************************************************************************80

void c4mat_nint ( int m, int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_NINT rounds the entries of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//
//    Input/output, complex <float> A[M*N], the matrix to be NINT'ed.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <float> ( 
        r4_nint ( real ( a[i+j*m] ) ), 
        r4_nint ( imag ( a[i+j*m] ) ) );
    }
  }
  return;
}
//****************************************************************************80

void c4mat_print ( int m, int n, complex <float> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT prints a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  c4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT_SOME prints some of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
  complex <float> c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
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
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= i4_min ( jhi, n ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << "     " << setw(10) << j << "     ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << ":";
//
//  Print out (up to) INCX entries in row I, that lie in the current strip.
//
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        cout << "  " << setw(8) << real ( c )
             << "  " << setw(8) << imag ( c );
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void c4mat_uniform_01 ( int m, int n, int *seed, complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
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
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C[M*N], the pseudorandom complex matrix.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return;
}
//****************************************************************************80

complex <float> *c4mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01_NEW returns a new unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
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
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <float> *c;
  int i;
  int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return c;
}
//****************************************************************************80

void c4vec_copy ( int n, complex <float> a1[], complex <float> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_COPY copies a C4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <float> A1[N], the vector to be copied.
//
//    Output, complex <float> A2[N], the copy of A1.
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

complex <float> *c4vec_copy_new ( int n, complex <float> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_COPY_NEW copies a C4VEC to a "new" C4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <float> A1[N], the vector to be copied.
//
//    Output, complex <float> C4VEC_COPY_NEW[N], the copy of A1.
//
{
  complex <float> *a2;
  int i;

  a2 = new complex <float>[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

complex <float> *c4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_INDICATOR_NEW sets a C4VEC to the indicator vector.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, complex <float> C4VEC_INDICATOR_NEW[N], the array.
//
{
  complex <float> *a;
  int i;

  a = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = complex <float> ( i+1, -i-1 );
  }

  return a;
}
//****************************************************************************80

float c4vec_norm_l2 ( int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_NORM_L2 returns the L2 norm of a C4VEC.
//
//  Discussion:
//
//    The vector L2 norm is defined as:
//
//      value = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, complex <float> A[N], the vector whose L2 norm is desired.
//
//    Output, float C4VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value 
          + real ( a[i] ) * real ( a[i] ) 
          + imag ( a[i] ) * imag ( a[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

void c4vec_print ( int n, complex <float> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT prints a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
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
//    Input, int N, the number of components of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
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
    cout << "  " << setw(8) << i
         << ": " << real ( a[i] )
         << "  " << imag ( a[i] ) << "\n";
  }

  return;
}
//****************************************************************************80

void c4vec_print_part ( int n, complex <float> a[], int max_print, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT_PART prints "part" of a C4VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines
//    to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
    cout << "  ........  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << "  " << setw(14) << real ( a[i] ) 
         << "  " << setw(14) << imag ( a[i] )
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

void c4vec_print_some ( int n, complex <float> a[], int i_lo, int i_hi, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT_SOME prints some of a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
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
//    Input, int N, the number of components of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
//
//    Input, int I_LO, I_HI, the first and last entries to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 0, i_lo ); i <= i4_min ( i_hi, n - 1 ); i++ )
  {
    cout << "  " << setw(6) << i
         << ": " << real ( a[i] )
         << "  " << imag ( a[i] ) << "\n";
  }

  return;
}
//****************************************************************************80

void c4vec_uniform_01 ( int n, int *seed, complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2006
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
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C[N], the pseudorandom 
//    complex vector.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return;
}
//****************************************************************************80

complex <float> *c4vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01_NEW returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2006
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
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
{
  complex <float> *c;
  int i;
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

complex <float> *c4vec_unity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNITY returns the N roots of unity in a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
//
//    X(1:N)^N = ( (1,0), (1,0), ..., (1,0) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, complex <float> C4VEC_UNITY[N], the N roots of unity.
//
{
  complex <float> *a;
  int i;
  float pi = 3.141592653589793;
  float theta;

  a = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( float ) ( 2 * i ) / ( float ) ( n );
    a[i] = complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return a;
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

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

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

float r4_max ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MAX returns the maximum of two R4's.
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
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R4_MAX, the maximum of X and Y.
//
{
  float value;

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

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
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
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }
  return value;
}
//****************************************************************************80

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a unit pseudorandom R4.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
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
