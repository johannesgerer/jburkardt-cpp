# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <cstring>

using namespace std;

# include "c8lib.hpp"

//****************************************************************************80

double c8_abs ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ABS returns the absolute value of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_ABS, the magnitude of X.
//
{
  double value;

  value = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

double c8_argument ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ARGUMENT returns the argument of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the value whose argument is desired.
//
//    Output, double C8_ARGUMENT, the argument of X.
//
{
  double value;

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

complex <double> c8_cube_root ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_CUBE_ROOT returns the principal cube root of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Input, complex <double> X, the number whose cube root is desired.
//
//    Output, complex <double> C8_CUBE_ROOT, the cube root of X.
//
{
  double argument;
  double magnitude;
  complex <double> value;

  argument = c8_argument ( x );
  magnitude = c8_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <double> ( 0.0, 0.0 );
  }
  else
  {
    value = pow ( magnitude, ( double ) ( 1.0 / 3.0 ) ) 
      * complex <double> ( cos ( argument / 3.0 ), sin ( argument / 3.0 ) );
  }

  return value;
}
//****************************************************************************80

complex <double> c8_i ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_I returns the value of the imaginary unit, i as a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Output, complex <double> C8_I, the value of complex i.
//
{
  complex <double> value;

  value = complex <double> ( 0.0, 1.0 );

  return value;
}
//****************************************************************************80

bool c8_le_l1 ( complex <double> x, complex <double> y )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_L1 := X <= Y for C8 values, and the L1 norm.
//
//  Discussion:
//
//    A C8 is a complex double precision value.
//
//    The L1 norm can be defined here as:
//
//      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
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
//    Input, complex <double> X, Y, the values to be compared.
//
//    Output, bool C8_LE_L1, is TRUE if X <= Y.
//
{
  bool value;

  if ( r8_abs ( real ( x ) ) + r8_abs ( imag ( x ) ) <= 
       r8_abs ( real ( y ) ) + r8_abs ( imag ( y ) ) )
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

bool c8_le_l2 ( complex <double> x, complex <double> y )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_L2 := X <= Y for C8 values, and the L2 norm.
//
//  Discussion:
//
//    A C8 is a complex double precision value.
//
//    The L2 norm can be defined here as:
//
//      C8_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
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
//    Input, complex <double> X, Y, the values to be compared.
//
//    Output, bool C8_LE_L2, is TRUE if X <= Y.
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

bool c8_le_li ( complex <double> x, complex <double> y )

//****************************************************************************80
//
//  Purpose:
//
//    C8_LE_LI := X <= Y for C8 values, and the L-oo norm.
//
//  Discussion:
//
//    A C8 is a complex double precision value.
//
//    The L-oo norm can be defined here as:
//
//      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
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
//    Input, complex <double> X, Y, the values to be compared.
//
//    Output, bool C8_LE_LI, is TRUE if X <= Y.
//
{
  bool value;

  if ( r8_max ( r8_abs ( real ( x ) ), r8_abs ( imag ( x ) ) ) <= 
       r8_max ( r8_abs ( real ( y ) ), r8_abs ( imag ( y ) ) ) )
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

double c8_magnitude ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_MAGNITUDE returns the magnitude of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_MAGNITUDE, the magnitude of X.
//
{
  double magnitude;

  magnitude = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return magnitude;
}
//****************************************************************************80

double c8_norm_l1 ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_L1 evaluates the L1 norm of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//    Numbers of equal norm lie along diamonds centered at (0,0).
//
//    The L1 norm can be defined here as:
//
//      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
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
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_NORM_L1, the norm of X.
//
{
  double value;

  value = r8_abs ( real ( x ) ) + r8_abs ( imag ( x ) );

  return value;
}
//****************************************************************************80

double c8_norm_l2 ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_L2 evaluates the L2 norm of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//    Numbers of equal norm lie on circles centered at (0,0).
//
//    The L2 norm can be defined here as:
//
//      C8_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_NORM_L2, the 2-norm of X.
//
{
  double value;

  value = sqrt ( pow ( real ( x ), 2 )
               + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

double c8_norm_li ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORM_LI evaluates the L-oo norm of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//    Numbers of equal norm lie along squares whose centers are at (0,0).
//
//    The L-oo norm can be defined here as:
//
//      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
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
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_NORM_LI, the L-oo norm of X.
//
{
  double value;

  value = r8_max ( r8_abs ( real ( x ) ), r8_abs ( imag ( x ) ) );

  return value;
}
//****************************************************************************80

complex <double> c8_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORMAL_01 returns a unit pseudonormal C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Output, complex <double> C8_NORMAL_01, a unit pseudornormal value.
//
{
  double pi = 3.141592653589793;
  double v1;
  double v2;
  double x_c;
  double x_r;
  complex <double> value;

  v1 = r8_uniform_01 ( seed );
  v2 = r8_uniform_01 ( seed );

  x_r = sqrt ( -2.0 * log ( v1 ) ) * cos ( 2.0 * pi * v2 );
  x_c = sqrt ( -2.0 * log ( v1 ) ) * sin ( 2.0 * pi * v2 );

  value = complex <double> ( x_r, x_c );

  return value;
}
//****************************************************************************80

complex <double> c8_one ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ONE returns the value of complex 1.
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
//    Output, complex <double> C8_ONE, the value of complex 1.
//
{
  complex <double> value;

  value = complex <double> ( 1.0, 0.0);

  return value;
}
//****************************************************************************80

void c8_print ( complex <double> a, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8_PRINT prints a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Input, complex <double> A, the value to be printed.
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

complex <double> c8_sqrt ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SQRT returns the principal square root of a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Input, complex <double> X, the number whose square root is desired.
//
//    Output, complex <double> C8_SQRT, the square root of X.
//
{
  double argument;
  double magnitude;
  complex <double> value;

  argument = c8_argument ( x );
  magnitude = c8_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <double> ( 0.0, 0.0 );
  }
  else
  {
    value = sqrt ( magnitude ) 
      * complex <double> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
  }

  return value;
}
//****************************************************************************80

void c8_swap ( complex <double> *x, complex <double> *y )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SWAP swaps two C8's.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, complex <double> *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  complex <double> z;

   z = *x; 
  *x = *y;
  *y =  z;

  return;
}
//****************************************************************************80

complex <double> c8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01 returns a unit pseudorandom C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8_UNIFORM_01, a pseudorandom complex value.
//
{
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;
  complex <double> value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = r * complex <double> ( cos ( theta ), sin ( theta ) );

  return value;
}
//****************************************************************************80

complex <double> c8_zero ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ZERO returns the value of 0 as a C8.
//
//  Discussion:
//
//    A C8 is a complex <double> value.
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
//    Output, complex <double> C8_ZERO, the value of complex 0.
//
{
  complex <double> value;

  value = complex <double> ( 0.0, 0.0 );

  return value;
}
//****************************************************************************80

void c8mat_copy ( int m, int n, complex <double> a1[], complex <double> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_COPY copies one C8MAT to another.
//
//  Discussion:
//
//    An C8MAT is a doubly dimensioned array of complex double precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <double> A1[M*N], the matrix to be copied.
//
//    Output, complex <double> A2[M*N], the copy of A1.
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

complex <double> *c8mat_copy_new ( int m, int n, complex <double> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_COPY_NEW copies one C8MAT to a "new" C8MAT.
//
//  Discussion:
//
//    An C8MAT is a doubly dimensioned array of complex double precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <double> A1[M*N], the matrix to be copied.
//
//    Output, complex <double> C8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  complex <double> *a2;
  int i;
  int j;

  a2 = new complex <double>[m*n];

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

complex <double> *c8mat_identity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_IDENTITY sets a C8MAT to the identity.
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
//    Output, complex <double> C8MAT_IDENTITY[N*N], the matrix.
//
{
  complex <double> *a;
  int i;
  int j;

  a = new complex <double> [n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = complex <double> ( 1.0, 0.0 );
      }
      else
      {
        a[i+j*n] = complex <double> ( 0.0, 0.0 );
      }
    }
  }
  return a;
}
//****************************************************************************80

complex <double> *c8mat_indicator_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_INDICATOR_NEW returns the C8MAT indicator matrix.
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
//    Output, complex <double> C8MAT_INDICATOR_NEW[M*N], the matrix.
//
{
  complex <double> *a;
  int i;
  int j;

  a = new complex <double> [m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <double> ( i, j );
    }
  }
  return a;
}
//****************************************************************************80

complex <double> *c8mat_mm_new ( int n1, int n2, int n3, complex <double> a[], 
  complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, complex <double> A[N1*N2], complex <double> B[N2*N3], 
//    the matrices to multiply.
//
//    Output, complex <double> C[N1*N3], the product matrix C = A * B.
//
{
  complex <double> *c;
  int i;
  int j;
  int k;

  c = new complex <double> [n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }
  return c;
}
//****************************************************************************80

void c8mat_nint ( int m, int n, complex <double> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_NINT rounds the entries of a C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
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
//    Input/output, complex <double> A[M*N], the matrix to be NINT'ed.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <double> ( 
        r8_nint ( real ( a[i+j*m] ) ), 
        r8_nint ( imag ( a[i+j*m] ) ) );
    }
  }
  return;
}
//****************************************************************************80

void c8mat_print ( int m, int n, complex <double> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_PRINT prints a C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <double> A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  c8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void c8mat_print_some ( int m, int n, complex <double> a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_PRINT_SOME prints some of a C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <double> A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
  complex <double> c;
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

void c8mat_uniform_01 ( int m, int n, int *seed, complex <double> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
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
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C[M*N], the pseudorandom complex matrix.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C8MAT_UNIFORM_01 - Fatal error!\n";
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

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <double> ( cos ( theta ), sin ( theta ) );
    }
  }

  return;
}
//****************************************************************************80

complex <double> *c8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01_NEW returns a unit pseudorandom C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
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
//    15 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
//    complex matrix.
//
{
  complex <double> *c;
  int i;
  int j;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  c = new complex <double> [m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <double> ( cos ( theta ), sin ( theta ) );
    }
  }
  return c;
}
//****************************************************************************80

void c8vec_copy ( int n, complex <double> a1[], complex <double> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_COPY copies a C8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <double> A1[N], the vector to be copied.
//
//    Output, complex <double> A2[N], the copy of A1.
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

complex <double> *c8vec_copy_new ( int n, complex <double> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_COPY_NEW copies a C8VEC to a "new" C8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <double> A1[N], the vector to be copied.
//
//    Output, complex <double> C8VEC_COPY_NEW[N], the copy of A1.
//
{
  complex <double> *a2;
  int i;

  a2 = new complex <double>[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

complex <double> *c8vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_INDICATOR_NEW sets a C8VEC to the indicator vector.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
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
//    Output, complex <double> C8VEC_INDICATOR_NEW[N], the array.
//
{
  complex <double> *a;
  int i;

  a = new complex <double> [n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = complex <double> ( i+1, -i-1 );
  }

  return a;
}
//****************************************************************************80

double c8vec_norm_l2 ( int n, complex <double> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_NORM_L2 returns the L2 norm of a C8VEC.
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
//    08 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, complex <double> A[N], the vector whose L2 norm is desired.
//
//    Output, double C8VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  double value;

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

void c8vec_print ( int n, complex <double> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_PRINT prints a C8VEC.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
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
//    Input, complex <double> A[N], the vector to be printed.
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

void c8vec_print_part ( int n, complex <double> a[], int max_print, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_PRINT_PART prints "part" of a C8VEC.
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
//    Input, complex <double> A[N], the vector to be printed.
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

void c8vec_print_some ( int n, complex <double> a[], int i_lo, int i_hi, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_PRINT_SOME prints some of a C8VEC.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, complex <double> A[N], the vector to be printed.
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

complex <double> *c8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNIFORM_01_NEW returns a unit pseudorandom C8VEC.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
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
//    15 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
{
  complex <double> *c;
  int i;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  c = new complex <double> [n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <double> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

complex <double> *c8vec_unity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNITY returns the N roots of unity in a C8VEC.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
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
//    Output, complex <double> C8VEC_UNITY[N], the N roots of unity.
//
{
  complex <double> *a;
  int i;
  double pi = 3.141592653589793;
  double theta;

  a = new complex <double> [n];

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( double ) ( 2 * i ) / ( double ) ( n );
    a[i] = complex <double> ( cos ( theta ), sin ( theta ) );
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

complex <double> r8_csqrt ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSQRT returns the complex square root of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose square root is desired.
//
//    Output, complex <double> R8_CSQRT, the square root of X:
//
{
  double argument;
  double magnitude;
  double pi = 3.141592653589793;
  complex <double> value;

  if ( 0.0 < x )
  {
    magnitude = x;
    argument = 0.0;
  }
  else if ( 0.0 == x )
  {
    magnitude = 0.0;
    argument = 0.0;
  }
  else if ( x < 0.0 )
  {
    magnitude = -x;
    argument = pi;
  }

  magnitude = sqrt ( magnitude );
  argument = argument / 2.0;

  value = magnitude * complex <double> ( cos ( argument ), sin ( argument ) );

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

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
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
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r8_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r8_abs ( x ) + 0.5 );
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

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
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
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
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
//    11 August 2004
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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8poly2_root ( double a, double b, double c, complex <double> *r1,
  complex <double> *r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY2_ROOT returns the two roots of a quadratic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X^2 + B * X + C = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Parameters:
//
//    Input, double A, B, C, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <double> *R1, *R2, the roots of the polynomial, which
//    might be real and distinct, real and equal, or complex conjugates.
//
{
  double disc;
  complex <double> q;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R8POLY2_ROOT - Fatal error!\n";
    cerr << "  The coefficient A is zero.\n";
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  q = -0.5 * ( b + r8_sign ( b ) * r8_csqrt ( disc ) );
  *r1 = q / a;
  *r2 = c / q;

  return;
}
//****************************************************************************80

void r8poly3_root ( double a, double b, double c, double d,
  complex <double> *r1, complex <double> *r2, complex <double> *r3 )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY3_ROOT returns the three roots of a cubic polynomial.
//
//  Discussion:
//
//    The polynomial has the form
//
//      A * X^3 + B * X^2 + C * X + D = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Parameters:
//
//    Input, double A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <double> *R1, *R2, *R3, the roots of the polynomial, which
//    will include at least one real root.
//
{
  complex <double> i;
  double pi = 3.141592653589793;
  double q;
  double r;
  double s1;
  double s2;
  double temp;
  double theta;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R8POLY3_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  i = complex <double> ( 0.0, 1.0 );

  q = ( pow ( b / a, 2 ) - 3.0 * ( c / a ) ) / 9.0;

  r = ( 2.0 * pow ( b / a, 3 ) - 9.0 * ( b / a ) * ( c / a )
      + 27.0 * ( d / a ) ) / 54.0;

  if ( r * r < q * q * q )
  {
    theta = acos ( r / sqrt ( pow ( q, 3 ) ) );
    *r1 = -2.0 * sqrt ( q ) * cos (   theta              / 3.0 );
    *r2 = -2.0 * sqrt ( q ) * cos ( ( theta + 2.0 * pi ) / 3.0 );
    *r3 = -2.0 * sqrt ( q ) * cos ( ( theta + 4.0 * pi ) / 3.0 );
  }
  else if ( q * q * q <= r * r )
  {
    temp = -r + sqrt ( r * r - q * q * q );
    s1 = r8_sign ( temp ) * pow ( r8_abs ( temp ), 1.0 / 3.0 );

    temp = -r - sqrt ( r * r - q * q * q );
    s2 = r8_sign ( temp ) * pow ( r8_abs ( temp ), 1.0 / 3.0 );

    *r1 = s1 + s2;
    *r2 = -0.5 * ( s1 + s2 ) + i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
    *r3 = -0.5 * ( s1 + s2 ) - i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
  }

  *r1 = *r1 - b / ( 3.0 * a );
  *r2 = *r2 - b / ( 3.0 * a );
  *r3 = *r3 - b / ( 3.0 * a );

  return;
}
//****************************************************************************80

void r8poly4_root ( double a, double b, double c, double d, double e,
  complex <double> *r1, complex <double> *r2, complex <double> *r3,
  complex <double> *r4 )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY4_ROOT returns the four roots of a quartic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2005
//
//  Parameters:
//
//    Input, double A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <double> *R1, *R2, *R3, *R4, the roots of the polynomial.
//
{
  double a3;
  double a4;
  double b3;
  double b4;
  double c3;
  double c4;
  double d3;
  double d4;
  complex <double> p;
  complex <double> q;
  complex <double> r;
  complex <double> zero;

  zero = 0.0;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R8POLY4_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  a4 = b / a;
  b4 = c / a;
  c4 = d / a;
  d4 = e / a;
//
//  Set the coefficients of the resolvent cubic equation.
//
  a3 = 1.0;
  b3 = -b4;
  c3 = a4 * c4 - 4.0 * d4;
  d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
//
//  Find the roots of the resolvent cubic.
//
  r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
//
//  Choose one root of the cubic, here R1.
//
//  Set R = sqrt ( 0.25 * A4^2 - B4 + R1 )
//
  r = c8_sqrt ( 0.25 * a4 * a4 - b4  + (*r1) );

  if ( real ( r ) != 0.0 || imag ( r ) != 0.0 )
  {
    p = c8_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        + 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );

    q = c8_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        - 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );
  }
  else
  {
    p = c8_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      + 2.0 * c8_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );

    q = c8_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      - 2.0 * c8_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );
  }
//
//  Set the roots.
//
  *r1 = -0.25 * a4 + 0.5 * r + 0.5 * p;
  *r2 = -0.25 * a4 + 0.5 * r - 0.5 * p;
  *r3 = -0.25 * a4 - 0.5 * r + 0.5 * q;
  *r4 = -0.25 * a4 - 0.5 * r - 0.5 * q;

  return;
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
