# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "mgs.hpp"

int main ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
int r4_nint ( float x );
void r4mat_delete ( float **a, int m, int n );
float **r4mat_new ( int m, int n );
void r4mat_uniform ( int m, int n, float b, float c, int *seed, float **r );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MGS_PRB.
//
//  Discussion:
//
//    MGS_PRB gives some test data to the MGS function.
//
//  Modified:
//
//    07 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  float **a;
  float **a_qr;
  float **a_save;
  float a_hi = 10.0;
  float a_lo = -10.0;
  float diff;
  int i;
  int j;
  int k;
  int m;
  int n;
  float **q;
  float **r;
  int seed = 123456789;
  int test;

  cout << "\n";
  cout << "MGS_PRB:\n";
  cout << "  Test cases for MGS.\n";
  cout << "\n";

  for ( test = 1; test <= 4; test++ )
  {
    m = i4_uniform ( 1, 20, &seed );
    n = i4_uniform ( 1, 20, &seed );

    a = r4mat_new ( m, n );
    a_qr = r4mat_new ( m, n );
    a_save = r4mat_new ( m, n );
    q = r4mat_new ( m, n );
    r = r4mat_new ( n, n );

    r4mat_uniform ( m, n, a_lo, a_hi, &seed, a );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a_save[i][j] = a[i][j];
      }
    }

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        q[i][j] = 0.0;
      }
    }

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i][j] = 0.0;
      }
    }

    mgs ( m, n, a, r, q );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a_qr[i][j] = 0.0;
        for ( k = 0; k < n; k++ )
        {
          a_qr[i][j] = a_qr[i][j] + q[i][k] * r[k][j];
        }
      }
    }
    diff = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        diff = diff + pow ( a_save[i][j] - a_qr[i][j], 2 );
      }
    }
    diff = diff / sqrt ( ( float ) ( m * n ) );

    cout << "  " << setw(2) << m
         << "  " << setw(2) << n
         << "  " << setw(14) << diff << "\n";

    r4mat_delete ( a, m, n );
    r4mat_delete ( a_qr, m, n );
    r4mat_delete ( a_save, m, n );
    r4mat_delete ( q, m, n );
    r4mat_delete ( r, n, n );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "MGS_PRB:\n";
  cout << "  Normal end of execution.\n";

  return 0;
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

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
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
//    in Handbook of Simulation,
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
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

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
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
}
//****************************************************************************80

void r4mat_delete ( float **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_DELETE frees memory associated with an R4MAT.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      float **a;
//      a = r4mat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float **A, the pointer to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

float **r4mat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_NEW allocates a new R4MAT.
//
//  Discussion:
//
//    A declaration of the form
//      float **a;
//    is necesary.  Then an assignment of the form:
//      a = r4mat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, float **R4MAT_NEW, a new M by N matrix.
//
{
  float **a;
  int i;

  a = new float *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R4MAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new float[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R4MAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}
//****************************************************************************80

void r4mat_uniform ( int m, int n, float b, float c, int *seed, float **r )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM returns a scaled pseudorandom R4MAT.
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
//    23 April 2008
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
//    Input, int M, N, the number of rows and columns.
//
//   Input, float B, C, the limits of the pseudorandom values.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, float **R, a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4MAT_UNIFORM - Fatal error!\n";
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

      r[i][j] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
