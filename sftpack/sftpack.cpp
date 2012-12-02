# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <complex>

using namespace std;

# include "sftpack.hpp"

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

complex <float> *c4mat_sftb ( int n1, int n2, complex <float> y[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_SFTB computes a "slow" backward Fourier transform of a C4MAT.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y,
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I1 <= N1 - 1, 
//        0 <= I2 <= N2 - 1,
//
//      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
//        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
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
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, complex <float> Y[N1*N2], the Fourier coefficients.
//
//    Output, complex <float> C4MAT_SFTB[N1*N2], the data.
//
{
  complex <float> cs1;
  complex <float> cs2;
  int i1;
  int i2;
  int j1;
  int j2;
  float pi = 3.141592653589793;
  float theta1;
  float theta2;
  complex <float> *x;

  x = new complex <float>[n1 * n2 ];

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      x[i1+i2*n1] = complex <float> ( 0.0, 0.0 );
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( j2 = 0; j2 < n2; j2++ )
    {
      theta2 = 2.0 * pi * ( float ) ( i2 * j2 ) / ( float ) ( n2 );
      cs2 = complex <float> ( cos ( theta2 ), - sin ( theta2 ) );
      for ( i1 = 0; i1 < n1; i1++ )
      {
        for ( j1 = 0; j1 < n1; j1++ )
        {
          theta1 = 2.0 * pi * ( float ) ( i1 * j1 ) / ( float ) ( n1 );
          cs1 = complex <float> ( cos ( theta1 ), - sin ( theta1 ) );
          x[i1+i2*n1] = x[i1+i2*n1] + y[j1+j2*n1] * cs1 * cs2;
        }
      }
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      x[i1+i2*n1] = x[i1+i2*n1] / ( float ) ( n1 * n2 );
    }
  }
  return x;
}
//****************************************************************************80

complex <float> *c4mat_sftf ( int n1, int n2, complex <float> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_SFTF computes a "slow" forward Fourier transform of a C4MAT.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y, 
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I1 <= N1 - 1, 
//        0 <= I2 <= N2 - 1,
//
//      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
//        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
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
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, complex <float> X[N1*N2], the data to be transformed.
//
//    Output, complex <float> C4MAT_SFTF[N1*N2], the Fourier coefficients.
//
{
  complex <float> cs1;
  complex <float> cs2;
  int i1;
  int i2;
  int j1;
  int j2;
  float pi = 3.141592653589793;
  float theta1;
  float theta2;
  complex <float> *y;

  y = new complex <float>[ n1 * n2 ];

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      y[i1+i2*n1] = complex <float> ( 0.0, 0.0 );
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( j2 = 0; j2 < n2; j2++ )
    {
      theta2 = - 2.0 * pi * ( float ) ( i2 * j2 ) / ( float ) ( n2 );
      cs2 = complex <float> ( cos ( theta2 ), - sin ( theta2 ) );
      for ( i1 = 0; i1 < n1; i1++ )
      {
        for ( j1 = 0; j1 < n1; j1++ )
        {
          theta1 = - 2.0 * pi * ( float ) ( i1 * j1 ) / ( float ) ( n1 );
          cs1 = complex <float> ( cos ( theta1 ), - sin ( theta1 ) );
          y[i1+i2*n1] = y[i1+i2*n1] + x[j1+j2*n1] * cs1 * cs2;
        }
      }
    }
  }
  return y;
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

complex <float> *c4vec_sftb ( int n, complex <float> y[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_SFTB computes a "slow" backward Fourier transform of a C4VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y,
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N - 1
//
//      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
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
//    Input, int N, the number of data values.
//
//    Input, complex <float> Y[N], the Fourier coefficients.
//
//    Output, complex <float> C4VEC_SFTB, the data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  complex <float> *x;

  x = new complex <float>[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( complex <float> ) ( 0.0, 0.0 );
    for ( j = 0; j < n; j++ )
    {
      theta = - 2.0 * pi * ( float ) ( i * j ) / ( float ) ( n );
      x[i] = x[i] + y[j] * complex <float> ( cos ( theta ), sin ( theta ) );
    }
    x[i] = x[i] / ( float ) ( n );
  }
  return x;
}
//****************************************************************************80

complex <float> *c4vec_sftf ( int n, complex <float> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_SFTF computes a "slow" forward Fourier transform of a C4VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y, 
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N - 1
//
//      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
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
//    Input, int N, the number of data values.
//
//    Input, complex <float> X[N], the data to be transformed.
//
//    Output, complex <float> C4VEC_SFTF[N], the Fourier coefficients.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  complex <float> *y;

  y = new complex <float>[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( complex <float> ) ( 0.0, 0.0 );
    for ( j = 0; j < n; j++ )
    {
      theta = - 2.0 * pi * ( float ) ( i * j ) / ( float ) ( n );
      y[i] = y[i] + x[j] * complex <float> ( cos ( theta ), - sin ( theta ) );
    }
  }
  return y;
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

complex <double> *c8mat_sftb ( int n1, int n2, complex <double> y[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_SFTB computes a "slow" backward Fourier transform of a C8MAT.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y,
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I1 <= N1 - 1, 
//        0 <= I2 <= N2 - 1,
//
//      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
//        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
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
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, complex <double> Y[N1*N2], the Fourier coefficients.
//
//    Output, complex <double> C8MAT_SFTB[N1*N2], the data.
//
{
  complex <double> cs1;
  complex <double> cs2;
  int i1;
  int i2;
  int j1;
  int j2;
  double pi = 3.141592653589793;
  double theta1;
  double theta2;
  complex <double> *x;

  x = new complex <double>[n1 * n2 ];

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      x[i1+i2*n1] = complex <double> ( 0.0, 0.0 );
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( j2 = 0; j2 < n2; j2++ )
    {
      theta2 = 2.0 * pi * ( double ) ( i2 * j2 ) / ( double ) ( n2 );
      cs2 = complex <double> ( cos ( theta2 ), - sin ( theta2 ) );
      for ( i1 = 0; i1 < n1; i1++ )
      {
        for ( j1 = 0; j1 < n1; j1++ )
        {
          theta1 = 2.0 * pi * ( double ) ( i1 * j1 ) / ( double ) ( n1 );
          cs1 = complex <double> ( cos ( theta1 ), - sin ( theta1 ) );
          x[i1+i2*n1] = x[i1+i2*n1] + y[j1+j2*n1] * cs1 * cs2;
        }
      }
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      x[i1+i2*n1] = x[i1+i2*n1] / ( double ) ( n1 * n2 );
    }
  }
  return x;
}
//****************************************************************************80

complex <double> *c8mat_sftf ( int n1, int n2, complex <double> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_SFTF computes a "slow" forward Fourier transform of a C8MAT.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y, 
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I1 <= N1 - 1, 
//        0 <= I2 <= N2 - 1,
//
//      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
//        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
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
//    Input, int N1, N2, the number of rows and columns of data.
//
//    Input, complex <double> X[N1*N2], the data to be transformed.
//
//    Output, complex <double> C8MAT_SFTF[N1*N2], the Fourier coefficients.
//
{
  complex <double> cs1;
  complex <double> cs2;
  int i1;
  int i2;
  int j1;
  int j2;
  double pi = 3.141592653589793;
  double theta1;
  double theta2;
  complex <double> *y;

  y = new complex <double>[ n1 * n2 ];

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( i1 = 0; i1 < n1; i1++ )
    {
      y[i1+i2*n1] = complex <double> ( 0.0, 0.0 );
    }
  }

  for ( i2 = 0; i2 < n2; i2++ )
  {
    for ( j2 = 0; j2 < n2; j2++ )
    {
      theta2 = - 2.0 * pi * ( double ) ( i2 * j2 ) / ( double ) ( n2 );
      cs2 = complex <double> ( cos ( theta2 ), - sin ( theta2 ) );
      for ( i1 = 0; i1 < n1; i1++ )
      {
        for ( j1 = 0; j1 < n1; j1++ )
        {
          theta1 = - 2.0 * pi * ( double ) ( i1 * j1 ) / ( double ) ( n1 );
          cs1 = complex <double> ( cos ( theta1 ), - sin ( theta1 ) );
          y[i1+i2*n1] = y[i1+i2*n1] + x[j1+j2*n1] * cs1 * cs2;
        }
      }
    }
  }
  return y;
}
//****************************************************************************80

complex <double> *c8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01_NEW returns a new unit pseudorandom C8MAT.
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
//    Output, complex <double> C8MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <double> *c;
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
    cerr << "C8MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <double>[m*n];

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

  return c;
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

complex <double> *c8vec_sftb ( int n, complex <double> y[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SFTB computes a "slow" backward Fourier transform of a C8VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y,
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N - 1
//
//      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, complex <double> Y[N], the Fourier coefficients.
//
//    Output, complex <double> C8VEC_SFTB, the data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  complex <double> *x;

  x = new complex <double>[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = complex <double> ( 0.0, 0.0 );
    for ( j = 0; j < n; j++ )
    {
      theta = - 2.0 * pi * ( double ) ( i * j ) / ( double ) ( n );
      x[i] = x[i] + y[j] * complex <double> ( cos ( theta ), sin ( theta ) );
    }
    x[i] = x[i] / ( double ) ( n );
  }
  return x;
}
//****************************************************************************80

complex <double> *c8vec_sftf ( int n, complex <double> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_SFTF computes a "slow" forward Fourier transform of a C8VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    X and apply SFTF to get Y, and then apply SFTB to Y, 
//    we should get back the original X.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N - 1
//
//      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, complex <double> X[N], the data to be transformed.
//
//    Output, complex <double> C8VEC_SFTF[N], the Fourier coefficients.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  complex <double> *y;

  y = new complex <double>[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = complex <double> ( 0.0, 0.0 );
    for ( j = 0; j < n; j++ )
    {
      theta = - 2.0 * pi * ( double ) ( i * j ) / ( double ) ( n );
      y[i] = y[i] + x[j] * complex <double> ( cos ( theta ), - sin ( theta ) );
    }
  }
  return y;
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
//    A C8VEC is a vector of C8's.
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
//    22 May 2006
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
  int i4_huge = 2147483647;
  int k;
  double pi = 3.141592653589793;
  double r;
  double theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <double> [n];

  for ( i = 0; i < n; i++ )
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

    c[i] = r * complex <double> ( cos ( theta ), sin ( theta ) );
  }

  return c;
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

void r4vec_print_part ( int n, float a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PRINT_PART prints "part" of an R4VEC.
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
//    09 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, float A[N], the vector to be printed.
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
           << "  " << setw(14) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    cout << "  ........  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i] << "\n";
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i]
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

float *r4vec_sct ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SCT computes a "slow" cosine transform of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//      Y(1) = Sum ( 1 <= J <= N ) X(J)
//
//      For 2 <= I <= N-1:
//
//        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
//          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
//
//      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
//
//    Applying the routine twice in succession should yield the original data,
//    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
//    and accuracy.
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
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SCT[N], the transformed data.
//
{
  float angle;
  int i;
  int j;
  float pi = 3.141592653589793;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[0] / 2.0;

    for ( j = 1; j < n - 1; j++ )
    {
      angle = pi * ( float ) ( ( i * j ) % ( 2 * ( n - 1 ) ) ) 
        / ( float ) ( n - 1 );
      y[i] = y[i] + x[j] * cos ( angle );
    }

    j = n - 1;

    angle = pi * ( float ) ( ( i * j ) % ( 2 * ( n - 1 ) ) ) 
      / ( float ) ( n - 1 );

    y[i] = y[i] + x[j] * cos ( angle ) / 2.0;
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = 2.0 * y[i] * sqrt ( ( float ) ( n ) / ( float ) ( n - 1 ) );
  }

  return y;
}
//****************************************************************************80

float *r4vec_sftb ( int n, float azero, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SFTB computes a "slow" backward Fourier transform of an R4VEC.
//
//  Discussion:
//
//    SFTB and SFTF are inverses of each other.  If we begin with data
//    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
//    resulting R vector, we should get back the original AZERO, A and B.
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
//    Input, int N, the number of data values.
//
//    Input, float AZERO, the constant Fourier coefficient.
//
//    Input, float A[N/2], B[N/2], the Fourier coefficients.
//
//    Output, float SFTB[N], the reconstructed data sequence.
//
{
  int i;
  int k;
  float pi = 3.141592653589793;
  float *r;
  float theta;

  r = new float[n];

  for ( i = 0; i < n; i++ )
  {
    r[i] = azero;
    for ( k = 0; k < ( n / 2 ); k++ )
    {
      theta = ( float ) ( ( k + 1 ) * i * 2 ) * pi / ( float ) ( n );
      r[i] = r[i] + a[k] * cos ( theta ) + b[k] * sin ( theta );
    }
  }

  return r;
}
//****************************************************************************80

void r4vec_sftf ( int n, float r[], float *azero, float a[], float b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SFTF computes a "slow" forward Fourier transform of an R4VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
//    A, and B, we should get back the original R.
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
//    Input, int N, the number of data values.
//
//    Input, float R[N], the data to be transformed.
//
//    Output, float *AZERO, = sum ( 1 <= I <= N ) R(I) / N.
//
//    Output, float A[N/2], B[N/2], the Fourier coefficients.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;

  *azero = 0.0;
  for ( i = 0; i < n; i++ )
  {
    *azero = *azero + r[i];
  }
  *azero = *azero / ( float ) ( n );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    a[i] = 0.0;
    b[i] = 0.0;

    for ( j = 0; j < n; j++ )
    {
      theta = ( float ) ( 2 * ( i + 1 ) * j ) * pi / ( float ) ( n );
      a[i] = a[i] + r[j] * cos ( theta );
      b[i] = b[i] + r[j] * sin ( theta );
    }

    a[i] = a[i] / ( float ) ( n );
    b[i] = b[i] / ( float ) ( n );

    if ( i != ( n / 2 - 1 ) )
    {
      a[i] = 2.0 * a[i];
      b[i] = 2.0 * b[i];
    }
  }
  return;
}
//****************************************************************************80

float *r4vec_sht ( int n, float a[]  )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SHT computes a "slow" Hartley transform of an R4VEC.
//
//  Discussion:
//
//    The discrete Hartley transform B of a set of data A is
//
//      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
//
//    Here, the data and coefficients are indexed from 0 to N-1.
//
//    With the above normalization factor of 1/sqrt(N), the Hartley
//    transform is its own inverse.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines.
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
//  Reference:
//
//    Ralph Hartley,
//    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
//    Proceedings of the Institute of Radio Engineers,
//    Volume 30, pages 144-150, 1942.
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float A[N], the data to be transformed.
//
//    Output, float SHT[N], the transformed data.
//
{
  float *b;
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;

  b = new float[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 2.0 * pi * ( float ) ( ( i * j ) % n ) / ( float ) ( n );
      b[i] = b[i] + a[j] * ( cos ( theta ) + sin ( theta ) );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    b[i] = b[i] / sqrt ( ( float ) ( n  ) );
  }

  return b;
}
//****************************************************************************80

float *r4vec_sqctb ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SQCTB computes a "slow" quarter cosine transform backward of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
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
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    LC: QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SQCTB[N], the transformed data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[0];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 1; j < n; j++ )
    {
      theta = 0.5 * pi * ( float ) ( j * ( 2 * i + 1 ) ) / ( float ) ( n );
      y[i] = y[i] + 2.0 * x[j] * cos ( theta  );
    }
  }
  return y;
}
//****************************************************************************80

float *r4vec_sqctf ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SQCTF computes a "slow" quarter cosine transform forward of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
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
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SQCTF[N], the transformed data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 0.5 * pi * ( float ) ( i * ( 2 * j + 1 ) ) / ( float ) ( n );
      y[i] = y[i] + x[j] * cos ( theta  );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = y[i] / ( float ) ( n );
  }

  return y;
}
//****************************************************************************80

float *r4vec_sqstb ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SQSTB computes a "slow" quarter sine transform backward of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
//             - X(N) * cos ( pi * I )
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
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SQSTB[N], the transformed data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n - 1; j++ )
    {
      theta = 0.5 * pi * ( float ) ( ( j + 1 ) * ( 2 * i + 1 ) ) 
        / ( float ) ( n );
      y[i] = y[i] - 2.0 * x[j] * sin ( theta  );
    }
    theta = pi * ( float ) ( i );
    y[i] = y[i] - x[n-1] * cos ( theta );

  }
  return y;
}
//****************************************************************************80

float *r4vec_sqstf ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SQSTF computes a "slow" quarter sine transform forward of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 1 <= I <= N,
//
//      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
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
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SQSTF{N], the transformed data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float theta;
  float *y;

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 0.5 * pi * ( float ) (  ( i + 1 ) * ( 2 * j + 1 ) ) 
        / ( float ) ( n );
      y[i] = y[i] + x[j] * sin ( theta );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = - y[i] / ( float ) ( n );
  }

  return y;
}
//****************************************************************************80

float *r4vec_sst ( int n, float x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_SST computes a "slow" sine transform of an R4VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 1 <= I <= N,
//
//      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
//
//    Applying the routine twice in succession should yield the original data,
//    multiplied by N / 2.  This is a good check for correctness and accuracy.
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
//    Input, int N, the number of data values.
//
//    Input, float X[N], the data sequence.
//
//    Output, float SST[N], the transformed data.
//
{
  int i;
  int j;
  float pi = 3.141592653589793;
  float *theta;
  float *y;

  theta = new float[n];

  for ( i = 0; i < n; i++ )
  {
    theta[i] = pi * ( float ) ( i + 1 ) / ( float ) ( n + 1 );
  }

  y = new float[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] + 2.0 * x[j] * sin ( ( float ) ( j + 1 ) * theta[i] );
    }
  }

  delete [] theta;

  return y;
}
//****************************************************************************80

float *r4vec_uniform_new ( int n, float b, float c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_NEW returns a scaled pseudorandom R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
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
//    Input, int N, the number of entries in the vector.
//
//    Input, float B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float R4VEC_UNIFORM_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new float[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec_print_part ( int n, double a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_PART prints "part" of an R8VEC.
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
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
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
           << "  " << setw(14) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i] << "\n";
    }
    cout << "  ........  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << "  " << setw(14) << a[i] << "\n";
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i] << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << "  " << setw(14) << a[i]
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_sct ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SCT computes a "slow" cosine transform of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//      Y(1) = Sum ( 1 <= J <= N ) X(J)
//
//      For 2 <= I <= N-1:
//
//        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
//          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
//
//      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
//
//    Applying the routine twice in succession should yield the original data,
//    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
//    and accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SCT[N], the transformed data.
//
{
  double angle;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[0] / 2.0;

    for ( j = 1; j < n - 1; j++ )
    {
      angle = pi * ( double ) ( ( i * j ) % ( 2 * ( n - 1 ) ) ) 
        / ( double ) ( n - 1 );
      y[i] = y[i] + x[j] * cos ( angle );
    }

    j = n - 1;

    angle = pi * ( double ) ( ( i * j ) % ( 2 * ( n - 1 ) ) ) 
      / ( double ) ( n - 1 );

    y[i] = y[i] + x[j] * cos ( angle ) / 2.0;
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = 2.0 * y[i] * sqrt ( ( double ) ( n ) / ( double ) ( n - 1 ) );
  }

  return y;
}
//****************************************************************************80

double *r8vec_sftb ( int n, double azero, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SFTB computes a "slow" backward Fourier transform of an R8VEC.
//
//  Discussion:
//
//    SFTB and SFTF are inverses of each other.  If we begin with data
//    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
//    resulting R vector, we should get back the original AZERO, A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double AZERO, the constant Fourier coefficient.
//
//    Input, double A[N/2], B[N/2], the Fourier coefficients.
//
//    Output, double SFTB[N], the reconstructed data sequence.
//
{
  int i;
  int k;
  double pi = 3.141592653589793;
  double *r;
  double theta;

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    r[i] = azero;
    for ( k = 0; k < ( n / 2 ); k++ )
    {
      theta = ( double ) ( ( k + 1 ) * i * 2 ) * pi / ( double ) ( n );
      r[i] = r[i] + a[k] * cos ( theta ) + b[k] * sin ( theta );
    }
  }

  return r;
}
//****************************************************************************80

void r8vec_sftf ( int n, double r[], double *azero, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SFTF computes a "slow" forward Fourier transform of an R8VEC.
//
//  Discussion:
//
//    SFTF and SFTB are inverses of each other.  If we begin with data
//    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
//    A, and B, we should get back the original R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double R[N], the data to be transformed.
//
//    Output, double *AZERO, = sum ( 1 <= I <= N ) R(I) / N.
//
//    Output, double A[N/2], B[N/2], the Fourier coefficients.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  *azero = 0.0;
  for ( i = 0; i < n; i++ )
  {
    *azero = *azero + r[i];
  }
  *azero = *azero / ( double ) ( n );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    a[i] = 0.0;
    b[i] = 0.0;

    for ( j = 0; j < n; j++ )
    {
      theta = ( double ) ( 2 * ( i + 1 ) * j ) * pi / ( double ) ( n );
      a[i] = a[i] + r[j] * cos ( theta );
      b[i] = b[i] + r[j] * sin ( theta );
    }

    a[i] = a[i] / ( double ) ( n );
    b[i] = b[i] / ( double ) ( n );

    if ( i != ( n / 2 - 1 ) )
    {
      a[i] = 2.0 * a[i];
      b[i] = 2.0 * b[i];
    }
  }
  return;
}
//****************************************************************************80

double *r8vec_sht ( int n, double a[]  )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SHT computes a "slow" Hartley transform of an R8VEC.
//
//  Discussion:
//
//    The discrete Hartley transform B of a set of data A is
//
//      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
//
//    Here, the data and coefficients are indexed from 0 to N-1.
//
//    With the above normalization factor of 1/sqrt(N), the Hartley
//    transform is its own inverse.
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ralph Hartley,
//    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
//    Proceedings of the Institute of Radio Engineers,
//    Volume 30, pages 144-150, 1942.
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double A[N], the data to be transformed.
//
//    Output, double SHT[N], the transformed data.
//
{
  double *b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 2.0 * pi * ( double ) ( ( i * j ) % n ) / ( double ) ( n );
      b[i] = b[i] + a[j] * ( cos ( theta ) + sin ( theta ) );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    b[i] = b[i] / sqrt ( ( double ) ( n  ) );
  }

  return b;
}
//****************************************************************************80

double *r8vec_sqctb ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SQCTB computes a "slow" quarter cosine transform backward of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    LC: QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SQCTB[N], the transformed data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[0];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 1; j < n; j++ )
    {
      theta = 0.5 * pi * ( double ) ( j * ( 2 * i + 1 ) ) / ( double ) ( n );
      y[i] = y[i] + 2.0 * x[j] * cos ( theta  );
    }
  }
  return y;
}
//****************************************************************************80

double *r8vec_sqctf ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SQCTF computes a "slow" quarter cosine transform forward of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SQCTF[N], the transformed data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 0.5 * pi * ( double ) ( i * ( 2 * j + 1 ) ) / ( double ) ( n );
      y[i] = y[i] + x[j] * cos ( theta  );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = y[i] / ( double ) ( n );
  }

  return y;
}
//****************************************************************************80

double *r8vec_sqstb ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SQSTB computes a "slow" quarter sine transform backward of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 0 <= I <= N-1,
//
//      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
//             - X(N) * cos ( pi * I )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SQSTB[N], the transformed data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n - 1; j++ )
    {
      theta = 0.5 * pi * ( double ) ( ( j + 1 ) * ( 2 * i + 1 ) ) 
        / ( double ) ( n );
      y[i] = y[i] - 2.0 * x[j] * sin ( theta  );
    }
    theta = pi * ( double ) ( i );
    y[i] = y[i] - x[n-1] * cos ( theta );

  }
  return y;
}
//****************************************************************************80

double *r8vec_sqstf ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SQSTF computes a "slow" quarter sine transform forward of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 1 <= I <= N,
//
//      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Briggs, Van Emden Henson,
//    The Discrete Fourier Transform,
//    SIAM, 1995,
//    QA403.5 B75
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SQSTF{N], the transformed data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta = 0.5 * pi * ( double ) (  ( i + 1 ) * ( 2 * j + 1 ) ) 
        / ( double ) ( n );
      y[i] = y[i] + x[j] * sin ( theta );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = - y[i] / ( double ) ( n );
  }

  return y;
}
//****************************************************************************80

double *r8vec_sst ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SST computes a "slow" sine transform of an R8VEC.
//
//  Discussion:
//
//    This routine is provided for illustration and testing.  It is inefficient
//    relative to optimized routines that use fast Fourier techniques.
//
//    For 1 <= I <= N,
//
//      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
//
//    Applying the routine twice in succession should yield the original data,
//    multiplied by N / 2.  This is a good check for correctness and accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double X[N], the data sequence.
//
//    Output, double SST[N], the transformed data.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double *theta;
  double *y;

  theta = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta[i] = pi * ( double ) ( i + 1 ) / ( double ) ( n + 1 );
  }

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] + 2.0 * x[j] * sin ( ( double ) ( j + 1 ) * theta[i] );
    }
  }

  delete [] theta;

  return y;
}
//****************************************************************************80

double *r8vec_uniform_new ( int n, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_NEW returns a scaled pseudorandom R8VEC.
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
//    30 January 2005
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
//    Input, double B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_NEW - Fatal error!\n";
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

    r[i] = b + ( c - b ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
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
