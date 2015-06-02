# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa053.hpp"

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

double r8_uniform_01 ( int &seed )

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
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
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

void r8pp_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT prints a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, string TITLE, a title.
//
{
  r8pp_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8pp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8PP_PRINT_SOME prints some of a R8PP matrix.
//
//  Discussion:
//
//    The R8PP storage format is appropriate for a symmetric positive
//    definite matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the R8PP matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          aij = a[i-1+(j*(j-1))/2];
        }
        else
        {
          aij = a[j-1+(i*(i-1))/2];
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8utp_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT prints a R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
//
//    Input, string TITLE, a title.
//
{
  r8utp_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r8utp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8UTP_PRINT_SOME prints some of an R8UTP matrix.
//
//  Discussion:
//
//    The R8UTP storage format is appropriate for an upper triangular
//    matrix.  Only the upper triangle of the matrix is stored,
//    by successive partial columns, in an array of length (N*(N+1))/2,
//    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[(N*(N+1))/2], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  double aij;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, n );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i <= j )
        {
          aij = a[i-1+(j*(j-1))/2];
        }
        else
        {
          aij = 0.0;
        }

        cout << setw(12) << aij << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void rnorm ( int &seed, double &u1, double &u2 )

//****************************************************************************80
//
//  Purpose:
//
//    RNORM returns two independent standard random normal deviates.
//
//  Discussion:
//
//    This routine sets U1 and U2 to two independent standardized 
//    random normal deviates.   This is a version of the 
//    method given in Knuth.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    Original FORTRAN77 version by William Smith, Ronald Hocking.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 2, Seminumerical Algorithms,
//    Third Edition,
//    Addison Wesley, 1997,
//    ISBN: 0201896842,
//    LC: QA76.6.K64.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double &U1, &U2, two standard random normal deviates.
//
{
  double s;
  double x;
  double y;

  for ( ; ; )
  {
    x = r8_uniform_01 ( seed );
    y = r8_uniform_01 ( seed );
    x = 2.0 * x - 1.0;
    y = 2.0 * y - 1.0;
    s = x * x + y * y;

    if ( s <= 1.0 )
    {
      s = sqrt ( - 2.0 * log ( s ) / s );
      u1 = x * s;
      u2 = y * s;
      break;
    }
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

double *wshrt ( double d[], int n, int np, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    WSHRT returns a random Wishart variate.
//
//  Discussion:
//
//    This routine is a Wishart variate generator.  
//
//    On output, SA is an upper-triangular matrix of size NP * NP,
//    written in linear form, column ordered, whose elements have a 
//    Wishart(N, SIGMA) distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    Original FORTRAN77 version by William Smith, Ronald Hocking.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    William Smith, Ronald Hocking,
//    Algorithm AS 53, Wishart Variate Generator,
//    Applied Statistics,
//    Volume 21, Number 3, pages 341-345, 1972.
//
//  Parameters:
//
//    Input, double D[NP*(NP+1)/2], the upper triangular array that
//    represents the Cholesky factor of the correlation matrix SIGMA.
//    D is stored in column-major form.
//
//    Input, int N, the number of degrees of freedom.
//    1 <= N <= NP.
//
//    Input, int NP, the size of variables.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double WSHART[NP*(NP+1)/2], a sample from the 
//    Wishart distribution.
//
{
  double c;
  double df;
  int i;
  int ii;
  int ip;
  int j;
  int k;
  int nnp;
  int nq;
  int nr;
  int ns;
  double rn;
  double *sa;
  double *sb;
  double u1;
  double u2;

  k = 0;
  nnp = ( np * ( np + 1 ) ) / 2;
//
//  Load SB with independent normal (0, 1) variates.
//
  sb = new double[nnp];

  while ( k < nnp )
  {
    rnorm ( seed, u1, u2 );

    sb[k] = u1;
    k = k + 1;

    if ( k < nnp )
    {
      sb[k] = u2;
      k = k + 1;
    }
  }
//
//  Load diagonal elements with square root of chi-square variates.
//
  ns = 0;

  for ( i = 1; i <= np; i++ )
  {
    df = ( double ) ( np - i + 1 );
    ns = ns + i;
    u1 = 2.0 / ( 9.0 * df );
    u2 = 1.0 - u1;
    u1 = sqrt ( u1 );
//
//  Wilson-Hilferty formula for approximating chi-square variates:
//  The original code did not take the absolute value!
//
    sb[ns-1] = sqrt ( df * fabs ( pow ( u2 + sb[ns-1] * u1, 3 ) ) );
  }

  sa = new double[nnp];

  rn = ( double ) ( n );
  nr = 1;

  for ( i = 1; i <= np; i++ )
  {
    nr = nr + i - 1;
    for ( j = i; j <= np; j++ )
    {
      ip = nr;
      nq = ( j * ( j - 1 ) ) / 2 + i - 1;
      c = 0.0;
      for ( k = i; k <= j; k++ )
      {
        ip = ip + k - 1;
        nq = nq + 1;
        c = c + sb[ip-1] * d[nq-1];
      }
      sa[ip-1] = c;
    }
  }

  for ( i = 1; i <= np; i++ )
  {
    ii = np - i + 1;
    nq = nnp - np;
    for ( j = 1; j <= i; j++ )
    {
      ip = ( ii * ( ii - 1 ) ) / 2;
      c = 0.0;
      for ( k = i; k <= np; k++ )
      {
        ip = ip + 1;
        nq = nq + 1;
        c = c + sa[ip-1] * sa[nq-1];
      }
      sa[nq-1] = c / rn;
      nq = nq - 2 * np + i + j - 1;
    }
  }

  delete [] sb;

  return sa;
}
