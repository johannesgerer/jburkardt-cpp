# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>

using namespace std;

# include "pdflib.hpp"
# include "rnglib.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] );
double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PDFLIB_PRB.
//
//  Discussion:
//
//    PDFLIB_PRB tests the PDFLIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PDFLIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PDFLIB library.\n";
//
//  Initialize the random number generator package.
//
  initialize ( );

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PDFLIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests R8MAT_POFAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double diff;
  int i;
  int j;
  int n = 5;
  double *r1;
  double *r2;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  R8MAT_POFAC computes the Cholesky factor R of a\n";
  cout << "  positive definite matrix A, so that A = R' * R.\n";
  cout << "\n";
  cout << "  Start with random R1;\n";
  cout << "  Compute A = R1' * R1.\n";
  cout << "  Call R8MAT_POFAC and see if you recover R2 = R1.\n";
//
//  Generate a random upper triangular matrix with positive diagonal.
//
  r1 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      if ( i == j )
      {
        r1[i+j*n] = fabs ( r8_uniform_01_sample ( ) );
      }
      else
      {
        r1[i+j*n] = r8_uniform_01_sample ( );
      }
    }
    for ( i = j + 1; i < n; i++ )
    {
      r1[i+j*n] = 0.0;
    }
  }
  r8mat_print ( n, n, r1, "  R1:" );
//
//  Compute a positive definite symmetric matrix A.
//
  a = r8mat_mtm_new ( n, n, n, r1, r1 );

  r8mat_print ( n, n, a, "  A:" );

  r2 = r8mat_pofac ( n, a );

  diff = r8mat_norm_fro_affine ( n, n, r1, r2 );

  cout << "\n";
  cout << "  Frobenius difference between R1 and R2 = " << diff << "\n";

  delete [] a;
  delete [] r1;
  delete [] r2;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8VEC_MULTINORMAL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double c_det;
  double *c_inv;
  double eps;
  int i;
  int j;
  double *mu;
  int n = 5;
  double pdf1;
  double pdf2;
  double pi = 3.141592653589793;
  double *r1;
  double *r2;
  double *x;
  double xcx;
  double *y;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R8VEC_MULTINORMAL_PDF evaluates the PDF for the\n";
  cout << "  multinormal distribution.\n";
  cout << "\n";
  cout << "  The covariance matrix is C.\n";
  cout << "  The definition uses the inverse of C;\n";
  cout << "  R8VEC_MULTINORMAL_PDF uses the Cholesky factor.\n";
  cout << "  Verify that the algorithms are equivalent.\n";
//
//  Generate a random upper triangular matrix with positive diagonal.
//
  r1 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      if ( i == j )
      {
        r1[i+j*n] = fabs ( r8_uniform_01_sample ( ) );
      }
      else
      {
        r1[i+j*n] = r8_uniform_01_sample ( );
      }
    }
    for ( i = j + 1; i < n; i++ )
    {
      r1[i+j*n] = 0.0;
    }
  }
  r8mat_print ( n, n, r1, "  R1:" );
//
//  Compute a positive definite symmetric matrix C.
//
  c = r8mat_mtm_new ( n, n, n, r1, r1 );
  r8mat_print ( n, n, c, "  C:" );
//
//  Compute the Cholesky factor.
//
  r2 = r8mat_pofac ( n, c );
  r8mat_print ( n, n, r2, "  R2:" );
//
//  Compute the determinant of C.
//
  c_det = r8mat_podet ( n, r2 );
  cout << "\n";
  cout << "  Determinant of C = " << c_det << "\n";
//
//  Compute the inverse of C.
//
  c_inv = r8mat_poinv ( n, r2 );
//
//  Compute a random set of means.
//
  mu = new double[n];
  for ( i = 0; i < n; i++ )
  {
    mu[i] = r8_normal_01_sample ( );
  }
//
//  Compute X as small variations from MU.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    eps = 0.01 * r8_normal_01_sample ( );
    x[i] = ( 1.0 + eps ) * mu[i];
  }
//
//  Compute PDF1 from the function.
//
  pdf1 = r8vec_multinormal_pdf ( n, mu, r2, c_det, x );
//
//  Compute PDF2 from the definition.
//
  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i] - mu[i];
  }

  xcx = r8mat_vtmv ( n, n, y, c_inv, y );

  pdf2 = 1.0 / sqrt ( pow ( 2.0 * pi, n ) ) 
    * 1.0 / sqrt ( c_det ) 
    * exp ( - 0.5 * xcx );

  cout << "\n";
  cout << "  PDF1 = " << pdf1 << "\n";
  cout << "  PDF2 = " << pdf2 << "\n";

  delete [] c;
  delete [] c_inv;
  delete [] mu;
  delete [] r1;
  delete [] r2;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 calls R8_CHI_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  double df;
  int g;
  int i;
  double u;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  R8_CHI_SAMPLE ( DF ) samples the Chi distribution with\n";
  cout << "  DF degrees of freedom.\n";
//
//  Initialize the package.
//
  cout << "\n";
  cout << "  INITIALIZE initializes the random number generator.\n";
  cout << "  It only needs to be called once before using the package.\n";

  initialize ( );
//
//  Set the current generator index to #2, which (this being C++) has index 1!.
//
  g = 1;
  cgn_set ( g );
  cout << "\n";
  cout << "  Current generator index = " << g << "\n";
//
//  Repeatedly call R8_CHI_SAMPLE ( DF ).
//
  cout << "\n";
  cout << "   I       DF       R8_CHI_SAMPLE ( DF )\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    df = 5.0 * r8_uniform_01_sample ( ) + 1.0;
    u = r8_chi_sample ( df );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << df
         << "  " << setw(14) << u << "\n";
  }

  return;
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

double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTM_NEW computes C = A' * B.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
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
//    07 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MTM_NEW[N1*N3], the product matrix C = A' * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[k+i*n2] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A1[M*N], A2[M,N], the matrice for whose difference the 
//    Frobenius norm is desired.
//
//    Output, double R8MAT_NORM_FRO_AFFINE, the Frobenius norm of A1 - A2.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a1[i+j*m] - a2[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
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
//    20 August 2010
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
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

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
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

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

double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_VTMV multiplies computes the scalar x' * A * y.
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
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of
//    the matrix.
//
//    Input, double X[N], the first vector factor.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, double Y[M], the second vector factor.
//
//    Output, double R8MAT_VTMV, the value of X' * A * Y.
//
{
  int i;
  int j;
  double vtmv;

  vtmv = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      vtmv = vtmv + x[i] * a[i+j*m] * y[j];
    }
  }
  return vtmv;
}
