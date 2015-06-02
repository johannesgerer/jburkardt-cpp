# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void r8cmat_delete ( double **a, int m, int n );
double **r8cmat_new ( int m, int n );
void r8rmat_delete ( double **a, int m, int n );
double **r8rmat_new ( int m, int n );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POINTERS.
//
//  Discussion:
//
//    This program demonstrates the use of pointers to define objects
//    that can be indexed like a two-dimensional array.  Either row-major
//    or column-major storage can be handled.  However, the column major
//    storage means that an expression like a[i][j] refers to the i-th
//    column, j-th row, in other words, what we would naturally refer
//    to as A(J,I) instead!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POINTERS\n";
  cout << "  C++ version\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POINTERS:\n";
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
//    TEST01 demonstrates R8RMAT_NEW and R8RMAT_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double **a;
  double **b;
  int c;
  int ma = 4;
  int mb;
  int na = 3;
  int nb;
  int r;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  R8RMAT_NEW creates row-major doubly dimensioned arrays.\n";
  cout << "  R8RMAT_DELETE deletes them.\n";
//
//  Create the matrix.
//
  a = r8rmat_new ( ma, na );
//
//  Put some stuff in it.
//
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      a[r][c] = 10.0 * r + c;
    }
  }
//
//  Print the matrix.
//
  cout << "\n";
  cout << "  Matrix A:\n";
  cout << "\n";
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      cout << "  " << setw(9) << a[r][c];
    }
    cout << "\n";
  }
//
//  Create the transpose.
//
  mb = na;
  nb = ma;

  b = r8rmat_new ( mb, nb );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      b[c][r] = a[r][c];
    }
  }
//
//  Print the transpose;
//
  cout << "\n";
  cout << "  Matrix B = transpose ( A ):\n";
  cout << "\n";
  for ( r = 0; r < mb; r++ )
  {
    for ( c = 0; c < nb; c++ )
    {
      cout << "  " << setw(9) << b[r][c];
    }
    cout << "\n";
  }
//
//  To print row R,
//  we need to access row R, column C, so realize that 
//  B[R] is a pointer to the first entry, so 
//  B[R]+C is a pointer to the entry B(R,C) and 
//  *(B[R]+C) is its value.
//  These are CONSECUTIVE memory locations.
//
  r = 1;
  cout << "\n";
  cout << "  Row " << r << " of matrix B:\n";
  cout << "  (consecutive memory locations)\n";
  cout << "\n";
  for ( c = 0; c < nb; c++ )
  {
    cout << "  " << setw(9) << *(b[r]+c);
  }
  cout << "\n";
//
//  To print column C, we end up using a similar formula,
//  but because B[R] changes on each step, these are NOT consecutive
//  memory locations.
//
  c = 2;
  cout << "\n";
  cout << "  Column " << c << " of matrix B:\n";
  cout << "\n";
  for ( r = 0; r < mb; r++ )
  {
    cout << "  " << setw(9) << *(b[r]+c);
  }
  cout << "\n";
//
//  Free memory.
//
  r8rmat_delete ( a, ma, na );
  r8rmat_delete ( b, mb, nb );

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates R8CMAT_NEW and R8CMAT_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double **a;
  double **b;
  int c;
  int ma = 4;
  int mb;
  int na = 3;
  int nb;
  int r;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R8CMAT_NEW creates column-major doubly dimensioned arrays.\n";
  cout << "  R8CMAT_DELETE deletes them.\n";
  cout << "\n";
  cout << "  Unfortunately, A(i,j) is referenced as A[J][I]!\n";
//
//  Create the matrix.
//
  a = r8cmat_new ( ma, na );
//
//  Put some stuff in it.
//
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      a[c][r] = 10.0 * r + c;
    }
  }
//
//  Print the matrix.
//
  cout << "\n";
  cout << "  Matrix A:\n";
  cout << "\n";
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      cout << "  " << setw(9) << a[c][r];
    }
    cout << "\n";
  }
//
//  Create the transpose.
//
  mb = na;
  nb = ma;

  b = r8cmat_new ( mb, nb );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      b[r][c] = a[c][r];
    }
  }
//
//  Print the transpose;
//
  cout << "\n";
  cout << "  Matrix B = transpose ( A ):\n";
  cout << "\n";
  for ( r = 0; r < mb; r++ )
  {
    for ( c = 0; c < nb; c++ )
    {
      cout << "  " << setw(9) << b[c][r];
    }
    cout << "\n";
  }
//
//  To print row R, B[C] changes on each step, these are NOT consecutive
//  memory locations.
//
  r = 1;
  cout << "\n";
  cout << "  Row " << r << " of matrix B:\n";
  cout << "\n";
  for ( c = 0; c < nb; c++ )
  {
    cout << "  " << setw(9) << *(b[c]+r);
  }
  cout << "\n";
//
//  To print column C,
//  we need to access row R, column C, so realize that 
//  B[C] is a pointer to the first entry of the column, so 
//  B[C]+R is a pointer to the entry B(R,C) and 
//  *(B[C]+R) is its value.
//  These are CONSECUTIVE memory locations.
//
  c = 2;
  cout << "\n";
  cout << "  Column " << c << " of matrix B:\n";
  cout << "  (consecutive memory locations)\n";
  cout << "\n";
  for ( r = 0; r < mb; r++ )
  {
    cout << "  " << setw(9) << *(b[c]+r);
  }
  cout << "\n";
//
//  Free memory.
//
  r8cmat_delete ( a, ma, na );
  r8cmat_delete ( b, mb, nb );

  return;
}
//****************************************************************************80

void r8cmat_delete ( double **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8CMAT_DELETE frees memory associated with an R8CMAT.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      double **a;
//      a = r8cmat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double **A, the pointer to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    delete [] a[j];
  }

  delete [] a;

  return;
}
//****************************************************************************80

double **r8cmat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8CMAT_NEW allocates a new R8CMAT.
//
//  Discussion:
//
//    An R8CMAT is a 2D double precision real column-major matrix.
//
//    A declaration of the form
//      double **a;
//    is necesary.  Then an assignment of the form:
//      a = r8cmat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//    Because of the column-major format, the interpretation of the
//    matrix indices is reversed, so that a[i][j] refers to column i,
//    row j!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, double **R8CMAT_NEW, a new matrix.
//
{
  double **a;
  int j;

  a = new double *[n];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8CMAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    a[j] = new double[m];
    if ( a[j] == NULL )
    {
      cerr << "\n";
      cerr << "R8CMAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}
//****************************************************************************80

void r8rmat_delete ( double **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_DELETE frees memory associated with an R8RMAT.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      double **a;
//      a = r8rmat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double **A, the pointer to the array.
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

double **r8rmat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_NEW allocates a new R8RMAT.
//
//  Discussion:
//
//    An R8RMAT is a 2D double precision real row-major matrix.
//
//    A declaration of the form
//      double **a;
//    is necesary.  Then an assignment of the form:
//      a = r8rmat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, double **R8RMAT_NEW, a new matrix.
//
{
  double **a;
  int i;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8RMAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8RMAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
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
