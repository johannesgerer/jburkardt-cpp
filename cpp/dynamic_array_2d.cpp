# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

int main ( int argc, char *argv[] );
void i4pp_delete ( int **a, int m, int n );
int **i4pp_new ( int m, int n );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    DYNAMIC_ARRAY_2D illustrates the creation of a dynamic 2D array.
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
{
  int **a;
  int **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  timestamp ( );
  cout << "\n";
  cout << "DYNAMIC_ARRAY_2D:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This program illustrates one way to create a dynamic 2D array\n";
  cout << "  in C++, that is, an array whose size and shape can be specified\n";
  cout << "  at run time, and whose entries can then be addressed using the\n";
  cout << "  notation \"a[i][j]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by I4PP_DELETE and new memory
//  requested by I4PP_NEW.
//
  m = 4;
  n = 5;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << m << " by " << n << ".\n";

  a = i4pp_new ( m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 10 * i + j;
    }
  }
//
//  Print A.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix A:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(2) << a[i][j];
    }
    cout << "\n";
  }
//
//  Create a new matrix B to store A' * A.
//
  b = i4pp_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
//
//  Print the matrix.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix B = A' * A:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << b[i][j];
    }
    cout << "\n";
  }
//
//  Free memory.
//
  i4pp_delete ( a, m, n );
  i4pp_delete ( b, n, n );
//
//  Terminate.
//
  cout << "\n";
  cout << "DYNAMIC_ARRAY_2D:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i4pp_delete ( int **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4PP_DELETE frees memory associated with an I4PP.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      int **a;
//      a = i4pp_new ( m, n );
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
//    Input, int **A, a pointer to the pointers to the array.
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

int **i4pp_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4PP_NEW allocates a new I4PP.
//
//  Discussion:
//
//    A declaration of the form
//      int **a;
//    is necesary.  Then an assignment of the form:
//      a = i4pp_new ( m, n );
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
//    Output, int **I4PP_NEW, a pointer to the pointers to the array.
//
{
  int **a;
  int i;

  a = new int *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "I4PP_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new int[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "I4PP_NEW - Fatal error!\n";
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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
