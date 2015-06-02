# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

int main ( );
void make_arrays ( int m, int **a, int n, int **b );
void i4vec_print ( int n, int a[], string title );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ARRAY_RETURN.
//
//  Discussion:
//
//    The correct form of this program was worked out with the somewhat
//    bemused assistance of Miro Stoyanov.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int *b;
  int m;
  int n;

  cout << "\n";
  cout << "ARRAY_RETURN:\n";
  cout << "  C++ version\n";
  cout << "  Create two arrays in a function, \n";
  cout << "  return them in the argument list.\n";
//
//  Specify the size of the arrays to be created.
//
  m = 10;
  n = 5;

  make_arrays ( m, &a, n, &b );
//
//  Verify that the arrays were created and transferred properly.
//
  i4vec_print ( m, a, "  A as received by main:" );
  i4vec_print ( n, b, "  B as received by main:" );
//
//  Free memory.
//
  delete [] a;
  delete [] b;
//
//  Terminate.
//
  cout << "\n";
  cout << "ARRAY_RETURN:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void make_arrays ( int m, int **a, int n, int **b )

//****************************************************************************80
//
//  Purpose:
//
//    MAKE_ARRAYS creates, sets, and returns two arrays using the argument list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the desired size of the first array.
//
//    Output, int **A, the first array.
//
//    Input, int N, the desired size of the second array.
//
//    Output, int **B, the second array.
//
{
  int i;
//
//  We create PA and PB simply as a convenience.
//  They make the code a little more readable.
//
  int *pa;
  int *pb;
  
  ( *a ) = new int[m];
  pa = *a;

  for ( i = 0; i < m; i++ )
  {
    pa[i] = 10 + i;
  }
  i4vec_print ( m, pa, "  A as defined in MAKE_ARRAYS:" );

  ( *b ) = new int[n];
  pb = *b;

  for ( i = 0; i < n; i++ )
  {
    pb[i] = 100 + 2 * i;
  }
  i4vec_print ( n, pb, "  B as defined in MAKE_ARRAYS:" );

  return;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
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
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
