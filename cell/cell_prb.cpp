# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cell.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN tests CELL.
//
//  Discussion:
//
//    An R8CVV is a "cell vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CELL_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the CELL library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CELL_PRB:\n";
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
//    TEST01 stores some of Pascal's triangle in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "cell array vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *ai;
  double aij;
  int col;
  int i;
  int in[4] = { 0, 1, 4, 4 };
  int j;
  int jn[4] = { 1, 2, 3, 7 };
  int m = 5;
  int mn;
  int nn;
  int nr[5] = { 4, 5, 6, 7, 8 };
  int nr_max;
  int nv;
  int *roff;
  int row;
  double *vn;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use a cell array (vector of vectors) to store rows 3:7\n";
  cout << "  of Pascal''s triangle.\n";

  i4vec_print ( m, nr, "  The row sizes:" );
//
//  From the NR information:
//  * determine the total size, MN
//
  mn = r8cvv_size ( m, nr );
  cout << "\n";
  cout << "  The storage for the cell array is " << mn << "\n";
//
//  Allocate the cell array.
//
  a = new double[mn];
//
//  Zero out the cell array.
//
  for ( i = 0; i < mn; i++ )
  {
    a[i] = 0.0;
  }
//
//  Allocate a vector big enough to hold any single row.
//
  nr_max = i4vec_max ( m, nr );
  ai = new double[nr_max];
//
//  From the NR information:
//  * determine the offsets.
//
  roff = r8cvv_offset ( m, nr );
  i4vec_print ( m + 1, roff, "  The row offsets:" );
//
//  Rows 1 through 5 of A will contain rows 3 through 7 of Pascal's triangle.
//  Set these values one row at a time.
//
  ai[0] = 1.0;

  for ( row = 1; row <= 7; row++ )
  {
    col = row + 1;
    ai[col-1] = ai[col-2];
    for ( col = row - 1; 1 <= col; col-- )
    {
      ai[col] = ai[col] + ai[col-1];
    }

    if ( 3 <= row )
    {
      i = row - 3;
      r8cvv_rset ( mn, a, m, roff, i, ai );
    }
  }
//
//  Print the cell array.
//
  r8cvv_print ( mn, a, m, roff, "  Rows 3:7 of Pascal's Triangle:" );
//
//  Retrieve the entry from cell array row 2, column 3:
//
  i = 2;
  j = 3;
  aij = r8cvv_iget ( mn, a, m, roff, i, j );
  cout << "\n";
  cout << "  A(" << i << "," << j << ") = " << aij << "\n";
//
//  Retrieve row 3:
//
  i = 3;
  ai = r8cvv_rget_new ( mn, a, m, roff, i );
  nv = roff[i+1] - roff[i];
  r8vec_transpose_print ( nv, ai, "  A(3,*):" );
//
//  Retrieve a list of entries.
//
  nn = 4;
  vn = r8cvv_nget_new ( mn, a, m, roff, nn, in, jn );
  cout << "\n";
  cout << "  Retrieve an arbitrary list of items:\n";
  cout << "\n";
  for ( i = 0; i < nn; i++ )
  {
    cout << "  A(" << in[i] << "," << jn[i] << ") = " << vn[i] << "\n";
  }
//
//  Free memory.
//
  delete [] a;
  delete [] ai;
  delete [] roff;
  delete [] vn;

  return;
}

