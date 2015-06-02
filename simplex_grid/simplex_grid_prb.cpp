# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "simplex_grid.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_TEST tests the SIMPLEX_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SIMPLEX_GRID_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the SIMPLEX_GRID library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SIMPLEX_GRID_TEST:\n";
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
//    TEST01 tests SIMPLEX_GRID_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int ng;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  SIMPLEX_GRID_SIZE counts the points in a regular grid\n";
  cout << "  with N+1 points on a side, in an M-dimensional simplex.\n";
  cout << "\n";
  cout << "        M: 0     1     2     3     4     5\n";
  cout << "    N:\n";
  for ( n = 0; n <= 10; n++ )
  {
    cout << "  " << setw(3) << n << ":";
    for ( m = 0; m <= 5; m++ )
    {
      ng = simplex_grid_size ( m, n );
      cout << setw(6) << ng;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SIMPLEX_GRID_INDEX_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *g;
  int i;
  int j;
  int m = 3;
  int n;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  SIMPLEX_GRID_INDEX_NEXT lists, one by one, the indices\n";
  cout << "  of a simplex grid that uses N+1 points on a side,\n";
  cout << "  in an M-dimensional simplex.\n";
  cout << "\n";
  cout << "   #:  1  2  3  (*)\n";
  cout << "\n";

  n = 3;

  j = 0;
  g = new int[m+1];
  for ( i = 0; i < m; i++ )
  {
    g[i] = 0;
  }
  g[m] = n;
  
  while ( 1 )
  {
    cout << "  " << setw(2) << j;
    for ( i = 0; i < m; i++ )
    {
      cout << setw(3) << g[i];
    }
    cout << " (" << setw(3) << g[m] << ")\n";

    if ( g[0] == n )
    {
      break;
    }

    simplex_grid_index_next ( m, n, g );

    j = j + 1;
  }

  delete [] g;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests SIMPLEX_GRID_INDEX_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *g;
  int i;
  int j;
  int m = 3;
  int n;
  int seed;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  SIMPLEX_GRID_INDEX_SAMPLE returns a randomly selected\n";
  cout << "  index of a simplex grid that uses N+1 points on a side,\n";
  cout << "  in an M-dimensional simplex.\n";
  cout << "\n";
  cout << "   #:  1  2  3  (*)\n";
  cout << "\n";

  n = 3;
  seed = 123456789;

  for ( j = 1; j <= 20; j++ )
  {
    g = simplex_grid_index_sample ( m, n, seed );

    cout << "  " << setw(2) << j << ":";
    for ( i = 0; i < m; i++ )
    {
      cout << setw(3) << g[i];
    }
    cout << " (" << setw(3) << g[m] << ")\n";

    delete [] g;
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests SIMPLEX_GRID_INDEX_TO_POINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *g;
  int i;
  int j;
  int m = 2;
  int n;
  int seed;
  double v[2*3] = {
    20.0,  0.0, 
    30.0, 40.0, 
    10.0, 20.0 };
  double *x;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point\n";
  cout << "  corresponding to a grid index of a simplex grid that\n";
  cout << "  that uses N+1 points on a side,\n";
  cout << "  in an M-dimensional simplex.\n";

  n = 5;

  r8mat_transpose_print ( m, m + 1, v, "  Simplex vertices:" );

  cout << "\n";
  cout << "  Choosing random simplex indices to convert:\n";
  cout << "   #:  1  2  3     X        Y\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 10; j++ )
  {
    g = simplex_grid_index_sample ( m, n, seed );
    x = simplex_grid_index_to_point ( m, n, 1, g, v );

    cout << "  " << setw(2) << j << ":";
    for ( i = 0; i <= m; i++ )
    {
      cout << setw(3) << g[i];
    }
    cout << "  ";
    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    cout << "\n";

    delete [] g;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests SIMPLEX_GRID_INDEX_ALL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *grid;
  int m;
  int n;
  int ng;

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  SIMPLEX_GRID_INDEX_ALL returns all the indices\n";
  cout << "  of a simplex grid that uses N+1 points on a side,\n";
  cout << "  in an M-dimensional simplex.\n";

  m = 3;
  n = 3;
  ng = simplex_grid_size ( m, n );

  grid = simplex_grid_index_all ( m, n, ng );

  i4mat_transpose_print ( m + 1, ng, grid, 
    "  Transposed Simplex Grid Index Matrix:" );

  delete [] grid;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests SIMPLEX_GRID_INDEX_TO_POINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *grid;
  int m = 2;
  int n;
  int ng;
  double v[2*3] = {
    20.0,  0.0, 
    30.0, 40.0, 
    10.0, 20.0 };
  double *x;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point\n";
  cout << "  corresponding to a grid index of a simplex grid that\n";
  cout << "  that uses N+1 points on a side,\n";
  cout << "  in an M-dimensional simplex.\n";

  n = 5;
  ng = simplex_grid_size ( m, n );

  r8mat_transpose_print ( m, m + 1, v, "  Simplex vertices:" );

  grid = simplex_grid_index_all ( m, n, ng );

  x = simplex_grid_index_to_point ( m, n, ng, grid, v );

  r8mat_transpose_print ( m, ng, x, "  Grid Point Coordinates:" );

  delete [] grid;
  delete [] x;

  return;
}

