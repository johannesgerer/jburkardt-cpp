# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "components.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMPONENTS_PRB.
//
//  Discussion:
//
//    COMPONENTS_PRB tests the COMPONENTS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "COMPONENTS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the COMPONENTS library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMPONENTS_PRB\n";
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
//    TEST01 tests I4VEC_COMPONENTS on a simple case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 28

  int a[N] = {
    0, 0, 1, 2, 4, 0, 0, 4, 0, 0,
    0, 8, 9, 9, 1, 2, 3, 0, 0, 5,
    0, 1, 6, 0, 0, 0, 4, 0 };
  int c[N];
  int component_num;
  int j;
  int n = N;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4VEC_COMPONENTS finds and labels connected\n";
  cout << "  components in a 1D integer vector.\n";

  cout << "\n";
  cout << "  A:\n";
  cout << "\n";
  cout << "    ";
  for ( j = 0; j < n; j++ )
  {
    cout << a[j];
  }
  cout << "\n";

  component_num = i4vec_components ( n, a, c );

  cout << "\n";
  cout << "  Number of components = " << component_num << "\n";
  cout << "\n";
  cout << "  C:\n";
  cout << "\n";
  cout << "    ";
  for ( j = 0; j < n; j++ )
  {
    cout << c[j];
  }
  cout << "\n";

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests I4MAT_COMPONENTS on a simple case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 9
# define N 17

  int a[M*N] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 1, 0, 1, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 1, 1, 1, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 0, 1, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 0, 1, 1, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 1, 0, 1, 1, 0,
    0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int *c;
  int component_num;
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  I4MAT_COMPONENTS finds and labels connected\n";
  cout << "  components in a 2D integer array.\n";

  cout << "\n";
  cout << "  A:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << "    ";
    for ( j = 0; j < n; j++ )
    {
      cout << a[i+j*m];
    }
    cout << "\n";
  }

  c = new int[m*n];

  component_num = i4mat_components ( m, n, a, c );

  cout << "\n";
  cout << "  Number of components = " << component_num << "\n";
  cout << "\n";
  cout << "  C:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << "    ";
    for ( j = 0; j < n; j++ )
    {
      cout << c[i+j*m];
    }
    cout << "\n";
  }

  delete [] c;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4BLOCK_COMPONENTS on a simple case.
//
//  Discussion:
//
//    This calculation is also done by a program called REGION.
//    The two programs differ in the number of components discovered
//    because REGION uses the full 3x3 block of pixels, resulting
//    in 26 potential neighbors, whereas I4BLOCK_COMPONENTS uses only
//    the north/south, east/west, up/down directions for 8 neighbors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2011
//
//  Author:
//
//    John Burkardt
//
{
# define L 64
# define M 64
# define N 26

  int a[L*M*N];
  int c[L*M*N];
  int component_num;
  string filename;
  int i;
  int i1;
  int *indices;
  int j;
  int j1;
  int k;
  int l = L;
  int m = M;
  int m1;
  int n = N;
  int n1;
  int *s;
  int s_total;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  I4BLOCK_COMPONENTS finds and labels connected\n";
  cout << "  components in a 3D integer block.\n";

  cout << "\n";
  cout << "  A is a 3D block of order " << l
       << " * " << m
       << " * " << n << "\n";

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        a[i+j*l+k*l*m] = 0;
      }
    }
  }
//
//  Retrieve the indices of nonzero data in A by reading a file.
//
  filename = "indices.txt";

  i4mat_header_read ( filename, &m1, &n1 );

  indices = i4mat_data_read ( filename, m1, n1 );

  for ( j1 = 0; j1 < n1; j1++ )
  {
    i = indices[0+j1*3] - 1;
    j = indices[1+j1*3] - 1;
    k = indices[2+j1*3] - 1;
    a[i+j*l+k*l*m] = 1;
  }

  delete [] indices;

  cout << "\n";
  cout << "  Number of nonzero A values is " << n1 << "\n";
//
//  Determine the components.
//
  component_num = i4block_components ( l, m, n, a, c );

  s = new int[component_num];

  for ( i = 0; i < component_num; i++ )
  {
    s[i] = 0;
  }

  cout << "\n";
  cout << "  Number of components = " << component_num << "\n";

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        if ( c[i+j*l+k*l*m] != 0 )
        {
          s[c[i+j*l+k*l*m]-1] = s[c[i+j*l+k*l*m]-1] + 1;
        }
      }
    }
  }

  cout << "\n";
  cout << "  Component  Size\n";
  cout << "\n";
  s_total = 0;
  for ( i = 0; i < component_num; i++ )
  {
    cout << "  " << setw(4) << i + 1
         << "  " << setw(8) << s[i] << "\n";
    s_total = s_total + s[i];
  }
  cout << "------  --------\n";
  cout << " Total  " << setw(8) << s_total << "\n";

  delete [] s;

  return;
# undef L
# undef M
# undef N
}
