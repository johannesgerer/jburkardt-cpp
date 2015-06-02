# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "index.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for INDEX_PRB.
//
//  Discussion:
//
//    INDEX_PRB tests the INDEX library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "INDEX_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the INDEX library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "INDEX_PRB:\n";
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
//    TEST01 tests INDEX0 and INDEX1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int value;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  INDEX0 indexes a 1D array with zero base,\n";
  cout << "  INDEX1 indexes a 1D array with  unit base.\n";
  cout << "\n";
  cout << "             Min Index   Max\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  cout << "  1D Index" 
       << "  " << setw(4) << i_min
       << "  " << setw(4) << i
       << "  " << setw(4) << i_max << "\n";

  value = index0 ( i_min, i, i_max );
  index_min = 0;
  index_max = index_min + i_max - i_min;
  cout << "  Index0  "
       << "  " << setw(4) << index_min
       << "  " << setw(4) << value
       << "  " << setw(4) << index_max << "\n";

  value = index1 ( i_min, i, i_max );
  index_min = 1;
  index_max = index_min + i_max - i_min; 
  cout << "  Index1  "
       << "  " << setw(4) << index_min
       << "  " << setw(4) << value
       << "  " << setw(4) << index_max << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests INDEX01, INDEX10, INDEX12 and INDEX21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int value;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a 2D array,\n";
  cout << "  INDEX01 column indexes with zero base,\n";
  cout << "  INDEX10 row indexes with zero base,\n";
  cout << "  INDEX12 column indexes with unit base,\n";
  cout << "  INDEX21 row indexes with unit base.\n";
  cout << "\n";
  cout << "                Min   Index     Max\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  cout << "  2D Index:"
       << "  " << setw(3) << i_min << setw(3) << j_min
       << "  " << setw(3) << i     << setw(3) << j
       << "  " << setw(3) << i_max << setw(3) << j_max << "\n";

  value = index01 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 0;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  cout << "  INDEX01: "
       << "  " << setw(6) << index_min
       << "  " << setw(6) << value
       << "  " << setw(6) << index_max << "\n";

  value = index10 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 0;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  cout << "  INDEX10: "
       << "  " << setw(6) << index_min
       << "  " << setw(6) << value
       << "  " << setw(6) << index_max << "\n";

  value = index12 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 1;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  cout << "  INDEX12: "
       << "  " << setw(6) << index_min
       << "  " << setw(6) << value
       << "  " << setw(6) << index_max << "\n";

  value = index21 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 1;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  cout << "  INDEX21: "
       << "  " << setw(6) << index_min
       << "  " << setw(6) << value
       << "  " << setw(6) << index_max << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests INDEX012, INDEX123, INDEX210, and INDEX321.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int m;
  int value;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a 3D array,\n";
  cout << "  INDEX012 column indexes with zero base,\n";
  cout << "  INDEX123 column indexes with unit base,\n";
  cout << "  INDEX210 row indexes with zero base.\n";
  cout << "  INDEX321 row indexes with unit base.\n";
  cout << "\n";
  cout << "                   Min      Index        Max\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;

  m = ( i_max - i_min + 1 ) 
    * ( j_max - j_min + 1 ) 
    * ( k_max - k_min + 1 );

  cout << "  3D Index:"
       << "  " << setw(3) << i_min << setw(3) << j_min << setw(3) << k_min
       << "  " << setw(3) << i     << setw(3) << j     << setw(3) << k
       << "  " << setw(3) << i_max << setw(3) << j_max << setw(3) << k_max << "\n";

  value = index012 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEX012:"
       << "  " << setw(9) << index_min
       << "  " << setw(9) << value
       << "  " << setw(9) << index_max << "\n";
 
  value = index123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEX123:"
       << "  " << setw(9) << index_min
       << "  " << setw(9) << value
       << "  " << setw(9) << index_max << "\n";

  value = index210 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEX210:"
       << "  " << setw(9) << index_min
       << "  " << setw(9) << value
       << "  " << setw(9) << index_max << "\n";

  value = index321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEX321:"
       << "  " << setw(9) << index_min
       << "  " << setw(9) << value
       << "  " << setw(9) << index_max << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests INDEX0123, INDEX1234, INDEX3210, and INDEX4321.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int l;
  int l_max;
  int l_min;
  int m;
  int value;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a 4D array,\n";
  cout << "  INDEX0123 column indexes with zero base,\n";
  cout << "  INDEX1234 column indexes with unit base,\n";
  cout << "  INDEX3210 row indexes with zero base,\n";
  cout << "  INDEX4321 row indexes with unit base.\n";
  cout << "\n";
  cout << "                       Min         Index           Max\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;

  m = ( i_max - i_min + 1 ) 
    * ( j_max - j_min + 1 ) 
    * ( k_max - k_min + 1 ) 
    * ( l_max - l_min + 1 );

  cout << "  4D Index:  "
       << "  " << setw(3) << i_min << setw(3) << j_min << setw(3) << k_min << setw(3) << l_min
       << "  " << setw(3) << i     << setw(3) << j     << setw(3) << k     << setw(3) << l
       << "  " << setw(3) << i_max << setw(3) << j_max << setw(3) << k_max << setw(3) << l_max << "\n";

  value = index0123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEX0123: "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = index1234 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEX1234: "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = index3210 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEX3210: "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = index4321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEX4321: "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests INDEX0N, INDEX1N, INDEXN0 and INDEXN1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//
{
  int i[4] = {3,2,1,2};
  int i_max[4] = {5,4,3,2};
  int i_min[4] = {1,1,1,1};
  int index_max;
  int index_min;
  int j;
  int m;
  int n = 4;
  int value;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For an N-dimensional array,\n";
  cout << "  INDEX0N column indexes with zero base,\n";
  cout << "  INDEX1N column indexes with unit base,\n";
  cout << "  INDEXN0 row indexes with zero base,\n";
  cout << "  INDEXN1 row indexes with unit base.\n";
  cout << "\n";
  cout << "                       Min         Index           Max\n";

  m = 1;
  for ( j = 0; j < n; j++ )
  {
    m = m * ( i_max[j] - i_min[j] + 1 );
  }

  cout << "  ND Index: "
       << "  " << setw(3) << i_min[0] << setw(3) << i_min[1] << setw(3) << i_min[2] << setw(3) << i_min[3]
       << "  " << setw(3) << i[0]     << setw(3) << i[1]     << setw(3) << i[2]     << setw(3) << i[3]
       << "  " << setw(3) << i_max[0] << setw(3) << i_max[1] << setw(3) << i_max[2] << setw(3) << i_max[3] << "\n";

  value = index0n ( n, i_min, i, i_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEX0N:  "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = index1n ( n, i_min, i, i_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEX1N:  "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = indexn0 ( n, i_min, i, i_max );
  index_min = 0;
  index_max = index_min + m - 1;
  cout << "  INDEXN0:  "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  value = indexn1 ( n, i_min, i, i_max );
  index_min = 1;
  index_max = index_min + m - 1;
  cout << "  INDEXN1:  "
       << "  " << setw(12) << index_min
       << "  " << setw(12) << value
       << "  " << setw(12) << index_max << "\n";

  return;
}

