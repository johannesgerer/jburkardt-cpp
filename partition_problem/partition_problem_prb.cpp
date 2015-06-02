# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "partition_problem.hpp"

int main ( );
void test01 ( int n, int w[] );
void test02 ( int n, int w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PARTITION_PROBLEM_PRB.
//
//  Discussion:
//
//    PARTITION_PROBLEM_PRB tests the PARTITION_PROBLEM library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int test;
  int test_num = 5;
  int *w;
  int w1[5] = { 19, 17, 13, 9, 6 };
  int w2[9] = { 484, 114, 205, 288, 506, 503, 201, 127, 410 };
  int w3[10] = { 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 };
  int w4[10] = { 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 };
  int w5[9] = { 3, 4, 3, 1, 3, 2, 3, 2, 1 };

  timestamp ( );
  cout << "\n";
  cout << "PARTITION_PROBLEM_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PARTITION_PROBLEM library.\n";
//
//  Find individual solutions.
//
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 5;
      w = i4vec_copy_new ( n, w1 );
    }
    else if ( test == 2 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w2 );
    }
    else if ( test == 3 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w3 );
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w4 );
    }
    else if ( test == 5 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w5 );
    }

    test01 ( n, w );

    delete [] w;
  }
//
//  Count solutions.
//
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 5;
      w = i4vec_copy_new ( n, w1 );
    }
    else if ( test == 2 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w2 );
    }
    else if ( test == 3 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w3 );
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w4 );
    }
    else if ( test == 5 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w5 );
    }

    test02 ( n, w );

    delete [] w;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PARTITION_PROBLEM_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int n, int w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PARTITION_BRUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, int W[N], a set of weights.
//
{
  int *c;
  int discrepancy;
  int i;
  int w0_sum;
  int w1_sum;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Partition a set of N integers W so that the subsets\n";
  cout << "  have equal sums.\n";

  c = new int[n];

  partition_brute ( n, w, c, discrepancy );

  cout << "\n";
  cout << "     I        W0        W1\n";
  cout << "\n";
  w0_sum = 0;
  w1_sum = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( c[i] == 0 )
    {
      w0_sum = w0_sum + w[i];
      cout << "  " << setw(4) << i
           << "  " << setw(8) << w[i] << "\n";
    }
    else
    {
      w1_sum = w1_sum + w[i];
      cout << "  " << setw(4) << i
           << "  " << "        "
           << "  " << setw(8) << w[i] << "\n";
    }
  }
  cout << "        --------  --------\n";
  cout << "  " << "    "
       << "  " << setw(8) << w0_sum
       << "  " << setw(8) << w1_sum << "\n";
  cout << "\n";
  cout << "  Discrepancy = " << setw(8) << discrepancy << "\n";

  delete [] c;

  return;
}
//****************************************************************************80

void test02 ( int n, int w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PARTITION_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, int W[N], a set of weights.
//
{
  int count;
  int i;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  PARTITION_COUNT counts the number of exact solutions\n";
  cout << "  of the partition problem.\n";

  count = partition_count ( n, w );

  cout << "\n";
  cout << "     I        W\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << w[i] << "\n";
  }
  cout << "\n";
  cout << "  Number of solutions = " << count << "\n";

  return;
}
