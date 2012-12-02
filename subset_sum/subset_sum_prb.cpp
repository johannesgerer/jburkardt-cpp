# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "subset_sum.hpp"

int main ( );
int test01 ( int n, int w[], int t, int ind_min, int ind_max );
void test02 ( int n, int w[], int t, int ind_min, int ind_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUBSET_SUM_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  int ind;
  int ind_max;
  int ind_min;
  int n;
  int t;
  int test;
  int test_num = 9;
  int *w;
  static int w01[8] = { 15, 22, 14, 26, 32, 9, 16, 8 };
  static int w02[8] = { 15, 22, 14, 26, 32, 9, 16, 8 };
  static int w03[8] = { 15, 22, 14, 26, 32, 9, 16, 8 };
  static int w04[10] = { 267,  493,  869,  961, 1000, 1153, 1246, 1598, 1766, 1922 };
  static int w05[21] = {518533, 1037066, 2074132, 1648264, 796528, 
             1593056,  686112, 1372224,  244448, 488896, 
              977792, 1955584, 1411168,  322336, 644672, 
             1289344,   78688,  157376,  314752, 629504, 
             1259008};
  static int w06[10] = { 41, 34, 21, 20,  8,  7,  7,  4,  3,  3 };
  static int w07[9] = { 81, 80, 43, 40, 30, 26, 12, 11, 9 };
  static int w08[6] = { 1, 2, 4, 8, 16, 32};
  static int w09[10] = { 25, 27, 3, 12, 6, 15, 9, 30, 21, 19 };

  timestamp ( );
  cout << "\n";
  cout << "SUBSET_SUM_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the SUBSET_SUM library.\n";
//
//  Find individual solutions.
//
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w01 );
      t = 53;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 2 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w02 );
      t = 53;
      ind_min = ind + 1;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 3 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w03 );
      t = 53;
      ind_min = ind + 1;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w04 );
      t = 5842;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 5 )
    {
      n = 21;
      w = i4vec_copy_new ( n, w05 );
      t = 2463098;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 6 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w06 );
      t = 50;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 7 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w07 );
      t = 100;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 8 )
    {
      n = 6;
      w = i4vec_copy_new ( n, w08 );
      t = 22;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 9 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w09 );
      t = 50;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }

    ind = test01 ( n, w, t, ind_min, ind_max );

    delete [] w;
  }
//
//  Simply count solutions.
//
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w01 );
      t = 53;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 2 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w02 );
      t = 53;
      ind_min = 68;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 3 )
    {
      n = 8;
      w = i4vec_copy_new ( n, w03 );
      t = 53;
      ind_min = 167;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w04 );
      t = 5842;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 5 )
    {
      n = 21;
      w = i4vec_copy_new ( n, w05 );
      t = 2463098;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 6 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w06 );
      t = 50;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 7 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w07 );
      t = 100;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 8 )
    {
      n = 6;
      w = i4vec_copy_new ( n, w08 );
      t = 22;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }
    else if ( test == 9 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w09 );
      t = 50;
      ind_min = 0;
      ind_max = i4_power ( 2, n ) - 1;
    }

    test02 ( n, w, t, ind_min, ind_max );

    delete [] w;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "SUBSET_SUM_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int test01 ( int n, int w[], int t, int ind_min, int ind_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 seeks a subset of a set that has a given sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, int W[N], a set of weights.  The length of this
//    array must be no more than 31.
//
//    Input, int T, the target value.
//
//    Input, int IND_MIN, IND_MAX, the lower and upper
//    limits to be searched.
//
//    Output, int TEST01, the index of a solution, if found,
//    or the value -1 otherwise.
//
{
  int *c;
  int i;
  int ind;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Seek a subset of W that sums to T.\n";
  cout << "\n";
  cout << "  Target value T = " << t << "\n";
  i4vec_print ( n, w, "  Weight vector W:" );

  c = subset_sum_find ( n, w, t, ind_min, ind_max, ind );

  if ( ind == -1 )
  {
    cout << "\n";
    cout << "  No solution was found.\n";
  }
  else
  {
    cout << "\n";
    cout << "  Solution index = " << ind << "\n";
    i4vec_print ( n, c, "  Solution:" );
  }
  delete [] c;

  return ind;
}
//****************************************************************************80

void test02 ( int n, int w[], int t, int ind_min, int ind_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 counts solutions to the subset sum problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, int W[N], a set of weights.  The length of this
//    array must be no more than 31.
//
//    Input, int T, the target value.
//
//    Input, int IND_MIN, IND_MAX, the lower and upper
//    limits to be searched.
//
{
  int i;
  int solution_num;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Count solutions to the subset sum problem.\n";
  cout << "\n";
  cout << "  Target value T = " << t << "\n";
  i4vec_print ( n, w, "  Weight vector W:" );

  solution_num = subset_sum_count ( n, w, t, ind_min, ind_max );

  cout << "\n";
  cout << "  Number of solutions = " << solution_num << "\n";

  return;
}
