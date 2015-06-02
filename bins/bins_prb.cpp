# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "bins.hpp"

int main ( );
void test0081 ( );
void test0835 ( );
void test0836 ( );
void test180 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BINS_PRB.
//
//  Discussion:
//
//    BINS_PRB tests the BINS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BINS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BINS library.\n";

  test0081 ( );
  test0835 ( );
  test0836 ( );
  test180 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BINS_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test0081 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0081 tests BIN_TO_R8_EVEN2, R8_TO_BIN_EVEN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  double a = 10.0;
  double b = 20.0;
  int bin;
  double c;
  double cmax;
  double cmin;
  int i;
  int nbin = 5;
  double rmax = 23.0;
  double rmin = 8.0;
  int seed;

  cout << "\n";
  cout << "TEST0676\n";
  cout << "  R8_TO_BIN_EVEN2 puts a number into a bin.\n";
  cout << "  BIN_TO_R8_EVEN2 returns the bin limits.\n";
  cout << "  The bins are equally spaced between A and B.\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  Total number of bins = " << nbin << "\n";
  cout << "\n";
  cout << "  Generate some random values C and put them in bins.\n";
  cout << "\n";
  cout << "       C      Bin   Bin_Min  Bin_Max\n";
  cout << "\n";

  seed = get_seed ( );

  for ( i = 1; i <= 30; i++ )
  {
    c = r8_uniform ( rmin, rmax, &seed );
    bin = r8_to_bin_even2 ( nbin, a, b, c );
    bin_to_r8_even2 ( nbin, bin, a, b, &cmin, &cmax );
    cout << "  " << setw(10) << c
         << "  " << setw(6)  << bin
         << "  " << setw(10) << cmin
         << "  " << setw(10) << cmax << "\n";
  }

  return;
}
//****************************************************************************80

void test0835 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0835 tests R82VEC_PART_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double a[N*2];
  double ahi[2] = { 10.0, 3.0 };
  double alo[2] = {  0.0, 2.0 };
  int i;
  int l;
  int r;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST0835\n";
  cout << "  R82VEC_PART_QUICK_A reorders an R2 vector\n";
  cout << "    as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  r82vec_uniform ( N, alo, ahi, &seed, a );

  r82vec_print ( N, a, "  Before rearrangment:" );

  r82vec_part_quick_a ( N, a, &l, &r );

  cout << "\n";
  cout << "  Rearranged array\n";
  cout << "  Left index =  " << l << "\n";
  cout << "  Key index =   " << l+1 << "\n";
  cout << "  Right index = " << r << "\n";

  r82vec_print ( l,     a,         "  Left half:" );
  r82vec_print ( 1,     a+2*l,     "  Key:" );
  r82vec_print ( N-l-1, a+2*(l+1), "  Right half:" );

  return;
# undef N
}
//****************************************************************************80

void test0836 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0836 tests R82VEC_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double a[N*2];
  double ahi[2] = { 10.0, 3.0 };
  double alo[2] = {  0.0, 2.0 };
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST0836\n";
  cout << "  D2VEC_SORT_QUICK_A sorts an R2 vector\n";
  cout << "    as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  r82vec_uniform ( N, alo, ahi, &seed, a );
//
//  For better testing, give a few elements the same first component.
//
  a[2*(3-1)+0] = a[2*(5-1)+0];
  a[2*(4-1)+0] = a[2*(12-1)+0];
//
//  Make two entries equal.
//
  a[2*(7-1)+0] = a[2*(11-1)+0];
  a[2*(7-1)+1] = a[2*(11-1)+1];

  r82vec_print ( N, a, "  Before sorting:" );

  r82vec_sort_quick_a ( N, a );

  r82vec_print ( N, a, "  Sorted array:" );

  return;
# undef N
}
//****************************************************************************80

void test180 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST180 tests SORT_HEAP_EXTERNAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int a[N];
  int i;
  int indx;
  int isgn;
  int itemp;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST180\n";
  cout << "  SORT_HEAP_EXTERNAL sorts objects externally.\n";

  indx = 0;
  i = 0;
  j = 0;
  isgn = 0;
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i] = i4_uniform ( 1, N, &seed );
  }

  i4vec_print ( N, a, "  Unsorted array:" );
//
//  Call the sort routine over and over.
//
  for ( ;; )
  {
    sort_heap_external ( N, &indx, &i, &j, isgn );
//
//  If the return value of INDX is negative, we're asked to compare
//  array elements I and J;
//
    if ( indx < 0 )
    {
      if ( a[i] <= a[j] )
      {
        isgn = -1;
      }
      else
      {
        isgn = 1;
      }

    }
//
//  ...and if the return value of INDX is positive, we're asked to switch
//  array elements I and J;
//
    else if ( 0 < indx )
    {
      i4_swap ( &a[i], &a[j] );
//
//  ...and if the return value of INDX is 0, we're done.
//
    }
    else
    {
      break;
    }

  }

  i4vec_print ( N, a, "  Sorted array:" );

  return;

# undef N
}
