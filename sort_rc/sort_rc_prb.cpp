# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "sort_rc.hpp"

int main ( );
void sort_rc_test ( );
void sort_safe_rc_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SORT_RC_PRB.
//
//  Discussion:
//
//    SORT_RC_PRB tests the SORT_RC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SORT_RC_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SORT_RC library.\n";

  sort_rc_test ( );
  sort_safe_rc_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SORT_RC_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void sort_rc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_RC_TEST tests SORT_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int i;
  int i4_hi;
  int i4_lo;
  int indx;
  int isgn;
  int j;
  int k;
  int n = 20;
  int seed;

  cout << "\n";
  cout << "SORT_RC_TEST\n";
  cout << "  SORT_RC sorts objects externally.\n";
  cout << "  This function relies on the use of persistent\n";
  cout << "  data stored internally.\n";
//
//  Generate some data to sort.
//
  i4_lo = 1;
  i4_hi = n;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, i4_lo, i4_hi, seed );
 
  i4vec_print ( n, a, "  Unsorted array:" );
//
//  Sort the data.
// 
  indx = 0;

  for ( ; ; )
  {
    sort_rc ( n, indx, i, j, isgn );
 
    if ( indx < 0 )
    {
      isgn = 1;
      if ( a[i-1] <= a[j-1] )
      {
        isgn = -1;
      }
    }
    else if ( 0 < indx )
    {
      k      = a[i-1];
      a[i-1] = a[j-1];
      a[j-1] = k;
    }
    else
    {
      break;
    }
  }
//
//  Display the sorted data.
//
  i4vec_print ( n, a, "  Sorted array:" );
 
  delete [] a;

  return;
}
//****************************************************************************80

void sort_safe_rc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_SAFE_RC_TEST tests SORT_SAFE_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int i;
  int i_save;
  int i4_hi;
  int i4_lo;
  int indx;
  int isgn;
  int j;
  int j_save;
  int k;
  int k_save;
  int l_save;
  int n = 20;
  int n_save;
  int seed;

  cout << "\n";
  cout << "SORT_SAFE_RC_TEST\n";
  cout << "  SORT_SAFE_RC sorts objects externally.\n";
  cout << "  This version of the algorithm does not rely on\n";
  cout << "  internally saved or 'persistent' or 'static' memory.\n";
//
//  Generate some data to sort.
//
  i4_lo = 1;
  i4_hi = n;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, i4_lo, i4_hi, seed );
 
  i4vec_print ( n, a, "  Unsorted array:" );
//
//  Sort the data.
//
  indx = 0;

  for ( ; ; )
  {
    sort_safe_rc ( n, indx, i, j, isgn, i_save, j_save, k_save, l_save, 
      n_save );
 
    if ( indx < 0 )
    {
      isgn = 1;
      if ( a[i-1] <= a[j-1] )
      {
        isgn = -1;
      }
    }
    else if ( 0 < indx )
    {
      k      = a[i-1];
      a[i-1] = a[j-1];
      a[j-1] = k;
    }
    else
    {
      break;
    }
  }
//
//  Display the sorted data.
//
  i4vec_print ( n, a, "  Sorted array:" );
 
  delete [] a;

  return;
}
