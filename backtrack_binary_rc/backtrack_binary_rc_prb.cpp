# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "backtrack_binary_rc.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BACKTRACK_BINARY_RC_PRB.
//
//  Discussion:
//
//    BACKTRACK_BINARY_RC_PRB tests BACKTRACK_BINARY_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BACKTRACK_BINARY_RC_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the BACKTRACK_BINARY_RC library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BACKTRACK_BINARY_RC_PRB:\n";
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
//    TEST01 seeks a selection of binary powers that have a given sum.
//
//  Discussion:
//
//    We consider the binary powers 1, 2, 4, ... 2^(n-1).
//
//    We wish to select some of these powers, so that the sum is equal
//    to a given target value.  We are actually simply seeking the binary
//    representation of an integer.
//
//    A partial solution is acceptable if it is less than the target value.
//
//    We list the powers in descending order, so that the bactracking
//    procedure makes the most significant choices first, thus quickly
//    eliminating many unsuitable choices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int call_num;
  int choice[8];
  int factor;
  int i;
  int n = 8;
  int n2;
  bool reject;
  int result;
  int target;
  int targets[3] = { 73, 299, -3 };
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use BACKBIN_RC to find the binary expansion of\n";
  cout << "  an integer between 0 and 255.\n";
  cout << "  The choices are 0/1 for the 8 digits.\n";

  for ( test = 0; test < test_num; test++ )
  {
    target = targets[test];
    cout << "\n";
    cout << "  TARGET = " << target << "\n";
    call_num = 0;
    n2 = -1;

    for ( ; ; )
    {
      backbin_rc ( n, reject, n2, choice );
      call_num = call_num + 1;

      if ( n2 == -1 )
      {
        cout << "  Termination without solution.\n";
        break;
      }
//
//  Evaluate the integer determined by the choices.
//
      factor = 1;
      for ( i = n; n2 < i; i-- )
      {
        factor = factor * 2;
      }

      result = 0;
      for ( i = 0; i < n2; i++ )
      {
        result = result * 2 + choice[i];
      }

      result = result * factor;
//
//  If the integer is too big, then we reject it, and
//  all the related integers formed by making additional choices.
//
      reject = ( target < result );
//
//  If we hit the target, then in this case, we can exit because
//  the solution is unique.
//
      if ( result == target )
      {
        break;
      }

    }

    cout << "  Number of calls = " << call_num << "\n";
    cout << "  Binary search space = " << i4_power ( 2, n ) << "\n";
    cout << "  ";
    for ( i = 0; i < n; i++ )
    {
      cout << setw(2) << choice[i];
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
//    TEST02 seeks a subset of a set of numbers which add to a given sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int call_num;
  int choice[8];
  int i;
  int n = 8;
  int n2;
  bool reject;
  int result;
  int target = 53;
  int test;
  int w[8] = { 15, 22, 14, 26, 32, 9, 16, 8 };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use BACKBIN_RC to seek subsets of a set W\n";
  cout << "  that sum to a given target value.\n";
  cout << "  The choices are 0/1 to select each element of W.\n";

  cout << "\n";
  cout << "  TARGET = " << target << "\n";
  cout << "\n";
  call_num = 0;
  n2 = -1;

  for ( ; ; )
  {
    backbin_rc ( n, reject, n2, choice );
    call_num = call_num + 1;

    if ( n2 == -1 )
    {
      break;
    }
//
//  Evaluate the partial sum.
//
    result = 0;
    for ( i = 0; i < n2; i++ )
    {
      result = result + choice[i] * w[i];
    }
//
//  If the sum is too big, then we reject it, and
//  all the related sums formed by making additional choices.
//
    reject = ( target < result );
//
//  If we hit the target, print out the information.
//
    if ( result == target && n2 == n )
    {
      cout << "  ";
      for ( i = 0; i < n; i++ )
      {
        cout << setw(2) << choice[i];
      }
      cout << "\n";
    }
  }

  cout << "\n";
  cout << "  Number of calls = " << call_num << "\n";
  cout << "  Binary search space = " << i4_power ( 2, n ) << "\n";

  return;
}
