# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "knapsack_01.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    KNAPSACK_01_TEST tests the KNAPSACK_01 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "KNAPSACK_01_TEST\n";
  cout << "  C++ version.\n";
  cout << "  Test the KNAPSACK_01 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "KNAPSACK_01_TEST\n";
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
//    TEST01 seeks a solution of the 0/1 Knapsack problem.
//
//  Discussion:
//
//    In the 0/1 knapsack problem, a knapsack of capacity C is given,
//    as well as N items, with the I-th item of weight W(I).
//
//    A selection is "acceptable" if the total weight is no greater than C.
//
//    It is desired to find an optimal acceptable selection, that is,
//    an acceptable selection such that there is no acceptable selection
//    of greater weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int i;
  int n = 6;
  int *s;
  int t;
  int w[6] = {
    16, 17, 23, 24, 39, 40 };

  c = 100;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Knapsack maximum capacity is " << c << "\n";
  cout << "  Come as close as possible to filling the knapsack.\n";

  s = knapsack_01 ( n, w, c );

  cout << "\n";
  cout << "   # 0/1  Weight\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << setw(4) << i << "  "
         << setw(1) << s[i] << "  "
         << setw(4) << w[i] << "\n";
  }
  t = 0;
  for ( i = 0; i < n; i++ )
  {
    t = t + s[i] * w[i];
  }
  cout << "\n";
  cout << "  Total:   " << t << "\n";

  delete [] s;

  return;
}
