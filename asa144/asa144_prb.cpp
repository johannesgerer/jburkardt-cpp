# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa144.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA144_PRB.
//
//  Discussion:
//
//    ASA144_PRB tests the ASA144 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA144_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA144 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA144_PRB\n";
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
//    TEST01 tests RCONT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2008
//
//  Author:
//
//    John Burkardt
//
{
# define NROW 5
# define NCOL 5

  int i;
  int ifault;
  bool key;
  int matrix[NROW*NCOL];
  int ncolt[NCOL] = { 2, 2, 2, 2, 1 };
  int nrowt[NROW] = { 3, 2, 2, 1, 1 };
  int nsubt[NCOL];
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  RCONT constructs a random matrix with\n";
  cout << "  given row and column sums.\n";

  i4vec_print ( NROW, nrowt, "  The rowsum vector:" );
  i4vec_print ( NCOL, ncolt, "  The columnsum vector: " );

  key = false;

  for ( test = 1; test <= test_num; test++ )
  {
    rcont ( NROW, NCOL, nrowt, ncolt, nsubt, matrix, &key, &ifault );

    if ( ifault != 0 )
    {
      cout << "\n";
      cout << "  RCONT returned IFAULT = " << ifault << "\n";
      return;
    }
    i4mat_print ( NROW, NCOL, matrix, "  The rowcolsum matrix:" );
  }

  return;
# undef NROW
# undef NCOL
}
