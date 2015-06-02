# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "asa314.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA314_PRB.
//
//  Discussion:
//
//    ASA314_PRB tests the ASA314 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA314_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA314 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA314_PRB:\n";
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
//    TEST01 tests INVMOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
{
  int cmod[3];
  int i;
  int ifault;
  int imat[3*3];
  int jmat[3*3] = { 1, 0, 0, 2, 1, 0, 1, 0, 1 };
  int mat[3*3] = { 1, 0, 0, 1, 1, 0, 2, 0, 1 };
  int nrow = 3;
  int rmod[3];

  for ( i = 0; i < nrow; i++ )
  {
    cmod[i] = 3;
  }
  for ( i = 0; i < nrow; i++ )
  {
    rmod[i] = 3;
  }

  cout << "\n";
  cout << "TEST01\n";
  cout << "  INVMOD computes the inverse of a matrix\n";
  cout << "  whose elements are subject to modulo arithmetic.\n";

  i4mat_print ( nrow, nrow, mat, "  The matrix to be inverted:" );

  invmod ( mat, imat, rmod, cmod, nrow, ifault );

  i4mat_print ( nrow, nrow, imat, "  The computed inverse:" );

  i4mat_print ( nrow, nrow, jmat, "  The correct inverse:" );

  return;
}

