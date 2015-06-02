# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "unicycle.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UNICYCLE_PRB.
//
//  Discussion:
//
//    UNICYCLE_PRB tests the UNICYCLE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "UNICYCLE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the UNICYCLE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "UNICYCLE_PRB\n";
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
//    TEST01 tests PERM_IS_UNICYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int *p;
  int seed;
  int test;
  int test_num = 10;
  int *u;
  bool value;

  seed = 123456789;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  PERM_IS_UNICYCLE determines whether a permutation\n";
  cout << "  is a unicyle\n";

  for ( test = 1; test <= test_num; test++ )
  {
    p = perm_random ( n, seed );

    value = perm_is_unicycle ( n, p );

    if ( value )
    {
      perm_print ( n, p, "  This permutation is a unicycle" );
      u = unicycle_index_to_sequence ( n, p );
      unicycle_print ( n, u, "  The permutation in sequence form" );
      delete [] u;
    }
    else
    {
      perm_print ( n, p, "  This permutation is NOT a unicycle" );
    }
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests UNICYCLE_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int n_max = 10;
  int num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  UNICYCLE_ENUM enumerates the unicycles of N objects.\n";
  cout << "\n";
  cout << "  N    Number\n";
  cout << "\n";

  for ( n = 0; n <= n_max; n++ )
  {
    num = unicycle_enum ( n );
    cout << "  " << setw(3) << n
         << "  " << setw(8) << num << "\n";
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests UNICYCLE_INVERSE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 7;
  int u[7] = { 1, 7, 6, 2, 4, 3, 5 };
  int *u_inverse;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  UNICYCLE_INVERSE inverts a unicycle;\n";

  unicycle_print ( n, u, "  The original unicycle:" );
 
  u_inverse = unicycle_inverse ( n, u );
 
  unicycle_print ( n, u_inverse, "  The inverse unicycle:" );
 
  delete [] u_inverse;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests UNICYCLE_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int rank;
  int *u;

  u = new int[n];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  UNICYCLE_NEXT generates unicycles in lex order.\n";
  cout << "\n";
  rank = -1;
 
  for ( ; ; )
  {
    unicycle_next ( n, u, rank );

    if ( rank == -1 )
    {
      break;
    }

    cout << "  " << setw(3) << rank << ":";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << u[i];
    }
    cout << "\n";
  }

  delete [] u;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests UNICYCLE_RANDOM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int seed;
  int *u;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  UNICYCLE_RANDOM produces a random unicyle\n";;
  cout << "  For this test, N = " << n << "\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    u = unicycle_random ( n, seed );
    unicycle_print ( n, u, " " );
    delete [] u;
  }
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests UNICYCLE_RANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int rank;
  int u[5] = { 1, 5, 2, 3, 4 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  UNICYCLE_RANK ranks a unicycle.\n";

  unicycle_print ( n, u, "  The unicycle:" );
 
  rank = unicycle_rank ( n, u );
 
  cout << "\n";
  cout << "  The rank is: " << rank << "\n";
 
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests UNICYCLE_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int rank;
  int *u;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  UNICYCLE_UNRANK, given a rank, computes the\n";
  cout << "  corresponding unicycle.\n";
  cout << "\n";
  rank = 6;
  cout << "  The requested rank is " << rank << "\n";
 
  u = unicycle_unrank ( n, rank );
 
  unicycle_print ( n, u, "  The unicycle:" );
 
  delete [] u;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests UNICYCLE_INDEX, UNICYCLE_INDEX_TO_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n = 6;
  int *u;
  int *u_index;
  int *u2;
  int seed;
  int test;
  int test_num = 5;

  seed = 123456789;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  UNICYCLE_INDEX converts a unicycle to index form.\n";
  cout << "  UNICYCLE_INDEX_TO_SEQUENCE converts an index to unicycle form.\n";

  for ( test = 1; test <= test_num; test++ )
  {
    u = unicycle_random ( n, seed );

    unicycle_print ( n, u, "  The unicycle:" );

    u_index = unicycle_index ( n, u );
    
    unicycle_index_print ( n, u_index, "  The index form:" );

    u2 = unicycle_index_to_sequence ( n, u_index );

    unicycle_print ( n, u2, "  The unicycle recovered:" );

    delete [] u;
    delete [] u_index;
    delete [] u2;
  }
  return;
}
