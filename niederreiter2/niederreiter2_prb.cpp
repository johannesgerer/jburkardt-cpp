# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "niederreiter2.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NIEDERREITER2_PRB.
//
//  Discussion:
//
//    NIEDERREITER2_PRB tests the NIEDERREITER2 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "NIEDERREITER2_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NIEDERREITER2 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "NIEDERREITER2_PRB\n";
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
//    TEST01 tests NIEDERREITER2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 4

  int dim_num;
  int i;
  int j;
  double r[DIM_MAX];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  NIEDERREITER2 computes the next element of \n";
  cout << "  a Niederreiter quasirandom sequence using base 2.\n";
  cout << "\n";
  cout << "  In this test, we call NIEDERREITER2 repeatedly.\n";

  for ( dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
  {
    seed = 0;

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "\n";
    cout << "  Seed  Seed   Niederreiter2\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 0; i <= 110; i++ )
    {
      seed_in = seed;
      niederreiter2 ( dim_num, &seed, r );
      seed_out = seed;
      if ( i <= 11 || 95 <= i )
      {
        cout << setw ( 6 ) << seed_in << "  ";
        cout << setw ( 6 ) << seed_out << " ";
        for ( j = 0; j < dim_num; j++ )
        {
          cout << setw(10) << r[j] << "  ";
        }
        cout << "\n";
      }
      else if ( i == 12 )
      {
        cout << "......................\n";
      }

    }

  }

  return;
# undef DIM_MAX
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests NIEDERREITER2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  int j;
  double r[DIM_NUM];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  NIEDERREITER2 computes the next element of \n";
  cout << "  a Niederreiter quasirandom sequence using base 2.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how the SEED can be\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";

  cout << "\n";
  cout << "  Using dimension DIM_NUM =   " << DIM_NUM << "\n";

  seed = 0;

  cout << "\n";
  cout << "  Seed  Seed   Niederreiter2\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    seed_in = seed;
    niederreiter2 ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw ( 6 ) << seed_in << "  ";
    cout << setw ( 6 ) << seed_out << " ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(10) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 100;

  cout << "\n";
  cout << "  Seed  Seed   Niederreiter2\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 0; i <= 5; i++ )
  {
    seed_in = seed;
    niederreiter2 ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw ( 6 ) << seed_in << "  ";
    cout << setw ( 6 ) << seed_out << " ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(10) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump back by decreasing SEED:\n";

  seed = 3;

  cout << "\n";
  cout << "  Seed  Seed   Niederreiter2\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    seed_in = seed;
    niederreiter2 ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw ( 6 ) << seed_in << "  ";
    cout << setw ( 6 ) << seed_out << " ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(10) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 98;

  cout << "\n";
  cout << "  Seed  Seed   Niederreiter2\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    niederreiter2 ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw ( 6 ) << seed_in << "  ";
    cout << setw ( 6 ) << seed_out << " ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(10) << r[j] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
