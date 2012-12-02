# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "niederreiter.hpp"

int main ( );
void test01 ( int base );
void test02 ( int base );
void test03 ( int base, int dim );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NIEDERREITER_PRB.
//
//  Discussion:
//
//    NIEDERREITER_PRB calls a set of problems for NIEDERREITER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int base;
  int dim_num;

  timestamp ( );

  cout << "\n";
  cout << "NIEDERREITER_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NIEDERREITER routines.\n";

  base = 2;
  test01 ( base );

  base = 3;
  test01 ( base );

  base = 13;
  test01 ( base );

  base = 2;
  test02 ( base );

  base = 3;
  test02 ( base );

  base = 2;
  dim_num = 20;
  test03 ( base, dim_num );

  base = 2;
  dim_num = 29;
  test03 ( base, dim_num );
//
//  Terminate.
//
  cout << "\n";
  cout << "NIEDERREITER_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int base )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests NIEDERREITER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BASE, the base to use in the computation.
//    BASE should be a prime, or a power of a prime.
//
{
  const int dim_max = 4;

  int dim;
  int dim_num;
  int i;
  double r[dim_max];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  NIEDERREITER computes the next element of \n";
  cout << "  a Niederreiter quasirandom sequence using base BASE.\n";
  cout << "\n";
  cout << "  In this test, we call NIEDERREITER repeatedly.\n";
  cout << "\n";
  cout << "  Using base BASE =      " << base << "\n";

  for ( dim_num = 2; dim_num <= dim_max; dim_num++ )
  {
    seed = 0;

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "\n";
    cout << "    Seed    Seed     Niederreiter\n";
    cout << "      In     Out\n";
    cout << "\n";
    for ( i = 0; i <= 110; i++ )
    {
      seed_in = seed;
      niederreiter ( dim_num, base, &seed, r );
      seed_out = seed;
      if ( i <= 11 || 95 <= i )
      {
        cout << "  " << setw(8) << seed_in
	       << "  " << setw(8) << seed_out;
        for ( dim = 0; dim < dim_num; dim++ )
        {
          cout << setw(10) << r[dim];
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
}
//****************************************************************************80

void test02 ( int base )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests NIEDERREITER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BASE, the base to use in the computation.
//    BASE should be a prime, or a power of a prime.
//
{
  const int dim_num = 3;

  int dim;
  int i;
  double r[dim_num];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  NIEDERREITER computes the next element of\n";
  cout << "  a Niederreiter quasirandom sequence using base BASE.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how the SEED can be\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";

  cout << "\n";
  cout << "  Using base BASE =           " << base << "\n";
  cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";

  seed = 0;

  cout << "\n";
  cout << "    Seed    Seed     Niederreiter\n";
  cout << "      In     Out\n";
  cout << "\n";
  for ( i = 0; i <= 10; i++ )
  {
    seed_in = seed;
    niederreiter ( dim_num, base, &seed, r );
    seed_out = seed;
    cout << "  " << setw(8) << seed_in
         << "  " << setw(8) << seed_out;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(10) << r[dim];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";
  cout << "\n";

  seed = 100;

  cout << "\n";
  cout << "    Seed    Seed     Niederreiter\n";
  cout << "      In     Out\n";
  cout << "\n";
  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    niederreiter ( dim_num, base, &seed, r );
    seed_out = seed;
    cout << "  " << setw(8) << seed_in
         << "  " << setw(8) << seed_out;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(10) << r[dim];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump back by decreasing SEED:\n";
  cout << "\n";

  seed = 3;

  cout << "\n";
  cout << "    Seed    Seed     Niederreiter\n";
  cout << "      In     Out\n";
  cout << "\n";
  for ( i = 0; i <= 10; i++ )
  {
    seed_in = seed;
    niederreiter ( dim_num, base, &seed, r );
    seed_out = seed;
    cout << "  " << setw(8) << seed_in
         << "  " << setw(8) << seed_out;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(10) << r[dim];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";
  cout << "\n";

  seed = 98;

  cout << "\n";
  cout << "    Seed    Seed     Niederreiter\n";
  cout << "      In     Out\n";
  cout << "\n";
  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    niederreiter ( dim_num, base, &seed, r );
    seed_out = seed;
    cout << "  " << setw(8) << seed_in
         << "  " << setw(8) << seed_out;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(10) << r[dim];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( int base, int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests NIEDERREITER.
//
//  Discussion:
//
//    Simply verify that a few terms of a sequence of given dimension
//    can be computed.  Most recently, the NIEDERREITER code was set
//    up to handle up to dimension 50...we think.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BASE, the base to use in the computation.
//    BASE should be a prime, or a power of a prime.
//
//    Input, int DIM, the spatial dimension.
//
{
  int dim;
  int i;
  double *r;
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  NIEDERREITER computes the next element of\n";
  cout << "  a Niederreiter quasirandom sequence using base BASE.\n";
  cout << "\n";
  cout << "  In this test, we simply generate ten elements in a given base\n";
  cout << "  and dimension.\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";

  cout << "\n";
  cout << "  Using base BASE =           " << base << "\n";
  cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";

  seed = 0;
  r = new double[dim_num];

  cout << "\n";
  cout << "    Seed    Seed     Niederreiter\n";
  cout << "      In     Out\n";
  cout << "\n";
  for ( i = 0; i <= 10; i++ )
  {
    seed_in = seed;
    niederreiter ( dim_num, base, &seed, r );
    seed_out = seed;
    cout << "  " << setw(8) << seed_in
         << "  " << setw(8) << seed_out;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(10) << r[dim];
      if ( ( ( dim + 1 ) % 5 == 0 ) && ( dim != dim_num ) )
      {
        cout << "\n";
        cout << "                    ";
      }
    }
    cout << "\n";
  }
  delete [] r;

  return;
}
