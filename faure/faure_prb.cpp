# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "faure.hpp"

int main ( );
void test005 ( );
void test006 ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FAURE_PRB.
//
//  Discussion:
//
//    FAURE_PRB tests the FAURE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FAURE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FAURE library.\n";
 
  test005 ( );
  test006 ( );
  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FAURE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests BINOMIAL_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *coef;
  int i;
  int j;
  int m = 10;
  int n = 7;
  int qs = 7;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  BINOMIAL_TABLE computes a table of binomial.\n";
  cout << "  coefficients mod QS.\n";
  cout << "\n";
  cout << "  Here, QS = " << qs << "\n";

  coef = binomial_table ( qs, m, n );

  cout << "\n";
  cout << "   I/J";
  for ( j = 0; j <= n; j++ )
  {
    cout << setw(8) << j;
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i <= m; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j <= n; j++ )
    {
      cout << setw(8) << coef[i+j*(m+1)];
    }
    cout << "\n";
  }

  delete [] coef;

  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests I4_LOG_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  int j4;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  I4_LOG_R8: whole part of log base B,\n";
  cout << "\n";
  cout << "        I4        J4 I4_LOG_J4\n";
  cout << "\n";

  for ( j4 = 2; j4 <= 5; j4++ )
  {
    for ( i4 = 0; i4 <= 10; i4++ )
    {
      cout << "  " << setw(8) << i4
           << "  " << setw(8) << j4
           << "  " << setw(8) << i4_log_i4 ( i4, j4 ) << "\n";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests FAURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 4

  int dim;
  int dim_num;
  int i;
  int qs;
  double r[DIM_MAX];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  FAURE computes the next element of a Faure sequence.\n";
  cout << "\n";
  cout << "  In this test, we call FAURE repeatedly.\n";

  for ( dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
  {

    seed = -1;
    qs = prime_ge ( dim_num );

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "  The underlying base is QS = " << qs << "\n";
    cout << "\n";
    cout << "  Seed  Seed   Faure\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      seed_in = seed;
      faure ( dim_num, &seed, r );
      seed_out = seed;
      cout << setw(6) << seed_in << "  "
           << setw(6) << seed_out << "  ";
      for ( dim = 0; dim < dim_num; dim++ )
      {
        cout << setw(10) << r[dim] << "  ";
      }
      cout << "\n";
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
//    TEST02 tests FAURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int dim;
  int i;
  int qs;
  double r[DIM_NUM];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  FAURE computes the next element of a Faure sequence.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how the SEED can be\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";

  qs = prime_ge ( DIM_NUM );

  cout << "\n";
  cout << "  Using dimension DIM_NUM =   " << DIM_NUM << "\n";
  cout << "  The underlying base is QS = " << qs << "\n";

  cout << "\n";
  cout << "  Note that on the first call to FAURE, if\n";
  cout << "  SEED is negative, it is reset to a value that\n";
  cout << "  is the recommended starting point:\n";
  
  seed = -1;
  cout << "\n";
  cout << "  Seed  Seed   Faure\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    faure ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( dim = 0; dim < DIM_NUM; dim++ )
    {
      cout << setw(10) << r[dim] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  However, if the input value of SEED is 0,\n";
  cout << "  then no initial skipping is done.\n";

  seed = 0;

  cout << "\n";
  cout << "  Seed  Seed   Faure\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    faure ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( dim = 0; dim < DIM_NUM; dim++ )
    {
      cout << setw(10) << r[dim] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";
  cout << "\n";

  seed = 100;

  cout << "\n";
  cout << "  Seed  Seed   Faure\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    faure ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( dim = 0; dim < DIM_NUM; dim++ )
    {
      cout << setw(10) << r[dim] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump back by decreasing SEED:\n";
  cout << "\n";

  seed = 3;

  cout << "\n";
  cout << "  Seed  Seed   Faure\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    faure ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( dim = 0; dim < DIM_NUM; dim++ )
    {
      cout << setw(10) << r[dim] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests FAURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_base = 10;
  int dim_num;
  int i;
  int qs;
  double *r;
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  FAURE computes the next element of a Faure sequence.\n";
  cout << "\n";
  cout << "  In this test, we try some large dimensions.\n";

  for ( dim_num = dim_base; dim_num <= 6 * dim_base; dim_num = dim_num + dim_base )
  {
    r = new double[dim_num];

    seed = -1;
    qs = prime_ge ( dim_num );

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "  The underlying base is QS = " << qs << "\n";
    cout << "\n";
    cout << "  Seed  Seed   Faure\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 1; i <= 2; i++ )
    {
      seed_in = seed;
      faure ( dim_num, &seed, r );
      seed_out = seed;
      cout << "  " << setw(8) << seed_in
           << "  " << setw(8) << seed_out << "\n";
      cout << "                    ";
      for ( dim = 0; dim < dim_num; dim++ )
      {
        cout << "  " << setw(10) << r[dim];
        if ( ( dim + 1 ) % 5 == 0 || dim + 1 == dim_num )
        {
          cout << "\n";
        }
        if ( ( dim + 1 ) % 5 == 0 && dim + 1 < dim_num )
        {
          cout << "                    ";
        }
      }
    }
    delete [] r;
  }

  return;
}
