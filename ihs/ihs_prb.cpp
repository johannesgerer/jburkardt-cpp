# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "ihs.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    IHS_PRB tests the improved hypercube sampling algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "IHS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the IHS library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "IHS_PRB\n";
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
//    TEST01 tests the improved distributed hypercube sampling algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  double average;
  double covc;
  int dim_num;
  int duplication = 5;
  int i;
  int j;
  double opt;
  int point_num = 10;
  int seed;
  double std;
  int *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  IHS implements the IHS Algorithm\n";
  cout << "  (Improved Distributed Hypercube Sampling)\n";
  cout << "\n";
  cout << "  Demonstrate the code for a fixed number of points\n";
  cout << "  and an increasing dimension.\n";

  for ( dim_num = 1; dim_num <= 4; dim_num++ )
  {
    x = new int [ dim_num * point_num ];

    seed = 17;

    opt = ( ( double ) point_num ) /
      pow ( ( ( double ) point_num ),
      ( double ) ( 1.0E+00 / ( ( double ) dim_num ) ) );

    cout << "\n";
    cout << "  Random number seed =       " << seed << "\n";
    cout << "  Spatial dimension =        " << dim_num << "\n";
    cout << "  Number of points =         " << point_num << "\n";
    cout << "  Duplication factor =       " << duplication << "\n";
    cout << "  Desired minimum distance = " << opt << "\n";
//
//  Get the points.
//
    ihs ( dim_num, point_num, duplication, &seed, x );
//
//  Compute the covariance.
//
    covariance ( dim_num, point_num, x, &average, &std, &covc );

    cout << "\n";
    cout << "  Average minimum distance " << average << "\n";
    cout << "  Standard deviation:      " << std << "\n";
    cout << "  Covariance:              " << covc << "\n";

    cout << "\n";

    for ( j = 0; j < point_num; j++ )
    {
      cout << setw(4) << j+1 << "    ";
      for ( i = 0; i < dim_num; i++ )
      {
        cout << setw(4) << x[i+j*dim_num] << "  ";
      }
      cout << "\n";
    }

    delete [] x;

  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests the improved distributed hypercube sampling algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2003
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define POINT_NUM 10

  double average;
  double covc;
  int duplication;
  int i;
  int j;
  double opt;
  int seed;
  double std;
  int x[ DIM_NUM * POINT_NUM ];
//
  cout << "\n";
  cout << "TEST02\n";
  cout << "  IHS implements the IHS Algorithm\n";
  cout << "  (Improved Distributed Hypercube Sampling)\n";
  cout << "\n";
  cout << "  Demonstrate the code for a fixed number of points\n";
  cout << "  and dimension, but vary the duplication value.\n";

  opt = ( ( double ) POINT_NUM ) /
    pow ( ( ( double ) POINT_NUM ),
    ( double ) ( 1.0E+00 / ( ( double ) DIM_NUM ) ) );

  cout << "\n";
  cout << "  Spatial dimension =        " << DIM_NUM << "\n";
  cout << "  Number of points =         " << POINT_NUM << "\n";
  cout << "  Desired minimum distance = " << opt << "\n";

  for ( duplication = 1; duplication <= 5; duplication++ )
  {
    seed = 17;

    cout << "\n";
    cout << "  Random number seed =       " << seed << "\n";
    cout << "  Duplication factor =       " << duplication << "\n";
//
//  Get the points.
//
    ihs ( DIM_NUM, POINT_NUM, duplication, &seed, x );
//
//  Compute the covariance.
//
    covariance ( DIM_NUM, POINT_NUM, x, &average, &std, &covc );

    cout << "\n";
    cout << "  Average minimum distance " << average << "\n";
    cout << "  Standard deviation:      " << std << "\n";
    cout << "  Covariance:              " << covc << "\n";

    cout << "\n";

    for ( j = 0; j < POINT_NUM; j++ )
    {
      cout << setw(4) << j+1 << "    ";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << setw(4) << x[i+j*DIM_NUM] << "  ";
      }
      cout << "\n";
    }

  }

  return;

# undef DIM_NUM
# undef POINT_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests the improved distributed hypercube sampling algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  double average;
  double covc;
  int dim_num = 2;
  int duplication = 5;
  int i;
  int j;
  int k;
  double opt;
  int point_num;
  int seed;
  double std;
  int *x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  IHS implements the IHS Algorithm\n";
  cout << "  (Improved Distributed Hypercube Sampling)\n";
  cout << "\n";
  cout << "  Demonstrate the code for a fixed dimension\n";
  cout << "  and duplication value, and increasing number of points.\n";

  cout << "\n";
  cout << "  Spatial dimension =        " << dim_num << "\n";
  cout << "  Duplication factor =       " << duplication << "\n";

  point_num = 5;

  for ( k = 1; k <= 5; k++ )
  {
    point_num = 2 * point_num;

    x = new int [ dim_num * point_num ];

    opt = ( ( double ) point_num ) /
      pow ( ( ( double ) point_num ),
      ( double ) ( 1.0E+00 / ( ( double ) dim_num ) ) );

    seed = 17;

    cout << "\n";
    cout << "  Random number seed =       " << seed << "\n";
    cout << "  Number of points =         " << point_num << "\n";
    cout << "  Desired minimum distance = " << opt << "\n";
//
//  Get the points.
//
    ihs ( dim_num, point_num, duplication, &seed, x );
//
//  Compute the covariance.
//
    covariance ( dim_num, point_num, x, &average, &std, &covc );

    cout << "\n";
    cout << "  Average minimum distance " << average << "\n";
    cout << "  Standard deviation:      " << std << "\n";
    cout << "  Covariance:              " << covc << "\n";

    cout << "\n";

    for ( j = 0; j < point_num; j++ )
    {
      if ( j <= 10 || point_num-10 <= j )
      {
        cout << setw(4) << j+1 << "    ";
        for ( i = 0; i < dim_num; i++ )
        {
          cout << setw(4) << x[i+j*dim_num] << "  ";
        }
        cout << "\n";
      }
      else if ( j == 11 )
      {
        cout << "....    ........\n";
      }
    }
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests the improved distributed hypercube sampling algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2003
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define POINT_NUM 10

  double average;
  double covc;
  int duplication = 5;
  int i;
  int j;
  int k;
  double opt;
  int seed;
  double std;
  int x[ DIM_NUM * POINT_NUM ];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  IHS implements the IHS Algorithm\n";
  cout << "  (Improved Distributed Hypercube Sampling)\n";
  cout << "\n";
  cout << "  Demonstrate the code for a fixed number of points,\n";
  cout << "  dimension, and duplication factor, but with a\n";
  cout << "  varying random number seed.\n";

  opt = ( ( double ) POINT_NUM ) /
    pow ( ( ( double ) POINT_NUM ),
    ( double ) ( 1.0 / ( ( double ) DIM_NUM ) ) );

  cout << "\n";
  cout << "  Spatial dimension =        " << DIM_NUM << "\n";
  cout << "  Number of points =         " << POINT_NUM << "\n";
  cout << "  Duplication factor =       " << duplication << "\n";
  cout << "  Desired minimum distance = " << opt << "\n";

  seed = 17;

  for ( k = 1; k <= 4; k++ )
  {
    cout << "\n";
    cout << "  Random number seed =       " << seed << "\n";
//
//  Get the points.
//
    ihs ( DIM_NUM, POINT_NUM, duplication, &seed, x );
//
//  Compute the covariance.
//
    covariance ( DIM_NUM, POINT_NUM, x, &average, &std, &covc );

    cout << "\n";
    cout << "  Average minimum distance " << average << "\n";
    cout << "  Standard deviation:      " << std << "\n";
    cout << "  Covariance:              " << covc << "\n";

    cout << "\n";

    for ( j = 0; j < POINT_NUM; j++ )
    {
      cout << setw(4) << j+1 << "    ";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << setw(4) << x[i+j*DIM_NUM] << "  ";
      }
      cout << "\n";
    }

  }

  return;
# undef DIM_NUM
# undef POINT_NUM
}
