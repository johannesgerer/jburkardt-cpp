# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "grid.hpp"

int main ( );
void test01 ( int center, int *seed );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GRID_PRB.
//
//  Discussion:
//
//    GRID_PRB tests the grid routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int center;
  int seed;

  timestamp ( );
  cout << "\n";
  cout << "GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the GRID library.\n";

  center = 1;
  seed = 123456789;

  test01 ( center, &seed );

  cout << "\n";
  cout << "  Repeat TEST01 with a different seed from the first run.\n";

  seed = 987654321;
  test01 ( center, &seed );

  cout << "\n";
  cout << "  Repeat TEST01 with the same seed as the first run.\n";

  seed = 123456789;
  test01 ( center, &seed );

  cout << "\n";
  cout << "  Repeat TEST01 with different centering values.\n";

  for ( center = 1; center <= 5; center ++ )
  {
    seed = 123456789;
    test01 ( center, &seed );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "GRID_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int center, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define POINT_NUM 10

  int i;
  int j;
  double x[DIM_NUM*POINT_NUM];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  GRID_GENERATE randomly chooses a given number of\n";
  cout << "  points on a uniform grid.\n";
  cout << "\n";
  cout << "  Spatial dimension =  " << DIM_NUM   << "\n";
  cout << "  Number of points =   " << POINT_NUM << "\n";
  cout << "  Random number SEED = " << *seed     << "\n";
  cout << "  Centering option =   " << center    << "\n";

  grid_generate ( DIM_NUM, POINT_NUM, center, seed, x );

  cout << "\n";
  cout << "  The grid points:\n";
  cout << "\n";

  for ( j = 0; j < POINT_NUM; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << x[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef POINT_NUM
}
