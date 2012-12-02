# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_interp_2d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN tests the TEST_INTERP_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_INTERP_2D_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INTERP_2D library.\n";
  cout << "  The test requires access to the R8LIB library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INTERP_2D_PRB\n";
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
//    TEST01 simply prints the title of each grid and function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  int f_num;
  int fi;
  string ft;
  int g_num;
  int gi;
  string gt;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For each grid and function, print the title.\n";

  g_num = g00_num ( );

  cout << "\n";
  cout << "  GRIDS:\n";
  cout << "  Index  Title\n";
  cout << "\n";

  for ( gi = 1; gi <= g_num; gi++ )
  {
    gt = g00_title ( gi );
    cout << "  " << setw(2) << gi
         << "  " << gt << "\n";
  }

  f_num = f00_num ( );

  cout << "\n";
  cout << "  FUNCTIONS:\n";
  cout << "  Index  Title\n";
  cout << "\n";

  for ( fi = 1; fi <= f_num; fi++ )
  {
    ft = f00_title ( fi );
    cout << "  " << setw(2) << fi
         << "  " << ft << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 samples each function using each grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *f;
  double f_ave;
  double f_max;
  double f_min;
  int f_num;
  int fi;
  string ft;
  int g_num;
  int gi;
  int gn;
  string gt;
  double *gx;
  double *gy;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Sample each function over each grid.\n";

  g_num = g00_num ( );
  f_num = f00_num ( );

  for ( fi = 1; fi <= f_num; fi++ )
  {
    ft = f00_title ( fi );
    cout << "\n";
    cout << "  " << setw(2) << fi
         << "  " << ft << "\n";
    cout << "        Grid Title                     ";
    cout << "Min(F)          Ave(F)           Max(F)\n";
    cout << "\n";

    for ( gi = 1; gi <= g_num; gi++ )
    {
      gt = g00_title ( gi );
      gn = g00_size ( gi );

      gx = new double[gn];
      gy = new double[gn];

      g00_xy ( gi, gn, gx, gy );

      f = new double[gn];

      f00_f0 ( fi, gn, gx, gy, f );

      f_max = r8vec_max ( gn, f );
      f_min = r8vec_min ( gn, f );
      f_ave = r8vec_sum ( gn, f );
      f_ave = f_ave / ( double ) ( gn );

      cout << "  " << setw(4) << gi
           << "  " << setw(25) << gt
           << "  " << setw(14) << f_min
           << "  " << setw(14) << f_ave
           << "  " << setw(14) << f_max << "\n";

      delete [] f;
      delete [] gx;
      delete [] gy;
    }
  }
  return;
}
