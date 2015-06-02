# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>

using namespace std;

# include "simplex_gm_rule.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SIMPLEX_GM_RULE_PRB.
//
//  Discussion:
//
//    SIMPLEX_GM_RULE_PRB tests the SIMPLEX_GM_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SIMPLEX_GM_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SIMPLEX_GM_RULE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SIMPLEX_GM_RULE_PRB\n";
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
//    TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int j;
  double *phy;
  double *phy_unit;
  int point_num = 10;
  double *ref;
  int seed = 123456789;
  double t[DIM_NUM*(DIM_NUM + 1 )] = {
    1.0, 1.0,
    3.0, 1.0,
    2.0, 5.0 };
  double t_unit[DIM_NUM*(DIM_NUM + 1 )] = {
    0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0 };
  int vertex_num = DIM_NUM + 1;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SIMPLEX_UNIT_TO_GENERAL\n";
  cout << "  maps points in the unit simplex to a general simplex.\n";
  cout << "\n";
  cout << "  Here we consider a simplex in 2D, a triangle.\n";
  cout << "\n";
  cout << "  The vertices of the general triangle are:\n";
  cout << "\n";
  for ( j = 0; j < vertex_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << t[dim+j*dim_num];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "   (  XSI     ETA )   ( X       Y  )\n";
  cout << "\n";

  phy_unit = simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit );

  for ( j = 0; j < dim_num + 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << t_unit[dim+j*dim_num];
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << phy_unit[dim+j*dim_num];
    }
    cout << "\n";
  }
  ref = simplex_unit_sample ( dim_num, point_num, &seed );

  phy = simplex_unit_to_general ( dim_num, point_num, t, ref );

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << ref[dim+j*dim_num];
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << phy[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] phy;
  delete [] phy_unit;
  delete [] ref;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SIMPLEX_UNIT_TO_GENERAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int dim;
  int dim_num = DIM_NUM;
  int vertex_num = DIM_NUM + 1;
  int j;
  double *phy;
  double *phy_unit;
  int point_num = 10;
  double *ref;
  int seed = 123456789;
  double t[DIM_NUM*(DIM_NUM + 1 )] = {
    1.0, 1.0, 1.0,
    3.0, 1.0, 1.0,
    1.0, 4.0, 1.0,
    1.0, 1.0, 5.0 };
  double t_unit[DIM_NUM*(DIM_NUM + 1 )] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  SIMPLEX_UNIT_TO_GENERAL\n";
  cout << "  maps points in the unit simplex to a general simplex.\n";
  cout << "\n";
  cout << "  Here we consider a simplex in 3D, a tetrahedron.\n";
  cout << "\n";
  cout << "  The vertices of the general tetrahedron are:\n";
  cout << "\n";
  for ( j = 0; j < vertex_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << t[dim+j*dim_num];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "   (  XSI     ETA     MU )    ( X       Y       Z )\n";
  cout << "\n";

  phy_unit = simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit );

  for ( j = 0; j < dim_num + 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << t_unit[dim+j*dim_num];
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << phy_unit[dim+j*dim_num];
    }
    cout << "\n";
  }

  ref = simplex_unit_sample ( dim_num, point_num, &seed );

  phy = simplex_unit_to_general ( dim_num, point_num, t, ref );

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << ref[dim+j*dim_num];
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(9) << phy[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] phy;
  delete [] phy_unit;
  delete [] ref;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests GM_RULE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  int dim_num;
  int dim_num_test[TEST_NUM] = { 2, 3, 5, 10 };
  int degree;
  int point_num;
  int rule;
  int test;
  int test_num = TEST_NUM;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  GM_RULE_SIZE returns POINT_NUM, the number of points\n";
  cout << "  associated with a Grundmann-Moeller quadrature rule\n";
  cout << "  for the unit simplex of dimension DIM_NUM\n";
  cout << "  with rule index RULE\n";
  cout << "  and degree of exactness DEGREE = 2*RULE+1.\n";

  cout << "\n";
  cout << "   DIM_NUM      RULE    DEGREE POINT_NUM\n";

  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_test[test];

    cout << "\n";

    for ( rule = 0; rule <= 5; rule++ )
    {
      point_num = gm_rule_size ( rule, dim_num );
      degree = 2 * rule + 1;

      cout << "  " << setw(8) << dim_num
           << "  " << setw(8) << rule
           << "  " << setw(8) << degree
           << "  " << setw(8) << point_num << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests GM_RULE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num;
  int point;
  int point_num;
  int rule;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  GM_RULE_SET determines the weights and abscissas\n";
  cout << "  of a Grundmann-Moeller quadrature rule for\n";
  cout << "  the DIM_NUM dimensional simplex,\n";
  cout << "  using a rule of in index RULE,\n";
  cout << "  which will have degree of exactness 2*RULE+1.\n";

  dim_num = 3;
  rule = 2;

  cout << "\n";
  cout << "  Here we use DIM_NUM = " << dim_num << "\n";
  cout << "  RULE = " << rule << "\n";
  cout << "  DEGREE = " << 2 * rule + 1 << "\n";

  point_num = gm_rule_size ( rule, dim_num );

  w = new double[point_num];
  x = new double[dim_num*point_num];

  gm_rule_set ( rule, dim_num, point_num, w, x );

  cout << "\n";
  cout << "     POINT        W             X             Y             Z\n";
  cout << "\n";

  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(8) << point + 1
         << "  " << setw(12) << w[point];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(12) << x[dim+point*dim_num];
    }
    cout << "\n";
  }

  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests GM_RULE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  int dim_num;
  int dim_num_test[TEST_NUM] = { 2, 3, 5, 10 };
  int point;
  int point_num;
  int rule;
  int test;
  int test_num = TEST_NUM;
  double *w;
  double w_sum;
  double *x;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  GM_RULE_SET determines the weights and abscissas\n";
  cout << "  of a Grundmann-Moeller quadrature rule for\n";
  cout << "  the DIM_NUM dimensional simplex,\n";
  cout << "  using a rule of in index RULE,\n";
  cout << "  which will have degree of exactness 2*RULE+1.\n";
  cout << "\n";
  cout << "  In this test, we compute various rules, and simply\n";
  cout << "  report the number of points, and the sum of weights.\n";

  cout << "\n";
  cout << "   DIM_NUM      RULE    POINT_NUM  WEIGHT SUM\n";

  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_test[test];

    cout << "\n";

    for ( rule = 0; rule <= 5; rule++ )
    {
      point_num = gm_rule_size ( rule, dim_num );

      w = new double[point_num];
      x = new double[dim_num*point_num];

      gm_rule_set ( rule, dim_num, point_num, w, x );

      w_sum = 0.0;
      for ( point = 0; point < point_num; point++ )
      {
        w_sum = w_sum + w[point];
      }

      cout << "  " << setw(8) << dim_num
           << "  " << setw(8) << rule
           << "  " << setw(8) << point_num
           << "  " << setprecision(16) << setw(24) << w_sum << "\n";

      delete [] w;
      delete [] x;
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests GM_RULE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int dim;
  int dim_num;
  int point;
  int point_num;
  int rule;
  double *w;
  char w_file[127];
  ofstream w_unit;
  double *x;
  char x_file[127];
  ofstream x_unit;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  GM_RULE_SET determines the weights and abscissas\n";
  cout << "  of a Grundmann-Moeller quadrature rule for\n";
  cout << "  the DIM_NUM dimensional simplex,\n";
  cout << "  using a rule of in index RULE,\n";
  cout << "  which will have degree of exactness 2*RULE+1.\n";
  cout << "\n";
  cout << "  In this test, we write a rule to a file.\n";

  dim_num = 3;
  rule = 2;

  cout << "\n";
  cout << "  Here we use DIM_NUM = " << dim_num << "\n";
  cout << "  RULE = " << rule << "\n";
  cout << "  DEGREE = " << 2 * rule + 1 << "\n";

  point_num = gm_rule_size ( rule, dim_num );

  w = new double[point_num];
  x = new double[dim_num*point_num];

  gm_rule_set ( rule, dim_num, point_num, w, x );

  sprintf ( w_file, "gm%d_%dd_w.txt", rule, dim_num );

  w_unit.open ( w_file );

  for ( point = 0; point < point_num; point++ )
  {
    w_unit << setprecision ( 16 ) << setw(20) << w[point] << "\n";
  }

  w_unit.close ( );

  sprintf ( x_file, "gm%d_%dd_x.txt", rule, dim_num );

  x_unit.open ( x_file );

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x_unit << setprecision ( 16 ) << setw(20) << x[dim+point*dim_num];
    }
    x_unit << "\n";
  }

  x_unit.close ( );

  cout << "\n";
  cout << "  Wrote rule " << rule
       << " to \"" << w_file
       << "\" and \"" << x_file << "\n";

  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests GM_RULE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int degree_max = 4;
  int dim_num = 5;
  int *expon;
  int h;
  int t;
  double *mono;
  bool more;
  int point;
  int point_num;
  double quad;
  double quad_error;
  int rule;
  int rule_max = 3;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  GM_RULE_SET determines the weights and abscissas\n";
  cout << "  of a Grundmann-Moeller quadrature rule for\n";
  cout << "  the DIM_NUM dimensional simplex,\n";
  cout << "  using a rule of in index RULE,\n";
  cout << "  which will have degree of exactness 2*RULE+1.\n";
  cout << "\n";
  cout << "  In this test, look at all the monomials up to\n";
  cout << "  some maximum degree, choose a few low order rules\n";
  cout << "  and determine the quadrature error for each.\n";
  cout << "\n";
  cout << "  Here we use DIM_NUM = " << dim_num << "\n";

  cout << "\n";
  cout << "      Rule     Order     Quad_Error\n";
  cout << "\n";

  expon = new int[dim_num];

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = false;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      cout << "\n";
      cout << "  F(X) = X1^" << expon[0]
           << " * X2^" << expon[1]
           << " * X3^" << expon[2]
           << " * X4^" << expon[3]
           << " * X5^" << expon[4] << "\n";
      cout << "\n";

      for ( rule = 0; rule <= rule_max; rule++ )
      {
        point_num = gm_rule_size ( rule, dim_num );

        mono = new double[point_num];
        w = new double[point_num];
        x = new double[dim_num*point_num];

        gm_rule_set ( rule, dim_num, point_num, w, x );

        quad_error = simplex_unit_monomial_quadrature ( dim_num, expon,
          point_num, x, w );

        cout << "  " << setw(8) << rule
             << "  " << setw(8) << point_num
             << "  " << setw(14) << quad_error << "\n";

        delete [] mono;
        delete [] w;
        delete [] x;
      }

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] expon;

  return;
}
