# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "wandzura.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//******************************************************************************

int main ( )

//******************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for WANDZURA_PRB.
//
//  Discussion:
//
//    WANDZURA_PRB calls the WANDZURA test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "WANDZURA_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the WANDZURA library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WANDZURA_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//******************************************************************************

void test01 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST01 tests WANDZURA_RULE_NUM, WANDZURA_DEGREE, WANDZURA_ORDER_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int order_num;
  int rule;
  int rule_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  WANDZURA_RULE_NUM returns the number of rules;\n";
  cout << "  WANDZURA_DEGREE returns the degree of a rule;\n";
  cout << "  WANDZURA_ORDER_NUM returns the order of a rule.\n";

  rule_num = wandzura_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Degree     Order\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = wandzura_order_num ( rule );
    degree = wandzura_degree ( rule );
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << degree
         << "  " << setw(8) << order_num << "\n";
  }

  return;
}
//******************************************************************************

void test02 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST02 tests WANDZURA_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int order;
  int order_num;
  int rule;
  int rule_num;
  double *w;
  double w_sum;
  double *xy;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  WANDZURA_RULE returns the points and weights\n";
  cout << "  of a Wandzura rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that the weights\n";
  cout << "  sum to 1.\n";

  rule_num = wandzura_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Sum of weights\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = wandzura_order_num ( rule );

    xy = new double[2*order_num];
    w = new double[order_num];

    wandzura_rule ( rule, order_num, xy, w );

    w_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      w_sum = w_sum + w[order];
    }

    cout << "  " << setw(8) << rule
         << "  " << setw(14) << w_sum << "\n";

    delete [] w;
    delete [] xy;
  }
  return;
}
//******************************************************************************

void test03 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST03 tests WANDZURA_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int rule;
  int rule_num;
  int suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
  double xyz_sum;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  WANDZURA_RULE returns the points and weights\n";
  cout << "  of a Wandzura rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that, for each\n";
  cout << "  quadrature point, the barycentric coordinates\n";
  cout << "  add up to 1.\n";

  rule_num = wandzura_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Suborder    Sum of coordinates\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = wandzura_suborder_num ( rule );

    suborder_xyz = new double[3*suborder_num];
    suborder_w = new double[suborder_num];

    wandzura_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

    cout << "\n";
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << suborder_num << "\n";

    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyz_sum = suborder_xyz[0+suborder*3]
              + suborder_xyz[1+suborder*3]
              + suborder_xyz[2+suborder*3];
     cout << "                    "
          << "  " << setprecision(16) << setw(25) << xyz_sum << "\n";
    }

    delete [] suborder_w;
    delete [] suborder_xyz;
  }
  return;
}
//******************************************************************************

void test04 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST04 tests WANDZURA_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double area = 0.5;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  int order_num;
  double quad;
  int rule;
  int rule_num;
  double value;
  double *w;
  double x;
  double *xy;
  double y;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  WANDZURA_RULE returns the points and weights of\n";
  cout << "  a Wandzura rule for the unit triangle.\n";
  cout << "\n";
  cout << "  This routine uses those rules to estimate the\n";
  cout << "  integral of monomomials in the unit triangle.\n";

  rule_num = wandzura_rule_num ( );

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
//
//  Multiplying X**A * Y**B by COEF will give us an integrand
//  whose integral is exactly 1.  This makes the error calculations easy.
//
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef * ( double ) ( a + i ) / ( double ) ( i );
      }

      cout << "\n";
      cout << "  Integrate " << coef
           << " * X^" << a
           << " * Y^" << b << "\n";
      cout << "\n";
      cout << "      Rule       QUAD           ERROR\n";
      cout << "\n";

      for ( rule = 1; rule <= rule_num; rule++ )
      {
        order_num = wandzura_order_num ( rule );

        xy = new double[2*order_num];
        w = new double[order_num];

        wandzura_rule ( rule, order_num, xy, w );

        quad = 0.0;

        for ( order = 0; order < order_num; order++ )
        {
          x = xy[0+order*2];
          y = xy[1+order*2];

          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( y, b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( x, a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( x, a ) * pow ( y, b );
          }
          quad = quad + w[order] * value;
        }
        quad = area * quad;

        exact = 1.0;
        err = fabs ( exact - quad );

        cout << "  " << setw(8)  << rule
             << "  " << setw(14) << quad
             << "  " << setw(14) << err << "\n";

        delete [] w;
        delete [] xy;
      }
    }
  }
  return;
}
//******************************************************************************

void test05 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST05 plots the Wandzura points in the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NODE_NUM 3

  char file_name[80];
  int i;
  int node_show = 1;
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0 };
  int order_num;
  int point_show = 2;
  int rule;
  int rule_num;
  double *w;
  double *xy;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  This routine creates an EPS plot of each\n";
  cout << "  set of Wandzura points.\n";
  cout << "\n";

  rule_num = wandzura_rule_num ( );

  strcpy ( file_name, "wandzura_rule_0.eps" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    file_name_inc ( file_name );

    order_num = wandzura_order_num ( rule );

    xy = new double[2*order_num];
    w = new double[order_num];

    wandzura_rule ( rule, order_num, xy, w );

    triangle_points_plot ( file_name, node_xy, node_show, order_num,
      xy, point_show );

    delete [] w;
    delete [] xy;

    cout << "  Rule " << rule << " plotted in \"" << file_name << "\".\n";
  }

  return;
# undef NODE_NUM
}
//******************************************************************************

void test06 ( )

//******************************************************************************
//
//  Purpose:
//
//    TEST06 demonstrates REFERENCE_TO_PHYSICAL_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NODE_NUM 3

  double area;
  double area2;
  int i;
  int node;
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0 };
  double node_xy2[2*NODE_NUM] = {
    1.0, 2.0,
    1.0, 1.0,
    3.0, 2.0 };
  int order;
  int order_num;
  int point_show = 2;
  int rule;
  double *w;
  double *xy;
  double *xy2;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  REFERENCE_TO_PHYSICAL_T3 transforms a rule\n";
  cout << "  on the unit (reference) triangle to a rule on \n";
  cout << "  an arbitrary (physical) triangle.\n";

  rule = 2;

  order_num = wandzura_order_num ( rule );

  xy = new double[2*order_num];
  xy2 = new double[2*order_num];
  w = new double[order_num];

  wandzura_rule ( rule, order_num, xy, w );
//
//  Here is the reference triangle, and its rule.
//
  cout << "\n";
  cout << "  The reference triangle:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xy[0+node*2]
         << "  " << setw(14) << node_xy[1+node*2] << "\n";
  }

  area = triangle_area ( node_xy );

  cout << "\n";
  cout << "  Rule " << rule << " for reference triangle\n";
  cout << "  with area = " << area << "\n";
  cout << "\n";
  cout << "                X               Y               W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xy[0+order*2]
         << "  " << setw(14) << xy[1+order*2]
         << "  " << setw(14) << w[order] << "\n";
  }
//
//  Transform the rule.
//
  reference_to_physical_t3 ( node_xy2, order_num, xy, xy2 );
//
//  Here is the physical triangle, and its transformed rule.
//
  cout << "\n";
  cout << "  The physical triangle:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xy2[0+node*2]
         << "  " << setw(14) << node_xy2[1+node*2] << "\n";
  }

  area2 = triangle_area ( node_xy2 );

  cout << "\n";
  cout << "  Rule " << rule << " for physical triangle\n";
  cout << "  with area = " << area2 << "\n";
  cout << "\n";
  cout << "                X               Y               W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xy2[0+order*2]
         << "  " << setw(14) << xy2[1+order*2]
         << "  " << setw(14) << w[order] << "\n";
  }

  delete [] w;
  delete [] xy;
  delete [] xy2;

  return;
# undef NODE_NUM
}
