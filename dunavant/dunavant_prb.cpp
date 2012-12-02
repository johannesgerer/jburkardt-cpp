# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "dunavant.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DUNAVANT_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Discussion:
//
//    DUNAVANT_PRB calls the DUNAVANT test routines.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "DUNAVANT_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the DUNAVANT library.\n";

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
  cout << "DUNAVANT_PRB:\n";
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
//    TEST01 tests DUNAVANT_RULE_NUM, DUNAVANT_DEGREE, DUNAVANT_ORDER_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
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
  cout << "  DUNAVANT_RULE_NUM returns the number of rules;\n";
  cout << "  DUNAVANT_DEGREE returns the degree of a rule;\n";
  cout << "  DUNAVANT_ORDER_NUM returns the order of a rule.\n";

  rule_num = dunavant_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Degree     Order\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = dunavant_order_num ( rule );
    degree = dunavant_degree ( rule );
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << degree
         << "  " << setw(8) << order_num << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DUNAVANT_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
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
  double *wtab;
  double wtab_sum;
  double *xytab;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  DUNAVANT_RULE returns the points and weights\n";
  cout << "  of a Dunavant rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that the weights\n";
  cout << "  sum to 1.\n";

  rule_num = dunavant_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";

  cout << "\n";
  cout << "      Rule      Order     Sum of weights\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = dunavant_order_num ( rule );

    xytab = new double[2*order_num];
    wtab = new double[order_num];

    dunavant_rule ( rule, order_num, xytab, wtab );

    wtab_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      wtab_sum = wtab_sum + wtab[order];
    }

    cout << "  " << setw(8) << rule
         << "  " << setw(8) << order_num
         << "  " << setw(14) << wtab_sum << "\n";

    delete [] wtab;
    delete [] xytab;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DUNAVANT_RULE.
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
  int rule;
  int rule_num;
  int suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
  double xyz_sum;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  DUNAVANT_RULE returns the points and weights\n";
  cout << "  of a Dunavant rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that, for each\n";
  cout << "  quadrature point, the barycentric coordinates\n";
  cout << "  add up to 1.\n";

  rule_num = dunavant_rule_num ( );

  cout << "\n";
  cout << "      Rule    Suborder    Sum of coordinates\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = dunavant_suborder_num ( rule );

    suborder_xyz = new double[3*suborder_num];
    suborder_w = new double[suborder_num];

    dunavant_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

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
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DUNAVANT_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
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
  int rule_max = 10;
  double value;
  double *wtab;
  double x;
  double *xytab;
  double y;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  DUNAVANT_RULE returns the points and weights of\n";
  cout << "  a Dunavant rule for the unit triangle.\n";
  cout << "\n";
  cout << "  This routine uses those rules to estimate the\n";
  cout << "  integral of monomomials in the unit triangle.\n";

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

      for ( rule = 1; rule <= rule_max; rule++ )
      {
        order_num = dunavant_order_num ( rule );

        xytab = new double[2*order_num];
        wtab = new double[order_num];

        dunavant_rule ( rule, order_num, xytab, wtab );

        quad = 0.0;

        for ( order = 0; order < order_num; order++ )
        {
          x = xytab[0+order*2];
          y = xytab[1+order*2];

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
          quad = quad + wtab[order] * value;
        }
        quad = area * quad;

        exact = 1.0;
        err = fabs ( exact - quad );

        cout << "  " << setw(8)  << rule
             << "  " << setw(14) << quad
             << "  " << setw(14) << err << "\n";

        delete [] wtab;
        delete [] xytab;
      }
    }
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 plots the Dunavant points in the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
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
  int rule_max = 10;
  double *wtab;
  double *xytab;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  This routine creates an EPS plot of each\n";
  cout << "  set of Dunavant points.\n";
  cout << "\n";

  strcpy ( file_name, "dunavant_rule_00.eps" );

  for ( rule = 1; rule <= rule_max; rule++ )
  {
    file_name_inc ( file_name );

    order_num = dunavant_order_num ( rule );

    xytab = new double[2*order_num];
    wtab = new double[order_num];

    dunavant_rule ( rule, order_num, xytab, wtab );

    triangle_points_plot ( file_name, node_xy, node_show, order_num,
      xytab, point_show );

    delete [] wtab;
    delete [] xytab;

    cout << "  Rule " << rule << " plotted in \"" << file_name << "\".\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
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

  order_num = dunavant_order_num ( rule );

  xy = new double[2*order_num];
  xy2 = new double[2*order_num];
  w = new double[order_num];

  dunavant_rule ( rule, order_num, xy, w );
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
