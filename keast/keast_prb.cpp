# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "keast.hpp"

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
//    MAIN is the main program for KEAST_PRB.
//
//  Discussion:
//
//    KEAST_PRB calls the KEAST test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "KEAST_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the KEAST library.\n";

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
  cout << "KEAST_PRB:\n";
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
//    TEST01 tests KEAST_RULE_NUM, KEAST_DEGREE, KEAST_ORDER_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
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
  cout << "  KEAST_RULE_NUM returns the number of rules;\n";
  cout << "  KEAST_DEGREE returns the degree of a rule;\n";
  cout << "  KEAST_ORDER_NUM returns the order of a rule.\n";

  rule_num = keast_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Degree     Order\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = keast_order_num ( rule );
    degree = keast_degree ( rule );
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
//    TEST02 tests KEAST_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
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
  double *xyztab;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  KEAST_RULE returns the points and weights\n";
  cout << "  of a Keast rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that the weights\n";
  cout << "  sum to 1.\n";

  rule_num = keast_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";

  cout << "\n";
  cout << "      Rule      Order     Sum of weights\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = keast_order_num ( rule );

    xyztab = new double[3*order_num];
    wtab = new double[order_num];

    keast_rule ( rule, order_num, xyztab, wtab );

    wtab_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      wtab_sum = wtab_sum + wtab[order];
    }

    cout << "  " << setw(8) << rule
         << "  " << setw(8) << order_num
         << "  " << setw(14) << wtab_sum << "\n";

    delete [] wtab;
    delete [] xyztab;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests KEAST_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
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
  double *suborder_xyzz;
  double xyzz_sum;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  KEAST_RULE returns the points and weights\n";
  cout << "  of a Keast rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that, for each\n";
  cout << "  quadrature point, the barycentric coordinates\n";
  cout << "  add up to 1.\n";

  rule_num = keast_rule_num ( );

  cout << "\n";
  cout << "      Rule    Suborder    Sum of coordinates\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = keast_suborder_num ( rule );

    suborder_xyzz = new double[4*suborder_num];
    suborder_w = new double[suborder_num];

    keast_subrule ( rule, suborder_num, suborder_xyzz, suborder_w );

    cout << "\n";
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << suborder_num << "\n";

    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyzz_sum = suborder_xyzz[0+suborder*4]
               + suborder_xyzz[1+suborder*4]
               + suborder_xyzz[2+suborder*4]
               + suborder_xyzz[3+suborder*4];
     cout << "                    "
          << "  " << setprecision(16) << setw(25) << xyzz_sum << "\n";
    }

    delete [] suborder_w;
    delete [] suborder_xyzz;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrates TETRAHEDRON_REFERENCE_TO_PHYSICAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NODE_NUM 4

  int i;
  int node;
  double node_xyz[3*NODE_NUM] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  double node_xyz2[3*NODE_NUM] = {
    1.0, 2.0, 3.0,
    2.0, 2.0, 3.0,
    1.0, 3.0, 3.0,
    1.0, 2.0, 9.0 };
  int order;
  int order_num;
  int rule;
  double volume;
  double volume2;
  double *w;
  double *xyz;
  double *xyz2;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  TETRAHEDRON_REFERENCE_TO_PHYSICAL transforms a rule\n";
  cout << "  on the unit (reference) tetrahedron to a rule on \n";
  cout << "  an arbitrary (physical) tetrahedron.\n";

  rule = 2;

  order_num = keast_order_num ( rule );

  xyz = new double[3*order_num];
  xyz2 = new double[3*order_num];
  w = new double[order_num];

  keast_rule ( rule, order_num, xyz, w );
//
//  Here is the reference tetrahedron, and its rule.
//
  cout << "\n";
  cout << "  The reference tetrahedron:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xyz[0+node*3]
         << "  " << setw(14) << node_xyz[1+node*3]
         << "  " << setw(14) << node_xyz[2+node*3] << "\n";
  }

  volume = tetrahedron_volume ( node_xyz );

  cout << "\n";
  cout << "  Rule " << rule << " for reference tetrahedron\n";
  cout << "  with volume = " << volume << "\n";
  cout << "\n";
  cout << "                X               Y               Z             W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xyz[0+order*3]
         << "  " << setw(14) << xyz[1+order*3]
         << "  " << setw(14) << xyz[2+order*3]
         << "  " << setw(14) << w[order] << "\n";
  }
//
//  Transform the rule.
//
  tetrahedron_reference_to_physical ( node_xyz2, order_num, xyz, xyz2 );
//
//  Here is the physical tetrahedron, and its transformed rule.
//
  cout << "\n";
  cout << "  The physical tetrahedron:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xyz2[0+node*3]
         << "  " << setw(14) << node_xyz2[1+node*3]
         << "  " << setw(14) << node_xyz2[2+node*3] << "\n";
  }

  volume2 = tetrahedron_volume ( node_xyz2 );

  cout << "\n";
  cout << "  Rule " << rule << " for physical tetrahedron\n";
  cout << "  with volume = " << volume2 << "\n";
  cout << "\n";
  cout << "                X               Y               Z             W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xyz2[0+order*3]
         << "  " << setw(14) << xyz2[1+order*3]
         << "  " << setw(14) << xyz2[2+order*3]
         << "  " << setw(14) << w[order] << "\n";
  }

  delete [] w;
  delete [] xyz;
  delete [] xyz2;

  return;
# undef NODE_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests KEAST_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int degree_max = 3;
  int dim_num = 3;
  int expon[3];
  int h;
  double *mono;
  bool more;
  int order;
  int order_num;
  double quad;
  int rule;
  int rule_num;
  int t;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Demonstrate the KEAST rules on a sequence of\n";
  cout << "  monomial integrands X^A Y^B Z^C\n";
  cout << "  on the unit tetrahedron.\n";

  rule_num = keast_rule_num ( );

  cout << "\n";
  cout << "      Rule     Order     Quad\n";
  cout << "\n";

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      cout << "\n";
      cout << "  F(X,Y,Z)"
           << " = X^" << expon[0]
           << " * Y^" << expon[1]
           << " * Z^" << expon[2] << "\n";
      cout << "\n";

      for ( rule = 1; rule <= rule_num; rule++ )
      {
        order_num = keast_order_num ( rule );

        xyz = new double[3*order_num];
        w = new double[order_num];

        keast_rule ( rule, order_num, xyz, w );

        mono = monomial_value ( dim_num, order_num, xyz, expon );

        quad = r8vec_dot ( order_num, w, mono );

        cout << "  " << setw(8) << rule
             << "  " << setw(8) << order_num
             << "  " << setw(14) << quad << "\n";

        delete [] mono;
        delete [] w;
        delete [] xyz;
      }
      if ( !more )
      {
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests KEAST_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int order;
  int order_num;
  int rule;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  KEAST_RULE returns the points and weights\n";
  cout << "  of a Keast rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply print a rule.\n";

  rule = 10;
  degree = keast_degree ( rule );
  order_num = keast_order_num ( rule );

  cout << "\n";
  cout << "  Rule =   " << rule << "\n";
  cout << "  Degree = " << degree << "\n";
  cout << "  Order =  " << order << "\n";

  cout << "\n";
  cout << "         I      W               X               Y               Z\n";
  cout << "\n";

  xyz = new double[3*order_num];
  w = new double[order_num];

  keast_rule ( rule, order_num, xyz, w );

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8) << order
         << "  " << setw(14) << w[order]
         << "  " << setw(14) << xyz[0+order*3]
         << "  " << setw(14) << xyz[1+order*3]
         << "  " << setw(14) << xyz[2+order*3] << "\n";
  }

  delete [] w;
  delete [] xyz;

  return;
}
