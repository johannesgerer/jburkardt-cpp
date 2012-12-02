# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "ncc_tetrahedron.hpp"

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
//    MAIN is the main program for NCC_TETRAHEDRON_PRB.
//
//  Discussion:
//
//    NCC_TETRAHEDRON_PRB calls the NCC_TETRAHEDRON test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "NCC_TETRAHEDRON_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the NCC_TETRAHEDRON library.\n";

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
  cout << "NCC_TETRAHEDRON_PRB:\n";
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
//    TEST01 tests NCC_TETRAHEDRON_RULE_NUM, NCC_TETRAHEDRON_DEGREE, NCC_TETRAHEDRON_ORDER_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
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
  cout << "  NCC_TETRAHEDRON_RULE_NUM returns the number of rules;\n";
  cout << "  NCC_TETRAHEDRON_DEGREE returns the degree of a rule;\n";
  cout << "  NCC_TETRAHEDRON_ORDER_NUM returns the order of a rule.\n";

  rule_num = ncc_tetrahedron_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Degree     Order\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = ncc_tetrahedron_order_num ( rule );
    degree = ncc_tetrahedron_degree ( rule );
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
//    TEST02 tests NCC_TETRAHEDRON_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num = 3;
  int order;
  int order_num;
  int rule;
  int rule_num;
  double *wtab;
  double wtab_sum;
  double *xyztab;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  NCC_TETRAHEDRON_RULE returns the points and weights\n";
  cout << "  of an NCC rule for the tetrahedron.\n";
  cout << "\n";
  cout << "  In this test, we simply check that the weights\n";
  cout << "  sum to 1.\n";

  rule_num = ncc_tetrahedron_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";

  cout << "\n";
  cout << "      Rule      Sum of weights\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = ncc_tetrahedron_order_num ( rule );

    xyztab = new double[dim_num*order_num];
    wtab = new double[order_num];

    ncc_tetrahedron_rule ( rule, order_num, xyztab, wtab );

    wtab_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      wtab_sum = wtab_sum + wtab[order];
    }

    cout << "  " << setw(8) << rule
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
//    TEST03 tests NCC_TETRAHEDRON_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
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
  cout << "  NCC_TETRAHEDRON_RULE returns the points and weights\n";
  cout << "  of an NCC rule for the tetrahedron.\n";
  cout << "\n";
  cout << "  In this test, we simply check that, for each\n";
  cout << "  quadrature point, the barycentric coordinates\n";
  cout << "  add up to 1.\n";

  rule_num = ncc_tetrahedron_rule_num ( );

  cout << "\n";
  cout << "      Rule    Suborder    Sum of coordinates\n";
  cout << "\n";

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = ncc_tetrahedron_suborder_num ( rule );

    suborder_xyz = new double[4*suborder_num];
    suborder_w = new double[suborder_num];

    ncc_tetrahedron_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

    cout << "\n";
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << suborder_num << "\n";

    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyz_sum = suborder_xyz[0+suborder*4]
              + suborder_xyz[1+suborder*4]
              + suborder_xyz[2+suborder*4]
              + suborder_xyz[3+suborder*4];
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
//    TEST04 tests NCC_TETRAHEDRON_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  double coef;
  int dim_num = 3;
  double err;
  double exact;
  int i;
  int order;
  int order_num;
  double quad;
  int rule;
  int rule_num;
  double value;
  double volume;
  double *wtab;
  double x;
  double *xyztab;
  double y;
  double z;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  NCC_TETRAHEDRON_RULE returns the points and weights of\n";
  cout << "  an NCC rule for the unit tetrahedron.\n";
  cout << "\n";
  cout << "  This routine uses those rules to estimate the\n";
  cout << "  integral of monomomials in the unit tetrahedron.\n";

  rule_num = ncc_tetrahedron_rule_num ( );

  volume = 1.0 / 6.0;

  for ( a = 0; a <= 6; a++ )
  {
    for ( b = 0; b <= 6 - a; b++ )
    {
      for ( c = 0; c <= 6 - a - b; c++ )
      {
//
//  Multiplying X**A * Y**B * Z**C by COEF will give us an integrand
//  whose integral is exactly 1.  This makes the error calculations easy.
//
        coef = 1.0;

//      for ( i = 1; i <= a; i++ )
//      {
//        coef = coef * i / i;
//      }
        for ( i = a + 1; i <= a + b; i++ )
        {
          coef = coef * ( double ) ( i ) / ( double ) ( i - a );
        }
        for ( i = a + b + 1; i <= a + b + c; i++ )
        {
          coef = coef * ( double ) ( i ) / ( double ) ( i - a - b );
        }
        for ( i = a + b + c + 1; i <= a + b + c + 3; i++ )
        {
          coef = coef * ( double ) ( i );
        }

        cout << "\n";
        cout << "  Integrate " << coef
             << " * X^" << a
             << " * Y^" << b
             << " * Z^" << c << "\n";
        cout << "\n";
        cout << "      Rule       QUAD           ERROR\n";
        cout << "\n";

        for ( rule = 1; rule <= rule_num; rule++ )
        {
          order_num = ncc_tetrahedron_order_num ( rule );

          xyztab = new double[dim_num*order_num];
          wtab = new double[order_num];

          ncc_tetrahedron_rule ( rule, order_num, xyztab, wtab );

          quad = 0.0;

          for ( order = 0; order < order_num; order++ )
          {
            x = xyztab[0+order*dim_num];
            y = xyztab[1+order*dim_num];
            z = xyztab[2+order*dim_num];
//
//  Some tedious calculations to avoid 0**0 complaints.
//
            value = coef;

            if ( a != 0 )
            {
              value = value * pow ( x, a );
            }

            if ( b != 0 )
            {
              value = value * pow ( y, b );
            }

            if ( c != 0 )
            {
              value = value * pow ( z, c );
            }

            quad = quad + wtab[order] * value;
          }
          quad = volume * quad;

          exact = 1.0;
          err = fabs ( exact - quad );

          cout << "  " << setw(8)  << rule
               << "  " << setw(14) << quad
               << "  " << setw(14) << err << "\n";

          delete [] wtab;
          delete [] xyztab;
        }
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
//    TEST05 demonstrates REFERENCE_TO_PHYSICAL_T4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define NODE_NUM 4

  int i;
  int node;
  double node_xyz[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  double node_xyz2[DIM_NUM*NODE_NUM] = {
    4.0, 5.0, 1.0,
    6.0, 5.0, 1.0,
    4.0, 8.0, 1.0,
    4.0, 5.0, 5.0 };
  int order;
  int order_num;
  int rule;
  double volume;
  double volume2;
  double *w;
  double *xyz;
  double *xyz2;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  REFERENCE_TO_PHYSICAL_T4 transforms a rule\n";
  cout << "  on the unit (reference) tetrahedron to a rule on \n";
  cout << "  an arbitrary (physical) tetrahedron.\n";

  rule = 3;

  order_num = ncc_tetrahedron_order_num ( rule );

  xyz = new double[DIM_NUM*order_num];
  xyz2 = new double[DIM_NUM*order_num];
  w = new double[order_num];

  ncc_tetrahedron_rule ( rule, order_num, xyz, w );
//
//  Here is the reference tetrahedron, and its rule.
//
  cout << "\n";
  cout << "  The reference tetrahedron:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xyz[0+node*DIM_NUM]
         << "  " << setw(14) << node_xyz[1+node*DIM_NUM]
         << "  " << setw(14) << node_xyz[2+node*DIM_NUM] << "\n";
  }

  volume = tetrahedron_volume ( node_xyz );

  cout << "\n";
  cout << "  Rule " << rule << " for reference tetrahedron\n";
  cout << "  with volume = " << volume << "\n";
  cout << "\n";
  cout << "                X               Y               Z               W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xyz[0+order*DIM_NUM]
         << "  " << setw(14) << xyz[1+order*DIM_NUM]
         << "  " << setw(14) << xyz[2+order*DIM_NUM]
         << "  " << setw(14) << w[order] << "\n";
  }
//
//  Transform the rule.
//
  reference_to_physical_t4 ( node_xyz2, order_num, xyz, xyz2 );
//
//  Here is the physical tetrahedron, and its transformed rule.
//
  cout << "\n";
  cout << "  The physical tetrahedron:\n";
  cout << "\n";

  for ( node = 0; node < NODE_NUM; node++ )
  {
    cout << "  " << setw(8)  << node+1
         << "  " << setw(14) << node_xyz2[0+node*DIM_NUM]
         << "  " << setw(14) << node_xyz2[1+node*DIM_NUM]
         << "  " << setw(14) << node_xyz2[2+node*DIM_NUM] << "\n";
  }

  volume2 = tetrahedron_volume ( node_xyz2 );

  cout << "\n";
  cout << "  Rule " << rule << " for physical tetrahedron\n";
  cout << "  with volume = " << volume2 << "\n";
  cout << "\n";
  cout << "                X               Y               Z               W\n";
  cout << "\n";

  for ( order = 0; order < order_num; order++ )
  {
    cout << "  " << setw(8)  << order
         << "  " << setw(14) << xyz2[0+order*DIM_NUM]
         << "  " << setw(14) << xyz2[1+order*DIM_NUM]
         << "  " << setw(14) << xyz2[2+order*DIM_NUM]
         << "  " << setw(14) << w[order] << "\n";
  }

  delete [] w;
  delete [] xyz;
  delete [] xyz2;

  return;
# undef DIM_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests NCC_TETRAHEDRON_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num = 3;
  int o;
  int order_num;
  int rule;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  NCC_TETRAHEDRON_RULE returns the points and weights\n";
  cout << "  of an NCC rule for the tetrahedron.\n";

  rule = 4;

  cout << "\n";
  cout << "  In this test, we simply print rule " << rule << "\n";

  suborder_num = ncc_tetrahedron_suborder_num ( rule );

  suborder = ncc_tetrahedron_suborder ( rule, suborder_num );

  suborder_w = new double[suborder_num];
  suborder_xyz = new double[4*suborder_num];

  ncc_tetrahedron_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

  cout << "\n";
  cout << "  The compressed rule:\n";
  cout << "\n";
  cout << "  Number of suborders = " << suborder_num << "\n";
  cout << "\n";
  cout << "     S   Sub     Weight     Xsi1      Xsi2      Xsi3      Xsi4\n";
  cout << "\n";

  for ( s = 0; s < suborder_num; s++ )
  {
    cout << "  " << setw(4) << s+1
         << "  " << setw(4) << suborder[s]
         << "  " << setw(8) << suborder_w[s]
         << "  " << setw(8) << suborder_xyz[0+s*4]
         << "  " << setw(8) << suborder_xyz[1+s*4]
         << "  " << setw(8) << suborder_xyz[2+s*4]
         << "  " << setw(8) << setprecision(4) << suborder_xyz[3+s*4] << "\n";
  }

  order_num = ncc_tetrahedron_order_num ( rule );

  xyz = new double[dim_num*order_num];
  w = new double[order_num];

  ncc_tetrahedron_rule ( rule, order_num, xyz, w );

  cout << "\n";
  cout << "  The full rule:\n";
  cout << "\n";
  cout << "  Order = " << order_num << "\n";
  cout << "\n";
  cout << "     O    Weight        X           Y           Z\n";
  cout << "\n";

  for ( o = o; o < order_num; o++ )
  {
    cout << "  " << setw(4) << o+1
         << "  " << setw(8) << w[o]
         << "  " << setw(8) << xyz[0+o*3]
         << "  " << setw(8) << xyz[1+o*3]
         << "  " << setw(8) << setprecision(4) << xyz[2+o*3] << "\n";
  }

  delete [] suborder;
  delete [] suborder_w;
  delete [] suborder_xyz;
  delete [] w;
  delete [] xyz;

  return;
}

