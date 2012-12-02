# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "lyness_rule.hpp"

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
//    MAIN is the main program for LYNESS_RULE_PRB.
//
//  Discussion:
//
//    LYNESS_RULE_PRB calls the LYNESS_RULE test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "LYNESS_RULE_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LYNESS_RULE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
// Terminate.
//
  cout << "\n";
  cout << "LYNESS_RULE_PRB:\n";
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
//    TEST01 tests LYNESS_RULE_NUM, LYNESS_DEGREE, LYNESS_ORDER_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int order;
  int precision;
  int rule;
  int rule_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  LYNESS_RULE_NUM returns the number of rules;\n";
  cout << "  LYNESS_DEGREE returns the degree of a rule;\n";
  cout << "  LYNESS_ORDER_NUM returns the order of a rule.\n";

  rule_num = lyness_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule     Order  Precision\n";
  cout << "\n";

  for ( rule = 0; rule <= rule_num; rule++ )
  {
    order = lyness_order ( rule );
    precision = lyness_precision ( rule );
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << order
         << "  " << setw(8) << precision << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 performs the weight sum test on Lyness rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int order;
  int rule;
  int rule_num;
  double *w;
  double w_sum;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  LYNESS_RULE returns the points and weights\n";
  cout << "  of a Lyness rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that the weights\n";
  cout << "  sum to 1.\n";

  rule_num = lyness_rule_num ( );

  cout << "\n";
  cout << "  Number of available rules = " << rule_num << "\n";
  cout << "\n";
  cout << "      Rule    Sum of weights\n";
  cout << "\n";

  for ( rule = 0; rule <= rule_num; rule++ )
  {
    order = lyness_order ( rule );

    x = new double[2*order];
    w = new double[order];

    lyness_rule ( rule, order, w, x );

    w_sum = r8vec_sum ( order, w );

    cout << "  " << setw(8) << rule
         << "  " << setw(25) << w_sum << "\n";

    delete [] w;
    delete [] x; 
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 performs the barycentric coordinate sum test on Lyness rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2010
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
  double *sub_w;
  double *sub_xyz;
  double xyz_sum;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  LYNESS_RULE returns the points and weights\n";
  cout << "  of a Lyness rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply check that, for each\n";
  cout << "  quadrature point, the barycentric coordinates\n";
  cout << "  sum to 1.\n";

  rule_num = lyness_rule_num ( );

  cout << "\n";
  cout << "      Rule   Suborder    Sum of coordinates\n";

  for ( rule = 0; rule < rule_num; rule++ )
  {
    suborder_num = lyness_suborder_num ( rule );

    sub_w = new double[suborder_num];
    sub_xyz = new double[3*suborder_num];

    lyness_subrule ( rule, suborder_num, sub_xyz, sub_w );

    cout << "\n";
    cout << "  " << setw(8) << rule
         << "  " << setw(8) << suborder_num << "\n";
    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyz_sum = r8vec_sum ( 3, sub_xyz+3*suborder );
      cout << "                           " << xyz_sum << "\n";
    }

    delete [] sub_w;
    delete [] sub_xyz; 
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 prints a rule generated by LYNESS_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int order;
  int precision;
  int rule;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  LYNESS_RULE returns the points and weights\n";
  cout << "  of a Lyness rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply print a rule.\n";

  rule = 18;
  order = lyness_order ( rule );
  precision = lyness_precision ( rule );

  cout << "\n";
  cout << "  Rule =      " << rule << "\n";
  cout << "  Order =     " << order << "\n";
  cout << "  Precision = " << precision << "\n";

  w = new double[order];
  x = new double[2*order];

  lyness_rule ( rule, order, w, x );

  cout << "\n";
  cout << "     I      W                   X                   Y\n";
  cout << "\n";

  for ( j = 0; j < order; j++ )
  {
    cout << "  " << setw(4) << j
         << "  " << setprecision(16) << setw(24) << w[j]
         << "  " << setprecision(16) << setw(24) << x[0+j*2]
         << "  " << setprecision(16) << setw(24) << x[1+j*2] << "\n";
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
//    TEST05 writes a rule created by LYNESS_RULE to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int order;
  int precision;
  double r[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 1.0 };
  int rule;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  LYNESS_RULE returns the points and weights\n";
  cout << "  of a Lyness rule for the triangle.\n";
  cout << "\n";
  cout << "  In this test, we simply print a rule.\n";

  rule = 18;
  order = lyness_order ( rule );
  precision = lyness_precision ( rule );

  cout << "\n";
  cout << "  Rule =      " << rule << "\n";
  cout << "  Order =     " << order << "\n";
  cout << "  Precision = " << precision << "\n";

  w = new double[order];
  x = new double[2*order];

  lyness_rule ( rule, order, w, x );

  r8mat_write ( "lyness_18_r.txt", 2, 3, r );
  cout << "\n";
  cout << "  Wrote the region file \"lyness_18_r.txt\".\n";
  r8mat_write ( "lyness_18_w.txt", 1, order, w );
  cout << "  Wrote the weight file \"lyness_18_w.txt\".\n";
  r8mat_write ( "lyness_18_x.txt", 2, order, x );
  cout << "  Wrote the point file \"lyness_18_x.txt\".\n";

  delete [] w;
  delete [] x;
    
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests the Lyness rules for exact integration of monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double area;
  int b;
  double coef;
  int degree;
  int degree_max = 10;
  double err;
  double exact;
  int i;
  int j;
  int order;
  double quad;
  int rule;
  int rule_num;
  double value;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  LYNESS_RULE returns the points and weights of\n";
  cout << "  a Lyness rule for the unit triangle.\n";
  cout << "\n";
  cout << "  This routine uses those rules to estimate the\n";
  cout << "  integral of monomomials in the unit triangle.\n";

  rule_num = lyness_rule_num ( );

  area = 0.5;

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    for ( a = 0; a <= degree; a++ )
    {
      b = degree - a;
//
//  Multiplying X^A * Y^B by COEF will give us an integrand
//  whose integral is exactly 1.  This makes the error calculations easy.
//
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef * ( double ) ( a + i ) / ( double ) ( i );
      }

      cout << "\n";
      cout << "  Integrate " << coef << " * X^" << a << " * Y^" << b << "\n";
      cout << "\n";
      cout << "      Rule       QUAD           ERROR\n";
      cout << "\n";

      for ( rule = 0; rule <= rule_num; rule++ )
      {
        order = lyness_order ( rule );

        w = new double[order];
        x = new double[2*order];

        lyness_rule ( rule, order, w, x );

        quad = 0.0;

        for ( j = 0; j < order; j++ )
        {
          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( x[1+j*2], b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( x[0+j*2], a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( x[0+j*2], a ) * pow ( x[1+j*2], b );
          }

          quad = quad + w[j] * value;
        }

        quad = area * quad;

        exact = 1.0;
        err = r8_abs ( exact - quad );

        cout << "  " << setw(8) << rule
             << "  " << setprecision(6) << setw(14) << quad
             << "  " << setprecision(2) << setw(10) << err << "\n";

        delete [] w;
        delete [] x;
      }
    }
  }
  return;
}
