# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "triangle_felippa_rule.hpp"

int main ( );
void triangle_unit_monomial_test ( int degree_max );
void triangle_unit_quad_test ( int degree_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_FELIPPA_RULE_PRB.
//
//  Discussion:
//
//    TRIANGLE_FELIPPA_RULE_PRB tests the TRIANGLE_FELIPPA_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_FELIPPA_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_FELIPPA_RULE library.\n";

  degree_max = 4;
  triangle_unit_monomial_test ( degree_max );

  degree_max = 7;
  triangle_unit_quad_test ( degree_max );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_FELIPPA_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void triangle_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_MONOMIAL_TEST tests TRIANGLE_UNIT_MONOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[2];
  double value;

  cout << "\n";
  cout << "TRIANGLE_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit triangle,\n";
  cout << "  TRIANGLE_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA\n";
  cout << "\n";
  cout << "  Volume = " << triangle_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;

      value = triangle_unit_monomial ( expon );

      cout << "  " << setw(8)  << expon[0]
           << "  " << setw(8)  << expon[1]
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void triangle_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_QUAD_TEST tests the rules for the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  bool more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *xy;

  cout << "\n";
  cout << "TRIANGLE_UNIT_QUAD_TEST\n";
  cout << "  For the unit triangle,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  TRIANGLE_UNIT_O01,\n";
  cout << "  TRIANGLE_UNIT_O03,\n";
  cout << "  TRIANGLE_UNIT_O03b,\n";
  cout << "  TRIANGLE_UNIT_O06,\n";
  cout << "  TRIANGLE_UNIT_O06b,\n";
  cout << "  TRIANGLE_UNIT_O07,\n";
  cout << "  TRIANGLE_UNIT_O12,\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    order = 1;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o01 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 3;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o03 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 3;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o03b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 6;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o06 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 6;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o06b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 7;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o07 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 12;
    w = new double[order];
    xy = new double[dim_num*order];
    triangle_unit_o12 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    cout << "\n";
    quad = triangle_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}