# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "pyramid_felippa_rule.hpp"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PYRAMID_FELIPPA_RULE_PRB.
//
//  Discussion:
//
//    PYRAMID_FELIPPA_RULE_PRB tests the PYRAMID_FELIPPA_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "PYRAMID_FELIPPA_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PYRAMID_FELIPPA_RULE library.\n";

  degree_max = 4;
  pyramid_unit_monomial_test ( degree_max );

  degree_max = 5;
  pyramid_unit_quad_test ( degree_max );
//
//  Terminate.
//
  cout << "\n";
  cout << "PYRAMID_FELIPPA_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void pyramid_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_MONOMIAL_TEST tests PYRAMID__UNIT_MONOMIAL.
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
  int expon[3];
  int gamma;
  double value;

  cout << "\n";
  cout << "PYRAMID__UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit pyramid,\n";
  cout << "  PYRAMID__UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << pyramid_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = pyramid_unit_monomial ( expon );

        cout << "  " << setw(8)  << expon[0]
             << "  " << setw(8)  << expon[1]
             << "  " << setw(8)  << expon[2]
             << "  " << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void pyramid_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_QUAD_TEST tests the rules for the unit pyramid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2008
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
# define DIM_NUM 3

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
  double *xyz;

  cout << "\n";
  cout << "PYRAMID__UNIT_QUAD_TEST\n";
  cout << "  For the unit pyramid,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  PYRAMID__UNIT_O01,\n";
  cout << "  PYRAMID__UNIT_O05,\n";
  cout << "  PYRAMID__UNIT_O06,\n";
  cout << "  PYRAMID__UNIT_O08,\n";
  cout << "  PYRAMID__UNIT_O08b,\n";
  cout << "  PYRAMID__UNIT_O09,\n";
  cout << "  PYRAMID__UNIT_O13,\n";
  cout << "  PYRAMID__UNIT_O18,\n";
  cout << "  PYRAMID__UNIT_O27,\n";
  cout << "  PYRAMID__UNIT_O48,\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[0] % 2 ) == 1 || ( expon[1] % 2 ) == 1 )
    {
      continue;
    }

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
    xyz = new double[dim_num*order];
    pyramid_unit_o01 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 5;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o05 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 6;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o06 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o08 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o08b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 9;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o09 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 13;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o13 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 18;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o18 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 27;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o27 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 48;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyramid_unit_o48 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    cout << "\n";
    quad = pyramid_unit_monomial ( expon );
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
