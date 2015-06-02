# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "wedge_felippa_rule.hpp"

int main ( );
void test01 ( int degree_max );
void test02 ( int degree_max );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WEDGE_FELIPPA_RULE_PRB.
//
//  Discussion:
//
//    FELIPPA_PRB tests the FELIPPA library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "WEDGE_FELIPPA_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WEDGE_FELIPPA_RULE library.\n";

  degree_max = 4;

  test01 ( degree_max );
  test02 ( degree_max );

  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WEDGE_FELIPPA_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests WEDGE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2014
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
  cout << "TEST01\n";
  cout << "  For the unit wedge,\n";
  cout << "  WEDGE_INTEGRAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << wedge_volume ( ) << "\n";
  cout << "\n";
  cout << "   ALPHA      BETA     GAMMA      INTEGRAL\n";
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
        value = wedge_integral ( expon );
        cout << setw(8) << expon[0] << "  "
             << setw(8) << expon[1] << "  "
             << setw(8) << expon[2] << "  "
             << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void test02 ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests the rules for the unit wedge.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2014
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
  int dim_num = 3;
  int expon[3];
  int h;
  int line_order;
  int line_order_array[7] = {
    1, 2, 2, 3, 2, 3, 4 };
  bool more;
  int order;
  double quad;
  int t;
  int test;
  int test_num = 7;
  int triangle_order;
  int triangle_order_index;
  int triangle_order_array[7] = {
    1, 3, -3, 6, -6, 7, 12 };
  double *v;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For the unit wedge,\n";
  cout << "  we approximate monomial integrals with WEDG_UNIT_RULE.\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[2] % 2 ) == 1 )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    cout << "\n";
    cout << "  Monomial exponents:   "
         << setw(2) << expon[0] << "  "
         << setw(2) << expon[1] << "  "
         << setw(2) << expon[2] << "\n";
    cout << "\n";

    for ( test = 0; test < test_num; test++ )
    {
      line_order = line_order_array[test];
      triangle_order = triangle_order_array[test];

      order = line_order * fabs ( triangle_order );

      w = new double[order];
      xyz = new double[dim_num * order];
      wedge_rule ( line_order, triangle_order, w, xyz );
      v = monomial_value ( dim_num, order, expon, xyz );
      quad = wedge_volume ( ) * r8vec_dot_product ( order, w, v );
      cout << setw(6) << triangle_order << "  "
           << setw(6) << line_order << "  "
           << setw(14) << quad << "\n";
      delete [] v;
      delete [] w;
      delete [] xyz;
    }

    cout << "\n";
    quad = wedge_integral ( expon );
    cout << "   Exact        "
         << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 writes out some rules for the unit wedge.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  int dim_num = 3;
  int line_order;
  int line_order_array[7] = {
    1, 2, 2, 3, 2, 3, 4 };
  int order;
  int rule;
  int rule_num = 7;
  int triangle_order;
  int triangle_order_array[7] = {
    1, 3, -3, 6, -6, 7, 12 };
  double *w;
  string w_filename;
  double *x;
  string x_filename;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For the unit wedge,\n";
  cout << "  write some rules to a file\n";
  cout << "\n";
  cout << "   Rule  Trig    Line   Total  W_File X_File\n";
  cout << "         Order   Order  Order\n";
  cout << "\n";

  for ( rule = 0; rule < rule_num; rule++ )
  {
    if ( rule == 0 )
    {
      w_filename = "wedge_felippa_1x1_w.txt";
      x_filename = "wedge_felippa_1x1_x.txt";
    }
    else if ( rule == 1 )
    {
      w_filename = "wedge_felippa_3x2_w.txt";
      x_filename = "wedge_felippa_3x2_x.txt";
    }
    else if ( rule == 2 )
    {
      w_filename = "wedge_felippa_3bx2_w.txt";
      x_filename = "wedge_felippa_3bx2_x.txt";
    }
    else if ( rule == 3 )
    {
      w_filename = "wedge_felippa_6x3_w.txt";
      x_filename = "wedge_felippa_6x3_x.txt";
    }
    else if ( rule == 4 )
    {
      w_filename = "wedge_felippa_6bx2_w.txt";
      x_filename = "wedge_felippa_6bx2_x.txt";
    }
    else if ( rule == 5 )
    {
      w_filename = "wedge_felippa_7x3_w.txt";
      x_filename = "wedge_felippa_7x3_x.txt";
    }
    else if ( rule == 6 )
    {
      w_filename = "wedge_felippa_12x4_w.txt";
      x_filename = "wedge_felippa_12x4_x.txt";
    }

    line_order = line_order_array[rule];
    triangle_order = triangle_order_array[rule];

    order = line_order * fabs ( triangle_order );

    w = new double[order];
    x = new double[dim_num * order];
    wedge_rule ( line_order, triangle_order, w, x );
    r8mat_write ( w_filename, 1, order, w );
    r8mat_write ( x_filename, dim_num, order, x );
    cout << setw(6) << rule << "  "
         << setw(6) << triangle_order << "  "
         << setw(6) << line_order << "  "
         << setw(6) << order << "  "
         << w_filename << "  "
         << x_filename << "\n";

    delete [] w;
    delete [] x;
  }

  return;
}
