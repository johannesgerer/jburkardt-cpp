# include <cstdlib>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hpp.hpp"

int main ( );
void hpp_test01 ( );
void hpp_test015 ( );
void hpp_test02 ( );
void hpp_test03 ( );
void hpp_test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HPP_PRB.
//
//  Discussion:
//
//    HPP_PRB tests the HERMITE_PRODUCT_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HPP_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the HERMITE_PRODUCT_POLYNOMIAL library.\n";

  hpp_test01 ( );
  hpp_test015 ( );
  hpp_test02 ( );
  hpp_test03 ( );
  hpp_test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HPP_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void hpp_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HPP_TEST01 tests routines for the GRLEX ordering of compositions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 2;
  int rank;
  int rank1;
  int rank2;
  int seed;
  int test;
  int *x;
  int x_sum;
  int x_sum_old;

  x = new int[k];

  cout << "\n";
  cout << "HPP_TEST01:\n";
  cout << "  COMP_NEXT_GRLEX is given a composition, and computes the \n";
  cout << "  next composition in grlex order.\n";

  cout << "\n";
  cout << "  Rank   Sum   Components\n";

  for ( i = 0; i < k; i++ )
  {
    x[i] = 0;
  }
  x_sum_old = -1;
  rank = 1;

  for ( ; ; )
  {
    x_sum = i4vec_sum ( k, x );

    if ( x_sum_old < x_sum )
    {
      x_sum_old = x_sum;
      cout << "\n";
    }

    cout << setw(6) << rank << "  "
         << setw(6) << x_sum;
    for ( i = 0; i < k; i++ )
    {
      cout << setw(4) << x[i];
    }
    cout << "\n";

    if ( 20 <= rank )
    {
      break;
    }

    comp_next_grlex ( k, x );
    rank = rank + 1;
  }
  delete [] x;

  cout << "\n";
  cout << "  COMP_UNRANK_GRLEX is given a rank and returns the\n";
  cout << "  corresponding set of multinomial exponents.\n";
  cout << "\n";
  cout << "  Rank   Sum   Components\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    rank = i4_uniform_ab ( 1, 20, seed );
    x = comp_unrank_grlex ( k, rank );
    x_sum = i4vec_sum ( k, x );
    cout << setw(6) << rank << "  "
         << setw(6) << x_sum;
    for ( i = 0; i < k; i++ )
    {
      cout << setw(4) << x[i];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  COMP_RANDOM_GRLEX randomly selects a composition\n";
  cout << "  between given lower and upper ranks.\n";
  cout << "\n";
  cout << "  Rank   Sum   Components\n";
  cout << "\n";

  seed = 123456789;
  rank1 = 5;
  rank2 = 20;

  for ( test = 1; test <= 5; test++ )
  {
    x = comp_random_grlex ( k, rank1, rank2, seed, rank );
    x_sum = i4vec_sum ( k, x );
    cout << setw(6) << rank << "  "
         << setw(6) << x_sum;
    for ( i = 0; i < k; i++ )
    {
      cout << setw(4) << x[i];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  COMP_RANK_GRLEX returns the rank of a given composition.\n";
  cout << "\n";
  cout << "  Rank   Sum   Components\n";
  cout << "\n";

  x = new int[k];

  x[0] = 4;
  x[1] = 0; 
  rank = comp_rank_grlex ( k, x );
  x_sum = i4vec_sum ( k, x );
  cout << setw(6) << rank << "  "
       << setw(6) << x_sum;
  for ( i = 0; i < k; i++ )
  {
    cout << setw(4) << x[i];
  }
  cout << "\n";;

  x[0] = 11;
  x[1] = 5; 
  rank = comp_rank_grlex ( k, x );
  x_sum = i4vec_sum ( k, x );
  cout << setw(6) << rank << "  "
       << setw(6) << x_sum;
  for ( i = 0; i < k; i++ )
  {
    cout << setw(4) << x[i];
  }
  cout << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void hpp_test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HPP_TEST015 tests HEP_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int *e;
  int *f;
  int l[1];
  int m;
  int n;
  int o;
  int o_max;
  ostringstream title;

  m = 1;

  cout << "\n";
  cout <<  "HPP_TEST015:\n";
  cout <<  "  HEP_COEFFICIENTS computes the coefficients and\n";
  cout <<  "  exponents of the Hermite polynomial He(n,x).\n";

  for ( n = 1; n <= 5; n++ )
  {
    o = ( n + 2 ) / 2;
    c = new double[o];
    e = new int[o];
    f = new int[o];

    hep_coefficients ( n, o, c, f );

    l[0] = n;
    o_max = o;

    hepp_to_polynomial ( m, l, o_max, o, c, e );

    cout << "\n";
    title << "  He(" << n << ",x) =";
    polynomial_print ( m, o, c, e, title.str ( ) );
    title.str ( "" );

    delete [] c;
    delete [] e;
    delete [] f;
  }

  return;
}
//****************************************************************************80

void hpp_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HPP_TEST02 tests HEP_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double e;
  int n;
  int n_data;
  int o;
  double x;
  double xvec[1];
  double fx1;
  double *fx2;

  n = 1;

  cout << "\n";
  cout << "HPP_TEST02:\n";
  cout << "  HEP_VALUES stores values of\n";
  cout << "  the Hermite polynomial He(o,x).\n";
  cout << "  HEP_VALUE evaluates a Hermite polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     O        X          He(O,X)                   He(O,X)";
  cout << "                   Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hep_values ( n_data, o, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    xvec[0] = x;

    fx2 = hep_value ( n, o, xvec );

    e = fx1 - fx2[0];

    cout << setw(6) << o << "  "
         << setw(12) << x << "  "
         << setw(24) << fx1 << "  "
         << setw(24) << fx2[0] << "  "
         << setw(8) << e << "\n";

    delete [] fx2;
  }

  return;
}
//****************************************************************************80

void hpp_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HPP_TEST03 tests HEPP_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int *e;
  int i;
  int *l;
  int m = 3;
  int n = 1;
  int o;
  int o_max;
  int rank;
  int seed;
  double *v1;
  double *v2;
  double *x;
  double xhi;
  double xlo;

  cout << "\n";
  cout << "HPP_TEST03:\n";
  cout << "  HEPP_VALUE evaluates a Hermite product polynomial.\n";
  cout << "  POLYNOMIAL_VALUE evaluates a polynomial.\n";

  xlo = -1.0;
  xhi = +1.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, xlo, xhi, seed );

  cout << "\n";
  cout << "  Evaluate at X = ";
  for ( i = 0; i < m; i++ )
  {
    cout << "  " << x[i+0*m];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Rank  I1  I2  I3:  He(I1,X1)*He(I2,X2)*He(I3,X3)    P(X1,X2,X3)\n";
  cout << "\n";

  for ( rank = 1; rank <= 20; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );
//
//  Evaluate the HePP directly.
//
    v1 = hepp_value ( m, n, l, x );
//
//  Convert the HePP to a polynomial.
//
    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = new double[o_max];
    e = new int[o_max];
 
    hepp_to_polynomial ( m, l, o_max, o, c, e );
//
//  Evaluate the polynomial.
//
    v2 = polynomial_value ( m, o, c, e, n, x );
//
//  Compare results.
//
    cout << setw(6) << rank << "  "
         << setw(2) << l[0] << "  "
         << setw(2) << l[1] << "  "
         << setw(2) << l[2] << "  "
         << setw(14) << v1[0] << "  "
         << setw(14) << v2[0] << "\n";
 
    delete [] c;
    delete [] e;
    delete [] l;
    delete [] v1;
    delete [] v2;
  }

  delete [] x;

  return;
}
//****************************************************************************80

void hpp_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HPP_TEST04 tests HEPP_TO_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int *e;
  int i;
  int *l;
  ostringstream label;
  int m = 2;
  int o;
  int o_max;
  int rank;

  cout << "\n";
  cout << "HPP_TEST04:\n";
  cout << "  HEPP_TO_POLYNOMIAL is given a Hermite product polynomial\n";
  cout << "  and determines its polynomial representation.\n";

  cout << "\n";
  cout << "  Using spatial dimension M = " << m << "\n";

  for ( rank = 1; rank <= 11; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );

    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = new double[o_max];
    e = new int[o_max];

    hepp_to_polynomial ( m, l, o_max, o, c, e );

    label << "  HePP #" << rank 
          << " = He(" << l[0]
          << ",X)*He(" << l[1] 
          << ",Y) =";

    cout << "\n";
    polynomial_print ( m, o, c, e, label.str ( ) );
    label.str ( "" );

    delete [] c;
    delete [] e;
    delete [] l;
  }

  return;
}

