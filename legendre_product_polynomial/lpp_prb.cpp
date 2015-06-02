# include <cstdlib>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "lpp.hpp"

int main ( );

void comp_enum_test ( );
void comp_next_grlex_test ( );
void comp_random_grlex_test ( );
void comp_rank_grlex_test ( );
void comp_unrank_grlex_test ( );

void i4_choose_test ( );
void i4_uniform_ab_test ( );

void i4vec_permute_test ( );
void i4vec_print_test ( );
void i4vec_sort_heap_index_a_test ( );
void i4vec_sum_test ( );
void i4vec_uniform_ab_new_test ( );

void mono_next_grlex_test ( );
void mono_print_test ( );
void mono_rank_grlex_test ( );
void mono_unrank_grlex_test ( );
void mono_upto_enum_test ( );
void mono_upto_next_grlex_test ( );
void mono_upto_random_test ( );

void lp_coefficients_test ( );
void lp_value_test ( );
void lp_values_test ( );

void lpp_to_polynomial_test ( );
void lpp_value_test ( );

void perm_uniform_test ( );

void polynomial_compress_test ( );
void polynomial_print_test ( );
void polynomial_sort_test ( );
void polynomial_value_test ( );

void r8mat_print_test ( );
void r8mat_print_some_test ( );
void r8mat_uniform_ab_new_test ( );

void r8vec_permute_test ( );
void r8vec_print_test ( );
void r8vec_uniform_ab_new_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LPP_PRB.
//
//  Discussion:
//
//    LPP_PRB tests the LEGENDRE_PRODUCT_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LPP_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LEGENDRE_PRODUCT_POLYNOMIAL library.\n";

  i4_choose_test ( );
  i4_uniform_ab_test ( );

  i4vec_permute_test ( );
  i4vec_print_test ( );
  i4vec_sort_heap_index_a_test ( );
  i4vec_sum_test ( );
  i4vec_uniform_ab_new_test ( );

  r8vec_permute_test ( );
  r8vec_print_test ( );
  r8vec_uniform_ab_new_test ( );

  r8mat_print_test ( );
  r8mat_print_some_test ( );
  r8mat_uniform_ab_new_test ( );

  perm_uniform_test ( );

  comp_enum_test ( );
  comp_next_grlex_test ( );
  comp_random_grlex_test ( );
  comp_rank_grlex_test ( );
  comp_unrank_grlex_test ( );

  mono_next_grlex_test ( );
  mono_print_test ( );
  mono_rank_grlex_test ( );
  mono_unrank_grlex_test ( );
  mono_upto_enum_test ( );
  mono_upto_next_grlex_test ( ); 
  mono_upto_random_test ( );

  polynomial_compress_test ( );
  polynomial_print_test ( );
  polynomial_sort_test ( );
  polynomial_value_test ( );

  lp_coefficients_test ( );
  lp_value_test ( );
  lp_values_test ( );

  lpp_to_polynomial_test ( );
  lpp_value_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LPP_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void comp_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_ENUM_TEST tests COMP_ENUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int num;
  int k;
  int n;

  cout << "\n";
  cout << "COMP_ENUM_TEST\n";
  cout << "  COMP_ENUM counts compositions;\n";
  cout << "\n";
  for ( n = 0; n <= 10; n++ )
  {
    for ( k = 1; k <= 10; k++ )
    {
      num = comp_enum ( n, k );
      cout << "  " << setw(6) << num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void comp_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT_GRLEX_TEST tests COMP_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int kc = 3;
  int nc;
  int rank;
  int xc[3];

  cout << "\n";
  cout << "COMP_NEXT_GRLEX_TEST\n";
  cout << "  A COMP is a composition of an integer N into K parts.\n";
  cout << "  Each part is nonnegative.  The order matters.\n";
  cout << "  COMP_NEXT_GRLEX determines the next COMP in\n";
  cout << "  graded lexicographic (grlex) order.\n";
  
  cout << "\n";
  cout << "  Rank:     NC       COMP\n";
  cout << "  ----:     --   ------------\n";

  for ( rank = 1; rank <= 71; rank++ )
  {
    if ( rank == 1 )
    {
      for ( j = 0; j < kc; j++ )
      {
        xc[j] = 0;
      }
    }
    else
    {
      comp_next_grlex ( kc, xc );
    }

    nc = i4vec_sum ( kc, xc );

    cout << "   " << setw(3) << rank << ": ";
    cout << "    " << setw(2) << nc << " = ";
    for ( j = 0; j < kc - 1; j++ )
    {
      cout << setw(2) << xc[j] << " + ";
    }
    cout << setw(2) << xc[kc-1] << "\n";
//
//  When XC(1) == NC, we have completed the compositions associated with
//  a particular integer, and are about to advance to the next integer.
//
    if ( xc[0] == nc )
    {
      cout << "  ----:     --   ------------\n";
    }
  }

  return;
}
//****************************************************************************80

void comp_random_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANDOM_GRLEX_TEST tests COMP_RANDOM_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int kc;
  int nc;
  int rank;
  int rank1;
  int rank2;
  int seed;
  int test;
  int *xc;

  cout << "\n";
  cout << "COMP_RANDOM_GRLEX_TEST\n";
  cout << "  A COMP is a composition of an integer N into K parts.\n";
  cout << "  Each part is nonnegative.  The order matters.\n";
  cout << "  COMP_RANDOM_GRLEX selects a random COMP in\n";
  cout << "  graded lexicographic (grlex) order between indices RANK1 and RANK2.\n";
  cout << "\n";

  kc = 3;
  rank1 = 20;
  rank2 = 60;
  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    xc = comp_random_grlex ( kc, rank1, rank2, seed, rank );
    nc = i4vec_sum ( kc, xc );

    cout << "   " << setw(3) << rank << ": ";
    cout << "    " << setw(2) << nc << " = ";
    for ( j = 0; j < kc - 1; j++ )
    {
      cout << setw(2) << xc[j] << " + ";
    }
    cout << setw(2) << xc[kc-1] << "\n";
    delete [] xc;
  }

  return;
}
//****************************************************************************80

void comp_rank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANK_GRLEX_TEST tests COMP_RANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int kc;
  int nc;
  int rank1;
  int rank2;
  int rank3;
  int rank4;
  int seed;
  int test;
  int *xc;

  cout << "\n";
  cout << "COMP_RANK_GRLEX_TEST\n";
  cout << "  A COMP is a composition of an integer N into K parts.\n";
  cout << "  Each part is nonnegative.  The order matters.\n";
  cout << "  COMP_RANK_GRLEX determines the rank of a COMP\n";
  cout << "  from its parts.\n";
  cout << "\n";
  cout << "        Actual  Inferred\n";
  cout << "  Test    Rank      Rank\n";
  cout << "\n";

  kc = 3;
  rank1 = 20;
  rank2 = 60;
  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    xc = comp_random_grlex ( kc, rank1, rank2, seed, rank3 );
    rank4 = comp_rank_grlex ( kc, xc );

    cout << "  " << setw(4) << test;
    cout << "  " << setw(6) << rank3;
    cout << "  " << setw(8) << rank4 << "\n";

    delete [] xc;
  }
  return;
}
//****************************************************************************80

void comp_unrank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_UNRANK_GRLEX_TEST tests COMP_UNRANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int kc = 3;
  int nc;
  int rank1;
  int rank2;
  int *xc;

  cout << "\n";
  cout << "COMP_UNRANK_GRLEX_TEST\n";
  cout << "  A COMP is a composition of an integer N into K parts.\n";
  cout << "  Each part is nonnegative.  The order matters.\n";
  cout << "  COMP_UNRANK_GRLEX determines the parts\n";
  cout << "  of a COMP from its rank.\n";
 
  cout << "\n";
  cout << "  Rank: ->  NC       COMP\n";
  cout << "  ----:     --   ------------\n";

  for ( rank1 = 1; rank1 <= 71; rank1++ )
  {
    xc = comp_unrank_grlex ( kc, rank1 );
    nc = i4vec_sum ( kc, xc );

    cout << "   " << setw(3) << rank1 << ": ";
    cout << "    " << setw(2) << nc << " = ";
    for ( j = 0; j < kc - 1; j++ )
    {
      cout << setw(2) << xc[j] << " + ";
    }
    cout << setw(2) << xc[kc-1] << "\n";
//
//  When XC(1) == NC, we have completed the compositions associated with
//  a particular integer, and are about to advance to the next integer.
//
    if ( xc[0] == nc )
    {
      cout << "  ----:     --   ------------\n";
    }
    delete [] xc;
  }
  return;
}
//****************************************************************************80

void i4_choose_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_TEST tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int cnk;
  int k;
  int n;

  cout << "\n";
  cout << "I4_CHOOSE_TEST\n";
  cout << "  I4_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "       N       K     CNK\n";

  for ( n = 0; n <= 4; n++ )
  {
    cout << "\n";
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      cout                   << "  "
           << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void i4_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed = 123456789;

  cout << "\n";
  cout << "I4_UNIFORM_AB_TEST\n";
  cout << "  I4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, seed );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";
  }

  return;
}
//****************************************************************************80

void i4vec_permute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PERMUTE_TEST tests I4VEC_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int b;
  int c;
  int n = 12;
  int *p;
  int seed;

  cout << "\n";
  cout << "I4VEC_PERMUTE_TEST\n";
  cout << "  I4VEC_PERMUTE reorders an integer vector\n";
  cout << "  according to a given permutation.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  b = 0;
  c = n;
  seed = 123456789;
  a = i4vec_uniform_ab_new ( n, b, c, seed );

  i4vec_print ( n, a, "  A, before rearrangement:" );

  p = perm_uniform_new ( n, seed );

  i4vec_print ( n, p, "  Permutation vector P:" );

  i4vec_permute ( n, p, a );

  i4vec_print ( n, a, "  A, after rearrangement:" );

  delete [] a;
  delete [] p;

  return;
}
//****************************************************************************80

void i4vec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT_TEST tests I4VEC_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n = 4;
  int v[4] = { 91, 92, 93, 94 };

  cout << "\n";
  cout << "I4VEC_PRINT_TEST\n";
  cout << "  I4VEC_PRINT prints an I4VEC\n";

  i4vec_print ( n, v, "  Here is the I4VEC:" );

  return;
}
//****************************************************************************80

void i4vec_sort_heap_index_a_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_INDEX_A_TEST tests I4VEC_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int b;
  int c;
  int i;
  int *indx;
  int n = 20;
  int seed;

  cout << "\n";
  cout << "I4VEC_SORT_HEAP_INDEX_A_TEST\n";
  cout << "  I4VEC_SORT_HEAP_INDEX_A creates an ascending\n";
  cout << "  sort index for an I4VEC.\n";

  b = 0;
  c = 3 * n;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, b, c, seed );

  i4vec_print ( n, a, "  Unsorted array:" );

  indx = i4vec_sort_heap_index_a ( n, a );

  i4vec_print ( n, indx, "  Sort vector INDX:" );

  cout << "\n";
  cout << "       I   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout                          << "  "
         << setw(8) <<        i   << "  "
         << setw(8) <<   indx[i]  << "  "
         << setw(8) << a[indx[i]] << "\n";
  }

  delete [] a;
  delete [] indx;

  return;
}
//****************************************************************************80

void i4vec_sum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM_TEST tests I4VEC_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int hi;
  int lo;
  int n;
  int s;
  int seed;

  cout << "\n";
  cout << "I4VEC_SUM_TEST\n";
  cout << "  I4VEC_SUM sums the entries of an I4VEC.\n";

  n = 5;
  lo = 0;
  hi = 10;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, lo, hi, seed );
  i4vec_print ( n, a, "  The vector:" );

  s = i4vec_sum ( n, a );
  cout << "\n";
  cout << "  The vector entries sum to " << s << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void i4vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int a = -100;
  int b = 200;
  int n = 20;
  int seed = 123456789;
  int *v;

  cout << "\n";
  cout << "I4VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  I4VEC_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  v = i4vec_uniform_ab_new ( n, a, b, seed );

  i4vec_print ( n, v, "  The random vector:" );

  delete [] v;

  return;
}
//****************************************************************************80

void lp_coefficients_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LP_COEFFICIENTS_TEST tests LP_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int *e;
  int *f;
  int i;
  ostringstream label;
  int m = 1;
  int n;
  int n_max = 10;
  int o;

  cout << "\n";
  cout << "LP_COEFFICIENTS_TEST\n";
  cout << "  LP_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).\n";
  cout << "\n";

  for ( n = 0; n <= n_max; n++ )
  {
    c = new double[n+1];
    f = new int[n+1];

    lp_coefficients ( n, o, c, f );

    e = new int[o];
    for ( i = 0; i < o; i++ )
    {
      e[i] = f[i] + 1;
    }
    
    label << "  P(" << n << ",x) = ";
    polynomial_print ( m, o, c, e, label.str ( ) );

    delete [] c;
    delete [] e;
    delete [] f;

    label.str ( "" );
   }

  return;
}
//****************************************************************************80

void lp_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LP_VALUE_TEST tests LP_VALUE.
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
  cout << "LP_VALUE_TEST:\n";
  cout << "  LP_VALUE evaluates a Legendre polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     O        X           L(O,X)                    L(O,X)";
  cout << "                   Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lp_values ( n_data, o, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    xvec[0] = x;

    fx2 = lp_value ( n, o, xvec );

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

void lp_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LP_VALUES_TEST tests LP_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  int o;
  double x;
  double fx;

  cout << "\n";
  cout << "LP_VALUES_TEST:\n";
  cout << "  LP_VALUES stores values of\n";
  cout << "  the Legendre polynomial P(o,x).\n";
  cout << "\n";
  cout << "                        Tabulated\n";
  cout << "     O        X           L(O,X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lp_values ( n_data, o, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    cout << setw(6) << o << "  "
         << setw(12) << x << "  "
         << setw(24) << fx << "\n";
  }

  return;
}
//****************************************************************************80

void lpp_to_polynomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LPP_TO_POLYNOMIAL_TEST tests LPP_TO_POLYNOMIAL.
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
  cout << "LPP_TO_POLYNOMIAL_TEST:\n";
  cout << "  LPP_TO_POLYNOMIAL is given a Legendre product polynomial\n";
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

    lpp_to_polynomial ( m, l, o_max, o, c, e );

    label << "  LPP #" << rank 
          << " = L(" << l[0]
          << ",X)*L(" << l[1] 
          << ",Y) = \n";

    cout << "\n";
    polynomial_print ( m, o, c, e, label.str ( ) );
    label.str ( "" );

    delete [] c;
    delete [] e;
    delete [] l;
  }

  return;
}
//****************************************************************************80

void lpp_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LPP_VALUE_TEST tests LPP_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2014
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
  cout << "LPP_VALUE_TEST:\n";
  cout << "  LPP_VALUE evaluates a Legendre product polynomial.\n";

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
  cout << "  Rank  I1  I2  I3:  L(I1,X1)*L(I2,X2)*L(I3,X3)    P(X1,X2,X3)\n";
  cout << "\n";

  for ( rank = 1; rank <= 20; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );
//
//  Evaluate the LPP directly.
//
    v1 = lpp_value ( m, n, l, x );
//
//  Convert the LPP to a polynomial.
//
    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = new double[o_max];
    e = new int[o_max];
 
    lpp_to_polynomial ( m, l, o_max, o, c, e );
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

void mono_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GRLEX_TEST tests MONO_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int m = 4;
  int i;
  int j;
  int k;
  int seed;
  int *x;

  cout << "\n";
  cout << "MONO_NEXT_GRLEX_TEST\n";
  cout << "  MONO_NEXT_GRLEX computes the next monomial\n";
  cout << "  in M variables, in grlex order.\n";
  cout << "\n";
  cout << "  Let M =  " << m << "\n";

  a = 0;
  b = 3;
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = i4vec_uniform_ab_new ( m, a, b, seed );
    cout << "\n";
    cout << "  ";
    for ( k = 0; k < m; k++ )
    {
      cout << setw(2) << x[k];
    }
    cout << "\n";

    for ( j = 1; j <= 5; j++ )
    {
      mono_next_grlex ( m, x );
      cout << "  ";
      for ( k = 0; k < m; k++ )
      {
        cout << setw(2) << x[k];
      }
      cout << "\n";
    }
    delete [] x;
  }

  return;
}
//****************************************************************************80

void mono_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_PRINT_TEST tests MONO_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int f1[1] = { 5 };
  int f2[1] = { -5 };
  int f3[4] = { 2, 1, 0, 3 };
  int f4[3] = { 17, -3, 199 };
  int m;

  cout << "\n";
  cout << "MONO_PRINT_TEST\n";
  cout << "  MONO_PRINT can print out a monomial.\n";
  cout << "\n";

  m = 1;
  mono_print ( m, f1, "  Monomial [5]:" );

  m = 1;
  mono_print ( m, f2, "  Monomial [5]:" );

  m = 4;
  mono_print ( m, f3, "  Monomial [2,1,0,3]:" );

  m = 3;
  mono_print ( m, f4, "  Monomial [17,-3,199]:" );

  return;
}
//****************************************************************************80

void mono_rank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_RANK_GRLEX_TEST tests MONO_RANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int test;
  int test_num = 8;
  int x[3];
  int x_test[3*8] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 0, 1, 
    0, 2, 0, 
    1, 0, 2, 
    0, 3, 1, 
    3, 2, 1, 
    5, 2, 1 };

  cout << "\n";
  cout << "MONO_RANK_GRLEX_TEST\n";
  cout << "  MONO_RANK_GRLEX returns the rank of a monomial in the sequence\n";
  cout << "  of all monomials in M dimensions, in grlex order.\n";

  cout << "\n";
  cout << "  Print a monomial sequence with ranks assigned.\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(3) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  cout << "\n";
  cout << "  Now, given a monomial, retrieve its rank in the sequence:\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    for ( j = 0; j < m; j++ )
    {
      x[j] = x_test[j+test*m];
    }
    rank = mono_rank_grlex ( m, x );

    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void mono_unrank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UNRANK_GRLEX_TEST tests MONO_UNRANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int rank_max;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_UNRANK_GRLEX_TEST\n";
  cout << "  MONO_UNRANK_GRLEX is given a rank, and returns the corresponding\n";
  cout << "  monomial in the sequence of all monomials in M dimensions\n";
  cout << "  in grlex order.\n";

  cout << "\n";
  cout << "  For reference, print a monomial sequence with ranks.\n";

  n = 4;
  rank_max = mono_upto_enum ( m, n );

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x = new int[n];
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(3) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  cout << "\n";
  cout << "  Now choose random ranks between 1 and " << rank_max << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    rank = i4_uniform_ab ( 1, rank_max, seed );
    x = mono_unrank_grlex ( m, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    free ( x );
  }

  return;
}
//****************************************************************************80

void mono_upto_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_ENUM_TEST tests MONO_UPTO_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int v;

  cout << "\n";
  cout << "MONO_UPTO_ENUM_TEST\n";
  cout << "  MONO_UPTO_ENUM can enumerate the number of monomials\n";
  cout << "  in M variables, of total degree 0 up to N.\n";

  cout << "\n";
  cout << "    N:\n";
  for ( n = 0; n <= 8; n++ )
  {
    cout << "  " << setw(4) << n;
  }
  cout << "\n";
  cout << "   m +------------------------------------------------------\n";
  for ( m = 1; m <= 8; m++ )
  {
    cout << "  " << setw(2) << m << "  |";
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_upto_enum ( m, n );
      cout << " " << setw(5) << v;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void mono_upto_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GRLEX_TEST tests MONO_UPTO_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  cout << "\n";
  cout << "MONO_UPTO_NEXT_GRLEX_TEST\n";
  cout << "  MONO_UPTO_NEXT_GRLEX can list the monomials\n";
  cout << "  in M variables, of total degree up to N,\n";
  cout << "  in grlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,0).\n";
  cout << "  The process ends with (N,0,...,0,0)\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_upto_random_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_RANDOM_TEST tests MONO_UPTO_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int j;
  int n;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_UPTO_RANDOM_TEST\n";
  cout << "  MONO_UPTO_RANDOM selects at random a monomial\n";
  cout << "  in M dimensions of total degree no greater than N.\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_upto_random ( m, n, seed, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  return;
}
//****************************************************************************80

void mono_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_VALUE_TEST tests MONO_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int *f;
  int j;
  int m = 3;
  int n;
  int nx = 2;
  int rank;
  int seed;
  int test;
  int test_num;
  double *v;
  double x[3*2] = {
     1.0, 2.0, 3.0, 
    -2.0, 4.0, 1.0 };

  cout << "\n";
  cout << "MONO_VALUE_TEST\n";
  cout << "  MONO_VALUE evaluates a monomial.\n";

  n = 6;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    f = mono_upto_random ( m, n, seed, rank );
    cout << "\n";
    mono_print ( m, f, "  M(X) = " );
    v = mono_value ( m, nx, f, x );
    for ( j = 0; j < nx; j++ )
    {
      cout << "  M(" << x[0+j*m] 
           << "," << x[1+j*m]
           << "," << x[2+j*m]
           << ") = " << v[j] << "\n";
    }
    delete [] f;
    delete [] v;
  }

  return;
}
//****************************************************************************80

void perm_uniform_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_TEST tests PERM_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  cout << "\n";
  cout << "PERM_UNIFORM_TEST\n";
  cout << "  PERM_UNIFORM randomly selects a permutation.\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm_uniform_new ( n, seed );
    cout << "  ";
    for ( i = 0; i < n; i++ )
    {
      cout << setw(4) << p[i];
    }
    cout << "\n";
    delete [] p;
  }
  return;
}
//****************************************************************************80

void polynomial_compress_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_COMPRESS_TEST tests POLYNOMIAL_COMPRESS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[10] = { 7.0, - 5.0, 5.0, 9.0, 11.0, 3.0, 6.0, 0.0, - 13.0, 1.0E-20 };
  double c2[10];
  int m = 3;
  int e[10] = { 1, 2, 2, 4, 5, 5, 5, 12, 33, 35 }; 
  int e2[10];
  int j;
  int nx = 2;
  int o = 10;
  int o2;
  string title;

  cout << "\n";
  cout << "POLYNOMIAL_COMPRESS_TEST\n";
  cout << "  POLYNOMIAL_COMPRESS compresses a polynomial.\n";

  cout << "\n";
  title = "  Uncompressed P(X) = ";
  polynomial_print ( m, o, c, e, title );

  polynomial_compress ( o, c, e, o2, c2, e2 );

  cout << "\n";
  title = "  Compressed P(X) = ";
  polynomial_print ( m, o2, c2, e2, title );

  return;
}
//****************************************************************************80

void polynomial_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_PRINT_TEST tests POLYNOMIAL_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = {7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 };
  int o = 6;
  string title = "  P1(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_PRINT_TEST\n";
  cout << "  POLYNOMIAL_PRINT prints a polynomial.\n";

  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_sort_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SORT_TEST tests POLYNOMIAL_SORT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 0.0, 9.0, -5.0, - 13.0, 7.0, 11.0 };
  int m = 3;
  int e[6] = { 12, 4, 2, 33, 1, 5 }; 
  int o = 6;
  string title;

  cout << "\n";
  cout << "POLYNOMIAL_SORT_TEST\n";
  cout << "  POLYNOMIAL_SORT sorts a polynomial by exponent index.\n";

  cout << "\n";
  title = "  Unsorted polynomial:";
  polynomial_print ( m, o, c, e, title );

  polynomial_sort ( o, c, e );

  cout << "\n";
  title = "  Sorted polynomial:";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_VALUE_TEST tests POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 }; 
  int j;
  int nx = 2;
  int o = 6;
  double *p;
  char title[] = "  P(X) =";
  double x[3*2] = {
     1.0, 2.0, 3.0,
    -2.0, 4.0, 1.0 };

  cout << "\n";
  cout << "POLYNOMIAL_VALUE_TEST\n";
  cout << "  POLYNOMIAL_VALUE evaluates a polynomial.\n";

  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  p = polynomial_value ( m, o, c, e, nx, x );

  cout << "\n";
  for ( j = 0; j < nx; j++ )
  {
    cout << "  P(" << x[0+j*m]
         << "," << x[1+j*m]
         << "," << x[2+j*m]
         << ") = " << p[j] << "\n";
  }

  delete [] p;

  return;
}
//****************************************************************************80

void r8mat_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_TEST tests R8MAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_TEST\n";
  cout << "  R8MAT_PRINT prints an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print ( m, n, a, "  The R8MAT:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_print_some_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_SOME_TEST\n";
  cout << "  R8MAT_PRINT_SOME prints some of an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print_some ( m, n, a, 2, 1, 4, 2, "  The R8MAT, rows 2:4, cols 1:2:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 4

  double *a;
  double b = 2.0E+00;
  double c = 10.0E+00;
  int seed = 123456789;

  cout << "\n";
  cout << "R8MAT_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8MAT_UNIFORM_AB_NEW returns a random R8MAT in [A,B].\n";
  cout << "\n";

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  The random R8MAT:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8vec_permute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PERMUTE_TEST tests R8VEC_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int p[5] = { 1, 3, 4, 0, 2 };
  double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  cout << "\n";
  cout << "R8VEC_PERMUTE_TEST\n";
  cout << "  R8VEC_PERMUTE permutes an R8VEC in place.\n";

  r8vec_print ( n, x, "  Original Array X[]:" );

  i4vec_print ( n, p, "  Permutation Vector P[]:" );

  r8vec_permute ( n, p, x );

  r8vec_print ( n, x, "  Permuted array X[P[]]:" );

  return;
}
//****************************************************************************80

void r8vec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_TEST tests R8VEC_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  cout << "\n";
  cout << "TEST1335\n";
  cout << "  R8VEC_PRINT prints an R8VEC.\n";

  r8vec_print ( n, a, "  The R8VEC:" );

  return;
}
//****************************************************************************80

void r8vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a = 10.0;
  double b = 20.0;
  int j;
  double *r;
  int seed;

  cout << "\n";
  cout << "R8VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM returns a random R8VEC\n";
  cout << "  with entries in a given range [ A, B ]\n";
  cout << "\n";
  cout << "  For this problem:\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    cout << "\n";
    cout << "  Input SEED = " << seed << "\n";
    cout << "\n";

    r = r8vec_uniform_ab_new ( N, a, b, seed );

    r8vec_print ( N, r, "  Random R8VEC:" );

    delete [] r;
  }

  return;
# undef N
}
