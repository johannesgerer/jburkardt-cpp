# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "combo.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );
void test40 ( );
void test41 ( );
void test42 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMBO_PRB.
//
//  Discussion:
//
//    COMBO_PRB tests the COMBO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "COMBO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the COMBO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMBO_PRB\n";
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
//    TEST01 tests BAL_SEQ_ENUM, BAL_SEQ_RANK, BAL_SEQ_SUCCESSOR, BAL_SEQ_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nseq;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Balanced sequences:\n";
  cout << "\n";
  cout << "  BAL_SEQ_ENUM enumerates,\n";
  cout << "  BAL_SEQ_RANK ranks,\n";
  cout << "  BAL_SEQ_SUCCESSOR lists,\n";
  cout << "  BAL_SEQ_UNRANK unranks.\n";
//
//  Enumerate.
//
  nseq = bal_seq_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of balanced sequences is " << nseq << "\n";
  cout << "\n";
//
//  List.
//
  t = new int[2*n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    bal_seq_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    cout << "  " << setw(4) << rank;
    for ( i = 0; i < 2*n; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nseq / 2;

  delete [] t;

  t = bal_seq_unrank ( rank, n );
  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( 2 * n, t, " " );
//
//  Rank.
//
  rank = bal_seq_rank ( n, t );

  i4vec_transpose_print ( 2 * n, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as: " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BAL_SEQ_TO_TABLEAU, TABLEAU_TO_BAL_SEQ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 4;
  int rank;
  int *t;
  int *tab;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  BAL_SEQ_TO_TABLEAU converts a balanced\n";
  cout << "  sequence to a tableau;\n";
  cout << "  TABLEAU_TO_BAL_SEQ converts a tableau\n";
  cout << "  to a balanced sequence.\n";
//
//  Pick a random balanced sequence.
//
  rank = 7;

  t = bal_seq_unrank ( rank, n );

  cout << "\n";
  cout << "  Random balanced sequence:\n";
  cout << "\n";
  i4vec_transpose_print ( 2 * n, t, " " );
//
//  Convert to a tableau.
//
  tab = bal_seq_to_tableau ( n, t );

  i4mat_print ( 2, n, tab, "  Corresponding tableau" );
//
//  Convert to a balanced sequence.
//
  delete [] t;

  t = tableau_to_bal_seq ( n, tab );

  i4vec_transpose_print ( 2 * n, t, "  Corresponding balanced sequence:" );

  delete [] t;
  delete [] tab;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BELL_NUMBERS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int bn;
  int n;
  int n_data;

  n_data = 0;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  BELL_NUMBERS computes Bell numbers.\n";
  cout << "\n";
  cout << "     N          BELL(N)      BELL_NUMBERS(N)\n";
  cout << "\n";
  for ( ; ; )
  {
    bell_values ( &n_data, &n, &bn );

    if ( n_data == 0 )
    {
      break;
    }
    b = bell_numbers ( n );
    cout << "  " << setw(8) << n
         << "  " << setw(12) << bn
         << "  " << setw(12) << b[n] << "\n";
    delete [] b;
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  I4_CHOOSE computes binomial coefficients.\n";

  for ( i = -1; i <= 5; i++ )
  {
    for ( j = - 1; j <= 5; j++ )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(4) << j
           << "  " << setw(12) << i4_choose ( i, j ) << "\n";
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
//    TEST05 tests CYCLE_TO_PERM, PERM_TO_CYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int jlo;
  int *index;
  int n = 7;
  int ncycle;
  int nperm;
  int *p;
  int rank;
  int *t;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  CYCLE_TO_PERM converts a permutation from\n";
  cout << "  cycle to array form;\n";
  cout << "  PERM_TO_CYCLE converts a permutation from\n";
  cout << "  array to cycle form.\n";
//
//  Enumerate.
//
  nperm = perm_enum ( n );
//
//  Choose a "random" permutation.
//
  rank = nperm / 2;

  p = perm_lex_unrank ( rank, n );

  perm_print ( n, p, "  Random permutation:" );
//
//  Convert the permutation to cycle form.
//
  t = new int[n];
  index = new int[n];

  perm_to_cycle ( n, p, ncycle, t, index );

  cout << "\n";
  cout << "  Corresponding cycle form:\n";
  cout << "  Number of cycles is " << ncycle << "\n";
  cout << "\n";
  jlo = 0;
  for ( i = 1; i <= ncycle; i++ )
  {
    for ( j = jlo + 1; j <= jlo + index[i-1]; j++ )
    {
      cout << "  " << setw(4) << t[j-1];
    }
    cout << "\n";
    jlo = jlo + index[i-1];
  }
//
//  Convert the set partition back to an RGF.
//
  delete [] p;

  p = cycle_to_perm ( n, ncycle, t, index );

  perm_print ( n, p, "  Corresponding permutation:" );

  delete [] index;
  delete [] p;
  delete [] t;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DIST_ENUM and DIST_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int idist;
  int k = 3;
  int m;
  bool more;
  int num_dist;
  int *q;

  k = 3;
  m = 5;

  num_dist = dist_enum ( k, m );

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a distribution of M indistinguishable\n";
  cout << "  objects among K distinguishable slots:\n";
  cout << "\n";
  cout << "  DIST_ENUM enumerates them;\n";
  cout << "  DIST_NEXT produces the \"next\" one.\n";
  cout << "\n";
  cout << "  Number of:\n";
  cout << "    indistinguishable objects = " << m << "\n";
  cout << "    distinguishable slots =     " << k << "\n";
  cout << "    distributions is            " << num_dist << "\n";
  cout << "\n";

  idist = 0;
  more = false;
  q = new int[k];

  for ( ; ; )
  {
    dist_next ( k, m, q, more );

    if ( !more )
    {
      break;
    }

    idist = idist + 1;
    cout << "  " << setw(4) << idist;
    for ( i = 0; i < k; i++ )
    {
      cout << "  " << setw(2) << q[i];
    }
    cout << "\n";
  }

  delete [] q;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests I4_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int fx;
  int fx2;
  int n;
  int x;

  cout << "\n";
  cout << "TEST07:\n";
  cout << "  I4_FACTORIAL evaluates the factorial function.\n";
  cout << "\n";
  cout << "     X       Exact F       FACTORIAL(X)\n";
  cout << "\n";

  n = 0;

  for ( ; ; )
  {
    i4_factorial_values ( n, x, fx );

    if ( n == 0 )
    {
      break;
    }

    if ( x <= 0.0 )
    {
      continue;
    }

    fx2 = i4_factorial ( x );

    cout << "  " << setw(4) << x
         << "  " << setw(12) << fx
         << "  " << setw(12) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests GRAY_CODE_*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int ngray;
  int rank;
  int rank_old;
  int *t;

  t = new int[n];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Gray codes:\n";
  cout << "\n";
  cout << "  GRAY_CODE_ENUM enumerates,\n";
  cout << "  GRAY_CODE_RANK ranks,\n";
  cout << "  GRAY_CODE_SUCCESSOR lists,\n";
  cout << "  GRAY_CODE_UNRANK unranks.\n";
//
//  Enumerate.
//
  ngray = gray_code_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of Gray code elements is " << ngray << "\n";
  cout << "\n";
//
//  List
//
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    gray_code_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = ngray / 2;

  delete [] t;

  t = gray_code_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << t[i];
  }
  cout << "\n";
//
//  Rank.
//
  rank = gray_code_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );

  cout << "\n";
  cout << "  Computed rank: " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests I4VEC_SEARCH_BINARY_A and I4VEC_SORT_INSERT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int a[10] = { 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 };
  int b;
  int index;
  int n = 10;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  Integer vectors:\n";
  cout << "\n";
  cout << "  I4VEC_SORT_INSERT_A ascending sorts;\n";
  cout << "  I4VEC_SEARCH_BINARY_A searches a ascending sorted vector.\n";

  i4vec_print ( n, a, "  Before ascending sort:" );

  i4vec_sort_insert_a ( n, a );

  i4vec_print ( n, a, "  After ascending sort:" );

  b = 5;

  cout << "\n";
  cout << "  Now search for an instance of the value " << b << "\n";

  index = i4vec_search_binary_a ( n, a, b );

  cout << "\n";
  if ( index == 0 )
  {
    cout << "  The value does not occur.\n";
  }
  else
  {
    cout << "  The value occurs at index = " << index << "\n";
  }
  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests I4VEC_SEARCH_BINARY_D and I4VEC_SORT_INSERT_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N] = { 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 };
  int b;
  int index;
  int n = N;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  Integer vectors:\n";
  cout << "\n";
  cout << "  I4VEC_SORT_INSERT_D descending sorts;\n";
  cout << "  I4VEC_SEARCH_BINARY_D searches a descending \n";
  cout << " sorted vector.\n";

  i4vec_print ( n, a, "  Before descending sort:" );

  i4vec_sort_insert_d ( n, a );

  i4vec_print ( n, a, "  After descending sort:" );

  b = 5;

  cout << "\n";
  cout << "  Now search for an instance of the value " << b << "\n";

  index = i4vec_search_binary_d ( n, a, b );

  cout << "\n";
  if ( index == 0 )
  {
    cout << "  The value does not occur.\n";
  }
  else
  {
    cout << "  The value occurs at index = " << index << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests KNAPSACK_REORDER and KNAPSACK_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  double mass;
  double mass_limit = 26.0;
  int n = N;
  double p[N] = { 24.0, 13.0, 23.0, 15.0, 16.0 };
  double profit;
  double w[N] = { 12.0,  7.0, 11.0,  8.0,  9.0 };
  double x[N];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  KNAPSACK_REORDER reorders the knapsack data.\n";
  cout << "  KNAPSACK_01 solves the 0/1 knapsack problem.\n";

  cout << "\n";
  cout << "  Object, Profit, Mass, Profit Density\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(7) << p[i]
         << "  " << setw(7) << w[i]
         << "  " << setw(7) << p[i] / w[i] << "\n";
  }

  knapsack_reorder ( n, p, w );

  cout << "\n";
  cout << "  After reordering by Profit Density:\n";
  cout << "\n";
  cout << "  Object, Profit, Mass, Profit Density\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(7) << p[i]
         << "  " << setw(7) << w[i]
         << "  " << setw(7) << p[i] / w[i] << "\n";
  }

  cout << "\n";
  cout << "  Total mass restriction is " << mass_limit << "\n";

  knapsack_01 ( n, mass_limit, p, w, x, mass, profit );

  cout << "\n";
  cout << "  Object, Density, Choice, Profit, Mass\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(7) << p[i] / w[i]
         << "  " << setw(7) << x[i]
         << "  " << setw(7) << x[i] * p[i]
         << "  " << setw(7) << x[i] * w[i] << "\n";
  }

  cout << "\n";
  cout << "  Total:            " << profit 
    << "  " << mass << "\n";

  return;
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests KNAPSACK_REORDER and KNAPSACK_RATIONAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  double mass;
  double mass_limit = 26.0;
  int n = N;
  double p[N] = { 24.0, 13.0, 23.0, 15.0, 16.0 };
  double profit;
  double w[N] = { 12.0,  7.0, 11.0,  8.0,  9.0 };
  double x[N];

  cout << "\n";
  cout << "TEST12\n";
  cout << "  KNAPSACK_REORDER reorders the knapsack data.\n";
  cout << "  KNAPSACK_RATIONAL solves the rational knapsack problem.\n";

  cout << "\n";
  cout << "  Object, Profit, Mass, Profit Density\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i + 1
         << "  " << setw(7) << p[i]
         << "  " << setw(7) << w[i]
         << "  " << setw(7) << p[i] / w[i] << "\n";
  }

  knapsack_reorder ( n, p, w );

  cout << "\n";
  cout << "  After reordering by Profit Density:\n";
  cout << "\n";
  cout << "  Object, Profit, Mass, Profit Density\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i + 1
         << "  " << setw(7) << p[i]
         << "  " << setw(7) << w[i]
         << "  " << setw(7) << p[i] / w[i] << "\n";
  }

  cout << "\n";
  cout << "  Total mass restriction is " << mass_limit << "\n";

  knapsack_rational ( n, mass_limit, p, w, x, mass, profit );

  cout << "\n";
  cout << "  Object, Density, Choice, Profit, Mass\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i + 1
         << "  " << setw(7) << p[i] / w[i]
         << "  " << setw(7) << x[i] * p[i]
         << "  " << setw(7) << x[i] * w[i] << "\n";
  }

  cout << "\n";
  cout << "  Total:            " << profit 
    << "  " << mass << "\n";

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests KSUBSET_COLEX_RANK, _SUCCESSOR, _UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  K-subsets of an N set,\n";
  cout << "  using the colexicographic ordering:\n";
  cout << "\n";
  cout << "  KSUBSET_COLEX_RANK ranks,\n";
  cout << "  KSUBSET_COLEX_SUCCESSOR lists,\n";
  cout << "  KSUBSET_COLEX_UNRANK unranks.\n";
  cout << "  KSUBSET_ENUM enumerates,\n";
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of K subsets is " << nksub << "\n";
  cout << "\n";
//
//  List
//
  t = new int[k];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_colex_successor ( k, n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < k; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nksub / 2;

  delete [] t;

  t = ksubset_colex_unrank ( rank, k, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_colex_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests KSUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  K-subsets of an N set,\n";
  cout << "  using the lexicographic ordering:\n";
  cout << "\n";
  cout << "  KSUBSET_ENUM enumerates,\n";
  cout << "  KSUBSET_LEX_RANK ranks,\n";
  cout << "  KSUBSET_LEX_SUCCESSOR lists,\n";
  cout << "  KSUBSET_LEX_UNRANK unranks.\n";
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of K subsets is " << nksub << "\n";
  cout << "\n";
//
//  List
//
  t = new int[k];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_lex_successor ( k, n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < k; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nksub / 2;

  delete [] t;

  t = ksubset_lex_unrank ( rank, k, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_lex_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests KSUBSET_ENUM, _REVDOOR_RANK, _REVDOOR_SUCCESSOR, _REVDOOR_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  K-subsets of an N set,\n";
  cout << "  using the revolving door ordering:\n";
  cout << "\n";
  cout << "  KSUBSET_ENUM enumerates,\n";
  cout << "  KSUBSET_REVDOOR_RANK ranks,\n";
  cout << "  KSUBSET_REVDOOR_SUCCESSOR lists,\n";
  cout << "  KSUBSET_REVDOOR_UNRANK unranks.\n";
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of K subsets is " << nksub << "\n";
  cout << "\n";
//
//  List
//
  t = new int[k];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_revdoor_successor ( k, n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < k; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nksub / 2;

  delete [] t;

  t = ksubset_revdoor_unrank ( rank, k, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_revdoor_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests MARRIAGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int *fiancee;
  int i;
  int n = N;
  int *next;
  int prefer[N*N] = {
    2, 1, 2, 1, 5, 
    5, 2, 3, 3, 3, 
    1, 3, 5, 2, 2, 
    3, 4, 4, 4, 1, 
    4, 5, 1, 5, 4 };
  int rank[N*N] = {
    2, 4, 1, 4, 5, 
    4, 3, 3, 2, 2, 
    5, 5, 4, 1, 3, 
    3, 1, 2, 3, 1, 
    1, 2, 5, 5, 4 };

  cout << "\n";
  cout << "TEST16\n";
  cout << "  MARRIAGE arranges a set of stable marriages\n";
  cout << "  given a set of preferences.\n";

  fiancee = new int[n];
  next = new int[n];

  marriage ( n, prefer, rank, fiancee, next );

  cout << "\n";
  cout << "  Man, Wife's rank, Wife\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << next[i-1]
         << "  " << setw(4) << prefer[i-1+(next[i-1]-1)*n] << "\n";
  }

  cout << "\n";
  cout << "  Woman, Husband's rank, Husband\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << rank[i-1+(fiancee[i-1]-1)*n]
         << "  " << setw(4) << fiancee[i-1] << "\n";
  }

  cout << "\n";
  cout << "  Correct result:\n";
  cout << "\n";
  cout << "  M:W 1  2  3  4  5\n";
  cout << "   1  +  .  .  .  .\n";
  cout << "   2  .  .  .  +  .\n";
  cout << "   3  .  .  .  .  +\n";
  cout << "   4  .  .  +  .  .\n";
  cout << "   5  .  +  .  .  .\n";

  delete [] fiancee;
  delete [] next;

  return;
# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests MOUNTAIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int x;
  int y;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  MOUNTAIN computes mountain numbers.\n";
  cout << "\n";
  cout << "  Y  MXY\n";
  cout << "\n";

  for ( y = 0; y <= n; y++ )
  {
    cout << "  " << setw(2) << y << "   ";

    for ( x = 0; x <= 2 * n; x++ )
    {
      cout << "  " << setw(4) << mountain ( n, x, y );
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests NPART_ENUM, _RSF_LEX_RANK, _RSF_LEX_SUCCESSOR, _RSF_LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int npart = 3;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  n = 12;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  Partitions of N with NPART parts\n";
  cout << "  in reverse standard form:\n";
  cout << "\n";
  cout << "  NPART_ENUM enumerates,\n";
  cout << "  NPART_RSF_LEX_RANK ranks,\n";
  cout << "  NPART_RSF_LEX_SUCCESSOR lists;\n";
  cout << "  NPART_RSF_LEX_UNRANK unranks.\n";
//
//  Enumerate.
//
  npartitions = npart_enum ( n, npart );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  and NPART = " << npart << "\n";
  cout << "  the number of partitions is " << npartitions << "\n";
  cout << "\n";
//
//  List.
//
  t = new int[npart];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    npart_rsf_lex_successor ( n, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = npartitions / 3;

  delete [] t;

  t = npart_rsf_lex_unrank ( rank, n, npart );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( npart, t, " " );
//
//  Rank.
//
  rank = npart_rsf_lex_rank ( n, npart, t );

  i4vec_transpose_print ( npart, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests NPART_RSF_LEX_RANDOM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 12;
  int npart = 3;
  int seed = 123456789;
  int *t;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  Partitions of N with NPART parts\n";
  cout << "  in reverse standard form:\n";
  cout << "\n";
  cout << "  NPART_RSF_LEX_RANDOM produces random examples.\n";

  for ( i = 1; i <= 10; i++ )
  {
    t = npart_rsf_lex_random ( n, npart, &seed );
    i4vec_transpose_print ( npart, t, " " );
    delete [] t;
  }

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests NPART_ENUM and NPART_SF_SUCCESSOR;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 12;
  int npart = 3;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  Partitions of N with NPART parts\n";
  cout << "  in standard form:\n";
  cout << "\n";
  cout << "  NPART_ENUM enumerates,\n";
  cout << "  NPART_SF_LEX_SUCCESSOR lists.\n";
//
//  Enumerate.
//
  npartitions = npart_enum ( n, npart );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  and NPART = " << npart << "\n";
  cout << "  the number of partitions is " << npartitions << "\n";
  cout << "\n";
//
//  List.
//
  t = new int[npart];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    npart_sf_lex_successor ( n, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests NPART_TABLE and PART_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int maxn = 10;
  int maxpart = 5;
  int *p;
  int *p2;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  NPART_TABLE tabulates partitions\n";
  cout << "  of N with NPART parts;\n";
  cout << "  PART_TABLE tabulates partitions of N.\n";

  p = npart_table ( maxn, maxpart );

  p2 = part_table ( maxn );

  cout << "\n";
  cout << "    I P(I)  P(I,0) P(I,1) P(I,2) P(I,3) P(I,4) P(I,5)\n";
  cout << "\n";

  for ( i = 0; i <= maxn; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(4) << p2[i];
    for ( j = 0; j <= maxpart; j++ )
    {
      cout << "  " << setw(4) << p[i+j*(maxn+1)];
    }
    cout << "\n";
  }

  delete [] p;
  delete [] p2;

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests PART_ENUM and PART_SUCCESSOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 8;
  int npart;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  PART_SUCCESSOR produces partitions of N,\n";
  cout << "  PART_ENUM enumerates.\n";
  cout << "\n";
  cout << "  Partitions of N = " << n << "\n";
//
//  Enumerate.
//
  npartitions = part_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of partitions is " << npartitions << "\n";
  cout << "\n";
//
//  List.
//
  t = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    part_successor ( n, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }

  delete [] t;

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests PART_SUCCESSOR and PART_SF_CONJUGATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int i;
  int n = 8;
  int npart;
  int npartb;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  PART_SUCCESSOR produces partitions of N,\n";
  cout << "  PART_SF_CONJUGATE produces the conjugate of a partition.\n";
  cout << "\n";
  cout << "  Partitions of N = " << n << "\n";
//
//  List.
//
  t = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    part_successor ( n, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";

    b = part_sf_conjugate ( n, npart, t, npartb );
    i4vec_transpose_print ( npartb, b, "  Con:" );
    delete [] b;
  }

  delete [] t;

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests PART_SF_MAJORIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 8

  int a[N] = { 2, 2, 2, 1, 1, 0, 0, 0 };
  int b[N] = { 3, 1, 1, 1, 1, 1, 0, 0 };
  int c[N] = { 2, 2, 1, 1, 1, 1, 0, 0 };
  int n = N;
  int nparta = 5;
  int npartb = 6;
  int npartc = 6;
  int result;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  PART_SF_MAJORIZE determines if one partition\n";
  cout << "  majorizes another.\n";
  cout << "\n";
  cout << "  Partitions of N = " << n << "\n";
  cout << "\n";
  i4vec_transpose_print ( nparta, a, "  A: " );
  i4vec_transpose_print ( npartb, b, "  B: " );
  i4vec_transpose_print ( npartc, c, "  C: " );

  result = part_sf_majorize ( n, nparta, a, npartb, b );
  cout << "\n";
  cout << "  A compare B: " << result << "\n";
  result = part_sf_majorize ( n, npartb, b, npartc, c );
  cout << "  B compare C: " << result << "\n";
  result = part_sf_majorize ( n, npartc, c, nparta, a );
  cout << "  C compare A: " << result << "\n";
  result = part_sf_majorize ( n, npartc, c, npartc, c );
  cout << "  C compare C: " << result << "\n";

  return;
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests PARTITION_GREEDY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a1[N] = { 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 };
  int a2[N] = { 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 };
  int i;
  int *indx;
  int n = N;
  int sums[2];

  cout << "\n";
  cout << "TEST25\n";
  cout << "  PARTITION_GREEDY partitions an integer vector into\n";
  cout << "  two subsets with nearly equal sum.\n";
  cout << "\n";

  indx = partition_greedy ( n, a1 );

  cout << "\n";
  cout << "\n";
  cout << "Data set #1 partitioned:\n";
  cout << "\n";
  sums[0] = 0;
  sums[1] = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( indx[i] == 1 ) 
    {
      sums[0] = sums[0] + a1[i];
      cout << "  " << setw(4) << a1[i] << "\n";
    }
    else
    {
      sums[1] = sums[1] + a1[i];
      cout << "  " << "    "
           << "  " << setw(4) << a1[i] << "\n";
    }
  }

  cout << "\n";
  cout << "Sums:\n";
  cout << "\n";
  cout << "  " << setw(4) << sums[0]
       << "  " << setw(4) << sums[1] << "\n";

  delete [] indx;

  indx = partition_greedy ( n, a2 );

  cout << "\n";
  cout << "\n";
  cout << "Data set #2 partitioned:\n";
  cout << "\n";

  sums[0] = 0;
  sums[1] = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( indx[i] == 1 ) 
    {
      sums[0] = sums[0] + a2[i];
      cout << "  " << setw(4) << a2[i] << "\n";
    }
    else
    {
      sums[1] = sums[1] + a2[i];
      cout << "  " << "    "
           << "  " << setw(4) << a2[i] << "\n";
    }
  }

  cout << "\n";
  cout << "Sums:\n";
  cout << "\n";
  cout << "  " << setw(4) << sums[0]
       << "  " << setw(4) << sums[1] << "\n";

  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests PARTN_ENUM, PARTN_SUCCESSOR and PART_SF_CONJUGATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int i;
  int n = 11;
  int nmax;
  int npart;
  int npart2;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  Partitions of N with maximum element NMAX:\n";
  cout << "\n";
  cout << "  PARTN_SUCCESSOR lists;\n";
  cout << "  PARTN_ENUM enumerates.\n";

  nmax = 4;
//
//  Enumerate.
//
  npartitions = partn_enum ( n, nmax );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  and NMAX = " << nmax << "\n";
  cout << "  the number of partitions is " << npartitions << "\n";
  cout << "\n";
//
//  List.
//
  t = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    partn_successor ( n, nmax, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  List conjugates.
//
  cout << "\n";
  cout << "  Repeat, but list RSF conjugated partitions.\n";
  cout << "\n";

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    partn_successor ( n, nmax, npart, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    b = part_sf_conjugate ( n, npart, t, npart2 );

    i4vec_reverse ( npart2, b );

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < npart2; i++ )
    {
      cout << "  " << setw(4) << b[i];
    }
    cout << "\n";
    delete [] b;
  }

  delete [] t;

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests PERM_INV and PERM_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 4;
  int nperm;
  int *p;
  int *q;
  int *r;
  int rank;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  Permutations of the integers:\n";
  cout << "\n";
  cout << "  PERM_INV computes an inverse permutation,\n";
  cout << "  PERM_MUL multiplies two permutations.\n";
//
//  Enumerate.
//
  nperm = perm_enum ( n );
//
//  Unrank.
//
  rank = nperm / 2;

  p = perm_lex_unrank ( rank, n );

  perm_print ( n, p, "  The permutation P:" );
//
//  Invert.
//
  q = perm_inv ( n, p );

  perm_print ( n, q, "  The inverse permutation Q:" );
//
//  Multiply.
//
  r = perm_mul ( n, p, q );

  perm_print ( n, r, "  The product R = P * Q:" );

  delete [] p;
  delete [] q;
  delete [] r;

  return;
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests PERM_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int nperm;
  int *pi;
  int rank;
  int rank_old;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  Permutations of the integers,\n";
  cout << "  using the lexicographic ordering:\n";
  cout << "\n";
  cout << "  PERM_ENUM enumerates,\n";
  cout << "  PERM_LEX_RANK ranks,\n";
  cout << "  PERM_LEX_SUCCESSOR lists,\n";
  cout << "  PERM_LEX_UNRANK unranks.\n";
//
//  Enumerate.
//
  nperm = perm_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of permutations is " << nperm << "\n";
  cout << "\n";
//
//  List
//
  pi = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    perm_lex_successor ( n, pi, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << pi[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nperm / 2;

  delete [] pi;
  pi = perm_lex_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";

  perm_print ( n, pi, " " );
//
//  Rank.
//
  rank = perm_lex_rank ( n, pi );

  perm_print ( n, pi, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] pi;

  return;
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests PERM_TJ_ENUM, _TJ_RANK, _TJ_SUCCESSOR, _TJ_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int nperm;
  int *pi;
  int rank;
  int rank_old;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  Permutations of the integers\n";
  cout << "  using the Trotter-Johnson ordering:\n";
  cout << "\n";
  cout << "  PERM_ENUM enumerates,\n";
  cout << "  PERM_TJ_RANK ranks,\n";
  cout << "  PERM_TJ_SUCCESSOR lists,\n";
  cout << "  PERM_TJ_UNRANK unranks.\n";
//
//  Enumerate.
//
  nperm = perm_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of permutations is " << nperm << "\n";
  cout << "\n";
//
//  List
//
  pi = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    perm_tj_successor ( n, pi, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << pi[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nperm / 2;

  delete [] pi;

  pi = perm_tj_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";

  perm_print ( n, pi, " " );
//
//  Rank.
//
  rank = perm_tj_rank ( n, pi );

  perm_print ( n, pi, "  The element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] pi;

  return;
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests PRUEFER_ENUM, PRUEFER_RANK, PRUEFER_SUCCESSOR, PRUEFER_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int ncode;
  int *p;
  int rank;
  int rank_old;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  Pruefer codes:\n";
  cout << "\n";
  cout << "  PRUEFER_ENUM enumerates,\n";
  cout << "  PRUEFER_RANK ranks,\n";
  cout << "  PRUEFER_SUCCESSOR lists,\n";
  cout << "  PRUEFER_UNRANK unranks.\n";
//
//  Enumerate.
//
  ncode = pruefer_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of Pruefer codes is " << ncode << "\n";
  cout << "\n";
//
//  List
//
  p = new int[n-2];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    pruefer_successor ( n, p, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n - 2; i++ )
    {
      cout << "  " << setw(4) << p[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = ncode / 2;

  delete [] p;

  p = pruefer_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( n - 2, p, " " );
//
//  Rank.
//
  rank = pruefer_rank ( n, p );

  i4vec_transpose_print ( n - 2, p, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] p;

  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests PRUEFER_TO_TREE and TREE_TO_PRUEFER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i4_hi;
  int i4_lo;
  int j;
  int n = 5;
  int *p;
  int pruefer_num;
  int rank;
  int seed = 123456789;
  int *t;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  PRUEFER_TO_TREE converts a Pruefer code to a tree;\n";
  cout << "  TREE_TO_PRUEFER converts a tree to a Pruefer code.\n";

  pruefer_num = pruefer_enum ( n );

  i4_lo = 0;
  i4_hi = pruefer_num - 1;

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Pick a "random" Pruefer code.
//
    rank = i4_uniform ( i4_lo, i4_hi, &seed );

    p = pruefer_unrank ( rank, n );

    cout << "\n";
    cout << "  Random Pruefer code of rank " << rank << "\n";
    i4vec_transpose_print ( n - 2, p, " " );
//
//  Convert the Pruefer code to a tree.
//
    t = pruefer_to_tree_new ( n, p );

    cout << "\n";
    cout << "  Edge list for the corresponding tree:\n";
    cout << "\n";
    for ( j = 0; j < n - 1; j++ )
    {
      cout << "  " << setw(2) << j
           << "  " << setw(4) << t[0+j*2]
           << "  " << setw(4) << t[1+j*2] << "\n";
    }
//
//  Convert the tree to a Pruefer code.
//
    delete [] p;

    p = tree_to_pruefer ( n, t );

    cout << "\n";
    i4vec_transpose_print ( n - 2, p, "  Pruefer code:" );

    delete [] p;
    delete [] t;
  }
  return;
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests QUEENS and BACKTRACK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 8

  int iarray[N];
  int indx;
  int istack[N*N];
  int k;
  int n = N;
  int maxstack = N * N;
  int nstack;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  QUEENS produces nonattacking queens\n";
  cout << "  on a chessboard.\n";
  cout << "  BACKTRACK supervises a backtrack search.\n";
  cout << "\n";

  indx = 0;

  for ( ; ; )
  {
    backtrack ( n, iarray, indx, k, nstack, istack, maxstack );

    if ( indx == 1 )
    {
      i4vec_transpose_print ( n, iarray, " " );
    }
    else if ( indx == 2 )
    {
      queens ( n, iarray, k, nstack, istack, maxstack );
    }
    else
    {
      break;
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests RGF_G_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *d;
  int i;
  int j;
  int m = 6;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  RGF_G_TABLE tabulates generalized restricted\n";
  cout << "  growth functions.\n";
  cout << "\n";

  d = rgf_g_table ( m );

  for ( i = 0; i <= m; i++ )
  {
    for ( j = 0; j <= m - i; j++ )
    {
      cout << "  " << setw(4) << d[i+j*(m+1)];
    }
    cout << "\n";
  }

  delete [] d;

  return;
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests RGF_ENUM, RGF_RANK, RGF_SUCCESSOR, RGF_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *f;
  int i;
  int m = 4;
  int nrgf;
  int rank;
  int rank_old;

  cout << "\n";
  cout << "TEST34\n";
  cout << "  Restricted growth functions:\n";
  cout << "\n";
  cout << "  RGF_ENUM enumerates,\n";
  cout << "  RGF_RANK ranks,\n";
  cout << "  RGF_SUCCESSOR lists;\n";
  cout << "  RGF_UNRANK unranks.\n";
//
//  Enumerate.
//
  nrgf = rgf_enum ( m );

  cout << "\n";
  cout << "  For M = " << m << "\n";
  cout << "  the number of RGF's is " << nrgf << "\n";
  cout << "\n";
//
//  List.
//
  f = new int[m];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    rgf_successor ( m, f, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(4) << f[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nrgf / 2;

  delete [] f;
  f = rgf_unrank ( rank, m );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( m, f, " " );
//
//  Rank.
//
  rank = rgf_rank ( m, f );

  i4vec_transpose_print ( m, f, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] f;

  return;
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests RGF_TO_SETPART and SETPART_TO_RGF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int jlo;
  int *f;
  int *index;
  int m = 8;
  int nsub;
  int rank;
  int *s;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  RGF_TO_SETPART converts a balanced\n";
  cout << "  sequence to a restricted growth function;\n";
  cout << "  SETPART_TO_RGF converts a restricted growth\n";
  cout << "  function to a balanced sequence.\n";
//
//  Choose a "random" RGF.
//
  rank = 7;
  f = rgf_unrank ( rank, m );

  cout << "\n";
  cout << "  Random restricted growth function:\n";
  cout << "\n";
  i4vec_transpose_print ( m, f, " " );
//
//  Convert the RGF to a set partition.
//
  s = new int[m];
  index = new int[m];

  rgf_to_setpart ( m, f, nsub, s, index );

  cout << "\n";
  cout << "  Corresponding set partition\n";
  cout << "\n";
  jlo = 1;
  for ( i = 1; i <= nsub; i++ )
  {
    for ( j = jlo; j <= index[i-1]; j++ )
    {
      cout << "  " << setw(4) << s[j-1];
    }
    cout << "\n";
    jlo = index[i-1] + 1;
  }
//
//  Convert the set partition back to an RGF.
//
  delete [] f;

  f = setpart_to_rgf ( m, nsub, s, index );

  i4vec_transpose_print ( m, f, "  Corresponding RGF:" );

  delete [] f;
  delete [] index;
  delete [] s;

  return;
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests SETPART_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int npart;

  cout << "\n";
  cout << "TEST36\n";
  cout << "  Set partitions:\n";
  cout << "\n";
  cout << "  SETPART_ENUM enumerates.\n";
  cout << "\n";
//
//  Enumerate.
//
  for ( n = 1; n <= 6; n++ )
  {
    npart = setpart_enum ( n );
    cout << "  " << setw(4) << n
         << "  " << setw(4) << npart << "\n";
  }

  return;
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests STIRLING_NUMBERS1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int maxm = 6;
  int maxn = 6;
  int *s;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  STIRLING_NUMBERS1 computes a table of Stirling\n";
  cout << "  numbers of the first kind.\n";

  s = stirling_numbers1 ( maxm, maxn );

  i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of first kind" ); 

  delete [] s;

  return;
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests STIRLING_NUMBERS2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int maxm = 6;
  int maxn = 6;
  int *s;

  cout << "\n";
  cout << "TEST38\n";
  cout << "  STIRLING_NUMBERS2 computes a table of Stirling\n";
  cout << "  numbers of the second kind.\n";

  s = stirling_numbers2 ( maxm, maxn );

  i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of second kind" ); 

  delete [] s;

  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests SUBSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST39\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the colexicographic ordering:\n";
  cout << "\n";
  cout << "  SUBSET_COLEX_RANK ranks,\n";
  cout << "  SUBSET_COLEX_SUCCESSOR lists,\n";
  cout << "  SUBSET_COLEX_UNRANK unranks.\n";
  cout << "  SUBSET_ENUM enumerates.\n";
//
//  Enumerate.
//
  nsub = subset_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  t = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    subset_colex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nsub / 3;

  delete [] t;
  t = subset_colex_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( n, t, " " );
//
//  Rank.
//
  rank = subset_colex_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests SUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  int *t;

  cout << "\n";
  cout << "TEST40\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the lexicographic ordering:\n";
  cout << "\n";
  cout << "  SUBSET_ENUM enumerates,\n";
  cout << "  SUBSET_LEX_RANK ranks,\n";
  cout << "  SUBSET_LEX_SUCCESSOR lists,\n";
  cout << "  SUBSET_LEX_UNRANK unranks.\n";
//
//  Enumerate.
//
  nsub = subset_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  t = new int[n];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    subset_lex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }

    cout << "  " << setw(4) << rank;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << t[i];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = nsub / 3;

  delete [] t;

  t = subset_lex_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";
  cout << "\n";
  i4vec_transpose_print ( n, t, " " );
//
//  Rank.
//
  rank = subset_lex_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests SUBSETSUM_SWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N] = { 12, 8, 11, 30, 8, 3, 7 };
  int i;
  int index[N];
  int n = N;
  int sum_achieved;
  int sum_desired = 17;

  cout << "\n";
  cout << "TEST41\n";
  cout << "  SUBSETSUM_SWAP seeks a solution of the subset\n";
  cout << "  sum problem using pair swapping.\n";
  cout << "\n";
  cout << "  The desired sum is " << sum_desired << "\n";

  sum_achieved = subsetsum_swap ( n, a, sum_desired, index );

  cout << "\n";
  cout << "    A(I), INDEX(I)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(5) << a[i]
         << "  " << setw(5) << index[i] << "\n";
  }

  cout << "\n";
  cout << "  The achieved sum is " << sum_achieved << "\n";

  return;
# undef N
}
//****************************************************************************80

void test42 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST42 tests TREE_ENUM, TREE_RANK, TREE_SUCCESSOR, TREE_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 4;
  int rank;
  int rank_old;
  int *t;
  int tree_num;

  cout << "\n";
  cout << "TEST42\n";
  cout << "  Trees:\n";
  cout << "\n";
  cout << "  TREE_ENUM enumerates,\n";
  cout << "  TREE_RANK ranks,\n";
  cout << "  TREE_SUCCESSOR lists,\n";
  cout << "  TREE_UNRANK unranks.\n";
//
//  Enumerate.
//
  tree_num = tree_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of trees is " << tree_num << "\n";
  cout << "\n";
//
//  List
//
  t = new int[2*(n-1)];

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    tree_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    cout << "  " << setw(4) << rank;
    for ( j = 0; j < n - 1; j++ )
    {
      cout << "  " << setw(4) << t[0+j*2];
    }
    cout << "\n";
    cout << "  " << "    ";
    for ( j = 0; j < n - 1; j++ )
    {
      cout << "  " << setw(4) << t[1+j*2];
    }
    cout << "\n";
  }
//
//  Unrank.
//
  rank = tree_num / 2;

  delete [] t;

  t = tree_unrank ( rank, n );

  cout << "\n";
  cout << "  The element of rank " << rank << "\n";

  i4mat_print ( 2, n - 1, t, " " );
//
//  Rank.
//
  rank = tree_rank ( n, t );

  i4mat_print ( 2, n - 1, t, "  Element to be ranked:" );
  cout << "\n";
  cout << "  Rank is computed as " << rank << "\n";

  delete [] t;

  return;
}
