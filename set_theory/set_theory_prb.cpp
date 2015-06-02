# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "set_theory.hpp"

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
//    MAIN is the main program for SET_THEORY_PRB.
//
//  Discussion:
//
//    SET_THEORY_PRB tests the SET_THEORY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SET_THEORY_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SET_THEORY library.\n";
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
  cout << "SET_THEORY_PRB\n";
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
//    TEST01 tests the B4SET routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int a_num;
  int *a_numeric;
  int b;
  int b_num = 16;
  int b_numeric[16] = {
    3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48 };
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int i;
  int n = 32;
  int u;
  int *u_numeric;
  int w;
  int w_num = 5;
  int w_numeric[5] = { 1, 11, 21, 31, 41 };
  int x;
  int y;
  int y_num = 4;
  int y_numeric[4] = { 4, 5, 6, 7 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the set theory functions\n";
  cout << "  with the B4SET representation of a set.\n";
//
//  Define the universal set.
//
  u_numeric = new int[n];

  for ( i = 0; i < n; i++ )
  {
    u_numeric[i] = i + 1;
  }
  u = i4vec_to_b4set ( n, u_numeric, n );
//
//  Define the set A by a numeric property.
//
  a_numeric = new int[n];
  a_num = 0;
  for ( i = 1; i <= n; i++ )
  {
    if ( ( i % 5 ) == 0 )
    {
      a_numeric[a_num] = i;
      a_num = a_num + 1;
    }
  }
  a = i4vec_to_b4set ( a_num, a_numeric, n );
  b4set_transpose_print ( n, a, "  A: " );
//
//  Define the set by starting with a numeric list of entries.
//
  b = i4vec_to_b4set ( b_num, b_numeric, n );
  b4set_transpose_print ( n, b, "  B: " );
//
//  C is the complement of B (with respect to the universal set).
//
  c = b4set_complement ( n, b );
  b4set_transpose_print ( n, c, "  C = ~ B:" );
//
//  D is the intersection of A and B.
//
  d = b4set_intersect ( n, a, b );
  b4set_transpose_print ( n, d, "  D = A intersect B:" );
//
//  E is the intersection of A and B.
//
  e = b4set_union ( n, a, b );
  b4set_transpose_print ( n, e, "  E = A union B:" );
//
//  F is the symmetric difference of A and B.
//
  f = b4set_xor ( n, a, b );
  b4set_transpose_print ( n, f, "  F = A xor B:" );
//
//  G is the complement of B with respect to A.
//  H is the complement of A with respect to B.
//
  g = b4set_complement_relative ( n, a, b );
  b4set_transpose_print ( n, g, "  G = A ~ B:" );

  h = b4set_complement_relative ( n, b, a );
  b4set_transpose_print ( n, h, "  H = B ~ A:" );
//
//  B4SET_IS_MEMBER checks if an element is in a set.
//
  cout << "\n";
  cout << "  B4SET_IS_MEMBER ( i, A ) reports whether i is a member of A\n";
  cout << "\n";

  for ( i = 10; i <= 20; i++ )
  {
    if ( b4set_is_member ( n, i, a ) )
    {
      cout << i << " is a member of A.\n";
    }
    else
    {
      cout << i << " is not a member of A.\n";
    }
  }
//
//  B4SET_IS_SUBSET checks whether a set is a subset.
//
  cout << "\n";
  cout << "  B4SET_IS_SUBSET ( D, A ) reports whether D is a subset of A\n";
  cout << "\n";

  d = b4set_intersect ( n, a, b );

  if ( b4set_is_subset ( n, d, a ) )
  {
    cout << "  ( A intersect B ) is a subset of A.\n";
  }
  else
  {
    cout << "  ( A intersect B)  is not a subset of A.\n";
  }
//
//  B4SET_INSERT adds an item to a set.
//
  w = i4vec_to_b4set ( w_num, w_numeric, n );
  b4set_transpose_print ( n, w, "  W: " );

  x = 6;
  w = b4set_insert ( n, x, w );
  b4set_transpose_print ( n, w, "  W := W + 6:" );

  x = 31;
  w = b4set_delete ( n, x, w );
  b4set_transpose_print ( n, w, "  W := W - 31:" );

  y = i4vec_to_b4set ( y_num, y_numeric, n );
  w = b4set_union ( n, w, y );

  b4set_transpose_print ( n, w, "  W := W union [ 4, 5, 6, 7 ]:" );

  delete [] a_numeric;
  delete [] u_numeric;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests B4SET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  string s;
  int t;
  string title;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the colexicographic ordering\n";
  cout << "  with the B4SET representation of a set.\n";
  cout << "\n";
  cout << "  B4SET_COLEX_RANK ranks,\n";
  cout << "  B4SET_COLEX_SUCCESSOR lists,\n";
  cout << "  B4SET_COLEX_UNRANK unranks.\n";
  cout << "  B4SET_ENUM enumerates.\n";
//
//  Enumerate.
//
  nsub = b4set_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    b4set_colex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    s = i4_to_s ( rank );
    title = "Rank: " + s;
    b4set_transpose_print ( n, t, title );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  t = b4set_colex_unrank ( rank, n );

  s = i4_to_s ( rank );
  title = "The element of rank " + s;
  b4set_transpose_print ( n, t, title );
//
//  Rank.
//
  rank = b4set_colex_rank ( n, t );

  cout << "\n";
  cout << "  The rank of this element is computed as " << rank << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests B4SET_LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  string s;
  int t;
  string title;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the lexicographic ordering,\n";
  cout << "  with the B4SET representation of a set.\n";
  cout << "\n";
  cout << "  B4SET_LEX_RANK ranks,\n";
  cout << "  B4SET_LEX_SUCCESSOR lists,\n";
  cout << "  B4SET_LEX_UNRANK unranks.\n";
  cout << "  B4SET_ENUM enumerates.\n";
//
//  Enumerate.
//
  nsub = b4set_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    b4set_lex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    s = i4_to_s ( rank );
    title = "Rank: " + s;
    b4set_transpose_print ( n, t, title );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  t = b4set_lex_unrank ( rank, n );

  s = i4_to_s ( rank );
  title = "The element of rank " + s;
  b4set_transpose_print ( n, t, title );
//
//  Rank.
//
  rank = b4set_lex_rank ( n, t );

  cout << "\n";
  cout << "  The rank of this element is computed as " << rank << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests the LSET routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  bool *a;
  bool *b;
  int b_num = 16;
  int b_numeric[16] = { 
    3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48 };
  bool *c;
  bool *d;
  bool *e;
  bool *f;
  bool *g;
  bool *h;
  int i;
  int n = 50;
  bool *u;
  int *u_numeric;
  bool *w;
  int w_num = 5;
  int w_numeric[5] = { 1, 11, 21, 31, 41 };
  int x;
  bool *y;
  int y_num = 4;
  int y_numeric[4] = { 16, 26, 36, 46 };

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test the set theory functions\n";
  cout << "  with the LSET representation of a set.\n";
//
//  Define the universal set.
//
  u_numeric = new int[n];
  u = new bool[n];

  for ( i = 0; i < n; i++ )
  {
    u_numeric[i] = i + 1;
    u[i] = true;
  }
//
//  Define the set A by a numeric property.
//
  a = new bool[n];
  for ( i = 0; i < n; i++ )
  {
    a[i] = ( ( u_numeric[i] % 5 ) == 0 );
  }
  lset_transpose_print ( n, a, "  A: " );
//
//  Define the set by starting with a numeric list of entries.
//
  b = i4vec_to_lset ( b_num, b_numeric, n );
  lset_transpose_print ( n, b, "  B: " );
//
//  C is the complement of B (with respect to the universal set).
//
  c = lset_complement ( n, b );
  lset_transpose_print ( n, c, "  C = ~ B:" );
//
//  D is the intersection of A and B.
//
  d = lset_intersect ( n, a, b );
  lset_transpose_print ( n, d, "  D = A intersect B:" );
//
//  E is the intersection of A and B.
//
  e = lset_union ( n, a, b );
  lset_transpose_print ( n, e, "  E = A union B:" );
//
//  F is the symmetric difference of A and B.
//
  f = lset_xor ( n, a, b );
  lset_transpose_print ( n, f, "  F = A xor B:" );
//
//  G is the complement of B with respect to A.
//  H is the complement of A with respect to B.
//
  g = lset_complement_relative ( n, a, b );
  lset_transpose_print ( n, g, "  G = A ~ B:" );

  h = lset_complement_relative ( n, b, a );
  lset_transpose_print ( n, h, "  H = B ~ A:" );
//
//  LSET_IS_MEMBER checks if an element is in a set.
//
  cout << "\n";
  cout << "  LSET_IS_MEMBER ( i, A ) reports whether i is a member of A\n";
  cout << "\n";

  for ( i = 10; i <= 20; i++ )
  {
    if ( lset_is_member ( n, i, a ) )
    {
      cout << i << " is a member of A.\n";
    }
    else
    {
      cout << i << " is not a member of A.\n";
    }
  }
//
//  LSET_IS_SUBSET checks whether a set is a subset.
//
  cout << "\n";
  cout << "  LSET_IS_SUBSET ( D, A ) reports whether D is a subset of A\n";
  cout << "\n";

  delete [] d;
  d = lset_intersect ( n, a, b );

  if ( lset_is_subset ( n, d, a ) )
  {
    cout << "  ( A intersect B ) is a subset of A.\n";
  }
  else
  {
    cout << "  ( A intersect B)  is not a subset of A.\n";
  }
//
//  LSET_INSERT adds an item to a set.
//
  w = i4vec_to_lset ( w_num, w_numeric, n );
  lset_transpose_print ( n, w, "  W:" );

  x = 6;
  lset_insert ( n, x, w );
  lset_transpose_print ( n, w, "  W := W + 6:" );

  x = 31;
  lset_delete ( n, x, w );
  lset_transpose_print ( n, w, "  W := W - 31:" );
  delete [] w;

  y = i4vec_to_lset ( y_num, y_numeric, n );
  w = lset_union ( n, w, y );

  lset_transpose_print ( n, w, "  W := W union [16, 26, 36, 46]:" );

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
  delete [] e;
  delete [] f;
  delete [] g;
  delete [] h;
  delete [] u;
  delete [] u_numeric;
  delete [] w;
  delete [] y;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests LSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2011
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
  string s;
  bool *t;
  string title;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the colexicographic ordering\n";
  cout << "  with the LSET representation of a set.\n";
  cout << "\n";
  cout << "  LSET_COLEX_RANK ranks,\n";
  cout << "  LSET_COLEX_SUCCESSOR lists,\n";
  cout << "  LSET_COLEX_UNRANK unranks.\n";
  cout << "  LSET_ENUM enumerates.\n";

  t = new bool[n];
//
//  Enumerate.
//
  nsub = lset_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  rank = -1;

  while ( true )
  {
    rank_old = rank;

    lset_colex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    s = i4_to_s ( rank );
    title = "Rank: " + s;
    lset_transpose_print ( n, t, title );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  delete [] t;

  t = lset_colex_unrank ( rank, n );

  s = i4_to_s ( rank );
  title = "  The element of rank " + s;
  lset_transpose_print ( n, t, title );
//
//  Rank.
//
  rank = lset_colex_rank ( n, t );

  cout << "\n";
  cout << "  The rank of this element is computed as " << rank << "\n";

  delete [] t;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests LSET_LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  string s;
  bool *t;
  string title;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  All subsets of a set,\n";
  cout << "  using the lexicographic ordering,\n";
  cout << "  with the LSET representation of a set.\n";
  cout << "\n";
  cout << "  LSET_LEX_RANK ranks,\n";
  cout << "  LSET_LEX_SUCCESSOR lists,\n";
  cout << "  LSET_LEX_UNRANK unranks.'\n";
  cout << "  LSET_ENUM enumerates.\n";

  t = new bool[n];
//
//  Enumerate.
//
  nsub = lset_enum ( n );

  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  the number of subsets is " << nsub << "\n";
  cout << "\n";
//
//  List
//
  rank = -1;

  while ( 1 )
  {
    rank_old = rank;

    lset_lex_successor ( n, t, rank );

    if ( rank <= rank_old )
    {
      break;
    }
    s = i4_to_s ( rank );
    title = "Rank: " + s;
    lset_transpose_print ( n, t, title );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  delete [] t;
  t = lset_lex_unrank ( rank, n );

  s = i4_to_s ( rank );
  title = "  The element of rank " + s;
  lset_transpose_print ( n, t, title );
//
//  Rank.
//
  rank = lset_lex_rank ( n, t );

  cout << "\n";
  cout << "  The rank of this element is computed as " << rank << "\n";

  delete [] t;

  return;
}
