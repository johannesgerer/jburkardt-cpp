# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "subset.hpp"

int main ( );

void test000 ( );
void test001 ( );
void test002 ( );
void test003 ( );
void test9001 ( );
void test9002 ( );
void test9003 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );

void test010 ( );
void test011 ( );
void test012 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test016 ( );
void test017 ( );
void test0174 ( );
void test0175 ( );
void test018 ( );
void test019 ( );

void test020 ( );
void test021 ( );
void test022 ( );
void test023 ( );
void test024 ( );
void test025 ( );
void test026 ( );
void test027 ( );
void test028 ( );
void test029 ( );

void test021a ( );
void test022a ( );
void test023a ( );
void test024a ( );
void test025a ( );
void test026a ( );
void test026b ( );
void test026c ( );
void test026d ( );
void test027a ( );
void test028a ( );
void test029a ( );
void test0295 ( );

void test0304 ( );
void test0305 ( );
void test031 ( );
void test032 ( );
void test0321 ( );
void test0322 ( );
void test03225 ( );
void test0323 ( );
void test0324 ( );
void test0325 ( );
void test0327 ( );
void test058 ( );
void test059 ( );
void test060 ( );
void test061 ( );
void test0615 ( );
void test062 ( );
void test06225 ( );
void test033 ( );
void test034 ( );
void test0625 ( );
void test035 ( );
void test0627 ( );
void test0364 ( );
void test036 ( );
void test037 ( );
void test038 ( );
void test039 ( );

void test040 ( );
void test041 ( );
void test042 ( );
void test043 ( );
void test044 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test048 ( );
void test049 ( );

void test050 ( );
void test051 ( );
void test052 ( );
void test053 ( );
void test054 ( );
void test055 ( );
void test056 ( );
void test057 ( );
void test063 ( );
void test064 ( );
void test065 ( );
void test066 ( );
void test067 ( );
void test0675 ( );
void test068 ( );
void test0683 ( );
void test0685 ( );
void test0686 ( );
void test0687 ( );
void test0688 ( );
void test0689 ( );
void test06895 ( );
void test069 ( );

void test070 ( );
void test071 ( );
void test072 ( );
void test073 ( );
void test074 ( );
void test075 ( );
void test076 ( );
void test077 ( );
void test0771 ( );
void test07715 ( );
void test0772 ( );
void test0773 ( );
void test078 ( );
void test079 ( );
void test0795 ( );

void test080 ( );
void test081 ( );
void test0813 ( );
void test0815 ( );
void test082 ( );
void test083 ( );
void test0835 ( );
void test084 ( );
void test085 ( );
void test086 ( );
void test087 ( );
void test088 ( );
void test089 ( );

void test090 ( );
void test091 ( );
void test092 ( );
void test093 ( );
void test094 ( );
void test095 ( );
void test0955 ( );
void test096 ( );
void test097 ( );
void test098 ( );
void test099 ( );

void test100 ( );
void test101 ( );
void test102 ( );
void test103 ( );
void test104 ( );
void test105 ( );
void test106 ( );
void test107 ( );
void test108 ( );
void test1085 ( );
void test109 ( );
void test110 ( );

void test111 ( );
void test112 ( );
void test113 ( );
void test114 ( );
void test115 ( );
void test030 ( );
void test1245 ( );
void test1155 ( );
void test116 ( );
void test1163 ( );
void test1165 ( );
void test117 ( );
void test118 ( );
void test119 ( );

void test120 ( );
void test121 ( );
void test1215 ( );
void test122 ( );
void test123 ( );
void test124 ( );
void test125 ( );
void test126 ( );
void test127 ( );
void test128 ( );
void test129 ( );

void test130 ( );
void test131 ( );
void test132 ( );
void test133 ( );
void test134 ( );
void test135 ( );
void test136 ( );
void test137 ( );
void test138 ( );
void test139 ( );
void test1395 ( );

void test140 ( );
void test141 ( );
void test142 ( );
void test143 ( );
void test1435 ( );
void test144 ( );
void test145 ( );
void test146 ( );
void test147 ( );
void test1475 ( );
void test1476 ( );
void test1477 ( );
void test1478 ( );
void test148 ( );
void test149 ( );

void test150 ( );
void test151 ( );
void test152 ( );
void test153 ( );
void test1531 ( );
void test0626 ( );
void test1535 ( );
void test1536 ( );
void test1537 ( );
void test154 ( );
void test155 ( );
void test156 ( );
void test1565 ( );
void test1566 ( );
void test1567 ( );
void test1568 ( );
void test1569 ( );
void test15695 ( );
void test15696 ( );
void test15698 ( );
void test157 ( );
void test158 ( );
void test159 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUBSET_PRB.
//
//  Discussion:
//
//    SUBSET_PRB calls the SUBSET test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "SUBSET_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SUBSET library.\n";

  test000 ( );
  test001 ( );
  test002 ( );
  test003 ( );
  test9001 ( );
  test9002 ( );
  test9003 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );

  test010 ( );
  test011 ( );
  test012 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test016 ( );
  test017 ( );
  test0174 ( );
  test0175 ( );
  test018 ( );
  test019 ( );

  test020 ( );
  test021 ( );
  test022 ( );
  test023 ( );
  test024 ( );
  test025 ( );
  test026 ( );
  test027 ( );
  test028 ( );
  test029 ( );

  test021a ( );
  test022a ( );
  test023a ( );
  test024a ( );
  test025a ( );
  test026a ( );
  test026b ( );
  test026c ( );
  test026d ( );
  test027a ( );
  test028a ( );
  test029a ( );
  test0295 ( );

  test0304 ( );
  test0305 ( );
  test031 ( );
  test032 ( );
  test0321 ( );
  test0322 ( );
  test03225 ( );
  test0323 ( );
  test0324 ( );
  test0325 ( );
  test0327 ( );
  test058 ( );
  test059 ( );
  test060 ( );
  test061 ( );
  test0615 ( );
  test062 ( );
  test06225 ( );
  test033 ( );
  test034 ( );
  test0625 ( );
  test035 ( );
  test0627 ( );
  test0364 ( );
  test036 ( );
  test037 ( );
  cout << "\n";
  cout << "SUBSET_PRB:\n";
  cout << "  TEST038 is imperfect because I4MAT_ROWCOLSUM2 is not ready.\n";
  test038 ( );
  test039 ( );

  test040 ( );
  test041 ( );
  test042 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test051 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test056 ( );
  test057 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test067 ( );
  test0675 ( );
  test068 ( );
  test0683 ( );
  test0685 ( );
  test0686 ( );
  test0687 ( );
  test0688 ( );
  test0689 ( );
  test06895 ( );
  test069 ( );

  test070 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test075 ( );
  test076 ( );
  test077 ( );
  test0771 ( );
  test07715 ( );
  test0772 ( );
  test0773 ( );
  test078 ( );
  test079 ( );
  test0795 ( );

  test080 ( );
  test081 ( );
  test0813 ( );
  test0815 ( );
  test083 ( );
  test0835 ( );
  test084 ( );
  test085 ( );
  test086 ( );
  test087 ( );
  test088 ( );
  test089 ( );

  test090 ( );
  test091 ( );
  test092 ( );
  test093 ( );
  test094 ( );
  test095 ( );
  test0955 ( );
  test096 ( );
  test097 ( );
  test098 ( );
  test099 ( );

  test100 ( );
  test101 ( );
  test102 ( );
  test103 ( );
  test104 ( );
  test105 ( );
  test106 ( );
  test107 ( );
  test108 ( );
  test1085 ( );
  test109 ( );
  test110 ( );

  test111 ( );
  test112 ( );
  test113 ( );
  test114 ( );
  test115 ( );
  test030 ( );
  test1245 ( );
  test1155 ( );
  test116 ( );
  test1163 ( );
  test1165 ( );
  test117 ( );
  test118 ( );
  test119 ( );

  test120 ( );
  test121 ( );
  test1215 ( );
  test122 ( );
  test123 ( );
  test124 ( );
  test125 ( );
  test126 ( );
  test127 ( );
  test128 ( );
  test129 ( );

  test130 ( );
  test131 ( );
  test132 ( );
  test133 ( );
  test134 ( );
  test135 ( );
  test136 ( );
  test137 ( );
  test138 ( );
  test139 ( );
  test1395 ( );

  test140 ( );
  test141 ( );
  test142 ( );
  test143 ( );
  test1435 ( );
  test144 ( );
  test145 ( );
  test146 ( );
  test147 ( );
  test1475 ( );
  test1476 ( );
  test1477 ( );
  test1478 ( );
  test148 ( );
  test149 ( );

  test150 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test1531 ( );
  test0626 ( );
  test1535 ( );
  test1536 ( );
  test1537 ( );
  test154 ( );
  test155 ( );
  test156 ( );
  test1565 ( );
  test1566 ( );
  test1567 ( );
  test1568 ( );
  test1569 ( );
  test15695 ( );
  test15696 ( );
  test15698 ( );
  test157 ( );
  test158 ( );
  test159 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SUBSET_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  exit ( 0 );
}
//****************************************************************************80

void test000 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST000 tests RANDOM_INITIALIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  unsigned long seed;

  cout << "\n";
  cout << "TEST000\n";
  cout << "  Call RANDOM_INITIALIZE to initialize the\n";
  cout << "  RANDOM random number generator.\n";

  seed = 0;
  seed = random_initialize ( seed );

  cout << "\n";
  cout << "  RANDOM_INITIALIZE returns SEED = " << seed << "\n";

  return;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests ASM_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  ASM_ENUM returns the number of alternating sign\n";
  cout << "  matrices of a given order.\n";

  cout << "\n";

  for ( n = 0; n <= 7; n++ )
  {
    cout << setw(4) << n << "  "
         << setw(6) << asm_enum ( n ) << "\n";
  }

  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests ASM_TRIANGLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 7

  int a[N_MAX+1];
  int i;
  int n;

  cout << "\n";
  cout << "TEST002\n";
  cout << "  ASM_TRIANGLE returns a row of the alternating sign\n";
  cout << "  matrix triangle.\n";
  cout << "\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    asm_triangle ( n, a );
    cout << setw(4) << n << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(8) << a[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests BELL and BELL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  BELL computes Bell numbers.\n";
  cout << "  BELL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    bell ( n, c2 );

    cout                 << "  "
         << setw(4) << n << "  "
         << setw(8) << c << "  "
         << setw(8) << c2[n] << "\n";

    delete [] c2;
  }

  return;
}
//****************************************************************************80

void test9001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST9001 tests BVEC_ADD and BVEC_SUB;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int bvec4[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST9001\n";
  cout << "  BVEC_ADD adds binary vectors representing integers;\n";
  cout << "  BVEC_SUB subtracts binary vectors representing integers;\n";
  cout << "\n";
  cout << "        I        J        K = I + J    L = I - J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( -100, 100, seed );
    j = i4_uniform ( -100, 100, seed );

    cout << "\n";
    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";

    k = i + j;
    l = i - j;

    cout << "  Directly:         "
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );

    bvec_add ( N, bvec1, bvec2, bvec3 );
    k = bvec_to_i4 ( N, bvec3 );

    bvec_sub ( N, bvec1, bvec2, bvec4 );
    l = bvec_to_i4 ( N, bvec4 );

    cout << "  BVEC_ADD, BVEC_SUB"
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test9002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST9002 tests BVEC_COMPLEMENT2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST9002\n";
  cout << "  BVEC_COMPLEMENT2 returns the two's complement\n";
  cout << "  of a (signed) binary vector;\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( -100, 100, seed );

    i4_to_bvec ( i, N, bvec1 );

    bvec_complement2 ( N, bvec1, bvec2 );

    j = bvec_to_i4 ( N, bvec2 );

    cout << "\n";
    cout << "  I = " << "  " << i << "\n";
    cout << "  J = " << "  " << j << "\n";
    bvec_print ( N, bvec1, " " );
    bvec_print ( N, bvec2, " " );

  }

  return;
# undef N
}
//****************************************************************************80

void test9003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST9003 tests BVEC_MUL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 15

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int i;
  int j;
  int k;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST9003\n";
  cout << "  BVEC_MUL multiplies binary vectors\n";
  cout << "  representing integers;\n";
  cout << "\n";
  cout << "        I        J        K = I * J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( -100, 100, seed );
    j = i4_uniform ( -100, 100, seed );

    cout << "\n";
    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";

    k = i * j;

    cout << "  Directly:         "
         << "  " << setw(8) << k<< "\n";

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );

    bvec_mul ( N, bvec1, bvec2, bvec3 );
    k = bvec_to_i4 ( N, bvec3 );

    cout << "  BVEC_MUL          "
         << "  " << setw(8) << k << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests CATALAN and CATALAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  CATALAN computes Catalan numbers.\n";
  cout << "  CATALAN_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    catalan ( n, c2 );

    cout                     << "  "
         << setw(4) << n     << "  "
         << setw(8) << c     << "  "
         << setw(8) << c2[n] << "\n";

    delete [] c2;

  }

  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests CATALAN_ROW_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  bool next;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  CATALAN_ROW_NEXT computes a row of the Catalan triangle.\n";
  cout << "\n";
  cout << "  First, compute row 7:\n";

  next = false;
  n = 7;
  catalan_row_next ( next, n, c );

  cout << setw(4) << n << "  ";
  for ( i = 0; i <= n; i++ )
  {
    cout << setw(8) << c[i] << "  ";
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now compute rows consecutively, one at a time:\n";
  cout << "\n";

  next = false;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = true;

    cout << setw(4) << n << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(6) << c[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests CFRAC_TO_RAT and RAT_TO_CFRAC.
//
//  Discussion:
//
//    Compute the continued fraction form of 4096/15625.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10

  int a[M];
  int bot = 15625;
  bool error;
  int i;
  int n;
  int p[M];
  int q[M];
  int top = 4096;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  RAT_TO_CFRAC fraction => continued fraction,\n";
  cout << "  CFRAC_TO_RAT continued fraction => fraction.\n";
  cout << "\n";
  cout << "  Regular fraction is " << top << "/" << bot << "\n";

  rat_to_cfrac ( top, bot, M, n, a, error );
 
  i4vec1_print ( n, a, "  Continued fraction coefficients:" );

  cfrac_to_rat ( n, a, p, q );
 
  cout << "\n";
  cout << "  The continued fraction convergents.\n";
  cout << "  The last row contains the value of the continued\n";
  cout << "  fraction, written as a common fraction.\n";
  cout << "\n";
  cout << "  I, P(I), Q(I), P(I)/Q(I)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << setw(3) <<   i  << "  "
         << setw(6) << p[i] << "  "
         << setw(6) << q[i] << "  "
         << setw(14) << ( double ) p[i] / ( double ) q[i] << "\n";
  }
 
  return;
# undef M
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests CFRAC_TO_RFRAC and RFRAC_TO_CFRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAXM 10

  bool error;
  double g[2*MAXM];
  double h[2*MAXM];
  int i;
  int m;
  double p[MAXM];
  double q[MAXM+1];

  m = 3;

  p[0] = 1.0;
  p[1] = 1.0;
  p[2] = 2.0;

  q[0] = 1.0;
  q[1] = 3.0;
  q[2] = 1.0;
  q[3] = 1.0;

  cout << "\n";
  cout << "TEST008\n";
  cout << "  CFRAC_TO_RFRAC: continued fraction to ratio;\n";
  cout << "  RFRAC_TO_CFRAC: ratio to continued fration.\n";

  cout << "\n";
  cout << "  Rational polynomial fraction coefficients:\n";
  cout << "\n";

  cout << "  P:  ";
  for ( i = 0; i < m; i++ )
  {
    cout << setw(12) << p[i];
  }
  cout << "\n";

  cout << "  Q:  ";
  for ( i = 0; i < m+1; i++ )
  {
    cout << setw(12) << q[i];
  }
  cout << "\n";

  rfrac_to_cfrac ( m, p, q, h, error );
 
  r8vec_print ( 2*m, h, "  Continued fraction coefficients:" );

  for ( i = 0; i < 2 * m; i++ )
  {
    g[i] = 1.0;
  }

  cfrac_to_rfrac ( 2*m, g, h, p, q );
 
  cout << "\n";
  cout << "  Recovered rational polynomial:\n";
  cout << "\n";
  
  cout << "  P:  ";
  for ( i = 0; i < m; i++ )
  {
    cout << setw(12) << p[i];
  }
  cout << "\n";

  cout << "  Q:  ";
  for ( i = 0; i < m+1; i++ )
  {
    cout << setw(12) << q[i];
  }
  cout << "\n";

  return;
# undef MAXM
}
//****************************************************************************80

void test009 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests CHANGE_GREEDY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define COIN_NUM 6

  int change[100];
  int change_num;
  int coin_value[COIN_NUM] = { 1, 5, 10, 25, 50, 100 };
  int i;
  int total;
  int total2;

  cout << "\n";
  cout << "TEST009\n";
  cout << "  CHANGE_GREEDY makes change using the biggest\n";
  cout << "  coins first.\n";

  total = 73;

  cout << "\n";
  cout << "  The total for which change is to be made: " << total << "\n";
  cout << "\n";
  cout << "  The available coins are:\n";
  cout << "\n";

  for ( i = 0; i < COIN_NUM; i++ )
  {
    cout << setw(6) << coin_value[i] << "\n";
  }

  change_greedy ( total, COIN_NUM, coin_value, change_num, change );

  cout << "\n";
  cout << "  The number of coins in change is: " << change_num << "\n";
  cout << "\n";

  cout << "        ";  
  for ( i = 0; i < change_num; i++ )
  {
    cout << setw(3) << change[i] << "  ";
  }
  cout << "\n";

  total2 = 0;
  for ( i = 0; i < change_num; i++ )
  {
    total2 = total2 + coin_value[change[i]];
  }

  cout << setw(6) << total2 << "  ";
  for ( i = 0; i < change_num; i++ )
  {
    cout << setw(3) << coin_value[change[i]] << "  ";
  }
  cout << "\n";

  return;
# undef COIN_NUM
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests CHANGE_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define COIN_NUM 6

  int change[100];
  int change_num;
  int coin_value[COIN_NUM] = { 1, 5, 10, 25, 50, 100 };
  bool done;
  int i;
  int j;
  int total;

  cout << "\n";
  cout << "TEST010\n";
  cout << "  CHANGE_NEXT displays the next possible way to make\n";
  cout << "  change for a given total\n";

  total = 50;

  cout << "\n";
  cout << "  The total for which change is to be made: " << total << "\n";
  cout << "\n";

  cout << "\n";
  cout << "  The available coins are:\n";
  cout << "\n";
  for ( i = 0; i < COIN_NUM; i++ )
  {
    cout << setw(6) << coin_value[i] << "\n";
  }

  done = true;
  i = 0;

  for ( ; ; )
  {
    change_next ( total, COIN_NUM, coin_value, change_num, change, done );

    if ( done || 9 < i )
    {
      break;
    }

    i = i + 1;
    cout << "\n";
    cout << setw(3) << i << "\n";
    for ( j = 0; j < change_num; j++ )
    {
      cout << setw(3) << coin_value[change[j]] << "  ";
    }
    cout << "\n";
  }

  return;
# undef COIN_NUM
}
//****************************************************************************80

void test011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests COMB_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  bool done;
  int i;
  int j;
  int k;
  int n = N;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  COMB_NEXT produces combinations.\n";

  for ( k = 1; k <= n; k++ )
  {
    cout << "\n";
    cout << "  Combinations of size K = " << k << "\n";
    cout << "\n";

    done = true;

    for ( ; ; )
    {
      comb_next ( n, k, a, done );

      if ( done )
      {
        break;
      }

      for ( i = 0; i < k; i++ )
      {
        cout << "  " << setw(4) << a[i];
      }
      cout << "\n";
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests COMB_ROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int c[N+1];
  int i;
  int j;
  bool next;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  COMB_ROW computes a row of the Pascal triangle.\n";
  cout << "\n";

  next = false;

  for ( i = 0; i <= N; i++ )
  {
    comb_row ( next, i, c );
    next = true;
    cout                 << "  "
         << setw(2) << i << "  ";
    for ( j = 0; j <= i; j++ )
    {
      cout << setw(5) << c[j];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test013 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests COMB_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int cnk;
  int i;
  int m = 10;
  int rank;

  cnk = i4_choose ( m, N );

  cout << "\n";
  cout << "TEST013\n";
  cout << "  COMB_UNRANK returns a combination of N things\n";
  cout << "  out of M, given the lexicographic rank.\n";
  cout << "\n";
  cout << "  The total set size is M = " << m << "\n";
  cout << "  The subset size is N =    " << N << "\n";
  cout << "  The number of combinations of N out of M is " << cnk << "\n";
  cout << "\n";
  cout << "   Rank	  Combination\n";
  cout << "\n";
 
  for ( rank = 1; rank <= 3; rank++ )
  {
    comb_unrank ( m, N, rank, a );
    cout                    << "  "
         << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  for ( rank = 6; rank <= 8; rank++ )
  {
    comb_unrank ( m, N, rank, a );
    cout << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  for ( rank = 250; rank <= 252; rank++ )
  {
    comb_unrank ( m, N, rank, a );
    cout                    << "  "
         << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test014 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST014 tests R8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cnk;
  int k;
  int n;

  cout << "\n";
  cout << "TEST014\n";
  cout << "  R8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "     N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );

      cout << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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
  cout << "TEST015\n";
  cout << "  I4_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "     N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      cout << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests COMP_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int h;
  int i;
  int index;
  bool more;
  int n = 6;
  int t;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  COMP_NEXT produces compositions.\n";
  cout << "\n";
  cout << "  Seeking all compositions of N = " << n << "\n";
  cout << "  using " << K << " parts.\n";
  cout << "\n";

  more = false;
  index = 0;

  for ( ; ; )
  {
    comp_next ( n, K, a, more, h, t );

    index = index + 1;
    cout << "  ";
    cout << "  " << setw(4) << index << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }
  }
 
  return;
# undef K
}
//****************************************************************************80

void test017 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST017 tests COMP_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define K 5

  int a[K];
  int i;
  int j;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "TEST017\n";
  cout << "  COMP_RANDOM produces compositions at random.\n";
  cout << "\n";
  cout << "  Seeking random compositions of N = " << n << "\n";
  cout << "  using " << K << " parts.\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 5; j++ )
  {
    comp_random ( n, K, seed, a );

    cout << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef K
}
//****************************************************************************80

void test0174 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0174 tests COMPNZ_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  bool more;
  int n = 6;

  cout << "\n";
  cout << "TEST0174\n";
  cout << "  COMPNZ_NEXT produces compositions using nonzero parts.\n";
  cout << "\n";
  cout << "  Seeking all compositions of N = " << n << "\n";
  cout << "  using " << K << " nonzero parts.\n";
  cout << "\n";

  more = false;
 
  for ( ; ; )
  {
    compnz_next ( n, K, a, more );

    cout << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef K
}
//****************************************************************************80

void test0175 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0175 tests COMPNZ_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define K 5

  int a[K];
  int i;
  int j;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "TEST0175\n";
  cout << "  COMPNZ_RANDOM produces compositions at random\n";
  cout << "  with only nonzero parts.\n";
  cout << "\n";
  cout << "  Seeking random compositions of N = " << n << "\n";
  cout << "  using " << K << " nonzero parts.\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 5; j++ )
  {
    compnz_random ( n, K, seed, a );

    cout << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef K
}
//****************************************************************************80

void test018 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST018 tests CONGRUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 20

  int a;
  int a_test[TEST_NUM] = {
     1027,   1027,  1027,   1027, -1027,
    -1027, - 1027, -1027,      6,     0,
        0,      0,     1,      1,     1,
     1024,      0,     0,      5,     2 };
  int b;
  int b_test[TEST_NUM] = {
      712,    712,  -712,   -712,   712,
      712,   -712,  -712,      8,     0,
        1,      1,     0,      0,     1,
   -15625,      0,     3,      0,     4 };
  int c;
  int c_test[TEST_NUM] = {
        7,     -7,     7,     -7,     7,
       -7,      7,    -7,     50,     0,
        0,      1,     0,      1,     0,
    11529,      1,    11,     19,     7 };
  bool error;
  int result;
  int test_i;
  int x;

  cout << "\n";
  cout << "TEST018\n";
  cout << "  CONGRUENCE solves a congruence equation:\n";
  cout << "    A * X = C mod ( B )\n";
  cout << "\n";
  cout << "   I        A         B         C         X     Mod ( A*X-C,B)\n";
  cout << "\n";

  for ( test_i = 1; test_i < TEST_NUM; test_i++ )
  {
    a = a_test[test_i];
    b = b_test[test_i];
    c = c_test[test_i];

    x = congruence ( a, b, c, error );

    if ( error )
    {
      cout                       << "  "
           << setw(2)  << test_i << "  "
           << setw(10) << a      << "  "
           << setw(10) << b      << "  "
           << setw(10) << c      << "  "
           << "(An error occurred)\n";
    }
    else
    {
      if ( b != 0 )
      {
        result = i4_modp ( a * x - c, b );
      }
      else
      {
        result = 0;
      }
      cout                        << "  "
           << setw(2)  << test_i  << "  "
           << setw(10) << a       << "  "
           << setw(10) << b       << "  "
           << setw(10) << c       << "  "
           << setw(10) << x       << "  "
           << setw(10) << result  << "\n";
    }

  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test019 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST019 tests COUNT_POSE_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int blocks[6];
  int goal;
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  COUNT_POSE_RANDOM poses a random problem for\n";
  cout << "  the game The Count is Good.\n";

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    count_pose_random ( seed, blocks, goal );

    cout << "\n";
    cout << "  Problem #" << i << "\n";
    cout << "\n";
    cout << "    The goal = " << goal << "\n";
    cout << "\n";
    cout << "    The available numbers are\n";
    cout << "\n";
    for ( j = 0; j < 6; j++ )
    {
      cout << setw(4) << blocks[j] << "  ";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test020 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST020 tests R8_TO_CFRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N+1];
  double error;
  int i;
  int p[N+2];
  int q[N+2];
  double r;
  double temp;

  cout << "\n";
  cout << "TEST020\n";
  cout << "  R8_TO_CFRAC converts a double precision number to\n";
  cout << "  a sequence of continued fraction convergents.\n";

  r = 2.0 * r8_pi ( );

  cout << "\n";
  cout << "  Use the real number R = " << r << "\n";

  r8_to_cfrac ( r, N, a, p, q );

  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    temp = ( double ) p[i+1] / ( double ) q[i+1];

    error = r - temp;

    cout                        << "  "
         << setw(12) << a[i]    << "  "
         << setw(12) << p[i+1]  << "  "
         << setw(12) << q[i+1]  << "  "
         << setw(14) << temp    << "  "
         << setw(14) << error   << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test021 ( )

//****************************************************************************80
//
//  Purpose:
//
//     TEST021 tests DEBRUIJN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NUM_TEST 3

  int i;
  int ihi;
  int m;
  int mtest[NUM_TEST] = { 2, 3, 2 };
  int n;
  int ntest[NUM_TEST] = { 3, 3, 4 };
  int string[28];
  int test;

  cout << "\n";
  cout << "TEST021\n";
  cout << "  DEBRUIJN computes a de Bruijn string.\n";

  for ( test = 0; test < NUM_TEST; test++ )
  {
    m = mtest[test];
    n = ntest[test];

    cout << "\n";
    cout << "  The alphabet size is M = " << m << "\n";
    cout << "  The string length is N = " << n << "\n";

    debruijn ( m, n, string );

    ihi = ( int ) pow ( ( double ) m, n );

    cout << "\n";
    cout << "  ";
    for ( i = 0; i < ihi; i++ )
    {
      cout << setw(1) << string[i];
    }
    cout << "\n";

  }

  return;
# undef NUM_TEST
}
//****************************************************************************80

void test022 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022 tests DEC_ADD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int abot;
  int atop;
  int bbot;
  int btop;
  int cbot;
  int ctop;
  int dec_digit;

  cout << "\n";
  cout << "TEST022\n";
  cout << "  DEC_ADD adds two decimals.\n";
  cout << "\n";

  dec_digit = 3;
  atop = 128;
  abot = -1;
  btop = 438;
  bbot = -2;

  dec_add ( atop, abot, btop, bbot, dec_digit, ctop, cbot );

  cout << "\n";
  cout << "  Number of decimal places is " << dec_digit << "\n";
  cout << "\n";

  cout << "  A = "         << atop << "*10^(" << abot << ")\n";
  cout << "  B = "         << btop << "*10^(" << bbot << ")\n";
  cout << "  C = A + B = " << ctop << "*10^(" << cbot << ")\n";
 
  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests DEC_DIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  int abot = -1;
  int atop = 523;
  int bbot = 2;
  int btop = 134;
  int cbot;
  int ctop;
  int dec_digit;
  bool error;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  DEC_DIV divides two decimals.\n";

  dec_digit = 3;

  dec_div ( atop, abot, btop, bbot, dec_digit, ctop, cbot, error );

  cout << "\n";
  cout << "  Number of decimal places is " << dec_digit << "\n";
  cout << "\n";

  cout << "  A = "         << atop << "*10^(" << abot << ")\n";
  cout << "  B = "         << btop << "*10^(" << bbot << ")\n";
  cout << "  C = A / B = " << ctop << "*10^(" << cbot << ")\n";

  return;
}
//****************************************************************************80

void test024 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024 tests DEC_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  int abot = -4;
  int atop = 14;
  int bbot = 2;
  int btop = 16;
  int cbot;
  int ctop;
  int dec_digit;

  cout << "\n";
  cout << "TEST024\n";
  cout << "  DEC_MUL multiplies two decimals.\n";

  dec_digit = 2;

  dec_mul ( atop, abot, btop, bbot, dec_digit, ctop, cbot );

  cout << "\n";
  cout << "  Number of decimal places is " << dec_digit << "\n";
  cout << "\n";

  cout << "  A = "         << atop << "*10^(" << abot << ")\n";
  cout << "  B = "         << btop << "*10^(" << bbot << ")\n";
  cout << "  C = A * B = " << ctop << "*10^(" << cbot << ")\n";

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests DEC_ROUND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 7

  int a;
  int a_test[N_TEST] = { 523, 523, 523, 523, 6340, 6340, 6340 };
  int b;
  int b_test[N_TEST] = { -1,  -1, -1, -1, 2, 2, 2 };
  int dec_digit;
  int r8_test[N_TEST] = { 1, 2, 3, 4, 2, 3, 4 };
  int i;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  DEC_ROUND rounds a decimal to a number of digits.\n";
  cout << "\n";
  cout << "           -----Before-------  -----After--------\n";
  cout << "  Digits   Mantissa  Exponent  Mantissa  Exponent\n";
  cout << "\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    dec_digit = r8_test[i];

    a = a_test[i];
    b = b_test[i];

    dec_round ( a, b, dec_digit, a, b );

    cout << setw(6) << r8_test[i] << "  "
         << setw(6) << a_test[i] << "  "
         << setw(6) << b_test[i] << "  "
         <<                         "  "
         << setw(6) << a         << "  "
         << setw(6) << b         << "\n";

  }

  return;
# undef N_TEST
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests DEC_TO_RAT and RAT_TO_DEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int exponent;
  int i;
  int mantissa;
  double r1;
  double r2;
  double r3;
  int rat_bot;
  int rat_bot2;
  int rat_top;
  int rat_top2;
  int seed;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  RAT_TO_DEC fraction => decimal,\n";
  cout << "  DEC_TO_RAT decimal => fraction.\n";
  cout << "\n";
  cout << "  In this test, choose the top and bottom\n";
  cout << "  of a rational at random, and compute the\n";
  cout << "  equivalent real number.\n";
  cout << "\n";
  cout << "  Then convert to decimal, and the equivalent real.\n";
  cout << "\n";
  cout << "  Then convert back to rational and the equivalent real.\n";
  
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    rat_top = i4_uniform ( -1000, 1000, seed );
    rat_bot = i4_uniform (     1, 1000, seed );

    r1 = ( double ) rat_top / ( double ) rat_bot;

    rat_to_dec ( rat_top, rat_bot, mantissa, exponent );
    r2 = ( double ) mantissa * pow ( 10.0, exponent );

    dec_to_rat ( mantissa, exponent, rat_top2, rat_bot2 );
    r3 = ( double ) rat_top2 / ( double ) rat_bot2;

    cout << "\n";
    cout << "  " << r1 << " = " << rat_top  << "/"     << rat_bot  << "\n";
    cout << "  " << r2 << " = " << mantissa << "*10^(" << exponent << ")\n";
    cout << "  " << r3 << " = " << rat_top2 << "/"     << rat_bot2 << "\n";
  }
 
  return;
}
//****************************************************************************80

void test027 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027 tests DEC_TO_S.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int exponent;
  int i;
  int mantissa;
  char *s;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  DEC_TO_S prints out a decimal.\n";
  cout << "\n";
  cout << "  Mantissa  Exponent  String\n";
  cout << "\n";

  mantissa = 523;
  exponent = -1;
  s = dec_to_s ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
                  << s        << "\n";

  mantissa = 134;
  exponent = 2;
  s = dec_to_s ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
                  << s        << "\n";

  mantissa = -134;
  exponent = 2;
  s = dec_to_s ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
                  << s        << "\n";
  mantissa = 0;
  exponent = 10;
  s = dec_to_s ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
                  << s        << "\n";

  for ( exponent = -8;exponent < 4; exponent++ )
  {
    mantissa = 123456;
    s = dec_to_s ( mantissa, exponent );
    cout << setw(8) << mantissa << "  "
         << setw(8) << exponent << "  "
                    << s        << "\n";
  }

  return;
}
//****************************************************************************80

void test028 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028 tests DEC_WIDTH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int exponent;
  int i;
  int mantissa;

  cout << "\n";
  cout << "TEST028\n";
  cout << "  DEC_WIDTH determines the \"width\" of a decimal.\n";
  cout << "\n";
  cout << "  Mantissa  Exponent  Width\n";
  cout << "\n";

  mantissa = 523;
  exponent = -1;
  i = dec_width ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
       << setw(8) << i        << "\n";

  mantissa = 134;
  exponent = 2;
  i = dec_width ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
       << setw(8) << i        << "\n";

  mantissa = -134;
  exponent = 2;
  i = dec_width ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
       << setw(8) << i        << "\n";

  mantissa = 0;
  exponent = 10;
  i = dec_width ( mantissa, exponent );
  cout << setw(8) << mantissa << "  "
       << setw(8) << exponent << "  "
       << setw(8) << i        << "\n";

  for ( exponent = -8; exponent < 4; exponent++ )
  {
    mantissa = 123456;
    i = dec_width ( mantissa, exponent );
    cout << setw(8) << mantissa << "  "
         << setw(8) << exponent << "  "
         << setw(8) << i        << "\n";
  }

  return;
}
//****************************************************************************80

void test029 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029 tests DECMAT_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N3 3
# define N4 4

  int a;
  int b;
  int a3[N3*N3];
  int a4[N4*N4];
  int b3[N3*N3];
  int b4[N4*N4];
  int i;
  int dbot;
  int dec_digit;
  int dtop;
  int j;
  int k;
  double r;

  cout << "\n";
  cout << "TEST029\n";
  cout << "  DECMAT_DET: determinant of a decimal matrix.\n";
  cout << "\n";
 
  dec_digit = 5;

  k = 0;
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      k = k + 1;
      a3[i+j*N3] = k;
    }
  }

  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      b3[i+j*N3] = 0;
    }
  }
 
  decmat_print ( N3, N3, a3, b3, "  The 123/456/789 matrix:" );

  decmat_det ( N3, a3, b3, dec_digit, dtop, dbot );
 
  cout << "\n";
  cout << "  Determinant of the 123/456/789 matrix = "
       << dtop << "* 10^(" 
       << dbot << ")\n"; 

  for ( i = 0; i < N4; i++ )
  {
    for ( j = 0; j < N4; j++ )
    {
      r = 1.0 / ( double ) ( i + j + 2 );
      
      r8_to_dec ( r, dec_digit, a, b );
      a4[i+j*N4] = a;
      b4[i+j*N4] = b;
    }
  }
 
  decmat_print ( N4, N4, a4, b4, "  The Hilbert matrix:" );
 
  decmat_det ( N4, a4, b4, dec_digit, dtop, dbot );
 
  cout << "\n";
  cout << "  Determinant of the Hilbert matrix = "
       << dtop << "* 10^(" 
       << dbot << ")\n"; 
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      if ( i == j )
      {
        a3[i+j*N3] = 2;
      }
      else if ( i == j+1 || i == j-1 )
      {
        a3[i+j*N3] = -1;
      }
      else
      {
        a3[i+j*N3] = 0;
      }
    }
  }
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      b3[i+j*N3] = 0;
    }
  }

  decmat_print ( N3, N3, a3, b3, "  The -1,2,-1 matrix:" );
 
  decmat_det ( N3, a3, b3, dec_digit, dtop, dbot );
 
  cout << "\n";
  cout << "  Determinant of the -1,2,-1 matrix = "
       << dtop << "* 10^("
       << dbot << ")\n"; 
 
  return;
# undef N3
# undef N4
}
//****************************************************************************80

void test021a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST021a tests DERANGE_ENUM, DERANGE_ENUM2 and DERANGE_ENUM3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int d[N+1];
  int i;

  cout << "\n";
  cout << "TEST021a\n";
  cout << "  DERANGE_ENUM counts derangements;\n";
  cout << "  DERANGE_ENUM2 counts derangements.\n";
  cout << "  DERANGE_ENUM3 counts derangements.\n";
  cout << "\n";
  cout << "       N    # of derangements\n";
  cout << "\n";

  derange_enum2 ( N, d );

  for ( i = 0; i<= N; i++ )
  {
    cout                                   << "  "
         << setw(8) << i                   << "  "
         << setw(8) << derange_enum ( i )  << "  "
         << setw(8) << d[i]                << "  "
         << setw(8) << derange_enum3 ( i ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test022a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022a tests DERANGE_BACK_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  bool more;
  int number;

  cout << "\n";
  cout << "TEST022a\n";
  cout << "  DERANGE_BACK_NEXT generates derangements\n";
  cout << "  using backtracking.\n";
  cout << "\n";

  more = false;
  number = 0;

  for ( ; ; )
  {
    derange_back_next ( N, a, more );

    if ( !more )
    {
      break;
    }

    number = number + 1;

    cout << setw(4) << number << "    ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test023a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023a tests DERANGE_WEED_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  int a[N];
  bool more;
  int number;

  cout << "\n";
  cout << "TEST023a\n";
  cout << "  DERANGE_WEED_NEXT generates derangements\n";
  cout << "  by generating ALL permutations, and weeding out\n";
  cout << "  the ones that are not derangements.\n";
  cout << "\n";

  more = false;
  number = 0;
 
  for ( ; ; )
  {
    derange_weed_next ( N, a, more );

    number = number + 1;

    cout << setw(4) << number << "    ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }
 
  }

  return;
# undef N
}
//****************************************************************************80

void test024a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024a calls DIGRAPH_ARC_EULER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NEDGE 7
# define NNODE 5

  int i;
  int in;
  int inode[NEDGE] = { 2, 1, 2, 1, 3, 5, 4 };
  int j;
  int jnode[NEDGE] = { 5, 4, 3, 2, 1, 1, 2 };
  int jp1;
  bool success;
  int trail[NEDGE];

  cout << "\n";
  cout << "TEST024a\n";
  cout << "  DIGRAPH_ARC_EULER finds an Euler circuit of a digraph.\n";

  digraph_arc_print ( NEDGE, inode, jnode, "  The arc list of the digraph:" );

  digraph_arc_euler ( NNODE, NEDGE, inode, jnode, success, trail );

  if ( success )
  {
    i4vec1_print ( NEDGE, trail, "  The edge list of the Euler circuit:" );

    cout << "\n";
    cout << "  The node list of the Euler circuit:\n";
    cout << "\n";
    cout << "	 I  Edge  Node\n";
    cout << "\n";

    for ( i = 0; i < NEDGE; i++ )
    {
      j = trail[i];

      if ( i+1 == NEDGE )
      {
        jp1 = trail[0];
      }
      else
      {
        jp1 = trail[i+1];
      }

      if ( jnode[j-1] == inode[jp1-1] )
      {
        in = jnode[j-1];
      }
      else
      {
        cout << "\n";
        cout << "The circuit has failed!\n";
        cout << "  JNODE[" << j-1 << "] = " << jnode[j-1] << "\n";
        cout << "  INODE[" << jp1-1 << "] = " << inode[jp1-1] << "\n";
        break;
      }

      cout << setw(6) << i << "  "
           << setw(6) << j << "  "
           << setw(6) << in << "\n";
    }
  }
  else
  {
    cout << "\n";
    cout << "  The digraph is not eulerian.\n";
    cout << "\n";
  }

  return;
# undef NEDGE
# undef NNODE
}
//****************************************************************************80

void test025a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025a tests DIOPHANTINE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 20

  int a;
  int a_test[TEST_NUM] = {
      1027,    1027,    1027,    1027,   -1027,
     -1027,   -1027,   -1027,       6,       0,
         0,       0,       1,       1,       1,
      1024,       0,       0,       5,       2 };
  int b;
  int b_test[TEST_NUM] = {
       712,     712,    -712,    -712,     712,
       712,    -712,    -712,       8,       0,
         1,       1,       0,       0,       1,
    -15625,       0,       3,       0,       4 };
  int c;
  int c_test[TEST_NUM] = {
         7,      -7,       7,      -7,       7,
        -7,       7,      -7,      50,       0,
         0,       1,       0,       1,       0,
     11529,       1,      11,      19,       7 };

  bool error;
  int r;
  int test_i;
  int x;
  int y;

  cout << "\n";
  cout << "TEST025a\n";
  cout << "  DIOPHANTINE solves a Diophantine equation:\n";
  cout << "    A * X + B * Y = C\n";
  cout << "\n";
  cout << "        A         B         C         X     Y     Residual\n";
  cout << "\n";

  for ( test_i = 0; test_i < TEST_NUM; test_i++ )
  {
    a = a_test[test_i];
    b = b_test[test_i];
    c = c_test[test_i];

    diophantine ( a, b, c, error, x, y );

    if ( error )
    {
      cout << setw(10) << a << "  "
           << setw(10) << b << "  "
           << setw(10) << c << "  "
           << "(Error occurred!)" << "\n";
    }
    else
    {
      r = a * x + b * y - c;
      cout << setw(10) << a << "  "
           << setw(10) << b << "  "
           << setw(10) << c << "  "
           << setw(10) << x << "  "
           << setw(10) << y << "  "
           << setw(10) << r << "\n";
    }

  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test026a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026a tests DIOPHANTINE_SOLUTION_MINIMIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a = 4096;
  int b = -15625;
  int c = 46116;
  int r;
  int x;
  int y;

  cout << "\n";
  cout << "TEST026a\n";
  cout << "  DIOPHANTINE_SOLUTION_MINIMIZE computes a minimal\n";
  cout << "  Euclidean norm solution of a Diophantine equation:\n";
  cout << "    A * X + B * Y = C\n";

  x = 665499996;
  y = 174456828;

  r = a * x + b * y - c;

  cout << "\n";
  cout << "  Coefficients:\n";
  cout << "    A = " << setw(12) << a << "\n";
  cout << "    B = " << setw(12) << b << "\n";
  cout << "    C = " << setw(12) << c << "\n";
  cout << "  Solution:\n";
  cout << "    X = " << setw(12) << x << "\n";
  cout << "    Y = " << setw(12) << y << "\n";
  cout << "  Residual R = A * X + B * Y - C:\n";
  cout << "    R = " << setw(12) << r << "\n";

  diophantine_solution_minimize ( a, b, x, y );

  r = a * x + b * y - c;

  cout << "\n";
  cout << "  DIOPHANTINE_SOLUTION_MINIMIZE returns\n";
  cout << "  the minimized solution:\n";
  cout << "    X = " << setw(12) << x << "\n";
  cout << "    Y = " << setw(12) << y << "\n";
  cout << "  Residual R = A * X + B * Y - C:\n";
  cout << "    R = " << setw(12) << r << "\n";

  x = 15621;
  y = 4092;

  r = a * x + b * y - c;

  cout << "\n";
  cout << "  Here is the minimal positive solution:\n";
  cout << "    X = " << setw(12) << x << "\n";
  cout << "    Y = " << setw(12) << y << "\n";
  cout << "  Residual R = A * X + B * Y - C:\n";
  cout << "    R = " << setw(12) << r << "\n";

  return;
}
//****************************************************************************80

void test026b ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026B tests DVEC_ADD and DVEC_SUB;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int dvec1[N];
  int dvec2[N];
  int dvec3[N];
  int dvec4[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST026B\n";
  cout << "  DVEC_ADD adds decimal vectors representing integers;\n";
  cout << "  DVEC_SUB subtracts decimal vectors representing integers;\n";
  cout << "\n";
  cout << "        I        J        K = I + J    L = I - J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( -100, 100, seed );
    j = i4_uniform ( -100, 100, seed );

    cout << "\n";
    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";

    k = i + j;
    l = i - j;

    cout << "  Directly:         "
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";

    i4_to_dvec ( i, N, dvec1 );
    i4_to_dvec ( j, N, dvec2 );

    dvec_add ( N, dvec1, dvec2, dvec3 );
    k = dvec_to_i4 ( N, dvec3 );

    dvec_sub ( N, dvec1, dvec2, dvec4 );
    l = dvec_to_i4 ( N, dvec4 );

    cout << "  DVEC_ADD, DVEC_SUB"
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test026c ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026C tests DVEC_COMPLEMENTX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int dvec1[N];
  int dvec2[N];
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST026C\n";
  cout << "  DVEC_COMPLEMENTX returns the ten's complement\n";
  cout << "  of a (signed) decimal vector;\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( -100, 100, seed );

    i4_to_dvec ( i, N, dvec1 );

    dvec_complementx ( N, dvec1, dvec2 );

    j = dvec_to_i4 ( N, dvec2 );

    cout << "\n";
    cout << "  I = " << "  " << i << "\n";
    cout << "  J = " << "  " << j << "\n";
    dvec_print ( N, dvec1, " " );
    dvec_print ( N, dvec2, " " );

  }

  return;
# undef N
}
//****************************************************************************80

void test026d ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026D tests DVEC_MUL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int dvec1[N];
  int dvec2[N];
  int dvec3[N];
  int i;
  int j;
  int k;
  int n2;
  int seed = 123456789;
  int test;
  int test_num = 10;
  int test2;
  int test2_num = 2;

  cout << "\n";
  cout << "TEST026D\n";
  cout << "  DVEC_MUL multiplies decimal vectors\n";
  cout << "  representing integers;\n";

  for ( test2 = 1; test2 <= test2_num; test2++ )
  {
    if ( test2 == 1 )
    {
      n2 = N;
    }
    else if ( test2 == 2 )
    {
      n2 = 6;

      cout << "\n";
      cout << "  NOW REPEAT THE TEST...\n";
      cout << "\n";
      cout << "  but use too few digits to represent big products.\n";
      cout << "  This corresponds to an \"overflow\".\n";
      cout << "  The result here should get the final decimal\n";
      cout << "  digits correctly, though.\n";
    }

    cout << "\n";
    cout << "        I        J        K = I * J\n";
    cout << "\n";

    for ( test = 1; test <= test_num; test++ )
    { 
      i = i4_uniform ( -1000, 1000, seed );
      j = i4_uniform ( -1000, 1000, seed );

      cout << "\n";
      cout << "  " << setw(8) << i
           << "  " << setw(8) << j << "\n";

      k = i * j;

      cout << "  Directly:         "
         << "  " << setw(8) << k << "\n";
      i4_to_dvec ( i, n2, dvec1 );
      i4_to_dvec ( j, n2, dvec2 );

      dvec_mul ( n2, dvec1, dvec2, dvec3 );
      k = dvec_to_i4 ( n2, dvec3 );

      cout << "  DVEC_MUL          "
           << "  " << setw(8) << k << "\n";
    }
  }
  return;
# undef N
}
//****************************************************************************80

void test027a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027a tests EQUIV_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int a[N];
  int i;
  int j;
  int jarray[N];
  bool more;
  int npart;
  int rank;

  cout << "\n";
  cout << "TEST027a\n";
  cout << "  EQUIV_NEXT generates all partitions of a set.\n";
  cout << "\n";
  cout << "  Rank//element:\n";
  cout << "\n";
  cout << "      ";
  for ( i = 1; i <= N; i++ )
  {
    cout << setw(2) << i << "  ";
  }
  cout << "\n";
  cout << "\n";
 
  rank = 0;
  more = false;
 
  for ( ; ; )
  {
    equiv_next ( N, npart, jarray, a, more );
 
    rank = rank + 1;

    cout                    << "  "
         << setw(2) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(2) << a[i] << "  ";
    }
    cout << "\n";
 
    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test028a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028a tests EQUIV_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int a[N];
  bool done;
  int i;
  int rank;

  cout << "\n";
  cout << "TEST028a\n";
  cout << "  EQUIV_NEXT2 generates all partitions of a set.\n";
  cout << "\n";
  cout << "  Rank//element:\n";
  cout << "\n";
  cout << "      ";
  for ( i = 1; i <= N; i++ )
  {
    cout << setw(2) << i << "  ";
  }
  cout << "\n";
  cout << "\n";
 
  rank = 0;
  done = true;
 
  for ( ; ; )
  {
    equiv_next2 ( done, a, N );

    if ( done )
    {
      break;
    }

    rank = rank + 1;

    cout                    << "  "
         << setw(2) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(2) << a[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test029a ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029a tests EQUIV_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int a[N];
  double b[N];
  int i;
  int npart;
  int seed;

  cout << "\n";
  cout << "TEST029a\n";
  cout << "  EQUIV_RANDOM selects a random set partition.\n";
 
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    equiv_random ( N, seed, npart, a, b );

    equiv_print ( N, a, "  The partition:" );
  }
  cout << "\n";
  cout << "  Now repeat, but print using EQUIV_PRINT2\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    equiv_random ( N, seed, npart, a, b );

    equiv_print2 ( N, a, "  The partition:" );
  }
  return;
# undef N
}
//****************************************************************************80

void test0295 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0295 tests EULER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 9

  int i;
  int ieuler[N_MAX+1];
  int n;

  cout << "\n";
  cout << "TEST0295\n";
  cout << "  EULER gets rows of the Euler triangle.\n";
  cout << "\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    euler ( n, ieuler );

    for ( i = 0; i <= n; i++ )
    {
      cout << setw(7) << ieuler[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N_MAX
}
//****************************************************************************80

void test030 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST030 tests R8_FALL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double x;

  cout << "\n";
  cout << "TEST030\n";
  cout << "  R8_FALL computes the falling factorial function.\n";
  cout << "  [X]_N = X * (X-1) * (X-2) * ... * ( X-N+1).\n";
  cout << "\n";
  cout << "    X          N  R8_FALL(X,N)\n";
  cout << "\n";

  x = 4.0;

  for ( n = -2; n <= 5; n++ )
  {
    cout << setw(10) << x << "  "
         << setw(8)  << n << "  "
         << setw(10) << r8_fall ( x, n ) << "\n";

  }

  return;
}
//****************************************************************************80

void test0304 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST0304 tests FROBENIUS_NUMBER_ORDER2 and FROBENIUS_NUMBER_ORDER2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  int c2;
  int f1;
  int f2;
  int n_data;

  cout << "\n";
  cout << "TEST0304\n";
  cout << "  FROBENIUS_NUMBER_ORDER2 computes Frobenius numbers of order 2.\n";
  cout << "  FROBENIUS_NUMBER_ORDER2_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "        C1        C1   exact F  comput F\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order2_values ( n_data, c1, c2, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = frobenius_number_order2 ( c1, c2 );

    cout << "  " << setw(8) << c1
         << "  " << setw(8) << c2
         << "  " << setw(8) << f1
         << "  " << setw(8) << f2 << "\n";
  }
  return;
}
//****************************************************************************80

void test0305 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0305 tests R8_GAMMA_LOG and GAMMA_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST0305:\n";
  cout << "  R8_GAMMA_LOG evaluates the logarithm of the Gamma function.\n";
  cout << "  GAMMA_LOG_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X       Exact F       GAMMA_LOG(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_gamma_log ( x );

    cout << setw(8) << x << "  "
         << setw(10) << fx << "  "
         << setw(10) << fx2 << "\n";

  }

  return;
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests GRAY_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int a[N];
  int change;
  int i;
  int k;

  cout << "\n";
  cout << "TEST031\n";
  cout << "  GRAY_NEXT returns the index of the single item\n";
  cout << "  to be changed in order to get the next Gray code.\n";

  cout << "\n";
  cout << "   K  Switch  Gray Code\n";
  cout << "\n";

  change = -N;
  k = 0;

  for ( ; ; )
  {
    gray_next ( N, change );

    if ( change == -N )
    {
      break;
    }
    else if ( change == 0 )
    {
      for ( i = 0; i < N; i++ )
      {
        a[i] = 0;
      }
    }
    else
    {
      a[abs(change)-1] = 1 - a[abs(change)-1];
    }

    cout                      << "  "
         << setw(2) << k      << "  "
         << setw(6) << change << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(1) << a[i];
    }
    cout << "\n";
    k = k + 1;
    
  }

  return;
# undef N
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests GRAY_RANK and GRAY_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int gray;
  int rank;
  int rank2;

  cout << "\n";
  cout << "TEST032\n";
  cout << "  GRAY_RANK ranks a Gray code;\n";
  cout << "  GRAY_UNRANK unranks a Gray code.\n";
  cout << "\n";
  cout << "    R  =                        RANK\n";
  cout << "    G  =            GRAY_UNRANK(RANK)\n";
  cout << "    R2 =  GRAY_RANK(GRAY_UNRANK(RANK))\n";
  cout << "\n";
  cout << "    R    G    R2\n";
  cout << "\n";
 
  for ( rank = 0; rank <= 24; rank++ )
  {
    gray = gray_unrank ( rank );

    rank2 = gray_rank ( gray );

    cout << setw(9) << rank << "  "
         << setw(9) << gray << "  "
         << setw(9) << rank2 << "\n";
  }
 
  return;
}
//****************************************************************************80

void test0321 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0321 tests GRAY_RANK2 and GRAY_UNRANK2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int gray;
  int rank;
  int rank2;

  cout << "\n";
  cout << "TEST0321\n";
  cout << "  GRAY_RANK2 ranks a Gray code;\n";
  cout << "  GRAY_UNRANK2 unranks a Gray code.\n";
  cout << "\n";
  cout << "    R  =                          RANK\n";
  cout << "    G  =             GRAY_UNRANK2(RANK)\n";
  cout << "    R2 =  GRAY_RANK2(GRAY_UNRANK2(RANK))\n";
  cout << "\n";
  cout << "    R    G    R2\n";
  cout << "\n";
 
  for ( rank = 0; rank <= 24; rank++ )
  {
    gray = gray_unrank2 ( rank );

    rank2 = gray_rank2 ( gray );

    cout << setw(9) << rank << "  "
         << setw(9) << gray << "  "
         << setw(9) << rank2 << "\n";
  }
 
  return;
}
//****************************************************************************80

void test0322 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0322 tests I4_BCLR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "TEST0322\n";
  cout << "  I4_BCLR sets a given bit to 0.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Working on I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit       I4_BCLR\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bclr ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test03225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03225 tests I4_BSET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "TEST03225\n";
  cout << "  I4_BSET sets a given bit to 0.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Working on I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit       I4_BSET\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bset ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0323 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0323 tests I4_BTEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "TEST0323\n";
  cout << "  I4_BTEST reports whether a given bit is 0 or 1.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Analyze the integer I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit  I4_BTEST\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_btest ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0324 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0324 tests I4_FACTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2005
//
{
# define FACTOR_MAX 10

  int factor[FACTOR_MAX];
  int factor_num;
  int i;
  int n;
  int nleft;
  int power[FACTOR_MAX];

  cout << "\n";
  cout << "TEST0324\n";
  cout << "  I4_FACTOR factors an integer,\n";

  n = 2 * 2 * 17 * 37;

  cout << "\n";
  cout << "  The integer is " << n << "\n";

  i4_factor ( n, FACTOR_MAX, factor_num, factor, power, nleft );

  cout << "\n";
  cout << "  Prime representation:\n";
  cout << "\n";
  cout << "  I  FACTOR(I)  POWER(I)\n";
  cout << "\n";

  if ( abs ( nleft ) != 1 )
  {
    cout << "  " << setw(6) << 0
         << "  " << setw(6) << nleft 
         << "  " << "(Unfactored portion)\n";
  }

  for ( i = 0; i < factor_num; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << factor[i]
         << "  " << setw(6) << power[i] << "\n";
  }
 
  return;
# undef FACTOR_MAX
}
//****************************************************************************80

void test0325 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0325 tests I4_GCD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int seed;

  cout << "\n";
  cout << "TEST0325\n";
  cout << "  I4_GCD computes the greatest common divisor\n";
  cout << "  of two integers.\n";

  cout << "\n";
  cout << "     I     J    I4_GCD(I,J)\n";
  cout << "\n";

  seed = 123456789;

  for ( k = 1; k <= 15; k++ )
  {
    i = i4_uniform ( -5, 15, seed );
    j = i4_uniform (  1, 15, seed );

    cout << setw(4) << i << "  "
         << setw(4) << j << "  "
         << setw(4) << i4_gcd ( i, j ) << "\n";
  }

  return;
}
//****************************************************************************80

void test0327 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0327 tests I4_LOG_10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 21

  int i;
  int x[N] = {
    0, 1, 2, 3, 9, 10, 11, 99, 100, 101, 999, 1000, 1001,
   -1, -2, -3, -9, -10, -11, -99, -101 };

  cout << "\n";
  cout << "TEST0327\n";
  cout << "  I4_LOG_10: whole part of log base 10,\n";
  cout << "\n";
  cout << "     X I4_LOG_10\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(6) << x[i] << "  "
         << setw(6) << i4_log_10 ( x[i] ) << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test058 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058 tests I4_PARTITION_CONJ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 14
# define NPART1 4

  int a1[NPART1] = { 2, 5, 1, 4 };
  int a2[N];
  int i;
  int mult1[NPART1] = { 1, 1, 3, 1 };
  int mult2[N];
  int npart2;

  cout << "\n";
  cout << "TEST058\n";
  cout << "  I4_PARTITION_CONJ conjugates an integer partition.\n";
  cout << "\n";
  cout << "  Original partition:\n";
  cout << "\n";

  i4_partition_print ( N, NPART1, a1, mult1 );

  i4_partition_conj ( N, a1, mult1, NPART1, a2, mult2, npart2 );

  cout << "\n";
  cout << "  Conjugate partition:\n";
  cout << "\n";

  i4_partition_print ( N, npart2, a2, mult2 );

  return;
# undef N
# undef NPART1
}

//****************************************************************************80

void test059 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST059 tests I4_PARTITION_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N];
  bool done;
  int i;
  int mult[N];
  int npart;
  int rank;

  cout << "\n";
  cout << "TEST059\n";
  cout << "  I4_PARTITION_NEXT generates partitions of an integer.\n";
  cout << "  Here N = " << N << "\n";
  cout << "\n";

  rank = 0;
  done = true;
 
  for ( ; ; )
  {
    i4_partition_next ( done, a, mult, N, npart );
 
    if ( done )
    {
      break;
    }

    rank = rank + 1;

    i4_partition_print ( N, npart, a, mult );

  }
 
  return;
# undef N
}
//****************************************************************************80

void test060 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST060 tests I4_PARTITION_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N];
  int i;
  bool more;
  int mult[N];
  int npart;

  cout << "\n";
  cout << "TEST060\n";
  cout << "  I4_PARTITION_NEXT2 produces partitions of an integer.\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    i4_partition_next2 ( N, a, mult, npart, more );

    i4_partition_print ( N, npart, a, mult );

    if ( !more )
    {
      break;
    }

  }
  
  return;
# undef N
}
//****************************************************************************80

void test061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST061 tests I4_PARTITION_COUNT and I4_PARTITION_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  int n;
  int n_data;
  int p;
  int p2[N_MAX+1];

  cout << "\n";
  cout << "TEST061\n";
  cout << "  I4_PARTITION_COUNT counts partitions of an integer.\n";
  cout << "  I4_PARTITION_COUNT_VALUES returns some exact values.\n";

  n_data = 0;

  cout << "\n";
  cout << "   N     Exact     Count\n";
  cout << "\n";

  for ( ; ; )
  {
    i4_partition_count_values ( n_data, n, p );

    if ( n_data == 0 )
    {
      break;
    }

    cout << setw(4)  << n << "  "
         << setw(10) << p << "  ";

    if ( n <= N_MAX )
    {
      i4_partition_count ( n, p2 );
      cout << setw(10) << p2[n];
    }

    cout << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0615 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0615 tests I4_PARTITION_COUNT2 and I4_PARTITION_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  int n;
  int n_data;
  int p;
  int *p2;

  cout << "\n";
  cout << "TEST0615\n";
  cout << "  I4_PARTITION_COUNT2 counts partitions of an integer.\n";
  cout << "  I4_PARTITION_COUNT_VALUES returns some exact values.\n";

  n_data = 0;

  cout << "\n";
  cout << "   N     Exact     Count\n";
  cout << "\n";

  for ( ; ; )
  {
    i4_partition_count_values ( n_data, n, p );

    if ( n_data == 0 )
    {
      break;
    }

    cout                  << "  "
         << setw(4)  << n << "  "
         << setw(10) << p << "  ";

    if ( n <= N_MAX )
    {
      p2 = i4_partition_count2 ( n );
      cout << "  "
           << "      "
           << setw(10) << p2[n];
      delete [] p2;
    }

    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests I4_PARTITION_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 8

  int a[N];
  int i;
  int j;
  int mult[N];
  int npart;
  int seed;
  int *table;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  I4_PARTITION_RANDOM generates a random partition.\n";
  cout << "  I4_PARTITION_COUNT2 sets up a necessary table.\n";
  cout << "\n";

  seed = 123456789;
//
//  Call once just to get the partition table.
//
  table = i4_partition_count2 ( N );

  cout << "\n";
  cout << "  The number of partitions of N.\n";
  cout << "\n";
  cout << "     N    Number of partitions\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    cout << setw(6) << j+1      << "  "
         << setw(6) << table[j] << "\n";
  }

  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    i4_partition_random ( N, table, seed, a, mult, npart );

    i4_partition_print ( N, npart, a, mult );

  }
 
  delete [] table;

  return;
# undef N
}
//****************************************************************************80

void test06225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06225 tests I4_PARTITIONS_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m[3];
  int msum;
  int s = 3;

  cout << "\n";
  cout << "TEST06225\n";
  cout << "  I4_PARTITIONS_NEXT produces the next\n";
  cout << "  nondecreasing partitions of an integer, and\n";
  cout << "  if necessary, increments the integer to keep on going.\n";

  i = 0;
  m[0] = 0;
  m[1] = 0;
  m[2] = 0;

  cout << "\n";
  cout << "   I Sum    Partition\n";
  cout << "\n";
  msum = i4vec_sum ( s, m );
  cout << "  " << setw(2) << i
       << "  " << setw(2) << msum  << "    ";
  for ( j = 0; j < s; j++ )
  {
    cout << setw(2) << m[j];
  }
  cout << "\n";

  for ( i = 1; i <= 15; i++ )
  {
    i4_partitions_next ( s, m );
    msum = i4vec_sum ( s, m );
    cout << "  " << setw(2) << i
         << "  " << setw(2) << msum  << "    ";
    for ( j = 0; j < s; j++ )
    {
      cout << setw(2) << m[j];
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "  You can start from any legal partition.\n";
  cout << "  Here, we restart at ( 2, 1, 0 ).\n";

  i = 0;
  m[0] = 2;
  m[1] = 1;
  m[2] = 0;

  cout << "\n";
  cout << "   I Sum    Partition\n";
  cout << "\n";
  msum = i4vec_sum ( s, m );
  cout << "  " << setw(2) << i
       << "  " << setw(2) << msum  << "    ";
  for ( j = 0; j < s; j++ )
  {
    cout << setw(2) << m[j];
  }
  cout << "\n";

  for ( i = 1; i <= 15; i++ )
  {
    i4_partitions_next ( s, m );
    msum = i4vec_sum ( s, m );
    cout << "  " << setw(2) << i
         << "  " << setw(2) << msum  << "    ";
    for ( j = 0; j < s; j++ )
    {
      cout << setw(2) << m[j];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests I4_SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int q;
  int r;

  cout << "\n";
  cout << "TEST033\n";
  cout << "  I4_SQRT computes the square root of an integer.\n";
  cout << "\n";
  cout << "       N  Sqrt(N) Remainder\n";
  cout << "\n";

  for ( n = -5; n <= 20; n++ )
  {
    i4_sqrt ( n, q, r );

    cout << setw(9) << n << "  "
         << setw(9) << q << "  "
         << setw(9) << r << "\n";
  }

  return;
}
//****************************************************************************80

void test034 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST034 tests I4_SQRT_CF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAX_TERM 100

  int b[MAX_TERM+1];
  int i;
  int n;
  int n_term;

  cout << "\n";
  cout << "TEST034\n";
  cout << "  I4_SQRT_CF computes the continued fraction form\n";
  cout << "  of the square root of an integer.\n";
  cout << "\n";
  cout << "   N  Period  Whole  Repeating Part\n";
  cout << "\n";

  for ( n = 1; n <= 20; n++ )
  {
    i4_sqrt_cf ( n, MAX_TERM, n_term, b );
    cout << setw(5) << n << "  "
         << setw(5) << n_term << "  ";
    for ( i = 0; i <= n_term; i++ )
    {
      cout << setw(5) << b[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef MAX_TERM
}
//****************************************************************************80

void test0625 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0625 tests I4_TO_BVEC and BVEC_TO_I4;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int bvec[N];
  int i;
  int i2;
  int j;

  cout << "\n";
  cout << "TEST0625\n";
  cout << "  I4_TO_BVEC converts an integer to a \n";
  cout << "    signed binary vector;\n";
  cout << "  BVEC_TO_I4 converts a signed binary vector\n";
  cout << "    to an integer;\n";
  cout << "\n";
  cout << "  I --> BVEC  -->  I\n";
  cout << "\n";

  for ( i = -3; i <= 10; i++ )
  {
    i4_to_bvec ( i, N, bvec );
    i2 = bvec_to_i4 ( N, bvec );

    cout << setw(3) << i << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(1) << bvec[j];
    }
    cout << "  ";
    cout << setw(3) << i2 << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 tests I4_TO_CHINESE and CHINESE_TO_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  int j;
  int j2;
  int m[N] = { 3, 4, 5, 7 };
  int r[N];

  cout << "\n";
  cout << "TEST035\n";
  cout << "  I4_TO_CHINESE computes the Chinese Remainder\n";
  cout << "    representation of an integer.\n";
  cout << "  CHINESE_TO_I4 computes an integer with the given\n";
  cout << "    Chinese Remainder representation.\n";

  i4vec1_print ( N, m, "  The moduli:" );

  j = 37;

  cout << "\n";
  cout << "  The number being analyzed is " << j << "\n";

  i4_to_chinese ( j, N, m, r );

  i4vec1_print ( N, r, "  The remainders:" );

  j2 = chinese_to_i4 ( N, m, r );

  cout << "\n";
  cout << "  The reconstructed number is " << j2 << "\n";

  i4_to_chinese ( j2, N, m, r );

  i4vec1_print ( N, r, "  The remainders of the reconstructed number are:" );

  return;
# undef N
}
//****************************************************************************80

void test0627 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0627 tests I4_TO_I4POLY and I4POLY_TO_I4;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DEGREE_MAX 5
# define TEST_NUM 9

  int a[DEGREE_MAX+1];
  int base;
  int base_test[TEST_NUM] = { 2, 2, 2, 3, 4, 5, 6, 23, 24 };
  int degree;
  int i;
  int intval;
  int intval2;
  int intval_test[TEST_NUM] = { 1, 6, 23, 23, 23, 23, 23, 23, 23 };
  int test;

  cout << "\n";
  cout << "TEST0627\n";
  cout << "  I4_TO_I4POLY converts an integer to a polynomial\n";
  cout << "    in a given base;\n";
  cout << "  I4POLY_TO_I4 evaluates an integer polynomial\n";
  cout << "    at a given point;\n";
  cout << "\n";
  cout << "       I    BASE  DEGREE  Coefficients\n";
  cout << "\n";
  for ( test = 0; test < TEST_NUM; test++ )
  {
    intval = intval_test[test];
    base = base_test[test];
    i4_to_i4poly ( intval, base, DEGREE_MAX, degree, a );
    cout                      << "  "
         << setw(6) << intval << "  "
         << setw(6) << base   << "  "
         << setw(6) << degree << "  ";
    for ( i = 0; i <= degree; i++ )
    {
      cout << setw(6) << a[i] << "  ";
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "  Now let I4_TO_I4POLY convert I to a polynomial,\n";
  cout << "  use I4POLY_TO_I4 to evaluate it, and compare.\n";
  cout << "\n";
  cout << "       I    I2\n";
  cout << "\n";
  for ( test = 0; test < TEST_NUM; test++ )
  {
    intval = intval_test[test];
    base = base_test[test];
    i4_to_i4poly ( intval, base, DEGREE_MAX, degree, a );
    intval2 = i4poly_to_i4 ( degree, a, base );
    cout                       << "  "
         << setw(6) << intval  << "  "
         << setw(6) << intval2 << "\n";
  }

  return;

# undef DEGREE_MAX
# undef TEST_NUM
}
//****************************************************************************80

void test0364 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0364 tests I4_TO_I4POLY and I4_TO_VAN_DER_CORPUT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DEGREE_MAX 10

  int base;
  int degree;
  double h1;
  double h2;
  int i;
  int ipoly[DEGREE_MAX+1];
  int j;
  int n_test = 10;
  int value;

  cout << "\n";
  cout << "TEST0364\n";
  cout << "  I4_TO_VAN_DER_CORPUT computes the elements\n";
  cout << "    of a van der Corput sequence.\n";
  cout << "\n";
  cout << "  I4_TO_I4POLY converts an integer to an integer\n";
  cout << "    polynomial in some base, and can be used to mimic\n";
  cout << "    the van der Corput calculation.\n";
  cout << "\n";

  for ( j = 1; j <= 3; j++ )
  {
    base = prime(j);
    cout << "\n";
    cout << "  BASE = " << base << "\n";
    cout << "\n";

    for ( i = 1; i <= n_test; i++ )
    {
      h1 = i4_to_van_der_corput ( i, base );

      i4_to_i4poly ( i, base, DEGREE_MAX, degree, ipoly );

      i4vec_reverse ( degree + 1, ipoly );

      value = i4poly_to_i4 ( degree, ipoly, base );

      h2 = ( double ) ( value ) / pow ( ( double ) base, ( double ) (degree+1) );

      cout                   << "  "
           << setw(6)  << i  << "  "
           << setw(14) << h1 << "  "
           << setw(14) << h2 << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 tests I4_TO_VAN_DER_CORPUT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int i;
  int j;
  int p;

  cout << "\n";
  cout << "TEST036\n";
  cout << "  I4_TO_VAN_DER_CORPUT computes the elements \n";
  cout << "  of a van der Corput sequence.\n";
  cout << "  The sequence depends on the prime number used\n";
  cout << "  as a base.\n";
  cout << "\n";
  cout << "Base: ";
  for ( j = 1; j <= 5; j++ )
  {
    p = prime ( j );
    cout << setw(10) << p << "  ";
  }
  cout << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << setw(4) << i << "  ";
    for ( j = 1; j <= 5; j++ )
    {
      p = prime ( j );
      h = i4_to_van_der_corput ( i, p );
      cout << setw(10) << h << "  ";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests I4MAT_01_ROWCOLSUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  int a[M*N];
  int c[N] = { 2, 2, 2, 2, 1 };
  bool error;
  int r[M] = { 3, 2, 2, 1, 1 };

  cout << "\n";
  cout << "TEST037\n";
  cout << "  I4MAT_01_ROWCOLSUM constructs a 01 matrix with\n";
  cout << "  given row and column sums.\n";
  
  i4vec1_print ( M, r, "  The rowsum vector:" );
  i4vec1_print ( N, c, "  The columnsum vector: " );

  i4mat_01_rowcolsum ( M, N, r, c, a, error );

  if ( error )
  {
    cout << "\n";
    cout << "  I4MAT_01_ROWCOLSUM returned error flag.\n";
  }
  else
  {
    i4mat_print ( M, N, a, "  The rowcolsum matrix:" );
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test038 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST038 tests I4MAT_01_ROWCOLSUM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  int a[M*N];
  int c[N] = { 2, 1, 2, 2, 2 };
  bool error;
  int r[M] = { 2, 1, 3, 1, 2 };

  cout << "\n";
  cout << "TEST038\n";
  cout << "  I4MAT_01_ROWCOLSUM2 constructs a 01 matrix with\n";
  cout << "  given row and column sums.\n";
  
  i4vec1_print ( M, r, "  The rowsum vector:" );
  i4vec1_print ( N, c, "  The columnsum vector: " );
//
//  This call commented out because the routine is not ready.
//
  cout << "\n";
  cout << "TEST038\n";
  cout << "  Cancel test!  The routine is not ready!\n";

  if ( 1 < 2 )
  {
    return;
  }

  i4mat_01_rowcolsum2 ( M, N, r, c, a, error );

  if ( error )
  {
    cout << "\n";
    cout << "  I4MAT_01_ROWCOLSUM2 returned an error flag.\n";
    cout << "  The matrix returned is not an exact solution.\n";
  }

  i4mat_print ( M, N, a, "  The rowcolsum matrix:" );

  cout << "\n";
  cout << "  Now repeat, with data for which there is not\n";
  cout << "  a solution.  The program will try its best anyway.\n";

  c[0] = 1;
  c[1] = 4;
  c[2] = 1;
  c[3] = 5;
  c[4] = 1;

  r[0] = 1;
  r[1] = 3;
  r[2] = 4;
  r[3] = 1;
  r[4] = 3;

  i4vec1_print ( M, r, "  The rowsum vector:" );
  i4vec1_print ( N, c, "  The columnsum vector: " );

//i4mat_01_rowcolsum2 ( M, N, r, c, a, error );

  if ( error )
  {
    cout << "\n";
    cout << "  I4MAT_01_ROWCOLSUM2 returned an error flag.\n";
    cout << "  The matrix returned is not an exact solution.\n";
  }

  i4mat_print ( M, N, a, "  The rowcolsum matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test039 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST039 tests I4MAT_01_ROWCOLSUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 13
# define N 17

  int a[M*N];
  int c[N] = {  4,  4, 11, 10, 10,  8,  9, 10,  8,  9,  3, 10,  4,  7,  9,  3,  3 };
  int c2[N];
  int c3[N];
  int *cperm;
  bool error;
  int i;
  int j;
  int r[M] = { 14, 13, 14, 10, 12,  2, 10,  1, 10, 11,  6,  2, 17 };
  int r2[M];
  int r3[M];
  int *rperm;

  cout << "\n";
  cout << "TEST039\n";
  cout << "  I4MAT_01_ROWCOLSUM constructs a 01 matrix with\n";
  cout << "  given row and column sums.\n";

  i4vec1_print ( M, r, "  The rowsum vector R:" );

  rperm = i4vec_sort_heap_index_d ( M, r );

  for ( i = 0; i < M; i++ )
  {
    r2[i] = r[rperm[i]];
  }

  i4vec1_print ( N, c, "  The columnsum vector C: " );

  cperm = i4vec_sort_heap_index_d ( N, c );

  for ( j = 0; j < N; j++ )
  {
    c2[j] = c[cperm[j]];
  }

  i4mat_01_rowcolsum ( M, N, r2, c2, a, error );

  if ( error )
  {
    cout << "\n";
    cout << "  I4MAT_01_ROWCOLSUM returned an error flag.\n";
    return;
  }
//
//  RPERM and CPERM are 0 based, but must be changed to 1 based.
//
  for ( i = 0; i < M; i++ )
  {
    rperm[i] = rperm[i] + 1;
  } 
  for ( i = 0; i < N; i++ )
  {
    cperm[i] = cperm[i] + 1;
  }   
//
//  Invert the R and C permutations.
//
  i4mat_perm2 ( M, N, a, rperm, cperm );

  i4mat_print ( M, N, a, "  The rowcolsum matrix:" );

  for ( i = 0; i < M; i++ )
  {
    r3[i] = 0;
    for ( j = 0; j < N; j++ )
    {
      r3[i] = r3[i] + a[i+j*M];
    }
  }

  i4vec1_print ( M, r3, "  Computed row sums" );

  for ( j = 0; j < N; j++ )
  {
    c3[j] = 0;
    for ( i = 0; i < M; i++ )
    {
      c3[j] = c3[j] + a[i+j*M];
    }
  }

  i4vec1_print ( N, c3, "  Computed column sums:" );

  delete [] cperm;
  delete [] rperm;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test040 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST040 tests I4MAT_U1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int a[N*N] = {
    1, 0, 0, 0, 0, 0, 
    1, 1, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 0, 1, 0, 
   75, 0, 0, 0, 1, 1 };
  int  b[N*N];

  cout << "\n";
  cout << "TEST040\n";
  cout << "  I4MAT_U1_INVERSE inverts a unit upper triangular matrix.\n";

  i4mat_print ( N, N, a, "  The input matrix:" );

  i4mat_u1_inverse ( N, a, b );

  i4mat_print ( N, N, b, "  The inverse matrix:" );

  return;
# undef N
}
//****************************************************************************80

void test041 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST041 tests I4MAT_PERM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  int a[N*N];
  int i;
  int j;
  int p[N] = { 2,3,9,6,7,8,5,4,1 };

  cout << "\n";
  cout << "TEST041\n";
  cout << "  I4MAT_PERM reorders an integer matrix in place.\n";
  cout << "  The rows and columns use the same permutation.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = (i+1) * 10 + (j+1);
    }
  }

  i4mat_print ( N, N, a, "  The input matrix:" );
 
  perm_print ( N, p, "  The row and column permutation:" );
 
  i4mat_perm ( N, a, p );
 
  i4mat_print ( N, N, a, "  The permuted matrix:" );
 
  return;
# undef N
}
//****************************************************************************80

void test042 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST042 tests I4MAT_PERM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 9
# define N 7

  int a[M*N];
  int i;
  int j;
  int p[M] = { 2,3,9,6,7,8,5,4,1 };
  int q[N] = { 3,4,5,6,7,1,2 };

  cout << "\n";
  cout << "TEST042\n";
  cout << "  I4MAT_PERM2 reorders an integer matrix in place.\n";
  cout << "  Rows and columns use different permutations.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = (i+1) * 10 + (j+1);
    }
  }
 
  i4mat_print ( M, N, a, "  The input matrix:" );
 
  perm_print ( M, p, "  The row permutation:" );

  perm_print ( N, q, "  The column permutation:" );
 
  i4mat_perm2 ( M, N, a, p, q );
 
  i4mat_print ( M, N, a, "  The permuted matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test043 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST043 tests INDEX_BOX_NEXT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  bool more;
  int n1 = 5;
  int n2 = 3;
  int n;

  cout << "\n";
  cout << "TEST043\n";
  cout << "  INDEX_BOX_NEXT_2D produces IJ indices that\n";
  cout << "  lie on the surface of a box in 2D.\n";
  cout << "\n";
  cout << "  The box has logical dimensions:\n";
  cout << setw(3) << n1 << "  "
       << setw(3) << n2 << "\n";
  cout << "\n";
  cout << "   #    I   J\n";
  cout << "\n";

  more = false;
  n = 0;

  for ( ; ; )
  {
    index_box_next_2d ( n1, n2, i, j, more );

    if ( !more )
    {
      break;
    }

    n = n + 1;
    cout << setw(3) << n << "  "
         << setw(3) << i << "  "
         << setw(3) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test044 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST044 tests INDEX_BOX_NEXT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  bool more;
  int n1 = 5;
  int n2 = 3;
  int n3 = 4;
  int n;

  cout << "\n";
  cout << "TEST044\n";
  cout << "  INDEX_BOX_NEXT_3D produces IJK indices that\n";
  cout << "  lie on the surface of a box.\n";
  cout << "\n";
  cout << "  The box has logical dimensions:\n";
  cout << setw(3) << n1 << "  "
       << setw(3) << n2 << "  "
       << setw(3) << n3 << "\n";
  cout << "\n";
  cout << "   #    I   J   K\n";
  cout << "\n";

  more = false;
  n = 0;

  for ( ; ; )
  {
    index_box_next_3d ( n1, n2, n3, i, j, k, more );

    if ( !more )
    {
      break;
    }

    n = n + 1;
    cout << setw(3) << n << "  "
         << setw(3) << i << "  "
         << setw(3) << j << "  "
         << setw(3) << k << "\n";

  }

  return;
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests INDEX_BOX2_NEXT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ic = 10;
  int j;
  int jc = 20;
  bool more;
  int n1 = 4;
  int n2 = 3;
  int n;

  cout << "\n";
  cout << "TEST045\n";
  cout << " INDEX_BOX2_NEXT_2D produces IJ indices that\n";
  cout << "  lie on the surface of a box2 in 2D.\n";
  cout << "\n";
  cout << "  The box has half-widths:\n";
  cout << setw(3) << n1 << "  "
       << setw(3) << n2 << "\n";
  cout << "\n";
  cout << "  and has center cell:\n";
  cout << setw(3) << ic << "  "
       << setw(3) << jc << "\n";
  cout << "\n";
  cout << "   #    I   J\n";
  cout << "\n";

  more = false;
  n = 0;

  for ( ; ; )
  {
    index_box2_next_2d ( n1, n2, ic, jc, i, j, more );

    if ( !more )
    {
      break;
    }

    n = n + 1;
    cout << setw(3) << n << "  "
         << setw(3) << i << "  "
         << setw(3) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test046 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST046 tests INDEX_BOX2_NEXT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ic = 10;
  int j;
  int jc = 20;
  int k;
  int kc = 30;
  bool more;
  int n1 = 5;
  int n2 = 3;
  int n3 = 4;
  int n;

  cout << "\n";
  cout << "TEST046\n";
  cout << "  INDEX_BOX2_NEXT_3D produces IJK indices that\n";
  cout << "  lie on the surface of a box.\n";
  cout << "\n";
  cout << "  The box has half widths:\n";
  cout << setw(3) << n1 << "  "
       << setw(3) << n2 << "  "
       << setw(3) << n3 << "\n";
  cout << "\n";
  cout << "  and central cell:\n";
  cout << setw(3) << ic << "  "
       << setw(3) << jc << "  "
       << setw(3) << kc << "\n";
  cout << "\n";
  cout << "  We will only print a PORTION of the data!\n";
  cout << "\n";
  cout << "   #    I   J   K\n";
  cout << "\n";

  more = false;
  n = 0;

  for ( ; ; )
  {
    index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more );

    if ( !more )
    {
      break;
    }

    n = n + 1;

    if ( n <= 10 || 370 <= n )
    {
      cout << setw(3) << n << "  "
           << setw(3) << i << "  "
           << setw(3) << j << "  "
           << setw(3) << k << "\n";
    }

  }

  return;
}
//****************************************************************************80

void test047 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST047 tests I4MAT_PERM2 and TRIANG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N*N] = {
    1,0,1,0,1,0,1,0,0,1, 
    0,1,0,0,1,0,0,0,0,0, 
    0,0,1,0,1,0,1,0,0,1, 
    0,1,1,1,1,1,1,1,0,1, 
    0,0,0,0,1,0,0,0,0,0, 
    0,1,0,0,1,1,1,0,0,0, 
    0,0,0,0,1,0,1,0,0,0, 
    0,1,0,0,1,1,1,1,0,1, 
    0,0,0,0,0,0,0,0,0,0, 
    0,0,0,0,1,0,1,0,0,1 };
  int i;
  int p[N];
  int j;

  cout << "\n";
  cout << "TEST047\n";
  cout << "  TRIANG relabels elements for a partial ordering,\n";
  cout << "  I4MAT_PERM2 reorders an integer matrix in place.\n";

  i4mat_print ( N, N, a, "  The input matrix:" );
 
  triang ( N, a, p );
 
  perm_print ( N, p, "  The new ordering:" );

  i4mat_perm2 ( N, N, a, p, p );
 
  i4mat_print ( N, N, a, "  The reordered matrix:" );
 
  return;
# undef N
}
//****************************************************************************80

void test048 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST048 tests INDEX_NEXT0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi = 3;
  int i;
  bool more;

  cout << "\n";
  cout << "TEST048\n";
  cout << "  INDEX_NEXT0 generates all indices of an\n";
  cout << "  array of given shape, with\n";
  cout << "  lower limit 1 and given upper limit.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";
  cout << "  Coordinate maximum HI =   " << hi << "\n";
 
  cout << "\n";
  cout << "  Index arrays:\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    index_next0 ( N, hi, a, more );

    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test049 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST049 tests INDEX_NEXT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi[N] = { 4, 2, 3 };
  int i;
  bool more;

  cout << "\n";
  cout << "TEST049\n";
  cout << "  INDEX_NEXT1 generates all indices of an\n";
  cout << "  array of given shape, with\n";
  cout << "  lower limit 1 and given upper limits.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";

  i4vec1_print ( N, hi, "  Coordinate maximum indices:" );
 
  cout << "\n";
  cout << "  Index arrays:\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    index_next1 ( N, hi, a, more );

    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test050 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST050 tests INDEX_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi[N] = { 11, -3, 1 };
  int i;
  int lo[N] = { 10, -5, 0 };
  bool more;

  cout << "\n";
  cout << "TEST050\n";
  cout << "  INDEX_NEXT2 generates all indices of an\n";
  cout << "  array of given shape with given\n";
  cout << "  lower and upper limits.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";
  cout << "\n";
  cout << "  Coordinate, Maximum Index\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(8) << i+1   << "  "
         << setw(8) << lo[i] << "  "
         << setw(8) << hi[i] << "\n";
  }
 
  cout << "\n";
  cout << "Index arrays:\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    index_next2 ( N, lo, hi, a, more );

    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test051 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST051 tests INDEX_RANK0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N] = { 3, 1, 2 };
  int hi = 3;
  int i;
  int rank;

  cout << "\n";
  cout << "TEST051\n";
  cout << "  INDEX_RANK0 ranks an index with\n";
  cout << "  lower limit 1 and given upper limit.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";
  cout << "\n";
  cout << "  Coordinate maximum Index = " << hi << "\n";
  cout << "\n";

  i4vec1_print ( N, a, "  The index array:" );

  rank = index_rank0 ( N, hi, a );

  cout << "\n";
  cout << "  The rank of this object is " << rank << "\n";

  return;
# undef N
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests INDEX_RANK1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N] = { 4, 1, 2 };
  int hi[N] = { 4, 2, 3 };
  int i;
  int rank;

  cout << "\n";
  cout << "TEST052\n";
  cout << "  INDEX_RANK1 ranks an index with\n";
  cout << "  lower limit 1 and given upper limits.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";
  cout << "\n";
  cout << "  Coordinate, Maximum Index\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(10) << i+1   << "  "
         << setw(10) << hi[i] << "\n";
  }
 
  i4vec1_print ( N, a, "  The index array:" );

  rank = index_rank1 ( N, hi, a );

  cout << "\n";
  cout << "  The rank of this object is " << rank << "\n";

  return;
# undef N
}
//****************************************************************************80

void test053 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST053 tests INDEX_RANK2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N] = { 1, 11, 5 };
  int hi[N] = { 2, 11, 6 };
  int i;
  int lo[N] = { 1, 10, 4 };
  int rank;

  cout << "\n";
  cout << "TEST053\n";
  cout << "  INDEX_RANK2 ranks an index with given\n";
  cout << "  lower and upper limits.\n";
  cout << "\n";
  cout << "  Number of index entries = " << N << "\n";
  cout << "\n";
  cout << "  Coordinate, Minimum index, Maximum Index\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(10) << i+1   << "  "
         << setw(10) << lo[i] << "  "
         << setw(10) << hi[i] << "\n";
  }
 
  i4vec1_print ( N, a, "  The index array:" );

  rank = index_rank2 ( N, lo, hi, a );

  cout << "\n";
  cout << "  The rank of this object is " << rank << "\n";

  return;
# undef N
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests INDEX_UNRANK0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi = 3;
  int i;
  int j;
  int maxrank;
  int rank;

  cout << "\n";
  cout << "TEST054\n";
  cout << "  INDEX_UNRANK0 unranks a multi-index.\n";
  cout << "\n";
  cout << "  The multi-index has dimension " << N << "\n";
  cout << "\n";
  cout << "  The upper limit is HI = " << hi << "\n";
  cout << "\n";
  cout << "  Rank, Multi-Index:\n";
  cout << "\n";
 
  maxrank = ( int ) pow ( ( double ) hi, N );

  for ( rank = 1; rank <= maxrank; rank++ )
  {
    index_unrank0 ( N, hi, rank, a );
    cout << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(6) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests INDEX_UNRANK1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi[N] = { 4, 2, 3 };
  int i;
  int j;
  int maxrank;
  int rank;

  cout << "\n";
  cout << "TEST055\n";
  cout << "  INDEX_UNRANK1 unranks a multi-index.\n";
  cout << "\n";
  cout << "  The multi-index has dimension " << N << "\n";

  i4vec1_print ( N, hi, "  The upper limits:" );

  cout << "\n";
  cout << "  Rank, Multi-Index:\n";
  cout << "\n";
 
  maxrank = i4vec_product ( N, hi );

  for ( rank = 1; rank <= maxrank; rank++ )
  {
    index_unrank1 ( N, hi, rank, a );
    cout << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(6) << a[i] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test056 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST056 tests INDEX_UNRANK2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int hi[N] = { 2, 11, 6 };
  int i;
  int j;
  int lo[N] = { 1, 10, 4 };
  int rank;

  cout << "\n";
  cout << "TEST056\n";
  cout << "  INDEX_UNRANK2 unranks a multi-index.\n";
  cout << "\n";
  cout << "  The multi-index has dimension " << N << "\n";
  cout << "\n";
  cout << "  The lower and upper limits are:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(10) << i     << "  "
         << setw(10) << lo[i] << "  "
         << setw(10) << hi[i] << "\n";
  }
  cout << "\n";
  cout << "  Rank, Multi-Index:\n";
  cout << "\n";
 
  rank = 7;

  index_unrank2 ( N, lo, hi, rank, a );
  cout << setw(3) << rank << "  ";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(6) << a[i] << "  ";
  }
  cout << "\n";
 
  return;
# undef N
}
//****************************************************************************80

void test057 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST057 tests INS_PERM and PERM_INS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
 
  int ins[N];
  int perm[N] = { 3, 5, 1, 4, 2 };
  int perm2[N];

  cout << "\n";
  cout << "TEST057\n";
  cout << "  PERM_INS computes the inversion sequence.\n";
  cout << "  INS_PERM recovers the permutation.\n";
  cout << "\n";

  i4vec1_print ( N, perm, "  The permutation:" );
  
  perm_ins ( N, perm, ins );

  i4vec1_print ( N, ins, "  The inversion sequence:" );

  ins_perm ( N, ins, perm2 );

  i4vec1_print ( N, perm2, "  The recovered permutation:" );
 
  return;
# undef N
}
//****************************************************************************80

void test063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST063 tests INVOLUTE_ENUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int s[N+1];

  cout << "\n";
  cout << "TEST063\n";
  cout << "  INVOLUTE_ENUM counts involutions;\n";
  cout << "\n";

  involute_enum ( N, s );

  cout << "\n";
  cout << "  N    # of involutions\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    cout << setw(10) << i    << "  "
         << setw(10) << s[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST064 test I4POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int i;
  int a[N];
  int iopt;
  int test;
  int val;
  int x0;

  cout << "\n";
  cout << "TEST064\n";
  cout << "  I4POLY converts between power sum, factorial\n";
  cout << "  and Taylor forms, and can evaluate a polynomial\n";
  cout << "\n";
 
  for ( test = 1; test <= 6; test++ )
  {
    if ( test == 1 )
    {
      iopt = -3;
    }
    else if ( test == 2 )
    {
      iopt = -2;
    }
    else if ( test == 3 )
    {
      iopt = -1;
      x0 = 2;
    }
    else if ( test == 4 )
    {
      iopt = 0;
      x0 = 2;
    }
    else if ( test == 5 )
    {
      iopt = 6;
      x0 = 2;
    }
    else if ( test == 6 )
    {
      iopt = 6;
      x0 = -2;
    }

    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
    a[4] = 0;
    a[5] = 1;

    if ( test == 1 )
    {
      i4vec1_print ( N, a, "  All calls have input A as follows:" );
    }
 
    i4poly ( N, a, x0, iopt, val );
 
    cout << "\n";
    cout << "  Option IOPT = " << iopt << "\n";

    if ( -1 <= iopt )
    {
      cout << "  X0 = " << x0 << "\n";
    }

    if ( iopt == -3 || iopt == -2 || iopt > 0 )
    {
      i4vec1_print ( N, a, "  Output array:" );
    }

    if ( iopt == -1 || iopt == 0 )
    {
      cout << "  Value = " << val << "\n";
    }
 
  }

  return;
# undef N
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests I4POLY_CYCLO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int phi[N_MAX+1];
  int n;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  I4POLY_CYCLO computes cyclotomic polynomials.\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    cout << "\n";
    cout << "  N = " << n << "\n";
    cout << "\n";

    i4poly_cyclo ( n, phi );

    i4poly_print ( n, phi, "  The cyclotomic polynomial:" );
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST066 tests I4POLY_DIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a[11];
  int b[11];
  int i;
  int na;
  int nb;
  int nq;
  int nr;
  int ntest = 2;
  int q[11];
  int r[11];
  int test;

  cout << "\n";
  cout << "TEST066\n";
  cout << "  I4POLY_DIV computes the quotient and\n";
  cout << "  remainder for polynomial division.\n";
  cout << "\n";
//
//  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
//     Quotient is 3+4*X+X**2, remainder is 0.
//
//  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
//     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
//
  for ( test = 1; test <= ntest; test++ )
  {
    if ( test == 1 )
    {
      na = 3;
      a[0] = -6;
      a[1] = -5;
      a[2] =  2;
      a[3] =  1;

      nb = 1;
      b[0] = -2;
      b[1] =  1;
    }
    else if ( test == 2 )
    {
      na = 4;
      a[0] = -2;
      a[1] =  5;
      a[2] =  2;
      a[3] =  3;
      a[4] =  1;
      nb = 2;
      b[0] = -3;
      b[1] =  1;
      b[2] =  1;
    }

    i4poly_print ( na, a, "  The polynomial to be divided, A:" );
    i4poly_print ( nb, b, "  The divisor polynomial, B:" );

    i4poly_div ( na, a, nb, b, nq, q, nr, r );
 
    i4poly_print ( nq, q, "  The quotient polynomial, Q:" );
    i4poly_print ( nr, r, "  The remainder polynomial, R:" );
  }

  return;
}
//****************************************************************************80

void test067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067 tests I4POLY_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAXN 5

  int a[MAXN+1];
  int b[MAXN+1];
  int c[MAXN+1];
  int na;
  int nb;
  int ntest = 2;
  int test;

  cout << "\n";
  cout << "TEST067\n";
  cout << "  I4POLY_MUL multiplies two polynomials.\n";
  cout << "\n";
//
//  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
//  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
//
  for ( test = 1; test <= ntest; test++ )
  {
    if ( test == 1 )
    {
      na = 1;
      a[0] =  1;
      a[1] =  1;
      nb = 1;
      b[0] =  1;
      b[1] = -1;
    }
    else if ( test == 2 )
    {
      na = 2;
      a[0] =  1;
      a[1] =  2;
      a[2] =  3;
      nb = 1;
      b[0] =  1;
      b[1] = -2;
    }

    i4poly_mul ( na, a, nb, b, c );

    i4poly_print ( na, a, "  The factor A:" );

    i4poly_print ( nb, b, "  The factor B:" );

    i4poly_print ( na+nb, c, "  The product C = A*B:" );

  }

  return;
# undef MAXN
}
//****************************************************************************80

void test0675 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0675 tests I4VEC_DESCENDS;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int *a;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST0675\n";
  cout << "  I4VEC_DESCENDS is true if an integer vector decreases.\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    a = i4vec_uniform ( N, 1, N, seed );

    i4vec1_print ( N, a, "  The integer array to search:" );
 
    if ( i4vec_descends ( N, a ) )
    {
      cout << "  The preceding vector is descending.\n";
    }
    else
    {
      cout << "  The preceding vector is not descending.\n";
    }
    delete [] a;
  }

  return;
# undef N
}
//****************************************************************************80

void test068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST068 tests I4VEC_FRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int afrac;
  int i;
  int k;
  int seed;

  cout << "\n";
  cout << "TEST068\n";
  cout << "  I4VEC_FRAC: K-th smallest integer vector entry.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 1, 2*N, seed );

  i4vec1_print ( N, a, "  The integer array to search:" );

  cout << "\n";
  cout << "     K   K-th smallest\n";
  cout << "\n";

  for ( k = 1; k <= N; k++ )
  {
    afrac = i4vec_frac ( N, a, k );

    cout << setw(6) << k        << "  "
         << setw(6) << afrac    << "\n";

  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test0683 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0683 tests I4VEC_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int aval;
  int first;
  int seed;

  cout << "\n";
  cout << "TEST0683\n";
  cout << "  I4VEC_INDEX returns the index of the first occurrence\n";
  cout << "  of a given value in an integer vector.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 1, N/2, seed );

  aval = a[N/2];

  i4vec1_print ( N, a, "  The integer array to search:" );

  first = i4vec_index ( N, a, aval );

  cout << "\n";
  cout << "  The value searched for is " << aval << "\n";
  cout << "  The index of first occurrence is " << first << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test0685 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0685 tests I4VEC_MAXLOC_LAST;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int last;
  int seed;

  cout << "\n";
  cout << "TEST0685\n";
  cout << "  I4VEC_MAXLOC_LAST: index of the last maximal\n";
  cout << "  entry in an integer vector.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 1, N/4, seed );

  i4vec1_print ( N, a, "  The integer array to search:" );
 
  last = i4vec_maxloc_last ( N, a );

  cout << "\n";
  cout << "  Index of last maximal entry is " << last << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test0686 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0686 tests I4VEC_PAIRWISE_PRIME;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int *a;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST0686\n";
  cout << "  I4VEC_PAIRWISE_PRIME is true if an integer vector\n";
  cout << "  is pairwise prime.\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    a = i4vec_uniform ( N, 1, N, seed );

    i4vec1_print ( N, a, "  The array to check:" );
 
    if ( i4vec_pairwise_prime ( N, a ) )
    {
      cout << "  The preceding vector is pairwise prime.\n";
    }
    else
    {
      cout << "  The preceding vector is not pairwise prime.\n";
    }
    delete [] a;
  }

  return;
# undef N
}
//****************************************************************************80

void test0687 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0687 tests I4VEC_REVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int *a;
  int seed;

  cout << "\n";
  cout << "TEST0687\n";
  cout << "  I4VEC_REVERSE reverses an integer vector.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 1, N, seed );

  i4vec1_print ( N, a, "  The integer array:" );

  i4vec_reverse ( N, a );

  i4vec1_print ( N, a, "  The reversed integer array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test0688 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0688 tests I4VEC_SORT_BUBBLE_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int seed;

  cout << "\n";
  cout << "TEST0688\n";
  cout << "  I4VEC_SORT_BUBBLE_A ascending sorts an integer vector\n";
  cout << "    using bubble sort.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 0, 3*N, seed );

  i4vec1_print ( N, a, "  Unsorted array:" );

  i4vec_sort_bubble_a ( N, a );

  i4vec1_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test0689 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0689 tests I4VEC_SORT_HEAP_INDEX_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "TEST0689\n";
  cout << "  I4VEC_SORT_HEAP_INDEX_D descending index-sorts\n";
  cout << "    an integer vector using heap sort.\n";
  cout << "\n";

  seed = 123456789;

  a = i4vec_uniform ( N, 0, 3*N, seed );

  i4vec1_print ( N, a, "  Unsorted array:" );

  indx = i4vec_sort_heap_index_d ( N, a );

  cout << "\n";
  cout << "     I  INDX[I]  A[INDX[I]-1]\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(6) << i          << "  "
         << setw(6) << indx[i]    << "  "
         << setw(6) << a[indx[i]] << "\n";
  }

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test06895 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06895 tests INVERSE_MOD_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2009
//
//  Author:
//
//    John Burkardt
//
{
  int b;
  int n;
  int y;
  int z;

  cout << "\n";
  cout << "TEST06895\n";
  cout << "  INVERSE_MOD_N seeks Y, the inverse of B mod N,\n";
  cout << "  so that mod ( B * Y, N ) = 1, but returns 0\n";
  cout << "  if the inverse does not exist.\n";

  cout << "\n";
  cout << "     B     N     Y     Z = ( ( B * Y ) % N )\n";
  cout << "\n";

  for ( n = 1; n <= 10;  n++ )
  {
    for ( b = 1; b < n; b++ )
    {
      y = inverse_mod_n ( b, n );
      z = ( ( b * y ) % n );
      cout << "  " << setw(2) << b
           << "  " << setw(2) << n
           << "  " << setw(2) << y
           << "  " << setw(2) << z << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test069 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST069 tests JFRAC_TO_RFRAC and RFRAC_TO_JFRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define MAXM 10

  int i;
  int m;
  double p[MAXM];
  double q[MAXM];
  double r[MAXM];
  double s[MAXM];
  int seed;
//
//  Generate the data, but force Q(M+1) to be 1.  
//  That will make it easier to see that the two operations are inverses
//  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
//  a scaling in which Q(M+1) is 1.
//
  seed = 123456789;
  m = 6;
  r8vec_uniform_01 ( m, seed, p );
  r8vec_uniform_01 ( m + 1, seed, q );

  for ( i = 0; i < m; i++ )
  {
    q[i] = q[i] / q[m];
  }
  q[m] = 1.0;

  cout << "\n";
  cout << "TEST069\n";
  cout << "  RFRAC_TO_JFRAC converts a rational polynomial\n";
  cout << "  fraction to a J fraction.\n";
  cout << "  JFRAC_TO_RFRAC converts a J fraction\n";
  cout << "  to a rational polynomial fraction.\n";
  cout << "\n";
  cout << "  The original rational polynomial coefficients:\n";
  cout << "\n";

  for ( i = 0; i < m; i++ )
  {
    cout << setw(14) << p[i] << "  ";
  }
  cout << "\n";

  for ( i = 0; i < m+1; i++ )
  {
    cout << setw(14) << q[i] << "  ";
  }
  cout << "\n";
 
  rfrac_to_jfrac ( m, p, q, r, s );
 
  cout << "\n";
  cout << "  The J fraction coefficients:\n";
  cout << "\n";

  for ( i = 0; i < m; i++ )
  {
    cout << setw(14) << r[i] << "  ";
  }
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << setw(14) << s[i] << "  ";
  }
  cout << "\n";
 
  jfrac_to_rfrac ( m, r, s, p, q );

  cout << "\n";
  cout << "  The recovered rational polynomial:\n";
  cout << "\n";

  for ( i = 0; i < m; i++ )
  {
    cout << setw(14) << p[i] << "  ";
  }
  cout << "\n";

  for ( i = 0; i < m+1; i++ )
  {
    cout << setw(14) << q[i] << "  ";
  }
  cout << "\n";

  return;
# undef MAXM
}
//****************************************************************************80

void test070 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST070 tests JOSEPHUS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int k;
  int m;
  int n;
  int x;

  cout << "\n";
  cout << "TEST070\n";
  cout << "  JOSEPHUS solves Josephus problems.\n";
  cout << "\n";
  cout << "    N    M    K	 X\n";
  cout << "\n";

  m = 3;
  n = 41;
  k = 41;
  x = josephus ( n, m, k );

  cout << setw(5) << n << "  "
       << setw(5) << m << "  "
       << setw(5) << k << "  "
       << setw(5) << x << "\n";

  m = -38;
  n = 41;
  k = 41;
  x = josephus ( n, m, k );

  cout << setw(5) << n << "  "
       << setw(5) << m << "  "
       << setw(5) << k << "  "
       << setw(5) << x << "\n";

  m = 3;
  n = 41;
  k = 40;
  x = josephus ( n, m, k );

  cout << setw(5) << n << "  "
       << setw(5) << m << "  "
       << setw(5) << k << "  "
       << setw(5) << x << "\n";

  m = 2;
  n = 64;
  k = 64;
  x = josephus ( n, m, k );

  cout << setw(5) << n << "  "
       << setw(5) << m << "  "
       << setw(5) << k << "  "
       << setw(5) << x << "\n";

  m = 2;
  n = 1000;
  k = 1000;
  x = josephus ( n, m, k );

  cout << setw(5) << n << "  "
       << setw(5) << m << "  "
       << setw(5) << k << "  "
       << setw(5) << x << "\n";

  return;
}
//****************************************************************************80

void test071 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST071 tests KSUB_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  bool more;
  int n = 5;
  int rank;

  for ( i = 0; i < K; i++ )
  {
    a[i] = 0;
  }

  cout << "\n";
  cout << "TEST071\n";
  cout << "  KSUB_NEXT generates all K subsets of an N set\n";
  cout << "  in lexicographic order.\n";
  cout << "\n";

  more = false;
  rank = 0;
 
  for ( ; ; )
  {
    ksub_next ( n, K, a, more );

    rank = rank + 1;

    cout << setw(4) << rank << "    ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef K
}
//****************************************************************************80

void test072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST072 tests KSUB_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  int i_in;
  int i_out;
  int j;
  bool more;
  int n = 5;
  int rank;

  cout << "\n";
  cout << "TEST072\n";
  cout << "  KSUB_NEXT2 generates the next K subset of an\n";
  cout << "  N set by the revolving door method.\n";
  cout << "\n";
  cout << "Rank  Subset  Added  Removed\n";
  cout << "\n";
//
//  KSUB_NEXT2 does not have a good way of stopping.  
//  We will save the starting subset, and stop when the
//  new subset is the same as the starting one.
//
  i_in = 0;
  i_out = 0;
  rank = 0;
 
  i4vec_indicator ( K, a );
 
  for ( ; ; )
  {
    rank = rank + 1;
    cout << setw(4) << rank << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(2) << a[i] << "  ";
    }
    cout << "   ";
    cout << setw(2) << i_in  << "  ";
    cout << setw(2) << i_out << "\n";
 
    ksub_next2 ( n, K, a, i_in, i_out );
 
    more = false;

    for ( i = 1; i <= K; i++ )
    {
      if ( a[i-1] != i )
      {
        more = true;
      }
    }

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef K
}
//****************************************************************************80

void test073 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST073 tests KSUB_NEXT3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  int i_in = 0;
  int i_out = 0;
  bool more;
  int n = 5;
  int rank;

  cout << "\n";
  cout << "TEST073\n";
  cout << "  KSUB_NEXT3 generates all K subsets of an N set\n";
  cout << "  using the revolving door method.\n";
  cout << "\n";
  cout << "Rank    Subset  Added Removed\n";
  cout << "\n";

  rank = 0;
  more = false;
 
  for ( ; ; )
  {
    ksub_next3 ( n, K, a, more, i_in, i_out );

    rank = rank + 1;
    cout << setw(4) << rank << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(2) << a[i] << "  ";
    }
    cout << "   ";
    cout << setw(2) << i_in  << "  ";
    cout << setw(2) << i_out << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;
# undef K
}
//****************************************************************************80

void test074 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST074 tests KSUB_NEXT4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  bool done;
  int i;
  int n = 5;
  int rank;

  cout << "\n";
  cout << "TEST074\n";
  cout << "  KSUB_NEXT4 generates K subsets of an N set.\n";
  cout << "  N = " << n << "\n";
  cout << "  K = " << K << "\n";
  cout << "\n";
  cout << "Rank    Subset\n";
  cout << "\n";

  done = true;
  rank = 0;
 
  for ( ; ; )
  {
    ksub_next4 ( n, K, a, done );
 
    if ( done )
    {
      break;
    }

    rank = rank + 1;
    cout << setw(4) << rank << "  ";
    cout << "  ";
    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

  }
 
  return;
# undef K
}
//****************************************************************************80

void test075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST075 tests KSUB_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  int j;
  int n = 5;
  int seed;

  cout << "\n";
  cout << "TEST075\n";
  cout << "  KSUB_RANDOM generates a random K subset of an N set.\n";
  cout << "  Set size is N =    " << n << "\n";
  cout << "  Subset size is K = " << K << "\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    ksub_random ( n, K, seed, a );
    for ( j = 0; j < K; j++ )
    {
      cout << "  " << setw(3) << a[j];
    }
    cout << "\n";
  }
 
  return;
# undef K
}
//****************************************************************************80

void test076 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests KSUB_RANDOM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int i;
  int j;
  int n = 5;
  int seed;

  cout << "\n";
  cout << "TEST076\n";
  cout << "  KSUB_RANDOM2 generates a random K subset of an N set.\n";
  cout << "  Set size is N =    " << n << "\n";
  cout << "  Subset size is K = " << K << "\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    ksub_random2 ( n, K, seed, a );
    for ( j = 0; j < K; j++ )
    {
      cout << "  " << setw(3) << a[j];
    }
    cout << "\n";
  }
 
  return;
# undef K
}
//****************************************************************************80

void test077 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST077 tests KSUB_RANDOM3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define K 3
# define N 5

  int a[N];
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST077\n";
  cout << "  KSUB_RANDOM3 generates a random K subset of an N set.\n";
  cout << "  Set size is N =    " << N << "\n";
  cout << "  Subset size is K = " << K << "\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    ksub_random3 ( N, K, seed, a );
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(3) << a[j];
    }
    cout << "\n";
  }
 
  return;
# undef K
# undef N
}
//****************************************************************************80

void test0771 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0771 tests KSUB_RANDOM4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define K 3
# define N 5

  int a[N];
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST0771\n";
  cout << "  KSUB_RANDOM4 generates a random K subset of an N set.\n";
  cout << "  Set size is N =    " << N << "\n";
  cout << "  Subset size is K = " << K << "\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    ksub_random4 ( N, K, seed, a );
    for ( j = 0; j < K; j++ )
    {
      cout << "  " << setw(3) << a[j];
    }
    cout << "\n";
  }
 
  return;
# undef K
# undef N
}
//****************************************************************************80

void test07715 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07715 tests KSUB_RANDOM5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int i;
  int j;
  int k = 5;
  int n = 52;
  int seed;

  cout << "\n";
  cout << "TEST07715\n";
  cout << "  KSUB_RANDOM5 generates a random K subset of an N set.\n";
  cout << "  Set size is N =    " << n << "\n";
  cout << "  Subset size is K = " << k << "\n";
  cout << "\n";
 
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    a = ksub_random5 ( n, k, seed );
    for ( j = 0; j < k; j++ )
    {
      cout << "  " << setw(3) << a[j];
    }
    cout << "\n";
    delete [] a;
  }
 
  return;
}
//****************************************************************************80

void test0772 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0772 tests KSUB_RANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
# define K 3

  int a[K] = { 1, 3, 5 };
  int i;
  int rank;

  cout << "\n";
  cout << "TEST0772\n";
  cout << "  KSUB_RANK: determine the rank of a K subset of an N set.\n";
  cout << "\n";
  cout << "  For N = " << N << "\n";
  cout << "  and K = " << K << "\n";
  cout << "  the subset is:\n";
  cout << "\n";

  for ( i = 0; i < K; i++ )
  {
    cout << setw(4) << a[i];
  }
  cout << "\n";

  ksub_rank ( K, a, rank );

  cout << "\n";
  cout << "  The computed rank is " << rank << "\n";

  return;

# undef N
# undef K
}
//****************************************************************************80

void test0773 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0773 tests KSUB_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int i;
  int a[K];
  int n = 5;
  int rank;

  rank = 8;

  cout << "\n";
  cout << "TEST0773\n";
  cout << "  KSUB_UNRANK: find the K-subset of an N set\n";
  cout << "  of a given rank.\n";
  cout << "\n";
  cout << "  For N = " << n << "\n";
  cout << "  and K = " << K << "\n";
  cout << "  and the desired rank is " << rank << "\n";

  ksub_unrank ( K, rank, a );

  cout << "\n";
  cout << "  The subset of the given rank is:\n";
  cout << "\n";

  for ( i = 0; i < K; i++ )
  {
    cout << setw(4) << a[i];
  }
  cout << "\n";

  return;
# undef K
}
//****************************************************************************80

void test078 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST078 tests MATRIX_PRODUCT_OPT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int cost;
  int i;
  int order[N-1];
  int rank[N+1] = { 4, 2, 3, 1, 2, 2, 3 };

  cout << "\n";
  cout << "TEST078\n";
  cout << "  MATRIX_PRODUCT_OPT seeks the optimal order\n";
  cout << "  for a chain of matrix products.\n";
  cout << "\n";
  cout << "  Matrix ranks:\n";
  cout << "\n";
  cout << "   I    R    C\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(5) << i         << "  "
         << setw(5) << rank[i]   << "  "
         << setw(5) << rank[i+1] << "\n";
  }

  matrix_product_opt ( N, rank, cost, order );

  cout << "\n";
  cout << "  Optimal cost is " << cost << "\n";

  i4vec1_print ( N-1, order, "  Ordering:" );

  return;
# undef N
}
//****************************************************************************80

void test079 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST079 tests MOEBIUS_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ih[N*N] = {
    0,0,1,1,0,0,0,0,0,0,0, 
    0,0,0,0,0,0,0,1,0,0,0, 
    0,1,0,0,0,0,0,0,0,0,0, 
    0,1,0,0,0,0,0,0,0,0,0, 
    0,0,0,1,0,0,0,0,0,0,0, 
    1,0,0,0,1,0,0,0,1,0,0, 
    0,0,0,0,0,1,0,0,0,1,1, 
    0,0,0,0,0,0,0,0,0,0,0, 
    0,0,1,1,0,0,0,0,0,0,0, 
    1,0,0,0,0,0,0,0,1,0,0, 
    0,0,0,0,0,0,0,0,1,0,0 };
  int j;
  int matrix[N*N];

  cout << "\n";
  cout << "TEST079\n";
  cout << "  MOEBIUS_MATRIX computes the Moebius matrix.\n";
 
  i4mat_print ( N, N, ih, "  The input matrix:" );

  moebius_matrix ( N, ih, matrix );
 
  i4mat_print ( N, N, matrix, "  The Moebius matrix:" );
 
  return;
# undef N
}
//****************************************************************************80

void test0795 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0795 tests MONOMIAL_COUNT and MONOMIAL_COUNTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int *counts;
  int degree;
  int degree_max = DEGREE_MAX;
  int dim;
  int total;
  int total2;

  cout << "\n";
  cout << "TEST0795\n";
  cout << "  MONOMIAL_COUNT counts the number of monomials of\n";
  cout << "  degrees 0 through DEGREE_MAX in a space of dimension DIM.\n";
  cout << "  MONOMIAL_COUNTS provides individual counts.\n";

  for ( dim = 1; dim <= 6; dim++ )
  {
    counts = monomial_counts ( degree_max, dim );
    total = monomial_count ( degree_max, dim );

    cout << "\n";
    cout << "  DIM = " << dim << "\n";
    cout << "\n";

    for ( degree = 0; degree <= degree_max; degree++ )
    {
      cout << "  " << setw(8) << degree
           << "  " << setw(8) << counts[degree] << "\n";
    }

    total2 = 0;
    for ( degree = 0; degree <= degree_max; degree++ )
    {
      total2 = total2 + counts[degree];
    }
    cout << "\n";
    cout << "     Total"
         << "  " << setw(8) << total2
         << "  " << setw(8) << total << "\n";

    delete [ ] counts;
  }

  return;
}
//****************************************************************************80

void test080 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST080 tests MORSE_THUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  int i;
  int ihi;
  int ilo;
  int s[N+1];

  cout << "\n";
  cout << "TEST080\n";
  cout << "  MORSE_THUE computes the Morse-Thue numbers.\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    s[i] = morse_thue ( i );
  }

  for ( ilo = 0; ilo <= N; ilo = ilo + 10 )
  {
    cout << "  ";
    ihi = i4_min ( ilo + 9, N );
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(1) << s[i];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test081 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST081 tests MULTINOMIAL_COEF1 and MULTINOMIAL_COEF2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXFACTOR 5

  int factor[MAXFACTOR];
  int i;
  int j;
  int n;
  int ncomb1;
  int ncomb2;
  int nfactor;

  cout << "\n";
  cout << "TEST081\n";
  cout << "  MULTINOMIAL_COEF1 computes multinomial\n";
  cout << "    coefficients using the Gamma function;\n";
  cout << "  MULTINOMIAL_COEF2 computes multinomial\n";
  cout << "    coefficients directly.\n";

  cout << "\n";
  cout << "  Line 10 of the BINOMIAL table:\n";
  cout << "\n";

  n = 10;
  nfactor = 2;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;
    factor[1] = n - i;

    ncomb1 = multinomial_coef1 ( nfactor, factor );

    ncomb2 = multinomial_coef2 ( nfactor, factor );

    cout << setw(4) << factor[0] << "  "
         << setw(4) << factor[1] << "  "
         << "    "
         << setw(5) << ncomb1 << "  "
         << setw(5) << ncomb2 << "\n";

  }

  cout << "\n";
  cout << "  Level 5 of the TRINOMIAL coefficients:\n";

  n = 5;
  nfactor = 3;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;

    cout << "\n";

    for ( j = 0; j <= n - factor[0]; j++ )
    {
      factor[1] = j;
      factor[2] = n - factor[0] - factor[1];

      ncomb1 = multinomial_coef1 ( nfactor, factor );

      ncomb2 = multinomial_coef2 ( nfactor, factor );

      cout << setw(4) << factor[0] << "  "
           << setw(4) << factor[1] << "  "
           << setw(4) << factor[2] << "  "
           << "    "
           << setw(5) << ncomb1 << "  "
           << setw(5) << ncomb2 << "\n";

    }

  }

  return;
# undef MAXFACTOR
}
//****************************************************************************80

void test0813 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0813 tests MULTIPERM_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int counts[N];
  int i;
  int k;
  int number;
  int seed = 123456789;
  int test;
  int test_num = 5;
  
  cout << "\n";
  cout << "TEST0813\n";
  cout << "  MULTIPERM_ENUM enumerates multipermutations.\n";
  cout << "\n";
  cout << "  N is the number of objects to be permuted.\n";
  cout << "  K is the number of distinct types of objects.\n";
  cout << "  COUNTS is the number of objects of each type.\n";
  cout << "  NUMBER is the number of multipermutations.\n";
  cout << "\n";
  cout << "  Number       N       K       Counts(1:K)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {  
    k = i4_uniform ( 1, N, seed );

    compnz_random ( N, k, seed, counts );

    number = multiperm_enum ( N, k, counts );

    cout << "  " << setw(6) << number
         << "  " << setw(6) << N
         << "  " << setw(6) << k;
    for ( i = 0; i < k; i++ )
    {
      cout << "  " << setw(4) << counts[i];
    }
    cout << "\n";
  }
  
  return;
# undef N
}
//****************************************************************************80

void test0815 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0815 tests MULTIPERM_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int a[N] = { 1, 2, 2, 3, 3, 3 };
  int i;
  bool more;
  int tally;
  
  cout << "\n";
  cout << "TEST0815\n";
  cout << "  MULTIPERM_NEXT computes multipermutations in\n";
  cout << "  lexical order.\n";
  cout << "\n";

  tally = 0;
  more = true;
 
  while ( more )
  {
    tally = tally + 1;

    cout << "  " << setw(4) << tally;
    for ( i = 0; i < N; i++ )
    {
      cout << "  " << setw(2) << a[i];
    }
    cout << "\n";

    multiperm_next ( N, a, more );
  }
 
  return;
# undef N
}
//****************************************************************************80

void test082 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST082 tests NETWORK_FLOW_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 6
# define NEDGE 20

  int cut[NNODE];
  int i;
  int icpflo[2*NEDGE] = {
    3,0,7,0,2,0,5,0,4,0,1,0,4,0,2,0,8,0,3,0, 
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  int iendpt[2*NEDGE] = {
    1,2, 1,3, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, 4,6, 5,6,
    2,1, 3,1, 3,2, 4,2, 5,2, 4,3, 5,3, 5,4, 6,4, 6,5 };
  int j;
  int node_flow[NNODE];
  int sink = 6;
  int source = 1;

  cout << "\n";
  cout << "TEST082\n";
  cout << "  NETWORK_FLOW_MAX finds the maximum flow on a network.\n";
  cout << "\n";
  cout << "  The source is node " << source << "\n";
  cout << "  The sink is node   " << sink   << "\n";
  cout << "\n";
  cout << "  Endpoint array:\n";
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << iendpt[0+j*2];
  }
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << iendpt[1+j*2];
  }
  cout << "\n";

  cout << "\n";
  cout << "  Input edge capacity array:\n";
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << icpflo[0+j*2];
  }
  cout << "\n";

  cout << "\n";
  cout << "TEST082\n";
  cout << "  Cancel this test, the routine is not ready!\n";
  return;

  network_flow_max ( NNODE, NEDGE, iendpt, icpflo, source, 
    sink, cut, node_flow );

  cout << "\n";
  cout << "  Reordered endpoint array:\n";
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << iendpt[0+j*2];
  }
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << iendpt[1+j*2];
  }
  cout << "\n";

  cout << "\n";
  cout << "  Output edge capacity//flow array:\n";
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << icpflo[0+j*2];
  }
  cout << "\n";
  for ( j = 0; j < NEDGE; j++ )
  {
    cout << setw(3) << icpflo[1+j*2];
  }
  cout << "\n";

  i4vec1_print ( NNODE, cut, "  Minimal node cut vector:" );

  i4vec1_print ( NNODE, node_flow, "  Nodal flow vector:" );

  return;
# undef NNODE
# undef NEDGE
}
//****************************************************************************80

void test083 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST083 tests NIM_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 32

  int i;
  unsigned int i1;
  int i1vec[N];
  int i2;
  int i2vec[N];
  unsigned int i3;
  int i3vec[N];
  int ihi = 1000;
  int ilo = 0;
  int ntest = 5;
  int seed;

  cout << "\n";
  cout << "TEST083\n";
  cout << "  NIM_SUM computes the Nim sum of two integers.\n";
  cout << "\n";
  cout << "    I    J    Nim(I+J)\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= ntest; i++ )
  {
    i1 = ( unsigned int ) i4_uniform ( ilo, ihi, seed );
    ui4_to_ubvec ( i1, N, i1vec );

    i2 = ( unsigned int ) i4_uniform ( ilo, ihi, seed );
    ui4_to_ubvec ( i2, N, i2vec );

    i3 = nim_sum ( i1, i2 );
    ui4_to_ubvec ( i3, N, i3vec );

    cout << "\n";
    cout << "  I1, I2, I3 in decimal:\n";
    cout << "\n";

    cout                  << "  "
         << setw(5) << i1 << "\n";
    cout                  << "  "
        << setw(5) << i2 << "\n";
    cout                  << "  "
        << setw(5) << i3 << "\n";

    cout << "\n";
    cout << "  I1, I2, I3 in binary:\n";
    cout << "\n";

    bvec_print ( N, i1vec, " " );
    bvec_print ( N, i2vec, " " );
    bvec_print ( N, i3vec, " " );

  }

  return;
# undef N
}
//****************************************************************************80

void test0835 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0835 tests PADOVAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int p[N];

  cout << "\n";
  cout << "TEST0835\n";
  cout << "  PADOVAN computes the Padovan numbers.\n";
  cout << "\n";

  padovan ( N, p );

  cout << "\n";
  cout << "   N    P(N)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                    << "  "
         << setw(4) << i    << "  "
         << setw(6) << p[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test084 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST084 tests PELL_BASIC and PELL_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int q;
  int r;
  int x0;
  int x1;
  int y0;
  int y1;

  cout << "\n";
  cout << "TEST084\n";
  cout << "  PELL_BASIC solves the basic Pell equation.\n";
  cout << "  PELL_NEXT computes the next solution.\n";
  cout << "\n";
  cout << "       D       X        Y         R\n";
  cout << "\n";

  for ( d = 2; d <= 20; d++ )
  {
    i4_sqrt ( d, q, r );

    if ( r != 0 )
    {
      pell_basic ( d, x0, y0 );

      r = x0 * x0 - d * y0 * y0;

      cout << setw(9) << d  << "  "
           << setw(9) << x0 << "  "
           << setw(9) << y0 << "  "
           << setw(9) << r  << "\n";

      pell_next ( d, x0, y0, x0, y0, x1, y1 );

      r = x1 * x1 - d * y1 * y1;

      cout << "         "   << "  "
           << setw(9) << x1 << "  "
           << setw(9) << y1 << "  "
           << setw(9) << r  << "\n";

    }

  }

  return;
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests PENT_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  PENT_ENUM counts points in pentagons.\n";
  cout << "\n";
  cout << "   N  Pent(N)\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    cout << setw(4) << i               << "  "
         << setw(6) << pent_enum ( i ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test086 ( )

//****************************************************************************80
// 
//  Purpose:
//
//    TEST086 tests PERM_BREAK_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int break_count;
  int p[N] = { 4, 5, 2, 1, 6, 3 };

  cout << "\n";
  cout << "TEST086\n";
  cout << "  PERM_BREAK_COUNT counts the breaks in a permutation.\n";
 
  perm_print ( N, p, "  The permutation:" );
 
  break_count = perm_break_count ( N, p );

  cout << "\n";
  cout << "  The number of breaks is " << break_count << "\n";

  return;
# undef N
}
//****************************************************************************80

void test087 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST087 tests PERM_CANON_TO_CYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int p1[N] = { 4, 5, 2, 1, 6, 3 };
  int p2[N];

  cout << "\n";
  cout << "TEST087\n";
  cout << "  PERM_CANON_TO_CYCLE converts a permutation from\n";
  cout << "  canonical to cycle form.\n";
 
  perm_print ( N, p1, "  The permutation in canonical form:");
 
  perm_canon_to_cycle ( N, p1, p2 );

  perm_print ( N, p2, "  The permutation in cycle form:" );
 
  return;
# undef N
}
//****************************************************************************80

void test088 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST088 tests PERM_CYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  int i;
  int iopt = 1;
  int isgn;
  int ncycle;
  int p[N] = { 2,3,9,6,7,8,5,4,1 };

  cout << "\n";
  cout << "TEST088\n";
  cout << "  PERM_CYCLE analyzes a permutation.\n";

  perm_print ( N, p, "  The permutation:" );

  perm_cycle ( N, p, isgn, ncycle, iopt );

  cout << "\n";
  cout << "  NCYCLE = " << ncycle << "\n";
  cout << "  ISGN =   " << isgn   << "\n";

  perm_print ( N, p, "  The permutation in cycle form:" );

  return;

# undef N
}
//****************************************************************************80

void test089 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST089 tests PERM_CYCLE_TO_CANON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int p1[N] = { -6, 3, 1, -5, 4, -2 };
  int p2[N];

  cout << "\n";
  cout << "TEST089\n";
  cout << "  PERM_CYCLE_TO_CANON converts a permutation from\n";
  cout << "  cycle to canonical form.\n";
 
  perm_print ( N, p1, "  The permutation in cycle form:" );
 
  perm_cycle_to_canon ( N, p1, p2 );

  perm_print ( N, p2, "  The permutation in canonical form:" );
 
  return;
# undef N
}
//****************************************************************************80

void test090 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST090 tests PERM_CYCLE_TO_INDEX and PERM_INDEX_TO_CYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  int p1[N] = { 2,3,9,6,7,8,5,4,1 };
  int p2[N];
  int p3[N];

  cout << "\n";
  cout << "TEST090\n";
  cout << "  PERM_CYCLE_TO_INDEX converts a permutation from\n";
  cout << "    cycle to standard index form.\n";
  cout << "  PERM_INDEX_TO_CYCLE converts a permutation from\n";
  cout << "    standard index to cycle form.\n";
 
  perm_print ( N, p1, "  The standard index form permutation:" );
 
  perm_index_to_cycle ( N, p1, p2 );

  perm_print ( N, p2, "  The permutation in cycle form:" );

  perm_cycle_to_index ( N, p2, p3 );
 
  perm_print ( N, p3, "  The standard index form permutation:" );
 
  return;
# undef N
}
//****************************************************************************80

void test091 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST091 tests PERM_DISTANCE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int k11;
  int k12;
  int k13;
  int k21;
  int k23;
  int p1[N];
  int p2[N];
  int p3[N];
  int seed;

  cout << "\n";
  cout << "TEST091\n";
  cout << "  PERM_DISTANCE computes the Ulam metric distance\n";
  cout << "  between two permutations.\n";

  seed = 123456789;

  perm_random3 ( N, seed, p1 );
  perm_random3 ( N, seed, p2 );
  perm_random3 ( N, seed, p3 );

  perm_print ( N, p1, "  Permutation P1" );
  perm_print ( N, p2, "  Permutation P2" );
  perm_print ( N, p3, "  Permutation P3" );

  k11 = perm_distance ( N, p1, p1 );
  k12 = perm_distance ( N, p1, p2 );
  k21 = perm_distance ( N, p2, p1 );
  k13 = perm_distance ( N, p1, p3 );
  k23 = perm_distance ( N, p2, p3 );

  cout << "\n";
  cout << "  K(P1,P1) should be 0.\n";
  cout << "\n";
  cout << "  K(P1,P1) = " << k11 << "\n";
  cout << "\n";
  cout << "  K(P1,P2) should equal K(P2,P1).\n";
  cout << "\n";
  cout << "  K(P1,P2) = " << k12 << "\n";
  cout << "  K(P2,P1) = " << k21 << "\n";
  cout << "\n";
  cout << "  K(P1,P2) + K(P2,P3) >= K(P1,P3).\n";
  cout << "\n";
  cout << "  K(P1,P3) = " << k13 << "\n";
  cout << "  K(P1,P2) = " << k12 << "\n";
  cout << "  K(P2,P3) = " << k23 << "\n";
  cout << "  K(P1,P2) + K(P2,P3) = " << k12 + k23 << "\n";

  return;
# undef N
}
//****************************************************************************80

void test092 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST092 tests PERM_FIXED_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int fnm;
  int m;
  int n = 10;

  cout << "\n";
  cout << "TEST092\n";
  cout << "  PERM_FIXED_ENUM enumerates the permutations\n";
  cout << "  of N objects that leave M unchanged.\n";
  cout << "\n";
  cout << "  For this test, N = " << n << "\n";
  cout << "\n";
  cout << "    M    F(N,M)\n";
  cout << "\n";

  for ( m = 0; m <= n; m++ )
  {
    fnm = perm_fixed_enum ( n, m );

    cout                   << "  "
         << setw(3) << m   << "  "
         << setw(8) << fnm << "\n";
  }

  return;
}
//****************************************************************************80

void test093 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST093 tests PERM_ASCEND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  int length;
  int p[N] = { 2,3,9,6,7,8,5,4,1 };
  int subseq[N];

  cout << "\n";
  cout << "TEST093\n";
  cout << "  PERM_ASCEND determines the length of the longest\n";
  cout << "  increasing subsequence in a permutation.\n";

  perm_print ( N, p, "  The permutation:" );

  perm_ascend ( N, p, length, subseq );

  cout << "\n";
  cout << "  The length of the longest increasing subsequence is " << 
    length << "\n";

  i4vec1_print ( length, subseq, "  A longest increasing subsequence:" );

  return;
# undef N
}
//****************************************************************************80

void test094 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST094 tests PERM_INVERSE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int p[N] = { 4, 3, 5, 1, 7, 6, 2 };

  cout << "\n";
  cout << "TEST094\n";
  cout << "  PERM_INVERSE inverts a permutation in place;\n";
  cout << "\n";

  perm_print ( N, p, "  The original permutation:" );
 
  perm_inverse ( N, p );
 
  perm_print ( N, p, "  The inverted permutation:" );
 
  return;
# undef N
}
//****************************************************************************80

void test095 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST095 tests PERM_INVERSE2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int p[N] = { 4, 3, 5, 1, 7, 6, 2 };

  cout << "\n";
  cout << "TEST095\n";
  cout << "  PERM_INVERSE2 inverts a permutation in place.\n";

  perm_print ( N, p, "  The original permutation:" );
 
  perm_inverse2 ( N, p );
 
  perm_print ( N, p, "  The inverted permutation:" );
 
  return;
# undef N
}
//****************************************************************************80

void test0955 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0955 tests PERM_INVERSE3_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int perm[N] = { 3, 2, 4, 0, 6, 5, 1 };
  int *perm_inv;

  cout << "\n";
  cout << "TEST0955\n";
  cout << "  PERM_INVERSE3_NEW inverts a permutation.\n";

  perm_print ( N, perm, "  The original permutation:" );
 
  perm_inv = perm_inverse3_new ( N, perm );
 
  perm_print ( N, perm_inv, "  The inverted permutation:" );
 
  delete [] perm_inv;

  return;
# undef N
}
//****************************************************************************80

void test096 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST096 tests PERM_LEX_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  bool more;
  int p[N];

  cout << "\n";
  cout << "TEST096\n";
  cout << "  PERM_LEX_NEXT generates permutations in order.\n";
  cout << "\n";

  more = false;
 
  for ( ; ; )
  {
    perm_lex_next ( N, p, more );

    if ( !more )
    {
      break;
    }

    perm_print ( N, p, " " );

  }
 
  return;
# undef N
}
//****************************************************************************80

void test097 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST097 tests PERM_LEX_NEXT and PERM_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  bool more;
  int p[N];
  int rank;
  int p_sign;

  cout << "\n";
  cout << "TEST097\n";
  cout << "  PERM_LEX_NEXT generates permutations in order.\n";
  cout << "  PERM_SIGN computes the sign of a permutation.\n";
  cout << "\n";
  cout << "  RANK  SIGN  Permutation\n";
  cout << "\n";

  more = false;
  rank = 0; 

  for ( ; ; )
  {
    perm_lex_next ( N, p, more );

    p_sign = perm_sign ( N, p );
    if ( !more )
    {
      break;
    }

    cout << setw(4) << rank   << "  "
         << setw(4) << p_sign << "  ";

    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << p[i] << "  ";
    }
    cout << "\n";

    rank = rank + 1;

  }
 
  return;
# undef N
}
//****************************************************************************80

void test098 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST098 tests PERM_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int p1[N];
  int p2[N];
  int p3[N];
  int seed;

  cout << "\n";
  cout << "TEST098\n";
  cout << "  PERM_MUL multiplies two permutations.\n";
  cout << "\n";

  seed = 123456789;

  perm_random ( N, seed, p1 );
  perm_random ( N, seed, p2 );

  perm_print ( N, p1, "  Permutation P1:" );

  perm_print ( N, p2, "  Permutation P2:" );

  perm_mul ( N, p1, p2, p3 );

  perm_print ( N, p3, "  Product permutation: P3" );

  return;
# undef N
}
//****************************************************************************80

void test099 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST099 tests PERM_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  bool even;
  bool more;
  int p[N];

  cout << "\n";
  cout << "TEST099\n";
  cout << "  PERM_NEXT generates permutations.\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    perm_next ( N, p, more, even );

    perm_print ( N, p, " " );

    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test100 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST100 tests PERM_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  bool done;
  int iactiv[N];
  int idir[N];
  int invers[N];
  int p[N];

  cout << "\n";
  cout << "TEST100\n";
  cout << "  PERM_NEXT2 generates permutations in order.\n";
  cout << "\n";

  done = true;
 
  for ( ; ; )
  {
    perm_next2 ( N, p, done );
 
    if ( done )
    {
      break;
    }

    perm_print ( N, p, " " );

  }
 
  return;
# undef N
}
//****************************************************************************80

void test101 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST101 tests PERM_NEXT3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  bool more;
  int p[N];

  cout << "\n";
  cout << "TEST101\n";
  cout << "  PERM_NEXT3 generates permutations in order.\n";
  cout << "\n";

  more = false;
 
  for ( ; ; )
  {
    perm_next3 ( N, p, more );

    perm_print ( N, p, " " );

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef N
}
//****************************************************************************80

void test102 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST102 tests PERM_RANDOM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  int p[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST102\n";
  cout << "  PERM_RANDOM produces a random permutation;\n";
  cout << "  For this test, N = " << N << "\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    perm_random ( N, seed, p );
    perm_print ( N, p, " " );
  }
 
  return;
# undef N
}
//****************************************************************************80

void test103 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST103 tests PERM_RANDOM2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  int p[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST103\n";
  cout << "  PERM_RANDOM2 produces a random permutation;\n";
  cout << "  For this test, N = " << N << "\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    p[0] = 101;
    p[1] = 202;
    p[2] = 303;
    p[3] = 404;

    perm_random2 ( N, seed, p );
    perm_print ( N, p, " " );
  }
 
  return;
# undef N
}
//****************************************************************************80

void test104 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST104 tests PERM_RANDOM3;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  int p[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST104\n";
  cout << "  PERM_RANDOM3 produces a random permutation;\n";
  cout << "  For this test, N = " << N << "\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    perm_random3 ( N, seed, p );
    perm_print ( N, p, " " );
  }
 
  return;
# undef N
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests PERM_RANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int invers[N];
  int p[N] = { 1, 4, 2, 3 };
  int rank;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  PERM_RANK ranks a permutation.\n";

  perm_print ( N, p, "  The permutation:" );
 
  rank = perm_rank ( N, p, invers );
 
  cout << "\n";
  cout << "  The rank is " << rank << "\n";
 
  return;
# undef N
}
//****************************************************************************80

void test106 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST106 tests PERM_TO_EQUIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  int a[N];
  int i;
  int jarray[N];
  int npart;
  int p[N] = { 2,3,9,6,7,8,5,4,1 };

  cout << "\n";
  cout << "TEST106\n";
  cout << "  PERM_TO_EQUIV returns the set partition\n";
  cout << "  or equivalence classes determined by a\n";
  cout << "  permutation.\n";

  perm_print ( N, p, "  The input permutation:" );
 
  perm_to_equiv ( N, p, npart, jarray, a );

  equiv_print ( N, a, "  The partition:" );
 
  return;
# undef N
}
//****************************************************************************80

void test107 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST107 tests PERM_TO_YTB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N];
  int lambda[N];
  int p[N] = { 7, 2, 4, 1, 5, 3, 6 };

  cout << "\n";
  cout << "TEST107\n";
  cout << "  PERM_TO_YTB converts a permutation to a\n";
  cout << "  Young tableau.\n";

  perm_print ( N, p, "  The permutation:" );
 
  perm_to_ytb ( N, p, lambda, a );

  ytb_print ( N, a, "  The Young tableau:" );
 
  return;
# undef N
}
//****************************************************************************80

void test108 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST108 tests PERM_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int p[N];
  int rank = 6;

  cout << "\n";
  cout << "TEST108\n";
  cout << "  PERM_UNRANK, given a rank, computes the\n";
  cout << "  corresponding permutation.\n";
  cout << "\n";
  cout << "  The requested rank is " << rank << "\n";
 
  perm_unrank ( N, rank, p );
 
  perm_print ( N, p, "  The permutation:" );
 
  return;
# undef N
}
//****************************************************************************80

void test1085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1085 tests PERRIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int p[N];

  cout << "\n";
  cout << "TEST1085\n";
  cout << "  PERRIN computes the Perrin numbers.\n";
  cout << "\n";

  perrin ( N, p );

  cout << "\n";
  cout << "   N    P(N)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                    << "  "
         << setw(4) << i    << "  "
         << setw(6) << p[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test109 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST109 tests POWER_MOD;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int m;
  int n;
  int x;

  cout << "\n";
  cout << "TEST109\n";
  cout << "  POWER_MOD computes the remainder of a power\n";
  cout << "  of an integer modulo another integer.\n";

  a = 7;
  n = 50;
  m = 11;

  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  N = " << n << "\n";
  cout << "  M = " << m << "\n";
  cout << "  mod ( A**N, M ) = " << power_mod ( a, n, m ) << "\n";

  a = 3;
  n = 118;
  m = 119;

  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  N = " << n << "\n";
  cout << "  M = " << m << "\n";
  cout << "  mod ( A**N, M ) = " << power_mod ( a, n, m ) << "\n";

  return;
}
//****************************************************************************80

void test110 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST110 tests POWER_SERIES1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N];
  double alpha;
  double b[N];
  int i;
  int n;

  cout << "\n";
  cout << "TEST110\n";
  cout << "  POWER_SERIES1 composes a power series;\n";

  alpha = 7.0;
 
  a[0] = 1.0;
  for ( i = 1; i < N; i++ )
  {
    a[i] = 0.0;
  }
 
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
  }

  cout << "\n";
  cout << "  Power series of G(x) = (1+F(x))**alpha\n";
  cout << "\n";
  cout << "  N = "     << N     << "\n";
  cout << "  ALPHA = " << alpha << "\n";

  r8vec_print ( N, a, "  Series for F(x):" );

  power_series1 ( N, alpha, a, b );
 
  r8vec_print ( N, b, "  Series for G(X):" );

  return;
# undef N
}
//****************************************************************************80

void test111 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST111 tests POWER_SERIES2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N];
  double b[N];
  int i;

  cout << "\n";
  cout << "TEST111\n";
  cout << "  POWER_SERIES2 composes a power series;\n";

  a[0] = -4.0;
  for ( i = 1; i < N; i++ )
  {
    a[i] = 0.0;
  } 

  cout << "\n";
  cout << "  Power series of G(x) = exp(F(x))-1\n";
  cout << "\n";
  cout << "  N = " << N << "\n";

  r8vec_print ( N, a, "  Series for F(X):" );

  power_series2 ( N, a, b );
 
  r8vec_print ( N, b, "  Series for G(X):" );

  return;
# undef N
}
//****************************************************************************80

void test112 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST112 tests POWER_SERIES3;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N];
  double b[N];
  double c[N];

  cout << "\n";
  cout << "TEST112\n";
  cout << "  POWER_SERIES3 composes two power series;\n";
 
  a[0] = 1.0;
  a[1] = 1.0;
  a[2] = 0.0;
  a[3] = 0.0;

  r8vec_print ( N, a, "  Series for F(X):" );

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 0.0;
  b[3] = 0.0;

  r8vec_print ( N, b, "  Series for G(X):" );

  power_series3 ( N, a, b, c );
 
  r8vec_print ( N, c, "  Series for H(X) = G(F(X)):" );

  return;
# undef N
}
//****************************************************************************80

void test113 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST113 tests POWER_SERIES4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
//
//  The order of arguments for POWER_SERIES4 is a shame.
//
  double a[N];
  double b[N];
  double c[N];
  int i;

  cout << "\n";
  cout << "TEST113\n";
  cout << "  POWER_SERIES4 composes a power series;\n";
  cout << "  Given power series for F(X) and G(X), we compute\n";
  cout << "  the power series of H(x) = G(1/F(x)).\n";

  for ( i = 0; i < N; i++ )
  {
    a[i] = 1.0 / ( double ) ( i + 1 );
  }
  r8vec_print ( N, a, "  Series for F(x):" );
 
  b[0] = 1.0;
  for ( i = 1; i < N; i++ )
  {
    b[i] = 0.0;
  }
  r8vec_print ( N, b, "  Series for G(x):" );

  power_series4 ( N, a, b, c );
 
  r8vec_print ( N, c, "  Series for H(x):" );

  return;
# undef N
}
//****************************************************************************80

void test114 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST114 tests PYTHAG_TRIPLE_NEXT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  int d;
  int e;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST114\n";
  cout << "  PYTHAG_TRIPLE_NEXT computes the next\n";
  cout << "    Pythagorean triple.\n";
  cout << "\n";
  cout << "     I       J       A       B       C A^2+B^2     C^2\n";
  cout << "\n";

  i = 0;
  j = 0;

  for ( k = 0; k <= 20; k++ )
  {
    pythag_triple_next ( i, j, a, b, c );

    d = a * a + b * b;
    e = c * c;

    cout << setw(6) << i << "  "
         << setw(6) << j << "  "
         << setw(6) << a << "  "
         << setw(6) << b << "  "
         << setw(6) << c << "  "
         << setw(6) << d << "  "
         << setw(6) << e << "\n";
  }

  return;
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests R8_AGM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int seed;
  double x;
  double y;
  double z;

  cout << "\n";
  cout << "TEST115\n";
  cout << "  R8_AGM computes the arithmetic-geometric mean (AGM)\n";
  cout << "  of two nonnegative real numbers.\n";

  cout << "\n";
  cout << "    X        Y    R8_AGM(X,Y)\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = ( double ) i4_uniform ( 1, 10, seed );

    y = ( double ) i4_uniform ( 1, 10, seed );

    z = r8_agm ( x, y );

    cout << setw(10) << x << "  "
         << setw(10) << y << "  "
         << setw(10) << z << "\n";
  }

  return;
}
//****************************************************************************

void test1155 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST1155 tests R8_EPSILON
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  double r;
  double s;

  cout << "\n";
  cout << "TEST1155\n";
  cout << "  R8_EPSILON produces the floating point machine precision.\n";
  cout << "\n";

  r = r8_epsilon ( );
  cout << "  R = R8_EPSILON() = " << setw(10) << r << "\n";

  s = ( double ) ( 1.0 + r ) - 1.0;
  cout << "  ( 1 + R ) - 1 = " << setw(10) << s << "\n";

  s = ( double ) ( 1.0 + ( double ) ( r / 2.0) ) - 1.0;
  cout << "  ( 1 + (R/2) ) - 1 = " << setw(10) << s << "\n";

  return;
}
//****************************************************************************80

void test116 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST116 tests R8_TO_CFRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N+1];
  double err;
  int i;
  int p[N+2];
  int q[N+2];
  double r;
  double temp;

  cout << "\n";
  cout << "TEST116\n";
  cout << "  R8_TO_CFRAC converts a real number to a sequence\n";
  cout << "  of continued fraction convergents.\n";

  r = 2.0 * r8_pi ( );

  cout << "\n";
  cout << "  Use the real number R = " << r << "\n";

  r8_to_cfrac ( r, N, a, p, q );

  cout << "\n";
  cout << " I        A[I]      P[I+1]      Q[I+1]   P[I]/Q[I]    Error\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    temp = ( double ) p[i+1] / ( double ) q[i+1];
    err = r - temp;
    cout                       << "  "
         << setw(2)  << i      << "  "
         << setw(10) << a[i]   << "  "
         << setw(10) << p[i+1] << "  "
         << setw(10) << q[i+1] << "  "
         << setw(14) << temp   << "  "
         << setw(14) << err    << "\n";
  }

  return;
# undef N
}
//****************************************************************************

void test1163 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST1163 tests R8_TO_DEC and DEC_TO_R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int dec_digit;
  int i;
  double r;
  double r2;
  int seed;

  cout << "\n";
  cout << "TEST1163\n";
  cout << "  R8_TO_DEC converts a real number to a decimal;\n";
  cout << "  DEC_TO_R8 converts a decimal to a real number.\n";

  dec_digit = 5;

  cout << "\n";
  cout << "  The maximum number of digits allowed is " << dec_digit << "\n";

  seed = 123456789;

  cout << "\n";
  cout << "     R   =>  A * 10^B  =>  R2\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_uniform_01 ( seed );
    r = 10.0 * ( r - 0.25 );

    r8_to_dec ( r, dec_digit, a, b );
    r2 = dec_to_r8 ( a, b );

    cout                   << "  "
         << setw(10) << r  << "  "
         << setw(6)  << a  << "  "
         << setw(6)  << b  << "  "
         << setw(10) << r2 << "\n";
  }

  return;
}
//****************************************************************************

void test1165 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST1165 tests R8_TO_RAT and RAT_TO_R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int i;
  int ndig = 4;
  double r;
  double r2;
  int seed;

  cout << "\n";
  cout << "TEST1165\n";
  cout << "  R8_TO_RAT converts a real number to a rational;\n";
  cout << "  RAT_TO_R8 converts a rational to a real number.\n";
  cout << "\n";
  cout << "  The maximum number of digits allowed is " << ndig << "\n";

  seed = 123456789;

  cout << "\n";
  cout << "     R   =>  A / B  =>  R2\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_uniform_01 ( seed );
    r = 10.0 * ( r - 0.25 );

    r8_to_rat ( r, ndig, a, b );
    r2 = rat_to_r8 ( a, b );

    cout                   << "  "
         << setw(10) << r  << "  "
         << setw(6)  << a  << "  "
         << setw(6)  << b  << "  "
         << setw(10) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void test117 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST117 tests RAT_ADD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int abot = 4;
  int atop = 3;
  int bbot = 7;
  int btop = 10;
  int cbot;
  int ctop;
  bool error;

  cout << "\n";
  cout << "TEST117\n";
  cout << "  RAT_ADD adds two rationals.\n";

  rat_add ( atop, abot, btop, bbot, ctop, cbot, error );

  cout << "\n";
  cout << "  A = "         << atop << "/" << abot << "\n";
  cout << "  B = "         << btop << "/" << bbot << "\n";
  cout << "  C = A + B = " << ctop << "/" << cbot << "\n";
 
  return;
}
//****************************************************************************80

void test118 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST118 tests RAT_DIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int abot = 4;
  int atop = 3;
  int bbot = 7;
  int btop = 10;
  int cbot;
  int ctop;
  bool error;

  cout << "\n";
  cout << "TEST118\n";
  cout << "  RAT_DIV divides two rationals.\n";

  rat_div ( atop, abot, btop, bbot, ctop, cbot, error );

  cout << "\n";
  cout << "  A = "         << atop << "/" << abot << "\n"; 
  cout << "  B = "         << btop << "/" << bbot << "\n"; 
  cout << "  C = A / B = " << ctop << "/" << cbot << "\n"; 
 
  return;
}
//****************************************************************************80

void test119 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST119 tests RAT_FAREY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define FRAC_MAX 20

  int a[FRAC_MAX];
  int b[FRAC_MAX];
  int frac_num;
  int i;
  int ihi;
  int ilo;
  int n;

  cout << "\n";
  cout << "TEST119\n";
  cout << "  RAT_FAREY computes a row of the Farey fraction table.\n";

  for ( n = 1; n <= 7; n++ )
  {
    rat_farey ( n, FRAC_MAX, frac_num, a, b );
 
    cout << "\n";
    cout << "  Row " << n << "\n";
    cout << "  Number of fractions: " << frac_num << "\n";
    cout << "\n";

    for ( ilo = 0; ilo < frac_num; ilo = ilo + 20 )
    {
      ihi = i4_min ( ilo+20-1, frac_num-1 );
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(3) << a[i];
      }
      cout << "\n";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(3) << b[i];
      }
      cout << "\n";
    }

  }

  return;
# undef FRAC_MAX
}
//****************************************************************************80

void test120 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST120 tests RAT_FAREY2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 4
# define TWO_POWER_MAX 16

  int a[TWO_POWER_MAX+1];
  int b[TWO_POWER_MAX+1];
  int i;
  int ihi;
  int ilo;
  int n;
  int num_frac;
  int two_power;

  cout << "\n";
  cout << "TEST120\n";
  cout << "  RAT_FAREY2 computes a row of the Farey fraction table.\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    rat_farey2 ( n, a, b );
 
    cout << "\n";
    cout << "  Row " << n+1 << "\n";

    two_power = ( int ) pow ( ( double ) 2, n );

    for ( ilo = 0; ilo <= two_power; ilo = ilo + 20 )

    {
      ihi = i4_min ( ilo+20-1, two_power );
      cout << "\n";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(3) << a[i];
      }
      cout << "\n";
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(3) << b[i];
      }
      cout << "\n";

    }

  }

  return;
# undef N_MAX
# undef TWO_POWER_MAX
}
//****************************************************************************80

void test121 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST121 tests RAT_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int abot = 4;
  int atop = 3;
  int bbot = 7;
  int btop = 10;
  int cbot;
  int ctop;
  bool error;

  cout << "\n";
  cout << "TEST121\n";
  cout << "  RAT_MUL multiplies two rationals.\n";

  rat_mul ( atop, abot, btop, bbot, ctop, cbot, error );

  cout << "\n";
  cout << "  A = "         << atop << "/" << abot << "\n";
  cout << "  B = "         << btop << "/" << bbot << "\n";
  cout << "  C = A * B = " << ctop << "/" << cbot << "\n";
 
  return;
}
//****************************************************************************80

void test1215 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1215 tests RAT_WIDTH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 17

  int a;
  int a_test[N_TEST] = {
    1000, 1000, 1000, 1000, 1000, 1, -1, -10, -100, -1000,
    1, 10, 100, 1000, 10000, 17, 4000000 };
  int b;
  int b_test[N_TEST] = {
    3, 40, 500, 6000, 70000, 1, 200, 200, 200, 200, 
   -200, -200, -200, -200, -200, 3000, 4000000 };
  int i;
  int width;

  cout << "\n";
  cout << "TEST1215\n";
  cout << "  RAT_WIDTH determines the \"width\" of a rational.\n";
  cout << "\n";
  cout << "  Top    Bottom  Width\n";
  cout << "\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    a = a_test[i];
    b = b_test[i];

    width = rat_width ( a, b );

    cout                     << "  "
         << setw(8) << a     << "  "
         << setw(8) << b     << "  "
         << setw(8) << width << "\n";
  }

  return;
# undef N_TEST
}
//****************************************************************************80

void test122 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST122 tests RATMAT_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int a[(N+1)*(N+1)];
  int b[(N+1)*(N+1)];

  cout << "\n";
  cout << "TEST122\n";
  cout << "  RAT_SUM_FORMULA computes the coefficients for the\n";
  cout << "  formulas for the sums of powers of integers.\n";
  
  rat_sum_formula ( N, a, b );

  ratmat_print ( N+1, N+1, a, b, "  Power Sum Coefficients:" );

  return;
# undef N
}
//****************************************************************************80

void test123 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST123 tests RATMAT_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N3 3

  int a3[N3*N3];
  int b3[N3*N3];
  bool error;
  int i;
  int idbot;
  int idtop;
  int j;
  int k;

  cout << "\n";
  cout << "TEST123\n";
  cout << "  RATMAT_DET: determinant of a rational matrix.\n";
  cout << "\n";
 
  k = 0;
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      k = k + 1;
      a3[i+j*N3] = k;
    }
  }

  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      b3[i+j*N3] = 1;
    }
  }
 
  ratmat_print ( N3, N3, a3, b3, "  The 123/456/789 matrix:" );

  ratmat_det ( N3, a3, b3, idtop, idbot, error );
 
  cout << "\n";
  cout << "  Determinant of the 123/456/789 matrix = "
       << idtop << "/" << idbot << "\n";
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      a3[i+j*N3] = 1;
    }
  }
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      b3[i+j*N3] = i + j + 2;
    }
  }
  ratmat_print ( N3, N3, a3, b3, "  The Hilbert matrix:" );

  ratmat_det ( N3, a3, b3, idtop, idbot, error );
 
  cout << "\n";
  cout << "  Determinant of the Hilbert matrix = "
       << idtop << "/" << idbot << "\n";
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      if ( i == j )
      {
        a3[i+j*N3] = 2;
      }
      else if ( i == j+1 || i == j-1 )
      {
        a3[i+j*N3] = -1;
      }
      else
      {
        a3[i+j*N3] = 0;
      }
    }
  }
 
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      b3[i+j*N3] = 1;
    }
  }

  ratmat_print ( N3, N3, a3, b3, "  The -1 2 -1 matrix:" );

  ratmat_det ( N3, a3, b3, idtop, idbot, error );
 
  cout << "\n";
  cout << "  Determinant of the -1,2,-1 matrix = "
       << idtop << "/" << idbot << "\n";
 
  return;
# undef N3
}
//****************************************************************************80

void test124 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST124 tests REGRO_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  bool done;
  int i;
  int rank;
  int v[N];
  int vmax[N];

  cout << "\n";
  cout << "TEST124\n";
  cout << "  REGRO_NEXT generates all restricted growth\n";
  cout << "  functions.\n";
  cout << "\n";

  rank = 0;

  done = true;
 
  for ( ; ; )
  {
    regro_next ( done, N, v, vmax );

    if ( done )
    {
      break;
    }

    rank = rank + 1;
    cout                    << "  "
         << setw(3) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(1) << v[i] << "  ";
    }
    cout << "\n";

  }
 
  return;
# undef N
}
//****************************************************************************

void test1245 ( )

//****************************************************************************
//
//  Purpose:
//
//    TEST1245 tests R8_RISE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double x;

  cout << "\n";
  cout << "TEST1245\n";
  cout << "  R8_RISE computes the rising factorial function.\n";
  cout << "  [X]_N = X * (X+1) * (X+2) * ... * ( X+N-1).\n";
  cout << "\n";
  cout << "    X          N  R8_RISE(X,N)\n";
  cout << "\n";

  x = 4.0;

  for ( n = -2; n <= 5; n++ )
  {
    cout << setw(10) << x             << "  "
         << setw(8)  << n             << "  "
         << setw(10) << r8_rise ( x, n ) << "\n";
  }

  return;
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests R8MAT_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N3 3
# define N4 4

  double a3[N3*N3];
  double a4[N4*N4];
  double det;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  R8MAT_DET: determinant of a real matrix.\n";
  cout << "\n";

  k = 0;
  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      k = k + 1;
      a3[i+j*N3] = ( double ) k;
    }
  }

  r8mat_print ( N3, N3, a3, "  The 123/456/789 matrix:" );

  det = r8mat_det ( N3, a3 );

  cout << "\n";
  cout << "  Determinant of the 123/456/789 matrix is " << det << "\n";

  for ( i = 0; i < N4; i++ )
  {
    for ( j = 0; j < N4; j++ )
    {
      a4[i+j*N4] = 1.0 / ( double ) ( 2 + i + j );
    }
  }

  r8mat_print ( N4, N4, a4, "  The Hilbert matrix:" );

  det = r8mat_det ( N4, a4 );

  cout << "\n";
  cout << "  Determinant of the Hilbert matrix is " << det << "\n";

  for ( i = 0; i < N3; i++ )
  {
    for ( j = 0; j < N3; j++ )
    {
      if ( i == j )
      {
        a3[i+j*N3] = 2.0;
      }
      else if ( i == j+1 || i == j-1 )
      {
        a3[i+j*N3] = -1.0;
      }
      else
      {
        a3[i+j*N3] = 0.0;
      }
    }
  }

  r8mat_print ( N3, N3, a3, "  The -1,2,-1 matrix:" );

  det = r8mat_det ( N3, a3 );

  cout << "\n";
  cout << "  Determinant of the -1,2,-1 matrix is " << det << "\n";

  return;
# undef N3
# undef N4
}
//****************************************************************************80

void test126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST126 tests R8MAT_PERM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  double a[N*N];
  int i;
  int p[N] = { 2,3,9,6,7,8,5,4,1 };
  int j;

  cout << "\n";
  cout << "TEST126\n";
  cout << "  R8MAT_PERM reorders a real matrix in place.\n";
  cout << "  The rows and columns use the same permutation.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( (i+1) * 10 + (j+1) );
    }
  }

  r8mat_print ( N, N, a, "  The original matrix" );

  perm_print ( N, p, "  The row and column permutation:" );

  r8mat_perm ( N, a, p );

  r8mat_print ( N, N, a, "  The permuted matrix" );

  return;
# undef N
}
//****************************************************************************80

void test127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST127 tests R8MAT_PERM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 9
# define N 7

  double a[M*N];
  int i;
  int j;
  int p[M] = { 2, 3, 9, 6, 7, 8, 5, 4, 1 };
  int q[N] = { 3, 4, 5, 6, 7, 1, 2 };

  cout << "\n";
  cout << "TEST127\n";
  cout << "  R8MAT_PERM2 reorders a real matrix in place.\n";
  cout << "  Rows and columns use different permutations.\n";
 
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
 
  r8mat_print ( M, N, a, "  The original matrix" );
 
  perm_print ( M, p, "  The row permutation:" );
 
  perm_print ( N, q, "  The column permutation:" );

  r8mat_perm2 ( M, N, a, p, q );
 
  r8mat_print ( M, N, a, "  The permuted matrix" );
 
  return;
# undef M
# undef N
}
//****************************************************************************80

void test128 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST128 tests R8MAT_PERMANENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int n;
  double perm;

  cout << "\n";
  cout << "TEST128\n";
  cout << "  R8MAT_PERMANENT: the matrix permanent function.\n";
  cout << "  We will analyze matrices with 0 diagonal and\n";
  cout << "  1 on all offdiagonals.\n";
  cout << "\n";
  cout << "  Order	    Permanent.\n";
  cout << "\n";
  cout << "DEBUG\n" << flush;
  for ( n = 2; n <= 12; n++ )
  {
    a = new double[n*n];

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( i == j ) 
        {
          a[i+j*n] = 0.0;
        }
        else
        {
          a[i+j*n] = 1.0;
        }
      }
    }

    perm = r8mat_permanent ( n, a );

    cout << setw(4)  << n    << "  "
         << setw(10) << perm << "\n";
    cout << "DEBUG, N = " << n << "\n" << flush;
    delete [] a;

  }

  return;
}
//****************************************************************************80

void test129 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST129 test R8POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int i;
  double a[N];
  int iopt;
  int test;
  double val;
  double x0;

  cout << "\n";
  cout << "TEST129\n";
  cout << "  R8POLY converts between power sum, factorial\n";
  cout << "  and Taylor forms, and can evaluate a polynomial\n";
  cout << "\n";
 
  for ( test = 1; test <= 6; test++ )
  {
    if ( test == 1 )
    {
      iopt = -3;
    }
    else if ( test == 2 )
    {
      iopt = -2;
    }
    else if ( test == 3 )
    {
      iopt = -1;
      x0 = 2.0;
    }
    else if ( test == 4 )
    {
      iopt = 0;
      x0 = 2.0;
    }
    else if ( test == 5 )
    {
      iopt = 6;
      x0 = 2.0;
    }
    else if ( test == 6 )
    {
      iopt = 6;
      x0 = -2.0;
    }

    for ( i = 0; i < N-1; i++ )
    {
      a[i] = 0.0;
    }
    a[N-1] = 1.0;

    if ( test == 1 )
    {
      r8vec_print ( N, a, "  All calls  have input A as follows:" );
    }
 
    r8poly ( N, a, x0, iopt, val );
 
    cout << "\n";
    cout << "  Option IOPT = " << iopt << "\n";

    if ( -1 <= iopt )
    {
      cout << "  X0 = " << x0 << "\n";
    }

    if ( iopt == -3 || iopt == -2 || 0 < iopt )
    {
      r8vec_print ( N, a, "  Output array:" );
    }

    if ( iopt == -1 || iopt == 0 )
    {
      cout << "  Value = " << val << "\n";
    }
 
  }

  return;
# undef N
}
//****************************************************************************80

void test130 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST130 tests R8POLY_DIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a[11];
  double b[11];
  int i;
  int na;
  int nb;
  int nq;
  int nr;
  int ntest = 2;
  double q[11];
  double r[11];
  int test;

  cout << "\n";
  cout << "TEST130\n";
  cout << "  R8POLY_DIV computes the quotient and\n";
  cout << "  remainder for polynomial division.\n";
  cout << "\n";
//
//  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
//     Quotient is 3+4*X+X**2, remainder is 0.
//
//  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
//     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
//
  for ( test = 1; test <= ntest; test++ )
  {
    if ( test == 1 )
    {
      na = 3;
      a[0] = -6.0;
      a[1] = -5.0;
      a[2] =  2.0;
      a[3] =  1.0;
      nb = 1;
      b[0] = -2.0;
      b[1] =  1.0;
    }
    else if ( test == 2 )
    {
      na = 4;
      a[0] = -2.0;
      a[1] =  5.0;
      a[2] =  2.0;
      a[3] =  3.0;
      a[4] =  1.0;
      nb = 2;
      b[0] = -3.0;
      b[1] =  1.0;
      b[2] =  1.0;
    }

    r8poly_print ( na, a, "  The polynomial to be divided, A:" );

    r8poly_print ( nb, b, "  The divisor polynomial, B:" );

    r8poly_div ( na, a, nb, b, nq, q, nr, r );

    r8poly_print ( nq, q, "  The quotient polynomial, Q:" );

    r8poly_print ( nr, r, "  The remainder polynomial, R:" );

  }
  return;
}
//****************************************************************************80

void test131 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST131 tests R8POLY_F2P and R8POLY_P2F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N];
  int i;

  cout << "\n";
  cout << "TEST131\n";
  cout << "  R8POLY_P2F: power sum => factorial;\n";
  cout << "  R8POLY_F2P: factorial => power sum.\n";

  r8vec_indicator ( N, a );

  r8poly_print ( N-1, a, "  The power sum polynomial:" );

  r8poly_p2f ( N, a );
 
  r8vec_print ( N, a, "  The factorial coefficients:" );

  r8poly_f2p ( N, a );
 
  r8poly_print ( N-1, a, "  The recovered power sum polynomial:" );

  return;
# undef N
}
//****************************************************************************80

void test132 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST132 tests R8POLY_FVAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N];
  int i;
  double val;
  double x;

  cout << "\n";
  cout << "TEST132\n";
  cout << "  R8POLY_FVAL evaluates a polynomial in factorial form.\n";

  r8vec_indicator ( N, a );
 
  r8vec_print ( N, a, "  The factorial coefficients:" );

  x = 2.0;

  val = r8poly_fval ( N, a, x );

  cout << "\n";
  cout << "  RPOLY (" << x << ") = " << val << "\n";
  cout << "  The correct value is 11.\n";
 
  return;
# undef N
}
//****************************************************************************80

void test133 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST133 tests R8POLY_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXN 5

  double a[MAXN+1];
  double b[MAXN+1];
  double c[MAXN+1];
  int na;
  int nb;
  int ntest = 2;
  int test;

  cout << "\n";
  cout << "TEST133\n";
  cout << "  R8POLY_MUL multiplies two polynomials.\n";
  cout << "\n";
//
//  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
//  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
//
  for ( test = 1; test <= ntest; test++ )
  {
    if ( test == 1 )
    {
      na = 1;
      a[0] = 1.0;
      a[1] = 1.0;
      nb = 1;
      b[0] =  1.0;
      b[1] = -1.0;
    }
    else if ( test == 2 )
    {
      na = 2;
      a[0] = 1.0;
      a[1] = 2.0;
      a[2] = 3.0;
      nb = 1;
      b[0] =  1.0;
      b[1] = -2.0;
    }

    r8poly_mul ( na, a, nb, b, c );

    r8poly_print ( na, a, "  The factor A:" );

    r8poly_print ( nb, b, "  The factor B:" );

    r8poly_print ( na+nb, c, "  The product C = A*B:" );

  }

  return;
# undef MAXN
}
//****************************************************************************80

void test134 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST134 tests R8POLY_N2P and R8POLY_P2N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  double a[N];
  double a2[N];
  double a3[N];

  cout << "\n";
  cout << "TEST134\n";
  cout << "  R8POLY_N2P: Newton => power sum;\n";
  cout << "  R8POLY_P2N: Power sum => Newton.\n";

  r8vec_indicator ( N, a );

  for ( i = 0; i < N; i++ )
  {
    a2[i] = 2.0 * a[i];
  }

  r8poly_print ( N-1, a, "  The power sum polynomial:" );

  r8poly_p2n ( N, a, a2 );
 
  r8vec_print ( N, a, "  Derived Newton form coefficients:" );

  r8vec_print ( N, a2, "  Newton form abscissas:" );

  r8poly_n2p ( N, a, a2 );
 
  r8poly_print ( N-1, a, "  The recovered power sum polynomial:" );

  return;
# undef N
}
//****************************************************************************80

void test135 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST135 tests R8POLY_NVAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N];
  double a2[N];
  int i;
  double val;
  double x;

  cout << "\n";
  cout << "TEST135\n";
  cout << "  R8POLY_NVAL evaluates a polynomial in Newton form.\n";

  r8vec_indicator ( N, a );

  for ( i = 0; i < N; i++ )
  {
    a2[i] = a[i] - 1.0;
  }
 
  r8vec_print ( N, a, "  Newton polynomial coefficients:" );

  r8vec_print ( N, a2, "  Newton polynomial abscissas:" );

  x = 2.0;
 
  val = r8poly_nval ( N, a, a2, x );
 

  cout << "\n";
  cout << "  RPOLY (" << x << ") = " << val << "\n";
  cout << "  The correct value is 11.\n";
 
  return;
# undef N
}
//****************************************************************************80

void test136 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST136 tests R8POLY_P2T and R8POLY_T2P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N+1];
  int i;
  double x = 2.0;

  cout << "\n";
  cout << "TEST136\n";
  cout << "  R8POLY_T2P: Taylor => Power sum;\n";
  cout << "  R8POLY_P2T: Power sum => Taylor.\n";
  cout << "  The Taylor form uses the base point X0 = " << x << "\n";

  r8vec_indicator ( N+1, a );

  r8vec_print ( N, a, "  Initial Taylor sum form:" );

  r8poly_t2p ( N, a, x );

  r8poly_print ( N, a, "  Power sum form:" );

  r8poly_p2t ( N, a, x );

  r8vec_print ( N, a, "  Recovered Taylor sum form:" );

  return;
# undef N
}
//****************************************************************************80

void test137 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST137 tests R8POLY_POWER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define LMAX 10

  double a[LMAX+1];
  double b[LMAX+1];
  int na;
  int p;

  cout << "\n";
  cout << "TEST137\n";
  cout << "  R8POLY_POWER takes a polynomial to a power.\n";
//
//  Cube (2-X).  Answer is 8-12*X+6*X**2-X**3.
//
  na = 1;
  a[0] =  2.0;
  a[1] = -1.0;
  p = 3;

  r8poly_print ( na, a, "  The polynomial A:" );
 
  r8poly_power ( na, a, p, b );
 
  r8poly_print ( p*na, b, "  Raised to the power 3:" );
//
//  Square X+X**2
//
  na = 2;
  a[0] =  0.0;
  a[1] =  1.0;
  a[2] =  1.0;
  p = 2;

  r8poly_print ( na, a, "  The polynomial A:" );
 
  r8poly_power ( na, a, p, b );
 
  r8poly_print ( p*na, b, "  Raised to the power 2:" );
 
  return;
# undef LMAX
}
//****************************************************************************80

void test138 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST138 tests R8POLY_PVAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  double a[N+1];
  double val;
  double x;

  r8vec_indicator ( N+1, a );

  cout << "\n";
  cout << "TEST138\n";
  cout << "  R8POLY_PVAL evaluates a polynomial\n";
  cout << "  in power sum form.\n";

  r8poly_print ( N, a, "  The polynomial to be evaluated:" );

  x = 2.0;
 
  val = r8poly_pval ( N, a, x );

  cout << "  At X = " << x << "\n";
  cout << "  Computed polynomial value is " << val << "\n";
  cout << "  Correct value is 129.\n";
 
  return;
# undef N
}
//****************************************************************************80

void test139 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST139 tests R8VEC_FRAC;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N];
  double ahi = 10.0;
  double alo = 0.0;
  int k;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST139\n";
  cout << "  R8VEC_FRAC: K-th smallest real vector entry;\n";

  r8vec_uniform ( N, alo, ahi, seed, a );

  r8vec_print ( N, a, "  The real array to search: " );

  cout << "\n";
  cout << "Frac     Value\n";
  cout << "\n";

  for ( k = 1; k <= N; k++ )
  {
    cout << setw(4)  << k                     << "  "
         << setw(10) << r8vec_frac ( N, a, k ) << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test1395 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1395 tests R8VEC_MIRROR_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N];
  bool done;

  cout << "\n";
  cout << "TEST1395\n";
  cout << "  R8VEC_MIRROR_NEXT generates all sign variations\n";
  cout << "    of a real vector.\n";

  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = 3.0;

  for ( ; ; )
  {
    r8vec_print ( N, a, "  Next vector:" );

    done = r8vec_mirror_next ( N, a );

    if ( done )
    {
      cout << "\n";
      cout << "  Done.\n";
      break;
    }
  }

  a[0] = 1.0;
  a[1] = 0.0;
  a[2] = 3.0;

  for ( ; ; )
  {
    r8vec_print ( N, a, "  Next vector:" );

    done = r8vec_mirror_next ( N, a );

    if ( done )
    {
      cout << "\n";
      cout << "  Done.\n";
      break;
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test140 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST140 tests SCHROEDER;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int s[N];

  cout << "\n";
  cout << "TEST140\n";
  cout << "  SCHROEDER computes the Schroeder numbers.\n";
  cout << "\n";

  schroeder ( N, s );

  cout << "\n";
  cout << "   N    S(N)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                    << "  "
         << setw(4) << i+1  << "  "
         << setw(6) << s[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test141 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST141 tests SORT_HEAP_EXTERNAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int i;
  int indx;
  int isgn;
  int itemp;
  int j;
  int seed;
  int temp;

  cout << "\n";
  cout << "TEST141\n";
  cout << "  SORT_HEAP_EXTERNAL sorts objects externally.\n";
  cout << "\n";

  indx = 0;
  i = 0;
  j = 0;
  isgn = 0;
  seed = 123456789;

  a = i4vec_uniform ( N, 1, N, seed );

  i4vec1_print ( N, a, "  Before sorting:" );
 
  for ( ; ; )
  {
    sort_heap_external ( N, indx, i, j, isgn );
 
    if ( indx < 0 )
    {

      if ( a[i-1] <= a[j-1] )
      {
        isgn = -1;
      }
      else
      {
        isgn = +1;
      }
    }
    else if ( 0 < indx )
    {
      temp = a[i-1];
      a[i-1] = a[j-1];
      a[j-1] = temp;
    }
    else
    {
      break;
    }

  }

  i4vec1_print ( N, a, "  After sorting:" ); 

  delete [] a;
 
  return;
# undef N
}
//****************************************************************************80

void test142 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST142 tests SUBSET_BY_SIZE_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  bool more;
  int rank;
  int size;

  cout << "\n";
  cout << "TEST142\n";
  cout << "  SUBSET_BY_SIZE_NEXT generates all subsets of an N set.\n";
  cout << "\n";

  more = false;
  rank = 0;

  for ( ; ; )
  {
    subset_by_size_next ( N, a, size, more );

    rank = rank + 1;

    cout << setw(4) << rank << "  ";

    if ( 0 < size )
    {
      for ( i = 0; i < size; i++ )
      {
        cout << setw(2) << a[i] << "  ";
      }
      cout << "\n";
    }
    else
    {
      cout << "The empty set\n";
    }

    if ( !more )
    {
      break;
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test143 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST143 tests SUBSET_LEX_NEXT without size restrictions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  int k;
  bool ltest;
  int ndim;
  int rank;

  cout << "\n";
  cout << "TEST143\n";
  cout << "  SUBSET_LEX_NEXT generates all subsets of an N set.\n";

  ndim = N;
  k = 0;
  rank = 0;
  ltest = false;
 
  for ( ; ; )
  {

    subset_lex_next ( N, ltest, ndim, k, a );
    rank = rank + 1;

    cout << setw(4) << rank << "  ";

    if ( k > 0 )
    {
      for ( i = 0; i < k; i++ )
      {
        cout << setw(2) << a[i] << "  ";
      }
      cout << "\n";
    }
    else
    {
      cout << "The empty set.\n";
    }
 
    if ( k == 0 )
    {
      break;
    }

  }
 
  return;
# undef NDIM
}
//****************************************************************************80

void test1435 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1435 tests SUBSET_LEX_NEXT with size restrictions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDIM 3

  int a[NDIM];
  int i;
  int k;
  bool ltest;
  int n = 5;

  cout << "\n";
  cout << "TEST1435\n";
  cout << "  SUBSET_LEX_NEXT generates all subsets of an N set.\n";
  cout << "  The user can impose a restriction on the\n";
  cout << "  maximum size of the subsets.\n";
  cout << "\n";
  cout << "  Here, we require the subsets to be no larger\n";
  cout << "  than NDIM = " << NDIM << "\n";

  k = 0;
 
  for ( ; ; )
  {
    ltest = ( k == NDIM );

    subset_lex_next ( n, ltest, NDIM, k, a );
 
    if ( 0 < k )
    {
      cout << "  ";
      for ( i = 0; i < k; i++ )
      {
        cout << setw(2) << a[i] << "  ";
      }
      cout << "\n";
    }
    else
    {
      cout << "  The empty set.\n";
    }
 
    if ( k == 0 )
    {
      break;
    }

  }
 
  return;
# undef NDIM
}
//****************************************************************************80

void test144 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST144 tests SUBSET_GRAY_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  int j;
  int iadd;
  bool more;
  int ncard;

  cout << "\n";
  cout << "TEST144\n";
  cout << "  SUBSET_GRAY_NEXT generates all subsets of an N set\n";
  cout << "  using the Gray code ordering:\n";
  cout << "  0 0 1 0 1 means the subset contains 3 and 5.\n";
  cout << "\n";
  cout << "  Gray code\n";
  cout << "\n";

  more = false;
  j = 0;

  for ( ; ; )
  {
    subset_gray_next ( N, a, more, ncard, iadd );

    j = j + 1;
    cout << setw(4) << j << "    ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(2) << a[i];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }

  return;

# undef N
}
//****************************************************************************80

void test145 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST145 tests SUBSET_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST145\n";
  cout << "  SUBSET_RANDOM picks a subset at random.\n";
  cout << "  The number of elements in the main set is " << N << "\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 5; j++ )
  {
    subset_random ( N, seed, a );

    cout << setw(4) << j << "    ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(2) << a[i];
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test146 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST146 tests SUBSET_GRAY_RANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N] = { 1, 0, 1, 1, 0 };
  int i;
  int rank;

  cout << "\n";
  cout << "TEST146\n";
  cout << "  SUBSET_GRAY_RANK returns rank of a subset of an N set\n";
  cout << "  using the Gray code ordering.\n";
  cout << "\n";
  cout << "  For N = " << N << ", the subset is:\n";

  cout << "  ";
  for ( i = 0; i < N; i++ )
  {
    cout << a[i] << " ";
  }
  cout << "\n";
 
  rank = subset_gray_rank ( N, a );
 
  cout << "\n";
  cout << "  The rank is " << rank << "\n";
 
  return;
# undef N
}
//****************************************************************************80

void test147 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST146 tests SUBSET_GRAY_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int i;
  int rank;

  cout << "\n";
  cout << "TEST147\n";
  cout << "  SUBSET_GRAY_UNRANK finds the subset of an N set\n";
  cout << "  of a given rank under the Gray code ordering.\n";
  cout << "\n";
  cout << "  N is " << N << "\n";
  cout << "\n";
  cout << "  Rank   Subset\n";
  cout << "\n";

  for ( rank = 1; rank <= 10; rank++ )
  {
    subset_gray_unrank ( rank, N, a );

    cout << "  "
         << setw(4) << rank << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(2) << a[i];
    }
    cout << "\n";

  }
 
  return;
# undef N
}
//****************************************************************************80

void test1475 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1475 tests SUBCOMP_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int count;
  int h;
  int i;
  bool more;
  int n = 6;
  int t;
  int total;

  cout << "\n";
  cout << "TEST1475\n";
  cout << "  SUBCOMP_NEXT generates subcompositions.\n";
  cout << "\n";
  cout << "  Seek all subcompositions of N = " << n << "\n";
  cout << "  using K = " << K << " parts.\n";
  cout << "\n";
  cout << "     #   Sum\n";
  cout << "\n";

  more = false;
  count = 0;

  for ( ; ; )
  {
    subcomp_next ( n, K, a, more, h, t );

    total = 0;
    for ( i = 0; i < K; i++ )
    {
      total = total + a[i];
    }
    count = count + 1;
    cout << "  " << setw(4) << count
         << "  " << setw(4) << total
         << "  ";

    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef K
}
//****************************************************************************80

void test1476 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1476 tests SUBCOMPNZ_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2005
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int count;
  int i;
  bool more;
  int n = 6;
  int total;

  cout << "\n";
  cout << "TEST1476\n";
  cout << "  SUBCOMPNZ_NEXT generates subcompositions using nonzero parts.\n";
  cout << "\n";
  cout << "  Seek all subcompositions of N = " << n << "\n";
  cout << "  using K = " << K << " nonzero parts.\n";
  cout << "\n";
  cout << "     #   Sum\n";
  cout << "\n";

  more = false;
  count = 0;

  for ( ; ; )
  {
    subcompnz_next ( n, K, a, more );

    total = 0;
    for ( i = 0; i < K; i++ )
    {
      total = total + a[i];
    }
    count = count + 1;
    cout << "  " << setw(4) << count
         << "  " << setw(4) << total
         << "  ";

    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef K
}
//****************************************************************************80

void test1477 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1477 tests SUBCOMPNZ2_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2005
//
//  Author:
//
//    John Burkardt
//
{
# define K 3

  int a[K];
  int count;
  int i;
  bool more;
  int n;
  int n_hi = 7;
  int n_lo = 5;

  cout << "\n";
  cout << "TEST1477\n";
  cout << "  SUBCOMPNZ2_NEXT generates subcompositions using nonzero parts.\n";
  cout << "\n";
  cout << "  Seek all subcompositions of N\n";
  cout << "  using K = " << K << " nonzero parts.\n";
  cout << "\n";
  cout << "  Here N is in the range " << n_lo << " <= N <= " << n_hi << "\n";
  cout << "\n";
  cout << "     #     N\n";
  cout << "\n";

  more = false;
  count = 0;

  for ( ; ; )
  {
    subcompnz2_next ( n_lo, n_hi, K, a, more );

    n = 0;
    for ( i = 0; i < K; i++ )
    {
      n = n + a[i];
    }
    count = count + 1;
    cout << "  " << setw(4) << count
         << "  " << setw(4) << n
         << "  ";

    for ( i = 0; i < K; i++ )
    {
      cout << setw(4) << a[i];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef K
}
//****************************************************************************80

void test1478 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1478 tests SUBTRIANGLE_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i1;
  int i2;
  int i3;
  int j1;
  int j2;
  int j3;
  bool more;
  int n;
  int rank;

  n = 4;
  rank = 0;

  more = false;
  i1 = 0;
  j1 = 0;
  i2 = 0;
  j2 = 0;
  i3 = 0;
  j3 = 0;

  cout << "\n";
  cout << "TEST1478\n";
  cout << "  SUBTRIANGLE_NEXT generates the indices of subtriangles\n";
  cout << "  in a triangle whose edges were divided into N subedges.\n";
  cout << "\n";
  cout << "  For this test, N = " << n << "\n";
  cout << "\n";
  cout << "  Rank    I1  J1    I2  J2    I3  J3\n";
  cout << "\n";

  for ( ; ; )
  {
    subtriangle_next ( n, more, i1, j1, i2, j2, i3, j3 );

    rank = rank + 1;

    cout << "  " << setw(4) << rank << "  "
         << "  " << setw(2) << i1
         << "  " << setw(2) << j1 << "  "
         << "  " << setw(2) << i2
         << "  " << setw(2) << j2 << "  "
         << "  " << setw(2) << i3
         << "  " << setw(2) << j3 << "\n"; 

    if ( !more )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void test148 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST148 tests THUE_BINARY_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 100

  int i;
  int j;
  int n;
  int thue[N_MAX];

  cout << "\n";
  cout << "TEST148\n";
  cout << "  THUE_BINARY_NEXT returns the next\n";
  cout << "  Thue binary sequence.\n";
  cout << "\n";

  n = 1;
  thue[0] = 0;
  cout << setw(4) << n << "    ";
  for ( i = 0; i < n; i++ )
  {
    cout << thue[i];
  }
  cout << "\n";

  for ( i = 1; i <= 6; i++ )
  {
    thue_binary_next ( n, thue );

    cout << setw(4) << n << "    ";
    for ( j = 0; j < n; j++ )
    {
      cout << thue[j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test149 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST149 tests THUE_TERNARY_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 100

  int i;
  int j;
  int n;
  int thue[N_MAX];

  cout << "\n";
  cout << "TEST149\n";
  cout << "  THUE_TERNARY_NEXT returns the next\n";
  cout << "  Thue ternary sequence.\n";
  cout << "\n";

  n = 1;
  thue[0] = 1;
  cout << setw(4) << n << "    ";
  for ( i = 0; i < n; i++ )
  {
    cout << thue[i];
  }
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    thue_ternary_next ( n, thue );

    cout << setw(4) << n << "    ";
    for ( j = 0; j < n; j++ )
    {
      cout << thue[j];
    }
    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test150 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST150 tests TUPLE_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int m1 = 2;
  int m2 = 4;
  int rank;
  int x[N];

  cout << "\n";
  cout << "TEST150\n";
  cout << "  TUPLE_NEXT returns the next \"tuple\", that is,\n";
  cout << "  a vector of N integers, each between M1 and M2.\n";
  cout << "\n";
  cout << "  M1 = " << m1 << "\n";
  cout << "  M2 = " << m2 << "\n";
  cout << "  N = " << N << "\n";
  cout << "\n";

  rank = 0;

  for ( ; ; )
  {
    tuple_next ( m1, m2, N, rank, x );

    if ( rank == 0 )
    {
      break;
    }

    cout << setw(4) << rank;
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << x[i] << "  ";
    }
    cout << "\n";

  }

  return;
}
//****************************************************************************80

void test151 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST151 tests TUPLE_NEXT_FAST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int m = 3;
  int rank;
  int rank_hi;
  int x[N];

  cout << "\n";
  cout << "TEST151\n";
  cout << "  TUPLE_NEXT_FAST returns the next \"tuple\", that is,\n";
  cout << "  a vector of N integers, each between 1 and M.\n";
  cout << "\n";
  cout << "  M = " << m << "\n";
  cout << "  N = " << N << "\n";
  cout << "\n";
//
//  Initialize.
//
  rank = -1;
  tuple_next_fast ( m, N, rank, x );

  rank_hi = ( int ) pow ( ( double ) m, N );

  for ( rank = 0; rank < rank_hi; rank++ )
  {
    tuple_next_fast ( m, N, rank, x );

    cout << setw(4) << rank;
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << x[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test152 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST152 tests TUPLE_NEXT_GE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int i;
  int m = 3;
  int rank;
  int x[N];

  cout << "\n";
  cout << "TEST152\n";
  cout << "  TUPLE_NEXT_GE returns the next nondecreasting \"tuple\",\n";
  cout << "  that is, a vector of N integers, each between 1 and M,\n";
  cout << "  with the additional property that the digits never decrease\n";
  cout << "  reading from left to right.\n";
  cout << "\n";
  cout << "  M = " << m << "\n";
  cout << "  N = " << N << "\n";
  cout << "\n";

  rank = 0;

  for ( ; ; )
  {
    tuple_next_ge ( m, N, rank, x );

    if ( rank == 0 )
    {
      break;
    }

    cout << setw(4) << rank;
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << x[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test153 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST153 tests TUPLE_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int i;
  int rank;
  int x[N];
  int xmin[N] = { 2, 3, 8 };
  int xmax[N] = { 4, 3, 5 };

  cout << "\n";
  cout << "TEST153\n";
  cout << "  TUPLE_NEXT_GE returns the next \"tuple\",\n";
  cout << "  that is, a vector of N integers.\n";
  cout << "  Each position in the vector has a separate min and max.\n";
  cout << "  reading from left to right.\n";
  cout << "\n";
  cout << "  N = " << N << "\n";
  cout << "\n";
  i4vec1_print ( N, xmin, "  The minimum values:" );
  i4vec1_print ( N, xmax, "  The maximum values:" );

  cout << "\n";
  cout << "\n";

  rank = 0;

  for ( ; ; )
  {
    tuple_next2 ( N, xmin, xmax, x, rank );

    if ( rank == 0 )
    {
      break;
    }

    cout << setw(4) << rank;
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << x[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test1531 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1531 tests UBVEC_ADD;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int i;
  int j;
  int k;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST1531\n";
  cout << "  UBVEC_ADD adds unsigned binary vectors \n";
  cout << "  representing unsigned integers;\n";
  cout << "\n";
  cout << "        I        J        K = I + J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( 0, 100, seed );
    j = i4_uniform ( 0, 100, seed );

    cout << "\n";
    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";

    k = i + j;

    cout << "  Directly:         "
         << "  " << setw(8) << k << "\n";

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );

    bvec_add ( N, bvec1, bvec2, bvec3 );
    k = bvec_to_i4 ( N, bvec3 );

    cout << "  BVEC_ADD          "
         << "  " << setw(8) << k << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test0626 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0626 tests UI4_TO_UBVEC and UBVEC_TO_UI4;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int bvec[N];
  unsigned int i;
  unsigned int i2;
  int j;

  cout << "\n";
  cout << "TEST0626\n";
  cout << "  UI4_TO_UBVEC converts an unsigned integer to an \n";
  cout << "    unsigned binary vector;\n";
  cout << "  UBVEC_TO_UI4 converts an unsigned binary vector\n";
  cout << "    to an unsigned integer;\n";
  cout << "\n";
  cout << "  I --> BVEC  -->  I\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    ui4_to_ubvec ( i, N, bvec );

    i2 = ubvec_to_ui4 ( N, bvec );

    cout << setw(3) << i << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(1) << bvec[j];
    }
    cout << "  ";
    cout << setw(3) << i2 << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test1535 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1535 tests VEC_COLEX_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int a[DIM_NUM];
  int base = 3;
  int i;
  bool more;

  cout << "\n";
  cout << "TEST1535\n";
  cout << "  VEC_COLEX_NEXT generates all DIM_NUM-vectors\n";
  cout << "  in colex order in a given base BASE.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << DIM_NUM << "\n";
  cout << "  The base BASE =                 " << base << "\n";

  cout << "\n";

  more = false;

  for ( ; ; )
  {
    vec_colex_next ( DIM_NUM, base, a, more );

    if ( !more ) 
    {
      break;
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1536 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1536 tests VEC_COLEX_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int a[DIM_NUM];
  int base[DIM_NUM] = { 2, 1, 3 };
  int i;
  bool more;

  cout << "\n";
  cout << "TEST1536\n";
  cout << "  VEC_COLEX_NEXT2 generates all DIM_NUM-vectors\n";
  cout << "  in colex order in a given base BASE.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << DIM_NUM << "\n";
  cout << "\n";
  cout << "  The base vector:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << setw(4) << base[i] << "  ";
  }
  cout << "\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    vec_colex_next2 ( DIM_NUM, base, a, more );

    if ( !more ) 
    {
      break;
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1537 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1537 tests VEC_COLEX_NEXT3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int a[DIM_NUM];
  int base[DIM_NUM] = { 2, 1, 3 };
  int i;
  bool more;

  cout << "\n";
  cout << "TEST1537\n";
  cout << "  VEC_COLEX_NEXT3 generates all DIM_NUM-vectors\n";
  cout << "  in colex order in a given base BASE.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << DIM_NUM << "\n";
  cout << "\n";
  cout << "  The base vector:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << setw(4) << base[i] << "  ";
  }
  cout << "\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    vec_colex_next3 ( DIM_NUM, base, a, more );

    if ( !more ) 
    {
      break;
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test155 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST155 tests VEC_GRAY_NEXT, VEC_GRAY_RANK and VEC_GRAY_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int a[N];
  int base[N] = { 2, 2, 1, 4 };
  int change;
  bool done;
  int i;
  int prod;
  int rank;

  prod = 1;
  for ( i = 0; i < N; i++ )
  {
    prod = prod * base[i];
  }

  cout << "\n";
  cout << "TEST155\n";
  cout << "  VEC_GRAY_NEXT generates product space elements.\n";
  cout << "  VEC_GRAY_RANK ranks them.\n";
  cout << "  VEC_GRAY_UNRANK unranks them.\n";
  cout << "\n";
  cout << "  The number of components is " << N << "\n";
  cout << "  The number of elements is " << prod << "\n";
  cout << "  Each component has its own number of degrees of\n";
  cout << "  freedom, which, for this example, are:\n";
  cout << "\n";
  cout << "  Rank Change     ";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(4) << base[i] << "  ";
  }
  cout << "\n";
  cout << "\n";

  rank = 0;
  done = true;

  for ( ; ; )
  {
    rank = rank + 1;

    vec_gray_next ( N, base, a, done, change );

    if ( done )
    {
      break;
    }

    cout << setw(4) << rank
         << setw(4) << change;
    for ( i = 0; i < N; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

  }

  for ( i = 0; i < N; i++ )
  {
    a[i] = base[i] / 2;
  }

  rank = vec_gray_rank ( N, base, a );

  cout << "\n";
  cout << "  VEC_GRAY_RANK reports the element\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(4) << a[i] << "  ";
  }
  cout << "\n";
  cout << "\n";
  cout << "  has rank " << rank << "\n";

  rank = 7;
  vec_gray_unrank ( N, base, rank, a );

  cout << "\n";
  cout << "  VEC_GRAY_UNRANK reports the element of rank " << rank << "  is:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(4) << a[i] << "  ";
  }
  cout << "\n";

  return;
# undef N
}
//****************************************************************************80

void test154 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST154 tests VEC_LEX_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int a[DIM_NUM];
  int base = 3;
  int i;
  bool more;

  cout << "\n";
  cout << "TEST154\n";
  cout << "  VEC_LEX_NEXT generates all DIM_NUM-vectors\n";
  cout << "  in a given base.  Here we use base " << base << "\n";
  cout << "\n";

  more = false;

  for ( ; ; )
  {
    vec_lex_next ( DIM_NUM, base, a, more );

    if ( !more ) 
    {
      break;
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(4) << a[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test156 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST156 tests VEC_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int a[N];
  int base;
  int i;
  int j;
  int seed;

  base = 3;
  seed = 123456789;

  cout << "\n";
  cout << "TEST156\n";
  cout << "  VEC_RANDOM generates a random N-vector\n";
  cout << "  in a given base.\n";
  cout << "  Here, we use base " << base << "\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    vec_random ( N, base, seed, a );

    cout << setw(4) << i << "    ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(4) << a[j] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test1565 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1565 tests VECTOR_CONSTRAINED_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int constraint;
  int i;
  int j;
  bool more;
  int x[N];
  int x_max[N] = { 4, 5, 3 };
  int x_min[N] = { 2, 2, 1 };
  int x_prod;

  cout << "\n";
  cout << "TEST1565\n";
  cout << "  VECTOR_CONSTRAINED_NEXT:\n";
  cout << "  Consider vectors:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    P = Product X_MAX(1:N)\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    sum ( (X(1:N)-1) * P / X_MAX(1:N) ) <= P\n";

  more = false;

  cout << "\n";
  cout << "  X_MIN:\n";
  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(4) << x_min[j];
  }
  cout << "\n";
  cout << "\n";
  cout << "  X_MAX:\n";
  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(4) << x_max[j];
  }
  cout << "\n";

  i = 0;

  x_prod = 1;
  for ( j = 0; j < N; j++ )
  {
    x_prod = x_prod * x_max[j];
  }

  cout << "\n";
  cout << "  Maximum allowed CONSTRAINT = P = " << x_prod << "\n";
  cout << "\n";

  for ( ; ; )
  {
    vector_constrained_next ( N, x_min, x_max, x, constraint, more );

    if ( !more )
    {
      break;
    }

    i = i + 1;
    cout << "  " << setw(8) << i;
    cout << "  " << setw(12) << constraint;
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(8) << x[j];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test1566 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1566 tests VECTOR_CONSTRAINED_NEXT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  int constraint;
  int i;
  int j;
  bool more;
  int n;
  int x[N_MAX];
  int x_max[N_MAX] = { 5, 6, 4 };
  int x_min[N_MAX] = { 1, 1, 1 };
  int x_prod;

  cout << "\n";
  cout << "TEST1566\n";
  cout << "  VECTOR_CONSTRAINED_NEXT2:\n";
  cout << "  Consider vectors:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    P = Product X_MAX(1:N)\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    sum ( X(1:N) * P / X_MAX(1:N) ) <= P\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "  X_MIN:\n";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_min[j];
    }
    cout << "\n";
    cout << "\n";
    cout << "  X_MAX:\n";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_max[j];
    }
    cout << "\n";

    i = 0;

    x_prod = 1;
    for ( j = 0; j < n; j++ )
    {
      x_prod = x_prod * x_max[j];
    }

    cout << "\n";
    cout << "  Maximum allowed CONSTRAINT = P = " << x_prod << "\n";
    cout << "\n";

    for ( ; ; )
    {
      vector_constrained_next2 ( n, x_min, x_max, x, constraint, more );

      if ( !more )
      {
        break;
      }

      i = i + 1;
      cout << "  " << setw(8) << i;
      cout << "  " << setw(12) << constraint;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test1567 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1567 tests VECTOR_CONSTRAINED_NEXT3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 April 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double constraint;
  int i;
  int j;
  bool more;
  int n;
  int x[N_MAX];
  int x_max[N_MAX] = { 5, 6, 4 };
  int x_min[N_MAX] = { 1, 1, 1 };

  cout << "\n";
  cout << "TEST1567\n";
  cout << "  VECTOR_CONSTRAINED_NEXT3:\n";
  cout << "  Consider vectors:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    CONSTRAINT = sum ( X(1:N) / X_MAX(1:N) )\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    CONSTRAINT <= 1\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "  X_MIN:\n";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_min[j];
    }
    cout << "\n";
    cout << "\n";
    cout << "  X_MAX:\n";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_max[j];
    }
    cout << "\n";
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      vector_constrained_next3 ( n, x_min, x_max, x, constraint, more );

      if ( !more )
      {
        break;
      }

      i = i + 1;
      cout << "  " << setw(8) << i;
      cout << "  " << setw(14) << constraint;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test1568 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1568 tests VECTOR_CONSTRAINED_NEXT4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double alpha[N_MAX] = { 4.0, 3.0, 5.0 };
  int i;
  int j;
  bool more;
  int n;
  double q = 20.0;
  double total;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };
  int x_min[N_MAX] = { 1, 0, 1 };

  cout << "\n";
  cout << "TEST1568\n";
  cout << "  VECTOR_CONSTRAINED_NEXT4:\n";
  cout << "  Consider vectors:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    TOTAL = sum ( ALPHA(1:N) * X(1:N) )\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    TOTAL <= Q\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "  ALPHA:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << alpha[j];
    }
    cout << "\n";
    cout << "  Q:    ";
    cout << "  " << setw(8) << q;
    cout << "\n";
    cout << "  X_MIN:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_min[j];
    }
    cout << "\n";
    cout << "  X_MAX:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_max[j];
    }
    cout << "\n";
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      vector_constrained_next4 ( n, alpha, x_min, x_max, x, q, more );

      if ( !more )
      {
        break;
      }

      total = 0.0;
      for ( j = 0; j < n; j++ )
      {
        total = total + alpha[j] * ( double ) x[j];
      }
      i = i + 1;
      cout << "  " << setw(8) << i;
      cout << "  " << setw(14) << total;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test1569 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1569 tests VECTOR_CONSTRAINED_NEXT5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int i;
  int j;
  bool more;
  int sum_max;
  int sum_min;
  int x[N];

  cout << "\n";
  cout << "TEST1569\n";
  cout << "  VECTOR_CONSTRAINED_NEXT5:\n";
  cout << "  Generate integer vectors X such that:\n";
  cout << "    SUM_MIN <= sum ( X(1:N) ) <= SUM_MAX,\n";
  cout << "  We require every X(I) to be at least 1.\n";

  more = false;
  sum_min = 5;
  sum_max = 7;

  cout << "\n";
  cout << "  N =       " << N << "\n";
  cout << "  SUM_MIN = " << sum_min << "\n";
  cout << "  SUM_MAX = " << sum_max << "\n";
  cout << "\n";
  cout << "         #        X(1)      X(2)      X(3)\n";
  cout << "\n";

  i = 0;

  for ( ; ; )
  {
    vector_constrained_next5 ( N, x, sum_min, sum_max, more );

    if ( !more )
    {
      break;
    }

    i = i + 1;
    cout << "  " << setw(8) << i;
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(8) << x[j];
    }
    cout << "\n";
  }
  return;
# undef N
}
//****************************************************************************80

void test15695 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15695 tests VECTOR_CONSTRAINED_NEXT6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double alpha[N_MAX] = { 4.0, 3.0, 5.0 };
  int i;
  int j;
  bool more;
  int n;
  double q_max = 20.0;
  double q_min = 16.0;
  double total;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };
  int x_min[N_MAX] = { 1, 0, 1 };

  cout << "\n";
  cout << "TEST15695\n";
  cout << "  VECTOR_CONSTRAINED_NEXT6:\n";
  cout << "  Consider vectors:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    TOTAL = sum ( ALPHA(1:N) * X(1:N) )\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    Q_MIN <= TOTAL <= Q_MAX\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "  ALPHA:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << alpha[j];
    }
    cout << "\n";
    cout << "  Q_MIN:";
    cout << "  " << setw(8) << q_min;
    cout << "\n";
    cout << "  Q_MAX:";
    cout << "  " << setw(8) << q_max;
    cout << "\n";
    cout << "  X_MIN:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_min[j];
    }
    cout << "\n";
    cout << "  X_MAX:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_max[j];
    }
    cout << "\n";
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      vector_constrained_next6 ( n, alpha, x_min, x_max, x, q_min, 
        q_max, more );

      if ( !more )
      {
        break;
      }

      total = 0.0;
      for ( j = 0; j < n; j++ )
      {
        total = total + alpha[j] * ( double ) x[j];
      }
      i = i + 1;
      cout << "  " << setw(8) << i;
      cout << "  " << setw(14) << total;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test15696 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15696 tests VECTOR_CONSTRAINED_NEXT7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double alpha[N_MAX] = { 4.0, 3.0, 5.0 };
  int i;
  int j;
  bool more;
  int n;
  double q_max = 20.0;
  double q_min = 16.0;
  double total;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };

  cout << "\n";
  cout << "TEST15696\n";
  cout << "  VECTOR_CONSTRAINED_NEXT7:\n";
  cout << "  Consider vectors:\n";
  cout << "    0 <= X(1:N) <= X_MAX(1:N),\n";
  cout << "  Set\n";
  cout << "    TOTAL = sum ( ALPHA(1:N) * X(1:N) )\n";
  cout << "  Accept only vectors for which:\n";
  cout << "    Q_MIN <= TOTAL <= Q_MAX\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "  ALPHA:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << alpha[j];
    }
    cout << "\n";
    cout << "  Q_MIN:";
    cout << "  " << setw(8) << q_min;
    cout << "\n";
    cout << "  Q_MAX:";
    cout << "  " << setw(8) << q_max;
    cout << "\n";
    cout << "  X_MAX:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << x_max[j];
    }
    cout << "\n";
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      vector_constrained_next7 ( n, alpha, x_max, x, q_min, 
        q_max, more );

      if ( !more )
      {
        break;
      }

      total = 0.0;
      for ( j = 0; j < n; j++ )
      {
        total = total + alpha[j] * ( double ) x[j];
      }
      i = i + 1;
      cout << "  " << setw(8) << i;
      cout << "  " << setw(14) << total;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test15698 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15698 tests VECTOR_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  int i;
  int j;
  bool more;
  int n;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };
  int x_min[N_MAX] = { 1, 4, 3 };

  cout << "\n";
  cout << "TEST15698\n";
  cout << "  VECTOR_NEXT:\n";
  cout << "  Generate all vectors X such that:\n";
  cout << "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n";

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = false;

    cout << "\n";
    cout << "    X_MIN:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << x_min[j];
    }
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      vector_next ( n, x_min, x_max, x, more );

      if ( !more )
      {
        break;
      }

      i = i + 1;
      cout << "  " << setw(8) << i;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(8) << x[j];
      }
      cout << "\n";
    }
    cout << "    X_MAX:";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << x_max[j];
    }
    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test157 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST157 tests YTB_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST157\n";
  cout << "  YTB_ENUM counts Young tableau.\n";
  cout << "\n";
  cout << "   N  YTB_ENUM(N)\n";
  cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    cout << setw(4)  << n              << "  "
         << setw(10) << ytb_enum ( n ) << "\n";
  }

  return;
}
//****************************************************************************80

void test158 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST158 tests YTB_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int i;
  int a[N];
  int lambda[N] = { 3, 2, 1, 0, 0, 0 };
  bool more;

  for ( i = 0; i < N; i++ )
  {
    a[i] = 0;
  }

  cout << "\n";
  cout << "TEST158\n";
  cout << "  YTB_NEXT generates Young tableaus.\n";
  cout << "\n";

  more = false;

  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    ytb_next ( N, lambda, a, more );

    ytb_print ( N, a, " " );

    if ( !more || 100 < i )
    {
      break;
    }

  }

  return;

# undef N
}
//****************************************************************************80

void test159 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST159 tests YTB_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int a[N];
  int i;
  int lambda[N] = { 3, 2, 1, 0, 0, 0 };
  int seed;

  cout << "\n";
  cout << "TEST159\n";
  cout << "  YTB_RANDOM generates a random Young tableau\n";

  seed = 123456789;

  for ( i = 1; i <=5; i++ )
  {
    ytb_random ( N, lambda, seed, a );

    ytb_print ( N, a, " " );

  }

  return;
# undef N
}
