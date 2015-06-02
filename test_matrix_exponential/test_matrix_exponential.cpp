# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_matrix_exponential.hpp"
# include "c8lib.hpp"
# include "r8lib.hpp"

//****************************************************************************80

complex <double> *c8mat_exp_a ( int test, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_EXP_A returns the matrix for a given complex test.
//
//  Discussion:
//
//    1) diagonal matrix, real.
//    2) diagonal matrix, pure imaginary.
//    3) diagonal matrix, complex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <double> C8MAT_EXP_A[N*N], the matrix.
//
{
  complex <double> *a;
  static complex <double> a01[2*2] = {
      1.0, 0.0, 
      0.0, 2.0 };
  static complex <double> a02[2*2] = {
      complex <double> ( 0.0, 3.0), 
      0.0,
      0.0,
      complex <double> ( 0.0, -4.0 ) };
  static complex <double> a03[2*2] = {
      complex <double> ( 5.0, 6.0 ), 
      0.0,
      0.0,
      complex <double> ( 7.0, -8.0 ) };
  int i;
  int j;

  if ( test == 1 )
  {
    a = c8mat_copy_new ( n, n, a01 );
  }
  else if ( test == 2 )
  {
    a = c8mat_copy_new ( n, n, a02 );
  }
  else if ( test == 3 )
  {
    a = c8mat_copy_new ( n, n, a03 );
  }
  else
  {
    cerr << "\n";
    cerr << "C8MAT_EXP_A - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return a;
}
//****************************************************************************80

complex <double> *c8mat_exp_expa ( int test, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_EXP_EXPA returns the "exact" exponential matrix for a given complex test.
//
//  Discussion:
//
//    1) diagonal matrix, real.
//    2) diagonal matrix, pure imaginary.
//    3) diagonal matrix, complex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <double> C8MAT_EXP_EXPA[N*N], the exponential of the test matrix.
//
{
  complex <double> *expa;
  static complex <double> expa01[2*2] = {
      complex <double> ( 2.718281828459046, 0.0 ),
      0.0,
      0.0,
      complex <double> ( 7.389056098930650, 0.0 ) };
  static complex <double> expa02[2*2] = {
      complex <double> ( -0.989992496600446, 0.141120008059867 ),
      0.0,
      0.0,
      complex <double> ( -0.653643620863612, 0.756802495307928 ) };
  static complex <double> expa03[2*2] = {
      complex <double> ( 142.501905518208, - 41.468936789923 ),
      0.0,
      0.0,
      complex <double> ( -159.560161626987, - 1084.963058811836 ) };

  if ( test == 1 )
  {
    expa = c8mat_copy_new ( n, n, expa01 );
  }
  else if ( test == 2 )
  {
    expa = c8mat_copy_new ( n, n, expa02 );
  }
  else if ( test == 3 )
  {
    expa = c8mat_copy_new ( n, n, expa03 );
  }
  else
  {
    cerr << "\n";
    cerr << "C8MAT_EXP_EXPA - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return expa;
}
//****************************************************************************80

int c8mat_exp_n ( int test )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_EXP_N returns the matrix order for a given complex test.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Output, int C8MAT_EXP_N, the order of the matrix.
//
{
  int n;

  if ( test == 1 )
  {
    n = 2;
  }
  else if ( test == 2 )
  {
    n = 2;
  }
  else if ( test == 3 )
  {
    n = 2;
  }
  else
  {
    cerr << "\n";
    cerr << "C8MAT_EXP_N - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }

  return n;
}
//****************************************************************************80

void c8mat_exp_story ( int test )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_EXP_STORY prints explanatory text for each complex problem.
//
//  Discussion:
//
//    1) diagonal matrix, real.
//    2) diagonal matrix, pure imaginary.
//    3) diagonal matrix, complex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
{
  if ( test == 1 )
  {
    cout << "\n";
    cout << "  This matrix is diagonal.\n";
    cout << "  The diagonal entries are real.\n";
  }
  else if ( test == 2 )
  {
    cout << "\n";
    cout << "  This matrix is diagonal.\n";
    cout << "  The diagonal entries are pure imaginary.\n";
  }
  else if ( test == 3 )
  {
    cout << "\n";
    cout << "  This matrix is diagonal.\n";
    cout << "  The diagonal entries are complex.\n";
  }
  else
  {
    cerr << "\n";
    cerr << "C8MAT_EXP_STORY - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int c8mat_exp_test_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_EXP_TEST_NUM returns the number of complex tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int C8MAT_EXP_TEST_NUM, the number of tests.
//
{
  int test_num;

  test_num = 3;

  return test_num;
}
//****************************************************************************80

double *r8mat_exp_a ( int test, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXP_A returns the matrix for a given real test.
//
//  Discussion:
//
//     1) Diagonal example
//     2) Symmetric example
//     3) Laub
//     4) Moler and Van Loan
//     5) Moler and Van Loan
//     6) Moler and Van Loan
//     7) Moler and Van Loan
//     8) Wikipedia example
//     9) NAG F01ECF
//    10) Ward #1
//    11) Ward #2
//    12) Ward #3
//    13) Ward #4
//    14) Moler example
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Laub,
//    Review of "Linear System Theory" by Joao Hespanha,
//    SIAM Review,
//    Volume 52, Number 4, December 2010, page 779-781.
//
//    Cleve Moler, Charles VanLoan,
//    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
//    Twenty-Five Years Later,
//    SIAM Review,
//    Volume 45, Number 1, March 2003, pages 3-49.
//
//    Cleve Moler,
//    Cleve's Corner: A Balancing Act for the Matrix Exponential,
//    July 23rd, 2012.
//
//    Robert Ward,
//    Numerical computation of the matrix exponential with accuracy estimate,
//    SIAM Journal on Numerical Analysis,
//    Volume 14, Number 4, September 1977, pages 600-610.
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Input, int N, the order of the matrix.
//
//    Output, double R8MAT_EXP_A[N*N], the matrix.
//
{
  double *a;
  static double a01[2*2] = {
      1.0, 0.0, 
      0.0, 2.0 };
  static double a02[2*2] = {
      1.0, 3.0, 
      3.0, 2.0 };
  static double a03[2*2] = {
      0.0, -39.0, 
      1.0, -40.0 };
  static double a04[2*2] = {
      -49.0, -64.0, 
       24.0,  31.0 };
  static double a05[4*4] = {
      0.0, 0.0, 0.0, 0.0, 
      6.0, 0.0, 0.0, 0.0, 
      0.0, 6.0, 0.0, 0.0, 
      0.0, 0.0, 6.0, 0.0 };
  static double a06[2*2] = {
      1.0, 0.0, 
      1.0, 1.0 };
  static double a08[3*3] = {
      21.0,  -5.0,   4.0, 
      17.0,  -1.0,   4.0, 
       6.0,  -6.0,  16.0 };
  static double a09[4*4] = {
      1.0, 3.0, 3.0, 3.0, 
      2.0, 1.0, 2.0, 3.0, 
      2.0, 1.0, 1.0, 3.0, 
      2.0, 2.0, 2.0, 1.0 };
  static double a10[3*3] = {
      4.0, 1.0, 1.0, 
      2.0, 4.0, 1.0, 
      0.0, 1.0, 4.0 };
  static double a11[3*3] = {
      29.87942128909879, 
       0.7815750847907159, 
      -2.289519314033932, 
       0.7815750847907159, 
      25.72656945571064, 
       8.680737820540137, 
      -2.289519314033932, 
       8.680737820540137, 
      34.39400925519054 };
  static double a12[3*3] = {
      -131.0, -390.0, -387.0, 
        19.0,   56.0,   57.0, 
        18.0,   54.0,   52.0 };
  int i;
  int j;

  if ( test == 1 )
  {
    a = r8mat_copy_new ( n, n, a01 );
  }
  else if ( test == 2 )
  {
    a = r8mat_copy_new ( n, n, a02 );
  }
  else if ( test == 3 )
  {
    a = r8mat_copy_new ( n, n, a03 );
  }
  else if ( test == 4 )
  {
    a = r8mat_copy_new ( n, n, a04 );
  }
  else if ( test == 5 )
  {
    a = r8mat_copy_new ( n, n, a05 );
  }
  else if ( test == 6 )
  {
    a = r8mat_copy_new ( n, n, a06 );
  }
  else if ( test == 7 )
  {
    a = new double[2*2];
    a[0+0*2] = 1.0 + r8_epsilon ( );
    a[1+0*2] = 0.0;
    a[0+1*2] = 0.0;
    a[1+1*2] = 1.0 - r8_epsilon ( );
  }
  else if ( test == 8 )
  {
    a = r8mat_copy_new ( n, n, a08 );
  }
  else if ( test == 9 )
  {
    a = r8mat_copy_new ( n, n, a09 );
  }
  else if ( test == 10 )
  {
    a = r8mat_copy_new ( n, n, a10 );
  }
  else if ( test == 11 )
  {
    a = r8mat_copy_new ( n, n, a11 );
  }
  else if ( test == 12 )
  {
    a = r8mat_copy_new ( n, n, a12 );
  }
  else if ( test == 13 )
  {
    a = new double[n*n];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        if ( j == i + 1 )
        {
          a[i+j*n] = 1.0;
        }
        else if ( i == n - 1 && j == 0 )
        {
          a[i+j*n] = 1.0E-10;
        }
        else
        {
          a[i+j*n] = 0.0;
        }
      }
    }
  }
  else if ( test == 14 )
  {
    a = new double[n*n];

    a[0+0*3] = 0.0;
    a[0+1*3] = 1.0E-08;
    a[0+2*3] = 0.0;
    a[1+0*3] = - 2.0E+10 - 2.0E+08 / 3.0;
    a[1+1*3] = - 3.0;
    a[1+2*3] = 2.0E+10;
    a[2+0*3] = 200.0E+00 / 3.0;
    a[2+1*3] = 0.0;
    a[2+2*3] = - 200.0E+00 / 3.0;
  }
  else
  {
    cerr << "\n";
    cerr << "R8MAT_EXP_A - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return a;
}
//****************************************************************************80

double *r8mat_exp_expa ( int test, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXP_EXPA returns the "exact" exponential matrix for a given real test.
//
//  Discussion:
//
//    In some cases, the "exact" value is given to six significant digits.
//
//     1) Diagonal example
//     2) Symmetric example
//     3) Laub
//     4) Moler and Van Loan
//     5) Moler and Van Loan
//     6) Moler and Van Loan
//     7) Moler and Van Loan
//     8) Wikipedia example
//     9) NAG F01ECF
//    10) Ward #1
//    11) Ward #2
//    12) Ward #3
//    13) Ward #4
//    14) Moler example
//
//    Thanks to Alex Griffing for correcting the value of matrix 3,
//    17 October 2012.
//
//    Thanks again to Alex Griffing for providing improved values for
//    matrices 4, 7 and 13, 03 September 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Laub,
//    Review of "Linear System Theory" by Joao Hespanha,
//    SIAM Review,
//    Volume 52, Number 4, December 2010, page 779-781.
//
//    Cleve Moler, Charles VanLoan,
//    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
//    Twenty-Five Years Later,
//    SIAM Review,
//    Volume 45, Number 1, March 2003, pages 3-49.
//
//    Cleve Moler,
//    Cleve's Corner: A Balancing Act for the Matrix Exponential,
//    July 23rd, 2012.
//
//    Robert Ward,
//    Numerical computation of the matrix exponential with accuracy estimate,
//    SIAM Journal on Numerical Analysis,
//    Volume 14, Number 4, September 1977, pages 600-610.
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Input, int N, the order of the matrix.
//
//    Output, double R8MAT_EXP_EXPA[N*N], the exponential of the test matrix.
//
{
  double exp16;
  double exp4;
  double *expa;
  static double expa01[2*2] = {
      2.718281828459046, 0.0, 
      0.0,               7.389056098930650 };
  static double expa02[2*2] = {
      39.322809708033859,  46.166301438885753, 
      46.166301438885768,  54.711576854329110 };
  static double expa03[2*2] = {
       0.37756048,  0.00968104,
      -0.37756048, -0.00968104 };
  static double expa04[2*2] = {
      -0.7357587581447531, -1.4715175990882605, 
       0.5518190996580977,  1.1036382407155727 };
  static double expa05[4*4] = {
      1.0,  0.0, 0.0, 0.0, 
      6.0,  1.0, 0.0, 0.0, 
     18.0,  6.0, 1.0, 0.0, 
     36.0, 18.0, 6.0, 1.0 };
  static double expa06[2*2] = {
      2.718281828459046, 0.0, 
      2.718281828459046, 2.718281828459046 };
  static double expa07[2*2] = {
      2.718281828459045235360287, 0.0, 
      2.718281828459045235360287, 2.718281828459045235360287 };
  static double expa09[4*4] = {
      740.7038, 731.2510, 823.7630, 998.4355, 
      610.8500, 603.5524, 679.4257, 823.7630, 
      542.2743, 535.0884, 603.5524, 731.2510, 
      549.1753, 542.2743, 610.8500, 740.7038 };
  static double expa10[3*3] = {
      147.8666224463699, 
      127.7810855231823, 
      127.7810855231824, 
      183.7651386463682, 
      183.7651386463682, 
      163.6796017231806, 
      71.79703239999647, 
      91.88256932318415, 
     111.9681062463718 };
  static double expa11[3*3] = {
     5.496313853692378E+15, 
    -1.823188097200899E+16, 
    -3.047577080858001E+16, 
    -1.823188097200898E+16, 
     6.060522870222108E+16, 
     1.012918429302482E+17, 
    -3.047577080858001E+16, 
     1.012918429302482E+17, 
     1.692944112408493E+17 };
  static double expa12[3*3] = {
    -1.509644158793135, 
    -5.632570799891469, 
    -4.934938326088363, 
     0.3678794391096522, 
     1.471517758499875, 
     1.103638317328798, 
     0.1353352811751005, 
     0.4060058435250609, 
     0.5413411267617766 };
  static double expa14[3*3] = {
    4.468494682831735E-01, 
   -5.743067779479621E+06, 
    4.477229778494929E-01, 
    1.540441573839520E-09, 
   -1.528300386868247E-02, 
    1.542704845195912E-09, 
    4.628114535587735E-01, 
   -4.526542712784168E+06, 
    4.634806488376499E-01 };
  int i;
  int j;
  int k;
  double value;

  if ( test == 1 )
  {
    expa = r8mat_copy_new ( n, n, expa01 );
  }
  else if ( test == 2 )
  {
    expa = r8mat_copy_new ( n, n, expa02 );
  }
  else if ( test == 3 )
  {
    expa = r8mat_copy_new ( n, n, expa03 );
  }
  else if ( test == 4 )
  {
    expa = r8mat_copy_new ( n, n, expa04 );
  }
  else if ( test == 5 )
  {
    expa = r8mat_copy_new ( n, n, expa05 );
  }
  else if ( test == 6 )
  {
    expa = r8mat_copy_new ( n, n, expa06 );
  }
  else if ( test == 7 )
  {
    expa = r8mat_copy_new ( n, n, expa07 );
  }
  else if ( test == 8 )
  {
    expa = new double[3*3];
    exp16 = exp ( 16.0 );
    exp4 = exp ( 4.0 );
    expa[0+0*3] = 0.25 * ( 13.0 * exp16 -       exp4 );
    expa[1+0*3] = 0.25 * ( -9.0 * exp16 +       exp4 );
    expa[2+0*3] = 0.25 * ( 16.0 * exp16 );
    expa[0+1*3] = 0.25 * ( 13.0 * exp16 - 5.0 * exp4 );
    expa[1+1*3] = 0.25 * ( -9.0 * exp16 + 5.0 * exp4 );
    expa[2+1*3] = 0.25 * ( 16.0 * exp16 );
    expa[0+2*3] = 0.25 * (  2.0 * exp16 - 2.0 * exp4 );
    expa[1+2*3] = 0.25 * ( -2.0 * exp16 + 2.0 * exp4 );
    expa[2+2*3] = 0.25 * (  4.0 * exp16 );
  }
  else if ( test == 9 )
  {
    expa = r8mat_copy_new ( n, n, expa09 );
  }
  else if ( test == 10 )
  {
    expa = r8mat_copy_new ( n, n, expa10 );
  }
  else if ( test == 11 )
  {
    expa = r8mat_copy_new ( n, n, expa11 );
  }
  else if ( test == 12 )
  {
    expa = r8mat_copy_new ( n, n, expa12 );
  }
  else if ( test == 13 )
  {
    expa = new double[n*n];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        expa[i+j*n] = 0.0;
      }
    }

    k = 0;
    for ( i = 1; i <= n; i++ )
    {
      expa[i-1+(i-1)*n] = 1.0;
    }

    value = 1.0;
    for ( k = 1; k < n; k++ )
    {
      value = value / ( double ) ( k );
      for ( i = 1; i <= n - k; i++ )
      {
        expa[i-1+(i+k-1)*n] = value;
      }
    }

    value = 1.0 / pow ( 10.0, n );
    for ( k = 1; k < n; k++ )
    {
      value = value / ( double ) ( k );
      for ( j = 1; j <= k; j++ )
      {
        expa[n+j-k-1+(j-1)*n] = value;
      }
    }

  }
  else if ( test == 14 )
  {
    expa = r8mat_copy_new ( n, n, expa14 );
  }
  else
  {
    cerr << "\n";
    cerr << "R8MAT_EXP_EXPA - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return expa;
}
//****************************************************************************80

int r8mat_exp_n ( int test )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXP_N returns the matrix order for a given real test.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
//    Output, int R8MAT_EXP_N, the order of the matrix.
//
{
  int n;

  if ( test == 1 )
  {
    n = 2;
  }
  else if ( test == 2 )
  {
    n = 2;
  }
  else if ( test == 3 )
  {
    n = 2;
  }
  else if ( test == 4 )
  {
    n = 2;
  }
  else if ( test == 5 )
  {
    n = 4;
  }
  else if ( test == 6 )
  {
    n = 2;
  }
  else if ( test == 7 )
  {
    n = 2;
  }
  else if ( test == 8 )
  {
    n = 3;
  }
  else if ( test == 9 )
  {
    n = 4;
  }
  else if ( test == 10 )
  {
    n = 3;
  }
  else if ( test == 11 )
  {
    n = 3;
  }
  else if ( test == 12 )
  {
    n = 3;
  }
  else if ( test == 13 )
  {
    n = 10;
  }
  else if ( test == 14 )
  {
    n = 3;
  }
  else
  {
    cerr << "\n";
    cerr << "R8MAT_EXP_N - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }

  return n;
}
//****************************************************************************80

void r8mat_exp_story ( int test )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXP_STORY prints explanatory text for each real problem.
//
//  Discussion:
//
//     1) Diagonal example
//     2) Symmetric example
//     3) Laub
//     4) Moler and Van Loan
//     5) Moler and Van Loan
//     6) Moler and Van Loan
//     7) Moler and Van Loan
//     8) Wikipedia example
//     9) NAG F01ECF
//    10) Ward #1
//    11) Ward #2
//    12) Ward #3
//    13) Ward #4
//    14) Moler example
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Laub,
//    Review of "Linear System Theory" by Joao Hespanha,
//    SIAM Review,
//    Volume 52, Number 4, December 2010, page 779-781.
//
//    Cleve Moler, Charles VanLoan,
//    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
//    Twenty-Five Years Later,
//    SIAM Review,
//    Volume 45, Number 1, March 2003, pages 3-49.
//
//    Cleve Moler,
//    Cleve's Corner: A Balancing Act for the Matrix Exponential,
//    July 23rd, 2012.
//
//    Robert Ward,
//    Numerical computation of the matrix exponential with accuracy estimate,
//    SIAM Journal on Numerical Analysis,
//    Volume 14, Number 4, September 1977, pages 600-610.
//
//  Parameters:
//
//    Input, int TEST, the index of the test case.
//
{
  if ( test == 1 )
  {
    cout << "\n";
    cout << "  This matrix is diagonal.\n";
    cout << "  The calculation of the matrix exponential is simple.\n";
  }
  else if ( test == 2 )
  {
    cout << "\n";
    cout << "  This matrix is symmetric.\n";
    cout << "  The calculation of the matrix exponential is straightforward.\n";
  }
  else if ( test == 3 )
  {
    cout << "\n";
    cout << "  This example is due to Laub.\n";
    cout << "  This matrix is ill-suited for the Taylor series approach.\n";
    cout << "  As powers of A are computed, the entries blow up too quickly.\n";
  }
  else if ( test == 4 )
  {
    cout << "\n";
    cout << "  This example is due to Moler and Van Loan.\n";
    cout << "  The example will cause problems for the series summation approach,\n";
    cout << "  as well as for diagonal Pade approximations.\n";
  }
  else if ( test == 5 )
  {
    cout << "\n";
    cout << "  This example is due to Moler and Van Loan.\n";
    cout << "  This matrix is strictly upper triangular\n";
    cout << "  All powers of A are zero beyond some (low) limit.\n";
    cout << "  This example will cause problems for Pade approximations.\n";
  }
  else if ( test == 6 )
  {
    cout << "\n";
    cout << "  This example is due to Moler and Van Loan.\n";
    cout << "  This matrix does not have a complete set of eigenvectors.\n";
    cout << "  That means the eigenvector approach will fail.\n";
  }
  else if ( test == 7 )
  {
    cout << "\n";
    cout << "  This example is due to Moler and Van Loan.\n";
    cout << "  This matrix is very close to example 5.\n";
    cout << "  Mathematically, it has a complete set of eigenvectors.\n";
    cout << "  Numerically, however, the calculation will be suspect.\n";
  }
  else if ( test == 8 )
  {
    cout << "\n";
    cout << "  This matrix was an example in Wikipedia.\n";
  }
  else if ( test == 9 )
  {
    cout << "\n";
    cout << "  This matrix is due to the NAG Library.\n";
    cout << "  It is an example for function F01ECF.\n";
  }
  else if ( test == 10 )
  {
    cout << "\n";
    cout << "  This is Ward's example #1.\n";
    cout << "  It is defective and nonderogatory.\n";
    cout << "  The eigenvalues are 3, 3 and 6.\n";
  }
  else if ( test == 11 )
  {
    cout << "\n";
    cout << "  This is Ward's example #2.\n";
    cout << "  It is a symmetric matrix.\n";
    cout << "  The eigenvalues are 20, 30, 40.\n";
  }
  else if ( test == 12 )
  {
    cout << "\n";
    cout << "  This is Ward's example #3.\n";
    cout << "  Ward's algorithm has difficulty estimating the accuracy\n";
    cout << "  of its results.  The eigenvalues are -1, -2, -20.\n";
  }
  else if ( test == 13 )
  {
    cout << "\n";
    cout << "  This is Ward's example #4.\n";
    cout << "  This is a version of the Forsythe matrix.\n";
    cout << "  The eigenvector problem is badly conditioned.\n";
    cout << "  Ward's algorithm has difficulty estimating the accuracy\n";
    cout << "  of its results for this problem.\n";
  }
  else if ( test == 14 )
  {
    cout << "\n";
    cout << "  This is Moler's example.\n";
    cout << "  This badly scaled matrix caused problems for MATLAB's expm().\n";
  }
  else
  {
    cerr << "\n";
    cerr << "R8MAT_EXP_STORY - Fatal error!\n";
    cerr << "  Illegal value of TEST = " << test << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int r8mat_exp_test_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXP_TEST_NUM returns the number of real tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int R8MAT_EXP_TEST_NUM, the number of tests.
//
{
  int test_num;

  test_num = 14;

  return test_num;
}
