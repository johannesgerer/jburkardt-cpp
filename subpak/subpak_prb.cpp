# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <cstring>

using namespace std;

# include "subpak.hpp"

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
void test225 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test29 ( );

void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test35 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUBPAK_PRB.
//
//  Discussion:
//
//    SUBPAK_PRB tests the SUBPAK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SUBPAK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SUBPAK library.\n";

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
  test225 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test35 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SUBPAK_PRB\n";
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
//    TEST01 tests ANGLE_SHIFT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double angle_hi;
  double angle_lo;
  double beta;
  double gamma;
  double pi = 3.141592653589793;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ANGLE_SHIFT shifts an angle by multiples of\n";
  cout << "  2 Pi until it lies between BETA and BETA+2Pi.\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA   BETA+2Pi\n";
  cout << "\n";

  angle_lo = -4.0 * pi;
  angle_hi = +4.0 * pi;

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    alpha = r8_uniform ( angle_lo, angle_hi, &seed );

    beta = r8_uniform ( angle_lo, angle_hi, &seed );

    gamma = angle_shift ( alpha, beta );

    cout << "  " << setw(8) << alpha
         << "  " << setw(8) << beta
         << "  " << setw(8) << gamma
         << "  " << setw(8) << beta + 2.0 * pi << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ANGLE_SHIFT_DEG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double angle_hi;
  double angle_lo;
  double beta;
  double gamma;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  ANGLE_SHIFT_DEG shifts an angle by multiples of\n";
  cout << "  360 until it lies between BETA and BETA+360.\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA   BETA+360\n";
  cout << "\n";

  angle_lo = -720.0;
  angle_hi = +720.0;

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    alpha = r8_uniform ( angle_lo, angle_hi, &seed );

    beta = r8_uniform ( angle_lo, angle_hi, &seed );

    gamma = angle_shift_deg ( alpha, beta );

    cout << "  " << setw(8) << alpha
         << "  " << setw(8) << beta
         << "  " << setw(8) << gamma
         << "  " << setw(8) << beta + 360.0 << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests ANGLE_TO_RGB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  double angle_lo = 0.0;
  double angle_hi = 360.0;
  int i;
  double *rgb;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  ANGLE_TO_RGB converts an angle into an RGB color.\n";
  cout << "\n";
  cout << "     ANGLE         R         G         B\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 0; test < test_num; test++ )
  {
    angle = r8_uniform ( angle_lo, angle_hi, &seed );

    rgb = angle_to_rgb ( angle );

    cout << "  " << setw(8) << angle
         << "  " << setw(8) << rgb[0]
         << "  " << setw(8) << rgb[1]
         << "  " << setw(8) << rgb[2] << "\n";

    delete [] rgb;
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests AXIS_LIMITS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int ndivs;
  int nticks;
  double pxdiv;
  double pxmax;
  double pxmin;
  double xmax;
  double xmin;

  xmin = 67.3;
  xmax = 114.7;
  ndivs = 6;

  axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

  cout << "\n";
  cout << "TEST04\n";
  cout << "  AXIS_LIMITS adjusts plot limits to \"nicer\" values.\n";
  cout << "\n";
  cout << "  Input XMIN =    " << xmin   << "\n";
  cout << "  Input XMAX =    " << xmax   << "\n";
  cout << "  Input NDIVS =   " << ndivs  << "\n";
  cout << "\n";
  cout << "  Output PXMIN =  " << pxmin  << "\n";
  cout << "  Output PXMAX =  " << pxmax  << "\n";
  cout << "  Output PXDIV =  " << pxdiv  << "\n";
  cout << "  Output NTICKS = " << nticks << "\n";

  xmin = -26.0;
  xmax = +26.0;
  ndivs = 10;

  axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

  cout << "\n";
  cout << "  Input XMIN =    " << xmin   << "\n";
  cout << "  Input XMAX =    " << xmax   << "\n";
  cout << "  Input NDIVS =   " << ndivs  << "\n";
  cout << "\n";
  cout << "  Output PXMIN =  " << pxmin  << "\n";
  cout << "  Output PXMAX =  " << pxmax  << "\n";
  cout << "  Output PXDIV =  " << pxdiv  << "\n";
  cout << "  Output NTICKS = " << nticks << "\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests AXIS_LIMITS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  int i;
  int ndivs;
  int nticks;
  double pxdiv;
  double pxmax;
  double pxmin;
  double test_max[TEST_NUM] = {
    9.0, 4.125, 193.75, 2000.250, 12.0 };
  double test_min[TEST_NUM] = {
    1.0, 1.003, 101.25, 2000.125, -7.0 };
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  AXIS_LIMITS computes \"nice\" limits for a graph\n";
  cout << "  that must include a given range.\n";

  ndivs = 5;

  cout << "\n";
  cout << "  All tests use NDIVS = " << ndivs << "\n";
  cout << "\n";
  cout << "          XMIN          XMAX         PXMIN"
       << "         PXMAX         PXDIV  NTICKS\n";
  cout << "\n";

  for ( i = 0; i < TEST_NUM; i++ )
  {

    xmin = test_min[i];
    xmax = test_max[i];

    axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

    cout                       << "  "
         << setw(12) << xmin   << "  "
         << setw(12) << xmax   << "  "
         << setw(12) << pxmin  << "  "
         << setw(12) << pxmax  << "  "
         << setw(12) << pxdiv  << "  "
         << setw(6)  << nticks << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests BAR_CHECK, BAR_CODE, BAR_DIGIT_CODE_LEFT, BAR_DIGIT_CODE_RIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  char *bar;
  int check;
  int digit[12];
  char *codel;
  char *coder;
  int i;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  BAR_CHECK checks digits for a barcode;\n";
  cout << "  BAR_CODE computes the barcode for a string of 11 digits;\n";
  cout << "  BAR_DIGIT_CODE_LEFT returns the left digit code.\n";
  cout << "  BAR_DIGIT_CODE_RIGHT returns the right digit code.\n";

  for ( i = 0; i <= 10; i++ )
  {
    digit[i] = ( i % 10 );
  }
 
  check = bar_check ( digit );
 
  cout << "\n";
  cout << "  The check digit is " << check << "\n";

  digit[11] = check;
 
  cout << "\n";
  cout << "  The left and right digit codes:\n";
  cout << "\n";

  for( i = 0; i <= 9; i++ )
  {
    codel = bar_digit_code_left ( i );
    coder = bar_digit_code_right ( i );
    cout << "  " << setw(2) << i
         << "  " << codel
         << "  " << coder << "\n";
    delete [] codel;
    delete [] coder;
  }
 
  bar = bar_code ( digit );
 
  cout << "\n";
  cout << "  Bar code:\n";
  cout << "\n";
  for ( i = 0; i < 9; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 9; i < 12; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 12; i < 19; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";

  for ( i = 19; i < 26; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 16; i < 33; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 33; i < 40; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 40; i < 47; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 47; i < 54; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 54; i < 59; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 59; i < 66; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 66; i < 73; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 73; i < 80; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 80; i < 87; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 87; i < 94; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 94; i < 101; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 101; i < 104; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";
  for ( i = 104; i < 113; i++ )
  {
    cout << bar[i];
  }
  cout << "\n";

  delete [] bar;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests BMI_ENGLISH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  double b;
  double bmi;
  double c;
  double h;
  double h_ft;
  double h_in;
  int seed;
  int test;
  double w;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  BMI_ENGLISH computes the Body Mass Index\n";
  cout << "  given body measurements in English Units.\n";
  cout << "\n";
  cout << "      Weight               Height            BMI\n";
  cout << "       (LB)          (FT          IN)\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = 100.0;
    c = 250.0;

    w = r8_uniform ( b, c, &seed );

    b = 4.0;
    c = 6.75;

    h = r8_uniform ( b, c, &seed );

    h_ft = ( int ) ( h );
    h_in = ( double ) ( ( int ) ( 12.0 * ( h - h_ft ) ) );
 
    bmi = bmi_english ( w, h_ft, h_in );
    cout << "  " << setw(10) << w
         << "  " << setw(10) << h_ft
         << "  " << setw(10) << h_in
         << "  " << setw(10) << bmi << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests FAC_DIV, FAC_GCD, FAC_LCM, FAC_MUL, FAC_TO_I4, I4_TO_FAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define PRIME_NUM 5

  int bot;
  int i1;
  int i2;
  int *npower1;
  int *npower2;
  int npower3[PRIME_NUM];
  int top;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For products of prime factors:\n";
  cout << "  FAC_DIV computes a quotient;\n";
  cout << "  FAC_MUL multiplies;\n";
  cout << "  FAC_LCM computes the LCM;\n";
  cout << "  FAC_GCD computes the GCD;\n";
  cout << "  I4_TO_FAC converts an integer;\n";
  cout << "  FAC_TO_I4 converts to an integer.\n";
  cout << "  FAC_TO_RAT converts to a ratio.\n";

  i1 = 720;
  i2 = 42;

  npower1 = i4_to_fac ( i1, PRIME_NUM );

  cout << "\n";
  cout << "  Representation of I1 = " << i1 << "\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower1 );

  npower2 = i4_to_fac ( i2, PRIME_NUM );

  cout << "\n";
  cout << "  Representation of I2 = " << i2 << "\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower2 );

  fac_lcm ( PRIME_NUM, npower1, npower2, npower3 );

  cout << "\n";
  cout << "  LCM of I1, I2:\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower3 );

  fac_gcd ( PRIME_NUM, npower1, npower2, npower3 );

  cout << "\n";
  cout << "  GCD of I1, I2:\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower3 );

  fac_mul ( PRIME_NUM, npower1, npower2, npower3 );

  cout << "\n";
  cout << "  Product of I1, I2:\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower3 );

  fac_div ( PRIME_NUM, npower2, npower1, npower3 );

  cout << "\n";
  cout << "  Quotient of I2 / I1:\n";
  cout << "\n";

  fac_print ( PRIME_NUM, npower3 );

  fac_to_rat ( PRIME_NUM, npower3, &top, &bot );

  cout << "\n";
  cout << "  Quotient as a rational: " << top << " / " << bot << "\n";

  delete [] npower1;
  delete [] npower2;

  return;
# undef PRIME_NUM
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests GAUSS_SUM
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
# define DIM_NUM 2
# define N 3

  double amplitude[N] = { 10.0, 5.0, -3.0 };
  double center[DIM_NUM*N] = { 2.0, 3.0,  5.0, 8.0,  7.0, 5.0 };
  double gxy;
  int i;
  int j;
  double width[N] = { 1.0, 2.0, 4.0 };
  double x[DIM_NUM];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  GAUSS_SUM evaluates a function which is the sum of\n";
  cout << "  Gaussian functions.\n";
  cout << "\n";
  cout << "  Number of component Gaussians = " << N << "\n";
  cout << "\n";
  cout << "          Center    Amplitude  Width\n";
  cout << "        X       Y\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                                   << "  "
         << setw(2) << j                   << "  "
         << setw(6) << center[0+j*DIM_NUM] << "  "
         << setw(6) << center[1+j*DIM_NUM] << "  "
         << setw(6) << amplitude[j]        << "  "
         << setw(6) << width[j]            << "\n";
  }

  cout << "\n";
  cout << "      X       Y        Gauss_Sum(X,Y)\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x[0] = ( double ) i;
    for ( j = 0; j <= 10; j++ )
    {
      x[1] = ( double ) j;
      gxy = gauss_sum ( DIM_NUM, N, amplitude, center, width, x );
      cout << "  "
           << setw(6) << x[0] << "  "
           << setw(6) << x[1] << "  "
           << setw(14) << gxy << "\n";
    }
  }

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests GET_SEED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  int seed;
  int seed_0;
  int seed_1;
  int seed_2;
  int seed_3;
  int test;
  double x;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  GET_SEED gets a seed for the random number\n";
  cout << "  generator.  These values are computed from\n";
  cout << "  the time and date.  Values computed nearby\n";
  cout << "  in time will be near to each other, and\n";
  cout << "  should be passed through a random number\n";
  cout << "  generator a few times before use.\n";
  cout << "\n";
  cout << "     I	     R(I)	 R2(I)        R3(I)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    seed = ( int ) get_seed ( );
    seed_0 = seed;
    x = r8_uniform_01 ( &seed );
    seed_1 = seed;
    x = r8_uniform_01 ( &seed );
    seed_2 = seed;
    x = r8_uniform_01 ( &seed );
    seed_3 = seed;

    cout << "  " << setw(12) << seed_0
         << "  " << setw(12) << seed_1
         << "  " << setw(12) << seed_2
         << "  " << setw(12) << seed_3 << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests GRID1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5
# define NSTEP 11

  int i;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  cout << "\n";
  cout << "TEST11\n";
  cout << "  GRID1 computes a 1D grid between\n";
  cout << "  two DIM_NUM dimensional points X1 and X2.\n";
  cout << "\n";
  cout << "  Here, we will use " << NSTEP << " steps\n";
  cout << "  going from: \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << x1[i];
  }
  cout << "\n";
  cout << "  to:\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << x2[i];
  }
  cout << "\n";
  cout << "\n";

  x = grid1 ( DIM_NUM, NSTEP, x1, x2 );

  r8mat_transpose_print ( DIM_NUM, NSTEP, x, "  The grid matrix:" );

  delete [] x;

  return;
# undef DIM_NUM
# undef NSTEP
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests GRID1N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5
# define NSTEP 11

  int i;
  int j;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  cout << "\n";
  cout << "TEST12\n";
  cout << "  GRID1N computes a 1D grid between\n";
  cout << "  two DIM_NUM dimensional points X1 and X2,\n";
  cout << "  one point at a time.\n";
  cout << "\n";
  cout << "  Here, we will use " << NSTEP << " steps\n";
  cout << "  going from \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << x1[i];
  }
  cout << "\n";
  cout << "  to\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << x2[i];
  }
  cout << "\n";
  cout << "\n";
 
  for ( j = 1; j <= NSTEP; j++ )
  {
    x = grid1n ( j, DIM_NUM, NSTEP, x1, x2 );

    cout << "  " << setw(6) << j  << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << x[i];
    }
    cout << "\n";
    delete [] x;
  }
 
  return;
# undef DIM_NUM
# undef NSTEP
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests GRID2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5
# define NSTEP 20

  int i;
  int j1;
  int j2;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  j1 = 3;
  j2 = 13;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  GRID2 computes a 1 D grid between\n";
  cout << "  two DIM_NUM dimensional points X1 and X2,\n";
  cout << "  computing X1 and X2 at user specified times.\n";
  cout << "\n";
  cout << "  Here, we will use " << NSTEP << " steps.\n";
  cout << "  and on step " << j1 << " we will compute\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "  and on step " << j2 << " we will compute\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "\n";
 
  x = grid2 ( j1, j2, DIM_NUM, NSTEP, x1, x2 );
 
  r8mat_print ( DIM_NUM, NSTEP, x, "  The grid matrix:" );

  delete [] x;

  return;
# undef DIM_NUM
# undef NSTEP
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests GRID2N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  j1 = 3;
  j2 = 13;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  GRID2N computes points from a 1D grid\n";
  cout << "  between two DIM_NUM dimensional points\n";
  cout << "  X1 and X2, one at a time, with X1 and X2\n";
  cout << "  having user specified J coordinates.\n";
  cout << "\n";
  cout << "  On step " << j1 << " we will compute\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "  and on step " << j2 << " we will compute\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "\n";
 
  for ( j = 1; j <= 20; j++ )
  {
    x = grid2n ( j, j1, j2, DIM_NUM, x1, x2 );

    cout << "  " << setw(6) << j;

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << x[i];
    }
    cout << "\n";

    delete [] x;
  }
 
  return;
# undef DIM_NUM
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests GRID3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5

  int i;
  int j;
  int k;
  int nstep1 = 3;
  int nstep2 = 6;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };
 
  cout << "\n";
  cout << "TEST15\n";
  cout << "  GRID3 computes a 2D grid in the plane\n";
  cout << "  containing the DIM_NUM-dimensional\n";
  cout << "  points X1, X2 and X3.\n";
  cout << "\n";
  cout << "  Here, we will use " << nstep1 << " steps\n";
  cout << "  going from \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "  to\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "  and " << nstep2 << " steps going to \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x3[i];
  }
  cout << "\n";
 
  x = grid3 ( DIM_NUM, nstep1, nstep2, x1, x2, x3 );
 
  for ( j = 1; j <= nstep1; j++ )
  {
    cout << "\n";
    for ( k = 1; k <= nstep2; k++ )
    { 
      cout << "  " << setw(3) << j
           << "  " << setw(3) << k;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(10) << x[i+(j-1)*DIM_NUM+(j-1)*(k-1)*DIM_NUM*nstep1];
      }
      cout << "\n";
    }
  }

  delete [] x;
 
  return;
# undef DIM_NUM
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests GRID3N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5

  int i;
  int j;
  int k;
  int nstep1 = 3;
  int nstep2 = 6;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  cout << "\n";
  cout << "TEST16\n";
  cout << "  GRID3N computes a point from a 2D\n";
  cout << "  grid in the plane containing the \n";
  cout << "  DIM_NUM-dimensional points X1, X2 and X3.\n";
  cout << "\n";
  cout << "  We use " << nstep1 << " steps from \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "  to\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "  and " << nstep2 << " steps going to \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x3[i];
  }
  cout << "\n";
 
  for ( j = 1; j <= nstep1; j++ )
  {
    cout << "\n";
    for ( k = 1; k <= nstep2; k++ )
    {
      x = grid3n ( j, k, DIM_NUM, nstep1, nstep2, x1, x2, x3 );
      cout << "  " << setw(3) << j
           << "  " << setw(3) << k;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(10) << x[i];
      }
      cout << "\n";
      delete [] x;
    }
  }
 
  return;
#  undef DIM_NUM
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests GRID4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  int k;
  int k1;
  int k2;
  int nstep1 = 6;
  int nstep2 = 10;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  j1 = 2;
  j2 = 5;
  k1 = 3;
  k2 = 9;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  GRID4 computes a 2D planar grid\n";
  cout << "  containing the DIM_NUM-dimensional\n";
  cout << "  points X1, X2 and X3.\n";
  cout << "\n";
  cout << "  We compute the points on the following steps:\n";
  cout << "\n";
  cout << "  X1 on step " << j1 << "  " << k1 << "\n";
  cout << "  X2 on step " << j2 << "  " << k1 << "\n";
  cout << "  X3 on step " << j1 << "  " << k2 << "\n";
  cout << "\n";
  cout << "  We use " << nstep1 << " steps in the J direction\n";
  cout << "  and " << nstep2 << " steps in the K direction.\n";
  cout << "\n";
  cout << "  The points X1, X2 and X3 are:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x3[i];
  }
  cout << "\n";
 
  x = grid4 ( j1, j2, k1, k2, DIM_NUM, nstep1, nstep2, x1, x2, x3 );
 
  for ( j = 1; j <= nstep1; j++ )
  {
    cout << "\n";
    for ( k = 1; k <= nstep2; k++ )
    { 
      cout << "  " << setw(3) << j
           << "  " << setw(3) << k;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(10) << x[i+j*DIM_NUM+k*DIM_NUM*nstep1];
      }
      cout << "\n"; 
    }
  }

  delete [] x;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests GRID4N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  int k;
  int k1;
  int k2;
  int nstep1 = 6;
  int nstep2 = 10;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  j1 = 2;
  j2 = 5;
  k1 = 3;
  k2 = 9;
  
  cout << "\n";
  cout << "TEST18\n";
  cout << "  GRID4N computes, one at a time, points\n";
  cout << "  on a 2D grid in the plane containing\n";
  cout << "  the DIM_NUM-dimensional points X1, X2 and X3.\n";
  cout << "\n";
  cout << "  We compute the points on the following steps:\n";
  cout << "\n";
  cout << "  X1 on step " << j1 << "  " << k1 << "\n";
  cout << "  X2 on step " << j2 << "  " << k1 << "\n";
  cout << "  X3 on step " << j1 << "  " << k2 << "\n";
  cout << "\n";
  cout << "  We use " << nstep1 << " steps in the J direction\n";
  cout << "  and " << nstep2 << " steps in the K direction.\n";
  cout << "\n";
  cout << "  The points X1, X2 and X3 are:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x1[i];
  }
  cout << "\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x2[i];
  }
  cout << "\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(10) << x3[i];
  }
  cout << "\n";

  for ( j = 1; j <= nstep1; j++ )
  {
    cout << "\n";
    for ( k = 1; k <= nstep2; k++ )
    {
      x = grid4n ( j, j1, j2, k, k1, k2, DIM_NUM, nstep1, nstep2, 
        x1, x2, x3 );
      cout << "  " << setw(6) << j
           << "  " << setw(6) << k;
      for ( i = 1; i <= DIM_NUM; i++ )
      {
        cout << "  " << setw(10) << x[i-1];
      }
      cout << "\n";
      delete [] x;
    }
  }
 
  return;
# undef DIM_NUM
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests INDEX1_COL, INDEX1_ROW, and related functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_max;
  int i_min;
  int in[4];
  int in_max[4];
  int in_min[4];
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int l;
  int l_max;
  int l_min;
  int n;
  int value;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  INDEX1_COL column indexes a 1D array,\n";
  cout << "  INDEX1_ROW row indexes a 1D array,\n";
  cout << "  and there are several more versions of these functions.\n";

  cout << "\n";
  cout << "  By COLS:\n";
  cout << "\n";
  cout << "  Imin     I  Imax  Xmin Index\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  index_min = 0;

  value = index1_col ( i_min, i, i_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "        INDEX1_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 1;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  index_min = 0;
  value = index2_col ( i_min, i, i_max, j_min, j, j_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "        INDEX2_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 2;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  index_min = 0;
  value = index3_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
   index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "  " << setw(4) << k_min << setw(4) << k << setw(4) << k_max << "\n";
  cout << "        INDEX3_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 3;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;
  index_min = 0;
  value = index4_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    l_min, l, l_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "  " << setw(4) << k_min << setw(4) << k << setw(4) << k_max << "\n";
  cout << "  " << setw(4) << l_min << setw(4) << l << setw(4) << l_max << "\n";
  cout << "        INDEX4_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 4;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  in_min[3] = 1;
  in[3] = 2;
  in_max[3] = 2;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_COL" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  cout << "\n";
  cout << "  By ROWS:\n";
  cout << "\n";
  cout << "  Imin     I  Imax  Xmin Index\n";
  cout << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  index_min = 0;
  value = index1_row ( i_min, i, i_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "        INDEX1_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 1;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  index_min = 0;
  value = index2_row ( i_min, i, i_max, j_min, j, j_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "        INDEX2_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 2;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  index_min = 0;
  value = index3_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "  " << setw(4) << k_min << setw(4) << k << setw(4) << k_max << "\n";
  cout << "        INDEX3_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 3;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;
  index_min = 0;
  value = index4_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    l_min, l, l_max, index_min );
  cout << "\n";
  cout << "  " << setw(4) << i_min << setw(4) << i << setw(4) << i_max << "\n";
  cout << "  " << setw(4) << j_min << setw(4) << j << setw(4) << j_max << "\n";
  cout << "  " << setw(4) << k_min << setw(4) << k << setw(4) << k_max << "\n";
  cout << "  " << setw(4) << l_min << setw(4) << l << setw(4) << l_max << "\n";
  cout << "        INDEX4_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  n = 4;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  in_min[3] = 1;
  in[3] = 2;
  in_max[3] = 2;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  cout << "        INDEXN_ROW" << "  " << setw(4) << index_min << "  " << setw(4) << value << "\n";

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests ISBN_CHECK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 8

  int check;
  char isbn[14];
  int test;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  ISBN_CHECK checks ISBN's.\n";
  cout << "\n";
  cout << "  A correct ISBN has a checksum of 0.\n";
  cout << "\n";
//
//  Sorry, but until I figure out a decent way to set up and
//  access an array of strings we'll have to do this the
//  bonehead way.
//
  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      strcpy ( isbn, "0-8493-9640-9" );
    }
    else if ( test == 1 ) 
    {
      strcpy ( isbn, "0-201-54275-7" );
    }
    else if ( test == 2 ) 
    {
      strcpy ( isbn, "0-521-35796-9" );
    }
    else if ( test == 3 ) 
    {
      strcpy ( isbn, "0-07-034025-0" );
    }
    else if ( test == 4 ) 
    {
      strcpy ( isbn, "0-7493-9640-9" );
    }
    else if ( test == 5 ) 
    {
      strcpy ( isbn, "0-201-54275-X" );
    }
    else if ( test == 6 ) 
    {
      strcpy ( isbn, "0-521-X5796-9" );
    }
    else if ( test == 7 ) 
    {
      strcpy ( isbn, "0-37-034025-0" );
    }
    check = isbn_check ( isbn );

    cout << "  " << isbn << "  " << check << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests ISBN_FILL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  int check;
  char isbn[20];
  char isbn_test[20];    
  int test;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  ISBN_FILL can fill in a single missing digit\n";
  cout << "  in an ISBN.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      strcpy ( isbn_test, "0-?493-9640-9" );
    }
    else if ( test == 1 )
    {
      strcpy ( isbn_test, "0-201-5427?-7" );
    }
    else if ( test == 2 )
    {
      strcpy ( isbn_test, "0-521-35796-?" );
    }
    else if ( test == 3 )
    {
      strcpy ( isbn_test, "?-07-034025-0" );
    }
    else if ( test == 4 )
    {
      strcpy ( isbn_test, "0-07-05?489-2" );
    }
    strcpy ( isbn, isbn_test );

    isbn_fill ( isbn );

    check = isbn_check ( isbn );

    cout << "  " << isbn_test
         << "  " << isbn
         << "  " << check << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests LCM_12N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  LCM_12N computes the least common multiple of the\n";
  cout << "  integers 1 through N.\n";
  cout << "\n";
  cout << "  N     LCM_12N ( N )\n";
  cout << "\n";
  for ( n = 1; n <= 12; n++ )
  {
    cout << "  " << setw(3) << n
         << "  " << lcm_12n ( n ) << "\n";
  }

  return;
}
//****************************************************************************80

void test225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST225 tests L4MAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  bool *a;
  int i;
  int j;
  int m = 20;
  int n = 50;

  a = new bool[m*n];

  cout << "\n";
  cout << "TEST225\n";
  cout << "  L4MAT_PRINT prints a logical matrix.\n";

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*m] = ( ( i + 1 ) % ( j + 1 ) == 0 );
    }
  }

  l4mat_print ( m, n, a, "  A(I,J) = I+1 is divisible by J+1" );

  delete [] a;

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests LUHN_CHECK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  int check_sum;
  int check_sum_test[TEST_NUM] = {
     6, 
    20, 
    40, 
    80 };
  int *digit;
  int digit_num;
  int digit_num_test[TEST_NUM] = {
     4, 
     4, 
     9, 
    15 };
  int digit_test[32] = {
    1, 1, 1, 1, 
    8, 7, 6, 3, 
    4, 4, 6, 6, 6, 7, 6, 5, 1, 
    3, 7, 7, 9, 5, 6, 5, 7, 0, 9, 4, 4, 7, 2, 6 };
  int i;
  int ilo;
  int test;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  LUHN_CHECK computes the Luhn checksum\n";
  cout << "  for a string of digits.\n";
  cout << "\n";
  cout << "  A correct string has a checksum divisible by 10.\n";

  ilo = 0;
  for ( test = 0; test < TEST_NUM; test++ )
  {
    digit_num = digit_num_test[test];

    digit = new int[digit_num];
    for ( i = 0; i < digit_num; i++ )
    {
      digit[i] = digit_test[i+ilo];
    }
    ilo = ilo + digit_num;

    check_sum = luhn_check ( digit_num, digit );

    cout << "\n";
    cout << "  Test number " << test << "\n";
    cout << "  Number of digits = " << digit_num << "\n";
    cout << "  Digits = ";
    for ( i = 0; i < digit_num; i++ )
    {
      cout << digit[i];
    }
    cout << "\n";
    cout << "  Computed check sum = " << check_sum << "\n";
    cout << "  Correct check sum =  " << check_sum_test[test] << "\n";

    delete [] digit;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests PERM_INVERSE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int p[N] = { 4, 3, 5, 1, 7, 6, 2 };

  cout << "\n";
  cout << "TEST24\n";
  cout << "  PERM_INVERSE inverts a permutation in place;\n";
  cout << "\n";

  perm_print ( N, p, "  The original permutation:" );

  perm_inverse ( N, p );

  perm_print ( N, p, "  The inverted permutation:" );

  return;
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests PRIME_GE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int p;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  PRIME_GE returns the smallest prime number greater\n";
  cout << "  than or equal to N.\n";
  cout << "\n";
  cout << "  N   PRIME_GE\n";
  cout << "\n";
  for ( n = 1; n <= 10; n++ )
  {
    p = prime_ge ( n );
    cout << "  " << setw ( 6 ) << n 
         << "  " << setw ( 6 ) << p << "\n";
  }

  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests RANDOM_INITIALIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double r1;
  double r2;
  double r3;
  unsigned long seed;
  unsigned long seed_in;
  unsigned long seed_out;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  RANDOM_INITIALIZE can make up a seed for the C\n";
  cout << "  random number generator RANDOM, or use a\n";
  cout << "  single SEED value from the user.\n";
  cout << "\n";
  cout << "  Calling RANDOM_INITIALIZE with a zero input value of SEED\n";
  cout << "  tells the routine to make up a seed.  And, at least for\n";
  cout << "  calls a few milliseconds apart, the output SEED should\n";
  cout << "  be different.\n";
  cout << "\n";
  cout << "  In any case, if RANDOM is restarted by calling\n";
  cout << "  RANDOM_INITIALIZE with a nonzero input SEED, then\n";
  cout << "  the random number sequence should repeat.\n";
  cout << "\n";
  cout << "  Call RANDOM_INITIALIZE 10 times, with a zero input SEED.\n";
  cout << "  Also, get the first three real random values.\n";
  cout << "\n";
  cout << "    SEED_IN         SEED_OUT     Random 1, 2, 3\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = 0;
    seed = seed_in;
    seed_out = random_initialize ( seed );

    r1 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r2 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r3 = ( double ) random ( ) / ( double ) ( RAND_MAX );

    cout << "  " << setw(12) << seed_in
         << "  " << setw(12) << seed_out
         << "  " << setw(12) << r1
         << "  " << setw(12) << r2
         << "  " << setw(12) << r3 << "\n";
  }

  cout << "\n";
  cout << "  Now call RANDOM_INITIALIZE with SEED = 5, 95, 5, 95.\n";
  cout << "  We promise the random numbers will repeat the second time.\n";
  cout << "\n";
  cout << "    SEED_IN         SEED_OUT     Random 1, 2, 3\n";
  cout << "\n";

  seed_in = 5;

  for ( i = 1; i <= 4; i++ )
  {
    seed = seed_in;
    seed_out = random_initialize ( seed );

    r1 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r2 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r3 = ( double ) random ( ) / ( double ) ( RAND_MAX );

    cout << "  " << setw(12) << seed_in
         << "  " << setw(12) << seed_out
         << "  " << setw(12) << r1
         << "  " << setw(12) << r2
         << "  " << setw(12) << r3 << "\n";

    seed_in = 100 - seed_in;
  }

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 1000

  int i;
  int n;
  int seed;
  double x[N_MAX];
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  RANDOM is an intrinsic C routine\n";
  cout << "  to computer uniform random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
//
//  Set the random number seed.
//
  seed = 123456789;
  srandom ( seed );
//
//  Test 1:
//  Simply call 5 times for 1 value, and print.
//
  cout << "\n";
  cout << "  Test #1: Call 5 times, 1 value each time.\n";
  cout << "\n";

  n = 1;
  for ( i = 1; i <= 5; i++ )
  {
    x[0] = ( double ) random ( ) / ( double ) ( RAND_MAX );
    cout << "  " << setw(6) << i 
         << "  " << setw(12) << x[0] << "\n";
  }
//
//  Test 2:
//  Restore the random number seed, and repeat.
//
  cout << "\n";
  cout << "  Test #2: Restore the random number seed.\n";
  cout << "  Call 5 times, 1 value each time.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  seed = 123456789;
  srandom ( seed );

  n = 1;
  for ( i = 1; i <= 5; i++ )
  {
    x[0] = ( double ) random ( ) / ( double ) ( RAND_MAX );
    cout << "  " << setw(6) << i 
         << "  " << setw(12) << x[0] << "\n";
  }
//
//  Test 5:
//  Determine the minimum, maximum, mean and variance.
//
  n = N_MAX;
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) random ( ) / ( double ) ( RAND_MAX );
  }
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );

  cout << "\n";
  cout << "  Test #5:\n";
  cout << "  Number of samples was " << n << "\n";
  cout << "  Minimum value was " << x_min << "\n";
  cout << "  Maximum value was " << x_max << "\n";
  cout << "  Average value was " << x_mean << "\n";
  cout << "  Variance was      " << x_var << "\n";
  cout << "  Expected average  " << 0.5 << "\n";
  cout << "  Expected variance " << 1.0 / 12.0 << "\n";

  return;
# undef N_MAX
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests RAT_FACTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define FACTOR_MAX 10

  int factor[FACTOR_MAX];
  int factor_num;
  int i;
  int m;
  int mleft;
  int n;
  int nleft;
  int power[FACTOR_MAX];

  cout << "\n";
  cout << "TEST29\n";
  cout << "  RAT_FACTOR factors a rational value.\n";

  m = 13 * 7 * 9 * 2;
  n = 12;

  cout << "\n";
  cout << "  Rational value is " << m << " / " << n << "\n";

  rat_factor ( m, n, FACTOR_MAX, &factor_num, factor, power, &mleft, &nleft );

  cout << "\n";
  cout << "  Prime representation:\n";
  cout << "\n";
  cout << "  I, FACTOR(I), POWER(I)\n";
  cout << "\n";

  if ( mleft != 1 || nleft != 1 )
  {
    cout << "  " << setw(6) << 0
         << "  " << setw(6) << mleft << " / "
         << "  " << setw(6) << nleft << " (UNFACTORED PORTION)\n";
    cout << "\n";
  }

  for ( i = 0; i < factor_num; i++ )
  {
    cout << "  " << setw(6) << i + 1
         << "  " << setw(6) << factor[i] << "   "
         << "  " << setw(6) << power[i] << "\n";
  }
 
  return;
# undef FACTOR_MAX
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests ROOTS_TO_R8POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *c;
  double *x;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  ROOTS_TO_R8POLY computes the coefficients of\n";
  cout << "  a polynomial from its roots.\n";
  cout << "  R8POLY_PRINT prints a polynomial.\n";

  x = r8vec_indicator_new ( N );

  r8vec_print ( N, x, "  Roots:" );

  c = roots_to_r8poly ( N, x );

  r8poly_print ( N, c, "  The polynomial" );

  delete [] c;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests SORT_HEAP_EXTERNAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int a[N];
  int b;
  int c;
  int i;
  int indx;
  int isgn;
  int itemp;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  SORT_HEAP_EXTERNAL sorts objects externally.\n";

  indx = 0;
  i = 0;
  j = 0;
  isgn = 0;

  b = 1;
  c = N;
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i] = i4_uniform ( b, c, &seed );
  }
 
  i4vec_print ( N, a, "  Unsorted array:" );
//
//  Call the sort routine over and over.
//
  for ( ;; )
  {
    sort_heap_external ( N, &indx, &i, &j, isgn );
//
//  If the return value of INDX is negative, we're asked to compare
//  array elements I and J;
//
    if ( indx < 0 )
    {
      if ( a[i] <= a[j] )
      {
        isgn = -1;
      }
      else
      {
        isgn = 1;
      }

    }
//
//  ...and if the return value of INDX is positive, we're asked to switch
//  array elements I and J;
//
    else if ( 0 < indx )
    {
      i4_swap ( &a[i], &a[j] );
//
//  ...and if the return value of INDX is 0, we're done.
//
    } 
    else
    {
      break;
    }

  }

  i4vec_print ( N, a, "  Sorted array:" );
 
  return;

# undef N
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int nt;
  double *t;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  For evenly spaced angles between 0 and 2*PI:\n";
  cout << "  TVEC_EVEN\n";
  cout << "  TVEC_EVEN2\n";
  cout << "  TVEC_EVEN3\n";

  nt = 4;

  t = tvec_even ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN" );

  delete [] t;

  nt = 4;

  t = tvec_even2 ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN2" );

  delete [] t;

  nt = 4;

  t = tvec_even3 ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN3" );

  delete [] t;

  return;
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2, TVEC_EVEN_BRACKET3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nt;
  double *t;
  double theta1;
  double theta2;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  For evenly spaced angles between THETA1 and THETA2:\n";
  cout << "  TVEC_EVEN_BRACKET\n";
  cout << "  TVEC_EVEN_BRACKET2.\n";
  cout << "  TVEC_EVEN_BRACKET3.\n";

  nt = 4;
  theta1 = 30.0;
  theta2 = 90.0;

  t = tvec_even_bracket ( nt, theta1, theta2 );

  cout << "\n";
  cout << "    NT = " << nt << "\n";
  cout << "    THETA1 = " << theta1 << "\n";
  cout << "    THETA2 = " << theta2 << "\n";

  r8vec_print ( nt, t, "  TVEC_BRACKET" );

  nt = 5;

  t = tvec_even_bracket2 ( nt, theta1, theta2 );

  cout << "\n";
  cout << "    NT = " << nt << "\n";
  cout << "    THETA1 = " << theta1 << "\n";
  cout << "    THETA2 = " << theta2 << "\n";

  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET2" );

  nt = 3;

  t = tvec_even_bracket3 ( nt, theta1, theta2 );

  cout << "\n";
  cout << "    NT = " << nt << "\n";
  cout << "    THETA1 = " << theta1 << "\n";
  cout << "    THETA2 = " << theta2 << "\n";

  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET3" );

  return;
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests UPC_CHECK_DIGIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2

  int c;
  int l;
  int l_test[TEST_NUM] = { 72890, 12345 };
  int p;
  int p_test[TEST_NUM] = { 0, 0 };
  int r;
  int r_test[TEST_NUM] = { 11, 67890 };
  int test;

  cout << "\n";
  cout << "TEST34\n";
  cout << "  UPC_CHECK_DIGIT determines the check digit for a UPC.\n";
  cout << "\n";
  cout << "  P-LLLLL-RRRRR-C\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];
    l = l_test[test];
    r = r_test[test];

    c = upc_check_digit ( p, l, r );

    cout << "  " << setw(1) << p 
         << "-"  << setfill ( '0' ) << setw(5) << l
         << "-"  << setfill ( '0' ) << setw(5) << r
         << "-"  << setw(1) << c << "\n";
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 calls VERSINE_PULSE.
//
//  Modified:
//
//    20 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double amp;
  int i;
  double t;
  double ta;
  double tb;
  double v;
  double v1;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  VERSINE_PULSE adds a versine pulse to a constant signal.\n";
  cout << "\n";

  ta = 2.0;
  tb = 4.0;
  v1 = 1.0;
  amp = 3.0;

  for ( i = 0; i <= 100; i++ )
  {
    t = ( double ) i / 10.0;
    v = versine_pulse ( t, ta, tb, v1, amp );
    cout << "  " << setw(4) << i
         << "  " << setw(10) << t
         << "  " << setw(10) << v << "\n";
  }
  return;
}
