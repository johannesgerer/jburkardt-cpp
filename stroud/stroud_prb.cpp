# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "stroud.hpp"
//
//  GLOBAL VARIABLES USED TO SELECT INTEGRATION FUNCTION.
//
int function_1d_index = 0;
int function_2d_index = 0;
int function_3d_index = 0;
int function_nd_index = 0;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test045 ( );
void test05 ( );
void test052 ( );
void test054 ( );
void test07 ( );
void test08 ( );
void test085 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test163 ( );
void cn_geg_test ( int n, double alpha, int expon[] );
void test165 ( );
void cn_jac_test ( int n, double alpha, double beta, int expon[] );
void test167 ( );
void cn_leg_test ( int n, int expon[] );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test205 ( );
void test207 ( );
void en_r2_test ( int n, int expon[] );
void test2075 ( );
void epn_glg_test ( int n, int expon[], double alpha );
void test208 ( );
void epn_lag_test ( int n, int expon[] );
void test21 ( );
void test215 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test255 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test322 ( );
void test324 ( );
void test326 ( );
void test33 ( );
void test335 ( );
void test34 ( );
void test345 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );
void test40 ( );
void test41 ( );
void test42 ( );
void test425 ( );
void test43 ( );
void test44 ( );
void test45 ( );
void test46 ( );
void test47 ( );
void test48 ( );
void test49 ( );
double fu18 ( double x );
double fl18 ( double x );
double fu28 ( double x, double y );
double fl28 ( double x, double y );

double function_1d ( double x );
void function_1d_name ( char *name );
int function_1d_num ( );

double function_2d ( double x, double y );
void function_2d_name ( char *name );
int function_2d_num ( );

double function_3d ( double x, double y, double z );
void function_3d_name ( char *name );
int function_3d_num ( );

double function_nd ( int n, double x[] );
void function_nd_name ( char *name );
int function_nd_num ( );

double f_1_2d ( double x, double y );
double f_x_2d ( double x, double y );
double f_r_2d ( double x, double y );
double mono_000_3d ( int n, double x[] );
double mono_111_3d ( int n, double x[] );
double mono_202_3d ( int n, double x[] );
double mono_422_3d ( int n, double x[] );
double *setsim ( int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for STROUD_PRB.
//
//  Discussion:
//
//    STROUD_PRB calls the STROUD tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "STROUD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the STROUD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test045 ( );
  test05 ( );
  test052 ( );
  test054 ( );
  test07 ( );
  test08 ( );
  test085 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test163 ( );
  test165 ( );
  test167 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test205 ( );
  test207 ( );
  test2075 ( );
  test208 ( );
  test21 ( );
  test215 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test255 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test322 ( );
  test324 ( );
  test326 ( );
  test33 ( );
  test335 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
  test425 ( );
  test43 ( );
  test44 ( );
  test45 ( );
  test46 ( );
  test47 ( );
  test48 ( );
  test49 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "STROUD_PRB\n";
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
//    TEST01 tests BALL_F1_ND, BALL_F3_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *center;
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result2;
  double r;

  r = 2.0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For integrals in a ball in ND:\n";
  cout << "  BALL_F1_ND approximates the integral;\n";
  cout << "  BALL_F3_ND approximates the integral.\n";

  for ( n = 2; n <= n_max; n++ )
  {
    center = new double[n];
    for ( i = 0; i < n; i++ )
    {
      center[i] = ( double ) ( i + 1 );
    }
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
    cout << "  Ball center:\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(10) << center[i];
    }
    cout << "\n";
    cout << "  Ball radius = " << r << "\n";
    cout << "  Ball volume = " << ball_volume_nd ( n, r ) << "\n";
    cout << "\n";
    cout << "    Rule:      F1          F3\n";
    cout << "    F(X)\n";
    cout << "\n";

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = ball_f1_nd ( function_nd, n, center, r );
      result2 = ball_f3_nd ( function_nd, n, center, r );
      cout << "  " << name
           << "  " << setw(14) << result1
           << "  " << setw(14) << result2 <<"\n";
    }
    delete [] center;
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BALL_MONOMIAL_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  int p[3];
  double result1;
  double result2;
  double r = 2.0;
  char string[11];
  int test;
  int test_num = 4;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For the integral of a monomial in a ball in ND:\n";
  cout << "  BALL_MONOMIAL_ND approximates the integral.\n";
  cout << "  BALL_F1_ND, which can handle general integrands,\n";
  cout << "    will be used for comparison.\n";
  cout << "\n";
  cout << "  Spatial dimension N = " << dim_num << "\n";
  cout << "  Ball radius = " << r << "\n";
  cout << "  Ball volume = " << ball_volume_nd ( dim_num, r ) << "\n";
  cout << "\n";
  cout << "    Rule:     MONOMIAL    F1\n";
  cout << "    F(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      strcpy ( string, "         1" );
      p[0] = 0;
      p[1] = 0;
      p[2] = 0;
      result2 = ball_f1_nd ( mono_000_3d, dim_num, center, r );
    }
    else if ( test == 2 )
    {
      strcpy ( string, "       xyz" );
      p[0] = 1;
      p[1] = 1;
      p[2] = 1;
      result2 = ball_f1_nd ( mono_111_3d, dim_num, center, r );
    }
    else if ( test == 3 )
    {
      strcpy ( string, "   x^2 z^2" );
      p[0] = 2;
      p[1] = 0;
      p[2] = 2;
      result2 = ball_f1_nd ( mono_202_3d, dim_num, center, r );
    }
    else if ( test == 4 )
    {
      strcpy ( string, " x^4y^2z^2" );
      p[0] = 4;
      p[1] = 2;
      p[2] = 2;
      result2 = ball_f1_nd ( mono_422_3d, dim_num, center, r );
    }

    result1 = ball_monomial_nd ( dim_num, p, r );

    cout << "  " << string
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2 << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BALL_UNIT_**_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For integrals in the unit ball in 3D:\n";
  cout << "  BALL_UNIT_07_3D uses a formula of degree 7;\n";
  cout << "  BALL_UNIT_14_3D uses a formula of degree 14;\n";
  cout << "  BALL_UNIT_15_3D uses a formula of degree 15.\n";
  cout << "\n";
  cout << "  Unit ball volume = " << ball_unit_volume_nd ( 3 ) << "\n";
  cout << "\n";
  cout << "    Rule:      #7             #14           #15\n";
  cout << "    F(X)\n";
  cout << "\n";
 
  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = ball_unit_07_3d ( function_3d );
    result2 = ball_unit_14_3d ( function_3d );
    result3 = ball_unit_15_3d ( function_3d );

    cout << "  " << name
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2
         << "  " << setw(14) << result3 << "\n";
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BALL_UNIT_F1_ND, BALL_UNIT_F3_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result2;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For integrals inside the unit ball in ND:\n";
  cout << "  BALL_UNIT_F1_ND approximates the integral;\n";
  cout << "  BALL_UNIT_F3_ND approximates the integral.\n";
  cout << "\n";
 
  for ( n = 2; n <= n_max; n++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
    cout << "  Unit ball volume = " << ball_unit_volume_nd ( n ) << "\n";
    cout << "\n";
    cout << "\n";
    cout << "    Rule:      F1          F3\n";
    cout << "    F(X)\n";
    cout << "\n";
 
    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = ball_unit_f1_nd ( function_nd, n );
      result2 = ball_unit_f3_nd ( function_nd, n );

      cout << "  " << name  
           << setw(14) << result1
           << setw(14) << result2 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests BALL_UNIT_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num = 3;

  cout << "\n";
  cout << "TEST045\n";
  cout << "  In 3 dimensions:\n";
  cout << "  BALL_UNIT_VOLUME_3D gets the volume of the unit ball.\n";
  cout << "  BALL_UNIT_VOLUME_ND will be called for comparison.\n";
  cout << "\n";
  cout << "    N    Volume    Method\n";
  cout << "\n";

  cout << "  " << setw(3) << dim_num 
       << "  " << setw(14) << ball_unit_volume_3d (   ) 
       << "  BALL_UNIT_VOLUME_3D\n";

  cout << "  " << setw(3) << dim_num 
       << "  " << setw(14) << ball_unit_volume_nd ( dim_num  ) 
       << "  BALL_UNIT_VOLUME_ND\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests BALL_UNIT_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  BALL_UNIT_VOLUME_ND computes the volume\n";
  cout << "    of the unit ball in ND.\n";
  cout << "\n";
  cout << "    N    Volume\n";
  cout << "\n";

  for ( dim_num = 2; dim_num <= 10; dim_num++ )
  {
    cout << "  " << setw(3) << dim_num
         << "  " <<  ball_unit_volume_nd ( dim_num ) << "\n";
  }
  return;
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests BALL_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 3;
  double r;

  cout << "\n";
  cout << "TEST052\n";
  cout << "  In 3 dimensions:\n";
  cout << "  BALL_VOLUME_3D computes the volume of a unit ball.\n";
  cout << "  BALL_VOLUME_ND will be called for comparison.\n";
  cout << "\n";
  cout << "    N    R      Volume    Method\n";
  cout << "\n";

  r = 1.0;

  for ( i = 1; i <= 3; i++ )
  {
    cout << "  " << setw(3) << n
         << "  " << setw(14) << r
         << "  " << ball_volume_3d ( r ) 
         << "  BALL_VOLUME_3D\n";

    cout << "  " << setw(3) << n
         << "  " << setw(14) << r
         << "  " << ball_volume_nd ( n, r )
         << "  " << "BALL_VOLUME_ND\n";

    r = r * 2.0;
  }
  return;
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests BALL_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double r;

  cout << "\n";
  cout << "TEST054\n";
  cout << "  BALL_UNIT_VOLUME_ND computes the volume of\n";
  cout << "    the unit ball in N dimensions.\n";
  cout << "\n";
  cout << "    N        R      Volume\n";
  cout << "\n";

  for ( n = 2; n <= 10; n++ )
  {
    r = 0.5;
    for ( i = 1; i <= 3; i++ )
    {
      cout << "  " << setw(3) << n
           << "  " << setw(14) << r
           << "  " << setw(14) <<  ball_volume_nd ( n, r ) << "\n";
      r = r * 2.0;
    }
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests CIRCLE_ANNULUS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double center[2];
  double center_test[2*2] = {
    0.0, 0.0, 
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int nr;
  double radius1;
  double radius1_test[2] = { 0.0, 1.0 };
  double radius2;
  double radius2_test[2] = { 1.0, 2.0 };
  double result;
  int test_num = 2;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  CIRCLE_ANNULUS estimates integrals\n";
  cout << "    in a circular annulus.\n";
  cout << "\n";
  cout << "        F       CENTER         Radius1   Radius2   NR  Result\n";
  cout << "\n";

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    cout << "\n";
    cout << "  " << "   Area"
                 << setw(10) << center[0]
                 << setw(10) << center[1]
                 << setw(10) << radius1
                 << setw(10) << radius2
         << "  " << setw(10) << area << "\n";

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;

      for ( nr = 1; nr <= 4; nr++ )
      {
        result = circle_annulus ( function_2d, center, radius1, radius2, nr );

        function_2d_name ( name );

        cout << "  " << name
                     << setw(10) << center[0]
                     << setw(10) << center[1]
                     << setw(10) << radius1
                     << setw(10) << radius2
                     << setw(2)  << nr
                     << setw(10) << result << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CIRCLE_ANNULUS, CIRCLE_RT_SET, CIRCLE_RT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double center[2];
  double center_test[2*3] = {
    0.0, 0.0, 
    0.0, 0.0, 
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int nc;
  int num;
  int nr;
  int nr2;
  int nt;
  double ra[5];
  double radius1;
  double radius1_test[3] = { 0.0, 1.0, 1.0 };
  double radius2;
  double radius2_test[3] = { 1.0, 2.0, 3.0 };
  double result1;
  double result2;
  double result3;
  double rw[5];
  int rule;
  double ta[20];
  int test_num = 3;
  double tw[20];
  double zw;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  CIRCLE_ANNULUS estimates integrals in a\n";
  cout << "    circular annulus.\n";
  cout << "  CIRCLE_RT_SET sets up a rule for a circle;\n";
  cout << "  CIRCLE_RT_SUM applies the rule.\n";
  cout << "\n";
  cout << "  RESULT1 = CIRCLE_ANNULUS result.\n";
  cout << "  RESULT2 = Difference of two CIRCLE_RT_SUM results.\n";
  cout << "\n";
  cout << "        F      CENTER       Radius1   Radius2   Result1 Result2\n";
  cout << "\n";

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    strcpy ( name, "   Area" );
    cout << "\n";
    cout << "  " << name
         << setw(11) << center[0]
         << setw(11) << center[1]
         << setw(11) << radius1
         << setw(11) << radius2
         << setw(11) << area << "\n";

    rule = 9;
    circle_rt_size ( rule, &nr2, &nt, &nc );
    circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;
      function_2d_name ( name );

      nr = 5;
      result1 = circle_annulus ( function_2d, center, radius1, radius2, nr );

      result2 = circle_rt_sum ( function_2d, center, radius1, nr2, ra, rw, nt, 
        ta, tw, zw );

      result3 = circle_rt_sum ( function_2d, center, radius2, nr2, ra, rw, nt, 
        ta, tw, zw );

      cout << "  " << name
           << setw(11) << center[0]
           << setw(11) << center[1]
           << setw(11) << radius1
           << setw(11) << radius2
           << setw(11) << result1
           << setw(11) << result3 - result2 << "\n";
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
//    TEST085 tests CIRCLE_ANNULUS_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double center[2];
  double center_test[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    3.0, 4.0 };
  int dim;
  int dim_num = 2;
  int i;
  int ntest = 3;
  double radius1;
  double radius1_test[3] = { 0.0, 1.0, 1.0 };
  double radius2;
  double radius2_test[3] = { 1.0, 2.0, 3.0 };

  cout << "\n";
  cout << "TEST085\n";
  cout << "  CIRCLE_ANNULUS_AREA_2D computes the area of a\n";
  cout << "    circular annulus.\n";
  cout << "\n";
  cout << "      CENTER       Radius1   Radius2   Area\n";
  cout << "\n";

  for ( i = 0; i < ntest; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }

    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    cout << "\n";
    cout << "  " << setw(9) << center[0]
         << "  " << setw(9) << center[1]
         << "  " << setw(9) << radius1
         << "  " << setw(9) << radius2
         << "  " << setw(9) << area << "\n";
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests CIRCLE_ANNULUS_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double as1;
  double as2;
  double as3;
  double as4;
  double center[2];
  int dim_num = 2;
  int j;
  int nc;
  char name[8];
  int num;
  int nr;
  int nr2;
  int nt;
  double pi = 3.141592653589793;
  double ra[5];
  double radius;
  double radius1a;
  double radius2a;
  double radius1b;
  double radius2b;
  double radius1c;
  double radius2c;
  double radius1d;
  double radius2d;
  double result1;
  double result2;
  int rule;
  double rw[5];
  double ta[20];
  int test_num = 4;
  double theta1a;
  double theta2a;
  double theta1b;
  double theta2b;
  double theta1c;
  double theta2c;
  double theta1d;
  double theta2d;
  double tw[20];
  double zw;

  nr = 5;

  rule = 9;
  circle_rt_size ( rule, &nr2, &nt, &nc );
  circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

  cout << "\n";
  cout << "TEST09\n";
  cout << "  CIRCLE_ANNULUS_SECTOR estimates an integral in a\n";
  cout << "    circular annulus sector.\n";
  cout << "  CIRCLE_RT_SET sets an integration rule in a circle.\n";
  cout << "  CIRCLE_RT_SUM uses an integration rule in a circle.\n";
  cout << "\n";
  cout << "  To test CIRCLE_ANNULUS_SECTOR, we estimate an integral\n";
  cout << "  over 4 annular sectors that make up the unit circle,\n";
  cout << "  and add to get RESULT1.\n";
  cout << "\n";
  cout << "  We will also estimate the integral over the unit circle\n";
  cout << "  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.\n";
  cout << "\n";
  cout << "  We will then compare RESULT1 and RESULT2.\n";
  cout << "\n";
  cout << "  CIRCLE_ANNULUS_SECTOR computations will use NR = " << nr << "\n";
  cout << "  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule " << rule << "\n";
  cout << "\n";
  cout << "  RESULT1 is the sum of Annulus Sector calculations.\n";
  cout << "  RESULT2 is for CIRCLE_RT_SET/CIRCLE_RT_SUM.\n";
  cout << "\n";

  center[0] = 0.0;
  center[1] = 0.0;
  radius = 1.0;

  radius1a = 0.0;
  radius2a = 0.25;
  theta1a = 0.0;
  theta2a = 0.5 * pi;

  radius1b = 0.0;
  radius2b = 0.25;
  theta1b = 0.5 * pi;
  theta2b = 2.0 * pi;

  radius1c = 0.25;
  radius2c = 1.0;
  theta1c = 0.0;
  theta2c = 0.25 * pi;

  radius1d = 0.25;
  radius2d = 1.0;
  theta1d = 0.25 * pi;
  theta2d = 2.0 * pi;

  cout << "\n";
  cout << "       F  Result1  Result2\n";
  cout << "\n";

  num = function_2d_num ( );

  for ( j = 1; j <= num; j++ )
  {
    function_2d_index = j;

    as1 = circle_annulus_sector ( function_2d, center, radius1a, radius2a, theta1a, 
      theta2a, nr );

    as2 = circle_annulus_sector ( function_2d, center, radius1b, radius2b, theta1b, 
      theta2b, nr );

    as3 = circle_annulus_sector ( function_2d, center, radius1c, radius2c, theta1c, 
      theta2c, nr );

    as4 = circle_annulus_sector ( function_2d, center, radius1d, radius2d, theta1d, 
      theta2d, nr );

    result1 = as1 + as2 + as3 + as4;

    result2 = circle_rt_sum ( function_2d, center, radius, nr2, ra, rw, nt, 
      ta, tw, zw );

    function_2d_name ( name );

    cout << "  " << name
         << setw(14) << result1
         << setw(14) << result2 << "\n";
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests CIRCLE_CUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int order;
  double r = 3.0;
 
  cout << "\n";
  cout << "TEST10\n";
  cout << "  CIRCLE_CUM approximates an integral over a circle.\n";
  cout << "\n";
  cout << "  We use radius R = " << r << "\n";
  cout << "  and center:\n";
  cout << "  CENTER = ( " << center[0] << ", " << center[1] << ").\n";
  cout << "\n";

  cout << "\n";
  cout << "    Order:      2             4              8            16\n";
  cout << "  F(X)\n";
  cout << "\n";
 
  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );

    cout << "  " << name;

    for ( j = 1; j <= 4; j++ )
    {
      order = i4_power ( 2, j );

      cout << setw(14) << circle_cum ( function_2d, center, r, order );
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests LENS_HALF_AREA_2D, CIRCLE_SECTOR_AREA_2D, CIRCLE_TRIANGLE_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area1;
  double area2;
  double area3;
  int i;
  double pi = 3.141592653589793;
  double r;
  double theta1;
  double theta2;

  r = 1.0;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  LENS_HALF_AREA_2D computes the area of a\n";
  cout << "    circular half lens, defined by joining the endpoints\n";
  cout << "    of a circular arc.\n";
  cout << "  CIRCLE_SECTOR_AREA_2D computes the area of a\n";
  cout << "    circular sector, defined by joining the endpoints\n";
  cout << "    of a circular arc to the center.\n";
  cout << "  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a\n";
  cout << "    triangle, defined by joining the endpoints\n";
  cout << "    of a circular arc and the center.\n";
  cout << "\n";
  cout << "      R       Theta1 Theta2        "
       << "Sector       Triangle     Half Lens\n";
  cout << "\n";

  for ( i = 0; i <= 12; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / 12.0;

    area1 = circle_sector_area_2d ( r, theta1, theta2 );

    area2 = circle_triangle_area_2d ( r, theta1, theta2 );

    area3 = lens_half_area_2d ( r, theta1, theta2 );

    cout << "  " << setw(9) << r
         << "  " << setw(9) << theta1
         << "  " << setw(14) << theta2
         << "  " << setw(14) << area1
         << "  " << setw(14) << area2
         << "  " << setw(14) << area3 << "\n";
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests LENS_HALF_AREA_2D, LENS_HALF_H_AREA_2D, LENS_HALF_W_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area1;
  double area2;
  double area3;
  double pi = 3.141592653589793;
  double h;
  int i;
  double r;
  double theta1;
  double theta2;
  double w;

  r = 50.0;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For the area of a circular half lens,\n";
  cout << "  LENS_HALF_AREA_2D uses two angles;\n";
  cout << "  LENS_HALF_H_AREA_2D works from the height;\n";
  cout << "  LENS_HALF_W_AREA_2D works from the width.\n";
  cout << "\n";
  cout << "  The circle has radius R = " << r << "\n";
  cout << "\n";
  cout << "  THETA1 THETA2  H     W  Area(THETA) Area(H)  Area(W)\n";
  cout << "\n";

  for ( i = 0; i <= 12; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / 12.0;
    w = 2.0 * r * sin ( 0.5 * ( theta2 - theta1 ) );
    h = r * ( 1.0 - cos ( 0.5 * ( theta2 - theta1 ) ) );

    area1 = lens_half_area_2d ( r, theta1, theta2 );

    area2 = lens_half_h_area_2d ( r, h );

    area3 = lens_half_w_area_2d ( r, w );

    cout << "  " << setw(6) << theta1
         << "  " << setw(6) << theta2
         << "  " << setw(6) << h
         << "  " << setw(6) << w
         << "  " << setw(10) << area1
         << "  " << setw(10) << area2
         << "  " << setw(10) << area3 << "\n";
  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests CIRCLE_SECTOR, CIRCLE_SECTOR_AREA_2D.
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
{
  double area;
  double center[2];
  double center_test[2*4] = {
    0.0, 0.0, 
    0.0, 0.0, 
    0.0, 0.0, 
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int nr;
  int nrhi = 5;
  int nrlo = 1;
  double pi = 3.141592653589793;
  double radius;
  double radius_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  double result;
  int test_num = 4;
  double theta1;
  double theta1_test[4] = { 0.0, 0.0, 0.0, 0.0 };
  double theta2;
  double theta2_test[4] = { 2.0, 1.0, 0.5, 0.25 };

  cout << "\n";
  cout << "TEST13\n";
  cout << "  CIRCLE_SECTOR_AREA_2D computes the area\n";
  cout << "    of a circular sector.\n";
  cout << "  CIRCLE_SECTOR estimates an integral \n";
  cout << "    in a circular sector.\n";
  cout << "\n";
  cout << "  The user can specify NR, the number of radial values\n";
  cout << "  used to approximated the integral.\n";
  cout << "\n";
  cout << "  In this test, computations will use values of NR\n";
  cout << "  from " << nrlo << "\n";
  cout << "  to   " << nrhi << "\n";
  cout << "\n";

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius = radius_test[i];
    theta1 = theta1_test[i] * pi;
    theta2 = theta2_test[i] * pi;

    area = circle_sector_area_2d ( radius, theta1, theta2 );

    cout << "\n";
    cout << "       CENTER      RADIUS  THETA1  THETA2  Area\n";
    cout << "\n";
    cout << setw(8) << center[0]
         << setw(8) << center[1]
         << setw(8) << radius
         << setw(8) << theta1
         << setw(8) << theta2
         << setw(8) << area << "\n";
    cout << "\n";
    cout << "       F   ";
    for ( nr = nrlo; nr <= nrhi; nr++ )
    {
      cout << "      " << setw(2) << nr << "      ";
    }
    cout << "\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;

      function_2d_name ( name );

      cout << "  " << name;

      for ( nr = nrlo; nr <= nrhi; nr++ )
      {
        result = circle_sector ( function_2d, center, radius, theta1, theta2, 
          nr );
        cout << setw(14) << result;
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests CIRCLE_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area1;
  double area2;
  double area3;
  double center[2];
  double center_test[2*4] = {
    0.0, 0.0, 
    0.0, 0.0, 
    0.0, 0.0, 
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int nc;
  int num;
  int nr;
  int nr2;
  int nt;
  double pi = 3.141592653589793;
  double ra[5];
  double radius;
  double radius_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  double result1;
  double result2;
  double resulta;
  double resultb;
  int rule;
  double rw[5];
  double ta[20];
  int test_num = 4;
  double theta1;
  double theta1_test[4] = { 0.0, 0.0, 0.0, 0.0 };
  double theta2;
  double theta2_test[4] = { 2.0, 1.0, 0.5, 0.25 };
  double theta3;
  double tw[20];
  double zw;

  nr = 5;

  rule = 9;
  circle_rt_size ( rule, &nr2, &nt, &nc );
  circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

  cout << "\n";
  cout << "TEST14\n";
  cout << "  CIRCLE_SECTOR estimates integrals in a circular sector.\n";
  cout << "  CIRCLE_RT_SET sets an integration rule in a circle.\n";
  cout << "  CIRCLE_RT_SUM uses an integration rule in a circle.\n";
  cout << "\n";
  cout << "  To test CIRCLE_SECTOR, we estimate an integral over\n";
  cout << "  a sector, and over its complement and add the results\n";
  cout << "  to get RESULT1.\n";
  cout << "\n";
  cout << "  We also estimate the integral over the whole circle\n";
  cout << "  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.\n";
  cout << "\n";
  cout << "  We will then compare RESULT1 and RESULT2.\n";
  cout << "\n";
  cout << "  CIRCLE_SECTOR computations will use NR = " << nr << "\n";
  cout << "  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule " << rule << "\n";
  cout << "\n";
  cout << "  'Sector1' and 'Sector2' are the CIRCLE_SECTOR\n";
  cout << "  computations\n";
  cout << "  for the sector and its complement.\n";
  cout << "  'Sum' is the sum of Sector1 and Sector2.\n";
  cout << "  'Circle' is the computation for \n";
  cout << "  CIRCLE_RT_SET + CIRCLE_RT_SUM.\n";
  cout << "\n";

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius = radius_test[i];

    theta1 = theta1_test[i] * pi;
    theta2 = theta2_test[i] * pi;
    theta3 = theta2 + 2.0 * pi - ( theta2 - theta1 );

    area1 = circle_sector_area_2d ( radius, theta1, theta2 );
    area2 = circle_sector_area_2d ( radius, theta2, theta3 );
    area3 = circle_area_2d ( radius );

    cout << "\n";
    cout << "      CENTER       RADIUS   THETA1   THETA2   Area1   Area2  Circle\n";
    cout << "\n";
    cout << setw(9) << center[0]
         << setw(9) << center[1]
         << setw(9) << radius
         << setw(9) << theta1
         << setw(9) << theta2
         << setw(9) << area1
         << setw(9) << area2
         << setw(9) << area3 << "\n";
    cout << "\n";
    cout << "       F   Sector1       Sector2         Sum         Circle\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;
      function_2d_name ( name );

      resulta = circle_sector ( function_2d, center, radius, theta1, 
        theta2, nr );

      resultb = circle_sector ( function_2d, center, radius, theta2, 
        theta3, nr );

      result1 = resulta + resultb;

      result2 = circle_rt_sum ( function_2d, center, radius, nr2, ra, rw, 
        nt, ta, tw, zw );

      cout << "  " << name
           << setw(14) << resulta
           << setw(14) << resultb
           << setw(14) << result1
           << setw(14) << result2 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests CIRCLE_RT_SET and CIRCLE_RT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 1.0, 1.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int nc;
  int num;
  int nr;
  int nt;
  double pi = 3.141592653589793;
  double r = 1.0;
  double *ra;
  double result;
  double *rw;
  int rule;
  int rule_max = 9;
  double *ta;
  double *tw;
  double zw;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For R, Theta product rules on the unit circle,\n";
  cout << "  CIRCLE_RT_SET sets a rule.\n";
  cout << "  CIRCLE_RT_SUM uses the rule in an arbitrary circle.\n";
  cout << "\n";
  cout << "  We use a radius " << r << "\n";
  cout << "  and center:\n";
  cout << "  CENTER = " << setw(14) << center[0] 
                        << setw(14) << center[1] << "\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    cout << "\n";
    cout << "  Rule:  ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << "       " << setw(7) << rule;
    }
    cout << "\n";
    cout << "Function\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        circle_rt_size ( rule, &nr, &nt, &nc );

        ra = new double[nr];
        rw = new double[nr];
        ta = new double[nt];
        tw = new double[nt];

        circle_rt_set ( rule, nr, nt, nc, ra, rw, ta, tw, &zw );

        result = circle_rt_sum ( function_2d, center, r, nr, ra, rw, nt, ta, 
          tw, zw );
        cout << setw(14) << result;

        delete [] ra;
        delete [] rw;
        delete [] ta;
        delete [] tw;
      }
      cout << "\n";      
    }
  }
  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests CIRCLE_XY_SET and CIRCLE_XY_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 1.0, 1.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double r;
  double result;
  int rule;
  int rule_max = 13;
  double *weight;
  double *xtab;
  double *ytab;

  r = 1.0;
  cout << "\n";
  cout << "TEST16\n";
  cout << "  CIRCLE_XY_SET sets a quadrature rule\n";
  cout << "    for the unit circle.\n";
  cout << "  CIRCLE_XY_SUM evaluates the quadrature rule\n";
  cout << "    in an arbitrary circle.\n";
  cout << "\n";
  cout << "  We use a radius " << r << "\n";
  cout << "  and center:\n";
  cout << "  CENTER = (" << center[0] << ", " << center[1] << ").\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    cout << "\n";
    cout << "  Rule:  ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << "       " << setw(7) << rule;
    }
    cout << "\n";
    cout << " 'Function\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );

      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = circle_xy_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        circle_xy_set ( rule, order, xtab, ytab, weight );

        result = circle_xy_sum ( function_2d, center, r, order, xtab, ytab,
          weight );

        cout << setw(14) << result;

        delete [] weight;
        delete [] xtab;
        delete [] ytab;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test163 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST163 tests the rules for CN with Gegenbauer weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  cout << "\n";
  cout << "TEST163\n";
  cout << "  Demonstrate the use of quadrature rules for the region\n";
  cout << "  CN_GEG, that is, the hypercube [-1,+1]^N, with the\n";
  cout << "  weight W(ALPHA;X) = product ( 1 <= I <= N )\n";
  cout << "    (1-X(I)^2)^ALPHA\n";
  cout << "\n";
  cout << "  We use the formulas to integrate various monomials of\n";
  cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  cout << "  and compare to the exact integral.\n";
  cout << "\n";
  cout << "  The precision of each formula is known, and we only use\n";
  cout << "  a formula if its precision indicates it should be able to\n";
  cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      cn_geg_test ( n, alpha, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_geg_test ( n, alpha, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];

        i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_geg_test ( n, alpha, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_geg_test ( n, alpha, expon );

    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void cn_geg_test ( int n, double alpha, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    cout << setw(4) << expon[i];
  }
  cout << "\n";
  d = i4vec_sum ( n, expon );
  cout << "  Degree = " << d << "\n";
  cout << "\n";

  exact = cn_geg_monomial_integral ( n, alpha, expon );

  p = 0;

  if ( d <= p )
  {
    o = cn_geg_00_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    cn_geg_00_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_GEG_00_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 1;

  if ( d <= p )
  {
    o = cn_geg_01_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    cn_geg_01_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_GEG_01_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_geg_02_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    cn_geg_02_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_GEG_02_XIU: "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / ( 2.0 * alpha + 3.0 );
    volume_1d = sqrt ( pi ) * r8_gamma ( alpha + 1.0 ) 
      / r8_gamma ( alpha + 1.5 );
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  GW_02_XIU:     "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o = cn_geg_03_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    cn_geg_03_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_GEG_03_XIU: "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  cout << "  EXACT                  "
       << "  " << setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test165 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST165 tests the rules for CN with Jacobi weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double alpha;
  double alpha_test[TEST_NUM] = { 0.0, 1.0, 0.0, 0.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.0, 0.0, 2.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  cout << "\n";
  cout << "TEST165\n";
  cout << "  Demonstrate the use of quadrature rules for the region\n";
  cout << "  CN_JAC, that is, the hypercube [-1,+1]^N, with the\n";
  cout << "  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )\n";
  cout << "    (1-X(I))^ALPHA (1+X(I))^BETA\n";
  cout << "\n";
  cout << "  We use the formulas to integrate various monomials of\n";
  cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  cout << "  and compare to the exact integral.\n";
  cout << "\n";
  cout << "  The precision of each formula is known, and we only use\n";
  cout << "  a formula if its precision indicates it should be able to\n";
  cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      cn_jac_test ( n, alpha, beta, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_jac_test ( n, alpha, beta, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        beta  = beta_test[test];

        i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_jac_test ( n, alpha, beta, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_jac_test ( n, alpha, beta, expon );

    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void cn_jac_test ( int n, double alpha, double beta, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA =  " << beta << "\n";
  cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    cout << setw(4) << expon[i];
  }
  cout << "\n";
  d = i4vec_sum ( n, expon );
  cout << "  Degree = " << d << "\n";
  cout << "\n";

  exact = cn_jac_monomial_integral ( n, alpha, beta, expon );

  p = 0;

  if ( d <= p )
  {
    o = cn_jac_00_1_size ( n, alpha, beta );
    x = new double[n*o];
    w = new double[o];
    cn_jac_00_1 ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_JAC_00_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 1;

  if ( d <= p )
  {
    o = cn_jac_01_1_size ( n, alpha, beta );
    x = new double[n*o];
    w = new double[o];
    cn_jac_01_1 ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_JAC_01_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_jac_02_xiu_size ( n, alpha, beta );
    x = new double[n*o];
    w = new double[o];
    cn_jac_02_xiu ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  CN_JAC_02_XIU: "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = ( alpha + beta + 2.0 ) / 2.0;
    delta0 = ( alpha - beta ) / 2.0;
    c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 ) 
      / ( alpha + beta + 2.0 );
    volume_1d = pow ( 2.0, alpha + beta + 1.0 ) * r8_gamma ( alpha + 1.0 ) 
      * r8_gamma ( beta + 1.0 ) / ( alpha + beta + 1.0 ) / r8_gamma ( alpha + beta + 1.0 );
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }
  cout << "  EXACT                  "
       << "  " << setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test167 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST167 tests the rules for CN with Legendre weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *expon;
  int i;
  int n;
  int test;

  std::cout << "\n";
  std::cout << "TEST167\n";
  std::cout << "  Demonstrate the use of quadrature rules for the region\n";
  std::cout << "  CN_LEG, that is, the hypercube [-1,+1]^N, with the\n";
  std::cout << "  Legendre weight W(X) = 1.\n";
  std::cout << "\n";
  std::cout << "  We use the formulas to integrate various monomials of\n";
  std::cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  std::cout << "  and compare to the exact integral.\n";
  std::cout << "\n";
  std::cout << "  The precision of each formula is known, and we only use\n";
  std::cout << "  a formula if its precision indicates it should be able to\n";
  std::cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    i4vec_zero ( n, expon );
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      cn_leg_test ( n, expon );
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[0] = 3;
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 4;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 3;
      expon[1] = 2;
      cn_leg_test ( n, expon );
    }

    delete [] expon;
  }

  return;
}
//****************************************************************************80

void cn_leg_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "  N = " << n << "\n";
  std::cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    std::cout << std::setw(4) << expon[i];
  }
  std::cout << "\n";
  d = i4vec_sum ( n, expon );
  std::cout << "  Degree = " << d << "\n";
  std::cout << "\n";

  exact = cn_leg_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o = cn_leg_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    cn_leg_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  CN_LEG_01_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_leg_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    cn_leg_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  CN_LEG_02_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / 3.0;
    volume_1d = 2.0;
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:     "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o = cn_leg_03_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    cn_leg_03_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  CN_LEG_03_1:   "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = cn_leg_03_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    cn_leg_03_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  CN_LEG_03_XIU: "
         << "  " << std::setw(6) << o
         << "  " << std::setw(14) << quad
         << "  " << std::setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 5;

  if ( d <= p )
  {
    if ( 4 <= n && n <= 6 )
    {
      o = cn_leg_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      option = 1;
      cn_leg_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_1(1):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( 4 <= n && n <= 5 )
    {
      o = cn_leg_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      option = 2;
      cn_leg_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_1(2):"
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    if ( 2 <= n )
    {
      o = cn_leg_05_2_size ( n );
      x = new double[n*o];
      w = new double[o];
      cn_leg_05_2 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      std::cout << "  CN_LEG_05_2:   "
           << "  " << std::setw(6) << o
           << "  " << std::setw(14) << quad
           << "  " << std::setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }

  std::cout << "  EXACT                  "
       << "  " << std::setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests CONE_UNIT_3D, CONE_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double h;
  int i;
  char name[8];
  int num;
  double r = 1.0;
  double result;

  h = 1.0;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  CONE_UNIT_3D approximates integrals in a unit cone.\n";
  cout << "\n";
  cout << "  Volume = " << cone_volume_3d ( r, h ) << "\n";
  cout << "\n";
  cout << "    F(X)    CONE_3D\n";
  cout << "\n";

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result = cone_unit_3d ( function_3d );

    cout << "  " << name
         << "  " << setw(14) << result << "\n";
  }
  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests CUBE_SHELL_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 4;
  char name[8];
  int num;
  double r1;
  double r1_test[2] = { 0.0, 1.0 };
  double r2;
  double r2_test[2] = { 1.0, 2.0 };
  double result;
  int test;
  int test_num = 2;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  CUBE_SHELL_ND approximates integrals in a\n";
  cout << "    cubical shell in ND.\n";

  for ( test = 0; test < test_num; test++ )
  {
    r1 = r1_test[test];
    r2 = r2_test[test];

    cout << "\n";
    cout << "  Inner radius = " << r1 << "\n";
    cout << "  Outer radius = " << r2 << "\n";
    cout << "\n";
 
    for ( n = 2; n <= n_max; n++ )
    {
      cout << "\n";
      cout << "  Spatial dimension N = " << n << "\n";
      cout << "  Volume = " << cube_shell_volume_nd ( n, r1, r2 ) << "\n";
      cout << "\n";
      cout << "    F(X)      CUBE_SHELL_ND\n";
      cout << "\n";
 
      num = function_nd_num ( );

      for ( i = 1; i <= num; i++ )
      {
        function_nd_index = i;
        function_nd_name ( name );

        result = cube_shell_nd ( function_nd, n, r1, r2 );

        cout << "  " << name 
             << setw(14) << result << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests CUBE_UNIT_3D, QMULT_3D, RECTANGLE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a1;
  double a[3];
  double b1;
  double b[3];
  int i;
  int n = 3;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;

  a1 = -1.0;
  b1 = +1.0;

  a[0]= -1.0;
  a[1] = -1.0;
  a[2] = -1.0;
  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  CUBE_UNIT_3D approximates integrals\n";
  cout << "    in the unit cube in 3D.\n";
  cout << "  QMULT_3D approximates triple integrals.\n";
  cout << "  RECTANGLE_3D approximates integrals\n";
  cout << "    in a rectangular block.\n";
  cout << "\n";
  cout << "    F(X)      CUBE_UNIT_3D    QMULT_3D        RECTANGLE_3D\n";
  cout << "\n";

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = cube_unit_3d ( function_3d );
    result2 = qmult_3d ( function_3d, a1, b1, fu18, fl18, fu28, fl28 );
    result3 = rectangle_3d ( function_3d, a, b );

    cout << "  " << name
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2
         << "  " << setw(14) << result3 << "\n";
  }
  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests CUBE_UNIT_ND.
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
{
  int i;
  int i_test;
  int k;
  int k2;
  int khi;
  int klo;
  int k_test[2] = { 10, 5 };
  int max_k = 10;
  int max_test = 2;
  int n;
  int n_test[2] = { 2, 3 };
  char name[8];
  int num;
  double qa[10];
  double qb[10];

  cout << "\n";
  cout << "TEST20\n";
  cout << "  CUBE_UNIT_ND approximates integrals inside \n";
  cout << "    the unit cube in ND.\n";

  for ( i_test = 0; i_test < max_test; i_test++ )
  {
    n = n_test[i_test];
    k = k_test[i_test];
    cout << "\n";
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
    cout << "  Value of K = " << k << "\n";
    cout << "\n";
    cout << "    F(X)    CUBE_UNIT_ND\n";
    cout << "\n";

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      cube_unit_nd ( function_nd, qa, qb, n, k );

      for ( klo = 0; klo <= k - 1; klo = klo + 5 )
      {
        khi = i4_min ( klo + 4, k-1 );
        if ( klo == 0 )
        {
          cout << "  " << name;
        }
        else
        {
          cout << "  " << "       ";
        }
        for ( k2 = klo; k2 <= khi; k2++ )
        {
          cout << setw(14) << qa[k2];
        }
        cout << "\n";
      }
      for ( klo = 0; klo <= k - 1; klo = klo + 5 )
      {
        khi = i4_min ( klo + 4, k - 1 );
        cout << "  " << "       ";
        for ( k2 = klo; k2 <= khi; k2++ )
        {
          cout << setw(14) << qb[k2];
        }
        cout << "\n";
      }                                    
    }
  }
  return;
}
//****************************************************************************80

void test205 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST205 tests ELLIPSE_AREA_2D, ELLIPSE_CIRCUMFERENCE_2D, ELLIPSE_ECCENTRICITY_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double e;
  int i;
  double p;
  double r1;
  double r2;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST205\n";
  cout << "  ELLIPSE_AREA_2D returns the area of an ellipse.\n";
  cout << "  ELLIPSE_ECCENTRICITY_2D returns the\n";
  cout << "    eccentricity of an ellipse.\n";
  cout << "  ELLIPSE_CIRCUMFERENCE_2D returns the \n";
  cout << "    circumference of an ellipse.\n";
  cout << "\n";
  cout << "        R1        R2        E         Circum    Area\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    if ( i == 1 )
    {
      r1 = 25.0;
      r2 = 20.0;
    }
    else
    {
      r1 = r8_uniform_01 ( &seed );
      r2 = r8_uniform_01 ( &seed );
    }

    e = ellipse_eccentricity_2d ( r1, r2 );
    p = ellipse_circumference_2d ( r1, r2 ) ;   
    area = ellipse_area_2d ( r1, r2 );

    cout << "  " << setw(10) << r1
         << "  " << setw(10) << r2
         << "  " << setw(10) << e
         << "  " << setw(10) << p
         << "  " << setw(10) << area << "\n";
  }
  cout << "\n";
  cout << "  (For the first example, \n";
  cout << "  the eccentricity should be 0.6,\n";
  cout << "  the circumference should be about 141.8).\n";

  return;
}
//****************************************************************************80

void test207 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST207 tests the Stroud EN_R2 rules on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *expon;
  int i;
  int n;

  cout << "\n";
  cout << "TEST207\n";
  cout << "  Demonstrate the use of Stroud rules for the region\n";
  cout << "  EN_R2, that is, all of N-dimensional space, with the\n";
  cout << "  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )\n";
  cout << "\n";
  cout << "  We use the formulas to integrate various monomials of\n";
  cout << "  the form X1^EXPON1 * X2^EXPON2 * ... XN^EXPONN\n";
  cout << "  and compare to the exact integral.\n";
  cout << "\n";
  cout << "  The precision of each formula is known, and we only use\n";
  cout << "  a formula if its precision indicates it should be able to\n";
  cout << "  produce an exact result.\n";

  for ( n = 1; n <= 7; n++ )
  {
    expon = new int[n];

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[0] = 2;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[1] = 4;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 3 % n;
    expon[i] = 6;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[0] = 2;
    expon[1] = 4;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 4 % n;
    expon[i] = 8;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 5 % n;
    expon[i] = 10;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = i + 1;
    }
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 2;
    }
    en_r2_test ( n, expon );

    delete [] expon;
  }

  return;
}
//****************************************************************************80

void en_r2_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_TEST tests the Stroud EN_R2 rules on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << expon[i];
  }
  cout << "\n";
  d = i4vec_sum ( n, expon );
  cout << "  Degree = " << d << "\n";
  cout << "\n";

  exact = en_r2_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o =  en_r2_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_01_1:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o =  en_r2_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_02_XIU:  "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = 2.0;
    delta0 = 0.0;
    c1 = 1.0;
    volume_1d = sqrt ( pi );
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  GW_02_XIU:     "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 3;

  if ( d <= p )
  {
    o =  en_r2_03_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_03_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_03_1:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o =  en_r2_03_2_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_03_2 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_03_2:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o =  en_r2_03_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_03_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_03_XIU:  "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 5;

  if ( d <= p )
  {
    if ( 2 <= n && n <= 7 )
    {
      option = 1;
      o =  en_r2_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_05_1(1): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    if ( n == 3 || n == 5 || n == 6 )
    {
      option = 2;
      o =  en_r2_05_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_05_1(2): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    o =  en_r2_05_2_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_05_2 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_05_2:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    if ( 3 <= n )
    {
      o =  en_r2_05_3_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_05_3 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_05_3:    "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    o =  en_r2_05_4_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_05_4 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_05_4:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o =  en_r2_05_5_size ( n );
    x = new double[n*o];
    w = new double[o];
    en_r2_05_5 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EN_R2_05_5:    "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    if ( 5 <= n )
    {
      o =  en_r2_05_6_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_05_6 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_05_6:    "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }
  p = 7;

  if ( d <= p )
  {
    if ( n == 3 || n == 4 || n == 6 || n == 7 )
    {
      option = 1;
      o =  en_r2_07_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_07_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_07_1(1): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( n == 3 || n == 4 )
    {
      option = 2;
      o =  en_r2_07_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_07_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_07_1(2): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( 3 <= n )
    {
      o =  en_r2_07_2_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_07_2 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_07_2:    "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }

    if ( 3 <= n && n <= 6 )
    {
      option = 1;
      o =  en_r2_07_3_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_07_3 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_07_3(1): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    if ( n == 3 || n == 4 )
    {
      option = 2;
      o =  en_r2_07_3_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_07_3 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_07_3(2): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }

  p = 9;

  if ( d <= p )
  {
    if ( 3 <= n && n <= 6 )
    {
      option = 1;
      o =  en_r2_09_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_09_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_09_1(1): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;

      option = 2;
      o =  en_r2_09_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_09_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_09_1(2): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }

  p = 11;

  if ( d <= p )
  {
    if ( 3 <= n && n <= 5 )
    {
      option = 1;
      o =  en_r2_11_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_11_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_11_1(1): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;

      option = 2;
      o =  en_r2_11_1_size ( n );
      x = new double[n*o];
      w = new double[o];
      en_r2_11_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      cout << "  EN_R2_11_1(2): "
           << "  " << setw(6) << o
           << "  " << setw(14) << quad
           << "  " << setw(14) << err << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
  }
  cout << "  EXACT                  "
       << "  " << setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test2075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2075 tests the rules for EPN with GLG on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 2.0 };
  int *expon;
  int i;
  int n;
  int test;

  cout << "\n";
  cout << "TEST2075\n";
  cout << "  Demonstrate the use of quadrature rules for the region\n";
  cout << "  EPN_GLG, that is, the positive half space [0,+oo)^N, with the\n";
  cout << "  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )\n";
  cout << "\n";
  cout << "  We use the formulas to integrate various monomials of\n";
  cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  cout << "  and compare to the exact integral.\n";
  cout << "\n";
  cout << "  The precision of each formula is known, and we only use\n";
  cout << "  a formula if its precision indicates it should be able to\n";
  cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    i4vec_zero ( n, expon );
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        epn_glg_test ( n, expon, alpha );
      }
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }
    delete [] expon;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void epn_glg_test ( int n, int expon[], double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    cout << setw(4) << expon[i];
  }
  cout << "\n";
  d = i4vec_sum ( n, expon );
  cout << "  Degree = " << d << "\n";
  cout << "\n";

  exact = epn_glg_monomial_integral ( n, expon, alpha );

  p = 0;

  if ( d <= p )
  {
    o =  epn_glg_00_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    epn_glg_00_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_GLG_00_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 1;

  if ( d <= p )
  {
    o =  epn_glg_01_1_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    epn_glg_01_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_GLG_01_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o =  epn_glg_02_xiu_size ( n, alpha );
    x = new double[n*o];
    w = new double[o];
    epn_glg_02_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_GLG_02_XIU: "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = alpha + 1.0;
    c1 = - alpha - 1.0;
    volume_1d = r8_gamma ( 1.0 + alpha );
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    std::cout << "  GW_02_XIU:      "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }
  cout << "  EXACT                   "
       << "  " << setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test208 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST208 tests the rules for EPN with Laguerre weight on monomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *expon;
  int i;
  int n;
  int test;

  cout << "\n";
  cout << "TEST208\n";
  cout << "  Demonstrate the use of quadrature rules for the region\n";
  cout << "  EPN_LAG, that is, the positive half space [0,+oo)^N, with the\n";
  cout << "  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )\n";
  cout << "\n";
  cout << "  We use the formulas to integrate various monomials of\n";
  cout << "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n";
  cout << "  and compare to the exact integral.\n";
  cout << "\n";
  cout << "  The precision of each formula is known, and we only use\n";
  cout << "  a formula if its precision indicates it should be able to\n";
  cout << "  produce an exact result.\n";

  for ( n = 1; n <= 6; n++ )
  {
    expon = new int[n];

    i4vec_zero ( n, expon );
    epn_lag_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    epn_lag_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      epn_lag_test ( n, expon );
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    epn_lag_test ( n, expon );

    delete [] expon;
  }

  return;
}
//****************************************************************************80

void epn_lag_test ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  EXPON = ";
  for ( i = 0; i < n; i++ )
  {
    cout << setw(4) << expon[i];
  }
  cout << "\n";
  d = i4vec_sum ( n, expon );
  cout << "  Degree = " << d << "\n";
  cout << "\n";

  exact = epn_lag_monomial_integral ( n, expon );

  p = 0;

  if ( d <= p )
  {
    o =  epn_lag_00_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    epn_lag_00_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_LAG_00_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 1;

  if ( d <= p )
  {
    o =  epn_lag_01_1_size ( n );
    x = new double[n*o];
    w = new double[o];
    epn_lag_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_LAG_01_1:   "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  p = 2;

  if ( d <= p )
  {
    o =  epn_lag_02_xiu_size ( n );
    x = new double[n*o];
    w = new double[o];
    epn_lag_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  EPN_LAG_02_XIU: "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    o = gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = 1.0;
    c1 = - 1.0;
    volume_1d = 1.0;
    x = new double[n*o];
    w = new double[o];
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    cout << "  GW_02_XIU:      "
         << "  " << setw(6) << o
         << "  " << setw(14) << quad
         << "  " << setw(14) << err << "\n";
    delete [] v;
    delete [] w;
    delete [] x;
  }

  cout << "  EXACT                   "
       << "  " << setw(14) << exact << "\n";

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST21 tests HEXAGON_UNIT_SET and HEXAGON_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double rad;
  double result;
  int rule;
  int rule_max = 4;
  double *weight;
  double *xtab;
  double *ytab;

  rad = 2.0;
 
  cout << "\n";
  cout << "TEST21\n";
  cout << "  HEXAGON_UNIT_SET sets a quadrature rule for the\n";
  cout << "    unit hexagon.\n";
  cout << "  HEXAGON_SUM evaluates the quadrature rule\n";
  cout << "    in an arbitrary hexagon.\n";
  cout << "\n";
  cout << "  We use a radius " << rad << "\n";
  cout << "  and center:\n";
  cout << "  CENTER = (" << center[0] << ", " << center[1] << ")\n"; 
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << setw(6) << rule;
    }
    cout << "\n";
    cout << "  Function\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = hexagon_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        hexagon_unit_set ( rule, order, xtab, ytab, weight );

        result = hexagon_sum ( function_2d, center, rad, order, xtab, ytab, 
          weight );

        cout << setw(14) << result;

        delete [] xtab;
        delete [] ytab;
        delete [] weight;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test215 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST215 tests LENS_HALF_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int i_max = 8;
  int order;
  double pi = 3.141592653589793;
  double r;
  double theta1;
  double theta2;
  double value;

  r = 1.0;

  cout << "\n";
  cout << "TEST215\n";
  cout << "  LENS_HALF_2D approximates an integral within a\n";
  cout << "    circular half lens, defined by joining the endpoints\n";
  cout << "    of a circular arc.\n";
  cout << "\n";
  cout << "  Integrate F(X,Y) = 1\n";
  cout << "\n";
  cout << "      R            Theta1      Theta2        "
       << "Area        Order Integral\n";
  cout << "\n";

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    cout << "\n";

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_1_2d, center, r, theta1, theta2, order );
      cout << setw(14) << r
           << setw(14) << theta1
           << setw(14) << theta2
           << setw(14) << area
           << setw(8)  << order
           << setw(14) << value << "\n";
    }
  }

  cout << "\n";
  cout << "\n";
  cout << "  Integrate F(X,Y) = X\n";
  cout << "\n";
  cout << "      R            Theta1      Theta2        "
       << "Area        Order Integral\n";
  cout << "\n";

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    cout << "\n";

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_x_2d, center, r, theta1, theta2, order );
      cout << setw(14) << r
           << setw(14) << theta1
           << setw(14) << theta2
           << setw(14) << area
           << setw(8)  << order
           << setw(14) << value << "\n";
    }
  }
  cout << "\n";
  cout << "\n";
  cout << "  Integrate F(X,Y) = R\n";
  cout << "\n";
  cout << "      R            Theta1      Theta2        "
       << "Area        Order Integral\n";
  cout << "\n";

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    cout << "\n";

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_r_2d, center, r, theta1, theta2, order );
      cout << setw(14) << r
           << setw(14) << theta1
           << setw(14) << theta2
           << setw(14) << area
           << setw(8)  << order
           << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests OCTAHEDRON_UNIT_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  OCTAHEDRON_UNIT_ND approximates integrals in a unit\n";
  cout << "    octahedron in N dimensions.\n";
  cout << "\n";
  cout << "\n";
  cout << "    F(X)    N = 1    N = 2   N = 3\n";
  cout << "\n";
 
  num = function_nd_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_nd_index = i;
    function_nd_name ( name );
    cout << "  " << name;

    for ( n = 1; n <= n_max; n++ )
    {
      result = octahedron_unit_nd ( function_nd, n );
      cout << setw(14) << result;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests PARALLELIPIPED_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  double *v;
  double volume;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  PARALLELIPIPED_VOLUME_ND computes the volume of a\n";
  cout << "    parallelipiped in N dimensions.\n";
  cout << "\n";

  for ( n = 2; n <= 4; n++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
//
//  Set the values of the parallelipiped.
//
    v = setsim ( n );

    cout << "\n";
    cout << "  Parallelipiped vertices:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        cout << setw(14) << v[i+j*n];
      }
      cout << "\n";
    }
    volume = parallelipiped_volume_nd ( n, v );

    cout << "\n";
    cout << "  Volume is " << volume << "\n";

    delete [] v;
  }
  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests POLYGON_**_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double result;
  int npts = 4;
  double x[4] = { 0.0, 1.0, 1.0, 0.0 };
  double y[4] = { 0.0, 0.0, 1.0, 1.0 };

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_1_2D integrates 1\n";
  cout << "  POLYGON_X_2D integrates X\n";
  cout << "  POLYGON_Y_2D integrates Y\n";
  cout << "  POLYGON_XX_2D integrates X*X\n";
  cout << "  POLYGON_XY_2D integrates X*Y\n";
  cout << "  POLYGON_YY_2D integrates Y*Y\n";
  cout << "\n";
  cout << "  F(X,Y)    Integral\n";
  cout << "\n";
 
  result = polygon_1_2d ( npts, x, y );
  cout << "     1    " << result << "\n";
 
  result = polygon_x_2d ( npts, x, y );
  cout << "     X    " << result << "\n";
 
  result = polygon_y_2d ( npts, x, y );
  cout << "     Y    " << result << "\n";
 
  result = polygon_xx_2d ( npts, x, y );
  cout << "   X*X    " << result << "\n";
 
  result = polygon_xy_2d ( npts, x, y );
  cout << "   X*Y    " << result << "\n";
 
  result = polygon_yy_2d ( npts, x, y );
  cout << "   Y*Y    " << result << "\n";
 
  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests PYRAMID_UNIT_O**_3D, PYRAMID_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  char name[8];
  int num;
  int order[10] = { 1, 5, 6, 8, 8, 9, 13, 18, 27, 48 };
  double result; 

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For the unit pyramid, we approximate integrals with:\n";
  cout << "  PYRAMID_UNIT_O01_3D, a 1 point rule.\n";
  cout << "  PYRAMID_UNIT_O05_3D, a 5 point rule.\n";
  cout << "  PYRAMID_UNIT_O06_3D, a 6 point rule.\n";
  cout << "  PYRAMID_UNIT_O08_3D, an 8 point rule.\n";
  cout << "  PYRAMID_UNIT_O08b_3D, an 8 point rule.\n";
  cout << "  PYRAMID_UNIT_O09_3D, a 9 point rule.\n";
  cout << "  PYRAMID_UNIT_O13_3D, a 13 point rule.\n";
  cout << "  PYRAMID_UNIT_O18_3D, a 18 point rule.\n";
  cout << "  PYRAMID_UNIT_O27_3D, a 27 point rule.\n";
  cout << "  PYRAMID_UNIT_O48_3D, a 48 point rule.\n";
  cout << "\n";
  cout << "  PYRAMID_UNIT_VOLUME_3D computes the volume of a unit pyramid.\n";
  cout << "\n";
  cout << "  Volume = " << pyramid_unit_volume_3d ( ) << "\n";
  cout << "\n";
    
  num = function_3d_num ( );

  for ( jlo = 1; jlo <= num; jlo = jlo + 5 )
  {
    jhi = i4_min ( jlo + 4, num );
    cout << "\n";
    cout << "  Order   ";
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      function_3d_name ( name );
      cout << "   " << name << "    ";
    }
    cout << "\n";
    cout << "\n";
    cout << "  " << setw(5) << order[0];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o01_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[1];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o05_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[2];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o06_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[3];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o08_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[4];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o08b_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[5];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o09_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[6];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o13_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[7];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o18_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[8];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o27_3d ( function_3d );
    }
    cout << "\n";
    cout << "  " << setw(5) << order[9];
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      cout << setw(14) << pyramid_unit_o48_3d ( function_3d );
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test255 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST255 tests PYRAMID_UNIT_MONOMIAL_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  int alpha;
  int beta;
  int degree_max = 4;
  int gamma;
  double value;

  cout << "\n";
  cout << "TEST255\n";
  cout << "  For the unit pyramid,\n";
  cout << "  PYRAMID_UNIT_MONOMIAL_3D returns the exact value of the\n";
  cout << " integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << pyramid_unit_volume_3d ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        value = pyramid_unit_monomial_3d ( alpha, beta, gamma );

        cout << "  " << setw(8)  << alpha
             << "  " << setw(8)  << beta
             << "  " << setw(8)  << gamma
             << "  " << setw(14) << value << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests QMULT_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  char name[8];
  int num;
  double result;

  a = -1.0;
  b = 1.0;
 
  cout << "\n";
  cout << "TEST26\n";
  cout << "  QMULT_1D approximates an integral on a\n";
  cout << "    one-dimensional interval.\n";
  cout << "\n";
  cout << "  We use the interval:\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "\n";
  cout << "    F(X)     QMULT_1D\n";
  cout << "\n";
 
  num = function_1d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_1d_index = i;
    function_1d_name ( name );

    result = qmult_1d ( function_1d, a, b );
    cout << "  " << name
         << "  " << setw(14) << result << "\n"; 
  }

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests SIMPLEX_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  char name[8];
  int num;
  double result;
  double *v;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  SIMPLEX_ND approximates integrals inside an\n";
  cout << "    arbitrary simplex in ND.\n";
  cout << "\n";
 
  for ( n = 2; n <= 4; n++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
//
//  Restore values of simplex.
//
    v = setsim ( n );

    cout << "\n";
    cout << "  Simplex vertices:\n";
    cout << "\n";
 
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        cout << setw(4) << v[i+j*n];
      }
      cout << "\n";
    }
 
    cout << "\n";
    cout << "  F(X)    SIMPLEX_ND\n";
    cout << "\n";

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result = simplex_nd ( function_nd, n, v );
      cout << "  " << name
           << "  " << setw(14) << result << "\n";
 
      v = setsim ( n );
    }
    delete [] v;
  }
  return;
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests SIMPLEX_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  double *v;
  double volume;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  SIMPLEX_VOLUME_ND computes the volume of a simplex\n";
  cout << "    in N dimensions.\n";
  cout << "\n";

  for ( n = 2; n <= 4; n++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
//
//  Set the values of the simplex.
//
    v = setsim ( n );

    cout << "\n";
    cout << "  Simplex vertices:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        cout << setw(14) << v[i+j*n];
      }
      cout << "\n";
    }
    volume = simplex_volume_nd ( n, v );

    cout << "\n";
    cout << "  Volume is " << volume << "\n";

    delete [] v;
  }
  return;
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests SIMPLEX_UNIT_**_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  char name[8];
  double result1;
  double result2;
  double result3;
  double result4;
  double volume;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  For integrals in the unit simplex in ND,\n";
  cout << "  SIMPLEX_UNIT_01_ND uses a formula of degree 1.\n";
  cout << "  SIMPLEX_UNIT_03_ND uses a formula of degree 3.\n";
  cout << "  SIMPLEX_UNIT_05_ND uses a formula of degree 5.\n";
  cout << "  SIMPLEX_UNIT_05_2_ND uses a formula of degree 5.\n";

  for ( i = 1; i <= 6; i++ )
  {
    function_nd_index = i;
    function_nd_name ( name );

    cout << "\n";
    cout << "  Check the integral of " << name << "\n";
    cout << "\n";
    cout << "  N     Volume         #1              #3"
         << "              #5              #5.2\n";
    cout << "\n";

    for ( n = 2; n <= 16; n++ )
    {
      result1 = simplex_unit_01_nd ( function_nd, n );
      result2 = simplex_unit_03_nd ( function_nd, n );
      result3 = simplex_unit_05_nd ( function_nd, n );
      result4 = simplex_unit_05_2_nd ( function_nd, n );

      volume = simplex_unit_volume_nd ( n );

      cout << "  " << setw(2) << n
           << "  " << setw(13) << volume
           << "  " << setw(13) << result1
           << "  " << setw(13) << result2
           << "  " << setw(13) << result3
           << "  " << setw(13) << result4 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests SPHERE_UNIT_**_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double result4;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  For integrals on the unit sphere in 3D:\n";
  cout << "  SPHERE_UNIT_07_3D uses a formula of degree 7.\n";
  cout << "  SPHERE_UNIT_11_3D uses a formula of degree 11.\n";
  cout << "  SPHERE_UNIT_14_3D uses a formula of degree 14.\n";
  cout << "  SPHERE_UNIT_15_3D uses a formula of degree 15.\n";
  cout << "\n";
  cout << "  Unit sphere area = " << sphere_unit_area_nd ( 3 ) << "\n";
  cout << "\n";
  cout << "    F(X)    S3S07        S3S11         S3S14 "
       << "         S3S15\n";
  cout << "\n";
 
  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = sphere_unit_07_3d ( function_3d );
    result2 = sphere_unit_11_3d ( function_3d );
    result3 = sphere_unit_14_3d ( function_3d );
    result4 = sphere_unit_15_3d ( function_3d );
 
    cout << "  " << name
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2
         << "  " << setw(14) << result3
         << "  " << setw(14) << result4 << "\n";
  }
  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests SPHERE_UNIT_**_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double result4;
  double result5;
  double result6;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  For integrals on the unit sphere in ND:\n";
  cout << "  SPHERE_UNIT_03_ND uses a formula of degree 3;\n";
  cout << "  SPHERE_UNIT_04_ND uses a formula of degree 4;\n";
  cout << "  SPHERE_UNIT_05_ND uses a formula of degree 5.\n";
  cout << "  SPHERE_UNIT_07_1_ND uses a formula of degree 7.\n";
  cout << "  SPHERE_UNIT_07_2_ND uses a formula of degree 7.\n";
  cout << "  SPHERE_UNIT_11_ND uses a formula of degree 11.\n";
  cout << "\n";

  for ( n = 3; n <= 10; n++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
    cout << "  Unit sphere area = " << sphere_unit_area_nd ( n ) << "\n";
    cout << "\n";
    cout << "    Rule:     #3            #4            #5\n";          
    cout << "              #7.1          #7.2          #11\n";
    cout << "    Function\n";
    cout << "\n";

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = sphere_unit_03_nd ( function_nd, n );
      result2 = sphere_unit_04_nd ( function_nd, n );
      result3 = sphere_unit_05_nd ( function_nd, n );
      result4 = sphere_unit_07_1_nd ( function_nd, n );
      result5 = sphere_unit_07_2_nd ( function_nd, n );
      result6 = sphere_unit_11_nd ( function_nd, n );

      cout << "  " << name
           << setw(14) << result1
           << setw(14) << result2
           << setw(14) << result3 << "\n";
      cout << "  " << "       "
           << setw(14) << result4
           << setw(14) << result5
           << setw(14) << result6 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests SPHERE_05_ND, SPHERE_07_1_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *center;
  int dim;
  int i;
  int n;
  char name[8];
  int num;
  double r;
  double result1;
  double result2;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  For integrals on a sphere in ND:\n";
  cout << "  SPHERE_05_ND uses a formula of degree 5.\n";
  cout << "  SPHERE_07_1_ND uses a formula of degree 7.\n";
  cout << "\n";
  r = 2.0;

  for ( n = 2; n <= 4; n++ )
  {
    center = new double[n];
    for ( dim = 0; dim < n; dim++ )
    {
      center[dim] = 1.0;
    }
    cout << "\n";
    cout << "  Spatial dimension N = " << n << "\n";
    cout << "  Sphere center = \n";
    for ( dim = 0; dim < n; dim++ )
    {
      cout << setw(14) << center[dim];
    }
    cout << "\n";
    cout << "  Sphere radius = " << r << "\n";
    cout << "  Sphere area = " << sphere_area_nd ( n, r ) << "\n";
    cout << "\n";
    cout << "    Rule:     #5           #7.1\n";
    cout << "    Function\n";
    cout << "\n";

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = sphere_05_nd ( function_nd, n, center, r );
      result2 = sphere_07_1_nd ( function_nd, n, center, r );

      cout << "  " << name
           << setw(14) << result1
           << setw(14) << result2 << "\n";
    }
    delete [] center;
  }
  return;
}
//****************************************************************************80

void test322 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST322 tests SPHERE_CAP_AREA_3D, SPHERE_CAP_AREA_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area1;
  double area2;
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  double h;
  int i;
  int ntest = 12;
  double r = 1.0;

  cout << "\n";
  cout << "TEST322\n";
  cout << "  SPHERE_CAP_AREA_3D computes the volume of a\n";
  cout << "    3D spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_AREA_ND computes the volume of an\n";
  cout << "    ND spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";

  area1 = sphere_area_3d ( r );

  cout << "\n";
  cout << "  Area of the total sphere in 3D = " << area1 << "\n";

  cout << "\n";
  cout << "        R           H           Cap         Cap\n";
  cout << "                                area_3d     area_nd\n";
  cout << "\n";

  for ( i= 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    area1 = sphere_cap_area_3d ( r, h );

    area2 = sphere_cap_area_nd ( dim_num, r, h );

    cout << "  " << setw(12) << r
         << "  " << setw(12) << h
         << "  " << setw(12) << area1
         << "  " << setw(12) << area2 << "\n";
  }
  return;
}
//****************************************************************************80

void test324 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST324 tests SPHERE_CAP_VOLUME_2D, SPHERE_CAP_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  double h;
  int i;
  int ntest = 12;
  double pi = 3.141592653589793;
  double r = 1.0;
  double volume1;
  double volume2;

  cout << "\n";
  cout << "TEST324\n";
  cout << "  SPHERE_CAP_VOLUME_2D computes the volume (area) of a\n";
  cout << "    spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_VOLUME_ND does the same operation,\n";
  cout << "    but in N dimensions.\n";
  cout << "\n";
  cout << "  Using a radius R = " << r << "\n";

  volume1 = sphere_volume_2d ( r );

  cout << "\n";
  cout << "  Volume of the total sphere in 2D = " << volume1 << "\n";

  cout << "\n";
  cout << "        H           Cap        Cap\n";
  cout << "                    vol_2d     vol_nd\n";
  cout << "\n";

  for ( i = 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    volume1 = sphere_cap_volume_2d ( r, h );

    volume2 = sphere_cap_volume_nd ( dim_num, r, h );

    cout << "  " << setw(12) << h
         << "  " << setw(12) << volume1
         << "  " << setw(12) << volume2 << "\n";
  }
  return;
}
//****************************************************************************80

void test326 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST326 tests SPHERE_CAP_VOLUME_3D, SPHERE_CAP_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  double h;
  int i;
  int ntest = 12;
  double r = 1.0;
  double volume1;
  double volume2;
  
  cout << "\n";
  cout << "TEST326\n";
  cout << "  SPHERE_CAP_VOLUME_3D computes the volume of a\n";
  cout << "    spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_VOLUME_ND does the same operation,\n";
  cout << "    but in N dimensions.\n";
  cout << "\n";
  cout << "  Using a radius R = " << r << "\n";

  volume1 = sphere_volume_3d ( r );

  cout << "\n";
  cout << "  Volume of the total sphere in 3D = " << volume1 << "\n";

  cout << "\n";
  cout << "        H           Cap        Cap\n";
  cout << "                    volume_3d  volume_nd\n";
  cout << "\n";

  for ( i = 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    volume1 = sphere_cap_volume_3d ( r, h );

    volume2 = sphere_cap_volume_nd ( dim_num, r, h );

    cout << "  " << setw(12) << h
         << "  " << setw(12) << volume1
         << "  " << setw(12) << volume2 << "\n";
  }

  return;
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests SPHERE_CAP_AREA_ND, SPHERE_CAP_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double h;
  int i;
  int n = 12;
  int dim_num;
  double pi = 3.141592653589793;
  double r;
  double volume;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  For a sphere in ND:\n";
  cout << "  SPHERE_CAP_AREA_ND computes the area\n";
  cout << "    of a spherical cap.\n";
  cout << "  SPHERE_CAP_VOLUME_ND computes the volume\n";
  cout << "    of a spherical cap.\n";
  cout << "\n";

  r = 1.0;

  for ( dim_num = 2; dim_num <= 5; dim_num++ )
  {
    cout << "\n";
    cout << "  Spatial dimension N = " << dim_num << "\n";
    cout << "  Radius =       " << r << "\n";
    cout << "  Area =         " << sphere_area_nd ( dim_num, r ) << "\n";
    volume = sphere_volume_nd ( dim_num, r );
    cout << "  Volume =       " << volume << "\n";

    cout << "\n";
    cout << "                 Sphere         Sphere\n";
    cout << "                 cap            cap\n";
    cout << "      H          area           volume\n";
    cout << "\n";

    for ( i = 0; i <= n + 1; i++ )
    {
      h = ( double ) ( 2 * i ) * r / ( double ) ( n );
      area = sphere_cap_area_nd ( dim_num, r, h );
      volume = sphere_cap_volume_nd ( dim_num, r, h );
      cout << "  " << setw(8) << h
           << "  " << setw(14) << area
           << "  " << setw(14) << volume << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test335 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST335 tests SPHERE_SHELL_03_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[3];
  int i;
  int j;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result3;
  double result4;
  double result5;
  double result6;
  double r1;
  double r2;

  cout << "\n";
  cout << "TEST335\n";
  cout << "  For integrals inside a spherical shell in ND:\n";
  cout << "  SPHERE_SHELL_03_ND approximates the integral.\n";
  cout << "\n";
  cout << "  We compare these results with those computed by\n";
  cout << "  from the difference of two ball integrals:\n";
  cout << "\n";
  cout << "  BALL_F1_ND approximates the integral;\n";
  cout << "  BALL_F3_ND approximates the integral\n";
  cout << "\n";

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 ) 
    {
      r1 = 0.0;
      r2 = 1.0;
      center[0] = 0.0;
      center[1] = 0.0;
      center[2] = 0.0;
    }
    else
    {
      r1 = 2.0;
      r2 = 3.0;
      center[0] =  1.0;
      center[1] = -1.0;
      center[2] =  2.0;
    }

    for ( n = 2; n <= n_max; n++ )
    {
      cout << "\n";
      cout << "  Spatial dimension N = " << n << "\n";
      cout << "  Sphere center:\n";
      for ( i = 0; i < n; i++ )
      {
        cout << setw(10) << center[i];
      }
      cout << "  Inner sphere radius = " << r1 << "\n";
      cout << "  Outer sphere radius = " << r2 << "\n";
      cout << "  Spherical shell volume = "
           << sphere_shell_volume_nd ( n, r1, r2 ) << "\n";
      cout << "\n";
      cout << "\n";
      cout << "    Rule:      #3       F1(R2)-F1(R1)  "      
           << "F3(R2)-F3(R1)\n";
      cout << "    F(X)\n";
      cout << "\n";

      num = function_nd_num ( );

      for ( i = 1; i <= num; i++ )
      {
        function_nd_index = i;
        function_nd_name ( name );
        cout << "  " << name;

        result1 = sphere_shell_03_nd ( function_nd, n, center, r1, r2 );

        result3 = ball_f1_nd ( function_nd, n, center, r1 );
        result4 = ball_f1_nd ( function_nd, n, center, r2 );
      
        result5 = ball_f3_nd ( function_nd, n, center, r1 );
        result6 = ball_f3_nd ( function_nd, n, center, r2 );

        cout << setw(14) << result1
             << setw(14) << result4 - result3
             << setw(14) << result6 - result5 << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests SPHERE_UNIT_AREA_ND, SPHERE_UNIT_AREA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double area2;
  int dim_num;
  int n_data;

  cout << "\n";
  cout << "TEST34:\n";
  cout << "  SPHERE_UNIT_AREA_ND evaluates the area of the unit\n";
  cout << "  sphere in N dimensions.\n";
  cout << "  SPHERE_UNIT_AREA_VALUES returns some test values.\n";
  cout << "\n";
  cout << "     dim_num    Exact          Computed\n";
  cout << "                Area           Area\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( &n_data, &dim_num, &area );

    if ( n_data == 0 )
    {
      break;
    }
    area2 = sphere_unit_area_nd ( dim_num );

    cout << "  " << setw(8) << dim_num
         << "  " << setw(10) << area
         << "  " << setw(10) << area2 << "\n";
  }
  return;
}
//****************************************************************************80

void test345 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST345 tests SPHERE_UNIT_VOLUME_ND, SPHERE_UNIT_VOLUME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int n_data;
  double volume;
  double volume2;

  cout << "\n";
  cout << "TEST345:\n";
  cout << "  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit\n";
  cout << "  sphere in N dimensions.\n";
  cout << "  SPHERE_UNIT_VOLUME_VALUES returns some test values.\n";
  cout << "\n";
  cout << "     dim_num    Exact          Computed\n";
  cout << "                Volume         Volume\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( &n_data, &dim_num, &volume );

    if ( n_data == 0 )
    {
      break;
    }

    volume2 = sphere_unit_volume_nd ( dim_num );

    cout << "  " << setw(8) << dim_num
         << "  " << setw(10) << volume
         << "  " << setw(10) << volume2 << "\n";
  }
  return;
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests SQUARE_UNIT_SET, RECTANGLE_SUB_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  char name[8];
  int num;
  int order;
  int nsub[2];
  double result;
  int rule;
  double *weight;
  double *xtab;
  double xval[2] = { 1.0, 3.0 };
  double *ytab;
  double yval[2] = { 2.0, 3.0 };

  cout << "\n";
  cout << "TEST35\n";
  cout << "  SQUARE_UNIT_SET sets up a quadrature rule \n";
  cout << "    on a unit square.\n";
  cout << "  RECTANGLE_SUB_2D applies it to subrectangles of an\n";
  cout << "    arbitrary rectangle.\n";
  cout << "\n";

  cout << "\n";
  cout << "  The corners of the rectangle are:\n";
  cout << "\n";
  cout << setw(14) << xval[0] << setw(14) << yval[0] << "\n";
  cout << setw(14) << xval[1] << setw(14) << yval[1] << "\n";
//
//  Get the quadrature abscissas and weights for a unit square.
//
  rule = 2;
  order = square_unit_size ( rule );

  xtab = new double[order];
  ytab = new double[order];
  weight = new double[order];

  square_unit_set ( rule, order, xtab, ytab, weight );

  cout << "\n";
  cout << "  Using unit square integration rule number " << rule << "\n";
  cout << "  Order of rule is " << order << "\n";
//
//  Set the function.
//
  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );
//
//  Try an increasing number of subdivisions.
//
    cout << "\n";
    cout << "    Function  Subdivisions  Integral\n";
    cout << "\n";

    for ( j = 1; j <= 5; j++ )
    {
      nsub[0] = j;
      nsub[1] = 2 * j;

      result = rectangle_sub_2d ( function_2d, xval, yval, nsub, order, xtab, 
        ytab, weight );

      cout << "  " << name
           << setw(4)  << nsub[0]
           << setw(4)  << nsub[1]
           << setw(14) << result << "\n";
    }
  }

  delete [] weight;
  delete [] xtab;
  delete [] ytab;

  return;
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests SQUARE_UNIT_SET and SQUARE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  double center[2] = { 2.0, 2.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double r;
  double result;
  int rule;
  int rule_max = 6;
  double *weight;
  double *xtab;
  double *ytab;

  r = 3.0;

  cout << "\n";
  cout << "TEST36\n";
  cout << "  SQUARE_UNIT_SET sets up quadrature on the unit square;\n";
  cout << "  SQUARE_SUM carries it out on an arbitrary square.\n";
  cout << "\n";
  cout << "  Square center:\n";
  cout << "  CENTER = ( " << center[0] << ", " << center[1] << ")\n";
  cout << "  Square radius is " << r << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    cout << "\n";
    cout << "  Rule:";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << "       " << setw(6) << rule;
    }
    cout << "\n";
    cout << "  Function \n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );

      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = square_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        square_unit_set ( rule, order, xtab, ytab, weight );

        result = square_sum ( function_2d, center, r, order, xtab, ytab, 
          weight );

        cout << setw(13) << result;

        delete [] weight;
        delete [] xtab;
        delete [] ytab;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests SQUARE_UNIT_SET and SQUARE_UNIT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 6;
  double *weight;
  double *xtab;
  double *ytab;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  SQUARE_UNIT_SET sets up quadrature on the unit square;\n";
  cout << "  SQUARE_UNIT_SUM carries it out on the unit square.\n";
  cout << "\n";
 
  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << setw(6) << rule;
    }
    cout << "\n";
    cout << "  Function\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = square_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        square_unit_set ( rule, order, xtab, ytab, weight );

        result = square_unit_sum ( function_2d, order, xtab, ytab, weight );

        cout << setw(13) << result;

        delete [] weight;
        delete [] xtab;
        delete [] ytab;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests TETRA_07, TETRA_TPRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int num;
  int order;
  int order_hi;
  int order_lo;
  int order_max = 9;
  double result;
  double x[4] = { 1.0, 4.0, 1.0, 1.0 };
  double y[4] = { 2.0, 2.0, 3.0, 2.0 };
  double z[4] = { 6.0, 6.0, 6.0, 8.0 };

  cout << "\n";
  cout << "TEST38\n";
  cout << "  For integrals inside an arbitrary tetrahedron:\n";
  cout << "  TETRA_07 uses a formula of degree 7;\n";
  cout << "  TETRA_TPRODUCT uses a triangular product formula\n";
  cout << "    of varying degree.\n";
  cout << "\n";
  cout << "  Tetrahedron vertices:\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(4) << x[i]
         << setw(4) << y[i]
         << setw(4) << z[i] << "\n";
  }
  cout << "\n";
  cout << "  Tetrahedron unit volume = " << tetra_unit_volume ( ) << "\n";
  cout << "  Tetrahedron Volume = " << tetra_volume ( x, y, z ) << "\n";
  cout << "\n";
  cout << "\n";
  cout << "  F(X)    TETRA_07\n";
  cout << "          TETRA_TPRODUCT(1:4)\n";
  cout << "          TETRA_TPRODUCT(5:8)\n";
  cout << "          TETRA_TPRODUCT(9)\n";
  cout << "\n";

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );
    cout << "  " << name;
    result = tetra_07 ( function_3d, x, y, z );
    cout << setw(14) << result << "\n";

    for ( order_lo = 1; order_lo <= order_max; order_lo = order_lo + 4 )
    {
      order_hi = i4_min ( order_lo + 3, order_max );
      cout << "  " << "       ";
      for ( order = order_lo; order <= order_hi; order++ )
      {
        result = tetra_tproduct ( function_3d, order, x, y, z );
        cout << setw(16) << result;
      }
      cout << "\n";
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests TETRA_UNIT_SET and TETRA_UNIT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double *weight;
  double *xtab;
  double *ytab;
  double *ztab;

  cout << "\n";
  cout << "TEST39\n";
  cout << "  TETRA_UNIT_SET sets quadrature rules\n";
  cout << "    for the unit tetrahedron;\n";
  cout << "  TETRA_UNIT_SUM applies them to the unit tetrahedron.\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << setw(6) << rule;
    }
    cout << "\n";
    cout << "Function\n";
    cout << "\n";

    num = function_3d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_3d_index = i;
      function_3d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = tetra_unit_size ( rule );

        weight = new double[order];
        xtab = new double[order];
        ytab = new double[order];
        ztab = new double[order];

        tetra_unit_set ( rule, order, xtab, ytab, ztab, weight );
 
        result = tetra_unit_sum ( function_3d, order, xtab, ytab, ztab,
          weight );

        cout << setw(14) << result;

        delete [] weight;
        delete [] xtab;
        delete [] ytab;
        delete [] ztab;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests TETRA_UNIT_SET and TETRA_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double value;
  double *weight;
  double x[4] = { 1.0, 4.0, 1.0, 1.0 };
  double *xtab;
  double y[4] = { 2.0, 2.0, 3.0, 2.0 };
  double *ytab;
  double z[4] = { 6.0, 6.0, 6.0, 8.0 };
  double *ztab;

  cout << "\n";
  cout << "TEST40\n";
  cout << "  TETRA_UNIT_SET sets quadrature rules\n";
  cout << "    for the unit tetrahedron;\n";
  cout << "  TETRA_SUM applies them to an arbitrary tetrahedron.\n";
  cout << "\n";
  cout << "  Tetrahedron vertices:\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(6) << x[i]
         << "  " << setw(6) << y[i]
         << "  " << setw(6) << z[i] << "\n";
  }


  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << "       " << setw(7) << rule;
    }
    cout << "\n";
    cout << "  Function\n";
    cout << "\n";

    num = function_3d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_3d_index = i;
      function_3d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = tetra_unit_size ( rule );

        weight = new double[order];
        xtab = new double[order];
        ytab = new double[order];
        ztab = new double[order];

        tetra_unit_set ( rule, order, xtab, ytab, ztab, weight );

        result = tetra_sum ( function_3d, x, y, z, order, xtab, ytab, ztab, 
          weight );

        cout << setw(14) << result;

        delete [] weight;
        delete [] xtab;
        delete [] ytab;
        delete [] ztab;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests TRIANGLE_UNIT_SET, TRIANGLE_SUB.
//
//  Discussion:
//
//    Break up the triangle into NSUB*NSUB equal subtriangles.  Approximate 
//    the integral over the triangle by the sum of the integrals over each
//    subtriangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int nsub;
  int num;
  int order;
  double result;
  int rule;
  double *weight;
  double *xtab;
  double xval[3] = { 0.0, 0.0, 1.0 };
  double *ytab;
  double yval[3] = { 0.0, 1.0, 0.0 };

  cout << "\n";
  cout << "TEST41\n";
  cout << "  TRIANGLE_UNIT_SET sets up a quadrature rule\n";
  cout << "    on a triangle.\n";
  cout << "  TRIANGLE_SUB applies it to subtriangles of an\n";
  cout << "    arbitrary triangle.\n";
  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "\n";
  cout << setw(14) << xval[0] << setw(14) << yval[0] << "\n";
  cout << setw(14) << xval[1] << setw(14) << yval[1] << "\n";
  cout << setw(14) << xval[2] << setw(14) << yval[2] << "\n";
//
//  Get the quadrature abscissas and weights for a unit triangle.
//
  rule = 3;
  order = triangle_unit_size ( rule );

  xtab = new double[order];
  ytab = new double[order];
  weight = new double[order];

  triangle_unit_set ( rule, order, xtab, ytab, weight );

  cout << "\n";
  cout << "  Using unit triangle quadrature rule " << rule << "\n";
  cout << "  Rule order = " << order << "\n";
  cout << "\n";
  cout << "  Function Nsub  Result\n";
  cout << "\n";
//
//  Set the function.
//
  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );
//
//  Try an increasing number of subdivisions.
//
    for ( nsub = 1; nsub <= 5; nsub++ )
    {
      result = triangle_sub ( function_2d, xval, yval, nsub, order, xtab, 
        ytab,  weight );
      cout << "  " << name << setw(4) << nsub << setw(14) << result << "\n";
    }
  }
  delete [] xtab;
  delete [] ytab;
  delete [] weight;

  return;
}
//****************************************************************************80

void test42 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST42 tests TRIANGLE_UNIT_SET and TRIANGLE_UNIT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 20;
  double *weight;
  double *xtab;
  double *ytab;

  cout << "\n";
  cout << "TEST42\n";
  cout << "  TRIANGLE_UNIT_SET sets up a quadrature\n";
  cout << "    in the unit triangle,\n";
  cout << "  TRIANGLE_UNIT_SUM applies it.\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << "       " << rule;
    }
    cout << "\n";
    cout << "Function\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        triangle_unit_set ( rule, order, xtab, ytab, weight );
 
        result = triangle_unit_sum ( function_2d, order, xtab, ytab, 
          weight );

        cout << setw(14) << result;

        delete [] xtab;
        delete [] ytab;
        delete [] weight;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test425 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST425 tests TRIANGLE_UNIT_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  double quad;
  int rule;
  int rule_max = 20;
  double value;
  double *weight;
  double *xtab;
  double *ytab;

  cout << "\n";
  cout << "TEST425\n";
  cout << "  TRIANGLE_UNIT_SET sets up a quadrature\n";
  cout << "    in the unit triangle,\n";
  cout << "\n";
  cout << "  Estimate integral of X^A * Y^B.\n";

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef *( double ) ( a + i ) / ( double ) ( i );
      }

      cout << "\n";
      cout << "  A = " << a 
           << "  B = " << b << "\n";
      cout << "\n";
      cout << "  Rule       QUAD           ERROR\n";
      cout << "\n";

      for ( rule = 1; rule <= rule_max; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];
        
        triangle_unit_set ( rule, order, xtab, ytab, weight );
 
        quad = 0.0;

        for ( i = 0; i < order; i++ )
        {
          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( ytab[i], b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( xtab[i], a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( xtab[i], a ) * pow ( ytab[i], b );
          }
          quad = quad + 0.5 * weight[i] * value;
        }
        exact = 1.0;
        err = r8_abs ( exact - quad );
        cout << "  " << setw(4) << rule
             << "  " << setw(14) << quad
             << "  " << setw(11) << err << "\n";

        delete [] xtab;
        delete [] ytab;
        delete [] weight;
      }
    }
  }
  return;
}
//****************************************************************************80

void test43 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST43 tests TRIANGLE_UNIT_PRODUCT_SET and TRIANGLE_UNIT_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double *weight;
  double *xtab;
  double *ytab;

  cout << "\n";
  cout << "TEST43\n";
  cout << "  TRIANGLE_UNIT_PRODUCT_SET sets up a product quadrature\n";
  cout << "    rule in the unit triangle,\n";
  cout << "  TRIANGLE_UNIT_SUM applies it.\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    cout << "\n";
    cout << "  Rule Order: ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << setw(6) << rule;
    }
    cout << "\n";
    cout << "Function\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule = rule++ )
      {
        order = triangle_unit_product_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        triangle_unit_product_set ( rule, order, xtab, ytab, weight );

        result = triangle_unit_sum ( function_2d, order, xtab, ytab, 
          weight );

        cout << setw(14) << result;

        delete [] xtab;
        delete [] ytab;
        delete [] weight;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test44 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST44 tests TRIANGLE_UNIT_SET and TRIANGLE_SUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 20;
  double *weight;
  double *xtab;
  double xval[3] = { 1.0, 3.0, 1.0 };
  double *ytab;
  double yval[3] = { 1.0, 1.0, 4.0 };

  cout << "\n";
  cout << "TEST44\n";
  cout << "  TRIANGLE_UNIT_SET sets up quadrature\n";
  cout << "    in the unit triangle,\n";
  cout << "  TRIANGLE_SUM applies it to an arbitrary triangle.\n";
  cout << "\n";

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    cout << "\n";
    cout << "  Rule:   ";
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      cout << setw(6) << rule;
    }
    cout << "\n";

    cout << "Function\n";
    cout << "\n";

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      cout << "  " << name;

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = new double[order];
        ytab = new double[order];
        weight = new double[order];

        triangle_unit_set ( rule, order, xtab, ytab, weight );
 
        result = triangle_sum ( function_2d, xval, yval, order, xtab, ytab, 
          weight );

        cout << setw(14) << result;

        delete [] xtab;
        delete [] ytab;
        delete [] weight;
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test45 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST45 tests TORUS_1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int j2;
  int n;
  char name[8];
  int num;
  double result;
  double r1;
  double r2;

  r1 = 0.5;
  r2 = 1.0;
  n = 10;
 
  cout << "\n";
  cout << "TEST45\n";
  cout << "  TORUS_1 approximates integrals on a torus.\n";
  cout << "\n";
  cout << "  The order N will be varied.\n";
  cout << "\n";
  cout << "  Inner radius = " << r1 << "\n";
  cout << "  Outer radius = " << r2 << "\n";
  cout << "  Area = " << torus_area_3d ( r1, r2 ) << "\n";
  cout << "\n";
  cout << "  " << "  F(X)  ";
  for ( j = 1; j <= 5; j++ )
  {
    j2 = 2 * ( j - 1 );
    cout << setw(14) << i4_power ( 2, j2 );
  }
  cout << "\n";
  cout << "\n";
 
  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    cout << "  " << name;

    for ( j = 1; j <= 5; j++ )
    {
      j2 = 2 * ( j - 1 );
      n = i4_power ( 2, j2 );
      result = torus_1 ( function_3d, r1, r2, n );
      cout << "  " << setw(14) << result;
    }
    cout << "\n";
   }
  return;
}
//****************************************************************************80

void test46 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST46 tests TORUS_5S2, TORUS_6S2 and TORUS_14S.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double r1;
  double r2;

  r1 = 0.5;
  r2 = 1.0;
 
  cout << "\n";
  cout << "TEST46\n";
  cout << "  For the interior of a torus,\n";
  cout << "  TORUS_5S2,\n";
  cout << "  TORUS_6S2, and\n";
  cout << "  TORUS_5S2 approximate integrals.\n";
  cout << "\n";
  cout << "  Inner radius = " << r1 << "\n";
  cout << "  Outer radius = " << r2 << "\n";
  cout << "  Volume = " << torus_volume_3d ( r1, r2 ) << "\n";
  cout << "\n";
  cout << "    Rule:        #5S2          #6S2          #14S\n";
  cout << "    F(X)\n";
  cout << "\n";
 
  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = torus_5s2 ( function_3d, r1, r2 );
    result2 = torus_6s2 ( function_3d, r1, r2 );
    result3 = torus_14s ( function_3d, r1, r2 );

    cout << "  " << name
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2
         << "  " << setw(14) << result3 << "\n";
  }
  return;
}
//****************************************************************************80

void test47 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST47 tests TORUS_SQUARE_5C2 and TORUS_SQUARE_14C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double r1;
  double r2;

  r1 = 1.0;
  r2 = 0.125;
 
  cout << "\n";
  cout << "TEST47\n";
  cout << "  For integrals inside a torus with square cross-section:\n";
  cout << "  TORUS_SQUARE_5C2 approximates the integral;\n";
  cout << "  TORUS_SQUARE_14C approximates the integral.\n";
  cout << "\n";
  cout << "  Inner radius = " << r1 << "\n";
  cout << "  Outer radius = " << r2 << "\n";
  cout << "  Volume = " << torus_square_volume_3d ( r1, r2 ) << "\n";
  cout << "\n";
  cout << "    F(X)    5C2           14C\n";
  cout << "\n";
 
  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = torus_square_5c2 ( function_3d, r1, r2 );
    result2 = torus_square_14c ( function_3d, r1, r2 );

    cout << "  " << name
         << "  " << setw(14) << result1
         << "  " << setw(14) << result2 << "\n";
  }
  return;
}
//****************************************************************************80

void test48 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST48 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int nt;
  double *t;

  cout << "\n";
  cout << "TEST48\n";
  cout << "  For evenly spaced angles between 0 and 2*PI:\n";
  cout << "  TVEC_EVEN\n";
  cout << "  TVEC_EVEN2\n";
  cout << "  TVEC_EVEN3\n";

  nt = 4;
  t = tvec_even ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN:" );
  delete [] t;

  nt = 4;
  t = tvec_even2 ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN2:" );
  delete [] t;

  nt = 4;
  t = tvec_even3 ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN3:" );
  delete [] t;

  return;
}
//****************************************************************************80

void test49 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST49 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2 and TVEC_EVEN_BRACKET3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int nt;
  double *t;
  double theta1;
  double theta2;

  cout << "\n";
  cout << "TEST49\n";
  cout << "  For evenly spaced angles between THETA1 and THETA2:\n";
  cout << "  TVEC_EVEN_BRACKET\n";
  cout << "  TVEC_EVEN_BRACKET2.\n";
  cout << "  TVEC_EVEN_BRACKET3.\n";

  theta1 = 30.0;
  theta2 = 90.0;

  cout << "\n";
  cout << "  THETA1 = " << theta1 << "\n";
  cout << "  THETA2 = " << theta2 << "\n";

  nt = 4;
  t = tvec_even_bracket ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET" );
  delete [] t;

  nt = 5;
  t = tvec_even_bracket2 ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET2" );
  delete [] t;

  nt = 3;
  t = tvec_even_bracket3 ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET3" );
  delete [] t;

  return;
}
//****************************************************************************80

double fu18 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FU18 is the upper limit of integration for x.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double fl18 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FL18 is the lower limit of integration for x.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double value;

  value = - 1.0;

  return value;
}
//****************************************************************************80

double fu28 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FU28 computes the upper limit of integration for (x,y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double fl28 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FL28 computes the lower limit of integration for (x,y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
{
  double value;

  value = -1.0;

  return value;
}
//****************************************************************************80

double function_1d ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_1D evaluates the current 1D function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_1D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double FUNCTION_1D, the value of the function.
//
{
  double value;

  if ( function_1d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_1d_index == 2 )
  {
    value = x;
  }
  else if ( function_1d_index == 3 )
  {
    value = x * x;
  }
  else if ( function_1d_index == 4 )
  {
    value = x * x * x * x;
  }
  else if ( function_1d_index == 5 )
  {
    value = x * x * x * x;
  }
  else if ( function_1d_index == 6 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_1d_index == 7 )
  {
    value = x * x * x * x * x * x;
  }
  else if ( function_1d_index == 8 )
  {
    value = r8_abs ( x );
  }
  else if ( function_1d_index == 9 )
  {
    value = sin ( x );
  }
  else if ( function_1d_index == 10 )
  {
    value = exp ( x );
  }
  else if ( function_1d_index == 11 )
  {
    value = 1.0 / ( 1.0 + r8_abs ( x ) );
  }
  else if ( function_1d_index == 12 )
  {
    value = sqrt ( r8_abs ( x ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void function_1d_name ( char *name )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_1D_NAME returns the name of the current 1D function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_1D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char NAME[8], the name of the current 1D function.
//
{
  if ( function_1d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_1d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_1d_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_1d_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_1d_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_1d_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_1d_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_1d_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_1d_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_1d_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_1d_index == 11 ) 
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_1d_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
//****************************************************************************80

int function_1d_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_1D_NUM returns the number of 1D functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int FUNCTION_1D_NUM, the number of 1D functions.
//
{
  int value = 12;

  return value;
}
//****************************************************************************80

double function_2d ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_2D evaluates the current 2D function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_2D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, double Y, the value of the variables.
//
//    Output, double FUNCTION_2D, the value of the function.
//
{
  double value;

  if ( function_2d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_2d_index == 2 )
  {
    value = x;
  }
  else if ( function_2d_index == 3 )
  {
    value = x * x;
  }
  else if ( function_2d_index == 4 )
  {
    value = x * x * x;
  }
  else if ( function_2d_index == 5 )
  {
    value = x * x * x * x;
  }
  else if ( function_2d_index == 6 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_2d_index == 7 )
  {
    value = x * x * x * x * x * x;
  }
  else if ( function_2d_index == 8 )
  {
    value = sqrt ( x * x + y * y );
  }
  else if ( function_2d_index == 9 )
  {
    value = sin ( x );
  }
  else if ( function_2d_index == 10 )
  {
    value = exp ( x );
  }
  else if ( function_2d_index == 11 )
  {
    value = 1.0 / ( 1.0 + sqrt ( x * x + y * y ) );
  }
  else if ( function_2d_index == 12 )
  {
    value = sqrt ( sqrt ( x * x + y * y ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void function_2d_name ( char *name )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_2D_NAME returns the name of the current 2D function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_2D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char NAME[8], the name of the current 2D function.
//
{
  if ( function_2d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_2d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_2d_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_2d_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_2d_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_2d_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_2d_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_2d_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_2d_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_2d_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_2d_index == 11 ) 
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_2d_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
//****************************************************************************80

int function_2d_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_2D_NUM returns the number of 2D functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int FUNCTION_2D_NUM, the number of 2D functions.
//
{
  int value = 12;

  return value;
}
//****************************************************************************80

double function_3d ( double x, double y, double z )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_3D evaluates a function F(X,Y,Z) of 3 variables.
//
///  Discussion:
//
//    This routine assumes that the global variable FUNCTION_3D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, Z, the value of the variables.
//
//    Output, double FUNC_3D, the value of the function.
//
{
  double value;

  if ( function_3d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_3d_index == 2 )
  {
    value = x;
  }
  else if ( function_3d_index == 3 )
  {
    value = y;
  }
  else if ( function_3d_index == 4 )
  {
    value = z;
  }
  else if ( function_3d_index == 5 )
  {
    value = x * x;
  }
  else if ( function_3d_index == 6 )
  {
    value = x * y;
  }
  else if ( function_3d_index == 7 )
  {
    value = x * z;
  }
  else if ( function_3d_index == 8 )
  {
    value = y * y;
  }
  else if ( function_3d_index == 9 )
  {
    value = y * z;
  }
  else if ( function_3d_index == 10 )
  {
    value = z * z;
  }
  else if ( function_3d_index == 11 )
  {
    value = x * x * x;
  }
  else if ( function_3d_index == 12 )
  {
    value = x * y * z;
  }
  else if ( function_3d_index == 13 )
  {
    value = z * z * z;
  }
  else if ( function_3d_index == 14 )
  {
    value = x * x * x * x;
  }
  else if ( function_3d_index == 15 )
  {
    value = x * x * z * z;
  }
  else if ( function_3d_index == 16 )
  {
    value = z * z * z * z;
  }
  else if ( function_3d_index == 17 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_3d_index == 18 )
  {
    value = pow ( x, 6 );
  }
  else if ( function_3d_index == 19 )
  {
    value = sqrt ( x * x + y * y + z * z );
  }
  else if ( function_3d_index == 20 )
  {
    value = sin ( x );
  }
  else if ( function_3d_index == 21 )
  {
    value = exp ( x );
  }
  else if ( function_3d_index == 22 )
  {
    value = 1.0 / sqrt ( 1.0 + x * x + y * y + z * z );
  }
  else if ( function_3d_index == 23 )
  {
    value = sqrt ( sqrt ( x * x + y * y + z * z ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void function_3d_name ( char *name )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_3D_NAME returns the name of the current 3D function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_3D_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char NAME[8], the name of the current 3D function.
//
{
  if ( function_3d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_3d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_3d_index == 3 )
  {
    strcpy ( name, "      Y" );
  }
  else if ( function_3d_index == 4 )
  {
    strcpy ( name, "      Z" );
  }
  else if ( function_3d_index == 5 )
  {
    strcpy ( name, "    X*X" );
  }
  else if ( function_3d_index == 6 )
  {
    strcpy ( name, "    X*Y" );
  }
  else if ( function_3d_index == 7 )
  {
    strcpy ( name, "    X*Z" );
  }
  else if ( function_3d_index == 8 )
  {
    strcpy ( name, "    Y*Y" );
  }
  else if ( function_3d_index == 9 )
  {
    strcpy ( name, "    Y*Z" );
  }
  else if ( function_3d_index == 10 )
  {
    strcpy ( name, "    Z*Z" );
  }
  else if ( function_3d_index == 11 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_3d_index == 12 )
  {
    strcpy ( name, "  X*Y*Z" );
  }
  else if ( function_3d_index == 13 )
  {
    strcpy ( name, "  Z*Z*Z" );
  }
  else if ( function_3d_index == 14 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_3d_index == 15 )
  {
    strcpy ( name, "X^2 Z^2" );
  }
  else if ( function_3d_index == 16 )
  {
    strcpy ( name, "    Z^4" );
  }
  else if ( function_3d_index == 17 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_3d_index == 18 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_3d_index == 19 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_3d_index == 20 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_3d_index == 21 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_3d_index == 22 )
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_3d_index == 23 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
//****************************************************************************80

int function_3d_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_3D_NUM returns the number of 3D functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int FUNCTION_3D_NUM, the number of 3D functions.
//
{
  int value = 23;

  return value;
}
//****************************************************************************80

double function_nd ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_ND evaluates the current ND function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_ND_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Input, double X[N], the value of the variables.
//
//    Output, double FUNCTION_ND, the value of the function.
//
{
  int i;
  double temp;
  double value;

  if ( function_nd_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_nd_index == 2 )
  {
    value = x[0];
  }
  else if ( function_nd_index == 3 )
  {
    value = pow ( x[0], 2 );
  }
  else if ( function_nd_index == 4 )
  {
    value = pow ( x[0], 3 );
  }
  else if ( function_nd_index == 5 )
  {
    value = pow ( x[0], 4 );
  }
  else if ( function_nd_index == 6 )
  {
    value = pow ( x[0], 5 );
  }
  else if ( function_nd_index == 7 )
  {
    value = pow ( x[0], 6 );
  }
  else if ( function_nd_index == 8 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = sqrt ( temp );
  }
  else if ( function_nd_index == 9 )
  {
    value = sin ( x[0] );
  }
  else if ( function_nd_index == 10 )
  {
    value = exp ( x[0] );
  }
  else if ( function_nd_index == 11 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = 1.0 / ( 1.0 + sqrt ( temp ) );
  }
  else if ( function_nd_index == 12 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = sqrt ( sqrt ( temp ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void function_nd_name ( char *name )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_ND_NAME returns the name of the current ND function.
//
//  Discussion:
//
//    This routine assumes that the global variable FUNCTION_ND_INDEX has been
//    set, and is accessible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char NAME[8], the name of the current ND function.
//
{
  if ( function_nd_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_nd_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_nd_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_nd_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_nd_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_nd_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_nd_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_nd_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_nd_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_nd_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_nd_index == 11 ) 
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_nd_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
//****************************************************************************80

int function_nd_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    FUNCTION_ND_NUM returns the number of ND functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int FUNCTION_ND_NUM, the number of ND functions.
//
{
  int value = 12;

  return value;
}
//****************************************************************************80

double f_1_2d ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    F_1_2D evaluates the function 1 in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the arguments.
//
//    Output, double F_1_2D, the value of the function.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double f_x_2d ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    F_X_2D evaluates the function X in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the arguments.
//
//    Output, double F_X_2D, the value of the function.
//
{
  double value;

  value = x;

  return value;
}
//****************************************************************************80

double f_r_2d ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    F_R_2D evaluates the function sqrt ( X**2 + Y**2 ) in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the arguments.
//
//    Output, double F_R_2D, the value of the function.
//
{
  double value;

  value = sqrt ( x * x + y * y );

  return value;
}
//****************************************************************************80

double mono_000_3d ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_000_3D evaluates X^0 Y^0 Z^0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension (which is 3 here).
//
//    Input, double X[N], the evaluation point.
//
//    Output, double MONO_000_3D, the value of the monomial at X.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double mono_111_3d ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_111_3D evaluates X^1 Y^1 Z^1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension (which is 3 here).
//
//    Input, double X[N], the evaluation point.
//
//    Output, double MONO_111_3D, the value of the monomial at X.
//
{
  double value;

  value = pow ( x[0], 1 ) 
        * pow ( x[1], 1 )
        * pow ( x[2], 1 );

  return value;
}
//****************************************************************************80

double mono_202_3d ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_202_3D evaluates X^2 Y^0 Z^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension (which is 3 here).
//
//    Input, double X[N], the evaluation point.
//
//    Output, double MONO_202_3D, the value of the monomial at X.
//
{
  double value;

  value = pow ( x[0], 2 ) 
        * pow ( x[1], 0 )
        * pow ( x[2], 2 );

  return value;
}
//****************************************************************************80

double mono_422_3d ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_422_3D evaluates X^4 Y^2 Z^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension (which is 3 here).
//
//    Input, double X[N], the evaluation point.
//
//    Output, double MONO_422_3D, the value of the monomial at X.
//
{
  double value;

  value = pow ( x[0], 4 ) 
        * pow ( x[1], 2 )
        * pow ( x[2], 2 );

  return value;
}
//****************************************************************************80

double *setsim ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SETSIM defines a unit simplex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, double V[N*(N+1)], the coordinates of the N+1 vertices.
//
{
  int i;
  int j;
  double *v;

  v = new double[n*(n+1)];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n + 1; j++ )
    {
      v[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    j = i + 1;
    v[i+j*n] = 1.0;
  }

  return v;
}
