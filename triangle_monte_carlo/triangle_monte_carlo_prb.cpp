# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "triangle_monte_carlo.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
double *triangle_integrand_user ( int p_num, double p[], int f_num );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_MONTE_CARLO_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_MONTE_CARLO library.\n";
//
//  Try each sampler on the unit triangle, integrating X^2, X*Y, Y^2.
//
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Try each sampler on a general triangle, integrating a selection of functions.
//
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_MONTE_CARLO_PRB\n";
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
//    TEST01 uses TRIANGLE_SAMPLE_01 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 3;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    1.0, 0.0,
    0.0, 1.0,
    0.0, 0.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_01\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_03\n";
  cout << "  Integration region is the unit triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a \"bad\" sampler.\n";

  seed = 123456789;

  cout << "\n";
  cout << "     P_NUM      X^2             X*Y             Y^2\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_01,
      triangle_integrand_03, &seed );

    cout << "  " << setw(8) << p_num
         << "  " << setw(14) << result[0]
         << "  " << setw(14) << result[1]
         << "  " << setw(14) << result[2] << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses TRIANGLE_SAMPLE_02 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 3;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    1.0, 0.0,
    0.0, 1.0,
    0.0, 0.0 };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_02\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_03\n";
  cout << "  Integration region is the unit triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a good sampler.\n";

  seed = 123456789;

  cout << "\n";
  cout << "     P_NUM      X^2             X*Y             Y^2\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_02,
      triangle_integrand_03, &seed );

    cout << "  " << setw(8) << p_num
         << "  " << setw(14) << result[0]
         << "  " << setw(14) << result[1]
         << "  " << setw(14) << result[2] << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses TRIANGLE_SAMPLE_03 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 3;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    1.0, 0.0,
    0.0, 1.0,
    0.0, 0.0 };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_03\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_03\n";
  cout << "  Integration region is the unit triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a good sampler.\n";

  seed = 123456789;

  cout << "\n";
  cout << "     P_NUM      X^2             X*Y             Y^2\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_03,
      triangle_integrand_03, &seed );

    cout << "  " << setw(8) << p_num
         << "  " << setw(14) << result[0]
         << "  " << setw(14) << result[1]
         << "  " << setw(14) << result[2] << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 uses TRIANGLE_SAMPLE_04 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 3;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    1.0, 0.0,
    0.0, 1.0,
    0.0, 0.0 };

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_04\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_03\n";
  cout << "  Integration region is the unit triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a good sampler.\n";

  seed = 123456789;

  cout << "\n";
  cout << "     P_NUM      X^2             X*Y             Y^2\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_04,
      triangle_integrand_03, &seed );

    cout << "  " << setw(8) << p_num
         << "  " << setw(14) << result[0]
         << "  " << setw(14) << result[1]
         << "  " << setw(14) << result[2] << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 uses TRIANGLE_SAMPLE_01 on a general triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 8;
  int i;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    4.0, 1.0,
    8.0, 3.0,
    0.0, 9.0 };

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_01\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_USER\n";
  cout << "  Integration region is over a general triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a \"bad\" sampler.\n";

  seed = 123456789;

  r8mat_transpose_print ( 2, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "     P_NUM\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_01,
      triangle_integrand_user, &seed );

    cout << "  " << setw(8) << p_num;
    for ( i = 0; i < f_num; i++ )
    {
      cout << "  " << setw(12) << result[i];
    }
    cout << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 uses TRIANGLE_SAMPLE_02 on a general triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 8;
  int i;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    4.0, 1.0,
    8.0, 3.0,
    0.0, 9.0 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_02\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_USER\n";
  cout << "  Integration region is over a general triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a \"good\" sampler.\n";

  seed = 123456789;

  r8mat_transpose_print ( 2, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "     P_NUM\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_02,
      triangle_integrand_user, &seed );

    cout << "  " << setw(8) << p_num;
    for ( i = 0; i < f_num; i++ )
    {
      cout << "  " << setw(12) << result[i];
    }
    cout << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 uses TRIANGLE_SAMPLE_03 on a general triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 8;
  int i;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    4.0, 1.0,
    8.0, 3.0,
    0.0, 9.0 };

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_03\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_USER\n";
  cout << "  Integration region is over a general triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a \"good\" sampler.\n";

  seed = 123456789;

  r8mat_transpose_print ( 2, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "     P_NUM\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_03,
      triangle_integrand_user, &seed );

    cout << "  " << setw(8) << p_num;
    for ( i = 0; i < f_num; i++ )
    {
      cout << "  " << setw(12) << result[i];
    }
    cout << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 uses TRIANGLE_SAMPLE_04 on a general triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int f_num = 8;
  int i;
  int p_num;
  double *result;
  int seed;
  double t[2*3] = {
    4.0, 1.0,
    8.0, 3.0,
    0.0, 9.0 };

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Sample using TRIANGLE_UNIT_SAMPLE_04\n";
  cout << "  Integrate TRIANGLE_UNIT_INTEGRAND_USER\n";
  cout << "  Integration region is over a general triangle.\n";
  cout << "\n";
  cout << "  Use an increasing number of points P_NUM.\n";
  cout << "  Note that the sample routine is a \"good\" sampler.\n";

  seed = 123456789;

  r8mat_transpose_print ( 2, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "     P_NUM\n";
  cout << "\n";

  p_num = 1;

  while ( p_num <= 65536 )
  {
    result = triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_04,
      triangle_integrand_user, &seed );

    cout << "  " << setw(8) << p_num;
    for ( i = 0; i < f_num; i++ )
    {
      cout << "  " << setw(12) << result[i];
    }
    cout << "\n";

    p_num = 2 * p_num;

    delete [] result;
  }

  return;
}
//****************************************************************************80

double *triangle_integrand_user ( int p_num, double p[], int f_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INTEGRAND_USER evaluates 8 integrand functions defined by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P(2,P_NUM), the evaluation points.
//
//    Input, int F_NUM, the number of integrands.
//
//    Output, double FP(F_NUM,P_NUM), the integrand values.
//
{
  double *fp;
  int j;

  fp = new double[f_num*p_num];

  for ( j = 0; j < p_num; j++ )
  {
    fp[0+j*f_num] = 1.0;
    fp[1+j*f_num] = p[0+j*2];
    fp[2+j*f_num] = p[1+j*2];
    fp[3+j*f_num] = p[0+j*2] * p[0+j*2];
    fp[4+j*f_num] = p[0+j*2] * p[1+j*2];
    fp[5+j*f_num] = p[1+j*2] * p[1+j*2];
    fp[6+j*f_num] = p[0+j*2] * p[0+j*2] * p[1+j*2];
    fp[7+j*f_num] = p[0+j*2] * p[0+j*2] * p[1+j*2] * p[1+j*2];
  }
  return fp;
}

