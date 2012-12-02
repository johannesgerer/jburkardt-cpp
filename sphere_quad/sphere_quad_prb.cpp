# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "sphere_quad.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void polyterm_exponent ( string s, int e[3] );
void polyterm_value_3d ( int n, double x[], double f[] );

//
//  Global data.
//
int e_save[3];

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_QUAD_PRB.
//
//  Discussion:
//
//    SPHERE_QUAD_PRB tests SPHERE_QUAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "SPHERE_QUAD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_QUAD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_QUAD_PRB\n";
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
//    TEST01 tests SPHERE01_QUAD_LL*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  double exact;
  double h;
  int h_test;
  int i;
  int n_llc;
  int n_llm;
  int n_llv;
  int n_mc;
  double result_llc;
  double result_llm;
  double result_llv;
  double result_mc;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Approximate the integral of a function on the unit sphere.\n";
  cout << "\n";
  cout << "  SPHERE01_QUAD_MC uses a Monte Carlo method.\n";
  cout << "  SPHERE01_QUAD_LLC uses centroids of spherical triangles.\n";
  cout << "  SPHERE01_QUAD_LLM uses midsides of spherical triangles.\n";
  cout << "  SPHERE01_QUAD_LLV uses vertices of spherical triangles.\n";
  cout << "\n";
  cout << "  H              QUAD_MC       QUAD_LLC      QUAD_LLM      QUAD_LLV         EXACT\n";

  for ( i = 0; i <= 17; i++ )
  {
    if ( i == 0 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    if ( i == 0 )
    {
      cout << "\n";
      cout << "Point counts per method:\n";
    }
    else
    {
      polyterm_exponent ( "PRINT", e );
    }

    for ( h_test = 1; h_test <= 3; h_test++ )
    {
      if ( h_test == 1 )
      {
        h = 1.0;
      }
      else if (  h_test == 2 )
      {
        h = 0.1;
      }
      else if (  h_test == 3 )
      {
        h = 0.01;
      }

      n_mc = sphere01_quad_mc_size ( h );

      result_mc = sphere01_quad_mc ( polyterm_value_3d, h, &seed, n_mc );

      result_llc = sphere01_quad_llc ( polyterm_value_3d, h, &n_llc );

      result_llm = sphere01_quad_llm ( polyterm_value_3d, h, &n_llm );

      result_llv = sphere01_quad_llv ( polyterm_value_3d, h, &n_llv );

      exact = sphere01_monomial_integral ( e );

      if ( i == 0 )
      {
        cout << "  " << setw(12) << h
             << "  " << setw(12) << n_mc
             << "  " << setw(12) << n_llc
             << "  " << setw(12) << n_llm
             << "  " << setw(12) << n_llv << "\n";
      }
      else
      {
        cout << "  " << setw(12) << h
             << "  " << setw(12) << result_mc
             << "  " << setw(12) << result_llc
             << "  " << setw(12) << result_llm
             << "  " << setw(12) << result_llv
             << "  " << setw(12) << exact << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SPHERE01_QUAD_ICOS1C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  double error;
  double exact;
  int factor;
  int factor_log;
  int i;
  int n;
  double result;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Approximate the integral of a function on the unit sphere.\n";
  cout << "  SPHERE01_QUAD_ICOS1C uses centroids of spherical triangles.\n";
  cout << "\n";
  cout << "FACTOR         N        QUAD          EXACT         ERROR\n";

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = 1;
    for ( factor_log = 0; factor_log <= 5; factor_log++ )
    {
      result = sphere01_quad_icos1c ( factor, polyterm_value_3d, &n );

      exact = sphere01_monomial_integral ( e );

      error = r8_abs ( exact - result );

      cout << "  " << setw(4) << factor
           << "  " << setw(8) << n
           << "  " << setw(14) << result
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";

      factor = factor * 2;
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
//    TEST03 tests SPHERE01_QUAD_ICOS1M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  double error;
  double exact;
  int factor;
  int factor_log;
  int i;
  int n;
  double result;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Approximate the integral of a function on the unit sphere.\n";
  cout << "  SPHERE01_QUAD_ICOS1M uses midpoints of sides of spherical triangles.\n";
  cout << "\n";
  cout << "FACTOR         N        QUAD          EXACT         ERROR\n";

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = 1;
    for ( factor_log = 0; factor_log <= 5; factor_log++ )
    {
      result = sphere01_quad_icos1m ( factor, polyterm_value_3d, &n );

      exact = sphere01_monomial_integral ( e );

      error = r8_abs ( exact - result );

      cout << "  " << setw(4) << factor
           << "  " << setw(8) << n
           << "  " << setw(14) << result
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";

      factor = factor * 2;
    }
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests SPHERE01_QUAD_ICOS1V.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  double error;
  double exact;
  int factor;
  int factor_log;
  int i;
  int n;
  double result;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Approximate the integral of a function on the unit sphere.\n";
  cout << "  SPHERE01_QUAD_ICOS1V uses vertices of spherical triangles.\n";
  cout << "\n";
  cout << "FACTOR         N        QUAD          EXACT         ERROR\n";

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = 1;
    for ( factor_log = 0; factor_log <= 5; factor_log++ )
    {
      result = sphere01_quad_icos1v ( factor, polyterm_value_3d, &n );

      exact = sphere01_monomial_integral ( e );

      error = r8_abs ( exact - result );

      cout << "  " << setw(4) << factor
           << "  " << setw(8) << n
           << "  " << setw(14) << result
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";

      factor = factor * 2;
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
//    TEST05 tests SPHERE01_QUAD_ICOS2V.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  double error;
  double exact;
  int factor;
  int factor_log;
  int i;
  int n;
  double result;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Approximate the integral of a function on the unit sphere.\n";
  cout << "  SPHERE01_QUAD_ICOS2V uses vertices of spherical triangles.\n";
  cout << "\n";
  cout << "FACTOR         N        QUAD          EXACT         ERROR\n";

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = 1;
    for ( factor_log = 0; factor_log <= 5; factor_log++ )
    {
      result = sphere01_quad_icos2v ( factor, polyterm_value_3d, &n );

      exact = sphere01_monomial_integral ( e );

      error = r8_abs ( exact - result );

      cout << "  " << setw(4) << factor
           << "  " << setw(8) << n
           << "  " << setw(14) << result
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";

      factor = factor * 2;
    }
  }
  return;
}
//****************************************************************************80

void polyterm_exponent ( string action, int e[3] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    'GET' asks the routine to return the current values in E.
//    'SET' asks the routine to set the current values to E.
//
//    Input/output, int E[3], storage used to set or get values.
//
{
  int i;

  if ( action[0] == 'G' )
  {
    for ( i = 0; i < 3; i++ )
    {
      e[i] = e_save[i];
    }
  }
  else if ( action[0] == 'P' )
  {
    cout << "\n";

    if ( e_save[0] == 0 && e_save[1] == 0 && e_save[2] == 0 )
    {
      cout << "P(X,Y,Z) = 1\n";
    }
    else
    {
      cout << "P(X,Y,Z) = ";

      if ( e_save[0] == 0 )
      {
      }
      else if ( e_save[0] == 1 )
      {
        cout << " X";
      }
      else
      {
        cout << " X^" << e_save[0];
      }
      if ( e_save[1] == 0 )
      {
      }
      else if ( e_save[1] == 1 )
      {
        cout << " Y";
      }
      else
      {
        cout << " Y^" << e_save[1];
      }
      if ( e_save[2] == 0 )
      {
      }
      else if ( e_save[2] == 1 )
      {
        cout << " Z";
      }
      else
      {
        cout << " Z^" << e_save[2];
      }
      cout << "\n";
    }
  } 
  else if ( action[0] == 'S' )
  {
    for ( i = 0; i < 3; i++ )
    {
      e_save[i] = e[i];
    }
  }

  return;
}
//****************************************************************************80

void polyterm_value_3d ( int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
//
//  Discussion:
//
//    The polynomial term has the form:
//
//      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
//
//    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[3*N], the points where the polynomial term 
//    is to be evaluated.
//
//    Output, double F[N], the value of the polynomial term.
//
{
  int e[3];
  int i;
  int j;

  polyterm_exponent ( "GET", e );

  for ( j = 0; j < n; j++ )
  {
    f[j] = 1.0;
  }
  for ( i = 0; i < 3; i++ )
  {
    if ( e[i] != 0 )
    {
      for ( j = 0; j < n; j++ )
      {
        f[j] = f[j] * pow ( x[i+j*3], e[i] );
      }
    }
  }
  
  return;
}
