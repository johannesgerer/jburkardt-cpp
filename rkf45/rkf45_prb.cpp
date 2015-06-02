# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "rkf45.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void r4_f1 ( float t, float y[], float yp[] );
float r4_y1x ( float t );
void r4_f2 ( float t, float y[], float yp[] );
void r8_f1 ( double t, double y[], double yp[] );
double r8_y1x ( double t );
void r8_f2 ( double t, double y[], double yp[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RKF45_PRB.
//
//  Discussion:
//
//    RKF45_PRB tests the RKF45 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "RKF45_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the RKF45 library.\n";

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
  cout << "RKF45_PRB\n";
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
//    TEST01 solves a scalar ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  float abserr;
  int flag;
  int i_step;
  int n_step;
  int neqn;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float *y;
  float *yp;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Solve a scalar equation using R4_RKF:\n";
  cout << "\n";
  cout << "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n";
  cout << "\n";

  neqn = 1;

  y = new float[neqn];
  yp = new float[neqn];

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r4_f1 ( t, y, yp );

  cout << "\n";
  cout << "FLAG             T          Y         Y'          Y_Exact         Error\n";
  cout << "\n";

  cout << setw(4)  << flag               << "  "
       << setw(12) << t                  << "  "
       << setw(12) << y[0]               << "  " 
       << setw(12) << yp[0]              << "  " 
       << setw(12) << r4_y1x ( t )        << "  "
       << setw(12) << y[0] - r4_y1x ( t ) << "\n";;

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start  
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start  
            + ( float ) (          i_step ) * t_stop )  
            / ( float ) ( n_step          );

    flag = r4_rkf45 ( r4_f1, neqn, y, yp, &t, t_out, &relerr, abserr, flag );

    cout << setw(4)  << flag               << "  "
         << setw(12) << t                  << "  "
         << setw(12) << y[0]               << "  "
         << setw(12) << yp[0]              << "  "  
         << setw(12) << r4_y1x ( t )        << "  "
         << setw(12) << y[0] - r4_y1x ( t ) << "\n";
  }

  delete [] y;
  delete [] yp;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 solves a vector ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  float abserr;
  int flag;
  int i_step;
  int n_step;
  int neqn;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float *y;
  float *yp;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Solve a vector equation using R4_RKF:\n";
  cout << "\n";
  cout << "  Y'(1) =  Y(2)\n";
  cout << "  Y'(2) = -Y(1)\n";
  cout << "\n";
  cout << "\n";
  cout << "  This system is equivalent to the following\n";
  cout << "  second order system:\n";
  cout << "\n";
  cout << "  Z\" = - Z.\n";

  neqn = 2;

  y = new float[neqn];
  yp = new float[neqn];

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 2.0 * 3.14159265;

  n_step = 12;

  t = 0.0;
  t_out = 0.0;

  y[0] = 1.0;
  y[1] = 0.0;

  cout << "\n";
  cout << "FLAG             T          Y(1)       Y(2)\n";
  cout << "\n";

  cout << setw(4)  << flag  << "  "
       << setw(12) << t     << "  "
       << setw(12) << y[0]  << "  "
       << setw(12) << y[1]  << "\n";

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start 
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start 
            + ( float ) (	   i_step ) * t_stop ) 
            / ( float ) ( n_step );

    flag = r4_rkf45 ( r4_f2, neqn, y, yp, &t, t_out, &relerr, abserr, flag );

    cout << setw(4)  << flag  << "  "
         << setw(12) << t     << "  "
         << setw(12) << y[0]  << "  " 
         << setw(12) << y[1]  << "\n";
  }

  delete [] y;
  delete [] yp;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 solves a scalar ODE using single step mode.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NEQN 1

  float abserr;
  int flag;
  int i_step;
  int n_step;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float y[NEQN];
  float yp[NEQN];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Solve a scalar equation using R4_RKF:\n";
  cout << "\n";
  cout << "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n";
  cout << "\n";
  cout << "  This routine uses the SINGLE STEP mode.\n";
  cout << "\n";

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = -1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r4_f1 ( t, y, yp );

  cout << "\n";
  cout << "FLAG             T          Y         Y'        Y_Exact         Error\n";
  cout << "\n";

  cout << setw(4)  << flag               << "  "
       << setw(12) << t                  << "  "
       << setw(12) << y[0]               << "  " 
       << setw(12) << yp[0]              << "  " 
       << setw(12) << r4_y1x ( t )        << "  "
       << setw(12) << y[0] - r4_y1x ( t ) << "\n";;

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start  
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start  
            + ( float ) (          i_step ) * t_stop )  
            / ( float ) ( n_step          );

    while ( flag < 0 )
    {
      flag = r4_rkf45 ( r4_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

      cout << setw(4)  << flag               << "  "
           << setw(12) << t                  << "  "
           << setw(12) << y[0]               << "  " 
           << setw(12) << yp[0]              << "  " 
           << setw(12) << r4_y1x ( t )        << "  "
           << setw(12) << y[0] - r4_y1x ( t ) << "\n";
    }
    flag = -2;
  }

  return;
# undef NEQN
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 solves a scalar ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NEQN 1

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Solve a scalar equation using R8_RKF:\n";
  cout << "\n";
  cout << "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n";
  cout << "\n";

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r8_f1 ( t, y, yp );

  cout << "\n";
  cout << "FLAG             T          Y         Y'          Y_Exact         Error\n";
  cout << "\n";

  cout << setw(4)  << flag               << "  "
       << setw(12) << t                  << "  "
       << setw(12) << y[0]               << "  " 
       << setw(12) << yp[0]              << "  " 
       << setw(12) << r8_y1x ( t )        << "  "
       << setw(12) << y[0] - r8_y1x ( t ) << "\n";;

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start  
            + ( double ) (          i_step ) * t_stop )  
            / ( double ) ( n_step          );

    flag = r8_rkf45 ( r8_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    cout << setw(4)  << flag               << "  "
         << setw(12) << t                  << "  "
         << setw(12) << y[0]               << "  "
         << setw(12) << yp[0]              << "  "  
         << setw(12) << r8_y1x ( t )        << "  "
         << setw(12) << y[0] - r8_y1x ( t ) << "\n";
  }

  return;
# undef NEQN
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 solves a vector ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NEQN 2

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Solve a vector equation using R8_RKF:\n";
  cout << "\n";
  cout << "  Y'(1) =  Y(2)\n";
  cout << "  Y'(2) = -Y(1)\n";
  cout << "\n";

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 2.0 * 3.14159265;

  n_step = 12;

  t = 0.0;
  t_out = 0.0;

  y[0] = 1.0;
  y[1] = 0.0;
  r8_f2 ( t, y,  yp );

  cout << "\n";
  cout << "FLAG             T          Y(1)       Y(2)\n";
  cout << "\n";

  cout << setw(4)  << flag  << "  "
       << setw(12) << t     << "  "
       << setw(12) << y[0]  << "  "
       << setw(12) << y[1]  << "\n";

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start 
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start 
            + ( double ) (	   i_step ) * t_stop ) 
            / ( double ) ( n_step );

    flag = r8_rkf45 ( r8_f2, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    cout << setw(4)  << flag  << "  "
         << setw(12) << t     << "  "
         << setw(12) << y[0]  << "  " 
         << setw(12) << y[1]  << "\n";
  }

  return;
# undef NEQN
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 solves a scalar ODE using single step mode.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NEQN 1

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Solve a scalar equation using R8_RKF:\n";
  cout << "\n";
  cout << "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n";
  cout << "\n";
  cout << "  This routine uses the SINGLE STEP mode.\n";
  cout << "\n";

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = -1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r8_f1 ( t, y, yp );

  cout << "\n";
  cout << "FLAG             T          Y         Y'        Y_Exact         Error\n";
  cout << "\n";

  cout << setw(4)  << flag               << "  "
       << setw(12) << t                  << "  "
       << setw(12) << y[0]               << "  " 
       << setw(12) << yp[0]              << "  " 
       << setw(12) << r8_y1x ( t )        << "  "
       << setw(12) << y[0] - r8_y1x ( t ) << "\n";;

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start  
            + ( double ) (          i_step ) * t_stop )  
            / ( double ) ( n_step          );

    while ( flag < 0 )
    {
      flag = r8_rkf45 ( r8_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

      cout << setw(4)  << flag               << "  "
           << setw(12) << t                  << "  "
           << setw(12) << y[0]               << "  " 
           << setw(12) << yp[0]              << "  " 
           << setw(12) << r8_y1x ( t )        << "  "
           << setw(12) << y[0] - r8_y1x ( t ) << "\n";
    }
    flag = -2;
  }

  return;
# undef NEQN
}
//****************************************************************************80

void r4_f1 ( float t, float y[], float yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4_F1 evaluates the derivative for the ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float T, the value of the independent variable.
//
//    Input, float Y[NEQN], the value of the dependent variable.
//
//    Output, float YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
//
{
  yp[0] = 0.25 * y[0] * ( 1.0 - y[0] / 20.0 );

  return;
}
//****************************************************************************80

float r4_y1x ( float t )

//****************************************************************************80
//
//  Purpose:
//
//    R4_Y1X evaluates the exact solution of the ODE.
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float T, the value of the independent variable.
//
//    Output, float R4_Y1X, the exact solution.
//
{
  float value;

  value = 20.0 / ( 1.0 + 19.0 * exp ( -0.25 * t ) );

  return value;
}
//****************************************************************************80

void r4_f2 ( float t, float y[], float yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4_F2 evaluates the derivative for the ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float T, the value of the independent variable.
//
//    Input, float Y(NEQN), the value of the dependent variable.
//
//    Output float YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
//
{
  yp[0] =  y[1];
  yp[1] = -y[0];

  return;
}
//****************************************************************************80

void r8_f1 ( double t, double y[], double yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8_F1 evaluates the derivative for the ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the value of the independent variable.
//
//    Input, double Y[NEQN], the value of the dependent variable.
//
//    Output, double YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
//
{
  yp[0] = 0.25 * y[0] * ( 1.0 - y[0] / 20.0 );

  return;
}
//****************************************************************************80

double r8_y1x ( double t )

//****************************************************************************80
//
//  Purpose:
//
//    R8_Y1X evaluates the exact solution of the ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the value of the independent variable.
//
//    Output, double R8_Y1X, the exact solution.
//
{
  double value;

  value = 20.0 / ( 1.0 + 19.0 * exp ( -0.25 * t ) );

  return value;
}
//****************************************************************************80

void r8_f2 ( double t, double y[], double yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8_F2 evaluates the derivative for the ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the value of the independent variable.
//
//    Input, double Y(NEQN), the value of the dependent variable.
//
//    Output double YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
//
{
  yp[0] =  y[1];
  yp[1] = -y[0];

  return;
}
